##  Chronic Insomnia Disorder (CID) — QSP Simulator (Shiny)
##  ============================================================================
##  Front-end for ins_mrgsolve_model.R (63-ODE chronic insomnia QSP model).
##
##  Run with:
##      shiny::runApp("ins_shiny_app.R")
##  from inside the chronic-insomnia-disorder/ directory, or
##      shiny::runApp("chronic-insomnia-disorder/ins_shiny_app.R")
##  from the repository root.
##
##  Requires: shiny, mrgsolve, dplyr, tidyr, ggplot2, DT
##
##  Design intent
##  -------------
##  Insomnia is the one sleep disorder where the patient's own behaviour is
##  part of the pathophysiology, so the app is laid out to keep the SCHEDULE
##  (time in bed, rise time, light) visible next to the pharmacology on every
##  tab. Two things the app is built to make obvious:
##
##    1. A single night and twelve weeks are different questions. The
##       "Two-process & switch" tab shows why a patient falls asleep tonight;
##       the "Perpetuating loop" and "Endpoints" tabs show why they will still
##       be an insomniac in three months. Hypnotics change the first plot and
##       not the second; CBT-I does the reverse.
##
##    2. Every arm is compared against its OWN untreated control, run from the
##       same lead-in state, because the phenotypes have different baselines.
##
##  Tabs
##    1  Patient profile        Patient & phenotype — arousal decomposition
##    2  Two-process & switch     Process S x C, wake drive, one-night sleep switch
##    3  Circadian & melatonin    Circadian phase, light, melatonin, temperature
##    4  Drug PK / target occupancy   Concentrations and receptor occupancy
##    5  Nightly sleep metrics       SOL / WASO / TST / SE trajectories
##    6  Sleep architecture            N3, REM, awakenings, hypnogram surrogate
##    7  Perpetuating loop (3-P)      COND / SEFF / DBAS / TIB and the CBT-I response
##    8  Clinical endpoints      ISI, subjective sleep, remission, safety
##    9  Scenario comparison        Several arms side by side with a summary table
##  ============================================================================

## NOTE: mrgsolve exports its own req(), which masks shiny::req() because it
## is attached later. Shiny's version is called with its namespace below.
library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

source("ins_mrgsolve_model.R")

## ---------------------------------------------------------------- helpers --
PHENOS <- c(
  "Good sleeper"                = "good_sleeper",
  "Chronic insomnia disorder"           = "chronic_insomnia",
  "Severe insomnia disorder"                     = "severe_insomnia",
  "Sleep-onset predominant type"              = "sleep_onset",
  "Sleep-maintenance predominant type"    = "sleep_maintenance",
  "Elderly"                             = "elderly",
  "Comorbid depression"            = "comorbid_depression",
  "Delayed sleep phase"               = "delayed_phase",
  "Shift worker"                  = "shift_worker"
)

DRUGS <- list(
  "None"                 = list(cmt = NA,     dose = 0),
  "Zolpidem"             = list(cmt = "ZOLD", dose = 10),
  "Eszopiclone"    = list(cmt = "ESZD", dose = 3),
  "Suvorexant"       = list(cmt = "SUVD", dose = 20),
  "Lemborexant"      = list(cmt = "LEMD", dose = 10),
  "Daridorexant"   = list(cmt = "DARD", dose = 50),
  "Ramelteon"          = list(cmt = "RAMD", dose = 8),
  "Melatonin"          = list(cmt = "MELD", dose = 2),
  "Low-dose doxepin"       = list(cmt = "DOXD", dose = 6),
  "Trazodone"          = list(cmt = "TRZD", dose = 50)
)

theme_ins <- function() {
  theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold"),
          legend.position = "bottom",
          plot.title = element_text(face = "bold"))
}

night_shade <- function(days, waket, tib) {
  bed <- (waket + 24 - tib) %% 24
  data.frame(
    xmin = sapply(0:(days - 1), function(d) d*24 + bed - ifelse(bed > waket, 0, 24)),
    xmax = sapply(0:(days - 1), function(d) d*24 + waket)
  )
}

## -------------------------------------------------------------------- UI ---
ui <- fluidPage(
  titlePanel("Chronic Insomnia Disorder QSP Simulator"),
  tags$p(style = "color:#666;margin-top:-8px",
         "Two-process (S x C) core · VLPO/arousal flip-flop · 24-h hyperarousal · ",
         "Spielman 3-P perpetuating loop · BzRA / DORA / melatonergic / H1 pharmacology"),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("① Patient"),
      selectInput("pheno", "Phenotype", choices = PHENOS,
                  selected = "chronic_insomnia"),
      sliderInput("A0", "Trait hyperarousal A0", 0.20, 1.20, 0.74, step = 0.02),
      sliderInput("TIB0", "Habitual time in bed TIB (h)", 5.0, 10.0, 8.8, step = 0.1),
      sliderInput("WAKET", "Rise time (clock h)", 4, 11, 7, step = 0.25),
      sliderInput("TSTNEED", "Individual sleep need (min)", 360, 540, 450, step = 10),

      h4("② Comorbidity"),
      sliderInput("PAIN", "Chronic pain burden (0-1)", 0, 1, 0, step = 0.05),
      sliderInput("VMS",  "Vasomotor symptom burden (0-1)", 0, 1, 0, step = 0.05),
      sliderInput("STR0", "Residual precipitating stress (0-1)", 0, 1, 0, step = 0.05),

      h4("③ Drug"),
      selectInput("drug", "Nighttime hypnotic", choices = names(DRUGS),
                  selected = "None"),
      numericInput("dose", "Dose per administration (mg)", 10, min = 0, max = 200, step = 0.5),
      sliderInput("txweeks", "Treatment duration (weeks)", 0, 12, 12, step = 1),
      sliderInput("bedhour", "Administration time (clock h)", 19, 24, 22, step = 0.25),
      sliderInput("CYP3A", "CYP3A4 activity multiplier", 0.2, 2.0, 1.0, step = 0.05),

      h4("④ Non-pharmacological"),
      checkboxInput("cbti", "Administer CBT-I", FALSE),
      sliderInput("cbtiday", "CBT-I start day", 1, 42, 1, step = 1),
      checkboxInput("blt", "Morning bright light therapy", FALSE),
      sliderInput("caf", "Caffeine (mg, at 16:00)", 0, 400, 0, step = 25),
      sliderInput("etoh", "Alcohol before bedtime (g)", 0, 60, 0, step = 5),
      sliderInput("luxeve", "Evening light level (lux)", 20, 600, 120, step = 10),

      hr(),
      sliderInput("weeks", "Observation period (weeks)", 2, 12, 12, step = 1),
      actionButton("run", "Run simulation", class = "btn-primary btn-block"),
      tags$small(style = "color:#888",
        "Each run brings the patient to a chronic equilibrium state via a 150-day untreated lead-in, then ",
        "displays only the results from the randomisation point onward. The control arm is an untreated ",
        "arm starting from the same lead-in.")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",

        tabPanel("1 Patient profile",
          fluidRow(
            column(6, plotOutput("p_arousal", height = 320)),
            column(6, plotOutput("p_hpa", height = 320))),
          fluidRow(column(12, h4("Arousal decomposition"),
                          DTOutput("t_profile")))),

        tabPanel("2 Two-process & switch",
          fluidRow(column(12, plotOutput("p_night", height = 380))),
          fluidRow(
            column(6, plotOutput("p_switch", height = 300)),
            column(6, plotOutput("p_drivecomp", height = 300))),
          tags$p(style = "color:#666",
            "Above: a representative night. Sleep probability SLP is sigmoid((S·wS + adenosine − C(t) − arousal − ",
            "ascending arousal drive + drug effect − θ)/k), and because arousal-system tone is then suppressed by ",
            "SLP, this switch is bistable. This is why sleep onset and waking do not transition smoothly but ",
            "'flip over'.")),

        tabPanel("3 Circadian & melatonin",
          fluidRow(
            column(6, plotOutput("p_phi", height = 300)),
            column(6, plotOutput("p_mel", height = 300))),
          fluidRow(
            column(6, plotOutput("p_temp", height = 280)),
            column(6, plotOutput("p_light", height = 280)))),

        tabPanel("4 Drug PK / target occupancy",
          fluidRow(column(12, plotOutput("p_pk", height = 330))),
          fluidRow(column(12, plotOutput("p_occ", height = 330)))),

        tabPanel("5 Nightly sleep metrics",
          fluidRow(column(12, plotOutput("p_metrics", height = 420))),
          fluidRow(column(12, DTOutput("t_weekly")))),

        tabPanel("6 Sleep architecture",
          fluidRow(
            column(6, plotOutput("p_stage", height = 320)),
            column(6, plotOutput("p_stagepct", height = 320))),
          fluidRow(column(12, plotOutput("p_naw", height = 260)))),

        tabPanel("7 Perpetuating loop (3-P)",
          fluidRow(column(12, plotOutput("p_loop", height = 400))),
          fluidRow(column(12, plotOutput("p_tib", height = 260))),
          tags$p(style = "color:#666",
            "This tab is the whole point of this model. Hypnotics barely touch these four state variables, ",
            "so the effect disappears when discontinued, whereas CBT-I directly targets these four variables, ",
            "so the effect remains even after treatment ends.")),

        tabPanel("8 Clinical endpoints",
          fluidRow(
            column(6, plotOutput("p_isi", height = 320)),
            column(6, plotOutput("p_subj", height = 320))),
          fluidRow(
            column(6, plotOutput("p_safety", height = 300)),
            column(6, plotOutput("p_conseq", height = 300))),
          fluidRow(column(12, DTOutput("t_endpoint")))),

        tabPanel("9 Scenario comparison",
          fluidRow(column(12,
            checkboxGroupInput("scen", "Scenarios to compare",
              choices = names(ins_scenarios()),
              selected = names(ins_scenarios())[c(1, 2, 7, 13)], inline = FALSE),
            actionButton("runscen", "Run selected scenarios", class = "btn-success"))),
          fluidRow(column(12, plotOutput("p_scen", height = 420))),
          fluidRow(column(12, DTOutput("t_scen"))))
      )
    )
  )
)

## ----------------------------------------------------------------- SERVER --
server <- function(input, output, session) {

  observeEvent(input$pheno, {
    p <- ins_phenotype(input$pheno)
    if (!is.null(p$A0))   updateSliderInput(session, "A0",   value = p$A0)
    if (!is.null(p$TIB0)) updateSliderInput(session, "TIB0", value = p$TIB0)
    if (!is.null(p$WAKET))updateSliderInput(session, "WAKET",value = p$WAKET)
  })

  observeEvent(input$drug, {
    d <- DRUGS[[input$drug]]
    updateNumericInput(session, "dose", value = d$dose)
  })

  build_events <- reactive({
    ev <- NULL
    d <- DRUGS[[input$drug]]
    nd <- input$txweeks*7
    if (!is.na(d$cmt) && input$dose > 0 && nd > 0)
      ev <- ev_nightly(d$cmt, input$dose, days = nd, bedhour = input$bedhour)
    if (input$caf > 0)
      ev <- ev %+% ev_daily("CAFD", input$caf, days = input$weeks*7, hour = 16)
    if (input$etoh > 0)
      ev <- ev %+% ev_nightly("ETHD", input$etoh, days = input$weeks*7,
                              bedhour = max(19, input$bedhour - 0.5))
    ev
  })

  user_param <- reactive({
    list(A0 = input$A0, TIB0 = input$TIB0, WAKET = input$WAKET,
         TSTNEED = input$TSTNEED, PAIN = input$PAIN, VMS = input$VMS,
         LUXEVE = input$luxeve)
  })

  treat_param <- reactive({
    p <- list(CYP3A = input$CYP3A)
    if (input$blt) p <- c(p, list(BLTON = 1, BLTLUX = 10000,
                                  BLTT0 = input$WAKET, BLTDUR = 0.75))
    p
  })

  sim <- eventReactive(input$run, {
    withProgress(message = "Simulating (150-day lead-in + observation period)...", value = 0.3, {
      act <- ins_simulate(phenotype = input$pheno, weeks = input$weeks,
                          param = user_param(), treat = treat_param(),
                          events = build_events(),
                          cbti_start_day = if (input$cbti) input$cbtiday else NA_real_)
      incProgress(0.5, detail = "Control arm")
      ctl <- ins_simulate(phenotype = input$pheno, weeks = input$weeks,
                          param = user_param(), treat = list(),
                          events = NULL, cbti_start_day = NA_real_)
      list(act = act, ctl = ctl,
           nact = ins_nightly(act), nctl = ins_nightly(ctl))
    })
  }, ignoreNULL = FALSE)

  ## a representative night, taken from the middle of the observation window
  night <- reactive({
    s <- sim()$act
    d0 <- floor(input$weeks*7/2)
    s[s$time >= (d0*24 + 16) & s$time <= (d0*24 + 40), ]
  })

  long_state <- function(df, vars, labs) {
    out <- df[, c("time", vars)]
    names(out) <- c("time", labs)
    tidyr::pivot_longer(out, -time, names_to = "state", values_to = "value")
  }

  ## ---- tab 1 --------------------------------------------------------------
  output$p_arousal <- renderPlot({
    s <- sim()$act
    d <- long_state(s, c("AROU", "COND", "SEFF", "DBAS"),
                    c("Hyperarousal AROU", "Conditioned arousal COND", "Sleep effort SEFF", "Beliefs DBAS"))
    ggplot(d, aes(time/24, value, colour = state)) + geom_line(linewidth = 0.6) +
      labs(x = "Day", y = "State value", colour = NULL,
           title = "Hyperarousal and its drivers") + theme_ins()
  })

  output$p_hpa <- renderPlot({
    s <- night()
    d <- long_state(s, c("CORT", "CRH", "ACTH"),
                    c("Cortisol (ug/dL)", "CRH", "ACTH"))
    ggplot(d, aes(time %% 24, value, colour = state)) + geom_line(linewidth = 0.7) +
      labs(x = "Time (clock h)", y = NULL, colour = NULL,
           title = "HPA axis — nocturnal cortisol nadir rise") + theme_ins()
  })

  output$t_profile <- renderDT({
    s <- sim()$act; n <- nrow(s)
    a <- sim()$nact; b <- sim()$nctl
    wk <- function(x) x[x$night > max(x$night) - 7, ]
    out <- data.frame(
      item = c("Hyperarousal AROU", "Conditioned arousal COND", "Sleep effort SEFF", "Dysfunctional beliefs DBAS",
               "Time in bed TIB (h)", "Circadian phase PHI (h)", "SOL (min)", "WASO (min)",
               "TST (min)", "SE (%)", "ISI"),
      active = round(c(s$AROU[n], s$COND[n], s$SEFF[n], s$DBAS[n], s$TIBS[n], s$PHI[n],
                       mean(wk(a)$SOL), mean(wk(a)$WASO), mean(wk(a)$TST),
                       mean(wk(a)$SE), mean(wk(a)$ISI)), 2),
      control = round(c(sim()$ctl$AROU[n], sim()$ctl$COND[n], sim()$ctl$SEFF[n],
                        sim()$ctl$DBAS[n], sim()$ctl$TIBS[n], sim()$ctl$PHI[n],
                        mean(wk(b)$SOL), mean(wk(b)$WASO), mean(wk(b)$TST),
                        mean(wk(b)$SE), mean(wk(b)$ISI)), 2)
    )
    names(out) <- c("Item", "Treatment arm", "Control arm")
    out
  }, options = list(dom = "t", pageLength = 12))

  ## ---- tab 2 --------------------------------------------------------------
  output$p_night <- renderPlot({
    s <- night()
    d <- long_state(s, c("S", "Cwake", "AROU", "Wdrv", "drive", "SLP"),
                    c("Process S", "C(t) arousal signal", "Hyperarousal A", "Ascending arousal drive",
                      "Net sleep drive", "Sleep probability SLP"))
    ggplot(d, aes(time %% 24, value, colour = state)) +
      geom_hline(yintercept = 0, colour = "grey70", linetype = 2) +
      geom_line(linewidth = 0.7) +
      scale_x_continuous(breaks = seq(0, 24, 3)) +
      labs(x = "Time (clock h)", y = NULL, colour = NULL,
           title = "A representative night — two processes, arousal, switch") + theme_ins()
  })

  output$p_switch <- renderPlot({
    s <- night()
    ggplot(s, aes(drive, SLP)) + geom_path(linewidth = 0.7, colour = "#3949ab") +
      labs(x = "Net sleep drive", y = "Sleep probability",
           title = "Switch trajectory (hysteresis)") + theme_ins()
  })

  output$p_drivecomp <- renderPlot({
    s <- night()
    d <- long_state(s, c("adoEff", "Ebz", "Emt", "oxBlk", "h1blk"),
                    c("Effective adenosine", "BzRA effect", "Melatonergic effect",
                      "Orexin blockade", "H1 blockade"))
    ggplot(d, aes(time %% 24, value, colour = state)) + geom_line(linewidth = 0.7) +
      labs(x = "Time (clock h)", y = NULL, colour = NULL,
           title = "Pharmacological drive components") + theme_ins()
  })

  ## ---- tab 3 --------------------------------------------------------------
  output$p_phi <- renderPlot({
    s <- sim()$act; c0 <- sim()$ctl
    d <- rbind(data.frame(time = s$time, PHI = s$PHI, arm = "Treatment arm"),
               data.frame(time = c0$time, PHI = c0$PHI, arm = "Control arm"))
    ggplot(d, aes(time/24, PHI, colour = arm)) + geom_line(linewidth = 0.6) +
      labs(x = "Day", y = "Phase shift Φ (h, + = delay)", colour = NULL,
           title = "Circadian phase — outcome of light-PRC entrainment") + theme_ins()
  })

  output$p_mel <- renderPlot({
    s <- night()
    ggplot(s, aes(time %% 24)) +
      geom_line(aes(y = MELP, colour = "Endogenous (pg/mL)"), linewidth = 0.8) +
      geom_line(aes(y = CMLX, colour = "Exogenous (pg/mL)"), linewidth = 0.8) +
      geom_line(aes(y = 100*occMT, colour = "MT occupancy (%)"), linetype = 2) +
      labs(x = "Time (clock h)", y = NULL, colour = NULL,
           title = "Melatonin — endogenous · exogenous · receptor occupancy") + theme_ins()
  })

  output$p_temp <- renderPlot({
    s <- night()
    ggplot(s, aes(time %% 24, TEMPC)) + geom_line(colour = "#c62828", linewidth = 0.8) +
      labs(x = "Time (clock h)", y = "Core body temperature (°C)",
           title = "Core body temperature — sleep maintenance is most stable near the nadir") + theme_ins()
  })

  output$p_light <- renderPlot({
    s <- night()
    ggplot(s, aes(time %% 24, lux + 1)) + geom_line(colour = "#f9a825", linewidth = 0.8) +
      scale_y_log10() +
      labs(x = "Time (clock h)", y = "Illuminance (lux, log)",
           title = "Light exposure profile") + theme_ins()
  })

  ## ---- tab 4 --------------------------------------------------------------
  output$p_pk <- renderPlot({
    s <- night()
    cs <- c("CZOL", "CESZ", "CSUV", "CLEM", "CDAR", "CRAM", "CRM2",
            "CDOX", "CTRZ", "CCAF", "CETH")
    keep <- cs[sapply(cs, function(v) max(s[[v]], na.rm = TRUE) > 1e-9)]
    if (!length(keep)) return(ggplot() + labs(title = "No drug administered") + theme_ins())
    d <- long_state(s, keep, keep)
    ggplot(d, aes(time %% 24, value, colour = state)) + geom_line(linewidth = 0.8) +
      labs(x = "Time (clock h)", y = "Concentration (mg/L; alcohol g/L)", colour = NULL,
           title = "Plasma concentration — a representative night") + theme_ins()
  })

  output$p_occ <- renderPlot({
    s <- night()
    d <- long_state(s, c("occBZ", "oxBlk", "occMT", "h1blk", "e5ht"),
                    c("GABA-A BZ site", "Orexin receptor blockade", "MT1/MT2 occupancy",
                      "H1 blockade", "5-HT2A blockade"))
    ggplot(d, aes(time %% 24, 100*value, colour = state)) + geom_line(linewidth = 0.8) +
      labs(x = "Time (clock h)", y = "Occupancy/blockade (%)", colour = NULL,
           title = "Target occupancy") + theme_ins()
  })

  ## ---- tab 5 --------------------------------------------------------------
  metrics_long <- reactive({
    a <- sim()$nact; b <- sim()$nctl
    f <- function(d, arm) data.frame(
      night = d$night, arm = arm,
      SOL = d$SOL, WASO = d$WASO, TST = d$TST, SE = d$SE)
    tidyr::pivot_longer(rbind(f(a, "Treatment arm"), f(b, "Control arm")),
                        c(SOL, WASO, TST, SE), names_to = "metric")
  })

  output$p_metrics <- renderPlot({
    ggplot(metrics_long(), aes(night, value, colour = arm)) +
      geom_line(linewidth = 0.6) +
      facet_wrap(~metric, scales = "free_y") +
      labs(x = "Night", y = NULL, colour = NULL,
           title = "Nightly sleep metrics — treatment arm vs untreated control arm") + theme_ins()
  })

  output$t_weekly <- renderDT({
    a <- sim()$nact; b <- sim()$nctl
    wk <- function(d) d %>% mutate(week = ceiling(night/7)) %>%
      group_by(week) %>% summarise(across(c(SOL, WASO, TST, SE, ISI), mean),
                                   .groups = "drop")
    A <- wk(a); B <- wk(b)
    out <- data.frame(week = A$week,
      SOL = round(A$SOL), SOLd = round(A$SOL - B$SOL),
      WASO = round(A$WASO), WASOd = round(A$WASO - B$WASO),
      TST = round(A$TST), TSTd = round(A$TST - B$TST),
      SEpct = round(A$SE, 1), ISI = round(A$ISI, 1),
      ISId = round(A$ISI - B$ISI, 1))
    names(out) <- c("Week", "SOL", "SOL \u0394", "WASO", "WASO \u0394",
                    "TST", "TST \u0394", "SE %", "ISI", "ISI \u0394")
    datatable(out, options = list(dom = "tp", pageLength = 12), rownames = FALSE)
  })

  ## ---- tab 6 --------------------------------------------------------------
  output$p_stage <- renderPlot({
    a <- sim()$nact
    d <- tidyr::pivot_longer(a[, c("night", "N3", "REM")], -night)
    ggplot(d, aes(night, value, colour = name)) + geom_line(linewidth = 0.7) +
      labs(x = "Night", y = "min/night", colour = NULL,
           title = "Sleep stages — absolute time") + theme_ins()
  })

  output$p_stagepct <- renderPlot({
    a <- sim()$nact
    d <- data.frame(night = a$night,
                    "N3 %" = 100*a$N3/pmax(a$TST, 1),
                    "REM %" = 100*a$REM/pmax(a$TST, 1), check.names = FALSE)
    d <- tidyr::pivot_longer(d, -night)
    ggplot(d, aes(night, value, colour = name)) + geom_line(linewidth = 0.7) +
      labs(x = "Night", y = "% of TST", colour = NULL,
           title = "Sleep stages — proportion of TST") + theme_ins()
  })

  output$p_naw <- renderPlot({
    a <- sim()$nact; b <- sim()$nctl
    d <- rbind(data.frame(night = a$night, NAW = a$NAW, arm = "Treatment arm"),
               data.frame(night = b$night, NAW = b$NAW, arm = "Control arm"))
    ggplot(d, aes(night, NAW, colour = arm)) + geom_line(linewidth = 0.6) +
      labs(x = "Night", y = "Number of awakenings", colour = NULL,
           title = "Nocturnal awakenings") + theme_ins()
  })

  ## ---- tab 7 --------------------------------------------------------------
  output$p_loop <- renderPlot({
    s <- sim()$act; c0 <- sim()$ctl
    mk <- function(df, arm) long_state(df, c("COND", "SEFF", "DBAS", "MISP"),
      c("Conditioned arousal COND", "Sleep effort SEFF", "Dysfunctional beliefs DBAS", "Sleep misperception MISP")) %>%
      mutate(arm = arm)
    d <- rbind(mk(s, "Treatment arm"), mk(c0, "Control arm"))
    ggplot(d, aes(time/24, value, colour = arm)) + geom_line(linewidth = 0.6) +
      facet_wrap(~state, scales = "free_y") +
      labs(x = "Day", y = NULL, colour = NULL,
           title = "3-P perpetuating loop state variables") + theme_ins()
  })

  output$p_tib <- renderPlot({
    s <- sim()$act; c0 <- sim()$ctl
    d <- rbind(data.frame(time = s$time, TIB = s$TIBS, SE = 100*s$SEBAR, arm = "Treatment arm"),
               data.frame(time = c0$time, TIB = c0$TIBS, SE = 100*c0$SEBAR, arm = "Control arm"))
    ggplot(d, aes(time/24)) +
      geom_line(aes(y = TIB, colour = arm), linewidth = 0.7) +
      geom_line(aes(y = SE/10, colour = arm), linetype = 2) +
      scale_y_continuous("Time in bed TIB (h)",
                         sec.axis = sec_axis(~.*10, name = "Sleep efficiency (%)")) +
      labs(x = "Day", colour = NULL,
           title = "Time-in-bed titration — solid line TIB, dashed line sleep efficiency") + theme_ins()
  })

  ## ---- tab 8 --------------------------------------------------------------
  output$p_isi <- renderPlot({
    s <- sim()$act; c0 <- sim()$ctl
    d <- rbind(data.frame(time = s$time, ISI = s$ISI, arm = "Treatment arm"),
               data.frame(time = c0$time, ISI = c0$ISI, arm = "Control arm"))
    ggplot(d, aes(time/24, ISI, colour = arm)) + geom_line(linewidth = 0.7) +
      geom_hline(yintercept = 8, linetype = 2, colour = "grey50") +
      annotate("text", x = 1, y = 8.8, label = "Remission threshold ISI<8",
               hjust = 0, size = 3.2, colour = "grey40") +
      labs(x = "Day", y = "ISI (0-28)", colour = NULL,
           title = "Insomnia Severity Index") + theme_ins()
  })

  output$p_subj <- renderPlot({
    a <- sim()$nact
    d <- data.frame(night = a$night, "objective TST" = a$TST,
                    "subjective TST" = a$sTST, check.names = FALSE)
    d <- tidyr::pivot_longer(d, -night)
    ggplot(d, aes(night, value, colour = name)) + geom_line(linewidth = 0.7) +
      labs(x = "Night", y = "min", colour = NULL,
           title = "Objective vs subjective sleep time — magnitude of misperception") + theme_ins()
  })

  output$p_safety <- renderPlot({
    s <- sim()$act
    d <- long_state(s, c("TOLB", "WDR", "RESBAR"),
                    c("BzRA tolerance", "Withdrawal/rebound", "Daytime residual sedation"))
    ggplot(d, aes(time/24, value, colour = state)) + geom_line(linewidth = 0.7) +
      labs(x = "Day", y = NULL, colour = NULL,
           title = "Tolerance · rebound · residual sedation") + theme_ins()
  })

  output$p_conseq <- renderPlot({
    s <- sim()$act
    d <- long_state(s, c("DEP", "IL6", "GLYM"),
                    c("Depression burden", "IL-6 (pg/mL)", "Glymphatic clearance debt"))
    ggplot(d, aes(time/24, value, colour = state)) + geom_line(linewidth = 0.7) +
      labs(x = "Day", y = NULL, colour = NULL,
           title = "Long-term outcome surrogates") + theme_ins()
  })

  output$t_endpoint <- renderDT({
    a <- sim()$nact; b <- sim()$nctl
    l <- function(d) d[d$night > max(d$night) - 7, ]
    A <- l(a); B <- l(b)
    out <- data.frame(
      endpoint = c("SOL (min)", "WASO (min)", "TST (min)", "SE (%)",
                     "Subjective TST (min)", "ISI", "N3 (min)", "REM (min)",
                     "Hangover index (0-100)", "Remission (ISI<8 & SE>85%)"),
      active = c(round(mean(A$SOL)), round(mean(A$WASO)), round(mean(A$TST)),
                 round(mean(A$SE), 1), round(mean(A$sTST)), round(mean(A$ISI), 1),
                 round(mean(A$N3)), round(mean(A$REM)), round(mean(A$HANG)),
                 ifelse(mean(A$ISI) < 8 && mean(A$SE) > 85, "Yes", "No")),
      control = c(round(mean(B$SOL)), round(mean(B$WASO)), round(mean(B$TST)),
                  round(mean(B$SE), 1), round(mean(B$sTST)), round(mean(B$ISI), 1),
                  round(mean(B$N3)), round(mean(B$REM)), round(mean(B$HANG)),
                  ifelse(mean(B$ISI) < 8 && mean(B$SE) > 85, "Yes", "No"))
    )
    names(out) <- c("Endpoint", "Treatment arm",
                    "Control arm")
    out
  }, options = list(dom = "t", pageLength = 12))

  ## ---- tab 9 --------------------------------------------------------------
  scen <- eventReactive(input$runscen, {
    shiny::req(length(input$scen) > 0)   # mrgsolve also exports req()
    withProgress(message = "Running scenarios...", value = 0, {
      out <- lapply(seq_along(input$scen), function(i) {
        incProgress(1/length(input$scen), detail = input$scen[i])
        r <- ins_run_scenario(input$scen[i], weeks = input$weeks)
        r$nightly$scenario <- input$scen[i]
        r$nightly
      })
      do.call(rbind, out)
    })
  })

  output$p_scen <- renderPlot({
    d <- scen()
    dl <- tidyr::pivot_longer(d[, c("night", "scenario", "SOL", "WASO", "TST", "ISI")],
                              c(SOL, WASO, TST, ISI), names_to = "metric")
    ggplot(dl, aes(night, value, colour = scenario)) + geom_line(linewidth = 0.6) +
      facet_wrap(~metric, scales = "free_y") +
      labs(x = "Night", y = NULL, colour = NULL, title = "Scenario comparison") +
      theme_ins() + theme(legend.text = element_text(size = 8))
  })

  output$t_scen <- renderDT({
    d <- scen()
    s <- d %>% group_by(scenario) %>%
      filter(night > max(night) - 7) %>%
      summarise(SOL = round(mean(SOL)), WASO = round(mean(WASO)),
                TST = round(mean(TST)), `SE %` = round(mean(SE), 1),
                ISI = round(mean(ISI), 1), sTST = round(mean(sTST)),
                N3 = round(mean(N3)), REM = round(mean(REM)),
                HANG = round(mean(HANG)), .groups = "drop")
    datatable(s, options = list(dom = "t"), rownames = FALSE)
  })
}

shinyApp(ui, server)
