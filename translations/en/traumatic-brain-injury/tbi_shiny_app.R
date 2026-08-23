## =============================================================================
##  Severe Traumatic Brain Injury -- QSP Shiny dashboard
##  Severe Traumatic Brain Injury QSP Dashboard
## =============================================================================
##
##  DESIGN PRINCIPLE
##  ----------------
##  This dashboard is built so that no number can be read in isolation, because
##  the single most consequential error in ICP management is reading the ICP
##  number by itself.  Three rules are enforced by the layout:
##
##   1. ICP is NEVER shown without the compensatory reserve beside it.  A green
##      ICP over a red reserve is the dangerous state, and it must look
##      dangerous.
##   2. Every ICP therapy is shown with its CURRENCY on the same panel -- the
##      oxygen it spends, the sodium it spends, the arterial pressure it spends.
##   3. Global monitors (ICP, CPP, SjvO2) and regional monitors (PbtO2, CBF in
##      the penumbra, LPR) are always plotted together, because the entire
##      point of the model is that they can disagree.
##
##  Run:  shiny::runApp("tbi_shiny_app.R")
##  Requires: shiny, mrgsolve, dplyr, tidyr, ggplot2, DT
##  Model:  tbi_mrgsolve_model.R must be in the same directory.
## =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

mod <- mread_cache("tbi", "tbi_mrgsolve_model.R")

THEME <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")

## thresholds that appear on more than one plot -------------------------------
THR <- list(ICP = 22, CPP = 60, PbtO2 = 20, LPR = 40, Na = 155, Osm = 320,
            PRx = 0.25)

band <- function(y, col = "#e74c3c", lty = 2)
  geom_hline(yintercept = y, colour = col, linetype = lty, linewidth = 0.4)

## ---------------------------------------------------------------------------
## Helper: build the event table for a scenario
## ---------------------------------------------------------------------------
build_events <- function(input) {
  ev_list <- list()
  if (input$osm_agent != "none" && input$osm_n > 0) {
    for (i in seq_len(input$osm_n)) {
      t0 <- input$osm_t1 * 60 + (i - 1) * input$osm_int * 60
      if (input$osm_agent == "hts3") {
        ev_list[[length(ev_list) + 1]] <-
          ev(time = t0, amt = input$osm_vol * 0.513, cmt = "Na_ecf", tinf = 30)
        ev_list[[length(ev_list) + 1]] <-
          ev(time = t0, amt = input$osm_vol, cmt = "V_ecf", tinf = 30)
      } else if (input$osm_agent == "hts23") {
        ev_list[[length(ev_list) + 1]] <-
          ev(time = t0, amt = input$osm_vol23 * 4.004, cmt = "Na_ecf", tinf = 12)
        ev_list[[length(ev_list) + 1]] <-
          ev(time = t0, amt = input$osm_vol23, cmt = "V_ecf", tinf = 12)
      } else if (input$osm_agent == "mann") {
        mmol <- input$mann_dose * input$wt * 1000 / 182.17
        ev_list[[length(ev_list) + 1]] <-
          ev(time = t0, amt = mmol, cmt = "Mann_c", tinf = 20)
        ev_list[[length(ev_list) + 1]] <-
          ev(time = t0, amt = input$mann_dose * input$wt * 5, cmt = "V_ecf",
             tinf = 20)
      }
    }
  }
  if (input$txa) {
    t0 <- input$txa_t * 60
    ev_list[[length(ev_list) + 1]] <- ev(time = t0, amt = 1000, cmt = "TXA_c",
                                         tinf = 10)
    ev_list[[length(ev_list) + 1]] <- ev(time = t0 + 10, amt = 1000,
                                         cmt = "TXA_c", tinf = 480)
  }
  # timed parameter changes: hyperventilation window, cooling, insult
  if (input$hv) {
    ev_list[[length(ev_list) + 1]] <- ev(time = input$hv_t1 * 60, amt = 0, cmt = 0, evid = 2, PACO2 = input$hv_co2)
    ev_list[[length(ev_list) + 1]] <- ev(time = input$hv_t2 * 60, amt = 0, cmt = 0, evid = 2, PACO2 = input$paco2)
  }
  if (input$cool) {
    ev_list[[length(ev_list) + 1]] <- ev(time = input$cool_t1 * 60, amt = 0, cmt = 0, evid = 2,
                                         T_TARGET = input$cool_temp)
    ev_list[[length(ev_list) + 1]] <- ev(time = input$cool_t2 * 60, amt = 0, cmt = 0, evid = 2, T_TARGET = 37)
  }
  if (input$insult) {
    ev_list[[length(ev_list) + 1]] <-
      ev(time = input$ins_t * 60, amt = 0, cmt = 0, evid = 2,
         MAP_BASE = input$ins_map, SAO2 = input$ins_sao2)
    ev_list[[length(ev_list) + 1]] <-
      ev(time = input$ins_t * 60 + input$ins_dur, amt = 0, cmt = 0, evid = 2,
         MAP_BASE = input$map, SAO2 = 0.98)
  }
  if (length(ev_list) == 0) return(NULL)
  bind_rows(ev_list) %>% arrange(time)
}

simulate_arm <- function(input, overrides = list()) {
  pp <- list(SEV = input$sev, AGE = input$age, WT = input$wt,
             MAP_BASE = input$map, PACO2 = input$paco2,
             PROP_RATE = input$prop, THIO_RATE = input$thio,
             NE_RATE = input$ne, EVD_OPEN = as.numeric(input$evd),
             EVD_SET = input$evd_set, AUTOREG = input$autoreg,
             CRANI_T = if (input$crani) input$crani_t * 60 else -1,
             PVI = input$pvi, Hb = input$hb, FLUID = input$fluid)
  pp <- modifyList(pp, overrides)
  m <- mod %>% param(pp)
  e <- build_events(input)
  if (!is.null(e)) m <- m %>% ev(e)
  m %>% mrgsim(end = input$tmax * 60, delta = 5, recsort = 3) %>%
    as_tibble() %>% mutate(hour = time / 60)
}

## outcome model -- mirrors gose_unfavourable() in the Python reference --------
burden <- function(d) {
  g <- diff(d$time)[1] / 60
  list(Fc   = tail(d$F_core, 1),
       Dicp = sum(pmax(0, d$ICP - 20)) * g,
       Dcpp = sum(pmax(0, 60 - d$CPP_o)) * g,
       Dpbt = sum(pmax(0, 20 - d$PbtO2)) * g,
       pupil = if (sum(d$ICP > 40) * g > 1) 2 else if (sum(d$ICP > 35) * g > 0.5) 1 else 0)
}
p_unfav <- function(d, age) {
  b <- burden(d)
  z <- -1.64 + 0.35 * (age - 40) / 10 + 6.00 * b$Fc +
    0.45 * log1p(b$Dicp / 20) + 0.30 * log1p(b$Dcpp / 40) +
    0.35 * log1p(b$Dpbt / 60) + 0.90 * b$pupil
  1 / (1 + exp(-z))
}
p_death <- function(d, age) {
  b <- burden(d)
  z <- -2.20 + 0.40 * (age - 40) / 10 + 5.50 * b$Fc +
    0.40 * log1p(b$Dicp / 20) + 0.90 * b$pupil
  1 / (1 + exp(-z))
}

## =============================================================================
##  UI
## =============================================================================
ui <- fluidPage(
  titlePanel("Severe Traumatic Brain Injury QSP Dashboard · Severe TBI QSP dashboard"),
  tags$p(style = "color:#666;margin-top:-8px;",
         "The skull is a closed box — every therapy removes one of four volumes, ",
         "and each pays in a different currency.  ",
         em("The skull is a closed box: every therapy removes one of four ",
            "volumes, and each pays in a different currency.")),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Patient"),
      sliderInput("sev", "Injury severity", 0, 1, 0.62, 0.01),
      sliderInput("age", "Age (yr)", 16, 85, 42, 1),
      sliderInput("wt", "Weight (kg)", 40, 130, 75, 1),
      sliderInput("pvi", "PVI craniospinal compliance (mL)", 12, 40, 26, 1),
      sliderInput("hb", "Hb (g/dL)", 6, 16, 12.5, 0.5),
      sliderInput("autoreg", "Autoregulatory gain (× normal)",
                  0, 1.5, 1.0, 0.05),
      helpText("0 = complete loss. This one value flips the sign of a vasopressor's effect on ICP."),
      hr(),
      h4("Basic Management · Tier 0"),
      sliderInput("map", "MAP target (mmHg)", 55, 120, 88, 1),
      sliderInput("paco2", "PaCO2 (mmHg)", 25, 50, 38, 1),
      sliderInput("prop", "Propofol (mg/kg/h)", 0, 8, 3, 0.5),
      sliderInput("ne", "Noradrenaline (µg/kg/min)", 0, 0.8, 0, 0.01),
      sliderInput("fluid", "Maintenance fluid (mL/min)", 0.5, 4, 1.6, 0.1),
      hr(),
      h4("Tier 1–2"),
      checkboxInput("evd", "External ventricular drain, EVD", FALSE),
      conditionalPanel("input.evd",
                       sliderInput("evd_set", "Drainage threshold (mmHg)", 5, 25, 12, 1)),
      selectInput("osm_agent", "Osmotherapy",
                  c("None" = "none", "3% NaCl" = "hts3",
                    "23.4% NaCl" = "hts23", "Mannitol" = "mann")),
      conditionalPanel("input.osm_agent == 'hts3'",
                       sliderInput("osm_vol", "3% NaCl (mL)", 100, 750, 250, 50)),
      conditionalPanel("input.osm_agent == 'hts23'",
                       sliderInput("osm_vol23", "23.4% NaCl (mL)", 10, 60, 30, 5)),
      conditionalPanel("input.osm_agent == 'mann'",
                       sliderInput("mann_dose", "Mannitol (g/kg)", 0.1, 1.5, 0.5, 0.05)),
      conditionalPanel("input.osm_agent != 'none'",
                       sliderInput("osm_n", "doses", 1, 8, 1, 1),
                       sliderInput("osm_t1", "Time of first dose (h)", 0, 48, 6, 0.5),
                       sliderInput("osm_int", "Dosing interval (h)", 2, 12, 4, 1)),
      checkboxInput("hv", "Hyperventilation", FALSE),
      conditionalPanel("input.hv",
                       sliderInput("hv_co2", "PaCO2 target", 22, 36, 30, 1),
                       sliderInput("hv_t1", "Start (h)", 0, 48, 6, 0.5),
                       sliderInput("hv_t2", "End (h)", 0, 96, 30, 0.5)),
      hr(),
      h4("Tier 3"),
      sliderInput("thio", "Thiopental (mg/kg/h)", 0, 8, 0, 0.5),
      checkboxInput("cool", "Hypothermia", FALSE),
      conditionalPanel("input.cool",
                       sliderInput("cool_temp", "Target temperature (°C)", 32, 36, 33.5, 0.5),
                       sliderInput("cool_t1", "Start (h)", 0, 48, 2, 1),
                       sliderInput("cool_t2", "Rewarming start (h)", 6, 96, 50, 1)),
      checkboxInput("crani", "Decompressive craniectomy", FALSE),
      conditionalPanel("input.crani",
                       sliderInput("crani_t", "Time of surgery (h)", 0, 72, 24, 1)),
      hr(),
      h4("Other"),
      checkboxInput("txa", "Tranexamic acid", FALSE),
      conditionalPanel("input.txa",
                       sliderInput("txa_t", "Administration time (h)", 0, 12, 1, 0.5)),
      checkboxInput("insult", "Secondary insult", FALSE),
      conditionalPanel("input.insult",
                       sliderInput("ins_t", "Time (h)", 0, 48, 6, 0.5),
                       sliderInput("ins_dur", "Duration (min)", 5, 120, 30, 5),
                       sliderInput("ins_map", "MAP drop (mmHg)", 40, 80, 55, 1),
                       sliderInput("ins_sao2", "SaO2", 0.6, 0.95, 0.75, 0.01)),
      hr(),
      sliderInput("tmax", "Simulation period (h)", 12, 168, 72, 12)
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs", type = "tabs",

        ## -- 1 ---------------------------------------------------------------
        tabPanel(
          "① Patient Profile",
          br(),
          fluidRow(
            column(3, wellPanel(h4("Peak ICP"), h2(textOutput("kpi_icp")),
                                span(style = "color:#888", "mmHg"))),
            column(3, wellPanel(h4("ICP dose"), h2(textOutput("kpi_dose")),
                                span(style = "color:#888", "mmHg·h > 20"))),
            column(3, wellPanel(h4("Final core"), h2(textOutput("kpi_core")),
                                span(style = "color:#888", "non-viable fraction"))),
            column(3, wellPanel(h4("P(GOS-E 1–4)"), h2(textOutput("kpi_gose")),
                                span(style = "color:#888", "unfavourable at 6 mo")))
          ),
          h4("The whole course on one axis"),
          plotOutput("p_overview", height = "520px"),
          helpText("Shading: red = ICP > 22 mmHg, blue = PbtO2 < 20 mmHg. ",
                   "The point of this model is the interval where the two shaded regions do not overlap.")
        ),

        ## -- 2 ---------------------------------------------------------------
        tabPanel(
          "② Pressure-Volume Reserve",
          br(),
          tags$div(style = "background:#fdf2e9;padding:10px;border-left:4px solid #e67e22;",
                   strong("ICP tells you almost nothing about the space you have left."),
                   " Because dP/dV ∝ P, the same 10 mL is invisible at ICP 8 and ",
                   "fatal at ICP 25. The two curves above must be read together."),
          br(),
          plotOutput("p_reserve", height = "440px"),
          h4("Where the buffer actually is"),
          plotOutput("p_buffers", height = "300px")
        ),

        ## -- 3 ---------------------------------------------------------------
        tabPanel(
          "③ Haemodynamics · Autoregulation",
          br(),
          plotOutput("p_haemo", height = "420px"),
          h4("Autoregulation curve (Lassen) · emergent, not coded"),
          plotOutput("p_lassen", height = "320px"),
          helpText("No autoregulation curve exists inside the model. ",
                   "It is the outcome produced by a single saturating sigmoid.")
        ),

        ## -- 4 ---------------------------------------------------------------
        tabPanel(
          "④ Oxygen Delivery: Global vs Local",
          br(),
          tags$div(style = "background:#eaf2f8;padding:10px;border-left:4px solid #2e86c1;",
                   strong("SjvO2 is a global average and PbtO2 is a local probe."),
                   " Because the dying tissue is a minority of compartments, jugular venous ",
                   "oxygen saturation can be normal while brain-tissue oxygen tension is at a dangerous level."),
          br(),
          plotOutput("p_oxygen", height = "460px"),
          h4("Hours the monitors disagree"),
          tableOutput("t_disagree")
        ),

        ## -- 5 ---------------------------------------------------------------
        tabPanel(
          "⑤ Osmotherapy",
          br(),
          tags$div(style = "background:#e8f5e9;padding:10px;border-left:4px solid #1b5e20;",
                   strong("Osmotherapy does not treat the lesion. It taxes the healthy brain."),
                   " In the injured region the blood-brain barrier is disrupted, so the reflection coefficient σ is close to 0, ",
                   "and therefore the osmotic gradient cannot do any work."),
          br(),
          plotOutput("p_osmo", height = "440px"),
          h4("Which region gave up the water"),
          plotOutput("p_water", height = "300px"),
          h4("Rebound and tachyphylaxis"),
          plotOutput("p_rebound", height = "300px")
        ),

        ## -- 6 ---------------------------------------------------------------
        tabPanel(
          "⑥ Currency",
          br(),
          tags$div(style = "background:#fce4ec;padding:10px;border-left:4px solid #880e4f;",
                   strong("The cranial cavity has only four volumes."),
                   " So every therapy can only do one of four things, ",
                   "and each pays in a different currency."),
          br(),
          actionButton("run_currency", "Price every therapy",
                       class = "btn-primary"),
          br(), br(),
          DTOutput("t_currency"),
          helpText("Result of applying each therapy for 1 hour, at the same time point in the same patient. ",
                   "Look not at the ICP column but at the columns to its right.")
        ),

        ## -- 7 ---------------------------------------------------------------
        tabPanel(
          "⑦ Metabolism · Microdialysis",
          br(),
          plotOutput("p_metab", height = "460px"),
          helpText("LPR > 40 indicates ischaemia or mitochondrial dysfunction. ",
                   "If PbtO2 is normal while LPR is high, it is the latter."),
          h4("Excitotoxic cascade"),
          plotOutput("p_excito", height = "320px")
        ),

        ## -- 8 ---------------------------------------------------------------
        tabPanel(
          "⑧ Oedema & Barrier",
          br(),
          plotOutput("p_oedema", height = "440px"),
          helpText("Blood-brain barrier opening is biphasic: ",
                   "mechanical injury (0–6 h, t½ ≈ 3 h) and inflammatory opening (24–72 h)."),
          h4("Haematoma and coagulopathy"),
          plotOutput("p_hem", height = "300px")
        ),

        ## -- 9 ---------------------------------------------------------------
        tabPanel(
          "⑨ Scenario Comparison",
          br(),
          checkboxGroupInput(
            "cmp", "Treatment arms to compare", inline = TRUE,
            choices = c("No therapy" = "none",
                        "Tier 0 (sedation+CPP)" = "t0",
                        "Tier 1 (+EVD+HTS)" = "t1",
                        "Tier 2 (+hyperventilation)" = "t2",
                        "Tier 3 (+barbiturate)" = "t3",
                        "Decompressive craniectomy" = "crani",
                        "Loss of autoregulation" = "noauto"),
            selected = c("none", "t1", "crani")),
          actionButton("run_cmp", "Run comparison", class = "btn-primary"),
          br(), br(),
          plotOutput("p_compare", height = "480px"),
          DTOutput("t_compare")
        ),

        ## -- 10 --------------------------------------------------------------
        tabPanel(
          "⑩ Biomarkers",
          br(),
          plotOutput("p_biomark", height = "420px"),
          helpText("All four markers are released from the same process (penumbra→core transition). ",
                   "What differs is only the half-life, and that alone is what makes the clinical meaning diverge: ",
                   "the fast marker reports instantaneous rate, the slow marker reports cumulative amount."),
          tableOutput("t_biomark")
        ),

        ## -- 11 --------------------------------------------------------------
        tabPanel(
          "⑪ Clinical Outcome",
          br(),
          fluidRow(
            column(6, h4("Burden measures"), tableOutput("t_burden")),
            column(6, h4("Predicted outcome"), tableOutput("t_outcome"))
          ),
          plotOutput("p_dose", height = "340px"),
          hr(),
          tags$div(style = "background:#f2f4f4;padding:12px;",
                   h4("⚠️ What this model cannot do"),
                   tags$ul(
                     tags$li("It cannot reproduce the harms of hypothermia, because pneumonia, coagulopathy, ",
                             "shivering, and arrhythmia are not in the model (see output R12)."),
                     tags$li("The prognosis logistic is fitted to three calibration points and ",
                             "should be read as an ordering, not as individual probabilities."),
                     tags$li("Surgical complications, cranial-defect syndrome, and rehabilitation effects are not modelled."),
                     tags$li("This is a single-patient deterministic model. There is no inter-individual variability.")))
        )
      )
    )
  )
)

## =============================================================================
##  SERVER
## =============================================================================
server <- function(input, output, session) {

  sim <- reactive(simulate_arm(input))

  ## ---- KPIs ---------------------------------------------------------------
  output$kpi_icp  <- renderText(sprintf("%.1f", max(sim()$ICP)))
  output$kpi_dose <- renderText(sprintf("%.0f", burden(sim())$Dicp))
  output$kpi_core <- renderText(sprintf("%.3f", tail(sim()$F_core, 1)))
  output$kpi_gose <- renderText(sprintf("%.0f%%", 100 * p_unfav(sim(), input$age)))

  ## ---- 1. overview --------------------------------------------------------
  output$p_overview <- renderPlot({
    d <- sim()
    long <- d %>%
      select(hour, ICP, CPP_o, PbtO2, CBF_pen, PRx, F_core) %>%
      pivot_longer(-hour)
    lv <- c(ICP = "ICP (mmHg)", CPP_o = "CPP (mmHg)",
            PbtO2 = "PbtO2 (mmHg)", CBF_pen = "CBF penumbra",
            PRx = "PRx", F_core = "non-viable fraction")
    long$name <- factor(lv[long$name], levels = lv)
    hl <- data.frame(name = factor(lv[c("ICP", "CPP_o", "PbtO2", "PRx")],
                                   levels = lv),
                     y = c(THR$ICP, THR$CPP, THR$PbtO2, THR$PRx))
    ggplot(long, aes(hour, value)) +
      geom_line(linewidth = 0.8, colour = "#2c3e50") +
      geom_hline(data = hl, aes(yintercept = y), colour = "#e74c3c",
                 linetype = 2, linewidth = 0.4) +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "hours post-injury", y = NULL) + THEME
  })

  ## ---- 2. reserve ---------------------------------------------------------
  output$p_reserve <- renderPlot({
    d <- sim()
    dd <- d %>% transmute(hour, ICP,
                          `compliance (mL/mmHg)` = denom,
                          `CSF buffer remaining (0-1)` = csf_avail,
                          `CSF volume (mL)` = V_csf) %>%
      pivot_longer(-hour)
    ggplot(dd, aes(hour, value)) +
      geom_line(linewidth = 0.9, colour = "#1f618d") +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "hours", y = NULL,
           title = "ICP is flat while reserve is gone") + THEME
  })

  output$p_buffers <- renderPlot({
    d <- sim()
    dd <- d %>% transmute(hour,
                          `CSF (mL)` = V_csf,
                          `venous blood (mL)` = Vv,
                          `arterial blood (mL)` = Va) %>%
      pivot_longer(-hour)
    ggplot(dd, aes(hour, value, colour = name)) +
      geom_line(linewidth = 0.9) +
      scale_colour_manual(values = c("#2e86c1", "#7d3c98", "#c0392b")) +
      labs(x = "hours", y = "mL", colour = NULL) + THEME
  })

  ## ---- 3. haemodynamics ---------------------------------------------------
  output$p_haemo <- renderPlot({
    d <- sim()
    dd <- d %>% transmute(hour, `MAP` = MAP_o, `ICP` = ICP, `CPP` = CPP_o,
                          `CBF global` = CBF100, `CBV (mL)` = CBV,
                          `Ca (mL/mmHg)` = Ca, `PRx` = PRx,
                          `Cushing (mmHg)` = cush) %>%
      pivot_longer(-hour)
    ggplot(dd, aes(hour, value)) + geom_line(linewidth = 0.8, colour = "#117a65") +
      facet_wrap(~name, scales = "free_y", ncol = 4) +
      labs(x = "hours", y = NULL) + THEME
  })

  output$p_lassen <- renderPlot({
    maps <- seq(40, 170, 10)
    res <- lapply(maps, function(m) {
      mod %>% param(SEV = 0, MAP_BASE = m, PROP_RATE = 0, PACO2 = 40) %>%
        mrgsim(end = 25, delta = 25) %>% as_tibble() %>% slice_tail(n = 1) %>%
        mutate(MAP_set = m)
    }) %>% bind_rows()
    ggplot(res, aes(CPP_o, CBF100)) +
      geom_line(linewidth = 1, colour = "#b7950b") +
      geom_point(size = 2, colour = "#b7950b") +
      annotate("rect", xmin = 60, xmax = 130, ymin = -Inf, ymax = Inf,
               alpha = 0.08, fill = "#f1c40f") +
      labs(x = "CPP (mmHg)", y = "CBF (mL/100 g/min)",
           title = "The autoregulation plateau is an output, not an input") + THEME
  })

  ## ---- 4. oxygen ----------------------------------------------------------
  output$p_oxygen <- renderPlot({
    d <- sim()
    dd <- d %>% transmute(hour,
                          `CBF global (mL/100g/min)` = CBF100,
                          `CBF penumbra (mL/100g/min)` = CBF_pen,
                          `SjvO2 (%) — GLOBAL` = SjvO2,
                          `PbtO2 (mmHg) — REGIONAL` = PbtO2,
                          `LPR` = LPR,
                          `supply adequacy r` = r_pen) %>%
      pivot_longer(-hour)
    ggplot(dd, aes(hour, value)) + geom_line(linewidth = 0.85, colour = "#af601a") +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "hours", y = NULL) + THEME
  })

  output$t_disagree <- renderTable({
    d <- sim(); g <- diff(d$time)[1] / 60
    data.frame(
      `Condition` = c("ICP < 22 but PbtO2 < 20",
                           "SjvO2 > 55% but PbtO2 < 20",
                           "PbtO2 < 20 (overall)",
                           "PbtO2 < 10 (critical)",
                           "LPR > 40 but PbtO2 > 20 (mitochondrial type)"),
      `Hours` = c(sum(d$ICP < 22 & d$PbtO2 < 20) * g,
                       sum(d$SjvO2 > 55 & d$PbtO2 < 20) * g,
                       sum(d$PbtO2 < 20) * g,
                       sum(d$PbtO2 < 10) * g,
                       sum(d$LPR > 40 & d$PbtO2 > 20) * g),
      check.names = FALSE)
  }, digits = 1)

  ## ---- 5. osmotherapy -----------------------------------------------------
  output$p_osmo <- renderPlot({
    d <- sim()
    dd <- d %>% transmute(hour, `ICP (mmHg)` = ICP,
                          `plasma Na (mmol/L)` = Na_p,
                          `plasma osmolality` = Osm_p,
                          `osmolar gap (mannitol)` = osm_gap,
                          `urine (mL/min)` = urine,
                          `GFR (L/min)` = gfr) %>%
      pivot_longer(-hour)
    hl <- data.frame(name = c("plasma Na (mmol/L)", "plasma osmolality"),
                     y = c(THR$Na, THR$Osm))
    ggplot(dd, aes(hour, value)) + geom_line(linewidth = 0.85, colour = "#1e8449") +
      geom_hline(data = hl, aes(yintercept = y), colour = "#e74c3c",
                 linetype = 2, linewidth = 0.4) +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "hours", y = NULL) + THEME
  })

  output$p_water <- renderPlot({
    d <- sim()
    dd <- d %>% transmute(hour,
                          `intact region (mL vs baseline)` = V_int - 950,
                          `injured region (mL vs baseline)` = V_inj - 170) %>%
      pivot_longer(-hour)
    ggplot(dd, aes(hour, value, colour = name)) + geom_line(linewidth = 0.95) +
      geom_hline(yintercept = 0, colour = "#999") +
      scale_colour_manual(values = c("#c0392b", "#2874a6")) +
      labs(x = "hours", y = "mL", colour = NULL,
           title = "Nearly all the water is lost from the healthy side") + THEME
  })

  output$p_rebound <- renderPlot({
    d <- sim()
    dd <- d %>% transmute(hour, `plasma mannitol (mmol/L)` = osm_gap,
                          `mannitol trapped in lesion` = Mann_inj_mM,
                          `barrier opening (0-1)` = BBB) %>%
      pivot_longer(-hour)
    ggplot(dd, aes(hour, value, colour = name)) + geom_line(linewidth = 0.9) +
      scale_colour_manual(values = c("#7d3c98", "#e67e22", "#16a085")) +
      labs(x = "hours", y = NULL, colour = NULL) + THEME
  })

  ## ---- 6. currency table --------------------------------------------------
  currency <- eventReactive(input$run_currency, {
    withProgress(message = "Calculating the price of each therapy...", value = 0, {
      base <- simulate_arm(input)
      # branch at the first crossing of ICP 22
      i <- which(base$ICP >= 22)[1]
      t0 <- if (is.na(i)) input$tmax * 60 / 3 else base$time[i]
      arms <- list(
        "nothing (control)"        = list(),
        "EVD to 12 mmHg"           = list(EVD_OPEN = 1, EVD_SET = 12),
        "PaCO2 -> 33"              = list(PACO2 = 33),
        "PaCO2 -> 28"              = list(PACO2 = 28),
        "propofol 6 mg/kg/h"       = list(PROP_RATE = 6),
        "thiopental 5 mg/kg/h"     = list(THIO_RATE = 5),
        "hypothermia 33.5 C"       = list(T_TARGET = 33.5),
        "noradrenaline +0.15"      = list(NE_RATE = input$ne + 0.15),
        "craniectomy"              = list(CRANI_T = t0))
      rows <- list()
      k <- 0
      for (nm in names(arms)) {
        k <- k + 1
        incProgress(1 / length(arms))
        d <- simulate_arm(input, arms[[nm]])
        j <- which.min(abs(d$time - (t0 + 60)))
        b <- which.min(abs(base$time - (t0 + 60)))
        rows[[nm]] <- data.frame(
          therapy = nm,
          ICP = round(d$ICP[j], 2),
          dICP = round(d$ICP[j] - base$ICP[b], 2),
          CPP = round(d$CPP_o[j], 1),
          MAP = round(d$MAP_o[j], 1),
          PbtO2 = round(d$PbtO2[j], 2),
          CBF_pen = round(d$CBF_pen[j], 2),
          Na = round(d$Na_p[j], 1),
          osm_gap = round(d$osm_gap[j], 1),
          urine = round(d$urine[j], 2))
      }
      bind_rows(rows)
    })
  })
  output$t_currency <- renderDT({
    datatable(currency(), rownames = FALSE,
              options = list(dom = "t", pageLength = 20)) %>%
      formatStyle("dICP", color = styleInterval(c(-0.2, 0.2),
                                                c("#1e8449", "#666", "#c0392b")))
  })

  ## ---- 7. metabolism ------------------------------------------------------
  output$p_metab <- renderPlot({
    d <- sim()
    dd <- d %>% transmute(hour, `CMRO2` = CMRO2_act, `LPR` = LPR,
                          `E_pump (0-1)` = E_pump,
                          `ionic failure (0-1)` = ionfail,
                          `brain lactate (mmol/L)` = Lac,
                          `brain glucose (mmol/L)` = Glc_br) %>%
      pivot_longer(-hour)
    ggplot(dd, aes(hour, value)) + geom_line(linewidth = 0.85, colour = "#1e8449") +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "hours", y = NULL) + THEME
  })
  output$p_excito <- renderPlot({
    d <- sim()
    dd <- d %>% transmute(hour, `glutamate (uM)` = Glu, `ECF K+ (mmol/L)` = K_ec,
                          `Ca_i` = Ca_i, `mitochondrial dysfunction` = MitoD,
                          `SD rate (/h)` = SD_rate) %>%
      pivot_longer(-hour)
    ggplot(dd, aes(hour, value)) + geom_line(linewidth = 0.85, colour = "#943126") +
      facet_wrap(~name, scales = "free_y", ncol = 5) +
      labs(x = "hours", y = NULL) + THEME
  })

  ## ---- 8. oedema ----------------------------------------------------------
  output$p_oedema <- renderPlot({
    d <- sim()
    dd <- d %>% transmute(hour, `excess brain water (mL)` = V_ed,
                          `barrier opening (0-1)` = BBB,
                          `intact osmolality` = Osm_int,
                          `lesion osmolality` = Osm_inj) %>%
      pivot_longer(-hour)
    ggplot(dd, aes(hour, value)) + geom_line(linewidth = 0.85, colour = "#0e6655") +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "hours", y = NULL) + THEME
  })
  output$p_hem <- renderPlot({
    d <- sim()
    dd <- d %>% transmute(hour, `haematoma (mL)` = V_hem,
                          `fibrinolysis` = Fib) %>% pivot_longer(-hour)
    ggplot(dd, aes(hour, value)) + geom_line(linewidth = 0.9, colour = "#566573") +
      facet_wrap(~name, scales = "free_y") + labs(x = "hours", y = NULL) + THEME
  })

  ## ---- 9. comparison ------------------------------------------------------
  ARMS <- list(
    none   = list(PROP_RATE = 3, EVD_OPEN = 0),
    t0     = list(PROP_RATE = 3, NE_RATE = 0.08),
    t1     = list(PROP_RATE = 3, NE_RATE = 0.08, EVD_OPEN = 1),
    t2     = list(PROP_RATE = 6, NE_RATE = 0.08, EVD_OPEN = 1, PACO2 = 33),
    t3     = list(PROP_RATE = 6, NE_RATE = 0.10, EVD_OPEN = 1, PACO2 = 33,
                  THIO_RATE = 4),
    crani  = list(PROP_RATE = 3, NE_RATE = 0.08, EVD_OPEN = 1, CRANI_T = 1440),
    noauto = list(PROP_RATE = 3, NE_RATE = 0.08, EVD_OPEN = 1, AUTOREG = 0))

  comparison <- eventReactive(input$run_cmp, {
    withProgress(message = "Running scenario comparison...", value = 0, {
      out <- lapply(input$cmp, function(a) {
        incProgress(1 / length(input$cmp))
        simulate_arm(input, ARMS[[a]]) %>% mutate(arm = a)
      })
      bind_rows(out)
    })
  })
  output$p_compare <- renderPlot({
    d <- comparison()
    dd <- d %>% select(hour, arm, ICP, CPP_o, PbtO2, F_core) %>%
      pivot_longer(c(-hour, -arm))
    ggplot(dd, aes(hour, value, colour = arm)) + geom_line(linewidth = 0.85) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "hours", y = NULL, colour = NULL) + THEME
  })
  output$t_compare <- renderDT({
    d <- comparison()
    d %>% group_by(arm) %>% group_modify(function(x, ...) {
      b <- burden(x)
      tibble(`peak ICP` = round(max(x$ICP), 1),
             `h ICP>22` = round(sum(x$ICP > 22) * diff(x$time)[1] / 60, 1),
             `ICP dose` = round(b$Dicp, 1),
             `CPP dose` = round(b$Dcpp, 1),
             `PbtO2 dose` = round(b$Dpbt, 1),
             `final core` = round(b$Fc, 4),
             `P(unfav)` = round(p_unfav(x, input$age), 3),
             `P(death)` = round(p_death(x, input$age), 3))
    }) %>% datatable(rownames = FALSE, options = list(dom = "t"))
  })

  ## ---- 10. biomarkers -----------------------------------------------------
  output$p_biomark <- renderPlot({
    d <- sim()
    dd <- d %>% transmute(hour, `S100B (t½ 1 h)` = S100B,
                          `UCH-L1 (t½ 7 h)` = UCHL1,
                          `GFAP (t½ 24 h)` = GFAP,
                          `NfL (t½ 21 d)` = NfL,
                          `kill rate (/min)` = kill) %>%
      pivot_longer(-hour)
    ggplot(dd, aes(hour, value)) + geom_line(linewidth = 0.9, colour = "#263238") +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "hours", y = NULL) + THEME
  })
  output$t_biomark <- renderTable({
    d <- sim()
    f <- function(v) d$hour[which.max(v)]
    data.frame(marker = c("S100B", "UCH-L1", "GFAP", "NfL"),
               `half-life` = c("~1 h", "~7 h", "~24 h", "~21 d"),
               `peak (h)` = c(f(d$S100B), f(d$UCHL1), f(d$GFAP), f(d$NfL)),
               reports = c("instantaneous rate", "recent rate",
                           "the integral so far", "cumulative axonal loss"),
               check.names = FALSE)
  }, digits = 1)

  ## ---- 11. outcome --------------------------------------------------------
  output$t_burden <- renderTable({
    b <- burden(sim())
    data.frame(quantity = c("final non-viable fraction",
                            "ICP dose (mmHg·h > 20)",
                            "CPP deficit dose (mmHg·h < 60)",
                            "PbtO2 deficit dose (mmHg·h < 20)",
                            "pupil score (herniation surrogate)"),
               value = c(b$Fc, b$Dicp, b$Dcpp, b$Dpbt, b$pupil))
  }, digits = 3)
  output$t_outcome <- renderTable({
    d <- sim()
    data.frame(endpoint = c("P(GOS-E 1–4) unfavourable at 6 months",
                            "P(death)"),
               probability = c(p_unfav(d, input$age), p_death(d, input$age)))
  }, digits = 3)
  output$p_dose <- renderPlot({
    d <- sim()
    dd <- d %>% transmute(hour, `cumulative ICP dose` = D_icp,
                          `cumulative CPP deficit dose` = D_cpp,
                          `non-viable fraction` = F_core,
                          `penumbra remaining` = F_pen) %>%
      pivot_longer(-hour)
    ggplot(dd, aes(hour, value)) + geom_line(linewidth = 0.9, colour = "#78909c") +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "hours", y = NULL) + THEME
  })
}

shinyApp(ui, server)
