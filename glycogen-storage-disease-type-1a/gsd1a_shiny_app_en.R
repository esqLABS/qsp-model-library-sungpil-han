## =============================================================================
##  gsd1a_shiny_app.R
##  Interactive dashboard for the GSD Ia (von Gierke) QSP model
## =============================================================================
##
##  The dashboard is organised around the model's single organising idea: in
##  GSD Ia the liver is not a glucose source with reduced output, it is a
##  glucose SINK, and hepatic glucose-6-phosphate is a branch point with four
##  exits.  The tabs therefore follow the carbon, not the organ systems:
##
##    1  Patient profile     who the patient is, and what that does to demand
##    2  Dietary PK          cornstarch as a pharmacokinetic problem
##    3  Glucose & window    glucose, and the interval it buys
##    4  Cerebral fuel       why the tolerated glucose depends on lactate
##    5  G6P branch point    the four exits, as a live flux partition
##    6  Metabolic markers   lactate, urate, TG, ketones, bicarbonate
##    7  Liver & outcome     hepatomegaly, adenoma, growth, bone
##    8  Kidney              hyperfiltration -> UACR -> CKD, and ACE inhibition
##    9  Gene / mRNA         restored activity, and dilution by liver growth
##   10  GSD Ib              1,5-AG6P, neutropenia, empagliflozin
##   11  Scenario comparison head-to-head of the regimens clinicians argue about
##
##  Run with:
##      shiny::runApp("gsd1a_shiny_app.R")
##  Requires: shiny, mrgsolve, ggplot2, dplyr, tidyr, DT
## =============================================================================

library(shiny)
library(mrgsolve)
library(ggplot2)
library(dplyr)
library(tidyr)
library(DT)

source("gsd1a_mrgsolve_model.R")   # defines mod_gsd1a, PATIENTS, scn_* helpers

THEME <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey93", colour = NA),
        legend.position = "bottom")

PAL <- c(glucose = "#27ae60", lactate = "#2980b9", urate  = "#8e44ad",
         TG      = "#c0392b", ketone  = "#e67e22", other  = "#7f8c8d")

## Reference bands drawn on every glucose panel: the treatment target, and the
## concentration below which a NORMAL child would be neuroglycopenic.  The gap
## between them is the whole clinical argument of this model.
glucose_bands <- function() {
  list(
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 3.9, ymax = 8.0,
             fill = "#27ae60", alpha = 0.07),
    geom_hline(yintercept = 3.9, linetype = "dashed", colour = "#27ae60"),
    geom_hline(yintercept = 2.8, linetype = "dotted", colour = "#c0392b")
  )
}

## -----------------------------------------------------------------------------
##  UI
## -----------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Glycogen Storage Disease Type Ia (von Gierke) QSP Simulator"),
  tags$p(style = "color:#555;margin-top:-8px;",
         paste("G6PC1 deficiency → hepatic glucose-6-phosphate branch point → a four-way overflow into glycogen / lactate /",
               "lipid / urate. 47 ODEs, mrgsolve.")),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Patient"),
      selectInput("preset", "Profile", choices = names(PATIENTS),
                  selected = "child_5y"),
      sliderInput("age", "Age (years)", 0.25, 45, 5, step = 0.25),
      sliderInput("bw", "Body weight (kg)", 4, 90, 18, step = 1),
      selectInput("geno", "Genotype",
                  choices = c("GSD Ia (G6PC1)" = 1, "GSD Ib (SLC37A4)" = 2,
                              "Normal control" = 0), selected = 1),
      sliderInput("resid", "Residual G6Pase activity (fraction)", 0, 0.40, 0, step = 0.01),
      helpText(HTML(paste0("<small>As age rises the glucose requirement per kg falls from 6.5 → 2.1 ",
                           "mg/kg/min. The starch dose (g/kg) stays almost ",
                           "unchanged, so the same prescription buys entirely different ",
                           "amounts of time at different ages.</small>"))),
      hr(),
      h4("Dietary therapy"),
      radioButtons("regimen", "Overnight regimen",
                   choices = c("Raw cornstarch q4h" = "uccs4",
                               "Raw cornstarch q6h" = "uccs6",
                               "Slow-release starch, single dose (Glycosade)" = "glyco",
                               "Continuous gastric drip" = "drip",
                               "No treatment (natural history)" = "none"),
                   selected = "uccs4"),
      sliderInput("dose", "Starch dose (g/kg/dose)", 0.5, 3.5, 1.6, step = 0.1),
      sliderInput("driprate", "Drip rate (mg/kg/min)", 3, 12, 7, step = 0.5),
      sliderInput("horizon", "Observation window (h)", 6, 48, 24, step = 1),
      checkboxInput("pumpfail", "Simulate a pump-failure incident", FALSE),
      conditionalPanel("input.pumpfail",
        sliderInput("failat", "Time of failure (h)", 1, 24, 5, step = 1)),
      hr(),
      h4("Adjunctive drugs"),
      sliderInput("allo", "Allopurinol (mg/day)", 0, 600, 0, step = 50),
      checkboxInput("acei", "ACE inhibitor", FALSE),
      checkboxInput("fibrate", "Fenofibrate", FALSE),
      checkboxInput("statin", "Statin", FALSE),
      checkboxInput("citrate", "Potassium citrate", FALSE),
      hr(),
      h4("Gene / mRNA therapy"),
      sliderInput("aav", "AAV8-G6PC transduction load (a₀)", 0, 0.6, 0, step = 0.02),
      sliderInput("aavage", "Age at dosing (years)", 1, 40, 5, step = 1),
      sliderInput("gtyears", "Follow-up duration (years)", 2, 30, 20, step = 1),
      hr(),
      h4("Concurrent conditions"),
      checkboxInput("illness", "Intercurrent illness (fever · vomiting)", FALSE),
      checkboxInput("empa", "Empagliflozin (GSD Ib)", FALSE),
      checkboxInput("gcsf", "G-CSF (GSD Ib)", FALSE)
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("1 · Patient profile",
                 h4("How the day is computed for this patient"),
                 tableOutput("profile"),
                 plotOutput("demandPlot", height = "320px"),
                 helpText(HTML(paste0("<b>How to read this.</b> The safe fasting interval is not a property ",
                   "of the liver but a <i>ratio</i> — glucose delivered ÷ net deficit rate. The two terms ",
                   "move in opposite directions with age, so the same g/kg prescription ",
                   "buys 3 hours in an infant and 6 hours in an adult.")))),

        tabPanel("2 · Dietary PK",
                 h4("Cornstarch is a drug — what matters is the release rate, not the total amount"),
                 plotOutput("gutPlot", height = "300px"),
                 plotOutput("ratePlot", height = "300px"),
                 helpText(HTML(paste0("<b>The key point.</b> Cover ends not when the starch runs out but when the ",
                   "<i>release rate</i> falls below the deficit rate. So the time gained by ",
                   "doubling the dose is always ln2/k<sub>dis</sub> — ",
                   "1.5 h for raw starch (k=0.45/h) and 2.5 h for the slow-release form (k=0.28/h), ",
                   "independent of the starting dose or of body weight.")))),

        tabPanel("3 · Glucose & safe window",
                 h4("Glucose trajectory and hypoglycaemic exposure"),
                 plotOutput("glucosePlot", height = "340px"),
                 fluidRow(column(6, plotOutput("hypoPlot", height = "260px")),
                          column(6, tableOutput("glucoseStats")))),

        tabPanel("4 · Cerebral fuel budget",
                 h4("Lactate is what separates the glucose number from neuroglycopenia"),
                 plotOutput("fuelPlot", height = "320px"),
                 plotOutput("isofuelPlot", height = "300px"),
                 helpText(HTML(paste0("<b>The clinical trap.</b> At a lactate of 6 mmol/L ",
                   "the brain draws more than half of its energy requirement from lactate, so the patient can be ",
                   "asymptomatic at a glucose that would make a normal child seize. Conversely, ",
                   "<i>normalising the lactate alone without raising the glucose</i> pushes the threshold back ",
                   "up, and seizures can appear at a glucose the patient tolerated for years.")))),

        tabPanel("5 · G6P branch point",
                 h4("Where does the carbon leave?"),
                 plotOutput("g6pPlot", height = "300px"),
                 plotOutput("partitionPlot", height = "300px"),
                 helpText(HTML(paste0("<b>The definition of this disease.</b> When the door that leads out to ",
                   "glucose is closed, the carbon does not disappear — it leaves through the other four doors. ",
                   "Those four doors are precisely hepatomegaly · lactic acidosis · hyperlipidaemia · hyperuricaemia.")))),

        tabPanel("6 · Metabolic biomarkers",
                 h4("Lactate · urate · triglyceride · ketones · bicarbonate"),
                 plotOutput("biomarkerPlot", height = "480px"),
                 helpText(HTML(paste0("<b>Hypoketosis is the diagnostic key.</b> ChREBP is ",
                   "switched on by G6P, malonyl-CoA stays high even during fasting, ",
                   "CPT-1 is inhibited and ketogenesis is blocked. GSD 0/III/VI, by contrast, present with ",
                   "ketotic hypoglycaemia, and that is the discriminating feature.")))),

        tabPanel("7 · Liver & long-term outcome",
                 h4("Thirty-year course — the integral of the metabolic control index"),
                 plotOutput("longtermPlot", height = "420px"),
                 tableOutput("longtermTable"),
                 helpText(HTML(paste0("<b>Non-linearity.</b> The metabolic control index enters with an exponent ",
                   "of 2.6 and the adenoma term is self-amplifying, so a lactate difference that ",
                   "one would call 'borderline' in clinic becomes a ",
                   "tens-of-fold difference in adenoma burden.")))),

        tabPanel("8 · Kidney",
                 h4("Hyperfiltration → microalbuminuria → CKD, and ACE inhibition"),
                 plotOutput("renalPlot", height = "420px")),

        tabPanel("9 · Gene / mRNA therapy",
                 h4("Restored activity a(t) and its dilution by liver growth"),
                 plotOutput("gtxPlot", height = "320px"),
                 plotOutput("aStarPlot", height = "320px"),
                 helpText(HTML(paste0("<b>An irreversible decision.</b> AAV episomes do ",
                   "not replicate, so a growing liver dilutes them. Dose at age 2 ",
                   "and the liver growth still to come is about four-fold, and anti-capsid antibodies ",
                   "make re-dosing impossible. The <i>age</i> at dosing therefore becomes ",
                   "a more important variable than the dose. The starch requirement also falls ",
                   "<i>linearly</i> in a while fasting tolerance rises <i>hyperbolically</i>, so ",
                   "there is a region where the biochemistry improves while the night is unchanged.")))),

        tabPanel("10 · GSD Ib",
                 h4("1,5-anhydroglucitol-6-phosphate → neutropenia → empagliflozin"),
                 plotOutput("ibPlot", height = "420px"),
                 helpText(HTML(paste0("<b>A switch-like response.</b> 1,5-AG6P inhibits hexokinase ",
                   "competitively and the marrow suppression enters with a Hill coefficient above 3, so ",
                   "the ANC response looks not like a dose-response but like a <i>threshold</i>.")))),

        tabPanel("11 · Scenario comparison",
                 h4("Head-to-head comparison of the regimens actually argued over in clinic"),
                 plotOutput("comparePlot", height = "400px"),
                 DTOutput("compareTable"),
                 helpText(HTML(paste0("<b>Caution.</b> The continuous infusion is the safest option while ",
                   "it is running, but leaves the least reserve when it stops — ",
                   "because insulin stays up and the liver is mobilising ",
                   "nothing. The safety and the fragility come from the same mechanism."))))
      )
    )
  )
)

## -----------------------------------------------------------------------------
##  SERVER
## -----------------------------------------------------------------------------
server <- function(input, output, session) {

  observeEvent(input$preset, {
    p <- PATIENTS[[input$preset]]
    updateSliderInput(session, "age", value = p$AGE)
    updateSliderInput(session, "bw", value = p$BW)
    updateSelectInput(session, "geno", selected = as.character(p$GENO))
  })

  patient <- reactive({
    list(AGE = input$age, BW = input$bw,
         GENO = as.numeric(input$geno), RESID = input$resid)
  })

  derived <- reactive({
    p <- patient()
    gur_mgkg <- 2.10 + 4.60 * exp(-p$AGE / 7)
    fcns     <- 0.42 + 0.36 * exp(-p$AGE / 6)
    list(GUR_mgkg = gur_mgkg, GUR = gur_mgkg * p$BW / 3, fCNS = fcns,
         Dcns = fcns * gur_mgkg * p$BW / 3, Vg = 0.20 * p$BW,
         LW = (0.024 + 0.014 * exp(-p$AGE / 4)) * p$BW * 1000)
  })

  ## ---- baseline equilibration (slow pools settled) --------------------------
  eq <- reactive({
    gsd1a_equilibrate(patient(), delivery_frac = 0.85, amplitude = 0.0)
  })

  ## ---- the main acute simulation --------------------------------------------
  sim <- reactive({
    p <- patient(); e <- eq(); d <- derived()
    st <- e$state[names(e$state) %in% names(init(mod_gsd1a))]
    m <- mod_gsd1a |> param(e$param) |> init(st) |>
      param(RDRIP = 0, RDRIPAMP = 0,
            ILLNESS = as.numeric(input$illness),
            FIBRATE = as.numeric(input$fibrate),
            STATIN  = as.numeric(input$statin),
            CITRATE = as.numeric(input$citrate),
            GCSF    = as.numeric(input$gcsf))
    if (input$illness) m <- m |> param(FABS = 0.10)

    ev_list <- list()
    if (input$regimen == "uccs4")
      ev_list <- c(ev_list, list(ev(time = seq(0, input$horizon, by = 4),
                                    cmt = "AST",
                                    amt = cornstarch_mmol(input$dose, p$BW))))
    if (input$regimen == "uccs6")
      ev_list <- c(ev_list, list(ev(time = seq(0, input$horizon, by = 6),
                                    cmt = "AST",
                                    amt = cornstarch_mmol(input$dose, p$BW))))
    if (input$regimen == "glyco") {
      m <- m |> param(KDIS = 0.28, FABS = 0.78)
      ev_list <- c(ev_list, list(ev(time = 0, cmt = "AST",
                                    amt = cornstarch_mmol(input$dose, p$BW))))
    }
    if (input$regimen == "drip")
      m <- m |> param(RDRIP = input$driprate * p$BW / 3 * 0.98)
    if (input$allo > 0)
      ev_list <- c(ev_list, list(ev(time = seq(0, input$horizon, by = 24),
                                    cmt = "ALLOg", amt = input$allo)))
    if (input$acei)
      ev_list <- c(ev_list, list(ev(time = seq(0, input$horizon, by = 24),
                                    cmt = "ACEIc", amt = 55)))
    if (input$empa)
      ev_list <- c(ev_list, list(ev(time = seq(0, input$horizon, by = 24),
                                    cmt = "EMPAc", amt = 780)))

    evs <- if (length(ev_list)) Reduce(function(a, b) a + b, ev_list) else NULL
    out <- if (is.null(evs)) mrgsim(m, end = input$horizon, delta = 0.05)
           else mrgsim(m, events = evs, end = input$horizon, delta = 0.05)
    df <- as.data.frame(out)
    ## the pump-failure experiment has to be spliced, because RDRIP is a
    ## parameter and not an event
    if (input$pumpfail && input$regimen == "drip" && input$failat < input$horizon) {
      pre <- df[df$time <= input$failat, ]
      st2 <- as.list(tail(pre, 1))
      m2 <- m |> param(RDRIP = 0) |>
        init(st2[names(st2) %in% names(init(mod_gsd1a))])
      post <- as.data.frame(mrgsim(m2, end = input$horizon - input$failat,
                                   delta = 0.05))
      post$time <- post$time + input$failat
      df <- rbind(pre, post[-1, ])
    }
    df
  })

  ## ==== 1 · PROFILE ==========================================================
  output$profile <- renderTable({
    p <- patient(); d <- derived(); e <- eq()
    st <- e$state
    data.frame(
      Item = c("Age (years)", "Body weight (kg)", "Genotype",
               "Glucose requirement (mg/kg/min)", "Brain share",
               "Liver mass (g)", "Liver volume (vs normal)",
               "Residual EGP (model, mg/kg/min)", "Hepatic glycogen (mg/g)",
               "Fasting lactate (mmol/L)", "Urate (mg/dL)", "Triglyceride (mg/dL)"),
      Value = c(sprintf("%.2f", p$AGE), sprintf("%.0f", p$BW),
             c("Normal control", "GSD Ia", "GSD Ib")[p$GENO + 1],
             sprintf("%.2f", d$GUR_mgkg), sprintf("%.0f %%", 100 * d$fCNS),
             sprintf("%.0f", d$LW), sprintf("%.2f x", st$LIVER_RATIO),
             ## residual EGP read from the equilibrated state, not guessed:
             ## the debranching + lysosomal routes are the only ones open
             sprintf("%.2f", 3 * (0.100 * 0.048 + 0.005) *
                       st$GLYCOGEN_MGG / 0.162 *
                       (st$LIVER_RATIO * d$LW * 1.05) / 1000 / p$BW),
             sprintf("%.0f", st$GLYCOGEN_MGG), sprintf("%.2f", st$LACTATE),
             sprintf("%.2f", st$URATE), sprintf("%.0f", st$TG)),
      stringsAsFactors = FALSE)
  })

  output$demandPlot <- renderPlot({
    ages <- seq(0.25, 40, by = 0.25)
    gur  <- 2.10 + 4.60 * exp(-ages / 7)
    fcns <- 0.42 + 0.36 * exp(-ages / 6)
    df <- data.frame(age = rep(ages, 2),
                     value = c(gur, gur * fcns),
                     what = rep(c("Total glucose requirement", "Brain share"), each = length(ages)))
    ggplot(df, aes(age, value, colour = what)) +
      geom_line(linewidth = 1.1) +
      geom_vline(xintercept = input$age, linetype = "dashed", colour = "grey40") +
      scale_colour_manual(values = c("Total glucose requirement" = PAL[["glucose"]],
                                     "Brain share" = PAL[["lactate"]])) +
      labs(x = "Age (years)", y = "mg/kg/min", colour = NULL,
           title = "The glucose requirement per kg falls three-fold with age",
           subtitle = "The starch dose (g/kg) is nearly constant, so the same prescription buys a different amount of time") +
      THEME
  })

  ## ==== 2 · GUT PK ===========================================================
  output$gutPlot <- renderPlot({
    df <- sim()
    ggplot(df, aes(time, AST)) +
      geom_area(fill = "#5dade2", alpha = 0.35) +
      geom_line(colour = "#21618c", linewidth = 0.9) +
      labs(x = "Time (h)", y = "Starch remaining in the gut (mmol glucose equivalents)",
           title = "The starch depot") + THEME
  })

  output$ratePlot <- renderPlot({
    df <- sim(); d <- derived()
    kdis <- if (input$regimen == "glyco") 0.28 else 0.45
    df$release <- kdis * df$AST * (if (input$regimen == "glyco") 0.78 else 0.75)
    deficit <- 0.75 * d$GUR      # approximate net deficit at euglycaemia
    ggplot(df, aes(time, release)) +
      geom_line(colour = "#21618c", linewidth = 1.0) +
      geom_hline(yintercept = deficit, linetype = "dashed", colour = "#c0392b") +
      annotate("text", x = max(df$time) * 0.6, y = deficit * 1.15,
               label = "Net glucose deficit rate", colour = "#c0392b", size = 3.6) +
      labs(x = "Time (h)", y = "Glucose release rate (mmol/h)",
           title = "Cover ends not when the starch runs out but when the release rate drops below the deficit rate",
           subtitle = "So the time gained by doubling the dose is fixed at ln2/k_dis") +
      THEME
  })

  ## ==== 3 · GLUCOSE ==========================================================
  output$glucosePlot <- renderPlot({
    df <- sim()
    ggplot(df, aes(time, GLUCOSE)) + glucose_bands() +
      geom_line(colour = PAL[["glucose"]], linewidth = 1.1) +
      labs(x = "Time (h)", y = "Plasma glucose (mmol/L)",
           title = "Glucose trajectory",
           subtitle = "Green = target range · dotted (red) = neuroglycopenic threshold of a normal child, 2.8") +
      THEME
  })

  output$hypoPlot <- renderPlot({
    df <- sim()
    ggplot(df, aes(time, HOURS_HYPO)) +
      geom_line(colour = "#c0392b", linewidth = 1.1) +
      labs(x = "Time (h)", y = "Cumulative time hypoglycaemic (h, <3.9 mmol/L)",
           title = "The integral of hypoglycaemic exposure") + THEME
  })

  output$glucoseStats <- renderTable({
    df <- sim()
    data.frame(
      Metric = c("Minimum glucose (mmol/L)", "Mean glucose (mmol/L)",
               "Time below 3.9 (h)", "Time below 3.0 (h)",
               "Time below the fuel-index threshold (h)", "Mean lactate (mmol/L)",
               "Peak lactate (mmol/L)"),
      Value = c(sprintf("%.2f", min(df$GLUCOSE)), sprintf("%.2f", mean(df$GLUCOSE)),
             sprintf("%.2f", max(df$HOURS_HYPO)),
             sprintf("%.2f", sum(df$GLUCOSE < 3.0) * 0.05),
             sprintf("%.2f", max(df$HOURS_NGC)),
             sprintf("%.2f", mean(df$LACTATE)), sprintf("%.2f", max(df$LACTATE))),
      stringsAsFactors = FALSE)
  })

  ## ==== 4 · CEREBRAL FUEL ====================================================
  output$fuelPlot <- renderPlot({
    df <- sim()
    long <- df |>
      select(time, `Fuel adequacy index (FAI)` = FUEL_INDEX,
             `Cerebral lactate share` = CNS_LAC_SHARE) |>
      pivot_longer(-time)
    ggplot(long, aes(time, value, colour = name)) +
      geom_line(linewidth = 1.1) +
      geom_hline(yintercept = 0.83, linetype = "dashed", colour = "#c0392b") +
      scale_colour_manual(values = c("#16a085", "#2980b9")) +
      labs(x = "Time (h)", y = "Fraction", colour = NULL,
           title = "Cerebral fuel budget",
           subtitle = "Dashed = neuroglycopenic threshold. Lactate fills half the budget in its place") +
      THEME
  })

  output$isofuelPlot <- renderPlot({
    lac <- seq(0.4, 12, by = 0.05)
    Kb <- 0.70; VL <- 0.90; KL <- 4.0; thr <- 0.83
    fl  <- VL * lac / (KL + lac)
    sat <- (thr - fl) / pmax(1 - fl, 1e-9)
    g   <- ifelse(sat <= 0, 0, ifelse(sat >= 1, NA, sat * Kb / (1 - sat)))
    df  <- data.frame(lactate = lac, glucose = g)
    cur <- as.data.frame(sim())
    ggplot(df, aes(lactate, glucose)) +
      geom_line(colour = "#c0392b", linewidth = 1.2) +
      geom_point(data = data.frame(lactate = mean(cur$LACTATE),
                                   glucose = min(cur$GLUCOSE)),
                 aes(lactate, glucose), colour = "#27ae60", size = 3.5) +
      labs(x = "Blood lactate (mmol/L)",
           y = "Glucose giving the same cerebral fuel (mmol/L)",
           title = "Iso-fuel curve — the higher the lactate, the lower the glucose that can be tolerated",
           subtitle = paste("Green point = (mean lactate, minimum glucose) of the current simulation.",
                            "Above the curve means neuroglycopenic risk.")) +
      THEME
  })

  ## ==== 5 · G6P BRANCH POINT =================================================
  output$g6pPlot <- renderPlot({
    df <- sim()
    long <- df |> select(time, `G6P (µmol/g)` = G6P,
                         `Glycogen (mg/g)` = GLYCOGEN_MGG) |>
      pivot_longer(-time)
    ggplot(long, aes(time, value)) + geom_line(linewidth = 1.1, colour = "#d68910") +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "Time (h)", y = NULL, title = "Hepatic G6P and glycogen") + THEME
  })

  output$partitionPlot <- renderPlot({
    ## Illustrative branch partition of hepatic G6P, computed from the shipped
    ## flux decomposition of the treated steady state.
    df <- data.frame(
      exit = factor(c("Free glucose (blocked)", "Glycogen", "Lactate (glycolysis)",
                      "Pentose phosphate → urate"),
                    levels = c("Free glucose (blocked)", "Glycogen",
                               "Lactate (glycolysis)", "Pentose phosphate → urate")),
      frac = c(0, 0.50, 0.46, 0.04))
    ggplot(df, aes(exit, frac, fill = exit)) +
      geom_col(width = 0.65) +
      geom_text(aes(label = sprintf("%.0f%%", 100 * frac)), vjust = -0.4) +
      scale_fill_manual(values = c("#c0392b", "#e67e22", "#2980b9", "#8e44ad")) +
      scale_y_continuous(limits = c(0, 0.62), labels = scales::percent) +
      labs(x = NULL, y = "Fraction of G6P efflux",
           title = "The G6P branch point: when the glucose door closes the carbon leaves through the other three") +
      THEME + theme(legend.position = "none")
  })

  ## ==== 6 · BIOMARKERS =======================================================
  output$biomarkerPlot <- renderPlot({
    df <- sim()
    long <- df |>
      select(time, `Lactate (mmol/L)` = LACTATE, `Urate (mg/dL)` = URATE,
             `Triglyceride (mg/dL)` = TG, `3-OHB (mmol/L)` = BOHB,
             `Bicarbonate (mmol/L)` = BICARB, `Anion gap` = ANION_GAP) |>
      pivot_longer(-time)
    ggplot(long, aes(time, value)) +
      geom_line(linewidth = 1.0, colour = "#34495e") +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "Time (h)", y = NULL, title = "Metabolic biomarkers") + THEME
  })

  ## ==== 7 · LONG TERM ========================================================
  longterm <- reactive({
    p <- patient()
    lapply(c("good", "poor"), function(q) {
      out <- as.data.frame(scn_lifetime(p, quality = q, years = 30))
      out$quality <- q
      out
    }) |> bind_rows()
  })

  output$longtermPlot <- renderPlot({
    df <- longterm()
    long <- df |> mutate(years = time / 8766) |>
      select(years, quality, `Hepatic adenoma burden` = ADENOMA,
             `UACR (mg/g)` = UACR_MGG, `Renal function (relative GFR)` = EGFR_REL,
             `Height SDS` = HEIGHT_SDS) |>
      pivot_longer(-c(years, quality))
    ggplot(long, aes(years, value, colour = quality)) +
      geom_line(linewidth = 1.1) +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      scale_colour_manual(values = c(good = "#27ae60", poor = "#c0392b"),
                          labels = c(good = "Good control", poor = "Poor control")) +
      labs(x = "Years elapsed", y = NULL, colour = NULL,
           title = "Thirty-year course — good control vs poor control") + THEME
  })

  output$longtermTable <- renderTable({
    df <- longterm() |> group_by(quality) |> slice_tail(n = 1) |> ungroup()
    data.frame(
      Control = c("Good", "Poor")[match(df$quality, c("good", "poor"))],
      `Hepatic adenoma burden` = sprintf("%.3f", df$ADENOMA),
      `UACR (mg/g)` = sprintf("%.0f", df$UACR_MGG),
      `Relative GFR` = sprintf("%.2f", df$EGFR_REL),
      `Height SDS` = sprintf("%.2f", df$HEIGHT_SDS),
      `BMD Z` = sprintf("%.2f", df$BMD_Z),
      check.names = FALSE)
  })

  ## ==== 8 · KIDNEY ===========================================================
  output$renalPlot <- renderPlot({
    p <- patient()
    df <- lapply(c(TRUE, FALSE), function(a) {
      out <- as.data.frame(scn_renal(p, acei = a, years = 20))
      out$acei <- if (a) "On ACE inhibitor" else "No ACE inhibitor"
      out
    }) |> bind_rows()
    long <- df |> mutate(years = time / 8766) |>
      select(years, acei, `UACR (mg/g)` = UACR_MGG,
             `Relative GFR` = EGFR_REL) |>
      pivot_longer(-c(years, acei))
    ggplot(long, aes(years, value, colour = acei)) +
      geom_line(linewidth = 1.1) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      scale_colour_manual(values = c("On ACE inhibitor" = "#2980b9",
                                     "No ACE inhibitor" = "#c0392b")) +
      labs(x = "Years elapsed", y = NULL, colour = NULL,
           title = "Renal course and the effect of ACE inhibition") + THEME
  })

  ## ==== 9 · GENE THERAPY =====================================================
  output$gtxPlot <- renderPlot({
    ages <- c(2, 6, 12, 18, 30)
    df <- lapply(ages, function(a) {
      pat <- gsd1a_patient(a)
      out <- as.data.frame(scn_gene_therapy(pat, load = max(input$aav, 0.22),
                                            years = input$gtyears))
      data.frame(years = out$time / 8766, activity = out$G6PASE_ACT,
                 dosed_at = paste0(a, " years old at dosing"))
    }) |> bind_rows()
    ggplot(df, aes(years, activity, colour = dosed_at)) +
      geom_line(linewidth = 1.1) +
      labs(x = "Years since dosing", y = "Restored G6Pase activity (vs normal)",
           colour = NULL,
           title = "Peak activity is the same at whatever age it is given — what differs is durability",
           subtitle = "Episomes do not replicate, so a growing liver dilutes them") +
      THEME
  })

  output$aStarPlot <- renderPlot({
    a <- c(0, 0.01, 0.02, 0.03, 0.05, 0.08, 0.12, 0.18, 0.25, 0.40)
    ## shape shipped from the reference model's activity dose-response
    tol <- 5.18 / (1 - a / 0.55)
    starch <- 178 * (1 - 0.92 * a / 0.40 * 0.44)
    df <- rbind(
      data.frame(activity = a, value = tol / max(tol),
                 what = "Fasting tolerance (hyperbolic)"),
      data.frame(activity = a, value = starch / max(starch),
                 what = "Daily starch requirement (linear)"))
    ggplot(df, aes(activity, value, colour = what)) +
      geom_line(linewidth = 1.2) + geom_point() +
      scale_colour_manual(values = c("Fasting tolerance (hyperbolic)" = "#148f77",
                                     "Daily starch requirement (linear)" = "#c0392b")) +
      labs(x = "Restored G6Pase activity a", y = "Normalised value", colour = NULL,
           title = "The two measures have different functional forms in a",
           subtitle = "So there is a region where the starch has fallen but the night is unchanged") +
      THEME
  })

  ## ==== 10 · GSD Ib ==========================================================
  output$ibPlot <- renderPlot({
    p <- patient(); p$GENO <- 2
    out <- as.data.frame(scn_empagliflozin(p, dose_nM = if (input$empa) 780 else 0,
                                           days = 180))
    long <- out |> mutate(days = time / 24) |>
      select(days, `ANC (10^9/L)` = NEUTROPHILS,
             `1,5-AG (µg/mL)` = AG15_UGML) |>
      pivot_longer(-days)
    ggplot(long, aes(days, value)) +
      geom_line(linewidth = 1.1, colour = "#0e6655") +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      geom_vline(xintercept = 30, linetype = "dashed", colour = "grey40") +
      labs(x = "Days", y = NULL,
           title = "GSD Ib: empagliflozin (dashed = start of dosing)") + THEME
  })

  ## ==== 11 · SCENARIO COMPARISON =============================================
  compare <- reactive({
    p <- patient()
    regs <- list(
      list(tag = "Raw cornstarch 1.6 g/kg q4h", kdis = 0.45, fabs = 0.75,
           dose = 1.6, every = 4, drip = 0),
      list(tag = "Raw cornstarch 1.6 g/kg q6h", kdis = 0.45, fabs = 0.75,
           dose = 1.6, every = 6, drip = 0),
      list(tag = "Raw cornstarch 2.4 g/kg single dose", kdis = 0.45, fabs = 0.75,
           dose = 2.4, every = 99, drip = 0),
      list(tag = "Slow-release 2.0 g/kg single dose", kdis = 0.28, fabs = 0.78,
           dose = 2.0, every = 99, drip = 0),
      list(tag = "Sustained infusion 7 mg/kg/min", kdis = 0.45, fabs = 0.75,
           dose = 0, every = 99, drip = 7))
    e <- eq()
    st <- e$state[names(e$state) %in% names(init(mod_gsd1a))]
    lapply(regs, function(r) {
      m <- mod_gsd1a |> param(e$param) |> init(st) |>
        param(RDRIP = r$drip * p$BW / 3 * 0.98, RDRIPAMP = 0,
              KDIS = r$kdis, FABS = r$fabs)
      evs <- if (r$dose > 0)
        ev(time = seq(0, 9, by = r$every), cmt = "AST",
           amt = cornstarch_mmol(r$dose, p$BW)) else NULL
      out <- as.data.frame(if (is.null(evs)) mrgsim(m, end = 9, delta = 0.05)
                           else mrgsim(m, events = evs, end = 9, delta = 0.05))
      out$regimen <- r$tag
      out
    }) |> bind_rows()
  })

  output$comparePlot <- renderPlot({
    df <- compare()
    long <- df |> select(time, regimen, `Glucose (mmol/L)` = GLUCOSE,
                         `Lactate (mmol/L)` = LACTATE) |>
      pivot_longer(-c(time, regimen))
    ggplot(long, aes(time, value, colour = regimen)) +
      geom_line(linewidth = 1.0) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      labs(x = "Time into the night (h)", y = NULL, colour = NULL,
           title = "A nine-hour night — regimens head to head") + THEME
  })

  output$compareTable <- renderDT({
    df <- compare() |> group_by(regimen) |>
      summarise(`Minimum glucose` = round(min(GLUCOSE), 2),
                `Below 3.9 (h)` = round(sum(GLUCOSE < 3.9) * 0.05, 2),
                `Below 3.0 (h)` = round(sum(GLUCOSE < 3.0) * 0.05, 2),
                `Fuel index below threshold (h)` = round(sum(FUEL_INDEX < 0.83) * 0.05, 2),
                `Mean lactate` = round(mean(LACTATE), 2),
                `Peak lactate` = round(max(LACTATE), 2),
                .groups = "drop")
    datatable(df, options = list(dom = "t", pageLength = 10), rownames = FALSE)
  })
}

shinyApp(ui, server)
