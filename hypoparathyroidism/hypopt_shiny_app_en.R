## ============================================================================
##  Chronic Hypoparathyroidism (HypoPT) — QSP Interactive Dashboard
##  ============================================================================
##  Front-end for hypopt_mrgsolve_model_en.R.
##
##  Run:
##      shiny::runApp("hypopt_shiny_app.R")
##  or from this directory:
##      R -e 'shiny::runApp(".", launch.browser = TRUE)'
##
##  Nine tabs:
##      1. Patient Profile   — aetiology, chief-cell mass, CaSR set-point, labs
##      2. PK              — PTH molecules, calcitriol, thiazide, encaleret
##      3. Calcium · Phosphate Homeostasis  — the homeostatic core (Ca, iCa, Pi, Mg, 1,25D, FGF23)
##      4. Renal Safety     — urine Ca, FE_Ca, TmP/GFR, nephrocalcinosis, eGFR
##      5. Bone              — turnover markers, osteoblast/osteoclast, BMD Z
##      6. Clinical Endpoints — symptoms, QTc, QoL, triple composite response
##      7. Scenario Comparison   — the 15 prebuilt scenarios side by side
##      8. Biomarker Summary — one table of every steady-state readout
##      9. References · Help     — model provenance and how to read it
##
##  DISCLAIMER: educational / research tool. Not for clinical use.
## ============================================================================

suppressPackageStartupMessages({
  library(shiny)
  library(mrgsolve)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(DT)
})

source("hypopt_mrgsolve_model_en.R", local = TRUE, chdir = TRUE)

MOD <- HYPOPT_build()

## ---------------------------------------------------------------------------
##  Plot furniture
## ---------------------------------------------------------------------------
theme_hp <- function() {
  theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          strip.text  = element_text(face = "bold"),
          legend.position = "bottom",
          plot.title  = element_text(face = "bold"))
}

PAL <- c("#2471A3", "#C0392B", "#1E8449", "#B9770E", "#7D3C98",
         "#117A65", "#CA6F1E", "#34495E", "#CB4335", "#5D6D7E")

## reference bands drawn on the calcium panels
band <- function(ymin, ymax, fill = "#1E8449", alpha = 0.08) {
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax,
           fill = fill, alpha = alpha)
}

lineplot <- function(d, y, ylab, title, xdays = TRUE, bands = NULL) {
  d$xx <- if (xdays) d$time/24 else d$time
  p <- ggplot(d, aes(x = xx, y = .data[[y]], colour = scenario))
  if (!is.null(bands)) p <- p + band(bands[1], bands[2])
  p + geom_line(linewidth = 0.7) +
    scale_colour_manual(values = rep(PAL, 4)) +
    labs(x = if (xdays) "Time (days)" else "Time (hours)", y = ylab,
         title = title, colour = NULL) +
    theme_hp()
}

## ---------------------------------------------------------------------------
##  UI
## ---------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Chronic Hypoparathyroidism — QSP Simulator"),
  tags$p(style = "color:#7F8C8D;margin-top:-8px;",
         "Educational/research quantitative systems pharmacology model. Not for use in clinical decision-making."),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      h4("Patient · Disease"),
      selectInput("aetiology", "Aetiology",
                  choices = c("Post-surgical"        = "surg",
                              "Autoimmune"          = "ai",
                              "ADH1 (CASR gain-of-function)"           = "adh1",
                              "Hypomagnesaemia (functional)"        = "mglow",
                              "Normal control"                    = "healthy"),
                  selected = "surg"),
      sliderInput("ptmass", "Residual chief-cell function (0-1)", 0, 1, 0.05, step = 0.01),
      sliderInput("setpt", "CaSR set-point (mmol/L)", 0.70, 1.40, 1.20, step = 0.01),
      sliderInput("alb", "Serum albumin (g/dL)", 2.0, 5.0, 4.0, step = 0.1),
      sliderInput("ph", "Arterial pH", 7.25, 7.60, 7.40, step = 0.01),
      sliderInput("gfr0", "Baseline eGFR (mL/min/1.73m²)", 20, 130, 120, step = 5),
      sliderInput("dietca", "Dietary calcium (mg/day)", 300, 1600, 1000, step = 50),
      sliderInput("dietpi", "Dietary phosphate (mg/day)", 400, 2000, 1200, step = 50),
      sliderInput("dietmg", "Dietary magnesium (mg/day)", 20, 500, 300, step = 10),

      hr(),
      h4("Conventional Treatment"),
      sliderInput("casup", "Oral calcium single dose (mg elemental calcium)", 0, 1500, 500, step = 50),
      selectInput("caii", "Calcium dosing interval", c("Once daily" = 24, "Twice daily" = 12,
                                              "Three times daily" = 8, "Four times daily" = 6),
                  selected = 8),
      sliderInput("ctr", "Calcitriol single dose (µg)", 0, 2, 0.25, step = 0.125),
      selectInput("ctrii", "Calcitriol interval", c("Once daily" = 24, "Twice daily" = 12),
                  selected = 12),
      sliderInput("d3", "Cholecalciferol (IU/day)", 0, 4000, 1000, step = 200),
      checkboxInput("thz", "Thiazide (HCTZ 25 mg BID)", FALSE),
      checkboxInput("binder", "Phosphate binder / low-phosphate diet", FALSE),
      sliderInput("mgsup", "Magnesium supplementation (mg/day)", 0, 1200, 0, step = 50),

      hr(),
      h4("PTH Replacement Therapy"),
      selectInput("pthdrug", "Drug",
                  choices = c("None"                                = "none",
                              "rhPTH(1-84) SC QD"                   = "rhpth84",
                              "Teriparatide PTH(1-34) SC"          = "terip",
                              "Palopegteriparatide (TransCon PTH)" = "transcon",
                              "Eneboparatide (RG-selective)"          = "eneb"),
                  selected = "none"),
      sliderInput("pthdose", "Single dose (µg)", 0, 120, 18, step = 2),
      selectInput("pthii", "Dosing interval", c("Once daily" = 24, "Twice daily" = 12,
                                          "Three times daily" = 8), selected = 24),

      hr(),
      h4("CaSR-Targeted Therapy"),
      sliderInput("enc", "Encaleret single dose (mg, BID)", 0, 120, 0, step = 10),

      hr(),
      sliderInput("days", "Simulation duration (days)", 14, 3650, 180, step = 14),
      actionButton("go", "Run simulation", class = "btn-primary btn-block")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",

        tabPanel("1. Patient Profile",
                 br(),
                 fluidRow(
                   column(6, h4("Baseline (pre-treatment) steady state"), tableOutput("tbl_base")),
                   column(6, h4("Post-treatment steady state"),        tableOutput("tbl_treat"))),
                 hr(),
                 h4("Calcium Set-Point Curve — This Patient's Ca–PTH Relationship"),
                 plotOutput("p_setpoint", height = "330px"),
                 helpText("The dashed line is the normal set-point (1.20 mmol/L); the solid line is the curve reflecting the current",
                          "set-point and residual chief-cell function.",
                          "In ADH1, the curve shifts leftward, so PTH looks",
                          "'inappropriately normal' even at low calcium.")),

        tabPanel("2. PK",
                 br(),
                 fluidRow(
                   column(6, plotOutput("p_pth", height = "300px")),
                   column(6, plotOutput("p_pthzoom", height = "300px"))),
                 fluidRow(
                   column(6, plotOutput("p_d125pk", height = "300px")),
                   column(6, plotOutput("p_drugpk", height = "300px"))),
                 helpText("The graph at top right is a zoom of the last 7 days.",
                          "The biphasic absorption of rhPTH(1-84) produces a spiky curve, while",
                          "the sustained release of TransCon PTH produces an almost flat curve.")),

        tabPanel("3. Calcium · Phosphate Homeostasis",
                 br(),
                 fluidRow(
                   column(6, plotOutput("p_ca",   height = "290px")),
                   column(6, plotOutput("p_ica",  height = "290px"))),
                 fluidRow(
                   column(6, plotOutput("p_pi",   height = "290px")),
                   column(6, plotOutput("p_mg",   height = "290px"))),
                 fluidRow(
                   column(6, plotOutput("p_d125", height = "290px")),
                   column(6, plotOutput("p_fgf",  height = "290px")))),

        tabPanel("4. Renal Safety",
                 br(),
                 fluidRow(
                   column(6, plotOutput("p_uca",  height = "290px")),
                   column(6, plotOutput("p_feca", height = "290px"))),
                 fluidRow(
                   column(6, plotOutput("p_tmp",  height = "290px")),
                   column(6, plotOutput("p_caxp", height = "290px"))),
                 fluidRow(
                   column(6, plotOutput("p_nc",   height = "290px")),
                   column(6, plotOutput("p_egfr", height = "290px"))),
                 helpText("This tab is the core of this model. Conventional treatment brings serum calcium",
                          "into target but fails to lower urinary calcium.",
                          "Only PTH replacement therapy lowers urinary calcium at the same serum calcium.")),

        tabPanel("5. Bone",
                 br(),
                 fluidRow(
                   column(6, plotOutput("p_turn", height = "290px")),
                   column(6, plotOutput("p_cells", height = "290px"))),
                 fluidRow(
                   column(6, plotOutput("p_markers", height = "290px")),
                   column(6, plotOutput("p_bmd",     height = "290px")))),

        tabPanel("6. Clinical Endpoints",
                 br(),
                 fluidRow(
                   column(6, plotOutput("p_sym", height = "290px")),
                   column(6, plotOutput("p_qtc", height = "290px"))),
                 fluidRow(
                   column(6, plotOutput("p_qol", height = "290px")),
                   column(6, plotOutput("p_bg",  height = "290px"))),
                 hr(),
                 h4("Triple Composite Response (PaTHway/PARALLAX Definition)"),
                 verbatimTextOutput("txt_resp3"),
                 helpText("Albumin-corrected calcium 8.0–10.6 mg/dL AND oral calcium ≤600 mg/day",
                          "AND discontinuation of active vitamin D — all three conditions must be met.")),

        tabPanel("7. Scenario Comparison",
                 br(),
                 checkboxGroupInput(
                   "scn", "Scenarios to compare",
                   choices = c(
                     "01 Normal control" = 1,  "02 No treatment" = 2,
                     "03 Conventional treatment" = 3,    "04 Conventional+thiazide" = 4,
                     "05 Overtreatment" = 5,       "06 rhPTH(1-84)" = 6,
                     "07 Teriparatide" = 7, "08 TransCon 18 µg" = 8,
                     "09 TransCon 30 µg" = 9, "10 ADH1 no treatment" = 10,
                     "11 ADH1+encaleret" = 11, "12 Low Mg → supplementation" = 12,
                     "13 Acute crisis (IV Ca)" = 13,
                     "14 10-year conventional treatment" = 14, "15 10-year TransCon" = 15),
                   selected = c(2, 3, 8), inline = TRUE),
                 actionButton("goscn", "Run scenarios", class = "btn-success"),
                 br(), br(),
                 fluidRow(
                   column(6, plotOutput("s_ca",  height = "300px")),
                   column(6, plotOutput("s_uca", height = "300px"))),
                 fluidRow(
                   column(6, plotOutput("s_pi",   height = "300px")),
                   column(6, plotOutput("s_egfr", height = "300px"))),
                 hr(),
                 h4("Scenario Summary"),
                 DT::dataTableOutput("s_tbl"),
                 hr(),
                 h4("Serum Calcium Swing Within the Dosing Interval (Ca swing)"),
                 DT::dataTableOutput("s_swing")),

        tabPanel("8. Biomarker Summary",
                 br(),
                 h4("Steady-State Biomarkers at Current Settings"),
                 DT::dataTableOutput("b_tbl"),
                 hr(),
                 h4("Conventional-Treatment Dose Grid — Renal Cost of Reaching Target Calcium"),
                 helpText("For each calcium/calcitriol combination, the serum calcium reached and, in exchange,",
                          "the resulting 24-hour urinary calcium is calculated (120-day steady state)."),
                 actionButton("gogrid", "Calculate grid (takes tens of seconds)"),
                 br(), br(),
                 DT::dataTableOutput("b_grid")),

        tabPanel("9. References · Help",
                 br(),
                 uiOutput("refs"))
      )
    )
  )
)

## ---------------------------------------------------------------------------
##  Server
## ---------------------------------------------------------------------------
server <- function(input, output, session) {

  ## keep the aetiology preset and the sliders in sync
  observeEvent(input$aetiology, {
    p <- switch(input$aetiology,
                surg    = list(pt = 0.05, sp = 1.20, mg = 300),
                ai      = list(pt = 0.10, sp = 1.20, mg = 300),
                adh1    = list(pt = 1.00, sp = 0.85, mg = 300),
                mglow   = list(pt = 1.00, sp = 1.20, mg = 40),
                healthy = list(pt = 1.00, sp = 1.20, mg = 300))
    updateSliderInput(session, "ptmass", value = p$pt)
    updateSliderInput(session, "setpt",  value = p$sp)
    updateSliderInput(session, "dietmg", value = p$mg)
  })

  ## ---- parameter block for the current UI state --------------------------
  base_pars <- reactive({
    list(PTMASS = input$ptmass, SETPT = input$setpt,
         ALB = input$alb, PHART = input$ph, GFRC0 = input$gfr0,
         DIETCA = input$dietca, DIETPI = input$dietpi, DIETMG = input$dietmg,
         MGSUP = 0, BINDER = 0)
  })

  treat_pars <- reactive({
    p <- base_pars()
    p$MGSUP  <- input$mgsup
    p$BINDER <- as.numeric(input$binder)
    p$DCASUP <- input$casup*(24/as.numeric(input$caii))
    p$DCTR   <- input$ctr*(24/as.numeric(input$ctrii))
    if (input$pthdrug != "none")
      p <- utils::modifyList(p, HYPOPT_pth_params[[input$pthdrug]])
    p
  })

  events <- reactive({
    d  <- input$days
    ev <- NULL
    add <- function(e) if (is.null(ev)) e else c(ev, e)
    if (input$casup > 0)
      ev <- add(hypopt_calcium(input$casup, ii = as.numeric(input$caii), days = d))
    if (input$ctr > 0)
      ev <- add(hypopt_calcitriol(input$ctr, ii = as.numeric(input$ctrii), days = d))
    if (input$d3 > 0)
      ev <- add(hypopt_vitd3(input$d3, ii = 24, days = d))
    if (isTRUE(input$thz))
      ev <- add(hypopt_thiazide(25, ii = 12, days = d))
    if (input$enc > 0)
      ev <- add(hypopt_encaleret(input$enc, ii = 12, days = d))
    if (input$pthdrug == "transcon" && input$pthdose > 0)
      ev <- add(hypopt_transcon(input$pthdose, ii = as.numeric(input$pthii), days = d))
    if (input$pthdrug %in% c("rhpth84", "terip", "eneb") && input$pthdose > 0)
      ev <- add(hypopt_pth_sc(input$pthdose, ii = as.numeric(input$pthii), days = d))
    ev
  })

  ## ---- the two runs: untreated baseline, then therapy --------------------
  baseline_state <- eventReactive(input$go, {
    withProgress(message = "Calculating disease steady state...", value = 0.3, {
      hypopt_baseline(MOD, base_pars(), days = 400)
    })
  }, ignoreNULL = FALSE)

  sim <- eventReactive(input$go, {
    withProgress(message = "Running simulation...", value = 0.6, {
      init <- hypopt_baseline(MOD, base_pars(), days = 400)
      m    <- update(MOD, param = treat_pars(), init = as.list(init))
      d    <- input$days
      delta <- if (d <= 60) 0.5 else if (d <= 400) 2 else 24*7
      ev   <- events()
      out  <- if (is.null(ev)) {
        mrgsim_df(m, end = d*24, delta = delta, hmax = 2)
      } else {
        mrgsim_df(m, events = ev, end = d*24, delta = delta, hmax = 2)
      }
      out$scenario <- "Current settings"
      out
    })
  }, ignoreNULL = FALSE)

  ss <- function(d, col, frac = 0.1) {
    n <- nrow(d); k <- max(1, floor(n*(1 - frac)))
    mean(d[[col]][k:n], na.rm = TRUE)
  }

  ## ---- Tab 1 -------------------------------------------------------------
  output$tbl_base <- renderTable({
    iv <- baseline_state()
    m  <- update(MOD, param = base_pars(), init = as.list(iv))
    d  <- mrgsim_df(m, end = 48, delta = 2)
    data.frame(
      Item = c("Albumin-corrected calcium (mg/dL)", "Ionised calcium (mmol/L)",
               "Serum phosphate (mg/dL)", "PTH (pg/mL)", "1,25(OH)₂D (pg/mL)",
               "24h urinary calcium (mg)", "Bone turnover index", "BMD Z-score"),
      Value = round(c(ss(d, "CACORRD"), ss(d, "CAIONMM"), ss(d, "PISER"),
                   ss(d, "PTHPG"), ss(d, "D125PG"), ss(d, "UCA24"),
                   ss(d, "TURNOV"), ss(d, "BMDZS")), 2))
  }, striped = TRUE)

  output$tbl_treat <- renderTable({
    d <- sim()
    data.frame(
      Item = c("Albumin-corrected calcium (mg/dL)", "Ionised calcium (mmol/L)",
               "Serum phosphate (mg/dL)", "PTH signal (pg/mL eq)", "1,25(OH)₂D (pg/mL)",
               "24h urinary calcium (mg)", "Bone turnover index", "BMD Z-score",
               "eGFR (mL/min/1.73m²)", "Symptom score (0-100)"),
      Value = round(c(ss(d, "CACORRD"), ss(d, "CAIONMM"), ss(d, "PISER"),
                   ss(d, "PTHPG"), ss(d, "D125PG"), ss(d, "UCA24"),
                   ss(d, "TURNOV"), ss(d, "BMDZS"), ss(d, "EGFR"),
                   ss(d, "SYMSCOR")), 2))
  }, striped = TRUE)

  output$p_setpoint <- renderPlot({
    ica  <- seq(0.60, 1.60, by = 0.005)
    curve <- function(sp, mass) {
      sec <- 0.10 + 0.90/(1 + (ica/sp)^6)
      mass*934*sec*0.823*0.835/10.4
    }
    df <- rbind(
      data.frame(ica = ica, pth = curve(1.20, 1.00), grp = "Normal (set-point 1.20)"),
      data.frame(ica = ica, pth = curve(input$setpt, input$ptmass),
                 grp = sprintf("Patient (set-point %.2f, chief-cell %.2f)",
                               input$setpt, input$ptmass)))
    ggplot(df, aes(ica, pth, colour = grp, linetype = grp)) +
      geom_line(linewidth = 0.9) +
      annotate("rect", xmin = 1.15, xmax = 1.30, ymin = -Inf, ymax = Inf,
               fill = "#1E8449", alpha = 0.08) +
      scale_colour_manual(values = c("#5D6D7E", "#C0392B")) +
      scale_linetype_manual(values = c("dashed", "solid")) +
      labs(x = "Ionised calcium (mmol/L)", y = "PTH (pg/mL)",
           title = "Ca–PTH Set-Point Curve", colour = NULL, linetype = NULL) +
      theme_hp()
  })

  ## ---- Tab 2 : PK --------------------------------------------------------
  output$p_pth <- renderPlot(
    lineplot(sim(), "PTHPG", "PTH signal (pg/mL eq)", "Total biological PTH exposure"))
  output$p_pthzoom <- renderPlot({
    d <- sim(); d <- d[d$time >= max(0, max(d$time) - 7*24), ]
    lineplot(d, "PTHPG", "PTH signal (pg/mL eq)", "Last-7-day zoom")
  })
  output$p_d125pk <- renderPlot(
    lineplot(sim(), "D125PG", "1,25(OH)₂D (pg/mL)", "Calcitriol exposure"))
  output$p_drugpk <- renderPlot({
    d <- sim() %>% select(time, PTHEXOG, D125PG, D25NG) %>%
      pivot_longer(-time)
    ggplot(d, aes(time/24, value, colour = name)) +
      geom_line(linewidth = 0.7) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      scale_colour_manual(values = PAL) +
      labs(x = "Time (days)", y = NULL, title = "Drug · Metabolite Concentrations", colour = NULL) +
      theme_hp() + theme(legend.position = "none")
  })

  ## ---- Tab 3 : homeostasis ----------------------------------------------
  output$p_ca <- renderPlot(
    lineplot(sim(), "CACORRD", "Corrected calcium (mg/dL)",
             "Albumin-corrected serum calcium", bands = c(8.0, 9.0)))
  output$p_ica <- renderPlot(
    lineplot(sim(), "CAIONMM", "Ionised calcium (mmol/L)",
             "Ionised calcium", bands = c(1.15, 1.30)))
  output$p_pi <- renderPlot(
    lineplot(sim(), "PISER", "Phosphate (mg/dL)", "Serum phosphate", bands = c(2.5, 4.5)))
  output$p_mg <- renderPlot(
    lineplot(sim(), "MGSER", "Magnesium (mg/dL)", "Serum magnesium",
             bands = c(1.7, 2.4)))
  output$p_d125 <- renderPlot(
    lineplot(sim(), "D125PG", "1,25(OH)₂D (pg/mL)", "Active vitamin D",
             bands = c(20, 60)))
  output$p_fgf <- renderPlot(
    lineplot(sim(), "FGF23RU", "FGF23 (relative units)", "FGF23"))

  ## ---- Tab 4 : renal -----------------------------------------------------
  output$p_uca <- renderPlot(
    lineplot(sim(), "UCA24", "24h urinary calcium (mg/day)",
             "Urinary calcium — the real cost of this disease", bands = c(0, 300)))
  output$p_feca <- renderPlot(
    lineplot(sim(), "FEXCA", "FE_Ca (%)", "Fractional excretion of calcium"))
  output$p_tmp <- renderPlot(
    lineplot(sim(), "TMPPERG", "TmP/GFR (mg/dL)", "Renal phosphate threshold",
             bands = c(2.5, 4.2)))
  output$p_caxp <- renderPlot(
    lineplot(sim(), "CAPROD", "Ca × P (mg²/dL²)", "Calcium-phosphate product",
             bands = c(0, 45)))
  output$p_nc <- renderPlot(
    lineplot(sim(), "NCBURD", "Nephrocalcinosis burden (units)", "Nephrocalcinosis accumulation"))
  output$p_egfr <- renderPlot(
    lineplot(sim(), "EGFR", "eGFR (mL/min/1.73m²)", "Renal function trajectory"))

  ## ---- Tab 5 : bone ------------------------------------------------------
  output$p_turn <- renderPlot(
    lineplot(sim(), "TURNOV", "Bone turnover index (1 = normal)", "Bone remodelling activity"))
  output$p_cells <- renderPlot({
    d <- sim() %>% select(time, OB, OC) %>% pivot_longer(-time)
    ggplot(d, aes(time/24, value, colour = name)) +
      geom_line(linewidth = 0.8) +
      scale_colour_manual(values = c(OB = "#7D3C98", OC = "#C0392B"),
                          labels = c(OB = "Osteoblast", OC = "Osteoclast")) +
      labs(x = "Time (days)", y = "Activity (normal = 1)",
           title = "Osteoblast · Osteoclast", colour = NULL) + theme_hp()
  })
  output$p_markers <- renderPlot({
    d <- sim() %>% select(time, P1NP, CTX, BSAP) %>% pivot_longer(-time)
    ggplot(d, aes(time/24, value, colour = name)) +
      geom_line(linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      scale_colour_manual(values = PAL) +
      labs(x = "Time (days)", y = NULL, title = "Bone Markers", colour = NULL) +
      theme_hp() + theme(legend.position = "none")
  })
  output$p_bmd <- renderPlot(
    lineplot(sim(), "BMDZS", "Lumbar BMD Z-score", "Bone density trajectory"))

  ## ---- Tab 6 : endpoints -------------------------------------------------
  output$p_sym <- renderPlot(
    lineplot(sim(), "SYMSCOR", "Symptom score (0-100)", "Neuromuscular symptoms"))
  output$p_qtc <- renderPlot(
    lineplot(sim(), "QTCMS", "QTc (ms)", "QTc interval", bands = c(350, 450)))
  output$p_qol <- renderPlot(
    lineplot(sim(), "QOLSC", "QoL score (0-100)", "Quality of life"))
  output$p_bg <- renderPlot(
    lineplot(sim(), "BGCALCI", "Basal ganglia calcification (units)", "Basal ganglia calcification accumulation"))

  output$txt_resp3 <- renderPrint({
    d  <- sim()
    ca <- ss(d, "CACORRD")
    tp <- treat_pars()
    cat(sprintf("Albumin-corrected calcium   : %.2f mg/dL  (%s)\n", ca,
                if (ca >= 8.0 && ca <= 10.6) "Met" else "Not met"))
    cat(sprintf("Oral calcium          : %.0f mg/day   (%s)\n", tp$DCASUP,
                if (tp$DCASUP <= 600) "Met" else "Not met"))
    cat(sprintf("Active vitamin D      : %.3f µg/day  (%s)\n", tp$DCTR,
                if (tp$DCTR <= 0) "Met" else "Not met"))
    cat(sprintf("\nTriple composite response     : %s\n",
                if (mean(d$RESP3) > 0.5) "Achieved" else "Not achieved"))
    cat(sprintf("Proportion of observation period achieved: %.0f%%\n", 100*mean(d$RESP3)))
  })

  ## ---- Tab 7 : scenarios -------------------------------------------------
  scn <- eventReactive(input$goscn, {
    withProgress(message = "Running scenarios (takes tens of seconds)...", value = 0.5, {
      HYPOPT_simulate_scenarios(MOD, which = as.integer(input$scn))
    })
  })

  output$s_ca   <- renderPlot(lineplot(scn(), "CACORRD", "Corrected calcium (mg/dL)",
                                       "Serum calcium", bands = c(8.0, 9.0)))
  output$s_uca  <- renderPlot(lineplot(scn(), "UCA24", "24h urinary calcium (mg/day)",
                                       "Urinary calcium", bands = c(0, 300)))
  output$s_pi   <- renderPlot(lineplot(scn(), "PISER", "Phosphate (mg/dL)", "Serum phosphate",
                                       bands = c(2.5, 4.5)))
  output$s_egfr <- renderPlot(lineplot(scn(), "EGFR", "eGFR", "Renal function"))

  output$s_tbl   <- DT::renderDataTable(
    DT::datatable(HYPOPT_summary(scn()), options = list(pageLength = 15,
                                                        scrollX = TRUE)))
  output$s_swing <- DT::renderDataTable({
    d <- scn()
    from <- max(0, floor(max(d$time)/24) - 30)
    DT::datatable(HYPOPT_swing(d, from_day = from),
                  options = list(pageLength = 15, scrollX = TRUE))
  })

  ## ---- Tab 8 : biomarkers ------------------------------------------------
  output$b_tbl <- DT::renderDataTable({
    d <- sim()
    cols <- c("CACORRD","CAIONMM","PISER","MGSER","PTHPG","PTHENDO","PTHEXOG",
              "D125PG","D25NG","FGF23RU","TMPPERG","FRCAPCT","FEXCA","UCA24",
              "UPI24","UCAMGKG","CAABS","FABSPCT","CAPROD","EGFR","NCBURD",
              "BMDZS","P1NP","CTX","BSAP","TURNOV","SYMSCOR","QTCMS","QOLSC")
    lab <- c("Corrected calcium (mg/dL)","Ionised calcium (mmol/L)","Phosphate (mg/dL)",
             "Magnesium (mg/dL)","PTH signal (pg/mL eq)","Endogenous PTH (pg/mL)",
             "Exogenous PTH (pg/mL)","1,25D (pg/mL)","25D (ng/mL)","FGF23 (RU)",
             "TmP/GFR (mg/dL)","Tubular Ca reabsorption (%)","FE_Ca (%)",
             "24h urinary Ca (mg)","24h urinary P (mg)","Urinary Ca (mg/kg/day)",
             "Absorbed calcium (mg/day)","Fractional calcium absorption (%)","Ca×P (mg²/dL²)",
             "eGFR","Nephrocalcinosis (units)","BMD Z-score","P1NP (µg/L)",
             "CTX (ng/mL)","BSAP (µg/L)","Bone turnover index","Symptom score","QTc (ms)",
             "QoL score")
    DT::datatable(
      data.frame(Metric = lab,
                 Steady_state = round(vapply(cols, function(c) ss(d, c), 0), 3)),
      options = list(pageLength = 30), rownames = FALSE)
  })

  grid <- eventReactive(input$gogrid, {
    withProgress(message = "Calculating dose grid...", value = 0.5, {
      HYPOPT_titrate_conventional(MOD)
    })
  })
  output$b_grid <- DT::renderDataTable(
    DT::datatable(grid(), options = list(pageLength = 16, scrollX = TRUE)))

  ## ---- Tab 9 : references ------------------------------------------------
  output$refs <- renderUI({
    f <- "hypopt_references_en.md"
    if (file.exists(f)) includeMarkdown(f) else p("Could not find the references file.")
  })
}

shinyApp(ui, server)
