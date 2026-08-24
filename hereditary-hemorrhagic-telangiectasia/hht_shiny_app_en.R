## ============================================================================
##  HHT (Hereditary Hemorrhagic Telangiectasia) QSP — Shiny Dashboard
##  Hereditary Hemorrhagic Telangiectasia · interactive QSP explorer
##
##  10 tabs:
##   1. Patient Profile     patient profile & natural history
##   2. Pharmacokinetics (PK)       five drugs, including the nasal route
##   3. Lesion Biology       pSMAD · VEGF · mural coverage · lesion number and calibre
##   4. ★ Shear Stress Feedback  the sign flip, the feeder ceiling, escape thresholds
##   5. Bleeding · Epistaxis      episodes, duration, rate, monthly duration
##   6. ★ Iron · Anaemia        the iron ceiling — where clinical benefit is threshold-like
##   7. Cardiac · Shunt        cardiac index, shunt flow, the anaemia feedback loop
##   8. Clinical Endpoints    ESS, HHT-QOL, transfusion burden, iron independence
##   9. Scenario Comparison      18 treatment scenarios side by side
##  10. Safety            drug-specific adverse effects
##
##  Run:
##    setwd("hereditary-hemorrhagic-telangiectasia"); shiny::runApp("hht_shiny_app.R")
##  Required packages: shiny, mrgsolve, ggplot2, dplyr, tidyr, DT
## ============================================================================

library(shiny)
library(mrgsolve)
library(ggplot2)
library(dplyr)
library(tidyr)

## the model, phenotypes, scenarios and helpers all live in the model file
source("hht_mrgsolve_model_en.R")

THEME <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey93", colour = NA),
        legend.position = "bottom")

PAL <- c("#2471a3", "#c0392b", "#1a7f37", "#8e44ad", "#e67e22",
         "#16a085", "#7f8c8d", "#d35400")

## ---------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("HHT (Hereditary Hemorrhagic Telangiectasia) QSP Dashboard"),
  tags$p(style = "color:#555;margin-top:-8px;",
         "ALK1/ENG/SMAD4 is the regulator of the endothelial shear-stress set point — lose it, and the sign of the feedback flips. ",
         "Clinical benefit is not continuous — it is decided by whether the iron-absorption ceiling is crossed."),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Patient"),
      selectInput("pheno", "Phenotype",
                  choices = c("Mild" = "mild",
                              "Moderate (PATH-HHT)" = "moderate",
                              "Severe · transfusion-dependent (Parambil)" = "severe",
                              "Hepatic AVM, high-output (Dupuis-Girod)" = "hepatic",
                              "HHT1 pulmonary AVM" = "hht1_pav",
                              "JP-HHT (SMAD4)" = "jphht"),
                  selected = "moderate"),
      sliderInput("age", "Age at treatment start (years)", 20, 75, 52, step = 1),
      sliderInput("wt", "Body weight (kg)", 45, 120, 70, step = 1),
      hr(),
      h4("Phenotype knobs"),
      sliderInput("sev",  "Somatic hit rate SEV (lesion count)", 0.1, 2.5, 1.05, step = 0.05),
      sliderInput("angf", "Angiogenic drive ANGF (lesion size)", 0.4, 1.6, 1.02, step = 0.02),
      sliderInput("gif",  "GI lesion burden GIF", 0, 50, 0.55, step = 0.25),
      sliderInput("hepf", "Hepatic involvement HEPF", 0, 1, 0.20, step = 0.05),
      sliderInput("pulf", "Pulmonary involvement PULF", 0, 1, 0.00, step = 0.05),
      hr(),
      h4("Treatment"),
      checkboxGroupInput(
        "drugs", NULL,
        choices = c("Bevacizumab IV (maintenance)"    = "bev",
                    "Bevacizumab nasal spray"      = "bevnas",
                    "Pazopanib"                = "paz",
                    "Pomalidomide"           = "pom",
                    "Tranexamic acid (oral)"        = "txa",
                    "ALK1 pathway restoration (hypothetical)"    = "alk1"),
        selected = character(0)),
      conditionalPanel("input.drugs.indexOf('paz') > -1",
        sliderInput("pazdose", "Pazopanib (mg/day)", 25, 800, 150, step = 25)),
      conditionalPanel("input.drugs.indexOf('pom') > -1",
        sliderInput("pomdose", "Pomalidomide (mg/day)", 1, 4, 4, step = 0.5)),
      conditionalPanel("input.drugs.indexOf('bevnas') > -1",
        sliderInput("bevnasdose", "Nasal spray dose per administration (mg)", 25, 1000, 75, step = 25)),
      selectInput("iron", "Iron support",
                  choices = c("None" = "none",
                              "Oral 65 mg/day" = "oral65",
                              "Oral 200 mg/day" = "oral200",
                              "Oral + IV 1 g/8 weeks" = "iv8",
                              "Oral + IV 1 g/4 weeks" = "iv4",
                              "Oral + IV 1 g/2 weeks" = "iv2"),
                  selected = "oral200"),
      sliderInput("dur", "Simulation duration (days)", 90, 1095, 365, step = 15),
      actionButton("go", "Run", class = "btn-primary", width = "100%")
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("1. Patient Profile", br(),
                 plotOutput("p_nat", height = "420px"),
                 h5("Natural history — lesions accumulate over decades and plateau in the 50s-60s."),
                 tableOutput("t_base")),
        tabPanel("2. Pharmacokinetics (PK)", br(),
                 plotOutput("p_pk", height = "440px"),
                 helpText("Bevacizumab is IV mg/kg; pazopanib, pomalidomide, and tranexamic acid are oral. ",
                          "Because of mucociliary clearance (t½ ~15 min), the nasal route has a high peak concentration, ",
                          "but the time above threshold is negligible.")),
        tabPanel("3. Lesion Biology", br(),
                 plotOutput("p_bio", height = "480px"),
                 helpText("When pSMAD disappears, VEGF and AKT rise and mural cell coverage collapses. ",
                          "Fragility FRAG = (1 − coverage)^1.55 determines bleeding frequency.")),
        tabPanel("4. ★ Shear Stress Feedback", br(),
                 fluidRow(column(6, plotOutput("p_wss", height = "360px")),
                          column(6, plotOutput("p_remod", height = "360px"))),
                 h5("The sign flip"),
                 helpText("REMOD = KSH · g(WSS) · (1 − 2·pSMAD).  If pSMAD ≈ 1, it is negative (inward remodelling, ",
                          "self-limiting). If pSMAD ≈ 0, it is positive (outward remodelling, self-amplifying). ",
                          "Shear stress does not increase monotonically with lumen size — once the feeder vessel dominates resistance, ",
                          "it falls again. That is why nasal lesions have a stable lumen and hepatic/pulmonary AVMs do not."),
                 tableOutput("t_escape")),
        tabPanel("5. Bleeding · Epistaxis", br(),
                 plotOutput("p_bleed", height = "460px"),
                 helpText("Tranexamic acid acts only on the clot, so it reduces duration but cannot change frequency ",
                          "(the dissociation observed in the ATERO trial — in this model this is a structural prediction, not a fit).")),
        tabPanel("6. ★ Iron · Anaemia", br(),
                 plotOutput("p_iron", height = "380px"),
                 h5("The iron-absorption ceiling — a clinical threshold"),
                 tableOutput("t_ceiling"),
                 helpText("Steady state requires (absorption − obligate loss) = 0.0347 × Hb × bleeding rate. ",
                          "Because hepcidin caps daily absorption at ~10 mg, oral iron 200 mg/day is no better ",
                          "than 65 mg/day. Whether this line is crossed decides transfusion dependence.")),
        tabPanel("7. Cardiac · Shunt", br(),
                 plotOutput("p_card", height = "440px"),
                 helpText("Anaemia is both an output and an input: Hb↓ → cardiac output↑ to defend oxygen delivery → ",
                          "perfusion pressure/shear stress↑ → lesion enlargement → bleeding↑. Because of this loop, ",
                          "correcting anaemia is not simple supportive care but an intervention that removes the mechanical drive for lesion growth.")),
        tabPanel("8. Clinical Endpoints", br(),
                 plotOutput("p_end", height = "440px"),
                 tableOutput("t_end"),
                 helpText("The MCID of the ESS is 0.71 points (PMID 26393959). ",
                          "But the continuous score cannot show the threshold event of iron independence.")),
        tabPanel("9. Scenario Comparison", br(),
                 actionButton("go_scen", "Run 18 scenarios", class = "btn-success"),
                 br(), br(),
                 plotOutput("p_scen", height = "520px"),
                 tableOutput("t_scen")),
        tabPanel("10. Safety", br(),
                 plotOutput("p_ae", height = "440px"),
                 helpText("Anti-VEGF: hypertension · proteinuria · delayed wound healing. Pazopanib: elevated liver enzymes. ",
                          "Pomalidomide: neutropenia · constipation · rash · venous thrombosis. ",
                          "The model predicts that at high pazopanib doses, mural cell loss from PDGFRβ blockade ",
                          "reverses the net effect."))
      )
    )
  )
)

## ---------------------------------------------------------------------------
server <- function(input, output, session) {

  observeEvent(input$pheno, {
    p <- PHENO[[input$pheno]]
    updateSliderInput(session, "sev",  value = p$SEV)
    updateSliderInput(session, "angf", value = p$ANGF)
    updateSliderInput(session, "gif",  value = p$GIF)
    updateSliderInput(session, "hepf", value = p$HEPF)
    updateSliderInput(session, "pulf", value = if (is.null(p$PULF)) 0 else p$PULF)
  })

  pars <- reactive({
    list(SEV = input$sev, ANGF = input$angf, GIF = input$gif,
         HEPF = input$hepf, PULF = input$pulf, WT = input$wt)
  })

  iron_par <- reactive({
    switch(input$iron,
           none    = list(FE_ORAL = 0,   FE_IV_RATE = 0),
           oral65  = list(FE_ORAL = 65,  FE_IV_RATE = 0),
           oral200 = list(FE_ORAL = 200, FE_IV_RATE = 0),
           iv8     = list(FE_ORAL = 200, FE_IV_RATE = 1000 / 56),
           iv4     = list(FE_ORAL = 200, FE_IV_RATE = 1000 / 28),
           iv2     = list(FE_ORAL = 200, FE_IV_RATE = 1000 / 14))
  })

  ## ---- natural history to the chosen age ---------------------------------
  natural <- eventReactive(input$go, {
    mod %>% param(pars()) %>%
      mrgsim(end = 365 * input$age, delta = 365, atol = 1e-8, rtol = 1e-6) %>%
      as.data.frame()
  }, ignoreNULL = FALSE)

  base_state <- reactive({
    df <- natural()
    as.list(df[nrow(df), names(init(mod))])
  })

  ## ---- treatment events ---------------------------------------------------
  events <- reactive({
    d <- input$drugs
    e <- NULL
    add <- function(x) if (is.null(e)) x else c(e, x)
    n <- input$dur
    if ("bev"    %in% d) e <- add(ev(amt = 5 * input$wt, cmt = 1, ii = 28, addl = floor(n / 28)))
    if ("bevnas" %in% d) e <- add(ev(amt = input$bevnasdose, cmt = 3, ii = 14, addl = floor(n / 14)))
    if ("paz"    %in% d) e <- add(ev(amt = input$pazdose, cmt = 5, ii = 1, addl = n - 1))
    if ("pom"    %in% d) e <- add(ev(amt = input$pomdose, cmt = 7, ii = 1, addl = n - 1))
    if ("txa"    %in% d) e <- add(ev(amt = 1000, cmt = 9, ii = 1/3, addl = 3 * n - 1))
    if (is.null(e)) ev(amt = 0, cmt = 1, time = 0) else e
  })

  sim <- eventReactive(input$go, {
    pp <- c(pars(), iron_par())
    if ("alk1" %in% input$drugs) pp$ALK1AG <- 0.60
    mod %>% param(pp) %>% init(base_state()) %>%
      mrgsim(events = events(), end = input$dur, delta = 1,
             atol = 1e-10, rtol = 1e-8) %>% as.data.frame()
  }, ignoreNULL = FALSE)

  ## ---- 1. natural history --------------------------------------------------
  output$p_nat <- renderPlot({
    df <- natural() %>% mutate(age = time / 365) %>%
      select(age, `Lesion count N_n` = N_n, `Mean lumen S_n` = S_n,
             `Haemoglobin Hb` = HB, `Bleeding rate mL/day` = BLRo,
             `ESS` = ESSo, `Cardiac index CI` = CIo) %>%
      pivot_longer(-age)
    ggplot(df, aes(age, value)) +
      geom_line(colour = PAL[1], linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "Age (years)", y = NULL) + THEME
  })

  output$t_base <- renderTable({
    df <- natural(); r <- df[nrow(df), ]
    data.frame(
      `Metric` = c("ESS", "Monthly epistaxis (min)", "Monthly epistaxis episodes", "Duration per event (min)",
                 "Total bleeding rate (mL/day)", "Haemoglobin (g/dL)", "Cardiac index", "SpO2"),
      `Value` = c(round(r$ESSo, 2), round(r$MEDo, 1), round(r$EPMo, 1),
               round(r$dur_evo, 1), round(r$BLRo, 1), round(r$HB, 2),
               round(r$CIo, 2), round(r$SPO2o, 3)),
      check.names = FALSE)
  })

  ## ---- 2. PK ---------------------------------------------------------------
  output$p_pk <- renderPlot({
    df <- sim() %>%
      select(time, `Bevacizumab (µg/mL)` = C_BEVo,
             `Pazopanib free (ng/mL)` = CU_PAZo,
             `Pomalidomide (ng/mL)` = C_POMo,
             `Tranexamic acid (µg/mL)` = C_TXAo,
             `Nasal submucosal (mg)` = BEVN_T) %>%
      pivot_longer(-time) %>% filter(value > 0)
    validate(need(nrow(df) > 0, "No drug is selected — choose a treatment on the left and run."))
    ggplot(df, aes(time, value, colour = name)) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = PAL, guide = "none") +
      labs(x = "Time (days)", y = "Concentration") + THEME
  })

  ## ---- 3. lesion biology ---------------------------------------------------
  output$p_bio <- renderPlot({
    df <- sim() %>%
      select(time, `pSMAD1/5 (nasal)` = PS_n, `ID1 transcript` = ID1_n,
             `Tissue VEGF-A (pg/mL)` = VEGF_n, `PI3K/AKT` = AKT_n,
             `PDGF-B` = PDGFB_n, `Mural coverage MU_n` = MU_n,
             `Lesion count N_n` = N_n, `Mean lumen S_n` = S_n) %>%
      pivot_longer(-time)
    ggplot(df, aes(time, value)) +
      geom_line(colour = PAL[3], linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 4) +
      labs(x = "Time (days)", y = NULL) + THEME
  })

  ## ---- 4. shear feedback ---------------------------------------------------
  output$p_wss <- renderPlot({
    p <- as.list(param(mod))
    S <- seq(0.5, 50, length.out = 400)
    beds <- data.frame(
      bed = rep(c("Nasal", "GI", "Hepatic", "Pulmonary"), each = length(S)),
      S = rep(S, 4),
      w0 = rep(c(p$WSS0_n, p$WSS0_g, p$WSS0_h, p$WSS0_p), each = length(S)),
      ss = rep(c(p$SSAT_n, p$SSAT_g, p$SSAT_h, p$SSAT_p), each = length(S)))
    beds$WSS <- beds$w0 * beds$S / (1 + (beds$S / beds$ss)^4)
    ggplot(beds, aes(S, WSS, colour = bed)) +
      geom_line(linewidth = 1) +
      scale_x_log10() +
      scale_colour_manual(values = PAL) +
      labs(title = "Shear stress does not increase monotonically with lumen size",
           subtitle = "WSS ∝ S/(1+(S/S_feeder)⁴) — it falls again once the feeder vessel dominates resistance",
           x = "Lumen index S (log)", y = "Wall shear stress") + THEME
  })

  output$p_remod <- renderPlot({
    p <- as.list(param(mod))
    ps <- seq(0, 1, length.out = 200)
    wss <- c(0.2, 0.5, 1.0, 3.0)
    df <- expand.grid(pSMAD = ps, WSS = wss)
    df$REMOD <- p$KSH * (df$WSS / (p$KWSS + df$WSS)) * (1 - 2 * df$pSMAD)
    df$WSSlab <- factor(paste0("WSS = ", df$WSS))
    ggplot(df, aes(pSMAD, REMOD, colour = WSSlab)) +
      geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
      geom_vline(xintercept = 0.5, linetype = 3, colour = "grey40") +
      geom_line(linewidth = 1) +
      annotate("text", x = 0.15, y = max(df$REMOD) * 0.85,
               label = "Outward\nself-amplifying", size = 3.4, colour = "#c0392b") +
      annotate("text", x = 0.85, y = min(df$REMOD) * 0.85,
               label = "Inward\nself-limiting", size = 3.4, colour = "#1a7f37") +
      scale_colour_manual(values = PAL, name = NULL) +
      labs(title = "The sign flip of the feedback",
           subtitle = "REMOD = KSH · g(WSS) · (1 − 2·pSMAD)",
           x = "Relative pSMAD1/5 activity", y = "Remodelling term") + THEME
  })

  output$t_escape <- renderTable({
    p <- as.list(param(mod))
    ps <- p$GD_LES * p$BMP9 / (p$KD9 + p$BMP9)
    beds <- list(n = c(p$WSS0_n, p$SSAT_n, 0.52), g = c(p$WSS0_g, p$SSAT_g, 0.52),
                 h = c(p$WSS0_h, p$SSAT_h, 0), p = c(p$WSS0_p, p$SSAT_p, 0),
                 c = c(p$WSS0_c, p$SSAT_c, 0))
    nm <- c(n = "Nasal", g = "GI", h = "Hepatic", p = "Pulmonary", c = "Cerebral")
    do.call(rbind, lapply(c(1.0, 1.45, 2.0, 2.6), function(co) {
      row <- sapply(names(beds), function(b) {
        v <- beds[[b]]
        Sp <- v[2] / 3^0.25                       # calibre where shear peaks
        w <- v[1] * co * Sp / (1 + (Sp / v[2])^4)
        remod <- p$KSH * (w / (p$KWSS + w)) * (1 - 2 * ps)
        sink <- p$KMAT * v[3] + p$KDEC
        sprintf("%s (%.2f)", ifelse(remod > sink, "Escape", "Suppressed"), remod / sink)
      })
      setNames(data.frame(sprintf("%.2f", co), t(row), stringsAsFactors = FALSE),
               c("CO_rel", nm[names(beds)]))
    }))
  }, caption = "When a fully effective antiangiogenic drug removes all drive, is the lesion sustained by shear stress alone? (ratio > 1 = drug-refractory)",
     caption.placement = "top")

  ## ---- 5. bleeding ---------------------------------------------------------
  output$p_bleed <- renderPlot({
    df <- sim() %>%
      select(time, `Monthly epistaxis (min)` = MEDo, `Monthly epistaxis episodes` = EPMo,
             `Duration per event (min)` = dur_evo, `Rate per event (mL/min)` = rate_evo,
             `Nasal bleeding rate (mL/day)` = BLR_No, `GI bleeding rate (mL/day)` = BLR_Go) %>%
      pivot_longer(-time)
    ggplot(df, aes(time, value)) +
      geom_line(colour = PAL[2], linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "Time (days)", y = NULL) + THEME
  })

  ## ---- 6. iron -------------------------------------------------------------
  output$p_iron <- renderPlot({
    df <- sim() %>%
      select(time, `Haemoglobin (g/dL)` = HB, `Iron stores (mg)` = FES,
             `Hepcidin` = HEP, `Bleeding iron loss (mg/day)` = FE_LOSSo,
             `Absorbed iron (mg/day)` = FE_ABSo,
             `Maximum sustainable bleeding (mL/day)` = FE_CEILo,
             `Actual bleeding rate (mL/day)` = BLRo) %>%
      pivot_longer(-time)
    ggplot(df, aes(time, value)) +
      geom_line(colour = PAL[5], linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "Time (days)", y = NULL) + THEME
  })

  output$t_ceiling <- renderTable({ iron_ceiling() },
    caption = "Maximum bleeding rate each iron-support strategy can sustain (at Hb 14.6 g/dL)",
    caption.placement = "top")

  ## ---- 7. cardiac ----------------------------------------------------------
  output$p_card <- renderPlot({
    df <- sim() %>%
      select(time, `Cardiac index (L/min/m²)` = CIo, `Cardiac output (L/min)` = CO,
             `Hepatic shunt index S_h` = S_h, `Pulmonary shunt index S_p` = S_p,
             `SpO2` = SPO2o, `Left ventricular remodelling` = LVR,
             `Right atrial pressure index` = RAP, `Haemoglobin` = HB) %>%
      pivot_longer(-time)
    ggplot(df, aes(time, value)) +
      geom_line(colour = PAL[4], linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 4) +
      labs(x = "Time (days)", y = NULL) + THEME
  })

  ## ---- 8. endpoints --------------------------------------------------------
  output$p_end <- renderPlot({
    df <- sim() %>%
      select(time, `ESS (0-10)` = ESSo, `HHT quality of life (0-16)` = QOL,
             `Cumulative bleeding (mL)` = CUM_BL, `Cumulative transfusion (unit-equivalents)` = CUM_TX,
             `Iron independence (1=yes)` = IRONINDo, `Haemoglobin` = HB) %>%
      pivot_longer(-time)
    ggplot(df, aes(time, value)) +
      geom_line(colour = PAL[1], linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "Time (days)", y = NULL) + THEME
  })

  output$t_end <- renderTable({
    d <- sim(); a <- d[1, ]; b <- d[nrow(d), ]
    data.frame(
      `Metric` = c("ESS", "HHT-QOL", "Monthly epistaxis (min)", "Total bleeding rate (mL/day)",
                 "Haemoglobin (g/dL)", "Cardiac index", "Iron independence"),
      `Baseline` = c(round(a$ESSo, 2), round(a$QOL, 2), round(a$MEDo, 1),
                 round(a$BLRo, 1), round(a$HB, 2), round(a$CIo, 2),
                 ifelse(a$IRONINDo > 0, "Yes", "No")),
      `End` = c(round(b$ESSo, 2), round(b$QOL, 2), round(b$MEDo, 1),
                 round(b$BLRo, 1), round(b$HB, 2), round(b$CIo, 2),
                 ifelse(b$IRONINDo > 0, "Yes", "No")),
      `Change` = c(round(b$ESSo - a$ESSo, 2), round(b$QOL - a$QOL, 2),
                 round(b$MEDo - a$MEDo, 1), round(b$BLRo - a$BLRo, 1),
                 round(b$HB - a$HB, 2), round(b$CIo - a$CIo, 2), ""),
      check.names = FALSE)
  })

  ## ---- 9. scenarios --------------------------------------------------------
  scen <- eventReactive(input$go_scen, {
    withProgress(message = "Running 18 scenarios…", value = 0.5, {
      run_scenarios(input$pheno, input$age, input$dur)
    })
  })

  output$p_scen <- renderPlot({
    df <- scen() %>% mutate(scenario = reorder(scenario, dESS))
    ggplot(df, aes(dESS, scenario, fill = dESS)) +
      geom_col() +
      geom_vline(xintercept = -0.71, linetype = 2, colour = "#c0392b") +
      annotate("text", x = -0.71, y = 1, label = "MCID 0.71", hjust = -0.08,
               size = 3.2, colour = "#c0392b") +
      scale_fill_gradient(low = "#1a7f37", high = "#c0392b", guide = "none") +
      labs(x = "ESS change (negative = improvement)", y = NULL,
           title = "Treatment scenario comparison", subtitle = "Red dashed line = minimal clinically important difference of 0.71 points") +
      THEME
  })

  output$t_scen <- renderTable({ scen() })

  ## ---- 10. safety ----------------------------------------------------------
  output$p_ae <- renderPlot({
    df <- sim() %>%
      select(time, `Systolic BP (mmHg)` = SBP, `Proteinuria (UPCR)` = UPCR,
             `ALT (U/L)` = ALT, `Neutrophils (10⁹/L)` = ANC,
             `Cumulative VTE risk` = VTEH, `Mural coverage MU_n` = MU_n) %>%
      pivot_longer(-time)
    ggplot(df, aes(time, value)) +
      geom_line(colour = PAL[8], linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "Time (days)", y = NULL) + THEME
  })
}

shinyApp(ui, server)
