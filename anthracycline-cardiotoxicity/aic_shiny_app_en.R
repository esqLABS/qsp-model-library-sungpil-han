## =============================================================================
## Anthracycline-Induced Cardiotoxicity (AIC) — Shiny dashboard
## Anthracycline-Induced Cardiotoxicity QSP Model Interactive Dashboard
##
## Eight tabs, each built around one thing the model claims:
##   1. Patient Profile (Patient)      risk phenotype -> which parameters move
##   2. PK · Myocardial Exposure (Exposure)    the two cardiac pools and their metrics
##   3. Injury Axes (Injury arms)        Top2b/DSB/p53 vs iron/ROS/mitochondria
##   4. Two Deficits and the Mask (Deficits)     MYO vs FUNC, and LVEF with/without the mask
##   5. Clinical Endpoints   LVEF, GLS, CTRCD, NYHA
##   6. Scenario Comparison (Scenarios)     same cumulative dose, different schedules
##   7. Biomarkers       troponin/NT-proBNP lead times
##   8. Window of Reversibility (Window)            when HF therapy is started
##
## Run:  shiny::runApp("aic_shiny_app.R")
## Requires: shiny, mrgsolve, dplyr, tidyr, ggplot2, DT
## =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

MODEL_FILE <- "aic_mrgsolve_model_en.R"
mod <- mread_cache("aic", MODEL_FILE)

BSA_DEFAULT <- 1.8
WT_DEFAULT  <- 70

## ---------------------------------------------------------------------------
## dosing helpers (mirror the R driver in aic_mrgsolve_model_en.R)
## ---------------------------------------------------------------------------
dox_ev <- function(dose_mg_m2, ncyc, interval, bsa, tinf = 0, lipo = FALSE,
                   start = 0) {
  ev(amt = dose_mg_m2 * bsa, cmt = if (lipo) 4 else 1, time = start,
     ii = interval, addl = ncyc - 1,
     rate = if (tinf > 0) dose_mg_m2 * bsa / tinf else 0)
}
dex_ev <- function(dose_mg_m2, ncyc, interval, bsa, ratio = 10) {
  ev(amt = ratio * dose_mg_m2 * bsa, cmt = 10, time = 0, ii = interval,
     addl = ncyc - 1)
}
tras_ev <- function(start, wt, weeks = 52) {
  c(ev(amt = 8 * wt, cmt = 12, time = start),
    ev(amt = 6 * wt, cmt = 12, time = start + 21, ii = 21,
       addl = floor(weeks / 3) - 2))
}

SCHEDULES <- c("IV push (q3w)"          = "push",
               "72-h continuous infusion" = "inf72",
               "Weekly fractionated"    = "weekly",
               "Pegylated liposomal (PLD)" = "pld")

theme_aic <- function() {
  theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank(),
          legend.position = "bottom",
          plot.title = element_text(face = "bold"))
}

## ---------------------------------------------------------------------------
## UI
## ---------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel(paste("QSP Model —",
                   "Anthracycline-Induced Cardiotoxicity")),
  helpText(paste("Cumulative dose is not an exposure metric: even at the same mg/m², changing the peak",
                 "concentration changes the injury. LVEF is a lagging indicator masked by compensatory",
                 "hypertrophy, and the chance of recovery is determined by the fibrosis clock.")),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Patient"),
      sliderInput("age", "Age (years)", 20, 85, 55, 1),
      sliderInput("bsa", "Body surface area BSA (m²)", 1.3, 2.3, BSA_DEFAULT, 0.05),
      sliderInput("ef0", "Baseline LVEF (%)", 45, 72, 62, 1),
      checkboxInput("htn", "Hypertension", FALSE),
      checkboxInput("rt", "Prior mediastinal radiation therapy (prior mediastinal RT)", FALSE),
      selectInput("cbr1", "CBR1/AKR1C3 metabolic phenotype (doxorubicinol FM)",
                  c("Slow (0.12)" = 0.12, "Intermediate (0.18)" = 0.18,
                    "Typical (0.25)" = 0.25, "Fast (0.34)" = 0.34,
                    "Very fast (0.45)" = 0.45), selected = 0.25),

      hr(), h4("Chemotherapy"),
      selectInput("sched", "Schedule · Formulation", SCHEDULES, selected = "push"),
      sliderInput("dose", "Dose per administration (mg/m²)", 10, 90, 60, 5),
      sliderInput("ncyc", "Number of cycles", 1, 24, 8, 1),
      sliderInput("ii", "Cycle interval (days)", 7, 28, 21, 7),
      helpText(textOutput("cumtxt")),

      hr(), h4("Cardioprotection"),
      checkboxInput("dex", "Dexrazoxane 10:1", FALSE),
      checkboxInput("sta", "Atorvastatin 40 mg", FALSE),
      sliderInput("sta_day", "Statin start day (day)", 0, 400, 0, 10),
      checkboxInput("ace", "ACE inhibitor / ARB", FALSE),
      checkboxInput("bb", "Beta-blocker (carvedilol)", FALSE),
      checkboxInput("arni", "Sacubitril/valsartan (ARNI)", FALSE),
      checkboxInput("sglt", "SGLT2 inhibitor", FALSE),
      sliderInput("hf_day", "ACEi/BB/ARNI start day (day)", 0, 540, 0, 30),

      hr(), h4("HER2-targeted therapy (Trastuzumab)"),
      selectInput("tras", "Trastuzumab",
                  c("None" = "none",
                    "Sequential (day 84)" = "seq",
                    "Concurrent (day 0)" = "conc"), selected = "none"),

      hr(), h4("Simulation"),
      sliderInput("tend", "Observation period (days)", 365, 1460, 730, 30),
      sliderInput("nid", "Number of virtual patients (virtual subjects)", 1, 500, 100, 1),
      actionButton("go", "Run", class = "btn-primary")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel(
          "1. Patient Profile",
          h4("Which parameters does the risk phenotype move"),
          tableOutput("risk_tbl"),
          helpText(paste("Older age · hypertension · radiation history are not a separate 'risk score' —",
                         "each one differently moves myocyte death susceptibility (KDROS), fibrotic tendency (KFIBIN),",
                         "regenerative capacity (KREG), and baseline reserve (EF0).")),
          h4("Compensatory reserve (masking reserve)"),
          plotOutput("reserve_plot", height = "300px"),
          helpText(paste("Because HYPMAX = 0.32, once the effective contractile-unit deficit exceeds about 20%,",
                         "compensation saturates and LVEF falls sharply.",
                         "If baseline LVEF is low, this margin is small."))
        ),

        tabPanel(
          "2. PK · Myocardial Exposure",
          h4("Plasma · Metabolite Concentration"),
          plotOutput("pk_plot", height = "260px"),
          h4("Two Pools Within the Myocardium — Different Exposure Metrics"),
          plotOutput("cardiac_plot", height = "300px"),
          tableOutput("exposure_tbl"),
          helpText(paste("The fast nuclear pool (CHF) tracks the peak, while the slow retained pool (CH,",
                         "CHM) integrates the AUC. Changing the schedule alters only the",
                         "former — this is why continuous infusion and liposomal formulations work."))
        ),

        tabPanel(
          "3. Injury Axes",
          fluidRow(
            column(6, h4("Top2b / DSB / p53 (peak-driven)"),
                   plotOutput("arm_nuc", height = "300px")),
            column(6, h4("Iron · ROS · Mitochondria (AUC-driven)"),
                   plotOutput("arm_redox", height = "300px"))),
          h4("Top2b Protein Level (Dexrazoxane Effect)"),
          plotOutput("t2b_plot", height = "220px"),
          helpText(paste("The two axes combine multiplicatively in the death rate:",
                         "kdeath ∝ Hill(ROS;n=3) × (1 + 1.2·p53).",
                         "So removing even just one side gives a benefit greater than its share."))
        ),

        tabPanel(
          "4. Two Deficits and the Mask",
          h4("Reversible Deficit FUNC vs Irreversible Myocyte Loss"),
          plotOutput("deficit_plot", height = "300px"),
          h4("Size of the Mask: Where Would LVEF Be Without Compensatory Hypertrophy"),
          plotOutput("mask_plot", height = "300px"),
          tableOutput("mask_tbl"),
          helpText(paste("The difference between LVEF (actual) and LVEF_NOMASK is the number of points",
                         "hidden by compensatory hypertrophy, and the difference from LVEF_IRREV is the",
                         "reversible component that can still be recovered."))
        ),

        tabPanel(
          "5. Clinical Endpoints",
          fluidRow(
            column(6, h4("LVEF (%)"), plotOutput("ef_plot", height = "300px")),
            column(6, h4("|GLS| (%)"), plotOutput("gls_plot", height = "300px"))),
          h4("Event incidence (virtual population)"),
          DTOutput("endpoint_tbl"),
          helpText(paste("CTRCD = LVEF drop of ≥10 points & <50% (ASE/ESC).",
                         "A relative 15% decrease in GLS captures the same injury at an earlier",
                         "time point."))
        ),

        tabPanel(
          "6. Scenario Comparison",
          h4("Same Cumulative Dose, Different Schedule · Formulation · Cardioprotectant"),
          actionButton("go_cmp", "Run scenario comparison", class = "btn-success"),
          plotOutput("cmp_plot", height = "340px"),
          DTOutput("cmp_tbl"),
          helpText(paste("Every row has the same cumulative dose. The difference comes entirely",
                         "from peak concentration and the Top2b axis — the fact that anti-tumour AUC",
                         "is preserved is the basis for the improved therapeutic index."))
        ),

        tabPanel(
          "7. Biomarkers",
          h4("hs-cTnI — a Marker of Injury Rate"),
          plotOutput("tni_plot", height = "280px"),
          h4("NT-proBNP — a Marker of Wall Stress"),
          plotOutput("bnp_plot", height = "250px"),
          h4("Lead Time"),
          tableOutput("lead_tbl"),
          helpText(paste("Because troponin's half-life is about 12 hours, it reads not cumulative damage but",
                         "'the rate of dying right now'. The fact that CTRCD is rare in patients who are",
                         "negative every cycle (high negative predictive value) is the",
                         "basis for the troponin-guided strategy."))
        ),

        tabPanel(
          "8. Window of Reversibility",
          h4("When to Start Heart-Failure Therapy"),
          actionButton("go_win", "Calculate window of reversibility", class = "btn-warning"),
          plotOutput("win_plot", height = "320px"),
          DTOutput("win_tbl"),
          helpText(paste("Recovery potential tracks not time itself but the level of fibrosis at that",
                         "time point: because the FUNC recovery rate slows as 1/(1+3·FIB).",
                         "What closes the window is the fibrosis clock."))
        )
      )
    )
  )
)

## ---------------------------------------------------------------------------
## server
## ---------------------------------------------------------------------------
server <- function(input, output, session) {

  output$cumtxt <- renderText({
    sprintf("Cumulative dose = %d mg/m²  (total %.0f mg, BSA %.2f m²)",
            input$dose * input$ncyc, input$dose * input$ncyc * input$bsa,
            input$bsa)
  })

  ## risk phenotype -> parameter multipliers
  risk_par <- reactive({
    p <- list(EF0 = input$ef0, FM = as.numeric(input$cbr1))
    kd <- 1; kf <- 1; kfib <- 1; kreg <- 1
    if (input$age >= 70) { kd <- kd * 1.6; kreg <- kreg * 0.4 }
    else if (input$age >= 60) { kd <- kd * 1.2; kreg <- kreg * 0.7 }
    if (input$age < 18) { kreg <- kreg * 1.4 }
    if (input$htn) { kfib <- kfib * 1.3; kf <- kf * 1.2 }
    if (input$rt)  { kfib <- kfib * 1.6; kd <- kd * 1.2 }
    p$KDROS  <- 4.8e-4 * kd
    p$KFIN   <- 1.0e-3 * kf
    p$KFIBIN <- 8.0e-3 * kfib
    p$KREG   <- 2.0e-5 * kreg
    p
  })

  drug_par <- reactive({
    p <- list()
    p$TGT_STA  <- if (input$sta)  0.75 else 0
    p$TON_STA  <- input$sta_day
    p$TGT_ACE  <- if (input$ace)  0.80 else 0
    p$TGT_BB   <- if (input$bb)   0.80 else 0
    p$TGT_ARNI <- if (input$arni) 0.90 else 0
    p$TGT_SGLT <- if (input$sglt) 0.70 else 0
    p$TON_ACE  <- input$hf_day
    p$TON_BB   <- input$hf_day
    p$TON_ARNI <- if (input$arni) input$hf_day else 1e6
    p$TON_SGLT <- if (input$sglt) 0 else 1e6
    p
  })

  regimen <- reactive({
    s <- input$sched
    e <- switch(
      s,
      push   = dox_ev(input$dose, input$ncyc, input$ii, input$bsa),
      inf72  = dox_ev(input$dose, input$ncyc, input$ii, input$bsa, tinf = 3),
      weekly = dox_ev(input$dose, input$ncyc, 7, input$bsa),
      pld    = dox_ev(input$dose, input$ncyc, 28, input$bsa,
                      lipo = TRUE))
    if (input$dex && s != "pld")
      e <- c(e, dex_ev(input$dose, input$ncyc, input$ii, input$bsa))
    if (input$tras == "seq")  e <- c(e, tras_ev(84, WT_DEFAULT))
    if (input$tras == "conc") e <- c(e, tras_ev(0,  WT_DEFAULT))
    e
  })

  sim <- eventReactive(input$go, {
    m <- param(mod, c(risk_par(), drug_par()))
    n <- input$nid
    out <- m %>% ev(regimen()) %>%
      mrgsim(end = input$tend, delta = 1, nid = n) %>% as_tibble()
    out
  }, ignoreNULL = FALSE)

  idx <- reactive({
    m <- param(mod, c(risk_par(), drug_par()))
    m %>% ev(regimen()) %>% mrgsim(end = input$tend, delta = 1) %>% as_tibble()
  })

  ## ---- tab 1 ------------------------------------------------------------
  output$risk_tbl <- renderTable({
    p <- risk_par()
    tibble(Parameter = c("KDROS (death susceptibility)", "KFIN (functional-injury susceptibility)",
                        "KFIBIN (fibrotic tendency)", "KREG (myocyte regeneration)",
                        "EF0 (baseline LVEF)", "FM (doxorubicinol formation fraction)"),
           Reference = c("4.8e-4", "1.0e-3", "8.0e-3", "2.0e-5", "62", "0.25"),
           Current = c(sprintf("%.2e", p$KDROS), sprintf("%.2e", p$KFIN),
                      sprintf("%.2e", p$KFIBIN), sprintf("%.2e", p$KREG),
                      sprintf("%.0f", p$EF0), sprintf("%.2f", p$FM)))
  })

  output$reserve_plot <- renderPlot({
    d <- tibble(deficit = seq(0, 0.5, 0.005)) %>%
      mutate(HYP = pmin(0.32, 1.6 * deficit),
             CONT = (1 - deficit) * (1 + HYP),
             LVEF = input$ef0 * CONT^0.9,
             LVEF_nomask = input$ef0 * (1 - deficit)^0.9)
    ggplot(d, aes(100 * deficit)) +
      geom_line(aes(y = LVEF, colour = "After compensation (observed LVEF)"), linewidth = 1.2) +
      geom_line(aes(y = LVEF_nomask, colour = "Absent compensation"),
                linewidth = 1.0, linetype = 2) +
      geom_vline(xintercept = 100 * 0.32 / 1.6, linetype = 3) +
      annotate("text", x = 20.5, y = input$ef0 * 0.55,
               label = "Reserve saturation (HYPMAX)", hjust = 0, size = 4) +
      labs(x = "Effective contractile-unit deficit (%)", y = "LVEF (%)", colour = NULL,
           title = "The algebraic form of the mask") +
      theme_aic()
  })

  ## ---- tab 2 ------------------------------------------------------------
  output$pk_plot <- renderPlot({
    d <- idx() %>% select(time, CPLASMA, CDOXOL) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 0.8) +
      scale_y_log10() +
      labs(x = "day", y = "concentration (mg/L, log)", colour = NULL,
           title = "Doxorubicin · Doxorubicinol") + theme_aic()
  })

  output$cardiac_plot <- renderPlot({
    d <- idx() %>% select(time, CHF, CH, CHM) %>% pivot_longer(-time)
    lab <- c(CHF = "CHF fast nuclear pool (tracks peak)",
             CH = "CH slow retained pool (integrates AUC)",
             CHM = "CHM myocardial doxorubicinol (accumulation)")
    ggplot(d, aes(time, value, colour = lab[name])) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~lab[name], scales = "free_y", ncol = 1) +
      labs(x = "day", y = "cardiac pool", colour = NULL) +
      theme_aic() + theme(legend.position = "none")
  })

  output$exposure_tbl <- renderTable({
    d <- idx()
    tibble(Metric = c("peak CHF (nuclear peak)", "peak CH (retained)",
                    "peak CHM (metabolite)", "AUC_CH (cumulative retained exposure)",
                    "peak p53 (genotoxic memory)"),
           Value = sprintf("%.3f", c(max(d$CHF), max(d$CH), max(d$CHM),
                                  max(d$AUCH), max(d$P53))))
  })

  ## ---- tab 3 ------------------------------------------------------------
  output$arm_nuc <- renderPlot({
    d <- idx() %>% select(time, DSB, P53) %>% pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      labs(x = "day", y = "level", colour = NULL) + theme_aic()
  })
  output$arm_redox <- renderPlot({
    d <- idx() %>% select(time, ROS, MITOD, LIP) %>% pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      labs(x = "day", y = "level", colour = NULL) + theme_aic()
  })
  output$t2b_plot <- renderPlot({
    ggplot(idx(), aes(time, TOP2B)) + geom_line(linewidth = 0.9) +
      ylim(0, 1.05) +
      labs(x = "day", y = "Top2b (1 = normal)",
           title = "Dexrazoxane degrades Top2b protein (resynthesis t½ 2 d)") +
      theme_aic()
  })

  ## ---- tab 4 ------------------------------------------------------------
  output$deficit_plot <- renderPlot({
    d <- idx() %>%
      transmute(time, `Irreversible myocyte loss (1-MYO)` = 1 - MYOCYTES,
                `Reversible functional deficit FUNC` = FUNCDEF,
                `Fibrosis FIB` = FIBROSIS, `Compensatory hypertrophy HYP` = MASK) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 1) +
      labs(x = "day", y = "fraction", colour = NULL) + theme_aic()
  })

  output$mask_plot <- renderPlot({
    d <- idx() %>% select(time, LVEF, LVEF_NOMASK, LVEF_IRREV) %>%
      pivot_longer(-time)
    lab <- c(LVEF = "LVEF (observed)", LVEF_NOMASK = "With compensatory hypertrophy removed",
             LVEF_IRREV = "With the reversible component recovered (ceiling)")
    ggplot(d, aes(time, value, colour = lab[name])) +
      geom_line(linewidth = 1) + geom_hline(yintercept = 50, linetype = 3) +
      labs(x = "day", y = "LVEF (%)", colour = NULL) + theme_aic()
  })

  output$mask_tbl <- renderTable({
    idx() %>% filter(time %in% c(0, 60, 120, 170, 240, 365, 540, 730)) %>%
      transmute(day = time, LVEF = round(LVEF, 2),
                `Mask score` = round(LVEF - LVEF_NOMASK, 2),
                `Reversible component score` = round(LVEF_IRREV - LVEF, 2),
                `Myocytes (%)` = round(100 * MYOCYTES, 1),
                `FIB` = round(FIBROSIS, 3))
  })

  ## ---- tab 5 ------------------------------------------------------------
  output$ef_plot <- renderPlot({
    d <- sim()
    q <- d %>% group_by(time) %>%
      summarise(med = median(LVEF), lo = quantile(LVEF, .05),
                hi = quantile(LVEF, .95), .groups = "drop")
    ggplot(q, aes(time)) +
      geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .2) +
      geom_line(aes(y = med), linewidth = 1.1) +
      geom_hline(yintercept = 50, linetype = 3) +
      geom_hline(yintercept = 40, linetype = 2, colour = "red") +
      labs(x = "day", y = "LVEF (%)",
           title = "Median and 90% interval") + theme_aic()
  })

  output$gls_plot <- renderPlot({
    d <- sim()
    q <- d %>% group_by(time) %>%
      summarise(med = median(GLS), lo = quantile(GLS, .05),
                hi = quantile(GLS, .95), .groups = "drop")
    ggplot(q, aes(time)) +
      geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .2) +
      geom_line(aes(y = med), linewidth = 1.1) +
      labs(x = "day", y = "|GLS| (%)",
           title = "GLS does not receive hypertrophic compensation") + theme_aic()
  })

  output$endpoint_tbl <- renderDT({
    d <- sim()
    s <- d %>% group_by(ID) %>%
      summarise(EF0 = first(LVEF), nadir = min(LVEF),
                CTRCD = any(LVEF < 50 & (first(LVEF) - LVEF) >= 10),
                HF = any(LVEF < 40),
                GLS15 = any((first(GLS) - GLS) / first(GLS) >= 0.15),
                .groups = "drop")
    datatable(tibble(
      Metric = c("CTRCD (LVEF ≥10-point↓ & <50%)", "Symptomatic heart-failure surrogate (LVEF <40%)",
               "GLS relative decrease ≥15%", "Mean LVEF nadir",
               "Mean LVEF decline (nadir)"),
      Value = c(sprintf("%.1f %%", 100 * mean(s$CTRCD)),
             sprintf("%.1f %%", 100 * mean(s$HF)),
             sprintf("%.1f %%", 100 * mean(s$GLS15)),
             sprintf("%.2f %%", mean(s$nadir)),
             sprintf("%.2f points", mean(s$EF0 - s$nadir)))),
      options = list(dom = "t"), rownames = FALSE)
  })

  ## ---- tab 6 ------------------------------------------------------------
  cmp <- eventReactive(input$go_cmp, {
    cum <- input$dose * input$ncyc
    base <- param(mod, c(risk_par(), list(TGT_STA = 0, TGT_ACE = 0,
                                          TGT_BB = 0, TGT_ARNI = 0,
                                          TGT_SGLT = 0)))
    defs <- list(
      list(k = "IV push", e = dox_ev(60, round(cum / 60), 21, input$bsa), p = list()),
      list(k = "72-h infusion",
           e = dox_ev(60, round(cum / 60), 21, input$bsa, tinf = 3), p = list()),
      list(k = "weekly 20 mg/m²",
           e = dox_ev(20, round(cum / 20), 7, input$bsa), p = list()),
      list(k = "liposomal (PLD)",
           e = dox_ev(40, round(cum / 40), 28, input$bsa, lipo = TRUE), p = list()),
      list(k = "push + dexrazoxane",
           e = c(dox_ev(60, round(cum / 60), 21, input$bsa),
                 dex_ev(60, round(cum / 60), 21, input$bsa)), p = list()),
      list(k = "push + statin",
           e = dox_ev(60, round(cum / 60), 21, input$bsa),
           p = list(TGT_STA = 0.75, TON_STA = 0)),
      list(k = "push + ACEi/BB",
           e = dox_ev(60, round(cum / 60), 21, input$bsa),
           p = list(TGT_ACE = 0.8, TON_ACE = 0, TGT_BB = 0.8, TON_BB = 0)))
    bind_rows(lapply(defs, function(d) {
      m <- if (length(d$p)) param(base, d$p) else base
      o <- m %>% ev(d$e) %>%
        mrgsim(end = input$tend, delta = 1, nid = input$nid) %>% as_tibble()
      i <- m %>% ev(d$e) %>% mrgsim(end = input$tend, delta = 1) %>%
        as_tibble()
      o %>% group_by(ID) %>%
        summarise(CTRCD = any(LVEF < 50 & (first(LVEF) - LVEF) >= 10),
                  HF = any(LVEF < 40), nadir = min(LVEF), .groups = "drop") %>%
        summarise(scenario = d$k, `CTRCD %` = 100 * mean(CTRCD),
                  `LVEF<40 %` = 100 * mean(HF),
                  `mean nadir` = mean(nadir)) %>%
        mutate(`peak nuclear CHF` = max(i$CHF), `peak p53` = max(i$P53),
               `AUC retained` = max(i$AUCH))
    }))
  })

  output$cmp_tbl <- renderDT({
    datatable(cmp() %>% mutate(across(where(is.numeric), ~round(.x, 2))),
              options = list(dom = "t", pageLength = 10), rownames = FALSE)
  })

  output$cmp_plot <- renderPlot({
    d <- cmp()
    ggplot(d, aes(reorder(scenario, `CTRCD %`), `CTRCD %`)) +
      geom_col(fill = "#c0392b", alpha = .85) +
      geom_text(aes(label = sprintf("%.1f%%", `CTRCD %`)), hjust = -0.15) +
      coord_flip() + ylim(0, max(d$`CTRCD %`) * 1.2) +
      labs(x = NULL, y = "CTRCD (%)",
           title = sprintf("Cumulative dose fixed at %d mg/m² — only schedule · formulation · cardioprotectant differ",
                           input$dose * input$ncyc)) +
      theme_aic()
  })

  ## ---- tab 7 ------------------------------------------------------------
  output$tni_plot <- renderPlot({
    d <- sim()
    q <- d %>% group_by(time) %>%
      summarise(med = median(CTNI), hi = quantile(CTNI, .95), .groups = "drop")
    ggplot(q, aes(time)) +
      geom_ribbon(aes(ymin = 0, ymax = hi), alpha = .15) +
      geom_line(aes(y = med), linewidth = 1) +
      geom_hline(yintercept = 14, linetype = 3) +
      labs(x = "day", y = "hs-cTnI (ng/L)",
           title = "Dotted line = 14 ng/L (99th percentile)") + theme_aic()
  })

  output$bnp_plot <- renderPlot({
    d <- sim()
    q <- d %>% group_by(time) %>%
      summarise(med = median(NTBNP), hi = quantile(NTBNP, .95),
                .groups = "drop")
    ggplot(q, aes(time)) +
      geom_ribbon(aes(ymin = 0, ymax = hi), alpha = .15) +
      geom_line(aes(y = med), linewidth = 1) +
      labs(x = "day", y = "NT-proBNP (pg/mL)") + theme_aic()
  })

  output$lead_tbl <- renderTable({
    d <- sim()
    f <- function(cond_expr) {
      d %>% group_by(ID) %>% filter({{ cond_expr }}) %>%
        summarise(day = min(time), .groups = "drop")
    }
    tni <- f(CTNI > 14)
    gls <- d %>% group_by(ID) %>%
      filter((first(GLS) - GLS) / first(GLS) >= 0.15) %>%
      summarise(day = min(time), .groups = "drop")
    ef <- d %>% group_by(ID) %>%
      filter(LVEF < 50, (first(LVEF) - LVEF) >= 10) %>%
      summarise(day = min(time), .groups = "drop")
    nid <- n_distinct(d$ID)
    tibble(Metric = c("hs-cTnI > 14 ng/L", "GLS relative 15% decrease", "CTRCD (LVEF)"),
           `Fraction of patients reaching` = sprintf("%.1f %%",
                                      100 * c(nrow(tni), nrow(gls),
                                              nrow(ef)) / nid),
           `Median day reached` = c(if (nrow(tni)) median(tni$day) else NA,
                             if (nrow(gls)) median(gls$day) else NA,
                             if (nrow(ef)) median(ef$day) else NA))
  })

  ## ---- tab 8 ------------------------------------------------------------
  win <- eventReactive(input$go_win, {
    starts <- c(60, 90, 120, 180, 270, 365, 540)
    bind_rows(lapply(starts, function(t0) {
      m <- param(mod, c(risk_par(), list(TGT_ACE = 0.8, TON_ACE = t0,
                                         TGT_BB = 0.8, TON_BB = t0,
                                         TGT_STA = 0, TGT_ARNI = 0,
                                         TGT_SGLT = 0)))
      o <- m %>% ev(regimen()) %>%
        mrgsim(end = t0 + 365, delta = 1, nid = input$nid) %>% as_tibble()
      o %>% group_by(ID) %>%
        summarise(had = any(LVEF < 50 & (first(LVEF) - LVEF) >= 10 &
                              time <= 365),
                  EF_base = first(LVEF),
                  EF_start = LVEF[which.min(abs(time - t0))],
                  FIB_start = FIBROSIS[which.min(abs(time - t0))],
                  EF_end = last(LVEF), .groups = "drop") %>%
        filter(had) %>%
        summarise(`Start day` = t0, n = n(),
                  `LVEF at start` = mean(EF_start),
                  `LVEF at 12 months` = mean(EF_end),
                  `Recovery rate %` = 100 * mean(EF_end >= 50 |
                                               (EF_base - EF_end) < 5),
                  `FIB at start` = mean(FIB_start))
    }))
  })

  output$win_tbl <- renderDT({
    datatable(win() %>% mutate(across(where(is.numeric), ~round(.x, 2))),
              options = list(dom = "t"), rownames = FALSE)
  })

  output$win_plot <- renderPlot({
    d <- win()
    ggplot(d, aes(`Start day`)) +
      geom_line(aes(y = `Recovery rate %`), linewidth = 1.2, colour = "#c0392b") +
      geom_point(aes(y = `Recovery rate %`), size = 3, colour = "#c0392b") +
      geom_line(aes(y = 100 * `FIB at start`), linewidth = 1,
                linetype = 2, colour = "#616161") +
      labs(x = "Heart-failure therapy start day (day)",
           y = "Recovery rate (%) / 100 × FIB",
           title = "Recovery tracks fibrosis, not time (dashed line = FIB)") +
      theme_aic()
  })
}

shinyApp(ui, server)
