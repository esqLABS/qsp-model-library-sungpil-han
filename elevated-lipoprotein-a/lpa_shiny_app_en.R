## =====================================================================
##  lpa_shiny_app.R
##  Elevated Lipoprotein(a) — QSP interactive dashboard
##  10 tabs · drives the 50-ODE model in lpa_mrgsolve_model_en.R
##
##  Elevated Lipoprotein(a) — QSP Dashboard
## ---------------------------------------------------------------------
##  The app is organised around the model's thesis rather than around the
##  compartment list.  Each tab answers ONE question:
##
##   1  Patient Profile   Who is this patient, and what does their KIV-2
##                      repeat number do to everything downstream?
##   2  Production · Assembly     Where in the production chain does each drug act?
##   3  Pharmacokinetics          PK of eight agents on one axis.
##   4  Lp(a) Response      The headline PD curve, four assays at once.
##   5  ★ The Measurement Trap   The tab the model exists for: the SAME patient
##                      reported four ways, and who gets missed.
##   6  Lipids · apoB     LDL-C contamination and the Dahlen correction.
##   7  Inflammation Axis         OxPL -> monocyte -> IL-6 -> CRP, and the loop gain.
##   8  Vascular · Risk     Plaque, vulnerability, and the MR-vs-trial split.
##   9  Aortic Valve      Why starting age, not potency, decides this one.
##  10  Scenario Comparison   All 20 scenarios side by side.
##
##  Requires: shiny, mrgsolve, dplyr, tidyr, ggplot2, DT
##  Run with:  shiny::runApp("lpa_shiny_app.R")
## =====================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

## Load the model from the companion file (keeps ONE source of truth for
## the equations; the app never redefines them).
LPA_MODEL_ONLY <- TRUE
source("lpa_mrgsolve_model_en.R", local = FALSE, echo = FALSE)
## `mod` is now available.

YR <- 365.25

DRUGS <- c("None"                 = "none",
           "High-intensity statin (Statin)"       = "statin",
           "Evolocumab"      = "evolocumab",
           "Pelacarsen ASO"   = "pelacarsen",
           "Olpasiran siRNA"    = "olpasiran",
           "Lepodisiran siRNA"= "lepodisiran",
           "Muvalaplin oral"  = "muvalaplin",
           "Niacin ER"           = "niacin",
           "Obicetrapib"    = "obicetrapib",
           "Ziltivekimab"     = "ziltivekimab",
           "Lipoprotein apheresis"  = "apheresis")

build_ev <- function(drugs, start_day = 0) {
  e <- NULL
  add <- function(x) if (is.null(e)) x else c(e, x)
  for (d in drugs) {
    e <- switch(d,
      statin       = add(ev(amt =   20, cmt = 12, ii = 1,   addl = 99999, time = start_day)),
      evolocumab   = add(ev(amt =  420, cmt =  9, ii = 28,  addl = 999,   time = start_day)),
      pelacarsen   = add(ev(amt =   80, cmt =  1, ii = 28,  addl = 999,   time = start_day)),
      olpasiran    = add(ev(amt =   75, cmt =  4, ii = 84,  addl = 999,   time = start_day)),
      lepodisiran  = add(ev(amt =  400, cmt =  4, ii = 168, addl = 999,   time = start_day)),
      muvalaplin   = add(ev(amt =  240, cmt =  7, ii = 1,   addl = 99999, time = start_day)),
      niacin       = add(ev(amt = 2000, cmt = 14, ii = 1,   addl = 99999, time = start_day)),
      obicetrapib  = add(ev(amt =   10, cmt = 17, ii = 1,   addl = 99999, time = start_day)),
      ziltivekimab = add(ev(amt =   30, cmt = 15, ii = 28,  addl = 999,   time = start_day)),
      e)
  }
  if (is.null(e)) ev(amt = 0, cmt = 1, time = 0) else e
}

THEME <- theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        legend.position = "bottom")

## =====================================================================
##  UI
## =====================================================================
ui <- fluidPage(
  titlePanel("Elevated Lipoprotein(a) QSP Dashboard"),
  tags$p(style = "color:#555;margin-top:-10px",
         "One particle · one rate-limiting step · three arms of action · two scales ",
         tags$em("(one particle, one rate-limiting step, three arms, two scales)")),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Patient"),
      sliderInput("LPA0", "Baseline Lp(a) particles (nmol/L)", 5, 500, 250, step = 5),
      sliderInput("NKIV2", "KIV-2 repeat number (dominantly expressed allele)", 5, 40, 12, step = 1),
      helpText("The smaller the repeat number, the higher the secretion efficiency and the higher the Lp(a), ",
               "and at the same time it is under-reported by mass assays."),
      sliderInput("LDLP0", "Baseline LDL particles (nmol/L)", 500, 2200, 1200, step = 50),
      sliderInput("ERISK", "Vascular risk burden (smoking · hypertension · diabetes)", 0, 1.5, 0, step = 0.1),
      sliderInput("IL6EXO", "Concomitant inflammation IL-6 (pg/mL)", 0, 40, 0, step = 1),
      selectInput("FLDLRFN", "LDL receptor function",
                  c("Normal" = 1, "Heterozygous FH" = 0.5, "Homozygous FH" = 0.05)),

      hr(), h4("Therapy"),
      checkboxGroupInput("drugs", NULL, choices = DRUGS, selected = "pelacarsen"),
      sliderInput("start_yr", "Treatment start (model year)", 0, 50, 0, step = 1),
      sliderInput("horizon", "Simulation duration (years)", 1, 62, 2, step = 1),

      hr(), h4("Assay assumptions"),
      sliderInput("NCAL", "Mass-assay calibrator isoform", 8, 40, 22, step = 1),
      sliderInput("RFREE", "Response per mole of free apo(a) in the traditional apo(a) assay", 1, 6, 1, step = 0.1),
      sliderInput("FCHOL", "Cholesterol fraction of Lp(a) mass", 0.15, 0.40, 0.30, step = 0.01),
      actionButton("go", "Run", class = "btn-primary btn-block")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("1 · Patient Profile",
                 br(), fluidRow(column(6, tableOutput("iso_tbl")),
                                column(6, plotOutput("iso_plot", height = 300))),
                 hr(), verbatimTextOutput("iso_note")),
        tabPanel("2 · Production · Assembly",  br(), plotOutput("prod_plot", height = 560)),
        tabPanel("3 · Pharmacokinetics (PK)",  br(), plotOutput("pk_plot",   height = 560)),
        tabPanel("4 · Lp(a) Response",   br(), plotOutput("pd_plot",   height = 420),
                 hr(), tableOutput("pd_tbl")),
        tabPanel("5 · ★ The Measurement Trap", br(),
                 h4("Reporting the same patient four different ways"),
                 plotOutput("assay_plot", height = 340),
                 hr(), h4("Same TRUE mass of 50 mg/dL, patients differing only in isoform"),
                 DTOutput("miss_tbl")),
        tabPanel("6 · Lipids · apoB",  br(), plotOutput("lipid_plot", height = 480),
                 hr(), verbatimTextOutput("lipid_note")),
        tabPanel("7 · Inflammation Axis",      br(), plotOutput("infl_plot", height = 480),
                 hr(), verbatimTextOutput("loop_note")),
        tabPanel("8 · Vascular · Risk",  br(), plotOutput("risk_plot", height = 480),
                 hr(), verbatimTextOutput("risk_note")),
        tabPanel("9 · Aortic Valve",   br(), plotOutput("valve_plot", height = 480),
                 hr(), verbatimTextOutput("valve_note")),
        tabPanel("10 · Scenario Comparison", br(), DTOutput("scen_tbl"))
      )
    )
  )
)

## =====================================================================
##  SERVER
## =====================================================================
server <- function(input, output, session) {

  pars <- reactive({
    list(LPA0 = input$LPA0, NKIV2 = input$NKIV2, LDLP0 = input$LDLP0,
         ERISK = input$ERISK, IL6EXO = input$IL6EXO,
         FLDLRFN = as.numeric(input$FLDLRFN), NCAL = input$NCAL,
         RFREE = input$RFREE, FCHOL = input$FCHOL,
         APHON = as.numeric("apheresis" %in% input$drugs))
  })

  sim <- eventReactive(input$go, {
    m  <- param(mod, pars())
    e  <- build_ev(setdiff(input$drugs, c("none", "apheresis")),
                   start_day = input$start_yr*YR)
    end <- max(input$horizon, input$start_yr + 1)*YR
    dl  <- if (end > 10*YR) 30 else 1
    mrgsim_df(m, events = e, end = end, delta = dl)
  }, ignoreNULL = FALSE)

  ## ---- tab 1 : patient / isoform ------------------------------------
  output$iso_tbl <- renderTable({
    d <- sim()[1, ]
    data.frame(
      Item = c("KIV-2 repeat number", "apo(a) molecular weight (kDa)", "Particle molecular weight (MDa)",
               "Secretion efficiency SECEFF", "Conversion factor (nmol/L per mg/dL)",
               "Mass-assay bias (x TRUE)", "TRUE mass (mg/dL)",
               "Mass-assay reported value (mg/dL)", "Particle concentration (nmol/L)"),
      Value = c(sprintf("%d", input$NKIV2),
             sprintf("%.0f", d$MWAPOAO/1000),
             sprintf("%.2f", (3.30e6 + d$MWAPOAO)/1e6),
             sprintf("%.3f", d$SECEFFO),
             sprintf("%.3f", d$CONVF),
             sprintf("%.3f", d$EPITB),
             sprintf("%.1f", d$LPA_MASS_T),
             sprintf("%.1f", d$ASSAY_MASS),
             sprintf("%.1f", d$LPA_NMOL))
    )
  })

  output$iso_plot <- renderPlot({
    nk <- 5:40
    conv <- 1e7/(3.30e6 + 14000*(nk + 10) + 35000)
    df <- data.frame(nk,
                     `Conversion factor (nmol/L per mg/dL)` = conv,
                     `Mass-assay bias` = (nk + 10)/(input$NCAL + 10),
                     `Secretion efficiency` = 22^3/(22^3 + nk^3), check.names = FALSE)
    pivot_longer(df, -nk) %>%
      ggplot(aes(nk, value)) +
      geom_line(linewidth = 1.1, colour = "#0D47A1") +
      geom_vline(xintercept = input$NKIV2, linetype = 2, colour = "#B71C1C") +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "KIV-2 repeat number", y = NULL,
           title = "Three things isoform size changes — and by how much") + THEME
  })

  output$iso_note <- renderText({
    d <- sim()[1, ]
    paste0(
      "The chemical conversion factor varies by only about 9%, from 2.52-2.79, across the entire isoform range.\n",
      "Over the same range, antibody bias varies by 150%, from 0.56-1.41.\n",
      "In other words, the problem is the antibody, not unit conversion.\n\n",
      sprintf("This patient: TRUE is %.1f mg/dL, but the mass assay reports %.1f mg/dL,\n",
              d$LPA_MASS_T, d$ASSAY_MASS),
      sprintf("and the particle concentration is %.0f nmol/L. 50 mg/dL threshold -> %s / 125 nmol/L threshold -> %s",
              d$LPA_NMOL,
              ifelse(d$ASSAY_MASS >= 50, "Positive", "Negative (missed)"),
              ifelse(d$LPA_NMOL  >= 125, "Positive", "Negative")))
  })

  ## ---- tab 2 : production / assembly chain --------------------------
  output$prod_plot <- renderPlot({
    sim() %>%
      select(time, MRNA, APOA_ER, APOA_FR, LPA_P, ASSEM, KNOCK, FMUV, TRANS) %>%
      pivot_longer(-time) %>%
      mutate(name = factor(name,
        levels = c("TRANS","MRNA","KNOCK","APOA_ER","APOA_FR","FMUV","ASSEM","LPA_P"),
        labels = c("Transcription TRANS","LPA mRNA","mRNA degradation multiplier KNOCK","ER apo(a) pool",
                   "Free apo(a)","Assembly residual fraction FMUV","Assembly flux ASSEM","Lp(a) particle"))) %>%
      ggplot(aes(time/YR, value)) +
      geom_line(linewidth = 1, colour = "#00695C") +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "Time (years)", y = NULL,
           title = "Drug action sites along the production chain") + THEME
  })

  ## ---- tab 3 : PK ----------------------------------------------------
  output$pk_plot <- renderPlot({
    sim() %>%
      select(time, CPELL, CMUV, CEVO, CSTA, CZIL, COBI, SIR_RISC) %>%
      pivot_longer(-time) %>% filter(value > 0) %>%
      mutate(name = recode(name,
        CPELL = "Pelacarsen liver (mg/L)", CMUV = "Muvalaplin plasma (mg/L)",
        CEVO = "Evolocumab plasma (mg/L)", CSTA = "Statin plasma (mg/L)",
        CZIL = "Ziltivekimab plasma (mg/L)", COBI = "Obicetrapib plasma (mg/L)",
        SIR_RISC = "siRNA RISC loading (mg)")) %>%
      ggplot(aes(time/YR, value)) +
      geom_line(linewidth = 0.9, colour = "#1B5E20") +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "Time (years)", y = NULL,
           title = "Pharmacokinetics — the RISC compartment's half-life creates PD persistence") + THEME
  })

  ## ---- tab 4 : Lp(a) response ---------------------------------------
  output$pd_plot <- renderPlot({
    sim() %>%
      ggplot(aes(time/YR, LPA_NMOL)) +
      geom_hline(yintercept = c(75, 125, 180), linetype = 3, colour = "grey50") +
      annotate("text", x = 0, y = c(75, 125, 180), hjust = 0, vjust = -0.4,
               size = 3, colour = "grey40",
               label = c("75 upper limit of normal", "125 at risk", "180 very high")) +
      geom_line(linewidth = 1.2, colour = "#B71C1C") +
      labs(x = "Time (years)", y = "Lp(a) particles (nmol/L)",
           title = "Lp(a) particle concentration") + THEME
  })

  output$pd_tbl <- renderTable({
    d <- sim(); b <- d[1, ]; e <- d[nrow(d), ]
    data.frame(
      Metric = c("Lp(a) (nmol/L)", "Lp(a) TRUE mass (mg/dL)",
               "Mass-assay reported value (mg/dL)", "Traditional apo(a) assay (nmol/L)",
               "Free apo(a) (nmol/L)", "Absolute decrease (nmol/L)"),
      Baseline = c(b$LPA_NMOL, b$LPA_MASS_T, b$ASSAY_MASS, b$ASSAY_APOA,
               b$ASSAY_APOA - b$LPA_NMOL, NA),
      End = c(e$LPA_NMOL, e$LPA_MASS_T, e$ASSAY_MASS, e$ASSAY_APOA,
               e$ASSAY_APOA - e$LPA_NMOL, NA),
      PctChange = c(100*(e$LPA_NMOL/b$LPA_NMOL - 1),
                 100*(e$LPA_MASS_T/b$LPA_MASS_T - 1),
                 100*(e$ASSAY_MASS/b$ASSAY_MASS - 1),
                 100*(e$ASSAY_APOA/b$ASSAY_APOA - 1),
                 100*((e$ASSAY_APOA - e$LPA_NMOL)/(b$ASSAY_APOA - b$LPA_NMOL) - 1),
                 e$LPA_NMOL - b$LPA_NMOL)
    )
  }, digits = 1)

  ## ---- tab 5 : the measurement trap ---------------------------------
  output$assay_plot <- renderPlot({
    sim() %>%
      transmute(time,
                `Particle molar assay (nmol/L)` = ASSAY_NMOL,
                `intact-Lp(a) assay (nmol/L)` = ASSAY_INTACT,
                `Traditional apo(a) assay (nmol/L)` = ASSAY_APOA,
                `Mass assay × conversion (nmol/L equiv.)` = ASSAY_MASS*CONVF) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time/YR, value, colour = name)) +
      geom_line(linewidth = 1.1) +
      scale_colour_manual(values = c("#FF6F00", "#00838F", "#6A1B9A", "#B71C1C")) +
      labs(x = "Time (years)", y = "Reported value (nmol/L basis)", colour = NULL,
           title = "Four assays, one patient — the divergence is largest with assembly inhibitors") +
      THEME
  })

  output$miss_tbl <- renderDT({
    nk <- c(6, 8, 12, 16, 22, 28, 35, 40)
    conv <- 1e7/(3.30e6 + 14000*(nk + 10) + 35000)
    true_nm <- 50*conv
    bias <- (nk + 10)/(input$NCAL + 10)
    data.frame(
      `KIV-2` = nk,
      `apo(a) MW (kDa)` = round((14000*(nk + 10) + 35000)/1000),
      `TRUE mass (mg/dL)` = 50,
      `TRUE particles (nmol/L)` = round(true_nm, 1),
      `Conversion factor` = round(conv, 3),
      `Antibody bias` = round(bias, 3),
      `Reported mass (mg/dL)` = round(50*bias, 1),
      `Mass 50 threshold` = ifelse(50*bias >= 50, "Positive", "Negative (missed)"),
      `Molar 125 threshold` = ifelse(true_nm >= 125, "Positive", "Negative"),
      check.names = FALSE
    )
  }, options = list(dom = "t", pageLength = 10), rownames = FALSE)

  ## ---- tab 6 : lipids / apoB ----------------------------------------
  output$lipid_plot <- renderPlot({
    sim() %>%
      select(time, LDLC_TRUE, LDLC_MEAS, LDLC_CORR, LPAC, APOB_TOT, LPA_APOB_PCT) %>%
      pivot_longer(-time) %>%
      mutate(name = recode(name,
        LDLC_TRUE = "TRUE LDL-C (mg/dL)", LDLC_MEAS = "Reported LDL-C (mg/dL)",
        LDLC_CORR = "Dahlen-corrected LDL-C (mg/dL)", LPAC = "Lp(a)-C (mg/dL)",
        APOB_TOT = "Total apoB (mg/dL)", LPA_APOB_PCT = "Lp(a) share of apoB particles (%)")) %>%
      ggplot(aes(time/YR, value)) +
      geom_line(linewidth = 1, colour = "#AD1457") +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "Time (years)", y = NULL,
           title = "What is actually inside 'LDL-C'") + THEME
  })

  output$lipid_note <- renderText({
    b <- sim()[1, ]
    paste0(
      sprintf("Of the reported LDL-C of %.1f mg/dL, %.1f mg/dL (%.0f%%) is not LDL but Lp(a) cholesterol.\n",
              b$LDLC_MEAS, b$LPAC, 100*b$LPAC/b$LDLC_MEAS),
      sprintf("Performing the Dahlen correction with the (biased) mass-assay value gives %.1f mg/dL,\n",
              b$LDLC_CORR),
      sprintf("which is off from the TRUE LDL-C of %.1f mg/dL by %+.1f mg/dL.\n", b$LDLC_TRUE,
              b$LDLC_CORR - b$LDLC_TRUE),
      sprintf("Lp(a) is only %.1f%% of apoB particles, yet accounts for %.1f%% of intimal retention flux.",
              b$LPA_APOB_PCT, b$RETSHARE))
  })

  ## ---- tab 7 : inflammation -----------------------------------------
  output$infl_plot <- renderPlot({
    sim() %>%
      select(time, OXPL, MONO, IL6, IL6EFF, CRP, FIL6) %>%
      pivot_longer(-time) %>%
      mutate(name = recode(name,
        OXPL = "OxPL-apoB (relative)", MONO = "Trained monocyte index",
        IL6 = "IL-6 (pg/mL)", IL6EFF = "Signalling-competent IL-6 (pg/mL)",
        CRP = "hsCRP (mg/L)", FIL6 = "IL-6 fold-effect on LPA transcription")) %>%
      ggplot(aes(time/YR, value)) +
      geom_line(linewidth = 1, colour = "#E65100") +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "Time (years)", y = NULL,
           title = "OxPL -> monocyte -> IL-6 -> LPA promoter") + THEME
  })

  output$loop_note <- renderText({
    paste0(
      "The feedforward loop is real but weak. The open-loop gain g is about 0.02 in the normal inflammatory state,\n",
      "and amplification is only 1/(1-g) = 1.02.\n\n",
      "What is not weak is the 'basal IL-6 tone': basal IL-6 alone accounts for about 25% of LPA transcription,\n",
      "which is why IL-6 blockade lowers Lp(a) by about 25% in non-inflammatory patients\n",
      "and by about 35-40% in rheumatoid arthritis — same equation, only IL-6 changed.")
  })

  ## ---- tab 8 : plaque / risk ----------------------------------------
  output$risk_plot <- renderPlot({
    sim() %>%
      select(time, INT_LPA, FOAM, NECRO, PLAQUE, VULN, HR_MACE, HR_SLOW, HR_FAST) %>%
      pivot_longer(-time) %>%
      mutate(name = recode(name,
        INT_LPA = "Intimal retained Lp(a)", FOAM = "Foam-cell burden",
        NECRO = "Necrotic core index", PLAQUE = "Plaque volume PAV (%)",
        VULN = "Vulnerability index", HR_MACE = "MACE hazard ratio HR",
        HR_SLOW = "HR — slow component (atheroma)", HR_FAST = "HR — fast component (vulnerability)")) %>%
      ggplot(aes(time/YR, value)) +
      geom_line(linewidth = 1, colour = "#3E2723") +
      facet_wrap(~name, scales = "free_y", ncol = 4) +
      labs(x = "Time (years)", y = NULL,
           title = "Slow component and fast component — why MR and trials diverge") + THEME
  })

  output$risk_note <- renderText({
    paste0(
      "Risk is carried by two components with different time constants.\n",
      "  Slow: atheroma burden, half-life about 14 years → about 70% of excess risk. What Mendelian randomisation sees.\n",
      "  Fast: OxPL/IL-6-mediated vulnerability, half-life about 2 months → about 30%. What a 5-year trial can move.\n\n",
      "As a result, the ratio of lifetime-exposure effect / 5-year-trial effect comes out 'automatically' at about 2.6-2.7-fold.\n",
      "This ratio is not a value entered as a correction factor but one derived from the two time constants.")
  })

  ## ---- tab 9 : aortic valve -----------------------------------------
  output$valve_plot <- renderPlot({
    sim() %>%
      select(time, ATXV, LYSOPA, VIC_OST, VCALC, AVA, MGRAD) %>%
      pivot_longer(-time) %>%
      mutate(name = recode(name,
        ATXV = "Valvular autotaxin", LYSOPA = "Valvular LysoPA",
        VIC_OST = "Osteogenic VIC fraction", VCALC = "Valvular calcium AVC (AU)",
        AVA = "Aortic valve area (cm2)", MGRAD = "Mean pressure gradient (mmHg)")) %>%
      ggplot(aes(time/YR, value)) +
      geom_line(linewidth = 1, colour = "#283593") +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "Time (years)", y = NULL,
           title = "The valve — the one endpoint where start timing matters more than potency") + THEME
  })

  output$valve_note <- renderText({
    paste0(
      "Valve calcification, once it crosses the threshold (Hill n=4), self-perpetuates regardless of Lp(a).\n",
      "In a 62-year simulation (Lp(a) 250 nmol/L, pelacarsen):\n",
      "  Start at age 30 → AVC about 206 AU, AVA 2.45 cm2   (87% of calcium prevented)\n",
      "  Start at age 60 → AVC about 1327 AU, AVA 0.93 cm2  (20 years of treatment reduces calcium by 15%)\n",
      "  No treatment    → AVC about 1558 AU, AVA 0.82 cm2\n\n",
      "This is why SEAS · ASTRONOMER · SALTIRE were all negative, and at the same time it is a\n",
      "prediction that the valve indication cannot succeed with a secondary-prevention design no matter how potent the drug.")
  })

  ## ---- tab 10 : scenario comparison ---------------------------------
  output$scen_tbl <- renderDT({
    withProgress(message = "Running 20 scenarios...", value = 0, {
      rows <- lapply(seq_along(sc), function(i) {
        incProgress(1/length(sc))
        s <- sc[[i]]
        m <- if (length(s$p)) param(mod, s$p) else mod
        o <- mrgsim_df(m, events = s$ev, end = 2*YR, delta = 7)
        b <- o[1, ]; e <- o[nrow(o), ]
        data.frame(
          Scenario = names(sc)[i],
          `Baseline nmol/L` = round(b$LPA_NMOL, 0),
          `2-year nmol/L` = round(e$LPA_NMOL, 1),
          `% change` = round(100*(e$LPA_NMOL/b$LPA_NMOL - 1), 1),
          `Absolute Δ nmol/L` = round(e$LPA_NMOL - b$LPA_NMOL, 1),
          `LDL-C (reported) %` = round(100*(e$LDLC_MEAS/b$LDLC_MEAS - 1), 1),
          `apoB %` = round(100*(e$APOB_TOT/b$APOB_TOT - 1), 1),
          `hsCRP %` = round(100*(e$CRP/b$CRP - 1), 1),
          check.names = FALSE)
      })
      do.call(rbind, rows)
    })
  }, options = list(dom = "t", pageLength = 25, scrollX = TRUE), rownames = FALSE)
}

shinyApp(ui, server)
