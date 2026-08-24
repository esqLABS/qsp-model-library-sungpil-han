## =====================================================================
##  rili_shiny_app.R
##  Radiation-Induced Lung Injury (RILI) — interactive QSP dashboard
##  Radiation-Induced Lung Injury QSP Dashboard
##
##  The app is built around the one thing this model does that a
##  mean-dose model cannot: the treatment plan is a DOSE-VOLUME
##  HISTOGRAM, so the user edits a vector of six volumes, not a dose.
##  Tab 2 exists to let you hold the mean lung dose fixed while you
##  reshape the histogram, and watch pneumonitis and fibrosis move in
##  opposite directions.
##
##  Ten tabs
##    1  Patient & Plan          DVH editor, MLD/V5/V20/V40
##    2  Iso-MLD experiment         same MLD, different shape
##    3  Pharmacokinetics            8 drugs
##    4  Fast Loop                DAMP → cytokine → oedema
##    5  Slow Loop                TGF-β → myofibroblast → collagen
##    6  Bistable Switch    phase plane + separatrix
##    7  Pulmonary Function          DLCO, FVC, per-bin map
##    8  Clinical Endpoints  NTCP, CTCAE grade, latency
##    9  Therapeutic Ratio           TCP vs NTCP vs UCP
##   10  Scenario Comparison   the 40 canned scenarios
##
##  RUN
##    shiny::runApp("rili_shiny_app.R")
##  The model file rili_mrgsolve_model_en.R must sit in the same directory.
## =====================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)

MODEL_DIR <- "."
mod <- mread_cache("rili_mrgsolve_model", MODEL_DIR)

## ---- plan library ---------------------------------------------------
PLANS <- list(
  "IMRT 60 Gy / 30 fx"       = list(DRX = 60, NFX = 30,
                                    v = c(.40, .20, .12, .11, .09, .08)),
  "3D-CRT 60 Gy / 30 fx"     = list(DRX = 60, NFX = 30,
                                    v = c(.52, .13, .08, .08, .09, .10)),
  "Proton 60 Gy / 30 fx"     = list(DRX = 60, NFX = 30,
                                    v = c(.58, .16, .09, .07, .055, .045)),
  "IMRT 74 Gy / 37 fx"       = list(DRX = 74, NFX = 37,
                                    v = c(.34, .20, .13, .12, .11, .10)),
  "Hypofx 60 Gy / 8 fx"      = list(DRX = 60, NFX = 8,
                                    v = c(.700, .130, .070, .045, .033,
                                          .022)),
  "SBRT 54 Gy / 3 fx"        = list(DRX = 54, NFX = 3,
                                    v = c(.860, .075, .032, .018, .010,
                                          .005)),
  "Low-dose bath (V5 high)"  = list(DRX = 60, NFX = 30,
                                    v = c(.06, .55, .29, .05, .03, .02)),
  "Focal hot spot (V40 high)" = list(DRX = 60, NFX = 30,
                                     v = c(.600, .110, .050, .050, .070,
                                           .120))
)
EDGEF <- c(0, .08, .24, .44, .66, .88, 1.05)
REPF  <- (EDGEF[-7] + EDGEF[-1]) / 2

mld_of <- function(v, DRX) sum(v / sum(v) * DRX * REPF)
vx_of <- function(v, DRX, x) {
  vv <- v / sum(v)
  tot <- 0
  for (j in seq_len(6)) {
    lo <- DRX * EDGEF[j]; hi <- DRX * EDGEF[j + 1]
    if (hi <= x) next
    tot <- tot + if (lo >= x) vv[j] else vv[j] * (hi - x) / (hi - lo)
  }
  tot
}
geud_of <- function(v, DRX, a) {
  vv <- v / sum(v); sum(vv * (DRX * REPF) ^ a) ^ (1 / a)
}

THEME <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 12))

## =====================================================================
##  UI
## =====================================================================
ui <- fluidPage(
  titlePanel(HTML(paste0(
    "<b>Radiation-Induced Lung Injury QSP Model</b> ",
    "<span style='font-size:70%'>Radiation-Induced Lung Injury &middot; ",
    "79 ODEs &middot; 6 DVH bins &times; 10 states</span>"))),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectInput("plan", "Plan preset",
                  choices = names(PLANS), selected = "IMRT 60 Gy / 30 fx"),
      fluidRow(
        column(6, numericInput("DRX", "Prescription dose D_RX (Gy)", 60, 20, 90,
                               step = 2)),
        column(6, numericInput("NFX", "Number of fractions NFX", 30, 1, 45, step = 1))
      ),
      tags$hr(),
      tags$b("Dose-volume histogram (DVH)"),
      helpText(HTML(paste0(
        "Fraction of lung volume per bin. Automatically normalised even if the total isn't 1.<br>",
        "<i>These six numbers are this model's 'prescription'.</i>"))),
      fluidRow(
        column(6, numericInput("v1", "bin1 (~4% D_RX)", .40, 0, 1,
                               step = .01)),
        column(6, numericInput("v2", "bin2 (~16%)", .20, 0, 1, step = .01))
      ),
      fluidRow(
        column(6, numericInput("v3", "bin3 (~34%)", .12, 0, 1, step = .01)),
        column(6, numericInput("v4", "bin4 (~55%)", .11, 0, 1, step = .01))
      ),
      fluidRow(
        column(6, numericInput("v5", "bin5 (~77%)", .09, 0, 1, step = .01)),
        column(6, numericInput("v6", "bin6 (~97%)", .08, 0, 1, step = .01))
      ),
      tags$hr(),
      tags$b("Patient reserve"),
      sliderInput("DLCO0", "Baseline DLCO (%pred)", 30, 110, 85, step = 5),
      sliderInput("FVC0", "Baseline FVC (%pred)", 40, 120, 95, step = 5),
      checkboxInput("emph", "Emphysema perfusion weighting",
                    FALSE),
      tags$hr(),
      tags$b("Drugs"),
      tags$div(style = "font-size:88%",
        tags$span(style = "color:#B8860B", "● Fast loop")),
      fluidRow(
        column(6, numericInput("rdex", "Prednisolone mg/d", 0, 0, 100,
                               step = 10)),
        column(6, numericInput("tdex", "Start day (day)", 56, 0, 400,
                               step = 7))
      ),
      checkboxInput("durv", "Durvalumab 1500 mg q4w (PACIFIC)", FALSE),
      tags$div(style = "font-size:88%",
        tags$span(style = "color:#9B4FA8", "● Slow loop")),
      fluidRow(
        column(6, numericInput("rpir", "Pirfenidone mg/d", 0, 0, 2403,
                               step = 801)),
        column(6, numericInput("tpir", "Start day (day)", 0, 0, 400,
                               step = 7))
      ),
      fluidRow(
        column(6, numericInput("rnin", "Nintedanib mg/d", 0, 0, 300,
                               step = 150)),
        column(6, numericInput("race", "Lisinopril mg/d", 0, 0, 40,
                               step = 10))
      ),
      tags$div(style = "font-size:88%",
        tags$span(style = "color:#E08B2E", "● Initial damage")),
      checkboxInput("ami", "Amifostine", FALSE),
      sliderInput("amidel", "Interval before irradiation (min before beam)", 0, 90, 20,
                  step = 5),
      helpText(HTML(paste0("WR-1065 half-life 8 min — the protection factor decays as ",
                           "exp(&minus;ln2&middot;&Delta;t/8)."))),
      tags$hr(),
      sliderInput("tend", "Simulation duration (days)", 200, 1095, 730,
                  step = 30)
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs", type = "tabs",
        tabPanel("1 Patient & Plan", br(),
                 fluidRow(column(7, plotOutput("dvhPlot", height = 300)),
                          column(5, tableOutput("planTab"))),
                 hr(), plotOutput("bedPlot", height = 260)),
        tabPanel("2 Iso-MLD experiment", br(),
                 htmlOutput("isoNote"),
                 plotOutput("isoPlot", height = 380),
                 hr(), tableOutput("isoTab")),
        tabPanel("3 Pharmacokinetics", br(),
                 plotOutput("pkPlot", height = 460),
                 hr(), htmlOutput("protNote")),
        tabPanel("4 Fast Loop", br(),
                 plotOutput("fastPlot", height = 500),
                 hr(), htmlOutput("gainNote")),
        tabPanel("5 Slow Loop", br(),
                 plotOutput("slowPlot", height = 500)),
        tabPanel("6 Bistable Switch", br(),
                 htmlOutput("switchNote"),
                 plotOutput("phasePlot", height = 420),
                 hr(), tableOutput("fpTab")),
        tabPanel("7 Pulmonary Function", br(),
                 plotOutput("pftPlot", height = 320),
                 hr(), plotOutput("binMap", height = 300)),
        tabPanel("8 Clinical Endpoints", br(),
                 fluidRow(column(6, plotOutput("ntcpPlot", height = 320)),
                          column(6, plotOutput("gradePlot", height = 320))),
                 hr(), tableOutput("endTab")),
        tabPanel("9 Therapeutic Ratio", br(),
                 plotOutput("trPlot", height = 400),
                 hr(), tableOutput("trTab")),
        tabPanel("10 Scenario Comparison", br(),
                 checkboxGroupInput(
                   "scen", "Scenarios to compare",
                   choices = names(PLANS),
                   selected = c("SBRT 54 Gy / 3 fx", "IMRT 60 Gy / 30 fx",
                                "IMRT 74 Gy / 37 fx"),
                   inline = TRUE),
                 plotOutput("scenPlot", height = 420),
                 hr(), tableOutput("scenTab"))
      )
    )
  )
)

## =====================================================================
##  SERVER
## =====================================================================
server <- function(input, output, session) {

  ## keep the DVH boxes in sync with the preset selector
  observeEvent(input$plan, {
    p <- PLANS[[input$plan]]
    updateNumericInput(session, "DRX", value = p$DRX)
    updateNumericInput(session, "NFX", value = p$NFX)
    for (j in 1:6)
      updateNumericInput(session, paste0("v", j), value = p$v[j])
  })

  vvec <- reactive({
    v <- c(input$v1, input$v2, input$v3, input$v4, input$v5, input$v6)
    v[is.na(v) | v < 0] <- 0
    if (sum(v) <= 0) v <- c(1, 0, 0, 0, 0, 0)
    v / sum(v)
  })

  parlist <- reactive({
    v <- vvec()
    pf <- if (isTRUE(input$emph)) c(1, 1, .8, .6, .5, .4) else rep(1, 6)
    lst <- list(DRX = input$DRX, NFX = input$NFX,
                V1 = v[1], V2 = v[2], V3 = v[3],
                V4 = v[4], V5 = v[5], V6 = v[6],
                PF1 = pf[1], PF2 = pf[2], PF3 = pf[3],
                PF4 = pf[4], PF5 = pf[5], PF6 = pf[6],
                DLCO0 = input$DLCO0, FVC0 = input$FVC0,
                RDEX = input$rdex, TDEX0 = if (input$rdex > 0)
                  input$tdex else 1e9,
                RPIR = input$rpir, TPIR0 = if (input$rpir > 0)
                  input$tpir else 1e9,
                TPIR1 = if (input$rpir > 0) input$tpir + 365 else 1e9,
                RNIN = input$rnin, TNIN0 = if (input$rnin > 0) 0 else 1e9,
                TNIN1 = if (input$rnin > 0) 365 else 1e9,
                RACE = input$race, TACE0 = if (input$race > 0) 0 else 1e9,
                TACE1 = if (input$race > 0) 730 else 1e9,
                AMION = as.numeric(isTRUE(input$ami)),
                AMIDEL = input$amidel)
    lst
  })

  dosing <- reactive({
    if (isTRUE(input$durv))
      ev(amt = 1500, cmt = "DURC", time = 56, ii = 28, addl = 12)
    else ev(amt = 0, cmt = "DURC", time = 0)
  })

  sim <- reactive({
    mod %>% param(parlist()) %>% ev(dosing()) %>%
      mrgsim(end = input$tend, delta = 1) %>% as_tibble()
  })

  simplan <- function(nm, extra = list()) {
    p <- PLANS[[nm]]
    v <- p$v / sum(p$v)
    pr <- c(list(DRX = p$DRX, NFX = p$NFX, V1 = v[1], V2 = v[2],
                 V3 = v[3], V4 = v[4], V5 = v[5], V6 = v[6],
                 DLCO0 = input$DLCO0, FVC0 = input$FVC0), extra)
    mod %>% param(pr) %>% mrgsim(end = input$tend, delta = 2) %>%
      as_tibble() %>% mutate(plan = nm)
  }

  tcourse <- reactive(input$NFX / 5 * 7)

  ## ---- 1  plan ------------------------------------------------------
  output$dvhPlot <- renderPlot({
    v <- vvec(); D <- input$DRX
    d <- tibble(bin = factor(1:6), dose = D * REPF, vol = 100 * v)
    ggplot(d, aes(dose, vol)) +
      geom_col(fill = "#4A90D9", width = D * .08) +
      geom_text(aes(label = sprintf("%.0f%%", vol)), vjust = -.4,
                size = 3.4) +
      geom_vline(xintercept = 42, linetype = "dashed",
                 colour = "#9B4FA8", linewidth = .9) +
      annotate("text", x = 42, y = max(100 * v) * .95,
               label = "Fibrotic threshold (42 Gy)",
               hjust = -.05, size = 3.2, colour = "#9B4FA8") +
      labs(title = "Dose-volume histogram (DVH)",
           subtitle = "Purple line = the fibrotic-conversion threshold this model computes",
           x = "Bin representative dose (Gy)", y = "Lung volume (%)") + THEME
  })

  output$planTab <- renderTable({
    v <- vvec(); D <- input$DRX
    tibble(
      `Metric` = c("Mean lung dose MLD (Gy)", "V5 (%)", "V20 (%)",
                          "V40 (%)", "gEUD a=1 (Gy)", "gEUD a=0.7 (Gy)",
                          "Fraction dose d (Gy)", "Treatment duration (days)",
                          "Tumour BED10 (Gy)", "Maximum bin lung BED3 (Gy)"),
      `Value` = c(
        sprintf("%.2f", mld_of(v, D)),
        sprintf("%.1f", 100 * vx_of(v, D, 5)),
        sprintf("%.1f", 100 * vx_of(v, D, 20)),
        sprintf("%.1f", 100 * vx_of(v, D, 40)),
        sprintf("%.2f", geud_of(v, D, 1)),
        sprintf("%.2f", geud_of(v, D, 0.7)),
        sprintf("%.2f", D / input$NFX),
        sprintf("%.0f", tcourse()),
        sprintf("%.1f", D * (1 + (D / input$NFX) / 10)),
        sprintf("%.1f", D * REPF[6] * (1 + (D * REPF[6] / input$NFX) / 3))))
  }, striped = TRUE)

  output$bedPlot <- renderPlot({
    v <- vvec(); D <- input$DRX; n <- input$NFX
    db <- D * REPF; dfx <- db / n
    d <- tibble(bin = factor(1:6), dose = db,
                BED3 = db * (1 + dfx / 3),
                BED10 = db * (1 + dfx / 10)) %>%
      pivot_longer(c(BED3, BED10))
    ggplot(d, aes(dose, value, colour = name)) +
      geom_line(linewidth = 1) + geom_point(size = 2.4) +
      scale_colour_manual(values = c(BED3 = "#C0392B", BED10 = "#4CA85F"),
                          labels = c("Lung BED₃ (α/β=3)",
                                     "Tumour BED₁₀ (α/β=10)")) +
      labs(title = "Same physical dose, different biologically equivalent dose",
           subtitle = paste0("Fraction dose ", sprintf("%.2f", D / n),
                             " Gy — the larger the fraction dose, the more the two curves diverge"),
           x = "Physical dose (Gy)", y = "BED (Gy)", colour = NULL) + THEME
  })

  ## ---- 2  iso-MLD experiment ---------------------------------------
  output$isoNote <- renderUI({
    HTML(paste0(
      "<div style='background:#F2F8FD;padding:10px;border-left:4px solid ",
      "#4A90D9;margin-bottom:8px'><b>This tab is this model's central claim.</b> ",
      "The two plans below have <b>almost the same mean lung dose</b> (15.5 vs 15.4 Gy) ",
      "yet V40 differs by 3.8-fold. Pneumonitis and fibrosis rank the two plans ",
      "<b>in exactly opposite order</b> — because pneumonitis is a volume-weighted ",
      "<i>sum</i> over a parallel organ, while fibrosis is a local <i>threshold</i>.</div>"))
  })

  isoSim <- reactive({
    bind_rows(simplan("Low-dose bath (V5 high)"),
              simplan("Focal hot spot (V40 high)"))
  })

  output$isoPlot <- renderPlot({
    d <- isoSim() %>%
      select(time, plan, NTCP, VFIB, PNIE, FIB) %>%
      mutate(NTCP = 100 * NTCP, VFIB = 100 * VFIB) %>%
      pivot_longer(c(NTCP, VFIB, PNIE, FIB))
    lab <- c(NTCP = "NTCP grade≥2 (%)", VFIB = "Fibrotic-conversion volume (%)",
             PNIE = "Pneumonitis index PNIe", FIB = "Fibrosis score FIB")
    ggplot(d, aes(time, value, colour = plan)) +
      geom_line(linewidth = 1) +
      facet_wrap(~ name, scales = "free_y",
                 labeller = labeller(name = lab)) +
      scale_colour_manual(values = c("#4A90D9", "#C0392B")) +
      labs(title = "Two plans with the same mean dose but different histogram shapes",
           x = "Time (days)", y = NULL, colour = NULL) + THEME
  })

  output$isoTab <- renderTable({
    isoSim() %>% group_by(plan) %>%
      summarise(MLD = first(MLD), V5 = 100 * first(V5G),
                V40 = 100 * first(V40G),
                `NTCP peak %` = 100 * max(NTCP),
                `Vfib final %` = 100 * last(VFIB),
                `DLCO final` = last(DLCO), .groups = "drop")
  }, digits = 2, striped = TRUE)

  ## ---- 3  PK --------------------------------------------------------
  output$pkPlot <- renderPlot({
    d <- sim() %>% select(time, CDEXo, CPIRo, CNINo, CDURo, CACEo) %>%
      pivot_longer(-time) %>% filter(value > 1e-9)
    if (!nrow(d)) {
      return(ggplot() + annotate("text", 0, 0, size = 5,
        label = "No drug selected") +
        theme_void())
    }
    lab <- c(CDEXo = "Prednisolone (mg/L)", CPIRo = "Pirfenidone (mg/L)",
             CNINo = "Nintedanib (mg/L)", CDURo = "Durvalumab (mg/L)",
             CACEo = "Lisinopril (mg/L)")
    ggplot(d, aes(time, value, colour = name)) +
      geom_line(linewidth = .9) +
      geom_vline(xintercept = tcourse(), linetype = "dotted") +
      facet_wrap(~ name, scales = "free_y",
                 labeller = labeller(name = lab)) +
      labs(title = "Pharmacokinetics (dotted line = end of radiotherapy)",
           x = "Time (days)", y = "Concentration", colour = NULL) +
      THEME + theme(legend.position = "none")
  })

  output$protNote <- renderUI({
    if (!isTRUE(input$ami))
      return(HTML("<i>Amifostine is turned off.</i>"))
    pf <- 0.65 * exp(-log(2) * input$amidel / 8)
    HTML(sprintf(paste0(
      "<div style='background:#FFF6EC;padding:10px;border-left:4px solid ",
      "#E08B2E'><b>Temporal coincidence.</b> Administered ",
      "%d min before irradiation → protection factor <b>%.3f</b> (maximum 0.650). ",
      "Because WR-1065's half-life is 8 minutes, protection nearly vanishes with a 30–60 minute delay. ",
      "This is this model's mechanistic explanation for why amifostine trial ",
      "results have been conflicting.</div>"), input$amidel, pf))
  })

  ## ---- 4  fast loop -------------------------------------------------
  output$fastPlot <- renderPlot({
    s <- sim()
    d <- s %>% select(time, DOM_6, CYT_6, PRM_6, SUR_6, CYTS, CYT_1) %>%
      pivot_longer(-time)
    lab <- c(DOM_6 = "Lethal-damage pool DOOM (bin6)",
             CYT_6 = "Local cytokine (bin6)",
             PRM_6 = "Alveolar oedema PERM (bin6)",
             SUR_6 = "Surfactant SURF (bin6)",
             CYTS  = "Systemic cytokine",
             CYT_1 = "Cytokine (bin1, low dose)")
    pk <- s$time[which.max(s$PNIE)]
    ggplot(d, aes(time, value)) +
      geom_line(linewidth = .9, colour = "#D6A81E") +
      geom_vline(xintercept = tcourse(), linetype = "dotted") +
      geom_vline(xintercept = pk, linetype = "dashed",
                 colour = "#C7538C") +
      facet_wrap(~ name, scales = "free_y",
                 labeller = labeller(name = lab)) +
      labs(title = "Fast loop — DAMP → NF-κB → cytokine → oedema",
           subtitle = paste0("Dotted line = end of RT (day ", round(tcourse()),
                             "), pink dashed line = pneumonitis index peak (day ",
                             round(pk), ")"),
           x = "Time (days)", y = NULL) + THEME
  })

  output$gainNote <- renderUI({
    g <- if (isTRUE(input$durv)) 1.35 else 1.00
    budget <- g * (0.025 + 2.0 * 0.015 / 1.0) + 0.10 * 0.05 / 0.30
    HTML(sprintf(paste0(
      "<div style='background:#FFFBEA;padding:10px;border-left:4px solid ",
      "#D6A81E'><b>Gain budget.</b> ",
      "GIMM(KAMP + KDAMP·KDCYT/KMCYT) + KSYSIN·KSOUT/KCLS = ",
      "<b>%.4f</b> &lt; KCE = 0.150. Because the fast loop is below the critical threshold, the pneumonitis ",
      "self-resolves. Durvalumab doesn't add a term but raises the gain GIMM by ",
      "%.2f-fold, extending the effective time constant from 12.8 days to %.1f days — ",
      "which is why the excess risk <i>grows</i> with mean lung dose.</div>"),
      budget, g, 1 / (0.15 - budget)))
  })

  ## ---- 5  slow loop -------------------------------------------------
  output$slowPlot <- renderPlot({
    d <- sim() %>%
      select(time, TGF_6, MFB_6, COL_6, XCL_6, TGF_4, COL_4) %>%
      pivot_longer(-time)
    lab <- c(TGF_6 = "Active TGF-β1 (bin6)", MFB_6 = "Myofibroblasts (bin6)",
             COL_6 = "Soluble collagen COL (bin6)",
             XCL_6 = "Cross-linked collagen XCOL (bin6)",
             TGF_4 = "Active TGF-β1 (bin4)",
             COL_4 = "Soluble collagen COL (bin4)")
    ggplot(d, aes(time, value)) +
      geom_line(linewidth = .9, colour = "#9B4FA8") +
      geom_hline(data = tibble(name = c("COL_6", "COL_4"), y = 1.586),
                 aes(yintercept = y), linetype = "dashed",
                 colour = "#C0392B") +
      facet_wrap(~ name, scales = "free_y",
                 labeller = labeller(name = lab)) +
      labs(title = "Slow loop — TGF-β1 → myofibroblasts → collagen → stiffness",
           subtitle = paste0("Red dashed line = separatrix COL 1.586. ",
                             "Crossing this line, there is no return."),
           x = "Time (days)", y = NULL) + THEME
  })

  ## ---- 6  bistable switch -------------------------------------------
  output$switchNote <- renderUI({
    HTML(paste0(
      "<div style='background:#FAF0FB;padding:10px;border-left:4px solid ",
      "#9B4FA8'><b>The slow loop is bistable.</b> The three fixed points are ",
      "computed algebraically from the parameters, and the simulation reproduces them. ",
      "Fibrosis is not a <i>continuous dose-response</i> but a <b>threshold</b> — ",
      "at 40 Gy the collagen excess is 0.09, at 42 Gy it is 2.85. ",
      "This is why, clinically, radiation fibrosis is <b>sharply bounded</b> ",
      "by the isodose line.</div>"))
  })

  output$phasePlot <- renderPlot({
    s <- sim()
    d <- bind_rows(
      tibble(bin = "bin6 (highest dose)", TGF = s$TGF_6,
             COLT = s$COL_6 + s$XCL_6, time = s$time),
      tibble(bin = "bin5", TGF = s$TGF_5, COLT = s$COL_5 + s$XCL_5,
             time = s$time),
      tibble(bin = "bin4", TGF = s$TGF_4, COLT = s$COL_4 + s$XCL_4,
             time = s$time),
      tibble(bin = "bin3", TGF = s$TGF_3, COLT = s$COL_3 + s$XCL_3,
             time = s$time))
    fp <- tibble(TGF = c(0.0401, 0.2115, 0.7767),
                 COLT = c(1.000, 1.586, 2.940),
                 lab = c("Normal (stable)", "Separatrix (UNSTABLE)",
                         "Fibrotic (stable)"))
    ggplot(d, aes(TGF, COLT, colour = bin)) +
      geom_path(linewidth = .9) +
      geom_point(data = fp, aes(TGF, COLT), inherit.aes = FALSE,
                 size = 4, shape = c(19, 1, 19), colour = "#7D3C98") +
      geom_text(data = fp, aes(TGF, COLT, label = lab), inherit.aes = FALSE,
                hjust = -.12, size = 3.3, colour = "#7D3C98") +
      geom_hline(yintercept = 1.586, linetype = "dashed",
                 colour = "#C0392B") +
      labs(title = "Phase plane — active TGF-β1 vs total collagen",
           subtitle = "Each trajectory is one dose bin. Only bins that cross the separatrix become fibrotic.",
           x = "Active TGF-β1", y = "Total collagen (soluble + cross-linked)",
           colour = NULL) + THEME
  })

  output$fpTab <- renderTable({
    tibble(
      `Fixed point` = c("Normal (healthy)", "Separatrix",
                                 "Fibrotic"),
      `Active TGF-β1` = c(0.0401, 0.2115, 0.7767),
      `Myofibroblasts` = c(0.00057, 0.3563, 0.9692),
      `Total collagen` = c(1.000, 1.586, 2.940),
      `Stability` = c("Stable", "UNSTABLE", "Stable"))
  }, digits = 4, striped = TRUE)

  ## ---- 7  pulmonary function ---------------------------------------
  output$pftPlot <- renderPlot({
    d <- sim() %>% select(time, DLCO, FVC) %>% pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) +
      geom_line(linewidth = 1) +
      geom_vline(xintercept = tcourse(), linetype = "dotted") +
      scale_colour_manual(values = c(DLCO = "#4CA85F", FVC = "#3FA3AD")) +
      labs(title = "Pulmonary function over time (DLCO · FVC)",
           subtitle = paste0("The permanent component comes from the capillary-loss ceiling ",
                             "EC_max = 1/(1+KECIRR·BED) — ",
                             "it exists even in bins that do not become fibrotic."),
           x = "Time (days)", y = "% predicted", colour = NULL) + THEME
  })

  output$binMap <- renderPlot({
    s <- tail(sim(), 1)
    v <- vvec()
    d <- tibble(
      bin = factor(1:6),
      dose = input$DRX * REPF,
      vol = 100 * v,
      AT2 = c(s$AT2_1, s$AT2_2, s$AT2_3, s$AT2_4, s$AT2_5, s$AT2_6),
      EC = c(s$EC_1, s$EC_2, s$EC_3, s$EC_4, s$EC_5, s$EC_6),
      COLT = c(s$COL_1 + s$XCL_1, s$COL_2 + s$XCL_2, s$COL_3 + s$XCL_3,
               s$COL_4 + s$XCL_4, s$COL_5 + s$XCL_5,
               s$COL_6 + s$XCL_6)) %>%
      pivot_longer(c(AT2, EC, COLT))
    lab <- c(AT2 = "Surviving AT2", EC = "Endothelial integrity EC",
             COLT = "Total collagen")
    ggplot(d, aes(dose, value, size = vol)) +
      geom_point(colour = "#6A75C9", alpha = .8) +
      facet_wrap(~ name, scales = "free_y",
                 labeller = labeller(name = lab)) +
      scale_size_continuous(range = c(2, 9), name = "Lung volume (%)") +
      labs(title = paste0("Final state by bin (day ", input$tend, ")"),
           subtitle = "Point size = lung volume occupied by that bin",
           x = "Bin representative dose (Gy)", y = NULL) + THEME
  })

  ## ---- 8  clinical endpoints ---------------------------------------
  output$ntcpPlot <- renderPlot({
    s <- sim()
    ggplot(s, aes(time, 100 * NTCP)) +
      geom_line(linewidth = 1, colour = "#C7538C") +
      geom_vline(xintercept = tcourse(), linetype = "dotted") +
      geom_hline(yintercept = 20, linetype = "dashed",
                 colour = "#888888") +
      annotate("text", x = input$tend * .7, y = 21,
               label = "QUANTEC 20% (MLD 20 Gy)", size = 3,
               colour = "#888888") +
      labs(title = "NTCP — risk of CTCAE grade ≥ 2 radiation pneumonitis",
           x = "Time (days)", y = "Risk (%)") + THEME
  })

  output$gradePlot <- renderPlot({
    s <- sim()
    ggplot(s, aes(time, GRADE)) +
      geom_step(linewidth = 1, colour = "#E06AA0") +
      geom_vline(xintercept = tcourse(), linetype = "dotted") +
      scale_y_continuous(breaks = 0:4, limits = c(0, 4)) +
      labs(title = "CTCAE grade (derived from the continuous state)",
           subtitle = "Grade is an output, not an input",
           x = "Time (days)", y = "grade") + THEME
  })

  output$endTab <- renderTable({
    s <- sim(); tc <- tcourse()
    pk <- s[which.max(s$PNIE), ]
    d365 <- s[which.min(abs(s$time - 365)), ]
    tibble(
      `Endpoint` = c(
        "Pneumonitis index peak PNIe", "Time of peak (days after end of RT)",
        "Time of peak (weeks)", "NTCP grade≥2 (%)", "Peak CTCAE grade",
        "DLCO nadir (%pred)", "DLCO change at 365 days (%)",
        "Fibrotic-conversion volume (%)", "Fibrosis score FIB at 365 days"),
      `Value` = c(
        sprintf("%.4f", max(s$PNIE)),
        sprintf("%.0f", pk$time - tc), sprintf("%.1f", (pk$time - tc) / 7),
        sprintf("%.1f", 100 * max(s$NTCP)),
        sprintf("%.0f", max(s$GRADE)),
        sprintf("%.2f", min(s$DLCO)),
        sprintf("%.2f", 100 * (d365$DLCO - input$DLCO0) / input$DLCO0),
        sprintf("%.1f", 100 * last(s$VFIB)),
        sprintf("%.4f", d365$FIB)))
  }, striped = TRUE)

  ## ---- 9  therapeutic ratio ----------------------------------------
  trData <- reactive({
    bind_rows(lapply(names(PLANS), function(nm) {
      s <- simplan(nm)
      tibble(plan = nm, MLD = first(s$MLD),
             dfx = PLANS[[nm]]$DRX / PLANS[[nm]]$NFX,
             TCP = max(s$TCP), NTCP = 100 * max(s$NTCP),
             UCP = max(s$TCP) / 100 * (1 - max(s$NTCP)) * 100)
    }))
  })

  output$trPlot <- renderPlot({
    d <- trData()
    ggplot(d, aes(NTCP, TCP, label = plan)) +
      geom_point(aes(size = dfx, colour = MLD)) +
      geom_text(hjust = -.08, size = 3.2) +
      scale_colour_gradient(low = "#4CA85F", high = "#C0392B",
                            name = "MLD (Gy)") +
      scale_size_continuous(range = c(3, 10), name = "Fraction dose (Gy)") +
      expand_limits(x = c(0, 45)) +
      labs(title = "Therapeutic ratio — tumour control probability vs pneumonitis risk",
           subtitle = paste0("Top-left is good. The reason SBRT wins on both axes ",
                             "is geometry, not radiobiology ",
                             "(MLD 17.8 → 4.3 Gy)."),
           x = "NTCP grade≥2 (%)", y = "TCP (%)") + THEME
  })

  output$trTab <- renderTable({
    trData() %>% arrange(desc(UCP)) %>%
      rename(`Plan` = plan, `Fraction dose (Gy)` = dfx,
             `MLD (Gy)` = MLD, `TCP (%)` = TCP,
             `NTCP (%)` = NTCP, `UCP (%)` = UCP)
  }, digits = 2, striped = TRUE)

  ## ---- 10  scenario comparison -------------------------------------
  scenData <- reactive({
    req(length(input$scen) > 0)
    bind_rows(lapply(input$scen, simplan))
  })

  output$scenPlot <- renderPlot({
    d <- scenData() %>%
      select(time, plan, DLCO, PNIE, NTCP, VFIB) %>%
      mutate(NTCP = 100 * NTCP, VFIB = 100 * VFIB) %>%
      pivot_longer(c(DLCO, PNIE, NTCP, VFIB))
    lab <- c(DLCO = "DLCO (%pred)", PNIE = "Pneumonitis index PNIe",
             NTCP = "NTCP grade≥2 (%)", VFIB = "Fibrotic-conversion volume (%)")
    ggplot(d, aes(time, value, colour = plan)) +
      geom_line(linewidth = .9) +
      facet_wrap(~ name, scales = "free_y",
                 labeller = labeller(name = lab)) +
      labs(title = "Scenario comparison", x = "Time (days)", y = NULL,
           colour = NULL) + THEME
  })

  output$scenTab <- renderTable({
    scenData() %>% group_by(plan) %>%
      summarise(`MLD (Gy)` = first(MLD), `V20 (%)` = 100 * first(V20G),
                `V40 (%)` = 100 * first(V40G),
                `PNIe peak` = max(PNIE),
                `NTCP (%)` = 100 * max(NTCP),
                `DLCO nadir` = min(DLCO),
                `Vfib (%)` = 100 * last(VFIB), .groups = "drop") %>%
      arrange(`MLD (Gy)`)
  }, digits = 2, striped = TRUE)
}

shinyApp(ui, server)
