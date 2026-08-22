## =============================================================================
##  CHRONIC HEPATITIS D (HDV) -- QSP DASHBOARD  (Shiny)
##  hdv_shiny_app.R      companion to hdv_mrgsolve_model.R
## =============================================================================
##
##  RUN
##      install.packages(c("shiny", "mrgsolve", "ggplot2", "dplyr", "tidyr", "DT"))
##      shiny::runApp("hdv_shiny_app.R")
##
##  WHAT THIS APP IS FOR
##  --------------------
##  Not "here are some curves".  Each tab answers one question that the model
##  exists to answer, and the tab titles say which:
##
##   1  Patient       who is this patient, and where does the floor sit for them
##   2  Drug PK       exposure and, for bulevirtide, DERIVED NTCP occupancy
##   3  Virology      serum HDV RNA + the intracellular pool, side by side --
##                    this is where the lonafarnib paradox becomes visible
##   4  Endpoints     HDV RNA response, ALT normalisation, combined response
##   5  Decoupling    why ALT and HDV RNA disagree, plotted against each other
##   6  The FLOOR     decomposition of the inflow maintaining Id, and the
##                    ceiling of the entire entry-inhibitor class
##   7  Bile acids    total bile acids as an NTCP-occupancy read-out, and the
##                    benefit-saturates-while-cost-does-not dose curve
##   8  Compare       any number of regimens overlaid on any output
##   9  Outcomes      fibrosis, platelets, HCC over 5 years
##
##  The numbers the app should reproduce are in hdv_model_report.txt.
## =============================================================================

library(shiny)
library(mrgsolve)
library(ggplot2)
library(dplyr)
library(tidyr)

MODEL_FILE <- "hdv_mrgsolve_model.R"
mod0 <- mread_cache("hdv", MODEL_FILE)

MW_BLV <- 5398.9

## ---------------------------------------------------------------- dosing ----
ev_blv <- function(mg, dur) if (mg <= 0 || dur <= 0) NULL else
  ev(time = 0, amt = mg * 1e-3 / MW_BLV * 1e9, cmt = "ASC", ii = 1, addl = dur - 1)
ev_ifn <- function(ug, dur) if (ug <= 0 || dur < 7) NULL else
  ev(time = 0, amt = ug, cmt = "ISC", ii = 7, addl = floor(dur / 7) - 1)
ev_lnf <- function(mg, rtv, dur) {
  if (mg <= 0 || dur <= 0) return(NULL)
  e <- ev(time = 0, amt = mg * 1e-3 / 638.6 * 1e6, cmt = "LGUT",
          ii = 0.5, addl = dur * 2 - 1)
  if (rtv > 0) e <- e + ev(time = 0, amt = 0.7 * rtv * 1e-3 / 720.9 * 1e6,
                           cmt = "RCEN", ii = 0.5, addl = dur * 2 - 1)
  e
}
ev_sir <- function(mg, dur) if (mg <= 0 || dur < 28) NULL else
  ev(time = 0, amt = mg * 1e-3 / 16000 * 1e9, cmt = "QCEN",
     ii = 28, addl = floor(dur / 28) - 1)

bind_ev <- function(...) {
  es <- Filter(Negate(is.null), list(...))
  if (!length(es)) return(NULL)
  Reduce(`+`, es)
}

## ------------------------------------------------------------ simulation ----
simulate_regimen <- function(pars, horizon_wk, blv, ifn, lam, lnf, rtv, sir,
                             tx_wk, occfix = NA) {
  dur <- tx_wk * 7
  m <- param(mod0, pars)
  if (!is.na(occfix)) m <- param(m, list(OCCFIX = occfix))
  m <- param(m, list(LAMBDA = if (lam) 1 else 0))
  e <- bind_ev(ev_blv(blv, dur), ev_ifn(ifn, dur),
               ev_lnf(lnf, rtv, dur), ev_sir(sir, dur))
  m <- update(m, end = horizon_wk * 7, delta = 1)
  d <- as.data.frame(if (is.null(e)) mrgsim(m) else mrgsim(m, events = e))
  d$week <- d$time / 7
  d
}

## per-patient parameter overrides from the sidebar
patient_pars <- function(input) {
  list(
    DIMM  = input$dimm,
    DHEP  = input$dhep,
    KCC   = 0.0137499 * input$ccmult,
    PHINT = input$phint,
    KCURE = input$kcure,
    ALTBASE = input$altbase,
    NUC   = if (input$nuc) 1 else 0
  )
}

## initial conditions from the sidebar (baseline severity)
patient_init <- function(m, input) {
  init(m, list(VD = 10^input$lgvd0, SSER = 10^input$lgsag0,
               ID = input$id0, ALT = input$alt0, FIB = input$fib0))
}

simulate_full <- function(input, blv, ifn, lam, lnf, rtv, sir, tx_wk,
                          horizon_wk, occfix = NA) {
  m <- param(mod0, patient_pars(input))
  if (!is.na(occfix)) m <- param(m, list(OCCFIX = occfix))
  m <- param(m, list(LAMBDA = if (lam) 1 else 0))
  m <- patient_init(m, input)
  dur <- tx_wk * 7
  e <- bind_ev(ev_blv(blv, dur), ev_ifn(ifn, dur),
               ev_lnf(lnf, rtv, dur), ev_sir(sir, dur))
  m <- update(m, end = horizon_wk * 7, delta = 1)
  d <- as.data.frame(if (is.null(e)) mrgsim(m) else mrgsim(m, events = e))
  d$week <- d$time / 7
  d$dlog <- d$LGVD - input$lgvd0
  d
}

THEME <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 10, colour = "grey30"),
        legend.position = "bottom")

## ==================================================================== UI ====
ui <- fluidPage(
  titlePanel("Chronic Hepatitis D — QSP dashboard"),
  tags$p(style = "color:#555;margin-top:-8px",
         HTML("A model that casts HDV as a <b>two-substrate assembly line</b>: the genome is supplied by host RNA Pol II, ",
              "the envelope by HBV cccDNA. Every drug is classified by <b>which flux it cuts</b> — ",
              "entry (bulevirtide) · envelopment (lonafarnib · siRNA) · intracellular replication and loss of infected cells (interferon).")),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("① Patient profile"),
      sliderInput("lgvd0", "Baseline HDV RNA (log10 IU/mL)", 3.0, 8.0, 5.5, 0.1),
      sliderInput("lgsag0", "Baseline HBsAg (log10 IU/mL)", 2.0, 5.0, 3.85, 0.05),
      sliderInput("id0", "Fraction of HDAg-positive hepatocytes", 0.01, 0.30, 0.10, 0.01),
      sliderInput("alt0", "Baseline ALT (U/L)", 40, 400, 110, 5),
      sliderInput("altbase", "Non-HDV ALT component (U/L)", 15, 80, 42, 1),
      sliderInput("fib0", "Baseline fibrosis (Ishak 0–6)", 0, 6, 2, 0.5),
      checkboxInput("nuc", "Concomitant nucleoside analogue (NUC)", TRUE),

      hr(), h4("② Host cell biology — the floor"),
      helpText(HTML("<small>These four set the limit of an entry inhibitor. ",
                    "Target affinity does not.</small>")),
      sliderInput("dimm", "Infected-cell killing rate DIMM (1/day)", 0.005, 0.08, 0.02636, 0.001),
      sliderInput("dhep", "Hepatocyte turnover rate DHEP (1/day)", 0.001, 0.012, 0.004, 0.0005),
      sliderInput("ccmult", "Cell-to-cell spread multiplier (× KCC)", 0.0, 3.0, 1.0, 0.1),
      sliderInput("phint", "NTCP-dependent fraction of cell-to-cell spread", 0.0, 1.0, 0.5, 0.05),
      sliderInput("kcure", "Non-cytolytic intracellular clearance (1/day)", 0.0, 0.01, 0.002, 0.0005),

      hr(), h4("③ Therapy"),
      sliderInput("blv", "Bulevirtide (mg SC qd)", 0, 20, 2, 0.5),
      sliderInput("ifn", "Peg-IFN (µg SC qw)", 0, 240, 0, 10),
      checkboxInput("lam", "Switch to Peg-IFN lambda (no haematological toxicity, hepatotoxicity present)", FALSE),
      sliderInput("lnf", "Lonafarnib (mg PO BID)", 0, 100, 0, 5),
      sliderInput("rtv", "Ritonavir booster (mg PO BID)", 0, 200, 100, 50),
      sliderInput("sir", "HBsAg siRNA (mg SC q4w)", 0, 400, 0, 50),
      sliderInput("txwk", "Treatment duration (weeks)", 4, 144, 48, 4),
      sliderInput("horiz", "Observation horizon (weeks)", 12, 260, 96, 4),
      checkboxInput("perfect", "Also show a hypothetical 'complete' entry inhibitor (100% occupancy)", TRUE)
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        type = "tabs",

        tabPanel("1 · Patient · the floor", br(),
          fluidRow(column(6, h5("Decomposition of the flux maintaining Id in this patient"),
                          plotOutput("floorbar", height = "300px")),
                   column(6, h5("Hepatocyte pools"), plotOutput("pools", height = "300px"))),
          hr(), verbatimTextOutput("patient_txt")),

        tabPanel("2 · Drug PK and target occupancy", br(),
          fluidRow(column(6, plotOutput("pk_blv", height = "280px")),
                   column(6, plotOutput("pk_occ", height = "280px"))),
          fluidRow(column(6, plotOutput("pk_ifn", height = "260px")),
                   column(6, plotOutput("pk_lnf", height = "260px"))),
          hr(), verbatimTextOutput("occ_txt")),

        tabPanel("3 · Virology (serum vs intracellular)", br(),
          plotOutput("vir_serum", height = "300px"),
          plotOutput("vir_intra", height = "280px"),
          helpText(HTML("<b>The lonafarnib paradox.</b> Block assembly and the genome is trapped inside the cell. ",
                        "Serum HDV RNA measures the <i>assembly flux</i>, not the size of the reservoir — ",
                        "which is why rebound on withdrawal is fast."))),

        tabPanel("4 · Clinical endpoints", br(),
          plotOutput("ep_plot", height = "320px"),
          hr(), h5("Summary at key time points"), tableOutput("ep_tab"),
          helpText(HTML("Response is defined as in MYR301: HDV RNA below the limit of detection or a fall of 2 log10 or more, ",
                        "together with ALT normalisation. Combined response requires both."))),

        tabPanel("5 · ALT–RNA decoupling", br(),
          fluidRow(column(7, plotOutput("dec_time", height = "320px")),
                   column(5, plotOutput("dec_scatter", height = "320px"))),
          hr(), verbatimTextOutput("dec_txt")),

        tabPanel("6 · The floor and the ceiling", br(),
          plotOutput("ceiling", height = "320px"),
          hr(), h5("Sensitivity to the floor parameters (log10 fall at week 48)"),
          tableOutput("sens_tab"),
          helpText(HTML("<b>The falsifiable central claim of this model:</b> the week-48 response is ",
                        "more sensitive to host cell biology (DHEP, KCC, PHINT) than to ",
                        "target affinity (KDNTCP)."))),

        tabPanel("7 · Bile acids = an occupancy meter", br(),
          fluidRow(column(6, plotOutput("ba_time", height = "300px")),
                   column(6, plotOutput("ba_dose", height = "300px"))),
          hr(), tableOutput("ba_tab"),
          helpText(HTML("NTCP is both the HDV receptor and a bile acid transporter. So a rise in total bile acids ",
                        "becomes a <b>direct read-out of target occupancy</b>. The benefit saturates but ",
                        "the cost does not."))),

        tabPanel("8 · Scenario comparison", br(),
          checkboxGroupInput("cmp", "Regimens to compare",
            choices = c("No treatment (NUC only)" = "none",
                        "BLV 2 mg" = "blv2",
                        "BLV 10 mg" = "blv10",
                        "Peg-IFN alfa 180 µg" = "ifn",
                        "Peg-IFN lambda 180 µg" = "lam",
                        "BLV 2 mg + Peg-IFN" = "blvifn",
                        "LNF 50 + RTV 100 BID" = "lnf",
                        "LNF/RTV + Peg-IFN" = "lnfifn",
                        "BLV 2 mg + HBsAg siRNA" = "blvsir",
                        "Total entry blockade (hypothetical)" = "perfect"),
            selected = c("none", "blv2", "blv10", "ifn", "blvifn", "lnf"),
            inline = TRUE),
          selectInput("cmp_y", "Variable to display",
            choices = c("Δlog10 HDV RNA" = "dlog", "ALT (U/L)" = "ALT",
                        "HDAg-positive hepatocytes (%)" = "HDAGPC",
                        "Intracellular HDV RNA (fold)" = "RGFOLD",
                        "Serum HBsAg (log10)" = "LGSAG",
                        "Total bile acids (fold)" = "TBAFLD",
                        "Exhaustion level" = "EXH",
                        "Fibrosis (Ishak)" = "FIB",
                        "Platelets (10^9/L)" = "PLT",
                        "Neutrophils (10^9/L)" = "NEU"),
            selected = "dlog"),
          plotOutput("cmp_plot", height = "420px"),
          hr(), tableOutput("cmp_tab")),

        tabPanel("9 · Long-term outcomes (5 years)", br(),
          fluidRow(column(6, plotOutput("out_fib", height = "290px")),
                   column(6, plotOutput("out_hcc", height = "290px"))),
          fluidRow(column(6, plotOutput("out_plt", height = "260px")),
                   column(6, plotOutput("out_alt", height = "260px"))),
          helpText(HTML("Fibrosis follows <b>ALT</b>, not HDV RNA. ",
                        "So the long-term ranking of regimens is decided not by virology but by ",
                        "normalisation of inflammation."))),

        tabPanel("ℹ Model information", br(),
          verbatimTextOutput("about"))
      )
    )
  )
)

## ================================================================ SERVER ====
server <- function(input, output, session) {

  cur <- reactive({
    simulate_full(input, input$blv, input$ifn, input$lam, input$lnf,
                  input$rtv, input$sir, input$txwk, input$horiz)
  })
  perf <- reactive({
    if (!isTRUE(input$perfect)) return(NULL)
    simulate_full(input, 0, 0, FALSE, 0, 0, 0, input$txwk, input$horiz,
                  occfix = 1.0)
  })
  untx <- reactive({
    simulate_full(input, 0, 0, FALSE, 0, 0, 0, input$txwk, input$horiz)
  })

  ## ---- 1. patient & floor ----
  output$floorbar <- renderPlot({
    d <- cur()
    p <- unlist(param(mod0))
    ID <- input$id0; IB <- 0.35
    occ <- d$OCCUP[which.min(abs(d$week - min(input$txwk, input$horiz)))]
    fentry <- p[["BETAD"]] * 10^input$lgvd0 * IB
    fcc    <- 0.0137499 * input$ccmult * ID * IB
    dthId  <- input$dhep + p[["DINN"]] + input$dimm
    Hd     <- (input$dhep + p[["DHBV"]]) * IB + dthId * ID
    fdiv   <- p[["LOCREN"]] * Hd * ID / (IB + ID)
    df <- data.frame(
      flux = factor(c("(E) entry/reinfection", "(C) cell-to-cell spread", "(D) division-mediated spread"),
                    levels = c("(E) entry/reinfection", "(C) cell-to-cell spread", "(D) division-mediated spread")),
      untreated = c(fentry, fcc, fdiv),
      treated = c(fentry * (1 - occ), fcc * (1 - input$phint * occ), fdiv)) %>%
      pivot_longer(-flux, names_to = "state", values_to = "v")
    df$state <- factor(df$state, levels = c("untreated", "treated"),
                       labels = c("untreated", sprintf("on treatment (occupancy %.2f)", occ)))
    ggplot(df, aes(flux, v, fill = state)) +
      geom_col(position = "dodge") +
      scale_fill_manual(values = c("#b04a4a", "#2f7a2f"), name = NULL) +
      labs(x = NULL, y = "Inflow flux (fraction/day)",
           title = "The three fluxes that maintain the infected-cell pool",
           subtitle = "(D) is NTCP-independent — the floor an entry inhibitor can never reach") +
      THEME
  })

  output$pools <- renderPlot({
    d <- cur() %>% select(week, T, IB, ID) %>%
      pivot_longer(-week, names_to = "pool", values_to = "v")
    d$pool <- factor(d$pool, levels = c("T", "IB", "ID"),
                     labels = c("T · HBsAg negative", "Ib · HBsAg+ HDV−(envelope supply)",
                                "Id · HBsAg+ HDV+"))
    ggplot(d, aes(week, v, colour = pool)) + geom_line(linewidth = .8) +
      scale_colour_manual(values = c("#556699", "#c07a30", "#b04a4a"), name = NULL) +
      labs(x = "Weeks", y = "Hepatocyte fraction", title = "Hepatocyte compartments") + THEME
  })

  output$patient_txt <- renderPrint({
    d <- cur(); u <- untx()
    i48 <- which.min(abs(d$week - 48))
    cat(sprintf(
"Baseline: HDV RNA %.2f log10 IU/mL · HBsAg %.2f log10 · HDAg+ hepatocytes %.0f%% · ALT %.0f U/L · Ishak %.1f
Week 48: HDV RNA change %+.2f log10 · ALT %.0f U/L (%s) · HDAg+ %.2f%% · intracellular RNA %.2f-fold
     NTCP occupancy %.3f · free NTCP %.1f%% · total bile acids %.1f-fold
Untreated control at week 48: HDV RNA change %+.2f log10 · ALT %.0f U/L
",
      input$lgvd0, input$lgsag0, 100 * input$id0, input$alt0, input$fib0,
      d$dlog[i48], d$ALT[i48],
      ifelse(d$ALT[i48] <= 40, "normalised", "abnormal"),
      d$HDAGPC[i48], d$RGFOLD[i48], d$OCCUP[i48], 100 * d$FREENT[i48],
      d$TBAFLD[i48], u$dlog[i48], u$ALT[i48]))
  })

  ## ---- 2. PK ----
  output$pk_blv <- renderPlot({
    d <- cur() %>% filter(week <= min(8, input$horiz))
    ggplot(d, aes(week * 7, CBLV)) + geom_line(colour = "#2f7a2f", linewidth = .8) +
      labs(x = "Days", y = "Bulevirtide (nM)", title = "Bulevirtide plasma concentration",
           subtitle = "Saturable target-mediated disposition → more-than-proportional exposure with dose") + THEME
  })
  output$pk_occ <- renderPlot({
    d <- cur()
    ggplot(d, aes(week, OCCUP)) + geom_line(colour = "#2f7a2f", linewidth = .9) +
      geom_hline(yintercept = 1, linetype = 3) +
      scale_y_continuous(limits = c(0, 1)) +
      labs(x = "Weeks", y = "NTCP occupancy", title = "NTCP occupancy (derived, not assumed)",
           subtitle = "Calculated through a Kd back-fitted from the two bile-acid fold-rise values") + THEME
  })
  output$pk_ifn <- renderPlot({
    d <- cur()
    ggplot(d, aes(week)) +
      geom_line(aes(y = CIFN, colour = "Peg-IFN (ng/mL)"), linewidth = .8) +
      geom_line(aes(y = 20 * ISG, colour = "ISG × 20"), linewidth = .8) +
      geom_line(aes(y = 20 * SOCS, colour = "SOCS1 × 20"), linewidth = .6, linetype = 2) +
      scale_colour_manual(values = c("#2f6a9a", "#3a8a8a", "#7d4a9a"), name = NULL) +
      labs(x = "Weeks", y = NULL, title = "Interferon exposure and ISG induction",
           subtitle = "SOCS1 negative feedback = tachyphylaxis") + THEME
  })
  output$pk_lnf <- renderPlot({
    d <- cur() %>% filter(week <= min(4, input$horiz))
    ggplot(d, aes(week * 7)) +
      geom_line(aes(y = CLNF, colour = "Lonafarnib (µM)"), linewidth = .8) +
      geom_line(aes(y = CRTV, colour = "Ritonavir (µM)"), linewidth = .8) +
      geom_line(aes(y = IFTC, colour = "FTase inhibition"), linewidth = .9, linetype = 2) +
      scale_colour_manual(values = c("#7d4a9a", "#a04a70", "#333333"), name = NULL) +
      labs(x = "Days", y = NULL, title = "Lonafarnib / ritonavir",
           subtitle = "CYP3A4-inhibition boosting raises exposure 4–5 fold") + THEME
  })
  output$occ_txt <- renderPrint({
    p <- unlist(param(mod0))
    css <- function(mg) {
      inp <- p[["FBLV"]] * (mg * 1e-3 / MW_BLV * 1e9)
      f <- function(c) p[["CLBLV"]] * c + p[["VMBLV"]] * c / (p[["KMBLV"]] + c) - inp
      uniroot(f, c(0, 1e5))$root
    }
    o <- function(mg) { c <- css(mg); c / (c + p[["KDNTCP"]]) }
    cat(sprintf(
"Fitted Kd_NTCP = %.3f nM  (back-fitted from the two bile-acid fold-rise values alone)
  2 mg : Css %.2f nM → occupancy %.4f → free NTCP %.1f%%  → bile acids %.1f-fold
 10 mg : Css %.2f nM → occupancy %.4f → free NTCP %.1f%%  → bile acids %.1f-fold
Ratio of residual entry flux (2 mg / 10 mg) = %.1f-fold

The key point: the approved 2 mg dose does not saturate NTCP. Even so, 10 mg barely raises the
combined response — because entry is only part of the inflow that maintains Id.
",
      p[["KDNTCP"]], css(2), o(2), 100 * (1 - o(2)),
      (1 + p[["ROATP"]]) / ((1 - o(2)) + p[["ROATP"]]),
      css(10), o(10), 100 * (1 - o(10)),
      (1 + p[["ROATP"]]) / ((1 - o(10)) + p[["ROATP"]]),
      (1 - o(2)) / (1 - o(10))))
  })

  ## ---- 3. virology ----
  output$vir_serum <- renderPlot({
    d <- cur(); u <- untx(); pf <- perf()
    df <- rbind(data.frame(week = d$week, v = d$dlog, arm = "Selected regimen"),
                data.frame(week = u$week, v = u$dlog, arm = "No treatment"))
    if (!is.null(pf)) df <- rbind(df, data.frame(week = pf$week, v = pf$dlog,
                                                 arm = "Total entry blockade (hypothetical)"))
    ggplot(df, aes(week, v, colour = arm)) + geom_line(linewidth = .9) +
      geom_hline(yintercept = -2, linetype = 2, colour = "grey40") +
      annotate("text", x = 2, y = -2.15, label = "2 log10 response threshold",
               hjust = 0, size = 3, colour = "grey35") +
      geom_vline(xintercept = input$txwk, linetype = 3) +
      scale_colour_manual(values = c("#b04a4a", "#2f7a2f", "#888888"), name = NULL) +
      labs(x = "Weeks", y = expression(Delta*log[10]*" HDV RNA"),
           title = "Serum HDV RNA",
           subtitle = "Dotted vertical line = end of treatment. An entry inhibitor has no phase 1 (first-week plunge)") + THEME
  })
  output$vir_intra <- renderPlot({
    d <- cur()
    ggplot(d, aes(week)) +
      geom_line(aes(y = RGFOLD, colour = "Intracellular HDV RNA (fold)"), linewidth = .9) +
      geom_line(aes(y = 10^(d$dlog), colour = "Serum HDV RNA (fold)"), linewidth = .9) +
      geom_hline(yintercept = 1, linetype = 3) +
      scale_y_log10() +
      scale_colour_manual(values = c("#c07a30", "#b04a4a"), name = NULL) +
      labs(x = "Weeks", y = "Fold vs baseline (log axis)",
           title = "Intracellular vs serum — they do not move in the same direction",
           subtitle = "Switch lonafarnib on and the two curves diverge (blocking assembly traps the genome)") + THEME
  })

  ## ---- 4. endpoints ----
  output$ep_plot <- renderPlot({
    d <- cur()
    df <- data.frame(week = rep(d$week, 3),
                     v = c(d$VR, d$ALTN, d$COMB),
                     ep = rep(c("HDV RNA response", "ALT normalisation", "Combined response"),
                              each = nrow(d)))
    ggplot(df, aes(week, v, colour = ep)) +
      geom_step(linewidth = .9) +
      geom_vline(xintercept = input$txwk, linetype = 3) +
      scale_y_continuous(breaks = c(0, 1), labels = c("not met", "met")) +
      scale_colour_manual(values = c("#2f6a9a", "#b04a4a", "#111111"), name = NULL) +
      labs(x = "Weeks", y = NULL, title = "Endpoints met (this patient)",
           subtitle = "Population response rates are computed in the virtual population of hdv_reference_model.py") + THEME
  })
  output$ep_tab <- renderTable({
    d <- cur()
    wks <- c(4, 12, 24, 48, 96, input$horiz)
    wks <- unique(wks[wks <= input$horiz])
    do.call(rbind, lapply(wks, function(wk) {
      i <- which.min(abs(d$week - wk))
      data.frame(Week = wk,
                 `Δlog10 HDV` = sprintf("%+.2f", d$dlog[i]),
                 `ALT (U/L)` = sprintf("%.0f", d$ALT[i]),
                 `HDAg+ (%)` = sprintf("%.2f", d$HDAGPC[i]),
                 `Intracellular RNA (fold)` = sprintf("%.2f", d$RGFOLD[i]),
                 `HBsAg (log10)` = sprintf("%.2f", d$LGSAG[i]),
                 `Bile acids (fold)` = sprintf("%.1f", d$TBAFLD[i]),
                 `RNA response` = ifelse(d$VR[i] > .5, "O", "-"),
                 `ALT normal` = ifelse(d$ALTN[i] > .5, "O", "-"),
                 `Combined` = ifelse(d$COMB[i] > .5, "O", "-"),
                 check.names = FALSE)
    }))
  })

  ## ---- 5. decoupling ----
  output$dec_time <- renderPlot({
    d <- cur()
    ggplot(d, aes(week)) +
      geom_line(aes(y = 100 * ALT / input$alt0, colour = "ALT (% of baseline)"),
                linewidth = .9) +
      geom_line(aes(y = 100 * 10^dlog, colour = "HDV RNA (% of baseline)"),
                linewidth = .9) +
      scale_colour_manual(values = c("#b04a4a", "#2f6a9a"), name = NULL) +
      labs(x = "Weeks", y = "% of baseline",
           title = "ALT falls first",
           subtitle = "ALT reads the 'inflow' of newly infected cells, HDV RNA reads the 'size' of the infected pool") +
      THEME
  })
  output$dec_scatter <- renderPlot({
    d <- cur()
    ggplot(d, aes(dlog, ALT, colour = week)) + geom_path(linewidth = .9) +
      geom_hline(yintercept = 40, linetype = 2) +
      geom_vline(xintercept = -2, linetype = 2) +
      annotate("text", x = -0.2, y = 38, label = "ALT normalisation only",
               hjust = 1, size = 3, colour = "grey30") +
      scale_colour_viridis_c(name = "Weeks") +
      labs(x = expression(Delta*log[10]*" HDV RNA"), y = "ALT (U/L)",
           title = "The trajectory of the two endpoints",
           subtitle = "Lower-left quadrant = combined response") + THEME
  })
  output$dec_txt <- renderPrint({
    d <- cur(); i4 <- which.min(abs(d$week - 4)); i48 <- which.min(abs(d$week - 48))
    cat(sprintf(
"Week 4:  ALT has already fallen %.0f%%, HDV RNA only %.0f%% (%+.2f log10)
Week 48: ALT down %.0f%%, HDV RNA down %.0f%% (%+.2f log10)

This decoupling is the model's 'hypothesis': the damage term contains kappa × (entry flux).
Falsification condition — even in virological non-responders ALT should fall in proportion to (1 − occupancy).
If the fall in ALT tracks the fall in HDV RNA patient by patient, this term is wrong.
",
      100 * (1 - d$ALT[i4] / input$alt0), 100 * (1 - 10^d$dlog[i4]), d$dlog[i4],
      100 * (1 - d$ALT[i48] / input$alt0), 100 * (1 - 10^d$dlog[i48]), d$dlog[i48]))
  })

  ## ---- 6. ceiling & sensitivity ----
  output$ceiling <- renderPlot({
    doses <- c(0, 0.5, 1, 2, 5, 10, 20)
    df <- do.call(rbind, lapply(doses, function(mg) {
      d <- simulate_full(input, mg, 0, FALSE, 0, 0, 0, 48, 48)
      data.frame(dose = mg, occ = tail(d$OCCUP, 1), dlog = tail(d$dlog, 1))
    }))
    pf <- simulate_full(input, 0, 0, FALSE, 0, 0, 0, 48, 48, occfix = 1)
    ggplot(df, aes(occ, dlog)) +
      geom_line(colour = "#2f7a2f", linewidth = .9) +
      geom_point(size = 2, colour = "#2f7a2f") +
      geom_hline(yintercept = tail(pf$dlog, 1), linetype = 2, colour = "#b04a4a") +
      annotate("text", x = 0.05, y = tail(pf$dlog, 1) + 0.12,
               label = "total blockade = the ceiling of this class", hjust = 0,
               size = 3.2, colour = "#b04a4a") +
      geom_text(aes(label = paste0(dose, " mg")), vjust = -0.8, size = 3) +
      labs(x = "NTCP occupancy", y = expression("week 48 "*Delta*log[10]*" HDV RNA"),
           title = "How far does raising occupancy get you",
           subtitle = "The ceiling is set by host cell biology, not by the drug") + THEME
  })
  output$sens_tab <- renderTable({
    base <- tail(simulate_full(input, 2, 0, FALSE, 0, 0, 0, 48, 48)$dlog, 1)
    pars <- c(DIMM = "dimm", DHEP = "dhep", KCC = "ccmult",
              PHINT = "phint", KCURE = "kcure")
    rows <- lapply(names(pars), function(pn) {
      key <- pars[[pn]]
      v0 <- input[[key]]
      got <- sapply(c(0.7, 1.3), function(k) {
        inp <- reactiveValuesToList(input)
        inp[[key]] <- v0 * k
        m <- param(mod0, patient_pars(inp))
        m <- patient_init(m, inp)
        m <- update(m, end = 336, delta = 14)
        d <- as.data.frame(mrgsim(m, events = ev_blv(2, 336)))
        tail(d$LGVD, 1) - input$lgvd0
      })
      data.frame(Parameter = pn, `-30%` = sprintf("%+.2f", got[1]),
                 Reference = sprintf("%+.2f", base),
                 `+30%` = sprintf("%+.2f", got[2]),
                 `Mean |effect|` = sprintf("%.3f", mean(abs(got - base))),
                 check.names = FALSE)
    })
    kd <- unlist(param(mod0))[["KDNTCP"]]
    gotkd <- sapply(c(0.7, 1.3), function(k) {
      m <- param(mod0, patient_pars(input)); m <- param(m, list(KDNTCP = kd * k))
      m <- patient_init(m, input); m <- update(m, end = 336, delta = 14)
      tail(as.data.frame(mrgsim(m, events = ev_blv(2, 336)))$LGVD, 1) - input$lgvd0
    })
    rows[[length(rows) + 1]] <- data.frame(
      Parameter = "KDNTCP (target affinity)", `-30%` = sprintf("%+.2f", gotkd[1]),
      Reference = sprintf("%+.2f", base), `+30%` = sprintf("%+.2f", gotkd[2]),
      `Mean |effect|` = sprintf("%.3f", mean(abs(gotkd - base))), check.names = FALSE)
    do.call(rbind, rows)
  })

  ## ---- 7. bile acids ----
  output$ba_time <- renderPlot({
    d <- cur()
    ggplot(d, aes(week, TBAFLD)) + geom_line(colour = "#a08a20", linewidth = .9) +
      geom_hline(yintercept = 1, linetype = 3) +
      geom_vline(xintercept = input$txwk, linetype = 3) +
      labs(x = "Weeks", y = "Total bile acids (fold vs baseline)",
           title = "Total bile acids — asymptomatic but informative",
           subtitle = "Rises without pruritus or cholestasis. A direct read-out of target occupancy") + THEME
  })
  output$ba_dose <- renderPlot({
    doses <- c(0.5, 1, 2, 5, 10, 20)
    df <- do.call(rbind, lapply(doses, function(mg) {
      d <- simulate_full(input, mg, 0, FALSE, 0, 0, 0, 48, 48)
      data.frame(dose = mg, ba = tail(d$TBAFLD, 1),
                 eff = -tail(d$dlog, 1))
    }))
    ggplot(df, aes(dose)) +
      geom_line(aes(y = eff, colour = "Efficacy (−Δlog10 HDV, week 48)"), linewidth = .9) +
      geom_point(aes(y = eff, colour = "Efficacy (−Δlog10 HDV, week 48)"), size = 2) +
      geom_line(aes(y = ba / 5, colour = "Cost (bile-acid fold ÷ 5)"), linewidth = .9) +
      geom_point(aes(y = ba / 5, colour = "Cost (bile-acid fold ÷ 5)"), size = 2) +
      geom_vline(xintercept = 2, linetype = 2, colour = "grey40") +
      annotate("text", x = 2.3, y = 0.3, label = "approved dose", hjust = 0, size = 3) +
      scale_x_log10(breaks = doses) +
      scale_colour_manual(values = c("#2f7a2f", "#a08a20"), name = NULL) +
      labs(x = "Bulevirtide dose (mg/day, log axis)", y = NULL,
           title = "The benefit saturates and the cost does not",
           subtitle = "2 mg sits at the knee — without saturating the target") + THEME
  })
  output$ba_tab <- renderTable({
    doses <- c(0.5, 1, 2, 5, 10, 20)
    do.call(rbind, lapply(doses, function(mg) {
      d <- simulate_full(input, mg, 0, FALSE, 0, 0, 0, 48, 48)
      i <- nrow(d)
      data.frame(`Dose (mg)` = mg,
                 `Css (nM)` = sprintf("%.2f", mean(tail(d$CBLV, 48))),
                 `NTCP occupancy` = sprintf("%.4f", d$OCCUP[i]),
                 `Free NTCP (%)` = sprintf("%.1f", 100 * d$FREENT[i]),
                 `Bile acids (fold)` = sprintf("%.1f", d$TBAFLD[i]),
                 `Week 48 Δlog10` = sprintf("%+.2f", d$dlog[i]),
                 `Week 48 ALT` = sprintf("%.0f", d$ALT[i]),
                 check.names = FALSE)
    }))
  })

  ## ---- 8. comparison ----
  cmp_data <- reactive({
    defs <- list(
      none    = list(0, 0, FALSE, 0, 0, 0, NA),
      blv2    = list(2, 0, FALSE, 0, 0, 0, NA),
      blv10   = list(10, 0, FALSE, 0, 0, 0, NA),
      ifn     = list(0, 180, FALSE, 0, 0, 0, NA),
      lam     = list(0, 180, TRUE, 0, 0, 0, NA),
      blvifn  = list(2, 180, FALSE, 0, 0, 0, NA),
      lnf     = list(0, 0, FALSE, 50, 100, 0, NA),
      lnfifn  = list(0, 180, FALSE, 50, 100, 0, NA),
      blvsir  = list(2, 0, FALSE, 0, 0, 200, NA),
      perfect = list(0, 0, FALSE, 0, 0, 0, 1.0)
    )
    labs <- c(none = "No treatment", blv2 = "BLV 2 mg", blv10 = "BLV 10 mg",
              ifn = "Peg-IFN alfa", lam = "Peg-IFN lambda",
              blvifn = "BLV 2 + Peg-IFN", lnf = "LNF/RTV",
              lnfifn = "LNF/RTV + Peg-IFN", blvsir = "BLV 2 + siRNA",
              perfect = "Total entry blockade")
    sel <- input$cmp
    if (!length(sel)) return(NULL)
    do.call(rbind, lapply(sel, function(k) {
      a <- defs[[k]]
      d <- simulate_full(input, a[[1]], a[[2]], a[[3]], a[[4]], a[[5]], a[[6]],
                         input$txwk, input$horiz, occfix = a[[7]])
      d$arm <- labs[[k]]
      d
    }))
  })

  output$cmp_plot <- renderPlot({
    d <- cmp_data(); if (is.null(d)) return(NULL)
    d$y <- d[[input$cmp_y]]
    ggplot(d, aes(week, y, colour = arm)) + geom_line(linewidth = .85) +
      geom_vline(xintercept = input$txwk, linetype = 3) +
      labs(x = "Weeks", y = names(which(c(
        "Δlog10 HDV RNA" = "dlog", "ALT (U/L)" = "ALT",
        "HDAg-positive hepatocytes (%)" = "HDAGPC", "Intracellular HDV RNA (fold)" = "RGFOLD",
        "Serum HBsAg (log10)" = "LGSAG", "Total bile acids (fold)" = "TBAFLD",
        "Exhaustion level" = "EXH", "Fibrosis (Ishak)" = "FIB",
        "Platelets (10^9/L)" = "PLT", "Neutrophils (10^9/L)" = "NEU") == input$cmp_y)),
        title = "Regimen comparison", subtitle = "Dotted line = end of treatment") +
      THEME + theme(legend.position = "right")
  })
  output$cmp_tab <- renderTable({
    d <- cmp_data(); if (is.null(d)) return(NULL)
    d %>% group_by(arm) %>%
      summarise(
        `Δlog10 at end of treatment` = sprintf("%+.2f", dlog[which.min(abs(week - input$txwk))]),
        `ALT at end of treatment` = sprintf("%.0f", ALT[which.min(abs(week - input$txwk))]),
        `Δlog10 at end of follow-up` = sprintf("%+.2f", dlog[n()]),
        `ALT at end of follow-up` = sprintf("%.0f", ALT[n()]),
        `Ishak at end of follow-up` = sprintf("%.2f", FIB[n()]),
        `Bile acids (fold)` = sprintf("%.1f", TBAFLD[which.min(abs(week - input$txwk))]),
        `Neutrophils` = sprintf("%.2f", NEU[which.min(abs(week - input$txwk))]),
        .groups = "drop") %>% as.data.frame()
  })

  ## ---- 9. long-term outcomes ----
  long5 <- reactive({
    defs <- list(
      list("No treatment", 0, 0, FALSE, 0, 0, 0),
      list("BLV 2 mg continuous", 2, 0, FALSE, 0, 0, 0),
      list("BLV 10 mg continuous", 10, 0, FALSE, 0, 0, 0),
      list("BLV 2 + Peg-IFN 48 weeks", 2, 180, FALSE, 0, 0, 0)
    )
    do.call(rbind, lapply(defs, function(a) {
      d <- simulate_full(input, a[[2]], a[[3]], a[[4]], a[[5]], a[[6]], a[[7]],
                         48, 260)
      d$arm <- a[[1]]; d
    }))
  })
  output$out_fib <- renderPlot({
    ggplot(long5(), aes(week / 52, FIB, colour = arm)) + geom_line(linewidth = .85) +
      geom_hline(yintercept = 5, linetype = 2, colour = "grey40") +
      annotate("text", x = 0.1, y = 5.15, label = "cirrhosis", hjust = 0, size = 3) +
      labs(x = "Years", y = "Ishak fibrosis stage", title = "Fibrosis progression/regression") + THEME
  })
  output$out_hcc <- renderPlot({
    ggplot(long5(), aes(week / 52, HCCINC, colour = arm)) +
      geom_line(linewidth = .85) +
      labs(x = "Years", y = "Cumulative HCC incidence (%)", title = "Hepatocellular carcinoma risk") + THEME
  })
  output$out_plt <- renderPlot({
    ggplot(long5(), aes(week / 52, PLT, colour = arm)) + geom_line(linewidth = .85) +
      labs(x = "Years", y = "Platelets (10^9/L)",
           title = "Platelets = a surrogate for portal pressure") + THEME
  })
  output$out_alt <- renderPlot({
    ggplot(long5(), aes(week / 52, ALT, colour = arm)) + geom_line(linewidth = .85) +
      geom_hline(yintercept = 40, linetype = 2, colour = "grey40") +
      labs(x = "Years", y = "ALT (U/L)",
           title = "What moves fibrosis is ALT") + THEME
  })

  output$about <- renderPrint({
    cat("
Chronic Hepatitis D (HDV) QSP model — 33 ODE compartments
==============================================================================
Central idea
  HDV is not 'a virus to be killed' but a two-substrate assembly line.
    Substrate 1  HDV genome     — made by host RNA Pol II (independent of HBV replication)
    Substrate 2  HBsAg envelope — made by HBV cccDNA (the same cell, or the surrounding HBsAg+ pool)
  So drugs are classified by 'which flux they cut'.
    (E) entry/reinfection          NTCP-dependent    bulevirtide
    (A) envelopment/release        FTase, HBsAg      lonafarnib, siRNA/NAP
    (C) intracellular + cell loss  Pol II, immunity  Peg-IFN alfa/lambda
  And the infected-cell pool is maintained by two further routes that have nothing to do with NTCP.
    division-mediated spread  — when an infected cell divides, both daughters are infected
    cell-to-cell spread       — only partly NTCP-dependent
  These two are the 'floor' an entry inhibitor cannot get past.

Only four parameters were fitted
  r_oatp, Kd_ntcp   the two total-bile-acid fold-rise values (2 mg ×3.2, 10 mg ×13.0)
                    → NTCP occupancy is then 'derived': 2 mg 0.710, 10 mg 0.954
  dth_Id_immune     MYR301 2 mg week-48 virological response rate 71% (fitted value 0.02636/day)
  ALT_base          MYR301 2 mg week-48 ALT normalisation rate 51%
Everything else is a literature value or back-calculated from the untreated steady-state condition.

Predictions (not used in the fit): the whole 10 mg arm, the untreated control, week 96, the on- and
off-treatment response to Peg-IFN, combination synergy, lonafarnib/ritonavir, Peg-IFN lambda,
the ceiling of the entry-inhibitor class, the bile-acid therapeutic-index curve, 5-year fibrosis and HCC estimates.

An honest discrepancy: tied to the bile-acid anchor, this model separates 2 mg from 10 mg on the
virological endpoints far more widely than the actual trial does. Three possible readings, and the
one among them that is testable, are in section A13 of hdv_model_report.txt.

Files
  hdv_qsp_model.dot/.svg/.png  mechanistic map (18 clusters, 157 nodes)
  hdv_mrgsolve_model.R         the ODE model this app uses (33 compartments, 14 scenarios)
  hdv_reference_model.py       dependency-free Python reference implementation + virtual population
  hdv_model_report.txt         the computed results for every number above (A0–A14)
  hdv_references.md            source per parameter (PubMed-confirmed entries marked)

For education and research only; it must not be used for clinical decision-making.
")
  })
}

shinyApp(ui, server)
