## =============================================================================
##  CHRONIC HYPERKALAEMIA -- QSP DASHBOARD (Shiny)
##  Companion to hk_mrgsolve_model.R
## =============================================================================
##
##  RUN
##      install.packages(c("shiny","mrgsolve","dplyr","tidyr","ggplot2","DT","gridExtra"))
##      shiny::runApp("hk_shiny_app.R")
##
##  DESIGN NOTE -- WHY THIS APP IS LAID OUT THE WAY IT IS
##  ------------------------------------------------------
##  Almost every hyperkalaemia tool ever built plots serum potassium and stops.
##  Serum potassium is the LEAST informative quantity in this disease, because
##  it is a ratio of two things that move on different timescales and for
##  different reasons.  So every tab here is built to show the reader something
##  the potassium number itself hides:
##
##    Tab 1  Patient          what you are simulating
##    Tab 2  The two pools    K_total vs LAMrel -- the decomposition
##    Tab 3  Kidney           the kaliuresis reserve (FE_K), the thing that
##                            was moving for years while serum K sat still
##    Tab 4  Drug PK          exposures, MR occupancy, binder capture fraction
##    Tab 5  Acute rescue     what each emergency drug actually does, with a
##                            running "mmol removed" counter next to it
##    Tab 6  The dilemma      RAASi dose vs potassium, clinician in the loop
##    Tab 7  Scenarios        side-by-side long-term comparison
##    Tab 8  Population       who a binder rescues, and who it cannot
##    Tab 9  Validation       every trial the model was checked against,
##                            including the ones it disagrees with
##
## =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

mod <- mread_cache("hk", "hk_mrgsolve_model.R")

RUN_IN <- 250

thm <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92", colour = NA),
        legend.position = "bottom")

## Reference bands used on every potassium plot
k_bands <- function() list(
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 3.5, ymax = 5.0,
           fill = "#2f7f4f", alpha = 0.07),
  geom_hline(yintercept = 5.5, linetype = "22", colour = "#c07000"),
  geom_hline(yintercept = 6.0, linetype = "22", colour = "#b03030"),
  geom_hline(yintercept = 3.5, linetype = "22", colour = "#1f6bb8"))

sim_ss <- function(pars, tend = RUN_IN) {
  mod %>% param(pars) %>%
    mrgsim(end = tend, delta = 1, atol = 1e-8, rtol = 1e-8) %>%
    as.data.frame()
}

# =============================================================================
#  UI
# =============================================================================
ui <- fluidPage(
  titlePanel("Chronic Hyperkalaemia QSP Dashboard"),
  tags$p(style = "color:#555;margin-top:-8px",
    HTML("serum K is a RATIO, not a pool: &nbsp;<code>K_total = Ce&middot;V_ECF +
          Ci0&middot;LAMrel&middot;(Ce/Ce0)<sup>&alpha;</sup>&middot;V_ICF</code>
          &nbsp;&mdash;&nbsp; every tab below is built to show what that number hides")),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Patient"),
      sliderInput("gfr", "eGFR (mL/min/1.73)", 8, 100, 25, step = 1),
      sliderInput("intake", "Dietary K (mmol/day)", 20, 180, 80, step = 5),
      sliderInput("hco3", "Serum HCO3 (mmol/L) — 0 = model-derived",
                  0, 28, 0, step = 1),
      sliderInput("bw", "Body weight (kg)", 40, 130, 70, step = 1),

      hr(), h4("Chronic therapy"),
      sliderInput("ace", "ACE inhibitor (mg/day, lisinopril-equivalent)",
                  0, 40, 20, step = 5),
      sliderInput("mra", "MRA (mg/day, spironolactone-equivalent)",
                  0, 50, 25, step = 5),
      radioButtons("binder", "K binder",
                   c("none", "patiromer", "SZC"), inline = TRUE),
      conditionalPanel("input.binder == 'patiromer'",
        sliderInput("pat", "patiromer (g/day)", 0, 33.6, 16.8, step = 4.2)),
      conditionalPanel("input.binder == 'SZC'",
        sliderInput("szc", "SZC (g/day)", 0, 30, 10, step = 5)),
      sliderInput("fur", "furosemide (mg/day)", 0, 160, 0, step = 20),
      sliderInput("bic", "oral alkali (mmol/day)", 0, 140, 0, step = 10),
      checkboxInput("sglt2i", "SGLT2 inhibitor", FALSE),

      hr(), h4("Long run"),
      checkboxInput("titrate", "Clinician titration loop", TRUE),
      checkboxInput("progress", "eGFR progression (let eGFR decline)", TRUE),
      sliderInput("years", "years", 1, 10, 5, step = 1),

      hr(),
      helpText(HTML("<b>How to read this.</b> Look not at the serum K column, but at the
        <b>RAASi dose</b> and <b>mmol removed</b> columns. Almost every regimen eventually brings K below 5.5
        — the difference is <i>what was lost to achieve it</i>."))
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",

        tabPanel("1. Patient",
          br(), fluidRow(
            column(4, wellPanel(h4("Steady state"), tableOutput("tbl_ss"))),
            column(8, plotOutput("p_runin", height = "380px"))),
          helpText("Values after running the model for 250 days to reach steady state.
                    Note that ALDO and RASDN take several days to stabilise —
                    that is exactly why MRA-induced hyperkalaemia appears days later,
                    not immediately.")),

        tabPanel("2. Two Pools",
          br(), fluidRow(
            column(6, plotOutput("p_pool", height = "330px")),
            column(6, plotOutput("p_buffer", height = "330px"))),
          hr(),
          fluidRow(column(12, wellPanel(
            h4("Same number, opposite pool"),
            tableOutput("tbl_two")))),
          helpText(HTML("The exchange ratio (mmol per mmol/L) is <b>not a fitted value</b> but is derived from
            &alpha;. The chronic value of ~224 mmol/(mmol/L) reproduces the classic
            deficit nomogram (serum K 3.0 → 200-400 mmol deficit), and
            the acute value of ~66 explains why 40 mmol of IV KCl is dangerous."))),

        tabPanel("3. Kidney Reserve",
          br(), plotOutput("p_fek", height = "420px"),
          hr(), fluidRow(
            column(6, plotOutput("p_thresh", height = "300px")),
            column(6, tableOutput("tbl_thresh"))),
          helpText("The FE_K curve is the whole of this disease. Serum K stays flat for a long
                    time even as eGFR falls, and what was moving all along was
                    the fractional excretion rate. Serum K is a
                    'lagging indicator' of already-exhausted reserve.")),

        tabPanel("4. Drug PK/PD",
          br(), fluidRow(
            column(6, plotOutput("p_pk", height = "300px")),
            column(6, plotOutput("p_occ", height = "300px"))),
          hr(), plotOutput("p_binder_dr", height = "320px"),
          helpText("The binder's phi (capture fraction) has a structural ceiling:
                    no dose can remove more than phi_max x f_abs x dietary intake.
                    For a patient whose balance has broken down below this ceiling, dialysis is
                    not a preference but a requirement of mass conservation.")),

        tabPanel("5. Acute Rescue",
          br(),
          fluidRow(column(3, numericInput("acuteK", "Starting serum K", 6.8, 5.5, 9, 0.1)),
                   column(3, numericInput("insU", "insulin (U)", 10, 0, 20, 1)),
                   column(3, numericInput("dex", "dextrose (g)", 25, 0, 50, 5)),
                   column(3, numericInput("salb", "salbutamol neb (mg)", 20, 0, 20, 5))),
          plotOutput("p_acute", height = "330px"),
          plotOutput("p_removed", height = "260px"),
          helpText(HTML("<b>You must look at these two figures together.</b> If the line that falls
            in the figure above stays at 0 in the figure below, that drug
            <i>only moved the potassium rather than removing it</i> — it will inevitably come back."))),

        tabPanel("6. RAASi Dilemma",
          br(), plotOutput("p_dilemma", height = "440px"),
          hr(), tableOutput("tbl_dilemma"),
          helpText("Not calcium or insulin — this tab is the real clinical problem in this disease.
                    Toggle the titration loop off and on, and you can see that most of the harm
                    from hyperkalaemia arises not from arrhythmia but through 'discontinued
                    renoprotective therapy'.")),

        tabPanel("7. Scenario Comparison",
          br(), checkboxGroupInput("scen", "Regimens to compare",
            choices = c("no binder", "patiromer 16.8 g", "SZC 10 g",
                        "low-K diet 50", "furosemide 40", "SGLT2i",
                        "alkali 70", "no MRA"),
            selected = c("no binder", "patiromer 16.8 g", "low-K diet 50", "no MRA"),
            inline = TRUE),
          plotOutput("p_scen", height = "460px"),
          DTOutput("dt_scen")),

        tabPanel("8. Virtual Population",
          br(), fluidRow(
            column(4, sliderInput("pop_n", "Grid resolution (grid per axis)", 4, 12, 8)),
            column(8, helpText("eGFR x dietary K x HCO3 grid. The patients a binder actually
                                rescues are those who exceed the threshold by 'less than the
                                binder's ceiling', and that band is a computable
                                quantity."))),
          plotOutput("p_pop", height = "440px"),
          tableOutput("tbl_pop")),

        tabPanel("9. Validation",
          br(), h4("What the model got right, and wrong"),
          DTOutput("dt_valid"),
          hr(),
          h4("Out of scope"),
          tags$ul(
            tags$li("AKI on CKD — the most common cause of real severe hyperkalaemia"),
            tags$li("Rhabdomyolysis / tumour lysis — release from the ICF pool"),
            tags$li("Pseudohyperkalaemia (haemolysis, leukocytosis/thrombocytosis)"),
            tags$li("Dialysis — an intermittent, ultra-high-clearance removal term"),
            tags$li("Exercise-induced transient rise"),
            tags$li("The hazard is based on observational data and is not causal —
                     it is likely to overestimate the benefit of 'lowering the number'.")))
      )
    )
  )
)

# =============================================================================
#  SERVER
# =============================================================================
server <- function(input, output, session) {

  base_pars <- reactive({
    p <- list(GFR_SET = input$gfr, INTAKE = input$intake, BW = input$bw,
              ACE_TGT = input$ace, MRA_TGT = input$mra,
              FUR_RATE = input$fur, BIC_RATE = input$bic,
              SGLT2I = as.numeric(input$sglt2i),
              HCO3_SET = input$hco3,
              PAT_RATE = if (input$binder == "patiromer") input$pat else 0,
              SZC_RATE = if (input$binder == "SZC") input$szc else 0)
    p
  })

  run_in <- reactive(sim_ss(base_pars()))
  ss     <- reactive(dplyr::filter(run_in(), time == max(time)))

  # ---- Tab 1 ---------------------------------------------------------------
  output$tbl_ss <- renderTable({
    s <- ss()
    data.frame(
      quantity = c("serum K (mmol/L)", "total body K (mmol)", "LAMrel",
                   "urine K (mmol/day)", "FE_K (%)", "colonic K (mmol/day)",
                   "aldosterone (ng/dL)", "MR occupancy", "HCO3 (mmol/L)",
                   "arterial pH", "binder capture phi", "QRS (ms)",
                   "composite hazard"),
      value = c(sprintf("%.2f", s$K), sprintf("%.0f", s$KTOT),
                sprintf("%.3f", s$LAMrel), sprintf("%.1f", s$UK),
                sprintf("%.1f", s$FEK), sprintf("%.1f", s$COLK_EX),
                sprintf("%.1f", s$ALDO), sprintf("%.2f", s$MROCC),
                sprintf("%.1f", s$HCO3ser), sprintf("%.3f", s$PHA),
                sprintf("%.3f", s$PHI), sprintf("%.0f", s$QRS),
                sprintf("%.2f", s$HZTOT)))
  }, colnames = FALSE)

  output$p_runin <- renderPlot({
    run_in() %>%
      select(time, K, ALDO, RASDN, UK, FEK) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) +
      geom_line(linewidth = 0.9, colour = "#1f6bb8") +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "days of run-in", y = NULL,
           title = "Approach to steady state — the ASDN transporter pool (RASDN) is the slowest") +
      thm
  })

  # ---- Tab 2 ---------------------------------------------------------------
  output$p_pool <- renderPlot({
    s <- ss()
    data.frame(pool = factor(c("ECF", "fast ICF", "slow ICF (muscle)"),
                             levels = c("ECF", "fast ICF", "slow ICF (muscle)")),
               mmol = c(s$KE, s$KIF, s$KIS)) %>%
      ggplot(aes(pool, mmol, fill = pool)) +
      geom_col(width = 0.6) +
      geom_text(aes(label = sprintf("%.0f\n(%.1f%%)", mmol, 100*mmol/s$KTOT)),
                vjust = -0.2, size = 3.6) +
      scale_fill_manual(values = c("#1f6bb8", "#7fb0dd", "#c9dcef"), guide = "none") +
      expand_limits(y = max(s$KIS)*1.15) +
      labs(title = "Whole-body potassium distribution", subtitle = "Serum only sees the left-hand bar",
           y = "mmol", x = NULL) + thm
  })

  output$p_buffer <- renderPlot({
    s <- ss()
    kk <- seq(2, 8, by = 0.05)
    a <- 0.25; ci0 <- 140; ce0 <- 4.2
    v_e <- 0.20*input$bw; v_i <- 0.36*input$bw; v_f <- 0.25*v_i
    ktot <- kk*v_e + ci0*(kk/ce0)^a*v_i
    kref <- ce0*v_e + ci0*v_i
    data.frame(K = kk, delta = ktot - kref) %>%
      ggplot(aes(K, delta)) +
      geom_hline(yintercept = 0, colour = "grey60") +
      geom_line(linewidth = 1.1, colour = "#a02a72") +
      geom_vline(xintercept = s$K, linetype = "22") +
      annotate("point", x = 3.0,
               y = 3.0*v_e + ci0*(3.0/ce0)^a*v_i - kref, size = 3, colour = "#1f6bb8") +
      annotate("text", x = 3.0, y = 3.0*v_e + ci0*(3.0/ce0)^a*v_i - kref,
               label = "  K 3.0: classic nomogram, 200-400 mmol deficit",
               hjust = 0, size = 3.4, colour = "#1f6bb8") +
      labs(title = "Serum K vs whole-body potassium excess/deficit",
           subtitle = "The slope of this curve is how many mmol 'one unit of the lab value' represents",
           x = "serum K (mmol/L)", y = "whole-body K, deviation from normal (mmol)") +
      thm
  })

  output$tbl_two <- renderTable({
    s <- ss()
    data.frame(
      patient = c("A — mass-balance type (CKD, binder is the only answer)",
                  "B — redistribution type (DKA, insulin is the answer and binder is harmful)"),
      `serum K` = c(sprintf("%.2f", s$K), "4.7"),
      LAMrel = c(sprintf("%.3f", s$LAMrel), "0.86"),
      `whole-body K` = c(sprintf("%+.0f mmol", s$KTOT - (4.2*0.20*input$bw + 140*0.36*input$bw)),
                         "-400 mmol"),
      `After treatment` = c("Slow decline with binder, does not rebound",
                     "Insulin drops K to 2.7 — that is cardiac arrest"),
      check.names = FALSE)
  })

  # ---- Tab 3 ---------------------------------------------------------------
  reserve <- reactive({
    gfrs <- c(90, 75, 60, 45, 35, 25, 20, 15, 12, 10)
    arms <- list("none" = list(ACE_TGT = 0, MRA_TGT = 0),
                 "ACEi" = list(ACE_TGT = 20, MRA_TGT = 0),
                 "ACEi+MRA" = list(ACE_TGT = 20, MRA_TGT = 25),
                 "ACEi+MRA+binder" = list(ACE_TGT = 20, MRA_TGT = 25, PAT_RATE = 16.8))
    do.call(rbind, lapply(names(arms), function(a)
      do.call(rbind, lapply(gfrs, function(g) {
        r <- dplyr::filter(sim_ss(c(base_pars()[c("INTAKE", "BW")], arms[[a]],
                                    list(GFR_SET = g))), time == max(time))
        data.frame(arm = a, eGFR = g, K = r$K, FEK = r$FEK, UK = r$UK,
                   SCAP = r$SCAP, COL = r$COLK_EX)
      }))))
  })

  output$p_fek <- renderPlot({
    d <- reserve()
    p1 <- ggplot(d, aes(eGFR, K, colour = arm)) + k_bands() +
      geom_line(linewidth = 1.1) + geom_point(size = 1.6) +
      scale_x_reverse() +
      labs(title = "Serum K stays flat until reserve is exhausted",
           y = "serum K (mmol/L)", colour = NULL) + thm
    p2 <- ggplot(dplyr::filter(d, arm == "none"), aes(eGFR, FEK)) +
      geom_line(linewidth = 1.1, colour = "#2f7f8c") + geom_point(size = 1.6) +
      scale_x_reverse() +
      labs(title = "…meanwhile, this is what was moving (fractional excretion of potassium)",
           y = "FE_K (%)") + thm
    gridExtra::grid.arrange(p1, p2, ncol = 1)
  })

  thresholds <- reactive({
    arms <- list("none" = list(ACE_TGT = 0, MRA_TGT = 0),
                 "ACEi" = list(ACE_TGT = 20, MRA_TGT = 0),
                 "ACEi+MRA" = list(ACE_TGT = 20, MRA_TGT = 25),
                 "ACEi+MRA+binder" = list(ACE_TGT = 20, MRA_TGT = 25, PAT_RATE = 16.8))
    data.frame(arm = names(arms), threshold = vapply(arms, function(a) {
      lo <- 5; hi <- 120
      for (i in 1:30) {
        mid <- (lo + hi)/2
        k <- dplyr::filter(sim_ss(c(base_pars()[c("INTAKE","BW")], a,
                                    list(GFR_SET = mid))), time == max(time))$K
        if (k > 5.5) lo <- mid else hi <- mid
      }
      (lo + hi)/2
    }, numeric(1)))
  })

  output$p_thresh <- renderPlot({
    ggplot(thresholds(), aes(reorder(arm, threshold), threshold, fill = arm)) +
      geom_col(width = 0.6) + coord_flip() +
      geom_text(aes(label = sprintf("%.1f", threshold)), hjust = -0.15, size = 4) +
      scale_fill_brewer(palette = "Blues", guide = "none") +
      expand_limits(y = max(thresholds()$threshold)*1.2) +
      labs(title = "eGFR threshold above which K exceeds 5.5",
           subtitle = "MRA does not raise K — it raises this threshold",
           x = NULL, y = "eGFR (mL/min/1.73)") + thm
  })
  output$tbl_thresh <- renderTable(thresholds(), digits = 1)

  # ---- Tab 4 ---------------------------------------------------------------
  output$p_pk <- renderPlot({
    run_in() %>% select(time, CACE, CMRA, PATC, SZCP) %>%
      pivot_longer(-time) %>% dplyr::filter(time <= 20) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.9, colour = "#7040a0") +
      facet_wrap(~name, scales = "free_y") +
      labs(title = "Drug exposure (first 20 days)", x = "days", y = NULL) + thm
  })

  output$p_occ <- renderPlot({
    run_in() %>% select(time, MROCC, RASDN, PHI) %>%
      pivot_longer(-time) %>% dplyr::filter(time <= 30) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 1.0) +
      labs(title = "MR occupancy · ASDN transporter pool · binder capture fraction",
           subtitle = "RASDN's t1/2 of 1.5 days creates the delay in MRA-induced hyperkalaemia",
           x = "days", y = NULL, colour = NULL) + thm
  })

  output$p_binder_dr <- renderPlot({
    doses <- seq(0, 34, by = 2)
    d <- do.call(rbind, lapply(c("patiromer", "SZC"), function(drug)
      do.call(rbind, lapply(doses, function(dz) {
        pp <- base_pars(); pp$PAT_RATE <- 0; pp$SZC_RATE <- 0
        if (drug == "patiromer") pp$PAT_RATE <- dz else pp$SZC_RATE <- dz
        r <- dplyr::filter(sim_ss(pp), time == max(time))
        data.frame(drug = drug, dose = dz, K = r$K, phi = r$PHI)
      }))))
    ggplot(d, aes(dose, K, colour = drug)) + k_bands() +
      geom_line(linewidth = 1.1) +
      labs(title = "Binder dose-response (steady state)",
           subtitle = sprintf("Note the structural ceiling at a dietary intake of %d mmol/day", input$intake),
           x = "dose (g/day)", y = "serum K (mmol/L)", colour = NULL) + thm
  })

  # ---- Tab 5 ---------------------------------------------------------------
  acute <- reactive({
    p <- base_pars()
    s <- ss()
    dK <- (input$acuteK - s$K)*s$BUF_CHR
    v_e <- 0.20*input$bw
    i0 <- list(KE = s$KE + dK*v_e/s$BUF_CHR,
               KIF = s$KIF + dK*(1 - v_e/s$BUF_CHR)*0.25,
               KIS = s$KIS + dK*(1 - v_e/s$BUF_CHR)*0.75)
    run <- function(lab, ev_ = NULL, extra = list()) {
      m <- mod %>% param(c(p, extra)) %>% init(i0)
      o <- if (is.null(ev_)) mrgsim(m, end = 2, delta = 1/288)
           else mrgsim(m, ev_, end = 2, delta = 1/288)
      as.data.frame(o) %>% mutate(arm = lab)
    }
    rbind(
      run("None"),
      run("calcium 1 g", ev(amt = 0.25, cmt = "CAE")),
      run(sprintf("insulin %g U + D50 %g g", input$insU, input$dex),
          ev(amt = input$insU*1000/12, cmt = "INS") +
          ev(amt = input$dex/180.15*1000/16, cmt = "GLU")),
      run(sprintf("salbutamol %g mg", input$salb),
          ev(amt = input$salb*0.20/150, cmt = "B2C")),
      run("SZC 10 g TID", NULL, list(SZC_RATE = 30)))
  })

  output$p_acute <- renderPlot({
    ggplot(acute(), aes(time*24, K, colour = arm)) + k_bands() +
      geom_line(linewidth = 1.0) +
      labs(title = "Serum potassium — looking only at this, insulin looks like the best drug",
           x = "hours", y = "serum K (mmol/L)", colour = NULL) + thm
  })

  output$p_removed <- renderPlot({
    ggplot(acute(), aes(time*24, KREM, colour = arm)) +
      geom_line(linewidth = 1.0) +
      labs(title = "…potassium actually removed from the body (mmol)",
           subtitle = "Calcium, insulin, and salbutamol stay at exactly 0",
           x = "hours", y = "cumulative K removed (mmol)", colour = NULL) + thm
  })

  # ---- Tab 6/7 -------------------------------------------------------------
  SCEN <- list("no binder"        = list(),
               "patiromer 16.8 g" = list(PAT_RATE = 16.8),
               "SZC 10 g"         = list(SZC_RATE = 10),
               "low-K diet 50"    = list(INTAKE = 50),
               "furosemide 40"    = list(FUR_RATE = 40),
               "SGLT2i"           = list(SGLT2I = 1),
               "alkali 70"        = list(BIC_RATE = 70),
               "no MRA"           = list(MRA_TGT = 0))

  long_run <- reactive({
    sel <- if (input$tabs == "6. RAASi Dilemma") names(SCEN) else input$scen
    do.call(rbind, lapply(sel, function(a) {
      pp <- base_pars()
      pp$PAT_RATE <- 0; pp$SZC_RATE <- 0
      pp <- modifyList(pp, SCEN[[a]])
      pp$TITRATE <- as.numeric(input$titrate)
      pp$PROGRESS <- as.numeric(input$progress)
      o <- as.data.frame(mrgsim(param(mod, pp), end = input$years*365, delta = 7,
                                atol = 1e-8, rtol = 1e-8))
      o$arm <- a; o
    }))
  })

  output$p_dilemma <- renderPlot({
    long_run() %>% select(time, arm, K, RDo, GFR, HZTOT) %>%
      pivot_longer(c(K, RDo, GFR, HZTOT)) %>%
      mutate(name = recode(name, K = "serum K (mmol/L)",
                           RDo = "prescribed RAASi dose (fraction)",
                           GFR = "eGFR (mL/min/1.73)",
                           HZTOT = "composite hazard ratio")) %>%
      ggplot(aes(time/365, value, colour = arm)) +
      geom_line(linewidth = 0.95) +
      facet_wrap(~name, scales = "free_y") +
      labs(title = "RAASi–potassium dilemma",
           subtitle = "Look at the RAASi-dose panel, not the K panel",
           x = "years", y = NULL, colour = NULL) + thm
  })

  output$tbl_dilemma <- renderTable({
    long_run() %>% group_by(arm) %>% dplyr::filter(time == max(time)) %>%
      transmute(arm, `K (5y)` = K, `RAASi dose` = RDo, `eGFR (5y)` = GFR,
                `days K>5.5` = T55, `mean hazard` = CHAZ/time,
                `K removed (mmol)` = KREM) %>% ungroup()
  }, digits = 2)

  output$p_scen <- renderPlot({
    long_run() %>%
      ggplot(aes(time/365, K, colour = arm)) + k_bands() +
      geom_line(linewidth = 1.0) +
      labs(x = "years", y = "serum K (mmol/L)", colour = NULL,
           title = "Scenario comparison") + thm
  })
  output$dt_scen <- renderDT({
    long_run() %>% group_by(arm) %>% dplyr::filter(time == max(time)) %>%
      transmute(arm, K = round(K, 2), RAASi = round(RDo, 2),
                eGFR = round(GFR, 1), days_K_over_5.5 = round(T55),
                mean_hazard = round(CHAZ/time, 3),
                K_removed_mmol = round(KREM)) %>% ungroup()
  }, options = list(dom = "t", pageLength = 10))

  # ---- Tab 8 ---------------------------------------------------------------
  output$p_pop <- renderPlot({
    n <- input$pop_n
    grid <- expand.grid(GFR_SET = seq(15, 60, length.out = n),
                        INTAKE = seq(40, 140, length.out = n),
                        HCO3_SET = seq(16, 26, length.out = 4))
    res <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
      pp <- as.list(grid[i, ])
      a <- dplyr::filter(sim_ss(c(pp, list(ACE_TGT = 20, MRA_TGT = 25,
                                           BW = input$bw))), time == max(time))$K
      b <- dplyr::filter(sim_ss(c(pp, list(ACE_TGT = 20, MRA_TGT = 25,
                                           BW = input$bw, PAT_RATE = 16.8))),
                         time == max(time))$K
      data.frame(grid[i, ], K_no = a, K_bind = b)
    }))
    res %>% mutate(status = case_when(
        K_no <= 5.5 ~ "Binder not needed (already < 5.5)",
        K_bind <= 5.5 ~ "Binder rescues (rescued)",
        TRUE ~ "Binder insufficient (still > 5.5)")) %>%
      ggplot(aes(GFR_SET, INTAKE, fill = status)) +
      geom_tile(colour = "white") +
      facet_wrap(~sprintf("HCO3 %.0f", HCO3_SET)) +
      scale_fill_manual(values = c("Binder not needed (already < 5.5)" = "#cfe6d2",
                                   "Binder rescues (rescued)" = "#1f6bb8",
                                   "Binder insufficient (still > 5.5)" = "#b03030")) +
      labs(title = "Who does a binder actually rescue?",
           x = "eGFR (mL/min/1.73)", y = "dietary K (mmol/day)", fill = NULL) + thm
  })

  output$tbl_pop <- renderTable({
    data.frame(note = c(
      "The benefit of a binder is not uniform.",
      "It is concentrated in patients who exceed the threshold by 'less than the binder's ceiling'.",
      "That ceiling is phi_max x f_abs x dietary intake, and it is computable.",
      "For a patient beyond this ceiling, dialysis is not a preference but a requirement of mass conservation."))
  }, colnames = FALSE)

  # ---- Tab 9 ---------------------------------------------------------------
  output$dt_valid <- renderDT({
    data.frame(
      target = c("eGFR 60 steady-state K", "eGFR 45 steady-state K", "eGFR 30 steady-state K",
                 "Whole-body deficit at serum K 3.0", "OPAL-HK patiromer 16.8 g, 4 weeks",
                 "HARMONIZE SZC 5/10/15 g, 28 days", "HARMONIZE SZC 10 g TID, 48 hours",
                 "RALES spironolactone 25 mg", "FIDELIO finerenone 20 mg",
                 "Insulin 10 U + glucose", "Salbutamol 20 mg nebulised",
                 "ACE inhibitor alone", "Bicarbonate, HCO3 18→24",
                 "Hypoglycaemia incidence after insulin/glucose"),
      observed = c("4.40", "4.55", "4.75", "-200 ~ -400 mmol", "-1.01 mmol/L",
                   "4.8 / 4.5 / 4.4", "about -1.1 mmol/L", "+0.30 mmol/L",
                   "+0.23 mmol/L", "-0.6 ~ -1.0 (30-60 min)", "-0.6 ~ -1.0",
                   "+0.1 ~ +0.4", "about -0.2 ~ -0.3", "15-20%"),
      model = c("4.47", "4.61", "4.80", "-302 mmol", "-0.95 mmol/L",
                "4.89 / 4.63 / 4.49", "-0.95 mmol/L", "+0.30 mmol/L",
                "+0.10 mmol/L", "-0.85 (about 50 min)", "-0.74 (2-4 hours)",
                "+0.21 ~ +0.23", "-0.28 (entirely via the renal route)", "0% — not reproduced"),
      verdict = c("Predicted (held-out)", "Predicted (held-out)", "Predicted (held-out)",
                  "Predicted — derived from α, not fitted",
                  "Match", "Match (consistently +0.1 higher)", "Match", "Fitting anchor",
                  "Mismatch — assumed MR load is low (72% needed on back-calculation)",
                  "Match", "Match", "Match", "Fitting anchor",
                  "Mismatch — structural limitation, model not reliable"),
      check.names = FALSE)
  }, options = list(dom = "t", pageLength = 20))
}

shinyApp(ui, server)
