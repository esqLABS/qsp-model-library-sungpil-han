# =============================================================================
#  Essential Tremor QSP — Shiny dashboard
#  ---------------------------------------------------------------------------
#  10 tabs.  The organising idea of the UI is the same as the model's: the user
#  does not set "tremor severity", they set LOOP GAIN, and then watch amplitude,
#  frequency and every clinical scale fall out of it.
#
#    1  Patient profile   Patient profile — G0, a_O, w_P, effector shares
#    2  Oscillator map     Oscillator — the bifurcation diagram, live
#    3  Drug PK         Plasma/brain concentrations, receptor occupancy
#    4  Gain decomposition       Gain decomposition: which phi is doing the work
#    5  Tremor amplitude/frequency Amplitude and frequency, with the waveform
#    6  Clinical scales       TETRAS / FTM / spiral / QUEST, and the log artefact
#    7  Scenario comparison   Side-by-side regimens
#    8  Surgery · neuromodulation   Lesion volume and DBS frequency sweeps
#    9  Differential diagnosis        Weight-loading test: ET vs enhanced physiological tremor
#   10  Organ systems · adverse effects HR, FEV1, sedation, ataxia, grip, HCO3, BMD
#
#  Run:  shiny::runApp("et_shiny_app.R")
#  Requires: shiny, mrgsolve, dplyr, tidyr, ggplot2, DT
#  The model code is sourced from et_mrgsolve_model_en.R (its scenario block is
#  skipped by the ET_SHINY guard at the bottom of that file's helpers).
# =============================================================================

library(shiny); library(mrgsolve); library(dplyr); library(tidyr)
library(ggplot2); library(DT)

# ---- load the model without executing the scenario block --------------------
src   <- readLines("et_mrgsolve_model_en.R")
cut   <- grep("^#  SCENARIOS", src)
if (length(cut)) src <- src[1:(cut[1] - 3)]
eval(parse(text = paste(src, collapse = "\n")))

W  <- 24 * 7
TH <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 13),
        legend.position = "bottom")
PAL <- c("#2c5aa0", "#c1121f", "#2d6a4f", "#e07b39", "#7b2cbf",
         "#0b7285", "#b08900", "#6a4c93")

# =============================================================================
#  UI
# =============================================================================
ui <- fluidPage(
  titlePanel("Essential Tremor QSP Model (oscillator formulation)"),
  tags$p(style = "color:#555;margin-top:-8px",
         HTML("Tremor is not a <b>level</b> but a <b>limit cycle</b>. ",
              "Amplitude is set by loop gain G, frequency by loop delay &tau;, and the ceiling on treatment is set by loop topology.")),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Patient"),
      sliderInput("G0", "Loop gain G0 (1.15 mild · 1.6 moderate · 6 severe)",
                  1.02, 12, 1.60, 0.02),
      sliderInput("AO", "Olivary fraction a_O  (ceiling for T-type blockers)", 0.05, 1.0, 0.35, 0.05),
      sliderInput("WP", "Peripheral loop weight w_P (upper limit of β₂ blockade)", 0.0, 0.8, 0.40, 0.05),
      sliderInput("AGE", "Age (changes frequency only)", 30, 90, 60, 1),
      sliderInput("HDG", "Cervical effector gain fraction HDG", 0.30, 1.20, 0.55, 0.05),
      sliderInput("STRESS", "Adrenergic drive multiplier (stress · caffeine)", 0, 9, 0, 0.5),
      checkboxInput("ASTHMA", "Comorbid asthma (computes β₂ contraindication)", FALSE),
      hr(),
      h4("Therapy"),
      selectInput("beta", "β blocker", c("None", "Propranolol", "Atenolol", "Nadolol")),
      conditionalPanel("input.beta == 'Propranolol'",
                       sliderInput("dprp", "Propranolol LA mg/day", 0, 320, 160, 20)),
      sliderInput("dprm", "Primidone mg/day (q8h)", 0, 750, 0, 25),
      sliderInput("dtop", "Topiramate mg/day (q12h)", 0, 400, 0, 25),
      sliderInput("dttb", "T-type calcium blocker mg/day", 0, 400, 0, 25),
      sliderInput("drinks", "Standard ethanol drinks (single occasion)", 0, 4, 0, 1),
      hr(),
      h4("Surgery · device"),
      sliderInput("vles", "Vim lesion volume mm³ (0 = none)", 0, 400, 0, 10),
      checkboxInput("dbs", "Turn on Vim DBS", FALSE),
      conditionalPanel("input.dbs == true",
                       sliderInput("fstim", "Stimulation frequency Hz", 10, 250, 130, 5),
                       sliderInput("vta", "Volume of tissue activated mm³", 30, 400, 250, 10)),
      sliderInput("btx", "Botulinum U (forearm)", 0, 150, 0, 10),
      sliderInput("fspill", "Diffusion fraction f_spill (0.15 guided · 0.45 unguided)",
                  0.05, 0.60, 0.15, 0.05),
      hr(),
      sliderInput("mload", "Wrist load mass kg (differential test)", 0, 1.0, 0, 0.1),
      sliderInput("weeks", "Simulation duration (weeks)", 1, 104, 24, 1),
      actionButton("go", "Run simulation", class = "btn-primary btn-block")
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("① Patient profile", br(),
                 fluidRow(column(6, tableOutput("prof")),
                          column(6, plotOutput("profplot", height = 300))),
                 hr(), htmlOutput("profnote")),
        tabPanel("② Oscillator map", br(),
                 plotOutput("bif", height = 380),
                 htmlOutput("bifnote"), hr(),
                 plotOutput("sqrtplot", height = 280)),
        tabPanel("③ Drug PK", br(),
                 plotOutput("pk", height = 340), hr(),
                 plotOutput("occ", height = 280)),
        tabPanel("④ Gain decomposition", br(),
                 plotOutput("phibar", height = 340),
                 htmlOutput("phinote"), hr(),
                 plotOutput("phitime", height = 280)),
        tabPanel("⑤ Amplitude · frequency", br(),
                 plotOutput("amp", height = 320), hr(),
                 fluidRow(column(6, plotOutput("freq", height = 280)),
                          column(6, plotOutput("wave", height = 280)))),
        tabPanel("⑥ Clinical scales", br(),
                 plotOutput("scales", height = 320), hr(),
                 htmlOutput("lognote"), plotOutput("logmap", height = 300)),
        tabPanel("⑦ Scenario comparison", br(),
                 DTOutput("cmp"), hr(), plotOutput("cmpplot", height = 340)),
        tabPanel("⑧ Surgery · neuromodulation", br(),
                 fluidRow(column(6, plotOutput("lesion", height = 320)),
                          column(6, plotOutput("dbsf", height = 320))),
                 hr(), htmlOutput("surgnote"), plotOutput("habit", height = 280)),
        tabPanel("⑨ Differential diagnosis", br(),
                 plotOutput("loadtest", height = 340),
                 htmlOutput("dxnote"), hr(), DTOutput("dxtab")),
        tabPanel("⑩ Organ systems · adverse effects", br(),
                 plotOutput("organs", height = 380), hr(),
                 fluidRow(column(6, plotOutput("btxwin", height = 300)),
                          column(6, tableOutput("adrtab"))))
      )
    )
  )
)

# =============================================================================
#  SERVER
# =============================================================================
server <- function(input, output, session) {

  pars <- reactive({
    p <- list(G0 = input$G0, AO = input$AO, AR = 1 - input$AO,
              WP = input$WP, WC = 1 - input$WP, AGE = input$AGE,
              HDG = input$HDG, STRESS = input$STRESS,
              ASTHMA = as.numeric(input$ASTHMA), MLOAD = input$mload,
              VLES = input$vles, LESION = as.numeric(input$vles > 0),
              DBSON = as.numeric(input$dbs), FSPILL = input$fspill)
    if (input$dbs) { p$FSTIM <- input$fstim; p$VTA <- input$vta }
    if (input$beta == "Propranolol") p$DPRP <- input$dprp
    p
  })

  evs <- reactive({
    e <- NULL; wk <- input$weeks
    add <- function(x) if (is.null(e)) x else c(e, x)
    if (input$beta == "Propranolol" && input$dprp > 0)
      e <- add(ev(amt = input$dprp, cmt = "A_PRPG", ii = 24, addl = wk*7 - 1))
    if (input$beta == "Atenolol")
      e <- add(ev(amt = 100, cmt = "A_ATNG", ii = 24, addl = wk*7 - 1))
    if (input$beta == "Nadolol")
      e <- add(ev(amt = 120, cmt = "A_NADG", ii = 24, addl = wk*7 - 1))
    if (input$dprm > 0)
      e <- add(ev(amt = input$dprm/3, cmt = "A_PRMG", ii = 8, addl = wk*21 - 1))
    if (input$dtop > 0)
      e <- add(ev(amt = input$dtop/2, cmt = "A_TOPG", ii = 12, addl = wk*14 - 1))
    if (input$dttb > 0)
      e <- add(ev(amt = input$dttb, cmt = "A_TTBG", ii = 24, addl = wk*7 - 1))
    if (input$drinks > 0)
      e <- add(ev(amt = 14*input$drinks, cmt = "A_ETHG", time = 24))
    if (input$btx > 0) {
      e <- add(ev(amt = input$btx*(1 - input$fspill), cmt = "A_BTXT"))
      e <- add(ev(amt = input$btx*input$fspill,       cmt = "A_BTXG"))
    }
    e
  })

  sim <- eventReactive(input$go, {
    m <- do.call(param, c(list(mod), pars()))
    dlt <- if (input$weeks <= 4) 0.25 else 2
    mrgsim(m, events = evs(), end = input$weeks*W, delta = dlt) %>% as.data.frame()
  }, ignoreNULL = FALSE)

  base_sim <- reactive({
    m <- do.call(param, c(list(mod), pars()))
    mrgsim(m, end = 8*W, delta = 4) %>% as.data.frame() %>% tail(1)
  })
  fin <- reactive({ d <- sim(); tail(d, min(25, nrow(d))) %>%
      summarise(across(where(is.numeric), mean)) })

  # ---- ① profile -----------------------------------------------------------
  output$prof <- renderTable({
    b <- base_sim()
    data.frame(Item = c("Loop gain G", "Bifurcation parameter μ = G−1", "Limit-cycle amplitude r*",
                        "Upper-limb amplitude A (cm)", "Physiological tremor floor (cm)",
                        "Central oscillator frequency (Hz)", "Limb mechanical resonance f₀ (Hz)",
                        "Observed dominant frequency (Hz)", "TETRAS performance (0-64)",
                        "Cervical μ_HD (negative = no head tremor)"),
               Value = sprintf("%.3f", c(b$GTOT, b$MU, sqrt(pmax(b$MU, 0)),
                                      b$A_UL, b$A_PHYS, b$FNEUR, b$F0M,
                                      b$F_OBS, b$TETRAS_PS, b$MU_HD)))
  }, striped = TRUE)

  output$profplot <- renderPlot({
    g <- seq(1.0, max(3, input$G0*1.3), length.out = 300)
    data.frame(G = g, A = ifelse(g > 1, sqrt(g - 1), 0)) %>%
      ggplot(aes(G, A)) + geom_line(linewidth = 1.2, colour = PAL[1]) +
      geom_vline(xintercept = 1, linetype = 2, colour = PAL[2]) +
      geom_point(data = data.frame(G = input$G0, A = sqrt(max(input$G0 - 1, 0))),
                 size = 4, colour = PAL[2]) +
      labs(title = "r* = √(G−1): where this patient sits on the curve",
           x = "Loop gain G", y = "Limit-cycle amplitude r*") + TH
  })
  output$profnote <- renderUI(HTML(
    "<div style='background:#eef3ff;padding:10px;border-left:4px solid #2c5aa0'>
     <b>G0</b> is this model's only severity axis. Because the curve has a <i>vertical tangent</i>
     at G=1, a patient near G=1 has their tremor <b>vanish</b> with the same drug, while a patient at
     <b>barely changes at all</b> at the far end. The responder/non-responder split repeatedly observed in trials is not two kinds of
     biology but a question of whether a single equation crosses a threshold.</div>"))

  # ---- ② bifurcation ------------------------------------------------------
  output$bif <- renderPlot({
    b <- base_sim(); f <- fin()
    g <- seq(0.85, max(3.5, input$G0*1.25), length.out = 400)
    df <- data.frame(G = g, r = ifelse(g > 1, sqrt(g - 1), 0))
    ggplot(df, aes(G, r)) +
      geom_line(linewidth = 1.2, colour = PAL[1]) +
      geom_hline(yintercept = 0, colour = "grey70") +
      geom_vline(xintercept = 1, linetype = 2, colour = PAL[2]) +
      annotate("text", x = 0.93, y = max(df$r)*0.9, label = "μ<0\nNo limit cycle",
               colour = PAL[2], size = 3.5) +
      geom_segment(aes(x = b$GTOT, xend = f$GTOT,
                       y = sqrt(max(b$MU, 0)), yend = sqrt(max(f$MU, 0))),
                   arrow = arrow(length = unit(0.22, "cm")), linewidth = 1,
                   colour = PAL[3]) +
      geom_point(data = data.frame(G = c(b$GTOT, f$GTOT),
                                   r = sqrt(pmax(c(b$MU, f$MU), 0)),
                                   w = c("Baseline", "On-treatment")),
                 aes(G, r, colour = w), size = 4) +
      scale_colour_manual(values = c("Baseline" = PAL[2], "On-treatment" = PAL[3]), name = NULL) +
      labs(title = "Supercritical Hopf bifurcation: treatment only pushes the patient left along the curve",
           x = "Loop gain G", y = "Amplitude envelope r*") + TH
  })
  output$bifnote <- renderUI({
    b <- base_sim(); f <- fin()
    HTML(sprintf("<div style='background:#f7f9fc;padding:8px'>G %.3f → %.3f
      (Δ %.1f%%) · μ %.3f → %.3f · amplitude %.3f → %.3f cm (<b>%+.1f%%</b>) ·
      TETRAS %.1f → %.1f (<b>%+.2f pts</b>)</div>",
      b$GTOT, f$GTOT, 100*(f$GTOT-b$GTOT)/b$GTOT, b$MU, f$MU,
      b$A_UL, f$A_UL, 100*(f$A_UL-b$A_UL)/b$A_UL,
      b$TETRAS_PS, f$TETRAS_PS, f$TETRAS_PS-b$TETRAS_PS))
  })
  output$sqrtplot <- renderPlot({
    # what the SAME gain reduction does across severities
    fr <- (fin()$GTOT/base_sim()$GTOT)
    g0 <- c(1.05, 1.15, 1.3, 1.6, 2.4, 4, 6, 9, 12)
    data.frame(G0 = g0,
               dA = 100*(sqrt(pmax(g0*fr - 1, 0))/sqrt(pmax(g0 - 1, 1e-9)) - 1)) %>%
      ggplot(aes(factor(G0), dA)) +
      geom_col(fill = PAL[1]) + geom_hline(yintercept = -50, linetype = 2, colour = PAL[2]) +
      labs(title = sprintf("Applying the gain ratio %.3f produced by the current prescription across severities", fr),
           x = "Patient's G0", y = "Amplitude change (%)") + TH
  })

  # ---- ③ PK ---------------------------------------------------------------
  output$pk <- renderPlot({
    sim() %>% select(time, C_PRP, C_PRM, C_PB, C_TOP, C_ETH) %>%
      pivot_longer(-time) %>% filter(value > 1e-9) %>%
      ggplot(aes(time/24, value, colour = name)) + geom_line(linewidth = 0.9) +
      scale_colour_manual(values = PAL, name = NULL) +
      facet_wrap(~name, scales = "free_y") +
      labs(title = "Plasma concentration (mg/L; ethanol g/L)", x = "Day", y = NULL) + TH
  })
  output$occ <- renderPlot({
    sim() %>% select(time, OCCB1, OCCB2, BLK_TT, P_EFF, REB) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time/24, value, colour = name)) + geom_line(linewidth = 0.9) +
      scale_colour_manual(values = PAL, name = NULL) +
      labs(title = "Receptor occupancy · blockade · GABA-A potentiation · rebound factor",
           subtitle = "OCCB2 alone produces both tremor efficacy and bronchoconstriction",
           x = "Day", y = "Fraction") + TH
  })

  # ---- ④ gain decomposition ----------------------------------------------
  output$phibar <- renderPlot({
    b <- base_sim(); f <- fin()
    data.frame(term = rep(c("φ_olive", "φ_cblthal", "φ_thal", "φ_ctx",
                            "φ_spindle", "φ_nmj"), 2),
               when = rep(c("Baseline", "On-treatment"), each = 6),
               v = c(b$PHI_OL, b$PHI_CBL, b$PHI_TH, 1, b$PHI_SPIN, b$PHI_NMJ,
                     f$PHI_OL, f$PHI_CBL, f$PHI_TH, 1, f$PHI_SPIN, f$PHI_NMJ)) %>%
      ggplot(aes(term, v, fill = when)) +
      geom_col(position = "dodge") + geom_hline(yintercept = 1, linetype = 2) +
      scale_fill_manual(values = c(PAL[2], PAL[3]), name = NULL) +
      labs(title = "Which term is actually doing the work", x = NULL, y = "Factor (1 = no change)") + TH
  })
  output$phinote <- renderUI(HTML(
    "<div style='background:#fff2f2;padding:10px;border-left:4px solid #c1121f'>
     <b>φ_olive</b> is a <i>parallel</i> branch, so even driving it to 0 cannot take Φ_C below a_R
     (the ceiling for T-type blockers). <b>φ_thal</b> is a <i>series</i> factor, so it multiplies the entire central loop
     (the reason surgery beats drugs is topology, not potency).</div>"))
  output$phitime <- renderPlot({
    sim() %>% select(time, PHI_OL, PHI_CBL, PHI_TH, PHI_SPIN, PHI_NMJ, GTOT) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time/24, value, colour = name)) + geom_line(linewidth = 0.9) +
      scale_colour_manual(values = PAL, name = NULL) +
      labs(title = "Time course of the gain terms", x = "Day", y = NULL) + TH
  })

  # ---- ⑤ amplitude / frequency / waveform --------------------------------
  output$amp <- renderPlot({
    sim() %>% select(time, A_UL, A_LC, A_PHYS, A_HD) %>% pivot_longer(-time) %>%
      ggplot(aes(time/24, value, colour = name)) + geom_line(linewidth = 0.9) +
      scale_colour_manual(values = PAL, name = NULL,
                          labels = c("Total observed amplitude (cm)", "Limit-cycle component (cm)",
                                     "Physiological floor (cm)", "Head (deg)")) +
      labs(title = "Amplitude: limit cycle + physiological floor (orthogonal sum)",
           x = "Day", y = NULL) + TH
  })
  output$freq <- renderPlot({
    sim() %>% select(time, FNEUR, F0M, F_OBS) %>% pivot_longer(-time) %>%
      ggplot(aes(time/24, value, colour = name)) + geom_line(linewidth = 1) +
      scale_colour_manual(values = PAL, name = NULL) +
      labs(title = "Frequency: the drug does not show up here",
           subtitle = "f is set by loop delay, and delay is set by age · Purkinje cell loss",
           x = "Day", y = "Hz") + TH
  })
  output$wave <- renderPlot({
    f <- fin(); fr <- f$F_OBS; amp <- f$A_UL
    t <- seq(0, 3, by = 0.002)
    data.frame(t = t, x = amp*sin(2*pi*fr*t)) %>%
      ggplot(aes(t, x)) + geom_line(colour = PAL[1]) +
      labs(title = sprintf("Waveform (%.2f Hz, %.2f cm)", fr, amp),
           x = "Time (s)", y = "cm") + TH
  })

  # ---- ⑥ scales ----------------------------------------------------------
  output$scales <- renderPlot({
    b <- base_sim(); f <- fin()
    data.frame(scale = rep(c("TETRAS-PS/64", "TETRAS-ADL/48", "FTM/144",
                             "Spiral/4", "Bain-Findley/10", "QUEST/100"), 2),
               when = rep(c("Baseline", "On-treatment"), each = 6),
               v = c(b$TETRAS_PS/64, b$TETRAS_ADL/48, b$FTM/144, b$SPIRAL/4,
                     b$BAINF/10, b$QUEST/100,
                     f$TETRAS_PS/64, f$TETRAS_ADL/48, f$FTM/144, f$SPIRAL/4,
                     f$BAINF/10, f$QUEST/100)) %>%
      ggplot(aes(scale, 100*v, fill = when)) + geom_col(position = "dodge") +
      scale_fill_manual(values = c(PAL[2], PAL[3]), name = NULL) +
      labs(title = "Clinical scales (% of each scale's maximum)", x = NULL, y = "%") + TH
  })
  output$lognote <- renderUI({
    b <- base_sim(); f <- fin()
    HTML(sprintf("<div style='background:#f2ecff;padding:10px;border-left:4px solid #5a189a'>
      The accelerometer reads <b>%+.1f%%</b>, yet TETRAS is <b>%+.2f points</b> (%+.1f%% of the baseline %.1f points).
      This is not a discrepancy but a <b>logarithm</b>. Because grade = 2 + 2·log₁₀(A), gaining 1 point requires
      the amplitude to shrink by <b>3.16-fold</b>.</div>",
      100*(f$A_UL-b$A_UL)/b$A_UL, f$TETRAS_PS-b$TETRAS_PS, b$TETRAS_PS,
      100*(f$TETRAS_PS-b$TETRAS_PS)/b$TETRAS_PS))
  })
  output$logmap <- renderPlot({
    a <- 10^seq(-1.5, 1.4, length.out = 300)
    data.frame(A = a, T = pmin(4, pmax(0, 2 + 2*log10(a)))) %>%
      ggplot(aes(A, T)) + geom_line(linewidth = 1.2, colour = PAL[5]) +
      scale_x_log10() +
      geom_point(data = data.frame(A = c(base_sim()$A_UL, fin()$A_UL),
                                   T = c(base_sim()$T_R, fin()$T_R),
                                   w = c("Baseline", "On-treatment")),
                 aes(colour = w), size = 4) +
      scale_colour_manual(values = c(PAL[2], PAL[3]), name = NULL) +
      labs(title = "Elble log transform: grade = 2 + 2·log₁₀(amplitude cm)",
           x = "Amplitude (cm, log axis)", y = "Item grade (0-4)") + TH
  })

  # ---- ⑦ scenario comparison ---------------------------------------------
  scen <- reactive({
    base_p <- pars(); base_p$DPRP <- 0
    mk <- function(nm, pp, e) {
      m <- do.call(param, c(list(mod), pp))
      d <- mrgsim(m, events = e, end = 8*W, delta = 4) %>% as.data.frame() %>% tail(7) %>%
        summarise(across(where(is.numeric), mean))
      data.frame(Regimen = nm, G = d$GTOT, `Amplitude_cm` = d$A_UL, TETRAS = d$TETRAS_PS,
                 QUEST = d$QUEST, HR = d$HR, Sedation = d$SED, Ataxia = d$ATAX,
                 `Grip strength` = d$GRIP, check.names = FALSE)
    }
    R <- function(a, c_, ii) ev(amt = a, cmt = c_, ii = ii, addl = ceiling(8*W/ii) - 1)
    bind_rows(
      mk("No treatment", base_p, NULL),
      mk("Propranolol 160", c(base_p, list(DPRP = 160)), R(160, "A_PRPG", 24)),
      mk("Atenolol 100", base_p, R(100, "A_ATNG", 24)),
      mk("Nadolol 120", base_p, R(120, "A_NADG", 24)),
      mk("Primidone 750", base_p, R(250, "A_PRMG", 8)),
      mk("Combination (propranolol+primidone)", c(base_p, list(DPRP = 160)),
         c(R(160, "A_PRPG", 24), R(250, "A_PRMG", 8))),
      mk("Topiramate 400", base_p, R(200, "A_TOPG", 12)),
      mk("T-type blocker 100", base_p, R(100, "A_TTBG", 24)),
      mk("MRgFUS 120 mm³", c(base_p, list(LESION = 1, VLES = 120)), NULL),
      mk("Vim DBS 130 Hz", c(base_p, list(DBSON = 1, FSTIM = 130)), NULL)
    ) %>% mutate(`AmplitudeChange_%` = 100*(Amplitude_cm/Amplitude_cm[1] - 1),
                 `TETRASChange` = TETRAS - TETRAS[1])
  })
  output$cmp <- renderDT(scen() %>% mutate(across(where(is.numeric), ~round(., 2))),
                         options = list(pageLength = 12, dom = "t"))
  output$cmpplot <- renderPlot({
    scen() %>% filter(Regimen != "No treatment") %>%
      ggplot(aes(reorder(Regimen, `AmplitudeChange_%`), `AmplitudeChange_%`)) +
      geom_col(fill = PAL[1]) + coord_flip() +
      geom_text(aes(label = sprintf("%+.1f%% / %+.1f pts", `AmplitudeChange_%`, `TETRASChange`)),
                hjust = 1.05, colour = "white", size = 3.2) +
      labs(title = "Amplitude change (%) and the TETRAS score change for the same outcome",
           x = NULL, y = "%") + TH
  })

  # ---- ⑧ surgery ----------------------------------------------------------
  output$lesion <- renderPlot({
    v <- seq(10, 400, by = 15)
    do.call(rbind, lapply(v, function(x) {
      d <- mrgsim(do.call(param, c(list(mod), pars(), list(LESION = 1, VLES = x))),
                  end = 4*W, delta = 24) %>% as.data.frame() %>% tail(1)
      data.frame(V = x, dA = 100*(d$A_UL/base_sim()$A_UL - 1), ATAX = d$ATAX)
    })) %>% pivot_longer(-V) %>%
      ggplot(aes(V, abs(value), colour = name)) + geom_line(linewidth = 1.1) +
      scale_colour_manual(values = c(PAL[3], PAL[2]), name = NULL,
                          labels = c("Ataxia index", "Amplitude reduction |%|")) +
      labs(title = "Lesion volume: efficacy saturates, ataxia does not",
           subtitle = "V50 efficacy 45 mm³ vs V50 ataxia 260 mm³ — the therapeutic window is derived",
           x = "Lesion volume (mm³)", y = NULL) + TH
  })
  output$dbsf <- renderPlot({
    f <- c(seq(10, 100, 10), 130, 160, 185, 220, 250)
    do.call(rbind, lapply(f, function(x) {
      d <- mrgsim(do.call(param, c(list(mod), pars(),
                                   list(DBSON = 1, FSTIM = x, VTA = input$vta))),
                  end = 4*W, delta = 24) %>% as.data.frame() %>% tail(1)
      data.frame(f = x, dA = 100*(d$A_UL/base_sim()$A_UL - 1))
    })) %>% ggplot(aes(f, dA)) + geom_line(linewidth = 1.1, colour = PAL[1]) +
      geom_hline(yintercept = 0, linetype = 2, colour = PAL[2]) +
      geom_vline(xintercept = 100, linetype = 3) +
      labs(title = "DBS frequency is a switch, not a dial",
           subtitle = "Hill 4, f₅₀ 80 Hz — at low frequency, entrainment actually worsens tremor",
           x = "Stimulation frequency (Hz)", y = "Amplitude change (%)") + TH
  })
  output$surgnote <- renderUI(HTML(
    "<div style='background:#eefaf1;padding:10px;border-left:4px solid #2d6a4f'>
     The same lesion volume produces both <b>therapeutic effect</b> and <b>gait ataxia</b>. Efficacy
     half-saturates at 45 mm³, but ataxia keeps rising up to 260 mm³, so enlarging the volume further yields
     almost no additional benefit while ataxia keeps accumulating. Habituation is a matter of <b>rewiring capacity</b>, not <b>duration</b>,
     and tremor returns only in patients who cross the threshold (and fairly abruptly when they do).</div>"))
  output$habit <- renderPlot({
    do.call(rbind, lapply(c(0.45, 0.75, 0.95, 1.00), function(rm) {
      do.call(rbind, lapply(c(0.25, 1, 2, 3, 4, 5), function(yr) {
        d <- mrgsim(do.call(param, c(list(mod), pars(),
                    list(LESION = 1, VLES = 120, REROUTE_MAX = rm))),
                    end = yr*8760, delta = 168) %>% as.data.frame() %>% tail(1)
        data.frame(rm = factor(rm), yr = yr,
                   dA = 100*(d$A_UL/base_sim()$A_UL - 1))
      }))
    })) %>% ggplot(aes(yr, dA, colour = rm)) +
      geom_line(linewidth = 1) + geom_point() +
      scale_colour_manual(values = PAL, name = "Rewiring capacity") +
      labs(title = "Habituation after MRgFUS: only patients who cross the threshold relapse",
           x = "Year", y = "Amplitude change (%)") + TH
  })

  # ---- ⑨ differential diagnosis -----------------------------------------
  dxd <- reactive({
    mk <- function(nm, pp, ml) {
      d <- mrgsim(do.call(param, c(list(mod), pars(), pp, list(MLOAD = ml))),
                  end = 4*W, delta = 24) %>% as.data.frame() %>% tail(1)
      data.frame(Condition = nm, Load_kg = ml, `f0_Hz` = d$F0M, `fcentral_Hz` = d$FNEUR,
                 `fobserved_Hz` = d$F_OBS, `A_cm` = d$A_UL,
                 `Limit cycle` = d$A_LC, Physiological = d$A_PHYS, check.names = FALSE)
    }
    bind_rows(
      mk("ET", list(), 0), mk("ET", list(), 0.5), mk("ET", list(), 1.0),
      mk("Enhanced physiological tremor (EPT)", list(G0 = 0.42, THYRO = 9), 0),
      mk("Enhanced physiological tremor (EPT)", list(G0 = 0.42, THYRO = 9), 0.5),
      mk("Enhanced physiological tremor (EPT)", list(G0 = 0.42, THYRO = 9), 1.0))
  })
  output$loadtest <- renderPlot({
    dxd() %>% pivot_longer(c(`f0_Hz`, `fobserved_Hz`)) %>%
      ggplot(aes(Load_kg, value, colour = Condition, linetype = name)) +
      geom_line(linewidth = 1.1) + geom_point(size = 2.5) +
      scale_colour_manual(values = c(PAL[1], PAL[2]), name = NULL) +
      scale_linetype_manual(values = c(2, 1), name = NULL) +
      labs(title = "Mass-loading test — the model reproduces the diagnostic test itself",
           subtitle = "The EPT peak falls as 1/√J because it is mechanical resonance; the ET peak barely moves because it is central delay",
           x = "Wrist load mass (kg)", y = "Frequency (Hz)") + TH
  })
  output$dxnote <- renderUI(HTML(
    "<div style='background:#f0f3f6;padding:10px;border-left:4px solid #495867'>
     One single equation produces both diseases. <b>ET</b> is a limit cycle with μ&gt;0, and its frequency is
     1/τ_loop; <b>EPT</b> has μ&lt;0, so there is no limit cycle, and because it is noise passing through
     mechanical resonance, its frequency is f₀. That is why <b>simply adding a weight</b> separates the two.</div>"))
  output$dxtab <- renderDT(dxd() %>% mutate(across(where(is.numeric), ~round(., 3))),
                           options = list(dom = "t"))

  # ---- ⑩ organs -----------------------------------------------------------
  output$organs <- renderPlot({
    sim() %>% select(time, HR, SBP, FEV1, SED, ATAX, GRIP, HCO3, COG) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time/24, value, colour = name)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 4) +
      scale_colour_manual(values = rep(PAL, 2), guide = "none") +
      labs(title = "Organ system indices", x = "Day", y = NULL) + TH
  })
  output$btxwin <- renderPlot({
    do.call(rbind, lapply(c(0.05, 0.15, 0.25, 0.35, 0.45, 0.60), function(fs) {
      d <- mrgsim(do.call(param, c(list(mod), pars(), list(FSPILL = fs))),
                  events = c(ev(amt = max(input$btx, 100)*(1 - fs), cmt = "A_BTXT"),
                             ev(amt = max(input$btx, 100)*fs,       cmt = "A_BTXG")),
                  end = 6*W, delta = 24) %>% as.data.frame()
      k <- which.min(d$A_UL)
      data.frame(f_spill = fs, `Tremor reduction` = -100*(d$A_UL[k]/base_sim()$A_UL - 1),
                 `Grip loss` = 100 - d$GRIP[k], check.names = FALSE)
    })) %>% pivot_longer(-f_spill) %>%
      ggplot(aes(f_spill, value, colour = name)) +
      geom_line(linewidth = 1.1) + geom_point() +
      scale_colour_manual(values = c(PAL[2], PAL[3]), name = NULL) +
      labs(title = "Botulinum: the therapeutic window comes from precision, not dose",
           subtitle = "Lowering the dose reduces both together; lowering f_spill separates them",
           x = "Diffusion fraction f_spill", y = "%") + TH
  })
  output$adrtab <- renderTable({
    f <- fin()
    data.frame(Metric = c("Heart rate (bpm)", "Systolic BP (mmHg)", "FEV₁ (L)",
                        "Sedation (0-100)", "Ataxia (0-100)", "Intoxication (0-100)",
                        "Grip strength (%)", "HCO₃⁻ (mmol/L)", "Body weight (kg)",
                        "Word-finding difficulty (0-100)", "Bone density index", "ALT (U/L)"),
               Value = sprintf("%.2f", c(f$HR, f$SBP, f$FEV1, f$SED, f$ATAX, f$INTOX,
                                      f$GRIP, f$HCO3, f$BW, f$COG, f$BMD, f$ALT)))
  }, striped = TRUE)
}

shinyApp(ui, server)
