## ===========================================================================
##  tap_shiny_app.R
##  Toxic alcohol poisoning (methanol / ethylene glycol) — interactive QSP
##  dashboard on top of tap_mrgsolve_model_en.R
##
##  Toxic alcohol poisoning (methanol · ethylene glycol) QSP dashboard
##
##  Run:
##      shiny::runApp("tap_shiny_app.R")
##  (tap_mrgsolve_model_en.R must sit in the same directory)
##
##  ---------------------------------------------------------------------------
##  WHAT THIS APP IS FOR
##  ---------------------------------------------------------------------------
##  Not "look at the curves".  Each tab is built to make ONE argument from the
##  model checkable by hand:
##
##   1  Patient · exposure       who they are, what they drank, and when they arrived
##   2  Parent alcohol PK        the alcohol is an osmole with a reservoir, not a poison
##   3  Acid metabolites         the poison is manufactured, and it lags
##   4  Two gap clocks           OG falls, AG rises, the RATIO reads the clock
##   5  Acid-base                bicarbonate is titrated, and the space expands
##   6  Ion trapping             why pH is in the DELIVERY term, drawn as f_HA
##   7  ADH blockade arithmetic  how much flux each antidote actually leaves running
##   8  Kidney · crystals        the only threshold in the model, and its rate-dependence
##   9  Antidote PK · dialysis   the margin C/K that decides whether dialysis un-blocks
##  10  Clinical endpoints       logMAR, putaminal injury, P(death), P(blindness)
##  11  Scenario comparison      25 pre-built scenarios, side by side
##  12  Biomarkers               everything a clinician would chart, in one grid
##
##  DISCLAIMER: educational / research model.  Not for clinical use.
## ===========================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)

source("tap_mrgsolve_model_en.R", local = TRUE)

MOD <- tap_model()
SCN <- tap_scenarios()

THEME <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "#eef2f7", colour = NA),
        strip.text = element_text(face = "bold", size = 10),
        legend.position = "bottom")

PAL <- c("#2b6ca3", "#b03a3a", "#2e7d5b", "#7a3b8f", "#c98b00", "#555555")

## ---------------------------------------------------------------------------
##  helpers
## ---------------------------------------------------------------------------
long_plot <- function(d, vars, labs = NULL, ylab = "", ncol = 2, logy = FALSE) {
  keep <- intersect(vars, names(d))
  dd <- d %>% select(all_of(c("time", keep))) %>%
    pivot_longer(-time, names_to = "v", values_to = "y")
  if (!is.null(labs)) dd$v <- factor(dd$v, levels = keep, labels = labs[keep])
  p <- ggplot(dd, aes(time, y)) +
    geom_line(colour = PAL[1], linewidth = 0.8) +
    facet_wrap(~ v, scales = "free_y", ncol = ncol) +
    labs(x = "Time since ingestion (h)", y = ylab) + THEME
  if (logy) p <- p + scale_y_log10()
  p
}

hd_bands <- function(p, s) {
  pr <- s$param
  for (w in list(c(pr$HD1ON, pr$HD1OFF), c(pr$HD2ON, pr$HD2OFF))) {
    if (w[1] < 1e5) p <- p + annotate("rect", xmin = w[1], xmax = w[2],
                                      ymin = -Inf, ymax = Inf,
                                      fill = "#2b6ca3", alpha = 0.10)
  }
  if (pr$CRON < 1e5) p <- p + annotate("rect", xmin = pr$CRON, xmax = pr$CROFF,
                                       ymin = -Inf, ymax = Inf,
                                       fill = "#2e7d5b", alpha = 0.08)
  p
}

## Build a scenario from the sidebar controls.
ui_scen <- function(i) {
  hd <- list(); if (i$use_hd) hd <- list(c(i$hd_start, i$hd_start + i$hd_len))
  scen(
    label = "custom", tend = i$tend, wt = i$wt, gfr = i$gfr,
    meoh_gkg = if (i$agent %in% c("methanol", "both")) i$meoh_dose else 0,
    eg_mL    = if (i$agent %in% c("ethylene glycol", "both")) i$eg_dose else 0,
    etoh_gkg = i$etoh_coing,
    fom_start = if (i$antidote == "fomepizole") i$ant_start else NA,
    fom_boost = i$fom_boost,
    eth_start = if (i$antidote == "ethanol") i$ant_start else NA,
    eth_rate_mgkgh = i$eth_rate, eth_hd_mult = i$eth_hd_mult,
    hd = hd, crrt = if (i$use_crrt) c(i$cr_start, i$cr_start + i$cr_len) else NULL,
    bic = if (i$use_bic) list(start = i$bic_start, stop = i$bic_start + i$bic_len,
                              rate = i$bic_rate) else NULL,
    fol = if (i$use_fol) c(i$ant_start, i$tend) else NULL,
    thi = if (i$use_thi) c(i$ant_start, i$tend) else NULL,
    pyr = if (i$use_thi) c(i$ant_start, i$tend) else NULL,
    ca  = if (i$use_ca) list(start = i$ant_start, stop = i$tend, rate = 3) else NULL,
    noacid = if (i$ketone) 1 else 0)
}

## ===========================================================================
##  UI
## ===========================================================================
ui <- fluidPage(
  titlePanel("Toxic Alcohol Poisoning QSP Simulator (Methanol / Ethylene Glycol)"),
  tags$p(style = "color:#666;margin-top:-10px",
         HTML("The parent alcohol is an osmole and <b>the poison is the metabolite</b>. Every tab in this app ",
              "is built so that one claim following from the single structure <i>injury = ∫(ADH flux)dt</i> ",
              "can be checked by hand, one tab at a time. ",
              "<b>Educational / research model; not for clinical use.</b>")),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Patient"),
      sliderInput("wt", "Body weight (kg)", 40, 130, 70, 1),
      sliderInput("gfr", "Baseline GFR (multiple of normal)", 0.15, 1.2, 1.0, 0.05),
      sliderInput("tend", "Observation period (h)", 24, 168, 96, 12),

      h4("Exposure"),
      selectInput("agent", "Substance ingested",
                  c("methanol", "ethylene glycol", "both"), "methanol"),
      conditionalPanel("input.agent != 'ethylene glycol'",
        sliderInput("meoh_dose", "Methanol (g/kg)", 0, 2.5, 0.7, 0.05)),
      conditionalPanel("input.agent != 'methanol'",
        sliderInput("eg_dose", "Ethylene glycol, neat volume (mL)", 0, 250, 90, 5)),
      sliderInput("etoh_coing", "Co-ingested ethanol (g/kg)", 0, 2, 0, 0.1),
      checkboxInput("ketone", "Control: metabolite is a ketone (isopropanol equivalent)", FALSE),

      h4("Antidote"),
      selectInput("antidote", "Antidote", c("none", "fomepizole", "ethanol"),
                  "fomepizole"),
      sliderInput("ant_start", "Time of first dose (h)", 0, 36, 8, 1),
      conditionalPanel("input.antidote == 'fomepizole'",
        checkboxInput("fom_boost", "Re-dose q4h during dialysis", TRUE)),
      conditionalPanel("input.antidote == 'ethanol'",
        sliderInput("eth_rate", "Ethanol maintenance rate (mg/kg/h)", 40, 200, 100, 10),
        sliderInput("eth_hd_mult", "Rate multiplier during dialysis", 1, 4, 2.5, 0.5)),

      h4("Extracorporeal removal"),
      checkboxInput("use_hd", "Intermittent haemodialysis (IHD)", TRUE),
      conditionalPanel("input.use_hd",
        sliderInput("hd_start", "IHD start (h)", 0, 48, 9, 1),
        sliderInput("hd_len", "IHD duration (h)", 2, 12, 6, 1)),
      checkboxInput("use_crrt", "CRRT (40 mL/min)", FALSE),
      conditionalPanel("input.use_crrt",
        sliderInput("cr_start", "CRRT start (h)", 0, 48, 9, 1),
        sliderInput("cr_len", "CRRT duration (h)", 6, 72, 30, 2)),

      h4("Adjunctive treatment"),
      checkboxInput("use_bic", "Sodium bicarbonate infusion", TRUE),
      conditionalPanel("input.use_bic",
        sliderInput("bic_rate", "Bicarbonate rate (mmol/h)", 0, 80, 30, 5),
        sliderInput("bic_start", "Start (h)", 0, 36, 8, 1),
        sliderInput("bic_len", "Duration (h)", 4, 48, 18, 2)),
      checkboxInput("use_fol", "Folinic acid (leucovorin)", FALSE),
      checkboxInput("use_thi", "Thiamine + pyridoxine", FALSE),
      checkboxInput("use_ca", "Calcium supplementation (empirical)", FALSE),
      hr(),
      helpText(HTML("Shading = during dialysis. <br>",
                    "If it runs slowly, reduce 'Observation period'."))
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs", type = "tabs",

        ## ---------------------------------------------------------------- 1
        tabPanel("1 · Patient · exposure",
          br(),
          fluidRow(column(6, h4("Ingested amount and initial distribution"), tableOutput("t_expose")),
                   column(6, h4("Clinical course the model predicts"), tableOutput("t_course"))),
          hr(),
          HTML("<p><b>What to check on this tab.</b> For the same <i>total ingested amount</i>,
                the <i>time of arrival</i> dominates the outcome. The 'fraction of dose already
                oxidised' in the table below expresses how little is left for an antidote to
                block, and the conclusion that a late-presenting patient needs dialysis rather
                than an antidote falls out of that one cell.</p>"),
          plotOutput("p_overview", height = "420px")),

        ## ---------------------------------------------------------------- 2
        tabPanel("2 · Parent alcohol PK",
          br(), plotOutput("p_parent", height = "520px"),
          hr(),
          HTML("<p><b>Blocking ADH does not make the poison disappear; it turns it into a reservoir.</b>
                Methanol has no renal escape route, so under blockade its half-life stretches into the 40-hour range,
                whereas 20–30% of ethylene glycol leaves in the urine unchanged, so blockade alone
                finishes the job — but only while renal function lasts. The peripheral compartment at bottom left
                is the source of the post-dialysis rebound.</p>"),
          tableOutput("t_halflife")),

        ## ---------------------------------------------------------------- 3
        tabPanel("3 · Acid metabolites",
          br(), plotOutput("p_acidmet", height = "520px"),
          hr(),
          HTML("<p><b>The poison is manufactured, so it arrives late.</b> The peak of formic acid (or
                glycolic acid) comes several hours after the peak of the parent alcohol, and
                co-ingested ethanol widens that interval further — the 12–24 hour asymptomatic
                latent period is not a separate phenomenon but the output of this delay.
                Once the folate pool (THF) is consumed, the capacity to clear formic acid itself falls.</p>")),

        ## ---------------------------------------------------------------- 4
        tabPanel("4 · Two gap clocks",
          br(), plotOutput("p_gaps", height = "420px"),
          hr(),
          HTML("<p><b>This is reading the same reaction twice.</b> The osmolal gap (OG) falls by
                as much as the parent alcohol has disappeared, and the anion gap (AG) rises by as much
                acid as has been made. So the <b>OG/ΔAG ratio reads elapsed time, not
                severity</b>. The <b>sum</b> of the two is an approximate conserved measure of the
                ingested dose, but it is not exactly conserved: OG is divided by the alcohol's volume of distribution
                (~42 L) and AG by the bicarbonate space (36→70 L and above, depending on HCO3),
                so they differ. The direction and size of that deviation is on tab 5.</p>"),
          tableOutput("t_clock")),

        ## ---------------------------------------------------------------- 5
        tabPanel("5 · Acid-base",
          br(), plotOutput("p_ab", height = "520px"),
          hr(),
          HTML("<p><b>Bicarbonate is titrated, not poured in.</b> In this model the
                bicarbonate infusion switches itself off according to a pH target. Exogenous
                NaHCO&#8323; also raises sodium, so <b>the anion gap barely changes under
                bicarbonate treatment</b> — the alkali fixes the pH without erasing the diagnostic
                clue. That the bicarbonate space grows as the acidosis deepens
                (0.40 + 2.6/[HCO3] L/kg) is why the fall in HCO3 decelerates.</p>"),
          tableOutput("t_space")),

        ## ---------------------------------------------------------------- 6
        tabPanel("6 · Ion trapping",
          br(),
          fluidRow(column(7, plotOutput("p_trap", height = "400px")),
                   column(5, h4("f_HA(pH) = 1/(1+10^(pH−pKa))"),
                          tableOutput("t_trap"))),
          hr(),
          HTML("<p><b>pH is a prognostic factor because pH sits inside the delivery term.</b>
                Only un-ionised formic acid crosses the membrane, and that fraction rises about
                four-fold from pH 7.40 → 6.80. At the same blood formate, the
                <i>influx rate constant</i> into the CNS is four times higher and the <i>equilibrium CNS/plasma ratio</i> 1.8 times higher,
                so the brain load differs by two- to three-fold. Bicarbonate reverses both terms
                at once and, on top of that, keeps ionised formate from being reabsorbed in alkaline urine,
                so it raises renal clearance too — one equation, three terms.</p>"),
          plotOutput("p_cns", height = "360px")),

        ## ---------------------------------------------------------------- 7
        tabPanel("7 · ADH blockade arithmetic",
          br(),
          h4("Competitive blockade is diluted by its own substrate"),
          tableOutput("t_block"),
          hr(),
          HTML("<p><b>‘Fomepizole is 1000 times ethanol’ is the wrong statement at the
                bedside.</b> The competitive inhibition factor is (1 + ΣS/Km + I/Ki)/(1 + ΣS/Km),
                so it shrinks as substrate accumulates. At methanol 100 mg/dL, ethanol at the target
                concentration still leaves about 18% of the acid-generating flux running, and
                fomepizole 10 µg/mL leaves 0.6% — the honest ratio is about 30-fold.
                At methanol 400 mg/dL, even fomepizole leaves 2.4%.</p>"),
          plotOutput("p_flux", height = "380px")),

        ## ---------------------------------------------------------------- 8
        tabPanel("8 · Kidney · crystals",
          br(), plotOutput("p_kidney", height = "520px"),
          hr(),
          HTML("<p><b>The only threshold in this model.</b> Calcium oxalate precipitates only
                when it exceeds the solubility product, so <b>delivering the same total amount
                slowly reduces the injury</b> — a property formic acid does not have. And because
                1 mmol of deposited crystal takes 1 mmol of calcium with it, ionised calcium is the
                <b>stoichiometric shadow</b> of the crystal load, and the QT reads that
                shadow.</p>"),
          tableOutput("t_ss")),

        ## ---------------------------------------------------------------- 9
        tabPanel("9 · Antidote PK · dialysis",
          br(), plotOutput("p_ant", height = "420px"),
          hr(),
          h4("Can dialysis actually undo the antidote's blockade — a question of margin"),
          tableOutput("t_margin"),
          HTML("<p>What matters is not how quickly dialysis removes the antidote but
                <b>the margin C/K between the therapeutic concentration and the inhibition constant</b>. For fomepizole it is
                about 800-fold, for ethanol about 22-fold. A 6-hour session cannot spend fomepizole's margin
                (16 continuous hours would be needed), but it spends most of ethanol's
                margin. In other words, <b>‘increasing the antidote during dialysis’ is a safety
                margin for fomepizole and a quantitative necessity for ethanol</b> — the opposite
                direction from what is usually taught. On the graph above, check <code>FOMKIRAT</code>
                (= C/Ki) for yourself.</p>")),

        ## ---------------------------------------------------------------- 10
        tabPanel("10 · Clinical endpoints",
          br(), plotOutput("p_end", height = "520px"),
          hr(), tableOutput("t_end"),
          HTML("<p><b>Vision loss is a threshold, not an integral.</b> The retina fails for the
                same reason as the basal ganglia (failure of oxidative phosphorylation), so the driving force for
                optic nerve injury in this model is not <i>vitreous formic acid</i> but the <i>retinal ATP
                deficit</i>. Written with a linear driving force, every patient exposed for long enough
                to a low concentration went blind, and that did not match the clinic.</p>")),

        ## ---------------------------------------------------------------- 11
        tabPanel("11 · Scenario comparison",
          br(),
          checkboxGroupInput("cmp", "Scenarios to compare",
                             choices = names(SCN),
                             selected = c("M1_untreated", "M2_fom_early",
                                          "M3_fom_late", "M4_fom_hd", "M9_full"),
                             inline = TRUE),
          selectInput("cmp_var", "Variable to display",
                      c("pHart", "HCO3o", "AG", "OG", "OGoverAG", "MEOH", "EG",
                        "FORMmM", "GLYCmM", "FCNS", "FVIT", "FOMug", "ETOH",
                        "LOGMAR", "PDEATH", "PBLIND", "GFRpct", "CAIONo",
                        "OXALuM", "CAOXKID", "FLUXM", "FOMKIRAT", "LACTGAP"),
                      "pHart"),
          plotOutput("p_cmp", height = "420px"),
          hr(), h4("Summary of the 25 scenarios"), tableOutput("t_all")),

        ## ---------------------------------------------------------------- 12
        tabPanel("12 · Biomarker dashboard",
          br(), plotOutput("p_bio", height = "760px"),
          hr(),
          HTML("<p><b>The lactate gap</b> is a free clue. Glycolic acid
                cross-reacts with some blood-gas lactate electrodes (lactate oxidase), so the
                POC lactate comes back far higher than the laboratory lactate. This model computes that
                artefact as <code>LACTPOC</code> — grounds for suspecting ethylene glycol,
                not grounds for treating a lactic acidosis.</p>"))
      )
    )
  )
)

## ===========================================================================
##  SERVER
## ===========================================================================
server <- function(input, output, session) {

  cur <- reactive({
    s <- ui_scen(input)
    list(s = s, d = tap_run(MOD, s, dt = 0.1))
  })

  ## ------------------------------------------------------------------- 1 --
  output$t_expose <- renderTable({
    i <- input
    data.frame(
      Item = c("Body weight (kg)", "Baseline GFR (multiple)", "Methanol (mmol)",
               "Ethylene glycol (mmol)", "Co-ingested ethanol (mmol)",
               "Maximum expected acid equivalent of methanol (mmol)", "Acid equivalent of ethylene glycol (mmol)"),
      Value = c(i$wt, i$gfr,
             round(if (i$agent != "ethylene glycol") mmol_meoh(i$meoh_dose, i$wt) else 0),
             round(if (i$agent != "methanol") mmol_eg(i$eg_dose) else 0),
             round(mmol_etoh(i$etoh_coing, i$wt)),
             round(if (i$agent != "ethylene glycol") mmol_meoh(i$meoh_dose, i$wt) else 0),
             round(1.25 * if (i$agent != "methanol") mmol_eg(i$eg_dose) else 0)))
  }, digits = 2)

  output$t_course <- renderTable({
    d <- cur()$d; fin <- d[nrow(d), ]
    dose <- max(1e-9, mmol_meoh(input$meoh_dose, input$wt) +
                      mmol_eg(if (input$agent != "methanol") input$eg_dose else 0))
    data.frame(
      Metric = c("pH nadir", "HCO3 nadir (mM)", "AG peak", "OG peak",
               "Formic acid peak (mM)", "Glycolic acid peak (mM)",
               "Fraction of dose already oxidised (%)", "GFR nadir (%)",
               "Final logMAR", "P(death)", "P(permanent visual loss)"),
      Value = c(round(min(d$pHart), 3), round(min(d$HCO3o), 1),
             round(max(d$AG), 1), round(max(d$OG), 1),
             round(max(d$FORMmM), 2), round(max(d$GLYCmM), 2),
             round(100 * (fin$MEOHOX + fin$EGOX) / dose, 1),
             round(min(d$GFRpct)), round(fin$LOGMAR, 2),
             round(fin$PDEATH, 3), round(fin$PBLIND, 3)))
  }, digits = 3)

  output$p_overview <- renderPlot({
    x <- cur()
    hd_bands(long_plot(x$d, c("MEOH", "EG", "FORMmM", "GLYCmM", "pHart", "AG"),
                       labs = c(MEOH = "Methanol (mg/dL)", EG = "EG (mg/dL)",
                                FORMmM = "Formic acid (mM)", GLYCmM = "Glycolic acid (mM)",
                                pHart = "Arterial pH", AG = "Anion gap"),
                       ncol = 3), x$s)
  })

  ## ------------------------------------------------------------------- 2 --
  output$p_parent <- renderPlot({
    x <- cur(); d <- x$d
    d$MEOH_periph <- d$A_MEOH2 / (0.25 * input$wt) * 32.042 / 10
    d$EG_periph   <- d$A_EG2   / (0.27 * input$wt) * 62.068 / 10
    hd_bands(long_plot(d, c("MEOH", "MEOH_periph", "EG", "EG_periph",
                            "ETOH", "OG"),
                       labs = c(MEOH = "Methanol, central (mg/dL)",
                                MEOH_periph = "Methanol, peripheral (mg/dL) — source of the rebound",
                                EG = "EG, central (mg/dL)",
                                EG_periph = "EG, peripheral (mg/dL)",
                                ETOH = "Ethanol (mg/dL)",
                                OG = "Osmolal gap (mOsm)"),
                       ncol = 2), x$s)
  })

  output$t_halflife <- renderTable({
    d <- cur()$d
    f <- function(v) {
      i0 <- which.max(v); m <- seq_along(v) > i0 + 40 & v > 2
      if (sum(m) < 10) return(NA_real_)
      -log(2) / coef(lm(log(v[m]) ~ d$time[m]))[2]
    }
    data.frame(
      Substance = c("Methanol", "Ethylene glycol"),
      `Peak (mg/dL)` = c(round(max(d$MEOH), 1), round(max(d$EG), 1)),
      `Observed elimination half-life (h)` = round(c(f(d$MEOH), f(d$EG)), 1),
      `Time to <20 mg/dL (h)` = c(
        suppressWarnings(min(d$time[d$MEOH < 20 & d$time > d$time[which.max(d$MEOH)]])),
        suppressWarnings(min(d$time[d$EG   < 20 & d$time > d$time[which.max(d$EG)]]))),
      check.names = FALSE)
  }, digits = 1)

  ## ------------------------------------------------------------------- 3 --
  output$p_acidmet <- renderPlot({
    x <- cur()
    hd_bands(long_plot(x$d, c("FORMmM", "FORMmgL", "FCNS", "FVIT",
                              "GLYCmM", "THFo"),
                       labs = c(FORMmM = "Plasma formic acid (mM)",
                                FORMmgL = "Plasma formic acid (mg/L)",
                                FCNS = "CNS formic acid (mM) — the injuring quantity",
                                FVIT = "Vitreous · retinal formic acid (mM)",
                                GLYCmM = "Plasma glycolic acid (mM)",
                                THFo = "Hepatic THF pool (relative)"),
                       ncol = 2), x$s)
  })

  ## ------------------------------------------------------------------- 4 --
  output$p_gaps <- renderPlot({
    x <- cur(); d <- x$d
    dd <- d %>% select(time, OG, DAG, GAPSUM, OGoverAG) %>%
      pivot_longer(-time, names_to = "v", values_to = "y")
    dd$v <- factor(dd$v, levels = c("OG", "DAG", "GAPSUM", "OGoverAG"),
                   labels = c("Osmolal gap OG (mOsm)",
                              "ΔAG (mEq) — same flux, opposite sign",
                              "OG + ΔAG — approximately conserved (dose estimate)",
                              "OG/ΔAG — the clock (elapsed time)"))
    p <- ggplot(dd, aes(time, y)) +
      geom_line(colour = PAL[2], linewidth = 0.9) +
      facet_wrap(~ v, scales = "free_y", ncol = 2) +
      labs(x = "Time (h)", y = NULL) + THEME
    hd_bands(p, x$s)
  })

  output$t_clock <- renderTable({
    d <- cur()$d
    tt <- c(0.5, 1, 2, 4, 6, 8, 12, 16, 20, 24, 36, 48)
    tt <- tt[tt <= max(d$time)]
    i <- vapply(tt, function(t) which.min(abs(d$time - t)), 1L)
    out <- d[i, c("time", "MEOH", "EG", "FORMmM", "GLYCmM", "OG", "DAG",
                  "OGoverAG", "GAPSUM", "pHart")]
    names(out) <- c("t (h)", "MeOH", "EG", "Formic acid", "Glycolic acid",
                    "OG", "ΔAG", "OG/ΔAG (clock)", "Sum", "pH")
    out
  }, digits = 2)

  ## ------------------------------------------------------------------- 5 --
  output$p_ab <- renderPlot({
    x <- cur()
    hd_bands(long_plot(x$d, c("pHart", "HCO3o", "PACO2o", "AG", "NAo", "LACTo"),
                       labs = c(pHart = "Arterial pH", HCO3o = "HCO3 (mM)",
                                PACO2o = "PaCO2 (mmHg)", AG = "Anion gap",
                                NAo = "Sodium (mM) — the price of the alkali load",
                                LACTo = "Lactate (mM)"),
                       ncol = 3), x$s)
  })

  output$t_space <- renderTable({
    h <- c(24, 20, 16, 12, 10, 8, 6, 4)
    data.frame(`HCO3 (mM)` = h,
               `Bicarbonate space (L)` = round((0.40 + 2.6 / h) * input$wt, 1),
               `Space / Vd(methanol)` = round((0.40 + 2.6 / h) * input$wt /
                                            (0.60 * input$wt), 3),
               check.names = FALSE)
  }, digits = 2)

  ## ------------------------------------------------------------------- 6 --
  output$p_trap <- renderPlot({
    tb <- tap_trapping_table()
    dd <- tb %>% select(pH, rel_entry, eq_CNS_over_plasma) %>%
      pivot_longer(-pH, names_to = "v", values_to = "y")
    dd$v <- factor(dd$v, levels = c("rel_entry", "eq_CNS_over_plasma"),
                   labels = c("CNS influx rate constant (pH 7.40 = 1)",
                              "Equilibrium CNS/plasma formic acid ratio"))
    ggplot(dd, aes(pH, y, colour = v)) +
      geom_line(linewidth = 1) + geom_point() +
      scale_x_reverse() + scale_colour_manual(values = PAL[1:2], name = NULL) +
      labs(x = "Arterial pH (acidic to the right)", y = "Relative value") + THEME
  })

  output$t_trap <- renderTable({
    tb <- tap_trapping_table()
    names(tb) <- c("pH", "f_HA plasma", "Relative influx rate", "f_HA brain", "Equilibrium CNS/plasma")
    tb
  }, digits = 4)

  output$p_cns <- renderPlot({
    x <- cur()
    hd_bands(long_plot(x$d, c("fHAratio", "FCNS", "FVIT", "ATPo"),
                       labs = c(fHAratio = "f_HA / f_HA(7.40) — delivery multiple",
                                FCNS = "CNS formic acid (mM)",
                                FVIT = "Vitreous formic acid (mM)",
                                ATPo = "CNS ATP (1 = normal)"),
                       ncol = 4), x$s)
  })

  ## ------------------------------------------------------------------- 7 --
  output$t_block <- renderTable({
    tb <- tap_blockade_table(MOD)
    names(tb) <- c("Methanol (mg/dL)", "Ethanol (mg/dL)", "Fomepizole (µg/mL)",
                   "Inhibition factor", "Flux remaining (%)", "Flux (mmol/h)")
    tb
  }, digits = 3)

  output$p_flux <- renderPlot({
    x <- cur()
    hd_bands(long_plot(x$d, c("FLUXM", "FLUXE", "INHFAC", "FLUXLEFT"),
                       labs = c(FLUXM = "Methanol oxidation flux (mmol/h)",
                                FLUXE = "EG oxidation flux (mmol/h)",
                                INHFAC = "Inhibition factor (fold)",
                                FLUXLEFT = "Flux remaining (%)"),
                       ncol = 2), x$s)
  })

  ## ------------------------------------------------------------------- 8 --
  output$p_kidney <- renderPlot({
    x <- cur()
    hd_bands(long_plot(x$d, c("OXALuM", "SSTUB", "CAOXKID", "PTINJo",
                              "GFRpct", "CAIONo", "QTCms", "SSPLAS"),
                       labs = c(OXALuM = "Plasma oxalate (µM)",
                                SSTUB = "Tubular supersaturation (multiple of Ksp)",
                                CAOXKID = "Renal CaOx deposition (mmol, cumulative)",
                                PTINJo = "Proximal tubular injury (0–1)",
                                GFRpct = "GFR (% of baseline)",
                                CAIONo = "Ionised calcium (mM)",
                                QTCms = "QTc (ms)",
                                SSPLAS = "Plasma supersaturation"),
                       ncol = 4), x$s)
  })

  output$t_ss <- renderTable({
    ox <- c(2, 5, 10, 20, 40, 80, 160)
    ksp <- 2.32e-3
    data.frame(`Oxalate (µM)` = ox,
               `Plasma SS` = round(1.10 * ox / 1000 / ksp, 2),
               `Tubular SS (×50)` = round(1.10 * ox / 1000 / ksp * 50, 0),
               Verdict = ifelse(1.10 * ox / 1000 / ksp > 4, "Tubule + tissue",
                       ifelse(1.10 * ox / 1000 / ksp * 50 > 8, "Tubule", "None")),
               check.names = FALSE)
  }, digits = 2)

  ## ------------------------------------------------------------------- 9 --
  output$p_ant <- renderPlot({
    x <- cur()
    hd_bands(long_plot(x$d, c("FOMug", "FOMKIRAT", "ETOH", "FLUXM"),
                       labs = c(FOMug = "Fomepizole (µg/mL)",
                                FOMKIRAT = "Fomepizole C/Ki (margin)",
                                ETOH = "Ethanol (mg/dL)",
                                FLUXM = "Acid-generating flux (mmol/h)"),
                       ncol = 2, logy = FALSE), x$s)
  })

  output$t_margin <- renderTable({
    tb <- tap_margin_table(MOD)
    names(tb) <- c("Antidote", "Therapeutic concentration (mM)", "Inhibition constant (mM)", "Dialysis CL (L/h)",
                   "Central volume (L)", "Margin C/K (fold)", "Dialysis k (1/h)",
                   "Dialysis time needed to spend the whole margin (h)")
    tb
  }, digits = 4)

  ## ------------------------------------------------------------------ 10 --
  output$p_end <- renderPlot({
    x <- cur()
    hd_bands(long_plot(x$d, c("LOGMAR", "OPTICo", "PUTo", "ATPo",
                              "COMA", "PDEATH", "PBLIND", "GFRpct"),
                       labs = c(LOGMAR = "Visual acuity logMAR (0 = 20/20)",
                                OPTICo = "Optic nerve injury index",
                                PUTo = "Basal ganglia (putamen) injury index",
                                ATPo = "CNS ATP",
                                COMA = "Depressed consciousness index",
                                PDEATH = "Cumulative P(death)",
                                PBLIND = "Cumulative P(permanent visual loss)",
                                GFRpct = "GFR (%)"),
                       ncol = 4), x$s)
  })

  output$t_end <- renderTable({
    d <- cur()$d; fin <- d[nrow(d), ]
    data.frame(
      Endpoint = c("Final logMAR", "logMAR interpretation", "Basal ganglia injury index",
                     "P(death)", "P(permanent visual loss)", "GFR nadir (%)",
                     "Unsurvivable pH reached", "Observation end time (h)"),
      Value = c(sprintf("%.2f", fin$LOGMAR),
             if (fin$LOGMAR < 0.1) "Close to normal" else
             if (fin$LOGMAR < 0.5) "Mild reduction" else
             if (fin$LOGMAR < 1.0) "Moderate reduction" else
             if (fin$LOGMAR < 1.6) "Within the legal blindness range" else "Loss of light perception level",
             sprintf("%.3f", fin$PUTo), sprintf("%.3f", fin$PDEATH),
             sprintf("%.3f", fin$PBLIND), sprintf("%.0f", min(d$GFRpct)),
             if (any(d$FATAL > 0)) "Yes" else "No",
             sprintf("%.1f", fin$time)))
  })

  ## ------------------------------------------------------------------ 11 --
  cmp_data <- reactive({
    req(length(input$cmp) > 0)
    do.call(rbind, lapply(input$cmp, function(k) {
      d <- tap_run(MOD, SCN[[k]], dt = 0.2)
      d$scenario <- k
      d
    }))
  })

  output$p_cmp <- renderPlot({
    d <- cmp_data(); v <- input$cmp_var
    ggplot(d, aes(time, .data[[v]], colour = scenario)) +
      geom_line(linewidth = 0.85) +
      labs(x = "Time (h)", y = v, colour = NULL) +
      scale_colour_manual(values = rep(PAL, 5)) + THEME
  })

  output$t_all <- renderTable({
    a <- tap_all(MOD, dt = 0.25)
    a$label <- NULL
    a
  }, digits = 3)

  ## ------------------------------------------------------------------ 12 --
  output$p_bio <- renderPlot({
    x <- cur()
    hd_bands(long_plot(x$d,
      c("MEOH", "EG", "FORMmM", "GLYCmM", "OXALuM", "FOMug", "ETOH",
        "pHart", "HCO3o", "PACO2o", "AG", "OG", "OGoverAG", "NAo",
        "CAIONo", "QTCms", "LACTo", "LACTPOC", "LACTGAP", "pHur",
        "GFRpct", "CAOXKID", "ATPo", "LOGMAR"),
      labs = c(MEOH = "Methanol mg/dL", EG = "EG mg/dL",
               FORMmM = "Formic acid mM", GLYCmM = "Glycolic acid mM",
               OXALuM = "Oxalate µM", FOMug = "Fomepizole µg/mL",
               ETOH = "Ethanol mg/dL", pHart = "pH", HCO3o = "HCO3",
               PACO2o = "PaCO2", AG = "AG", OG = "OG",
               OGoverAG = "OG/ΔAG clock", NAo = "Na",
               CAIONo = "Ionised Ca", QTCms = "QTc",
               LACTo = "Lactate (laboratory)", LACTPOC = "Lactate (POC, with interference)",
               LACTGAP = "Lactate gap", pHur = "Urine pH",
               GFRpct = "GFR %", CAOXKID = "Renal CaOx mmol",
               ATPo = "CNS ATP", LOGMAR = "logMAR"),
      ncol = 4), x$s)
  })
}

shinyApp(ui, server)
