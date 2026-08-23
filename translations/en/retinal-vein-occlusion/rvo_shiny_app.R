# =============================================================================
#  rvo_shiny_app.R
#  Retinal vein occlusion QSP model — interactive Shiny dashboard
#  Retinal Vein Occlusion QSP Model — Interactive Dashboard
# =============================================================================
#
#  WHAT THIS APP IS FOR
#  --------------------
#  The app exists to make ONE question answerable by moving sliders:
#
#      is this eye's residual oedema a PHARMACOLOGICAL problem
#      or a PRESSURE problem?
#
#  Everything is organised around that. The "Pressure vs permeability" tab shows
#  P_c against the critical capillary pressure P_c*, and the verdict banner on the
#  patient tab says, in words, whether any anti-VEGF dose can dry this macula.
#  Increasing the dose, shortening the interval, or switching the agent moves the
#  permeability arm; only the blood-pressure slider and the collateral parameters
#  move the pressure arm. An app in which every control does the same thing would
#  hide the point of the model.
#
#  USAGE
#  -----
#      Rscript -e 'shiny::runApp("rvo_shiny_app.R", port = 8080)'
#
#  This file sources rvo_mrgsolve_model.R from the same directory for the model
#  specification, the drug library, the phenotypes and the treat-and-extend
#  controller, so the app and the batch analyses cannot drift apart.
# =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)

# -----------------------------------------------------------------------------
#  Load the model. Everything below assumes rvo_mrgsolve_model.R has defined:
#  mod, DRUGS, PHENOTYPES, molar_dose(), drug_params(), run_case(),
#  run_reactive(), q4w(), suppression_duration(), critical_pressure(),
#  critical_rv(), BASE_START, HORIZON
# -----------------------------------------------------------------------------
model_file <- "rvo_mrgsolve_model.R"
if (!file.exists(model_file)) {
  stop("rvo_mrgsolve_model.R must be in the working directory.")
}
local({
  # source everything except the interactive() demonstration block
  src <- readLines(model_file)
  cut <- grep("^if \\(interactive\\(\\)\\) \\{", src)
  if (length(cut)) src <- src[seq_len(cut[1] - 1)]
  eval(parse(text = paste(src, collapse = "\n")), envir = globalenv())
})

THEME <- theme_bw(base_size = 12) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95"))

PAL <- c("#1565c0", "#c62828", "#2e7d32", "#ef6c00", "#6a1b9a", "#00838f",
         "#ad1457", "#455a64")

PHENO_LABELS <- c(
  "crvo_nonisch" = "Non-ischaemic CRVO",
  "crvo_isch"    = "Ischaemic CRVO",
  "brvo"         = "Branch retinal vein occlusion (BRVO)",
  "crvo_mild"    = "Mild occlusion")

REGIMEN_LABELS <- c(
  "none" = "Untreated (natural history)",
  "tae"  = "Treat-and-extend",
  "prn"  = "PRN (monthly PRN)",
  "q4w"  = "Fixed q4w (fixed monthly)",
  "q8w"  = "Fixed q8w (fixed 8-weekly)",
  "dex"  = "Dexamethasone implant q6mo",
  "perf" = "Infinite blockade experiment")

# =============================================================================
#  UI
# =============================================================================
ui <- fluidPage(
  titlePanel(HTML(paste0(
    "<b>Retinal Vein Occlusion QSP Simulator</b> &nbsp;",
    "<span style='font-size:15px;color:#555'>Retinal Vein Occlusion QSP simulator ",
    "&mdash; J<sub>v</sub> = L<sub>p</sub>&middot;S&middot;[(P<sub>c</sub>&minus;P<sub>t</sub>) ",
    "&minus; &sigma;(&pi;<sub>c</sub>&minus;&pi;<sub>t</sub>)]</span>"))),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      h4("1 · The Eye"),
      selectInput("pheno", "Phenotype",
                  choices = setNames(names(PHENO_LABELS), PHENO_LABELS),
                  selected = "crvo_nonisch"),
      sliderInput("occ0", "Occlusion severity, OCC0",
                  min = 1, max = 45, value = 10, step = 0.5),
      sliderInput("occres", "Non-recanalisable residual, OCC_RES",
                  min = 0, max = 0.95, value = 0.55, step = 0.05),
      sliderInput("collmax", "Collateral ceiling, COLL_MAX",
                  min = 0.1, max = 5, value = 0.7, step = 0.1),

      h4("2 · Systemic"),
      sliderInput("pa", "Retinal arteriolar inflow pressure, P_a (mmHg)",
                  min = 28, max = 50, value = 40, step = 0.5),
      helpText(HTML(paste0("<span style='font-size:11px'>P<sub>a</sub> is the only ",
        "slider that moves the <b>pressure arm</b>. Everything under 'treatment' ",
        "moves the permeability arm.</span>"))),
      sliderInput("iop0", "Baseline intraocular pressure, IOP (mmHg)", min = 8, max = 30,
                  value = 15, step = 1),

      h4("3 · Treatment"),
      selectInput("regimen", "Regimen",
                  choices = setNames(names(REGIMEN_LABELS), REGIMEN_LABELS),
                  selected = "tae"),
      selectInput("drug", "Agent",
                  choices = setNames(DRUGS$drug, DRUGS$label),
                  selected = "aflibercept"),
      sliderInput("start", "Day of first treatment",
                  min = 7, max = 420, value = 42, step = 7),
      sliderInput("loadn", "Loading injections",
                  min = 0, max = 9, value = 6, step = 1),
      sliderInput("maxinj", "Injection cap",
                  min = 0, max = 40, value = 40, step = 1),
      sliderInput("cstthr", "Retreatment threshold, CST (µm)", min = 260, max = 450,
                  value = 310, step = 10),
      checkboxInput("dexadd", "Add dexamethasone implant (3·9 months)", FALSE),
      checkboxInput("prp", "Panretinal photocoagulation (PRP, 3 months)", FALSE),

      h4("4 · Horizon"),
      sliderInput("end", "Observation period (months)", min = 6, max = 36,
                  value = 36, step = 1),
      actionButton("go", "Simulate", class = "btn-primary",
                   width = "100%")
    ),

    mainPanel(
      width = 9,
      htmlOutput("verdict"),
      tabsetPanel(
        id = "tabs",

        # ------------------------------------------------------------------
        tabPanel(
          "① Patient Profile",
          h4("The pressure arm at a glance"),
          tableOutput("profile"),
          plotOutput("profile_plot", height = "330px"),
          helpText(HTML(paste0("<b>P<sub>c</sub>* = P<sub>t</sub> + ",
            "&sigma;<sub>max</sub>(&pi;<sub>c</sub> &minus; &pi;<sub>t0</sub>)</b>. ",
            "Above this capillary pressure the Starling bracket stays positive ",
            "even with the barrier fully restored, so the macula cannot be dried ",
            "by any anti-VEGF dose. The critical venous resistance R<sub>v</sub>* ",
            "diverges as P<sub>a</sub> falls towards P<sub>c</sub>*: if arterial ",
            "inlet pressure is below P<sub>c</sub>*, no occlusion can create a ",
            "permanent floor.")))
        ),

        # ------------------------------------------------------------------
        tabPanel(
          "② Pressure Arm vs Permeability Arm",
          h4("Decomposing the two arms"),
          plotOutput("arms_plot", height = "420px"),
          plotOutput("bracket_plot", height = "300px"),
          helpText(HTML(paste0("Top: P<sub>c</sub> against the intact-barrier ",
            "P<sub>c</sub>*. A crossing means the eye becomes dryable. ",
            "Bottom: the Starling bracket, the permeability multiplier ",
            "L<sub>p,rel</sub>, and the irreducible flux J<sub>v,floor</sub> ",
            "obtained by forcing the permeability arm to its floor.")))
        ),

        # ------------------------------------------------------------------
        tabPanel(
          "③ Drug PK",
          h4("Vitreous · Effect Site · Plasma"),
          plotOutput("pk_plot", height = "420px"),
          h4("Suppression-duration decomposition"),
          tableOutput("supp_table"),
          helpText(HTML(paste0("t<sub>sup</sub> = (t<sub>&frac12;</sub>/ln2)&middot;",
            "ln[C<sub>0,eff</sub> / ((R&minus;1)&middot;&alpha;&middot;K<sub>D</sub>)]. ",
            "Duration is <b>logarithmic</b> in molar dose and affinity but ",
            "<b>linear</b> in half-life, which is why a 94-fold affinity gap ",
            "between ranibizumab and aflibercept becomes only a ~4-fold gap in ",
            "time &mdash; and why SCORE2 found bevacizumab non-inferior to ",
            "aflibercept.")))
        ),

        # ------------------------------------------------------------------
        tabPanel(
          "④ Barrier and Fluid",
          h4("Junction · Permeability · Filtration/Drainage Balance"),
          plotOutput("barrier_plot", height = "420px"),
          plotOutput("flux_plot", height = "300px"),
          helpText(HTML(paste0("The clearance term J<sub>out</sub> = ",
            "0.10&middot;W + 32&middot;W/(45+W) &micro;m/d <b>saturates</b>. ",
            "Once J<sub>v</sub> exceeds the pump's ceiling, thickness is limited ",
            "only by the small linear term, which is why oedema is so much more ",
            "sensitive to the bracket at high flux than at low flux.")))
        ),

        # ------------------------------------------------------------------
        tabPanel(
          "⑤ Clinical Endpoints",
          fluidRow(
            column(6, plotOutput("cst_plot", height = "300px")),
            column(6, plotOutput("bcva_plot", height = "300px"))),
          h4("Where the letters went"),
          plotOutput("letters_plot", height = "300px"),
          tableOutput("endpoint_table"),
          helpText(HTML(paste0("L<sub>ed</sub> is a function of the current ",
            "state and comes back when the macula dries. L<sub>pr</sub> is the ",
            "time-integral of thickening &times; ischaemia and does not. ",
            "That asymmetry is the entire argument for treating early.")))
        ),

        # ------------------------------------------------------------------
        tabPanel(
          "⑥ Ischaemia · Neovascularisation",
          plotOutput("isch_plot", height = "420px"),
          tableOutput("nvg_table"),
          helpText(HTML(paste0("Anti-VEGF suppresses the <i>signalling</i> drive ",
            "for neovascularisation without removing the <i>ischaemic</i> drive. ",
            "Stopping while non-perfusion persists therefore produces rebound: ",
            "run the same eye with a low injection cap to see it.")))
        ),

        # ------------------------------------------------------------------
        tabPanel(
          "⑦ Scenario Comparison",
          h4("Compare regimens on the same eye"),
          checkboxGroupInput(
            "compare", NULL, inline = TRUE,
            choices = c("None" = "none",
                        "Ranibizumab PRN" = "rbz_prn",
                        "Aflibercept T&E" = "afl_tae",
                        "Faricimab T&E" = "far_tae",
                        "Bevacizumab PRN" = "bev_prn",
                        "Dexamethasone q6mo" = "dex",
                        "Blood-pressure lowering + aflibercept" = "bp_afl",
                        "Infinite blockade (floor)" = "perfect"),
            selected = c("none", "afl_tae", "perfect")),
          plotOutput("compare_plot", height = "480px"),
          tableOutput("compare_table")
        ),

        # ------------------------------------------------------------------
        tabPanel(
          "⑧ Biomarkers",
          h4("Tone versus delivered signal"),
          plotOutput("biomarker_plot", height = "420px"),
          h4("Systemic free-VEGF suppression"),
          plotOutput("sysvegf_plot", height = "260px"),
          helpText(HTML(paste0("Aqueous VEGF and IL-6 are shown in the units the ",
            "clinical literature reports (pg/mL). Note that a drug moves the ",
            "<i>delivered</i> signal, not the <i>tone</i>: total mediator ",
            "production is unchanged, which is why the tone rebounds the moment ",
            "the drug falls below its apparent IC<sub>50</sub>.")))
        ),

        # ------------------------------------------------------------------
        tabPanel(
          "⑨ The Cost of Delay",
          h4("Permanent letters per month of delay"),
          actionButton("run_delay", "Compute delay curve"),
          plotOutput("delay_plot", height = "380px"),
          tableOutput("delay_table")
        ),

        # ------------------------------------------------------------------
        tabPanel(
          "⑩ Virtual Population",
          h4("Virtual Population"),
          sliderInput("npop", "Number of patients (n)", min = 20, max = 300,
                      value = 60, step = 20),
          actionButton("run_pop", "Simulate population"),
          plotOutput("pop_plot", height = "380px"),
          tableOutput("pop_table"),
          helpText(HTML(paste0("The key population read-out is the split between ",
            "eyes that are wet <i>with</i> a pressure floor (switching agent will ",
            "not help) and eyes that are wet <i>without</i> one (a genuine ",
            "pharmacological miss).")))
        )
      )
    )
  )
)

# =============================================================================
#  SERVER
# =============================================================================
server <- function(input, output, session) {

  # keep phenotype defaults in sync with the sliders
  observeEvent(input$pheno, {
    ph <- PHENOTYPES[[input$pheno]]
    updateSliderInput(session, "occ0", value = ph$OCC0)
    updateSliderInput(session, "occres", value = ph$OCCRES)
    updateSliderInput(session, "collmax", value = ph$COLLMAX)
  })

  extra_params <- reactive({
    list(OCC0 = input$occ0, OCCRES = input$occres, COLLMAX = input$collmax,
         PA = input$pa, IOP0 = input$iop0)
  })

  # ---- the main simulation ---------------------------------------------------
  sim <- eventReactive(input$go, {
    endd <- input$end * 30.4
    dex <- if (isTRUE(input$dexadd)) input$start + c(91, 274) else numeric(0)
    prp <- if (isTRUE(input$prp)) input$start + 91 else NULL
    ex <- extra_params()

    out <- switch(
      input$regimen,
      none = run_case(input$pheno, input$drug, extra = ex, end = endd,
                      dex_times = dex, prp_time = prp),
      perf = run_case(input$pheno, input$drug, ivt_times = input$start,
                      perfect = TRUE, pfrom = input$start, extra = ex, end = endd),
      dex  = run_case(input$pheno, input$drug, extra = ex, end = endd,
                      dex_times = input$start + 182.5 * (0:5), prp_time = prp),
      q4w  = run_case(input$pheno, input$drug, extra = ex, end = endd,
                      ivt_times = q4w(min(input$maxinj,
                                          floor((endd - input$start) / 28) + 1),
                                      input$start, 28),
                      dex_times = dex, prp_time = prp),
      q8w  = run_case(input$pheno, input$drug, extra = ex, end = endd,
                      ivt_times = c(q4w(min(3, input$maxinj), input$start, 28),
                                    q4w(max(0, input$maxinj - 3),
                                        input$start + 84, 56)),
                      dex_times = dex, prp_time = prp),
      tae  = run_reactive(input$pheno, input$drug, mode = "tae",
                          load_n = input$loadn, start = input$start,
                          cst_thr = input$cstthr, max_inj = input$maxinj,
                          dex_times = dex, extra = ex, end = endd),
      prn  = run_reactive(input$pheno, input$drug, mode = "prn",
                          load_n = input$loadn, start = input$start,
                          cst_thr = input$cstthr, max_inj = input$maxinj,
                          dex_times = dex, extra = ex, end = endd))
    out
  }, ignoreNULL = FALSE)

  n_inj <- reactive({
    a <- attr(sim(), "n_inj")
    if (is.null(a)) NA_integer_ else a
  })

  # ---- verdict banner -------------------------------------------------------
  output$verdict <- renderUI({
    d <- sim()
    last <- tail(d, 1)
    pcs <- last$PCSTAR
    floor_present <- last$PCt > pcs
    wet <- last$CSTt > 300
    if (floor_present && wet) {
      col <- "#c62828"
      msg <- paste0(
        "<b>PRESSURE-ARM EYE.</b> P<sub>c</sub> = ",
        sprintf("%.1f", last$PCt), " mmHg &gt; P<sub>c</sub>* = ",
        sprintf("%.1f", pcs), " mmHg, and the macula is wet (CST ",
        sprintf("%.0f", last$CSTt), " µm). ",
        "Switching agent, raising the dose or shortening the interval cannot ",
        "dry this eye &mdash; only the pressure arm can. Try the P<sub>a</sub> ",
        "slider.")
    } else if (floor_present && !wet) {
      col <- "#ef6c00"
      msg <- paste0(
        "<b>BORDERLINE.</b> P<sub>c</sub> = ", sprintf("%.1f", last$PCt),
        " mmHg is above P<sub>c</sub>* = ", sprintf("%.1f", pcs),
        " mmHg, so a floor exists, but the pump is still clearing it (CST ",
        sprintf("%.0f", last$CSTt), " µm). Expect recurrence on withdrawal.")
    } else if (wet) {
      col <- "#1565c0"
      msg <- paste0(
        "<b>PHARMACOLOGICAL MISS.</b> P<sub>c</sub> = ",
        sprintf("%.1f", last$PCt), " mmHg is BELOW P<sub>c</sub>* = ",
        sprintf("%.1f", pcs), " mmHg, so this macula is dryable &mdash; and it ",
        "is wet (CST ", sprintf("%.0f", last$CSTt),
        " µm). More drug, or more often, should work.")
    } else {
      col <- "#2e7d32"
      msg <- paste0("<b>DRY.</b> CST ", sprintf("%.0f", last$CSTt),
                    " µm, BCVA ", sprintf("%.1f", last$BCVA), " letters, ",
                    ifelse(is.na(n_inj()), "no", n_inj()), " injections.")
    }
    HTML(paste0("<div style='padding:9px 12px;margin-bottom:9px;border-left:6px ",
                "solid ", col, ";background:#fafafa;font-size:14px'>", msg,
                "</div>"))
  })

  # ---- ① patient profile ----------------------------------------------------
  output$profile <- renderTable({
    d <- sim(); last <- tail(d, 1)
    pcs <- critical_pressure()
    rvstar <- critical_rv(pcs, input$pa)
    tibble(
      Item = c("Capillary pressure, P_c (mmHg)", "Critical pressure, P_c* (mmHg)",
               "Venous resistance, R_v", "Critical venous resistance, R_v*",
               "Pressure transmission coefficient, dP_c/dP_a", "Relative blood flow, Q_rel",
               "Non-perfusion (disc areas)", "Cumulative oedema exposure (µm·d)",
               "Number of injections"),
      Value = c(sprintf("%.1f", last$PCt), sprintf("%.2f", pcs),
             sprintf("%.2f", last$RVREt),
             ifelse(is.na(rvstar), "unreachable", sprintf("%.2f", rvstar)),
             sprintf("%.3f", last$DPCDPA), sprintf("%.3f", last$QRELt),
             sprintf("%.1f", last$NP), sprintf("%.0f", last$CUMED),
             ifelse(is.na(n_inj()), "-", as.character(n_inj()))))
  }, striped = TRUE, spacing = "xs")

  output$profile_plot <- renderPlot({
    pcs <- critical_pressure()
    tibble(Pa = seq(31.4, 50, by = 0.2)) %>%
      mutate(Rv_star = critical_rv(pcs, Pa)) %>%
      ggplot(aes(Pa, Rv_star)) +
      geom_line(linewidth = 1, colour = PAL[1]) +
      geom_vline(xintercept = input$pa, linetype = "dashed", colour = PAL[2]) +
      geom_hline(yintercept = 0.25, linetype = "dotted") +
      scale_y_log10() +
      annotate("text", x = 32.5, y = 0.28, hjust = 0, size = 3.3,
               label = "healthy R_v = 0.25") +
      labs(x = expression("arteriolar inlet pressure P"[a]*" (mmHg)"),
           y = expression("critical venous resistance R"[v]*"*  (log scale)"),
           title = "How hard the vein must be squeezed before a permanent oedema floor appears",
           subtitle = paste0("dashed line = this patient's P_a. The curve diverges as P_a falls to P_c* = ",
                             sprintf("%.1f", pcs), " mmHg")) + THEME
  })

  # ---- ② the two arms -------------------------------------------------------
  output$arms_plot <- renderPlot({
    d <- sim()
    d %>% select(time, PCt, PCSTAR, PCCRIt) %>%
      pivot_longer(-time) %>%
      mutate(name = recode(name,
        PCt = "P_c (capillary pressure)",
        PCSTAR = "P_c* (intact barrier)",
        PCCRIt = "P_c* (current barrier)")) %>%
      ggplot(aes(time / 30.4, value, colour = name)) +
      geom_line(linewidth = 0.9) +
      scale_colour_manual(values = PAL) +
      labs(x = "months from occlusion", y = "mmHg",
           title = "The crossing that decides everything") + THEME
  })

  output$bracket_plot <- renderPlot({
    d <- sim()
    d %>% select(time, BRACKt, LPRELt, JVt, JVFLRt) %>%
      pivot_longer(-time) %>%
      mutate(name = recode(name,
        BRACKt = "bracket (mmHg)", LPRELt = "L_p,rel (-)",
        JVt = "J_v (µm/d)", JVFLRt = "J_v,floor (µm/d)")) %>%
      ggplot(aes(time / 30.4, value, colour = name)) +
      geom_line(linewidth = 0.85) + geom_hline(yintercept = 0, colour = "grey60") +
      facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = PAL, guide = "none") +
      labs(x = "months", y = NULL) + THEME
  })

  # ---- ③ PK -----------------------------------------------------------------
  output$pk_plot <- renderPlot({
    d <- sim()
    d %>% select(time, CRETt, SUPPRt, CPUG, CVITUG) %>%
      pivot_longer(-time) %>%
      mutate(name = recode(name,
        CRETt = "effect-site concentration (pM)",
        SUPPRt = "VEGF suppression ratio (-)",
        CPUG = "plasma concentration (µg/mL)",
        CVITUG = "vitreous concentration (µg/mL)"),
        value = pmax(value, 1e-6)) %>%
      ggplot(aes(time / 30.4, value)) +
      geom_line(linewidth = 0.85, colour = PAL[1]) +
      geom_hline(data = tibble(name = "VEGF suppression ratio (-)", y = 4.9),
                 aes(yintercept = y), linetype = "dashed", colour = PAL[2]) +
      facet_wrap(~name, scales = "free_y") + scale_y_log10() +
      labs(x = "months", y = NULL,
           title = "Drug kinetics on a log scale; dashed line = the adequate-suppression threshold") +
      THEME
  })

  output$supp_table <- renderTable({
    bind_rows(lapply(DRUGS$drug, function(dd)
      suppression_duration(dd, target_ratio = 4.9))) %>%
      transmute(Agent = label,
                `Molar dose (nmol)` = sprintf("%.1f", nmol),
                `C0 vitreous (µM)` = sprintf("%.2f", C0_vit_uM),
                `K_D (pM)` = sprintf("%.2f", KD_pM),
                `Threshold (pM)` = sprintf("%.0f", threshold_pM),
                `t½ (d)` = sprintf("%.2f", t_half),
                `t_sup (d)` = sprintf("%.1f", t_sup),
                `Reservoir term` = sprintf("%.2f", reservoir_term),
                `Affinity term` = sprintf("%.2f", affinity_term),
                `t½/ln2` = sprintf("%.2f", halflife_mult))
  }, striped = TRUE, spacing = "xs")

  # ---- ④ barrier ------------------------------------------------------------
  output$barrier_plot <- renderPlot({
    d <- sim()
    d %>% select(time, TJ, LPRELt, SIGMAt, CHRON) %>%
      pivot_longer(-time) %>%
      mutate(name = recode(name,
        TJ = "TJ — junction integrity", LPRELt = "L_p,rel",
        SIGMAt = "σ — reflection coefficient",
        CHRON = "CHRON — chronic remodelling")) %>%
      ggplot(aes(time / 30.4, value, colour = name)) +
      geom_line(linewidth = 0.85) + facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = PAL, guide = "none") +
      labs(x = "months", y = NULL,
           title = "The permeability arm, and the ceiling chronicity puts on its recovery") +
      THEME
  })

  output$flux_plot <- renderPlot({
    d <- sim()
    jout <- 0.10 * pmax(d$W, 0) + 32 * pmax(d$W, 0) / (45 + pmax(d$W, 0))
    tibble(time = d$time, `J_v (in)` = d$JVt, `J_out (out)` = jout,
           `net dW/dt` = d$JVt - jout) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time / 30.4, value, colour = name)) +
      geom_line(linewidth = 0.85) + geom_hline(yintercept = 0, colour = "grey60") +
      scale_colour_manual(values = PAL) +
      labs(x = "months", y = "µm/day",
           title = "Filtration against a saturable pump") + THEME
  })

  # ---- ⑤ endpoints ----------------------------------------------------------
  output$cst_plot <- renderPlot({
    d <- sim()
    ggplot(d, aes(time / 30.4, CSTt)) +
      geom_line(linewidth = 0.95, colour = PAL[1]) +
      geom_hline(yintercept = c(250, 300), linetype = c("solid", "dashed"),
                 colour = "grey55") +
      labs(x = "months", y = "central subfield thickness (µm)",
           title = "CST") + THEME
  })

  output$bcva_plot <- renderPlot({
    d <- sim()
    ggplot(d, aes(time / 30.4, BCVA)) +
      geom_line(linewidth = 0.95, colour = PAL[3]) +
      labs(x = "months", y = "ETDRS letters", title = "BCVA") + THEME
  })

  output$letters_plot <- renderPlot({
    d <- sim()
    d %>% select(time, LEDt, LPRt, LISCt) %>%
      pivot_longer(-time) %>%
      mutate(name = recode(name,
        LEDt = "L_ed — reversible (thickening)",
        LPRt = "L_pr — irreversible (photoreceptors)",
        LISCt = "L_isch — macular ischaemia")) %>%
      ggplot(aes(time / 30.4, value, fill = name)) +
      geom_area(alpha = 0.85) + scale_fill_manual(values = PAL[c(1, 2, 4)]) +
      labs(x = "months", y = "letters lost",
           title = "Only the blue band comes back") + THEME
  })

  output$endpoint_table <- renderTable({
    d <- sim()
    pick <- function(day, col) d[[col]][which.min(abs(d$time - day))]
    days <- c(input$start, input$start + 182, input$start + 365, 730,
              max(d$time))
    days <- unique(days[days <= max(d$time)])
    tibble(`Time point (day)` = sprintf("%.0f", days),
           `CST (µm)` = sprintf("%.0f", sapply(days, pick, "CSTt")),
           `BCVA (letters)` = sprintf("%.1f", sapply(days, pick, "BCVA")),
           `Snellen 20/x` = sprintf("%.0f", sapply(days, pick, "SNELL")),
           `EZ` = sprintf("%.3f", sapply(days, pick, "EZ")),
           `P_c` = sprintf("%.1f", sapply(days, pick, "PCt")))
  }, striped = TRUE, spacing = "xs")

  # ---- ⑥ ischaemia ----------------------------------------------------------
  output$isch_plot <- renderPlot({
    d <- sim()
    d %>% select(time, NP, MI, NVI, IOP, QRELt, VACTt) %>%
      pivot_longer(-time) %>%
      mutate(name = recode(name,
        NP = "NP — non-perfusion (disc areas)", MI = "MI — macular ischaemia",
        NVI = "NVI — iris/angle neovascularisation", IOP = "IOP (mmHg)",
        QRELt = "Q_rel — relative flow", VACTt = "delivered VEGF (pM)")) %>%
      ggplot(aes(time / 30.4, value, colour = name)) +
      geom_line(linewidth = 0.85) + facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = PAL, guide = "none") +
      labs(x = "months", y = NULL) + THEME
  })

  output$nvg_table <- renderTable({
    d <- sim(); last <- tail(d, 1)
    tibble(Item = c("Non-perfusion, NP (DA)", "Macular ischaemia, MI", "NVI", "IOP (mmHg)",
                    "Neovascular glaucoma", "Delivered VEGF (pM)"),
           Value = c(sprintf("%.1f", last$NP), sprintf("%.2f", last$MI),
                  sprintf("%.3f", last$NVI), sprintf("%.1f", last$IOP),
                  ifelse(last$NVG > 0.5, "Yes", "No"),
                  sprintf("%.2f", last$VACTt)))
  }, striped = TRUE, spacing = "xs")

  # ---- ⑦ comparison ---------------------------------------------------------
  comparison <- eventReactive(list(input$go, input$compare), {
    endd <- input$end * 30.4
    ex <- extra_params()
    make <- function(k) switch(
      k,
      none     = run_case(input$pheno, "aflibercept", extra = ex, end = endd),
      rbz_prn  = run_reactive(input$pheno, "ranibizumab", mode = "prn",
                              load_n = 6, start = input$start, extra = ex,
                              end = endd),
      afl_tae  = run_reactive(input$pheno, "aflibercept", mode = "tae",
                              load_n = 6, start = input$start, extra = ex,
                              end = endd),
      far_tae  = run_reactive(input$pheno, "faricimab", mode = "tae",
                              load_n = 6, start = input$start, extra = ex,
                              end = endd),
      bev_prn  = run_reactive(input$pheno, "bevacizumab", mode = "prn",
                              load_n = 6, start = input$start, extra = ex,
                              end = endd),
      dex      = run_case(input$pheno, "aflibercept", extra = ex, end = endd,
                          dex_times = input$start + 182.5 * (0:5)),
      bp_afl   = run_reactive(input$pheno, "aflibercept", mode = "tae",
                              load_n = 6, start = input$start,
                              extra = modifyList(ex, list(PA = ex$PA - 8)),
                              end = endd),
      perfect  = run_case(input$pheno, "aflibercept", ivt_times = input$start,
                          perfect = TRUE, pfrom = input$start, extra = ex,
                          end = endd))
    sel <- input$compare
    if (!length(sel)) sel <- "afl_tae"
    setNames(lapply(sel, function(k) {
      o <- make(k); o$arm <- k
      o$n_inj <- if (is.null(attr(o, "n_inj"))) NA else attr(o, "n_inj")
      o
    }), sel)
  })

  output$compare_plot <- renderPlot({
    bind_rows(comparison()) %>%
      select(time, arm, CSTt, BCVA, PCt, VACTt) %>%
      pivot_longer(c(CSTt, BCVA, PCt, VACTt)) %>%
      mutate(name = recode(name, CSTt = "CST (µm)", BCVA = "BCVA (letters)",
                           PCt = "P_c (mmHg)", VACTt = "delivered VEGF (pM)")) %>%
      ggplot(aes(time / 30.4, value, colour = arm)) +
      geom_line(linewidth = 0.85) + facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = PAL) +
      labs(x = "months", y = NULL) + THEME
  })

  output$compare_table <- renderTable({
    bind_rows(lapply(comparison(), function(d) {
      pick <- function(day, col) d[[col]][which.min(abs(d$time - day))]
      b <- pick(input$start, "BCVA")
      tibble(Regimen = d$arm[1],
             `Injections (inj)` = ifelse(is.na(d$n_inj[1]), "-",
                                   as.character(d$n_inj[1])),
             `CST start` = sprintf("%.0f", pick(input$start, "CSTt")),
             `CST final` = sprintf("%.0f", tail(d$CSTt, 1)),
             `BCVA start` = sprintf("%.1f", b),
             `BCVA final` = sprintf("%.1f", tail(d$BCVA, 1)),
             `Change` = sprintf("%+.1f", tail(d$BCVA, 1) - b),
             `Days dry` = sum(d$CSTt < 300),
             `Days suppressed` = sprintf("%.0f", tail(d$TSUP, 1)))
    }))
  }, striped = TRUE, spacing = "xs")

  # ---- ⑧ biomarkers ---------------------------------------------------------
  output$biomarker_plot <- renderPlot({
    d <- sim()
    tibble(time = d$time,
           `VEGF tone (pg/mL)` = d$VTONE / 0.0222,
           `VEGF delivered (pg/mL)` = d$VACTt / 0.0222,
           `IL-6 (pg/mL)` = d$IL6,
           `Ang-2 delivered (pM)` = d$AACTt,
           `Ang-2 tone (pM)` = d$ATONE,
           `HIF-1α (0-1)` = d$HIF) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time / 30.4, pmax(value, 1e-3), colour = name)) +
      geom_line(linewidth = 0.85) + scale_y_log10() +
      scale_colour_manual(values = PAL) +
      labs(x = "months", y = "log scale",
           title = "A drug moves the delivered signal, not the tone") + THEME
  })

  output$sysvegf_plot <- renderPlot({
    d <- sim()
    dr <- DRUGS[DRUGS$drug == input$drug, ]
    supp <- 1 - 1 / (1 + d$CPUG / dr$mw * 1e6 * 1e3 / dr$kdv)
    tibble(time = d$time, `systemic free-VEGF suppression` = supp) %>%
      ggplot(aes(time / 30.4, 100 * `systemic free-VEGF suppression`)) +
      geom_line(linewidth = 0.9, colour = PAL[5]) +
      labs(x = "months", y = "% suppression of plasma free VEGF",
           title = paste0(dr$label,
                          " — the systemic cost of an intraocular injection")) +
      THEME
  })

  # ---- ⑨ delay curve --------------------------------------------------------
  delay <- eventReactive(input$run_delay, {
    ex <- extra_params()
    bind_rows(lapply(c(0, 0.5, 1, 2, 3, 4, 6, 9, 12), function(mo) {
      d <- run_reactive(input$pheno, input$drug, mode = "tae", load_n = 6,
                        start = 42 + mo * 30.4, extra = ex, end = 1095)
      late <- d$BCVA[d$time >= 730]
      tibble(delay_months = mo,
             bcva_start = d$BCVA[which.min(abs(d$time - (42 + mo * 30.4)))],
             peak = max(d$BCVA), late_mean = mean(late),
             ez = tail(d$EZ, 1), n_inj = attr(d, "n_inj"))
    })) %>% mutate(permanent_loss = late_mean - late_mean[1])
  })

  output$delay_plot <- renderPlot({
    delay() %>% select(delay_months, late_mean, ez, permanent_loss) %>%
      pivot_longer(-delay_months) %>%
      mutate(name = recode(name,
        late_mean = "mean BCVA, months 24-36 (letters)",
        ez = "ellipsoid zone integrity at month 36",
        permanent_loss = "permanent letters forfeited")) %>%
      ggplot(aes(delay_months, value, colour = name)) +
      geom_line(linewidth = 0.9) + geom_point(size = 2) +
      facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = PAL, guide = "none") +
      labs(x = "months of delay before the first injection", y = NULL) + THEME
  })

  output$delay_table <- renderTable({
    delay() %>% transmute(`Delay (months)` = sprintf("%.1f", delay_months),
                          `Starting BCVA` = sprintf("%.1f", bcva_start),
                          `Peak BCVA` = sprintf("%.1f", peak),
                          `Mean m24-36` = sprintf("%.1f", late_mean),
                          `EZ m36` = sprintf("%.3f", ez),
                          `Permanent loss` = sprintf("%+.1f", permanent_loss),
                          `Injections` = n_inj)
  }, striped = TRUE, spacing = "xs")

  # ---- ⑩ virtual population -------------------------------------------------
  pop <- eventReactive(input$run_pop, {
    withProgress(message = "simulating population", {
      virtual_population(n = input$npop, end = 730)
    })
  })

  output$pop_plot <- renderPlot({
    p <- pop()
    p %>% mutate(group = case_when(
      floor == 1 & dry == 0 ~ "wet, pressure floor",
      floor == 1 & dry == 1 ~ "dry despite a floor",
      floor == 0 & dry == 0 ~ "wet, no floor (pharmacological miss)",
      TRUE ~ "dry, no floor")) %>%
      ggplot(aes(pc, cst, colour = group, size = n_inj)) +
      geom_point(alpha = 0.75) +
      geom_vline(xintercept = critical_pressure(), linetype = "dashed") +
      geom_hline(yintercept = 300, linetype = "dotted") +
      scale_colour_manual(values = PAL) + scale_size(range = c(1, 4)) +
      labs(x = expression("capillary pressure P"[c]*" at month 24 (mmHg)"),
           y = "residual CST at month 24 (µm)",
           title = "Residual thickness tracks capillary pressure, not injection count",
           subtitle = paste0("dashed line = P_c* = ",
                             sprintf("%.1f", critical_pressure()), " mmHg")) +
      THEME
  })

  output$pop_table <- renderTable({
    summarise_population(pop()) %>%
      pivot_longer(everything(), names_to = "Metric", values_to = "Value") %>%
      mutate(Value = sprintf("%.2f", Value))
  }, striped = TRUE, spacing = "xs")
}

shinyApp(ui, server)
