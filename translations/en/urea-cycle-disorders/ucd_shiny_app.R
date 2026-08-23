## ============================================================================
## Urea Cycle Disorders (UCD) QSP — Shiny Dashboard
## ----------------------------------------------------------------------------
## 8 tabs:
##   1 Patient profile   — build the virtual UCD patient (locus, residual
##                         activity, diet) and read the analytically derived
##                         untreated baseline BEFORE any drug is given
##   2 Drug PK           — phenylbutyrate -> phenylacetate -> PAGN and the
##                         benzoate -> hippurate axis; NaPBA vs GPB profiles
##   3 Ammonia & N flux  — plasma ammonia, glutamine and the complete
##                         nitrogen-disposal ledger (urea / PAGN / hippurate /
##                         renal), in umol N/h and as % of the daily load
##   4 Brain             — brain ammonia, the astrocytic glutamine trap,
##                         myo-inositol depletion, brain water, ICP, HE stage
##   5 Clinical endpoints— coma hours, cumulative injury, modelled IQ and the
##                         natural-protein tolerance that patients actually feel
##   6 Scenario compare  — all 14 prebuilt arms side by side + summary table
##   7 Biomarkers        — U-PAGN, PAA:PAGN ratio, orotate, citrulline,
##                         arginine, glycine, BCAA: the monitoring panel
##   8 Safety            — phenylacetate neurotoxicity, sodium load, BCAA and
##                         glycine depletion, dialysis rebound
##
## Dependencies: shiny, shinydashboard, mrgsolve, dplyr, tidyr, ggplot2, DT
## Run with:  shiny::runApp("ucd_shiny_app.R")
##            (keep ucd_mrgsolve_model.R in the same directory)
##
## EDUCATIONAL / RESEARCH USE ONLY — not validated for clinical decisions.
## ----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(shiny)
  library(shinydashboard)
  library(mrgsolve)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(DT)
})

## ---------------------------------------------------------------- model load
UCD_SPEC <- "ucd_mrgsolve_model.R"
if (!file.exists(UCD_SPEC)) {
  stop("ucd_mrgsolve_model.R not found — keep it in the same directory as the app.")
}
source(UCD_SPEC, local = FALSE)
MOD <- UCD_build()

THEME <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.position  = "bottom",
        legend.title     = element_blank(),
        plot.title       = element_text(face = "bold", size = 12))

PAL <- c("#1b6ca8", "#c0392b", "#27865c", "#b8860b", "#7d3c98",
         "#16a085", "#d35400", "#2c3e50")

## band helper — draws a reference band behind a line plot
band <- function(lo, hi, fill = "#2ecc71", alpha = 0.10) {
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = lo, ymax = hi,
           fill = fill, alpha = alpha)
}

## ============================================================================
## UI
## ============================================================================
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "UCD QSP Dashboard", titleWidth = 300),

  dashboardSidebar(
    width = 300,
    sidebarMenu(
      id = "tabs",
      menuItem("1. Patient Profile",      tabName = "patient",  icon = icon("user")),
      menuItem("2. Drug PK",            tabName = "pk",       icon = icon("pills")),
      menuItem("3. Ammonia · Nitrogen Balance", tabName = "nitrogen", icon = icon("scale-balanced")),
      menuItem("4. Brain (Astrocytes)",    tabName = "brain",    icon = icon("brain")),
      menuItem("5. Clinical Endpoints",    tabName = "endpoint", icon = icon("stethoscope")),
      menuItem("6. Scenario Comparison",      tabName = "compare",  icon = icon("layer-group")),
      menuItem("7. Biomarkers",         tabName = "biomark",  icon = icon("vial")),
      menuItem("8. Safety",             tabName = "safety",   icon = icon("triangle-exclamation"))
    ),
    hr(),
    h4("Patient", style = "padding-left:15px;"),
    selectInput("locus", "Enzyme block",
                choices = c("OTC (X-linked)"          = "OTC",
                            "CPS1"                    = "CPS1",
                            "NAGS"                    = "NAGS",
                            "ASS1 / ASL (distal arm)" = "DISTAL"),
                selected = "OTC"),
    sliderInput("act",  "Residual enzyme activity (%)",
                min = 0, max = 100, value = 32, step = 1),
    sliderInput("wt",   "Weight (kg)", min = 3, max = 100, value = 70, step = 0.5),
    sliderInput("vbr",  "Brain volume (L)", min = 0.35, max = 1.5,
                value = 1.35, step = 0.05),
    sliderInput("prot", "Natural protein intake (g/kg/day)",
                min = 0, max = 2.5, value = 0.60, step = 0.05),

    hr(),
    h4("Therapy", style = "padding-left:15px;"),
    sliderInput("gpb",  "Glycerol phenylbutyrate GPB (mL/day)",
                min = 0, max = 25, value = 0, step = 0.5),
    sliderInput("napba","Sodium phenylbutyrate NaPBA (g/day)",
                min = 0, max = 30, value = 0, step = 0.5),
    sliderInput("bz",   "Sodium benzoate Na-benzoate (mg/kg/day)",
                min = 0, max = 500, value = 0, step = 25),
    sliderInput("cit",  "L-citrulline (mg/kg/day)",
                min = 0, max = 250, value = 0, step = 10),
    sliderInput("ncg",  "Carglumic acid (mg/kg/day)",
                min = 0, max = 250, value = 0, step = 10),

    hr(),
    h4("Stress / Rescue", style = "padding-left:15px;"),
    sliderInput("cat",  "Catabolic amplitude",
                min = 0, max = 4, value = 0, step = 0.1),
    sliderInput("catt", "Catabolic stress onset (h)", min = 0, max = 168, value = 24, step = 2),
    sliderInput("catd", "Catabolic stress duration (h)", min = 6, max = 168, value = 48, step = 6),
    sliderInput("anab", "Anabolic rescue (0-1)",
                min = 0, max = 1, value = 0, step = 0.1),
    sliderInput("hd",   "Dialysis clearance CL (L/h)",
                min = 0, max = 25, value = 0, step = 0.5),
    sliderInput("hdt",  "Dialysis start (h)", min = 0, max = 168, value = 36, step = 2),
    sliderInput("hdd",  "Dialysis duration (h)", min = 1, max = 48, value = 6, step = 1),

    hr(),
    sliderInput("days", "Simulation duration Horizon (days)",
                min = 2, max = 180, value = 14, step = 1),
    actionButton("go", "Run simulation", icon = icon("play"),
                 class = "btn-primary", style = "margin:10px 15px;width:80%;")
  ),

  dashboardBody(
    tags$head(tags$style(HTML(".content-wrapper{background-color:#f7f9fb;}"))),
    tabItems(

      ## ---------------------------------------------------------- 1 patient
      tabItem(
        "patient",
        fluidRow(
          box(width = 12, title = "What This Shows",
              status = "primary", solidHeader = TRUE,
              HTML("The model analytically computes the pre-treatment steady state from the patient
                    descriptors (enzyme block · residual activity · weight · protein intake), and uses
                    <b>this computed steady state</b> as the initial value for every compartment.
                    So if no drug is given, the curve is completely flat, and every change seen in
                    the simulation <b>comes from the mechanism, not an error in the initial value.</b><br>
                    The untreated steady state is solved from the patient
                    descriptors and used as the initial condition, so an
                    untreated arm is exactly stationary."))),
        fluidRow(
          valueBoxOutput("vb_nh3",  width = 3),
          valueBoxOutput("vb_gln",  width = 3),
          valueBoxOutput("vb_tol",  width = 3),
          valueBoxOutput("vb_cap",  width = 3)),
        fluidRow(
          box(width = 6, title = "Residual Activity vs. Baseline Ammonia (Activity–Ammonia Cliff)",
              status = "info", solidHeader = TRUE, plotOutput("p_cliff", height = 330)),
          box(width = 6, title = "Protein Intake vs. Baseline Ammonia (Protein–Ammonia Cliff)",
              status = "info", solidHeader = TRUE, plotOutput("p_protcliff", height = 330))),
        fluidRow(
          box(width = 12, title = "Derived Baseline State",
              status = "info", solidHeader = TRUE, DTOutput("t_base")))),

      ## --------------------------------------------------------------- 2 PK
      tabItem(
        "pk",
        fluidRow(
          box(width = 6, title = "Phenylbutyrate / Phenylacetate (PBA · PAA)",
              status = "primary", solidHeader = TRUE, plotOutput("p_pba", height = 320)),
          box(width = 6, title = "Phenylacetylglutamine (PAGN)",
              status = "primary", solidHeader = TRUE, plotOutput("p_pagn", height = 320))),
        fluidRow(
          box(width = 6, title = "Benzoate · Hippurate",
              status = "info", solidHeader = TRUE, plotOutput("p_bz", height = 300)),
          box(width = 6, title = "Urinary Excretion (PAGN & Hippurate)",
              status = "info", solidHeader = TRUE, plotOutput("p_uexc", height = 300))),
        fluidRow(
          box(width = 12, title = "NaPBA vs GPB — Same PBA Molar Dose, Different Profile",
              status = "warning", solidHeader = TRUE,
              plotOutput("p_pbavsgpb", height = 320),
              HTML("<small>Even when the same PBA molar dose is given (NaPBA 20 g/day ≈ GPB 17.5 mL/day ≈ 107-109 mmol/day),
                    GPB has a <b>lower Cmax and a flatter 24-hour
                    profile</b> because of its pancreatic-lipase-dependent sustained release — the result of the crossover design in Diaz 2013 Hepatology 57:2171.</small>")))),

      ## --------------------------------------------------------- 3 nitrogen
      tabItem(
        "nitrogen",
        fluidRow(
          box(width = 8, title = "Plasma Ammonia",
              status = "primary", solidHeader = TRUE, plotOutput("p_nh3", height = 340)),
          box(width = 4, title = "Nitrogen Disposal Share",
              status = "primary", solidHeader = TRUE, plotOutput("p_nshare", height = 340))),
        fluidRow(
          box(width = 6, title = "Plasma Glutamine — an Early-Warning Indicator",
              status = "info", solidHeader = TRUE, plotOutput("p_gln", height = 300),
              HTML("<small>Glutamine rises <b>before</b> ammonia does. Exceeding 1000 µmol/L is
                    treated clinically as a crisis warning.</small>")),
          box(width = 6, title = "Deep Tissue Glutamine Reservoir",
              status = "info", solidHeader = TRUE, plotOutput("p_glndeep", height = 300),
              HTML("<small>The source of the post-dialysis rebound.</small>"))),
        fluidRow(
          box(width = 12, title = "Nitrogen Balance (Nitrogen Ledger, µmol N/h)",
              status = "info", solidHeader = TRUE, plotOutput("p_nledger", height = 320)))),

      ## ------------------------------------------------------------ 4 brain
      tabItem(
        "brain",
        fluidRow(
          box(width = 12, title = "Osmotic Hypothesis (the Osmotic Cascade)",
              status = "primary", solidHeader = TRUE,
              HTML("Blood NH<sub>3</sub> → diffusion across the BBB → <b>astrocyte glutamine synthetase</b> →
                    glutamine accumulation → acts as an osmolyte → cell swelling. The defence mechanism is <b>myo-inositol efflux</b>, but
                    <b>it is already depleted in chronic hyperammonaemia</b>, so there is no buffer left for an acute rise.
                    This interaction is the single most important counter-intuitive behaviour in UCD."))),
        fluidRow(
          box(width = 6, title = "Brain Ammonia",
              status = "info", solidHeader = TRUE, plotOutput("p_nh3b", height = 300)),
          box(width = 6, title = "Brain Glutamine (MRS 2.05-2.45 ppm)",
              status = "info", solidHeader = TRUE, plotOutput("p_glnb", height = 300))),
        fluidRow(
          box(width = 6, title = "Myo-inositol (Osmolyte Reserve Capacity)",
              status = "warning", solidHeader = TRUE, plotOutput("p_mins", height = 300)),
          box(width = 6, title = "Brain Oedema & Intracranial Pressure (Brain Water · ICP)",
              status = "danger", solidHeader = TRUE, plotOutput("p_bw", height = 300))),
        fluidRow(
          box(width = 12, title = "Encephalopathy Stage (0-4)",
              status = "danger", solidHeader = TRUE, plotOutput("p_he", height = 300)))),

      ## -------------------------------------------------------- 5 endpoints
      tabItem(
        "endpoint",
        fluidRow(
          valueBoxOutput("vb_coma", width = 3),
          valueBoxOutput("vb_iq",   width = 3),
          valueBoxOutput("vb_peak", width = 3),
          valueBoxOutput("vb_auc",  width = 3)),
        fluidRow(
          box(width = 6, title = "Cumulative Coma-Range Hours",
              status = "danger", solidHeader = TRUE, plotOutput("p_coma", height = 300)),
          box(width = 6, title = "Cumulative Neurological Injury Index & Predicted IQ",
              status = "danger", solidHeader = TRUE, plotOutput("p_iq", height = 300))),
        fluidRow(
          box(width = 6, title = "Ammonia AUC>100 (µmol/L·h)",
              status = "warning", solidHeader = TRUE, plotOutput("p_auc", height = 300)),
          box(width = 6, title = "Natural-Protein Tolerance",
              status = "success", solidHeader = TRUE, plotOutput("p_tol", height = 300),
              HTML("<small>The treatment effect the patient actually experiences —
                    ‘what and how much they can eat’.</small>")))),

      ## --------------------------------------------------------- 6 compare
      tabItem(
        "compare",
        fluidRow(
          box(width = 12, title = "Predefined Scenarios (14 Prebuilt Arms)",
              status = "primary", solidHeader = TRUE,
              checkboxGroupInput("scen", NULL, inline = TRUE,
                                 choices  = setNames(1:14, paste0("S", 1:14)),
                                 selected = c(1, 3, 4, 5, 6, 12)),
              actionButton("goscen", "Run scenarios",
                           icon = icon("play"), class = "btn-primary"))),
        fluidRow(
          box(width = 6, title = "Ammonia", status = "info", solidHeader = TRUE,
              plotOutput("c_nh3", height = 330)),
          box(width = 6, title = "Glutamine", status = "info", solidHeader = TRUE,
              plotOutput("c_gln", height = 330))),
        fluidRow(
          box(width = 6, title = "Protein Tolerance", status = "info", solidHeader = TRUE,
              plotOutput("c_tol", height = 330)),
          box(width = 6, title = "Brain Glutamine", status = "info", solidHeader = TRUE,
              plotOutput("c_glnb", height = 330))),
        fluidRow(
          box(width = 12, title = "Scenario Summary",
              status = "primary", solidHeader = TRUE, DTOutput("c_tab")))),

      ## -------------------------------------------------------- 7 biomarker
      tabItem(
        "biomark",
        fluidRow(
          box(width = 6, title = "Urinary PAGN — Dose Adequacy · Adherence",
              status = "primary", solidHeader = TRUE, plotOutput("b_upagn", height = 300)),
          box(width = 6, title = "PAA : PAGN Ratio — an Indicator of Conjugation Saturation",
              status = "warning", solidHeader = TRUE, plotOutput("b_ratio", height = 300),
              HTML("<small>Exceeding 2.5 on a µg/mL basis = conjugation saturation, a warning that PAA will rise.</small>"))),
        fluidRow(
          box(width = 4, title = "Orotic Acid", status = "info",
              solidHeader = TRUE, plotOutput("b_orot", height = 280),
              HTML("<small>Low citrulline + <b>high orotic acid</b> = OTC deficiency;
                    low citrulline + normal orotic acid = CPS1/NAGS deficiency.</small>")),
          box(width = 4, title = "Citrulline", status = "info",
              solidHeader = TRUE, plotOutput("b_cit", height = 280)),
          box(width = 4, title = "Arginine", status = "info",
              solidHeader = TRUE, plotOutput("b_arg", height = 280))),
        fluidRow(
          box(width = 6, title = "Glycine — Consumed by Benzoate",
              status = "info", solidHeader = TRUE, plotOutput("b_gly", height = 280)),
          box(width = 6, title = "Branched-Chain Amino Acids (BCAA) — Consumed by Phenylbutyrate",
              status = "info", solidHeader = TRUE, plotOutput("b_bcaa", height = 280)))),

      ## ----------------------------------------------------------- 8 safety
      tabItem(
        "safety",
        fluidRow(
          box(width = 12, title = "Safety Signals",
              status = "danger", solidHeader = TRUE,
              HTML("① <b>Phenylacetate neurotoxicity</b> — drowsiness, nausea, headache above 500 µg/mL.
                    Risk increases with renal failure, hepatic failure, and high doses.<br>
                    ② <b>Sodium load</b> — NaPBA and sodium benzoate each deliver a substantial amount of Na along with the drug.<br>
                    ③ <b>BCAA and glycine depletion</b> — characteristic side effects of phenylbutyrate and benzoate, respectively.<br>
                    ④ <b>Post-dialysis rebound</b> — nitrogen redistributes from the deep glutamine reservoir.<br>
                    ⑤ <b>Valproic acid</b> inhibits NAGS/CPS1, so it is close to contraindicated in UCD."))),
        fluidRow(
          box(width = 6, title = "Plasma PAA vs. the 500 µg/mL Toxicity Threshold",
              status = "danger", solidHeader = TRUE, plotOutput("s_paa", height = 320)),
          box(width = 6, title = "Cumulative Excess PAA Exposure (µg/mL·h)",
              status = "danger", solidHeader = TRUE, plotOutput("s_paatox", height = 320))),
        fluidRow(
          box(width = 6, title = "Estimated Sodium Load (mmol/day)",
              status = "warning", solidHeader = TRUE, plotOutput("s_na", height = 300)),
          box(width = 6, title = "Post-Dialysis Rebound",
              status = "warning", solidHeader = TRUE, plotOutput("s_rebound", height = 300))),
        fluidRow(
          box(width = 12, title = "Disclaimer", status = "danger",
              solidHeader = TRUE,
              HTML("<b>This dashboard is a semi-quantitative QSP model for educational and research purposes.</b>
                    It must not be used for clinical decision-making, prescribing, or regulatory submission."))))
    )
  )
)

## ============================================================================
## Server
## ============================================================================
server <- function(input, output, session) {

  ## ---------------------------------------------------- parameter assembly
  pars <- reactive({
    a <- input$act/100
    p <- list(WT = input$wt, VBR = input$vbr, PROT = input$prot,
              OTCACT = 1, CPS1A = 1, NAGSA = 1, DISTALA = 1,
              CATAMP = input$cat, CATT0 = input$catt, CATDUR = input$catd,
              ANAB = input$anab,
              CLHD = input$hd, HDT0 = ifelse(input$hd > 0, input$hdt, 1e9),
              HDDUR = input$hdd,
              BSA = round(0.024265 * (input$wt^0.5378) * (170^0.3964), 3))
    switch(input$locus,
           OTC    = { p$OTCACT  <- a },
           CPS1   = { p$CPS1A   <- a },
           NAGS   = { p$NAGSA   <- a },
           DISTAL = { p$DISTALA <- a })
    ## citrulline is diagnostic: low in proximal blocks, high in distal ones
    p$CITBASE <- if (input$locus %in% c("OTC", "CPS1", "NAGS")) 12 else 180
    p
  })

  evts <- reactive({
    d  <- input$days
    e  <- NULL
    add <- function(x) if (is.null(e)) e <<- x else e <<- c(e, x)
    if (input$gpb   > 0) add(ucd_gpb(input$gpb,  n = 3, ii = 8, addl = ucd_addl(d, 8)))
    if (input$napba > 0) add(ucd_napba(input$napba, n = 3, ii = 8, addl = ucd_addl(d, 8)))
    if (input$bz    > 0) add(ucd_benzoate(input$bz, input$wt, n = 3, ii = 8,
                                          addl = ucd_addl(d, 8)))
    if (input$cit   > 0) add(ucd_citrulline(input$cit, input$wt, n = 3, ii = 8,
                                            addl = ucd_addl(d, 8)))
    if (input$ncg   > 0) add(ucd_carglumic(input$ncg, input$wt, n = 2, ii = 12,
                                           addl = ucd_addl(d, 12)))
    e
  })

  sim <- eventReactive(input$go, {
    m <- do.call(mrgsolve::param, c(list(MOD), pars()))
    dl <- max(0.1, input$days*24/1500)
    as.data.frame(mrgsolve::mrgsim(m, events = evts(), end = input$days*24,
                                   delta = dl, hmax = 0.25))
  }, ignoreNULL = FALSE)

  base <- reactive({
    m <- do.call(mrgsolve::param, c(list(MOD), pars()))
    as.data.frame(mrgsolve::mrgsim(m, end = 1, delta = 1))[1, ]
  })

  tt <- function(d) d$time/24

  ## -------------------------------------------------------------- 1 patient
  output$vb_nh3 <- renderValueBox({
    v <- round(base()$NH3, 1)
    valueBox(paste0(v, " µmol/L"), "Baseline plasma ammonia (NH3)",
             icon = icon("droplet"),
             color = if (v < 50) "green" else if (v < 100) "yellow" else "red")
  })
  output$vb_gln <- renderValueBox({
    v <- round(base()$GLNP)
    valueBox(paste0(v, " µmol/L"), "Baseline plasma glutamine (Gln)",
             icon = icon("flask"),
             color = if (v < 800) "green" else if (v < 1000) "yellow" else "red")
  })
  output$vb_tol <- renderValueBox({
    v <- round(base()$PROTTOL, 2)
    valueBox(paste0(v, " g/kg/d"), "Natural protein tolerance",
             icon = icon("drumstick-bite"),
             color = if (v > 1.0) "green" else if (v > 0.6) "yellow" else "red")
  })
  output$vb_cap <- renderValueBox({
    v <- round(100*base()$ACTEFF, 1)
    valueBox(paste0(v, " %"), "Effective urea cycle activity",
             icon = icon("gauge-high"),
             color = if (v > 50) "green" else if (v > 20) "yellow" else "red")
  })

  output$p_cliff <- renderPlot({
    p <- pars()
    grid <- seq(0.01, 1, by = 0.01)
    y <- vapply(grid, function(a) {
      q <- p
      switch(input$locus, OTC = {q$OTCACT <- a}, CPS1 = {q$CPS1A <- a},
             NAGS = {q$NAGSA <- a}, DISTAL = {q$DISTALA <- a})
      m <- do.call(mrgsolve::param, c(list(MOD), q))
      as.data.frame(mrgsolve::mrgsim(m, end = 1, delta = 1))$NH3[1]
    }, numeric(1))
    ggplot(data.frame(a = 100*grid, y = y), aes(a, y)) +
      band(0, 35) + geom_line(colour = PAL[1], linewidth = 1.1) +
      geom_vline(xintercept = input$act, linetype = 2, colour = PAL[2]) +
      geom_hline(yintercept = 100, linetype = 3) +
      scale_y_log10() +
      labs(x = "Residual enzyme activity (%)", y = "Baseline ammonia (µmol/L, log)",
           title = "Why This Disease Is a Cliff, Not a Slope") + THEME
  })

  output$p_protcliff <- renderPlot({
    p <- pars()
    grid <- seq(0.1, 2.5, by = 0.05)
    y <- vapply(grid, function(g) {
      q <- p; q$PROT <- g
      m <- do.call(mrgsolve::param, c(list(MOD), q))
      as.data.frame(mrgsolve::mrgsim(m, end = 1, delta = 1))$NH3[1]
    }, numeric(1))
    ggplot(data.frame(g = grid, y = y), aes(g, y)) +
      band(0, 35) + geom_line(colour = PAL[3], linewidth = 1.1) +
      geom_vline(xintercept = input$prot, linetype = 2, colour = PAL[2]) +
      geom_hline(yintercept = 100, linetype = 3) +
      scale_y_log10() +
      labs(x = "Natural protein intake (g/kg/day)", y = "Baseline ammonia (µmol/L, log)",
           title = "Why Protein Prescribing Is So Difficult") + THEME
  })

  output$t_base <- renderDT({
    b <- base()
    out <- data.frame(
      item = c("Plasma ammonia (µmol/L)", "Plasma glutamine (µmol/L)",
               "Brain ammonia (µmol/L)", "Brain glutamine (mmol/L)",
               "Brain myo-inositol (mmol/L)", "Excess brain water (%)",
               "Intracranial pressure (mmHg)", "Plasma citrulline (µmol/L)",
               "Plasma arginine (µmol/L)", "Urea nitrogen disposal (µmol N/h)",
               "Total nitrogen disposal (µmol N/h)", "Natural protein tolerance (g/kg/d)"),
      value = round(c(b$NH3, b$GLNP, b$NH3BR, b$GLNBRAIN, b$MINSMM, b$BW, b$ICP,
                      b$CITP, b$ARGP, b$NUREA, b$NTOTAL, b$PROTTOL), 2),
      stringsAsFactors = FALSE)
    colnames(out) <- c("Item", "Value")
    out
  }, options = list(dom = "t", pageLength = 12), rownames = FALSE)

  ## ------------------------------------------------------------------ 2 PK
  output$p_pba <- renderPlot({
    d <- sim()
    dd <- d |> select(time, PBAUG, PAAUG) |>
      pivot_longer(-time, names_to = "k", values_to = "v")
    ggplot(dd, aes(time/24, v, colour = k)) + geom_line(linewidth = 0.9) +
      scale_colour_manual(values = PAL[c(1, 2)],
                          labels = c("PAA (active form)", "PBA (precursor)")) +
      labs(x = "Days", y = "µg/mL", title = "PBA → PAA") + THEME
  })
  output$p_pagn <- renderPlot({
    ggplot(sim(), aes(time/24, PGNUG)) +
      geom_line(colour = PAL[3], linewidth = 0.9) +
      labs(x = "Days", y = "PAGN (µg/mL)",
           title = "Phenylacetylglutamine — Carries Out Two Nitrogens") + THEME
  })
  output$p_bz <- renderPlot({
    dd <- sim() |> select(time, BZUG, HIPUG) |>
      pivot_longer(-time, names_to = "k", values_to = "v")
    ggplot(dd, aes(time/24, v, colour = k)) + geom_line(linewidth = 0.9) +
      scale_colour_manual(values = PAL[c(4, 5)],
                          labels = c("Benzoate", "Hippurate")) +
      labs(x = "Days", y = "µg/mL", title = "Benzoate → Hippurate (One Nitrogen)") + THEME
  })
  output$p_uexc <- renderPlot({
    dd <- sim() |> select(time, UPAGN, UHIPP) |>
      pivot_longer(-time, names_to = "k", values_to = "v")
    ggplot(dd, aes(time/24, v, colour = k)) + geom_line(linewidth = 0.9) +
      scale_colour_manual(values = PAL[c(3, 5)],
                          labels = c("Urinary PAGN", "Urinary hippurate")) +
      labs(x = "Days", y = "mmol/day", title = "Urinary Excretion Rate") + THEME
  })
  output$p_pbavsgpb <- renderPlot({
    m <- do.call(mrgsolve::param, c(list(MOD), pars()))
    a <- as.data.frame(mrgsolve::mrgsim(
      m, events = ucd_napba(20, n = 3, ii = 8, addl = 8), end = 72,
      delta = 0.1, hmax = 0.25)); a$arm <- "NaPBA 20 g/day TID"
    b <- as.data.frame(mrgsolve::mrgsim(
      m, events = ucd_gpb(17.5, n = 3, ii = 8, addl = 8), end = 72,
      delta = 0.1, hmax = 0.25)); b$arm <- "GPB 17.5 mL/day TID"
    ggplot(rbind(a, b), aes(time, PAAUG, colour = arm)) +
      geom_line(linewidth = 0.9) +
      scale_colour_manual(values = PAL[c(2, 1)]) +
      labs(x = "Hours", y = "Plasma PAA (µg/mL)",
           title = "Same PBA Molar Dose, Different Cmax") + THEME
  })

  ## ------------------------------------------------------------ 3 nitrogen
  output$p_nh3 <- renderPlot({
    ggplot(sim(), aes(time/24, NH3)) +
      band(0, 35) +
      geom_hline(yintercept = c(100, 200, 360), linetype = c(3, 2, 1),
                 colour = c("#b8860b", "#d35400", "#c0392b")) +
      geom_line(colour = PAL[1], linewidth = 1.0) +
      labs(x = "Days", y = "Plasma ammonia (µmol/L)",
           title = "100 = Treat, 200 = Coma Risk, 360 = Consider Dialysis") + THEME
  })
  output$p_nshare <- renderPlot({
    d <- sim(); k <- d[nrow(d), ]
    df <- data.frame(
      route = c("Urea", "PAGN", "Hippurate", "Kidney", "Citrulline pathway"),
      v = c(k$NUREA, k$NPAGN, k$NHIPP, k$NRENL, k$NCITR))
    df$route <- factor(df$route, levels = df$route)
    ggplot(df, aes(route, v, fill = route)) +
      geom_col(width = 0.7) + scale_fill_manual(values = PAL[1:5]) +
      labs(x = NULL, y = "µmol N/h", title = "Nitrogen Disposal Routes at the Final Time Point") +
      THEME + theme(legend.position = "none",
                    axis.text.x = element_text(angle = 20, hjust = 1))
  })
  output$p_gln <- renderPlot({
    ggplot(sim(), aes(time/24, GLNP)) +
      band(450, 750) + geom_hline(yintercept = 1000, linetype = 2,
                                  colour = PAL[2]) +
      geom_line(colour = PAL[3], linewidth = 1.0) +
      labs(x = "Days", y = "Plasma glutamine (µmol/L)", title = NULL) + THEME
  })
  output$p_glndeep <- renderPlot({
    dd <- sim() |> select(time, GLNP, GLNDEEP) |>
      pivot_longer(-time, names_to = "k", values_to = "v")
    ggplot(dd, aes(time/24, v, colour = k)) + geom_line(linewidth = 0.9) +
      scale_colour_manual(values = PAL[c(6, 3)],
                          labels = c("Deep reservoir", "Plasma")) +
      labs(x = "Days", y = "µmol/L", title = NULL) + THEME
  })
  output$p_nledger <- renderPlot({
    dd <- sim() |> select(time, NUREA, NPAGN, NHIPP, NRENL, NCITR) |>
      pivot_longer(-time, names_to = "k", values_to = "v")
    dd$k <- factor(dd$k, levels = c("NUREA", "NPAGN", "NHIPP", "NCITR", "NRENL"))
    ggplot(dd, aes(time/24, v, fill = k)) +
      geom_area(alpha = 0.85) +
      scale_fill_manual(values = PAL[1:5],
                        labels = c("Urea", "PAGN", "Hippurate", "Citrulline", "Kidney")) +
      labs(x = "Days", y = "µmol N/h", title = "Nitrogen Disposal Routes Over Time") + THEME
  })

  ## --------------------------------------------------------------- 4 brain
  output$p_nh3b <- renderPlot({
    ggplot(sim(), aes(time/24, NH3BR)) +
      geom_line(colour = PAL[1], linewidth = 1.0) +
      labs(x = "Days", y = "Brain ammonia (µmol/L)", title = NULL) + THEME
  })
  output$p_glnb <- renderPlot({
    ggplot(sim(), aes(time/24, GLNBRAIN)) +
      band(4.0, 6.0) + geom_line(colour = PAL[2], linewidth = 1.0) +
      labs(x = "Days", y = "Brain glutamine (mmol/L)", title = NULL) + THEME
  })
  output$p_mins <- renderPlot({
    ggplot(sim(), aes(time/24, MINSMM)) +
      band(4.0, 6.0) + geom_line(colour = PAL[4], linewidth = 1.0) +
      labs(x = "Days", y = "Brain myo-inositol (mmol/L)",
           title = "Exhausted Buffer Capacity = Why the Next Crisis Is Fatal") + THEME
  })
  output$p_bw <- renderPlot({
    d <- sim()
    ggplot(d, aes(time/24)) +
      geom_line(aes(y = BW), colour = PAL[2], linewidth = 1.0) +
      geom_line(aes(y = ICP/25), colour = PAL[1], linewidth = 0.8, linetype = 2) +
      scale_y_continuous(name = "Excess brain water (%)",
                         sec.axis = sec_axis(~.*25, name = "ICP (mmHg)")) +
      labs(x = "Days", title = "Solid = Brain Water, Dashed = ICP") + THEME
  })
  output$p_he <- renderPlot({
    ggplot(sim(), aes(time/24, HESTAGE)) +
      geom_line(colour = PAL[2], linewidth = 1.1) +
      scale_y_continuous(limits = c(0, 4)) +
      labs(x = "Days", y = "Encephalopathy stage (0-4)", title = NULL) + THEME
  })

  ## ----------------------------------------------------------- 5 endpoints
  output$vb_coma <- renderValueBox({
    v <- round(sim()$COMAH[nrow(sim())], 1)
    valueBox(paste0(v, " h"), "Cumulative time in coma range (>200 µmol/L)", icon = icon("bed"),
             color = if (v < 1) "green" else if (v < 24) "yellow" else "red")
  })
  output$vb_iq <- renderValueBox({
    v <- round(sim()$IQEST[nrow(sim())], 1)
    valueBox(v, "Model-predicted full-scale IQ", icon = icon("brain"),
             color = if (v > 90) "green" else if (v > 75) "yellow" else "red")
  })
  output$vb_peak <- renderValueBox({
    v <- round(max(sim()$NH3), 0)
    valueBox(paste0(v, " µmol/L"), "Peak ammonia", icon = icon("arrow-up"),
             color = if (v < 100) "green" else if (v < 360) "yellow" else "red")
  })
  output$vb_auc <- renderValueBox({
    v <- round(sim()$NAUC[nrow(sim())], 0)
    valueBox(format(v, big.mark = ","), "AUC>100 (µmol/L·h)",
             icon = icon("chart-area"),
             color = if (v < 500) "green" else if (v < 5000) "yellow" else "red")
  })
  output$p_coma <- renderPlot({
    ggplot(sim(), aes(time/24, COMAH)) +
      geom_line(colour = PAL[2], linewidth = 1.0) +
      labs(x = "Days", y = "Cumulative time (h)", title = NULL) + THEME
  })
  output$p_iq <- renderPlot({
    d <- sim()
    ggplot(d, aes(time/24)) +
      geom_line(aes(y = IQEST), colour = PAL[1], linewidth = 1.0) +
      geom_line(aes(y = 100 - NEURO), colour = PAL[2], linewidth = 0.7,
                linetype = 2) +
      labs(x = "Days", y = "IQ (solid) / 100−injury index (dashed)", title = NULL) + THEME
  })
  output$p_auc <- renderPlot({
    ggplot(sim(), aes(time/24, NAUC)) +
      geom_line(colour = PAL[4], linewidth = 1.0) +
      labs(x = "Days", y = "µmol/L·h", title = NULL) + THEME
  })
  output$p_tol <- renderPlot({
    ggplot(sim(), aes(time/24, PROTTOL)) +
      geom_line(colour = PAL[3], linewidth = 1.0) +
      geom_hline(yintercept = input$prot, linetype = 2, colour = PAL[2]) +
      labs(x = "Days", y = "g/kg/day",
           title = "Dashed Line = Currently Prescribed Intake") + THEME
  })

  ## ------------------------------------------------------------- 6 compare
  scen <- eventReactive(input$goscen, {
    req(length(input$scen) > 0)
    UCD_simulate_scenarios(MOD, which = as.integer(input$scen))
  }, ignoreNULL = FALSE)

  cplot <- function(v, lab, logy = FALSE) {
    renderPlot({
      d <- scen(); req(nrow(d) > 0)
      g <- ggplot(d, aes(time/24, .data[[v]], colour = scenario)) +
        geom_line(linewidth = 0.8) +
        labs(x = "Days", y = lab, title = NULL) + THEME +
        theme(legend.text = element_text(size = 7))
      if (logy) g <- g + scale_y_log10()
      g
    })
  }
  output$c_nh3  <- cplot("NH3",      "Plasma ammonia (µmol/L)", TRUE)
  output$c_gln  <- cplot("GLNP",     "Plasma glutamine (µmol/L)")
  output$c_tol  <- cplot("PROTTOL",  "Protein tolerance (g/kg/day)")
  output$c_glnb <- cplot("GLNBRAIN", "Brain glutamine (mmol/L)")
  output$c_tab  <- renderDT({
    d <- scen(); req(nrow(d) > 0)
    datatable(UCD_summarise(d), options = list(scrollX = TRUE, pageLength = 14),
              rownames = FALSE)
  })

  ## ----------------------------------------------------------- 7 biomarker
  simple <- function(v, lab, hline = NULL, bandlo = NULL, bandhi = NULL,
                     col = PAL[1]) {
    renderPlot({
      g <- ggplot(sim(), aes(time/24, .data[[v]]))
      if (!is.null(bandlo)) g <- g + band(bandlo, bandhi)
      g <- g + geom_line(colour = col, linewidth = 0.9) +
        labs(x = "Days", y = lab, title = NULL) + THEME
      if (!is.null(hline)) g <- g + geom_hline(yintercept = hline, linetype = 2,
                                               colour = PAL[2])
      g
    })
  }
  output$b_upagn <- simple("UPAGN", "Urinary PAGN (mmol/day)", col = PAL[3])
  output$b_ratio <- simple("PARATIO", "PAA : PAGN (on a µg/mL basis)", hline = 2.5,
                           col = PAL[4])
  output$b_orot  <- simple("OROTC", "Orotic acid pool (µmol)", col = PAL[5])
  output$b_cit   <- simple("CITP", "Citrulline (µmol/L)", bandlo = 20, bandhi = 45)
  output$b_arg   <- simple("ARGP", "Arginine (µmol/L)", bandlo = 50, bandhi = 130,
                           col = PAL[6])
  output$b_gly   <- simple("GLYC", "Glycine (µmol/L)", bandlo = 180, bandhi = 320,
                           col = PAL[7])
  output$b_bcaa  <- simple("BCAC", "BCAA (µmol/L)", bandlo = 300, bandhi = 500,
                           col = PAL[8])

  ## -------------------------------------------------------------- 8 safety
  output$s_paa <- renderPlot({
    ggplot(sim(), aes(time/24, PAAUG)) +
      band(0, 100) +
      geom_hline(yintercept = 500, linetype = 1, colour = PAL[2], linewidth = 0.8) +
      geom_line(colour = PAL[1], linewidth = 1.0) +
      labs(x = "Days", y = "Plasma PAA (µg/mL)",
           title = "Red Line = 500 µg/mL Neurotoxicity Threshold") + THEME
  })
  output$s_paatox <- renderPlot({
    ggplot(sim(), aes(time/24, PAATOX)) +
      geom_line(colour = PAL[2], linewidth = 1.0) +
      labs(x = "Days", y = "µg/mL·h (excess above 500)", title = NULL) + THEME
  })
  output$s_na <- renderPlot({
    ## NaPBA 186.2 g/mol and Na-benzoate 144.1 g/mol each carry 1 mol Na
    na <- input$napba/186.2*1000 + input$bz*input$wt/1000/144.11*1000
    df <- data.frame(source = c("NaPBA", "Na-benzoate"),
                     mmol = c(input$napba/186.2*1000,
                              input$bz*input$wt/1000/144.11*1000))
    ggplot(df, aes(source, mmol, fill = source)) +
      geom_col(width = 0.6) +
      geom_hline(yintercept = 100, linetype = 2, colour = PAL[2]) +
      scale_fill_manual(values = PAL[c(1, 4)]) +
      labs(x = NULL, y = "mmol Na/day",
           title = paste0("Total ", round(na), " mmol Na/day (dashed = 100)")) +
      THEME + theme(legend.position = "none")
  })
  output$s_rebound <- renderPlot({
    d <- sim()
    g <- ggplot(d, aes(time/24, NH3)) +
      geom_line(colour = PAL[1], linewidth = 1.0) +
      labs(x = "Days", y = "Plasma ammonia (µmol/L)",
           title = "After Dialysis Ends, Nitrogen Redistributes from the Deep Glutamine Reservoir") + THEME
    if (input$hd > 0) {
      g <- g + annotate("rect", xmin = input$hdt/24,
                        xmax = (input$hdt + input$hdd)/24,
                        ymin = -Inf, ymax = Inf, fill = PAL[3], alpha = 0.15)
    }
    g
  })
}

shinyApp(ui, server)
