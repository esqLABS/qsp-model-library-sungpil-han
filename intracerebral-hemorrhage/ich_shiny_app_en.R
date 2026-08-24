## =============================================================================
##  ich_shiny_app.R
##  Interactive dashboard for the spontaneous intracerebral haemorrhage (ICH)
##  QSP model.
##
##  Run with:
##      shiny::runApp("ich_shiny_app.R")
##  or from this directory:
##      R -e 'shiny::runApp(".")'
##
##  The app sources ich_mrgsolve_model.R for the compiled model, the scenario
##  library and the dosing helpers, so the two files can never drift apart.
##
##  EIGHT TABS
##   1. Patient & Lesion           — covariates, the volume budget
##   2. Drug PK                    — every drug on one time axis
##   3. Haemostasis                — the clot-competence factor
##   4. Haematoma · bleeding flux  — the product that drives dV/dt
##   5. ICP & Perfusion            — the two appearances of MAP
##   6. Iron & Inflammation        — the delayed chemical channel
##   7. Clinical Endpoints         — NIHSS, GCS, mRS, mortality
##   8. Scenario Comparison        — all 12 arms + the U-shape surface
## =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)

## ---------------------------------------------------------------------------
##  Load the model, scenarios and dosing helpers
## ---------------------------------------------------------------------------
.here <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) ".")
if (is.null(.here) || !nzchar(.here)) .here <- "."
source(file.path(.here, "ich_mrgsolve_model.R"), local = FALSE)

THEME <- theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "grey92"),
        plot.title = element_text(face = "bold", size = 12))

PAL <- c("#3060b0", "#c05050", "#2e8b57", "#8000a0", "#e08050",
         "#7fa8b8", "#b8845f", "#555555")

## helper: long-format selected columns
gather_vars <- function(d, vars) {
  d %>% select(time, all_of(vars)) %>%
    pivot_longer(-time, names_to = "variable", values_to = "value") %>%
    mutate(variable = factor(variable, levels = vars))
}

## ===========================================================================
##  UI
## ===========================================================================
ui <- fluidPage(
  titlePanel("Spontaneous Intracerebral Haemorrhage (ICH) — QSP Model"),
  tags$p(style = "color:#666; margin-top:-10px; font-size:90%;",
         "Bleeding flux = driving pressure x haemostatic incompetence x number of open bleeding sites. ",
         "Outcome = the mass-effect channel + the haem·iron channel. ",
         "MAP appears in two places with opposite signs (CPP = MAP - ICP)."),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      selectInput("scenario", "Scenario",
                  choices = names(scenarios),
                  selected = "3_bp_intensive_140"),
      helpText(textOutput("scen_note")),
      tags$hr(),

      h4("Patient"),
      sliderInput("VHEM0", "Initial haematoma volume (mL)", 5, 90, 30, step = 1),
      sliderInput("AGE", "Age (years)", 40, 95, 68, step = 1),
      sliderInput("WT", "Body weight (kg)", 45, 130, 75, step = 1),
      selectInput("LOC", "Location",
                  choices = c("Deep" = 0, "Lobar" = 1,
                              "Infratentorial" = 2),
                  selected = 0),
      sliderInput("SVD", "Small-vessel disease burden (0-1)", 0, 1, 0.4, step = 0.05),
      sliderInput("FIVH", "Fraction with intraventricular extension", 0, 0.6, 0.10, step = 0.02),
      sliderInput("SBPBASE", "Untreated SBP set point (mmHg)", 130, 240, 185, step = 5),
      sliderInput("HPG", "Haptoglobin scavenging capacity (Hp1-1=1.0, Hp2-2=0.6)",
                  0.4, 1.2, 1.0, step = 0.05),

      tags$hr(),
      h4("Antithrombotic on board"),
      sliderInput("FII0", "Prothrombin activity (%)  [20% ~ INR 3.0]",
                  10, 100, 100, step = 5),
      sliderInput("APIX0", "Apixaban amount in the body (mg)", 0, 8, 0, step = 0.2),
      sliderInput("DABI0", "Dabigatran amount in the body (mg)", 0, 30, 0, step = 1),
      sliderInput("IASA", "Aspirin platelet-function loss", 0, 0.6, 0, step = 0.05),
      sliderInput("IP2Y12", "P2Y12 inhibitor function loss", 0, 0.6, 0, step = 0.05),

      tags$hr(),
      h4("Interventions"),
      checkboxInput("EVDON", "External ventricular drain (EVD)", FALSE),
      checkboxInput("SURGON", "Haematoma evacuation", FALSE),
      conditionalPanel("input.SURGON",
        sliderInput("FEVAC", "Fraction evacuated", 0.2, 0.95, 0.70, step = 0.05),
        sliderInput("TSURG", "Time of surgery (h)", 4, 96, 30, step = 2)),
      checkboxInput("DCON", "Decompressive craniectomy", FALSE),
      checkboxInput("INSON", "Insulin protocol", FALSE),
      checkboxInput("COOLON", "Fever management / cooling", FALSE),
      checkboxInput("VITK", "Vitamin K 10 mg IV", FALSE),
      sliderInput("REHAB", "Rehabilitation intensity (0-1.5)", 0, 1.5, 1.0, step = 0.1),

      tags$hr(),
      sliderInput("tmax", "Simulation duration (days)", 2, 90, 14, step = 1),
      actionButton("run", "Run", class = "btn-primary")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",

        ## ---- 1 --------------------------------------------------------
        tabPanel("1. Patient · lesion",
          h4("The Monro-Kellie volume budget"),
          p("ICP is not a state variable but the residual of this budget. The curve is flat until the CSF reserve is exhausted, and steep from then on."),
          plotOutput("p_budget", height = "300px"),
          fluidRow(
            column(6, h4("Haematoma · oedema · intraventricular blood"), plotOutput("p_vols", height = "260px")),
            column(6, h4("Pressure-volume relation (the realised trajectory)"), plotOutput("p_pv", height = "260px"))
          ),
          h4("Key numbers"), tableOutput("t_summary")
        ),

        ## ---- 2 --------------------------------------------------------
        tabPanel("2. Drug PK",
          h4("Concentration and effect site of the drugs given"),
          p("An empty panel means that drug was not given in that scenario."),
          plotOutput("p_pk", height = "420px"),
          fluidRow(
            column(6, h4("Reversal-agent binding (andexanet / idarucizumab)"),
                      plotOutput("p_reversal", height = "280px")),
            column(6, h4("Free anticoagulant activity"),
                      plotOutput("p_anticoag", height = "280px"))
          ),
          p(strong("Rebound is not coded as a separate term."),
            " The rise in anti-Xa after andexanet is stopped emerges from complex dissociation (KOFFA) and peripheral redistribution alone. Set KOFFA = 0 and it disappears.")
        ),

        ## ---- 3 --------------------------------------------------------
        tabPanel("3. Haemostasis",
          h4("Clot competence — the second factor of the bleeding flux"),
          plotOutput("p_clot", height = "320px"),
          fluidRow(
            column(6, h4("Coagulation factors · platelets"), plotOutput("p_coag", height = "280px")),
            column(6, h4("Thrombin generation · fibrinolysis"), plotOutput("p_thr", height = "280px"))
          ),
          p("CLOT is built from ", strong("two arms"), ": a platelet arm (WPLT, unaffected by VKAs·DOACs) and a coagulation arm (WFIB). Because of this separation, platelet transfusion and DDAVP do nothing for warfarin bleeding, and an INR of 3.0 increases bleeding about 2-fold rather than 3.5-fold.")
        ),

        ## ---- 4 --------------------------------------------------------
        tabPanel("4. Haematoma · bleeding flux",
          h4("dV/dt = KBLEED x NOPEN x (1 - CLOT) x max(0, MAP - PTISS)"),
          p("The four factors are shown together. If any one of them reaches zero the whole flux disappears — this is why a single-factor intervention has a small effect and combinations are synergistic."),
          plotOutput("p_flux", height = "420px"),
          fluidRow(
            column(6, h4("Haematoma volume and expansion"), plotOutput("p_vhem", height = "280px")),
            column(6, h4("Tamponade: the tissue counter-pressure"), plotOutput("p_tamp", height = "280px"))
          )
        ),

        ## ---- 5 --------------------------------------------------------
        tabPanel("5. ICP · perfusion",
          h4("The two faces of MAP"),
          p("The same MAP is both the driving pressure that pushes the bleeding out (lower is better) and the numerator of CPP (lower is worse). The U-shaped curve is coded into none of the equations."),
          plotOutput("p_icp", height = "380px"),
          fluidRow(
            column(6, h4("Autoregulation · relative perfusion"), plotOutput("p_auto", height = "280px")),
            column(6, h4("Cumulative exposure (observed quantities)"), plotOutput("p_accum", height = "280px"))
          )
        ),

        ## ---- 6 --------------------------------------------------------
        tabPanel("6. Iron · inflammation",
          h4("The delayed chemical channel: erythrocyte lysis → haem → HO-1 → free Fe2+ → lipid peroxidation"),
          plotOutput("p_iron", height = "380px"),
          fluidRow(
            column(6, h4("Neuroinflammation"), plotOutput("p_inflam", height = "280px")),
            column(6, h4("Tissue survival (neurons · white matter)"), plotOutput("p_tissue", height = "280px"))
          ),
          p("Deferoxamine has ", strong("not a single edge leading to VHEM"),
            " — it acts on this channel only. So the 24-hour haematoma volume is identical to control while only the 90-day outcome moves (the i-DEF pattern).")
        ),

        ## ---- 7 --------------------------------------------------------
        tabPanel("7. Clinical endpoints",
          h4("NIHSS · GCS · 90-day mRS · probability of death"),
          plotOutput("p_clin", height = "380px"),
          fluidRow(
            column(6, h4("Components of the ICH score"), tableOutput("t_ichscore")),
            column(6, h4("Outcome probabilities"), plotOutput("p_outcome", height = "260px"))
          ),
          p(strong("Note: "), "The self-fulfilling prophecy created by DNR orders and withdrawal of care (Becker 2001) confounds the interpretation of every ICH trial. The model includes it as a node on the map but does not implement it as an ODE — the mortality here reflects the biological trajectory alone.")
        ),

        ## ---- 8 --------------------------------------------------------
        tabPanel("8. Scenario comparison",
          h4("Twelve scenarios"),
          p("The table below computes the measure each trial actually reported, so that it can be compared directly with the literature anchors. The computation takes 30-90 seconds."),
          actionButton("run_all", "Run all 12 scenarios", class = "btn-warning"),
          tags$br(), tags$br(),
          tableOutput("t_all"),
          tags$hr(),
          h4("The emergent U-shape"),
          p("Sweeps how deeply the blood pressure is lowered, from 0 to 1. Haematoma expansion falls monotonically but exposure to CPP<60 rises monotonically, so the outcome is best in the middle."),
          actionButton("run_surf", "Compute the response surface", class = "btn-warning"),
          plotOutput("p_surface", height = "340px"),
          tableOutput("t_surface")
        )
      )
    )
  )
)

## ===========================================================================
##  SERVER
## ===========================================================================
server <- function(input, output, session) {

  output$scen_note <- renderText({
    s <- scenarios[[input$scenario]]
    if (is.null(s)) "" else s$note
  })

  ## when a scenario is chosen, adopt its parameter overrides in the sidebar
  observeEvent(input$scenario, {
    p <- scenarios[[input$scenario]]$param
    if (!is.null(p$VHEM0))   updateSliderInput(session, "VHEM0",   value = p$VHEM0)
    if (!is.null(p$SBPBASE)) updateSliderInput(session, "SBPBASE", value = p$SBPBASE)
    if (!is.null(p$FII0))    updateSliderInput(session, "FII0",    value = p$FII0)
    if (!is.null(p$APIX0))   updateSliderInput(session, "APIX0",   value = p$APIX0)
    if (!is.null(p$FIVH))    updateSliderInput(session, "FIVH",    value = p$FIVH)
    updateCheckboxInput(session, "EVDON",  value = isTRUE(p$EVDON  == 1))
    updateCheckboxInput(session, "SURGON", value = isTRUE(p$SURGON == 1))
    updateCheckboxInput(session, "INSON",  value = isTRUE(p$INSON  == 1))
    updateCheckboxInput(session, "COOLON", value = isTRUE(p$COOLON == 1))
    updateCheckboxInput(session, "VITK",   value = isTRUE(p$VITK   == 1))
  })

  ## ---- the single simulation the first seven tabs share -----------------
  sim <- eventReactive(input$run, {
    s   <- scenarios[[input$scenario]]
    end <- input$tmax * 24

    pars <- c(
      s$param,
      list(VHEM0 = input$VHEM0, AGE = input$AGE, WT = input$WT,
           LOC = as.numeric(input$LOC), SVD = input$SVD, FIVH = input$FIVH,
           SBPBASE = input$SBPBASE, HPG = input$HPG,
           FII0 = input$FII0, APIX0 = input$APIX0, DABI0 = input$DABI0,
           IASA = input$IASA, IP2Y12 = input$IP2Y12,
           EVDON = as.numeric(input$EVDON), SURGON = as.numeric(input$SURGON),
           DCON = as.numeric(input$DCON), INSON = as.numeric(input$INSON),
           COOLON = as.numeric(input$COOLON), VITK = as.numeric(input$VITK),
           REHAB = input$REHAB)
    )
    ## sidebar wins over the scenario's own overrides for the shared covariates
    pars <- pars[!duplicated(names(pars), fromLast = TRUE)]
    if (input$SURGON) {
      pars$FEVAC <- input$FEVAC
      pars$TSURG <- input$TSURG
    }

    m <- param(ich, pars)
    withProgress(message = "Integrating...", value = 0.5, {
      out <- if (is.null(s$events)) mrgsim(m, end = end, delta = 0.25)
             else mrgsim(m, events = s$events, end = end, delta = 0.25)
    })
    as_tibble(out) %>% mutate(day = time / 24)
  }, ignoreNULL = FALSE)

  ## ================= TAB 1 : patient & lesion ==========================
  output$p_budget <- renderPlot({
    d <- sim()
    d %>% transmute(day,
                    `Haematoma VHEM` = VHEM,
                    `Early oedema OEDE` = OEDE,
                    `Late oedema OEDL` = OEDL,
                    `IVH blood VIVH` = VIVH,
                    `CSF displacement (negative = compensated)` = VCSF - 140) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value, fill = name)) +
      geom_area(alpha = 0.85) +
      geom_line(data = transmute(d, day, value = VADD, name = "Uncompensated volume VADD"),
                aes(day, value), inherit.aes = FALSE,
                colour = "black", linewidth = 1.1, linetype = "dashed") +
      scale_fill_manual(values = PAL) +
      labs(x = "Days", y = "Volume (mL)", fill = NULL,
           subtitle = "Dashed = the residual volume not absorbed by compensation (it sets ICP)") + THEME
  })

  output$p_vols <- renderPlot({
    gather_vars(sim(), c("VHEM", "OEDE", "OEDL", "VIVH")) %>%
      mutate(day = time / 24) %>%
      ggplot(aes(day, value, colour = variable)) +
      geom_line(linewidth = 0.8) + scale_colour_manual(values = PAL) +
      labs(x = "Days", y = "mL", colour = NULL) + THEME
  })

  output$p_pv <- renderPlot({
    ggplot(sim(), aes(VADD, ICP, colour = day)) +
      geom_path(linewidth = 1) +
      scale_colour_viridis_c(option = "plasma") +
      labs(x = "Uncompensated volume VADD (mL)", y = "ICP (mmHg)", colour = "Days",
           subtitle = "The slope differs before and after the CSF reserve is exhausted") + THEME
  })

  output$t_summary <- renderTable({
    d <- sim(); e <- max(d$time)
    at <- function(v, t) d[[v]][which.min(abs(d$time - t))]
    tibble(
      `Measure` = c("Initial haematoma (mL)", "24-hour haematoma (mL)", "Expansion (mL)", "Expansion (%)",
               "Peak PHE (mL)", "Peak PHE (day)", "Peak ICP (mmHg)",
               "Lowest CPP (mmHg)", "Cumulative CPP<60 (h)", "Cumulative ICP>20 (h)",
               "Peak midline shift (mm)", "NIHSS (24 h)", "NIHSS (final)",
               "Lowest GCS", "ICH score (24 h)", "P(mRS 0-2)", "P(death)"),
      `Value` = c(round(at("VHEM", 0), 1), round(at("VHEM", 24), 1),
             round(at("VHEM", 24) - at("VHEM", 0), 2),
             round(100 * (at("VHEM", 24) - at("VHEM", 0)) / at("VHEM", 0), 1),
             round(max(d$PHE), 1), round(d$time[which.max(d$PHE)] / 24, 2),
             round(max(d$ICP), 1), round(min(d$CPP), 1),
             round(at("TCPP", e), 1), round(at("TICP", e), 1),
             round(max(d$MIDLINE), 1),
             round(at("NIH", 24), 1), round(at("NIH", e), 1),
             round(min(d$GCSE), 1), round(at("ICHSC", 24), 0),
             round(at("PMRS02", e), 3), round(at("PMORT", e), 3))
    )
  }, digits = 3)

  ## ================= TAB 2 : PK ========================================
  output$p_pk <- renderPlot({
    d <- sim()
    d %>% transmute(day,
                    `Nicardipine central (mg)` = NICA1,
                    `Nicardipine effect site` = NICE,
                    `Labetalol (mg)` = LABE,
                    `TXA central (mg)` = TXAC,
                    `4F-PCC remaining (IU)` = PCCA,
                    `Deferoxamine (mg)` = DFOC,
                    `Mannitol (g)` = MANN,
                    `Intraventricular alteplase (mg)` = IVTPA) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value)) + geom_line(colour = PAL[1], linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y", ncol = 4) +
      labs(x = "Days", y = NULL) + THEME
  })

  output$p_reversal <- renderPlot({
    sim() %>% transmute(day,
                        `Free andexanet (NEQ)` = ANDXE,
                        `Andexanet-apixaban complex` = CPLXA,
                        `Free idarucizumab (NEQ)` = IDAE,
                        `Ida-dabi complex` = CPLXD) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value, colour = name)) + geom_line(linewidth = 0.8) +
      scale_colour_manual(values = PAL) +
      labs(x = "Days", y = "mg-equivalents", colour = NULL) + THEME
  })

  output$p_anticoag <- renderPlot({
    sim() %>% transmute(day, `Free apixaban anti-Xa (ng/mL)` = AXA,
                        `Free dabigatran (ng/mL)` = DTT, INR = INR) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value)) + geom_line(colour = PAL[2], linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      labs(x = "Days", y = NULL) + THEME
  })

  ## ================= TAB 3 : haemostasis ===============================
  output$p_clot <- renderPlot({
    sim() %>% transmute(day, `Clot competence CLOT` = CLOT,
                        `Open bleeding sites NOPEN` = NOPEN,
                        `Thrombin generation THRGEN` = THRGEN) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value, colour = name)) + geom_line(linewidth = 1) +
      scale_colour_manual(values = PAL) +
      labs(x = "Days", y = "Relative value", colour = NULL) + THEME
  })

  output$p_coag <- renderPlot({
    sim() %>% transmute(day, `Prothrombin activity (%)` = FII,
                        `Fibrinogen (mg/dL)` = FIB,
                        `Platelet function (0-1)` = PLTF) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value)) + geom_line(colour = PAL[3], linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      labs(x = "Days", y = NULL) + THEME
  })

  output$p_thr <- renderPlot({
    sim() %>% transmute(day, `Plasmin activity` = PLAS,
                        `Tissue thrombin THRT` = THRT,
                        `INR` = INR) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value)) + geom_line(colour = PAL[4], linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      labs(x = "Days", y = NULL) + THEME
  })

  ## ================= TAB 4 : bleeding flux =============================
  output$p_flux <- renderPlot({
    sim() %>% transmute(day,
                        `Bleeding flux (mL/h)` = BLEEDR,
                        `Driving pressure PDRIVE (mmHg)` = PDRIVE,
                        `Haemostatic incompetence (1-CLOT)` = 1 - CLOT,
                        `Open bleeding sites NOPEN` = NOPEN) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value)) +
      geom_line(colour = PAL[2], linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "Days", y = NULL,
           subtitle = "Top left = the product of the three factors") + THEME
  })

  output$p_vhem <- renderPlot({
    d <- sim(); v0 <- d$VHEM[1]
    ggplot(d, aes(day, VHEM)) +
      geom_line(colour = PAL[2], linewidth = 1) +
      geom_hline(yintercept = v0, linetype = "dotted") +
      geom_hline(yintercept = v0 * 1.33, linetype = "dashed", colour = "grey40") +
      annotate("text", x = max(d$day) * 0.7, y = v0 * 1.36,
               label = "Definition of expansion (+33%)", size = 3.2, colour = "grey30") +
      labs(x = "Days", y = "Haematoma volume (mL)") + THEME
  })

  output$p_tamp <- renderPlot({
    sim() %>% transmute(day, `MAP` = MAP, `Tissue counter-pressure PTISS` = PTISS,
                        `ICP` = ICP) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value, colour = name)) + geom_line(linewidth = 0.9) +
      scale_colour_manual(values = PAL) +
      labs(x = "Days", y = "mmHg", colour = NULL,
           subtitle = "When MAP and PTISS meet, the bleeding stops of its own accord") + THEME
  })

  ## ================= TAB 5 : ICP & perfusion ===========================
  output$p_icp <- renderPlot({
    sim() %>% transmute(day, `SBP` = SBP, `MAP` = MAP, `ICP` = ICP,
                        `CPP` = CPP) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value, colour = name)) + geom_line(linewidth = 0.9) +
      geom_hline(yintercept = 60, linetype = "dashed", colour = "grey50") +
      geom_hline(yintercept = 20, linetype = "dotted", colour = "grey50") +
      scale_colour_manual(values = PAL) +
      labs(x = "Days", y = "mmHg", colour = NULL,
           subtitle = "Dashed = the CPP 60 threshold, dotted = the ICP 20 threshold") + THEME
  })

  output$p_auto <- renderPlot({
    sim() %>% transmute(day, `Autoregulatory capacity AUTOR` = AUTOR,
                        `Relative perfusion CBFP` = CBFP,
                        `Ischaemic stress ISCH` = ISCH) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value, colour = name)) + geom_line(linewidth = 0.9) +
      scale_colour_manual(values = PAL) +
      labs(x = "Days", y = "Relative value", colour = NULL) + THEME
  })

  output$p_accum <- renderPlot({
    sim() %>% transmute(day, `SBP>140 AUC` = ASBP, `Hours CPP<60` = TCPP,
                        `Hours ICP>20` = TICP, `Fe2+ AUC` = AFE,
                        `IL-6 AUC` = AIL6) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value)) + geom_line(colour = PAL[6], linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "Days", y = NULL,
           subtitle = "These are not parameters but integrals, that is, observed quantities") + THEME
  })

  ## ================= TAB 6 : iron & inflammation =======================
  output$p_iron <- renderPlot({
    sim() %>% transmute(day, `Free haem HEME` = HEME, `Free Fe2+ FEII` = FEII,
                        `Ferritin sequestration FERR` = FERR,
                        `Lipid peroxidation LPO` = LPO,
                        `Feroxamine (excreted)` = FERX) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value, colour = name)) + geom_line(linewidth = 0.9) +
      scale_colour_manual(values = PAL) +
      labs(x = "Days", y = "Relative value", colour = NULL) + THEME
  })

  output$p_inflam <- renderPlot({
    sim() %>% transmute(day, `Pro-inflammatory microglia MG1` = MG1,
                        `Repair-type MG2` = MG2, `Neutrophils NEU` = NEU,
                        `IL-6` = IL6, `IL-10` = IL10,
                        `MMP-9` = MMP9, `BBB permeability` = BBBP) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value, colour = name)) + geom_line(linewidth = 0.8) +
      scale_colour_manual(values = c(PAL, "#999999")) +
      labs(x = "Days", y = "Relative value", colour = NULL) + THEME
  })

  output$p_tissue <- renderPlot({
    sim() %>% transmute(day, `Surviving neuron fraction NEUR` = NEUR,
                        `White matter integrity WMI` = WMI,
                        `Plasticity reserve PLAST` = PLAST,
                        `Mechanical strain STRAIN` = STRAIN) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value, colour = name)) + geom_line(linewidth = 0.9) +
      scale_colour_manual(values = PAL) +
      labs(x = "Days", y = "0-1 (STRAIN is an index)", colour = NULL) + THEME
  })

  ## ================= TAB 7 : clinical ==================================
  output$p_clin <- renderPlot({
    sim() %>% transmute(day, `NIHSS` = NIH, `Estimated GCS` = GCSE,
                        `Temperature (degC)` = TEMP, `Glucose (mmol/L)` = GLU,
                        `Midline shift (mm)` = MIDLINE,
                        `ICH score` = ICHSC) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value)) + geom_line(colour = PAL[1], linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "Days", y = NULL) + THEME
  })

  output$t_ichscore <- renderTable({
    d <- sim(); at <- function(v, t) d[[v]][which.min(abs(d$time - t))]
    g <- at("GCSE", 24); v <- at("VHEM", 24); iv <- at("VIVH", 24)
    tibble(
      `Component` = c("GCS", "Volume >= 30 mL", "Intraventricular haemorrhage", "Infratentorial", "Age >= 80", "Total"),
      `Value` = c(sprintf("%.1f", g), sprintf("%.1f mL", v), sprintf("%.1f mL", iv),
             ifelse(as.numeric(input$LOC) > 1.5, "Yes", "No"),
             sprintf("%d years", input$AGE), ""),
      `Points` = c(ifelse(g <= 4, 2, ifelse(g <= 12, 1, 0)),
               ifelse(v >= 30, 1, 0), ifelse(iv > 0.5, 1, 0),
               ifelse(as.numeric(input$LOC) > 1.5, 1, 0),
               ifelse(input$AGE >= 80, 1, 0),
               at("ICHSC", 24))
    )
  })

  output$p_outcome <- renderPlot({
    sim() %>% transmute(day, `P(mRS 0-2)` = PMRS02, `P(death)` = PMORT,
                        `Utility-weighted mRS` = UWMRS) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value, colour = name)) + geom_line(linewidth = 1) +
      scale_colour_manual(values = PAL) + ylim(0, 1) +
      labs(x = "Days", y = "Probability / utility", colour = NULL) + THEME
  })

  ## ================= TAB 8 : scenario comparison =======================
  all_sim <- eventReactive(input$run_all, {
    withProgress(message = "Running the 12 scenarios...", value = 0.3, {
      summarise_scenarios()
    })
  })

  output$t_all <- renderTable({
    req(all_sim())
    all_sim() %>%
      transmute(`Scenario` = scenario,
                `Expansion (mL)` = expansion_mL,
                `Expansion (%)` = expansion_pct,
                `SBP nadir` = SBP_nadir_mmHg,
                `Peak PHE (mL)` = PHE_peak_mL,
                `Peak ICP` = ICP_peak_mmHg,
                `CPP<60 (h)` = hours_CPP_lt60,
                `Fe AUC (14d)` = AUC_Fe2_14d,
                `NIHSS 90d` = NIHSS_90d,
                `P(mRS 0-2)` = P_mRS02_90d,
                `P(death)` = P_death)
  }, digits = 3)

  surf <- eventReactive(input$run_surf, {
    withProgress(message = "Computing the response surface...", value = 0.3, {
      bp_response_surface(speeds = c(1, 4), end = 1440)
    })
  })

  output$p_surface <- renderPlot({
    req(surf())
    s <- surf()
    best <- s %>% group_by(ttt) %>% slice_max(P_mRS02, n = 1)
    ggplot(s, aes(SBP_nadir, P_mRS02, colour = factor(ttt))) +
      geom_line(linewidth = 1) + geom_point(size = 2) +
      geom_point(data = best, size = 5, shape = 21, stroke = 1.4,
                 fill = NA, colour = "black") +
      scale_x_reverse() +
      scale_colour_manual(values = PAL[c(1, 2)]) +
      labs(x = "SBP nadir achieved (mmHg) — deeper lowering towards the right",
           y = "P(mRS 0-2) at 90 d", colour = "Time to target (h)",
           subtitle = "Black circle = the optimum. The U-shaped curve is not coded, it emerges.") +
      THEME
  })

  output$t_surface <- renderTable({
    req(surf())
    surf() %>% filter(ttt == 1) %>%
      transmute(`Lowering intensity` = intensity, `SBP nadir` = SBP_nadir,
                `CPP nadir` = CPP_min, `Expansion (mL)` = expansion,
                `CPP<60 (h)` = hrs_CPP60, `Neuron survival` = NEUR_end,
                `NIHSS` = NIHSS_end, `P(mRS 0-2)` = P_mRS02)
  }, digits = 3)
}

shinyApp(ui, server)
