# =============================================================================
#  nhb_shiny_app.R
#  Neonatal hyperbilirubinaemia QSP model — interactive dashboard
#  Neonatal Hyperbilirubinaemia QSP Model · Shiny Dashboard
#  ---------------------------------------------------------------------------
#  Run:
#     shiny::runApp("nhb_shiny_app.R")
#  Requires nhb_mrgsolve_model.R in the same directory (it is sourced for the
#  compiled model object `mod`, the genotype table and the dosing helpers).
#
#  The app is built around the model's central claim: the number on the chart
#  (TSB) is not the number that injures (free bilirubin).  Tab 3 exists to make
#  that visible, and every other tab reports both.
#
#  NOT FOR CLINICAL USE.  Thresholds are analytic approximations of the AAP
#  2022 curves and must not be used to make care decisions.
# =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)

source("nhb_mrgsolve_model.R", local = TRUE)   # provides mod, GENOTYPES, doses

THEME <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92"),
        legend.position = "bottom")

# -----------------------------------------------------------------------------
#  UI
# -----------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel(
    div(
      h3("Neonatal Hyperbilirubinaemia QSP Simulator"),
      h5(em(paste("34-ODE mrgsolve model · haem catabolism → albumin binding →",
                  "UGT1A1 → enterohepatic shunt → photochemistry → BBB"))),
      h6(strong("For education and research purposes only — do not use for clinical judgement"),
         style = "color:#b03a2e")
    )
  ),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      h4("1 · Patient Profile / Infant"),
      sliderInput("ga", "Gestational age (wk)", 28, 42, 40, 1),
      sliderInput("bw", "Birth weight (kg)", 0.8, 4.8, 3.40, 0.05),
      sliderInput("alb", "Serum albumin (g/dL)", 1.8, 4.5, 3.40, 0.1),
      sliderInput("hb", "Cord haemoglobin (g/dL)", 8, 22, 17, 0.5),
      selectInput("geno", "UGT1A1 genotype",
                  choices = c("wild type (*1/*1)"        = "wild",
                              "*28 heterozygote"          = "UGT1A1_28_het",
                              "Gilbert (*28/*28)"         = "Gilbert",
                              "*6/*6 (East Asian, G71R)"  = "UGT1A1_6",
                              "Crigler-Najjar II"         = "CN2",
                              "Crigler-Najjar I"          = "CN1"),
                  selected = "wild"),

      h4("2 · Haemolytic Load"),
      sliderInput("abmat", "Alloantibody load", 0, 0.6, 0, 0.02),
      checkboxInput("g6pd", "G6PD deficient", FALSE),
      conditionalPanel(
        "input.g6pd == true",
        sliderInput("oxstart", "Oxidant challenge start (h)",
                    0, 168, 60, 6),
        sliderInput("oxdur", "Duration (h)", 0, 48, 24, 6)),
      sliderInput("cephal", "Cephalohaematoma blood volume (extravasated Hb, g)", 0, 20, 0, 1),

      h4("3 · Binding & Barrier Modifying Factors"),
      sliderInput("facid", "Acidaemia factor on KA", 0.5, 1.0, 1.0, 0.05),
      sliderInput("fdisp", "Displacer factor", 0.5, 1.0, 1.0, 0.05),
      sliderInput("fbbb", "BBB permeability (fold)", 1.0, 3.0, 1.0, 0.1),

      h4("4 · Feeding"),
      radioButtons("breast", "Feeding method",
                   c("Human milk" = 1, "Formula" = 0), 1, inline = TRUE),
      sliderInput("fiearly", "Early intake adequacy",
                  0.1, 1.0, 1.0, 0.05),
      sliderInput("fisw", "Intake improves at (h)", 0, 240, 96, 12),

      h4("5 · Therapy"),
      radioButtons("ptmode", "Phototherapy mode",
                   c("None"                       = "none",
                     "Threshold-driven"     = "auto",
                     "Fixed window"           = "fixed"),
                   "auto"),
      conditionalPanel(
        "input.ptmode != 'none'",
        sliderInput("irr", "Irradiance (µW/cm²/nm)", 0, 60, 30, 1),
        sliderInput("fbsa", "Exposed BSA fraction", 0.2, 1.0, 0.80, 0.05)),
      conditionalPanel(
        "input.ptmode == 'fixed'",
        sliderInput("ptwin", "Phototherapy window (h)", 0, 336, c(24, 96), 6)),
      conditionalPanel(
        "input.ptmode == 'auto'",
        sliderInput("ptduty", "Max hours/day", 6, 24, 24, 2)),

      checkboxInput("ivig", "IVIG 1 g/kg", FALSE),
      conditionalPanel("input.ivig == true",
                       sliderInput("tivig", "Administration time (h)", 0, 168, 20, 2)),
      checkboxInput("et", "Double-volume exchange transfusion", FALSE),
      conditionalPanel("input.et == true",
                       sliderInput("tet", "Start (h)", 6, 240, 30, 2),
                       sliderInput("etdur", "Duration (h)", 1.5, 5, 3, 0.5)),
      checkboxInput("snmp", "Stannsoporfin 4.5 mg/kg IM", FALSE),
      checkboxInput("pheno", "Phenobarbital 5 mg/kg/d", FALSE),
      checkboxInput("udca", "UDCA 10 mg/kg q12h", FALSE),
      checkboxInput("agar", "Oral agar 250 mg/kg/d", FALSE),
      checkboxInput("aav", "AAV8-hUGT1A1 gene transfer", FALSE),
      conditionalPanel("input.aav == true",
                       sliderInput("taav", "Administration time (days)", 1, 90, 60, 1),
                       sliderInput("tgxmax", "Transgene activity (% adult)",
                                   1, 40, 10, 1)),

      hr(),
      sliderInput("tend", "Simulation horizon (days)", 3, 60, 14, 1),
      actionButton("run", "Simulate",
                   class = "btn-primary", width = "100%")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs", type = "tabs",

        # ------------------------------------------------------------- TAB 1
        tabPanel(
          "① Bilirubin trajectory · Bilirubin & thresholds",
          br(),
          plotOutput("p_traj", height = "430px"),
          fluidRow(
            column(6, h5("Decision summary"),
                   tableOutput("t_decision")),
            column(6, h5("Weight & growth"), plotOutput("p_wt", height = "220px"))
          ),
          helpText(paste("Dashed line = approximate AAP 2022 phototherapy threshold, dotted line = exchange-transfusion threshold.",
                         "Thresholds are a function of gestational age, time, and risk factors."))
        ),

        # ------------------------------------------------------------- TAB 2
        tabPanel(
          "② Free Bilirubin",
          br(),
          fluidRow(
            column(6, plotOutput("p_bf", height = "300px")),
            column(6, plotOutput("p_bf_vs_tsb", height = "300px"))
          ),
          fluidRow(
            column(6, plotOutput("p_ba", height = "260px")),
            column(6, h5("Iso-risk (iso-Bf) contour"),
                   plotOutput("p_iso", height = "260px"))
          ),
          helpText(paste("This tab is the heart of the model. The relationship between the measured value (TSB) and",
                         "the value that causes injury (Bf) is a saturating binding isotherm, and when albumin is low or affinity",
                         "is reduced, the same TSB means an entirely different risk."))
        ),

        # ------------------------------------------------------------- TAB 3
        tabPanel(
          "③ Flux Balance · Production vs Clearance",
          br(),
          plotOutput("p_flux", height = "340px"),
          fluidRow(
            column(6, plotOutput("p_prod", height = "260px")),
            column(6, plotOutput("p_ugt", height = "260px"))
          ),
          helpText(paste("Because TSB is the integral of the difference between two fluxes, it cannot by itself",
                         "distinguish the cause. ETCOc (end-tidal CO) reads the production flux",
                         "directly, separating overproduction from underclearance."))
        ),

        # ------------------------------------------------------------- TAB 4
        tabPanel(
          "④ Phototherapy Dose-Response",
          br(),
          fluidRow(
            column(7, plotOutput("p_ptgrid", height = "340px")),
            column(5, plotOutput("p_photoflux", height = "340px"))
          ),
          h5("Photoisomers are assayed as 'bilirubin'"),
          plotOutput("p_isomer", height = "300px"),
          helpText(paste("Increasing exposed body-surface area from 0.35 → 0.80 is far more effective",
                         "than raising irradiance from 30 → 50 µW/cm²/nm (optical saturation)."))
        ),

        # ------------------------------------------------------------- TAB 5
        tabPanel(
          "⑤ Exchange Transfusion",
          br(),
          plotOutput("p_et", height = "420px"),
          tableOutput("t_et"),
          helpText(paste("Exchange transfusion is not 'fast phototherapy'. Because donor plasma resets",
                         "albumin and removes maternal antibody, it acts simultaneously on load, binding capacity, and production,",
                         "so the fall in Bf is larger than the fall in TSB."))
        ),

        # ------------------------------------------------------------- TAB 6
        tabPanel(
          "⑥ Drug PK/PD",
          br(),
          fluidRow(
            column(6, plotOutput("p_pk1", height = "300px")),
            column(6, plotOutput("p_pk2", height = "300px"))
          ),
          plotOutput("p_pdeff", height = "280px"),
          helpText(paste("Because phenobarbital induces UGT1A1 transcription, it is by definition ineffective",
                         "in CN-I (GENO = 0). You cannot induce a gene product that does not exist."))
        ),

        # ------------------------------------------------------------- TAB 7
        tabPanel(
          "⑦ Neurotoxicity",
          br(),
          fluidRow(
            column(6, plotOutput("p_brain", height = "300px")),
            column(6, plotOutput("p_inj", height = "300px"))
          ),
          fluidRow(
            column(6, plotOutput("p_abr", height = "260px")),
            column(6, h5("Risk summary"), tableOutput("t_risk"))
          ),
          helpText(paste("ABR deficits are modelled as reversible; cumulative damage is modelled as irreversible.",
                         "This corresponds to the clinical observation that auditory brainstem responses normalise after treatment."))
        ),

        # ------------------------------------------------------------- TAB 8
        tabPanel(
          "⑧ Biomarkers",
          br(),
          fluidRow(
            column(6, plotOutput("p_etco", height = "290px")),
            column(6, plotOutput("p_hb", height = "290px"))
          ),
          fluidRow(
            column(6, plotOutput("p_tcb", height = "270px")),
            column(6, plotOutput("p_stool", height = "270px"))
          ),
          helpText(paste("ETCOc is a direct indicator of the heme-breakdown flux (normal ~1.4 ppm, haemolysis 2-4 ppm),",
                         "while TcB reflects the skin compartment and so dissociates from serum during phototherapy."))
        ),

        # ------------------------------------------------------------- TAB 9
        tabPanel(
          "⑨ Scenario Comparison · Scenario Library",
          br(),
          checkboxGroupInput("scen", "Scenarios to compare",
                             choices = names(scenarios),
                             selected = c("S1", "S3", "S5", "S6", "S7"),
                             inline = TRUE),
          actionButton("runscen", "Run library", class = "btn-info"),
          br(), br(),
          plotOutput("p_scen", height = "420px"),
          tableOutput("t_scen")
        ),

        # ------------------------------------------------------------ TAB 10
        tabPanel(
          "⑩ Crigler-Najjar Growth Analysis · Growth Decomposition",
          br(),
          h5(paste("Why does the phototherapy ceiling rise as the child grows —",
                   "the surface-area/weight term is offset by a fall in production")),
          plotOutput("p_cn", height = "380px"),
          tableOutput("t_cn"),
          helpText(paste("[G] geometry only / [G+P] + per-kg decline in production / [G+P+O] + skin optics.",
                         "BSA/W falls to 0.42-fold from birth to age 18, but bilirubin production per kg",
                         "also falls to 0.45-fold, so the two cancel out, and geometry alone cannot explain",
                         "the phototherapy failures observed clinically."))
        ),

        # ------------------------------------------------------------ TAB 11
        tabPanel(
          "⑪ Model Card",
          br(),
          h4("34 State Variables · 34 ODEs"),
          verbatimTextOutput("t_cmt"),
          h4("Current parameter set"),
          verbatimTextOutput("t_par"),
          h4("Known biases"),
          tags$ul(
            tags$li(paste("Double-volume exchange transfusion removes 40-45% of total body load",
                          "(the classic teaching figure is ~25%). This results from treating the extravascular",
                          "compartment as a single homogeneous compartment, and for the same reason the CN-I phototherapy ceiling is optimistic.")),
            tags$li(paste("TAUUGT lumps UGT1A1 protein maturation and OATP1B1/ligandin maturation into a single",
                          "time constant, fitted to the TSB trajectory.")),
            tags$li("The AAP 2022 thresholds are an analytical approximation, not the tabulated values."),
            tags$li(paste("The direction of the skin-optics terms (FOPTMIN, TAUOPT) is well established, but",
                          "the values are illustrative."))
          )
        )
      )
    )
  )
)

# -----------------------------------------------------------------------------
#  SERVER
# -----------------------------------------------------------------------------
server <- function(input, output, session) {

  # -- assemble the parameter set from the sidebar ---------------------------
  pset <- reactive({
    p <- list(
      GA = input$ga, W0 = input$bw, ALB0 = input$alb, ALBSET = max(input$alb, 3.0),
      HB0 = input$hb, GENO = GENOTYPES[[input$geno]],
      ABMAT0 = input$abmat, BLEX0 = input$cephal,
      G6PD = as.numeric(input$g6pd),
      FACID = input$facid, FDISP = input$fdisp, FBBB = input$fbbb,
      BREAST = as.numeric(input$breast),
      FIEARLY = input$fiearly, FISW = input$fisw, FILATE = 1.0,
      # gestational age drives RBC lifespan, binding affinity and the barrier
      LRBC  = 50 + 30*pmin(1, pmax(0, (input$ga - 32)/8)),
      FMATK = 0.70 + 0.30*pmin(1, pmax(0, (input$ga - 30)/8)),
      # AAP neurotoxicity risk factors, as the guideline defines them
      RF = as.numeric(input$ga < 38 || input$alb < 3.0 || input$abmat > 0 ||
                      input$g6pd || input$facid < 1.0)
    )
    if (input$g6pd) {
      p$OXSTART <- input$oxstart; p$OXDUR <- input$oxdur; p$OXSET <- 1.0
    }
    if (input$ptmode == "auto") {
      p$AUTOPT <- 1; p$IRRSET <- input$irr; p$FBSASET <- input$fbsa
      p$PTOFFM <- 2.0
    } else if (input$ptmode == "fixed") {
      p$AUTOPT <- 0; p$IRRSET <- input$irr; p$FBSASET <- input$fbsa
      p$PTSTART <- input$ptwin[1]; p$PTSTOP <- input$ptwin[2]
    }
    if (input$et) {
      p$ETSTART <- input$tet; p$ETDUR <- input$etdur; p$ETRATE <- 170/input$etdur
    }
    if (input$aav) {
      p$GTSTART <- input$taav*24; p$TGXMAX <- input$tgxmax/100
    }
    p
  })

  eset <- reactive({
    wt <- input$bw
    e <- NULL
    add <- function(a, b) if (is.null(a)) b else c(a, b)
    if (input$ivig)  e <- add(e, dose_ivig(wt, input$tivig))
    if (input$snmp)  e <- add(e, dose_snmp(wt, 24))
    if (input$pheno) e <- add(e, dose_pheno(wt, 0, input$tend))
    if (input$udca)  e <- add(e, dose_udca(wt, 0, input$tend))
    if (input$agar)  e <- add(e, dose_agar(wt, 0, input$tend))
    e
  })

  sim <- eventReactive(input$run, {
    m <- param(mod, pset())
    end <- input$tend*24
    e <- eset()
    out <- if (is.null(e)) mrgsim(m, end = end, delta = 0.5, hmax = 0.25)
           else mrgsim(m, events = e, end = end, delta = 0.5, hmax = 0.25)
    as_tibble(out)
  }, ignoreNULL = FALSE)

  # ------------------------------------------------------------------ TAB 1
  output$p_traj <- renderPlot({
    d <- sim()
    shade <- d %>% mutate(on = IRRAD > 0)
    ggplot(d, aes(time/24)) +
      geom_area(data = filter(shade, on), aes(y = Inf), fill = "#5dade2",
                alpha = 0.12) +
      geom_line(aes(y = THRESHPT, colour = "AAP phototherapy threshold"),
                linetype = 2, linewidth = 0.7) +
      geom_line(aes(y = THRESHET, colour = "AAP exchange threshold"),
                linetype = 3, linewidth = 0.7) +
      geom_line(aes(y = TSBOUT, colour = "TSB (reported by the assay)"),
                linewidth = 1.1) +
      geom_line(aes(y = UCBNAT, colour = "native (4Z,15Z) pigment"),
                linewidth = 0.8) +
      geom_line(aes(y = TCB, colour = "transcutaneous TcB"), linewidth = 0.6) +
      scale_colour_manual(values = c(
        "TSB (reported by the assay)" = "#c0392b",
        "native (4Z,15Z) pigment"     = "#8e44ad",
        "transcutaneous TcB"          = "#27ae60",
        "AAP phototherapy threshold"  = "grey35",
        "AAP exchange threshold"      = "grey10")) +
      labs(x = "postnatal age (days)", y = "bilirubin (mg/dL)", colour = NULL,
           title = "Total serum bilirubin against this infant's own thresholds",
           subtitle = "shaded = lamps on") + THEME
  })

  output$p_wt <- renderPlot({
    ggplot(sim(), aes(time/24, WTOUT)) + geom_line(linewidth = 1) +
      labs(x = "days", y = "weight (kg)") + THEME
  })

  output$t_decision <- renderTable({
    d <- sim(); dt <- d$time[2] - d$time[1]
    tibble(
      quantity = c("peak TSB (mg/dL)", "time of peak (h)",
                   "peak free bilirubin (nM)", "phototherapy hours",
                   "hours above phototherapy threshold",
                   "hours above exchange threshold",
                   "AUC above threshold (mg·h/dL)",
                   "exchange volume given (mL/kg)"),
      value = c(sprintf("%.2f", max(d$TSBOUT)),
                sprintf("%.0f", d$time[which.max(d$TSBOUT)]),
                sprintf("%.1f", max(d$BF)),
                sprintf("%.0f", max(d$PTH)),
                sprintf("%.0f", sum(d$OVERPT)*dt),
                sprintf("%.0f", sum(d$OVERET)*dt),
                sprintf("%.0f", max(d$AUCX)),
                sprintf("%.0f", max(d$ETV))))
  })

  # ------------------------------------------------------------------ TAB 2
  output$p_bf <- renderPlot({
    ggplot(sim(), aes(time/24, BF)) +
      geom_hline(yintercept = 30, linetype = 2, colour = "#c0392b") +
      annotate("text", x = 0.2, y = 32, hjust = 0, size = 3.4,
               label = "~30 nM: injury accrual begins in this model") +
      geom_line(linewidth = 1.1, colour = "#2471a3") +
      labs(x = "postnatal age (days)", y = "free bilirubin (nM)",
           title = "The species that crosses the barrier") + THEME
  })

  output$p_bf_vs_tsb <- renderPlot({
    ggplot(sim(), aes(TSBOUT, BF, colour = time/24)) +
      geom_path(linewidth = 1) +
      geom_hline(yintercept = 30, linetype = 2) +
      scale_colour_viridis_c(name = "days") +
      labs(x = "TSB, what is measured (mg/dL)", y = "Bf, what injures (nM)",
           title = "A convex map, not a proportionality") + THEME
  })

  output$p_ba <- renderPlot({
    ggplot(sim(), aes(time/24)) +
      geom_line(aes(y = BARATIO), linewidth = 1, colour = "#7d3c98") +
      geom_hline(yintercept = 0.8, linetype = 2) +
      labs(x = "days", y = "bilirubin / albumin molar ratio",
           title = "B/A: albumin-free, but NOT affinity-free") + THEME
  })

  output$p_iso <- renderPlot({
    grid <- expand.grid(tsb = seq(2, 30, 0.5),
                        alb = c(2.0, 2.5, 3.0, 3.5, 4.0))
    grid$bf <- free_bilirubin(grid$tsb, grid$alb,
                              fk = input$facid*input$fdisp)
    ggplot(grid, aes(tsb, bf, colour = factor(alb))) +
      geom_line(linewidth = 0.9) +
      geom_hline(yintercept = 30, linetype = 2) +
      coord_cartesian(ylim = c(0, 150)) +
      labs(x = "TSB (mg/dL)", y = "free bilirubin (nM)",
           colour = "albumin\n(g/dL)",
           title = "Iso-risk: the TSB that carries equal Bf") + THEME
  })

  # ------------------------------------------------------------------ TAB 3
  output$p_flux <- renderPlot({
    sim() %>%
      select(time, production = PROD24) %>%
      ggplot(aes(time/24, production)) + geom_line(linewidth = 1) +
      labs(x = "days", y = "bilirubin production (mg/kg/day)",
           title = "The production arm of the flux balance",
           subtitle = "term reference ~7.7 mg/kg/day, roughly twice the adult 3.8") +
      THEME
  })

  output$p_prod <- renderPlot({
    ggplot(sim(), aes(time/24, ETCOC)) + geom_line(linewidth = 1, colour = "#c0392b") +
      geom_hline(yintercept = 1.7, linetype = 2) +
      labs(x = "days", y = "ETCOc (ppm)",
           title = "End-tidal CO reads the production arm directly") + THEME
  })

  output$p_ugt <- renderPlot({
    ggplot(sim(), aes(time/24, UGTPCT)) + geom_line(linewidth = 1, colour = "#1e8449") +
      labs(x = "days", y = "UGT1A1 activity (% of adult)",
           title = "The clearance arm: ontogeny × genotype × induction") + THEME
  })

  # ------------------------------------------------------------------ TAB 4
  output$p_ptgrid <- renderPlot({
    grid <- expand.grid(IRRSET = c(0, 8, 15, 30, 50),
                        FBSASET = c(0.35, 0.60, 0.80, 1.00))
    res <- bind_rows(lapply(seq_len(nrow(grid)), function(i) {
      m <- param(mod, c(pset(), list(AUTOPT = 0, PTSTART = 24, PTSTOP = 1e9,
                                     IRRSET = grid$IRRSET[i],
                                     FBSASET = grid$FBSASET[i])))
      o <- as_tibble(mrgsim(m, end = 96, delta = 2, hmax = 0.25))
      tibble(irr = grid$IRRSET[i], fbsa = grid$FBSASET[i],
             change = 100*(o$TSBOUT[o$time == 48]/o$TSBOUT[o$time == 24] - 1))
    }))
    ggplot(res, aes(irr, change, colour = factor(fbsa))) +
      geom_line(linewidth = 1) + geom_point() +
      labs(x = "spectral irradiance (µW/cm²/nm)",
           y = "change in reported TSB over 24 h (%)",
           colour = "exposed\nBSA fraction",
           title = "Area beats irradiance once the skin saturates") + THEME
  })

  output$p_photoflux <- renderPlot({
    ggplot(sim(), aes(time/24, PHOTOFLX)) + geom_line(linewidth = 1) +
      labs(x = "days", y = "irreversible photochemical removal (mg/h)",
           title = "Lumirubin + E-isomer into bile",
           subtitle = "the only bilirubin exit that does not need UGT1A1") + THEME
  })

  output$p_isomer <- renderPlot({
    sim() %>%
      select(time, `reported TSB` = TSBOUT, `native pigment` = UCBNAT,
             photoisomers = ISOMER) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time/24, value, colour = name)) + geom_line(linewidth = 1) +
      labs(x = "days", y = "mg/dL", colour = NULL,
           title = "Photoisomers are counted as bilirubin by the assay") + THEME
  })

  # ------------------------------------------------------------------ TAB 5
  output$p_et <- renderPlot({
    d <- sim()
    d %>% select(time, TSB = TSBOUT, `free bilirubin (nM)` = BF,
                 `albumin (g/dL)` = ALB, `haemoglobin (g/dL)` = HB,
                 `alloantibody load` = ABMAT) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "postnatal age (h)", y = NULL,
           title = "Exchange transfusion acts on load, capacity and production") +
      THEME
  })

  output$t_et <- renderTable({
    d <- sim()
    if (max(d$ETV) == 0) return(tibble(note = "no exchange transfusion in this run"))
    t0 <- input$tet; t1 <- input$tet + input$etdur
    pre  <- d$TSBOUT[which.min(abs(d$time - t0))]
    post <- d$TSBOUT[which.min(abs(d$time - t1))]
    reb  <- d$TSBOUT[which.min(abs(d$time - (t1 + 6)))]
    bfpre  <- d$BF[which.min(abs(d$time - t0))]
    bfpost <- d$BF[which.min(abs(d$time - t1))]
    tibble(quantity = c("pre-exchange TSB", "immediately post",
                        "removed (%)", "rebound at +6 h (% of pre)",
                        "pre-exchange Bf (nM)", "post-exchange Bf (nM)",
                        "Bf removed (%)"),
           value = sprintf("%.2f", c(pre, post, 100*(1-post/pre),
                                     100*reb/pre, bfpre, bfpost,
                                     100*(1-bfpost/bfpre))))
  })

  # ------------------------------------------------------------------ TAB 6
  output$p_pk1 <- renderPlot({
    sim() %>% select(time, `IVIG central (g)` = IGGC,
                     `stannsoporfin (µg/mL)` = CSNMP) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time/24, value)) + geom_line(linewidth = 1) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      labs(x = "days", y = NULL, title = "Production-arm drugs") + THEME
  })

  output$p_pk2 <- renderPlot({
    sim() %>% select(time, `phenobarbital (mg/L)` = CPB,
                     `UDCA (µmol/L)` = CUDCA) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time/24, value)) + geom_line(linewidth = 1) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      labs(x = "days", y = NULL, title = "Clearance-arm drugs") + THEME
  })

  output$p_pdeff <- renderPlot({
    sim() %>% select(time, `UGT1A1 (% adult)` = UGTPCT,
                     `transgene (fraction)` = TGX,
                     `luminal binder (mg)` = GBIND) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time/24, value)) + geom_line(linewidth = 1) +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "days", y = NULL, title = "Pharmacodynamic states") + THEME
  })

  # ------------------------------------------------------------------ TAB 7
  output$p_brain <- renderPlot({
    ggplot(sim(), aes(time/24, BBR)) +
      geom_hline(yintercept = 35, linetype = 2, colour = "#c0392b") +
      geom_line(linewidth = 1.1, colour = "#943126") +
      labs(x = "days", y = "brain bilirubin (nM-equivalent)",
           title = "Basal-ganglia bilirubin",
           subtitle = "dashed = injury threshold") + THEME
  })

  output$p_inj <- renderPlot({
    sim() %>% select(time, `cumulative injury` = INJ,
                     `P(kernicterus)` = KERN) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time/24, value, colour = name)) + geom_line(linewidth = 1.1) +
      coord_cartesian(ylim = c(0, 1)) +
      labs(x = "days", y = NULL, colour = NULL,
           title = "Irreversible accrual and its logistic read-out") + THEME
  })

  output$p_abr <- renderPlot({
    ggplot(sim(), aes(time/24, ABRD)) + geom_line(linewidth = 1) +
      labs(x = "days", y = "ABR deficit (relative)",
           title = "Auditory brainstem response — a reversible marker") + THEME
  })

  output$t_risk <- renderTable({
    d <- sim()
    tibble(quantity = c("peak brain bilirubin (nM-eq)",
                        "hours above the injury threshold",
                        "final cumulative injury index",
                        "P(kernicterus)",
                        "peak ABR deficit",
                        "AAP neurotoxicity risk factor flag"),
           value = c(sprintf("%.1f", max(d$BBR)),
                     sprintf("%.0f", sum(d$BBR > 35)*(d$time[2]-d$time[1])),
                     sprintf("%.3f", last(d$INJ)),
                     sprintf("%.3f", last(d$KERN)),
                     sprintf("%.3f", max(d$ABRD)),
                     ifelse(pset()$RF > 0, "present", "absent")))
  })

  # ------------------------------------------------------------------ TAB 8
  output$p_etco <- renderPlot({
    ggplot(sim(), aes(time/24, ETCOC)) + geom_line(linewidth = 1) +
      geom_hline(yintercept = c(1.3, 1.7), linetype = 3) +
      labs(x = "days", y = "ETCOc (ppm)",
           title = "End-tidal CO", subtitle = "dotted = normal band") + THEME
  })
  output$p_hb <- renderPlot({
    sim() %>% select(time, `haemoglobin (g/dL)` = HB,
                     `reticulocytes (%)` = RET) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time/24, value)) + geom_line(linewidth = 1) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      labs(x = "days", y = NULL, title = "Erythron") + THEME
  })
  output$p_tcb <- renderPlot({
    ggplot(sim(), aes(TSBOUT, TCB, colour = IRRAD > 0)) + geom_point(size = 0.7) +
      geom_abline(slope = 1, intercept = 0, linetype = 2) +
      labs(x = "TSB (mg/dL)", y = "TcB (mg/dL)", colour = "lamps on",
           title = "TcB dissociates from TSB under the lamps") + THEME
  })
  output$p_stool <- renderPlot({
    sim() %>% select(time, `cumulative output (mg)` = BSTL,
                     `gut unconjugated (mg)` = BGU) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time/24, value)) + geom_line(linewidth = 1) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      labs(x = "days", y = NULL, title = "The enterohepatic shunt") + THEME
  })

  # ------------------------------------------------------------------ TAB 9
  scen_sim <- eventReactive(input$runscen, {
    sel <- input$scen
    if (length(sel) == 0) return(NULL)
    bind_rows(lapply(sel, function(s) run_scenario(scenarios[[s]], end = 336)))
  })

  output$p_scen <- renderPlot({
    d <- scen_sim(); if (is.null(d)) return(NULL)
    ggplot(d, aes(time/24, TSBOUT, colour = scenario)) +
      geom_line(aes(y = THRESHPT), linetype = 2, colour = "grey55",
                show.legend = FALSE) +
      geom_line(linewidth = 1) +
      facet_wrap(~scenario, ncol = 2) +
      labs(x = "postnatal age (days)", y = "TSB (mg/dL)") +
      THEME + theme(legend.position = "none")
  })

  output$t_scen <- renderTable({
    d <- scen_sim(); if (is.null(d)) return(NULL)
    summarise_scenarios(d)
  })

  # ----------------------------------------------------------------- TAB 10
  cn <- eventReactive(input$run, { cn_ceiling(days = 21) })

  output$p_cn <- renderPlot({
    cn() %>%
      pivot_longer(c(TSB_G, TSB_GP, TSB_GPO)) %>%
      mutate(name = recode(name,
                           TSB_G   = "[G] geometry only",
                           TSB_GP  = "[G+P] + production ontogeny",
                           TSB_GPO = "[G+P+O] + skin optics")) %>%
      ggplot(aes(age_y, value, colour = name)) +
      geom_line(linewidth = 1.1) + geom_point() +
      labs(x = "age (years)", y = "phototherapy ceiling TSB (mg/dL)",
           colour = NULL,
           title = "Crigler-Najjar type I on 12 h/day intensive phototherapy",
           subtitle = paste("geometry alone predicts a rising ceiling; adding the",
                            "fall in production per kg cancels it")) + THEME
  })

  output$t_cn <- renderTable({ cn() })

  # ----------------------------------------------------------------- TAB 11
  output$t_cmt <- renderText(paste(names(init(mod)), collapse = " · "))
  output$t_par <- renderPrint({ str(pset()) })
}

shinyApp(ui, server)

# =============================================================================
#  END
# =============================================================================
