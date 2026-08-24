## =============================================================================
##  abm_shiny_app.R
##  Acute Bacterial Meningitis (pneumococcal) — QSP dashboard
##  Acute bacterial meningitis QSP interactive dashboard (12 tabs)
##
##  Run with:
##      shiny::runApp("abm_shiny_app.R")
##  (abm_mrgsolve_model_en.R must be in the same directory)
##
##  This app sets out to show one thing: **the injury flux is a product, and
##  that product is maximal at the first antibiotic dose.**  So the "dexamethasone
##  timing" slider is the app's central control, and every other tab shows where that manipulation propagates.
##
##      injury flux = k_kill(t) · N(t) · Y_lysis · (1 − E_dex(t))
##
##  The second thing it sets out to show: the blood-CSF barrier is a two-way door.  In the
##  "Barrier and penetration" tab, switching dexamethasone on brings inflammation
##  (and oedema·ICP) down and the vancomycin CSF concentration down with it — harmless in a susceptible strain, harmful in a resistant one.
## =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

source("abm_mrgsolve_model_en.R")

theme_abm <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")

## Literature reference ranges (the shaded bands in the diagnostic panel tab)
ref_bands <- list(
  PMN      = c(1000, 5000),
  PROT_CSF = c(100, 500),
  GLC_CSF  = c(0, 40),
  LAC_CSF  = c(3.5, 12),
  QALB     = c(30, 100),
  TNF      = c(100, 1000),
  MMP9     = c(100, 1000),
  ICP      = c(20, 40)
)

## =============================================================================
##  UI
## =============================================================================
ui <- fluidPage(
  titlePanel("Acute Bacterial Meningitis (pneumococcal) — QSP dashboard · 63-ODE model"),
  tags$p(style = "color:#444;",
         HTML("injury flux = <b>k_kill(t) · N(t) · Y_lysis · (1 − E_dex(t))</b> — ",
              "this product is maximal <b>immediately after the first antibiotic dose</b>. ",
              "So the steroid has to be on already <b>before</b> that peak.")),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Patient at presentation"),
      sliderInput("logN0", "CSF bacterial density log₁₀(CFU/mL)", 4, 9, 7, step = 0.25),
      sliderInput("logNb0", "Bacteraemia log₁₀(CFU/mL)", 0, 6, 3, step = 0.5),
      sliderInput("host", "Host defence index", 0.3, 1.0, 1.0, step = 0.05),
      sliderInput("wt", "Body weight (kg)", 40, 120, 70, step = 5),

      hr(), h4("Antibiotics"),
      checkboxInput("use_cef", "Ceftriaxone 2 g q12h", TRUE),
      checkboxInput("cef_ci", "  └ as continuous infusion (4 g/day)", FALSE),
      sliderInput("cef_delay", "Antibiotic dosing delay (h)", 0, 24, 0, step = 0.5),
      selectInput("mic_cef", "Pneumococcal ceftriaxone MIC (mg/L)",
                  c("0.03 (susceptible)" = 0.03, "0.5 (intermediate)" = 0.5,
                    "2 (resistant)" = 2, "4 (highly resistant)" = 4), selected = 0.03),
      checkboxInput("use_van", "Vancomycin 15 mg/kg q6h", FALSE),
      checkboxInput("use_rif", "Rifampicin 600 mg q12h", FALSE),
      sliderInput("rif_lead", "  └ Rifampicin lead time (h)", 0, 6, 0, step = 0.5),

      hr(), h4("Dexamethasone — the central control of this app"),
      checkboxInput("use_dex", "Dexamethasone 0.15 mg/kg q6h", TRUE),
      sliderInput("dex_time", "Dose time (relative to antibiotic, h)", -2, 24, -0.33, step = 0.33),
      sliderInput("dex_days", "Dosing duration (days)", 1, 4, 4, step = 1),

      hr(), h4("Adjunctive therapy"),
      checkboxInput("use_mann", "Mannitol 0.5 g/kg q6h", FALSE),
      checkboxInput("use_gly", "Glycerol 1.5 g/kg q6h", FALSE),
      sliderInput("drain", "EVD drainage (mL/h)", 0, 24, 0, step = 2),
      sliderInput("antipyr", "Antipyretic effect", 0, 1, 0, step = 0.1),
      sliderInput("anticonv", "Anticonvulsant effect", 0, 1, 0, step = 0.1),

      hr(),
      sliderInput("tend", "Simulation duration (h)", 48, 336, 336, step = 24),
      actionButton("go", "Recalculate", class = "btn-primary")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs", type = "tabs",

        ## ------------------------------------------------------------------ 1
        tabPanel("1 · Patient profile",
                 h4("State at presentation and the immediately predicted outcome"),
                 fluidRow(column(6, tableOutput("tbl_profile")),
                          column(6, tableOutput("tbl_endpoints"))),
                 hr(),
                 h4("Treatment timeline"),
                 plotOutput("p_timeline", height = "220px"),
                 helpText("t = 0 is not the onset of infection but arrival at hospital (presentation). ",
                          "Symptoms have already been running for 12–36 hours, so the initial conditions ",
                          "are matched to the reported CSF panel of pneumococcal meningitis.")),

        ## ------------------------------------------------------------------ 2
        tabPanel("2 · Antibiotic PK",
                 h4("Plasma and CSF — the delivery step sets the therapeutic index"),
                 plotOutput("p_pk", height = "440px"),
                 fluidRow(column(6, h5("CSF exposure summary"), tableOutput("tbl_pk")),
                          column(6, h5("Penetration ratio (CSF / total plasma concentration)"),
                                 plotOutput("p_pen", height = "220px"))),
                 helpText("Ceftriaxone protein binding is saturable — the free fraction rises ",
                          "from 0.076 at 30 mg/L to 0.167 at 250 mg/L. ",
                          "Only free drug crosses into the CSF.")),

        ## ------------------------------------------------------------------ 3
        tabPanel("3 · Bacterial kinetics",
                 h4("Free CSF bacteria · adherent subpopulation · bacteraemia"),
                 plotOutput("p_bact", height = "340px"),
                 plotOutput("p_kill", height = "240px"),
                 tableOutput("tbl_ster"),
                 helpText("Sterilisation is judged on free CSF bacteria N_c < 10 CFU/mL (the culture ",
                          "detection limit). The adherent·sequestered subpopulation receives only 35 % of the kill effect, so ",
                          "it always disappears later — the reason a full treatment course is needed ",
                          "is the gap between these two curves.")),

        ## ------------------------------------------------------------------ 4
        tabPanel("4 · Injury flux (the product)",
                 h4("k_kill × N × Y × (1 − E_dex) — the central claim of this model"),
                 plotOutput("p_product", height = "420px"),
                 fluidRow(column(6, plotOutput("p_cargo", height = "260px")),
                          column(6, plotOutput("p_burst_zoom", height = "260px"))),
                 helpText("k_kill jumps 0 → E_max at the first dose and N is maximal at exactly ",
                          "that moment. So the peak of the injury flux lies in the first few ",
                          "hours of treatment. Delay the antibiotic and N has grown exponentially, so ",
                          "the peak grows — delay is bad not because 'there are many bacteria' ",
                          "but because 'there is much to kill'.")),

        ## ------------------------------------------------------------------ 5
        tabPanel("5 · Inflammation",
                 h4("CSF cytokines · neutrophils · proteolysis · oxidative burden"),
                 plotOutput("p_cyto", height = "420px"),
                 plotOutput("p_pmn", height = "260px"),
                 helpText("The shaded bands are the CSF ranges reported in pneumococcal meningitis. ",
                          "The dexamethasone effect compartment (IDEX) is asymmetric — on t½ 1.5 h, ",
                          "off t½ 17 h. That asymmetry is what makes the dose timing matter.")),

        ## ------------------------------------------------------------------ 6
        tabPanel("6 · Barrier and penetration (two doors)",
                 h4("The same door that inflammation opens and the steroid closes"),
                 plotOutput("p_barrier", height = "400px"),
                 h5("Where the sign flips"),
                 tableOutput("tbl_signflip"),
                 helpText("When the barrier Pb opens, hydrophilic antibiotic comes in and at the ",
                          "same time albumin·neutrophils·water come in. Dexamethasone reduces both. ",
                          "Ceftriaxone has 200-fold C/MIC headroom so the net effect is a gain, but ",
                          "vancomycin against a cephalosporin-resistant strain has no headroom, so ",
                          "the same manoeuvre delays sterilisation 2-fold.")),

        ## ------------------------------------------------------------------ 7
        tabPanel("7 · CSF diagnostic panel",
                 h4("The values actually measured — against the literature ranges"),
                 plotOutput("p_panel", height = "520px"),
                 tableOutput("tbl_panel"),
                 helpText("Glucose competition: 4,000/µL neutrophils eat 20 mg/dL/h and ",
                          "10⁷ CFU/mL bacteria eat 2 mg/dL/h. A ten-fold difference. ",
                          "So the low glucose persists after the CSF has turned culture-negative — ",
                          "glucose cannot be an early marker of sterilisation.")),

        ## ------------------------------------------------------------------ 8
        tabPanel("8 · Intracranial pressure and perfusion",
                 h4("CSF dynamics · oedema · CPP · ischaemia"),
                 plotOutput("p_icp", height = "420px"),
                 plotOutput("p_perf", height = "260px"),
                 helpText("An exponential compliance model (Marmarou PVI 25 mL) that at steady ",
                          "state reduces to ICP = P_ss + Q_f · R_out. ",
                          "The potency of osmotherapy is proportional to 1/Pb — the barrier damage ",
                          "that made the oedema also lets the osmotic agent leak.")),

        ## ------------------------------------------------------------------ 9
        tabPanel("9 · Clinical endpoints",
                 h4("Death · hearing · cognition · focal deficit"),
                 fluidRow(column(5, tableOutput("tbl_out")),
                          column(7, plotOutput("p_hazard", height = "300px"))),
                 plotOutput("p_injury", height = "280px"),
                 helpText("The hazard function integrates only the acute physiological derangement and ",
                          "evaluates structural damage once, at the end. Multiply permanent damage by time ",
                          "and even a recovered patient accrues unbounded hazard.")),

        ## ------------------------------------------------------------------ 10
        tabPanel("10 · Dexamethasone timing sweep",
                 h4("When to raise the shield"),
                 actionButton("run_sweep", "Run sweep", class = "btn-warning"),
                 plotOutput("p_sweep", height = "420px"),
                 tableOutput("tbl_sweep"),
                 helpText("The cytokine burst passes its peak 1–4 hours after the first antibiotic ",
                          "dose. The transcriptional effect takes t½ 1.5 hours to switch on, so ",
                          "a dose 4 hours late misses the burst. This is the mechanistic content of ",
                          "the recommendation 'before or with the antibiotic'.")),

        ## ------------------------------------------------------------------ 11
        tabPanel("11 · Scenario comparison",
                 h4("26 treatment scenarios"),
                 actionButton("run_all", "Run all scenarios (takes several minutes)",
                              class = "btn-warning"),
                 tableOutput("tbl_all"),
                 hr(),
                 h5("Virtual cohort of 10 × (DEX on/off)"),
                 actionButton("run_coh", "Run cohort", class = "btn-warning"),
                 tableOutput("tbl_coh"),
                 helpText("Reference: de Gans & van de Beek NEJM 2002 pneumococcal subgroup — ",
                          "death 34 % → 14 %, unfavourable outcome 52 % → 26 %.")),

        ## ------------------------------------------------------------------ 12
        tabPanel("12 · Methods and verification",
                 h4("How this model was verified"),
                 verbatimTextOutput("txt_checks"),
                 hr(),
                 h4("Defects caught by numerical verification"),
                 htmlOutput("txt_bugs"))
      )
    )
  )
)

## =============================================================================
##  Server
## =============================================================================
server <- function(input, output, session) {

  dosing <- reactive({
    wt <- input$wt
    e <- NULL
    add <- function(a, b) if (is.null(a)) b else a + b
    if (input$use_cef) {
      e <- add(e, if (input$cef_ci) ev_cef_ci(input$cef_delay)
                  else ev_cef(input$cef_delay))
    }
    if (input$use_van) {
      e <- add(e, ev(time = input$cef_delay, amt = 15 * wt, cmt = 4,
                     rate = 15 * wt, ii = 6, addl = 55))
    }
    if (input$use_rif) {
      e <- add(e, ev_rif(input$cef_delay - input$rif_lead))
    }
    if (input$use_dex) {
      e <- add(e, ev(time = input$cef_delay + input$dex_time,
                     amt = 0.15 * wt, cmt = 9, rate = 0.15 * wt / 0.25,
                     ii = 6, addl = 4 * input$dex_days - 1))
    }
    if (input$use_mann) {
      e <- add(e, ev(time = 1, amt = 0.5 * wt * 1000, cmt = 12,
                     rate = 0.5 * wt * 1000 / 0.5, ii = 6, addl = 11))
    }
    if (input$use_gly) {
      e <- add(e, ev(time = 0, amt = 1.5 * wt * 1000, cmt = 13, ii = 6, addl = 15))
    }
    e
  })

  sim <- eventReactive(input$go, ignoreNULL = FALSE, {
    withProgress(message = "Integrating...", {
      run_abm(dosing(), end = input$tend, delta = 0.1,
              N0 = 10^input$logN0, NB0 = 10^input$logNb0,
              HOST_DEF = input$host,
              MU_SCALE = 1 + 0.5 * (1 - input$host),
              MIC_cef = as.numeric(input$mic_cef),
              ANTIPYR = input$antipyr, ANTICONV = input$anticonv,
              CSF_DRAIN = input$drain)
    })
  })

  D <- reactive(as.data.frame(sim()))
  E <- reactive(endpoints_abm(sim()))

  long <- function(d, vars) {
    d %>% select(time, all_of(vars)) %>%
      pivot_longer(-time, names_to = "var", values_to = "val")
  }

  band <- function(v) {
    b <- ref_bands[[v]]
    if (is.null(b)) return(NULL)
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = b[1], ymax = b[2],
             fill = "steelblue", alpha = 0.10)
  }

  ## ---------------------------------------------------------------- tab 1
  output$tbl_profile <- renderTable({
    d <- D()[1, ]
    data.frame(
      Variable = c("CSF bacterial density (CFU/mL)", "Bacteraemia (CFU/mL)", "CSF white cells (/µL)",
               "CSF protein (mg/dL)", "CSF glucose (mg/dL)", "CSF glucose ratio",
               "CSF lactate (mmol/L)", "Albumin index Q_alb", "Barrier permeability Pb",
               "ICP (mmHg)", "CPP (mmHg)", "Temperature (°C)", "SOFA"),
      At_presentation = c(sprintf("%.1e", 10^input$logN0), sprintf("%.1e", 10^input$logNb0),
                 sprintf("%.0f", d$PMN), sprintf("%.0f", d$PROT_CSF),
                 sprintf("%.1f", d$GLC_CSF), sprintf("%.2f", d$GLCR),
                 sprintf("%.1f", d$LAC_CSF), sprintf("%.0f", d$QALB),
                 sprintf("%.1f", d$PB), sprintf("%.1f", d$ICP),
                 sprintf("%.1f", d$CPPo), sprintf("%.1f", d$TEMP),
                 sprintf("%.1f", d$SOFA)))
  }, striped = TRUE)

  output$tbl_endpoints <- renderTable({
    e <- E()
    data.frame(
      Endpoint = c("CSF sterilisation (h)", "Adherent subpopulation cleared (h)", "Peak ICP (mmHg)",
                     "Minimum CPP (mmHg)", "Hearing threshold shift (dB)", "Cognition z",
                     "Death probability (%)", "Acute hazard integral"),
      Value = c(ifelse(is.na(e$sterile_h), "not reached", sprintf("%.1f", e$sterile_h)),
             ifelse(is.na(e$adh_clear_h), "not reached", sprintf("%.1f", e$adh_clear_h)),
             sprintf("%.1f", e$peak_ICP), sprintf("%.1f", e$min_CPP),
             sprintf("%.1f", e$hear_dB), sprintf("%.2f", e$cog_z),
             sprintf("%.1f", 100 * e$p_death), sprintf("%.3f", e$haz_acute)))
  }, striped = TRUE)

  output$p_timeline <- renderPlot({
    ev <- dosing()
    if (is.null(ev)) return(NULL)
    df <- as.data.frame(ev@data)
    nm <- c("1" = "Ceftriaxone", "4" = "Vancomycin", "7" = "Rifampicin",
            "9" = "Dexamethasone", "12" = "Mannitol", "13" = "Glycerol")
    df$drug <- nm[as.character(df$cmt)]
    df <- df %>% mutate(n = ifelse(is.na(addl), 1, addl + 1),
                        ii2 = ifelse(is.na(ii), 0, ii))
    ex <- do.call(rbind, lapply(seq_len(nrow(df)), function(i)
      data.frame(drug = df$drug[i],
                 time = df$time[i] + (0:(df$n[i] - 1)) * df$ii2[i])))
    ex <- ex[ex$time <= min(input$tend, 96), ]
    ggplot(ex, aes(time, drug, colour = drug)) +
      geom_point(size = 3, shape = 25, fill = "white") +
      labs(x = "Time (h)", y = NULL, title = "Dose times in the first 96 hours") +
      theme_abm + theme(legend.position = "none")
  })

  ## ---------------------------------------------------------------- tab 2
  output$p_pk <- renderPlot({
    d <- D()
    d$CEF_PL <- d$CEF_C / 8
    d$VAN_PL <- d$VAN_C / 20
    d$RIF_PL <- d$RIF_C / 50
    d$DEX_PL <- d$DEX_C / 70
    p1 <- long(d, c("CEF_PL", "VAN_PL", "RIF_PL", "DEX_PL")) %>%
      ggplot(aes(time, val, colour = var)) + geom_line(linewidth = 0.7) +
      scale_y_log10() + labs(x = NULL, y = "Total plasma concentration (mg/L)", colour = NULL,
                             title = "Plasma") + theme_abm
    p2 <- long(d, c("CEF_CSF", "VAN_CSF", "RIF_CSF", "DEX_CSF")) %>%
      ggplot(aes(time, val, colour = var)) + geom_line(linewidth = 0.7) +
      scale_y_log10() + labs(x = "Time (h)", y = "CSF concentration (mg/L)", colour = NULL,
                             title = "CSF — this is where the bacteria are") + theme_abm
    gridExtra::grid.arrange(p1, p2, ncol = 1)
  })

  output$p_pen <- renderPlot({
    long(D(), c("PEN_CEF", "PEN_VAN")) %>%
      ggplot(aes(time, val, colour = var)) + geom_line(linewidth = 0.8) +
      labs(x = "Time (h)", y = "CSF / total plasma concentration", colour = NULL) +
      theme_abm
  })

  output$tbl_pk <- renderTable({
    e <- E()
    mic <- as.numeric(input$mic_cef)
    data.frame(
      Metric = c("CEF CSF AUC (mg·h/L)", "VAN CSF AUC (mg·h/L)",
               sprintf("CEF T>4×MIC (h)  [4×MIC = %.2f]", 4 * mic),
               "VAN T>4×MIC (h)  [4 mg/L]", "Peak Pb"),
      Value = c(sprintf("%.0f", D()$AUC_CEF[nrow(D())]),
             sprintf("%.0f", e$AUC_van), sprintf("%.1f", e$T_cef_4MIC),
             sprintf("%.1f", e$T_van_4MIC), sprintf("%.1f", e$peak_Pb)))
  }, striped = TRUE)

  ## ---------------------------------------------------------------- tab 3
  output$p_bact <- renderPlot({
    d <- D()
    d$log_NC <- log10(pmax(d$NC, 1e-2))
    d$log_NADH <- log10(pmax(d$NADH, 1e-2))
    d$log_NB <- log10(pmax(d$NB, 1e-2))
    long(d, c("log_NC", "log_NADH", "log_NB")) %>%
      ggplot(aes(time, val, colour = var)) + geom_line(linewidth = 0.8) +
      geom_hline(yintercept = 1, linetype = 2, colour = "grey40") +
      annotate("text", x = Inf, y = 1.3, hjust = 1.05, size = 3,
               label = "culture detection limit 10 CFU/mL") +
      labs(x = "Time (h)", y = "log₁₀ CFU/mL", colour = NULL) + theme_abm
  })

  output$p_kill <- renderPlot({
    long(D(), c("KCEF", "KVAN", "KRIF", "KTOT")) %>%
      ggplot(aes(time, val, colour = var)) + geom_line(linewidth = 0.7) +
      coord_cartesian(xlim = c(0, min(72, input$tend))) +
      labs(x = "Time (h)", y = "Kill rate (1/h)", colour = NULL,
           title = "Emax kill rate — as a function of the C/MIC ratio") + theme_abm
  })

  output$tbl_ster <- renderTable({
    e <- E()
    data.frame(
      Variable = c("Free CSF bacteria sterilised", "Adherent subpopulation cleared", "Gap"),
      Time_h = c(ifelse(is.na(e$sterile_h), "not reached", sprintf("%.1f", e$sterile_h)),
                 ifelse(is.na(e$adh_clear_h), "not reached", sprintf("%.1f", e$adh_clear_h)),
                 ifelse(is.na(e$adh_clear_h) | is.na(e$sterile_h), "—",
                        sprintf("%.1f", e$adh_clear_h - e$sterile_h))))
  })

  ## ---------------------------------------------------------------- tab 4
  output$p_product <- renderPlot({
    d <- D()
    d$factor_kkill <- d$KTOT / max(d$KTOT, 1e-9)
    d$factor_N <- d$NC / max(d$NC, 1e-9)
    d$factor_shield <- 1 - d$IDEX
    d$product <- d$LYSFLUX / max(d$LYSFLUX, 1e-9)
    long(d, c("factor_kkill", "factor_N", "factor_shield", "product")) %>%
      mutate(var = factor(var,
        levels = c("factor_kkill", "factor_N", "factor_shield", "product"),
        labels = c("k_kill (normalised)", "N (normalised)", "1 − E_dex (shield)",
                   "injury flux (normalised)"))) %>%
      ggplot(aes(time, val, colour = var)) +
      geom_line(linewidth = 0.9) +
      coord_cartesian(xlim = c(0, min(48, input$tend))) +
      labs(x = "Time (h)", y = "Normalised value", colour = NULL,
           title = "The three factors of the product, and the product") + theme_abm
  })

  output$p_cargo <- renderPlot({
    long(D(), c("CW", "PLY")) %>%
      ggplot(aes(time, val, colour = var)) + geom_line(linewidth = 0.8) +
      coord_cartesian(xlim = c(0, min(96, input$tend))) +
      labs(x = "Time (h)", y = "Concentration", colour = NULL,
           title = "Released cargo — cell wall and pneumolysin") + theme_abm
  })

  output$p_burst_zoom <- renderPlot({
    d <- D()
    ggplot(d, aes(time, LYSFLUX)) +
      geom_area(fill = "firebrick", alpha = 0.35) +
      geom_line(colour = "firebrick", linewidth = 0.8) +
      coord_cartesian(xlim = c(0, 12)) +
      labs(x = "Time (h)", y = "Lysis flux (CWU/mL/h)",
           title = "The first 12 hours — the peak is here") + theme_abm
  })

  ## ---------------------------------------------------------------- tab 5
  output$p_cyto <- renderPlot({
    d <- D()
    vars <- c("TNF", "IL1", "IL6", "IL10", "CXCL8", "MMP9")
    long(d, vars) %>%
      ggplot(aes(time, val)) + geom_line(colour = "firebrick", linewidth = 0.8) +
      facet_wrap(~var, scales = "free_y", ncol = 3) +
      labs(x = "Time (h)", y = "CSF concentration (pg/mL, MMP9 in ng/mL)") + theme_abm
  })

  output$p_pmn <- renderPlot({
    d <- D()
    p1 <- ggplot(d, aes(time, PMN)) + band("PMN") +
      geom_line(colour = "steelblue4", linewidth = 0.9) +
      labs(x = "Time (h)", y = "CSF neutrophils (/µL)",
           title = "Neutrophils (shading = literature range)") + theme_abm
    p2 <- long(d, c("MG", "ROS", "IDEX")) %>%
      ggplot(aes(time, val, colour = var)) + geom_line(linewidth = 0.8) +
      labs(x = "Time (h)", y = NULL, colour = NULL,
           title = "Macrophage activation · oxidative burden · steroid effect") + theme_abm
    gridExtra::grid.arrange(p1, p2, ncol = 2)
  })

  ## ---------------------------------------------------------------- tab 6
  output$p_barrier <- renderPlot({
    d <- D()
    p1 <- ggplot(d, aes(time, PB)) +
      geom_line(colour = "darkcyan", linewidth = 1) +
      labs(x = NULL, y = "Barrier permeability Pb", title = "The door") + theme_abm
    p2 <- ggplot(d, aes(time, QALB)) + band("QALB") +
      geom_line(colour = "darkcyan", linewidth = 0.9) +
      labs(x = NULL, y = "Q_alb (×10⁻³)", title = "The door as measured — albumin index") +
      theme_abm
    p3 <- ggplot(d, aes(time, VAN_CSF)) +
      geom_line(colour = "purple4", linewidth = 0.9) +
      geom_hline(yintercept = 4, linetype = 2) +
      annotate("text", x = Inf, y = 4.4, hjust = 1.05, size = 3,
               label = "4×MIC (bactericidal threshold)") +
      labs(x = "Time (h)", y = "Vancomycin CSF (mg/L)",
           title = "Close the door and this dries up first") + theme_abm
    p4 <- ggplot(d, aes(time, VBR)) +
      geom_line(colour = "sienna", linewidth = 0.9) +
      labs(x = "Time (h)", y = "Excess brain water (mL)",
           title = "Close the door and this falls too") + theme_abm
    gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2)
  })

  output$tbl_signflip <- renderTable({
    withProgress(message = "Comparing four arms...", {
      as.data.frame(claim3_signflip()) %>%
        transmute(Arm = arm,
                  Sterilised_h = round(sterile_h, 1),
                  VAN_CSF_AUC = round(AUC_van, 0),
                  Peak_Pb = round(peak_Pb, 1),
                  Peak_ICP = round(peak_ICP, 1),
                  Hearing_dB = round(hear_dB, 1),
                  Death_pct = round(100 * p_death, 1))
    })
  }, striped = TRUE)

  ## ---------------------------------------------------------------- tab 7
  output$p_panel <- renderPlot({
    d <- D()
    mk <- function(v, ylab, ttl) {
      g <- ggplot(d, aes(time, .data[[v]]))
      b <- band(v); if (!is.null(b)) g <- g + b
      g + geom_line(colour = "grey20", linewidth = 0.9) +
        labs(x = NULL, y = ylab, title = ttl) + theme_abm
    }
    gridExtra::grid.arrange(
      mk("PMN", "/µL", "CSF white cells"),
      mk("PROT_CSF", "mg/dL", "CSF protein"),
      mk("GLC_CSF", "mg/dL", "CSF glucose"),
      mk("LAC_CSF", "mmol/L", "CSF lactate"),
      mk("QALB", "×10⁻³", "Albumin index"),
      ggplot(d, aes(time, GLCR)) +
        geom_hline(yintercept = 0.4, linetype = 2, colour = "firebrick") +
        geom_line(colour = "grey20", linewidth = 0.9) +
        labs(x = NULL, y = "Ratio", title = "CSF/plasma glucose ratio (diagnostic threshold 0.4)") +
        theme_abm,
      ncol = 3)
  })

  output$tbl_panel <- renderTable({
    e <- E()
    data.frame(
      Metric = c("Peak white cells (/µL)", "Peak protein (mg/dL)", "Minimum glucose (mg/dL)",
               "Peak lactate (mmol/L)", "Peak Q_alb"),
      Model = c(sprintf("%.0f", e$peak_PMN), sprintf("%.0f", e$peak_Prot),
               sprintf("%.1f", e$min_Glc), sprintf("%.1f", e$peak_Lac),
               sprintf("%.0f", e$peak_Qalb)),
      Literature_range = c("1,000–5,000", "100–500", "<40 (ratio <0.4)", ">3.5 (commonly 6–12)",
                   "30–100+"))
  }, striped = TRUE)

  ## ---------------------------------------------------------------- tab 8
  output$p_icp <- renderPlot({
    d <- D()
    p1 <- ggplot(d, aes(time, ICP)) + band("ICP") +
      geom_line(colour = "darkred", linewidth = 1) +
      geom_hline(yintercept = 25, linetype = 2) +
      labs(x = NULL, y = "ICP (mmHg)", title = "Intracranial pressure") + theme_abm
    p2 <- ggplot(d, aes(time, R_OUT)) +
      geom_line(colour = "chocolate4", linewidth = 0.9) +
      labs(x = NULL, y = "mmHg/(mL/h)", title = "CSF outflow resistance") + theme_abm
    p3 <- ggplot(d, aes(time, VBR)) +
      geom_line(colour = "sienna", linewidth = 0.9) +
      labs(x = "Time (h)", y = "mL", title = "Excess brain water") + theme_abm
    p4 <- long(d, c("MAP", "CPPo")) %>%
      ggplot(aes(time, val, colour = var)) + geom_line(linewidth = 0.9) +
      geom_hline(yintercept = 60, linetype = 2, colour = "firebrick") +
      labs(x = "Time (h)", y = "mmHg", colour = NULL,
           title = "MAP and CPP (dashed = CPP target 60)") + theme_abm
    gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2)
  })

  output$p_perf <- renderPlot({
    d <- D()
    p1 <- ggplot(d, aes(time, CBFo)) +
      geom_hline(yintercept = 27.5, linetype = 2, colour = "firebrick") +
      geom_line(colour = "navy", linewidth = 0.9) +
      annotate("text", x = Inf, y = 29, hjust = 1.05, size = 3,
               label = "ischaemia threshold") +
      labs(x = "Time (h)", y = "CBF (mL/100 g/min)", title = "Cerebral blood flow") + theme_abm
    p2 <- long(d, c("AUTOR", "ISCHo")) %>%
      ggplot(aes(time, val, colour = var)) + geom_line(linewidth = 0.9) +
      labs(x = "Time (h)", y = NULL, colour = NULL,
           title = "Autoregulation preservation and ischaemia index") + theme_abm
    gridExtra::grid.arrange(p1, p2, ncol = 2)
  })

  ## ---------------------------------------------------------------- tab 9
  output$tbl_out <- renderTable({
    e <- E()
    d <- D()[nrow(D()), ]
    data.frame(
      Endpoint = c("Death probability (%)", "Hearing threshold shift (dB)", "Any hearing loss (>25 dB)",
                     "Severe hearing loss (>60 dB)", "Cognition z", "Focal deficit probability",
                     "Cortical surviving fraction", "Dentate gyrus surviving fraction", "Hair cell surviving fraction",
                     "Labyrinthitis ossificans"),
      Value = c(sprintf("%.1f", 100 * e$p_death), sprintf("%.1f", e$hear_dB),
             ifelse(e$hear_dB > 25, "Yes", "No"),
             ifelse(e$hear_dB > 60, "Yes", "No"),
             sprintf("%.2f", e$cog_z), sprintf("%.2f", d$FOCAL),
             sprintf("%.3f", d$NCORT), sprintf("%.3f", d$NDG),
             sprintf("%.3f", d$HC), sprintf("%.3f", d$OSS)))
  }, striped = TRUE)

  output$p_hazard <- renderPlot({
    d <- D()[nrow(D()), ]
    p <- as.list(param(mod))
    df <- data.frame(
      Component = c("ICP", "CPP", "SOFA", "Ischaemia", "Bacteraemia", "CSF infection", "Baseline", "Structural damage"),
      Contribution = c(p$h_icp * d$I_ICP, p$h_cpp * d$I_CPP, p$h_sofa * d$I_SOFA,
               p$h_isch * d$I_ISCH, p$h_bact * d$I_BACT, p$h_ncsf * d$I_NCSF,
               p$h0 * d$time, p$h_cort_final * (1 - d$NCORT)))
    ggplot(df, aes(reorder(Component, Contribution), Contribution)) +
      geom_col(fill = "firebrick", alpha = 0.8) + coord_flip() +
      labs(x = NULL, y = "Cumulative hazard contribution", title = "What actually kills") +
      theme_abm
  })

  output$p_injury <- renderPlot({
    long(D(), c("NCORT", "NDG", "HC", "OSS")) %>%
      mutate(var = factor(var, levels = c("NCORT", "NDG", "HC", "OSS"),
                          labels = c("Cortical neurones", "Hippocampal dentate gyrus", "Cochlear hair cells",
                                     "Labyrinthitis ossificans"))) %>%
      ggplot(aes(time, val, colour = var)) + geom_line(linewidth = 0.9) +
      labs(x = "Time (h)", y = "Surviving fraction (for ossificans, progression fraction)", colour = NULL) +
      theme_abm
  })

  ## ---------------------------------------------------------------- tab 10
  sweep <- eventReactive(input$run_sweep, {
    withProgress(message = "Timing sweep...", claim2_timing())
  })
  output$p_sweep <- renderPlot({
    s <- sweep()
    s2 <- s %>% pivot_longer(-dex_time_h)
    ggplot(s2, aes(dex_time_h, value)) +
      geom_line(colour = "purple4", linewidth = 0.9) +
      geom_point(colour = "purple4") +
      geom_vline(xintercept = 0, linetype = 2) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "Dexamethasone dose time (relative to antibiotic, h)", y = NULL) + theme_abm
  })
  output$tbl_sweep <- renderTable(sweep(), striped = TRUE, digits = 3)

  ## ---------------------------------------------------------------- tab 11
  allsc <- eventReactive(input$run_all, {
    withProgress(message = "26 scenarios...", run_all_scenarios())
  })
  output$tbl_all <- renderTable({
    allsc() %>%
      transmute(id, Scenario = label,
                Sterilised_h = round(sterile_h, 1), Adherent_cleared_h = round(adh_clear_h, 1),
                Peak_ICP = round(peak_ICP, 1), Min_CPP = round(min_CPP, 1),
                Min_glucose = round(min_Glc, 1), Peak_Lac = round(peak_Lac, 1),
                Hearing_dB = round(hear_dB, 1), Cognition_z = round(cog_z, 2),
                Death_pct = round(100 * p_death, 1))
  }, striped = TRUE)

  coh <- eventReactive(input$run_coh, {
    withProgress(message = "20 cohort integrations...", run_cohort())
  })
  output$tbl_coh <- renderTable({
    c0 <- coh()
    rbind(
      c0 %>% transmute(Patient = as.character(id), Bact_load = round(logN0, 1),
                       Delay_h = delay_h, Host_defence = host_def,
                       Death_no_DEX = round(100 * death_no_dex, 1),
                       Death_DEX = round(100 * death_dex, 1),
                       Hearing_no_DEX = round(dB_no_dex, 1),
                       Hearing_DEX = round(dB_dex, 1)),
      data.frame(Patient = "Mean", Bact_load = NA, Delay_h = NA, Host_defence = NA,
                 Death_no_DEX = round(100 * mean(c0$death_no_dex), 1),
                 Death_DEX = round(100 * mean(c0$death_dex), 1),
                 Hearing_no_DEX = round(mean(c0$dB_no_dex), 1),
                 Hearing_DEX = round(mean(c0$dB_dex), 1)))
  }, striped = TRUE)

  ## ---------------------------------------------------------------- tab 12
  output$txt_checks <- renderPrint(structural_checks())

  output$txt_bugs <- renderUI({
    HTML(paste0(
      "<p>Because the equations of this model were written in an environment with no R runtime, ",
      "they were first independently implemented in a dependency-free Python RK4 (<code>abm_reference_python_en.py</code>) ",
      "and actually integrated to verify them. That process exposed <b>32 real defects</b>",
      ". The full list is in that file's header, [F1]–[F32], and ",
      "the complete run log is <code>abm_reference_output_en.txt</code>. Carrying over only ",
      "the ones that changed the structure:</p><ul>",
      "<li><b>F1</b> the glucose consumption term had no substrate dependence, so consumption ",
      "continued even after CSF glucose had reached 0 → lactate 243 mmol/L (measured 6–12). Lactate comes only ",
      "from glucose.</li>",
      "<li><b>F14</b> blood bacteria had no carrying capacity, so N_b diverged to 2×10¹⁰ CFU/mL and ",
      "that back-contaminated the CSF.</li>",
      "<li><b>F15</b> the osmolality conversion used a molecular weight of 18.2 instead of 182, so mannitol 0.5 g/kg ",
      "raised plasma osmolality by 113 mOsm/kg (actually 11.3). A 10-fold unit error.</li>",
      "<li><b>F16</b> the hazard function integrated permanent damage across all 14 days, so even a recovered patient ",
      "reached a death probability of 99 %, and 22 of the 26 scenarios were indistinguishable.</li>",
      "<li><b>F17</b> the autoregulation formula made autoregulation <i>harmful</i> — ",
      "at CPP 38 a patient with autoregulation intact was more ischaemic than a pressure-passive one.</li>",
      "<li><b>F28</b> if the ICP·CPP hazard terms are linear, 'CPP 45 for 100 h' and 'CPP 5 for 10 h' ",
      "carry the same hazard. Mild sustained derangement dominated the hazard function.</li>",
      "<li><b>F29</b> writing the bacterial extinction floor inside the derivative block makes it a discontinuous switch that ",
      "stalls the integrator. The mrgsolve version leaves the state untouched and judges ",
      "sterilisation solely by the reporting threshold.</li>",
      "<li><b>F30</b> instead of matching the hazard coefficients by eye, the per-component integrals were added ",
      "as state variables and the trial targets (pneumococcal death 34 %/14 %) were solved <i>arithmetically</i>. ",
      "The solution was a single scale factor of 0.30 on all coefficients.</li>",
      "</ul>"))
  })
}

shinyApp(ui, server)
