## =====================================================================
##  PATENT DUCTUS ARTERIOSUS OF PREMATURITY -- QSP EXPLORER (Shiny)
## =====================================================================
##
##  A dashboard over pda_mrgsolve_model.R.  The app is organised around
##  the question the field actually argues about, which is NOT "does the
##  drug close the duct" (it does) but "does closing it help" (mostly it
##  has not).  So the tabs deliberately separate three things that the
##  clinical literature keeps fusing:
##
##    * DUCTAL RESPONSE      -- tone, diameter, closure          (tab 3)
##    * HAEMODYNAMIC BURDEN  -- integral of significant shunt     (tab 4)
##    * OUTCOME              -- death / BPD / NEC / IVH           (tab 5)
##
##  Tab 8 is the one to read first: it shows the untreated duct closing on
##  its own, which is the competing process that caps every treatment
##  effect and which the negative trials were really measuring.
##
##  Run:  shiny::runApp("pda_shiny_app.R")
##  Requires: shiny, mrgsolve, dplyr, ggplot2, tidyr, DT
## =====================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)

have_DT <- requireNamespace("DT", quietly = TRUE)

## ---------------------------------------------------------------------
## Load the model.  Sourcing pda_mrgsolve_model.R defines `mod`, the
## rx_* regimen builders, SCENARIOS and the helper functions, so the app
## and the batch model can never diverge.
## ---------------------------------------------------------------------
.here <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) ".")
if (is.null(.here) || !nzchar(.here)) .here <- "."
model_file <- file.path(.here, "pda_mrgsolve_model.R")
if (!file.exists(model_file)) model_file <- "pda_mrgsolve_model.R"
source(model_file, local = FALSE, chdir = TRUE)

THEME <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "#eceff1", colour = NA),
        legend.position = "bottom")

PAL <- c("ibuprofen" = "#00796b", "indomethacin" = "#00838f",
         "acetaminophen" = "#558b2f", "expectant" = "#90a4ae")

## =====================================================================
##  UI
## =====================================================================
ui <- fluidPage(
  titlePanel("Patent Ductus Arteriosus of Prematurity (PDA) · QSP simulator"),
  tags$p(style = "color:#555; margin-top:-8px;",
    HTML("A model for education and research. It must not be used for clinical decision-making. &nbsp;|&nbsp;
          Educational model — not for clinical use.")),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Patient"),
      sliderInput("GA", "Gestational age (wk)", 23, 40, 26, step = 0.5),
      sliderInput("BW", "Birth weight (kg)", 0.4, 3.5, 0.80, step = 0.02),
      sliderInput("PAO2", "PaO₂ (mmHg)", 30, 90, 55, step = 1),
      checkboxInput("ANTESTER", "Antenatal steroids completed", TRUE),
      checkboxInput("SEPSIS",
        "Chorioamnionitis / early-onset sepsis — raised ductal wall peroxide", FALSE),
      checkboxInput("HCORT", "Concomitant early hydrocortisone", FALSE),

      hr(),
      h4("Treatment"),
      selectInput("drug", "Drug", c("No treatment (expectant)" = "none",
                                    "Ibuprofen IV" = "ibu",
                                    "Indomethacin IV" = "ind",
                                    "Acetaminophen IV" = "apap",
                                    "Ibuprofen + acetaminophen" = "combo")),
      sliderInput("start_d", "Treatment start (postnatal day)", 0.25, 21, 2, step = 0.25),
      conditionalPanel("input.drug == 'ibu' || input.drug == 'combo'",
        radioButtons("ibu_mode", "Ibuprofen regimen",
          c("Standard 10-5-5 mg/kg q24h" = "std",
            "High dose 20-10-10 mg/kg q24h" = "high",
            "Continuous infusion" = "inf"), selected = "std")),
      conditionalPanel("input.drug == 'ind'",
        sliderInput("ind_dose", "Indomethacin loading dose (mg/kg)",
                    0.05, 0.30, 0.20, step = 0.05),
        sliderInput("ind_n", "Number of doses", 1, 6, 3, step = 1)),
      conditionalPanel("input.drug == 'apap' || input.drug == 'combo'",
        sliderInput("apap_dose", "Acetaminophen dose per administration (mg/kg)",
                    7.5, 25, 15, step = 2.5),
        sliderInput("apap_days", "Duration of dosing (days)", 1, 7, 3, step = 1)),
      checkboxInput("second_course", "Second course on reopening (d+5)", FALSE),

      hr(),
      sliderInput("horizon", "Observation horizon (days)", 20, 120, 90, step = 5),
      actionButton("go", "Run simulation", class = "btn-primary",
                   icon = icon("play")),
      hr(),
      htmlOutput("verdict")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs", type = "tabs",

        ## ---------------------------------------------------------- 1
        tabPanel("① Patient profile",
          br(),
          fluidRow(
            column(6, h4("What is physically possible at this gestational age"),
                   verbatimTextOutput("profile_txt")),
            column(6, h4("Occlusion ceiling"),
                   plotOutput("plot_ceiling", height = 300))),
          hr(),
          HTML("<p><b>The core structure:</b> the maximum attainable constriction (TMAXGA) is a
            function of gestational age, and if the residual lumen it leaves is larger than the
            closure criterion (0.30 mm) then <b>no drug exposure whatsoever will close that duct.</b>
            At 24 weeks more than 0.4 mm remains even at the ceiling. The steep gradient of
            treatment success with gestational age that is seen clinically is <i>derived</i> from
            this structure; it was not <i>fitted</i> to it.</p>")),

        ## ---------------------------------------------------------- 2
        tabPanel("② Drug PK · COX inhibition",
          br(),
          fluidRow(column(6, plotOutput("plot_pk", height = 300)),
                   column(6, plotOutput("plot_cox", height = 300))),
          hr(),
          fluidRow(column(6, h4("Exposure summary"), tableOutput("tbl_pk")),
                   column(6, h4("Two enzyme sites"),
                     HTML("<p>Ibuprofen and indomethacin compete with the substrate at the
                     <b>arachidonic acid channel</b>, while acetaminophen reduces the
                     ferryl-protoporphyrin radical at the physically separate
                     <b>peroxidase site</b>. Acetaminophen therefore competes not with the
                     substrate but <b>with peroxide</b>, and as ductal wall peroxide tone
                     rises its apparent IC₅₀ rises too. Tick chorioamnionitis in the
                     sidebar and compare the two drugs — this is the falsifiable
                     prediction of this model.</p>")))),

        ## ---------------------------------------------------------- 3
        tabPanel("③ Ductal response (tone · lumen · remodelling)",
          br(),
          fluidRow(column(6, plotOutput("plot_duct", height = 300)),
                   column(6, plotOutput("plot_wall", height = 300))),
          hr(),
          fluidRow(column(12, plotOutput("plot_drive", height = 260))),
          HTML("<p><b>Constriction is not closure.</b> Permanent closure requires hypoxia of the
            ductal <i>wall</i> — when the lumen closes, oxygen diffusion from the lumen is cut
            off, the media becomes hypoxic, and HIF-1α/VEGF/TGF-β1 build the neointimal
            cushion. In a preterm infant the ductal wall is <i>thin</i>, so oxygen still reaches
            it by diffusion even during constriction, the remodelling signal never switches on,
            and once the drug is cleared the duct <b>reopens</b>.</p>")),

        ## ---------------------------------------------------------- 4
        tabPanel("④ Haemodynamics · shunt burden",
          br(),
          fluidRow(column(6, plotOutput("plot_shunt", height = 300)),
                   column(6, plotOutput("plot_press", height = 300))),
          hr(),
          fluidRow(column(6, plotOutput("plot_organ", height = 280)),
                   column(6, plotOutput("plot_burden", height = 280))),
          HTML("<p>Immediately after birth pulmonary vascular resistance (PVR) is high, so even a
            large duct shunts little. As PVR falls over days 2–5 the shunt is <b>revealed</b> —
            this is why the 'large PDA on day 3' appears, and why prophylaxis within 24 hours
            ends up treating a duct that is not yet a problem.
            Because resistance varies inversely with the <b>fourth power</b> of the lumen,
            1.5 mm is not an arbitrary cut-off but a real threshold.</p>")),

        ## ---------------------------------------------------------- 5
        tabPanel("⑤ Clinical endpoints",
          br(),
          fluidRow(column(7, plotOutput("plot_outcome", height = 320)),
                   column(5, h4("At 36 weeks corrected age"),
                          tableOutput("tbl_outcome"))),
          hr(),
          HTML("<p>What drives outcome risk is not 'closed by day 7' but the
            <b>PDA burden</b> — the time integral of significant shunt. This is how the model
            reproduces the Baby-OSCAR result in which the duct closed better and yet death/BPD
            did not improve (69.4% vs 63.5%, if anything worse).
            Gestational age itself dominates BPD risk, and the burden component a drug can
            reduce is small beside it.</p>")),

        ## ---------------------------------------------------------- 6
        tabPanel("⑥ Scenario comparison",
          br(),
          fluidRow(column(4, actionButton("run_scen", "Run all 16 scenarios",
                                          class = "btn-warning")),
                   column(8, helpText("This takes a little while to compute."))),
          br(),
          if (have_DT) DT::dataTableOutput("tbl_scen") else tableOutput("tbl_scen_plain"),
          hr(),
          plotOutput("plot_scen", height = 340)),

        ## ---------------------------------------------------------- 7
        tabPanel("⑦ Biomarkers · safety",
          br(),
          fluidRow(column(6, plotOutput("plot_renal", height = 300)),
                   column(6, plotOutput("plot_safety", height = 300))),
          hr(),
          fluidRow(column(12, tableOutput("tbl_safety"))),
          HTML("<p><b>Efficacy and toxicity are the same molecular event in different organs.</b>
            COX appears four times in this model — duct, kidney, intestinal mucosa, platelets.
            What separates the two NSAIDs is not COX potency but the
            <b>non-COX vasoconstriction</b> that indomethacin has and ibuprofen does not
            (cerebral blood flow −20~40% vs almost no change).</p>")),

        ## ---------------------------------------------------------- 8
        tabPanel("⑧ Spontaneous closure (competing risk)",
          br(),
          fluidRow(column(4, actionButton("run_spont", "Run gestational-age sweep",
                                          class = "btn-warning")),
                   column(8, helpText("Simulates gestational ages 24–38 weeks without treatment."))),
          br(),
          fluidRow(column(6, plotOutput("plot_spont", height = 320)),
                   column(6, tableOutput("tbl_spont"))),
          hr(),
          HTML("<p>Semberova (2017) reported that under a non-intervention policy 94–98% of
            infants below 27 weeks closed spontaneously before discharge, with a median of
            71 days below 26 weeks. If almost every duct closes eventually, then <b>the drug
            effect on a fixed-time closure endpoint is bounded above by the fraction still
            unclosed at that time.</b> It is this ceiling — not drug potency — that the
            negative trials were actually measuring.</p>")),

        ## ---------------------------------------------------------- 9
        tabPanel("⑨ Trial reproduction",
          br(),
          h4("What the model was fitted to and what it predicted"),
          tableOutput("tbl_trials"),
          hr(),
          HTML("<p><b>Fitted:</b> the spontaneous closure times of 8 gestational ages →
            NET50 · NETW · TAUSYN0 · KTAUSYN · KINVGA. Seven-day closure rates by drug →
            Kᵢ·IC₅₀. Both arms of the Baby-OSCAR primary outcome → 3 risk coefficients.<br>
            <b>Predicted (no parameter adjusted after that point):</b>
            the reopening rate and its gestational-age dependence · loss of acetaminophen
            potency when peroxide rises · the need for high doses in late rescue therapy ·
            the BeNeDuctus composite outcome · the IVH reduction without outcome benefit in
            TIPP · separation of renal/cerebral/intestinal flow between the two NSAIDs · the whole targeted-treatment strategy.</p>"))
      )
    )
  )
)

## =====================================================================
##  SERVER
## =====================================================================
server <- function(input, output, session) {

  build_rx <- reactive({
    s <- input$start_d * 24
    e <- NULL
    if (input$drug %in% c("ibu", "combo")) {
      e <- if (input$ibu_mode == "inf") rx_ibuprofen_inf(s)
           else rx_ibuprofen(s, high = (input$ibu_mode == "high"))
    }
    if (input$drug == "ind")
      e <- rx_indomethacin(s, dose = input$ind_dose,
                           maint = min(input$ind_dose, 0.1),
                           n = input$ind_n)
    if (input$drug %in% c("apap", "combo")) {
      a <- rx_acetaminophen(s, dose = input$apap_dose, days = input$apap_days)
      e <- if (is.null(e)) a else ev_bind(e, a)
    }
    if (!is.null(e) && isTRUE(input$second_course)) {
      s2 <- s + 5 * 24
      e2 <- if (input$drug == "ind")
              rx_indomethacin(s2, dose = input$ind_dose,
                              maint = min(input$ind_dose, 0.1), n = input$ind_n)
            else if (input$drug == "apap")
              rx_acetaminophen(s2, dose = input$apap_dose, days = input$apap_days)
            else rx_ibuprofen(s2, high = TRUE)
      e <- ev_bind(e, e2)
    }
    e
  })

  sim <- eventReactive(input$go, ignoreNULL = FALSE, {
    withProgress(message = "Simulating...", value = 0.4, {
      as.data.frame(
        sim_pda(doses = build_rx(), days = input$horizon, delta = 1,
                GA = input$GA, BW = input$BW, PAO2 = input$PAO2,
                SEPSIS = as.numeric(input$SEPSIS),
                HCORT = as.numeric(input$HCORT),
                ANTESTER = as.numeric(input$ANTESTER)))
    })
  })

  ## reference: identical patient, no drug -- every "effect" is a difference
  sim_ref <- eventReactive(input$go, ignoreNULL = FALSE, {
    as.data.frame(
      sim_pda(doses = NULL, days = input$horizon, delta = 1,
              GA = input$GA, BW = input$BW, PAO2 = input$PAO2,
              SEPSIS = as.numeric(input$SEPSIS),
              HCORT = as.numeric(input$HCORT),
              ANTESTER = as.numeric(input$ANTESTER)))
  })

  cd <- function(d) { i <- which(d$DDUCT < 0.30); if (length(i)) d$DAY[i[1]] else NA }

  ## -------------------------------------------------------------- verdict
  output$verdict <- renderUI({
    d <- sim(); r <- sim_ref()
    c1 <- cd(d); c0 <- cd(r)
    reop <- { i <- which(d$DDUCT < 0.30)
              if (length(i)) any(d$DDUCT[i[1]:nrow(d)] > 0.48) else FALSE }
    b1 <- d$PDABUR[nrow(d)]; b0 <- r$PDABUR[nrow(r)]
    p1 <- 100 * d$PCOMP[nrow(d)]; p0 <- 100 * r$PCOMP[nrow(r)]
    HTML(sprintf(
      "<div style='font-size:12px;line-height:1.5'>
       <b>Closure day</b>: %s (untreated %s)<br>
       <b>Acceleration</b>: %s days<br>
       <b>Reopening</b>: %s<br>
       <b>PDA burden</b>: %.0f (untreated %.0f, %+.0f%%)<br>
       <b>Death/BPD</b>: %.1f%% (untreated %.1f%%)</div>",
      ifelse(is.na(c1), "not closed", sprintf("d%.1f", c1)),
      ifelse(is.na(c0), "not closed", sprintf("d%.1f", c0)),
      ifelse(is.na(c1) || is.na(c0), "—", sprintf("%.1f", c0 - c1)),
      ifelse(reop, "<span style='color:#c62828'>yes</span>", "no"),
      b1, b0, ifelse(b0 > 0, 100 * (b1 - b0) / b0, 0), p1, p0))
  })

  ## -------------------------------------------------------------- tab 1
  output$profile_txt <- renderPrint({
    ga <- input$GA
    tmax <- 1 / (1 + exp(-(ga - 21.5) * 0.55))
    dmax <- max(0.6, -1.35 + 0.135 * ga)
    dmin <- dmax * (1 - 0.965 * tmax)
    p50  <- 30 * exp(0.130 * (28 - ga))
    wall <- max(0.05, -1.20 + 0.075 * ga)
    cat(sprintf("Gestational age                         %.1f wk\n", ga))
    cat(sprintf("Unconstricted lumen d_max               %.2f mm\n", dmax))
    cat(sprintf("Maximum attainable constriction TMAXGA  %.3f\n", tmax))
    cat(sprintf("Residual lumen at the ceiling           %.3f mm\n", dmin))
    cat(sprintf("Closure criterion                       %.2f mm\n", 0.30))
    cat(sprintf("→ closable by drug?                     %s\n",
                ifelse(dmin < 0.30, "yes", "no (structural limit)")))
    cat(sprintf("\nO2 sensitivity P50                      %.1f mmHg  (PaO2 %.0f)\n",
                p50, input$PAO2))
    cat(sprintf("Ductal wall thickness index             %.2f\n", wall))
    cat(sprintf("→ wall hypoxia during constriction?     %s\n",
                ifelse(wall > 0.75, "yes → permanent remodelling", "weak → risk of reopening")))
  })

  output$plot_ceiling <- renderPlot({
    ga <- seq(23, 40, 0.1)
    tmax <- 1 / (1 + exp(-(ga - 21.5) * 0.55))
    dmax <- pmax(0.6, -1.35 + 0.135 * ga)
    df <- data.frame(ga, resid = dmax * (1 - 0.965 * tmax))
    ggplot(df, aes(ga, resid)) +
      geom_hline(yintercept = 0.30, linetype = 2, colour = "#c62828") +
      geom_line(linewidth = 1.1, colour = "#3949ab") +
      geom_vline(xintercept = input$GA, colour = "#00796b", linetype = 3) +
      annotate("text", x = 36, y = 0.33, label = "closure criterion 0.30 mm",
               colour = "#c62828", size = 3.4) +
      labs(x = "Gestational age (wk)", y = "Residual lumen at maximal constriction (mm)",
           title = "Structural closure limit") + THEME
  })

  ## -------------------------------------------------------------- tab 2
  output$plot_pk <- renderPlot({
    d <- sim()
    df <- d %>% select(DAY, CIBU, CIND, CAPAP) %>%
      pivot_longer(-DAY) %>%
      mutate(name = recode(name, CIBU = "ibuprofen", CIND = "indomethacin",
                           CAPAP = "acetaminophen")) %>%
      filter(value > 1e-6)
    if (!nrow(df)) return(ggplot() + annotate("text", 1, 1, label = "no treatment") +
                          THEME + theme_void())
    ggplot(df, aes(DAY, value, colour = name)) +
      geom_line(linewidth = 0.9) + scale_colour_manual(values = PAL) +
      coord_cartesian(xlim = c(0, min(30, input$horizon))) +
      labs(x = "Postnatal day", y = "Total plasma concentration (mg/L)", colour = NULL,
           title = "Drug concentrations") + THEME
  })

  output$plot_cox <- renderPlot({
    d <- sim()
    df <- d %>% select(DAY, ICOXD, ICHAN, IPEROX, ICOXK, ICOXG, ICOXPLT) %>%
      pivot_longer(-DAY) %>%
      mutate(name = factor(name,
        levels = c("ICOXD", "ICHAN", "IPEROX", "ICOXK", "ICOXG", "ICOXPLT"),
        labels = c("Duct (net)", "Channel component", "Peroxidase component",
                   "Kidney", "Intestinal mucosa", "Platelets")))
    ggplot(df, aes(DAY, 100 * value, colour = name)) +
      geom_line(linewidth = 0.85) +
      coord_cartesian(xlim = c(0, min(30, input$horizon))) +
      labs(x = "Postnatal day", y = "COX inhibition (%)", colour = NULL,
           title = "One enzyme, four organs") + THEME
  })

  output$tbl_pk <- renderTable({
    d <- sim()
    f <- function(x, nm, u) {
      if (max(x) < 1e-6) return(NULL)
      data.frame(Item = nm, Cmax = sprintf("%.2f", max(x)),
                 Units = u, stringsAsFactors = FALSE)
    }
    rows <- rbind(f(d$CIBU, "Ibuprofen (total)", "mg/L"),
                  f(d$UIBU, "Ibuprofen (free)", "uM"),
                  f(d$CIND, "Indomethacin (total)", "mg/L"),
                  f(d$UIND, "Indomethacin (free)", "uM"),
                  f(d$CAPAP, "Acetaminophen (total)", "mg/L"),
                  f(d$UAPAP, "Acetaminophen (free)", "uM"))
    if (is.null(rows)) rows <- data.frame(Item = "no treatment", Cmax = "-", Units = "-")
    rbind(rows, data.frame(Item = "Maximum ductal COX inhibition",
                           Cmax = sprintf("%.1f", 100 * max(d$ICOXD)),
                           Units = "%"))
  }, digits = 3)

  ## -------------------------------------------------------------- tab 3
  output$plot_duct <- renderPlot({
    d <- sim(); r <- sim_ref()
    df <- rbind(data.frame(DAY = d$DAY, v = d$DDUCT, arm = "Treated"),
                data.frame(DAY = r$DAY, v = r$DDUCT, arm = "Untreated"))
    ggplot(df, aes(DAY, v, colour = arm)) +
      geom_hline(yintercept = 0.30, linetype = 2, colour = "#c62828") +
      geom_line(linewidth = 1) +
      scale_colour_manual(values = c("Treated" = "#00796b", "Untreated" = "#90a4ae")) +
      labs(x = "Postnatal day", y = "Ductal lumen (mm)", colour = NULL,
           title = "Ductal lumen — treatment vs natural course") + THEME
  })

  output$plot_wall <- renderPlot({
    d <- sim()
    df <- d %>% select(DAY, TONE, REMOD, WALLO2) %>%
      mutate(WALLO2 = WALLO2 / 60) %>% pivot_longer(-DAY) %>%
      mutate(name = recode(name, TONE = "Ductal tone (0-1)",
                           REMOD = "Neointimal cushion (0-1)",
                           WALLO2 = "Ductal wall PO2 / 60 mmHg"))
    ggplot(df, aes(DAY, value, colour = name)) +
      geom_line(linewidth = 0.95) +
      coord_cartesian(xlim = c(0, min(40, input$horizon))) +
      labs(x = "Postnatal day", y = NULL, colour = NULL,
           title = "Constriction · wall oxygen · permanent remodelling") + THEME
  })

  output$plot_drive <- renderPlot({
    d <- sim()
    df <- d %>% select(DAY, RELAXP, RELAXN, O2GAIN, NETD, TONETGT) %>%
      pivot_longer(-DAY) %>%
      mutate(name = recode(name, RELAXP = "PGE2/EP4 relaxation",
                           RELAXN = "NO relaxation", O2GAIN = "O2 constriction gain",
                           NETD = "Net constriction drive",
                           TONETGT = "Target tone"))
    ggplot(df, aes(DAY, value, colour = name)) +
      geom_line(linewidth = 0.85) +
      coord_cartesian(xlim = c(0, min(40, input$horizon))) +
      labs(x = "Postnatal day", y = NULL, colour = NULL,
           title = "Tone is a sigmoid of the net drive, not a product of gains") + THEME
  })

  ## -------------------------------------------------------------- tab 4
  output$plot_shunt <- renderPlot({
    d <- sim(); r <- sim_ref()
    df <- rbind(data.frame(DAY = d$DAY, v = d$QSH, arm = "Treated"),
                data.frame(DAY = r$DAY, v = r$QSH, arm = "Untreated"))
    ggplot(df, aes(DAY, v, colour = arm)) +
      geom_hline(yintercept = 100, linetype = 2, colour = "#c62828") +
      geom_line(linewidth = 1) +
      scale_colour_manual(values = c("Treated" = "#c62828", "Untreated" = "#90a4ae")) +
      coord_cartesian(xlim = c(0, min(40, input$horizon))) +
      annotate("text", x = 30, y = 112, size = 3.2, colour = "#c62828",
               label = "significance threshold") +
      labs(x = "Postnatal day", y = "Left-to-right shunt (mL/min/kg)", colour = NULL,
           title = "The shunt — revealed as PVR falls") + THEME
  })

  output$plot_press <- renderPlot({
    d <- sim()
    df <- d %>% select(DAY, PAOM, PPA, PDIA, QPQS) %>%
      mutate(QPQS = QPQS * 10) %>% pivot_longer(-DAY) %>%
      mutate(name = recode(name, PAOM = "Mean aortic pressure (mmHg)",
                           PPA = "Mean pulmonary artery pressure (mmHg)",
                           PDIA = "Diastolic pressure (mmHg)",
                           QPQS = "Qp:Qs x 10"))
    ggplot(df, aes(DAY, value, colour = name)) + geom_line(linewidth = 0.85) +
      coord_cartesian(xlim = c(0, min(40, input$horizon))) +
      labs(x = "Postnatal day", y = NULL, colour = NULL, title = "Pressures and flow ratio") +
      THEME
  })

  output$plot_organ <- renderPlot({
    d <- sim()
    df <- d %>% select(DAY, QCERREL, QMESREL, QRENREL) %>%
      pivot_longer(-DAY) %>%
      mutate(name = recode(name, QCERREL = "Brain", QMESREL = "Intestine",
                           QRENREL = "Kidney"))
    ggplot(df, aes(DAY, 100 * value, colour = name)) +
      geom_line(linewidth = 0.9) +
      coord_cartesian(xlim = c(0, min(20, input$horizon))) +
      labs(x = "Postnatal day", y = "Blood flow vs baseline (%)", colour = NULL,
           title = "Diastolic runoff + drug vasoconstriction") + THEME
  })

  output$plot_burden <- renderPlot({
    d <- sim(); r <- sim_ref()
    df <- rbind(data.frame(DAY = d$DAY, v = d$PDABUR, arm = "Treated"),
                data.frame(DAY = r$DAY, v = r$PDABUR, arm = "Untreated"))
    ggplot(df, aes(DAY, v, colour = arm)) + geom_line(linewidth = 1) +
      scale_colour_manual(values = c("Treated" = "#37474f", "Untreated" = "#90a4ae")) +
      labs(x = "Postnatal day", y = "PDA burden (cumulative)", colour = NULL,
           title = "What outcome sees is burden, not closure") + THEME
  })

  ## -------------------------------------------------------------- tab 5
  output$plot_outcome <- renderPlot({
    d <- sim(); r <- sim_ref()
    mk <- function(x, arm) x %>% select(DAY, PBPD, PNEC, PIVH, PDTH, PCOMP) %>%
      pivot_longer(-DAY) %>% mutate(arm = arm)
    df <- rbind(mk(d, "Treated"), mk(r, "Untreated")) %>%
      mutate(name = recode(name, PBPD = "Moderate-severe BPD", PNEC = "NEC >= stage 2",
                           PIVH = "Severe IVH", PDTH = "Death",
                           PCOMP = "Death or BPD"))
    ggplot(df, aes(DAY, 100 * value, colour = name, linetype = arm)) +
      geom_line(linewidth = 0.85) +
      labs(x = "Postnatal day", y = "Cumulative probability (%)", colour = NULL, linetype = NULL,
           title = "Outcome risk") + THEME
  })

  output$tbl_outcome <- renderTable({
    d <- sim(); r <- sim_ref()
    g <- function(x, nm) x[[nm]][nrow(x)] * 100
    data.frame(
      Endpoint = c("Death or moderate-severe BPD", "Moderate-severe BPD", "Death",
                     "NEC >= stage 2", "Severe IVH", "Spontaneous intestinal perforation", "PDA burden"),
      Treated = c(g(d, "PCOMP"), g(d, "PBPD"), g(d, "PDTH"), g(d, "PNEC"),
               g(d, "PIVH"), g(d, "PSIP"), d$PDABUR[nrow(d)]),
      Untreated = c(g(r, "PCOMP"), g(r, "PBPD"), g(r, "PDTH"), g(r, "PNEC"),
                 g(r, "PIVH"), g(r, "PSIP"), r$PDABUR[nrow(r)]),
      stringsAsFactors = FALSE)
  }, digits = 1)

  ## -------------------------------------------------------------- tab 6
  scen <- eventReactive(input$run_scen, {
    withProgress(message = "Running 16 scenarios...", value = 0.2,
                 run_scenarios(days = 90))
  })
  if (have_DT) {
    output$tbl_scen <- DT::renderDataTable({
      DT::datatable(scen(), options = list(pageLength = 16, dom = "t",
                                          scrollX = TRUE), rownames = FALSE)
    })
  } else {
    output$tbl_scen_plain <- renderTable({ scen() }, digits = 2)
  }

  output$plot_scen <- renderPlot({
    s <- scen()
    s$scenario <- factor(s$id, levels = rev(s$id))
    ggplot(s, aes(P_comp, scenario)) +
      geom_col(aes(fill = burden), width = 0.72) +
      scale_fill_gradient(low = "#b2dfdb", high = "#c62828") +
      labs(x = "Death or moderate-severe BPD (%)", y = NULL, fill = "PDA burden",
           title = "The low-burden scenarios have the better outcomes — not closure itself") +
      THEME
  })

  ## -------------------------------------------------------------- tab 7
  output$plot_renal <- renderPlot({
    d <- sim(); r <- sim_ref()
    mk <- function(x, arm) data.frame(DAY = x$DAY, SCr = x$SCR, UO = x$UO,
                                      GFR = x$GFR, arm = arm)
    df <- rbind(mk(d, "Treated"), mk(r, "Untreated")) %>% pivot_longer(-c(DAY, arm))
    ggplot(df, aes(DAY, value, colour = name, linetype = arm)) +
      geom_line(linewidth = 0.9) +
      coord_cartesian(xlim = c(0, min(20, input$horizon))) +
      labs(x = "Postnatal day", y = NULL, colour = NULL, linetype = NULL,
           title = "Kidney: SCr (mg/dL) · urine output (mL/kg/h) · GFR") + THEME
  })

  output$plot_safety <- renderPlot({
    d <- sim()
    df <- d %>% select(DAY, ALT, TXA2, BTIME, BFREE) %>%
      mutate(ALT = ALT / 20) %>% pivot_longer(-DAY) %>%
      mutate(name = recode(name, ALT = "ALT / 20 U/L",
                           TXA2 = "Platelet TXA2 (fraction)",
                           BTIME = "Bleeding time index",
                           BFREE = "Free bilirubin index"))
    ggplot(df, aes(DAY, value, colour = name)) + geom_line(linewidth = 0.9) +
      coord_cartesian(xlim = c(0, min(20, input$horizon))) +
      labs(x = "Postnatal day", y = NULL, colour = NULL, title = "Safety indices") +
      THEME
  })

  output$tbl_safety <- renderTable({
    d <- sim(); r <- sim_ref()
    e <- d[d$DAY <= 12, ]; e0 <- r[r$DAY <= 12, ]
    data.frame(
      Metric = c("Peak SCr (mg/dL)", "Lowest urine output (mL/kg/h)", "Lowest GFR (mL/min/kg)",
               "Lowest cerebral flow (%)", "Lowest intestinal flow (%)", "Peak ALT (U/L)",
               "Peak bleeding time index", "Peak free bilirubin index"),
      Treated = c(max(e$SCR), min(e$UO), min(e$GFR), 100 * min(e$QCERREL),
               100 * min(e$QMESREL), max(e$ALT), max(e$BTIME), max(e$BFREE)),
      Untreated = c(max(e0$SCR), min(e0$UO), min(e0$GFR), 100 * min(e0$QCERREL),
                 100 * min(e0$QMESREL), max(e0$ALT), max(e0$BTIME),
                 max(e0$BFREE)), stringsAsFactors = FALSE)
  }, digits = 2)

  ## -------------------------------------------------------------- tab 8
  spont <- eventReactive(input$run_spont, {
    withProgress(message = "Gestational-age sweep...", value = 0.2,
                 spontaneous_closure_table())
  })
  output$tbl_spont <- renderTable({ spont() }, digits = 2)
  output$plot_spont <- renderPlot({
    s <- spont()
    tgt <- data.frame(GA = c(24, 25, 26, 27, 28, 30, 32, 38),
                      obs = c(85, 71, 48, 32, 21, 10, 5, 1.5))
    ggplot(s, aes(GA, closure_day)) +
      geom_line(linewidth = 1.1, colour = "#3949ab") +
      geom_point(size = 2, colour = "#3949ab") +
      geom_point(data = tgt, aes(GA, obs), colour = "#c62828",
                 shape = 4, size = 3, stroke = 1.2) +
      scale_y_log10() +
      labs(x = "Gestational age (wk)", y = "Spontaneous closure day (log)",
           title = "Model (line) vs observed target (×)") + THEME
  })

  ## -------------------------------------------------------------- tab 9
  output$tbl_trials <- renderTable({
    data.frame(
      Trial = c("Semberova 2017", "Baby-OSCAR 2024", "BeNeDuctus 2023",
               "TIPP 2001", "PDA-TOLERATE 2019", "Cochrane / NMA"),
      Observation = c("<27 wk 94-98% spontaneous closure, <26 wk median 71 days",
               "Early ibuprofen vs expectant: death/BPD 69.4% vs 63.5%",
               "Expectant management non-inferior: composite 46.5% vs 63.5%",
               "Prophylactic indomethacin: severe IVH 9% vs 13%, no change in 18-month outcome",
               "No benefit from early treatment at 6-14 days",
               "Closure rates similar for all three drugs; ibuprofen has less NEC and AKI"),
      Role = c("Fitted (involution parameters)", "Fitted (3 risk coefficients)",
               "Predicted", "Predicted", "Predicted", "Fitted (Ki/IC50) + predicted (toxicity separation)"),
      stringsAsFactors = FALSE)
  })
}

shinyApp(ui, server)
