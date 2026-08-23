## =============================================================================
##  Familial Chylomicronemia Syndrome (FCS) — QSP Shiny dashboard
##  fcs_shiny_app.R
## -----------------------------------------------------------------------------
##  Run with:
##      setwd("familial-chylomicronemia-syndrome")
##      shiny::runApp("fcs_shiny_app.R")
##
##  Requires: shiny, mrgsolve, ggplot2, DT  (dplyr/tidyr optional)
##  The app sources fcs_mrgsolve_model.R, so keep both files together.
##
##  The dashboard is organised around the model's central claim: TRL clearance
##  is a SUM of two limbs, the genotype multiplies the first by zero, and every
##  therapeutic decision in FCS is therefore a decision about the second limb
##  or about the input.  Each tab makes one part of that argument visible.
## =============================================================================

library(shiny)
library(mrgsolve)
library(ggplot2)

if(!exists("mod")) source("fcs_mrgsolve_model.R")

THEME <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 13),
        legend.position = "bottom")

PAL <- c(placebo = "#90a4ae", diet = "#7cb342", fibrate = "#546e7a",
         evinacumab = "#6d4c41", volanesorsen = "#00897b",
         olezarsen = "#43a047", plozasiran = "#0288d1")

## =============================================================================
##  UI
## =============================================================================
ui <- fluidPage(
  titlePanel(
    div(style = "line-height:1.25",
        strong("Familial Chylomicronemia Syndrome (FCS) — QSP dashboard"),
        div(style = "font-size:13px;color:#555;font-weight:normal",
            "FCS as a saturated clearance capacity · ",
            "CL_total = CL_LPL(genotype ≈ 0) + CL_independent(apoC-III, saturable)"))),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      h4("① Patient"),
      selectInput("geno", "Genotype",
                  choices = c("LPL-null FCS"            = "fcs_null",
                              "APOC2-deficient FCS" = "fcs_apoc2",
                              "GPIHBP1-deficient FCS" = "fcs_gpihbp1",
                              "Residual LPL 5% (partial deficiency)" = "partial_05",
                              "Residual LPL 20% (mild)" = "partial_20",
                              "Multifactorial chylomicronemia syndrome (MCS)" = "mcs",
                              "Normal control (healthy)" = "healthy"),
                  selected = "fcs_null"),
      sliderInput("fat", "Dietary long-chain fat, LCT (g/day)", 5, 100, 20, step = 5),
      sliderInput("fmct", "MCT substitution fraction", 0, 0.8, 0, step = 0.1),
      checkboxInput("alc",  "Alcohol", FALSE),
      checkboxInput("est",  "Oral oestrogen", FALSE),
      checkboxInput("preg", "Pregnancy", FALSE),

      hr(),
      h4("② Therapy"),
      selectInput("drug", "Drug",
                  choices = c("None (diet only)" = "none",
                              "Fenofibrate + ω-3" = "fibrate",
                              "Evinacumab 15 mg/kg q4w" = "evinacumab",
                              "Volanesorsen 300 mg qw" = "volanesorsen",
                              "Olezarsen 80 mg q4w" = "olezarsen",
                              "Plozasiran 25 mg q12w" = "plozasiran"),
                  selected = "plozasiran"),
      conditionalPanel("input.drug == 'olezarsen'",
                       sliderInput("ole_mg", "Olezarsen dose (mg)", 50, 80, 80, 10)),
      conditionalPanel("input.drug == 'plozasiran'",
                       sliderInput("plo_mg", "Plozasiran dose (mg)", 10, 50, 25, 5)),
      sliderInput("days", "Observation period (days)", 90, 730, 365, step = 30),
      checkboxInput("binge", "60 g fat binge every 30 days (non-adherence)", FALSE),

      hr(),
      h4("③ Model parameters"),
      sliderInput("vmax", "Vmax(limb 2, LPL-independent) mg/h", 600, 2600, 1300, 100),
      sliderInput("imax", "apoC-III maximum amplification, IMAX_C3", 0.5, 4, 2.2, 0.1),
      sliderInput("hill", "AP risk convexity, HILL_AP", 1.0, 3.0, 1.7, 0.1),
      actionButton("go", "Run simulation", class = "btn-primary btn-block"),

      hr(),
      div(style = "font-size:11px;color:#777",
          "This model is for education and research purposes only. It must not be used for clinical decision-making.",
          br(), "Research and education only — not for clinical use.")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs", type = "tabs",

        ## ---------------------------------------------------------------- 1
        tabPanel("① Patient profile",
          br(),
          fluidRow(
            column(4, wellPanel(h4(textOutput("kpi_tg")),  p("Mean TG (mg/dL)"))),
            column(4, wellPanel(h4(textOutput("kpi_c3")),  p("Plasma apoC-III (mg/dL)"))),
            column(4, wellPanel(h4(textOutput("kpi_ap")),  p("Pancreatitis risk (events/py)")))),
          fluidRow(
            column(4, wellPanel(h4(textOutput("kpi_tat")), p("Days TG > 880 / year"))),
            column(4, wellPanel(h4(textOutput("kpi_plt")), p("Platelet nadir (10⁹/L)"))),
            column(4, wellPanel(h4(textOutput("kpi_lpl")), p("LPL clearance rate (dL/h)")))),
          hr(),
          h4("Patient summary"),
          verbatimTextOutput("profile_txt")),

        ## ---------------------------------------------------------------- 2
        tabPanel("② Lipid time course (PD)",
          br(), plotOutput("p_tg", height = "330px"),
          plotOutput("p_species", height = "300px"),
          helpText("Dashed line = 880 mg/dL (10 mmol/L) acute-pancreatitis threshold. ",
                   "In FCS a truly fasting state never exists, and postprandial fluctuation ",
                   "accounts for much of the risk.")),

        ## ---------------------------------------------------------------- 3
        tabPanel("③ Saturation cliff",
          br(), plotOutput("p_sat", height = "420px"),
          DT::dataTableOutput("t_sat"),
          helpText("In normal physiology, fat intake and TG are nearly linear, but ",
                   "in LPL-null FCS it becomes hyperbolic. The reason '<20 g/day' is ",
                   "a cliff edge rather than a slope.")),

        ## ---------------------------------------------------------------- 4
        tabPanel("④ Limb decomposition",
          br(), plotOutput("p_limb", height = "330px"),
          DT::dataTableOutput("t_limb"),
          helpText("limb 1 (LPL) is exactly 0 mg/h in every LPL-null scenario, and ",
                   "no drug can move it. 0 × anything = 0.")),

        ## ---------------------------------------------------------------- 5
        tabPanel("⑤ Drug PK",
          br(), plotOutput("p_pk", height = "330px"),
          plotOutput("p_target", height = "300px"),
          helpText("GalNAc conjugation did not change the target, but reduced systemic ",
                   "exposure 20-30-fold. Volanesorsen's platelet signal was not ",
                   "pharmacology but a PK problem.")),

        ## ---------------------------------------------------------------- 6
        tabPanel("⑥ Clinical endpoints",
          br(),
          fluidRow(column(6, plotOutput("p_haz",  height = "300px")),
                   column(6, plotOutput("p_prob", height = "300px"))),
          plotOutput("p_tat", height = "260px"),
          helpText("Because the hazard function is convex in TG (Hill 1.7), a ",
                   "reduction in 'time above threshold' explains the actual event ",
                   "reduction better than a % change in mean TG.")),

        ## ---------------------------------------------------------------- 7
        tabPanel("⑦ Jensen gap (postprandial fluctuation)",
          br(), plotOutput("p_jensen", height = "330px"),
          DT::dataTableOutput("t_jensen"),
          helpText("E[λ(TG)] > λ(E[TG]). Assessing risk from a single fasting-TG ",
                   "value systematically underestimates it.")),

        ## ---------------------------------------------------------------- 8
        tabPanel("⑧ Safety",
          br(),
          fluidRow(column(6, plotOutput("p_plt", height = "300px")),
                   column(6, plotOutput("p_alt", height = "300px"))),
          DT::dataTableOutput("t_safety"),
          helpText("Platelet reduction is proportional to systemic PS-ASO tissue ",
                   "exposure; it essentially disappears with GalNAc conjugates.")),

        ## ---------------------------------------------------------------- 9
        tabPanel("⑨ Scenario comparison",
          br(), plotOutput("p_cmp", height = "360px"),
          DT::dataTableOutput("t_ledger"),
          helpText("Compare dTG (%) and AP event reduction (%) side by side. ",
                   "The latter is always larger — the dividend of convexity.")),

        ## ---------------------------------------------------------------- 10
        tabPanel("⑩ Biomarkers · calibration",
          br(), plotOutput("p_bio", height = "330px"),
          h4("Literature calibration table"),
          DT::dataTableOutput("t_calib"),
          helpText("See fcs_references.md for the source of each target."))
      )
    )
  )
)

## =============================================================================
##  SERVER
## =============================================================================
server <- function(input, output, session) {

  pars <- reactive({
    p <- geno[[input$geno]]
    p$FMCT     <- input$fmct
    p$ALC      <- as.numeric(input$alc || isTRUE(p$ALC == 1))
    p$EST      <- as.numeric(input$est)
    p$PREG     <- as.numeric(input$preg)
    p$VMAX_IND <- input$vmax
    p$IMAX_C3  <- input$imax
    p$HILL_AP  <- input$hill
    if(input$drug == "fibrate") p$OM3_ON <- 1
    p
  })

  drug_events <- reactive({
    d <- input$days
    switch(input$drug,
      none         = NULL,
      fibrate      = fenofibrate(0, d),
      evinacumab   = evinacumab(0, ceiling(d/28)),
      volanesorsen = volanesorsen(0, ceiling(d/7)),
      olezarsen    = olezarsen(0, ceiling(d/28), mg = input$ole_mg),
      plozasiran   = plozasiran(0, ceiling(d/90), mg = input$plo_mg))
  })

  sim <- eventReactive(input$go, {
    e <- drug_events()
    if(input$binge) {
      b <- do.call(c, lapply(seq(30, input$days, by = 30),
                             function(k) fat_binge(60, k)))
      e <- if(is.null(e)) b else e + b
    }
    withProgress(message = "Running simulation...", value = 0.5,
      run_arm(pars(), events = e, days = input$days, fat = input$fat, delta = 2))
  }, ignoreNULL = FALSE)

  on_treat <- reactive({ d <- sim(); subset(d, day >= 0) })

  ## ------------------------------------------------------------------ KPIs
  output$kpi_tg  <- renderText(sprintf("%.0f", mean(on_treat()$TG)))
  output$kpi_c3  <- renderText(sprintf("%.1f", mean(on_treat()$APOC3)))
  output$kpi_ap  <- renderText({
    d <- on_treat(); sprintf("%.3f", (max(d$CUMHAZ)-min(d$CUMHAZ))/(max(d$day)/365))})
  output$kpi_tat <- renderText({
    d <- on_treat()
    sprintf("%.0f", (max(d$TAT880)-min(d$TAT880))/24/(max(d$day)/365))})
  output$kpi_plt <- renderText(sprintf("%.0f", min(on_treat()$PLT)))
  output$kpi_lpl <- renderText(sprintf("%.2f", mean(on_treat()$CLLPL)))

  output$profile_txt <- renderPrint({
    d <- on_treat(); b <- sim()
    base <- mean(b$TG[b$day > -30 & b$day <= 0])
    cat("Genotype:", input$geno, "\n")
    cat("Dietary LCT:", input$fat, "g/day  (MCT substitution", input$fmct, ")\n")
    cat("Therapy:", input$drug, "\n")
    cat("Pre-treatment mean TG:", sprintf("%.0f mg/dL (%.1f mmol/L)",
                                          base, base/88.57), "\n")
    cat("Post-treatment mean TG:", sprintf("%.0f mg/dL  (%+.1f %%)",
                                          mean(d$TG), 100*(mean(d$TG)/base-1)), "\n")
    cat("Chylomicron fraction:", sprintf("%.0f %%", 100*mean(d$CM_C)/mean(d$TG)), "\n")
    cat("Plasma viscosity index:", sprintf("%.2f", mean(d$VISCI)), "\n")
    cat("Eruptive xanthoma score:", sprintf("%.2f / 10", tail(d$XANTH, 1)), "\n")
    cat("Hepatosplenomegaly index:", sprintf("%.2f / 10", tail(d$HSM, 1)), "\n")
    cat("Cognitive fog score:", sprintf("%.2f / 10", tail(d$FOG, 1)), "\n")
    cat("1-year pancreatitis probability:", sprintf("%.1f %%", 100*max(d$PROB_AP)), "\n")
  })

  ## ------------------------------------------------------------------ tab 2
  output$p_tg <- renderPlot({
    d <- sim()
    ggplot(d, aes(day, TG)) +
      geom_hline(yintercept = 880, linetype = 2, colour = "grey40") +
      geom_vline(xintercept = 0, linetype = 3, colour = "grey60") +
      geom_line(colour = "#d84315", linewidth = 0.6) +
      labs(title = "Total plasma triglyceride",
           x = "Days since treatment start", y = "TG (mg/dL)") + THEME
  })

  output$p_species <- renderPlot({
    d <- sim()
    long <- rbind(
      data.frame(day = d$day, conc = d$CM_C,   species = "chylomicron-TG"),
      data.frame(day = d$day, conc = d$VLDL_C, species = "VLDL-TG"),
      data.frame(day = d$day, conc = d$REM_C,  species = "remnant-TG"))
    ggplot(long, aes(day, conc, colour = species)) +
      geom_line(linewidth = 0.6) +
      scale_colour_manual(values = c("chylomicron-TG" = "#0277bd",
                                     "VLDL-TG" = "#ef6c00",
                                     "remnant-TG" = "#7e57c2")) +
      labs(title = "Breakdown by lipoprotein species (chylomicrons dominate in FCS)",
           x = "day", y = "mg/dL", colour = NULL) + THEME
  })

  ## ------------------------------------------------------------------ tab 3
  sat_tab <- reactive({
    fats <- c(5,10,15,20,25,30,40,50,60,80,100)
    do.call(rbind, lapply(c("fcs_null","partial_05","partial_20","healthy"),
      function(g) do.call(rbind, lapply(fats, function(f) {
        d  <- run_arm(geno[[g]], days = 30, burn = 120, fat = f, delta = 4)
        ss <- subset(d, day > 18)
        data.frame(genotype = g, fat = f, TG = mean(ss$TG),
                   peak = max(ss$TG),
                   pct_above = 100*mean(ss$TG > 880))
      }))))
  })

  output$p_sat <- renderPlot({
    ggplot(sat_tab(), aes(fat, TG, colour = genotype)) +
      geom_hline(yintercept = 880, linetype = 2, colour = "grey40") +
      geom_line(linewidth = 0.9) + geom_point(size = 2) +
      scale_colour_manual(values = c(fcs_null = "#d84315", partial_05 = "#f9a825",
                                     partial_20 = "#43a047", healthy = "#1565c0")) +
      labs(title = "Saturation cliff: steady-state TG vs dietary fat",
           x = "Dietary long-chain fat (g/day)", y = "Steady-state TG (mg/dL)",
           colour = NULL) + THEME
  })
  output$t_sat <- DT::renderDataTable(
    DT::datatable(sat_tab(), options = list(pageLength = 8), rownames = FALSE) |>
      DT::formatRound(c("TG","peak","pct_above"), 1))

  ## ------------------------------------------------------------------ tab 4
  output$p_limb <- renderPlot({
    d <- sim()
    long <- rbind(
      data.frame(day = d$day, flux = d$FLUX1, limb = "limb 1 · LPL"),
      data.frame(day = d$day, flux = d$FLUX2, limb = "limb 2 · LPL-independent"),
      data.frame(day = d$day, flux = d$FLUX3, limb = "residual scavenging"))
    ggplot(long, aes(day, flux, colour = limb)) +
      geom_line(linewidth = 0.7) +
      scale_colour_manual(values = c("limb 1 · LPL" = "#d81b60",
                                     "limb 2 · LPL-independent" = "#7e57c2",
                                     "residual scavenging" = "#90a4ae")) +
      labs(title = "Flux by clearance pathway (mg/h)", x = "day", y = "flux (mg/h)",
           colour = NULL) + THEME
  })
  output$t_limb <- DT::renderDataTable({
    d <- on_treat()
    tb <- data.frame(
      limb = c("limb 1 (LPL)", "limb 2 (LPL-independent)", "residual"),
      mean_flux_mg_h = c(mean(d$FLUX1), mean(d$FLUX2), mean(d$FLUX3)))
    tb$percent <- 100*tb$mean_flux_mg_h/sum(tb$mean_flux_mg_h)
    DT::datatable(tb, rownames = FALSE, options = list(dom = "t")) |>
      DT::formatRound(c("mean_flux_mg_h","percent"), 2)
  })

  ## ------------------------------------------------------------------ tab 5
  output$p_pk <- renderPlot({
    d <- sim()
    long <- rbind(
      data.frame(day = d$day, amt = d$VOL_LIV,  cmt = "volanesorsen · liver"),
      data.frame(day = d$day, amt = d$VOL_SYS,  cmt = "volanesorsen · systemic"),
      data.frame(day = d$day, amt = d$OLE_LIV,  cmt = "olezarsen · liver"),
      data.frame(day = d$day, amt = d$OLE_SYS,  cmt = "olezarsen · systemic"),
      data.frame(day = d$day, amt = d$PLO_RISC, cmt = "plozasiran · RISC"),
      data.frame(day = d$day, amt = d$CEVI2*10, cmt = "evinacumab · plasma x10"))
    long <- long[long$amt > 1e-6, ]
    if(!nrow(long)) return(ggplot() + labs(title = "No drug") + THEME)
    ggplot(long, aes(day, amt, colour = cmt)) + geom_line(linewidth = 0.7) +
      labs(title = "Drug tissue/effect compartments (mg or AU)", x = "day", y = "amount",
           colour = NULL) + THEME
  })
  output$p_target <- renderPlot({
    d <- sim()
    ggplot(d, aes(day, APOC3)) + geom_line(colour = "#00897b", linewidth = 0.8) +
      labs(title = "Target engagement: plasma apoC-III", x = "day",
           y = "apoC-III (mg/dL)") + THEME
  })

  ## ------------------------------------------------------------------ tab 6
  output$p_haz <- renderPlot({
    d <- sim()
    ggplot(d, aes(day, HAZ_YR)) + geom_line(colour = "#bf360c", linewidth = 0.7) +
      labs(title = "Acute pancreatitis instantaneous hazard", x = "day",
           y = "events / patient-year") + THEME
  })
  output$p_prob <- renderPlot({
    d <- on_treat()
    ggplot(d, aes(day, 100*PROB_AP)) +
      geom_line(colour = "#5e35b1", linewidth = 0.8) +
      labs(title = "Cumulative pancreatitis probability", x = "day", y = "%") + THEME
  })
  output$p_tat <- renderPlot({
    d <- on_treat()
    ggplot(d, aes(day, TAT880/24)) +
      geom_line(colour = "#4527a0", linewidth = 0.8) +
      labs(title = "TG > 880 mg/dL cumulative time above threshold",
           x = "day", y = "cumulative days") + THEME
  })

  ## ------------------------------------------------------------------ tab 7
  jen <- eventReactive(input$go, { FCS_jensen_gap() })
  output$t_jensen <- DT::renderDataTable(
    DT::datatable(jen(), rownames = FALSE, options = list(dom = "t")) |>
      DT::formatRound(2:9, 3))
  output$p_jensen <- renderPlot({
    j <- jen()
    long <- rbind(
      data.frame(fat = j$fat_g_day, h = j$haz_of_fasting, k = "λ(fasting TG)"),
      data.frame(fat = j$fat_g_day, h = j$haz_of_mean,    k = "λ(mean TG)"),
      data.frame(fat = j$fat_g_day, h = j$mean_of_haz,    k = "mean λ(TG) — true risk"))
    ggplot(long, aes(factor(fat), h, fill = k)) +
      geom_col(position = "dodge") +
      scale_fill_manual(values = c("#90a4ae", "#f9a825", "#bf360c")) +
      labs(title = "Jensen gap: convex hazard function + postprandial fluctuation",
           x = "Dietary fat (g/day)", y = "events / patient-year", fill = NULL) +
      THEME
  })

  ## ------------------------------------------------------------------ tab 8
  output$p_plt <- renderPlot({
    d <- sim()
    ggplot(d, aes(day, PLT)) +
      geom_hline(yintercept = 100, linetype = 2, colour = "#d84315") +
      geom_hline(yintercept = 50,  linetype = 3, colour = "#b71c1c") +
      geom_line(colour = "#5e35b1", linewidth = 0.8) +
      labs(title = "Platelet count", x = "day", y = "10⁹/L") + THEME
  })
  output$p_alt <- renderPlot({
    d <- sim()
    ggplot(d, aes(day, ALT)) + geom_line(colour = "#6d4c41", linewidth = 0.8) +
      labs(title = "Serum ALT", x = "day", y = "U/L") + THEME
  })
  output$t_safety <- DT::renderDataTable({
    d <- on_treat()
    tb <- data.frame(
      metric = c("Platelet nadir (10⁹/L)", "Time with platelets <100 (%)",
                 "ALT peak (U/L)", "Systemic ASO tissue amount (mg)",
                 "Hepatic ASO tissue amount (mg)"),
      value = c(min(d$PLT), 100*mean(d$PLT < 100), max(d$ALT),
                max(d$VOL_SYS + d$OLE_SYS), max(d$VOL_LIV + d$OLE_LIV)))
    DT::datatable(tb, rownames = FALSE, options = list(dom = "t")) |>
      DT::formatRound("value", 2)
  })

  ## ------------------------------------------------------------------ tab 9
  ledger <- eventReactive(input$go, { FCS_trial_ledger(days = input$days) })
  output$t_ledger <- DT::renderDataTable(
    DT::datatable(ledger(), rownames = FALSE,
                  options = list(pageLength = 8, scrollX = TRUE)) |>
      DT::formatRound(2:12, 2))
  output$p_cmp <- renderPlot({
    L <- ledger()
    long <- rbind(
      data.frame(arm = L$arm, v = -L$dTG_m12_pct,     k = "TG reduction (%)"),
      data.frame(arm = L$arm, v = L$AP_reduction_pct, k = "Pancreatitis event reduction (%)"))
    ggplot(long, aes(reorder(arm, v), v, fill = k)) +
      geom_col(position = "dodge") + coord_flip() +
      scale_fill_manual(values = c("#0288d1", "#bf360c")) +
      labs(title = "% TG change vs actual event reduction — the dividend of convexity",
           x = NULL, y = "%", fill = NULL) + THEME
  })

  ## ------------------------------------------------------------------ tab 10
  output$p_bio <- renderPlot({
    d <- sim()
    long <- rbind(
      data.frame(day = d$day, v = d$XANTH, k = "Eruptive xanthoma"),
      data.frame(day = d$day, v = d$HSM,   k = "Hepatosplenomegaly index"),
      data.frame(day = d$day, v = d$FOG,   k = "Cognitive fog"),
      data.frame(day = d$day, v = d$VISCI, k = "Plasma viscosity index"))
    ggplot(long, aes(day, v, colour = k)) + geom_line(linewidth = 0.7) +
      facet_wrap(~k, scales = "free_y") +
      labs(title = "Non-pancreatic disease-burden biomarkers", x = "day", y = NULL) +
      THEME + theme(legend.position = "none")
  })

  output$t_calib <- DT::renderDataTable({
    tb <- data.frame(
      target = c("Normal fasting TG",
                 "FCS · 10 g fat/day",
                 "FCS · 20 g fat/day",
                 "FCS · 60 g fat/day",
                 "Residual LPL 5%",
                 "Volanesorsen apoC-III",
                 "Volanesorsen TG (3 months)",
                 "Volanesorsen platelets <100",
                 "Olezarsen TG (12 months, open-label extension)",
                 "Olezarsen pancreatitis events",
                 "Plozasiran TG (10 months)",
                 "Plozasiran pancreatitis events",
                 "Fibrate+ω-3 (LPL-null)",
                 "Placebo-arm pancreatitis incidence"),
      published = c("< 150 mg/dL", "~600-800 mg/dL", "1000-2000 mg/dL",
                    "> 5000 mg/dL", "300-500 mg/dL", "-76 ± 8 %", "-77 %",
                    "about 47 % of patients", "-73.7 %", "1 event vs 11 with placebo", "about -80 %",
                    "-83 %", "-10 ~ -20 %", "0.20-0.30 /patient-year"),
      source = c("General reference", "Williams 2018 / clinical experience", "Brahm 2015",
                 "Gaudet 2014", "Rahalkar 2009", "Witztum NEJM 2019 (APPROACH)",
                 "Witztum NEJM 2019", "Witztum NEJM 2019",
                 "Stroes NEJM 2024 (Balance)", "Stroes NEJM 2024",
                 "Watts NEJM 2024/25 (PALISADE)", "Watts NEJM 2024/25",
                 "Falko 2018 / Gaudet 2014", "Balance·PALISADE placebo arm"))
    DT::datatable(tb, rownames = FALSE, options = list(pageLength = 14, dom = "t"))
  })
}

shinyApp(ui, server)
