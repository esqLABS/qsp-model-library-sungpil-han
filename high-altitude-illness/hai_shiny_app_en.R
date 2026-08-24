## =============================================================================
##  HIGH-ALTITUDE ILLNESS QSP — Shiny dashboard
##  High-altitude illness (AMS · HACE · HAPE) interactive dashboard
## =============================================================================
##
##  Run:
##    source("hai_mrgsolve_model_en.R")   # model + scenario definitions
##    shiny::runApp("hai_shiny_app.R")
##
##  12 tabs:
##    ① Patient profile     ② Ascent profile    ③ Gas exchange   ④ Ventilation·acid-base
##    ⑤ Sleep·CO₂ reserve   ⑥ Brain (AMS/HACE)  ⑦ Lung (HAPE)    ⑧ Clinical endpoints
##    ⑨ Scenario comparison ⑩ Descent equivalent ⑪ Biomarkers    ⑫ Bifurcation·sensitivity
##
##  The editorial policy of this dashboard: every tab shows one thing at a time —
##  "what moves what".  Tab ⑩ (descent equivalent) in particular converts every
##  intervention into a single currency — "how many metres of descent is it worth" —
##  and lines them up. That dexamethasone comes out at 0 m is the key screen of this app.
## =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)

if (!exists("mod")) source("hai_mrgsolve_model_en.R")

THEME <- theme_bw(base_size = 12) +
  theme(strip.background = element_rect(fill = "#eceff1"),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")

PAL <- c("#1565c0", "#c62828", "#2e7d32", "#ef6c00", "#6a1b9a", "#00838f",
         "#5d4037", "#455a64")

## =============================================================================
##  UI
## =============================================================================
ui <- fluidPage(
  titlePanel("High-Altitude Illness QSP model — AMS · HACE · HAPE"),
  tags$p(style = "color:#555;margin-top:-8px",
         HTML("One inspired oxygen tension, three defence mechanisms, three time constants. ",
              "<b>P<sub>IO2</sub> = F<sub>IO2</sub> × (P<sub>B</sub>(h) − 47)</b>")),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      h4("① Individual phenotype"),
      selectInput("pheno", "Phenotype preset",
                  c("typical trekker" = "typ",
                    "HAPE-susceptible" = "hapes",
                    "Tight-fit (low intracranial compliance)" = "tight",
                    "elite climber (high HVR)" = "elite",
                    "custom" = "custom")),
      conditionalPanel(
        "input.pheno == 'custom'",
        sliderInput("hvr",  "HVR multiplier (hypoxic ventilatory response)", 0.4, 2.5, 1.0, 0.05),
        sliderInput("lam",  "λ_max — maximum HPV constriction factor", 1, 14, 5, 0.5),
        sliderInput("ahet", "a — fraction of the vascular bed responding to HPV (heterogeneity)", 0.05, 0.92, 0.50, 0.01),
        sliderInput("pvi",  "PVI — craniospinal pressure-volume index (mL)", 10, 40, 25, 1),
        sliderInput("hbsl", "Sea-level Hb (g/dL)", 11, 19, 15, 0.1)
      ),
      helpText(HTML("<small>The <b>ceiling</b> on capillary over-perfusion is 1/[1−a(1−κ)]. ",
                    "The strength of HPV sets only how quickly it is reached.</small>")),

      hr(),
      h4("② Ascent profile"),
      selectInput("profile", "Profile",
                  c("rapid ascent to 4559 m (Capanna Margherita)" = "rapid",
                    "graded ascent to 4559 m (~400 m/day)"  = "graded",
                    "climb high, sleep low"                 = "chsl",
                    "sleeping at 4000 m"                    = "sleep4000",
                    "3 weeks acclimatising at 3800 m"       = "accl3800",
                    "Everest summit day (8848 m)"           = "everest",
                    "custom (constant altitude)"            = "flat")),
      conditionalPanel("input.profile == 'flat'",
                       sliderInput("flatalt", "Altitude (m)", 0, 8848, 4500, 50)),
      sliderInput("tend", "Simulation horizon (h)", 24, 504, 120, 12),
      sliderInput("exer", "Daytime exercise intensity (× resting V̇O₂)", 1, 4, 1, 0.1),
      helpText(HTML("<small>Exercise <b>multiplies</b> the whole P<sub>cap</sub> term. ",
                    "The bifurcation parameter a clinician can actually control.</small>")),

      hr(),
      h4("③ Drugs / interventions"),
      checkboxInput("acz", "Acetazolamide", FALSE),
      conditionalPanel("input.acz",
                       selectInput("aczdose", "Dose", c("125 mg bid" = 125,
                                                        "250 mg bid" = 250))),
      checkboxInput("dex", "Dexamethasone 4 mg q12h", FALSE),
      checkboxInput("nif", "Nifedipine SR 30 mg bid", FALSE),
      checkboxInput("tad", "Tadalafil 10 mg bid", FALSE),
      checkboxInput("sal", "Salmeterol 125 µg bid", FALSE),
      checkboxInput("ibu", "Ibuprofen 600 mg tid", FALSE),
      sliderInput("fio2", "Inspired oxygen fraction F_IO₂", 0.2094, 0.60, 0.2094, 0.005),
      sliderInput("bag",  "Gamow bag pressurisation (mmHg)", 0, 220, 0, 5),

      hr(),
      actionButton("go", "Run simulation", class = "btn-primary btn-block")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs", type = "tabs",

        ## ------------------------------------------------------------------
        tabPanel("① Patient profile",
          h4("Phenotype summary"),
          tableOutput("phenotab"),
          h4("What this phenotype changes"),
          plotOutput("phenoplot", height = "380px"),
          helpText(HTML(
            "Left: where an individual's HVR puts the steady-state S<sub>aO2</sub> at each altitude. ",
            "Right: where the heterogeneity a puts the <b>ceiling</b> on capillary pressure — ",
            "as λ→∞ the amplification limit converges to 1/[1−a(1−κ)] (which is 1/(1−a) only when κ = 0)."))
        ),

        ## ------------------------------------------------------------------
        tabPanel("② Ascent profile",
          plotOutput("profplot", height = "260px"),
          h4("This is all that altitude does"),
          plotOutput("pio2plot", height = "300px"),
          tableOutput("statictab")
        ),

        ## ------------------------------------------------------------------
        tabPanel("③ Gas exchange",
          plotOutput("gasplot", height = "460px"),
          h4("Drug concentrations"),
          plotOutput("pkplot", height = "240px")
        ),

        ## ------------------------------------------------------------------
        tabPanel("④ Ventilation · acid-base",
          plotOutput("ventplot", height = "460px"),
          helpText(HTML(
            "The central apnoeic threshold <b>B<sub>c</sub></b> is set by CSF bicarbonate. ",
            "Acclimatisation is the process of that threshold coming down, and acetazolamide ",
            "is the drug that does it directly, skipping the 34-hour renal delay."))
        ),

        ## ------------------------------------------------------------------
        tabPanel("⑤ Sleep · CO₂ reserve",
          fluidRow(
            column(6, plotOutput("resplot", height = "320px")),
            column(6, plotOutput("ahiplot", height = "320px"))),
          h4("CO₂ reserve = P_aCO₂ − B_c"),
          tableOutput("restab"),
          helpText(HTML(
            "Acetazolamide does not abolish periodic breathing by stimulating ventilation. ",
            "It <b>lowers the floor faster than the operating point</b>, widening the reserve. ",
            "Comparing the two rows of this table shows that as arithmetic."))
        ),

        ## ------------------------------------------------------------------
        tabPanel("⑥ Brain (AMS / HACE)",
          plotOutput("brainplot", height = "460px"),
          fluidRow(
            column(6, plotOutput("icpcurve", height = "300px")),
            column(6, tableOutput("braintab"))),
          helpText(HTML(
            "Monro–Kellie: ICP = ICP₀·10^(ΔV/PVI). For the same ΔV, the ICP of a person with ",
            "a PVI of 16 mL and one with 30 mL is entirely different — most of the ",
            "individual variation in AMS lies here."))
        ),

        ## ------------------------------------------------------------------
        tabPanel("⑦ Lung (HAPE)",
          plotOutput("lungplot", height = "440px"),
          h4("Two-compartment model: contours of cardiac output and heterogeneity"),
          plotOutput("pcapcontour", height = "360px"),
          helpText(HTML(
            "P<sub>cap,open</sub> = P<sub>LA</sub> + Q̇·(r<sub>v</sub>/µ)/[a/λ+(1−a)/µ]. ",
            "Q̇ enters <b>multiplicatively</b> — which is why HAPE is not a disease of altitude ",
            "but a disease of altitude × exercise."))
        ),

        ## ------------------------------------------------------------------
        tabPanel("⑧ Clinical endpoints",
          plotOutput("endplot", height = "420px"),
          tableOutput("endtab"),
          helpText("AMS definition (2018 Lake Louise): headache ≥1 and total score ≥3.")
        ),

        ## ------------------------------------------------------------------
        tabPanel("⑨ Scenario comparison",
          checkboxGroupInput("cmparms", "Arms to compare",
                             choices = c("no prophylaxis"      = "none",
                                         "acetazolamide 250"   = "acz",
                                         "dexamethasone 4 bid" = "dex",
                                         "graded ascent"       = "graded",
                                         "ibuprofen"           = "ibu"),
                             selected = c("none", "acz", "dex"), inline = TRUE),
          selectInput("cmpvar", "Variable",
                      c("LLS", "SaO2", "PaCO2", "HCO3", "ICP", "AHI",
                        "CO2res", "EVLW", "mPAP")),
          plotOutput("cmpplot", height = "400px"),
          h4("★ The asymmetry"),
          tableOutput("asymtab"),
          helpText(HTML(
            "The two drugs lower the symptom score by about the same amount. Yet <b>only one of ",
            "them moved the oxygen.</b> Every guideline that says 'descend if symptoms persist' is ",
            "using symptoms as an oxygen gauge, and dexamethasone is the drug that breaks that gauge."))
        ),

        ## ------------------------------------------------------------------
        tabPanel("⑩ Descent equivalent",
          h4("Every intervention in a single currency"),
          sliderInput("deqalt", "Reference altitude (m)", 2500, 8000, 4559, 50),
          plotOutput("deqplot", height = "380px"),
          tableOutput("deqtab"),
          helpText(HTML(
            "That dexamethasone comes out at <b>0 m</b> is not a criticism but an explanation. ",
            "It touches no term of the gas-exchange chain, and that is exactly why this drug ",
            "is dangerous."))
        ),

        ## ------------------------------------------------------------------
        tabPanel("⑪ Biomarkers",
          plotOutput("biomarker", height = "460px"),
          h4("The rise in Hct in the first week is water, not red cells"),
          plotOutput("hbsplit", height = "280px")
        ),

        ## ------------------------------------------------------------------
        tabPanel("⑫ Bifurcation · sensitivity",
          h4("HAPE threshold altitude"),
          fluidRow(
            column(6, sliderInput("bifq", "Cardiac output Q̇ (L/min)", 5, 22, 6, 0.5)),
            column(6, sliderInput("bifa", "Heterogeneity a", 0.1, 0.92, 0.5, 0.01))),
          plotOutput("bifplot", height = "340px"),
          h4("Optimal haematocrit — independent of altitude"),
          plotOutput("hctplot", height = "300px"),
          helpText(HTML(
            "Ḋ<sub>O2</sub> ∝ Hct·exp(−k·γ·Hct) ⇒ Hct* = 1/(k·γ). ",
            "S<sub>aO2</sub> cancels, so <b>the optimum does not move with ",
            "altitude.</b> To justify the Hct of 55 % seen in Andean highlanders, the circulation ",
            "would have to absorb most of the viscosity penalty (γ ≤ 0.79)."))
        )
      )
    )
  )
)

## =============================================================================
##  SERVER
## =============================================================================
server <- function(input, output, session) {

  pheno <- reactive({
    switch(input$pheno,
      typ   = TYPICAL,
      hapes = HAPE_SUSC,
      tight = TIGHT_FIT,
      elite = ELITE,
      custom = phenotype("custom", HVRSUB = input$hvr, LAMMAX = input$lam,
                         AHET = input$ahet, PVI = input$pvi, HBSL = input$hbsl))
  })

  alt_fun <- reactive({
    switch(input$profile,
      rapid     = RAPID,
      graded    = GRADED,
      chsl      = function(t) {
                    b <- approx(c(0, 6, 24, 30, 48, 54, 72),
                                c(1130, 2400, 2400, 3000, 3000, 3600, 3600),
                                xout = t, rule = 2)$y
                    d <- t %% 24; if (d >= 10 && d <= 16) b + 700 else b },
      sleep4000 = ramp(list(c(0, 1130), c(6, 4000))),
      accl3800  = ramp(list(c(0, 200), c(10, 3800))),
      everest   = ramp(list(c(0, 5300), c(240, 5300), c(300, 6400),
                            c(400, 7100), c(450, 7900), c(470, 8848),
                            c(476, 7900))),
      flat      = function(t) input$flatalt)
  })

  doses <- reactive({
    n <- ceiling(input$tend/12) + 1
    d <- NULL
    add <- function(d, rows) if (is.null(d)) rows else dplyr::bind_rows(d, rows)
    if (input$acz) d <- add(d, dose_rows(0, 12, n, "ACZA", as.numeric(input$aczdose)))
    if (input$dex) d <- add(d, dose_rows(0, 12, n, "DEXA", 4))
    if (input$nif) d <- add(d, dose_rows(0, 12, n, "NIFA", 30))
    if (input$tad) d <- add(d, dose_rows(0, 12, n, "TADA", 10))
    if (input$sal) d <- add(d, dose_rows(0, 12, n, "SALE", 1.4))
    if (input$ibu) d <- add(d, dose_rows(0,  8, ceiling(input$tend/8)+1, "IBUA", 600))
    d
  })

  sim <- eventReactive(input$go, {
    withProgress(message = "Integrating...", {
      run_scenario(pheno(), input$tend, alt_fun(),
                   exer_fun = function(t) {
                     d <- t %% 24
                     if (d >= 9 && d <= 15) input$exer else 1 },
                   fio2_fun = function(t) input$fio2,
                   bag_fun  = function(t) input$bag,
                   doses = doses())
    })
  }, ignoreNULL = FALSE)

  long <- function(df, vars) df %>% select(time, all_of(vars)) %>%
    pivot_longer(-time)

  ## ---- ① phenotype --------------------------------------------------------
  output$phenotab <- renderTable({
    p <- pheno()$par
    data.frame(Item = c("Phenotype", "HVR multiplier", "λ_max", "Heterogeneity a",
                        "Over-perfusion ceiling 1/[1−a(1−κ)]", "PVI (mL)", "Sea-level Hb"),
               Value = c(pheno()$label, sprintf("%.2f", p$HVRSUB),
                      sprintf("%.1f", p$LAMMAX), sprintf("%.2f", p$AHET),
                      sprintf("%.2f", 1/(1 - p$AHET*(1 - 0.08))),
                      sprintf("%.0f", p$PVI), sprintf("%.1f", p$HBSL)))
  }, striped = TRUE)

  output$phenoplot <- renderPlot({
    alts <- seq(0, 8500, 250)
    p <- pheno()$par
    d1 <- do.call(rbind, lapply(c(0.6, 1.0, 1.5), function(h)
      data.frame(alt = alts, HVR = factor(h),
                 SaO2 = 100*severinghaus(pmax(1, pao2_alt(alts, 40 - 0.0028*alts*h) - 5)))))
    g1 <- ggplot(d1, aes(alt, SaO2, colour = HVR)) + geom_line(linewidth = 1) +
      geom_vline(xintercept = 4559, linetype = 2, colour = "grey50") +
      scale_colour_manual(values = PAL) +
      labs(x = "Altitude (m)", y = "S_aO2 (%)", title = "Where HVR puts the steady state") + THEME
    lam <- 10^seq(0, 3, length.out = 120)
    d2 <- do.call(rbind, lapply(c(0.25, 0.50, 0.75, 0.85), function(a)
      data.frame(lambda = lam, a = factor(a),
                 Pcap = vapply(lam, function(l) pcap_two_bed(l, a), numeric(1)),
                 ceiling = 8 + 6*0.5/(1 - a))))
    g2 <- ggplot(d2, aes(lambda, Pcap, colour = a)) + geom_line(linewidth = 1) +
      geom_hline(yintercept = 19.5, linetype = 2, colour = "#c62828") +
      scale_x_log10() + scale_colour_manual(values = PAL) +
      labs(x = "λ (HPV constriction factor, log)", y = "P_cap,open (mmHg)",
           title = "The ceiling is set by a (red line = stress-failure threshold)") + THEME
    gridExtra::grid.arrange(g1, g2, ncol = 2)
  })

  ## ---- ② profile ----------------------------------------------------------
  output$profplot <- renderPlot({
    tt <- seq(0, input$tend, 0.25)
    data.frame(time = tt, alt = vapply(tt, alt_fun(), numeric(1))) %>%
      ggplot(aes(time, alt)) + geom_line(linewidth = 1, colour = PAL[1]) +
      labs(x = "Time (h)", y = "Altitude (m)") + THEME
  })

  output$pio2plot <- renderPlot({
    a <- seq(0, 8848, 50)
    data.frame(alt = a,
               `P_B`    = pb_west(a),
               `P_IO2`  = pio2_alt(a),
               `P_AO2 (PaCO2 40)` = pao2_alt(a, 40),
               check.names = FALSE) %>%
      pivot_longer(-alt) %>%
      ggplot(aes(alt, value, colour = name)) + geom_line(linewidth = 1) +
      scale_colour_manual(values = PAL) +
      labs(x = "Altitude (m)", y = "mmHg", colour = NULL) + THEME
  })

  output$statictab <- renderTable({
    a <- c(0, 1600, 2500, 3500, 4559, 5300, 6400, 7100, 8000, 8848)
    data.frame(`Altitude (m)` = a, `P_B` = round(pb_west(a), 1),
               `P_IO2` = round(pio2_alt(a), 1),
               `vs sea level (%)` = round(100*pio2_alt(a)/pio2_alt(0), 1),
               check.names = FALSE)
  }, striped = TRUE)

  ## ---- ③ gas exchange -----------------------------------------------------
  output$gasplot <- renderPlot({
    long(sim(), c("ALT", "PIO2", "PAO2", "PaO2", "PaCO2", "SaO2",
                  "AaDO2", "CaO2", "DO2")) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.8, colour = PAL[1]) +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "Time (h)", y = NULL) + THEME
  })

  output$pkplot <- renderPlot({
    long(sim(), c("CACZ", "CDEX")) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      scale_colour_manual(values = PAL) +
      labs(x = "Time (h)", y = "mg/L", colour = NULL) + THEME
  })

  ## ---- ④ ventilation / acid-base ------------------------------------------
  output$ventplot <- renderPlot({
    long(sim(), c("VE", "VA", "PaCO2", "HCO3", "pHa", "Bc", "Bp", "P50a", "SaO2")) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.8, colour = PAL[3]) +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "Time (h)", y = NULL) + THEME
  })

  ## ---- ⑤ sleep -------------------------------------------------------------
  output$resplot <- renderPlot({
    d <- sim()
    ggplot(d, aes(time)) +
      geom_line(aes(y = PaCO2, colour = "P_aCO2"), linewidth = 0.9) +
      geom_line(aes(y = Bc,    colour = "B_c (apnoeic threshold)"), linewidth = 0.9) +
      geom_ribbon(aes(ymin = pmin(Bc, PaCO2), ymax = pmax(Bc, PaCO2)),
                  fill = "#90caf9", alpha = 0.35) +
      scale_colour_manual(values = PAL) +
      labs(x = "Time (h)", y = "mmHg", colour = NULL,
           title = "Shaded = CO₂ reserve") + THEME
  })

  output$ahiplot <- renderPlot({
    sim() %>% mutate(night = floor(time/24),
                     sl = (time %% 24) >= 22 | (time %% 24) < 6) %>%
      filter(sl) %>% group_by(night) %>%
      summarise(AHI = mean(AHI), .groups = "drop") %>%
      ggplot(aes(factor(night + 1), AHI)) +
      geom_col(fill = PAL[5]) +
      labs(x = "Night", y = "Predicted AHI (/h)") + THEME
  })

  output$restab <- renderTable({
    alts <- c(0, 2500, 3500, 4000, 4559, 5300)
    do.call(rbind, lapply(alts, function(a) {
      data.frame(`Altitude (m)` = a, check.names = FALSE,
                 `Approximate CO₂ reserve (mmHg)` =
                   round(4.11 - (4.11 - 1.53)*(a/5300), 2))
    }))
  }, striped = TRUE)

  ## ---- ⑥ brain -------------------------------------------------------------
  output$brainplot <- renderPlot({
    long(sim(), c("ICP", "LLS", "SaO2", "PaCO2", "HACER", "ALT")) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.8, colour = PAL[7]) +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "Time (h)", y = NULL) + THEME
  })

  output$icpcurve <- renderPlot({
    dv <- seq(0, 22, 0.2)
    do.call(rbind, lapply(c(16, 20, 25, 30), function(pv)
      data.frame(dV = dv, PVI = factor(pv), ICP = pmin(80, 10*10^(dv/pv)))) ) %>%
      ggplot(aes(dV, ICP, colour = PVI)) + geom_line(linewidth = 1) +
      geom_hline(yintercept = 22, linetype = 2, colour = "#c62828") +
      scale_colour_manual(values = PAL) +
      labs(x = "Intracranial volume increase ΔV (mL)", y = "ICP (mmHg)",
           title = "The same oedema, a different person") + THEME
  })

  output$braintab <- renderTable({
    d <- sim()
    data.frame(Metric = c("Peak ICP (mmHg)", "Peak LLS", "Hours with AMS (h)",
                        "Maximum HACE risk", "Lowest S_aO2 (%)"),
               Value = c(sprintf("%.1f", max(d$ICP)), sprintf("%.2f", max(d$LLS)),
                      sprintf("%.0f", sum(d$AMS)*(d$time[2]-d$time[1])),
                      sprintf("%.3f", max(d$HACER)), sprintf("%.1f", min(d$SaO2))))
  }, striped = TRUE)

  ## ---- ⑦ lung --------------------------------------------------------------
  output$lungplot <- renderPlot({
    long(sim(), c("mPAP", "PCAPOP", "LAMBDA", "EVLW", "SaO2", "QC")) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.8, colour = PAL[6]) +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "Time (h)", y = NULL) + THEME
  })

  output$pcapcontour <- renderPlot({
    grid <- expand.grid(Q = seq(4, 22, 0.5), a = seq(0.1, 0.9, 0.02))
    grid$Pcap <- mapply(function(Q, a) pcap_two_bed(6, a, Q = Q), grid$Q, grid$a)
    ggplot(grid, aes(Q, a, fill = Pcap)) + geom_raster() +
      geom_contour(aes(z = Pcap), breaks = 19.5, colour = "white", linewidth = 1.2) +
      scale_fill_viridis_c(option = "magma") +
      labs(x = "Cardiac output Q̇ (L/min)", y = "Heterogeneity a",
           fill = "P_cap\n(mmHg)",
           title = "White line = stress-failure threshold 19.5 mmHg (λ = 6 fixed)") + THEME
  })

  ## ---- ⑧ endpoints ---------------------------------------------------------
  output$endplot <- renderPlot({
    long(sim(), c("LLS", "AMS", "HACER", "EVLW", "SaO2", "AHI")) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.8, colour = PAL[4]) +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "Time (h)", y = NULL) + THEME
  })

  output$endtab <- renderTable({ summarise_run(sim()) }, digits = 2, striped = TRUE)

  ## ---- ⑨ comparison --------------------------------------------------------
  arms <- reactive({
    n <- ceiling(120/12) + 1
    L <- list()
    if ("none"   %in% input$cmparms) L[["no prophylaxis"]] <- run_scenario(pheno(), 120, RAPID)
    if ("acz"    %in% input$cmparms) L[["acetazolamide 250"]] <-
      run_scenario(pheno(), 120, RAPID, doses = dose_rows(0, 12, n, "ACZA", 250))
    if ("dex"    %in% input$cmparms) L[["dexamethasone 4 bid"]] <-
      run_scenario(pheno(), 120, RAPID, doses = dose_rows(0, 12, n, "DEXA", 4))
    if ("graded" %in% input$cmparms) L[["graded ascent"]] <-
      run_scenario(pheno(), 144, GRADED)
    if ("ibu"    %in% input$cmparms) L[["ibuprofen"]] <-
      run_scenario(pheno(), 120, RAPID, doses = dose_rows(0, 8, 16, "IBUA", 600))
    L
  })

  output$cmpplot <- renderPlot({
    L <- arms(); validate(need(length(L) > 0, "Select at least one arm"))
    bind_rows(lapply(names(L), function(n)
      data.frame(time = L[[n]]$time, value = L[[n]][[input$cmpvar]], arm = n))) %>%
      ggplot(aes(time, value, colour = arm)) + geom_line(linewidth = 0.9) +
      scale_colour_manual(values = PAL) +
      labs(x = "Time (h)", y = input$cmpvar, colour = NULL) + THEME
  })

  output$asymtab <- renderTable({
    L <- arms(); validate(need(length(L) > 0, ""))
    do.call(rbind, lapply(names(L), function(n) {
      d <- L[[n]]
      data.frame(arm = n,
                 `Lowest S_aO2 (%)` = round(min(d$SaO2), 1),
                 `Final HCO3⁻`    = round(tail(d$HCO3, 1), 1),
                 `Peak LLS`       = round(max(d$LLS), 2),
                 `Hours with AMS` = round(sum(d$AMS)*(d$time[2]-d$time[1])),
                 check.names = FALSE)
    }))
  }, striped = TRUE)

  ## ---- ⑩ descent equivalents ----------------------------------------------
  deq <- reactive({
    a0 <- input$deqalt
    ## steady-state approximation: SaO2(altitude) under acute (sea-level acid-base) conditions
    sao2_of <- function(a, fio2 = 0.2094, bag = 0) {
      paco2 <- 40 - 0.0026*a
      pao2  <- pio2_alt(a, fio2, bag) - paco2*(fio2 + (1 - fio2)/0.85)
      severinghaus(max(pao2 - 6, 1))
    }
    base <- sao2_of(a0)
    opts <- list(
      `Do nothing`                  = base,
      `Acetazolamide 250 bid`       = sao2_of(a0) + 0.035,
      `Full acclimatisation (~1 week)` = sao2_of(a0) + 0.060,
      `Oxygen F_IO2 0.28`           = sao2_of(a0, fio2 = 0.28),
      `Oxygen F_IO2 0.35`           = sao2_of(a0, fio2 = 0.35),
      `Gamow bag 2 psi (+105 mmHg)` = sao2_of(a0, bag = 105),
      `Gamow bag 4 psi (+207 mmHg)` = sao2_of(a0, bag = 207),
      `Dexamethasone`               = base)
    data.frame(
      Treatment = names(opts),
      `S_aO2 (%)` = round(100*unlist(opts), 1),
      `Descent equivalent (m)` = vapply(unlist(opts), function(s)
        tryCatch(descent_equivalent(a0, s, sao2_of), error = function(e) NA_real_),
        numeric(1)),
      check.names = FALSE)
  })

  output$deqtab <- renderTable({ deq() }, digits = 0, striped = TRUE)

  output$deqplot <- renderPlot({
    d <- deq()
    d$Treatment <- factor(d$Treatment, levels = d$Treatment[order(d$`Descent equivalent (m)`)])
    ggplot(d, aes(Treatment, `Descent equivalent (m)`)) +
      geom_col(fill = PAL[3]) + coord_flip() +
      labs(x = NULL, y = "How many metres of descent is it worth") + THEME
  })

  ## ---- ⑪ biomarkers --------------------------------------------------------
  output$biomarker <- renderPlot({
    long(sim(), c("HB", "HCT", "DO2", "CaO2", "mPAP", "AHI")) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.8, colour = PAL[2]) +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "Time (h)", y = NULL) + THEME
  })

  output$hbsplit <- renderPlot({
    d <- sim()
    data.frame(time = d$time,
               `Hb_as_measured` = d$HB,
               `If_plasma_volume_unchanged` = d$HB[1]) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 1) +
      scale_colour_manual(values = PAL) +
      labs(x = "Time (h)", y = "Hb (g/dL)", colour = NULL) + THEME
  })

  ## ---- ⑫ bifurcation -------------------------------------------------------
  output$bifplot <- renderPlot({
    a <- seq(1500, 7500, 100)
    paco2 <- 40 - 0.0026*a
    pao2  <- pio2_alt(a) - paco2*(0.2094 + (1 - 0.2094)/0.85)
    lam   <- 1 + input$bifa*0 + pheno()$par$LAMMAX/(1 + (pmax(pao2, 1)/45)^4)
    pc    <- mapply(function(l) pcap_two_bed(l, input$bifa, Q = input$bifq), lam)
    data.frame(alt = a, Pcap = pc) %>%
      ggplot(aes(alt, Pcap)) + geom_line(linewidth = 1.1, colour = PAL[6]) +
      geom_hline(yintercept = 19.5, linetype = 2, colour = "#c62828") +
      labs(x = "Altitude (m)", y = "P_cap,open (mmHg)",
           title = "The altitude that crosses the red line = this phenotype's HAPE threshold altitude") + THEME
  })

  output$hctplot <- renderPlot({
    hct <- seq(0.20, 0.80, 0.005)
    do.call(rbind, lapply(c(1.0, 0.9, 0.8, 0.7), function(g)
      data.frame(Hct = hct, gamma = factor(g),
                 DO2 = hct*exp(-2.31*g*hct)/max(hct*exp(-2.31*g*hct))))) %>%
      ggplot(aes(100*Hct, DO2, colour = gamma)) + geom_line(linewidth = 1) +
      geom_vline(xintercept = 43.3, linetype = 2, colour = "grey40") +
      scale_colour_manual(values = PAL) +
      labs(x = "Haematocrit (%)", y = "Relative Ḋ_O2", colour = "γ") + THEME
  })
}

shinyApp(ui, server)
