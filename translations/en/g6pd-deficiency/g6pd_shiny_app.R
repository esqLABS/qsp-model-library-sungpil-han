## =====================================================================
##  G6PD DEFICIENCY — interactive QSP dashboard (Shiny)
##  Glucose-6-phosphate dehydrogenase deficiency · interactive simulator
##
##  Run with:
##      shiny::runApp("g6pd_shiny_app.R")
##  (the app sources g6pd_mrgsolve_model.R from the same directory)
##
##  The app is organised around ONE control that matters more than all the
##  others put together: the variant selector, which sets exactly two
##  numbers, E0 and TAU. Tab 1 draws the curve those two numbers make.
##  Every other tab is a consequence of it.
##
##  Tabs
##    1  Patient & variant              the age-activity curve, and a*
##    2  Pharmacokinetics               drug concentrations, oxidant flux
##    3  Erythrocyte age structure      the histogram with a* drawn on it
##    4  Haematology                    Hb, reticulocytes, EPO
##    5  Scenario comparison            side-by-side, incl. the 1-vs-2 case
##    6  Biomarkers                     MetHb, bilirubin, haptoglobin, AKI
##    7  The assay trap                 what the laboratory would report
##    8  Neonatal jaundice              the G6PD x UGT1A1 product
## =====================================================================

library(shiny)
library(mrgsolve)
library(ggplot2)
library(dplyr)
library(tidyr)

## NOTE ON MASKING. mrgsolve exports req(), param(), plot() and filter(), and
## it is attached AFTER shiny, so a bare req() in the server would resolve to
## mrgsolve::req and fail at runtime (a parse check does not catch this).
## dplyr is attached after mrgsolve, so filter() correctly resolves to dplyr.
## Anything ambiguous below is namespace-qualified.

source("g6pd_mrgsolve_model.R", local = TRUE)

THEME <- theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey93"),
        legend.position = "bottom")

VARIANT_CHOICES <- setNames(names(VARIANTS),
                            vapply(VARIANTS, function(v) v$label, ""))

DRUG_CHOICES <- c("none"                             = "none",
                  "primaquine daily"                 = "pq",
                  "primaquine weekly"                = "pqw",
                  "tafenoquine single dose"          = "tq",
                  "dapsone"                          = "dap",
                  "rasburicase (TLS)"                = "rbx",
                  "fava bean meal"                   = "fava")

## ---------------------------------------------------------------------
##  UI
## ---------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("G6PD deficiency QSP simulator — the enzyme is not a number but a function of cell age"),
  tags$p(style = "color:#666; margin-top:-8px;",
         "Glucose-6-phosphate dehydrogenase deficiency · age-structured ",
         "erythrocyte redox model · for education and research only, not for clinical use"),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      h4("1. Genotype (the whole disease)"),
      selectInput("variant", "G6PD variant", VARIANT_CHOICES, selected = "Aminus"),
      helpText(HTML("What the variant changes is <b>exactly two numbers</b>: the activity E0
                     of a new erythrocyte, and TAU, the half-life over which the enzyme
                     disappears inside the cell. Everything else is derived from these.")),
      checkboxInput("custom", "Specify E0 / TAU directly", FALSE),
      conditionalPanel("input.custom",
        sliderInput("E0",  "E0 (activity of a new erythrocyte, normal=1)", 0.005, 1.0, 0.55, 0.005),
        sliderInput("TAU", "TAU (intracellular enzyme half-life, days)", 3, 80, 13, 1)),

      hr(),
      h4("2. Patient"),
      numericInput("wt", "Body weight (kg)", 70, 3, 150, 1),
      selectInput("ugt", "UGT1A1 promoter (TA)n", c("6/6", "6/7", "7/7"), "6/6"),
      sliderInput("cyp2d6", "CYP2D6 activity score", 0, 2, 1, 0.25),
      selectInput("nat2", "NAT2 acetylation", c("fast", "slow"), "fast"),
      sliderInput("marrow", "Marrow reserve (1=normal, 0.15=aplastic crisis)",
                  0.05, 1.5, 1.0, 0.05),
      checkboxInput("splenx", "Splenectomised", FALSE),
      checkboxInput("infect", "Concurrent infection (adds oxidative stress)", FALSE),

      hr(),
      h4("3. Exposure"),
      selectInput("drug", "Drug / trigger", DRUG_CHOICES, selected = "pq"),
      conditionalPanel("input.drug == 'pq'",
        sliderInput("pq_dose", "Primaquine (mg/day)", 7.5, 60, 30, 7.5),
        sliderInput("pq_days", "Duration of dosing (days)", 1, 60, 60, 1)),
      conditionalPanel("input.drug == 'pqw'",
        sliderInput("pqw_dose", "Primaquine (mg, once weekly)", 15, 60, 45, 7.5),
        sliderInput("pqw_wk", "Number of weeks", 1, 12, 8, 1)),
      conditionalPanel("input.drug == 'tq'",
        sliderInput("tq_dose", "Tafenoquine single dose (mg)", 50, 600, 300, 50)),
      conditionalPanel("input.drug == 'dap'",
        sliderInput("dap_dose", "Dapsone (mg/day)", 25, 300, 100, 25),
        checkboxInput("mb", "Methylene blue 2 mg/kg given (day 30)", FALSE)),
      conditionalPanel("input.drug == 'rbx'",
        sliderInput("urate0", "Baseline urate (mg/dL)", 4, 30, 20, 1),
        sliderInput("rbx_days", "Number of rasburicase dosing days", 1, 5, 3, 1)),
      conditionalPanel("input.drug == 'fava'",
        sliderInput("fava_amt", "Absorbable aglycone (mg)", 100, 3000, 1000, 100)),

      hr(),
      sliderInput("days", "Observation period (days)", 20, 180, 60, 5),
      actionButton("go", "Run simulation", class = "btn-primary btn-block")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs", type = "tabs",

        ## ---------------------------------------------------------------
        tabPanel("1 · Patient & variant",
          br(),
          fluidRow(
            column(7, h4("E(a) = E0 · exp(−ln2 · a / TAU)"),
                      plotOutput("p_curve", height = 330)),
            column(5, h4("Critical age a*"),
                      plotOutput("p_astar", height = 330))),
          hr(),
          h4("These two curves are the whole app"),
          HTML("<p>Left: mature erythrocytes have no nucleus and no ribosomes, so they cannot
                make new enzyme. It is loaded once and thereafter only <b>declines</b>. Activity
                is therefore an exponential function of cell age, and a variant is that <b>slope</b>.</p>
                <p>Right: given an oxidant load OX, a cell survives only while
                <code>VMAXOX·E(a) &gt; OX</code>. That boundary is the critical
                age <b>a*</b>, and it solves in closed form:</p>
                <p style='margin-left:2em'><code>a* = (TAU / ln2) ·
                ln(E0 · VMAXOX / OX<sub>eff</sub>)</code>&nbsp;&nbsp;&nbsp;
                <small>(OX<sub>eff</sub> = OX − KPIT·HZ50, that is OX minus the small load
                that can be borne without triggering removal)</small></p>
                <p>The point is that <b>there is no fitted haemolysis rate constant</b> in
                this expression. What it contains is the genotype (E0, TAU), the biochemical
                capacity (VMAXOX), and the exposure (OX) — nothing else.</p>
                <p>Haemolysis is not a rate constant multiplying a pool but a <b>blade
                passing through the age histogram</b>. In A− the slope is steep, so a*
                falls <i>inside</i> the distribution and leaves young, resistant survivors —
                which is why recovery happens even while the drug is continued. In
                Mediterranean the curve is flat and low, and there are no survivors to leave.</p>"),
          hr(),
          tableOutput("t_variants")),

        ## ---------------------------------------------------------------
        tabPanel("2 · Pharmacokinetics",
          br(),
          fluidRow(column(6, plotOutput("p_pk", height = 300)),
                   column(6, plotOutput("p_ox", height = 300))),
          hr(),
          HTML("<p><b>The half-life is the toxicology.</b> Primaquine has a half-life of about
                7 hours, so stopping it actually works once haemolysis begins. Tafenoquine has
                one of about 15 days, so there is <i>nothing to stop</i> — and the reticulocyte
                rescue call takes 4–7 days to arrive. This is why the same total exposure gives
                a different haemoglobin nadir depending on whether it is delivered by a drug
                that switches off quickly or by one that switches off slowly, and why the FDA
                requires a quantitative assay for tafenoquine alone.</p>"),
          plotOutput("p_astar_t", height = 280)),

        ## ---------------------------------------------------------------
        tabPanel("3 · Erythrocyte age structure",
          br(),
          sliderInput("snap", "Time point observed (days)", 0, 60, 8, 1, width = "70%"),
          plotOutput("p_hist", height = 340),
          hr(),
          plotOutput("p_heat", height = 300),
          HTML("<p>Above: the erythrocyte age histogram at that time point. The red vertical
                line is a*, and the cells to its right are the ones being destroyed now.
                Below: how the cell count in each age bin moves on the time × age plane. In A−
                only the older bins empty while the young end refills.</p>")),

        ## ---------------------------------------------------------------
        tabPanel("4 · Haematology",
          br(),
          fluidRow(column(6, plotOutput("p_hb",   height = 300)),
                   column(6, plotOutput("p_ret",  height = 300))),
          fluidRow(column(6, plotOutput("p_epo",  height = 280)),
                   column(6, plotOutput("p_lys",  height = 280))),
          hr(), tableOutput("t_summary")),

        ## ---------------------------------------------------------------
        tabPanel("5 · Scenario comparison",
          br(),
          checkboxGroupInput("cmp", "Variants to compare (the same exposure applied to each)",
            choices = VARIANT_CHOICES,
            selected = c("normal", "Aminus", "Mediterranean"), inline = TRUE),
          actionButton("go_cmp", "Run comparison", class = "btn-primary"),
          br(), br(),
          plotOutput("p_cmp", height = 380),
          hr(),
          tableOutput("t_cmp"),
          HTML("<p><b>This tab is the model's central result.</b> The same drug, the same
                dose, and on both the same label of 'severe deficiency'. A− hits bottom and
                then climbs back <i>while the drug is still being taken</i>; Mediterranean does not.
                The only parameter that makes the difference is the decay constant.</p>")),

        ## ---------------------------------------------------------------
        tabPanel("6 · Biomarkers",
          br(),
          fluidRow(column(6, plotOutput("p_met",  height = 280)),
                   column(6, plotOutput("p_bili", height = 280))),
          fluidRow(column(6, plotOutput("p_hapto",height = 280)),
                   column(6, plotOutput("p_aki",  height = 280))),
          hr(),
          HTML("<p><b>Why methylene blue is contraindicated is not a fact to memorise but
                the result of a calculation.</b> One NADPH pool feeds two consumers:
                glutathione reductase (oxidative defence) and NADPH-methaemoglobin reductase
                (the route converting methylene blue to leucomethylene blue). In deficiency
                the second does not turn, so methylene blue <i>fails as an antidote and
                becomes a new oxidant at once</i>. Turn it on in the dapsone tab and compare normal with A−.</p>")),

        ## ---------------------------------------------------------------
        tabPanel("7 · The assay trap",
          br(),
          plotOutput("p_assay", height = 360),
          hr(),
          HTML("<p>A quantitative G6PD assay reports the <b>age-weighted mean over all</b>
                circulating erythrocytes. But haemolysis has just erased the old, enzyme-free
                cells, and their place has been taken by reticulocytes carrying
                <b>1.5 times the activity of a young cell</b>. So the assay value
                <b>rises</b> during and immediately after haemolysis.</p>
                <p>In the plot above a genuine A− patient reads within the normal range near
                the nadir. A normal result during an acute episode therefore cannot exclude
                deficiency, and the retest must be done <b>three months later</b> — one
                erythrocyte lifespan.</p>"),
          tableOutput("t_assay")),

        ## ---------------------------------------------------------------
        tabPanel("8 · Neonatal jaundice",
          br(),
          p("A 3 kg term neonate, from 12 hours of age for 12 days. G6PD and UGT1A1 act on",
            "production (the numerator) and conjugation (the denominator) respectively, so their",
            "effects multiply rather than add. All four combinations are run together."),
          actionButton("go_neo", "Simulate the four neonatal groups", class = "btn-primary"),
          br(), br(),
          plotOutput("p_neo", height = 380),
          hr(),
          tableOutput("t_neo"),
          HTML("<p>With only one of the two present the curve stops near the phototherapy
                threshold; only when both are present does it cross into exchange-transfusion
                territory. G6PD deficiency is a leading cause of kernicterus worldwide not
                because the deficiency itself is especially powerful but because it <b>multiplies with a common second defect</b>.</p>"))
      )
    )
  )
)

## ---------------------------------------------------------------------
##  SERVER
## ---------------------------------------------------------------------
server <- function(input, output, session) {

  patient <- reactive({
    p <- g6pd_patient(variant = input$variant, wt = input$wt, ugt = input$ugt,
                      cyp2d6 = input$cyp2d6, nat2 = input$nat2,
                      spleen = if (input$splenx) 0 else 1,
                      marrow = input$marrow)
    if (isTRUE(input$custom)) { p$E0 <- input$E0; p$TAU <- input$TAU }
    if (isTRUE(input$infect)) p$INFON <- 1
    p
  })

  regimen <- reactive({
    switch(input$drug,
      none  = NULL,
      pq    = rx_primaquine(input$pq_dose, days = input$pq_days),
      pqw   = rx_pq_weekly(input$pqw_dose, weeks = input$pqw_wk),
      tq    = rx_tafenoquine(input$tq_dose),
      dap   = if (isTRUE(input$mb))
                c(rx_dapsone(input$dap_dose, days = input$days),
                  rx_methyleneblue(2 * input$wt, start = 30))
              else rx_dapsone(input$dap_dose, days = input$days),
      rbx   = rx_rasburicase(0.2 * input$wt, days = input$rbx_days),
      fava  = rx_fava(input$fava_amt))
  })

  extra <- reactive({
    if (input$drug == "rbx") {
      ## urate mg/dL -> mmol in a 0.5 L/kg volume
      vur <- 0.5 * input$wt
      list(VUR = vur, URATE_0 = input$urate0 * 10 * vur / 168.1,
           KGENUR = 3 * input$wt / 70 * (input$urate0 / 6))
    } else list()
  })

  sim <- eventReactive(input$go, {
    withProgress(message = "Simulating...", {
      run_g6pd(patient = patient(), rx = regimen(), days = input$days,
               extra = extra(), delta = 0.05,
               name = VARIANTS[[input$variant]]$label)
    })
  }, ignoreNULL = FALSE)

  observe({
    updateSliderInput(session, "snap", max = input$days,
                      value = min(8, input$days))
  })

  ## ---- tab 1 -------------------------------------------------------
  output$p_curve <- renderPlot({
    a <- seq(0, 120, 1)
    d <- do.call(rbind, lapply(names(VARIANTS), function(v) {
      p <- VARIANTS[[v]]
      data.frame(age = a, E = p$E0 * exp(-log(2) * a / p$TAU),
                 variant = p$label, sel = (v == input$variant))
    }))
    ggplot(d, aes(age, E, colour = variant, linewidth = sel)) +
      geom_line() +
      scale_linewidth_manual(values = c(`FALSE` = 0.5, `TRUE` = 1.8),
                             guide = "none") +
      scale_y_log10(labels = scales::percent) +
      labs(x = "Erythrocyte age (days)", y = "G6PD activity (vs a normal young cell, log)",
           colour = NULL,
           title = "Variant = the slope of this curve") +
      THEME
  })

  output$p_astar <- renderPlot({
    p  <- patient()
    ox <- 10^seq(-3, 2, length.out = 400)
    ## the same closed form the model reports as ASTAR, offset included
    oxe <- ox - 0.5 * 0.10
    astar <- ifelse(oxe <= 0, 120,
                    pmin(120, pmax(0, (p$TAU / log(2)) *
                                     log(p$E0 * 60 / pmax(oxe, 1e-12)))))
    ggplot(data.frame(ox, astar), aes(ox, astar)) +
      geom_line(linewidth = 1.2, colour = "#a02c2c") +
      geom_hline(yintercept = 120, linetype = 3) +
      annotate("text", x = 0.0015, y = 116, hjust = 0, size = 3.4,
               label = "No cells at risk") +
      scale_x_log10() +
      labs(x = "Oxidant load OX (mmol H2O2 / L RBC / day, log)",
           y = "Critical age a* (days)",
           title = "Cells older than a* haemolyse") +
      THEME
  })

  output$t_variants <- renderTable(g6pd_age_table(), digits = 3)

  ## ---- tab 2 -------------------------------------------------------
  output$p_pk <- renderPlot({
    d <- sim() %>% filter(time >= 0) %>%
      select(time, CPQ_o, CTQ_o, CDP_o) %>%
      pivot_longer(-time) %>% filter(value > 1e-8)
    if (!nrow(d)) return(ggplot() + THEME + labs(title = "No drug"))
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      scale_colour_manual(values = c(CPQ_o = "#2f6fb5", CTQ_o = "#a02c2c",
                                     CDP_o = "#2c7a2c"),
        labels = c(CPQ_o = "Primaquine", CTQ_o = "Tafenoquine", CDP_o = "Dapsone")) +
      labs(x = "Time (days)", y = "Plasma concentration (mg/L)", colour = NULL,
           title = "Pharmacokinetics") + THEME
  })

  output$p_ox <- renderPlot({
    ggplot(filter(sim(), time >= 0), aes(time, OXFLUX)) +
      geom_line(colour = "#b5651d", linewidth = 1) +
      labs(x = "Time (days)", y = "OX (mmol H2O2 / L RBC / day)",
           title = "Total oxidant load — the one point where every trigger meets") + THEME
  })

  output$p_astar_t <- renderPlot({
    ggplot(filter(sim(), time >= 0), aes(time, ASTAR)) +
      geom_line(colour = "#a02c2c", linewidth = 1.1) +
      geom_hline(yintercept = 120, linetype = 3) +
      labs(x = "Time (days)", y = "Critical age a* (days)",
           title = "How deeply a* sweeps through the age distribution") +
      ylim(0, 125) + THEME
  })

  ## ---- tab 3 -------------------------------------------------------
  age_long <- reactive({
    s <- sim()
    rn <- paste0("R", 1:12)
    s %>% select(time, all_of(rn)) %>%
      pivot_longer(-time, names_to = "bin", values_to = "cells") %>%
      mutate(bin = as.integer(sub("R", "", bin)),
             age = 10 * (bin - 0.5))
  })

  output$p_hist <- renderPlot({
    s  <- sim()
    tt <- s$time[which.min(abs(s$time - input$snap))]
    d  <- filter(age_long(), abs(time - tt) < 1e-6)
    as <- s$ASTAR[which.min(abs(s$time - tt))]
    ggplot(d, aes(age, cells)) +
      geom_col(fill = "#a9cdee", colour = "#25588c", width = 9) +
      geom_vline(xintercept = as, colour = "#a02c2c", linewidth = 1.3) +
      annotate("text", x = as, y = Inf, vjust = 1.6, hjust = -0.08,
               colour = "#a02c2c", size = 4.2,
               label = sprintf("a* = %.0f days", as)) +
      labs(x = "Erythrocyte age (days)", y = "Cell count (10^12/L)",
           title = sprintf("Age histogram, day %.0f", tt)) + THEME
  })

  output$p_heat <- renderPlot({
    ggplot(filter(age_long(), time >= 0), aes(time, age, fill = cells)) +
      geom_raster() +
      geom_line(data = filter(sim(), time >= 0),
                aes(time, ASTAR), inherit.aes = FALSE,
                colour = "white", linewidth = 1.1) +
      scale_fill_viridis_c(option = "mako", direction = -1) +
      labs(x = "Time (days)", y = "Erythrocyte age (days)", fill = "10^12/L",
           title = "The white line is a* — the emptying above it is haemolysis") + THEME
  })

  ## ---- tab 4 -------------------------------------------------------
  output$p_hb <- renderPlot({
    s <- sim(); b <- s$HB[which.min(abs(s$time))]
    ggplot(s, aes(time, HB)) +
      geom_vline(xintercept = 0, linetype = 3) +
      geom_hline(yintercept = b * 0.75, linetype = 2, colour = "#a02c2c") +
      annotate("text", x = max(s$time), y = b * 0.75, vjust = -0.5, hjust = 1,
               size = 3.3, colour = "#a02c2c", label = "-25% from baseline") +
      geom_line(linewidth = 1.1, colour = "#25588c") +
      labs(x = "Time (days)", y = "Haemoglobin (g/dL)", title = "Haemoglobin") + THEME
  })

  output$p_ret <- renderPlot({
    ggplot(sim(), aes(time, RETPCT)) +
      geom_line(linewidth = 1.1, colour = "#2c7a2c") +
      labs(x = "Time (days)", y = "Reticulocytes (%)",
           title = "The rescue call — 4–7 days to arrive") + THEME
  })

  output$p_epo <- renderPlot({
    ggplot(sim(), aes(time, EPOOUT)) +
      geom_line(linewidth = 1, colour = "#6b4a99") +
      scale_y_log10() +
      labs(x = "Time (days)", y = "EPO (mIU/mL, log)", title = "Erythropoietin") +
      THEME
  })

  output$p_lys <- renderPlot({
    ggplot(filter(sim(), time >= 0), aes(time, ATRISK)) +
      geom_area(fill = "#f2b3ae", colour = "#963029") +
      labs(x = "Time (days)", y = "Erythrocytes older than a* (%)",
           title = "The erythrocyte mass now at risk") + THEME
  })

  output$t_summary <- renderTable(summarise_g6pd(sim()))

  ## ---- tab 5 -------------------------------------------------------
  cmp <- eventReactive(input$go_cmp, {
    shiny::req(length(input$cmp) > 0)
    withProgress(message = "Comparing variants...", {
      do.call(rbind, lapply(input$cmp, function(v) {
        p <- g6pd_patient(v, wt = input$wt, ugt = input$ugt,
                          cyp2d6 = input$cyp2d6, nat2 = input$nat2,
                          spleen = if (input$splenx) 0 else 1,
                          marrow = input$marrow)
        if (isTRUE(input$infect)) p$INFON <- 1
        run_g6pd(patient = p, rx = regimen(), days = input$days,
                 extra = extra(), delta = 0.05,
                 name = VARIANTS[[v]]$label)
      }))
    })
  })

  output$p_cmp <- renderPlot({
    d <- cmp()
    d %>% group_by(scenario) %>%
      mutate(pct = 100 * HB / HB[which.min(abs(time))]) %>%
      ggplot(aes(time, pct, colour = scenario)) +
      geom_hline(yintercept = 100, linetype = 3) +
      geom_hline(yintercept = 75, linetype = 2, colour = "#a02c2c") +
      geom_line(linewidth = 1.1) +
      labs(x = "Time (days)", y = "Haemoglobin (% of baseline)", colour = NULL,
           title = "Same drug, same dose — different decay constant") + THEME
  })

  output$t_cmp <- renderTable({
    d <- cmp()
    do.call(rbind, lapply(split(d, d$scenario), summarise_g6pd))
  })

  ## ---- tab 6 -------------------------------------------------------
  output$p_met <- renderPlot({
    ggplot(sim(), aes(time, METPCT)) +
      geom_hline(yintercept = c(15, 30), linetype = 2, colour = "#a02c2c") +
      geom_line(linewidth = 1, colour = "#2f4a91") +
      labs(x = "Time (days)", y = "Methaemoglobin (%)",
           title = "MetHb — 15% cyanosis, 30% symptoms") + THEME
  })

  output$p_bili <- renderPlot({
    ggplot(sim(), aes(time, TSB)) +
      geom_line(linewidth = 1, colour = "#8f7212") +
      labs(x = "Time (days)", y = "Total bilirubin (mg/dL)", title = "Bilirubin") + THEME
  })

  output$p_hapto <- renderPlot({
    s <- sim()
    ggplot(s, aes(time)) +
      geom_line(aes(y = HPOUT, colour = "Haptoglobin (g/L)"), linewidth = 1) +
      geom_line(aes(y = FHBOUT * 5, colour = "Free Hb x5 (g/L)"), linewidth = 1) +
      scale_colour_manual(values = c("#2b6f68", "#963029")) +
      labs(x = "Time (days)", y = "Concentration", colour = NULL,
           title = "Intravascular haemolysis — haptoglobin reaches zero first") + THEME
  })

  output$p_aki <- renderPlot({
    ggplot(sim(), aes(time, CREA)) +
      geom_line(linewidth = 1, colour = "#963029") +
      labs(x = "Time (days)", y = "Serum creatinine (mg/dL)",
           title = "Tubular injury from filtered free haemoglobin") + THEME
  })

  ## ---- tab 7 -------------------------------------------------------
  output$p_assay <- renderPlot({
    s <- sim()
    ggplot(s, aes(time)) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = 70, ymax = Inf,
               fill = "#e8f6e8") +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = 30, ymax = 70,
               fill = "#fdf6ee") +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 30,
               fill = "#fbe0e0") +
      geom_hline(yintercept = c(30, 70), linetype = 2) +
      geom_line(aes(y = G6PDPCT), linewidth = 1.3, colour = "#8c4c11") +
      geom_line(aes(y = 100 * HB / max(HB)), linewidth = 0.8,
                linetype = 2, colour = "#25588c") +
      annotate("text", x = max(s$time), y = 74, hjust = 1, size = 3.4,
               label = "70% — tafenoquine threshold") +
      annotate("text", x = max(s$time), y = 33, hjust = 1, size = 3.4,
               label = "30% — daily primaquine threshold") +
      labs(x = "Time (days)",
           y = "G6PD activity the assay would report (% of normal) · dashed is haemoglobin (relative)",
           title = "The solid line is the assay value — it rises during haemolysis") + THEME
  })

  output$t_assay <- renderTable({
    s <- filter(sim(), time >= 0)
    nadir_t <- s$time[which.min(s$HB)]
    j14 <- which.min(abs(s$time - (nadir_t + 14)))
    ## column names are given as STRINGS, not identifiers: R will not parse a
    ## Hangul identifier under a C locale, which is what a bare container has
    data.frame(
      "Timepoint" = c("Before exposure", "Haemoglobin nadir", "Nadir +14 days", "End of observation"),
      "Day"       = round(c(0, nadir_t, nadir_t + 14, max(s$time)), 1),
      "Assay (%)" = round(c(s$G6PDPCT[1], s$G6PDPCT[which.min(s$HB)],
                              s$G6PDPCT[j14], tail(s$G6PDPCT, 1)), 1),
      "Hb"        = round(c(s$HB[1], min(s$HB), s$HB[j14],
                              tail(s$HB, 1)), 2),
      check.names = FALSE)
  })

  ## ---- tab 8 -------------------------------------------------------
  neo <- eventReactive(input$go_neo, {
    withProgress(message = "Simulating the four neonatal groups...", {
      grid <- list(c("normal", "6/6", "Neither"),
                   c("Mediterranean", "6/6", "G6PD deficiency only"),
                   c("normal", "7/7", "UGT1A1 7/7 only"),
                   c("Mediterranean", "7/7", "Both"))
      do.call(rbind, lapply(grid, function(g)
        run_g6pd(patient = g6pd_patient(g[1], wt = 3, ugt = g[2],
                                        neonate = TRUE, pna = 0.5),
                 rx = NULL, days = 12, lead_in = 0, delta = 0.02,
                 name = g[3])))
    })
  })

  output$p_neo <- renderPlot({
    ggplot(neo(), aes(time, TSB, colour = scenario)) +
      geom_hline(yintercept = 15, linetype = 2, colour = "#c9a227") +
      geom_hline(yintercept = 20, linetype = 2, colour = "#a02c2c") +
      annotate("text", x = 0.2, y = 15.6, hjust = 0, size = 3.4,
               colour = "#8f7212", label = "Phototherapy threshold (approx.)") +
      annotate("text", x = 0.2, y = 20.6, hjust = 0, size = 3.4,
               colour = "#a02c2c", label = "Exchange-transfusion threshold (approx.)") +
      geom_line(linewidth = 1.2) +
      labs(x = "Days of age", y = "Total serum bilirubin (mg/dL)", colour = NULL,
           title = "G6PD × UGT1A1 — they multiply, they do not add") + THEME
  })

  output$t_neo <- renderTable({
    d <- neo()
    do.call(rbind, lapply(split(d, d$scenario), function(x)
      data.frame("Group"              = x$scenario[1],
                 "Peak TSB (mg/dL)"   = round(max(x$TSB), 2),
                 "Peak day"           = round(x$time[which.max(x$TSB)], 1),
                 "Peak free bilirubin (nM)" = round(max(x$BFREE), 1),
                 "Nadir Hb (g/dL)"    = round(min(x$HB), 2),
                 check.names = FALSE)))
  })
}

shinyApp(ui, server)
