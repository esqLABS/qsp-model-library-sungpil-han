## =============================================================================
##  nb_shiny_app.R — High-Risk Neuroblastoma QSP explorer
##
##  Nine tabs, organised so that the model's two structural claims are the first
##  things a user can interrogate:
##
##    1  Patient & regimen                who is being treated, and with what
##    2  Pharmacokinetics                 four drug classes, four PK shapes
##    3  The two Fc arms                  the core: ADCC n=1 vs CDC n=2
##    4  Delivery barriers                permeability vs dose, side by side
##    5  Tumour burden & PD               four pools, four drug sensitivities
##    6  Clinical endpoints               relapse, pain, ANC/PLT, hearing
##    7  Scenario comparison              the nine derived comparisons
##    8  MIBG dosimetry                   transporter competition + dosimetry
##    9  Virtual cohort (EFS)             ANBL0032 and HR-NBL1 reproduction
##
##  Run:  shiny::runApp("nb_shiny_app.R")
##  Requires: shiny, mrgsolve, dplyr, tidyr, ggplot2, DT
##  The model itself lives in nb_mrgsolve_model.R and is sourced below.
## =============================================================================

suppressPackageStartupMessages({
  library(shiny); library(mrgsolve); library(dplyr); library(tidyr)
  library(ggplot2); library(DT)
})

source("nb_mrgsolve_model.R")

THEME <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")

## a colour per target compartment — used consistently in every tab so that
## "marrow", "nerve" and "solid tumour" always mean the same colour
CMP <- c("solid tumour" = "#C0392B", "bone marrow" = "#1E8449",
         "DRG / nerve" = "#CA6F1E", "plasma" = "#2874A6")

# =============================================================================
ui <- fluidPage(
  titlePanel("High-Risk Neuroblastoma QSP model"),
  tags$p(style = "color:#555;margin-top:-10px",
         HTML("Two Fc effector arms (ADCC ∝ bound IgG¹ · CDC ∝ bound IgG²) and three
               permeabilities (solid tumour : marrow : dorsal root ganglion ≈ 1 : 15 : 200) are the whole of this model.
               <b>The whole model is two Fc arms and three permeabilities.</b>")),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Patient"),
      sliderInput("BSA", "Body surface area BSA (m²)", 0.4, 1.6, 0.8, 0.05),
      sliderInput("WT", "Body weight (kg)", 8, 60, 20, 1),
      sliderInput("tp0", "Initial solid burden (g)", 1, 400, 100, 1),
      sliderInput("tm0", "Initial marrow MRD (10⁹ cells)", 0, 30, 5, 0.5),
      checkboxInput("MYCN", "MYCN amplified", TRUE),
      selectInput("FCG", "FcγR genotype (FCGR3A/2A)",
                  c("low affinity F/F (1.00)" = 1.00, "heterozygous V/F (1.35)" = 1.35,
                    "high affinity V/V (1.80)" = 1.80), selected = 1.35),
      sliderInput("GD2", "GD2 density (10⁶ /cell)", 1, 20, 8, 0.5),
      sliderInput("FMES", "MES (antigen-low) fraction", 0, 0.6, 0.05, 0.01),

      hr(), h4("Regimen"),
      checkboxInput("immuno", "Anti-GD2 antibody (dinutuximab)", TRUE),
      selectInput("fc", "Fc variant",
                  c("ch14.18 (C1q normal)" = 1.00,
                    "hu14.18K322A (C1q ×0.1)" = 0.10), selected = 1.00),
      sliderInput("abdose", "Antibody dose (mg/m²/day)", 0.5, 70, 17.5, 0.5),
      sliderInput("abhrs", "Infusion duration (h)", 0.5, 24, 10, 0.5),
      checkboxInput("use_gm", "GM-CSF (cycles 1,3,5)", TRUE),
      checkboxInput("use_il2", "IL-2 (cycles 2,4 — COG)", TRUE),
      checkboxInput("use_ra", "Isotretinoin", TRUE),
      selectInput("frel", "Isotretinoin administration",
                  c("capsule swallowed whole (F 1.00)" = 1.00,
                    "capsule opened, mixed with food (F 0.60)" = 0.60), selected = 1.00),
      checkboxInput("ra_conc", "Isotretinoin given concurrently with induction", FALSE),
      selectInput("hdct", "High-dose consolidation", c("CEM", "BuMel")),
      checkboxInput("tandem", "Tandem ASCT (ANBL0532)", FALSE),
      sliderInput("ests", "Sodium thiosulfate inner-ear protection", 0, 0.75, 0, 0.05),
      sliderInput("ftumsts", "  …assumed tumour protection (the unwanted effect)", 0, 0.30, 0, 0.05),

      hr(), h4("Delivery"),
      sliderInput("psgtu", "Solid-tumour permeability PSG_TU ×", 0.1, 25, 1, 0.1),
      sliderInput("vifp", "IFP half-effect volume VIFP (L)", 0.005, 0.30, 0.033, 0.005),

      hr(),
      sliderInput("tend", "Simulation horizon (days)", 400, 1500, 1200, 50),
      actionButton("go", "Run", class = "btn-primary btn-block")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs", type = "tabs",

        ## ------------------------------------------------- 1. patient & regimen
        tabPanel("1 · Patient & regimen",
          h4("Dosing schedule as built"),
          DTOutput("schedule"),
          h4("Derived quantities at these settings"),
          tableOutput("derived"),
          tags$small(HTML(
            "The GD2 antigen load is computed as <b>moles per gram</b>: 8×10⁶ copies/cell →
             13.3 nmol/g. One course of dinutuximab is 373 nmol, so the antibody is not
             the scarce resource — what decides the outcome is <b>vascular permeability</b>."))),

        ## ------------------------------------------------------------- 2. PK
        tabPanel("2 · Pharmacokinetics (PK)",
          fluidRow(
            column(6, plotOutput("pk_ab", height = 260)),
            column(6, plotOutput("pk_ra", height = 260))),
          fluidRow(
            column(6, plotOutput("pk_ct", height = 260)),
            column(6, plotOutput("pk_cyt", height = 260))),
          tags$small("The antibody is shown as concentrations in four places: plasma and the
                      interstitium of solid tumour, marrow and dorsal root ganglion. The three tissues differ only in permeability.")),

        ## --------------------------------------------------- 3. the two Fc arms
        tabPanel("3 · The two Fc arms",
          fluidRow(
            column(6, plotOutput("bpc", height = 280)),
            column(6, plotOutput("arms", height = 280))),
          h4("Dose-response geometry"),
          plotOutput("ti_curve", height = 300),
          tags$small(HTML(
            "ADCC is <b>first order</b> in bound IgG per cell (one Fc binds one
             FcγR), while CDC is <b>second order</b> (C1q must bridge <b>two</b>
             adjacent Fc). Because the two arms sit in different compartments and
             those compartments differ in permeability, the optimal dose depends on
             <b>where the lesion is</b>."))),

        ## --------------------------------------------- 4. delivery vs dose
        tabPanel("4 · Delivery barriers",
          plotOutput("perm_vs_dose", height = 330),
          h4("Fate of the injected antibody"),
          tableOutput("fate"),
          tags$small("Dose and permeability enter the model only as a product, so for the
                      solid-tumour arm they are completely interchangeable. Raising permeability,
                      however, does not touch the pain compartment.")),

        ## ------------------------------------------------------ 5. tumour PD
        tabPanel("5 · Tumour burden · PD",
          plotOutput("burden", height = 300),
          fluidRow(
            column(6, plotOutput("pools", height = 260)),
            column(6, plotOutput("killrates", height = 260))),
          tags$small("The four tumour pools have different drug sensitivities:
                      TP (cytotoxic + ADCC + radiation), TQ (G0 — resistant to cytotoxics),
                      TD (differentiated — unresponsive to cytotoxics), TM (marrow — where
                      ADCC actually works).")),

        ## ------------------------------------------------------- 6. endpoints
        tabPanel("6 · Clinical endpoints",
          fluidRow(
            column(4, wellPanel(h4("Time to relapse"), textOutput("t_relapse"))),
            column(4, wellPanel(h4("Nadir burden"), textOutput("t_nadir"))),
            column(4, wellPanel(h4("INRC response"), textOutput("t_resp")))),
          fluidRow(
            column(6, plotOutput("pain", height = 260)),
            column(6, plotOutput("heme", height = 260))),
          fluidRow(
            column(6, plotOutput("oto", height = 240)),
            column(6, tableOutput("tox_table")))),

        ## ------------------------------------------------------- 7. scenarios
        tabPanel("7 · Scenario comparison",
          checkboxGroupInput("scen", "Axes to run",
            choices = c("Antibody dose (both Fc arms)" = "dose",
                        "Permeability vs dose" = "perm",
                        "K322A Fc variant" = "k322a",
                        "Tumour burden dependence" = "burden",
                        "IL-2 / GM-CSF" = "cyto",
                        "Isotretinoin schedule" = "retinoid",
                        "Consolidation intensity" = "consol",
                        "STS ototoxicity trade-off" = "sts"),
            selected = c("dose", "k322a"), inline = TRUE),
          actionButton("runscen", "Run scenarios", class = "btn-warning"),
          hr(), uiOutput("scen_out")),

        ## ------------------------------------------------------------ 8. MIBG
        tabPanel("8 · MIBG dosimetry",
          fluidRow(
            column(4, sliderInput("mibg_mci", "¹³¹I-MIBG (mCi/kg)", 0, 24, 18, 1)),
            column(4, selectInput("mibg_sa", "Specific activity",
                     c("no-carrier-added NCA (1.1e5 MBq/µmol)" = 1.10e5,
                       "carrier-added ×10 (1.1e4)" = 1.10e4,
                       "carrier-added ×100 (1.1e3)" = 1.10e3))),
            column(4, selectInput("netx", "NET-blocking drug (failed washout)",
                     c("none (1.00)" = 1.00, "tricyclic antidepressant (0.55)" = 0.55,
                       "labetalol (0.35)" = 0.35)))),
          fluidRow(
            column(4, sliderInput("kexcr", "Whole-body biological clearance (1/d)",
                                  0.08, 0.40, 0.19, 0.01)),
            column(4, checkboxInput("ki", "Potassium iodide (KI) thyroid blockade", TRUE)),
            column(4, br(), actionButton("runmibg", "Run MIBG",
                                         class = "btn-info"))),
          hr(),
          fluidRow(
            column(6, plotOutput("mibg_dose", height = 280)),
            column(6, plotOutput("mibg_act", height = 280))),
          tableOutput("mibg_tab"),
          tags$small(HTML(
            "The whole-body dose coefficient is a <b>calibrated</b> value (0.20 mGy/MBq), so it is
             not a prediction. What is predicted is how the <b>administered activity</b> needed to
             meet a 2 Gy whole-body dose target changes when biological clearance changes."))),

        ## ---------------------------------------------------------- 9. cohort
        tabPanel("9 · Virtual cohort (EFS)",
          fluidRow(
            column(3, numericInput("ncoh", "Number of patients", 20, 5, 200, 5)),
            column(3, numericInput("seed", "seed", 20260730, 1, 99999999, 1)),
            column(3, checkboxGroupInput("coh_arms", "Arms to compare",
                     c("Isotretinoin alone" = "ctl", "COG full" = "cog",
                       "IL-2 omitted (SIOPEN)" = "noil2"),
                     selected = c("ctl", "cog"))),
            column(3, br(), actionButton("runcoh", "Run cohort",
                                         class = "btn-danger"))),
          hr(),
          plotOutput("efs_curve", height = 340),
          tableOutput("efs_tab"),
          tags$small("Calibration target: ANBL0032 2-year EFS 66% (immunotherapy) vs 46% (control).
                      In a small cohort the confidence interval is wide, so it is shown in the table too."))
      )
    )
  )
)

# =============================================================================
server <- function(input, output, session) {

  ## ------------------------------------------------ parameters from the sidebar
  pset <- reactive({
    list(BSA = input$BSA, WT = input$WT,
         MYCN = as.numeric(input$MYCN),
         FCG = as.numeric(input$FCG),
         GD2DENS = input$GD2 * 1e6,
         FMES = input$FMES,
         C1QEFF = as.numeric(input$fc),
         PSG_TU = 0.020 * input$psgtu,
         VIFP = input$vifp,
         ESTS = input$ests, FTUMSTS = input$ftumsts,
         IMMUNO = as.numeric(input$immuno),
         USE_RA = as.numeric(input$use_ra))
  })

  regimen <- reactive({
    nb_regimen(immuno = input$immuno, use_il2 = input$use_il2,
               use_gm = input$use_gm, use_ra = input$use_ra,
               hdct = input$hdct, tandem = input$tandem,
               ra_concurrent = input$ra_conc, ab_mgm2 = input$abdose,
               ab_hours = input$abhrs, BSA = input$BSA,
               FREL = as.numeric(input$frel))
  })

  sim <- eventReactive(input$go, ignoreNULL = FALSE, {
    withProgress(message = "Integrating…", {
      nb_run(e = regimen(), tp0 = input$tp0, tm0 = input$tm0,
             end = input$tend, delta = 0.25, param = pset())
    })
  })

  ## ------------------------------------------------------------ 1. schedule
  output$schedule <- renderDT({
    regimen() %>% as_tibble() %>%
      mutate(cmt = as.character(cmt)) %>%
      arrange(time) %>% head(200) %>%
      datatable(options = list(pageLength = 8, dom = "tp"), rownames = FALSE)
  })

  output$derived <- renderTable({
    p <- pset()
    ab_nmol_day <- input$abdose * input$BSA / 1000 / MW_AB * 1e9
    ag_per_g <- 1e9 * p$GD2DENS / 6.022e14
    tibble(
      quantity = c("dinutuximab, nmol/day", "dinutuximab, nmol/course (4 days)",
                   "GD2 antigen, nmol/g tumour",
                   "Total solid-tumour antigen (nmol)",
                   "PS product — solid tumour : marrow : nerve (L/d)",
                   "ADCC50 / CDC50 (molecules/cell)"),
      value = c(sprintf("%.1f", ab_nmol_day),
                sprintf("%.1f", 4 * ab_nmol_day),
                sprintf("%.2f", ag_per_g),
                sprintf("%.0f", ag_per_g * 1.25 * input$tp0),
                sprintf("%.4f : %.3f : %.3f",
                        p$PSG_TU * 1.25 * input$tp0 / 1000, 0.300 * 0.300,
                        4.000 * 0.0020),
                "2.0e4 / 5.0e4"))
  })

  ## ------------------------------------------------------------------- 2. PK
  output$pk_ab <- renderPlot({
    o <- sim()
    o %>% transmute(time,
                    `plasma` = CPn,
                    `solid tumour` = ABT / (0.35 * pmax(VTU, 1e-9)),
                    `bone marrow` = ABM / (0.35 * 0.300),
                    `DRG / nerve` = ABN / (0.35 * 0.0020)) %>%
      pivot_longer(-time) %>%
      filter(time >= 195, time <= 380) %>%
      ggplot(aes(time, pmax(value, 1e-4), colour = name)) +
      geom_line(linewidth = 0.6) + scale_y_log10() +
      scale_colour_manual(values = CMP, name = NULL) +
      labs(title = "Anti-GD2 antibody: free concentration in four compartments",
           subtitle = "immunotherapy window only · the only difference is permeability",
           x = "Day", y = "Free antibody (nM, log)") + THEME
  })

  output$pk_ra <- renderPlot({
    sim() %>% filter(time >= 195) %>%
      ggplot(aes(time, CRAo)) + geom_line(colour = "#B9770E", linewidth = 0.6) +
      geom_hline(yintercept = 2, linetype = 2, colour = "grey40") +
      annotate("text", x = Inf, y = 2, label = "2 µM target", hjust = 1.05,
               vjust = -0.5, size = 3, colour = "grey30") +
      labs(title = "Isotretinoin (13-cis-RA)",
           subtitle = "CYP26A1 autoinduction makes exposure fall within a course",
           x = "Day", y = "Plasma concentration (µM)") + THEME
  })

  output$pk_ct <- renderPlot({
    sim() %>% ggplot(aes(time, CCTo)) +
      geom_line(colour = "#5D6D7E", linewidth = 0.5) +
      labs(title = "Composite cytotoxic drug",
           subtitle = "six induction cycles → high-dose consolidation → ASCT", x = "Day",
           y = "Concentration (mg/L)") + THEME
  })

  output$pk_cyt <- renderPlot({
    sim() %>% filter(time >= 195) %>%
      transmute(time, `NK (blood)` = NKB, Treg = TREG,
                `ANC ×50` = ANC * 50) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.6) +
      labs(title = "Effector supply",
           subtitle = "IL-2 expands Treg at a lower EC50 (high-affinity CD25)",
           x = "Day", y = "cells/µL", colour = NULL) + THEME
  })

  ## --------------------------------------------------------- 3. the Fc arms
  output$bpc <- renderPlot({
    sim() %>% filter(time >= 195) %>%
      transmute(time, `solid tumour` = BPCTo, `bone marrow` = BPCMo,
                `DRG / nerve` = BPCNo) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, pmax(value, 1), colour = name)) +
      geom_line(linewidth = 0.7) + scale_y_log10() +
      geom_hline(yintercept = 2e4, linetype = 2, colour = "#1E8449") +
      geom_hline(yintercept = 5e4, linetype = 3, colour = "#CA6F1E") +
      scale_colour_manual(values = CMP, name = NULL) +
      labs(title = "Bound IgG per cell — the actual PD driver",
           subtitle = "dashed = ADCC50 (2×10⁴) · dotted = CDC50 (5×10⁴)",
           x = "Day", y = "Molecules/cell (log)") + THEME
  })

  output$arms <- renderPlot({
    o <- sim() %>% filter(time >= 195)
    a <- 0.90 * (o$BPCTo / (o$BPCTo + 2e4))
    m <- 0.90 * 1.25 * (o$BPCMo / (o$BPCMo + 2e4))
    r <- (o$BPCNo / 5e4)^2
    tibble(time = o$time, `ADCC · solid tumour` = a, `ADCC · bone marrow` = m,
           `CDC · nociceptor` = r / (1 + r) * o$CPL) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.7) +
      labs(title = "Time course of the two Fc arms", x = "Day",
           subtitle = "CDC decays within a course through consumption of complement (CPL)",
           y = "Relative drive", colour = NULL) + THEME
  })

  output$ti_curve <- renderPlot({
    ## algebraic, not simulated — instantaneous geometry at plasma steady state
    d <- 10^seq(log10(0.05), log10(300), length.out = 120)
    scale <- input$abdose
    bT <- 1.42e4 * d / scale; bM <- 3.01e6 * d / scale
    bMs <- 3.73e6 * bM / (bM + 3.73e6)            # marrow saturates on antigen
    bN <- 2.00e5 * d / scale; bNs <- 4.2e5 * bN / (bN + 4.2e5)
    adcc_t <- bT / (bT + 2e4); adcc_m <- bMs / (bMs + 2e4)
    r <- (bNs / 5e4)^2; cdc <- as.numeric(input$fc) * r / (1 + r)
    tibble(dose = d, `ADCC · solid tumour` = adcc_t,
           `ADCC · bone marrow` = adcc_m, `CDC · nociceptor` = cdc) %>%
      pivot_longer(-dose) %>%
      ggplot(aes(dose, value, colour = name)) + geom_line(linewidth = 0.8) +
      geom_vline(xintercept = 17.5, linetype = 2, colour = "grey40") +
      annotate("text", x = 17.5, y = 0.05, label = "Clinical dose", angle = 90,
               vjust = -0.4, size = 3, colour = "grey30") +
      scale_x_log10() +
      labs(title = "Dose-response geometry: the marrow arm is already saturated, the solid-tumour arm still linear",
           subtitle = "so the optimal dose depends on where the lesion is",
           x = "dinutuximab (mg/m²/day, log)", y = "Drive (0-1)",
           colour = NULL) + THEME
  })

  ## ------------------------------------------------ 4. permeability vs dose
  output$perm_vs_dose <- renderPlot({
    f <- c(0.1, 0.25, 0.5, 1, 2, 5, 10, 25)
    ## the model enters dose and permeability as a product, so one curve serves
    tibble(fold = f,
           `Permeability ×N (pain unchanged)` = 1.42e4 * f,
           `Dose ×N (pain increases)` = 1.42e4 * f) %>%
      pivot_longer(-fold) %>%
      mutate(drive = value / (value + 2e4)) %>%
      ggplot(aes(fold, drive, colour = name, linetype = name)) +
      geom_line(linewidth = 0.9) + geom_point(size = 1.8) + scale_x_log10() +
      labs(title = "Solid-tumour ADCC: permeability and dose are completely interchangeable",
           subtitle = "because they enter only as a product. But permeability does not touch the pain compartment",
           x = "Fold", y = "ADCC drive", colour = NULL,
           linetype = NULL) + THEME
  })

  output$fate <- renderTable({
    o <- sim(); dt <- diff(o$time)
    intg <- function(x) sum(dt * head(x, -1))
    dose <- 4 * 5 * input$abdose * input$BSA / 1000 / MW_AB * 1e9
    tibble(compartment = c("Total injected (nmol)", "Bound in solid tumour (nmol)",
                           "Bound in marrow (nmol)", "Bound in dorsal root ganglion (nmol)"),
           amount = c(dose, 0.030 * intg(o$BNDT), 0.030 * intg(o$BNDM),
                      0.030 * intg(o$BNDN)),
           `% of dose` = c(100, 100 * 0.030 * intg(o$BNDT) / dose,
                           100 * 0.030 * intg(o$BNDM) / dose,
                           100 * 0.030 * intg(o$BNDN) / dose))
  }, digits = 6)

  ## ---------------------------------------------------------- 5. tumour PD
  output$burden <- renderPlot({
    o <- sim()
    ggplot(o, aes(time, pmax(BURDEN, 1e-9))) +
      geom_line(linewidth = 0.7, colour = "#922B21") + scale_y_log10() +
      geom_hline(yintercept = 0.5, linetype = 2, colour = "grey40") +
      annotate("rect", xmin = 200, xmax = 368, ymin = 1e-9, ymax = Inf,
               alpha = 0.06, fill = "#1E8449") +
      annotate("text", x = 284, y = 1e-8, label = "Immunotherapy window", size = 3,
               colour = "#145A32") +
      labs(title = "Total tumour burden", subtitle = "dashed = relapse detection threshold (0.5 g)",
           x = "Day", y = "10⁹ cells (log)") + THEME
  })

  output$pools <- renderPlot({
    sim() %>% transmute(time, TP, TQ, TD, TM) %>% pivot_longer(-time) %>%
      ggplot(aes(time, pmax(value, 1e-9), colour = name)) +
      geom_line(linewidth = 0.6) + scale_y_log10() +
      labs(title = "The four tumour pools", x = "Day", y = "10⁹ cells (log)",
           colour = NULL) + THEME
  })

  output$killrates <- renderPlot({
    o <- sim() %>% filter(time >= 195)
    tibble(time = o$time,
           `ADCC (marrow)` = 0.90 * 1.25 * o$BPCMo / (o$BPCMo + 2e4),
           `ADCC (solid)` = 0.90 * o$BPCTo / (o$BPCTo + 2e4),
           `cytotoxic` = 1.95 * o$CCTo / (o$CCTo + 1.10)) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.6) +
      labs(title = "Kill rates compared (1/day)", x = "Day", y = "kill rate (1/d)",
           colour = NULL) + THEME
  })

  ## ---------------------------------------------------------- 6. endpoints
  output$t_relapse <- renderText({
    t <- nb_relapse(sim())
    if (is.infinite(t)) "No relapse (complete eradication in the model)" else sprintf("%.0f days", t)
  })
  output$t_nadir <- renderText(sprintf("%.3g × 10⁹ cells", min(sim()$BURDEN)))
  output$t_resp <- renderText({
    o <- sim(); b0 <- o$BURDEN[1]; bn <- min(o$BURDEN)
    r <- 1 - bn / b0
    if (bn < 1e-6) "CR (below the model's detection limit)" else
      if (r > 0.9) "VGPR / PR" else if (r > 0.3) "MR" else "SD / PD"
  })

  output$pain <- renderPlot({
    sim() %>% filter(time >= 195) %>%
      ggplot(aes(time, PAIN)) +
      geom_line(colour = "#CA6F1E", linewidth = 0.7) +
      labs(title = "Pain (allodynia) intensity — complement mediated",
           subtitle = "the within-course decay comes from consumption of the complement pool",
           x = "Day", y = "Relative intensity") + THEME
  })

  output$heme <- renderPlot({
    sim() %>% transmute(time, `ANC (10⁹/L)` = ANC,
                        `PLT /100 (10⁹/L)` = PLT / 100) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.6) +
      geom_hline(yintercept = 0.5, linetype = 2, colour = "grey50") +
      labs(title = "Myelosuppression (Friberg structure)", x = "Day", y = NULL,
           colour = NULL) + THEME
  })

  output$oto <- renderPlot({
    sim() %>% ggplot(aes(time, OTO)) +
      geom_line(colour = "#943126", linewidth = 0.6) +
      labs(title = "Cumulative ototoxicity (cumulative platinum exposure)", x = "Day",
           y = "Relative units") + THEME
  })

  output$tox_table <- renderTable({
    o <- sim()
    tibble(endpoint = c("ANC nadir", "PLT nadir", "Peak pain", "Complement nadir (CPL)",
                        "Cumulative ototoxicity", "Cumulative platinum AUC"),
           value = c(min(o$ANC), min(o$PLT), max(o$PAIN), min(o$CPL),
                     max(o$OTO), max(o$PLATAUC)))
  }, digits = 3)

  ## ---------------------------------------------------------- 7. scenarios
  scen_res <- eventReactive(input$runscen, {
    out <- list()
    withProgress(message = "Running scenarios…", {
      if ("dose" %in% input$scen)     out$dose <- nb_scen_dose()
      if ("perm" %in% input$scen)     out$perm <- nb_scen_perm()
      if ("k322a" %in% input$scen)    out$k322a <- nb_scen_k322a()
      if ("burden" %in% input$scen)   out$burden <- nb_scen_burden()
      if ("cyto" %in% input$scen)     out$cyto <- nb_scen_cytokine()
      if ("retinoid" %in% input$scen) out$retinoid <- nb_scen_retinoid()
      if ("sts" %in% input$scen)      out$sts <- nb_scen_sts()
    })
    out
  })

  output$scen_out <- renderUI({
    r <- scen_res()
    if (!length(r)) return(helpText("Choose an axis and press 'Run scenarios'."))
    do.call(tagList, lapply(names(r), function(nm) {
      tagList(h4(nm), renderTable(r[[nm]], digits = 4)())
    }))
  })

  ## --------------------------------------------------------------- 8. MIBG
  mibg <- eventReactive(input$runmibg, {
    MBq <- input$mibg_mci * 37 * input$WT
    withProgress(message = "Integrating MIBG…", {
      nb_run(e = ev_mibg(0, MBq, as.numeric(input$mibg_sa)),
             tp0 = input$tp0, tm0 = input$tm0, end = 60, delta = 0.05,
             param = c(pset(), list(NETX = as.numeric(input$netx),
                                    KEXCR = input$kexcr,
                                    KIBLOCK = if (input$ki) 0.10 else 1.00)))
    })
  })

  output$mibg_dose <- renderPlot({
    mibg() %>% transmute(time, `Whole body` = DWB,
                         `Tumour` = DTU) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.7) +
      geom_hline(yintercept = 2, linetype = 2, colour = "grey40") +
      labs(title = "Cumulative absorbed dose", subtitle = "dashed = 2 Gy whole-body dose (the stem-cell-rescue threshold)",
           x = "Day", y = "Gy", colour = NULL) + THEME
  })

  output$mibg_act <- renderPlot({
    mibg() %>% transmute(time, `Blood` = MBC, `Tumour` = MBT,
                         `Thyroid` = THY) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, pmax(value, 1e-2), colour = name)) +
      geom_line(linewidth = 0.7) + scale_y_log10() +
      labs(title = "Residual activity", x = "Day", y = "MBq (log)", colour = NULL) +
      THEME
  })

  output$mibg_tab <- renderTable({
    o <- mibg(); MBq <- input$mibg_mci * 37 * input$WT
    tibble(quantity = c("Administered activity (MBq)", "Whole-body dose (Gy)",
                        "mGy / MBq", "Tumour dose (Gy)", "ANC nadir",
                        "Peak thyroid activity (MBq)"),
           value = c(MBq, max(o$DWB), 1000 * max(o$DWB) / MBq, max(o$DTU),
                     min(o$ANC), max(o$THY)))
  }, digits = 3)

  ## -------------------------------------------------------------- 9. cohort
  coh <- eventReactive(input$runcoh, {
    co <- nb_cohort(input$ncoh, input$seed)
    res <- list()
    withProgress(message = "Running the cohort (this takes some minutes)…", {
      if ("ctl" %in% input$coh_arms)
        res[["Isotretinoin alone"]] <- nb_efs(co, immuno = FALSE)
      if ("cog" %in% input$coh_arms)
        res[["COG full (Ab+GM+IL2)"]] <- nb_efs(co, immuno = TRUE)
      if ("noil2" %in% input$coh_arms)
        res[["IL-2 omitted (SIOPEN)"]] <- nb_efs(co, immuno = TRUE, use_il2 = FALSE)
    })
    res
  })

  output$efs_curve <- renderPlot({
    r <- coh()
    grid <- seq(0, 1100, 10)
    purrr::imap_dfr(r, function(t, nm)
      tibble(arm = nm, time = grid,
             efs = sapply(grid, function(g) mean(t > g)))) %>%
      ggplot(aes(time, efs, colour = arm)) + geom_step(linewidth = 0.9) +
      geom_vline(xintercept = 730, linetype = 2, colour = "grey40") +
      scale_y_continuous(limits = c(0, 1)) +
      labs(title = "Event-free survival",
           subtitle = "dashed = 2 years (the ANBL0032 endpoint)", x = "Day",
           y = "EFS", colour = NULL) + THEME
  })

  output$efs_tab <- renderTable({
    r <- coh()
    purrr::imap_dfr(r, function(t, nm) {
      p2 <- mean(t > 730); n <- length(t)
      se <- sqrt(p2 * (1 - p2) / n)
      tibble(arm = nm, n = n,
             `1-year EFS` = mean(t > 365), `2-year EFS` = p2,
             `2-year 95% CI` = sprintf("%.2f – %.2f", max(0, p2 - 1.96 * se),
                                    min(1, p2 + 1.96 * se)),
             `3-year EFS` = mean(t > 1095))
    })
  }, digits = 3)
}

shinyApp(ui, server)
