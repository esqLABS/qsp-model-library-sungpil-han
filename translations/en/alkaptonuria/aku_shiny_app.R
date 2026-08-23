## ============================================================================
##  Alkaptonuria (AKU) QSP — Shiny dashboard
##  12 tabs.  The app is organised around the three balances that the model is
##  built on, because those are what a user needs to be able to interrogate:
##
##    Balance 1  the flux is conserved; nitisinone changes its EXIT, not its
##               size, so DOSE sets HGA and DIET sets tyrosine
##    Balance 2  the toxic branch is ~1.5e-5 of the flux and is one-way, so
##               benefit is bounded by the integral not yet accrued
##    Balance 3  HGA clearance sits at the renal-plasma-flow ceiling, so the
##               trial endpoint (urinary HGA) and the causal quantity (serum
##               HGA) diverge on treatment
##
##  Run with:  shiny::runApp("aku_shiny_app.R")
##  Requires:  shiny, mrgsolve, dplyr, tidyr, ggplot2, DT
##  The model file must sit alongside this one.
## ============================================================================

suppressPackageStartupMessages({
  library(shiny); library(mrgsolve); library(dplyr); library(tidyr)
  library(ggplot2); library(DT)
})

source("aku_mrgsolve_model.R", local = TRUE)

AKU_MOD <- aku_model()

THEME <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(colour = "grey35", size = 10),
        legend.position = "bottom")

PAL <- c("#e63946", "#457b9d", "#2a9d8f", "#e9c46a", "#7b2cbf",
         "#fb8500", "#8d99ae", "#0288d1")

## ---------------------------------------------------------------------------
## UI
## ---------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel(
    div(
      h3("Alkaptonuria QSP Dashboard"),
      p(style = "color:#555;font-size:13px;margin-top:-6px;",
        HTML("HGD deficiency · ochronosis as a one-way integral · nitisinone as a ",
             "<b>flux-redirection</b> agent. ",
             "The disease is the <b>integral</b> of an irreversible branch ",
             "amounting to 1.5e-5 of ",
             "the conserved flux."))
    )
  ),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Patient"),
      sliderInput("bmi", "BMI", 18, 40, 25, 0.5),
      sliderInput("occup", "Occupational load", 0.8, 1.8, 1.0, 0.1),
      radioButtons("sexm", "Sex", c("Male" = 1, "Female" = 0),
                   selected = 1, inline = TRUE),
      sliderInput("resact", "Residual HGD activity %", 0, 5, 0, 0.5),
      sliderInput("ckdx", "Renal decline factor (x)", 1, 4, 1, 0.5),

      hr(), h4("Treatment"),
      sliderInput("dose", "Nitisinone (mg/day)", 0, 20, 10, 1),
      sliderInput("start", "Age at initiation (y)", 2, 65, 25, 1),
      sliderInput("stop", "Age at stop (y; 70 = continue)", 5, 70, 70, 1),
      sliderInput("adh", "Adherence", 0.4, 1, 1, 0.05),

      hr(), h4("Diet — the only control point for tyrosine"),
      sliderInput("prot", "Protein intake (g/day)", 35, 120, 70, 1),
      sliderInput("aasupp", "Phe/Tyr-free supplement (g/day)", 0, 40, 0, 1),

      hr(), h4("Adjuncts"),
      checkboxInput("nsaid", "NSAID", FALSE),
      checkboxInput("physio", "Physiotherapy · weight management", FALSE),
      checkboxInput("alkali", "Urine alkalinisation · hydration", FALSE),
      sliderInput("casc", "Plasma ascorbate (umol/L)", 30, 150, 50, 5),
      checkboxInput("ideal", "IDEAL comparator: HGA blocked with no tyrosine rise", FALSE),

      hr(),
      sliderInput("endage", "Simulation end age", 40, 85, 70, 1),
      actionButton("go", "Run", class = "btn-primary btn-block"),
      helpText(style = "font-size:11px;",
               "The 70-year integration takes a few seconds. The matched untreated control is always computed alongside it.")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs", type = "tabs",

        tabPanel("① Patient profile",
          br(), fluidRow(column(12, uiOutput("kpis"))),
          br(), plotOutput("p_overview", height = "420px"),
          br(), h5("Summary of the current settings"), DTOutput("t_profile")),

        tabPanel("② Nitisinone PK",
          br(), plotOutput("p_pk", height = "300px"),
          plotOutput("p_pkss", height = "260px"),
          helpText(HTML("Because the half-life is about 54 hours, <b>once-daily and twice-daily ",
                        "dosing are almost indistinguishable</b>. The half-life is far longer than the dosing interval.")),
          br(), DTOutput("t_pk")),

        tabPanel("③ Pathway flux (Balance 1)",
          br(), h5("What comes in must go out — what changes is the composition of the exits"),
          plotOutput("p_flux", height = "330px"),
          plotOutput("p_exit", height = "300px"),
          br(), DTOutput("t_mass")),

        tabPanel("④ HGA: urine vs serum (Balance 3)",
          br(), h5("The divergence between the trial endpoint and the causal quantity"),
          plotOutput("p_hga", height = "330px"),
          plotOutput("p_gap", height = "300px"),
          helpText(HTML("What pigments cartilage is not urinary HGA but <b>plasma HGA</b>. ",
                        "On treatment HPPA and HPLA rise more than 14-fold and competitively occupy the ",
                        "same organic anion secretion route, so HGA clearance actually falls, and as a result ",
                        "the plasma concentration comes down far less than urinary excretion does."))),

        tabPanel("⑤ Tyrosinaemia (dose vs diet)",
          br(), h5("2x2: dose sets HGA, diet sets tyrosine"),
          plotOutput("p_tyr", height = "320px"),
          plotOutput("p_2x2", height = "330px"),
          br(), DTOutput("t_tyr")),

        tabPanel("⑥ Pigment integral (Balance 2)",
          br(), plotOutput("p_pig", height = "330px"),
          plotOutput("p_pigdist", height = "300px"),
          helpText(HTML("Pigment in articular cartilage, the discs and the aortic valve has <b>no loss term</b>, ",
                        "because that collagen is not replaced. A loss term matching the collagen turnover rate ",
                        "was given only to skin, sclera and ear cartilage, and the model predicts that the ",
                        "reported reversal of pigmentation is confined to those tissues."))),

        tabPanel("⑦ Joints · spine",
          br(), plotOutput("p_joint", height = "330px"),
          plotOutput("p_spine", height = "300px"),
          br(), DTOutput("t_struct")),

        tabPanel("⑧ Heart · kidney · other organs",
          br(), plotOutput("p_valve", height = "300px"),
          plotOutput("p_organ", height = "320px")),

        tabPanel("⑨ Clinical endpoints · AKUSSI",
          br(), plotOutput("p_akussi", height = "340px"),
          plotOutput("p_hazard", height = "300px"),
          br(), DTOutput("t_end")),

        tabPanel("⑩ Safety: keratopathy",
          br(), plotOutput("p_kerato", height = "330px"),
          helpText(HTML("Management thresholds: tighten protein restriction above a tyrosine of 700 umol/L, ",
                        "add a Phe/Tyr-free amino acid supplement above 900 umol/L.")),
          plotOutput("p_safety", height = "290px")),

        tabPanel("⑪ Age at initiation and headroom",
          br(), h5("Same drug, same dose — only the remaining integral differs"),
          actionButton("go_hr", "Run the headroom scan (13 ages at initiation)",
                       class = "btn-warning"),
          br(), br(), plotOutput("p_headroom", height = "340px"),
          plotOutput("p_headroom2", height = "300px"),
          br(), DTOutput("t_headroom")),

        tabPanel("⑫ Scenario comparison · validation",
          br(), h5("The 24 predefined scenarios and the held-out validation"),
          actionButton("go_scn", "Run the scenario set", class = "btn-warning"),
          actionButton("go_val", "Run the held-out validation", class = "btn-info"),
          br(), br(), plotOutput("p_scn", height = "340px"),
          br(), h5("Held-out validation (values not used in the calibration)"),
          DTOutput("t_val"),
          br(), h5("Scenario summary at end of follow-up"), DTOutput("t_scn"))
      )
    )
  )
)

## ---------------------------------------------------------------------------
## Server
## ---------------------------------------------------------------------------
server <- function(input, output, session) {

  pars <- reactive({
    list(BMI = input$bmi, OCCUP = input$occup, SEXM = as.numeric(input$sexm),
         RESACT = input$resact/100, CKDX = input$ckdx,
         PROT = input$prot, AASUPP = input$aasupp,
         FBIO = input$adh, CASC = input$casc,
         NSAIDON = as.numeric(input$nsaid), PHYSIO = as.numeric(input$physio),
         ALKALI = as.numeric(input$alkali), IDEALDRG = as.numeric(input$ideal))
  })

  ## treated arm and its matched untreated control -- always both, so that no
  ## conclusion on this dashboard can come from an unmatched comparison
  sims <- eventReactive(input$go, {
    withProgress(message = "Integrating 70 years...", value = 0.2, {
      tx <- sim_aku(AKU_MOD, dose_mg = input$dose, start_age = input$start,
                    stop_age = if (input$stop >= 70) Inf else input$stop,
                    end_age = input$endage, pars = pars(), label = "treated")
      incProgress(0.5)
      ct <- sim_aku(AKU_MOD, dose_mg = 0, start_age = input$start,
                    end_age = input$endage, pars = pars(), label = "control")
      list(tx = tx, ct = ct, both = bind_rows(tx, ct))
    })
  }, ignoreNULL = FALSE)

  fin <- reactive(tail(sims()$tx, 1))
  finc <- reactive(tail(sims()$ct, 1))

  ## ---------------------------------------------------------------- ① profile
  output$kpis <- renderUI({
    f <- fin(); c0 <- finc()
    box <- function(lab, val, sub = "", col = "#457b9d")
      div(style = paste0("display:inline-block;width:15.5%;margin:3px;padding:9px;",
                         "border-left:4px solid ", col, ";background:#fafafa;"),
          div(style = "font-size:10.5px;color:#666;", lab),
          div(style = "font-size:19px;font-weight:bold;", val),
          div(style = "font-size:10px;color:#888;", sub))
    tagList(
      box("Serum HGA", sprintf("%.2f", f$CHGAo), "umol/L (vs untreated)", "#e63946"),
      box("Urinary HGA", sprintf("%.0f", f$UHGA24), "umol/day", "#c1121f"),
      box("Plasma tyrosine", sprintf("%.0f", f$CTYRo), "umol/L", "#e9c46a"),
      box("cAKUSSI", sprintf("%.1f", f$CAKUSSI), sprintf("untreated %.1f", c0$CAKUSSI), "#212121"),
      box("Pigment spared", sprintf("%.0f%%", 100*f$PIGSPARE), "vs counterfactual", "#7b2cbf"),
      box("Joint replacement probability", sprintf("%.0f%%", 100*f$PJR), sprintf("untreated %.0f%%", 100*c0$PJR), "#2a9d8f")
    )
  })

  output$p_overview <- renderPlot({
    d <- sims()$both %>%
      select(AGEYo, scenario, CHGAo, CTYRo, CART, DISCH, CAKUSSI, PAINVAS) %>%
      pivot_longer(-c(AGEYo, scenario))
    lab <- c(CHGAo = "Serum HGA (umol/L)", CTYRo = "Plasma Tyr (umol/L)",
             CART = "Intact cartilage (fraction)", DISCH = "Disc height (fraction)",
             CAKUSSI = "cAKUSSI", PAINVAS = "Pain VAS")
    d$name <- factor(lab[d$name], levels = lab)
    ggplot(d, aes(AGEYo, value, colour = scenario)) +
      geom_line(linewidth = 0.85) +
      geom_vline(xintercept = input$start, linetype = 3, colour = "grey40") +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      scale_colour_manual(values = c(control = "#8d99ae", treated = "#e63946")) +
      labs(title = "The treated arm and its matched untreated control",
           subtitle = paste0("Dotted = start of treatment (age ", input$start, "). ",
                             "Every other setting is identical in the untreated arm."),
           x = "Age (years)", y = NULL, colour = NULL) + THEME
  })

  output$t_profile <- renderDT({
    p <- pars()
    datatable(data.frame(parameter = names(p), value = unlist(p)),
              options = list(dom = "t", pageLength = 20), rownames = FALSE)
  })

  ## --------------------------------------------------------------------- ② PK
  output$p_pk <- renderPlot({
    st <- input$start*365.25
    mo <- aku_init(AKU_MOD) %>% param(pars())
    e  <- ev(time = st, amt = input$dose*NT_MG_TO_UMOL, ii = 1, addl = 60)
    o  <- mo %>% ev(e) %>%
      mrgsim(start = st, end = st + 62, delta = 0.05,
             atol = 1e-12, rtol = 1e-10, maxsteps = 2e6) %>% as_tibble()
    ggplot(o, aes((TIME-st), CNTo)) + geom_line(colour = PAL[2], linewidth = 0.8) +
      labs(title = "Nitisinone plasma concentration — 60 days of accumulation",
           subtitle = "With a t1/2 of about 54 hours it takes 1-2 weeks to reach steady state",
           x = "Days after start of treatment", y = "Nitisinone (umol/L)") + THEME
  })

  output$p_pkss <- renderPlot({
    d <- bind_rows(lapply(c(1, 2, 4, 8, 10, 20), function(dd) {
      b <- probe_biochem(AKU_MOD, dd, 0.5, input$start, pars())
      tibble(dose = dd, CNT = b$CNTo, uHGA = b$UHGA24, sHGA = b$CHGAo, Tyr = b$CTYRo)
    }))
    ggplot(d, aes(dose, CNT)) +
      geom_line(colour = PAL[2]) + geom_point(size = 2.6, colour = PAL[2]) +
      labs(title = "Dose–steady-state concentration is linear",
           subtitle = "The PK is linear, yet the HGA suppression that follows it is not (next tab)",
           x = "Nitisinone (mg/day)", y = "Steady-state nitisinone (umol/L)") + THEME
  })

  output$t_pk <- renderDT({
    d <- bind_rows(lapply(c(0, 1, 2, 4, 8, 10, 20), function(dd) {
      b <- probe_biochem(AKU_MOD, dd, 1, input$start, pars())
      tibble(`Dose mg/day` = dd, `Nitisinone umol/L` = round(b$CNTo, 2),
             `Urinary HGA umol/day` = round(b$UHGA24),
             `Serum HGA umol/L` = round(b$CHGAo, 2),
             `Plasma Tyr umol/L` = round(b$CTYRo))
    }))
    datatable(d, options = list(dom = "t"), rownames = FALSE)
  })

  ## ------------------------------------------------------------------- ③ flux
  output$p_flux <- renderPlot({
    d <- sims()$both %>% select(AGEYo, scenario, FLUXIN, EXITSUM) %>%
      pivot_longer(-c(AGEYo, scenario))
    ggplot(d, aes(AGEYo, value/1000, colour = name, linetype = scenario)) +
      geom_line(linewidth = 0.85) +
      scale_colour_manual(values = c(FLUXIN = "#2a9d8f", EXITSUM = "#e63946"),
                          labels = c(EXITSUM = "Measured sum of exits", FLUXIN = "Dietary Phe+Tyr inflow")) +
      labs(title = "Mass conservation check: inflow and the sum of the exits",
           subtitle = "Nitisinone does not change this total. It changes only the composition of the exits.",
           x = "Age (years)", y = "mmol/day", colour = NULL, linetype = NULL) + THEME
  })

  output$p_exit <- renderPlot({
    d <- sims()$both %>%
      select(AGEYo, scenario, UHGA24, UTYR24, UHPP24, UHPL24, UCON24) %>%
      pivot_longer(-c(AGEYo, scenario))
    lab <- c(UHGA24 = "Urinary HGA", UTYR24 = "Urinary tyrosine", UHPP24 = "Urinary HPPA",
             UHPL24 = "Urinary HPLA", UCON24 = "Tyrosine conjugates")
    d$name <- factor(lab[d$name], levels = lab)
    ggplot(d, aes(AGEYo, value/1000, fill = name)) +
      geom_area(position = "stack") + facet_wrap(~scenario) +
      scale_fill_manual(values = PAL[1:5]) +
      labs(title = "The composition of the exits — what nitisinone actually does",
           subtitle = "Close the HGA exit and the same quantity leaves as HPPA · HPLA · tyrosine · conjugates",
           x = "Age (years)", y = "mmol/day", fill = NULL) + THEME
  })

  output$t_mass <- renderDT({
    f <- fin(); c0 <- finc()
    d <- tibble(
      `Item` = c("Dietary Phe+Tyr inflow", "Urinary HGA", "Urinary tyrosine", "Urinary HPPA",
               "Urinary HPLA", "Tyrosine conjugates", "Sum of exits", "Exit/inflow"),
      `Treated` = c(f$FLUXIN, f$UHGA24, f$UTYR24, f$UHPP24, f$UHPL24, f$UCON24,
               f$EXITSUM, f$EXITSUM/f$FLUXIN),
      `Untreated` = c(c0$FLUXIN, c0$UHGA24, c0$UTYR24, c0$UHPP24, c0$UHPL24,
                 c0$UCON24, c0$EXITSUM, c0$EXITSUM/c0$FLUXIN))
    d$`Treated` <- ifelse(d$`Item` == "Exit/inflow", round(d$`Treated`, 3), round(d$`Treated`))
    d$`Untreated` <- ifelse(d$`Item` == "Exit/inflow", round(d$`Untreated`, 3), round(d$`Untreated`))
    datatable(d, options = list(dom = "t"), rownames = FALSE)
  })

  ## ------------------------------------------------------------------- ④ HGA
  output$p_hga <- renderPlot({
    d <- sims()$both %>% select(AGEYo, scenario, CHGAo, UHGA24)
    ggplot(d) +
      geom_line(aes(AGEYo, CHGAo, colour = scenario), linewidth = 0.9) +
      geom_line(aes(AGEYo, UHGA24/1000, colour = scenario), linetype = 2) +
      scale_colour_manual(values = c(control = "#8d99ae", treated = "#e63946")) +
      labs(title = "Serum HGA (solid, umol/L) vs urinary HGA (dashed, mmol/day)",
           subtitle = paste("That the two curves fall by different amounts is the central claim of this model.",
                            "The rate of pigment deposition follows the solid line."),
           x = "Age (years)", y = NULL, colour = NULL) + THEME
  })

  output$p_gap <- renderPlot({
    g <- analysis_endpoint_gap(AKU_MOD, pars())
    d <- g %>% select(dose_mg, uHGA_pct_drop, sHGA_pct_drop) %>%
      pivot_longer(-dose_mg)
    ggplot(d, aes(dose_mg, value, colour = name)) +
      geom_line(linewidth = 0.9) + geom_point(size = 2.4) +
      scale_colour_manual(values = c(uHGA_pct_drop = "#c1121f", sHGA_pct_drop = "#457b9d"),
                          labels = c(sHGA_pct_drop = "Serum HGA reduction %",
                                     uHGA_pct_drop = "Urinary HGA reduction % (the trial endpoint)")) +
      labs(title = "The endpoint looks better than the causal quantity",
           subtitle = "The urinary reduction systematically outruns the serum reduction",
           x = "Nitisinone (mg/day)", y = "Reduction vs untreated (%)", colour = NULL) + THEME
  })

  ## ------------------------------------------------------------------- ⑤ Tyr
  output$p_tyr <- renderPlot({
    d <- sims()$both %>% select(AGEYo, scenario, CTYRo, CPHEo, CHPPo, CHPLAo) %>%
      pivot_longer(-c(AGEYo, scenario))
    lab <- c(CTYRo = "Tyrosine", CPHEo = "Phenylalanine", CHPPo = "HPPA", CHPLAo = "HPLA")
    d$name <- factor(lab[d$name], levels = lab)
    ggplot(d, aes(AGEYo, value, colour = name, linetype = scenario)) +
      geom_line(linewidth = 0.8) +
      geom_hline(yintercept = c(700, 900), linetype = 3, colour = "grey50") +
      scale_colour_manual(values = PAL[c(4,3,2,5)]) +
      labs(title = "Back-pressure in the upstream metabolites",
           subtitle = "Dotted = the tyrosine 700 / 900 umol/L management thresholds",
           x = "Age (years)", y = "umol/L", colour = NULL, linetype = NULL) + THEME
  })

  output$p_2x2 <- renderPlot({
    d <- analysis_dose_vs_diet(AKU_MOD, pars())
    ggplot(d, aes(dose_mg, sTYR, colour = factor(protein_g))) +
      geom_line(linewidth = 0.9) + geom_point(size = 2) +
      geom_hline(yintercept = 700, linetype = 3) +
      scale_colour_manual(values = PAL[1:5], name = "Protein g/day") +
      labs(title = "Tyrosine is set by diet, not by dose",
           subtitle = paste("Move the horizontal axis 20-fold and the curves stay almost flat;",
                            "the spacing between the curves (diet) is everything"),
           x = "Nitisinone (mg/day)", y = "Plasma tyrosine (umol/L)") + THEME
  })

  output$t_tyr <- renderDT({
    d <- analysis_dose_vs_diet(AKU_MOD, pars()) %>%
      mutate(across(c(sHGA, sTYR, uHGA24, pigment_rate), ~round(.x, 3)))
    datatable(d, options = list(pageLength = 12), rownames = FALSE)
  })

  ## ---------------------------------------------------------------- ⑥ pigment
  output$p_pig <- renderPlot({
    d <- sims()$both %>% select(AGEYo, scenario, PIGTOT, CFPIG, PDCo, BRCo) %>%
      pivot_longer(-c(AGEYo, scenario))
    lab <- c(PIGTOT = "Total pigment (umol)", CFPIG = "Untreated counterfactual (cartilage, umol)",
             PDCo = "Cartilage pigment density", BRCo = "Cartilage embrittlement (0-1)")
    d$name <- factor(lab[d$name], levels = lab)
    ggplot(d, aes(AGEYo, value, colour = scenario)) +
      geom_line(linewidth = 0.85) + facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = c(control = "#8d99ae", treated = "#e63946")) +
      labs(title = "Pigment is an integral, and treatment touches only the future part of it",
           subtitle = "Embrittlement is a steep Hill function of pigment density, so it appears abruptly in the third decade",
           x = "Age (years)", y = NULL, colour = NULL) + THEME
  })

  output$p_pigdist <- renderPlot({
    f <- fin()
    d <- tibble(depot = c("Articular cartilage", "Discs", "Aortic valve", "Tendon",
                          "Sclera", "Ear", "Skin", "Prostate · kidney · other"),
                umol = c(f$PCART, f$PDISC, f$PVALV, f$PTEND,
                         f$PSCL, f$PEAR, f$PSKIN, f$POTH),
                turnover = c("None", "None", "None", "None",
                             "Present", "Present", "Present", "None"))
    ggplot(d, aes(reorder(depot, umol), umol, fill = turnover)) +
      geom_col() + coord_flip() +
      scale_fill_manual(values = c("None" = "#212529", "Present" = "#457b9d"),
                        name = "Collagen turnover") +
      labs(title = paste0("Distribution of pigment deposition (at age ", round(f$AGEYo), ")"),
           subtitle = "Pigment can fade only in tissues that turn over — articular cartilage has no loss term",
           x = NULL, y = "umol HGA-equivalents") + THEME
  })

  ## ------------------------------------------------------------------ ⑦ joint
  output$p_joint <- renderPlot({
    d <- sims()$both %>% select(AGEYo, scenario, CART, FRAGC, SYN, SUBCH, OSTEO) %>%
      pivot_longer(-c(AGEYo, scenario))
    lab <- c(CART = "Intact cartilage", FRAGC = "Intra-articular pigment fragments", SYN = "Synovitis",
             SUBCH = "Subchondral sclerosis", OSTEO = "Osteophytes")
    d$name <- factor(lab[d$name], levels = lab)
    ggplot(d, aes(AGEYo, value, colour = scenario)) +
      geom_line(linewidth = 0.85) + facet_wrap(~name, scales = "free_y", ncol = 3) +
      scale_colour_manual(values = c(control = "#8d99ae", treated = "#e63946")) +
      labs(title = "The amplifying loop: fragments → synovitis → MMP → more cartilage loss",
           subtitle = "The acceleration comes out of this loop, not out of a parameter",
           x = "Age (years)", y = NULL, colour = NULL) + THEME
  })

  output$p_spine <- renderPlot({
    d <- sims()$both %>% select(AGEYo, scenario, DISCH, DCALC, ANKY) %>%
      pivot_longer(-c(AGEYo, scenario))
    lab <- c(DISCH = "Disc height", DCALC = "Disc calcification", ANKY = "Spinal ankylosis index")
    d$name <- factor(lab[d$name], levels = lab)
    ggplot(d, aes(AGEYo, value, colour = scenario)) +
      geom_line(linewidth = 0.85) + facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = c(control = "#8d99ae", treated = "#e63946")) +
      labs(title = "The spine — the earliest radiographic trace of the pigment",
           x = "Age (years)", y = NULL, colour = NULL) + THEME
  })

  output$t_struct <- renderDT({
    f <- fin(); c0 <- finc()
    datatable(tibble(
      `Metric` = c("Intact cartilage", "Synovitis", "Disc height", "Disc calcification",
               "Spinal ankylosis", "Tendon integrity", "Bone mineral density"),
      `Treated` = round(c(f$CART, f$SYN, f$DISCH, f$DCALC, f$ANKY, f$TENDI, f$BMD), 3),
      `Untreated` = round(c(c0$CART, c0$SYN, c0$DISCH, c0$DCALC, c0$ANKY, c0$TENDI, c0$BMD), 3)),
      options = list(dom = "t"), rownames = FALSE)
  })

  ## ------------------------------------------------------------------ ⑧ organ
  output$p_valve <- renderPlot({
    d <- sims()$both %>% select(AGEYo, scenario, PMAXS, AVA, LVMI) %>%
      pivot_longer(-c(AGEYo, scenario))
    lab <- c(PMAXS = "Peak pressure gradient (mmHg)", AVA = "Valve area (cm2)", LVMI = "Excess LV mass")
    d$name <- factor(lab[d$name], levels = lab)
    ggplot(d, aes(AGEYo, value, colour = scenario)) +
      geom_line(linewidth = 0.85) + facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = c(control = "#8d99ae", treated = "#e63946")) +
      labs(title = "The aortic valve — a route deliberately left out of the reversible channel",
           subtitle = "In SONIA 2 nitisinone did not slow valve progression (p=0.53)",
           x = "Age (years)", y = NULL, colour = NULL) + THEME
  })

  output$p_organ <- renderPlot({
    d <- sims()$both %>% select(AGEYo, scenario, STONE, PSTONE, RFUN, HEAR, SKINP, TENDI) %>%
      pivot_longer(-c(AGEYo, scenario))
    lab <- c(STONE = "Renal stone burden", PSTONE = "Prostatic stones", RFUN = "Relative renal plasma flow",
             HEAR = "Hearing threshold shift (dB)", SKINP = "Visible pigment", TENDI = "Tendon integrity")
    d$name <- factor(lab[d$name], levels = lab)
    ggplot(d, aes(AGEYo, value, colour = scenario)) +
      geom_line(linewidth = 0.85) + facet_wrap(~name, scales = "free_y", ncol = 3) +
      scale_colour_manual(values = c(control = "#8d99ae", treated = "#e63946")) +
      labs(title = "Outcomes by organ — declining renal function feeds back by raising plasma HGA",
           x = "Age (years)", y = NULL, colour = NULL) + THEME
  })

  ## ----------------------------------------------------------------- ⑨ AKUSSI
  output$p_akussi <- renderPlot({
    d <- sims()$both %>% select(AGEYo, scenario, CAKUSSI, AKJ, AKS, AKC, PAINVAS) %>%
      pivot_longer(-c(AGEYo, scenario))
    lab <- c(CAKUSSI = "cAKUSSI (total)", AKJ = "Joint domain", AKS = "Spine domain",
             AKC = "Clinical domain", PAINVAS = "Pain VAS")
    d$name <- factor(lab[d$name], levels = lab)
    ggplot(d, aes(AGEYo, value, colour = scenario)) +
      geom_line(linewidth = 0.85) + facet_wrap(~name, scales = "free_y", ncol = 3) +
      geom_vline(xintercept = input$start, linetype = 3, colour = "grey40") +
      scale_colour_manual(values = c(control = "#8d99ae", treated = "#e63946")) +
      labs(title = "cAKUSSI and its component domains",
           subtitle = "The absolute values cannot be compared directly with the clinical scale — compare slopes only",
           x = "Age (years)", y = NULL, colour = NULL) + THEME
  })

  output$p_hazard <- renderPlot({
    d <- sims()$both %>% select(AGEYo, scenario, PJR, PAVR, PRUP, PKER) %>%
      pivot_longer(-c(AGEYo, scenario))
    lab <- c(PJR = "Joint replacement", PAVR = "Aortic valve replacement", PRUP = "Tendon rupture",
             PKER = "Keratopathy")
    d$name <- factor(lab[d$name], levels = lab)
    ggplot(d, aes(AGEYo, value, colour = name, linetype = scenario)) +
      geom_line(linewidth = 0.85) +
      scale_colour_manual(values = PAL[c(1,2,3,6)]) +
      labs(title = "Cumulative event probabilities",
           subtitle = "Only keratopathy rises in the treated arm — everything else falls",
           x = "Age (years)", y = "Cumulative probability", colour = NULL, linetype = NULL) + THEME
  })

  output$t_end <- renderDT({
    f <- fin(); c0 <- finc()
    datatable(tibble(
      `Endpoint` = c("cAKUSSI", "Pain VAS", "Joint replacement probability", "Valve replacement probability",
                     "Tendon rupture probability", "Keratopathy probability", "Fraction of pigment spared",
                     "HGA exposure avoided (umol/L*day)"),
      `Treated` = round(c(f$CAKUSSI, f$PAINVAS, f$PJR, f$PAVR, f$PRUP, f$PKER,
                     f$PIGSPARE, f$AVOIDI), 3),
      `Untreated` = round(c(c0$CAKUSSI, c0$PAINVAS, c0$PJR, c0$PAVR, c0$PRUP,
                       c0$PKER, 0, 0), 3)),
      options = list(dom = "t"), rownames = FALSE)
  })

  ## ---------------------------------------------------------------- ⑩ kerato
  output$p_kerato <- renderPlot({
    d <- sims()$both %>% select(AGEYo, scenario, CTYRo, CORTYR, PKER) %>%
      pivot_longer(-c(AGEYo, scenario))
    lab <- c(CTYRo = "Plasma tyrosine (umol/L)", CORTYR = "Corneal crystal burden",
             PKER = "Cumulative keratopathy probability")
    d$name <- factor(lab[d$name], levels = lab)
    ggplot(d, aes(AGEYo, value, colour = scenario)) +
      geom_line(linewidth = 0.85) + facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = c(control = "#8d99ae", treated = "#e63946")) +
      labs(title = "Keratopathy is a question of dietary management, not of dose",
           subtitle = "Lower the protein intake slider on the left and all three panels come down",
           x = "Age (years)", y = NULL, colour = NULL) + THEME
  })

  output$p_safety <- renderPlot({
    grid <- expand.grid(dose = c(2, 10, 20), prot = c(42, 56, 70, 84, 105))
    d <- bind_rows(Map(function(dd, pp) {
      b <- probe_biochem(AKU_MOD, dd, 3, input$start,
                         modifyList(pars(), list(PROT = pp)))
      tibble(dose = dd, protein = pp, Tyr = b$CTYRo, crystal = b$CORTYR)
    }, grid$dose, grid$prot))
    ggplot(d, aes(protein, Tyr, colour = factor(dose))) +
      geom_line(linewidth = 0.9) + geom_point(size = 2.2) +
      geom_hline(yintercept = c(700, 900), linetype = 3) +
      scale_colour_manual(values = PAL[1:3], name = "mg/day") +
      labs(title = "Protein intake vs tyrosine, at three doses",
           subtitle = "The three lines almost coincide — this is the most testable prediction of the model",
           x = "Dietary protein (g/day)", y = "Plasma tyrosine (umol/L)") + THEME
  })

  ## --------------------------------------------------------------- ⑪ headroom
  HR <- eventReactive(input$go_hr, {
    withProgress(message = "Scanning 13 ages at initiation...", value = 0.3,
                 analysis_headroom(AKU_MOD, pars()))
  })

  output$p_headroom <- renderPlot({
    d <- HR() %>% select(init_age, spared_frac, cart_70, p_joint_repl) %>%
      pivot_longer(-init_age)
    lab <- c(spared_frac = "Fraction of pigment spared", cart_70 = "Intact cartilage at 70",
             p_joint_repl = "Joint replacement probability at 70")
    d$name <- factor(lab[d$name], levels = lab)
    ggplot(d, aes(init_age, value, colour = name)) +
      geom_line(linewidth = 1) + geom_point(size = 2.4) +
      scale_colour_manual(values = PAL[c(5,3,1)]) +
      labs(title = "Same drug, same dose — only the remaining integral differs",
           subtitle = "The slope of every curve is the 'cost of waiting'",
           x = "Age at start of treatment (years)", y = NULL, colour = NULL) + THEME
  })

  output$p_headroom2 <- renderPlot({
    ggplot(HR(), aes(init_age, cAKUSSI_70)) +
      geom_line(colour = PAL[1], linewidth = 1) + geom_point(size = 2.4, colour = PAL[1]) +
      labs(title = "cAKUSSI at 70 vs age at initiation",
           subtitle = "The further right along the horizontal axis, the closer it converges on the untreated trajectory",
           x = "Age at start of treatment (years)", y = "cAKUSSI at 70") + THEME
  })

  output$t_headroom <- renderDT({
    datatable(HR() %>% mutate(across(where(is.numeric), ~round(.x, 3))),
              options = list(pageLength = 13), rownames = FALSE)
  })

  ## -------------------------------------------------------------- ⑫ scenarios
  SC <- eventReactive(input$go_scn, {
    withProgress(message = "Running 24 scenarios...", value = 0.2,
                 run_scenarios(AKU_MOD))
  })
  VAL <- eventReactive(input$go_val, {
    withProgress(message = "Running the held-out validation...", value = 0.3,
                 validate_aku(AKU_MOD))
  })

  output$p_scn <- renderPlot({
    d <- SC() %>% filter(id %in% c("S01","S09","S10","S11","S12","S13","S23","S24"))
    ggplot(d, aes(AGEYo, CAKUSSI, colour = id)) +
      geom_line(linewidth = 0.9) +
      labs(title = "Age at initiation · stopping · ideal comparator scenarios",
           subtitle = "S01 untreated · S09-S13 initiation at ages 5/15/25/40/55 · S23 stop at 45 · S24 IDEAL",
           x = "Age (years)", y = "cAKUSSI", colour = NULL) + THEME
  })

  output$t_val <- renderDT({
    datatable(VAL() %>% mutate(across(where(is.numeric), ~signif(.x, 3))),
              options = list(pageLength = 21), rownames = FALSE)
  })

  output$t_scn <- renderDT({
    d <- SC() %>% group_by(id, scenario) %>% slice_tail(n = 1) %>% ungroup() %>%
      select(id, scenario, AGEYo, CTYRo, CHGAo, UHGA24, CAKUSSI, PIGTOT, PJR, PKER) %>%
      mutate(across(where(is.numeric), ~signif(.x, 3)))
    datatable(d, options = list(pageLength = 24, scrollX = TRUE), rownames = FALSE)
  })
}

shinyApp(ui, server)
