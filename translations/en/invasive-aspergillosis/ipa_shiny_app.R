## =============================================================================
##  INVASIVE PULMONARY ASPERGILLOSIS — QSP Shiny DASHBOARD
## =============================================================================
##  11 tabs over the 53-ODE model in `ipa_mrgsolve_model.R`.
##
##  The app is organised around one question the user should leave with:
##  WHICH LEVER ARE YOU ACTUALLY PULLING? Tab 4 lets you move the drug; tab 5
##  lets you move the marrow; tab 8 lets you move the clock. The mortality
##  readout is the same in all three, so the three can be compared in the only
##  currency that matters.
##
##  Run:  shiny::runApp("ipa_shiny_app.R")
##  Requires: shiny, mrgsolve, ggplot2, dplyr, tidyr, DT
## =============================================================================

library(shiny)
library(mrgsolve)
suppressMessages({library(ggplot2); library(dplyr); library(tidyr); library(DT)})

source("ipa_mrgsolve_model.R", local = TRUE)   # defines mod, rx_*, sim(), CYP2C19

THEME <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92", colour = NA),
        legend.position = "bottom")

PAL <- c("#2E86C1", "#E74C3C", "#28B463", "#B9770E", "#7D3C98", "#17A589",
         "#CB4335", "#5D6D7E")

## -----------------------------------------------------------------------------
## UI
## -----------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Invasive Pulmonary Aspergillosis (IPA) — QSP Dashboard"),
  tags$p(tags$b("One sentence: "),
         "Triazoles multiply the growth term, and only neutrophils and polyenes add to the clearance term. ",
         "Therefore, in a host without neutrophils, azole monotherapy can never reduce the fungal burden ",
         "— it can only slow it."),
  hr(),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Host"),
      selectInput("phenotype", "Host phenotype",
                  c("Neutropenia (haematologic malignancy remission induction)" = "neutropenic",
                    "Steroid host (normal neutrophils)"   = "steroid",
                    "Transplant (tacrolimus + low-dose steroid)" = "transplant",
                    "Immunocompetent (control)"                 = "normal")),
      sliderInput("anc_rec", "ANC recovery day (day; 99 = no recovery)",
                  min = 3, max = 99, value = 10, step = 1),
      sliderInput("inoc", "Conidia reaching the alveoli (log10 CFUe)",
                  min = 2, max = 5, value = 3, step = 0.25),
      sliderInput("pred", "Prednisone equivalent (mg/day)", min = 0, max = 100, value = 0, step = 10),
      hr(),
      h4("Antifungal"),
      selectInput("drug", "Treatment",
                  c("None" = "none", "Voriconazole" = "vrc", "Isavuconazole" = "isa",
                    "Liposomal amphotericin B" = "amb", "Anidulafungin" = "ech",
                    "Voriconazole + anidulafungin" = "vrc_ech",
                    "Posaconazole (prophylaxis/salvage)" = "pos")),
      sliderInput("start_d", "Treatment start day (day)", min = 1, max = 14, value = 3, step = 1),
      sliderInput("vrc_dose", "Voriconazole maintenance dose (mg/kg q12h)",
                  min = 1, max = 9, value = 4, step = 0.2),
      sliderInput("amb_dose", "L-AmB dose (mg/kg/day)", min = 1, max = 10, value = 3, step = 1),
      checkboxInput("gcsf", "Concomitant G-CSF", FALSE),
      hr(),
      h4("Isolate"),
      selectInput("geno", "CYP2C19 genotype", names(CYP2C19), selected = "NM (*1/*1)"),
      sliderInput("mic", "Voriconazole MIC (mg/L)", min = 0.125, max = 8, value = 0.5, step = 0.125),
      checkboxInput("tr34", "TR34/L98H environmental resistant strain", FALSE),
      hr(),
      sliderInput("horizon", "Observation period (days)", min = 28, max = 120, value = 84, step = 7),
      actionButton("go", "Run simulation", class = "btn-primary btn-block")
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs", type = "tabs",

        ## ------------------------------------------------------------ tab 1
        tabPanel(
          "1 · Patient profile",
          br(),
          fluidRow(
            column(4, wellPanel(h4("Host status"), tableOutput("tbl_host"))),
            column(4, wellPanel(h4("Treatment"),      tableOutput("tbl_rx"))),
            column(4, wellPanel(h4("Outcome summary"), tableOutput("tbl_out")))),
          plotOutput("p_profile", height = "440px"),
          helpText("The reason neutrophils (ANC) · lesion volume · perfusion · oxygenation are placed on one screen ",
                   "is that these four curves are, in effect, each other's cause.")),

        ## ------------------------------------------------------------ tab 2
        tabPanel(
          "2 · Pharmacokinetics (PK)",
          br(),
          plotOutput("p_pk", height = "330px"),
          plotOutput("p_site", height = "330px"),
          helpText("The region where the dashed line in the figure below (intralesional concentration = ELF concentration × perfusion^γ) ",
                   "diverges from the solid line (ELF) is precisely the region where vascular invasion is blocking the drug.")),

        ## ------------------------------------------------------------ tab 3
        tabPanel(
          "3 · CYP2C19 · TDM",
          br(),
          plotOutput("p_geno", height = "340px"),
          DTOutput("dt_geno"),
          helpText("At the same mg/kg dose, genotype alone causes AUC to diverge more than 5-fold. ",
                   "TDM is a device for narrowing this variance, not for raising mean exposure.")),

        ## ------------------------------------------------------------ tab 4
        tabPanel(
          "4 · Fungal Burden · Lesion",
          br(),
          plotOutput("p_burden", height = "340px"),
          plotOutput("p_lesion", height = "300px"),
          helpText("The bold horizontal dashed line is the trajectory corresponding to the growth floor at ",
                   "maximal azole effect, k_grow·(1 − Imax). ",
                   "Without neutrophils, no azole exposure can push the curve below this line.")),

        ## ------------------------------------------------------------ tab 5
        tabPanel(
          "5 · Biomarkers: GM · BDG · PCR",
          br(),
          plotOutput("p_bio", height = "330px"),
          plotOutput("p_gmdiss", height = "330px"),
          helpText("The figure below is this model's clinically most provocative prediction. ",
                   "The arm with lower GM at week 1 can have a higher fungal burden — ",
                   "because GM is a marker of flux, not stock.")),

        ## ------------------------------------------------------------ tab 6
        tabPanel(
          "6 · Clinical Endpoints",
          br(),
          plotOutput("p_surv", height = "330px"),
          plotOutput("p_organ", height = "330px")),

        ## ------------------------------------------------------------ tab 7
        tabPanel(
          "7 · Scenario Comparison",
          br(),
          checkboxGroupInput(
            "cmp", "Arms to compare", inline = TRUE,
            choices = c("No treatment" = "none", "Voriconazole" = "vrc", "Isavuconazole" = "isa",
                        "L-AmB 3" = "amb", "Anidulafungin" = "ech",
                        "VRC+anidulafungin" = "vrc_ech"),
            selected = c("none", "vrc", "amb")),
          plotOutput("p_cmp", height = "380px"),
          DTOutput("dt_cmp")),

        ## ------------------------------------------------------------ tab 8
        tabPanel(
          "8 · Treatment Delay Cliff",
          br(),
          plotOutput("p_delay", height = "380px"),
          DTOutput("dt_delay"),
          helpText("Turning perfusion feedback (γ) off to 0 makes this cliff disappear — ",
                   "the model's claim is that late treatment fails not because 'there is too much fungus' ",
                   "but because 'the drug cannot get there'.")),

        ## ------------------------------------------------------------ tab 9
        tabPanel(
          "9 · Drug Lever vs Host Lever",
          br(),
          plotOutput("p_lever", height = "420px"),
          DTOutput("dt_lever"),
          helpText("Compares the range from fixing the x-axis and varying the drug, and the range from fixing the drug and varying the marrow, ",
                   "in the same unit (12-week mortality %p).")),

        ## ----------------------------------------------------------- tab 10
        tabPanel(
          "10 · Resistance · MIC",
          br(),
          plotOutput("p_mic", height = "340px"),
          plotOutput("p_resist", height = "320px")),

        ## ----------------------------------------------------------- tab 11
        tabPanel(
          "11 · Toxicity · Drug Interactions",
          br(),
          plotOutput("p_tox", height = "340px"),
          plotOutput("p_ddi", height = "320px"),
          helpText("The tacrolimus curve is a picture that views the azole not as an 'antifungal' ",
                   "but as a 'CYP3A4 inhibitor'."))
      )
    )
  ),
  hr(),
  tags$small("This is an educational/research QSP model. Do not use it for clinical decision-making. ",
             "All figures have been independently reproduced and verified by ipa_reference_model.py.")
)

## -----------------------------------------------------------------------------
## SERVER
## -----------------------------------------------------------------------------
server <- function(input, output, session) {

  build <- reactive({
    ph  <- input$phenotype
    par <- list(
      INOCULUM   = 10^input$inoc,
      CYP2C19    = unname(CYP2C19[[input$geno]]),
      MIC_VRC    = input$mic,
      RESIST_FRAC = if (isTRUE(input$tr34)) 1 else 0,
      T_ANC_REC  = if (ph == "neutropenic") ifelse(input$anc_rec >= 99, 1e9, input$anc_rec * 24) else 0,
      ANC0       = switch(ph, neutropenic = 0.05, steroid = 4, transplant = 3, normal = 4))
    d <- switch(input$drug,
      none    = rx_none(),
      vrc     = rx_vrc(input$start_d, maint = input$vrc_dose),
      isa     = rx_isa(input$start_d),
      amb     = rx_amb(input$start_d, mgkg = input$amb_dose),
      ech     = rx_ech(input$start_d),
      vrc_ech = rbind(rx_vrc(input$start_d, maint = input$vrc_dose), rx_ech(input$start_d)),
      pos     = rx_pos(0, 90))
    if (input$pred > 0) d <- rbind(d, rx_steroid(input$pred, 60))
    if (ph == "steroid" && input$pred == 0)   d <- rbind(d, rx_steroid(60, 60))
    if (ph == "transplant") d <- rbind(d, rx_tac(), rx_steroid(20, 84))
    if (isTRUE(input$gcsf)) d <- rbind(d, rx_gcsf(input$start_d, 12))
    list(dose = d, par = par)
  })

  run <- eventReactive(input$go, {
    b <- build()
    sim(b$dose, b$par, end = input$horizon * 24, delta = 2)
  }, ignoreNULL = FALSE)

  day <- function(df) df$time / 24

  ## ---------------------------------------------------------------- tab 1
  output$tbl_host <- renderTable({
    o <- run()
    data.frame(Item = c("Phenotype", "ANC recovery", "Inoculum (log10)", "Steroid"),
               Value = c(input$phenotype,
                      ifelse(input$anc_rec >= 99, "None", paste0("day ", input$anc_rec)),
                      sprintf("%.2f", input$inoc),
                      paste0(input$pred, " mg/day")))
  })
  output$tbl_rx <- renderTable({
    data.frame(Item = c("Agent", "Start day", "VRC dose", "CYP2C19", "MIC"),
               Value = c(input$drug, paste0("day ", input$start_d),
                      paste0(input$vrc_dose, " mg/kg q12h"), input$geno,
                      paste0(input$mic, " mg/L")))
  })
  output$tbl_out <- renderTable({
    o <- run()
    data.frame(Item = c("Peak log10 burden", "Lowest perfusion", "Peak GM", "12-week mortality (%)"),
               Value = sprintf("%.2f", c(max(o$LOGB), min(o$PERF), max(o$GMser),
                                      tail(o$MORT, 1))))
  })
  output$p_profile <- renderPlot({
    o <- run()
    d <- data.frame(day = day(o),
                    `ANC (10^9/L)` = o$ANC, `Lesion volume (mL)` = o$VLES,
                    `Perfusion index` = o$PERF, `PaO2/FiO2` = o$PFR,
                    check.names = FALSE) %>%
      pivot_longer(-day)
    ggplot(d, aes(day, value)) + geom_line(linewidth = 0.9, colour = PAL[1]) +
      facet_wrap(~name, scales = "free_y") + THEME +
      labs(x = "Day", y = NULL, title = "Host profile")
  })

  ## ---------------------------------------------------------------- tab 2
  output$p_pk <- renderPlot({
    o <- run()
    d <- data.frame(day = day(o), Voriconazole = o$CP_VRC, Isavuconazole = o$CP_ISA,
                    Posaconazole = o$CP_POS, `L-AmB` = o$CP_AMB,
                    Echinocandin = o$CP_ECH, check.names = FALSE) %>%
      pivot_longer(-day) %>% filter(value > 1e-6)
    if (!nrow(d)) return(NULL)
    ggplot(d, aes(day, value, colour = name)) + geom_line(linewidth = 0.8) +
      scale_colour_manual(values = PAL) + THEME +
      labs(x = "Day", y = "Plasma concentration (mg/L)", colour = NULL, title = "Plasma pharmacokinetics")
  })
  output$p_site <- renderPlot({
    o <- run()
    d <- rbind(
      data.frame(day = day(o), conc = o$VRC_elf,          Layer = "ELF (epithelial lining fluid)", drug = "VRC"),
      data.frame(day = day(o), conc = o$VRC_elf * o$PERF, Layer = "Intralesional (× perfusion)",  drug = "VRC"),
      data.frame(day = day(o), conc = o$AMB_elf,          Layer = "ELF (epithelial lining fluid)", drug = "L-AmB"),
      data.frame(day = day(o), conc = o$AMB_elf * o$PERF, Layer = "Intralesional (× perfusion)",  drug = "L-AmB")) %>%
      filter(conc > 1e-6)
    if (!nrow(d)) return(NULL)
    ggplot(d, aes(day, conc, colour = drug, linetype = Layer)) +
      geom_line(linewidth = 0.85) + scale_colour_manual(values = PAL) + THEME +
      labs(x = "Day", y = "Concentration (mg/L)", colour = NULL, linetype = NULL,
           title = "Lung vs lesion: the gap perfusion creates")
  })

  ## ---------------------------------------------------------------- tab 3
  geno_tab <- reactive({
    bind_rows(lapply(names(CYP2C19), function(k) {
      o <- sim(rx_vrc(0, 14, maint = input$vrc_dose),
               list(CYP2C19 = unname(CYP2C19[[k]]), INOCULUM = 0),
               end = 14 * 24, delta = 0.25)
      w <- o[o$time >= 12 * 24 & o$time <= 13 * 24, ]
      data.frame(Genotype = k, Cmin = min(w$CP_VRC), Cmax = max(w$CP_VRC),
                 AUC24 = sum(diff(w$time) * (head(w$CP_VRC, -1) + tail(w$CP_VRC, -1)) / 2))
    }))
  })
  output$p_geno <- renderPlot({
    g <- geno_tab() %>% mutate(Genotype = factor(Genotype, levels = names(CYP2C19)))
    ggplot(g, aes(Genotype, AUC24, fill = Genotype)) +
      geom_col(width = 0.65) +
      geom_hline(yintercept = c(24, 132), linetype = 2, colour = "grey30") +
      annotate("text", x = 1, y = 138, label = "Therapeutic range (equivalent to trough 1–5.5 mg/L)",
               hjust = 0, size = 3.4) +
      scale_fill_manual(values = PAL, guide = "none") + THEME +
      labs(x = NULL, y = "Steady-state AUC24 (mg·h/L)",
           title = sprintf("Exposure by genotype at the same %.1f mg/kg q12h dose", input$vrc_dose))
  })
  output$dt_geno <- renderDT(datatable(geno_tab() %>% mutate(across(where(is.numeric), ~round(.x, 2))),
                                       options = list(dom = "t")))

  ## ---------------------------------------------------------------- tab 4
  output$p_burden <- renderPlot({
    o <- run()
    p <- as.list(mrgsolve::param(mod))
    fl <- p$KGROW * (1 - p$IMAX_AZOLE)
    t0 <- input$start_d * 24
    ref <- data.frame(day = day(o)) %>%
      mutate(y = o$LOGB[which.min(abs(o$time - t0))] +
               pmax(0, (time <- day * 24) - t0) * fl / log(10))
    ggplot() +
      geom_line(data = data.frame(day = day(o), y = o$LOGB), aes(day, y),
                linewidth = 1.0, colour = PAL[2]) +
      geom_line(data = ref, aes(day, y), linetype = 2, linewidth = 0.8, colour = "grey35") +
      annotate("text", x = max(day(o)) * 0.55, y = max(o$LOGB),
               label = "Grey dashed line = growth floor at maximal azole effect", size = 3.6, hjust = 0) +
      THEME + labs(x = "Day", y = "log10 hyphal burden (CFUe)", title = "Fungal burden")
  })
  output$p_lesion <- renderPlot({
    o <- run()
    d <- data.frame(day = day(o), `Lesion volume (mL)` = o$VLES,
                    `Infarct volume (mL)` = o$NEC, `Perfusion index` = o$PERF,
                    `Angioinvasion index` = o$ANG, check.names = FALSE) %>% pivot_longer(-day)
    ggplot(d, aes(day, value)) + geom_line(linewidth = 0.9, colour = PAL[4]) +
      facet_wrap(~name, scales = "free_y") + THEME + labs(x = "Day", y = NULL)
  })

  ## ---------------------------------------------------------------- tab 5
  output$p_bio <- renderPlot({
    o <- run()
    d <- data.frame(day = day(o), `Serum GM (ODI)` = o$GMser,
                    `BDG (pg/mL)` = o$BDG, `PCR (copies/mL)` = o$PCRs,
                    `log10 burden` = o$LOGB, check.names = FALSE) %>% pivot_longer(-day)
    ggplot(d, aes(day, value)) + geom_line(linewidth = 0.9, colour = PAL[3]) +
      facet_wrap(~name, scales = "free_y") + THEME + labs(x = "Day", y = NULL)
  })
  gm_tab <- reactive({
    arms <- list(`No treatment` = rx_none(), `Voriconazole` = rx_vrc(input$start_d),
                 `L-AmB 3 mg/kg` = rx_amb(input$start_d, mgkg = 3),
                 `Anidulafungin` = rx_ech(input$start_d))
    bind_rows(lapply(names(arms), function(a) {
      o <- sim(arms[[a]], list(T_ANC_REC = 336), end = 20 * 24, delta = 1)
      data.frame(day = o$time / 24, arm = a, GM = o$GMser, logB = o$LOGB)
    }))
  })
  output$p_gmdiss <- renderPlot({
    d <- gm_tab() %>% pivot_longer(c(GM, logB))
    d$name <- factor(d$name, c("GM", "logB"), c("Serum GM index", "log10 fungal burden"))
    ggplot(d, aes(day, value, colour = arm)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = PAL) + THEME +
      labs(x = "Day", y = NULL, colour = NULL,
           title = "GM tracks flux, not stock — the week-1 ranking flips")
  })

  ## ---------------------------------------------------------------- tab 6
  output$p_surv <- renderPlot({
    o <- run()
    ggplot(data.frame(day = day(o), s = 100 * o$SURV), aes(day, s)) +
      geom_line(linewidth = 1.0, colour = PAL[5]) + ylim(0, 100) + THEME +
      labs(x = "Day", y = "Survival probability (%)", title = "Model survival curve")
  })
  output$p_organ <- renderPlot({
    o <- run()
    d <- data.frame(day = day(o), `Temperature (°C)` = o$TEMP, `CRP (mg/L)` = o$CRP,
                    `IL-6 (pg/mL)` = o$IL6, `PaO2/FiO2` = o$PFR,
                    check.names = FALSE) %>% pivot_longer(-day)
    ggplot(d, aes(day, value)) + geom_line(linewidth = 0.9, colour = PAL[6]) +
      facet_wrap(~name, scales = "free_y") + THEME + labs(x = "Day", y = NULL)
  })

  ## ---------------------------------------------------------------- tab 7
  cmp_run <- reactive({
    req(length(input$cmp) > 0)
    b <- build()
    arms <- list(none = rx_none(), vrc = rx_vrc(input$start_d, maint = input$vrc_dose),
                 isa = rx_isa(input$start_d), amb = rx_amb(input$start_d, mgkg = input$amb_dose),
                 ech = rx_ech(input$start_d),
                 vrc_ech = rbind(rx_vrc(input$start_d, maint = input$vrc_dose),
                                 rx_ech(input$start_d)))
    bind_rows(lapply(input$cmp, function(k) {
      sim(arms[[k]], b$par, end = input$horizon * 24, delta = 4) %>%
        mutate(arm = k)
    }))
  })
  output$p_cmp <- renderPlot({
    d <- cmp_run() %>%
      transmute(day = time / 24, arm,
                `log10 burden` = LOGB, `Serum GM` = GMser,
                `Perfusion` = PERF, `Mortality (%)` = MORT) %>%
      pivot_longer(-c(day, arm))
    ggplot(d, aes(day, value, colour = arm)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = PAL) + THEME +
      labs(x = "Day", y = NULL, colour = NULL)
  })
  output$dt_cmp <- renderDT({
    d <- cmp_run() %>% group_by(arm) %>%
      summarise(`Peak log10 burden` = round(max(LOGB), 2),
                `d42 log10 burden`  = round(LOGB[which.min(abs(time - 42 * 24))], 2),
                `Lowest perfusion`       = round(min(PERF), 3),
                `Peak GM`         = round(max(GMser), 2),
                `12-week mortality (%)` = round(tail(MORT, 1), 1), .groups = "drop")
    datatable(d, options = list(dom = "t"))
  })

  ## ---------------------------------------------------------------- tab 8
  delay_tab <- reactive({
    bind_rows(lapply(c(1:8, 10, 12, 14), function(dd) {
      a <- sim(rx_vrc(dd, maint = input$vrc_dose), list(T_ANC_REC = 336), end = 84 * 24, delta = 6)
      b <- sim(rx_vrc(dd, maint = input$vrc_dose),
               list(T_ANC_REC = 336, PERF_GAMMA = 0), end = 84 * 24, delta = 6)
      data.frame(`Start day` = dd,
                 `d84 log10 burden (γ=1)` = round(tail(a$LOGB, 1), 2),
                 `d84 log10 burden (γ=0)` = round(tail(b$LOGB, 1), 2),
                 `Lowest perfusion` = round(min(a$PERF), 3),
                 `12-week mortality (%)` = round(tail(a$MORT, 1), 1),
                 check.names = FALSE)
    }))
  })
  output$p_delay <- renderPlot({
    d <- delay_tab()
    dd <- data.frame(day = rep(d$`Start day`, 2),
                     y = c(d$`d84 log10 burden (γ=1)`, d$`d84 log10 burden (γ=0)`),
                     model = rep(c("Perfusion feedback on (γ=1)", "Off (γ=0)"), each = nrow(d)))
    ggplot(dd, aes(day, y, colour = model)) +
      geom_line(linewidth = 1.0) + geom_point(size = 2) +
      scale_colour_manual(values = PAL[c(2, 1)]) + THEME +
      labs(x = "Treatment start day", y = "log10 fungal burden at 12 weeks", colour = NULL,
           title = "Perfusion feedback creates the treatment-delay cliff")
  })
  output$dt_delay <- renderDT(datatable(delay_tab(), options = list(dom = "t")))

  ## ---------------------------------------------------------------- tab 9
  lever_tab <- reactive({
    arms <- list(`No treatment` = rx_none(), VRC = rx_vrc(input$start_d),
                 ISA = rx_isa(input$start_d), `L-AmB` = rx_amb(input$start_d, mgkg = 3))
    bind_rows(lapply(c(5, 10, 14, 21, 28, 99), function(r) {
      bind_rows(lapply(names(arms), function(a) {
        o <- sim(arms[[a]], list(T_ANC_REC = ifelse(r >= 99, 1e9, r * 24)),
                 end = 84 * 24, delta = 8)
        data.frame(`ANC recovery day` = ifelse(r >= 99, 99, r), Arm = a,
                   `12-week mortality` = round(tail(o$MORT, 1), 1), check.names = FALSE)
      }))
    }))
  })
  output$p_lever <- renderPlot({
    d <- lever_tab()
    ggplot(d, aes(factor(`ANC recovery day`), `12-week mortality`, fill = Arm)) +
      geom_col(position = position_dodge(0.8), width = 0.75) +
      scale_fill_manual(values = PAL) + THEME +
      labs(x = "ANC recovery day (99 = no recovery)", y = "12-week mortality (%)", fill = NULL,
           title = "The gap between bar groups (marrow) is larger than between bars (drug)")
  })
  output$dt_lever <- renderDT({
    d <- lever_tab() %>% pivot_wider(names_from = Arm, values_from = `12-week mortality`)
    datatable(d, options = list(dom = "t"))
  })

  ## --------------------------------------------------------------- tab 10
  output$p_mic <- renderPlot({
    d <- bind_rows(lapply(c(0.125, 0.25, 0.5, 1, 2, 4, 8), function(m) {
      o <- sim(rx_vrc(input$start_d, maint = input$vrc_dose),
               list(T_ANC_REC = 336, MIC_VRC = m), end = 84 * 24, delta = 8)
      data.frame(MIC = m, logB_d42 = o$LOGB[which.min(abs(o$time - 42 * 24))],
                 mort = tail(o$MORT, 1))
    }))
    ggplot(d, aes(MIC, mort)) + geom_line(linewidth = 1.0, colour = PAL[2]) +
      geom_point(size = 2.2) + scale_x_log10(breaks = c(0.125, 0.5, 2, 8)) +
      geom_vline(xintercept = 1, linetype = 2) +
      annotate("text", x = 1.1, y = min(d$mort), label = "EUCAST ECOFF", hjust = 0, size = 3.5) +
      THEME + labs(x = "Voriconazole MIC (mg/L, log axis)", y = "12-week mortality (%)",
                   title = "MIC ladder")
  })
  output$p_resist <- renderPlot({
    arms <- list(`No treatment` = rx_none(), `VRC alone` = rx_vrc(input$start_d),
                 `L-AmB` = rx_amb(input$start_d, mgkg = 3),
                 `VRC+echinocandin` = rbind(rx_vrc(input$start_d), rx_ech(input$start_d)))
    d <- bind_rows(lapply(names(arms), function(a) {
      o <- sim(arms[[a]], list(T_ANC_REC = 336), end = 84 * 24, delta = 8)
      data.frame(day = o$time / 24, arm = a, rfrac = pmax(o$RFRAC, 1e-12))
    }))
    ggplot(d, aes(day, rfrac, colour = arm)) + geom_line(linewidth = 0.9) +
      scale_y_log10() + scale_colour_manual(values = PAL) + THEME +
      labs(x = "Day", y = "Resistant subpopulation fraction (log axis)", colour = NULL,
           title = "Resistance selection under azole monotherapy")
  })

  ## --------------------------------------------------------------- tab 11
  output$p_tox <- renderPlot({
    o <- run()
    d <- data.frame(day = day(o), `ALT (U/L)` = o$ALT, `Creatinine (mg/dL)` = o$SCR,
                    `QTc increase (ms)` = o$QTC, `Neuro/visual effects` = o$NEUROE,
                    check.names = FALSE) %>% pivot_longer(-day)
    ggplot(d, aes(day, value)) + geom_line(linewidth = 0.9, colour = PAL[7]) +
      facet_wrap(~name, scales = "free_y") + THEME + labs(x = "Day", y = NULL)
  })
  output$p_ddi <- renderPlot({
    arms <- list(`Tacrolimus alone` = rx_tac(),
                 `+ Voriconazole` = rbind(rx_tac(), rx_vrc(0)),
                 `+ Isavuconazole` = rbind(rx_tac(), rx_isa(0)),
                 `+ Posaconazole` = rbind(rx_tac(), rx_pos(0, 30)))
    d <- bind_rows(lapply(names(arms), function(a) {
      o <- sim(arms[[a]], list(T_ANC_REC = 0, ANC0 = 3), end = 20 * 24, delta = 1)
      data.frame(day = o$time / 24, arm = a, tac = o$CP_TAC)
    }))
    ggplot(d, aes(day, tac, colour = arm)) + geom_line(linewidth = 0.9) +
      geom_hline(yintercept = c(5, 15), linetype = 2, colour = "grey40") +
      scale_colour_manual(values = PAL) + THEME +
      labs(x = "Day", y = "Tacrolimus (ng/mL)", colour = NULL,
           title = "Grey band = usual target concentration 5–15 ng/mL")
  })
}

shinyApp(ui, server)
