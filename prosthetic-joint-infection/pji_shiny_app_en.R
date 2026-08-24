## =============================================================================
##  pji_shiny_app.R — Prosthetic joint infection (PJI) QSP dashboard
##  Prosthetic Joint Infection — interactive QSP dashboard
##
##  Run:
##      # from the same directory
##      shiny::runApp("pji_shiny_app.R")
##  Required packages: shiny, mrgsolve, ggplot2, dplyr, tidyr, DT
##
##  The app loads the model from pji_mrgsolve_model_en.R unchanged and shows it in 13 tabs.
##  Tab 3 (bone concentration / MBEC ratio) and tab 7 (mutant supply) are this model's
##  argument — the other tabs support why those two pictures shape the clinical guidance.
## =============================================================================

library(shiny)
library(mrgsolve)
library(ggplot2)
library(dplyr)
library(tidyr)
library(DT)

`%||%` <- function(a, b) if (is.null(a)) b else a

## ---- Load the model ---------------------------------------------------------
## pji_mrgsolve_model_en.R builds the model object pji_mod and the scenario list scen.
if (!exists("pji_mod")) {
  src <- file.path(dirname(sys.frame(1)$ofile %||% "."), "pji_mrgsolve_model_en.R")
  if (!file.exists(src)) src <- "pji_mrgsolve_model_en.R"
  source(src, local = FALSE)
}

DAY <- 24
THEME <- theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        legend.position  = "bottom",
        plot.title       = element_text(face = "bold", size = 14),
        strip.text       = element_text(face = "bold"))

DRUGCOL <- c(VAN = "#2E6DA4", RIF = "#C0392B", LVX = "#27AE60",
             DAP = "#8E44AD", LZD = "#E67E22", CFZ = "#16A085",
             LOCAL = "#7F8C8D")
POPCOL  <- c("Planktonic"                = "#3498DB",
             "Biofilm"                   = "#E67E22",
             "Persister"                 = "#C0392B",
             "SCV"                       = "#8E44AD",
             "Intracellular"             = "#16A085",
             "rpoB resistant"            = "#000000",
             "gyrA resistant"            = "#7F8C8D")

## =============================================================================
##  UI
## =============================================================================
ui <- fluidPage(
  titlePanel("Prosthetic Joint Infection (PJI) QSP dashboard"),
  tags$p(style = "color:#555;margin-top:-10px;",
         "Implant-associated Staphylococcus aureus osteomyelitis · 53-ODE QSP model · for education and research"),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      h4("① Patient & lesion"),
      radioButtons("organism", "Organism", c("MRSA" = 0, "MSSA" = 1), selected = 0, inline = TRUE),
      sliderInput("imsup", "Host immune coefficient (1 = normal)", 0.2, 1.2, 1.0, step = 0.1),
      numericInput("inoc", "Inoculum (CFU)", 100, min = 1, max = 1e8),
      checkboxInput("implant", "Implant present (foreign body)", TRUE),
      sliderInput("tdx", "Time of diagnosis / surgery (days)", 7, 240, 21, step = 7),

      hr(), h4("② Surgical strategy"),
      selectInput("surg", "Primary surgery",
                  c("none (antibiotics alone)"      = "none",
                    "DAIR (debridement, implant retained)" = "dair",
                    "one-stage exchange"            = "one",
                    "two-stage exchange (+ spacer)" = "two")),
      sliderInput("logk", "Bacterial reduction at surgery (log10)", 0, 7, 2.5, step = 0.5),
      conditionalPanel("input.surg == 'two'",
        sliderInput("spacer", "Elutable vancomycin in the spacer (mg)", 0, 3000, 1900, step = 100),
        sliderInput("reimp", "Reimplantation (weeks after surgery)", 4, 16, 8, step = 1)),

      hr(), h4("③ Antibiotics"),
      checkboxGroupInput("drugs", "Regimen",
        c("Vancomycin IV 1 g q12h"       = "VAN",
          "Rifampicin PO 450 mg q12h"    = "RIF",
          "Levofloxacin PO 750 mg qd"    = "LVX",
          "Daptomycin IV 8 mg/kg qd"     = "DAP",
          "Linezolid PO 600 mg q12h"     = "LZD",
          "Cefazolin IV 2 g q8h (MSSA)"  = "CFZ"),
        selected = c("LVX", "RIF")),
      sliderInput("wk_backbone", "Backbone duration (weeks)", 2, 26, 6, step = 1),
      sliderInput("wk_rif",      "Rifampicin duration (weeks)", 0, 26, 12, step = 1),

      hr(), h4("④ Simulation"),
      sliderInput("tend", "Observation period (days)", 90, 730, 365, step = 30),
      actionButton("go", "Run simulation", class = "btn-primary btn-block"),
      br(),
      tags$small(style = "color:#777;",
        "In summary: apart from rifampicin, the free bone concentration reachable with ",
        "systemic antibiotics is only 1~2% of the MBEC. Surgery does not circumvent that ",
        "fact; it is the step that lowers the residual burden N, dropping the probability ",
        "1−e^(−μN) that an rpoB mutant exists and so making rifampicin usable.")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs", type = "tabs",

        tabPanel("1. Overview",
          br(), fluidRow(
            column(3, wellPanel(h5("Final burden (log10 CFU)"), h3(textOutput("kpi_burden")))),
            column(3, wellPanel(h5("Probability of cure"),      h3(textOutput("kpi_cure")))),
            column(3, wellPanel(h5("P(pre-existing rpoB mutant)"), h3(textOutput("kpi_prpob")))),
            column(3, wellPanel(h5("2018 ICM score"),           h3(textOutput("kpi_icm"))))),
          plotOutput("p_overview", height = "460px"),
          tags$hr(), htmlOutput("narrative")),

        tabPanel("2. Antibiotic PK (plasma vs bone)",
          br(), plotOutput("p_pk", height = "560px"),
          tags$small("Dashed = total plasma concentration, solid = free bone concentration. ",
                     "The gap between the two lines is (bone penetration x free fraction).")),

        tabPanel("3. ★ Bone concentration / MBEC ratio",
          br(), h4("This tab is the argument of this model"),
          plotOutput("p_ratio", height = "380px"),
          br(), DT::dataTableOutput("t_ratio"),
          tags$small("If the ratio never passes 1, that agent cannot eradicate biofilm. ",
                     "Only rifampicin comes close to 1 by the systemic route, and the only thing that ",
                     "passes 1 is the local concentration from antibiotic-loaded cement.")),

        tabPanel("4. Bacterial subpopulations",
          br(), plotOutput("p_pops", height = "520px"),
          checkboxInput("logscale", "log10 axis", TRUE)),

        tabPanel("5. Biofilm · maturation",
          br(), plotOutput("p_biofilm", height = "520px"),
          tags$small("EPS maturity is a state variable that literally multiplies the biofilm ",
                     "tolerance factor (MBEC/MIC). This is why DAIR only holds 'within three weeks'.")),

        tabPanel("6. Local elution (spacer)",
          br(), plotOutput("p_local", height = "480px"),
          tags$small("Antibiotic-loaded PMMA exceeds the MBEC only over the first few days. ",
                     "After that it falls into the same range as systemic therapy.")),

        tabPanel("7. ★ Emergence of resistance / mutant supply",
          br(), plotOutput("p_res", height = "380px"),
          br(), plotOutput("p_mutsupply", height = "300px"),
          tags$small("The figure below is not a simulation but arithmetic: ",
                     "P = 1 − exp(−μN), μ = 1e−8. N before and after debridement is marked.")),

        tabPanel("8. Host immunity · inflammation",
          br(), plotOutput("p_immune", height = "560px")),

        tabPanel("9. Osteolysis · implant loosening",
          br(), plotOutput("p_bone", height = "560px")),

        tabPanel("10. Biomarkers · endpoints",
          br(), plotOutput("p_biomarker", height = "520px"),
          tags$small("Dashed = 2018 ICM diagnostic thresholds (CRP 10 mg/L, synovial WBC 3,000/uL, ",
                     "PMN 80%, alpha-defensin 1.0 S/CO)")),

        tabPanel("11. Scenario comparison",
          br(), actionButton("runall", "Run the 21 standard scenarios", class = "btn-success"),
          br(), br(), DT::dataTableOutput("t_scen"),
          br(), plotOutput("p_scen", height = "460px")),

        tabPanel("12. Toxicity · safety",
          br(), plotOutput("p_tox", height = "520px"),
          tags$small("Nephrotoxicity above vancomycin AUC24 > 600 mg·h/L, ",
                     "thrombocytopenia with prolonged linezolid, raised ALT with rifampicin.")),

        tabPanel("13. Parameters · references",
          br(), DT::dataTableOutput("t_par"),
          br(), htmlOutput("reflink"))
      )
    )
  )
)

## =============================================================================
##  SERVER
## =============================================================================
server <- function(input, output, session) {

  ## ---- Regimen → mrgsolve events --------------------------------------------
  build <- reactive({
    tdx <- input$tdx*DAY
    e <- NULL
    add <- function(e, x) if (is.null(e)) x else e + x
    if ("VAN" %in% input$drugs) e <- add(e, ev_van(tdx, input$wk_backbone))
    if ("LVX" %in% input$drugs) e <- add(e, ev_lvx(tdx, input$wk_backbone))
    if ("DAP" %in% input$drugs) e <- add(e, ev_dap(tdx, input$wk_backbone))
    if ("LZD" %in% input$drugs) e <- add(e, ev_lzd(tdx, input$wk_backbone))
    if ("CFZ" %in% input$drugs) e <- add(e, ev_cfz(tdx, input$wk_backbone))
    if ("RIF" %in% input$drugs && input$wk_rif > 0) e <- add(e, ev_rif(tdx, input$wk_rif))
    if (is.null(e)) e <- ev(time = 0, cmt = "NP", amt = 0)

    p <- list(MSSA = as.numeric(input$organism),
              IMSUP = input$imsup,
              IMPL = as.numeric(input$implant),
              INOC0 = input$inoc)
    if (input$surg != "none") {
      p$TSURG1 <- tdx
      p$LOGK1  <- input$logk
      p$EXCH1  <- ifelse(input$surg == "dair", 0, 1)
    }
    if (input$surg == "two") {
      p$TSURG2 <- tdx + input$reimp*7*DAY
      p$LOGK2  <- 1.5
      p$EXCH2  <- 1
      p$SPCF0  <- 0.63*input$spacer
      p$SPCS0  <- 0.37*input$spacer
    }
    list(ev = e, par = p, tdx = tdx)
  })

  sim <- eventReactive(input$go, {
    b <- build()
    m <- mrgsolve::param(pji_mod, b$par)
    out <- mrgsolve::mrgsim(m, events = b$ev, end = input$tend*DAY, delta = 4) %>%
      as.data.frame()
    out$day <- out$time/DAY
    out
  }, ignoreNULL = FALSE)

  ## ---- KPI -------------------------------------------------------------------
  output$kpi_burden <- renderText(sprintf("%.2f", tail(sim()$LOGNTOT, 1)))
  output$kpi_cure   <- renderText(sprintf("%.1f%%", 100*tail(sim()$PCURE, 1)))
  output$kpi_prpob  <- renderText({
    d <- sim(); sprintf("%.3f", 1 - exp(-1e-8*min(10^d$LOGNTOT))) })
  output$kpi_icm    <- renderText(sprintf("%.0f / 10", max(sim()$ICMSC)))

  output$narrative <- renderUI({
    d <- sim(); b <- build()
    nadir <- min(10^d$LOGNTOT); pr <- 1 - exp(-1e-8*nadir)
    cure  <- tail(d$PCURE, 1)
    HTML(sprintf(
      "<div style='background:#F7F9FB;padding:12px;border-left:4px solid #2E6DA4'>
       <b>What the model says</b><br>
       At a nadir burden of <b>%.2e CFU</b> the probability that an rpoB mutant already exists is <b>%.1f%%</b>.
       The largest (free bone concentration / MBEC) ratio achieved was <b>%.3f</b>,
       and the final probability of cure is <b>%.1f%%</b>.<br>%s</div>",
      nadir, 100*pr, max(d$RATMAX), 100*cure,
      if (!("RIF" %in% input$drugs))
        "<span style='color:#C0392B'>Rifampicin is missing — no systemic agent passes 2% of the MBEC, so the biofilm is not eradicated.</span>"
      else if (length(setdiff(input$drugs, "RIF")) == 0)
        "<span style='color:#C0392B'>This is rifampicin monotherapy — in the fraction of patients given above an rpoB mutant is selected and treatment fails.</span>"
      else if (input$surg == "none")
        "<span style='color:#C0392B'>There is no debridement — at N ~ 1e10 the probability that a mutant exists is effectively 1, and rifampicin is already powerless at the outset.</span>"
      else "<span style='color:#2E7D32'>Debridement plus rifampicin — the only combination this model regards as viable.</span>"))
  })

  ## ---- 1. Overview -----------------------------------------------------------
  output$p_overview <- renderPlot({
    d <- sim()
    dd <- d %>% select(day, `Burden log10 CFU` = LOGNTOT, `CRP (mg/L)` = CRP,
                       `Synovial WBC (/uL)` = SYNWBC, `Bone mass (% of baseline)` = BONEV,
                       `Loosening index` = LOOSEN, `Probability of cure` = PCURE) %>%
      pivot_longer(-day)
    ggplot(dd, aes(day, value)) +
      geom_line(linewidth = 0.9, colour = "#2E6DA4") +
      geom_vline(xintercept = input$tdx, linetype = 2, colour = "#C0392B") +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "Day", y = NULL,
           title = "PJI course summary (red dashed = start of surgery / therapy)") + THEME
  })

  ## ---- 2. PK -----------------------------------------------------------------
  output$p_pk <- renderPlot({
    d <- sim()
    pl <- d %>% select(day, VAN = CPVAN, RIF = CPRIF, LVX = CPLVX,
                       DAP = CPDAP, LZD = CPLZD, CFZ = CPCFZ) %>%
      pivot_longer(-day, names_to = "drug", values_to = "conc") %>% mutate(site = "Plasma (total)")
    bn <- d %>% select(day, VAN = VANB, RIF = RIFB, LVX = LVXB,
                       DAP = DAPB, LZD = LZDB, CFZ = CFZB) %>%
      pivot_longer(-day, names_to = "drug", values_to = "conc") %>% mutate(site = "Bone (free)")
    keep <- union(input$drugs, "RIF")
    ggplot(bind_rows(pl, bn) %>% filter(drug %in% keep, conc > 1e-4),
           aes(day, conc, colour = drug, linetype = site)) +
      geom_line(linewidth = 0.8) + scale_y_log10() +
      scale_colour_manual(values = DRUGCOL) +
      scale_linetype_manual(values = c("Plasma (total)" = 2, "Bone (free)" = 1)) +
      labs(x = "Day", y = "Concentration (mg/L, log)", colour = "Agent", linetype = "Site",
           title = "Total plasma concentration vs free bone concentration") + THEME
  })

  ## ---- 3. Bone concentration / MBEC ------------------------------------------
  output$p_ratio <- renderPlot({
    d <- sim()
    r <- d %>% select(day, VAN = RATVAN, RIF = RATRIF, LVX = RATLVX,
                      DAP = RATDAP, LZD = RATLZD, CFZ = RATCFZ) %>%
      pivot_longer(-day, names_to = "drug", values_to = "ratio") %>%
      filter(drug %in% union(input$drugs, "RIF"), ratio > 1e-6)
    ggplot(r, aes(day, ratio, colour = drug)) +
      geom_hline(yintercept = 1, linetype = 2, colour = "#C0392B", linewidth = 0.9) +
      annotate("text", x = Inf, y = 1.35, hjust = 1.05, label = "MBEC (= 1)",
               colour = "#C0392B", size = 4) +
      geom_line(linewidth = 0.9) + scale_y_log10() +
      scale_colour_manual(values = DRUGCOL) +
      labs(x = "Day", y = "C_bone,free / MBEC (log)", colour = "Agent",
           title = "Achievable free bone concentration / biofilm-eradication concentration") + THEME
  })

  output$t_ratio <- DT::renderDataTable({
    p <- as.list(mrgsolve::param(pji_mod))
    reg <- data.frame(
      Drug  = c("Vancomycin 1 g q12h", "Rifampicin 450 mg q12h (after induction)",
                "Levofloxacin 750 mg qd", "Daptomycin 8 mg/kg qd",
                "Linezolid 600 mg q12h", "Cefazolin 2 g q8h"),
      AUC24 = c(2*1000/p$CLVAN, 2*450*p$FRIF/(2*p$CLRIF0), 750*p$FLVX/p$CLLVX,
                560/p$CLDAP, 2*600*p$FLZD/p$CLLZD0, 3*2000/p$CLCFZ),
      f_bone= c(p$PENVAN*p$FUVAN, p$PENRIF*p$FURIF, p$PENLVX*p$FULVX,
                p$PENDAP*p$FUDAP, p$PENLZD*p$FULZD, p$PENCFZ*p$FUCFZ),
      MIC   = c(p$MICVAN, p$MICRIF, p$MICLVX, p$MICDAP, p$MICLZD, p$MICCFZ),
      MBEC  = c(p$MBCVAN, p$MBCRIF, p$MBCLVX, p$MBCDAP, p$MBCLZD, p$MBCCFZ))
    reg$Cbone_free <- round(reg$AUC24/24*reg$f_bone, 3)
    reg$`Ratio / MIC`  <- round(reg$Cbone_free/reg$MIC, 1)
    reg$`Ratio / MBEC` <- round(reg$Cbone_free/reg$MBEC, 4)
    reg$AUC24 <- round(reg$AUC24, 0); reg$f_bone <- round(reg$f_bone, 3)
    DT::datatable(reg, rownames = FALSE, options = list(dom = "t", pageLength = 10)) %>%
      DT::formatStyle("Ratio / MBEC",
        backgroundColor = DT::styleInterval(c(0.02, 0.10), c("#FDECEC", "#FDF3E7", "#EAF6F0")))
  })

  ## ---- 4. Subpopulations -----------------------------------------------------
  output$p_pops <- renderPlot({
    d <- sim()
    dd <- d %>% transmute(day,
      `Planktonic`    = NP,  `Biofilm`              = NB,
      `Persister`     = NPER,`SCV`                  = NSCV,
      `Intracellular` = NIC,
      `rpoB resistant` = RP + RB,      `gyrA resistant` = QP + QB) %>%
      pivot_longer(-day, names_to = "pop", values_to = "n") %>%
      mutate(n = pmax(n, 1e-3))
    g <- ggplot(dd, aes(day, n, colour = pop)) +
      geom_line(linewidth = 0.85) +
      geom_vline(xintercept = input$tdx, linetype = 2, colour = "#555") +
      scale_colour_manual(values = POPCOL) +
      labs(x = "Day", y = "CFU", colour = NULL,
           title = "Bacterial subpopulation dynamics (dashed = surgery)") + THEME
    if (input$logscale) g <- g + scale_y_log10(limits = c(1e-2, 1e11))
    g
  })

  ## ---- 5. Biofilm ------------------------------------------------------------
  output$p_biofilm <- renderPlot({
    d <- sim()
    dd <- d %>% select(day, `EPS matrix (0-1)` = EPS, `EPS maturity` = EPSMAT,
                       `agr signal AIP` = AIP, `Sequestrum SEQ` = SEQ,
                       `MDSC` = MDSC, `Biofilm burden log10` = LOGNB) %>%
      pivot_longer(-day)
    ggplot(dd, aes(day, value)) + geom_line(linewidth = 0.9, colour = "#4E9E7F") +
      geom_vline(xintercept = input$tdx, linetype = 2, colour = "#C0392B") +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "Day", y = NULL, title = "Biofilm matrix maturation and its by-products") + THEME
  })

  ## ---- 6. Local elution ------------------------------------------------------
  output$p_local <- renderPlot({
    d <- sim(); p <- as.list(mrgsolve::param(pji_mod))
    ggplot(d, aes(day, pmax(LOCJ, 1e-3))) +
      geom_line(linewidth = 1, colour = DRUGCOL[["LOCAL"]]) +
      geom_hline(yintercept = p$MBCVAN, linetype = 2, colour = "#C0392B") +
      annotate("text", x = Inf, y = p$MBCVAN*1.6, hjust = 1.05,
               label = "Vancomycin MBEC 512 mg/L", colour = "#C0392B") +
      geom_hline(yintercept = p$MICVAN, linetype = 3, colour = "#2E7D32") +
      scale_y_log10() +
      labs(x = "Day", y = "Synovial fluid antibiotic concentration (mg/L, log)",
           title = "Local elution from an antibiotic-loaded cement spacer") + THEME
  })

  ## ---- 7. Resistance --------------------------------------------------------
  output$p_res <- renderPlot({
    d <- sim()
    dd <- d %>% transmute(day,
      `Susceptible total` = pmax(NP + NB + NPER + NSCV + NIC, 1e-3),
      `rpoB mutant`   = pmax(RP + RB, 1e-3),
      `gyrA mutant`   = pmax(QP + QB, 1e-3)) %>%
      pivot_longer(-day, names_to = "clone", values_to = "n")
    ggplot(dd, aes(day, n, colour = clone)) +
      geom_line(linewidth = 0.95) +
      geom_hline(yintercept = 1, linetype = 3, colour = "#888") +
      geom_vline(xintercept = input$tdx, linetype = 2, colour = "#555") +
      scale_y_log10(limits = c(1e-2, 1e11)) +
      scale_colour_manual(values = c("Susceptible total" = "#2E6DA4",
                                     "rpoB mutant" = "#C0392B", "gyrA mutant" = "#7F8C8D")) +
      labs(x = "Day", y = "CFU (log)", colour = NULL,
           title = "Clone dynamics — resistance is selected, not created") + THEME
  })

  output$p_mutsupply <- renderPlot({
    d <- sim()
    nad <- min(10^d$LOGNTOT); pk <- max(10^d$LOGNTOT)
    N <- 10^seq(3, 11, by = 0.05)
    df <- data.frame(N = N, P = 1 - exp(-1e-8*N))
    ggplot(df, aes(N, P)) + geom_line(linewidth = 1, colour = "#B08040") +
      geom_vline(xintercept = pk,  linetype = 2, colour = "#C0392B") +
      geom_vline(xintercept = max(nad, 1), linetype = 2, colour = "#2E7D32") +
      annotate("text", x = pk, y = 0.15, hjust = 1.1,
               label = sprintf("Peak burden %.1e", pk), colour = "#C0392B") +
      annotate("text", x = max(nad, 1), y = 0.85, hjust = -0.05,
               label = sprintf("Nadir burden %.1e", nad), colour = "#2E7D32") +
      scale_x_log10() +
      labs(x = "Burden N (CFU, log)", y = "P(pre-existing rpoB mutant) = 1 − e^(−μN)",
           title = "μ = 1e−8 — what surgery actually changes") + THEME
  })

  ## ---- 8. Immunity ----------------------------------------------------------
  output$p_immune <- renderPlot({
    d <- sim()
    dd <- d %>% select(day, `Neutrophils (/uL)` = PMN, `Monocytes (/uL)` = MONO,
                       `MDSC` = MDSC, `Macrophages` = MAC,
                       `IL-1beta` = IL1B, `TNF-alpha` = TNFA,
                       `IL-6` = IL6, `IL-10` = IL10) %>% pivot_longer(-day)
    ggplot(dd, aes(day, value)) + geom_line(linewidth = 0.85, colour = "#5E86BF") +
      geom_vline(xintercept = input$tdx, linetype = 2, colour = "#C0392B") +
      facet_wrap(~name, scales = "free_y", ncol = 4) +
      labs(x = "Day", y = NULL, title = "Innate immune cells and cytokines") + THEME
  })

  ## ---- 9. Bone -------------------------------------------------------------
  output$p_bone <- renderPlot({
    d <- sim()
    dd <- d %>% select(day, `RANKL` = RANKL, `OPG` = OPG,
                       `Osteoblasts` = OBL, `Osteoclasts` = OCL,
                       `Peri-implant bone mass (%)` = BONEV, `Loosening index (0-100)` = LOOSEN) %>%
      pivot_longer(-day)
    ggplot(dd, aes(day, value)) + geom_line(linewidth = 0.9, colour = "#A98A5F") +
      geom_vline(xintercept = input$tdx, linetype = 2, colour = "#C0392B") +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "Day", y = NULL, title = "The RANKL/OPG axis and peri-implant osteolysis") + THEME
  })

  ## ---- 10. Biomarkers -------------------------------------------------------
  output$p_biomarker <- renderPlot({
    d <- sim()
    dd <- d %>% select(day, `CRP (mg/L)` = CRP, `ESR (mm/h)` = ESR,
                       `Synovial WBC (/uL)` = SYNWBC, `PMN (%)` = PMNPCT,
                       `alpha-defensin (S/CO)` = ADEF, `ICM score` = ICMSC) %>%
      pivot_longer(-day)
    thr <- data.frame(name = c("CRP (mg/L)", "ESR (mm/h)", "Synovial WBC (/uL)",
                               "PMN (%)", "alpha-defensin (S/CO)", "ICM score"),
                      y = c(10, 30, 3000, 80, 1.0, 6))
    ggplot(dd, aes(day, value)) + geom_line(linewidth = 0.9, colour = "#7592AE") +
      geom_hline(data = thr, aes(yintercept = y), linetype = 2, colour = "#C0392B") +
      geom_vline(xintercept = input$tdx, linetype = 3, colour = "#555") +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "Day", y = NULL, title = "Diagnostic biomarkers and the 2018 ICM score") + THEME
  })

  ## ---- 11. Scenario comparison ----------------------------------------------
  allout <- eventReactive(input$runall, {
    withProgress(message = "Running the 21 scenarios...", value = 0, {
      res <- lapply(seq_along(scen), function(i) {
        incProgress(1/length(scen), detail = names(scen)[i])
        d <- run_scen(pji_mod, scen[[i]], end = input$tend*DAY)
        d$scenario <- names(scen)[i]; d
      })
      bind_rows(res)
    })
  })

  output$t_scen <- DT::renderDataTable({
    s <- summarise_scen(allout())
    DT::datatable(s, rownames = FALSE,
                  options = list(pageLength = 21, scrollX = TRUE)) %>%
      DT::formatStyle("P_cure",
        backgroundColor = DT::styleInterval(c(0.05, 0.80), c("#FDECEC", "#FDF3E7", "#EAF6F0")))
  })

  output$p_scen <- renderPlot({
    d <- allout()
    ggplot(d, aes(time/DAY, LOGNTOT, colour = scenario)) +
      geom_line(linewidth = 0.7, show.legend = FALSE) +
      facet_wrap(~scenario, ncol = 4) +
      labs(x = "Day", y = "Total burden (log10 CFU)",
           title = "Burden course of the 21 standard scenarios") + THEME +
      theme(strip.text = element_text(size = 8))
  })

  ## ---- 12. Toxicity ---------------------------------------------------------
  output$p_tox <- renderPlot({
    d <- sim()
    dd <- d %>% select(day, `Serum creatinine (mg/dL)` = SCR,
                       `VAN AUC24 (mg*h/L)` = AUCF,
                       `Platelets (1e9/L)` = PLT, `ALT (U/L)` = ALT) %>% pivot_longer(-day)
    thr <- data.frame(name = c("Serum creatinine (mg/dL)", "VAN AUC24 (mg*h/L)",
                               "Platelets (1e9/L)", "ALT (U/L)"),
                      y = c(1.35, 600, 100, 75))
    ggplot(dd, aes(day, value)) + geom_line(linewidth = 0.9, colour = "#CB8B62") +
      geom_hline(data = thr, aes(yintercept = y), linetype = 2, colour = "#C0392B") +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "Day", y = NULL, title = "Drug toxicity measures (dashed = warning threshold)") + THEME
  })

  ## ---- 13. Parameters -------------------------------------------------------
  output$t_par <- DT::renderDataTable({
    p <- as.list(mrgsolve::param(pji_mod))
    DT::datatable(data.frame(parameter = names(p), value = unlist(p)),
                  rownames = FALSE, options = list(pageLength = 25, scrollX = TRUE))
  })

  output$reflink <- renderUI(HTML(
    "<p>For the provenance of the parameters and the verification anchors see <code>pji_references_en.md</code> (68 papers, every PMID checked) and
     the VERIFICATION section in the header of <code>pji_mrgsolve_model_en.R</code>.
     The mechanistic map is <code>pji_qsp_model_en.svg</code>.</p>
     <p style='color:#C0392B'><b>Disclaimer</b>: this is a semi-quantitative model for education and research.
     Do not use it for clinical decisions, prescribing, or regulatory submission.</p>"))
}

shinyApp(ui, server)
