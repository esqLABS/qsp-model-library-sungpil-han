## ============================================================================
##  Acne vulgaris (ACN) — QSP Interactive Dashboard
##  ============================================================================
##  Front-end for acn_mrgsolve_model.R.
##
##  Run:
##      shiny::runApp("acn_shiny_app.R")
##  or from this directory:
##      R -e 'shiny::runApp(".", launch.browser = TRUE)'
##
##  Ten tabs:
##      1. Patient profile           — phenotype, androgen/metabolic set-point, exposome
##      2. Regimen design            — build a regimen (topical / systemic / hormonal)
##      3. PK                — every drug concentration on one time axis
##      4. Sebaceous gland·androgen  — SGM, SER, AR signal, SHBG/FAI, IGF-1, mTORC1
##      5. Microbiome·antibiotic resistance — C. acnes planktonic/biofilm and the resistant
##                             fraction; the stewardship argument in one plot
##      6. Inflammatory cascade      — TLR2 -> IL-1beta -> IL-8 -> neutrophils -> MMP
##      7. Lesion counts·endpoints   — the lesion-transit chain and IGA
##      8. Scenario comparison       — the 17 prebuilt scenarios side by side
##      9. Safety                    — TG, ALT, K+, mucocutaneous, cumulative mg/kg
##     10. Scarring·pigmentation     — PIH and atrophic scarring, and what prevents them
##
##  DISCLAIMER: educational / research tool. Not for clinical use.
## ============================================================================

suppressPackageStartupMessages({
  library(shiny)
  library(mrgsolve)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(DT)
})

source("acn_mrgsolve_model.R", local = TRUE, chdir = TRUE)

MOD <- ACN_build()

## ---------------------------------------------------------------------------
##  Plot furniture
## ---------------------------------------------------------------------------
theme_acn <- function() {
  theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          strip.text      = element_text(face = "bold"),
          legend.position = "bottom",
          legend.title    = element_blank(),
          plot.title      = element_text(face = "bold"))
}

PAL <- c("#C0392B", "#2874A6", "#1E8449", "#B9770E", "#7D3C98",
         "#117A65", "#CA6F1E", "#5D6D7E", "#A93226", "#1F618D")

long_plot <- function(d, vars, labels = vars, ylab = "", title = "",
                      hline = NULL, logy = FALSE) {
  dd <- d[, c("week", vars), drop = FALSE]
  names(dd) <- c("week", labels)
  dd <- tidyr::pivot_longer(dd, -week, names_to = "v", values_to = "y")
  dd$v <- factor(dd$v, levels = labels)
  p <- ggplot(dd, aes(week, y, colour = v)) +
    geom_line(linewidth = 0.9) +
    scale_colour_manual(values = rep(PAL, 4)) +
    labs(x = "Weeks", y = ylab, title = title) +
    theme_acn()
  if (!is.null(hline)) p <- p + geom_hline(yintercept = hline, linetype = 2,
                                           colour = "grey40")
  if (logy) p <- p + scale_y_log10()
  p
}

## ---------------------------------------------------------------------------
##  UI
## ---------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Acne vulgaris QSP dashboard"),
  tags$p(style = "color:#7B241C;font-weight:bold;",
         "This is an educational and research model. Do not use it for clinical decisions. / ",
         "Educational model. Not for clinical use."),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      h4("1. Patient phenotype"),
      selectInput("pheno", "Phenotype",
                  choices = c("mild comedonal"              = "mild_comedonal",
                              "moderate"                    = "moderate",
                              "moderate male"               = "moderate_male",
                              "severe nodular"              = "severe_nodular",
                              "adult female"                = "adult_female",
                              "with PCOS"                   = "pcos",
                              "skin of colour"              = "skin_of_colour"),
                  selected = "moderate"),
      sliderInput("SEVX",   "Constitutional sebum drive SEVX", 1.0, 2.4, 1.55, 0.05),
      sliderInput("NODPROP","Nodular propensity NODPROP", 0.1, 3.0, 1.00, 0.1),
      sliderInput("SKINPIG","Pigmentation propensity (Fitzpatrick)", 0.3, 2.4, 1.0, 0.1),

      h4("2. Exposome"),
      sliderInput("GLYLOAD","Glycaemic load (0=low GI, 1=high GI)", 0, 1, 0.5, 0.05),
      sliderInput("DAIRY",  "Dairy / whey intake",         0, 1, 0.3, 0.05),
      sliderInput("STRESS", "Psychological stress",        0, 1, 0.3, 0.05),
      sliderInput("UVX",    "UV / oxidative exposure",     0, 1, 0.2, 0.05),
      sliderInput("PHYLIA", "C. acnes IA1 predominance",   0, 1, 0.75, 0.05),

      h4("3. Topical therapy"),
      selectInput("ret", "Topical retinoid",
                  choices = c("none" = "none", names(ACN_RETPOT)),
                  selected = "adapalene_0.1"),
      checkboxInput("bpo",   "Benzoyl peroxide 2.5% QD", FALSE),
      checkboxInput("cli",   "Clindamycin 1% BID",       FALSE),
      checkboxInput("aze",   "Azelaic acid 15% BID",     FALSE),
      checkboxInput("clasco","Clascoterone 1% BID",      FALSE),
      checkboxInput("dap",   "Dapsone 7.5% QD",          FALSE),
      sliderInput("ADHERE",  "Adherence", 0.2, 1.0, 1.0, 0.05),

      h4("4. Systemic therapy"),
      selectInput("tet", "Tetracycline class",
                  choices = c("none"                          = "none",
                              "doxycycline 100 mg QD"         = "doxy100",
                              "doxycycline 40 mg MR (sub-antimicrobial)" = "sub40",
                              "minocycline 100 mg QD"         = "mino100",
                              "sarecycline 1.5 mg/kg QD"      = "sare"),
                  selected = "none"),
      numericInput("tet_wk", "Antibiotic duration (weeks)", 12, 0, 52, 1),
      sliderInput("iso",  "Isotretinoin (mg/kg/day)", 0, 1.2, 0, 0.05),
      numericInput("iso_wk", "Isotretinoin duration (weeks)", 24, 0, 60, 1),
      checkboxInput("food", "Taken with a high-fat meal", TRUE),
      checkboxInput("lidose", "Lidose / micronised formulation", FALSE),
      sliderInput("spiro","Spironolactone (mg/day)", 0, 200, 0, 25),
      sliderInput("coc",  "COC ethinylestradiol (µg)", 0, 35, 0, 5),

      h4("5. Simulation"),
      numericInput("weeks", "Simulation duration (weeks)", 52, 4, 160, 4),
      actionButton("go", "Run", class = "btn-primary btn-block")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",

        tabPanel("1. Patient profile",
          h4("Baseline steady state (untreated)"),
          p("Every simulation starts from this patient's own steady state. ",
            "The lesion count is not an input but a value that emerges from the ",
            "androgen, metabolic and exposome settings."),
          DT::dataTableOutput("tbl_base"),
          hr(),
          h4("The four pathogenic pillars"),
          plotOutput("p_pillars", height = "320px"),
          hr(),
          verbatimTextOutput("txt_pheno")),

        tabPanel("2. Regimen design",
          h4("The selected regimen"),
          verbatimTextOutput("txt_regimen"),
          hr(),
          h4("Guideline sanity checks"),
          htmlOutput("txt_flags")),

        tabPanel("3. PK",
          plotOutput("p_pk_sys", height = "300px"),
          plotOutput("p_pk_top", height = "300px"),
          plotOutput("p_pk_eff", height = "280px")),

        tabPanel("4. Sebaceous gland·androgen",
          plotOutput("p_seb",  height = "300px"),
          plotOutput("p_horm", height = "300px"),
          plotOutput("p_met",  height = "260px")),

        tabPanel("5. Microbiome·resistance",
          p(strong("This tab is the antibiotic stewardship argument."),
            " Under antibiotic monotherapy the resistant fraction (RESF) rises, the ",
            "bacterial load comes back up and lesions rebound. Adding BPO clears the ",
            "resistant fraction and the rebound disappears."),
          plotOutput("p_micro", height = "320px"),
          plotOutput("p_res",   height = "300px")),

        tabPanel("6. Inflammatory cascade",
          plotOutput("p_infl", height = "320px"),
          plotOutput("p_ker",  height = "300px")),

        tabPanel("7. Lesion counts·endpoints",
          plotOutput("p_les",  height = "340px"),
          plotOutput("p_iga",  height = "260px"),
          hr(),
          h4("Endpoint summary"),
          DT::dataTableOutput("tbl_end")),

        tabPanel("8. Scenario comparison",
          checkboxGroupInput("scn", "Select scenarios", inline = TRUE,
            choices = setNames(1:17, paste0("S", 1:17)),
            selected = c(2, 3, 5, 8, 9, 11)),
          numericInput("scn_wk", "Assessment time point (weeks)", 12, 1, 80, 1),
          actionButton("go_scn", "Run scenarios", class = "btn-primary"),
          hr(),
          DT::dataTableOutput("tbl_scn"),
          plotOutput("p_scn", height = "360px")),

        tabPanel("9. Safety",
          plotOutput("p_saf", height = "320px"),
          plotOutput("p_cum", height = "260px"),
          hr(),
          htmlOutput("txt_saf")),

        tabPanel("10. Scarring·pigmentation",
          plotOutput("p_scar", height = "320px"),
          hr(),
          p("Atrophic scarring accumulates as the time integral of nodule burden x MMP ",
            "activity and cannot be reversed within the model. The only way to reduce ",
            "scarring is to shorten the time nodules exist — that is, early adequate treatment."),
          DT::dataTableOutput("tbl_scar")),

        tabPanel("11. References·help",
          h4("Model structure"),
          verbatimTextOutput("txt_help"))
      )
    )
  )
)

## ---------------------------------------------------------------------------
##  Server
## ---------------------------------------------------------------------------
server <- function(input, output, session) {

  ## ---- assemble the parameter overrides and the event table ---------------
  regimen <- reactive({
    ev  <- NULL
    lab <- character(0)
    prm <- list(ADHERE = input$ADHERE)
    wk  <- input$weeks

    if (input$ret != "none") {
      ev <- c(ev, acn_retinoid(input$ret, days = wk * 7))
      lab <- c(lab, paste("Topical retinoid:", input$ret))
    }
    if (isTRUE(input$bpo)) {
      ev <- c(ev, acn_bpo(2.5, days = wk * 7)); lab <- c(lab, "BPO 2.5% QD")
    }
    if (isTRUE(input$cli)) {
      ev <- c(ev, acn_clindamycin(1.0, days = wk * 7)); lab <- c(lab, "Clindamycin 1% BID")
    }
    if (isTRUE(input$aze)) {
      ev <- c(ev, acn_azelaic(15, days = wk * 7)); lab <- c(lab, "Azelaic acid 15% BID")
    }
    if (isTRUE(input$clasco)) {
      ev <- c(ev, acn_clascoterone(1.0, days = wk * 7)); lab <- c(lab, "Clascoterone 1% BID")
    }
    if (isTRUE(input$dap)) {
      ev <- c(ev, acn_dapsone(7.5, days = wk * 7)); lab <- c(lab, "Dapsone 7.5% QD")
    }

    if (input$tet != "none" && input$tet_wk > 0) {
      spec <- switch(input$tet,
        doxy100 = list(dose = 100, p = list()),
        sub40   = list(dose = 40,  p = list(SUBANTI = 1)),
        mino100 = list(dose = 100, p = list(CLTET = 1.6, VTET = 80)),
        sare    = list(dose = round(1.5 * 62), p = list(NARROW = 1, CLTET = 1.75, VTET = 55)))
      ev  <- c(ev, acn_tetracycline(spec$dose, days = input$tet_wk * 7))
      prm <- utils::modifyList(prm, spec$p)
      lab <- c(lab, sprintf("%s x %d weeks", input$tet, input$tet_wk))
    }

    if (input$iso > 0 && input$iso_wk > 0) {
      ev  <- c(ev, acn_isotretinoin(input$iso, wt = 62, days = input$iso_wk * 7))
      prm <- utils::modifyList(prm, list(FOOD = as.numeric(isTRUE(input$food)),
                                         LIDOSE = as.numeric(isTRUE(input$lidose))))
      lab <- c(lab, sprintf("Isotretinoin %.2f mg/kg/day x %d weeks (cumulative %.0f mg/kg)",
                            input$iso, input$iso_wk, input$iso * input$iso_wk * 7))
    }
    if (input$spiro > 0) {
      ev <- c(ev, acn_spironolactone(input$spiro, days = wk * 7))
      lab <- c(lab, sprintf("Spironolactone %d mg/day", input$spiro))
    }
    if (input$coc > 0) {
      ev <- c(ev, acn_coc(input$coc, cycles = ceiling(wk / 4)))
      lab <- c(lab, sprintf("COC EE %d µg", input$coc))
    }
    if (!length(lab)) lab <- "No treatment (natural history)"
    list(ev = ev, param = prm, label = lab)
  })

  base <- eventReactive(input$go, {
    extra <- list(SEVX = input$SEVX, NODPROP = input$NODPROP,
                  SKINPIG = input$SKINPIG, GLYLOAD = input$GLYLOAD,
                  DAIRY = input$DAIRY, STRESS = input$STRESS,
                  UVX = input$UVX, PHYLIA = input$PHYLIA)
    acn_baseline(MOD, input$pheno, extra = extra)
  }, ignoreNULL = FALSE)

  sim <- eventReactive(input$go, {
    b <- base()
    r <- regimen()
    p <- utils::modifyList(b$param,
                           list(PBOMAX = as.numeric(mrgsolve::param(MOD)$PBOMAX)))
    p <- utils::modifyList(p, r$param)
    m <- mrgsolve::update(MOD, param = p, init = b$init)
    d <- if (is.null(r$ev)) {
      mrgsolve::mrgsim_df(m, end = input$weeks * 168, delta = 12, hmax = 4)
    } else {
      mrgsolve::mrgsim_df(m, events = r$ev, end = input$weeks * 168,
                          delta = 12, hmax = 4)
    }
    d$week <- d$time / 168
    d
  }, ignoreNULL = FALSE)

  ## ---- 1. patient profile -------------------------------------------------
  output$tbl_base <- DT::renderDataTable({
    s <- base()$summary
    data.frame(
      Metric = c("Inflammatory lesion count", "Non-inflammatory lesion count", "Total lesion count", "IGA (0-4)",
               "Microcomedone reservoir MC", "Sebum excretion rate SER (µg/cm²/min)",
               "C. acnes load (relative)", "Resistant fraction", "Keratinisation index KER",
               "AR signal ARS", "Free androgen index FAI", "IGF-1 (ng/mL)",
               "PIH index", "Scar index"),
      Value = round(c(s$INFLAM, s$NONINF, s$TOTLES, s$IGA, s$MC, s$SERO,
                   s$CACNT, s$RESF, s$KER, s$ARS, s$FAIO, s$IGF1,
                   s$PIH, s$SCAR), 2))
  }, options = list(dom = "t", pageLength = 20), rownames = FALSE)

  output$p_pillars <- renderPlot({
    s <- base()$summary
    d <- data.frame(
      pillar = factor(c("① Seborrhoea\n(SER/ref)", "② Hyperkeratinisation\n(KER)",
                        "③ C. acnes\n(load/ref)", "④ Inflammation\n(IL-1β)"),
                      levels = c("① Seborrhoea\n(SER/ref)", "② Hyperkeratinisation\n(KER)",
                                 "③ C. acnes\n(load/ref)", "④ Inflammation\n(IL-1β)")),
      value  = c(s$SERO / 0.9, s$KER, s$CACNT / 1.5, s$IL1B))
    ggplot(d, aes(pillar, value, fill = pillar)) +
      geom_col(width = 0.6) +
      geom_hline(yintercept = 1, linetype = 2) +
      scale_fill_manual(values = PAL[c(4, 3, 6, 1)], guide = "none") +
      labs(x = NULL, y = "Fold of healthy reference",
           title = "The four pathogenic pillars — 1.0 = a healthy pilosebaceous unit") +
      theme_acn()
  })

  output$txt_pheno <- renderPrint({
    p <- base()$param
    cat("Phenotype parameters\n")
    str(p[order(names(p))], give.attr = FALSE)
  })

  ## ---- 2. regimen ---------------------------------------------------------
  output$txt_regimen <- renderPrint({ cat(paste(regimen()$label, collapse = "\n")) })

  output$txt_flags <- renderUI({
    r  <- regimen(); msg <- character(0)
    abx_top <- isTRUE(input$cli)
    abx_sys <- input$tet != "none" && input$tet_wk > 0 && input$tet != "sub40"
    has_bpo <- isTRUE(input$bpo)
    has_ret <- input$ret != "none"

    if ((abx_top || abx_sys) && !has_bpo)
      msg <- c(msg, paste("⚠️ An antibiotic is being used without BPO. The resistant",
                          "fraction rises and the effect is lost (see tab 5). The guidelines",
                          "require BPO or a retinoid alongside any antibiotic."))
    if (abx_sys && input$tet_wk > 16)
      msg <- c(msg, "⚠️ Systemic antibiotic beyond 16 weeks. The guideline recommendation is within 3–4 months.")
    if ((abx_top || abx_sys) && !has_ret)
      msg <- c(msg, paste("⚠️ Anti-inflammatory therapy is being used without a retinoid.",
                          "The microcomedone reservoir remains, so relapse follows withdrawal."))
    if (input$iso > 0 && input$iso * input$iso_wk * 7 < 120)
      msg <- c(msg, sprintf(paste("⚠️ Projected cumulative isotretinoin dose %.0f mg/kg —",
                                  "below the 120–150 mg/kg target, so relapse risk is high."),
                            input$iso * input$iso_wk * 7))
    if (input$iso > 0 && abx_sys)
      msg <- c(msg, "⛔ Tetracycline + isotretinoin is contraindicated because of the risk of idiopathic intracranial hypertension.")
    if (input$iso > 0)
      msg <- c(msg, "⛔ Women of childbearing potential: teratogenicity — two forms of contraception and a pregnancy testing programme are mandatory.")
    if (!length(msg)) msg <- "✅ No particular warnings."
    HTML(paste0("<ul>", paste0("<li>", msg, "</li>", collapse = ""), "</ul>"))
  })

  ## ---- 3. PK --------------------------------------------------------------
  output$p_pk_sys <- renderPlot({
    long_plot(sim(), c("CTETO", "CISOO", "COXOO", "CSPIO"),
              c("Tetracycline", "Isotretinoin", "4-oxo metabolite", "Canrenone"),
              "Concentration (mg/L)", "Systemic drug concentrations")
  })
  output$p_pk_top <- renderPlot({
    long_plot(sim(), c("BPOS", "RETS", "CLIS", "AZES", "CLAS", "DAPS"),
              c("BPO", "Retinoid", "Clindamycin", "Azelaic acid",
                "Clascoterone", "Dapsone"),
              "Relative intrafollicular concentration", "Topical drug follicular reservoir")
  })
  output$p_pk_eff <- renderPlot({
    long_plot(sim(), c("ISOEFFT", "AIEFFT", "CEEO"),
              c("Isotretinoin effect (0-1)", "Tetracycline anti-inflammatory effect (0-1)",
                "EE concentration (pg/mL)"),
              "", "Pharmacodynamic effect measures")
  })

  ## ---- 4. sebaceous / hormonal -------------------------------------------
  output$p_seb <- renderPlot({
    long_plot(sim(), c("SER", "SGM", "LIP", "DURAB"),
              c("Sebum excretion rate SER", "Sebaceous gland mass SGM", "Lipogenesis drive LIP",
                "Durable gland shrinkage DURAB"),
              "", "The sebaceous gland — the target of isotretinoin", hline = 1)
  })
  output$p_horm <- renderPlot({
    long_plot(sim(), c("TT", "SHBG", "FAIO"),
              c("Total testosterone (ng/dL)", "SHBG (nmol/L)", "Free androgen index"),
              "", "The androgen axis — target of COC / spironolactone / clascoterone")
  })
  output$p_met <- renderPlot({
    long_plot(sim(), c("IGF1", "INS", "ARS"),
              c("IGF-1 (ng/mL)", "Insulin (µU/mL)", "AR signal (relative)"),
              "", "Metabolic and receptor signalling")
  })

  ## ---- 5. microbiology ----------------------------------------------------
  output$p_micro <- renderPlot({
    long_plot(sim(), c("CAP", "CAB", "CACNT", "FFA", "PORP"),
              c("Planktonic C. acnes", "Biofilm", "Total load",
                "Free fatty acids", "Porphyrins"),
              "Relative (healthy reference = 1)", "C. acnes ecology")
  })
  output$p_res <- renderPlot({
    d <- sim()
    ggplot(d, aes(week)) +
      geom_line(aes(y = RESF, colour = "Antibiotic-resistant fraction RESF"), linewidth = 1.1) +
      geom_line(aes(y = CACNT / max(1e-6, max(d$CACNT)),
                    colour = "Total bacterial load (normalised)"), linewidth = 0.9) +
      geom_line(aes(y = INFLAM / max(1e-6, max(d$INFLAM)),
                    colour = "Inflammatory lesions (normalised)"), linewidth = 0.9) +
      scale_colour_manual(values = PAL[c(1, 3, 2)]) +
      labs(x = "Weeks", y = "0–1", title = "Resistance selection pressure and lesion rebound") +
      theme_acn()
  })

  ## ---- 6. inflammation ----------------------------------------------------
  output$p_infl <- renderPlot({
    long_plot(sim(), c("TLR", "IL1B", "IL8", "TNF", "IL17", "NEU", "MMP"),
              c("TLR2 signal", "IL-1β", "IL-8", "TNF-α", "IL-17A",
                "Neutrophils", "MMP"),
              "Relative", "Innate → effector inflammatory cascade", hline = 1)
  })
  output$p_ker <- renderPlot({
    long_plot(sim(), c("KER", "IL1A", "LA", "SQOX"),
              c("Keratinisation index KER", "IL-1α", "Linoleic acid (dilution)",
                "Squalene peroxide"),
              "Relative", "Infundibular hyperkeratinisation — the retinoid target", hline = 1)
  })

  ## ---- 7. lesions ---------------------------------------------------------
  output$p_les <- renderPlot({
    long_plot(sim(), c("MC", "CC", "OC", "PAP", "PUS", "NOD"),
              c("Microcomedone (MC, clinically undetectable)", "Closed comedone", "Open comedone",
                "Papule", "Pustule", "Nodule"),
              "Lesion count", "The lesion-transit chain — MC is the reservoir of relapse")
  })
  output$p_iga <- renderPlot({
    long_plot(sim(), c("INFLAM", "NONINF", "IGA"),
              c("Inflammatory lesion count", "Non-inflammatory lesion count", "IGA (0-4)"),
              "", "Clinical endpoints")
  })
  output$tbl_end <- DT::renderDataTable({
    d <- sim(); b <- d[1, ]; e <- d[nrow(d), ]
    data.frame(
      Metric = c("Inflammatory lesions", "Non-inflammatory lesions", "Total lesions", "IGA",
               "Sebum excretion rate", "C. acnes total load", "Resistant fraction",
               "PIH index", "Scar index"),
      Baseline = round(c(b$INFLAM, b$NONINF, b$TOTLES, b$IGA, b$SERO,
                     b$CACNT, b$RESF, b$PIH, b$SCAR), 2),
      Final = round(c(e$INFLAM, e$NONINF, e$TOTLES, e$IGA, e$SERO,
                     e$CACNT, e$RESF, e$PIH, e$SCAR), 2),
      Change_pct = round(100 * (c(e$INFLAM, e$NONINF, e$TOTLES, e$IGA, e$SERO,
                                  e$CACNT, e$RESF, e$PIH, e$SCAR) -
                                c(b$INFLAM, b$NONINF, b$TOTLES, b$IGA, b$SERO,
                                  b$CACNT, b$RESF, b$PIH, b$SCAR)) /
                        pmax(1e-6, c(b$INFLAM, b$NONINF, b$TOTLES, b$IGA, b$SERO,
                                     b$CACNT, b$RESF, b$PIH, b$SCAR)), 1))
  }, options = list(dom = "t", pageLength = 12), rownames = FALSE)

  ## ---- 8. scenarios -------------------------------------------------------
  scn_sim <- eventReactive(input$go_scn, {
    ids <- as.integer(input$scn)
    if (!length(ids)) return(NULL)
    ACN_simulate(MOD, which = ids, delta = 24)
  })
  output$tbl_scn <- DT::renderDataTable({
    s <- scn_sim(); if (is.null(s)) return(NULL)
    ACN_summary(s, at_week = input$scn_wk)
  }, options = list(scrollX = TRUE, pageLength = 20), rownames = FALSE)
  output$p_scn <- renderPlot({
    s <- scn_sim(); if (is.null(s)) return(NULL)
    ggplot(s, aes(week, INFLAM, colour = factor(sid))) +
      geom_line(linewidth = 0.9) +
      scale_colour_manual(values = rep(PAL, 3),
                          labels = unique(paste0("S", s$sid, " ", s$scenario))) +
      labs(x = "Weeks", y = "Inflammatory lesion count",
           title = "Inflammatory lesion course by scenario") +
      theme_acn() + theme(legend.text = element_text(size = 8))
  })

  ## ---- 9. safety ----------------------------------------------------------
  output$p_saf <- renderPlot({
    long_plot(sim(), c("TG", "ALT", "KSER", "MUCO"),
              c("Triglycerides (mg/dL)", "ALT (U/L)", "Serum potassium (mEq/L)",
                "Mucocutaneous toxicity (0-10)"),
              "", "Safety tracking")
  })
  output$p_cum <- renderPlot({
    d <- sim()
    ggplot(d, aes(week, CUMISO)) +
      geom_line(linewidth = 1.1, colour = PAL[1]) +
      geom_hline(yintercept = c(120, 150), linetype = 2, colour = "grey35") +
      annotate("text", x = max(d$week) * 0.05, y = 135,
               label = "Target cumulative 120–150 mg/kg", hjust = 0, size = 3.5) +
      labs(x = "Weeks", y = "Cumulative isotretinoin (mg/kg)",
           title = "Cumulative dose — the variable that determines the relapse rate") +
      theme_acn()
  })
  output$txt_saf <- renderUI({
    d <- sim(); e <- d[nrow(d), ]
    HTML(sprintf(paste0(
      "<ul><li>Peak triglycerides <b>%.0f mg/dL</b> (baseline %.0f)</li>",
      "<li>Peak ALT <b>%.0f U/L</b></li>",
      "<li>Peak serum potassium <b>%.2f mEq/L</b></li>",
      "<li>Peak mucocutaneous toxicity score <b>%.1f / 10</b> — cheilitis occurs in virtually every case and ",
      "is used as a surrogate for confirming adherence.</li>",
      "<li>Total cumulative isotretinoin <b>%.0f mg/kg</b></li></ul>"),
      max(d$TG), d$TG[1], max(d$ALT), max(d$KSER), max(d$MUCO), max(d$CUMISO)))
  })

  ## ---- 10. scarring -------------------------------------------------------
  output$p_scar <- renderPlot({
    long_plot(sim(), c("NOD", "MMP", "SCAR", "PIH"),
              c("Nodule count", "MMP activity", "Atrophic scar index (irreversible)",
                "PIH index"),
              "", "Accumulation of scarring and pigmentation")
  })
  output$tbl_scar <- DT::renderDataTable({
    d <- sim(); b <- d[1, ]; e <- d[nrow(d), ]
    data.frame(
      Metric = c("Cumulative nodule-time (nodule-weeks)", "Final scar index",
               "Final PIH index", "Estimated PIH half-life (weeks)"),
      Value = round(c(sum(d$NOD) * (d$week[2] - d$week[1]),
                   e$SCAR, e$PIH, log(2) / (1.2e-3 * 168)), 2))
  }, options = list(dom = "t"), rownames = FALSE)

  ## ---- 11. help -----------------------------------------------------------
  output$txt_help <- renderPrint({
    cat("Acne vulgaris QSP model\n")
    cat("  - 55 ODE compartments, 213 parameters, 17 scenarios\n")
    cat("  - mechanistic map : acn_qsp_model_en.dot / .svg / .png\n")
    cat("  - model code      : acn_mrgsolve_model.R\n")
    cat("  - references      : acn_references_en.md\n\n")
    cat("Reading the model in one sentence:\n")
    cat("  every visible lesion descends from an invisible microcomedone (MC),\n")
    cat("  so anti-inflammatory therapy empties the downstream compartments fast\n")
    cat("  and rebounds, while retinoids drain the reservoir slowly and hold.\n\n")
    cat("Compartments:\n")
    print(names(mrgsolve::init(MOD)))
  })
}

shinyApp(ui, server)
