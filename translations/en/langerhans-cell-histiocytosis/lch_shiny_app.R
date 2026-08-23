# =============================================================================
# lch_shiny_app.R
# Interactive dashboard for the Langerhans Cell Histiocytosis QSP model
# -----------------------------------------------------------------------------
# Run with:  shiny::runApp("lch_shiny_app.R")
# Requires:  shiny, mrgsolve, dplyr, tidyr, ggplot2, DT
#
# The app sources lch_mrgsolve_model.R for the model, phenotypes, regimen
# builders and the run_scenario() driver, so there is exactly one definition of
# the model in this directory.
#
# Eleven tabs:
#   1  Patient profile        Patient / cell-of-origin profile
#   2  Drug PK              Drug pharmacokinetics
#   3  MAPK signalling            pERK, cyclin D1, BCL2A1, SASP
#   4  Cell compartments · Lesions     Precursor reservoir and lesional burden
#   5  Cytokines · Secretome  Secretome
#   6  Biomarkers           cfDNA, sCD163, CRP, ferritin
#   7  Clinical endpoints      DAS, response, reactivation
#   8  Permanent sequelae          Irreversible sequelae (CDI, ND, biliary, lung, bone)
#   9  Toxicity                 ANC (Friberg), skin, LVEF
#  10  Scenario comparison        Scenario comparison
#  11  Structural-choice validation        Falsification / kill-switch panel
# =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

source("lch_mrgsolve_model.R", local = TRUE)

THEME <- theme_bw(base_size = 12) +
  theme(strip.background = element_rect(fill = "grey92"),
        legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 12))

PHENO_LABEL <- c(
  SSb      = "SS-b single-system bone lesion (tissue-restricted origin, FSR 0.10)",
  MSROneg  = "MS RO-negative multisystem (risk organs not involved)",
  MSROpos  = "MS RO-positive multisystem (HSC origin, risk organs involved)",
  CNSrisk  = "Craniofacial CNS-risk + pituitary involvement",
  PLCH     = "Adult pulmonary LCH (smoker)"
)

SCEN_LABEL <- setNames(
  vapply(scenarios, function(s) s$label, character(1)),
  names(scenarios)
)

long_of <- function(df, cols) {
  df %>% select(time, all_of(cols)) %>% pivot_longer(-time)
}

# =============================================================================
# UI
# =============================================================================
ui <- fluidPage(
  titlePanel("Langerhans Cell Histiocytosis (LCH) — QSP Model Dashboard"),
  p(tags$em(paste(
    "62 ODEs · Cell-of-origin partitioning · ERK-driven secretome · Cytostatic MAPK inhibition, and",
    "relapse after discontinuation · permanent sequelae as a time integral"))),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectInput("scen", "Scenario",
                  choices = names(SCEN_LABEL), selected = "S5_vem_continuous"),
      helpText(textOutput("scen_desc")),
      hr(),
      h4("Patient · Cell of Origin"),
      selectInput("geno", "Driver genotype",
                  choices = c("BRAF V600E" = 1, "MAP2K1 / ARAF / negative" = 2),
                  selected = 1),
      sliderInput("wt", "Weight (kg)", 5, 80, 12, step = 1),
      sliderInput("fsr", "FSR — cell-of-origin self-renewal capacity", 0, 1, 1.0, step = 0.05),
      sliderInput("precm0", "PRECM0 — initial mutant precursor reservoir",
                  0, 1.5, 0.45, step = 0.05),
      sliderInput("thr", "THR — risk-organ seeding fraction", 0, 0.7, 0.34,
                  step = 0.02),
      sliderInput("thp", "THP — pituitary seeding fraction", 0, 0.5, 0.08,
                  step = 0.02),
      checkboxInput("smoke", "Smoking exposure (drives pulmonary LCH)", FALSE),
      hr(),
      h4("Pharmacodynamics · Structural Switches"),
      sliderInput("ic50vem", "IC50 vemurafenib (free, mg/L)",
                  0.005, 0.08, 0.018, step = 0.001),
      sliderInput("kprol", "KPROL — intralesional proliferation rate (1/day)",
                  0, 0.06, 0.022, step = 0.002),
      helpText(HTML(paste(
        "If KPROL exceeds <b>0.0288</b> (the minimum immune clearance rate), the lesion becomes",
        "self-sustaining without any precursor supply, contradicting the low Ki-67 (≈9%)."))),
      sliderInput("kill", "SL_MAPKI_KILL — falsification switch",
                  0, 0.4, 0.0, step = 0.05),
      helpText(HTML(paste(
        "If nonzero, the MAPK inhibitor becomes <b>cytocidal</b>, and relapse after discontinuation",
        "disappears — contradicting the literature."))),
      sliderInput("fnapo", "FN_APO — niche protection from apoptosis",
                  0, 1, 0.25, step = 0.05),
      sliderInput("pniche", "PNICHE — clonal extinction threshold",
                  0.0005, 0.02, 0.003, step = 0.0005),
      sliderInput("kavp", "KAVP — AVP neuron loss rate (1/day)",
                  0, 0.04, 0.010, step = 0.002),
      hr(),
      actionButton("go", "Recalculate (Simulate)", class = "btn-primary"),
      br(), br(),
      helpText("Sliders override the phenotype defaults. If you change the scenario,",
               "refresh the page to reset to that phenotype's values.")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("1. Patient Profile",
                 h4("Cell-of-origin partitioning — the model's first structural choice"),
                 plotOutput("p_partition", height = "300px"),
                 h4("Status at diagnosis (at end of run-in)"),
                 DTOutput("t_baseline")),
        tabPanel("2. Drug PK",
                 plotOutput("p_pk", height = "620px"),
                 helpText("PK parameters are referenced to 70 kg and are corrected for body weight",
                          "using allometric scaling (CL^0.75, V^1.0).")),
        tabPanel("3. MAPK Signalling",
                 plotOutput("p_mapk", height = "600px"),
                 helpText(HTML(paste(
                   "pERK must fall to <b>less than 20%</b> of baseline (i.e., more than 80% inhibition) for",
                   "cyclin D1 to drop sharply — this is the threshold created by a Hill coefficient of 2",
                   "(Bollag 2010, PMID 20823850)。")))),
        tabPanel("4. Cell Compartments · Lesions",
                 plotOutput("p_cells", height = "620px"),
                 helpText(HTML(paste(
                   "Key contrast: lesional mass disappears, yet the <b>bone-marrow precursor reservoir</b>",
                   "remains. That remaining reservoir is the source of relapse after discontinuation.")))),
        tabPanel("5. Cytokines · Secretome",
                 plotOutput("p_cyto", height = "620px")),
        tabPanel("6. Biomarkers",
                 plotOutput("p_bio", height = "560px"),
                 helpText("cfDNA tracks the death flux of mutant cells, so even when mass",
                          "remains, a cytostatic drug produces a gentle plateau",
                          "instead.")),
        tabPanel("7. Clinical Endpoints",
                 plotOutput("p_das", height = "420px"),
                 h4("Response · Reactivation Summary"),
                 DTOutput("t_response")),
        tabPanel("8. Permanent Sequelae",
                 plotOutput("p_seq", height = "620px"),
                 helpText(HTML(paste(
                   "The five pools (AVPN·ANTPIT·NEUR·BILF·LUNGC) are <b>monotone",
                   "(monotone)</b> — no therapy can reverse them.",
                   "Therefore TTET (days of active disease) explains sequelae",
                   "better than drug efficacy does.")))),
        tabPanel("9. Toxicity",
                 plotOutput("p_tox", height = "620px")),
        tabPanel("10. Scenario Comparison",
                 checkboxGroupInput(
                   "cmp", "Scenarios to compare",
                   choices = names(SCEN_LABEL),
                   selected = c("S5_vem_continuous", "S6_vem_stop",
                                "S4_cladarac_upfront",
                                "S8_bridge_consolidate"),
                   inline = TRUE),
                 actionButton("go_cmp", "Run comparison", class = "btn-primary"),
                 plotOutput("p_cmp", height = "560px"),
                 DTOutput("t_cmp")),
        tabPanel("11. Structural-Choice Validation",
                 h4("Falsification Panel"),
                 p(HTML(paste(
                   "Each structural choice has a <b>kill switch</b>.",
                   "When the switch is turned on, the model must show behaviour that contradicts the literature",
                   "— that is the evidence that the choice is actually doing work."))),
                 actionButton("go_kill", "Run all 3 kill switches",
                              class = "btn-primary"),
                 plotOutput("p_kill", height = "620px"),
                 verbatimTextOutput("t_kill"))
      )
    )
  )
)

# =============================================================================
# SERVER
# =============================================================================
server <- function(input, output, session) {

  output$scen_desc <- renderText(SCEN_LABEL[[input$scen]])

  overrides <- reactive({
    list(GENO = as.numeric(input$geno), WT = input$wt, FSR = input$fsr,
         PRECM0 = input$precm0, THR = input$thr, THP = input$thp,
         SMOKE = as.numeric(input$smoke), IC50_VEM = input$ic50vem,
         KPROL = input$kprol, SL_MAPKI_KILL = input$kill,
         FN_APO = input$fnapo, PNICHE = input$pniche, KAVP = input$kavp)
  })

  sim <- eventReactive(input$go, {
    withProgress(message = "Integrating (62 ODEs)...", value = 0.5, {
      run_scenario(input$scen, overrides())
    })
  }, ignoreNULL = FALSE)

  # ---- 1. patient profile -------------------------------------------------
  output$p_partition <- renderPlot({
    p <- modifyList(pheno[[scenarios[[input$scen]]$kind]], overrides())
    th <- c(Bone = p$THB, Skin = p$THS, `Risk organs` = p$THR,
            Pituitary = p$THP, CNS = p$THC, Lung = p$THL)
    th[is.na(th)] <- 0
    data.frame(site = names(th), theta = as.numeric(th)) %>%
      ggplot(aes(reorder(site, theta), theta, fill = site)) +
      geom_col(show.legend = FALSE) +
      coord_flip() +
      labs(title = paste0("Seeding partition  (FSR = ", p$FSR,
                          ",  PRECM0 = ", p$PRECM0, ")"),
           x = NULL, y = "theta (fraction of circulating progeny)") +
      THEME
  })

  output$t_baseline <- renderDT({
    d <- sim()
    b <- d[1, ]
    data.frame(
      quantity = c("DAS at diagnosis", "Total lesional burden (units)",
                   "Risk-organ burden", "Pituitary burden", "CNS burden",
                   "cfDNA signal", "AVP neuron pool", "Neuron pool",
                   "Biliary fibrosis", "ANC (1e9/L)"),
      value = round(c(b$DAS, b$LTOT_units, b$LRO, b$LPIT, b$LCNS, b$CFDNA,
                      b$AVPN, b$NEUR, b$BILF, b$ANC), 4)
    ) %>% datatable(rownames = FALSE, options = list(dom = "t"))
  })

  # ---- 2. PK --------------------------------------------------------------
  output$p_pk <- renderPlot({
    long_of(sim(), c("CVEM_mgL", "CDAB_mgL", "CTRA_ngmL", "CVBL_mgL",
                     "CPRE_mgL", "ARATP", "CLATP", "VEMI")) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.5) +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(title = "Drug concentration · Intracellular active metabolite · CYP3A4 autoinduction",
           x = "Days since diagnosis", y = NULL) + THEME
  })

  # ---- 3. MAPK ------------------------------------------------------------
  output$p_mapk <- renderPlot({
    d <- sim()
    long_of(d, c("pERK_pct", "CCND", "BCL", "SASP", "AUCERK")) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.6) +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      geom_hline(data = data.frame(name = "pERK_pct", y = 20),
                 aes(yintercept = y), linetype = 2, colour = "red") +
      labs(title = "MAPK signalling output (red line = 80% inhibition threshold)",
           x = "Days since diagnosis", y = NULL) + THEME
  })

  # ---- 4. cells -----------------------------------------------------------
  output$p_cells <- renderPlot({
    d <- sim()
    long_of(d, c("PRECM", "CIRC", "LBONE", "LSKIN", "LRO", "LPIT", "LCNS",
                 "LLUNG", "TREG", "OCL")) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.6) +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(title = "Precursor reservoir and organ-specific lesional burden",
           x = "Days since diagnosis", y = "burden units") + THEME
  })

  # ---- 5. secretome -------------------------------------------------------
  output$p_cyto <- renderPlot({
    long_of(sim(), c("IL1B", "TNFA", "IL6", "OSM", "MMP9", "RANKL")) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.6) +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(title = "ERK/SASP-driven secretome", x = "Days since diagnosis", y = NULL) +
      THEME
  })

  # ---- 6. biomarkers ------------------------------------------------------
  output$p_bio <- renderPlot({
    d <- sim()
    lod <- as.numeric(param(mod)$CF_LOD)
    long_of(d, c("CFDNA", "SCD163", "CRP", "FERR")) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.6) +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      geom_hline(data = data.frame(name = "CFDNA", y = lod),
                 aes(yintercept = y), linetype = 2, colour = "red") +
      labs(title = "Biomarkers (red line = ddPCR limit of detection)",
           x = "Days since diagnosis", y = NULL) + THEME
  })

  # ---- 7. endpoints -------------------------------------------------------
  output$p_das <- renderPlot({
    d <- sim()
    ggplot(d, aes(time, DAS)) +
      geom_line(linewidth = 0.7) +
      geom_hline(yintercept = as.numeric(param(mod)$DAS_ACTIVE),
                 linetype = 2, colour = "red") +
      labs(title = "Disease Activity Score (red line = active disease threshold)",
           x = "Days since diagnosis", y = "DAS") + THEME
  })

  output$t_response <- renderDT({
    d <- sim()
    d0 <- d$DAS[1]
    at_t <- function(tt) d[which.min(abs(d$time - tt)), ]
    w6 <- at_t(42); w12 <- at_t(84); end <- d[nrow(d), ]
    resp <- function(x) {
      if (x$DAS < 0.1 * d0) "NAD / better"
      else if (x$DAS < d0) "intermediate"
      else "worse"
    }
    data.frame(
      timepoint = c("baseline", "week 6", "week 12", "end of follow-up"),
      DAS = round(c(d0, w6$DAS, w12$DAS, end$DAS), 2),
      cfDNA = round(c(d$CFDNA[1], w6$CFDNA, w12$CFDNA, end$CFDNA), 4),
      reservoir = round(c(d$PRECM[1], w6$PRECM, w12$PRECM, end$PRECM), 4),
      category = c("-", resp(w6), resp(w12), resp(end)),
      TTET_days = round(c(0, w6$TTET, w12$TTET, end$TTET), 0)
    ) %>% datatable(rownames = FALSE, options = list(dom = "t"))
  })

  # ---- 8. sequelae --------------------------------------------------------
  output$p_seq <- renderPlot({
    d <- sim()
    long_of(d, c("AVPN", "ANTPIT", "NEUR", "BILF", "LUNGC", "BVOL",
                 "TTET", "CUMDAS")) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.6) +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(title = "Irreversible organ pools and cumulative active disease",
           x = "Days since diagnosis", y = NULL) + THEME
  })

  # ---- 9. toxicity -------------------------------------------------------
  output$p_tox <- renderPlot({
    d <- sim()
    long_of(d, c("ANC", "PROL", "SKTOX", "LVEF")) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.6) +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      geom_hline(data = data.frame(name = "ANC", y = 0.5),
                 aes(yintercept = y), linetype = 2, colour = "red") +
      labs(title = "Myelosuppression (Friberg) · Paradoxical skin toxicity · LVEF",
           x = "Days since diagnosis", y = NULL) + THEME
  })

  # ---- 10. scenario comparison -------------------------------------------
  cmp <- eventReactive(input$go_cmp, {
    req(length(input$cmp) > 0)
    withProgress(message = "Comparing scenarios...", value = 0.4, {
      bind_rows(lapply(input$cmp, function(s) run_scenario(s, overrides())))
    })
  })

  output$p_cmp <- renderPlot({
    cmp() %>%
      select(time, scenario, DAS, CFDNA, PRECM, AVPN, NEUR, BVOL) %>%
      pivot_longer(-c(time, scenario)) %>%
      ggplot(aes(time, value, colour = scenario)) +
      geom_line(linewidth = 0.6) +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(title = "Scenario comparison", x = "Days since diagnosis", y = NULL) + THEME
  })

  output$t_cmp <- renderDT({
    cmp() %>% group_by(scenario) %>% slice_tail(n = 1) %>%
      transmute(scenario,
                DAS = round(DAS, 2), cfDNA = round(CFDNA, 4),
                reservoir = round(PRECM, 5),
                AVPN = round(AVPN, 3), NEUR = round(NEUR, 3),
                BILF = round(BILF, 2), LUNGC = round(LUNGC, 2),
                CDI, ND, TTET = round(TTET, 0)) %>%
      datatable(rownames = FALSE, options = list(dom = "t", scrollX = TRUE))
  })

  # ---- 11. falsification panel -------------------------------------------
  kill <- eventReactive(input$go_kill, {
    withProgress(message = "Running all 3 kill switches...", value = 0.3, {
      list(
        c1_ss   = run_scenario("S1_observation"),
        c1_swap = run_scenario("S1_observation", pheno$MSROpos),
        c2_base = run_scenario("S6_vem_stop"),
        c2_kill = run_scenario("S6_vem_stop", list(SL_MAPKI_KILL = 0.25)),
        c3_late = run_scenario("S9_delayed_dx"),
        c3_early = run_scenario("S9b_early_dx")
      )
    })
  })

  output$p_kill <- renderPlot({
    k <- kill()
    bind_rows(
      k$c1_ss    %>% transmute(time, y = LRO,   panel = "① Partitioning: SS-b origin descriptor",         arm = "baseline"),
      k$c1_swap  %>% transmute(time, y = LRO,   panel = "① Partitioning: swapped to MS RO+ origin descriptor", arm = "switch"),
      k$c2_base  %>% transmute(time, y = CFDNA, panel = "② Cytostasis: SL_MAPKI_KILL = 0",     arm = "baseline"),
      k$c2_kill  %>% transmute(time, y = CFDNA, panel = "② Cytostasis: SL_MAPKI_KILL = 0.25",  arm = "switch"),
      k$c3_early %>% transmute(time, y = AVPN,  panel = "③ Time integral: treatment started at day 14",        arm = "baseline"),
      k$c3_late  %>% transmute(time, y = AVPN,  panel = "③ Time integral: 180-day delay",              arm = "switch")
    ) %>%
      ggplot(aes(time, y, colour = arm)) + geom_line(linewidth = 0.7) +
      facet_wrap(~panel, scales = "free_y", ncol = 2) +
      labs(title = "Kill switches for the three structural choices",
           x = "Days since diagnosis", y = NULL) + THEME
  })

  output$t_kill <- renderPrint({
    k <- kill()
    last <- function(d, v) d[[v]][nrow(d)]
    cat("① Only the cell-of-origin partitioning swapped (identical lesional rate constants)\n")
    cat(sprintf("   Risk-organ burden d730:  SS-b = %.4f   ->  MS RO+ descriptor = %.4f\n",
                last(k$c1_ss, "LRO"), last(k$c1_swap, "LRO")))
    cat("\n② Giving the MAPK inhibitor a cytocidal term (SL_MAPKI_KILL = 0.25)\n")
    cat(sprintf("   Peak DAS after discontinuation:  0 = %.2f   ->  0.25 = %.2f  (relapse abolished)\n",
                max(k$c2_base$DAS[k$c2_base$time > 365]),
                max(k$c2_kill$DAS[k$c2_kill$time > 365])))
    cat(sprintf("   Reservoir d365:         0 = %.4f   ->  0.25 = %.6f\n",
                k$c2_base$PRECM[which.min(abs(k$c2_base$time - 365))],
                k$c2_kill$PRECM[which.min(abs(k$c2_kill$time - 365))]))
    cat("\n③ Only the treatment start time changed (day 14 vs day 180)\n")
    cat(sprintf("   AVP neuron pool d730:    day 14 = %.3f  ->  day 180 = %.3f\n",
                last(k$c3_early, "AVPN"), last(k$c3_late, "AVPN")))
    cat(sprintf("   Central diabetes insipidus:       day 14 = %d      ->  day 180 = %d\n",
                as.integer(last(k$c3_early, "CDI")),
                as.integer(last(k$c3_late, "CDI"))))
    cat(sprintf("   Days of active disease, TTET: day 14 = %.0f     ->  day 180 = %.0f\n",
                last(k$c3_early, "TTET"), last(k$c3_late, "TTET")))
  })
}

shinyApp(ui, server)
