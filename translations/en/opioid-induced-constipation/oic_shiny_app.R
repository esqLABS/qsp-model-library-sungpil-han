## =============================================================================
##  oic_shiny_app.R
##  Opioid-Induced Constipation (OIC) QSP — interactive dashboard
##  Opioid-Induced Constipation QSP — Interactive Dashboard
##
##  10 tabs:
##    1  Patient profile   Patient profile & regimen
##    2  PK              Opioid and antagonist pharmacokinetics
##    3  Receptor occupancy   Receptor occupancy — the two compartments side by side
##    4  Selectivity index, SI     Selectivity index & the P-gp drug-interaction explorer
##    5  Colonic transit       Colonic transit, segment loads and hydration
##    6  Clinical endpoints     SBM / CSBM / Bristol / straining / PAC-SYM
##    7  Scenario comparison   Scenario comparison (16 regimens)
##    8  Brake decomposition   Which brake carries the drug response?
##    9  Tolerance titration       Tolerance titration over 24 weeks
##   10  Biomarkers      Biomarker panel & safety
##
##  REQUIRES oic_mrgsolve_model.R in the same directory.
##
##  NOTE ON PROVENANCE: like the model file, this app has never been executed --
##  the build container had no R.  It is written against the mrgsolve model's
##  documented output columns.  See README.md.
## =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)   # grid.arrange, used on tabs 5 and 6

source("oic_mrgsolve_model.R", local = TRUE)

THEME <- theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

DRUGS <- c("None (none)"                    = "none",
           "Naloxegol"             = "naloxegol",
           "Naldemedine"           = "naldemedine",
           "Methylnaltrexone SC" = "methylnaltrexone_sc",
           "Methylnaltrexone PO 450 mg"          = "methylnaltrexone_po",
           "Naloxone PO (control)"      = "naloxone_po")

DEFAULT_DOSE <- c(none = 0, naloxegol = 25, naldemedine = 0.2,
                  methylnaltrexone_sc = 12, methylnaltrexone_po = 450,
                  naloxone_po = 20)

## =============================================================================
##  UI
## =============================================================================
ui <- fluidPage(
  titlePanel("Opioid-Induced Constipation (OIC) — QSP Simulator"),
  tags$p(style = "color:#555;margin-top:-8px",
         "One receptor, two compartments, one barrier. ",
         tags$b("For education and research purposes only — not for use in clinical decision-making.")),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Opioid"),
      sliderInput("op_dose", "Oxycodone dose per administration (mg, q12h)",
                  min = 0, max = 120, value = 30, step = 5),
      checkboxInput("methadone", "Substitute methadone (direct ClC-2 block)", FALSE),

      h4("Antagonist"),
      selectInput("drug", "Drug", choices = DRUGS, selected = "naloxegol"),
      numericInput("pam_dose", "Dose per administration (mg)", value = 25, min = 0, step = 0.1),
      numericInput("pam_int", "Dosing interval (h)", value = 24, min = 1, step = 1),

      h4("Adjuncts"),
      numericInput("peg_g",  "PEG 3350 (g/day)",     value = 0, min = 0, step = 1),
      numericInput("lac_g",  "Lactulose (g/day)",       value = 0, min = 0, step = 5),
      numericInput("lub_ug", "Lubiprostone (µg, bid)", value = 0, min = 0, step = 8),
      numericInput("lin_ug", "Linaclotide (µg/day)", value = 0, min = 0, step = 72),
      numericInput("pro_mg", "Prucalopride (mg/day)", value = 0, min = 0, step = 1),

      h4("Drug interaction & patient factors"),
      sliderInput("pgp", "P-gp inhibition (Kp,uu multiple)", min = 1, max = 20,
                  value = 1, step = 1),
      sliderInput("cyp", "CYP3A4 inhibition (AUC multiple)", min = 1, max = 6,
                  value = 1, step = 0.2),
      sliderInput("fluid", "Fluid intake multiple", min = 0.6, max = 1.6,
                  value = 1, step = 0.05),
      sliderInput("absx", "Colonic water-absorption capacity multiple", min = 0.7, max = 1.5,
                  value = 1, step = 0.05),
      checkboxInput("rescue", "Allow rescue laxative (trial protocol)", TRUE),

      sliderInput("days", "Simulation period (days)", min = 14, max = 168,
                  value = 84, step = 7),
      actionButton("go", "Run simulation", class = "btn-primary")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",

        tabPanel("1 · Patient profile",
          br(),
          fluidRow(
            column(4, wellPanel(h4("SBM / week"), h2(textOutput("kpi_sbm")),
                                tags$small("Responder threshold ≥3"))),
            column(4, wellPanel(h4("Bristol"), h2(textOutput("kpi_bsfs")))),
            column(4, wellPanel(h4("PAC-SYM"), h2(textOutput("kpi_pacsym"))))
          ),
          fluidRow(
            column(4, wellPanel(h4("Pain NRS"), h2(textOutput("kpi_pain")))),
            column(4, wellPanel(h4("COWS withdrawal"), h2(textOutput("kpi_cows")),
                                tags$small("Safety endpoint"))),
            column(4, wellPanel(h4("Rescue medication / week"), h2(textOutput("kpi_resc"))))
          ),
          hr(),
          h4("Summary"), tableOutput("tbl_profile")),

        tabPanel("2 · PK",
          br(), plotOutput("p_pk", height = "620px"),
          tags$p(tags$small(
            "Top: opioid plasma & free-brain concentration. Bottom: antagonist plasma & free-brain concentration. ",
            "The ratio of the two brain curves is the whole of selectivity."))),

        tabPanel("3 · Receptor occupancy",
          br(), plotOutput("p_occ", height = "620px"),
          tags$p(tags$small(
            "Same receptor, two compartments. The gap between the intestinal curve and the central curve is the therapeutic window."))),

        tabPanel("4 · Selectivity index · Drug interaction",
          br(),
          h4("Selectivity index SI = intestinal antagonist occupancy / central antagonist occupancy"),
          tableOutput("tbl_si"), hr(),
          h4("P-gp inhibitors change the 'ratio', not the exposure"),
          tableOutput("tbl_ddi"),
          tags$p(tags$small(
            "CYP3A4 inhibition is a pure exposure shift, so dose reduction restores most of it. ",
            "P-gp inhibition moves Kp,uu itself, so dose reduction barely recovers safety while ",
            "efficacy alone is lost — a case where the label's dose-reduction instruction targets the wrong variable."))),

        tabPanel("5 · Colonic transit",
          br(), plotOutput("p_transit", height = "620px"),
          tags$p(tags$small(
            "Solids (g) and water fraction by segment. The lower bound of water fraction (≈0.62) is bound water in the solid phase, and ",
            "without this lower bound the model diverges even in healthy patients."))),

        tabPanel("6 · Clinical endpoints",
          br(), plotOutput("p_end", height = "620px"),
          tags$p(tags$small(
            "SBM is defined as 'a bowel movement without rescue-laxative use within 24 hours'. ",
            "In this model a bowel movement induced by rescue medication empties the colon but is not counted as an SBM."))),

        tabPanel("7 · Scenario comparison",
          br(), actionButton("go_scen", "Run 16 scenarios", class = "btn-default"),
          br(), br(), plotOutput("p_scen", height = "560px"),
          tableOutput("tbl_scen")),

        tabPanel("8 · Brake decomposition",
          br(), actionButton("go_brake", "Run brake decomposition", class = "btn-default"),
          br(), br(), tableOutput("tbl_brake"),
          tags$p(tags$small(
            "Each brake is switched off one at a time and the same 0 → 25 mg naloxegol comparison is re-measured. ",
            "The term whose multiple approaches 1 is the term carrying the drug effect. In the reference implementation, ",
            "even completely removing the transit–water feedback loop, the multiple went from 2.89 → 2.91, ",
            "barely changing — that loop is not an amplifier."))),

        tabPanel("9 · Tolerance titration",
          br(), actionButton("go_tol", "Run 24-week titration", class = "btn-default"),
          br(), br(), plotOutput("p_tol", height = "520px"),
          tableOutput("tbl_tol"),
          tags$p(tags$small(
            "Dose is adjusted weekly to keep pain at NRS 4.0. ",
            "Central receptor availability collapses while intestinal availability is preserved."))),

        tabPanel("10 · Biomarkers · Safety",
          br(), plotOutput("p_bio", height = "620px"),
          tableOutput("tbl_bio"))
      )
    )
  )
)

## =============================================================================
##  SERVER
## =============================================================================
server <- function(input, output, session) {

  observeEvent(input$drug, {
    updateNumericInput(session, "pam_dose",
                       value = unname(DEFAULT_DOSE[input$drug]))
    updateNumericInput(session, "pam_int",
                       value = if (input$drug == "naloxone_po") 8 else 24)
  })

  sim <- eventReactive(input$go, {
    pars <- list(PGPINH = input$pgp, CYP3A4I = input$cyp,
                 FLUIDX = input$fluid, ABSX = input$absx,
                 METHADON = as.numeric(input$methadone),
                 RESCUEON = as.numeric(input$rescue))
    run_oic(drug = input$drug, days = input$days, params = pars,
            op_dose = input$op_dose, op_int = 12,
            pam_dose = input$pam_dose, pam_int = input$pam_int,
            pro_dose = input$pro_mg, peg_g = input$peg_g, lac_g = input$lac_g,
            lub_ug = input$lub_ug, lin_ug = input$lin_ug)
  }, ignoreNULL = FALSE)

  smry <- reactive(summarise_oic(sim()))

  fmt <- function(x, d = 2) formatC(x, digits = d, format = "f")
  output$kpi_sbm    <- renderText(fmt(smry()$SBM_wk))
  output$kpi_bsfs   <- renderText(fmt(smry()$BSFS))
  output$kpi_pacsym <- renderText(fmt(smry()$PACSYM))
  output$kpi_pain   <- renderText(fmt(smry()$PAIN))
  output$kpi_cows   <- renderText(fmt(smry()$COWS))
  output$kpi_resc   <- renderText(fmt(smry()$RESC_wk))

  output$tbl_profile <- renderTable({
    s <- smry()
    tibble(
      Metric = c("SBM/week", "CSBM/week", "Rescue medication/week", "Bristol", "Straining", "Bloating",
               "PAC-SYM", "PAC-QOL", "Pain NRS", "COWS",
               "Intestinal agonist occupancy", "Intestinal antagonist occupancy",
               "Central antagonist occupancy", "Selectivity index, SI",
               "Rectal water fraction", "Colonic transit time (h)", "Impaction-risk integral"),
      Value = c(s$SBM_wk, s$CSBM_wk, s$RESC_wk, s$BSFS, s$STRAIN, s$DIST,
             s$PACSYM, s$QOL, s$PAIN, s$COWS, s$OCCG_AG, s$OCCG_ANT,
             s$OCCC_ANT, s$SI, s$W4, s$CTT_H, s$IMPACT))
  }, digits = 4)

  ## ---- 2 PK ----------------------------------------------------------------
  output$p_pk <- renderPlot({
    d <- sim() %>% filter(time >= max(time) - 168)
    d %>%
      transmute(time = time - min(time),
                `Opioid plasma (mg/L)` = CP_OP,
                `Opioid free-brain (nM)` = COP_BR,
                `Antagonist plasma (mg/L)` = CP_PAM,
                `Antagonist free-brain (nM)` = CPAM_BR) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) +
      geom_line(linewidth = 0.8, colour = "#4f88a8") +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "Time in final week (h)", y = NULL,
           title = "PK — final 7 days") + THEME
  })

  ## ---- 3 occupancy ---------------------------------------------------------
  output$p_occ <- renderPlot({
    sim() %>%
      transmute(day = time / 24,
                `Intestinal MOR — agonist` = OCCG_AG,
                `Intestinal MOR — antagonist` = OCCG_ANT,
                `Central MOR — agonist` = OCCC_AG,
                `Central MOR — antagonist` = OCCC_ANT) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value, colour = name)) +
      geom_line(linewidth = 0.8) +
      scale_y_continuous(limits = c(0, 1)) +
      scale_colour_manual(values = c("#c0392b", "#2e8b57", "#8b7bb8", "#c9a227")) +
      labs(x = "Day", y = "Occupancy", colour = NULL,
           title = "Same receptor, two compartments") + THEME
  })

  ## ---- 4 selectivity & DDI -------------------------------------------------
  output$tbl_si  <- renderTable(selectivity_table(days = 56), digits = 5)
  output$tbl_ddi <- renderTable(ddi_table(days = 56), digits = 4)

  ## ---- 5 transit -----------------------------------------------------------
  output$p_transit <- renderPlot({
    d <- sim() %>% mutate(day = time / 24)
    p1 <- d %>%
      select(day, `Ascending` = S1, `Transverse` = S2, `Descending` = S3, `Rectosigmoid` = S4) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value, colour = name)) + geom_line(linewidth = 0.8) +
      labs(x = NULL, y = "Solids (g)", colour = NULL,
           title = "Solid load by segment") + THEME
    p2 <- d %>%
      transmute(day,
                `Ascending` = W1/(W1+S1), `Transverse` = W2/(W2+S2),
                `Descending` = W3/(W3+S3), `Rectosigmoid` = W4/(W4+S4)) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value, colour = name)) + geom_line(linewidth = 0.8) +
      geom_hline(yintercept = 0.615, linetype = "dashed", colour = "#c0392b") +
      labs(x = "Day", y = "Water fraction", colour = NULL,
           title = "Water fraction by segment (red line = bound-water lower bound)") + THEME
    gridExtra::grid.arrange(p1, p2, ncol = 1)
  })

  ## ---- 6 endpoints ---------------------------------------------------------
  output$p_end <- renderPlot({
    d <- sim() %>% mutate(day = time / 24)
    wk <- d %>% mutate(week = floor(day / 7)) %>% group_by(week) %>%
      summarise(SBM = max(CUM_SBM) - min(CUM_SBM),
                CSBM = max(CUM_CSBM) - min(CUM_CSBM),
                RESC = max(CUM_RESC) - min(CUM_RESC), .groups = "drop")
    p1 <- wk %>% pivot_longer(-week) %>%
      ggplot(aes(week, value, fill = name)) +
      geom_col(position = "dodge") +
      geom_hline(yintercept = 3, linetype = "dashed", colour = "#c0392b") +
      labs(x = NULL, y = "Count per week", fill = NULL,
           title = "Weekly bowel movements — red line = responder threshold") + THEME
    p2 <- d %>% select(day, Bristol = BSFS, `Straining` = STRAIN,
                       `Bloating` = DIST, `PAC-SYM` = PACSYM) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value, colour = name)) + geom_line(linewidth = 0.8) +
      labs(x = "Day", y = "Score", colour = NULL, title = "Symptom scores") + THEME
    gridExtra::grid.arrange(p1, p2, ncol = 1)
  })

  ## ---- 7 scenarios ---------------------------------------------------------
  scen <- eventReactive(input$go_scen, run_all_scenarios(days = 56))
  output$p_scen   <- renderPlot(plot_scenarios(scen()))
  output$tbl_scen <- renderTable(scen(), digits = 3)

  ## ---- 8 brake decomposition ----------------------------------------------
  brake <- eventReactive(input$go_brake, brake_decomposition(days = 56))
  output$tbl_brake <- renderTable(brake(), digits = 3)

  ## ---- 9 tolerance ---------------------------------------------------------
  tol <- eventReactive(input$go_tol, tolerance_titration(weeks = 24))
  output$p_tol <- renderPlot({
    tol() %>%
      select(week, `Oxycodone mg/day` = oxy_mg_day,
             `Central receptor availability` = OCCG_AG, `SBM/week` = SBM_wk,
             `Pain NRS` = PAIN) %>%
      pivot_longer(-week) %>%
      ggplot(aes(week, value)) + geom_line(linewidth = 0.9, colour = "#8b7bb8") +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "Week", y = NULL, title = "24-week titration — attempt to maintain analgesia") + THEME
  })
  output$tbl_tol <- renderTable(tol(), digits = 3)

  ## ---- 10 biomarkers -------------------------------------------------------
  output$p_bio <- renderPlot({
    sim() %>% mutate(day = time / 24) %>%
      select(day, `Intestinal MOR availability` = RG_AV, `Central MOR availability` = RC_AV,
             `Intestinal cAMP` = CAMP, `ACh release` = ACH,
             `HAPC (1/h)` = HAPC, `Segmental tone` = TONE,
             `ClC-2 activity` = CLC2, `Impaction risk` = IMPACT) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value)) + geom_line(linewidth = 0.8, colour = "#4f9a9a") +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "Day", y = NULL, title = "Mechanistic biomarkers") + THEME
  })
  output$tbl_bio <- renderTable({
    d <- tail(sim(), 1)
    tibble(Biomarker = c("Intestinal MOR availability", "Central MOR availability", "Intestinal cAMP",
                          "HAPC (1/h)", "Segmental tone", "Total colonic solids (g)",
                          "Colonic transit time (h)", "Impaction-risk integral"),
           Value = c(d$RG_AV, d$RC_AV, d$CAMP, d$HAPC, d$TONE, d$STOT,
                  d$CTT_H, d$IMPACT))
  }, digits = 4)
}

shinyApp(ui, server)
