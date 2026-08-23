##############################################################################
# Fabry Disease QSP — Shiny Interactive Dashboard
# Fabry Disease Quantitative Systems Pharmacology Model — Shiny Dashboard
# 8 Tabs: Patient Profile · PK/Enzyme · Gb3 Dynamics · Renal Function · Cardiac Function ·
#         Treatment Scenario Comparison · Biomarker Panel · Virtual Patient Population
##############################################################################

library(shiny)
library(shinydashboard)
library(ggplot2)
library(dplyr)
library(tidyr)
library(DT)

# ── Colour palette ────────────────────────────────────────────────────────
scenario_colors <- c(
  "Natural history"       = "#616161",
  "Agalsidase beta"       = "#1565C0",
  "Agalsidase alfa"       = "#0288D1",
  "Migalastat"            = "#2E7D32",
  "Pegunigalsidase alfa"  = "#7B1FA2",
  "ERT+Lucerastat"        = "#E65100"
)

##############################################################################
# CORE SIMULATION FUNCTION (simplified analytical / semi-ODE)
##############################################################################
simulate_fabry <- function(
    phenotype      = "classic",   # "classic" | "late"
    is_amenable    = FALSE,       # migalastat amenable mutation
    treatment      = "none",      # "none"|"agaB"|"agaA"|"migalastat"|"pegu"|"combo"
    agaB_dose      = 1.0,         # mg/kg
    agaA_dose      = 0.2,         # mg/kg
    mig_dose       = 150,         # mg QOD
    luc_dose       = 1000,        # mg TID
    sim_years      = 5,
    age_start      = 30,
    sex            = "male",
    eGFR0          = 90
) {
  t <- seq(0, sim_years * 365, by = 7)  # weekly steps
  n <- length(t)

  # Baseline enzyme activity
  E0 <- if (phenotype == "classic") 0.05 else 2.5   # nmol/h/mg

  # Treatment-specific enzyme restoration
  E_max_treat <- switch(treatment,
    "none"       = 0,
    "agaB"       = 70 * agaB_dose,   # ~70 for 1 mg/kg
    "agaA"       = 15 * agaA_dose / 0.2, # ~15 for 0.2 mg/kg
    "migalastat" = if (is_amenable) 6.0 else 0,
    "pegu"       = 80 * agaB_dose,   # similar to agaB but longer duration
    "combo"      = 70 * agaB_dose,   # ERT component
    0
  )
  E_total <- E0 + E_max_treat

  # GCS inhibition (lucerastat component)
  GCS_inhib_frac <- if (treatment %in% c("combo")) 0.40 else 0.0

  # Enzyme activity over time (approach to steady-state, weeks)
  t_ss <- 4 * 7  # 4 weeks to reach ~63% of SS enzyme effect
  E_t  <- E_total - (E_total - E0) * exp(-t / t_ss)

  # Lyso-Gb3 dynamics
  LGB3_0 <- if (phenotype == "classic") 12.0 else 4.5   # μg/L
  LGB3_ss_treat <- LGB3_0 * (1 - min(0.85, (E_t[n] - E0) / 40)) *
                   (1 - GCS_inhib_frac * 0.35)
  LGB3_t <- LGB3_0 - (LGB3_0 - pmax(LGB3_ss_treat, 1.5)) *
            (1 - exp(-t / (90 * 7)))  # 90 weeks to reach new SS

  # Gb3 plasma
  Gb3_0 <- if (phenotype == "classic") 0.8 else 0.4
  Gb3_ss <- Gb3_0 * (1 - min(0.75, (E_total - E0) / 50))
  Gb3_t  <- Gb3_0 - (Gb3_0 - pmax(Gb3_ss, 0.05)) * (1 - exp(-t / (78 * 7)))

  # eGFR decline
  # Natural history: ~ -3 to -6 mL/min/1.73m²/yr (classic male)
  eGFR_annual_loss_base <- if (sex == "male" & phenotype == "classic") -4.5
                           else if (sex == "female") -1.5 else -2.0

  eGFR_treat_benefit <- switch(treatment,
    "agaB"  = 0.65,  # 65% reduction in GFR loss (Banikazemi 2007)
    "agaA"  = 0.50,
    "pegu"  = 0.70,
    "migalastat" = if (is_amenable) 0.60 else 0,
    "combo" = 0.75,
    "none"  = 0
  )
  eGFR_annual_loss <- eGFR_annual_loss_base * (1 - eGFR_treat_benefit)
  eGFR_t <- pmax(0, eGFR0 + eGFR_annual_loss * t / 365)

  # LVMi
  LVMi0 <- if (phenotype == "classic" & sex == "male") 148 else 120
  LVMi_rate <- if (phenotype == "classic") 3.0 else 1.5  # g/m²/yr natural
  LVMi_benefit <- switch(treatment,
    "agaB"  = 0.80, "agaA" = 0.65, "pegu" = 0.85,
    "migalastat" = if (is_amenable) 0.60 else 0,
    "combo" = 0.90, "none" = 0
  )
  LVMi_t <- LVMi0 + (LVMi_rate * (1 - LVMi_benefit)) * t / 365

  # Neuropathic Pain (BPI-SF)
  Pain0 <- 6.5
  Pain_benefit <- switch(treatment,
    "agaB"  = 0.35, "agaA" = 0.30, "pegu" = 0.40,
    "migalastat" = if (is_amenable) 0.45 else 0,
    "combo" = 0.45, "none" = 0
  )
  Pain_t <- Pain0 * (1 - Pain_benefit * (1 - exp(-t / (52 * 7))))

  # UPCR
  UPCR0 <- if (phenotype == "classic") 350 else 150
  UPCR_benefit <- switch(treatment,
    "agaB" = 0.40, "agaA" = 0.30, "pegu" = 0.45,
    "migalastat" = if (is_amenable) 0.50 else 0,
    "combo" = 0.55, "none" = 0
  )
  UPCR_t <- UPCR0 * (1 - UPCR_benefit * (1 - exp(-t / (24 * 7))))

  data.frame(
    time     = t,
    Year     = t / 365,
    LysoGb3  = pmax(0.5, LGB3_t),
    Gb3_plasma = pmax(0.02, Gb3_t),
    eGFR     = eGFR_t,
    LVMi     = LVMi_t,
    Pain     = pmax(1, Pain_t),
    UPCR     = pmax(30, UPCR_t),
    EnzymeActivity = E_t,
    Inflammation = 1.8 - 0.9 * (E_t - E0) / (E_total - E0 + 0.1)
  )
}

##############################################################################
# UI
##############################################################################
ui <- dashboardPage(
  skin = "blue",

  dashboardHeader(title = "Fabry Disease QSP Dashboard"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("1. Patient Profile",      tabName = "patient",   icon = icon("user")),
      menuItem("2. PK / Enzyme Activity",     tabName = "pk",        icon = icon("flask")),
      menuItem("3. Gb3 Dynamics",         tabName = "gb3",       icon = icon("dna")),
      menuItem("4. Renal Function",          tabName = "renal",     icon = icon("filter")),
      menuItem("5. Cardiac Function",          tabName = "cardiac",   icon = icon("heart")),
      menuItem("6. Scenario Comparison",      tabName = "scenarios", icon = icon("chart-line")),
      menuItem("7. Biomarker Panel",    tabName = "biomarkers",icon = icon("vial")),
      menuItem("8. Virtual Patient Population",     tabName = "population",icon = icon("users"))
    )
  ),

  dashboardBody(
    tabItems(

      ##########################################################################
      # TAB 1: PATIENT PROFILE
      ##########################################################################
      tabItem("patient",
        fluidRow(
          box(title = "Fabry Disease Patient Settings", width = 4, status = "primary", solidHeader = TRUE,
            selectInput("phenotype", "Phenotype",
                        choices = c("Classic" = "classic",
                                    "Late-onset" = "late")),
            radioButtons("sex", "Sex", choices = c("Male" = "male", "Female" = "female"), inline = TRUE),
            sliderInput("age_start", "Age at diagnosis", min = 5, max = 70, value = 25, step = 1),
            sliderInput("eGFR0", "Baseline eGFR (mL/min/1.73m²)", min = 30, max = 120, value = 90),
            sliderInput("sim_years", "Simulation duration (years)", min = 1, max = 15, value = 5)
          ),
          box(title = "Mutation and Treatment Eligibility", width = 4, status = "warning", solidHeader = TRUE,
            checkboxInput("is_amenable", "Carries a migalastat-amenable mutation (~40%)", value = FALSE),
            helpText("Migalastat (Galafold) is effective only in mutations whose amenability has been\nconfirmed by the HEK293 cell assay."),
            selectInput("treatment", "Select treatment",
              choices = c(
                "No treatment (natural history)" = "none",
                "Agalsidase beta (Fabrazyme)" = "agaB",
                "Agalsidase alfa (Replagal)" = "agaA",
                "Migalastat (Galafold)" = "migalastat",
                "Pegunigalsidase alfa (Elfabrio)" = "pegu",
                "ERT + lucerastat (combination)" = "combo"
              ))
          ),
          box(title = "Dose Adjustment", width = 4, status = "info", solidHeader = TRUE,
            sliderInput("agaB_dose", "Agalsidase beta (mg/kg)", min = 0.3, max = 2.0, value = 1.0, step = 0.1),
            sliderInput("agaA_dose", "Agalsidase alfa (mg/kg)", min = 0.1, max = 0.6, value = 0.2, step = 0.05),
            sliderInput("mig_dose",  "Migalastat (mg, QOD)", min = 50, max = 300, value = 150, step = 25),
            sliderInput("luc_dose",  "Lucerastat (mg, TID)",   min = 250, max = 1500, value = 1000, step = 250)
          )
        ),
        fluidRow(
          box(title = "Fabry Disease Overview — Treatment Target Pathway", width = 12, status = "success",
            p(strong("Fabry disease"), " is an X-linked recessive genetic disorder in which",
              strong("GLA gene mutation"), " causes deficiency of the enzyme α-galactosidase A (α-Gal A)."),
            p("Due to the deficient enzyme, glycosphingolipids (Gb3 and lyso-Gb3)
              progressively accumulate in the kidneys, heart, CNS, skin, and other organs, causing multi-organ damage."),
            tags$ul(
              tags$li(strong("Enzyme replacement therapy (ERT):"), " Agalsidase beta (Fabrazyme) 1 mg/kg Q2W,
                      agalsidase alfa (Replagal) 0.2 mg/kg Q2W,
                      pegunigalsidase alfa (Elfabrio) 1 mg/kg Q4W"),
              tags$li(strong("Chaperone therapy:"), " Migalastat (Galafold) 150 mg QOD —
                      effective only in amenable mutations (~40%), oral administration"),
              tags$li(strong("Substrate reduction therapy (SRT):"), " Lucerastat 1000 mg TID —
                      reduces Gb3 precursor via GCS inhibition (MODIFY trial)")
            )
          )
        )
      ),

      ##########################################################################
      # TAB 2: PK / ENZYME ACTIVITY
      ##########################################################################
      tabItem("pk",
        fluidRow(
          box(title = "α-galactosidase A Enzyme Activity Over Time", width = 6, status = "primary", solidHeader = TRUE,
            plotOutput("plot_enzyme")
          ),
          box(title = "Mechanism of Enzyme Activity Restoration by Treatment", width = 6, status = "info",
            h4("ERT (enzyme replacement therapy)"),
            p("After IV infusion, the enzyme is delivered to the lysosome via the M6P receptor. Agalsidase beta has
              a short plasma half-life of t½ ~45 minutes, while pegunigalsidase alfa (PEGylation)
              extends this to t½ ~80 hours."),
            h4("Migalastat (chaperone)"),
            p("As a small-molecule chaperone, it prevents misfolding of the mutant α-Gal A protein and
              improves lysosomal delivery. It is effective only in amenable mutations (GalaxIES assay)."),
            h4("Lucerastat (SRT)"),
            p("By inhibiting glucosylceramide synthase (GCS/UGCG), it reduces synthesis of the Gb3 precursor.
              The MODIFY trial confirmed a reduction in Gb3.")
          )
        ),
        fluidRow(
          box(title = "Summary of Treatment PK Parameters", width = 12,
            DTOutput("pk_param_table")
          )
        )
      ),

      ##########################################################################
      # TAB 3: Gb3 DYNAMICS
      ##########################################################################
      tabItem("gb3",
        fluidRow(
          box(title = "Plasma Lyso-Gb3 (a sensitive biomarker of treatment response)", width = 6,
              status = "danger", solidHeader = TRUE,
            plotOutput("plot_lysoGb3")
          ),
          box(title = "Plasma Gb3 Over Time", width = 6, status = "warning", solidHeader = TRUE,
            plotOutput("plot_plasma_Gb3")
          )
        ),
        fluidRow(
          box(title = "Clinical Significance of Lyso-Gb3", width = 12,
            p(strong("Plasma Lyso-Gb3"), " is currently the most sensitive biomarker for Fabry disease,
              reflecting baseline disease activity and treatment response."),
            p("Normal reference range: men <2 μg/L, women <2 μg/L (some institutions use <1.4 μg/L)"),
            p("Untreated classic-type men: mean 30–80 μg/L"),
            tags$ul(
              tags$li("Aerts JM et al. (2008) PNAS: discovery of Lyso-Gb3 and its clinical significance"),
              tags$li("Smid BE et al. (2014) Mol Genet Metab: sharp decrease after ERT"),
              tags$li("Nowak A et al. (2017) Kidney Int: correlation with eGFR slope")
            )
          )
        )
      ),

      ##########################################################################
      # TAB 4: RENAL FUNCTION
      ##########################################################################
      tabItem("renal",
        fluidRow(
          box(title = "eGFR Over Time (renal function preservation)", width = 6,
              status = "danger", solidHeader = TRUE,
            plotOutput("plot_eGFR")
          ),
          box(title = "UPCR Over Time (proteinuria)", width = 6, status = "warning", solidHeader = TRUE,
            plotOutput("plot_UPCR")
          )
        ),
        fluidRow(
          box(title = "Staging of Fabry Nephropathy", width = 12,
            p(strong("Natural history:"), " without treatment, eGFR declines by about -3~-12 mL/min/1.73m²
              per year in classic-type men (Warnock 2012)."),
            p(strong("ERT renal protection:"), " Banikazemi et al. (2007) Annals of Internal Medicine
              reduced the incidence of the composite renal event endpoint by 61%."),
            p(strong("Treatment goal:"), " improve the eGFR slope to the -1 mL/min/1.73m²/yr level typical of normal aging."),
            h5("CKD staging criteria (KDIGO 2012):"),
            tags$ul(
              tags$li("G1: ≥90 mL/min/1.73m²"),
              tags$li("G2: 60–89"),
              tags$li("G3a: 45–59"),
              tags$li("G3b: 30–44"),
              tags$li("G4: 15–29"),
              tags$li("G5 (ESRD): <15")
            )
          )
        )
      ),

      ##########################################################################
      # TAB 5: CARDIAC FUNCTION
      ##########################################################################
      tabItem("cardiac",
        fluidRow(
          box(title = "LVMi Over Time (reversal of cardiac hypertrophy)", width = 7,
              status = "primary", solidHeader = TRUE,
            plotOutput("plot_LVMi")
          ),
          box(title = "Neuropathic Pain (BPI-SF)", width = 5,
              status = "warning", solidHeader = TRUE,
            plotOutput("plot_Pain")
          )
        ),
        fluidRow(
          box(title = "Fabry Cardiomyopathy — Key Points", width = 12,
            p(strong("Fabry cardiomyopathy"), " includes left ventricular hypertrophy (LVH), diastolic dysfunction, conduction abnormalities,
              arrhythmia, and myocardial fibrosis (LGE)."),
            p("Upper limit of normal LVMi: men 115 g/m², women 95 g/m²"),
            tags$ul(
              tags$li("Weidemann F et al. (2003) Circulation: LVMi decrease after ERT"),
              tags$li("Desnick RJ (2010) J Inherit Metab Dis: cardiac Gb3 responsiveness to ERT"),
              tags$li("Linhart A (2020) Eur Heart J: natural history of Fabry cardiomyopathy")
            )
          )
        )
      ),

      ##########################################################################
      # TAB 6: SCENARIO COMPARISON
      ##########################################################################
      tabItem("scenarios",
        fluidRow(
          box(title = "Comparison of 6 Treatment Scenarios — Lyso-Gb3", width = 6,
              status = "primary", solidHeader = TRUE,
            plotOutput("plot_scen_lgb3")
          ),
          box(title = "Comparison of 6 Treatment Scenarios — eGFR", width = 6,
              status = "danger", solidHeader = TRUE,
            plotOutput("plot_scen_eGFR")
          )
        ),
        fluidRow(
          box(title = "Scenario Comparison — LVMi", width = 6,
              status = "warning", solidHeader = TRUE,
            plotOutput("plot_scen_LVMi")
          ),
          box(title = "5-Year Projected Outcomes by Treatment Scenario", width = 6,
            DTOutput("outcome_table")
          )
        )
      ),

      ##########################################################################
      # TAB 7: BIOMARKER PANEL
      ##########################################################################
      tabItem("biomarkers",
        fluidRow(
          box(title = "Biomarker Dashboard", width = 12, status = "success",
            fluidRow(
              valueBoxOutput("vb_lysoGb3"),
              valueBoxOutput("vb_eGFR"),
              valueBoxOutput("vb_LVMi")
            ),
            fluidRow(
              valueBoxOutput("vb_UPCR"),
              valueBoxOutput("vb_pain"),
              valueBoxOutput("vb_enzyme")
            )
          )
        ),
        fluidRow(
          box(title = "Biomarker Interpretation Guide", width = 12,
            DTOutput("biomarker_guide_table")
          )
        )
      ),

      ##########################################################################
      # TAB 8: VIRTUAL POPULATION
      ##########################################################################
      tabItem("population",
        fluidRow(
          box(title = "Virtual Patient Population Settings (Monte Carlo)", width = 4,
              status = "info", solidHeader = TRUE,
            sliderInput("n_patients", "Number of simulated patients", min = 20, max = 200, value = 50),
            sliderInput("cv_enzyme", "Enzyme activity CV% (inter-individual variability)", min = 10, max = 80, value = 40),
            sliderInput("cv_eGFR", "Baseline eGFR CV%", min = 10, max = 50, value = 25),
            selectInput("pop_treatment", "Treatment",
              choices = c("none", "agaB", "agaA", "migalastat", "pegu", "combo")),
            actionButton("run_pop", "Run population simulation", class = "btn-primary")
          ),
          box(title = "eGFR Population Distribution (5-year projection)", width = 8, status = "success",
            plotOutput("plot_pop_eGFR")
          )
        ),
        fluidRow(
          box(title = "Lyso-Gb3 Response Distribution", width = 6,
            plotOutput("plot_pop_lgb3")
          ),
          box(title = "Proportion of ERT Responders vs Non-Responders", width = 6,
            plotOutput("plot_responder")
          )
        )
      )
    )
  )
)

##############################################################################
# SERVER
##############################################################################
server <- function(input, output, session) {

  # ── Reactive simulation (single patient) ─────────────────────────────
  sim_data <- reactive({
    simulate_fabry(
      phenotype   = input$phenotype,
      is_amenable = input$is_amenable,
      treatment   = input$treatment,
      agaB_dose   = input$agaB_dose,
      agaA_dose   = input$agaA_dose,
      mig_dose    = input$mig_dose,
      luc_dose    = input$luc_dose,
      sim_years   = input$sim_years,
      age_start   = input$age_start,
      sex         = input$sex,
      eGFR0       = input$eGFR0
    )
  })

  # ── All 6 scenarios for comparison ───────────────────────────────────
  all_scen_data <- reactive({
    treatments <- c("none", "agaB", "agaA", "migalastat", "pegu", "combo")
    labels     <- c("Natural history", "Agalsidase beta", "Agalsidase alfa",
                    "Migalastat", "Pegunigalsidase alfa", "ERT+Lucerastat")
    purrr::map2_dfr(treatments, labels, function(trt, lbl) {
      simulate_fabry(
        phenotype   = input$phenotype,
        is_amenable = TRUE,  # assume amenable for migalastat comparison
        treatment   = trt,
        sim_years   = input$sim_years,
        eGFR0       = input$eGFR0
      ) %>% mutate(Scenario = lbl)
    })
  })

  # ── TAB 2: Enzyme Activity ────────────────────────────────────────────
  output$plot_enzyme <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(Year, EnzymeActivity)) +
      geom_line(color = "#1565C0", size = 1.3) +
      geom_hline(yintercept = 8.0, linetype = "dashed", color = "gray50") +
      annotate("text", x = max(d$Year) * 0.8, y = 8.5,
               label = "Normal range (≥8 nmol/h/mg)", size = 3.5) +
      labs(title = "α-galactosidase A Activity Over Time",
           x = "Time (years)", y = "α-Gal A (nmol/h/mg)") +
      theme_bw(base_size = 13)
  })

  output$pk_param_table <- renderDT({
    data.frame(
      Drug = c("Agalsidase beta (Fabrazyme)", "Agalsidase alfa (Replagal)",
               "Pegunigalsidase alfa (Elfabrio)", "Migalastat (Galafold)", "Lucerastat"),
      Dosing = c("1 mg/kg IV Q2W", "0.2 mg/kg IV Q2W", "1 mg/kg IV Q4W",
               "150 mg PO QOD", "1000 mg PO TID"),
      `t½` = c("~45 min", "~45–110 min", "~80 h", "~3.5 h", "~8 h"),
      `F(%)` = c("100 (IV)", "100 (IV)", "100 (IV)", "75%", "65%"),
      MOA = c("M6P receptor → lysosomal delivery", "M6P receptor → lysosomal delivery",
                   "PEGylation → extended t½", "α-Gal A chaperone (amenable mutations)", "GCS inhibition (SRT)"),
      `EC50/IC50` = c("about 1.5 ng/mL (Km)", "about 1.5 ng/mL (Km)", "about 0.8 ng/mL", "0.25 μg/mL", "0.18 μg/mL"),
      `Emax_Gb3_decrease(%)` = c("~80%", "~70%", "~85%", "~60% (amenable)", "~40%"),
      Trial = c("FABRY-001 (Eng 2001 NEJM)", "Schiffmann 2001 Ann Intern Med",
                   "BRIGHT (Schiffmann 2021 JAMA)", "ATTRACT (Germain 2016 NEJM)",
                   "MODIFY (Lenders 2022 Lancet DE)")
    ) %>% datatable(options = list(pageLength = 5, scrollX = TRUE), rownames = FALSE)
  })

  # ── TAB 3: Gb3 Dynamics ───────────────────────────────────────────────
  output$plot_lysoGb3 <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(Year, LysoGb3)) +
      geom_line(color = "#C62828", size = 1.3) +
      geom_hline(yintercept = 2.0, linetype = "dashed", color = "forestgreen") +
      annotate("text", x = max(d$Year) * 0.75, y = 2.5,
               label = "Upper limit of normal (<2 μg/L)", color = "forestgreen", size = 3.5) +
      labs(title = "Plasma Lyso-Gb3 Over Time",
           x = "Time (years)", y = "Lyso-Gb3 (μg/L)") +
      theme_bw(base_size = 13)
  })

  output$plot_plasma_Gb3 <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(Year, Gb3_plasma)) +
      geom_line(color = "#E65100", size = 1.3) +
      geom_hline(yintercept = 0.1, linetype = "dashed", color = "gray50") +
      labs(title = "Plasma Gb3 (Globotriaosylceramide)",
           x = "Time (years)", y = "Gb3 (μg/mL)") +
      theme_bw(base_size = 13)
  })

  # ── TAB 4: Renal ──────────────────────────────────────────────────────
  output$plot_eGFR <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(Year, eGFR)) +
      geom_line(color = "#AD1457", size = 1.3) +
      geom_hline(yintercept = 60, linetype = "dashed", color = "red") +
      annotate("text", x = max(d$Year) * 0.7, y = 63, label = "CKD G3 threshold (60)",
               color = "red", size = 3.5) +
      labs(title = "eGFR Over Time (CKD-EPI)",
           x = "Time (years)", y = "eGFR (mL/min/1.73m²)") +
      ylim(0, 120) + theme_bw(base_size = 13)
  })

  output$plot_UPCR <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(Year, UPCR)) +
      geom_line(color = "#6A1B9A", size = 1.3) +
      geom_hline(yintercept = 300, linetype = "dashed", color = "orange") +
      labs(title = "UPCR (proteinuria)",
           x = "Time (years)", y = "UPCR (mg/g Cr)") +
      theme_bw(base_size = 13)
  })

  # ── TAB 5: Cardiac ────────────────────────────────────────────────────
  output$plot_LVMi <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(Year, LVMi)) +
      geom_line(color = "#1565C0", size = 1.3) +
      geom_hline(yintercept = 115, linetype = "dashed", color = "blue") +
      annotate("text", x = max(d$Year) * 0.7, y = 118, label = "Upper limit of normal in men (115 g/m²)",
               color = "blue", size = 3.5) +
      labs(title = "LVMi (Left Ventricular Mass Index)",
           x = "Time (years)", y = "LVMi (g/m²)") +
      theme_bw(base_size = 13)
  })

  output$plot_Pain <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(Year, Pain)) +
      geom_line(color = "#FF6F00", size = 1.3) +
      ylim(0, 10) +
      labs(title = "Neuropathic Pain (BPI-SF)",
           x = "Time (years)", y = "Pain score (0–10)") +
      theme_bw(base_size = 13)
  })

  # ── TAB 6: Scenario Comparison ────────────────────────────────────────
  output$plot_scen_lgb3 <- renderPlot({
    d <- all_scen_data()
    ggplot(d, aes(Year, LysoGb3, color = Scenario)) +
      geom_line(size = 1.1) +
      scale_color_manual(values = scenario_colors) +
      geom_hline(yintercept = 2.0, linetype = "dashed") +
      labs(title = "Lyso-Gb3 Scenario Comparison",
           x = "Time (years)", y = "Lyso-Gb3 (μg/L)") +
      theme_bw(base_size = 12) +
      theme(legend.position = "bottom", legend.title = element_blank())
  })

  output$plot_scen_eGFR <- renderPlot({
    d <- all_scen_data()
    ggplot(d, aes(Year, eGFR, color = Scenario)) +
      geom_line(size = 1.1) +
      scale_color_manual(values = scenario_colors) +
      geom_hline(yintercept = 60, linetype = "dashed", color = "red") +
      labs(title = "eGFR Scenario Comparison",
           x = "Time (years)", y = "eGFR (mL/min/1.73m²)") +
      theme_bw(base_size = 12) +
      theme(legend.position = "bottom", legend.title = element_blank())
  })

  output$plot_scen_LVMi <- renderPlot({
    d <- all_scen_data()
    ggplot(d, aes(Year, LVMi, color = Scenario)) +
      geom_line(size = 1.1) +
      scale_color_manual(values = scenario_colors) +
      labs(title = "LVMi Scenario Comparison",
           x = "Time (years)", y = "LVMi (g/m²)") +
      theme_bw(base_size = 12) +
      theme(legend.position = "bottom", legend.title = element_blank())
  })

  output$outcome_table <- renderDT({
    d <- all_scen_data()
    last_yr <- max(d$Year)
    d %>% filter(abs(Year - last_yr) < 0.15) %>%
      group_by(Scenario) %>%
      summarise(
        `Lyso-Gb3 (μg/L)` = round(mean(LysoGb3), 1),
        `eGFR`             = round(mean(eGFR), 1),
        `LVMi (g/m²)`      = round(mean(LVMi), 1),
        `Pain BPI`          = round(mean(Pain), 1),
        .groups = "drop"
      ) %>%
      datatable(rownames = FALSE,
                options = list(pageLength = 8, dom = "t"))
  })

  # ── TAB 7: Biomarker Valuebox ─────────────────────────────────────────
  output$vb_lysoGb3 <- renderValueBox({
    d <- tail(sim_data(), 1)
    valueBox(round(d$LysoGb3, 1), "Lyso-Gb3 (μg/L)", icon = icon("vial"),
             color = if (d$LysoGb3 < 2) "green" else if (d$LysoGb3 < 10) "yellow" else "red")
  })
  output$vb_eGFR <- renderValueBox({
    d <- tail(sim_data(), 1)
    valueBox(round(d$eGFR, 1), "eGFR (mL/min/1.73m²)", icon = icon("filter"),
             color = if (d$eGFR >= 60) "green" else if (d$eGFR >= 30) "yellow" else "red")
  })
  output$vb_LVMi <- renderValueBox({
    d <- tail(sim_data(), 1)
    valueBox(round(d$LVMi, 1), "LVMi (g/m²)", icon = icon("heart"),
             color = if (d$LVMi <= 115) "green" else if (d$LVMi <= 140) "yellow" else "red")
  })
  output$vb_UPCR <- renderValueBox({
    d <- tail(sim_data(), 1)
    valueBox(round(d$UPCR, 0), "UPCR (mg/g)", icon = icon("tint"),
             color = if (d$UPCR < 150) "green" else if (d$UPCR < 500) "yellow" else "red")
  })
  output$vb_pain <- renderValueBox({
    d <- tail(sim_data(), 1)
    valueBox(round(d$Pain, 1), "Pain BPI-SF (0–10)", icon = icon("bolt"),
             color = if (d$Pain < 3) "green" else if (d$Pain < 6) "yellow" else "red")
  })
  output$vb_enzyme <- renderValueBox({
    d <- tail(sim_data(), 1)
    valueBox(round(d$EnzymeActivity, 2), "α-Gal A (nmol/h/mg)", icon = icon("flask"),
             color = if (d$EnzymeActivity >= 5) "green" else if (d$EnzymeActivity >= 1) "yellow" else "red")
  })

  output$biomarker_guide_table <- renderDT({
    data.frame(
      Biomarker = c("Plasma Lyso-Gb3", "Plasma Gb3", "Leukocyte α-Gal A activity",
                     "Urinary Gb3", "eGFR (CKD-EPI)", "LVMi", "UPCR"),
      Unit = c("μg/L", "μg/mL", "nmol/h/mg", "nmol/mg Cr",
               "mL/min/1.73m²", "g/m²", "mg/g Cr"),
      NormalRange = c("<2.0", "<0.1", "≥8 (men)", "Low", "≥60 (G1-2)", "≤115 (men)", "<150"),
      ClinicalSignificance = c("Sensitive treatment response marker (ERT/chaperone)", "Indirect marker of long-term Gb3 accumulation",
                   "Diagnostic criterion (men)", "Direct measurement of renal Gb3 deposition", "Renal function preservation target",
                   "Cardiac hypertrophy reversal target", "Glomerular protein filtration barrier"),
      ResponseSpeed = c("Fast (weeks)", "Intermediate (months)", "Intermediate", "Fast", "Slow (1–2 years)",
                     "Slow (1–2 years)", "Intermediate (6 months)")
    ) %>% datatable(rownames = FALSE, options = list(pageLength = 10, dom = "t"))
  })

  # ── TAB 8: Population ─────────────────────────────────────────────────
  pop_data <- eventReactive(input$run_pop, {
    n <- input$n_patients
    cv_e  <- input$cv_enzyme / 100
    cv_gfr <- input$cv_eGFR / 100
    set.seed(42)

    purrr::map_dfr(seq_len(n), function(i) {
      this_eGFR <- rnorm(1, input$eGFR0, input$eGFR0 * cv_gfr) %>% pmax(20)
      simulate_fabry(
        phenotype   = input$phenotype,
        is_amenable = TRUE,
        treatment   = input$pop_treatment,
        sim_years   = input$sim_years,
        eGFR0       = this_eGFR
      ) %>% mutate(PatientID = i)
    })
  })

  output$plot_pop_eGFR <- renderPlot({
    req(pop_data())
    d <- pop_data()
    ggplot(d, aes(Year, eGFR, group = PatientID)) +
      geom_line(alpha = 0.25, color = "#1565C0") +
      stat_summary(aes(group = 1), fun = median, geom = "line",
                   color = "navy", size = 1.5, linetype = "solid") +
      geom_hline(yintercept = 60, linetype = "dashed", color = "red") +
      labs(title = "Virtual Patient Population eGFR Distribution (median = bold line)",
           x = "Time (years)", y = "eGFR (mL/min/1.73m²)") +
      theme_bw(base_size = 13)
  })

  output$plot_pop_lgb3 <- renderPlot({
    req(pop_data())
    d <- pop_data()
    ggplot(d, aes(Year, LysoGb3, group = PatientID)) +
      geom_line(alpha = 0.25, color = "#C62828") +
      stat_summary(aes(group = 1), fun = median, geom = "line",
                   color = "darkred", size = 1.5) +
      geom_hline(yintercept = 2.0, linetype = "dashed") +
      labs(title = "Virtual Patient Population Lyso-Gb3 Distribution",
           x = "Time (years)", y = "Lyso-Gb3 (μg/L)") +
      theme_bw(base_size = 13)
  })

  output$plot_responder <- renderPlot({
    req(pop_data())
    d <- pop_data()
    resp <- d %>%
      filter(abs(Year - max(Year)) < 0.2) %>%
      group_by(PatientID) %>%
      summarise(LysoGb3_end = mean(LysoGb3)) %>%
      mutate(Response = ifelse(LysoGb3_end < 5.0, "Responder (<5 μg/L)", "Non-responder"))
    ggplot(resp, aes(Response, fill = Response)) +
      geom_bar() +
      scale_fill_manual(values = c("Responder (<5 μg/L)" = "#2E7D32", "Non-responder" = "#C62828")) +
      labs(title = "Proportion of Responders at End of Treatment (Lyso-Gb3 <5 μg/L threshold)",
           x = "", y = "Number of patients") +
      theme_bw(base_size = 13) + theme(legend.position = "none")
  })
}

##############################################################################
# RUN APP
##############################################################################
shinyApp(ui, server)
