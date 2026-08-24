################################################################################
# Cushing's Syndrome QSP Shiny App
# Cushing's Syndrome QSP interactive dashboard
#
# 8 tabs: Patient profile · HPA/PK · Steroidogenesis · Clinical indices ·
#         Scenario comparison · Biomarkers · Metabolic complications · Virtual population
################################################################################

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(tidyverse)
library(plotly)
library(DT)

# ============================================================
# Model Definition (inline)
# ============================================================
cs_code <- '
[PARAM] @annotated
k_CRH_syn=0.06:CRH synthesis
k_CRH_deg=0.20:CRH degradation
k_ACTH_syn=0.12:ACTH synthesis
k_ACTH_deg=0.25:ACTH clearance
CRH_ss=0.30:Normal CRH ss
circ_amp=0.55:Circadian amplitude
circ_peak=8.0:Peak hour
CS_type=1.0:Disease type
tumor_fold=3.5:Tumor ACTH fold
k_F_syn=0.30:Max cortisol synthesis
Km_ACTH=18.0:Km ACTH
k_F_pl_cl=0.15:Cortisol clearance
F_ss=12.0:Normal cortisol ss
IC50_F_ACTH=22.0:IC50 cortisol feedback ACTH
IC50_F_CRH=16.0:IC50 cortisol feedback CRH
nHill_HPA=2.0:Hill coefficient
CBG_tot=700.0:CBG capacity
Kd_CBG=10.0:CBG Kd
GR_tot=100.0:Total GR
kon_GR=0.08:GR-F association
koff_GR=0.40:GR-F dissociation
k_GRnuc_in=0.20:GR nuclear in
k_GRnuc_out=0.10:GR nuclear out
Gluc_base=5.2:Baseline glucose
k_GR_gluc=0.025:GR glucose production
k_gluc_cl=0.12:Glucose clearance
k_ins_sec=0.18:Insulin secretion
k_ins_cl=0.22:Insulin clearance
Ins_base=8.0:Baseline insulin
VAT_base=5.0:Baseline VAT
k_VAT_acc=0.001:VAT accumulation
k_VAT_cl=0.002:VAT clearance
Musc_base=30.0:Baseline muscle
k_musc_loss=0.0008:Muscle loss
k_musc_syn=0.002:Muscle synthesis
BMD_base=0.0:Baseline BMD
k_BMD_loss=0.0003:BMD loss
k_BMD_syn=0.0001:BMD synthesis
BP_base=120.0:Baseline BP
k_BP_F=0.06:Cortisol BP effect
k_BP_cl=0.025:BP clearance
ka_pas=0.22:Pasireotide SC absorption
CL_pas=8.5:Pasireotide CL
V1_pas=28.0:Pasireotide V1
Q_pas=5.5:Pasireotide Q
V2_pas=55.0:Pasireotide V2
F_pas=0.88:Pasireotide F
Emax_pas=0.62:Pasireotide Emax
EC50_pas=0.50:Pasireotide EC50
nH_pas=1.8:Pasireotide Hill
ka_keto=0.85:Ketoconazole absorption
CL_keto=4.8:Ketoconazole CL
V1_keto=38.0:Ketoconazole V1
Emax_keto=0.72:Ketoconazole Emax
EC50_keto=3.2:Ketoconazole EC50
nH_keto=2.2:Ketoconazole Hill
ka_mety=1.10:Metyrapone absorption
CL_mety=9.0:Metyrapone CL
V1_mety=32.0:Metyrapone V1
Emax_mety=0.82:Metyrapone Emax
EC50_mety=2.2:Metyrapone EC50
ka_osilo=0.95:Osilodrostat absorption
CL_osilo=12.0:Osilodrostat CL
V1_osilo=95.0:Osilodrostat V1
Emax_osilo=0.85:Osilodrostat Emax
EC50_osilo=0.15:Osilodrostat EC50
nH_osilo=1.5:Osilodrostat Hill
ka_mife=0.55:Mifepristone absorption
CL_mife=3.2:Mifepristone CL
V1_mife=115.0:Mifepristone V1
Emax_mife=0.82:Mifepristone Emax
EC50_mife=0.45:Mifepristone EC50
UFC_coef=0.012:UFC coefficient

[CMT]
CRH ACTH_PIT ACTH_PL F_ADR F_PL GR_FREE GR_BOUND GR_NUC
GLUCOSE INSULIN VAT MUSCLE BMD BP UFC_ACC
A_PAS_C A_PAS_P A_KETO A_METY A_OSILO A_MIFE

[MAIN]
double t_hr = fmod(TIME, 24.0);
double circ = 1.0 + circ_amp * cos(2.0 * M_PI * (t_hr - circ_peak) / 24.0);
double tumor_ACTH_add = 0.0;
double tumor_F_add = 0.0;
if(CS_type == 1.0) tumor_ACTH_add = (tumor_fold - 1.0) * k_ACTH_syn * CRH_ss;
if(CS_type == 2.0) tumor_ACTH_add = 7.0 * k_ACTH_syn * CRH_ss;
if(CS_type == 3.0) tumor_F_add = 2.5 * k_F_syn;
double E_pas  = (A_PAS_C > 1e-6) ? Emax_pas  * pow(A_PAS_C, nH_pas)   / (pow(EC50_pas, nH_pas)   + pow(A_PAS_C, nH_pas))   : 0.0;
double E_keto = (A_KETO  > 1e-6) ? Emax_keto * pow(A_KETO,  nH_keto)  / (pow(EC50_keto, nH_keto)  + pow(A_KETO,  nH_keto))  : 0.0;
double E_mety = (A_METY  > 1e-6) ? Emax_mety * A_METY  / (EC50_mety  + A_METY)  : 0.0;
double E_osilo= (A_OSILO > 1e-6) ? Emax_osilo * pow(A_OSILO, nH_osilo) / (pow(EC50_osilo, nH_osilo) + pow(A_OSILO, nH_osilo)) : 0.0;
double E_mife = (A_MIFE  > 1e-6) ? Emax_mife * A_MIFE  / (EC50_mife  + A_MIFE)  : 0.0;
double E_CYP11B1 = 1.0 - (1.0 - E_keto) * (1.0 - E_mety) * (1.0 - E_osilo);
double steroid_inh = (E_CYP11B1 > 0.95) ? 0.05 : (1.0 - E_CYP11B1);
double GR_eff = GR_NUC * (1.0 - E_mife);
double FB_CRH  = 1.0 / (1.0 + pow(GR_eff / IC50_F_CRH,  nHill_HPA));
double FB_ACTH = 1.0 / (1.0 + pow(GR_eff / IC50_F_ACTH, nHill_HPA));
double F_syn_rate = k_F_syn * ACTH_PL / (Km_ACTH + ACTH_PL) * steroid_inh + tumor_F_add;

[ODE]
dxdt_CRH     = k_CRH_syn * circ * FB_CRH - k_CRH_deg * CRH;
dxdt_ACTH_PIT = k_ACTH_syn * CRH * FB_ACTH + tumor_ACTH_add - k_ACTH_deg * ACTH_PIT;
dxdt_ACTH_PL  = ACTH_PIT * (1.0 - E_pas) - k_ACTH_deg * ACTH_PL;
dxdt_F_ADR   = F_syn_rate - k_F_pl_cl * F_ADR;
dxdt_F_PL    = k_F_pl_cl * F_ADR - k_F_pl_cl * F_PL;
dxdt_GR_FREE  = -kon_GR * GR_FREE * F_PL + koff_GR * GR_BOUND + k_GRnuc_out * GR_NUC;
dxdt_GR_BOUND =  kon_GR * GR_FREE * F_PL - koff_GR * GR_BOUND - k_GRnuc_in * GR_BOUND;
dxdt_GR_NUC   =  k_GRnuc_in * GR_BOUND - k_GRnuc_out * GR_NUC;
double gluc_prod = k_GR_gluc * GR_eff + 0.02;
double gluc_cl   = k_gluc_cl * (GLUCOSE / Gluc_base);
dxdt_GLUCOSE = gluc_prod - gluc_cl;
dxdt_INSULIN = k_ins_sec * (GLUCOSE - Gluc_base) - k_ins_cl * (INSULIN - Ins_base);
dxdt_VAT     = k_VAT_acc * F_PL - k_VAT_cl * INSULIN * VAT / Ins_base;
dxdt_MUSCLE  = k_musc_syn - k_musc_loss * GR_eff * MUSCLE;
dxdt_BMD     = k_BMD_syn - k_BMD_loss * GR_eff;
dxdt_BP      = k_BP_F * (F_PL - F_ss) - k_BP_cl * (BP - BP_base);
dxdt_UFC_ACC = UFC_coef * F_PL;
dxdt_A_PAS_C = -CL_pas/V1_pas * A_PAS_C - Q_pas/V1_pas * A_PAS_C + Q_pas/V2_pas * A_PAS_P;
dxdt_A_PAS_P =  Q_pas/V1_pas * A_PAS_C - Q_pas/V2_pas * A_PAS_P;
dxdt_A_KETO  = -CL_keto/V1_keto * A_KETO;
dxdt_A_METY  = -CL_mety/V1_mety * A_METY;
dxdt_A_OSILO = -CL_osilo/V1_osilo * A_OSILO;
dxdt_A_MIFE  = -CL_mife/V1_mife * A_MIFE;

[TABLE]
capture cortisol_free   = F_PL;
capture ACTH_pl         = ACTH_PL;
capture glucose_mmol    = GLUCOSE;
capture insulin_uU      = INSULIN;
capture VAT_kg          = VAT;
capture muscle_kg       = MUSCLE;
capture BMD_T           = BMD;
capture BP_sys          = BP;
capture UFC_cumul       = UFC_ACC;
capture GR_nuc_pct      = GR_NUC;
capture E_pas_out       = E_pas;
capture E_keto_out      = E_keto;
capture E_mety_out      = E_mety;
capture E_osilo_out     = E_osilo;
capture E_mife_out      = E_mife;
capture pas_C           = A_PAS_C;
capture keto_C          = A_KETO;
capture mety_C          = A_METY;
capture osilo_C         = A_OSILO;
capture mife_C          = A_MIFE;
capture LNSC_nmol       = F_PL * 27.6 * 0.85;
'

mod_cs <- mcode("cs_shiny", cs_code, quiet = TRUE)

# ============================================================
# Helper: run simulation
# ============================================================
run_simulation <- function(cs_type, tumor_fold, duration_days,
                           drug, dose, freq_h,
                           baseline_cortisol = 25.0,
                           baseline_acth = 55.0,
                           baseline_glucose = 7.8) {
  cmt_map <- list(
    "pasireotide"  = "A_PAS_C",
    "ketoconazole" = "A_KETO",
    "metyrapone"   = "A_METY",
    "osilodrostat" = "A_OSILO",
    "mifepristone" = "A_MIFE"
  )

  init_vals <- list(
    CRH = 0.35, ACTH_PIT = 4.5, ACTH_PL = baseline_acth,
    F_ADR = baseline_cortisol * 1.1, F_PL = baseline_cortisol,
    GR_FREE = 58, GR_BOUND = 22, GR_NUC = 18,
    GLUCOSE = baseline_glucose, INSULIN = 20,
    VAT = 8.0, MUSCLE = 24, BMD = -1.2, BP = 145,
    UFC_ACC = 0,
    A_PAS_C = 0, A_PAS_P = 0, A_KETO = 0,
    A_METY = 0, A_OSILO = 0, A_MIFE = 0
  )

  m <- mod_cs %>%
    param(CS_type = cs_type, tumor_fold = tumor_fold)

  addl_doses <- ceiling(duration_days * 24 / freq_h) - 1

  if (!is.null(drug) && drug != "none" && dose > 0) {
    cmt <- cmt_map[[drug]]
    e <- ev(amt = dose, cmt = cmt, ii = freq_h, addl = addl_doses)
    out <- m %>% init(init_vals) %>%
      ev(e) %>%
      mrgsim(end = duration_days * 24, delta = 2.0)
  } else {
    out <- m %>% init(init_vals) %>%
      mrgsim(end = duration_days * 24, delta = 2.0)
  }
  as_tibble(out)
}

# ============================================================
# UI
# ============================================================
ui <- dashboardPage(
  skin = "red",
  dashboardHeader(title = "Cushing's Syndrome QSP Dashboard"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient profile",       tabName = "patient",   icon = icon("user-md")),
      menuItem("HPA axis / PK dynamics", tabName = "hpa_pk",    icon = icon("chart-line")),
      menuItem("Steroidogenesis",       tabName = "steroid",   icon = icon("flask")),
      menuItem("Clinical indices",      tabName = "clinical",  icon = icon("stethoscope")),
      menuItem("Scenario comparison",   tabName = "scenario",  icon = icon("sliders-h")),
      menuItem("Biomarker panel",       tabName = "biomarker", icon = icon("vial")),
      menuItem("Metabolic complications", tabName = "metabolic", icon = icon("heartbeat")),
      menuItem("Virtual population",    tabName = "virtual",   icon = icon("users"))
    ),
    hr(),
    # Global Controls
    h5("Global settings", style = "padding-left:15px; color:#ccc"),
    selectInput("cs_type", "Cushing's syndrome subtype",
                choices = c("Cushing's disease (pituitary adenoma)" = "1",
                            "Ectopic ACTH syndrome"                = "2",
                            "Adrenal adenoma (ACTH-independent)"   = "3"),
                selected = "1"),
    sliderInput("tumor_fold", "Tumour ACTH fold (Cushing's disease / ectopic)",
                min = 1.5, max = 12.0, value = 3.5, step = 0.5),
    sliderInput("duration", "Simulation duration (days)",
                min = 7, max = 365, value = 30, step = 7),
    hr(),
    h5("Drug selection", style = "padding-left:15px; color:#ccc"),
    selectInput("drug", "Treatment drug",
                choices = c("No treatment" = "none",
                            "Pasireotide"  = "pasireotide",
                            "Ketoconazole" = "ketoconazole",
                            "Metyrapone"   = "metyrapone",
                            "Osilodrostat" = "osilodrostat",
                            "Mifepristone" = "mifepristone"),
                selected = "none"),
    conditionalPanel(
      "input.drug != 'none'",
      numericInput("dose", "Dose", value = 5),
      selectInput("freq", "Dosing interval",
                  choices = c("Once daily (QD)" = "24",
                              "Twice daily (BID)" = "12",
                              "Three times daily (TID)" = "8"),
                  selected = "12")
    )
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .main-header .logo { font-size: 14px; }
      .box-header { font-weight: bold; }
      .info-box { min-height: 80px; }
    "))),
    tabItems(

      # ===================================================
      # TAB 1: Patient profile
      # ===================================================
      tabItem(tabName = "patient",
        fluidRow(
          box(title = "Patient baseline characteristics", status = "danger", solidHeader = TRUE, width = 4,
              sliderInput("pt_cortisol", "Baseline plasma cortisol (μg/dL)", 15, 60, 25, 1),
              sliderInput("pt_acth", "Baseline plasma ACTH (pg/mL)", 20, 250, 55, 5),
              sliderInput("pt_glucose", "Baseline blood glucose (mmol/L)", 4.5, 15.0, 7.8, 0.1),
              selectInput("pt_sex", "Sex", choices = c("Female" = "F", "Male" = "M")),
              sliderInput("pt_age", "Age (years)", 18, 75, 42, 1),
              actionButton("run_sim", "Run simulation",
                           class = "btn btn-danger btn-block")
          ),
          box(title = "Cushing's syndrome diagnostic criteria", status = "warning", solidHeader = TRUE, width = 4,
              tableOutput("diag_table")
          ),
          box(title = "Patient profile summary", status = "primary", solidHeader = TRUE, width = 4,
              uiOutput("pt_summary")
          )
        ),
        fluidRow(
          infoBoxOutput("box_cortisol"), infoBoxOutput("box_acth"),
          infoBoxOutput("box_ufc"),      infoBoxOutput("box_lnsc")
        ),
        fluidRow(
          box(title = "Pathophysiology by aetiology", status = "info", solidHeader = TRUE, width = 12,
              tableOutput("etiology_table"))
        )
      ),

      # ===================================================
      # TAB 2: HPA axis / PK dynamics
      # ===================================================
      tabItem(tabName = "hpa_pk",
        fluidRow(
          box(title = "Circadian CRH and cortisol rhythm (24h)", status = "primary",
              solidHeader = TRUE, width = 6,
              plotlyOutput("plot_circadian", height = 350)),
          box(title = "ACTH dynamics", status = "danger",
              solidHeader = TRUE, width = 6,
              plotlyOutput("plot_acth", height = 350))
        ),
        fluidRow(
          box(title = "Plasma cortisol time course", status = "warning",
              solidHeader = TRUE, width = 6,
              plotlyOutput("plot_cortisol", height = 350)),
          box(title = "Drug PK profile", status = "success",
              solidHeader = TRUE, width = 6,
              plotlyOutput("plot_drug_pk", height = 350))
        ),
        fluidRow(
          box(title = "HPA axis feedback strength", status = "info", solidHeader = TRUE, width = 12,
              plotlyOutput("plot_feedback", height = 280))
        )
      ),

      # ===================================================
      # TAB 3: Steroidogenesis
      # ===================================================
      tabItem(tabName = "steroid",
        fluidRow(
          box(title = "Steroidogenic pathway enzyme inhibition", status = "warning",
              solidHeader = TRUE, width = 6,
              plotlyOutput("plot_enzyme_inh", height = 350)),
          box(title = "Cortisol synthesis inhibition dynamics", status = "danger",
              solidHeader = TRUE, width = 6,
              plotlyOutput("plot_synth_inh", height = 350))
        ),
        fluidRow(
          box(title = "Overview of the main steroidogenic pathway", status = "info",
              solidHeader = TRUE, width = 12,
              tableOutput("steroid_pathway_table"))
        )
      ),

      # ===================================================
      # TAB 4: Clinical indices
      # ===================================================
      tabItem(tabName = "clinical",
        fluidRow(
          box(title = "24h UFC (urinary free cortisol)", status = "danger",
              solidHeader = TRUE, width = 6,
              plotlyOutput("plot_ufc", height = 300)),
          box(title = "Late-night salivary cortisol (LNSC)", status = "warning",
              solidHeader = TRUE, width = 6,
              plotlyOutput("plot_lnsc", height = 300))
        ),
        fluidRow(
          box(title = "Time at which remission criteria are met", status = "success",
              solidHeader = TRUE, width = 6,
              verbatimTextOutput("remission_check")),
          box(title = "Clinical diagnostic algorithm", status = "info",
              solidHeader = TRUE, width = 6,
              tableOutput("diag_algorithm"))
        )
      ),

      # ===================================================
      # TAB 5: Scenario comparison
      # ===================================================
      tabItem(tabName = "scenario",
        fluidRow(
          box(title = "Setup for comparing 5 treatment scenarios", status = "primary",
              solidHeader = TRUE, width = 3,
              checkboxGroupInput("scen_drugs", "Treatments to compare",
                choices = c("No treatment"           = "none",
                            "Pasireotide 0.6mg BID"  = "pasireotide",
                            "Ketoconazole 400mg BID" = "ketoconazole",
                            "Osilodrostat 5mg BID"   = "osilodrostat",
                            "Mifepristone 600mg QD"  = "mifepristone"),
                selected = c("none", "ketoconazole", "osilodrostat")),
              actionButton("run_compare", "Run comparison", class = "btn btn-primary btn-block")
          ),
          box(title = "Cortisol scenario comparison", status = "danger",
              solidHeader = TRUE, width = 9,
              plotlyOutput("plot_scen_cortisol", height = 350))
        ),
        fluidRow(
          box(title = "ACTH scenario comparison", status = "warning",
              solidHeader = TRUE, width = 6,
              plotlyOutput("plot_scen_acth", height = 280)),
          box(title = "Blood glucose scenario comparison", status = "success",
              solidHeader = TRUE, width = 6,
              plotlyOutput("plot_scen_glucose", height = 280))
        ),
        fluidRow(
          box(title = "Key values by scenario (final time point)", status = "info",
              solidHeader = TRUE, width = 12,
              DTOutput("scen_summary_table"))
        )
      ),

      # ===================================================
      # TAB 6: Biomarkers
      # ===================================================
      tabItem(tabName = "biomarker",
        fluidRow(
          box(title = "GR nuclear activation dynamics", status = "primary",
              solidHeader = TRUE, width = 6,
              plotlyOutput("plot_GR_nuc", height = 300)),
          box(title = "ACTH / cortisol negative feedback efficiency", status = "danger",
              solidHeader = TRUE, width = 6,
              plotlyOutput("plot_feedback_loop", height = 300))
        ),
        fluidRow(
          box(title = "Biomarker reference value panel", status = "warning",
              solidHeader = TRUE, width = 6,
              DTOutput("biomarker_table")),
          box(title = "Dexamethasone suppression test simulation", status = "success",
              solidHeader = TRUE, width = 6,
              plotlyOutput("plot_dex_test", height = 300))
        )
      ),

      # ===================================================
      # TAB 7: Metabolic complications
      # ===================================================
      tabItem(tabName = "metabolic",
        fluidRow(
          box(title = "Blood glucose and insulin resistance", status = "danger",
              solidHeader = TRUE, width = 6,
              plotlyOutput("plot_glucose_ins", height = 300)),
          box(title = "Visceral fat and skeletal muscle change", status = "warning",
              solidHeader = TRUE, width = 6,
              plotlyOutput("plot_body_comp", height = 300))
        ),
        fluidRow(
          box(title = "Bone mineral density change (BMD T-score)", status = "primary",
              solidHeader = TRUE, width = 6,
              plotlyOutput("plot_BMD", height = 300)),
          box(title = "Blood pressure change", status = "success",
              solidHeader = TRUE, width = 6,
              plotlyOutput("plot_BP", height = 300))
        ),
        fluidRow(
          box(title = "Metabolic complication risk score", status = "info",
              solidHeader = TRUE, width = 12,
              uiOutput("metabolic_risk"))
        )
      ),

      # ===================================================
      # TAB 8: Virtual population analysis
      # ===================================================
      tabItem(tabName = "virtual",
        fluidRow(
          box(title = "Population simulation settings", status = "primary",
              solidHeader = TRUE, width = 3,
              numericInput("n_virtual", "Number of virtual patients", 100, 10, 500, 10),
              sliderInput("var_tumor", "Tumour ACTH fold variability (CV%)", 5, 50, 20, 5),
              sliderInput("var_enz", "Enzyme activity variability (CV%)", 5, 40, 15, 5),
              selectInput("virt_drug", "Treatment simulation",
                choices = c("No treatment" = "none",
                            "Ketoconazole 400mg BID" = "ketoconazole",
                            "Osilodrostat 5mg BID" = "osilodrostat"),
                selected = "osilodrostat"),
              actionButton("run_virtual", "Population simulation", class = "btn btn-primary btn-block")
          ),
          box(title = "UFC distribution (simulated population)", status = "danger",
              solidHeader = TRUE, width = 9,
              plotlyOutput("plot_virtual_ufc", height = 350))
        ),
        fluidRow(
          box(title = "Treatment response rate (UFC normalisation criterion)", status = "success",
              solidHeader = TRUE, width = 6,
              plotlyOutput("plot_virtual_response", height = 280)),
          box(title = "Population cortisol distribution (Day 30)", status = "warning",
              solidHeader = TRUE, width = 6,
              plotlyOutput("plot_virtual_dist", height = 280))
        )
      )

    ) # end tabItems
  )
)

# ============================================================
# SERVER
# ============================================================
server <- function(input, output, session) {

  # Reactive simulation result
  sim_result <- eventReactive(input$run_sim, {
    run_simulation(
      cs_type        = as.numeric(input$cs_type),
      tumor_fold     = input$tumor_fold,
      duration_days  = input$duration,
      drug           = if(input$drug == "none") NULL else input$drug,
      dose           = if(input$drug == "none") 0 else input$dose,
      freq_h         = as.numeric(input$freq),
      baseline_cortisol = input$pt_cortisol,
      baseline_acth  = input$pt_acth,
      baseline_glucose = input$pt_glucose
    )
  }, ignoreNULL = FALSE)

  # Auto-trigger on load
  observe({ input$run_sim })

  # ===========================
  # TAB 1: Patient Profile
  # ===========================
  output$diag_table <- renderTable({
    tibble(
      Test = c("24h UFC", "LNSC", "1mg DST", "Plasma ACTH", "High-dose DST"),
      `Normal range` = c("<50 μg/24h", "<4 nmol/L", "<1.8 μg/dL", "10-46 pg/mL", ">50% suppression"),
      `Cushing's syndrome` = c(">150", ">10", ">1.8", "Variable", "Not suppressed")
    )
  })

  output$etiology_table <- renderTable({
    tibble(
      Type = c("Cushing's disease", "Ectopic ACTH", "Adrenal adenoma", "Adrenal carcinoma"),
      Frequency = c("70%", "10%", "15%", "5%"),
      ACTH = c("↑ (mild-moderate)", "↑↑↑ (marked)", "↓ (suppressed)", "↓ (suppressed)"),
      Features = c("Pituitary adenoma, USP8 mutation ~50%", "Small-cell lung cancer, carcinoid",
               "Autonomous cortisol secretion", "Cortisol + androgen excess"),
      Treatment = c("Transsphenoidal surgery", "Treatment of the primary tumour", "Adrenalectomy", "Adrenalectomy + mitotane")
    )
  }, striped = TRUE, hover = TRUE)

  output$box_cortisol <- renderInfoBox({
    infoBox("Plasma cortisol", paste0(input$pt_cortisol, " μg/dL"),
            color = if(input$pt_cortisol > 23) "red" else "green",
            icon = icon("tint"), fill = TRUE)
  })
  output$box_acth <- renderInfoBox({
    infoBox("Plasma ACTH", paste0(input$pt_acth, " pg/mL"),
            color = if(input$pt_acth > 46) "red" else "green",
            icon = icon("flask"), fill = TRUE)
  })
  output$box_ufc <- renderInfoBox({
    est_ufc <- input$pt_cortisol * 0.012 * 24
    infoBox("Estimated UFC", paste0(round(est_ufc, 0), " μg/24h"),
            color = if(est_ufc > 50) "red" else "green",
            icon = icon("vial"), fill = TRUE)
  })
  output$box_lnsc <- renderInfoBox({
    lnsc <- input$pt_cortisol * 27.6 * 0.85
    infoBox("Estimated LNSC", paste0(round(lnsc, 1), " nmol/L"),
            color = if(lnsc > 4) "red" else "green",
            icon = icon("moon"), fill = TRUE)
  })

  output$pt_summary <- renderUI({
    cs_name <- switch(input$cs_type,
      "1" = "Cushing's disease (pituitary ACTH adenoma)",
      "2" = "Ectopic ACTH syndrome",
      "3" = "ACTH-independent adrenal adenoma"
    )
    tagList(
      h4(cs_name), hr(),
      tags$b("Sex: "), input$pt_sex, br(),
      tags$b("Age: "), input$pt_age, "years", br(),
      tags$b("Baseline cortisol: "), input$pt_cortisol, "μg/dL", br(),
      tags$b("Baseline ACTH: "), input$pt_acth, "pg/mL", br(),
      tags$b("Baseline blood glucose: "), input$pt_glucose, "mmol/L", br(),
      tags$b("Tumour fold: "), input$tumor_fold, "×"
    )
  })

  # ===========================
  # TAB 2: HPA/PK
  # ===========================
  output$plot_circadian <- renderPlotly({
    df <- sim_result() %>% filter(time <= 48)
    p <- plot_ly(df, x = ~time, y = ~cortisol_free, type = "scatter",
                 mode = "lines", name = "Free cortisol (μg/dL)",
                 line = list(color = "#E53935", width = 2)) %>%
      add_trace(y = ~LNSC_nmol / 10, name = "LNSC/10 (nmol/L)",
                line = list(color = "#1E88E5", width = 2, dash = "dash")) %>%
      layout(title = "Circadian rhythm (48 hours)",
             xaxis = list(title = "Time (h)"),
             yaxis = list(title = "Cortisol (μg/dL)"))
    p
  })

  output$plot_acth <- renderPlotly({
    df <- sim_result()
    plot_ly(df, x = ~time / 24, y = ~ACTH_pl, type = "scatter", mode = "lines",
            line = list(color = "#E65100", width = 2)) %>%
      add_segments(x = 0, xend = max(df$time) / 24, y = 46, yend = 46,
                   line = list(color = "gray", dash = "dash")) %>%
      layout(title = "Plasma ACTH dynamics",
             xaxis = list(title = "Time (days)"),
             yaxis = list(title = "ACTH (pg/mL)"))
  })

  output$plot_cortisol <- renderPlotly({
    df <- sim_result()
    plot_ly(df, x = ~time / 24, y = ~cortisol_free, type = "scatter", mode = "lines",
            name = "Free cortisol", line = list(color = "#F9A825", width = 2)) %>%
      layout(title = "Plasma free cortisol",
             xaxis = list(title = "Time (days)"),
             yaxis = list(title = "Cortisol (μg/dL)"))
  })

  output$plot_drug_pk <- renderPlotly({
    df <- sim_result() %>% filter(time <= 120)
    p <- plot_ly(df, x = ~time / 24) %>%
      add_trace(y = ~pas_C, name = "Pasireotide (ng/mL)", type = "scatter",
                mode = "lines", line = list(color = "#1E88E5")) %>%
      add_trace(y = ~keto_C, name = "Ketoconazole (μg/mL)", type = "scatter",
                mode = "lines", line = list(color = "#43A047")) %>%
      add_trace(y = ~osilo_C * 10, name = "Osilodrostat×10 (μg/mL)", type = "scatter",
                mode = "lines", line = list(color = "#8E24AA")) %>%
      add_trace(y = ~mife_C, name = "Mifepristone (μg/mL)", type = "scatter",
                mode = "lines", line = list(color = "#F4511E")) %>%
      layout(title = "Drug PK profile (5 days)",
             xaxis = list(title = "Time (days)"),
             yaxis = list(title = "Plasma concentration"))
    p
  })

  output$plot_feedback <- renderPlotly({
    df <- sim_result()
    fb_acth <- 1 / (1 + (df$GR_nuc_pct / 22)^2)
    fb_crh  <- 1 / (1 + (df$GR_nuc_pct / 16)^2)
    plot_ly(df, x = ~time / 24) %>%
      add_trace(y = ~fb_acth, name = "ACTH feedback strength",
                type = "scatter", mode = "lines", line = list(color = "#E53935")) %>%
      add_trace(y = ~fb_crh, name = "CRH feedback strength",
                type = "scatter", mode = "lines", line = list(color = "#1E88E5")) %>%
      layout(title = "HPA negative feedback efficiency (1 = maximal suppression)",
             xaxis = list(title = "Time (days)"),
             yaxis = list(title = "Feedback strength (0-1)"))
  })

  # ===========================
  # TAB 3: Steroidogenesis
  # ===========================
  output$plot_enzyme_inh <- renderPlotly({
    df <- sim_result()
    plot_ly(df, x = ~time / 24) %>%
      add_trace(y = ~E_keto_out,  name = "Ketoconazole (CYP17A1+11B1)",
                type = "scatter", mode = "lines", line = list(color = "#43A047")) %>%
      add_trace(y = ~E_mety_out,  name = "Metyrapone (CYP11B1)",
                type = "scatter", mode = "lines", line = list(color = "#FFA000")) %>%
      add_trace(y = ~E_osilo_out, name = "Osilodrostat (CYP11B1/B2)",
                type = "scatter", mode = "lines", line = list(color = "#8E24AA")) %>%
      layout(title = "Enzyme inhibition effect",
             xaxis = list(title = "Time (days)"),
             yaxis = list(title = "Inhibition (0-1)", range = c(0, 1)))
  })

  output$plot_synth_inh <- renderPlotly({
    df <- sim_result()
    plot_ly(df, x = ~time / 24, y = ~cortisol_free,
            type = "scatter", mode = "lines",
            line = list(color = "#E53935", width = 2)) %>%
      add_segments(x = 0, xend = max(df$time)/24, y = 23, yend = 23,
                   line = list(color = "gray", dash = "dot"), name = "Upper limit of normal") %>%
      layout(title = "Cortisol synthesis inhibition dynamics",
             xaxis = list(title = "Time (days)"),
             yaxis = list(title = "Free cortisol (μg/dL)"))
  })

  output$steroid_pathway_table <- renderTable({
    tibble(
      Enzyme = c("StAR", "CYP11A1 (P450scc)", "CYP17A1 (17α-OHase)",
               "HSD3B2", "CYP21A2", "CYP11B1 (11β-OHase)", "CYP11B2 (Aldosterone synthase)"),
      Location = c("Inner mitochondrial membrane", "Mitochondria", "Endoplasmic reticulum", "Endoplasmic reticulum", "Endoplasmic reticulum", "Mitochondria", "Mitochondria"),
      Substrate = c("Cholesterol transport", "Cholesterol", "Pregnenolone/progesterone", "Pregnenolone", "17-OHP",
               "11-deoxycortisol", "Corticosterone"),
      Product = c("(transport facilitation)", "Pregnenolone", "17-OHP", "Progesterone", "11-deoxycortisol",
               "Cortisol", "Aldosterone"),
      Inhibitor = c("-", "Mitotane", "Ketoconazole, LCI699", "Etomidate", "-",
                 "Metyrapone, osilodrostat, ketoconazole", "Osilodrostat")
    )
  }, striped = TRUE, hover = TRUE)

  # ===========================
  # TAB 4: Clinical Endpoints
  # ===========================
  output$plot_ufc <- renderPlotly({
    df <- sim_result()
    plot_ly(df, x = ~time / 24, y = ~UFC_cumul,
            type = "scatter", mode = "lines",
            line = list(color = "#E53935", width = 2)) %>%
      add_segments(x = 0, xend = max(df$time)/24, y = 50, yend = 50,
                   line = list(color = "orange", dash = "dash"), name = "Upper limit of normal (50 μg/24h)") %>%
      layout(title = "UFC (cumulative urinary free cortisol)",
             xaxis = list(title = "Time (days)"),
             yaxis = list(title = "UFC (μg/24h)"))
  })

  output$plot_lnsc <- renderPlotly({
    df <- sim_result() %>%
      filter(time %% 24 >= 23 | time %% 24 < 1)  # midnight samples
    plot_ly(df, x = ~time / 24, y = ~LNSC_nmol,
            type = "scatter", mode = "markers+lines",
            line = list(color = "#1565C0", width = 1.5)) %>%
      add_segments(x = 0, xend = max(df$time)/24, y = 4, yend = 4,
                   line = list(color = "orange", dash = "dash")) %>%
      layout(title = "Late-night salivary cortisol (LNSC)",
             xaxis = list(title = "Time (days)"),
             yaxis = list(title = "LNSC (nmol/L)"))
  })

  output$remission_check <- renderText({
    df <- sim_result()
    last <- tail(df, 1)
    ufc_ok  <- last$UFC_cumul < 50
    lnsc_ok <- last$LNSC_nmol < 4
    cortisol_ok <- last$cortisol_free < 23
    cat_lines <- c(
      paste("UFC:", ifelse(ufc_ok, "✅ Normalisation achieved", "❌ Not achieved"), "(< 50 μg/24h)"),
      paste("LNSC:", ifelse(lnsc_ok, "✅ Normalisation achieved", "❌ Not achieved"), "(< 4 nmol/L)"),
      paste("Cortisol:", ifelse(cortisol_ok, "✅ Normalisation achieved", "❌ Not achieved"), "(< 23 μg/dL)"),
      "",
      if(ufc_ok && lnsc_ok && cortisol_ok) "→ Biochemical remission achieved" else "→ Further treatment required"
    )
    paste(cat_lines, collapse = "\n")
  })

  output$diag_algorithm <- renderTable({
    tibble(
      Step = c("1st: screening", "2nd: confirmation", "3rd: aetiology", "4th: localisation"),
      Test = c("24h UFC, LNSC measurement",
               "1mg DST, late-night serum cortisol",
               "Plasma ACTH, CRH stimulation test, 8mg DST",
               "Pituitary MRI, IPSS, abdominal CT/MRI"),
      Meaning = c("Confirms biochemical hypercortisolism",
               "Confirms loss of circadian rhythm and failure to suppress",
               "Distinguishes ACTH-dependent from ACTH-independent",
               "Localises the lesion")
    )
  }, striped = TRUE)

  # ===========================
  # TAB 5: Scenario Comparison
  # ===========================
  scen_data <- eventReactive(input$run_compare, {
    req(input$scen_drugs)
    drug_config <- list(
      "none"         = list(drug = NULL, dose = 0, freq = 12, label = "No treatment"),
      "pasireotide"  = list(drug = "pasireotide", dose = 600, freq = 12, label = "Pasireotide 0.6mg BID"),
      "ketoconazole" = list(drug = "ketoconazole", dose = 400, freq = 12, label = "Ketoconazole 400mg BID"),
      "osilodrostat" = list(drug = "osilodrostat", dose = 5, freq = 12, label = "Osilodrostat 5mg BID"),
      "mifepristone" = list(drug = "mifepristone", dose = 600, freq = 24, label = "Mifepristone 600mg QD")
    )
    bind_rows(lapply(input$scen_drugs, function(d) {
      cfg <- drug_config[[d]]
      tryCatch(
        run_simulation(as.numeric(input$cs_type), input$tumor_fold,
                       input$duration, cfg$drug, cfg$dose, cfg$freq,
                       input$pt_cortisol, input$pt_acth, input$pt_glucose) %>%
          mutate(scenario = cfg$label),
        error = function(e) NULL
      )
    }))
  }, ignoreNULL = FALSE)

  output$plot_scen_cortisol <- renderPlotly({
    df <- scen_data()
    req(nrow(df) > 0)
    p <- plot_ly()
    for (s in unique(df$scenario)) {
      d <- df %>% filter(scenario == s)
      p <- p %>% add_trace(data = d, x = ~time/24, y = ~cortisol_free,
                           type = "scatter", mode = "lines", name = s)
    }
    p %>% add_segments(x = 0, xend = max(df$time)/24, y = 23, yend = 23,
                       line = list(color = "gray", dash = "dot"), name = "Upper limit of normal") %>%
      layout(title = "Cortisol comparison by scenario",
             xaxis = list(title = "Time (days)"),
             yaxis = list(title = "Cortisol (μg/dL)"))
  })

  output$plot_scen_acth <- renderPlotly({
    df <- scen_data()
    req(nrow(df) > 0)
    p <- plot_ly()
    for (s in unique(df$scenario)) {
      d <- df %>% filter(scenario == s)
      p <- p %>% add_trace(data = d, x = ~time/24, y = ~ACTH_pl,
                           type = "scatter", mode = "lines", name = s)
    }
    p %>% layout(title = "ACTH scenario comparison",
                 xaxis = list(title = "Time (days)"),
                 yaxis = list(title = "ACTH (pg/mL)"))
  })

  output$plot_scen_glucose <- renderPlotly({
    df <- scen_data()
    req(nrow(df) > 0)
    p <- plot_ly()
    for (s in unique(df$scenario)) {
      d <- df %>% filter(scenario == s)
      p <- p %>% add_trace(data = d, x = ~time/24, y = ~glucose_mmol,
                           type = "scatter", mode = "lines", name = s)
    }
    p %>% layout(title = "Blood glucose scenario comparison",
                 xaxis = list(title = "Time (days)"),
                 yaxis = list(title = "Blood glucose (mmol/L)"))
  })

  output$scen_summary_table <- renderDT({
    df <- scen_data()
    req(nrow(df) > 0)
    df %>%
      group_by(scenario) %>%
      filter(time == max(time)) %>%
      summarise(
        `Cortisol (μg/dL)`  = round(mean(cortisol_free), 1),
        `ACTH (pg/mL)`      = round(mean(ACTH_pl), 1),
        `Blood glucose (mmol/L)` = round(mean(glucose_mmol), 1),
        `Blood pressure (mmHg)`  = round(mean(BP_sys), 0),
        `Muscle (kg)`       = round(mean(muscle_kg), 1),
        `BMD T-score`       = round(mean(BMD_T), 2),
        .groups = "drop"
      ) %>%
      datatable(rownames = FALSE, options = list(pageLength = 10))
  })

  # ===========================
  # TAB 6: Biomarker
  # ===========================
  output$plot_GR_nuc <- renderPlotly({
    df <- sim_result()
    plot_ly(df, x = ~time/24, y = ~GR_nuc_pct, type = "scatter", mode = "lines",
            line = list(color = "#3949AB", width = 2)) %>%
      layout(title = "GR nuclear activation (% total GR)",
             xaxis = list(title = "Time (days)"),
             yaxis = list(title = "Nuclear GR (%)"))
  })

  output$plot_feedback_loop <- renderPlotly({
    df <- sim_result()
    plot_ly(df, x = ~cortisol_free, y = ~ACTH_pl,
            type = "scatter", mode = "markers",
            marker = list(color = ~time, colorscale = "RdYlBu", size = 4,
                          colorbar = list(title = "Time (h)"))) %>%
      layout(title = "Cortisol-ACTH feedback loop",
             xaxis = list(title = "Cortisol (μg/dL)"),
             yaxis = list(title = "ACTH (pg/mL)"))
  })

  output$biomarker_table <- renderDT({
    tibble(
      Biomarker = c("UFC 24h", "LNSC", "Morning serum cortisol", "Plasma ACTH",
                    "1mg DST cortisol", "Fasting glucose", "BMD T-score",
                    "Systolic blood pressure", "HDRS-17", "CRP"),
      `Normal range` = c("<50 μg/24h", "<4 nmol/L", "6-23 μg/dL", "10-46 pg/mL",
                 "<1.8 μg/dL", "<5.6 mmol/L", ">-1.0", "<140 mmHg", "<7 points", "<1 mg/L"),
      `Cushing's syndrome` = c(">150", ">10", ">23 (or cycling)", "↑ (pit/ectopic)",
                    ">1.8 (unsuppressed)", "5.6-11.1+", "<-1.0 (often below -2)", ">140",
                    "7-17 points (mild-moderate)", "1-3 (raised in some)"),
      `Clinical meaning` = c("Primary index of hypercortisolism", "Most sensitive screening test", "Diagnostic when the circadian rhythm is lost",
                  "Dependent vs independent", "Confirms loss of feedback", "Insulin resistance", "Osteoporosis risk",
                  "Cardiovascular risk", "Psychiatric complications", "Systemic inflammation")
    ) %>%
      datatable(rownames = FALSE)
  })

  output$plot_dex_test <- renderPlotly({
    # Simulate 1mg dexamethasone suppression test (administered at 23h, measured at 8h)
    dex_times <- c(0, 23, 31)
    dex_cortisol_normal <- c(15, 14.5, 0.8)
    dex_cortisol_cd     <- c(25, 24,   18.5)
    dex_cortisol_ectopic <- c(35, 34.5, 32.0)
    df_dex <- tibble(
      time = rep(dex_times, 3),
      cortisol = c(dex_cortisol_normal, dex_cortisol_cd, dex_cortisol_ectopic),
      group = rep(c("Normal", "Cushing's disease", "Ectopic ACTH"), each = 3)
    )
    plot_ly(df_dex, x = ~time, y = ~cortisol, color = ~group,
            type = "scatter", mode = "lines+markers") %>%
      add_segments(x = 23, xend = 23, y = 0, yend = 40,
                   line = list(color = "gray", dash = "dot"), name = "Dexamethasone administration") %>%
      add_segments(x = 0, xend = 31, y = 1.8, yend = 1.8,
                   line = list(color = "orange", dash = "dash"), name = "Cut-off 1.8 μg/dL") %>%
      layout(title = "1mg dexamethasone suppression test simulation",
             xaxis = list(title = "Time (h)"),
             yaxis = list(title = "Cortisol (μg/dL)"))
  })

  # ===========================
  # TAB 7: Metabolic
  # ===========================
  output$plot_glucose_ins <- renderPlotly({
    df <- sim_result()
    plot_ly(df, x = ~time/24) %>%
      add_trace(y = ~glucose_mmol, name = "Blood glucose (mmol/L)",
                type = "scatter", mode = "lines", line = list(color = "#E53935")) %>%
      add_trace(y = ~insulin_uU / 10, name = "Insulin/10 (μU/mL)",
                type = "scatter", mode = "lines", line = list(color = "#1E88E5")) %>%
      layout(title = "Blood glucose and insulin dynamics",
             xaxis = list(title = "Time (days)"),
             yaxis = list(title = "Value"))
  })

  output$plot_body_comp <- renderPlotly({
    df <- sim_result()
    plot_ly(df, x = ~time/24) %>%
      add_trace(y = ~VAT_kg, name = "Visceral fat (kg)",
                type = "scatter", mode = "lines", line = list(color = "#F4511E")) %>%
      add_trace(y = ~muscle_kg, name = "Skeletal muscle (kg)",
                type = "scatter", mode = "lines", line = list(color = "#43A047")) %>%
      layout(title = "Body composition change",
             xaxis = list(title = "Time (days)"),
             yaxis = list(title = "Mass (kg)"))
  })

  output$plot_BMD <- renderPlotly({
    df <- sim_result()
    plot_ly(df, x = ~time/24, y = ~BMD_T, type = "scatter", mode = "lines",
            line = list(color = "#795548", width = 2)) %>%
      add_segments(x = 0, xend = max(df$time)/24, y = -1.0, yend = -1.0,
                   line = list(color = "orange", dash = "dash")) %>%
      add_segments(x = 0, xend = max(df$time)/24, y = -2.5, yend = -2.5,
                   line = list(color = "red", dash = "dot")) %>%
      layout(title = "Bone mineral density change (BMD T-score)",
             xaxis = list(title = "Time (days)"),
             yaxis = list(title = "T-score"))
  })

  output$plot_BP <- renderPlotly({
    df <- sim_result()
    plot_ly(df, x = ~time/24, y = ~BP_sys, type = "scatter", mode = "lines",
            line = list(color = "#1565C0", width = 2)) %>%
      add_segments(x = 0, xend = max(df$time)/24, y = 140, yend = 140,
                   line = list(color = "red", dash = "dash")) %>%
      layout(title = "Systolic blood pressure dynamics",
             xaxis = list(title = "Time (days)"),
             yaxis = list(title = "Blood pressure (mmHg)"))
  })

  output$metabolic_risk <- renderUI({
    df <- tail(sim_result(), 1)
    risks <- c(
      if(df$glucose_mmol > 7.0) "Hyperglycaemia / diabetes risk" else NULL,
      if(df$BMD_T < -1.0) "Osteopenia / osteoporosis risk" else NULL,
      if(df$BP_sys > 140) "Hypertension risk" else NULL,
      if(df$muscle_kg < 22) "Sarcopenia risk" else NULL
    )
    if (length(risks) == 0) risks <- "No major metabolic risk within the current simulation period"
    tagList(
      h4("Current metabolic risk assessment"),
      tags$ul(lapply(risks, tags$li))
    )
  })

  # ===========================
  # TAB 8: Virtual Population
  # ===========================
  virt_data <- eventReactive(input$run_virtual, {
    n <- input$n_virtual
    set.seed(42)
    results <- lapply(seq_len(n), function(i) {
      tf_i   <- rnorm(1, input$tumor_fold, input$tumor_fold * input$var_tumor / 100)
      tf_i   <- max(1.5, tf_i)
      kf_i   <- rnorm(1, 0.30, 0.30 * input$var_enz / 100)
      kf_i   <- max(0.10, kf_i)

      drug_i <- if(input$virt_drug == "none") NULL else input$virt_drug
      dose_i <- switch(input$virt_drug,
                       "none" = 0, "ketoconazole" = 400, "osilodrostat" = 5)
      freq_i <- 12

      tryCatch({
        run_simulation(1.0, tf_i, 30, drug_i, dose_i, freq_i, 25, 55, 7.8) %>%
          filter(time == max(time)) %>%
          mutate(subject_id = i, tumor_fold_i = tf_i)
      }, error = function(e) NULL)
    })
    bind_rows(Filter(Negate(is.null), results))
  })

  output$plot_virtual_ufc <- renderPlotly({
    df <- virt_data()
    req(nrow(df) > 0)
    plot_ly(df, x = ~subject_id, y = ~UFC_cumul, type = "bar",
            marker = list(color = ~ifelse(UFC_cumul < 50, "#43A047", "#E53935"))) %>%
      add_segments(x = 0, xend = nrow(df), y = 50, yend = 50,
                   line = list(color = "orange", dash = "dash")) %>%
      layout(title = paste0("UFC distribution (N=", nrow(df), " virtual patients)"),
             xaxis = list(title = "Patient ID"),
             yaxis = list(title = "UFC (μg/24h)"))
  })

  output$plot_virtual_response <- renderPlotly({
    df <- virt_data()
    req(nrow(df) > 0)
    resp_rate <- mean(df$UFC_cumul < 50) * 100
    plot_ly() %>%
      add_pie(labels = c("Responder (UFC<50)", "Non-responder"),
              values = c(resp_rate, 100 - resp_rate),
              marker = list(colors = c("#43A047", "#E53935"))) %>%
      layout(title = paste0("Treatment response rate: ", round(resp_rate, 1), "%\n(UFC normalisation criterion)"))
  })

  output$plot_virtual_dist <- renderPlotly({
    df <- virt_data()
    req(nrow(df) > 0)
    plot_ly(df, x = ~cortisol_free, type = "histogram",
            marker = list(color = "#1565C0", line = list(color = "white", width = 0.5))) %>%
      add_segments(x = 23, xend = 23, y = 0, yend = nrow(df)/5,
                   line = list(color = "red", dash = "dash")) %>%
      layout(title = "Day 30 cortisol distribution (virtual population)",
             xaxis = list(title = "Free cortisol (μg/dL)"),
             yaxis = list(title = "Number of patients"))
  })
}

shinyApp(ui = ui, server = server)
