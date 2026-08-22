## ============================================================
## Wilson's Disease QSP — Interactive Shiny Application
## File: wd_shiny_app.R
## Tabs:
##   1. Patient Profile
##   2. Drug Pharmacokinetics
##   3. Copper Kinetics
##   4. Hepatic Outcomes
##   5. Neuro/Ophthalmic Outcomes
##   6. Scenario Comparison
##   7. Biomarker Explorer
##   8. Model Information
## ============================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(DT)
library(shinythemes)
library(plotly)

## ============================================================
## Embed mrgsolve model code
## ============================================================
wd_code <- '
$PARAM @annotated
Cu_intake   : 1.2   : Dietary Cu (mg/day)
f_abs_Cu    : 0.60  : GI Cu absorption fraction
ATP7B_func  : 0.05  : Residual ATP7B function (0-1)
k_bil_WT    : 3.0   : Biliary Cu excretion rate (WT)
k_Cp_synth  : 0.15  : Cp synthesis rate
k_Cp_deg    : 0.03  : Cp degradation rate
Cp_baseline : 30.0  : Baseline Cp (WT)
k_hep_NCBC  : 0.08  : Hepatic Cu overflow rate
k_NCBC_uri  : 0.15  : NCBC renal CL
k_NCBC_brain: 0.005 : NCBC brain transfer
k_NCBC_kid  : 0.003 : NCBC kidney transfer
k_NCBC_corn : 0.0005: NCBC cornea transfer
k_brain_out : 0.001 : Brain Cu clearance
k_kidney_out: 0.002 : Kidney Cu clearance
k_cornea_out: 0.0002: Cornea Cu clearance
k_ROS_gen   : 0.4   : ROS generation
k_ROS_scav  : 0.2   : ROS scavenging
k_ALT_prod  : 5.0   : ALT release rate
k_ALT_elim  : 0.1   : ALT elimination
ALT_base    : 20.0  : Baseline ALT
k_fib_prog  : 0.002 : Fibrosis progression
k_fib_reg   : 0.001 : Fibrosis regression
Fib_max     : 4.0   : Maximum Metavir score
k_neuro_prog: 0.003 : Neurodegeneration progression
k_neuro_reg : 0.001 : Neurological recovery
MT_max      : 50.0  : MT max capacity
MT_Km       : 25.0  : MT half-saturation
k_MT_degrad : 0.05  : MT degradation rate
ka_DPA      : 2.0   : DPA absorption ka
F_DPA       : 0.55  : DPA bioavailability
Vd_DPA      : 0.5   : DPA Vd (L/kg)
CL_DPA      : 3.5   : DPA clearance (L/hr)
Kchel_DPA   : 0.0   : DPA chelation potency
ka_Zn       : 1.5   : Zinc ka
F_Zn        : 0.15  : Zinc bioavailability
Vd_Zn       : 0.25  : Zinc Vd
CL_Zn       : 1.0   : Zinc clearance
IC50_Zn     : 2.0   : Zinc IC50 for Cu absorption
ka_TRI      : 1.8   : Trientine ka
F_TRI       : 0.45  : Trientine bioavailability
Vd_TRI      : 0.4   : Trientine Vd
CL_TRI      : 4.0   : Trientine clearance
Kchel_TRI   : 0.0   : Trientine chelation
ka_TTM      : 0.8   : ALXN1840 ka
F_TTM       : 0.30  : ALXN1840 bioavailability
Vd_TTM      : 1.0   : ALXN1840 Vd
CL_TTM      : 0.7   : ALXN1840 clearance
Emax_TTM    : 0.0   : TTM Emax (0=off)
EC50_TTM    : 0.5   : TTM EC50

$INIT
GUT_DPA:0 CENT_DPA:0 GUT_ZN:0 CENT_ZN:0
GUT_TRI:0 CENT_TRI:0 GUT_TTM:0 CENT_TTM:0
CU_GI:0.5 CU_HEP:80 MT_HEP:20 CU_NCBC:25 CP_SERUM:15
CU_URINE:0 CU_BRAIN:5 CU_KIDNEY:2 CU_CORNEA:1
ROS_HEP:5 ALT_SERUM:45 FIBROSIS:0.5 NEURODEGENERATION:2

$ODE
double Vd_DPA_tot = Vd_DPA * 70;
double Cc_DPA = CENT_DPA / Vd_DPA_tot;
double Vd_ZN_tot = Vd_Zn * 70;
double Cc_ZN = CENT_ZN / Vd_ZN_tot;
double Zn_eff = Cc_ZN / (IC50_Zn + Cc_ZN);
double Vd_TRI_tot = Vd_TRI * 70;
double Cc_TRI = CENT_TRI / Vd_TRI_tot;
double Vd_TTM_tot = Vd_TTM * 70;
double Cc_TTM = CENT_TTM / Vd_TTM_tot;
double TTM_eff = Emax_TTM * Cc_TTM / (EC50_TTM + Cc_TTM);

dxdt_GUT_DPA  = -ka_DPA * GUT_DPA;
dxdt_CENT_DPA = ka_DPA * GUT_DPA - (CL_DPA / Vd_DPA_tot) * CENT_DPA;
dxdt_GUT_ZN   = -ka_Zn * GUT_ZN;
dxdt_CENT_ZN  = ka_Zn * GUT_ZN - (CL_Zn / Vd_ZN_tot) * CENT_ZN;
dxdt_GUT_TRI  = -ka_TRI * GUT_TRI;
dxdt_CENT_TRI = ka_TRI * GUT_TRI - (CL_TRI / Vd_TRI_tot) * CENT_TRI;
dxdt_GUT_TTM  = -ka_TTM * GUT_TTM;
dxdt_CENT_TTM = ka_TTM * GUT_TTM - (CL_TTM / Vd_TTM_tot) * CENT_TTM;

double f_abs_eff = f_abs_Cu * (1.0 - Zn_eff * 0.8);
double Cu_absorb_rate = Cu_intake / 24.0;
dxdt_CU_GI = Cu_absorb_rate - f_abs_eff * CU_GI;

double k_bil_eff = k_bil_WT * ATP7B_func;
double MT_saturation = CU_HEP / (MT_Km + CU_HEP);
double MT_binding = k_MT_degrad * MT_HEP;
double Cu_to_MT = (1.0 - MT_saturation) * CU_HEP * 0.1;
double DPA_chel = Kchel_DPA * Cc_DPA;
double TRI_chel = Kchel_TRI * Cc_TRI;
double total_chelation = (DPA_chel + TRI_chel) * CU_HEP / (CU_HEP + 10.0);

dxdt_CU_HEP = f_abs_eff * CU_GI * 24.0
             - k_bil_eff * CU_HEP
             - Cu_to_MT + MT_binding * 0.5
             - total_chelation;
dxdt_MT_HEP = Cu_to_MT - MT_binding - MT_HEP * 0.01;

double Cp_synth = k_Cp_synth * CU_HEP * ATP7B_func;
dxdt_CP_SERUM = Cp_synth - k_Cp_deg * CP_SERUM;

double k_NCBC_gen = k_hep_NCBC * MT_saturation * CU_HEP;
double NCBC_chel_out = (DPA_chel*0.6 + TRI_chel*0.4) * CU_NCBC/(CU_NCBC+5.0);
double NCBC_TTM_out  = TTM_eff * CU_NCBC;
dxdt_CU_NCBC = k_NCBC_gen - k_NCBC_uri*CU_NCBC
              - k_NCBC_brain*CU_NCBC - k_NCBC_kid*CU_NCBC
              - k_NCBC_corn*CU_NCBC - NCBC_chel_out - NCBC_TTM_out;

dxdt_CU_URINE = k_NCBC_uri*CU_NCBC + DPA_chel*CU_HEP*0.01 + TRI_chel*CU_HEP*0.005;
dxdt_CU_BRAIN  = k_NCBC_brain*CU_NCBC - k_brain_out*CU_BRAIN - TTM_eff*0.01*CU_BRAIN;
dxdt_CU_KIDNEY = k_NCBC_kid*CU_NCBC - k_kidney_out*CU_KIDNEY;
dxdt_CU_CORNEA = k_NCBC_corn*CU_NCBC - k_cornea_out*CU_CORNEA;

double ROS_gen  = k_ROS_gen  * (CU_HEP/100.0) * (CU_HEP/100.0);
double ROS_scav = k_ROS_scav * ROS_HEP;
dxdt_ROS_HEP = ROS_gen - ROS_scav;

double ALT_gen = k_ALT_prod * (ROS_HEP/20.0);
dxdt_ALT_SERUM = ALT_gen - k_ALT_elim*(ALT_SERUM - ALT_base);

double fib_prog = k_fib_prog * ROS_HEP * FIBROSIS * (1.0 - FIBROSIS/Fib_max);
double fib_reg  = k_fib_reg  * (DPA_chel + TRI_chel + TTM_eff) * FIBROSIS;
dxdt_FIBROSIS = fib_prog - fib_reg;
if(FIBROSIS > Fib_max) dxdt_FIBROSIS = 0.0;
if(FIBROSIS < 0.0)     dxdt_FIBROSIS = 0.0;

double neuro_prog = k_neuro_prog * CU_BRAIN;
double neuro_reg  = k_neuro_reg  * (TTM_eff*0.5 + DPA_chel*0.1);
dxdt_NEURODEGENERATION = neuro_prog - neuro_reg;
if(NEURODEGENERATION < 0.0) dxdt_NEURODEGENERATION = 0.0;

$TABLE
capture Cc_DPA  = CENT_DPA / (Vd_DPA * 70);
capture Cc_ZN   = CENT_ZN  / (Vd_Zn  * 70);
capture Cc_TRI  = CENT_TRI / (Vd_TRI * 70);
capture Cc_TTM  = CENT_TTM / (Vd_TTM * 70);
capture NCBC    = CU_NCBC;
capture Cp_mg   = CP_SERUM;
capture Cu_hep  = CU_HEP;
capture Cu_brain_idx = CU_BRAIN;
capture Cu_corn_idx  = CU_CORNEA;
capture Cu_uri_rate  = CU_URINE;
capture ROS_idx = ROS_HEP;
capture ALT_val = ALT_SERUM;
capture Fib_val = FIBROSIS;
capture Neuro_val = NEURODEGENERATION;
capture Cu_total_serum = CP_SERUM * 3.15 + CU_NCBC;
capture UWDRS_proxy = NEURODEGENERATION * 5.0;
'

## Compile model once at startup
wd_mod <- mcode("WD_Shiny", wd_code)

## ============================================================
## Simulation helper
## ============================================================
run_sim <- function(mod,
                    atp7b_func = 0.05,
                    years = 3,
                    drug = "none",
                    dose_mg = 500,
                    interval_h = 8,
                    cu_intake = 1.2,
                    init_CU_HEP = 80,
                    init_NCBC = 25,
                    init_Neuro = 2) {

  hrs <- years * 365 * 24
  n_doses <- years * 365 * (24 / interval_h)

  init_vals <- list(
    CU_HEP = init_CU_HEP,
    CU_NCBC = init_NCBC,
    NEURODEGENERATION = init_Neuro
  )

  p <- list(
    ATP7B_func = atp7b_func,
    Cu_intake  = cu_intake
  )

  e <- ev(time = 0, amt = 0, cmt = "GUT_DPA")  # dummy

  if (drug == "DPA") {
    p$Kchel_DPA <- 0.3
    e <- ev(amt = dose_mg * 0.55, cmt = "GUT_DPA",
            ii = interval_h, addl = round(n_doses) - 1, time = 0)
  } else if (drug == "Zinc") {
    p$IC50_Zn <- 2.0
    e <- ev(amt = dose_mg * 0.15, cmt = "GUT_ZN",
            ii = interval_h, addl = round(n_doses) - 1, time = 0)
  } else if (drug == "Trientine") {
    p$Kchel_TRI <- 0.15
    e <- ev(amt = dose_mg * 0.45, cmt = "GUT_TRI",
            ii = interval_h, addl = round(n_doses) - 1, time = 0)
  } else if (drug == "ALXN1840") {
    p$Emax_TTM <- 0.98
    e <- ev(amt = 15 * 0.30, cmt = "GUT_TTM",
            ii = 24, addl = years * 365 - 1, time = 0)
  }

  mod %>%
    param(p) %>%
    init(init_vals) %>%
    mrgsim(ev = e, end = hrs, delta = 24) %>%
    as.data.frame() %>%
    mutate(Day = time / 24, Year = time / (365 * 24))
}

## ============================================================
## UI
## ============================================================
ui <- fluidPage(
  theme = shinytheme("darkly"),
  tags$head(tags$style(HTML("
    .navbar-brand { font-weight: bold; font-size: 18px; }
    .well { background-color: #2c3e50; border: 1px solid #3d5166; }
    h4 { color: #3498db; }
    .tab-content { padding-top: 15px; }
    .summary-box { background: #34495e; padding: 10px; border-radius: 5px;
                   margin: 5px; text-align: center; }
    .summary-value { font-size: 24px; font-weight: bold; color: #3498db; }
    .summary-label { font-size: 11px; color: #bdc3c7; }
  "))),

  navbarPage(
    title = "Wilson's Disease QSP Model",
    id = "main_nav",

    ## ---- TAB 1: Patient profile ----
    tabPanel("1. Patient profile",
      sidebarLayout(
        sidebarPanel(width = 3,
          h4("Patient settings"),
          selectInput("ptype", "WD phenotype",
            choices = c("Hepatic" = "hepatic",
                        "Neuropsychiatric" = "neuro",
                        "Mixed" = "mixed",
                        "Asymptomatic" = "asymp")),
          selectInput("mutation", "ATP7B mutation",
            choices = c("p.His1069Gln (Europe, most common)" = "H1069Q",
                        "p.Arg778Leu (Asian type)" = "R778L",
                        "p.Gly710Ser (US type)" = "G710S",
                        "Other truncating mutations" = "trunc",
                        "Normal WT" = "WT")),
          sliderInput("atp7b_f", "Residual ATP7B function (%)",
                      min = 0, max = 100, value = 5, step = 1),
          sliderInput("cu_diet", "Copper intake (mg/day)",
                      min = 0.5, max = 3.0, value = 1.2, step = 0.1),
          hr(),
          h4("Initial disease state"),
          sliderInput("init_hep_cu", "Initial hepatic copper (μg/g dw)",
                      min = 10, max = 500, value = 80),
          sliderInput("init_ncbc", "Initial NCBC (μg/dL)",
                      min = 5, max = 80, value = 25),
          sliderInput("init_neuro", "Initial neurological damage index",
                      min = 0, max = 20, value = 2)
        ),
        mainPanel(width = 9,
          fluidRow(
            column(3, div(class="summary-box",
              div("Leipzig Score", class="summary-label"),
              div(textOutput("Leipzig_out"), class="summary-value")
            )),
            column(3, div(class="summary-box",
              div("Expected NCBC", class="summary-label"),
              div(textOutput("NCBC_init_out"), class="summary-value")
            )),
            column(3, div(class="summary-box",
              div("Expected Cp (mg/dL)", class="summary-label"),
              div(textOutput("Cp_init_out"), class="summary-value")
            )),
            column(3, div(class="summary-box",
              div("Hepatic copper status", class="summary-label"),
              div(textOutput("HepCu_stage_out"), class="summary-value")
            ))
          ),
          hr(),
          fluidRow(
            column(6,
              h4("WD diagnostic criteria (Leipzig Score)"),
              DTOutput("leipzig_table")
            ),
            column(6,
              h4("Disease stage definitions"),
              DTOutput("staging_table")
            )
          ),
          hr(),
          h4("Genotype-phenotype correlation"),
          plotOutput("genotype_plot", height = "280px")
        )
      )
    ),

    ## ---- TAB 2: Drug PK ----
    tabPanel("2. Drug PK",
      sidebarLayout(
        sidebarPanel(width = 3,
          h4("Drug selection"),
          checkboxGroupInput("pk_drugs", "Drugs to display",
            choices = c("D-Penicillamine" = "DPA",
                        "Zinc Acetate" = "Zinc",
                        "Trientine" = "Trientine",
                        "ALXN1840 (TTM)" = "ALXN1840"),
            selected = c("DPA", "ALXN1840")),
          hr(),
          h4("DPA settings"),
          numericInput("dpa_dose", "DPA dose (mg/dose)", value = 500, min = 125, max = 1500),
          selectInput("dpa_freq", "Dosing interval",
                      choices = c("QD (once daily)" = "24", "BID (twice daily)" = "12",
                                  "TID (three times daily)" = "8", "QID (four times daily)" = "6"),
                      selected = "8"),
          hr(),
          h4("Simulation"),
          sliderInput("pk_days", "Duration (days)", min = 1, max = 90, value = 30),
          actionButton("run_pk", "Run PK simulation", class = "btn-primary btn-block")
        ),
        mainPanel(width = 9,
          tabsetPanel(
            tabPanel("PK profile",
              plotlyOutput("pk_plot", height = "400px"),
              br(),
              h4("PK parameter summary"),
              DTOutput("pk_params_table")
            ),
            tabPanel("NCBC vs time",
              plotlyOutput("ncbc_time_plot", height = "400px"),
              p("NCBC = Non-Ceruloplasmin Bound Copper. Target: <15 μg/dL")
            ),
            tabPanel("PK parameters",
              fluidRow(
                column(12,
                  h4("PK parameters by drug"),
                  DTOutput("full_pk_params")
                )
              )
            )
          )
        )
      )
    ),

    ## ---- TAB 3: Copper kinetics ----
    tabPanel("3. Copper kinetics",
      sidebarLayout(
        sidebarPanel(width = 3,
          h4("Treatment scenario"),
          selectInput("cu_drug", "Treatment drug",
            choices = c("Untreated" = "none",
                        "D-Penicillamine" = "DPA",
                        "Zinc Acetate" = "Zinc",
                        "Trientine" = "Trientine",
                        "ALXN1840 (TTM)" = "ALXN1840")),
          sliderInput("cu_years", "Simulation duration (years)", min = 1, max = 10, value = 5),
          hr(),
          h4("Treatment targets"),
          div(class="summary-box",
            p("Treatment target values:", class="summary-label"),
            tags$ul(
              tags$li("NCBC < 15 μg/dL"),
              tags$li("Ceruloplasmin > 20 mg/dL"),
              tags$li("24h urinary Cu < 100 μg/day (when stable)"),
              tags$li("Hepatic Cu < 50 μg/g dw")
            )
          ),
          actionButton("run_cu", "Run copper simulation", class = "btn-success btn-block")
        ),
        mainPanel(width = 9,
          tabsetPanel(
            tabPanel("Serum copper",
              plotlyOutput("serum_cu_plot", height = "380px"),
              fluidRow(
                column(4, div(class="summary-box",
                  div("NCBC (1 year of treatment)", class="summary-label"),
                  div(textOutput("ncbc_1yr"), class="summary-value")
                )),
                column(4, div(class="summary-box",
                  div("Ceruloplasmin (1 year)", class="summary-label"),
                  div(textOutput("cp_1yr"), class="summary-value")
                )),
                column(4, div(class="summary-box",
                  div("Hepatic Cu (1 year, μg/g)", class="summary-label"),
                  div(textOutput("hepcu_1yr"), class="summary-value")
                ))
              )
            ),
            tabPanel("Organ copper distribution",
              plotlyOutput("organ_cu_plot", height = "400px")
            ),
            tabPanel("Urinary copper",
              plotlyOutput("urine_cu_plot", height = "380px"),
              p("A steep rise in urinary copper during early treatment indicates copper mobilisation.")
            )
          )
        )
      )
    ),

    ## ---- TAB 4: Hepatic outcomes ----
    tabPanel("4. Hepatic outcomes",
      sidebarLayout(
        sidebarPanel(width = 3,
          h4("Hepatic disease parameters"),
          selectInput("hep_drug", "Treatment drug",
            choices = c("Untreated" = "none",
                        "D-Penicillamine" = "DPA",
                        "Zinc Acetate" = "Zinc",
                        "Trientine" = "Trientine",
                        "ALXN1840" = "ALXN1840")),
          sliderInput("hep_atp7b", "Residual ATP7B function (%)", min = 0, max = 30, value = 5),
          sliderInput("hep_years", "Follow-up duration (years)", min = 1, max = 10, value = 5),
          sliderInput("fib_init", "Initial fibrosis score (Metavir)", min = 0, max = 4, value = 1, step = 0.5),
          hr(),
          h4("Liver transplant risk index"),
          sliderInput("meld_score", "MELD Score", min = 6, max = 40, value = 10),
          div(textOutput("lt_risk"), style = "color: #e74c3c; font-weight: bold;"),
          actionButton("run_hep", "Run hepatic outcome simulation", class = "btn-warning btn-block")
        ),
        mainPanel(width = 9,
          tabsetPanel(
            tabPanel("ALT & fibrosis",
              plotlyOutput("hep_altscore_plot", height = "400px")
            ),
            tabPanel("ROS & oxidative stress",
              plotlyOutput("ros_plot", height = "380px"),
              p("Oxidative stress index: a marker of the ROS generated by the Fenton reaction driven by copper accumulation")
            ),
            tabPanel("Hepatic pathophysiology summary",
              fluidRow(
                column(12,
                  h4("WD hepatic staging"),
                  DTOutput("hep_staging_dt")
                )
              )
            )
          )
        )
      )
    ),

    ## ---- TAB 5: Neuro/ophthalmic outcomes ----
    tabPanel("5. Neurological/ophthalmic outcomes",
      sidebarLayout(
        sidebarPanel(width = 3,
          h4("Neurological disease parameters"),
          selectInput("neuro_drug", "Treatment drug",
            choices = c("Untreated" = "none",
                        "D-Penicillamine (caution: paradoxical worsening)" = "DPA",
                        "ALXN1840 (preferred in the neurological form)" = "ALXN1840",
                        "Trientine" = "Trientine")),
          sliderInput("neuro_years", "Follow-up duration (years)", min = 1, max = 10, value = 5),
          sliderInput("init_brain_cu", "Initial brain copper index",
                      min = 1, max = 20, value = 5),
          sliderInput("init_neuro_dmg", "Initial neurological damage index",
                      min = 0, max = 20, value = 3),
          hr(),
          div(class="summary-box",
            p("Risk of paradoxical neurological worsening on DPA:", class="summary-label"),
            p("in neurological WD, on initial DPA dosing", class="summary-label"),
            p("worsening of neurological symptoms is reported in ~50%", class="summary-label"),
            p("(a transient rise in brain Cu from the Cu mobilisation effect)", class="summary-label")
          ),
          actionButton("run_neuro", "Run neurological simulation", class = "btn-info btn-block")
        ),
        mainPanel(width = 9,
          tabsetPanel(
            tabPanel("Brain copper & neurodegeneration",
              plotlyOutput("brain_cu_plot", height = "400px")
            ),
            tabPanel("KF Ring & UWDRS",
              plotlyOutput("kf_uwdrs_plot", height = "380px"),
              p("KF Ring: can disappear more than 5 years after treatment starts (falling corneal Cu)")
            ),
            tabPanel("Neurological symptom progression",
              fluidRow(
                column(12,
                  h4("WD neurological symptom classification"),
                  DTOutput("neuro_symptoms_dt")
                )
              )
            )
          )
        )
      )
    ),

    ## ---- TAB 6: Scenario comparison ----
    tabPanel("6. Scenario comparison",
      sidebarLayout(
        sidebarPanel(width = 3,
          h4("Comparison scenarios"),
          checkboxGroupInput("compare_scenarios",
            "Select scenarios to compare",
            choices = c(
              "S1: untreated WD" = "S1",
              "S2: DPA 500mg TID" = "S2",
              "S3: Zinc 50mg TID" = "S3",
              "S4: Trientine 500mg TID" = "S4",
              "S5: ALXN1840 15mg QD" = "S5",
              "S6: DPA→Zinc switch (1 year)" = "S6",
              "S7: ALXN1840+Trientine combination" = "S7",
              "S8: normal WT control" = "S8"
            ),
            selected = c("S1", "S2", "S5", "S8")
          ),
          sliderInput("comp_years", "Simulation duration (years)", min = 1, max = 10, value = 5),
          hr(),
          h4("Comparison metric"),
          selectInput("comp_endpoint",
            "Primary comparison metric",
            choices = c(
              "NCBC (free Cu)" = "NCBC",
              "Serum ALT" = "ALT",
              "Hepatic copper" = "Cu_hep",
              "Fibrosis score" = "Fib",
              "Neurodegeneration index" = "Neuro",
              "Brain copper" = "Cu_brain",
              "KF Ring index" = "KF"
            )
          ),
          actionButton("run_compare", "Run comparison", class = "btn-primary btn-block")
        ),
        mainPanel(width = 9,
          tabsetPanel(
            tabPanel("Comparison plot",
              plotlyOutput("compare_plot", height = "450px")
            ),
            tabPanel("Results summary table",
              h4("Summary of results at the end of the simulation"),
              DTOutput("compare_table")
            ),
            tabPanel("ATLAS trial comparison",
              h4("ATLAS Trial (ALXN1840) vs simulation"),
              fluidRow(
                column(6,
                  div(class="summary-box",
                    h4("ATLAS Trial results (NEJM 2022)"),
                    tags$ul(
                      tags$li("NCBC ↓ 98.6% (vs placebo)"),
                      tags$li("Neurological WD patients enrolled preferentially"),
                      tags$li("Significant improvement in UWDRS"),
                      tags$li("No serious adverse events"),
                      tags$li("Indication: neurological WD, first-line treatment")
                    )
                  )
                ),
                column(6,
                  div(class="summary-box",
                    h4("DPA vs ALXN1840 — the core comparison"),
                    tags$ul(
                      tags$li("DPA: potent in hepatic WD, paradoxical worsening in the neurological form"),
                      tags$li("ALXN1840: superior in neurological WD"),
                      tags$li("Zinc: safe, maintenance therapy / pregnancy"),
                      tags$li("Trientine: 2nd-line when DPA causes adverse effects")
                    )
                  )
                )
              )
            )
          )
        )
      )
    ),

    ## ---- TAB 7: Biomarker explorer ----
    tabPanel("7. Biomarkers",
      sidebarLayout(
        sidebarPanel(width = 3,
          h4("Copper biomarkers"),
          selectInput("bm_x", "X-axis biomarker",
            choices = c("Day", "NCBC", "Cu_hep", "Cu_brain_idx",
                        "Cu_corn_idx", "Cu_total_serum", "Cp_mg",
                        "Cu_uri_rate", "ALT_val", "Fib_val")),
          selectInput("bm_y", "Y-axis biomarker",
            choices = c("NCBC", "ALT_val", "Fib_val", "Neuro_val",
                        "Cu_brain_idx", "Cu_corn_idx", "UWDRS_proxy",
                        "Cp_mg", "Cu_total_serum"),
            selected = "ALT_val"),
          checkboxGroupInput("bm_drugs", "Drugs to display",
            choices = c("none", "DPA", "Zinc", "Trientine", "ALXN1840"),
            selected = c("none", "DPA", "ALXN1840")),
          sliderInput("bm_years", "Duration (years)", min = 1, max = 10, value = 3),
          hr(),
          h4("Diagnostic thresholds"),
          div(class="summary-box",
            tags$ul(
              tags$li("NCBC >20 μg/dL: WD diagnostic criterion"),
              tags$li("Cp <20 mg/dL: WD suspected"),
              tags$li("24h urinary Cu >100 μg: significant copper excretion"),
              tags$li("Hepatic Cu >250 μg/g dw: WD diagnosis"),
              tags$li("NCBC <15 μg/dL: treatment target")
            )
          ),
          actionButton("run_bm", "Analyse biomarkers", class = "btn-success btn-block")
        ),
        mainPanel(width = 9,
          tabsetPanel(
            tabPanel("Biomarker correlation",
              plotlyOutput("bm_scatter", height = "400px")
            ),
            tabPanel("Biomarker time series",
              plotlyOutput("bm_timeseries", height = "400px")
            ),
            tabPanel("Diagnostic panel",
              h4("WD diagnostic biomarker panel"),
              DTOutput("diagnostic_panel")
            )
          )
        )
      )
    ),

    ## ---- TAB 8: Model information ----
    tabPanel("8. Model information",
      fluidRow(
        column(6,
          h4("Model structure"),
          div(class="summary-box",
            h5("24 ODE compartments"),
            tags$ul(
              tags$li("Drug PK: DPA, Zinc, Trientine, ALXN1840 (8 compartments)"),
              tags$li("Copper kinetics: GI, liver, MT, NCBC, CP, urine (6 compartments)"),
              tags$li("Organ distribution: brain, kidney, cornea (3 compartments)"),
              tags$li("Hepatic pathology: ROS, ALT, fibrosis (3 compartments)"),
              tags$li("Neuropathology: neurodegeneration (1 compartment)"),
              tags$li("24-compartment ODE system in total")
            )
          ),
          br(),
          h4("Sources of the main parameters"),
          DTOutput("param_ref_table")
        ),
        column(6,
          h4("Clinical trial evidence"),
          div(class="summary-box",
            h5("ATLAS Trial (ALXN1840, NEJM 2022)"),
            p("N=173, randomised double-blind placebo-controlled"),
            p("Primary result: NCBC ↓98.6% (p<0.001)"),
            p("Significant improvement in UWDRS score"),
            hr(),
            h5("RCT Wilson (Trientine, 2011)"),
            p("Trientine vs DPA: similar efficacy"),
            p("Fewer DPA adverse effects"),
            hr(),
            h5("Zinc Acetate (Brewer 1998)"),
            p("15-year long-term efficacy confirmed as maintenance therapy"),
            p("Excellent safety profile")
          ),
          br(),
          h4("Mechanistic map"),
          tags$img(src = "wd_qsp_model.png",
                   style = "width:100%; border-radius:5px;",
                   alt = "WD QSP Mechanistic Map")
        )
      )
    )
  )
)

## ============================================================
## Server
## ============================================================
server <- function(input, output, session) {

  ## ---- Tab 1: Patient Profile ----
  output$Leipzig_out <- renderText({
    score <- 2  # KF rings
    if (input$init_ncbc > 20) score <- score + 2
    if (input$init_ncbc > 5 && input$init_ncbc <= 20) score <- score + 1
    if (input$atp7b_f < 20) score <- score + 2
    if (input$init_hep_cu > 250) score <- score + 2
    paste0(score, " pts")
  })
  output$NCBC_init_out  <- renderText({ paste0(input$init_ncbc, " μg/dL") })
  output$Cp_init_out    <- renderText({
    cp <- round(15 * (input$atp7b_f / 5) * (0.9 + 0.1), 1)
    paste0(min(cp, 40), " mg/dL")
  })
  output$HepCu_stage_out <- renderText({
    cu <- input$init_hep_cu
    if (cu > 250) "Severe"
    else if (cu > 100) "Moderate"
    else if (cu > 50)  "Mild"
    else "Normal range"
  })

  output$leipzig_table <- renderDT({
    df <- data.frame(
      Item = c("KF Ring present", "NCBC >25 μg/dL", "NCBC 5-25 μg/dL",
               "ATP7B pathogenic variant (biallelic)", "Hepatic Cu >250 μg/g dw",
               "Serum Cp <20 mg/dL", "Neurological symptoms / brain MRI change",
               "Coombs-negative haemolytic anaemia"),
      Score = c(2, 2, 1, 4, 2, 1, 1, 1),
      Interpretation = c("Slit-lamp examination of both eyes", "Direct assay", "Direct assay",
               "Two or more pathogenic variants", "Fresh biopsy tissue",
               "Immunoturbidimetric assay", "Brain MRI abnormality", "Autoimmune-negative")
    )
    datatable(df, options = list(dom = 't', pageLength = 10), rownames = FALSE,
              class = "table-sm")
  })

  output$staging_table <- renderDT({
    df <- data.frame(
      Stage = c("I (preclinical)", "II (hepatic)", "III (neurological)", "IV (mixed)"),
      Features = c("Asymptomatic, hepatic Cu accumulation",
               "Acute/chronic hepatitis, ALF possible",
               "Predominantly neurological/psychiatric symptoms",
               "Simultaneous hepatic + neurological involvement"),
      Leipzig_score = c("≥4", "≥4", "≥4", "≥4"),
      Recommended_treatment = c("Zinc (prophylaxis)", "DPA/Trientine", "ALXN1840 or Trientine", "ALXN1840")
    )
    datatable(df, options = list(dom = 't'), rownames = FALSE,
              class = "table-sm")
  })

  output$genotype_plot <- renderPlot({
    df <- data.frame(
      Mutation = c("p.His1069Gln", "p.Arg778Leu", "p.Met645Arg",
                   "p.Gly710Ser", "p.Ala874Val", "Others"),
      Frequency_pct = c(35, 20, 8, 6, 5, 26),
      Population = c("Europe", "Asia", "Europe", "USA", "Middle East", "Global")
    )
    ggplot(df, aes(x = reorder(Mutation, Frequency_pct), y = Frequency_pct,
                  fill = Population)) +
      geom_col() +
      coord_flip() +
      scale_fill_brewer(palette = "Set2") +
      labs(x = "", y = "Approximate frequency (%)", title = "Frequency of the main ATP7B mutations") +
      theme_minimal(base_size = 13) +
      theme(plot.background = element_rect(fill = "#2c3e50", color = NA),
            panel.background = element_rect(fill = "#2c3e50"),
            text = element_text(color = "white"),
            axis.text = element_text(color = "white"),
            legend.background = element_rect(fill = "#2c3e50"))
  })

  ## ---- Tab 2: Drug PK ----
  pk_data <- eventReactive(input$run_pk, {
    days <- input$pk_days
    hrs <- days * 24
    results <- list()

    if ("DPA" %in% input$pk_drugs) {
      interval_h <- as.numeric(input$dpa_freq)
      e <- ev(amt = input$dpa_dose * 0.55, cmt = "GUT_DPA",
              ii = interval_h, addl = ceiling(hrs / interval_h) - 1, time = 0)
      res <- wd_mod %>%
        param(list(Kchel_DPA = 0.3)) %>%
        mrgsim(ev = e, end = hrs, delta = 1) %>%
        as.data.frame() %>%
        mutate(Drug = "DPA", Day = time / 24, Conc = Cc_DPA * 1000)
      results[["DPA"]] <- res
    }
    if ("Zinc" %in% input$pk_drugs) {
      e <- ev(amt = 50 * 0.15, cmt = "GUT_ZN", ii = 8,
              addl = ceiling(hrs / 8) - 1, time = 0)
      res <- wd_mod %>%
        mrgsim(ev = e, end = hrs, delta = 1) %>%
        as.data.frame() %>%
        mutate(Drug = "Zinc", Day = time / 24, Conc = Cc_ZN * 1000)
      results[["Zinc"]] <- res
    }
    if ("Trientine" %in% input$pk_drugs) {
      e <- ev(amt = 500 * 0.45, cmt = "GUT_TRI", ii = 8,
              addl = ceiling(hrs / 8) - 1, time = 0)
      res <- wd_mod %>%
        param(list(Kchel_TRI = 0.15)) %>%
        mrgsim(ev = e, end = hrs, delta = 1) %>%
        as.data.frame() %>%
        mutate(Drug = "Trientine", Day = time / 24, Conc = Cc_TRI * 1000)
      results[["Trientine"]] <- res
    }
    if ("ALXN1840" %in% input$pk_drugs) {
      e <- ev(amt = 15 * 0.30, cmt = "GUT_TTM", ii = 24,
              addl = days - 1, time = 0)
      res <- wd_mod %>%
        param(list(Emax_TTM = 0.98)) %>%
        mrgsim(ev = e, end = hrs, delta = 1) %>%
        as.data.frame() %>%
        mutate(Drug = "ALXN1840", Day = time / 24, Conc = Cc_TTM * 1000)
      results[["ALXN1840"]] <- res
    }
    bind_rows(results)
  })

  output$pk_plot <- renderPlotly({
    req(pk_data())
    p <- ggplot(pk_data(), aes(x = Day, y = Conc, color = Drug)) +
      geom_line(size = 1) +
      labs(title = "Drug Plasma Concentration vs Time",
           x = "Day", y = "Plasma Conc (μg/L)") +
      scale_color_brewer(palette = "Set1") +
      theme_minimal(base_size = 12)
    ggplotly(p)
  })

  output$pk_params_table <- renderDT({
    df <- data.frame(
      Drug = c("D-Penicillamine", "Zinc", "Trientine", "ALXN1840"),
      `F (%)` = c(50, 15, 45, 30),
      `Tmax (hr)` = c("1-2", "1-2", "1-2", "2"),
      `t½ (hr)` = c("1.7", "2", "2-3", "19"),
      `Vd (L)` = c(35, 17.5, 28, 70),
      `CL (L/hr)` = c(3.5, 1.0, 4.0, 0.7),
      Mechanism = c("Chelation (NCBC, liver)", "MT induction", "Chelation", "TTM tripartite complex"),
      check.names = FALSE
    )
    datatable(df, options = list(dom = 't'), rownames = FALSE, class = "table-sm")
  })

  output$ncbc_time_plot <- renderPlotly({
    req(pk_data())
    p <- ggplot(pk_data(), aes(x = Day, y = NCBC, color = Drug)) +
      geom_line(size = 1) +
      geom_hline(yintercept = 15, linetype = "dashed", color = "yellow") +
      annotate("text", x = 5, y = 16, label = "treatment target <15", color = "yellow") +
      labs(title = "NCBC (Free Copper) Response",
           x = "Day", y = "NCBC (μg/dL)") +
      theme_minimal()
    ggplotly(p)
  })

  output$full_pk_params <- renderDT({
    df <- data.frame(
      Parameter = c("DPA ka", "DPA F", "DPA Vd", "DPA CL", "DPA Kchel",
                    "Zinc ka", "Zinc F", "Zinc IC50 (Cu absorption)",
                    "Trientine ka", "Trientine Kchel",
                    "ALXN1840 Emax (NCBC)", "ALXN1840 EC50", "ALXN1840 t½"),
      Value = c("2.0/hr", "55%", "35 L (0.5 L/kg)", "3.5 L/hr", "0.3/hr·μM",
                "1.5/hr", "15%", "2 μg/mL",
                "1.8/hr", "0.15/hr·μg/mL",
                "98%", "0.5 μg/mL", "19 hr"),
      Source = c(rep("Weiss 2010", 5),
                 rep("Brewer 1998", 3),
                 rep("Weiss 2011", 2),
                 rep("Schilsky 2022 (ATLAS)", 3))
    )
    datatable(df, options = list(dom = 't', pageLength = 15), rownames = FALSE,
              class = "table-sm")
  })

  ## ---- Tab 3: Copper Kinetics ----
  cu_sim <- eventReactive(input$run_cu, {
    run_sim(mod = wd_mod,
            atp7b_func = input$atp7b_f / 100,
            years = input$cu_years,
            drug = input$cu_drug,
            cu_intake = input$cu_diet,
            init_CU_HEP = input$init_hep_cu,
            init_NCBC = input$init_ncbc,
            init_Neuro = input$init_neuro)
  })

  output$serum_cu_plot <- renderPlotly({
    req(cu_sim())
    df <- cu_sim() %>%
      pivot_longer(c(NCBC, Cp_mg, Cu_total_serum), names_to = "Marker", values_to = "Value") %>%
      mutate(Marker = recode(Marker,
        NCBC = "NCBC (μg/dL)", Cp_mg = "Ceruloplasmin (mg/dL)",
        Cu_total_serum = "Total Serum Cu (μg/dL)"))
    p <- ggplot(df, aes(x = Year, y = Value, color = Marker)) +
      geom_line(size = 1) + facet_wrap(~Marker, scales = "free_y") +
      labs(x = "Year", title = "Serum copper markers over time") + theme_minimal()
    ggplotly(p)
  })

  output$ncbc_1yr <- renderText({
    req(cu_sim())
    v <- cu_sim() %>% filter(abs(Year - 1) < 0.05) %>% slice(1) %>% pull(NCBC)
    paste0(round(v, 1), " μg/dL")
  })
  output$cp_1yr <- renderText({
    req(cu_sim())
    v <- cu_sim() %>% filter(abs(Year - 1) < 0.05) %>% slice(1) %>% pull(Cp_mg)
    paste0(round(v, 1), " mg/dL")
  })
  output$hepcu_1yr <- renderText({
    req(cu_sim())
    v <- cu_sim() %>% filter(abs(Year - 1) < 0.05) %>% slice(1) %>% pull(Cu_hep)
    paste0(round(v, 1))
  })

  output$organ_cu_plot <- renderPlotly({
    req(cu_sim())
    df <- cu_sim() %>%
      pivot_longer(c(Cu_hep, Cu_brain_idx, CU_KIDNEY, CU_CORNEA),
                   names_to = "Organ", values_to = "Cu") %>%
      mutate(Organ = recode(Organ,
        Cu_hep = "Liver", Cu_brain_idx = "Brain",
        CU_KIDNEY = "Kidney", CU_CORNEA = "Cornea"))
    p <- ggplot(df, aes(x = Year, y = Cu, color = Organ)) +
      geom_line(size = 1) + facet_wrap(~Organ, scales = "free_y") +
      labs(x = "Year", y = "Cu Index", title = "Copper accumulation by organ") + theme_minimal()
    ggplotly(p)
  })

  output$urine_cu_plot <- renderPlotly({
    req(cu_sim())
    p <- ggplot(cu_sim(), aes(x = Year, y = Cu_uri_rate)) +
      geom_line(size = 1, color = "#f39c12") +
      geom_hline(yintercept = 100, linetype = "dashed", color = "red") +
      annotate("text", x = 0.5, y = 110, label = "diagnostic threshold 100 μg/day",
               color = "red", size = 3) +
      labs(x = "Year", y = "μg/day", title = "24-hour urinary copper excretion rate") +
      theme_minimal()
    ggplotly(p)
  })

  ## ---- Tab 4: Hepatic Outcomes ----
  hep_sim <- eventReactive(input$run_hep, {
    run_sim(mod = wd_mod,
            atp7b_func = input$hep_atp7b / 100,
            years = input$hep_years,
            drug = input$hep_drug)
  })

  output$hep_altscore_plot <- renderPlotly({
    req(hep_sim())
    df <- hep_sim() %>%
      pivot_longer(c(ALT_val, Fib_val), names_to = "Marker", values_to = "Value") %>%
      mutate(Marker = recode(Marker,
        ALT_val = "Serum ALT (U/L)", Fib_val = "Fibrosis Score (Metavir)"))
    p <- ggplot(df, aes(x = Year, y = Value, color = Marker)) +
      geom_line(size = 1) + facet_wrap(~Marker, scales = "free_y") +
      labs(x = "Year", title = "Course of hepatic injury and fibrosis") + theme_minimal()
    ggplotly(p)
  })

  output$ros_plot <- renderPlotly({
    req(hep_sim())
    p <- ggplot(hep_sim(), aes(x = Year, y = ROS_idx)) +
      geom_line(size = 1, color = "#e74c3c") +
      labs(x = "Year", y = "ROS Index", title = "Course of hepatic oxidative stress") +
      theme_minimal()
    ggplotly(p)
  })

  output$lt_risk <- renderText({
    meld <- input$meld_score
    if (meld >= 25) "Urgent indication for liver transplant (MELD ≥25)"
    else if (meld >= 15) "Liver transplant to be actively considered (MELD 15-24)"
    else "Manageable with drug therapy"
  })

  output$hep_staging_dt <- renderDT({
    df <- data.frame(
      Stage = c("WD hepatitis (mild)", "WD hepatitis (moderate)", "WD cirrhosis (compensated)",
               "WD cirrhosis (decompensated)", "Acute liver failure (ALF-WD)"),
      ALT = c("<100", "100-500", "100-300", ">300 or normalised", "Markedly raised"),
      Fibrosis = c("F0-F1", "F1-F2", "F3", "F4", "F4 + ascites/bleeding"),
      `MELD` = c("<10", "10-14", "15-19", "20-24", "≥25"),
      Treatment = c("DPA or Zinc", "DPA + planned switch to Zinc", "DPA + consider TIPS",
               "Liver transplant listing", "Emergency liver transplant (Cu plasma exchange)")
    )
    datatable(df, options = list(dom = 't'), rownames = FALSE, class = "table-sm")
  })

  ## ---- Tab 5: Neuro Outcomes ----
  neuro_sim <- eventReactive(input$run_neuro, {
    run_sim(mod = wd_mod,
            atp7b_func = input$atp7b_f / 100,
            years = input$neuro_years,
            drug = input$neuro_drug,
            init_CU_HEP = 80,
            init_NCBC = 30,
            init_Neuro = input$init_neuro_dmg)
  })

  output$brain_cu_plot <- renderPlotly({
    req(neuro_sim())
    df <- neuro_sim() %>%
      pivot_longer(c(Cu_brain_idx, Neuro_val),
                   names_to = "Marker", values_to = "Value") %>%
      mutate(Marker = recode(Marker,
        Cu_brain_idx = "Brain copper index",
        Neuro_val    = "Neurodegeneration index"))
    p <- ggplot(df, aes(x = Year, y = Value, color = Marker)) +
      geom_line(size = 1) + facet_wrap(~Marker, scales = "free_y") +
      labs(x = "Year", title = "Course of brain copper and neurodegeneration") + theme_minimal()
    ggplotly(p)
  })

  output$kf_uwdrs_plot <- renderPlotly({
    req(neuro_sim())
    df <- neuro_sim() %>%
      pivot_longer(c(Cu_corn_idx, UWDRS_proxy),
                   names_to = "Marker", values_to = "Value") %>%
      mutate(Marker = recode(Marker,
        Cu_corn_idx  = "KF Ring index (corneal Cu)",
        UWDRS_proxy  = "UWDRS estimate"))
    p <- ggplot(df, aes(x = Year, y = Value, color = Marker)) +
      geom_line(size = 1) + facet_wrap(~Marker, scales = "free_y") +
      labs(x = "Year", title = "Course of KF Ring and UWDRS") + theme_minimal()
    ggplotly(p)
  })

  output$neuro_symptoms_dt <- renderDT({
    df <- data.frame(
      Symptom = c("Wing-beating tremor", "Dysarthria",
               "Dysphagia", "Dystonia",
               "Parkinsonism", "Cerebellar ataxia",
               "Behavioural change", "Depression / psychosis"),
      Main_brain_region = c("Cerebello-thalamic circuit", "Brainstem/cerebral cortex",
                     "Basal ganglia/medulla", "Basal ganglia (putamen)",
                     "Nigrostriatal", "Cerebellar dentate nucleus",
                     "Frontal lobe/limbic system", "Frontal lobe/dopaminergic system"),
      UWDRS_contribution = c("High", "High", "Moderate", "High", "Moderate", "Moderate", "Moderate", "Low"),
      ALXN1840_response = c("Excellent", "Excellent", "Moderate", "Excellent", "Excellent", "Moderate", "Excellent", "Fair")
    )
    datatable(df, options = list(dom = 't'), rownames = FALSE, class = "table-sm")
  })

  ## ---- Tab 6: Scenario Comparison ----
  comp_data <- eventReactive(input$run_compare, {
    scenario_map <- list(
      S1 = list(drug = "none", label = "S1: untreated WD"),
      S2 = list(drug = "DPA",  label = "S2: DPA 500mg TID"),
      S3 = list(drug = "Zinc", label = "S3: Zinc 50mg TID"),
      S4 = list(drug = "Trientine", label = "S4: Trientine 500mg TID"),
      S5 = list(drug = "ALXN1840",  label = "S5: ALXN1840 15mg QD"),
      S7 = list(drug = "ALXN1840",  label = "S7: ALXN1840+Trientine"),
      S8 = list(drug = "none",      label = "S8: normal WT control")
    )

    results <- list()
    for (sid in input$compare_scenarios) {
      if (!sid %in% names(scenario_map)) next
      sc <- scenario_map[[sid]]
      atp7b <- if (sid == "S8") 1.0 else input$atp7b_f / 100
      init_hep <- if (sid == "S8") 15 else input$init_hep_cu
      init_ncbc <- if (sid == "S8") 5  else input$init_ncbc
      init_neuro <- if (sid == "S8") 0 else input$init_neuro

      res <- run_sim(wd_mod, atp7b_func = atp7b, years = input$comp_years,
                     drug = sc$drug, init_CU_HEP = init_hep,
                     init_NCBC = init_ncbc, init_Neuro = init_neuro) %>%
        mutate(Scenario = sc$label)
      results[[sid]] <- res
    }
    bind_rows(results)
  })

  output$compare_plot <- renderPlotly({
    req(comp_data())
    yvar <- switch(input$comp_endpoint,
      NCBC = "NCBC", ALT = "ALT_val", Cu_hep = "Cu_hep",
      Fib = "Fib_val", Neuro = "Neuro_val",
      Cu_brain = "Cu_brain_idx", KF = "Cu_corn_idx"
    )
    ytitle <- switch(input$comp_endpoint,
      NCBC = "NCBC (μg/dL)", ALT = "ALT (U/L)", Cu_hep = "Hepatic Cu (μg/g)",
      Fib = "Fibrosis Score", Neuro = "Neuro Index",
      Cu_brain = "Brain Cu Index", KF = "KF Ring Index"
    )
    p <- ggplot(comp_data(), aes_string(x = "Year", y = yvar, color = "Scenario")) +
      geom_line(size = 1) +
      labs(x = "Year", y = ytitle,
           title = paste("Comparison:", input$comp_endpoint)) +
      scale_color_brewer(palette = "Set1") +
      theme_minimal()
    ggplotly(p)
  })

  output$compare_table <- renderDT({
    req(comp_data())
    end_year <- max(comp_data()$Year)
    comp_data() %>%
      filter(abs(Year - end_year) < 0.1) %>%
      group_by(Scenario) %>%
      summarise(
        `NCBC (μg/dL)` = round(mean(NCBC), 1),
        `ALT (U/L)` = round(mean(ALT_val), 1),
        `Hepatic Cu (μg/g)` = round(mean(Cu_hep), 1),
        `Fibrosis` = round(mean(Fib_val), 2),
        `Neuro Index` = round(mean(Neuro_val), 2),
        `UWDRS Proxy` = round(mean(UWDRS_proxy), 1),
        .groups = "drop"
      ) %>%
      datatable(options = list(dom = 't'), rownames = FALSE, class = "table-sm")
  })

  ## ---- Tab 7: Biomarker Explorer ----
  bm_data <- eventReactive(input$run_bm, {
    results <- list()
    for (drug in input$bm_drugs) {
      res <- run_sim(wd_mod, atp7b_func = input$atp7b_f / 100,
                     years = input$bm_years, drug = drug,
                     init_CU_HEP = input$init_hep_cu,
                     init_NCBC = input$init_ncbc,
                     init_Neuro = input$init_neuro) %>%
        mutate(Drug = drug)
      results[[drug]] <- res
    }
    bind_rows(results)
  })

  output$bm_scatter <- renderPlotly({
    req(bm_data())
    p <- ggplot(bm_data(), aes_string(x = input$bm_x, y = input$bm_y, color = "Drug")) +
      geom_path(alpha = 0.7, size = 0.8) +
      labs(title = paste(input$bm_x, "vs", input$bm_y)) +
      scale_color_brewer(palette = "Set2") + theme_minimal()
    ggplotly(p)
  })

  output$bm_timeseries <- renderPlotly({
    req(bm_data())
    df <- bm_data() %>%
      pivot_longer(c(NCBC, Cp_mg, ALT_val, Fib_val, Neuro_val, Cu_brain_idx),
                   names_to = "Marker", values_to = "Value") %>%
      mutate(Marker = recode(Marker,
        NCBC = "NCBC", Cp_mg = "Ceruloplasmin", ALT_val = "ALT",
        Fib_val = "Fibrosis", Neuro_val = "Neurodegeneration", Cu_brain_idx = "Brain Cu"))
    p <- ggplot(df, aes(x = Year, y = Value, color = Drug)) +
      geom_line(alpha = 0.8) +
      facet_wrap(~Marker, scales = "free_y", ncol = 3) +
      labs(title = "WD biomarker panel time series") + theme_minimal(base_size = 10)
    ggplotly(p)
  })

  output$diagnostic_panel <- renderDT({
    df <- data.frame(
      Biomarker = c("NCBC (free Cu)", "Serum Ceruloplasmin", "24h urinary Cu",
                     "Hepatic Cu (biopsy)", "Serum Cp+NCBC (total)", "UWDRS score",
                     "KF Ring (slit lamp)", "ALP/AST ratio"),
      `Normal value` = c("<10 μg/dL", "20-60 mg/dL", "<100 μg/day",
                    "<50 μg/g dw", "70-140 μg/dL", "0 (normal)", "None", ">4"),
      `WD abnormal value` = c(">20 μg/dL", "<20 mg/dL", ">100 μg/day",
                       ">250 μg/g dw", "<70 μg/dL", ">10 (severe cases)", "Present", "<4"),
      `Treatment target` = c("<15 μg/dL", ">20 mg/dL", "100-500 (early treatment)",
                       "<50 μg/g dw", ">70 μg/dL", "<5 improvement", "Disappearance (long term)", ">4"),
      Sensitivity = c("97%", "75%", "90%", "99%", "80%", "Neurological form", "Neurological form 95%", "80%"),
      check.names = FALSE
    )
    datatable(df, options = list(dom = 't', pageLength = 10), rownames = FALSE,
              class = "table-sm")
  })

  ## ---- Tab 8: Model Info ----
  output$param_ref_table <- renderDT({
    df <- data.frame(
      Parameter = c("ATP7B_func (WD)", "k_bil_WT", "NCBC normal range",
                    "Cp synthesis rate", "DPA F", "DPA CL",
                    "Zinc IC50", "ALXN1840 Emax"),
      Value = c("0.05 (5%)", "3.0/day", "5-10 μg/dL",
               "0.15 mg/dL/day", "55%", "3.5 L/hr",
               "2 μg/mL", "98%"),
      Source = c("Ferenci 2019", "Medici 2010",
                "Członkowska 2018 Nat Rev",
                "Gitlin 1988", "Weiss 2010",
                "Weiss 2010", "Brewer 1998",
                "Schilsky ATLAS 2022")
    )
    datatable(df, options = list(dom = 't'), rownames = FALSE, class = "table-sm")
  })
}

## ============================================================
## Run App
## ============================================================
shinyApp(ui = ui, server = server)
