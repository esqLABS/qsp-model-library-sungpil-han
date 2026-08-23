##############################################################################
# Chronic Hypothyroidism QSP — Interactive Shiny Dashboard
# Chronic Hypothyroidism QSP Interactive Dashboard
#
# Tabs
# ----
# 1. Patient Profile   — Patient profile (disease type, severity, weight)
# 2. HPT Axis          — Hypothalamic-pituitary-thyroid axis simulation
# 3. Drug PK           — Levothyroxine/liothyronine pharmacokinetics
# 4. Organ Effects     — Organ effects (cardiovascular, metabolic, neurological)
# 5. Scenario Compare  — Comparison of 7 treatment scenarios
# 6. Biomarkers & Dose — Biomarker tracking & dose optimisation
##############################################################################

library(shiny)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(shinydashboard)
library(DT)

# ============================================================
# mrgsolve MODEL (inline)
# ============================================================
hypo_code <- '
$PARAM
TRH0=5.0 TSH0=2.0 TT40=100.0 TT30=1.8 rT30=0.45
ksyn_TRH=5.0 kdeg_TRH=1.04 EC50_TRH=5.4 n_TRH=2.0
ksyn_TSH=2.08 kdeg_TSH=0.693
EC50_T4_fb=20.0 EC50_T3_fb=5.4 n_fb=2.0 n_TRH_stim=1.5
ksec_T4=1.4 ksec_T3=0.052 thyroid_cap=1.0
kconv_T4T3=0.0043 kconv_T4rT3=0.0019
kel_T4=0.0041 kel_T3=0.029 kel_rT3=0.145
ff_T4=0.0002 ff_T3=0.003
Ka_LT4=0.35 F_LT4=0.70 V1_LT4=15.0 Vss_LT4=700.0
k12_LT4=0.020 k21_LT4=0.004 CL_LT4=0.50 MW_T4=776.87
Ka_LT3=0.80 F_LT3=0.95 V1_LT3=15.0 CL_LT3=2.5 MW_T3=650.97
HR0=70.0 kout_HR=0.10 LDL0=3.2 kout_LDL=0.005
BMR0=100.0 kout_BMR=0.008 Emax_Sym=95.0 EC50_Sym=5.4 kout_Sym=0.02
kout_BMD=0.0001

$CMT A_TRH A_TSH A_TT4 A_TT3 A_rT3
     A_LT4_gut A_LT4_c A_LT4_p A_LT3_gut A_LT3_c
     Eff_HR Eff_LDL Eff_BMR Eff_Sym Eff_BMD

$INIT A_TRH=5.0 A_TSH=2.0 A_TT4=100.0 A_TT3=1.8 A_rT3=0.45
      A_LT4_gut=0 A_LT4_c=0 A_LT4_p=0 A_LT3_gut=0 A_LT3_c=0
      Eff_HR=70.0 Eff_LDL=3.2 Eff_BMR=100.0 Eff_Sym=5.0 Eff_BMD=0.0

$MAIN
double fT4 = A_TT4 * ff_T4 * 1000.0;
double fT3 = A_TT3 * ff_T3 * 1000.0;
double fT3_base = TT30 * ff_T3 * 1000.0;
double TSH_rel = A_TSH / TSH0;
double T4_sec = ksec_T4 * TSH_rel * thyroid_cap;
double T3_sec = ksec_T3 * TSH_rel * thyroid_cap;
double HR_target = HR0 * (fT3 / fT3_base);
if(HR_target < 35) HR_target = 35;
if(HR_target > 130) HR_target = 130;
double LDL_target = LDL0 * (1.0 + 0.60 * (1.0 - fT3/(5.4+fT3)) / (1.0 - fT3_base/(5.4+fT3_base)));
if(LDL_target < 0.5) LDL_target = 0.5;
double BMR_target = BMR0 * (fT3 / fT3_base);
if(BMR_target < 50) BMR_target = 50;
if(BMR_target > 130) BMR_target = 130;
double fT3_safe = (fT3 < 0.01) ? 0.01 : fT3;
double Sym_target = Emax_Sym*(1.0 - fT3_safe/(EC50_Sym+fT3_safe))/(1.0 - fT3_base/(EC50_Sym+fT3_base));
if(Sym_target < 0) Sym_target = 0;
if(Sym_target > 100) Sym_target = 100;
double BMD_target = -0.5*(fT3/fT3_base - 1.0)*12.0;
if(BMD_target < -8) BMD_target = -8;
if(BMD_target > 2)  BMD_target = 2;

$ODE
double fT4_ = A_TT4*ff_T4*1000.0;
double fT3_ = A_TT3*ff_T3*1000.0;
double Inh_TRH = pow(EC50_TRH,n_TRH)/(pow(EC50_TRH,n_TRH)+pow(fT3_,n_TRH));
dxdt_A_TRH = ksyn_TRH*Inh_TRH - kdeg_TRH*A_TRH;
double Stim_TRH_ = pow(A_TRH/TRH0, n_TRH_stim);
double Inh_T4_ = pow(EC50_T4_fb,n_fb)/(pow(EC50_T4_fb,n_fb)+pow(fT4_,n_fb));
double Inh_T3_ = pow(EC50_T3_fb,n_fb)/(pow(EC50_T3_fb,n_fb)+pow(fT3_,n_fb));
dxdt_A_TSH = ksyn_TSH*Stim_TRH_*Inh_T4_*Inh_T3_ - kdeg_TSH*A_TSH;
double LT4_abs = Ka_LT4*A_LT4_gut*F_LT4/Vss_LT4/MW_T4*1.0e6;
double LT3_abs = Ka_LT3*A_LT3_gut*F_LT3/V1_LT3/MW_T3*1.0e6;
dxdt_A_TT4 = ksec_T4*(A_TSH/TSH0)*thyroid_cap + LT4_abs
             - kconv_T4T3*A_TT4 - kconv_T4rT3*A_TT4 - kel_T4*A_TT4;
dxdt_A_TT3 = ksec_T3*(A_TSH/TSH0)*thyroid_cap + kconv_T4T3*A_TT4 + LT3_abs - kel_T3*A_TT3;
dxdt_A_rT3 = kconv_T4rT3*A_TT4 - kel_rT3*A_rT3;
dxdt_A_LT4_gut = -Ka_LT4*A_LT4_gut;
double k10_ = CL_LT4/V1_LT4;
dxdt_A_LT4_c = Ka_LT4*A_LT4_gut*F_LT4 - (k10_+k12_LT4)*A_LT4_c + k21_LT4*A_LT4_p;
dxdt_A_LT4_p = k12_LT4*A_LT4_c - k21_LT4*A_LT4_p;
dxdt_A_LT3_gut = -Ka_LT3*A_LT3_gut;
dxdt_A_LT3_c = Ka_LT3*A_LT3_gut*F_LT3 - (CL_LT3/V1_LT3)*A_LT3_c;
dxdt_Eff_HR  = kout_HR *(HR_target  - Eff_HR);
dxdt_Eff_LDL = kout_LDL*(LDL_target - Eff_LDL);
dxdt_Eff_BMR = kout_BMR*(BMR_target - Eff_BMR);
dxdt_Eff_Sym = kout_Sym*(Sym_target - Eff_Sym);
dxdt_Eff_BMD = kout_BMD*(BMD_target - Eff_BMD);

$TABLE
double fT4_pmol = A_TT4*ff_T4*1000.0;
double fT3_pmol = A_TT3*ff_T3*1000.0;
double LT4_conc = A_LT4_c/V1_LT4;
double WeightEst = (100.0-Eff_BMR)*0.15;

$CAPTURE fT4_pmol fT3_pmol A_TSH A_TT4 A_TT3 A_rT3
         LT4_conc Eff_HR Eff_LDL Eff_BMR Eff_Sym Eff_BMD WeightEst
'

mod <- mcode("hypo_shiny", hypo_code, quiet = TRUE)

# ============================================================
# UTILITY FUNCTIONS
# ============================================================
make_events <- function(lt4, lt3, n_days) {
  evs <- list()
  if (lt4 > 0) evs[[1]] <- ev(amt = lt4, cmt = "A_LT4_gut", ii = 24, addl = n_days-1, time=0)
  if (lt3 > 0) {
    evs[[2]] <- ev(amt = lt3/2, cmt = "A_LT3_gut", ii = 24, addl = n_days-1, time=8)
    evs[[3]] <- ev(amt = lt3/2, cmt = "A_LT3_gut", ii = 24, addl = n_days-1, time=20)
  }
  if (length(evs) == 0) return(ev(amt=0, cmt=1, time=1e6))
  do.call(c, evs)
}

run_model <- function(lt4, lt3, cap, n_days=365, delta=6) {
  evs <- make_events(lt4, lt3, n_days)
  mrgsim(mod, ev = evs, end = n_days*24, delta = delta,
         param = list(thyroid_cap = cap)) %>%
    as.data.frame() %>%
    mutate(time_weeks = time / 168)
}

theme_hypo <- function() {
  theme_minimal(base_size = 12) +
    theme(legend.position = "bottom",
          plot.title = element_text(face = "bold", size = 11),
          panel.grid.minor = element_blank())
}

ref_band <- function(lo, hi) {
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = lo, ymax = hi,
           fill = "#4caf50", alpha = 0.08)
}

# ============================================================
# UI
# ============================================================
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "Hypothyroidism QSP", titleWidth = 300),

  dashboardSidebar(
    width = 300,
    sidebarMenu(
      menuItem("Patient profile",   tabName = "patient",   icon = icon("user-md")),
      menuItem("HPT axis",          tabName = "hpt",       icon = icon("project-diagram")),
      menuItem("Drug PK",         tabName = "pk",        icon = icon("pills")),
      menuItem("Organ effects",       tabName = "organs",    icon = icon("heartbeat")),
      menuItem("Scenario comparison",   tabName = "scenarios", icon = icon("chart-line")),
      menuItem("Biomarkers · dose", tabName = "biomarker", icon = icon("flask"))
    ),

    hr(),
    h5("  Treatment settings", style = "color:#aaa; padding-left:15px"),
    sliderInput("lt4_dose",  "LT4 dose (μg/day)", 0, 250, 100, step=12.5),
    sliderInput("lt3_dose",  "LT3 dose (μg/day)", 0,  50,   0, step=2.5),
    sliderInput("thyroid_cap", "Residual thyroid function (%)", 0, 100, 0, step=5),
    sliderInput("sim_days",  "Simulation period (days)",  30, 730, 365, step=30),
    actionButton("run_btn", "Run simulation", icon=icon("play"),
                 class="btn-primary btn-block")
  ),

  dashboardBody(
    tabItems(

      # ------ TAB 1: PATIENT PROFILE ------
      tabItem("patient",
        fluidRow(
          box(title="Select disease type", width=6, status="primary",
            selectInput("dis_type", "Type of hypothyroidism",
              choices = c("Hashimoto's thyroiditis (Hashimoto's)" = "hashimoto",
                          "Post-thyroidectomy (Post-Thyroidectomy)" = "post_tx",
                          "Post-radioactive-iodine treatment (Post-RAI)" = "post_rai",
                          "Central hypothyroidism (Central Hypothyroidism)" = "central",
                          "Subclinical hypothyroidism (Subclinical)" = "subclinical",
                          "Drug-induced (Drug-Induced)" = "drug_ind"),
              selected = "post_tx"),
            numericInput("bw", "Weight (kg)", 70, min=40, max=150),
            numericInput("age", "Age (years)", 45, min=18, max=85),
            selectInput("sex", "Sex", choices=c("Female"="F","Male"="M"), selected="F")
          ),
          box(title="Enter initial laboratory results", width=6, status="info",
            numericInput("init_tsh", "Initial TSH (mIU/L)", 25.0, min=0.01, max=200),
            numericInput("init_ft4", "Initial FT4 (pmol/L)", 6.0, min=0.5, max=40),
            numericInput("init_ft3", "Initial FT3 (pmol/L)", 2.5, min=0.5, max=15),
            actionButton("profile_btn", "Apply profile", class="btn-info")
          )
        ),
        fluidRow(
          box(title="Disease-specific pathophysiology summary", width=12, status="warning",
            uiOutput("disease_info")
          )
        ),
        fluidRow(
          valueBoxOutput("vb_tsh"),
          valueBoxOutput("vb_ft4"),
          valueBoxOutput("vb_ft3")
        )
      ),

      # ------ TAB 2: HPT AXIS ------
      tabItem("hpt",
        fluidRow(
          box(title="Hypothalamic-Pituitary-Thyroid (HPT) Axis", width=12, status="primary",
            p("Shows HPT axis hormone changes after LT4/LT3 treatment. Dashed lines are the normal reference range."),
            fluidRow(
              column(6, plotOutput("plot_tsh",   height=300)),
              column(6, plotOutput("plot_ft4",   height=300))
            ),
            fluidRow(
              column(6, plotOutput("plot_ft3",   height=300)),
              column(6, plotOutput("plot_rt3",   height=300))
            )
          )
        ),
        fluidRow(
          box(title="Expected time to reach TSH target", width=6, status="info",
            tableOutput("tsh_target_table")
          ),
          box(title="HPT axis normalisation indicators", width=6, status="success",
            tableOutput("hpt_ss_table")
          )
        )
      ),

      # ------ TAB 3: DRUG PK ------
      tabItem("pk",
        fluidRow(
          box(title="Levothyroxine (LT4) pharmacokinetics", width=6, status="primary",
            plotOutput("plot_lt4_pk", height=350),
            p(strong("LT4 PK parameters: "), "F=0.70 (fasting), Ka=0.35 h⁻¹, Vss=700L, t½=7 days")
          ),
          box(title="Liothyronine (LT3) pharmacokinetics", width=6, status="warning",
            plotOutput("plot_lt3_pk", height=350),
            p(strong("LT3 PK parameters: "), "F=0.95, Ka=0.80 h⁻¹, V1=15L, t½=24h")
          )
        ),
        fluidRow(
          box(title="T4/T3 profile by absorption phase (first 72 hours)", width=6, status="info",
            plotOutput("plot_early_pk", height=300)
          ),
          box(title="LT4 drug-interaction effects", width=6, status="danger",
            tableOutput("drug_interaction_table"),
            p("Ca²⁺: F ↓ 20-40% | Fe²⁺: F ↓ 30-40% | PPI: F ↓ 10-20%",
              style="font-size:11px; color:gray")
          )
        )
      ),

      # ------ TAB 4: ORGAN EFFECTS ------
      tabItem("organs",
        fluidRow(
          box(title="Cardiovascular effects", width=6, status="danger",
            plotOutput("plot_hr",  height=250),
            plotOutput("plot_ldl", height=250)
          ),
          box(title="Metabolic · neurological effects", width=6, status="warning",
            plotOutput("plot_bmr", height=250),
            plotOutput("plot_sym", height=250)
          )
        ),
        fluidRow(
          box(title="Estimated weight change", width=6, status="info",
            plotOutput("plot_weight", height=280)
          ),
          box(title="Bone mineral density tracking (long-term)", width=6, status="success",
            plotOutput("plot_bmd", height=280)
          )
        )
      ),

      # ------ TAB 5: SCENARIO COMPARISON ------
      tabItem("scenarios",
        fluidRow(
          box(title="Comparison of 7 treatment scenarios", width=12, status="primary",
            p("Each scenario shares the same HPT axis and PD parameters, differing only in treatment method."),
            plotOutput("plot_scenarios_tsh", height=350),
            plotOutput("plot_scenarios_ft4", height=350)
          )
        ),
        fluidRow(
          box(title="Steady-state comparison table (after 1 year)", width=12, status="info",
            DTOutput("scenario_table")
          )
        )
      ),

      # ------ TAB 6: BIOMARKERS & DOSE OPTIMIZATION ------
      tabItem("biomarker",
        fluidRow(
          box(title="Dose-response relationship (TSH vs LT4 dose)", width=6, status="primary",
            plotOutput("plot_dose_response", height=350),
            p("Steady-state TSH by LT4 dose in a post-thyroidectomy patient")
          ),
          box(title="Estimated LT4 dose to reach target TSH", width=6, status="success",
            sliderInput("target_tsh", "Target TSH (mIU/L)", 0.1, 4.5, 2.0, step=0.1),
            numericInput("bw_opt", "Weight (kg)", 70, min=40, max=150),
            actionButton("calc_dose", "Calculate optimal dose", class="btn-success"),
            verbatimTextOutput("optimal_dose_text"),
            hr(),
            h4("Recommended monitoring schedule by treatment period"),
            tableOutput("monitoring_schedule")
          )
        ),
        fluidRow(
          box(title="Biomarker panel trend", width=12, status="info",
            plotOutput("plot_biomarker_panel", height=400)
          )
        )
      )
    )
  )
)

# ============================================================
# SERVER
# ============================================================
server <- function(input, output, session) {

  # Reactive simulation
  sim_data <- eventReactive(input$run_btn, {
    cap <- input$thyroid_cap / 100
    run_model(input$lt4_dose, input$lt3_dose, cap, input$sim_days)
  }, ignoreNULL = FALSE)

  # Disease info text
  output$disease_info <- renderUI({
    info <- switch(input$dis_type,
      hashimoto   = "Hashimoto's thyroiditis: autoimmune destruction of thyroid follicles by anti-TPO/anti-Tg antibodies. Progressive thyroid failure. Goitre (early) → atrophy (late).",
      post_tx     = "After thyroidectomy: loss of thyroid tissue prevents endogenous T4/T3 synthesis. Complete levothyroxine replacement is required. Residual function varies with the extent of surgery.",
      post_rai    = "After radioactive iodine (I-131) treatment: radiation destruction of thyroid follicular cells. Hypothyroidism develops within 3-6 months of treatment. Progresses to permanent hypothyroidism.",
      central     = "Central hypothyroidism: insufficient pituitary TSH secretion or hypothalamic TRH deficiency. Low T4 despite normal/low serum TSH. Use FT4, not TSH, to guide LT4 dose adjustment.",
      subclinical = "Subclinical hypothyroidism: elevated TSH (4.5-10 mIU/L), normal FT4. Mostly early Hashimoto's. Symptoms are mild. Consider treatment if TSH>10 or symptoms are present.",
      drug_ind    = "Drug-induced hypothyroidism: amiodarone (iodine excess · D1 inhibition), lithium (inhibits iodination), IFN-α, anti-PD-1 immunotherapy. Reversible after stopping the causative drug."
    )
    div(class="alert alert-info", p(strong(info)))
  })

  # Value boxes
  output$vb_tsh <- renderValueBox({
    tsh <- isolate(input$init_tsh)
    color <- if(tsh < 0.5) "yellow" else if(tsh > 4.5) "red" else "green"
    valueBox(paste0(tsh, " mIU/L"), "Serum TSH", icon=icon("vials"), color=color)
  })
  output$vb_ft4 <- renderValueBox({
    ft4 <- isolate(input$init_ft4)
    color <- if(ft4 < 12) "red" else if(ft4 > 22) "yellow" else "green"
    valueBox(paste0(ft4, " pmol/L"), "Free T4 (FT4)", icon=icon("flask"), color=color)
  })
  output$vb_ft3 <- renderValueBox({
    ft3 <- isolate(input$init_ft3)
    color <- if(ft3 < 3.5) "red" else if(ft3 > 6.5) "yellow" else "green"
    valueBox(paste0(ft3, " pmol/L"), "Free T3 (FT3)", icon=icon("flask"), color=color)
  })

  # ---- HPT AXIS PLOTS ----
  output$plot_tsh <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(x=time_weeks, y=A_TSH)) +
      ref_band(0.5, 4.5) +
      geom_line(color="#d32f2f", size=1.1) +
      geom_hline(yintercept=c(0.5, 4.5), linetype="dashed", color="grey50") +
      labs(title="Serum TSH", x="Time (weeks)", y="TSH (mIU/L)") +
      theme_hypo()
  })
  output$plot_ft4 <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(x=time_weeks, y=fT4_pmol)) +
      ref_band(12, 22) +
      geom_line(color="#1565c0", size=1.1) +
      geom_hline(yintercept=c(12, 22), linetype="dashed", color="grey50") +
      labs(title="Free T4", x="Time (weeks)", y="FT4 (pmol/L)") +
      theme_hypo()
  })
  output$plot_ft3 <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(x=time_weeks, y=fT3_pmol)) +
      ref_band(3.5, 6.5) +
      geom_line(color="#2e7d32", size=1.1) +
      geom_hline(yintercept=c(3.5, 6.5), linetype="dashed", color="grey50") +
      labs(title="Free T3", x="Time (weeks)", y="FT3 (pmol/L)") +
      theme_hypo()
  })
  output$plot_rt3 <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(x=time_weeks, y=A_rT3)) +
      geom_line(color="#6a1b9a", size=1.1) +
      labs(title="Reverse T3", x="Time (weeks)", y="rT3 (nmol/L)") +
      theme_hypo()
  })

  output$tsh_target_table <- renderTable({
    df <- sim_data()
    wk_cross <- df %>% filter(A_TSH <= 4.5 & A_TSH >= 0.5) %>% head(1)
    data.frame(
      Item = c("Reach TSH target (≤4.5 mIU/L)", "LT4 half-life (reach SS)", "Estimated SS TSH"),
      Value   = c(if(nrow(wk_cross)>0) paste0(round(wk_cross$time_weeks[1],1), " wk") else "Not reached within period",
               "6-8 weeks",
               paste0(round(tail(df$A_TSH,14) %>% mean(), 2), " mIU/L"))
    )
  })

  output$hpt_ss_table <- renderTable({
    df <- sim_data() %>% tail(14*4)  # last 2 weeks
    data.frame(
      Metric    = c("TSH (mIU/L)", "FT4 (pmol/L)", "FT3 (pmol/L)", "rT3 (nmol/L)"),
      NormalRange = c("0.5–4.5", "12–22", "3.5–6.5", "0.2–0.5"),
      Simulated = c(round(mean(df$A_TSH), 2),
                   round(mean(df$fT4_pmol), 1),
                   round(mean(df$fT3_pmol), 1),
                   round(mean(df$A_rT3), 2))
    )
  })

  # ---- PK PLOTS ----
  output$plot_lt4_pk <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(x=time_weeks, y=LT4_conc)) +
      geom_line(color="#0277bd", size=1) +
      labs(title="LT4 central compartment concentration", x="Time (weeks)", y="LT4 (μg/L)") +
      theme_hypo()
  })
  output$plot_lt3_pk <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(x=time_weeks)) +
      geom_line(aes(y=fT3_pmol), color="#f57f17", size=1) +
      labs(title="Free T3 (including LT3 contribution)", x="Time (weeks)", y="FT3 (pmol/L)") +
      theme_hypo()
  })
  output$plot_early_pk <- renderPlot({
    df <- sim_data() %>% filter(time <= 72)
    ggplot(df) +
      geom_line(aes(x=time, y=fT4_pmol, color="FT4"), size=1) +
      geom_line(aes(x=time, y=fT3_pmol*3, color="FT3×3"), size=1) +
      scale_color_manual(values=c(FT4="#1565c0", "FT3×3"="#2e7d32")) +
      labs(title="T4/T3 absorption over the first 72 hours", x="Time (h)", y="pmol/L", color=NULL) +
      theme_hypo()
  })
  output$drug_interaction_table <- renderTable({
    data.frame(
      Drug = c("Calcium carbonate (CaCO₃)", "Ferrous sulfate (FeSO₄)", "Aluminium hydroxide", "Cholestyramine",
               "PPI (omeprazole)", "Sucralfate", "Food (regular meal)"),
      AbsorptionReduction = c("20-40%", "30-40%", "15-20%", "30-40%", "10-20%", "15-20%", "20%"),
      RecommendedInterval = c("≥4 hours", "≥4 hours", "≥4 hours", "≥6 hours", "co-administration OK", "≥4 hours", "take fasting")
    )
  })

  # ---- ORGAN EFFECT PLOTS ----
  output$plot_hr <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(x=time_weeks, y=Eff_HR)) +
      ref_band(60, 100) +
      geom_line(color="#c62828", size=1) +
      labs(title="Heart Rate", x="Time (weeks)", y="HR (bpm)") +
      theme_hypo()
  })
  output$plot_ldl <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(x=time_weeks, y=Eff_LDL)) +
      geom_hline(yintercept=3.4, linetype="dashed", color="grey50") +
      geom_line(color="#e65100", size=1) +
      labs(title="LDL cholesterol", x="Time (weeks)", y="LDL (mmol/L)") +
      theme_hypo()
  })
  output$plot_bmr <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(x=time_weeks, y=Eff_BMR)) +
      ref_band(85, 115) +
      geom_line(color="#6a1b9a", size=1) +
      labs(title="Basal metabolic rate (BMR)", x="Time (weeks)", y="BMR (% normal)") +
      theme_hypo()
  })
  output$plot_sym <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(x=time_weeks, y=Eff_Sym)) +
      geom_line(color="#1b5e20", size=1) +
      labs(title="Hypothyroid symptom score (lower is better)", x="Time (weeks)", y="Symptom score (0-100)") +
      theme_hypo()
  })
  output$plot_weight <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(x=time_weeks, y=WeightEst)) +
      geom_hline(yintercept=0, linetype="dashed") +
      geom_line(color="#0277bd", size=1) +
      labs(title="Estimated weight change (kg)", x="Time (weeks)", y="Weight change (kg)") +
      theme_hypo()
  })
  output$plot_bmd <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(x=time_weeks, y=Eff_BMD)) +
      geom_hline(yintercept=0, linetype="dashed") +
      geom_line(color="#004d40", size=1) +
      labs(title="Bone mineral density change (%)", x="Time (weeks)", y="BMD change (%)") +
      theme_hypo()
  })

  # ---- SCENARIO COMPARISON ----
  scenario_results <- reactive({
    scenarios <- list(
      list(lt4=0,   lt3=0,  cap=1.0, lbl="1. Euthyroid"),
      list(lt4=0,   lt3=0,  cap=0.0, lbl="2. Untreated hypothyroidism"),
      list(lt4=100, lt3=0,  cap=0.0, lbl="3. LT4 100 μg/day"),
      list(lt4=50,  lt3=0,  cap=0.3, lbl="4. LT4 50 μg (subclinical)"),
      list(lt4=175, lt3=0,  cap=0.0, lbl="5. LT4 175 μg (overtreatment)"),
      list(lt4=100, lt3=10, cap=0.0, lbl="6. LT4+LT3 combination"),
      list(lt4=125, lt3=0,  cap=0.0, lbl="7. Post-Tx LT4 125 μg")
    )
    bind_rows(lapply(scenarios, function(sc) {
      df <- run_model(sc$lt4, sc$lt3, sc$cap, n_days=365, delta=24)
      df$scenario <- sc$lbl
      df
    }))
  })

  output$plot_scenarios_tsh <- renderPlot({
    df <- scenario_results()
    cols <- c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02","#a6761d")
    ggplot(df, aes(x=time_weeks, y=A_TSH, color=scenario)) +
      ref_band(0.5, 4.5) +
      geom_line(size=0.9) +
      scale_color_manual(values=cols) +
      labs(title="TSH change by scenario", x="Time (weeks)", y="TSH (mIU/L)", color=NULL) +
      coord_cartesian(ylim=c(0, 30)) +
      theme_hypo()
  })
  output$plot_scenarios_ft4 <- renderPlot({
    df <- scenario_results()
    cols <- c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02","#a6761d")
    ggplot(df, aes(x=time_weeks, y=fT4_pmol, color=scenario)) +
      ref_band(12, 22) +
      geom_line(size=0.9) +
      scale_color_manual(values=cols) +
      labs(title="FT4 change by scenario", x="Time (weeks)", y="FT4 (pmol/L)", color=NULL) +
      theme_hypo()
  })
  output$scenario_table <- renderDT({
    df <- scenario_results() %>%
      filter(time > 0.85*365*24) %>%
      group_by(scenario) %>%
      summarise(
        TSH       = round(mean(A_TSH), 2),
        FT4_pmol  = round(mean(fT4_pmol), 1),
        FT3_pmol  = round(mean(fT3_pmol), 1),
        HR_bpm    = round(mean(Eff_HR), 0),
        LDL_mmol  = round(mean(Eff_LDL), 2),
        BMR_pct   = round(mean(Eff_BMR), 1),
        Sym_score = round(mean(Eff_Sym), 1),
        .groups="drop"
      )
    datatable(df, options=list(pageLength=10, dom='t'), rownames=FALSE) %>%
      formatStyle("TSH", backgroundColor = styleInterval(c(0.5, 4.5), c("#ffeb3b","#c8e6c9","#ffcdd2")))
  })

  # ---- BIOMARKER & DOSE ----
  output$plot_dose_response <- renderPlot({
    doses <- seq(25, 225, by=25)
    dr <- lapply(doses, function(d) {
      ev_obj <- ev(amt=d, cmt="A_LT4_gut", ii=24, addl=364)
      df_ss <- mrgsim(mod, ev=ev_obj, end=8760, delta=24,
                      param=list(thyroid_cap=0)) %>%
        as.data.frame() %>% tail(14) %>%
        summarise(dose=d, TSH=mean(A_TSH), FT4=mean(fT4_pmol))
    })
    dr_df <- bind_rows(dr)
    ggplot(dr_df, aes(x=dose)) +
      geom_line(aes(y=TSH), color="#d32f2f", size=1.2) +
      geom_point(aes(y=TSH), color="#d32f2f", size=3) +
      geom_hline(yintercept=c(0.5, 2.5), linetype="dashed", color="grey40") +
      annotate("rect", xmin=-Inf, xmax=Inf, ymin=0.5, ymax=2.5, fill="#4caf50", alpha=0.1) +
      labs(title="LT4 dose-TSH response curve (post-thyroidectomy)",
           x="LT4 daily dose (μg/day)", y="Steady-state TSH (mIU/L)") +
      theme_hypo()
  })

  observeEvent(input$calc_dose, {
    target <- input$target_tsh
    bw_kg  <- input$bw_opt
    doses  <- seq(10, 250, by=5)
    tsh_vals <- sapply(doses, function(d) {
      ev_obj <- ev(amt=d, cmt="A_LT4_gut", ii=24, addl=364)
      mean(tail(mrgsim(mod, ev=ev_obj, end=8760, delta=48,
                       param=list(thyroid_cap=0)) %>%
                  as.data.frame() %>% pull(A_TSH), 7))
    })
    idx <- which.min(abs(tsh_vals - target))
    opt_dose <- doses[idx]
    opt_dose_kg <- round(opt_dose / bw_kg, 2)
    output$optimal_dose_text <- renderText({
      paste0("Target TSH: ", target, " mIU/L\n",
             "Estimated LT4 dose: ", opt_dose, " μg/day\n",
             "Weight-based: ", opt_dose_kg, " μg/kg/day\n",
             "(reference range: 1.2-1.8 μg/kg/day)")
    })
  })

  output$monitoring_schedule <- renderTable({
    data.frame(
      Period          = c("6-8 weeks after starting treatment", "6 months after TSH normalisation", "1 year after stabilisation", "Long-term monitoring"),
      TestItems       = c("TSH, FT4", "TSH, FT4, lipids", "TSH, FT4, FT3, lipids, BMD", "TSH, FT4, lipids"),
      DoseAdjustmentCriteria   = c("TSH>4.5→increase, TSH<0.5→decrease", "Check FT3 if symptoms persist", "Consider adding LT3 if FT4 normal + symptomatic", "Reassess annually"),
      stringsAsFactors = FALSE
    )
  })

  output$plot_biomarker_panel <- renderPlot({
    df <- sim_data()
    df_long <- df %>%
      select(time_weeks, A_TSH, fT4_pmol, fT3_pmol, Eff_HR, Eff_LDL, Eff_BMR, Eff_Sym) %>%
      pivot_longer(-time_weeks, names_to="biomarker", values_to="value") %>%
      mutate(biomarker = recode(biomarker,
        A_TSH="TSH (mIU/L)", fT4_pmol="FT4 (pmol/L)", fT3_pmol="FT3 (pmol/L)",
        Eff_HR="HR (bpm)", Eff_LDL="LDL (mmol/L)", Eff_BMR="BMR (%)",
        Eff_Sym="Symptom score"))
    ggplot(df_long, aes(x=time_weeks, y=value, color=biomarker)) +
      geom_line(size=0.9) +
      facet_wrap(~biomarker, scales="free_y", ncol=4) +
      labs(title="Biomarker panel trend", x="Time (weeks)", y=NULL) +
      theme_hypo() +
      theme(legend.position="none")
  })
}

# ============================================================
# RUN APPLICATION
# ============================================================
shinyApp(ui = ui, server = server)
