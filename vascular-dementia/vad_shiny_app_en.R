## ============================================================
## Vascular Dementia (VaD) — Shiny Interactive QSP Dashboard
## Vascular Dementia Quantitative Systems Pharmacology Interactive Dashboard
## ============================================================
## How to run: shiny::runApp("vad_shiny_app.R")
## ============================================================

library(shiny)
library(shinydashboard)
library(ggplot2)
library(dplyr)
library(tidyr)
library(DT)

## ── Simulation engine (self-contained Euler integrator) ────────────
simulate_VaD <- function(
    MMSE0, WMH0, CBF0, BP0, LDL0, ACh0, Syn0,
    use_AHT, use_APT, use_STATIN, use_AChEI, use_MEM, use_CIL,
    dose_AHT, dose_APT, dose_STATIN, dose_AChEI, dose_MEM, dose_CIL,
    n_days = 730, dt = 1
) {
  # PK steady-state approximations (Cp_ss = F*Dose/CL)
  Cp_AHT    <- if (use_AHT    == 1) dose_AHT    * 0.65 / 12.0 / 24 else 0
  Cp_APT    <- if (use_APT    == 1) dose_APT    * 0.80 / 8.0  / 24 else 0
  Cp_ST     <- if (use_STATIN == 1) dose_STATIN * 0.12 / 200.0/ 24 else 0
  Cp_AChEI  <- if (use_AChEI  == 1) dose_AChEI  * 0.90 / 3.5  / 24 * 15 else 0  # *Kp_brain
  Cp_MEM    <- if (use_MEM    == 1) dose_MEM    * 0.85 / 4.5  / 24 * 8  else 0
  Cp_CIL    <- if (use_CIL    == 1) dose_CIL    * 0.90 / 35.0 / 24 else 0

  # PD effects at steady-state
  E_AHT_BP  <- use_AHT    * 25.0 * (Cp_AHT^1.5)  / (0.8^1.5  + Cp_AHT^1.5)
  E_ST_LDL  <- use_STATIN * 0.50 * Cp_ST          / (1.2      + Cp_ST)
  E_APT_SVD <- use_APT    * 0.35 * Cp_APT         / (0.5      + Cp_APT)
  E_CIL_CBF <- use_CIL    * 0.18 * Cp_CIL         / (0.6      + Cp_CIL)
  E_ST_CBF  <- use_STATIN * 0.10 * Cp_ST          / (2.0      + Cp_ST)
  AChE_inhib<- use_AChEI  * 0.75 * Cp_AChEI       / (0.05     + Cp_AChEI)
  E_MEM     <- use_MEM    * 0.60 * Cp_MEM         / (0.10     + Cp_MEM)
  E_ST_inf  <- use_STATIN * 0.30 * Cp_ST          / (2.0      + Cp_ST)
  E_NIC_CBF <- 0  # not simulated separately here

  # State variables
  BP   <- BP0; LDL  <- LDL0; CBF  <- CBF0
  WMH  <- WMH0; Inf  <- 5.0
  MG   <- 0.25; CYT  <- 3.0; ROS  <- 2.5
  ACh  <- ACh0; Syn  <- Syn0; MMSE <- MMSE0

  # Adjust BP/LDL for drug effect (immediate)
  BP  <- BP  - E_AHT_BP
  LDL <- LDL * (1 - E_ST_LDL)

  # CBF with drug effects
  CBF <- CBF * (1 + E_CIL_CBF + E_ST_CBF)

  # Storage
  out <- data.frame(
    Day = numeric(n_days + 1),
    BP = numeric(n_days + 1),   LDL = numeric(n_days + 1),
    CBF = numeric(n_days + 1),  WMH = numeric(n_days + 1),
    Inf = numeric(n_days + 1),  MG  = numeric(n_days + 1),
    CYT = numeric(n_days + 1),  ROS = numeric(n_days + 1),
    ACh = numeric(n_days + 1),  Syn = numeric(n_days + 1),
    MMSE = numeric(n_days + 1)
  )
  out[1, ] <- c(0, BP, LDL, CBF, WMH, Inf, MG, CYT, ROS, ACh, Syn, MMSE)

  for (i in seq_len(n_days)) {
    # WMH
    svd_rate <- 0.0018 * (1 + 0.0012 * max(BP - 130, 0)) *
                           (1 + 0.0006 * max(LDL - 100, 0)) *
                           (1 - E_APT_SVD) * (1 - E_AHT_BP / 25 * 0.4)
    WMH <- WMH + svd_rate * dt

    # Microinfarcts
    Inf <- Inf + 0.0005 * (WMH / WMH0) * (1 - E_APT_SVD * 0.5) * dt

    # CBF dynamics
    CBF_target <- CBF0 * (1 - 0.15 * WMH / (WMH + 15)) *
                          (1 - 0.1 * max(BP - 130, 0) / 30) *
                          (1 + E_CIL_CBF + E_ST_CBF)
    dCBF <- 0.005 * (CBF_target - CBF)
    CBF <- CBF + dCBF * dt

    # Microglia
    hypox <- max(0.5 * (1 - CBF / 55), 0)
    dMG_in  <- 0.003 * (hypox + 0.2 * max(ROS / 2.5 - 1, 0)) * (1 - MG)
    dMG_out <- 0.002 * MG * (1 + E_ST_inf)
    MG <- max(0, min(1, MG + (dMG_in - dMG_out) * dt))

    # Cytokine
    CYT <- max(0.5, CYT + (0.5 * MG - 0.003 * CYT) * dt)

    # ROS
    ros_prod  <- 0.004 * (1 + 0.5 * (CYT / 3 - 1)) * max(1 - CBF / 60, 0.01)
    ros_clear <- 0.003 * ROS * (1 + E_ST_inf * 0.3)
    ROS <- max(0.5, ROS + (ros_prod - ros_clear) * dt)

    # ACh
    ach_loss <- 0.0001 * (1 + 0.15 * max(1 - CBF / 55, 0)) * ACh
    ach_rest <- AChE_inhib * 0.30 * (1 - ACh)
    ACh <- max(0.1, min(1.0, ACh + (-ach_loss + ach_rest) * dt))

    # Synaptic density
    syn_loss <- 0.00015 * (1 + 0.10 * max(ROS / 2.5 - 1, 0) +
                                0.08 * max(CYT / 3 - 1, 0)) *
                            (1 - E_MEM * 0.3) * Syn
    syn_rest <- 0.00005 * (1 + E_ST_CBF * 0.5 + AChE_inhib * 0.2) *
                            max(0.9 - Syn, 0)
    Syn <- max(0.1, min(1.0, Syn + (-syn_loss + syn_rest) * dt))

    # MMSE
    MMSE_exp <- 8.0 * ACh + 12.0 * Syn -
                0.3 * max(WMH - WMH0, 0) - 0.4 * max(Inf - 5, 0)
    MMSE <- max(0, min(30, MMSE + 0.002 * (MMSE_exp - MMSE) * dt))

    out[i + 1, ] <- c(i, BP, LDL, CBF, WMH, Inf, MG, CYT, ROS, ACh, Syn, MMSE)
  }
  out
}

## ── UI ──────────────────────────────────────────────────────────────

ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "VaD QSP Dashboard — Vascular Dementia"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("1. Patient Profile",     tabName = "profile",  icon = icon("user")),
      menuItem("2. Drug PK/PD",        tabName = "pk",       icon = icon("pills")),
      menuItem("3. Vascular · Cerebral Perfusion",       tabName = "vascular", icon = icon("heart")),
      menuItem("4. Neurobiological Mechanisms", tabName = "neuro",    icon = icon("brain")),
      menuItem("5. Clinical Endpoints",   tabName = "outcomes", icon = icon("chart-line")),
      menuItem("6. Scenario Comparison",     tabName = "scenario", icon = icon("sliders-h")),
      menuItem("7. Biomarkers",        tabName = "biomarker",icon = icon("flask"))
    )
  ),

  dashboardBody(
    tabItems(

      ## ── Tab 1: Patient Profile ────────────────────────────────────
      tabItem(tabName = "profile",
        fluidRow(
          box(title = "Baseline Patient Characteristics", width = 4, status = "primary", solidHeader = TRUE,
            sliderInput("MMSE0",  "Baseline MMSE score (0-30)",      22, min = 10, max = 28, step = 1),
            sliderInput("WMH0",   "Baseline WMH volume (mL)",          8, min = 1,  max = 30, step = 0.5),
            sliderInput("BP0",    "Baseline systolic blood pressure (mmHg)",    145, min = 110,max = 185, step = 5),
            sliderInput("CBF0",   "Baseline CBF (mL/100g/min)",      50, min = 30, max = 70, step = 2),
            sliderInput("LDL0",   "Baseline LDL-C (mg/dL)",        130, min = 60, max = 220, step = 5),
            sliderInput("ACh0",   "Baseline ACh tone (rel, 0-1)",   0.60, min = 0.2, max = 1.0, step = 0.05),
            sliderInput("Syn0",   "Baseline synaptic density (rel, 0-1)",0.75, min = 0.3, max = 1.0, step = 0.05)
          ),
          box(title = "Risk Factor Summary", width = 4, status = "warning", solidHeader = TRUE,
            htmlOutput("riskSummary")
          ),
          box(title = "VaD Diagnostic Criteria (VASCOG 2014)", width = 4, status = "info",
            h5("Key diagnostic features:"),
            tags$ul(
              tags$li("Evidence of neurocognitive impairment (MMSE / MoCA / NTB)"),
              tags$li("Imaging evidence of cerebrovascular disease (MRI)"),
              tags$li("Clinical causal relationship (temporal and spatial association)"),
              tags$li("Exclusion of other causes of dementia")
            ),
            h5("Key subtypes:"),
            tags$ul(
              tags$li(strong("Small vessel disease type"), " — white matter degeneration, lacunar infarcts"),
              tags$li(strong("Strategic infarct type"), " — thalamus · hippocampus · angular gyrus"),
              tags$li(strong("Cortical infarct type"), " — multiple cerebral infarcts"),
              tags$li(strong("Mixed type"), " — Alzheimer's + vascular")
            )
          )
        )
      ),

      ## ── Tab 2: Drug PK/PD ─────────────────────────────────────────
      tabItem(tabName = "pk",
        fluidRow(
          box(title = "Drug Selection & Dosing", width = 4, status = "success", solidHeader = TRUE,
            h5("Vascular risk factor control"),
            checkboxInput("use_AHT",    "Antihypertensive",  TRUE),
            sliderInput("dose_AHT",     "Dose (mg/day)", 10, min = 2.5, max = 20, step = 2.5),
            checkboxInput("use_APT",    "Antiplatelet (Aspirin)",  TRUE),
            sliderInput("dose_APT",     "Dose (mg/day)", 100, min = 75, max = 325, step = 25),
            checkboxInput("use_STATIN", "Statin", TRUE),
            sliderInput("dose_STATIN",  "Dose (mg/day)",  40, min = 10, max = 80, step = 10),
            hr(),
            h5("Symptomatic treatment"),
            checkboxInput("use_AChEI",  "Acetylcholinesterase inhibitor (AChEI)", FALSE),
            sliderInput("dose_AChEI",   "Dose (mg/day)", 10, min = 5, max = 23, step = 5),
            checkboxInput("use_MEM",    "Memantine", FALSE),
            sliderInput("dose_MEM",     "Dose (mg/day)", 20, min = 5, max = 20, step = 5),
            hr(),
            h5("Other"),
            checkboxInput("use_CIL",    "Cilostazol", FALSE),
            sliderInput("dose_CIL",     "Dose (mg/day)", 200, min = 50, max = 200, step = 50)
          ),
          box(title = "Pharmacokinetics: Blood Concentration (Steady State)", width = 8, status = "success",
            plotOutput("pkPlot", height = "450px")
          )
        ),
        fluidRow(
          box(title = "Pharmacodynamic Effect Summary", width = 12, status = "primary",
            tableOutput("pdSummary")
          )
        )
      ),

      ## ── Tab 3: Vascular / Cerebral Perfusion ──────────────────────
      tabItem(tabName = "vascular",
        fluidRow(
          box(title = "Cerebral Perfusion & Vascular Mechanism Simulation", width = 12, status = "danger",
            sliderInput("sim_duration", "Simulation duration (days)", 730, min = 90, max = 1460, step = 90),
            actionButton("run_sim", "Run simulation", class = "btn-primary btn-lg")
          )
        ),
        fluidRow(
          box(title = "Systolic Blood Pressure (mmHg)", width = 6, status = "danger",
            plotOutput("bpPlot", height = "280px")),
          box(title = "Cerebral Blood Flow CBF (mL/100g/min)", width = 6, status = "warning",
            plotOutput("cbfPlot", height = "280px"))
        ),
        fluidRow(
          box(title = "White Matter Hyperintensity (WMH) Progression (mL)", width = 6, status = "warning",
            plotOutput("wmhPlot", height = "280px")),
          box(title = "Cumulative Microinfarcts (n)", width = 6, status = "danger",
            plotOutput("infPlot", height = "280px"))
        )
      ),

      ## ── Tab 4: Neurobiological Mechanisms ─────────────────────────
      tabItem(tabName = "neuro",
        fluidRow(
          box(title = "Neurobiological Mechanism Dynamics", width = 12, status = "info",
            p("Visualises the key neurobiological mechanisms of vascular dementia.")
          )
        ),
        fluidRow(
          box(title = "Microglial (M1) Activation & Cytokine Index", width = 6, status = "warning",
            plotOutput("microglia_plot", height = "280px")),
          box(title = "Oxidative Stress Index (ROS)", width = 6, status = "success",
            plotOutput("rosPlot", height = "280px"))
        ),
        fluidRow(
          box(title = "Acetylcholine Tone & AChE Inhibition Effect", width = 6, status = "info",
            plotOutput("achPlot", height = "280px")),
          box(title = "Synaptic Density Change", width = 6, status = "info",
            plotOutput("synPlot", height = "280px"))
        )
      ),

      ## ── Tab 5: Clinical Endpoints ──────────────────────────────────
      tabItem(tabName = "outcomes",
        fluidRow(
          box(title = "MMSE Score Trajectory", width = 8, status = "success", solidHeader = TRUE,
            plotOutput("mmsePlot", height = "400px")),
          box(title = "Outcome Summary (6·12·24 months)", width = 4, status = "primary",
            tableOutput("outcomeTable"))
        ),
        fluidRow(
          box(title = "Contribution by Cognitive Domain", width = 6, status = "info",
            plotOutput("cogDomainPlot", height = "300px")),
          box(title = "Vascular Dementia Stage Classification", width = 6, status = "warning",
            htmlOutput("stageHTML"))
        )
      ),

      ## ── Tab 6: Scenario Comparison ────────────────────────────────
      tabItem(tabName = "scenario",
        fluidRow(
          box(title = "Comparison of 6 Treatment Scenarios", width = 12, status = "primary", solidHeader = TRUE,
            p("Compares 2-year outcomes across 6 treatment strategies for the baseline patient profile.")
          )
        ),
        fluidRow(
          box(title = "MMSE Score Comparison (2 years)", width = 6, status = "success",
            plotOutput("scenMMSE", height = "380px")),
          box(title = "WMH Progression Comparison (2 years)", width = 6, status = "warning",
            plotOutput("scenWMH", height = "380px"))
        ),
        fluidRow(
          box(title = "Table of 2-Year Outcomes by Scenario", width = 12, status = "primary",
            DTOutput("scenTable"))
        )
      ),

      ## ── Tab 7: Biomarker Panel ────────────────────────────────────
      tabItem(tabName = "biomarker",
        fluidRow(
          box(title = "Biomarker Panel & Treatment Response Prediction", width = 12, status = "info",
            p("Visualises pre- and post-treatment changes in key vascular dementia biomarkers.")
          )
        ),
        fluidRow(
          box(title = "Brain MRI Biomarkers", width = 6, status = "primary",
            plotOutput("mri_biomarker", height = "320px")),
          box(title = "Blood/CSF Biomarker Changes", width = 6, status = "info",
            plotOutput("csf_biomarker", height = "320px"))
        ),
        fluidRow(
          box(title = "Biomarker Normal Reference Range & VaD Range", width = 12, status = "warning",
            tableOutput("biomarkerRef"))
        )
      )

    ) # end tabItems
  ) # end dashboardBody
)

## ── Server ──────────────────────────────────────────────────────────

server <- function(input, output, session) {

  ## Reactive simulation
  sim_result <- eventReactive(input$run_sim, {
    simulate_VaD(
      MMSE0 = input$MMSE0, WMH0 = input$WMH0, CBF0 = input$CBF0,
      BP0 = input$BP0, LDL0 = input$LDL0,
      ACh0 = input$ACh0, Syn0 = input$Syn0,
      use_AHT = as.integer(input$use_AHT),
      use_APT = as.integer(input$use_APT),
      use_STATIN = as.integer(input$use_STATIN),
      use_AChEI = as.integer(input$use_AChEI),
      use_MEM = as.integer(input$use_MEM),
      use_CIL = as.integer(input$use_CIL),
      dose_AHT = input$dose_AHT, dose_APT = input$dose_APT,
      dose_STATIN = input$dose_STATIN, dose_AChEI = input$dose_AChEI,
      dose_MEM = input$dose_MEM, dose_CIL = input$dose_CIL,
      n_days = input$sim_duration
    )
  }, ignoreNULL = FALSE)

  ## ── Tab 1: Risk summary ──────────────────────────────────────────
  output$riskSummary <- renderUI({
    bp_risk  <- if (input$BP0 >= 160) "Very high" else if (input$BP0 >= 140) "High" else "Moderate"
    mmse_cat <- if (input$MMSE0 >= 24) "Minimal" else if (input$MMSE0 >= 18) "Mild" else "Moderate"
    wmh_cat  <- if (input$WMH0 < 5) "Mild (Grade 1)" else if (input$WMH0 < 15) "Moderate (Grade 2-3)" else "Severe (Grade 4)"
    HTML(paste0(
      "<b>Blood pressure:</b> ", input$BP0, " mmHg (", bp_risk, ")<br>",
      "<b>LDL-C:</b> ", input$LDL0, " mg/dL<br>",
      "<b>Cognitive stage:</b> ", mmse_cat, " (MMSE=", input$MMSE0, ")<br>",
      "<b>WMH grade:</b> ", wmh_cat, "<br>",
      "<b>CBF:</b> ", input$CBF0, " mL/100g/min<br>",
      "<b>ACh tone:</b> ", input$ACh0, " (normal = 1.0)<br>",
      "<b>Synaptic density:</b> ", input$Syn0, " (normal = 1.0)"
    ))
  })

  ## ── Tab 2: PK Plot ───────────────────────────────────────────────
  output$pkPlot <- renderPlot({
    t <- seq(0, 48, 0.5)  # 48h PK profile
    drugs <- list()
    if (input$use_AHT)    drugs[["Antihypertensive (ARB/ACEi)"]] <-
      input$dose_AHT * 0.65 / 12 * exp(-0.25 * t) + input$dose_AHT * 0.65 / 12 * 0.3 * exp(-0.08 * t)
    if (input$use_APT)    drugs[["Antiplatelet (ASA)"]] <-
      input$dose_APT * 0.80 / 8  * exp(-0.35 * t)
    if (input$use_STATIN) drugs[["Statin"]] <-
      input$dose_STATIN * 0.12 / 200 * exp(-0.20 * t)
    if (input$use_AChEI)  drugs[["AChEI (donepezil)"]] <-
      input$dose_AChEI * 0.90 / 3.5 * exp(-0.029 * t)
    if (input$use_MEM)    drugs[["Memantine"]] <-
      input$dose_MEM   * 0.85 / 4.5 * exp(-0.033 * t)
    if (input$use_CIL)    drugs[["Cilostazol"]] <-
      input$dose_CIL   * 0.90 / 35  * exp(-0.15 * t)
    if (length(drugs) == 0) {
      plot(0, type = "n", xlab = "Time (h)", ylab = "Cp (μg/mL)", main = "No drug selected")
      return()
    }
    df <- do.call(rbind, lapply(names(drugs), function(nm)
      data.frame(t = t, Cp = pmax(drugs[[nm]], 0), Drug = nm)))
    ggplot(df, aes(x = t, y = Cp, color = Drug)) +
      geom_line(size = 1.2) +
      labs(x = "Time (h)", y = "Cp (μg/mL)", title = "Pharmacokinetics: Plasma Concentration-Time Profile (48h)") +
      theme_bw() + theme(legend.position = "bottom", legend.text = element_text(size = 9))
  })

  output$pdSummary <- renderTable({
    Cp_AHT    <- if (input$use_AHT)    input$dose_AHT    * 0.65 / 12.0 / 24 else 0
    Cp_APT    <- if (input$use_APT)    input$dose_APT    * 0.80 / 8.0  / 24 else 0
    Cp_ST     <- if (input$use_STATIN) input$dose_STATIN * 0.12 / 200.0/ 24 else 0
    Cp_AChEI  <- if (input$use_AChEI)  input$dose_AChEI  * 0.90 / 3.5  / 24 * 15 else 0
    Cp_MEM    <- if (input$use_MEM)    input$dose_MEM    * 0.85 / 4.5  / 24 * 8  else 0
    Cp_CIL    <- if (input$use_CIL)    input$dose_CIL    * 0.90 / 35.0 / 24 else 0

    E_AHT_BP <- round(input$use_AHT * 25 * Cp_AHT^1.5 / (0.8^1.5 + Cp_AHT^1.5), 1)
    E_LDL    <- round(input$use_STATIN * 100 * 0.50 * Cp_ST / (1.2 + Cp_ST), 1)
    E_AChE   <- round(input$use_AChEI * 100 * 0.75 * Cp_AChEI / (0.05 + Cp_AChEI), 1)
    E_NMDA   <- round(input$use_MEM * 100 * 0.60 * Cp_MEM / (0.10 + Cp_MEM), 1)
    E_CBF    <- round((input$use_CIL * 18 * Cp_CIL / (0.6 + Cp_CIL)) +
                      (input$use_STATIN * 10 * Cp_ST / (2.0 + Cp_ST)), 1)

    data.frame(
      Drug = c("Antihypertensive (ARB/ACEi)", "Statin", "AChEI (donepezil)", "Memantine", "Cilostazol"),
      `Used` = c(input$use_AHT, input$use_STATIN, input$use_AChEI, input$use_MEM, input$use_CIL),
      `Effect_Metric` = c("SBP reduction (mmHg)", "LDL reduction (%)", "AChE inhibition (%)", "NMDA blockade (%)", "CBF increase (%)"),
      `Predicted_Effect` = c(E_AHT_BP, E_LDL, E_AChE, E_NMDA, E_CBF),
      stringsAsFactors = FALSE
    )
  })

  ## ── Vascular Tab ─────────────────────────────────────────────────
  output$bpPlot <- renderPlot({
    df <- sim_result()
    ggplot(df, aes(x = Day, y = BP)) +
      geom_line(color = "#C62828", size = 1.2) +
      geom_hline(yintercept = 130, linetype = "dashed", color = "steelblue") +
      annotate("text", x = 5, y = 131.5, label = "Target: 130 mmHg",
               hjust = 0, color = "steelblue", size = 3.5) +
      labs(x = "Day", y = "SBP (mmHg)", title = "Systolic Blood Pressure Trend") +
      theme_bw()
  })

  output$cbfPlot <- renderPlot({
    df <- sim_result()
    ggplot(df, aes(x = Day, y = CBF)) +
      geom_line(color = "#1565C0", size = 1.2) +
      geom_hline(yintercept = input$CBF0, linetype = "dashed", color = "gray50") +
      labs(x = "Day", y = "CBF (mL/100g/min)", title = "Change in Cerebral Blood Flow") +
      theme_bw()
  })

  output$wmhPlot <- renderPlot({
    df <- sim_result()
    ggplot(df, aes(x = Day, y = WMH)) +
      geom_area(fill = "#FF8A65", alpha = 0.4) +
      geom_line(color = "#BF360C", size = 1.2) +
      labs(x = "Day", y = "WMH Volume (mL)", title = "White Matter Hyperintensity (WMH) Progression") +
      theme_bw()
  })

  output$infPlot <- renderPlot({
    df <- sim_result()
    ggplot(df, aes(x = Day, y = Inf)) +
      geom_line(color = "#6A1B9A", size = 1.2) +
      labs(x = "Day", y = "Microinfarct Count (cumulative)", title = "Cumulative Cortical Microinfarcts") +
      theme_bw()
  })

  ## ── Neuro Tab ────────────────────────────────────────────────────
  output$microglia_plot <- renderPlot({
    df <- sim_result()
    df2 <- df %>% select(Day, MG, CYT) %>%
      pivot_longer(-Day, names_to = "Marker", values_to = "Value") %>%
      mutate(Marker = recode(Marker, "MG" = "M1 microglia (0-1)", "CYT" = "Cytokine index (AU)"))
    ggplot(df2, aes(x = Day, y = Value, color = Marker)) +
      geom_line(size = 1.1) +
      facet_wrap(~Marker, scales = "free_y", ncol = 1) +
      labs(x = "Day", y = "State value", title = "Neuroinflammatory Mechanism") +
      theme_bw() + theme(legend.position = "none")
  })

  output$rosPlot <- renderPlot({
    df <- sim_result()
    ggplot(df, aes(x = Day, y = ROS)) +
      geom_line(color = "#2E7D32", size = 1.2) +
      geom_hline(yintercept = 1.0, linetype = "dashed", color = "gray50") +
      annotate("text", x = 5, y = 1.05, label = "Normal = 1.0", hjust = 0, size = 3.5) +
      labs(x = "Day", y = "ROS Index (AU)", title = "Oxidative Stress Index") +
      theme_bw()
  })

  output$achPlot <- renderPlot({
    df <- sim_result()
    ggplot(df, aes(x = Day, y = ACh)) +
      geom_line(color = "#1565C0", size = 1.2) +
      geom_hline(yintercept = 1.0, linetype = "dashed", color = "gray50") +
      scale_y_continuous(limits = c(0, 1.05)) +
      labs(x = "Day", y = "ACh Tone (rel, normal=1)", title = "Acetylcholine Tone") +
      theme_bw()
  })

  output$synPlot <- renderPlot({
    df <- sim_result()
    ggplot(df, aes(x = Day, y = Syn)) +
      geom_area(fill = "#7B1FA2", alpha = 0.3) +
      geom_line(color = "#4A148C", size = 1.2) +
      scale_y_continuous(limits = c(0, 1.05)) +
      labs(x = "Day", y = "Synaptic Density (rel)", title = "Synaptic Density Change") +
      theme_bw()
  })

  ## ── Outcomes Tab ─────────────────────────────────────────────────
  output$mmsePlot <- renderPlot({
    df <- sim_result()
    ggplot(df, aes(x = Day, y = MMSE)) +
      geom_ribbon(aes(ymin = MMSE - 1.5, ymax = MMSE + 1.5), fill = "#A5D6A7", alpha = 0.3) +
      geom_line(color = "#1B5E20", size = 1.5) +
      geom_hline(yintercept = input$MMSE0, linetype = "dashed", color = "gray60") +
      geom_hline(yintercept = 10, linetype = "dotted", color = "red", size = 0.8) +
      annotate("text", x = 5, y = 10.5, label = "Moderate dementia threshold (10)", hjust = 0, color = "red", size = 3) +
      scale_y_continuous(limits = c(0, 32)) +
      labs(x = "Day", y = "MMSE Score (0-30)", title = "MMSE Score Trajectory",
           subtitle = paste0("Baseline: ", input$MMSE0)) +
      theme_bw(base_size = 14)
  })

  output$outcomeTable <- renderTable({
    df <- sim_result()
    tpts <- c(180, 365, 730)
    df %>% filter(Day %in% tpts) %>%
      mutate(
        `Timepoint` = paste0(round(Day / 30.4, 0), " months"),
        `MMSE` = round(MMSE, 1),
        `Change` = round(MMSE - input$MMSE0, 1),
        `WMH(mL)` = round(WMH, 1),
        `CBF` = round(CBF, 1)
      ) %>%
      select(`Timepoint`, `MMSE`, `Change`, `WMH(mL)`, `CBF`)
  })

  output$cogDomainPlot <- renderPlot({
    df <- sim_result()
    last_row <- tail(df, 1)
    domains <- data.frame(
      Domain = c("Episodic memory", "Executive function", "Processing speed", "Language function", "Visuospatial function"),
      Score  = c(
        last_row$ACh * 0.9,
        last_row$Syn * 0.8 * (last_row$CBF / input$CBF0),
        last_row$CBF / input$CBF0 * 0.95,
        last_row$Syn * 0.85,
        last_row$Syn * 0.80 * (last_row$CBF / input$CBF0)
      )
    )
    ggplot(domains, aes(x = Domain, y = Score, fill = Domain)) +
      geom_col() +
      geom_hline(yintercept = 1.0, linetype = "dashed", color = "gray40") +
      scale_y_continuous(limits = c(0, 1.1)) +
      labs(x = "", y = "Relative Function (normal=1)", title = "Cognitive Domain Levels at 2 Years") +
      theme_bw() + theme(legend.position = "none",
                         axis.text.x = element_text(angle = 20, hjust = 1, size = 9))
  })

  output$stageHTML <- renderUI({
    df <- sim_result()
    last_mmse <- tail(df$MMSE, 1)
    stage <- if (last_mmse >= 24) "MCI / minimal"
             else if (last_mmse >= 18) "Mild dementia"
             else if (last_mmse >= 10) "Moderate dementia"
             else "Severe dementia"
    color <- if (last_mmse >= 24) "success" else if (last_mmse >= 18) "info"
              else if (last_mmse >= 10) "warning" else "danger"
    HTML(paste0(
      "<div class='alert alert-", color, "'>",
      "<strong>Stage at 2 years: ", stage, "</strong><br>",
      "MMSE: ", round(last_mmse, 1), " (baseline: ", input$MMSE0, ")<br>",
      "Change: ", round(last_mmse - input$MMSE0, 1), " points<br>",
      "WMH: ", round(tail(df$WMH, 1), 1), " mL<br>",
      "CBF: ", round(tail(df$CBF, 1), 1), " mL/100g/min",
      "</div>"
    ))
  })

  ## ── Scenario Tab ─────────────────────────────────────────────────
  scenario_data <- reactive({
    base_args <- list(
      MMSE0 = input$MMSE0, WMH0 = input$WMH0, CBF0 = input$CBF0,
      BP0 = input$BP0, LDL0 = input$LDL0,
      ACh0 = input$ACh0, Syn0 = input$Syn0, n_days = 730
    )
    scen_list <- list(
      list(label = "1. No treatment", use_AHT=0,use_APT=0,use_STATIN=0,use_AChEI=0,use_MEM=0,use_CIL=0),
      list(label = "2. Antihypertensive alone", use_AHT=1,use_APT=0,use_STATIN=0,use_AChEI=0,use_MEM=0,use_CIL=0),
      list(label = "3. Vascular combination\n(antihypertensive+antiplatelet+statin)", use_AHT=1,use_APT=1,use_STATIN=1,use_AChEI=0,use_MEM=0,use_CIL=0),
      list(label = "4. Symptomatic treatment\n(AChEI+memantine)", use_AHT=0,use_APT=0,use_STATIN=0,use_AChEI=1,use_MEM=1,use_CIL=0),
      list(label = "5. Comprehensive treatment\n(vascular+symptomatic+cilostazol)", use_AHT=1,use_APT=1,use_STATIN=1,use_AChEI=1,use_MEM=1,use_CIL=1),
      list(label = "6. Optimal+\n(high-dose statin)", use_AHT=1,use_APT=1,use_STATIN=1,use_AChEI=1,use_MEM=1,use_CIL=1)
    )
    dose_args <- list(
      dose_AHT=10, dose_APT=100, dose_STATIN=40, dose_AChEI=10, dose_MEM=20, dose_CIL=200
    )

    bind_rows(lapply(scen_list, function(sc) {
      args <- c(base_args, sc[setdiff(names(sc), "label")], dose_args)
      if (sc$label == "6. Optimal+\n(high-dose statin)") args$dose_STATIN <- 80
      df <- do.call(simulate_VaD, args)
      df$Scenario <- sc$label
      df
    }))
  })

  output$scenMMSE <- renderPlot({
    df <- scenario_data()
    ggplot(df, aes(x = Day, y = MMSE, color = Scenario)) +
      geom_line(size = 1.1) +
      scale_y_continuous(limits = c(max(0, input$MMSE0 - 10), input$MMSE0 + 1)) +
      labs(x = "Day", y = "MMSE Score", title = "MMSE Trend by Scenario") +
      theme_bw() + theme(legend.position = "bottom", legend.text = element_text(size = 8))
  })

  output$scenWMH <- renderPlot({
    df <- scenario_data()
    ggplot(df, aes(x = Day, y = WMH, color = Scenario)) +
      geom_line(size = 1.1) +
      labs(x = "Day", y = "WMH Volume (mL)", title = "WMH Progression by Scenario") +
      theme_bw() + theme(legend.position = "bottom", legend.text = element_text(size = 8))
  })

  output$scenTable <- renderDT({
    df <- scenario_data() %>%
      filter(Day == 730) %>%
      mutate(
        `MMSE (24 months)` = round(MMSE, 1),
        `MMSE change`     = round(MMSE - input$MMSE0, 1),
        `WMH (24 months, mL)` = round(WMH, 1),
        `WMH change`      = round(WMH - input$WMH0, 1),
        `CBF (24 months)`  = round(CBF, 1),
        `SBP`           = round(BP, 1)
      ) %>%
      select(Scenario, `MMSE (24 months)`, `MMSE change`,
             `WMH (24 months, mL)`, `WMH change`, `CBF (24 months)`, `SBP`)
    datatable(df, options = list(pageLength = 10, dom = 't'), rownames = FALSE)
  })

  ## ── Biomarker Tab ────────────────────────────────────────────────
  output$mri_biomarker <- renderPlot({
    df <- sim_result()
    last <- tail(df, 1)
    bio <- data.frame(
      Marker = c("WMH Volume (mL)", "Microinfarct count\n(cumulative)", "CBF\n(mL/100g/min)"),
      Baseline = c(input$WMH0, 5, input$CBF0),
      End_of_study = c(last$WMH, last$Inf, last$CBF)
    ) %>%
      pivot_longer(-Marker, names_to = "Timepoint", values_to = "Value") %>%
      mutate(Timepoint = recode(Timepoint,
                                "Baseline" = "Baseline", "End_of_study" = "24 months"))
    ggplot(bio, aes(x = Marker, y = Value, fill = Timepoint)) +
      geom_col(position = "dodge") +
      scale_fill_manual(values = c("Baseline" = "#90CAF9", "24 months" = "#1565C0")) +
      labs(x = "", y = "Value", title = "MRI Imaging Biomarkers (baseline vs 24 months)") +
      theme_bw() + theme(axis.text.x = element_text(angle = 10, hjust = 1))
  })

  output$csf_biomarker <- renderPlot({
    df <- sim_result()
    last <- tail(df, 1)
    bio <- data.frame(
      Marker = c("ACh tone\n(rel)", "Synaptic density\n(rel)", "M1 microglia\n(0-1)", "ROS index\n(AU)"),
      Baseline = c(input$ACh0, input$Syn0, 0.25, 2.5),
      End_of_study = c(last$ACh, last$Syn, last$MG, last$ROS)
    ) %>%
      pivot_longer(-Marker, names_to = "Timepoint", values_to = "Value") %>%
      mutate(Timepoint = recode(Timepoint,
                                "Baseline" = "Baseline", "End_of_study" = "24 months"))
    ggplot(bio, aes(x = Marker, y = Value, fill = Timepoint)) +
      geom_col(position = "dodge") +
      scale_fill_manual(values = c("Baseline" = "#C8E6C9", "24 months" = "#2E7D32")) +
      labs(x = "", y = "Value", title = "Neurobiological Biomarkers (baseline vs 24 months)") +
      theme_bw() + theme(axis.text.x = element_text(angle = 10, hjust = 1))
  })

  output$biomarkerRef <- renderTable({
    data.frame(
      `Biomarker` = c(
        "WMH Volume (MRI)",
        "Cerebral Blood Flow (CBF, ASL-MRI)",
        "CSF p-Tau 181",
        "CSF Aβ42",
        "Serum NfL",
        "hsCRP",
        "MMSE",
        "CDR-SB"
      ),
      `Normal Reference Range` = c(
        "< 2 mL (Fazekas grade 1-2)",
        "50-80 mL/100g/min",
        "< 80 pg/mL",
        "> 550 pg/mL",
        "< 10 pg/mL (adult)",
        "< 1.0 mg/L",
        "27-30 (normal cognition)",
        "0 (normal CDR)"
      ),
      `VaD Range` = c(
        "5-40+ mL (severe)",
        "30-55 mL/100g/min",
        "80-250 pg/mL",
        "200-500 pg/mL",
        "20-80+ pg/mL",
        "2-5 mg/L",
        "10-24 (mild-moderate)",
        "2-12 (mild-moderate)"
      ),
      `Clinical Significance` = c(
        "WMH ↑ = small vessel disease, mortality risk ↑",
        "CBF ↓ = reduced cerebral blood flow, cognition ↓",
        "p-Tau ↑ = tau pathology, differentiates mixed type",
        "Aβ42 ↓ = amyloid deposition, mixed type",
        "NfL ↑ = marker of axonal injury",
        "CRP ↑ = systemic inflammation, cerebrovascular risk",
        "MMSE ↓ = degree of cognitive decline",
        "CDR-SB ↑ = severity of functional decline"
      ),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })

} # end server

shinyApp(ui, server)
