## =============================================================================
## PTSD QSP Shiny App — skeleton
## Requires: shiny, mrgsolve, dplyr, tidyr, ggplot2, DT
## Companion model: ptsd_mrgsolve_model.R
## =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

mod <- mread("ptsd_mrgsolve_model.R")

SCENARIOS <- c(
  "1. Natural History - Moderate Trauma (Untreated)"                       = "natural_moderate",
  "2. Natural History - Severe/Repeated Trauma + FKBP5 + Dissociative Subtype (Worst Prognosis)"  = "natural_severe_dissoc",
  "3. Sertraline 100-200mg/day"                                  = "sertraline",
  "4. Paroxetine 20-50mg/day"                                    = "paroxetine",
  "5. Prazosin 1-15mg qhs (Nightmare Add-on Therapy)"                        = "prazosin_addon",
  "6. Trauma-Focused Psychotherapy Weekly x12 (PE/CPT/EMDR)"                  = "trauma_psychotherapy",
  "7. Ketamine 0.5mg/kg IV x6 (Over 2 Weeks)"                            = "ketamine",
  "8. MDMA-Assisted Psychotherapy (3 Sessions + 12 Preparation/Integration Sessions)"             = "mdma_therapy",
  "9. High Resilience + Mild Trauma (Spontaneous Remission)"                     = "high_resilience",
  "10. Combination: SSRI + Weekly Psychotherapy + Prazosin"                     = "combination_soc"
)

TRAUMA_TYPES <- c(
  "Mild/Single Trauma (SEVERITY=0.6)" = "0.6",
  "Standard/Moderate (Combat, Violence) (SEVERITY=1.0)" = "1.0",
  "Severe/Repeated (Childhood Abuse + Re-traumatization) (SEVERITY=1.5)" = "1.5"
)

ui <- navbarPage(
  title = "PTSD QSP Explorer",

  ## ---- Tab 1: Patient Profile ----
  tabPanel("Patient Profile",
    sidebarLayout(
      sidebarPanel(
        h4("Patient/Risk Factor Parameters"),
        selectInput("trauma_severity", "Trauma Severity", choices = TRAUMA_TYPES, selected = "1.0"),
        checkboxInput("fkbp5_risk", "FKBP5 Risk Allele x Childhood Trauma", value = FALSE),
        checkboxInput("dissociation", "Dissociative Subtype (Peritraumatic Dissociation)", value = FALSE),
        numericInput("resilience", "Resilience Index (0.5=Vulnerable ~ 1.5=High Resilience)", value = 1.0, min = 0.5, max = 1.5, step = 0.1),
        numericInput("age_yr", "Age (years)", value = 32, min = 18, max = 80),
        numericInput("horizon_wk", "Simulation Horizon (weeks)", value = 52, min = 4, max = 260),
        hr(),
        actionButton("simulate", "Run Simulation", class = "btn-primary")
      ),
      mainPanel(
        h4("Patient Summary"),
        verbatimTextOutput("patient_summary"),
        plotOutput("baseline_caps_plot")
      )
    )
  ),

  ## ---- Tab 2: PK (Pharmacokinetics) ----
  tabPanel("PK (Pharmacokinetics)",
    sidebarLayout(
      sidebarPanel(
        selectInput("pk_drug", "Select Drug",
          choices = c("Sertraline", "Paroxetine", "Prazosin", "Ketamine", "MDMA")),
        numericInput("dose", "Dose (mg)", value = 100),
        numericInput("interval_h", "Dosing Interval (h)", value = 24)
      ),
      mainPanel(plotOutput("pk_plot"))
    )
  ),

  ## ---- Tab 3: PD Key Metrics (Fear Circuit/HPA) ----
  tabPanel("PD Key Metrics (Fear Circuit/HPA)",
    sidebarLayout(
      sidebarPanel(
        selectInput("scenario_pd", "Treatment Scenario", choices = SCENARIOS)
      ),
      mainPanel(
        plotOutput("circuit_plot"),
        plotOutput("cortisol_ne_plot")
      )
    )
  ),

  ## ---- Tab 4: Clinical Endpoints ----
  tabPanel("Clinical Endpoints",
    sidebarLayout(
      sidebarPanel(
        selectInput("scenario_clin", "Treatment Scenario", choices = SCENARIOS)
      ),
      mainPanel(
        plotOutput("caps5_trajectory_plot"),
        plotOutput("symptom_clusters_plot"),
        verbatimTextOutput("remission_summary")
      )
    )
  ),

  ## ---- Tab 5: Scenario Comparison ----
  tabPanel("Scenario Comparison",
    sidebarLayout(
      sidebarPanel(
        checkboxGroupInput("compare_scenarios", "Scenarios to Compare (max 5)",
          choices = SCENARIOS, selected = c("natural_moderate", "sertraline", "trauma_psychotherapy", "mdma_therapy", "combination_soc"))
      ),
      mainPanel(
        plotOutput("compare_caps5_plot"),
        DTOutput("compare_table")
      )
    )
  ),

  ## ---- Tab 6: Biomarkers ----
  tabPanel("Biomarkers",
    sidebarLayout(
      sidebarPanel(
        selectInput("biomarker", "Biomarker",
          choices = c("Cortisol", "Locus Coeruleus/NE Tone", "Amygdala Reactivity Index", "vmPFC Tone", "Extinction Memory Strength", "Sleep Disturbance Index"))
      ),
      mainPanel(plotOutput("biomarker_plot"))
    )
  ),

  ## ---- Tab 7: Sleep/Nocturnal Symptoms ----
  tabPanel("Sleep/Nocturnal Symptoms",
    sidebarLayout(
      sidebarPanel(
        selectInput("scenario_sleep", "Treatment Scenario", choices = SCENARIOS)
      ),
      mainPanel(
        plotOutput("sleep_plot")
      )
    )
  ),

  ## ---- Tab 8: References ----
  tabPanel("References",
    mainPanel(
      includeMarkdown("ptsd_references_en.md")
    )
  )
)

server <- function(input, output, session) {

  sim_data <- eventReactive(input$simulate, {
    mod %>%
      param(TRAUMA_SEVERITY = as.numeric(input$trauma_severity),
            FKBP5_RISK = as.numeric(input$fkbp5_risk),
            DISSOCIATION = as.numeric(input$dissociation),
            RESILIENCE = input$resilience) %>%
      mrgsim(end = 24 * 7 * input$horizon_wk, delta = 24) %>%
      as_tibble()
  })

  output$patient_summary <- renderPrint({
    cat("Trauma severity:", input$trauma_severity, "\n")
    cat("FKBP5 risk allele x childhood trauma:", input$fkbp5_risk, "\n")
    cat("Dissociative subtype:", input$dissociation, "\n")
    cat("Resilience index:", input$resilience, "\n")
    cat("Age (yr):", input$age_yr, "\n")
    cat("Simulation horizon (wk):", input$horizon_wk, "\n")
  })

  output$baseline_caps_plot <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(FX_WEEKS, CAPS5)) + geom_line(color = "#2874A6") +
      labs(x = "Weeks", y = "CAPS-5-like total score", title = "Symptom Severity — Natural History")
  })

  output$pk_plot <- renderPlot({
    cmt <- switch(input$pk_drug,
      "Sertraline" = "SERT_GUT", "Paroxetine" = "PAROX_GUT", "Prazosin" = "PRAZ_GUT",
      "Ketamine" = "KET_CENT", "MDMA" = "MDMA_GUT")
    e <- ev(amt = input$dose, cmt = cmt, ii = input$interval_h, addl = 30)
    out <- mod %>% ev(e) %>% mrgsim(end = 24 * 14, delta = 1) %>% as_tibble()
    ycol <- switch(input$pk_drug,
      "Sertraline" = "SERT_CENT", "Paroxetine" = "PAROX_CENT", "Prazosin" = "PRAZ_CENT",
      "Ketamine" = "KET_CENT", "MDMA" = "MDMA_CENT")
    ggplot(out, aes(time, .data[[ycol]])) + geom_line(color = "#B03A2E") +
      labs(x = "Time (h)", y = "Plasma amount", title = paste(input$pk_drug, "PK profile"))
  })

  output$circuit_plot <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(FX_WEEKS)) +
      geom_line(aes(y = AMYG_REACT, color = "Amygdala reactivity")) +
      geom_line(aes(y = VMPFC_TONE, color = "vmPFC inhibitory tone")) +
      geom_line(aes(y = EXT_MEM, color = "Extinction-memory strength")) +
      geom_line(aes(y = FEAR_MEM, color = "Fear-memory strength")) +
      labs(x = "Weeks", y = "Index (0-1)", color = "", title = "Fear-Extinction Circuit Trajectory")
  })

  output$cortisol_ne_plot <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(FX_WEEKS)) +
      geom_line(aes(y = CORTISOL, color = "Serum cortisol (ug/dL)")) +
      geom_line(aes(y = NE_TONE * 20, color = "NE/LC tone (x20 scaled)")) +
      labs(x = "Weeks", y = "Value", color = "", title = "HPA Axis & Noradrenergic Tone")
  })

  output$caps5_trajectory_plot <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(FX_WEEKS, CAPS5)) + geom_line(color = "#212F3C") +
      geom_hline(yintercept = 20, linetype = "dashed", color = "darkgreen") +
      labs(x = "Weeks", y = "CAPS-5-like total score", title = "CAPS-5 Trajectory (remission threshold <=20)")
  })

  output$symptom_clusters_plot <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(FX_WEEKS)) +
      geom_line(aes(y = INTRUSION, color = "B: Intrusion")) +
      geom_line(aes(y = AVOIDANCE, color = "C: Avoidance")) +
      geom_line(aes(y = NEGCOG, color = "D: Negative cognition/mood")) +
      geom_line(aes(y = HYPERAROUSE, color = "E: Hyperarousal")) +
      labs(x = "Weeks", y = "Cluster severity (0-40)", color = "", title = "DSM-5 Symptom Cluster Trajectories")
  })

  output$remission_summary <- renderPrint({
    df <- sim_data()
    rem_row <- df[df$CAPS5 <= 20, ]
    if (nrow(rem_row) > 0) {
      cat("Time to remission (CAPS-5<=20):", round(rem_row$FX_WEEKS[1], 1), "weeks\n")
    } else {
      cat("Remission not reached within simulation horizon\n")
    }
  })

  output$compare_caps5_plot <- renderPlot({
    plot.new()
    title("Scenario comparison placeholder — loop mod %>% param(...) %>% ev(...) per scenario")
  })

  output$compare_table <- renderDT({
    datatable(data.frame(Scenario = names(SCENARIOS), Time_to_Remission_wk = NA))
  })

  output$biomarker_plot <- renderPlot({
    df <- sim_data()
    col <- switch(input$biomarker,
      "Cortisol" = "CORTISOL", "Locus Coeruleus/NE Tone" = "NE_TONE",
      "Amygdala Reactivity Index" = "AMYG_REACT", "vmPFC Tone" = "VMPFC_TONE",
      "Extinction Memory Strength" = "EXT_MEM", "Sleep Disturbance Index" = "SLEEP_DIST")
    ggplot(df, aes(FX_WEEKS, .data[[col]])) + geom_line(color = "#1B4F72") +
      labs(x = "Weeks", y = col, title = paste(input$biomarker, "Trajectory"))
  })

  output$sleep_plot <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(FX_WEEKS, SLEEP_DIST)) + geom_line(color = "#78281F") +
      labs(x = "Weeks", y = "Sleep disturbance index (0-1)", title = "Nightmare/Sleep Disturbance Trajectory")
  })
}

shinyApp(ui, server)
