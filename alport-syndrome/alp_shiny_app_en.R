## =============================================================================
## Alport Syndrome (AS) QSP Shiny App — skeleton
## Requires: shiny, mrgsolve, dplyr, tidyr, ggplot2, DT
## Companion model: alp_mrgsolve_model.R
## =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

mod <- mread("alp_mrgsolve_model.R")

SCENARIOS <- c(
  "1. Natural History - XLAS Male"                          = "natural_xlas_m",
  "2. Natural History - ARAS (Autosomal Recessive)"                = "natural_aras",
  "3. Natural History - ADAS/Heterozygous Thin-BM"               = "natural_adas",
  "4. Ramipril (ACEi), Started After Proteinuria Onset"           = "ramipril_late",
  "5. Ramipril Started Early (Asymptomatic Phase) (EARLY PRO-TECT)"  = "ramipril_early",
  "6. Losartan (ARB) Monotherapy"                            = "losartan",
  "7. Sparsentan (dual ETA/AT1) Monotherapy"                 = "sparsentan",
  "8. Bardoxolone Methyl (Nrf2 Activator)"             = "bardoxolone",
  "9. Lademirsen (RG-012, anti-miR-21 ASO)"           = "lademirsen",
  "10. Combination: Maximal RAAS Blockade + Dapagliflozin + Sparsentan" = "combo_max"
)

GENOTYPES <- c(
  "ADAS - Heterozygous Thin-BM (Mild)"           = "0.4",
  "XLAS Female - Mosaic (Mild-Moderate)"        = "0.55",
  "XLAS Male / ARAS (Standard Severity)"           = "1.0",
  "ARAS - Truncating Variant (Most Severe)"       = "1.15"
)

ui <- navbarPage(
  title = "Alport Syndrome (AS) QSP Explorer",

  ## ---- Tab 1: Patient Profile ----
  tabPanel("Patient Profile",
    sidebarLayout(
      sidebarPanel(
        h4("Patient/Genotype Parameters"),
        selectInput("genotype", "Genotype (Severity)", choices = GENOTYPES, selected = "1.0"),
        numericInput("age_yr", "Age (years)", value = 10, min = 0, max = 60),
        numericInput("horizon_yr", "Simulation Horizon (years)", value = 25, min = 1, max = 50),
        hr(),
        actionButton("simulate", "Run Simulation", class = "btn-primary")
      ),
      mainPanel(
        h4("Patient Summary"),
        verbatimTextOutput("patient_summary"),
        plotOutput("baseline_gbm_plot")
      )
    )
  ),

  ## ---- Tab 2: PK (Drug Concentration) ----
  tabPanel("PK (Pharmacokinetics)",
    sidebarLayout(
      sidebarPanel(
        selectInput("pk_drug", "Select Drug",
          choices = c("Ramipril/Ramiprilat", "Losartan/E-3174", "Sparsentan",
                      "Bardoxolone methyl", "Lademirsen", "Dapagliflozin", "Finerenone")),
        numericInput("dose", "Dose (mg or relative units)", value = 5),
        numericInput("interval_h", "Dosing Interval (h)", value = 24)
      ),
      mainPanel(plotOutput("pk_plot"))
    )
  ),

  ## ---- Tab 3: PD Key Metrics (Kidney) ----
  tabPanel("PD Key Metrics (Kidney)",
    sidebarLayout(
      sidebarPanel(
        selectInput("scenario_pd", "Treatment Scenario", choices = SCENARIOS)
      ),
      mainPanel(
        plotOutput("gbm_podocyte_plot"),
        plotOutput("fibrosis_hemodynamics_plot")
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
        plotOutput("egfr_trajectory_plot"),
        plotOutput("uacr_plot"),
        verbatimTextOutput("esrd_time_summary")
      )
    )
  ),

  ## ---- Tab 5: Scenario Comparison ----
  tabPanel("Scenario Comparison",
    sidebarLayout(
      sidebarPanel(
        checkboxGroupInput("compare_scenarios", "Scenarios to Compare (max 5)",
          choices = SCENARIOS, selected = c("natural_xlas_m", "ramipril_late", "ramipril_early", "combo_max"))
      ),
      mainPanel(
        plotOutput("compare_egfr_plot"),
        DTOutput("compare_table")
      )
    )
  ),

  ## ---- Tab 6: Biomarkers ----
  tabPanel("Biomarkers",
    sidebarLayout(
      sidebarPanel(
        selectInput("biomarker", "Biomarker",
          choices = c("UACR", "GBM Structural Index", "miR-21 Activity", "Fibrosis (IFTA) Index"))
      ),
      mainPanel(plotOutput("biomarker_plot"))
    )
  ),

  ## ---- Tab 7: Extrarenal Lesions (Cochlea/Eye) ----
  tabPanel("Extrarenal Lesions (Cochlea/Eye)",
    sidebarLayout(
      sidebarPanel(
        selectInput("scenario_extra", "Treatment Scenario", choices = SCENARIOS)
      ),
      mainPanel(
        plotOutput("hearing_loss_plot"),
        plotOutput("ocular_score_plot")
      )
    )
  ),

  ## ---- Tab 8: References ----
  tabPanel("References",
    mainPanel(
      includeMarkdown("alp_references_en.md")
    )
  )
)

server <- function(input, output, session) {

  sim_data <- eventReactive(input$simulate, {
    mod %>%
      param(SEVERITY = as.numeric(input$genotype)) %>%
      mrgsim(end = 24 * 365 * input$horizon_yr, delta = 24 * 7) %>%
      as_tibble()
  })

  output$patient_summary <- renderPrint({
    cat("Genotype severity:", input$genotype, "\n")
    cat("Age (yr):", input$age_yr, "\n")
    cat("Simulation horizon (yr):", input$horizon_yr, "\n")
  })

  output$baseline_gbm_plot <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(FX_YEARS, GBM_INTEG)) + geom_line(color = "#2874A6") +
      labs(x = "Years", y = "GBM structural integrity", title = "GBM Integrity — Natural History")
  })

  output$pk_plot <- renderPlot({
    e <- ev(amt = input$dose, cmt = "RAM_GUT", ii = input$interval_h, addl = 30)
    out <- mod %>% ev(e) %>% mrgsim(end = 24 * 40, delta = 1) %>% as_tibble()
    ggplot(out, aes(time, RAM_CENT)) + geom_line(color = "#B03A2E") +
      labs(x = "Time (h)", y = "Plasma amount", title = paste(input$pk_drug, "PK profile"))
  })

  output$gbm_podocyte_plot <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(FX_YEARS)) +
      geom_line(aes(y = GBM_INTEG, color = "GBM integrity")) +
      geom_line(aes(y = PODO_FRAC, color = "Podocyte fraction")) +
      labs(x = "Years", y = "Fraction", color = "", title = "GBM & Podocyte Trajectory")
  })

  output$fibrosis_hemodynamics_plot <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(FX_YEARS, FIBROSIS)) + geom_line(color = "#943126") +
      labs(x = "Years", y = "Fibrosis (IFTA) index", title = "Fibrotic Progression")
  })

  output$egfr_trajectory_plot <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(FX_YEARS, EGFR)) + geom_line(color = "#212F3C") +
      geom_hline(yintercept = 15, linetype = "dashed", color = "red") +
      labs(x = "Years", y = "eGFR (mL/min/1.73m2)", title = "eGFR Trajectory (ESRD threshold = 15)")
  })

  output$uacr_plot <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(FX_YEARS, UACR)) + geom_line(color = "#6E2C00") +
      labs(x = "Years", y = "UACR (mg/g)", title = "Proteinuria Progression")
  })

  output$esrd_time_summary <- renderPrint({
    df <- sim_data()
    esrd_row <- df[df$EGFR <= 15, ]
    if (nrow(esrd_row) > 0) {
      cat("Time to ESRD (eGFR<=15):", round(esrd_row$FX_YEARS[1], 1), "years\n")
    } else {
      cat("ESRD not reached within simulation horizon\n")
    }
  })

  output$compare_egfr_plot <- renderPlot({
    plot.new()
    title("Scenario comparison placeholder — loop mod %>% param(...) %>% ev(...) per scenario")
  })

  output$compare_table <- renderDT({
    datatable(data.frame(Scenario = names(SCENARIOS), Time_to_ESRD_yr = NA))
  })

  output$biomarker_plot <- renderPlot({
    df <- sim_data()
    col <- switch(input$biomarker,
      "UACR" = "UACR", "GBM Structural Index" = "GBM_INTEG",
      "miR-21 Activity" = "MIR21", "Fibrosis (IFTA) Index" = "FIBROSIS")
    ggplot(df, aes(FX_YEARS, .data[[col]])) + geom_line(color = "#1B4F72") +
      labs(x = "Years", y = col, title = paste(input$biomarker, "Trajectory"))
  })

  output$hearing_loss_plot <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(FX_YEARS, HEARING_LOSS)) + geom_line(color = "#117A65") +
      labs(x = "Years", y = "Hearing threshold shift (dB HL)", title = "Progressive Sensorineural Hearing Loss")
  })

  output$ocular_score_plot <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(FX_YEARS, OCULAR_SCORE)) + geom_line(color = "#CA6F1E") +
      labs(x = "Years", y = "Ocular severity index", title = "Lenticonus / Retinopathy Progression")
  })
}

shinyApp(ui, server)
