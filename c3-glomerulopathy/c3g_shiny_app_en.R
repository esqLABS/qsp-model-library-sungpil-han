## =============================================================================
## C3 Glomerulopathy (C3G: DDD & C3GN) QSP Shiny App — skeleton
## Requires: shiny, mrgsolve, dplyr, tidyr, ggplot2, DT
## Companion model: c3g_mrgsolve_model.R
## =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

mod <- mread("c3g_mrgsolve_model.R")

SCENARIOS <- c(
  "1. Natural History - C3GN Phenotype"                              = "natural_c3gn",
  "2. Natural History - DDD Phenotype (High-Titre C3NeF)"                 = "natural_ddd",
  "3. Iptacopan 200mg BID (Oral Factor B Inhibitor)"           = "iptacopan",
  "4. Pegcetacoplan 1080mg SC Twice Weekly (C3/C3b Inhibitor)"        = "pegcetacoplan",
  "5. Eculizumab 900mg IV q2w (C5 Inhibitor, Off-label)"           = "eculizumab",
  "6. Ravulizumab Weight-Based IV q8w (Next-Generation C5 Inhibitor)"        = "ravulizumab",
  "7. Danicopan 150mg TID Add-on (Eculizumab + Factor D Inhibition)" = "danicopan_addon",
  "8. CFH Loss-of-Function Variant + Iptacopan"                        = "cfh_lof_iptacopan",
  "9. Post-Transplant Recurrence + Pegcetacoplan Prophylaxis"                 = "transplant_peg",
  "10. Iptacopan Discontinued After 2 Years -> Relapse"                        = "iptacopan_withdrawal"
)

PHENOTYPES <- c(
  "C3GN - Mild (SEVERITY=0.7)"                    = "0.7",
  "C3GN - Standard (SEVERITY=1.0)"                     = "1.0",
  "DDD - Severe/Intramembranous Ribbon-like Deposits (SEVERITY=1.35)"     = "1.35"
)

ui <- navbarPage(
  title = "C3 Glomerulopathy (C3G) QSP Explorer",

  ## ---- Tab 1: Patient Profile ----
  tabPanel("Patient Profile",
    sidebarLayout(
      sidebarPanel(
        h4("Patient/Genotype Parameters"),
        selectInput("phenotype", "Phenotype (Severity)", choices = PHENOTYPES, selected = "1.0"),
        numericInput("c3nef_titer", "C3NeF Titre (0=Negative, 1=Standard Positive, 2=High-Titre Persistent)", value = 1.0, min = 0, max = 2, step = 0.1),
        checkboxInput("cfh_lof", "Carries CFH Loss-of-Function Variant", value = FALSE),
        numericInput("age_yr", "Age (years)", value = 25, min = 2, max = 80),
        numericInput("horizon_yr", "Simulation Horizon (years)", value = 10, min = 1, max = 30),
        hr(),
        actionButton("simulate", "Run Simulation", class = "btn-primary")
      ),
      mainPanel(
        h4("Patient Summary"),
        verbatimTextOutput("patient_summary"),
        plotOutput("baseline_ap_plot")
      )
    )
  ),

  ## ---- Tab 2: PK (Pharmacokinetics) ----
  tabPanel("PK (Pharmacokinetics)",
    sidebarLayout(
      sidebarPanel(
        selectInput("pk_drug", "Select Drug",
          choices = c("Iptacopan", "Pegcetacoplan", "Eculizumab", "Ravulizumab", "Danicopan")),
        numericInput("dose", "Dose (mg)", value = 200),
        numericInput("interval_h", "Dosing Interval (h)", value = 12)
      ),
      mainPanel(plotOutput("pk_plot"))
    )
  ),

  ## ---- Tab 3: PD Key Metrics (Complement System) ----
  tabPanel("PD Key Metrics (Complement System)",
    sidebarLayout(
      sidebarPanel(
        selectInput("scenario_pd", "Treatment Scenario", choices = SCENARIOS)
      ),
      mainPanel(
        plotOutput("ap_activity_plot"),
        plotOutput("c3_sc5b9_plot")
      )
    )
  ),

  ## ---- Tab 4: Clinical Endpoints (Kidney) ----
  tabPanel("Clinical Endpoints",
    sidebarLayout(
      sidebarPanel(
        selectInput("scenario_clin", "Treatment Scenario", choices = SCENARIOS)
      ),
      mainPanel(
        plotOutput("egfr_trajectory_plot"),
        plotOutput("upcr_plot"),
        verbatimTextOutput("eskd_time_summary")
      )
    )
  ),

  ## ---- Tab 5: Glomerular Histology ----
  tabPanel("Glomerular Histology",
    sidebarLayout(
      sidebarPanel(
        selectInput("scenario_histo", "Treatment Scenario", choices = SCENARIOS)
      ),
      mainPanel(
        plotOutput("deposit_mesangial_plot"),
        plotOutput("podocyte_fibrosis_plot")
      )
    )
  ),

  ## ---- Tab 6: Scenario Comparison ----
  tabPanel("Scenario Comparison",
    sidebarLayout(
      sidebarPanel(
        checkboxGroupInput("compare_scenarios", "Scenarios to Compare (max 5)",
          choices = SCENARIOS, selected = c("natural_c3gn", "natural_ddd", "iptacopan", "pegcetacoplan", "eculizumab"))
      ),
      mainPanel(
        plotOutput("compare_egfr_plot"),
        DTOutput("compare_table")
      )
    )
  ),

  ## ---- Tab 7: Biomarkers ----
  tabPanel("Biomarkers",
    sidebarLayout(
      sidebarPanel(
        selectInput("biomarker", "Biomarker",
          choices = c("Serum C3", "sC5b-9", "AP Activity Index", "Glomerular Deposit Burden", "Mesangial Activity Index"))
      ),
      mainPanel(plotOutput("biomarker_plot"))
    )
  ),

  ## ---- Tab 8: References ----
  tabPanel("References",
    mainPanel(
      includeMarkdown("c3g_references_en.md")
    )
  )
)

server <- function(input, output, session) {

  sim_data <- eventReactive(input$simulate, {
    mod %>%
      param(SEVERITY = as.numeric(input$phenotype),
            C3NEF_TITER = input$c3nef_titer,
            CFH_LOF = as.numeric(input$cfh_lof)) %>%
      mrgsim(end = 24 * 365 * input$horizon_yr, delta = 24 * 7) %>%
      as_tibble()
  })

  output$patient_summary <- renderPrint({
    cat("Phenotype severity:", input$phenotype, "\n")
    cat("C3NeF titer:", input$c3nef_titer, "\n")
    cat("CFH LOF variant:", input$cfh_lof, "\n")
    cat("Age (yr):", input$age_yr, "\n")
    cat("Simulation horizon (yr):", input$horizon_yr, "\n")
  })

  output$baseline_ap_plot <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(FX_YEARS, AP_ACTIVITY)) + geom_line(color = "#2874A6") +
      labs(x = "Years", y = "AP convertase activity index", title = "Alternative-Pathway Activity — Natural History")
  })

  output$pk_plot <- renderPlot({
    cmt <- switch(input$pk_drug,
      "Iptacopan" = "IPTA_GUT", "Pegcetacoplan" = "PEG_SC",
      "Eculizumab" = "ECU_CENT", "Ravulizumab" = "RAVU_CENT", "Danicopan" = "DANI_GUT")
    e <- ev(amt = input$dose, cmt = cmt, ii = input$interval_h, addl = 30)
    out <- mod %>% ev(e) %>% mrgsim(end = 24 * 60, delta = 1) %>% as_tibble()
    ycol <- switch(input$pk_drug,
      "Iptacopan" = "IPTA_CENT", "Pegcetacoplan" = "PEG_CENT",
      "Eculizumab" = "ECU_CENT", "Ravulizumab" = "RAVU_CENT", "Danicopan" = "DANI_CENT")
    ggplot(out, aes(time, .data[[ycol]])) + geom_line(color = "#B03A2E") +
      labs(x = "Time (h)", y = "Plasma amount", title = paste(input$pk_drug, "PK profile"))
  })

  output$ap_activity_plot <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(FX_YEARS, AP_ACTIVITY)) + geom_line(color = "#1B4F72") +
      labs(x = "Years", y = "AP activity index", title = "Alternative-Pathway Activity Trajectory")
  })

  output$c3_sc5b9_plot <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(FX_YEARS)) +
      geom_line(aes(y = C3_LEVEL, color = "Serum C3 (mg/dL)")) +
      geom_line(aes(y = SC5B9 / 10, color = "sC5b-9 (ng/mL / 10)")) +
      labs(x = "Years", y = "Value", color = "", title = "Serum C3 & sC5b-9 Biomarkers")
  })

  output$egfr_trajectory_plot <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(FX_YEARS, EGFR)) + geom_line(color = "#212F3C") +
      geom_hline(yintercept = 15, linetype = "dashed", color = "red") +
      labs(x = "Years", y = "eGFR (mL/min/1.73m2)", title = "eGFR Trajectory (ESKD threshold = 15)")
  })

  output$upcr_plot <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(FX_YEARS, UPCR)) + geom_line(color = "#6E2C00") +
      labs(x = "Years", y = "UPCR (mg/g)", title = "Proteinuria Progression")
  })

  output$eskd_time_summary <- renderPrint({
    df <- sim_data()
    eskd_row <- df[df$EGFR <= 15, ]
    if (nrow(eskd_row) > 0) {
      cat("Time to ESKD (eGFR<=15):", round(eskd_row$FX_YEARS[1], 1), "years\n")
    } else {
      cat("ESKD not reached within simulation horizon\n")
    }
  })

  output$deposit_mesangial_plot <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(FX_YEARS)) +
      geom_line(aes(y = GLOM_DEPOSIT, color = "Glomerular deposit burden")) +
      geom_line(aes(y = MESANGIAL, color = "Mesangial activity index")) +
      labs(x = "Years", y = "Index (0-1)", color = "", title = "Glomerular Deposition & Mesangial Activity")
  })

  output$podocyte_fibrosis_plot <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(FX_YEARS)) +
      geom_line(aes(y = PODO_FRAC, color = "Podocyte viable fraction")) +
      geom_line(aes(y = FIBROSIS, color = "Fibrosis (IF-TA) index")) +
      labs(x = "Years", y = "Fraction / Index", color = "", title = "Podocyte Loss & Fibrotic Progression")
  })

  output$compare_egfr_plot <- renderPlot({
    plot.new()
    title("Scenario comparison placeholder — loop mod %>% param(...) %>% ev(...) per scenario")
  })

  output$compare_table <- renderDT({
    datatable(data.frame(Scenario = names(SCENARIOS), Time_to_ESKD_yr = NA))
  })

  output$biomarker_plot <- renderPlot({
    df <- sim_data()
    col <- switch(input$biomarker,
      "Serum C3" = "C3_LEVEL", "sC5b-9" = "SC5B9",
      "AP Activity Index" = "AP_ACTIVITY", "Glomerular Deposit Burden" = "GLOM_DEPOSIT",
      "Mesangial Activity Index" = "MESANGIAL")
    ggplot(df, aes(FX_YEARS, .data[[col]])) + geom_line(color = "#1B4F72") +
      labs(x = "Years", y = col, title = paste(input$biomarker, "Trajectory"))
  })
}

shinyApp(ui, server)
