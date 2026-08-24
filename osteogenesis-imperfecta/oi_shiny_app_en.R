## =============================================================================
## Osteogenesis Imperfecta (OI) QSP Shiny App — skeleton
## Requires: shiny, mrgsolve, dplyr, tidyr, ggplot2, DT
## Companion model: oi_mrgsolve_model.R
## =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

mod <- mread("oi_mrgsolve_model.R")

SCENARIOS <- c(
  "1. Natural History (Untreated, type III/IV)"          = "untreated",
  "2. IV Pamidronate Cyclical Therapy (Glorieux 1998)" = "pamidronate",
  "3. IV Zoledronic Acid Twice Yearly"                = "zoledronic",
  "4. Denosumab SC (severe recessive OI)"       = "denosumab",
  "5. Teriparatide SC (Adult Type I Only)"       = "teriparatide",
  "6. Setrusumab (Anti-Sclerostin, Investigational)"     = "setrusumab",
  "7. Fresolimumab (Anti-TGF-b, Investigational)"          = "fresolimumab"
)

ui <- navbarPage(
  title = "Osteogenesis Imperfecta (OI) QSP Explorer",

  ## ---- Tab 1: Patient Profile ----
  tabPanel("Patient Profile",
    sidebarLayout(
      sidebarPanel(
        h4("Patient/Disease Parameters"),
        selectInput("oi_type", "Sillence Classification",
          choices = c("Type I (Mild)" = "0.4", "Type III/IV (Moderate-Severe)" = "1.0",
                      "Type V-VIII (Rare Recessive Forms)" = "1.4")),
        numericInput("age_yr", "Age (years)", value = 8, min = 0, max = 80),
        numericInput("weight_kg", "Weight (kg)", value = 25, min = 3, max = 100),
        checkboxInput("growth_plate_open", "Growth Plates Open (Paediatric)", TRUE),
        hr(),
        actionButton("simulate", "Run Simulation", class = "btn-primary")
      ),
      mainPanel(
        h4("Patient Summary"),
        verbatimTextOutput("patient_summary"),
        plotOutput("baseline_bmd_plot")
      )
    )
  ),

  ## ---- Tab 2: PK ----
  tabPanel("PK (Pharmacokinetics)",
    sidebarLayout(
      sidebarPanel(
        selectInput("pk_scenario", "Scenario", choices = SCENARIOS),
        checkboxGroupInput("pk_drugs", "Displayed Drugs",
          choices = c("Pamidronate" = "PAM_Cp", "Zoledronic acid" = "ZOL_Cp",
                      "Denosumab" = "DMAB_Cp", "Teriparatide" = "TPTD_Cp",
                      "Setrusumab" = "SETRU_Cp", "Fresolimumab" = "FRESO_Cp"),
          selected = c("PAM_Cp", "ZOL_Cp"))
      ),
      mainPanel(plotOutput("pk_plot", height = 450), DTOutput("pk_table"))
    )
  ),

  ## ---- Tab 3: PD Key Metrics (Bone Remodeling Pathway) ----
  tabPanel("PD Key Metrics",
    sidebarLayout(
      sidebarPanel(selectInput("pd_scenario", "Scenario", choices = SCENARIOS)),
      mainPanel(
        plotOutput("sost_rankl_plot", height = 300),
        plotOutput("ob_oc_plot", height = 300)
      )
    )
  ),

  ## ---- Tab 4: Clinical Endpoints ----
  tabPanel("Clinical Endpoints",
    sidebarLayout(
      sidebarPanel(selectInput("clin_scenario", "Scenario", choices = SCENARIOS)),
      mainPanel(
        plotOutput("bmd_plot", height = 300),
        plotOutput("fracture_plot", height = 300),
        verbatimTextOutput("clin_summary")
      )
    )
  ),

  ## ---- Tab 5: Scenario Comparison ----
  tabPanel("Scenario Comparison",
    sidebarLayout(
      sidebarPanel(
        checkboxGroupInput("compare_scenarios", "Scenarios to Compare",
          choices = SCENARIOS, selected = c("untreated", "pamidronate", "zoledronic", "setrusumab"))
      ),
      mainPanel(
        plotOutput("compare_bmd_plot", height = 350),
        plotOutput("compare_fracture_plot", height = 350),
        DTOutput("compare_summary_table")
      )
    )
  ),

  ## ---- Tab 6: Biomarkers (Bone Turnover Markers) ----
  tabPanel("Biomarkers",
    sidebarLayout(
      sidebarPanel(selectInput("bio_scenario", "Scenario", choices = SCENARIOS)),
      mainPanel(
        plotOutput("p1np_ctx_plot", height = 350),
        helpText("P1NP: bone formation marker, CTX: bone resorption marker. OI typically shows a high-turnover pattern.")
      )
    )
  ),

  ## ---- Tab 7: Growth/Skeleton ----
  tabPanel("Growth/Skeleton",
    sidebarLayout(
      sidebarPanel(selectInput("growth_scenario", "Scenario", choices = SCENARIOS)),
      mainPanel(
        plotOutput("height_plot", height = 300),
        helpText("Teriparatide is contraindicated in children with open growth plates due to osteosarcoma risk (a warning is shown if administered).")
      )
    )
  ),

  ## ---- Tab 8: References ----
  tabPanel("References",
    mainPanel(
      h4("Key References"),
      includeMarkdown("oi_references_en.md")
    )
  )
)

server <- function(input, output, session) {

  run_sim <- function(scenario_key, end_time = 24*365*4) {
    mod2 <- mod %>% param(SEVERITY = as.numeric(input$oi_type))

    dosing <- switch(scenario_key,
      untreated    = ev(time = 0, amt = 0, cmt = "PAM_CENT"),
      pamidronate  = ev(time = seq(0, end_time, 24*90), amt = input$weight_kg*1,
                         cmt = "PAM_CENT", rate = input$weight_kg*1/2, ii = 24*90, addl = 15),
      zoledronic   = ev(time = seq(0, end_time, 24*182), amt = input$weight_kg*0.05,
                         cmt = "ZOL_CENT", rate = input$weight_kg*0.05/1, ii = 24*182, addl = 7),
      denosumab    = ev(time = seq(0, end_time, 24*90), amt = input$weight_kg*1,
                         cmt = "DMAB_DEPOT", ii = 24*90, addl = 15),
      teriparatide = ev(time = seq(0, end_time, 24), amt = 20, cmt = "TPTD_DEPOT",
                         ii = 24, addl = 1459),
      setrusumab   = ev(time = seq(0, end_time, 24*28), amt = input$weight_kg*20,
                         cmt = "SETRU_CENT", rate = input$weight_kg*20/1, ii = 24*28, addl = 51),
      fresolimumab = ev(time = seq(0, end_time, 24*28), amt = input$weight_kg*1,
                         cmt = "FRESO_CENT", rate = input$weight_kg*1/1, ii = 24*28, addl = 51)
    )

    mod2 %>% ev(dosing) %>% mrgsim(end = end_time, delta = 24) %>% as_tibble()
  }

  output$patient_summary <- renderText({
    paste0("Sillence Classification Severity: ", input$oi_type, "\nAge: ", input$age_yr, " years\n",
           "Weight: ", input$weight_kg, "kg\nGrowth Plates Open: ", ifelse(input$growth_plate_open, "Yes", "No"))
  })

  output$baseline_bmd_plot <- renderPlot({
    df <- run_sim("untreated")
    ggplot(df, aes(time/365, BMC)) + geom_line(color = "firebrick") +
      labs(x = "Year", y = "BMC Index (Change from Baseline)", title = "Untreated Natural History") + theme_minimal()
  })

  output$pk_plot <- renderPlot({
    df <- run_sim(input$pk_scenario)
    df_long <- df %>% select(time, all_of(input$pk_drugs)) %>%
      pivot_longer(-time, names_to = "analyte", values_to = "conc")
    ggplot(df_long, aes(time/24, conc, color = analyte)) + geom_line() +
      labs(x = "Day", y = "Concentration (mg/L or ug/L)", color = "") + theme_minimal()
  })
  output$pk_table <- renderDT({ run_sim(input$pk_scenario) %>% select(time, all_of(input$pk_drugs)) %>% head(50) })

  output$sost_rankl_plot <- renderPlot({
    df <- run_sim(input$pd_scenario)
    ggplot(df, aes(time/365)) +
      geom_line(aes(y = SOST_eff, color = "Sclerostin (effective)")) +
      geom_line(aes(y = RANKL_free, color = "RANKL (free)")) +
      labs(x = "Year", y = "Relative Level", color = "") + theme_minimal()
  })
  output$ob_oc_plot <- renderPlot({
    df <- run_sim(input$pd_scenario)
    ggplot(df, aes(time/365)) +
      geom_line(aes(y = OB, color = "Osteoblast Activity")) +
      geom_line(aes(y = OC, color = "Osteoclast Activity")) +
      labs(x = "Year", y = "Activity Index", color = "") + theme_minimal()
  })

  output$bmd_plot <- renderPlot({
    df <- run_sim(input$clin_scenario)
    ggplot(df, aes(time/365, BMC)) + geom_line(color = "steelblue") +
      labs(x = "Year", y = "BMC Index (Z-score Change)") + theme_minimal()
  })
  output$fracture_plot <- renderPlot({
    df <- run_sim(input$clin_scenario)
    ggplot(df, aes(time/365, FX_CUM)) + geom_line(color = "darkorange") +
      labs(x = "Year", y = "Cumulative Fracture Count") + theme_minimal()
  })
  output$clin_summary <- renderText({
    df <- run_sim(input$clin_scenario)
    paste0("BMC Change at Year 4: ", round(tail(df$BMC, 1), 2),
           "\nCumulative Fractures at Year 4: ", round(tail(df$FX_CUM, 1), 2))
  })

  output$compare_bmd_plot <- renderPlot({
    dfs <- lapply(input$compare_scenarios, function(s) run_sim(s) %>% mutate(scenario = s))
    df_all <- bind_rows(dfs)
    ggplot(df_all, aes(time/365, BMC, color = scenario)) + geom_line() +
      labs(x = "Year", y = "BMC Index", color = "Scenario") + theme_minimal()
  })
  output$compare_fracture_plot <- renderPlot({
    dfs <- lapply(input$compare_scenarios, function(s) run_sim(s) %>% mutate(scenario = s))
    df_all <- bind_rows(dfs)
    ggplot(df_all, aes(time/365, FX_CUM, color = scenario)) + geom_line() +
      labs(x = "Year", y = "Cumulative Fracture Count", color = "Scenario") + theme_minimal()
  })
  output$compare_summary_table <- renderDT({
    dfs <- lapply(input$compare_scenarios, function(s) {
      d <- run_sim(s); tibble(scenario = s, final_BMC = tail(d$BMC, 1), final_fracture = tail(d$FX_CUM, 1))
    })
    bind_rows(dfs)
  })

  output$p1np_ctx_plot <- renderPlot({
    df <- run_sim(input$bio_scenario)
    ggplot(df, aes(time/365)) +
      geom_line(aes(y = P1NP, color = "P1NP (formation)")) +
      geom_line(aes(y = CTX*50, color = "CTX x50 (resorption)")) +
      labs(x = "Year", y = "Marker Level", color = "") + theme_minimal()
  })

  output$height_plot <- renderPlot({
    df <- run_sim(input$growth_scenario)
    ggplot(df, aes(time/365, HEIGHT_Z)) + geom_line(color = "darkgreen") +
      labs(x = "Year", y = "Height Z-score") + theme_minimal()
  })
}

shinyApp(ui, server)
