## =============================================================================
## Trigeminal Neuralgia (TN) QSP Shiny App — skeleton
## Requires: shiny, mrgsolve, dplyr, tidyr, ggplot2, DT, plotly (optional)
## Companion model: tn_mrgsolve_model.R
## =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

mod <- mread(here_or_local <- "tn_mrgsolve_model.R")

SCENARIOS <- c(
  "1. Untreated natural history"                 = "untreated",
  "2. Carbamazepine monotherapy (titrated)"       = "cbz_mono",
  "3. Oxcarbazepine monotherapy"                  = "oxc_mono",
  "4. CBZ + Baclofen combination"                 = "cbz_baclofen",
  "5. MVD (microvascular decompression), day 14"  = "mvd_post",
  "6. CBZ-intolerant -> Gabapentin + Pregabalin"  = "gbp_pgb_switch",
  "7. Percutaneous RF rhizotomy, day 30"          = "rf_rhizotomy"
)

ui <- navbarPage(
  title = "Trigeminal Neuralgia (TN) QSP Explorer",

  ## ---- Tab 1: Patient Profile ----
  tabPanel("Patient Profile",
    sidebarLayout(
      sidebarPanel(
        h4("Patient/Disease Parameters"),
        sliderInput("nvc_severity", "Neurovascular Compression (NVC) Severity (0-2)", 0, 2, 1.2, step = 0.1),
        selectInput("branch", "Affected Branch", choices = c("V2 Maxillary", "V3 Mandibular", "V2+V3 Combined", "V1 Ophthalmic (rare)")),
        checkboxInput("secondary_ms", "Secondary TN (with multiple sclerosis)", FALSE),
        hr(),
        actionButton("simulate", "Run Simulation", class = "btn-primary")
      ),
      mainPanel(
        h4("Patient Summary"),
        verbatimTextOutput("patient_summary"),
        plotOutput("baseline_severity_plot")
      )
    )
  ),

  ## ---- Tab 2: PK ----
  tabPanel("PK (Pharmacokinetics)",
    sidebarLayout(
      sidebarPanel(
        selectInput("pk_scenario", "Scenario", choices = SCENARIOS),
        checkboxGroupInput("pk_drugs", "Displayed Drugs",
          choices = c("Carbamazepine" = "CBZ_conc", "CBZ-epoxide" = "CBZ_epox_c",
                      "Oxcarbazepine-MHD" = "OXC_mhd_c", "Baclofen" = "BAC_conc",
                      "Gabapentin" = "GBP_conc", "Pregabalin" = "PGB_conc"),
          selected = c("CBZ_conc", "CBZ_epox_c"))
      ),
      mainPanel(plotOutput("pk_plot", height = 450), DTOutput("pk_table"))
    )
  ),

  ## ---- Tab 3: Pathway PD (Nav channelopathy / central sensitization) ----
  tabPanel("Pathway PD",
    sidebarLayout(
      sidebarPanel(selectInput("pd_scenario", "Scenario", choices = SCENARIOS)),
      mainPanel(
        plotOutput("navblock_plot", height = 300),
        plotOutput("centsens_plot", height = 300)
      )
    )
  ),

  ## ---- Tab 4: Clinical Endpoints ----
  tabPanel("Clinical Endpoints",
    sidebarLayout(
      sidebarPanel(selectInput("clin_scenario", "Scenario", choices = SCENARIOS)),
      mainPanel(
        plotOutput("paroxysm_plot", height = 300),
        plotOutput("pain_nrs_plot", height = 300),
        verbatimTextOutput("bni_summary")
      )
    )
  ),

  ## ---- Tab 5: Scenario Comparison ----
  tabPanel("Scenario Comparison",
    sidebarLayout(
      sidebarPanel(
        checkboxGroupInput("compare_scenarios", "Scenarios to Compare",
          choices = SCENARIOS, selected = c("untreated", "cbz_mono", "mvd_post"))
      ),
      mainPanel(
        plotOutput("compare_paroxysm_plot", height = 350),
        plotOutput("compare_pain_plot", height = 350),
        DTOutput("compare_summary_table")
      )
    )
  ),

  ## ---- Tab 6: Biomarkers/Imaging ----
  tabPanel("Biomarkers",
    mainPanel(
      h4("MRI FIESTA/CISS Neurovascular Compression Findings (Conceptual Indicator)"),
      plotOutput("nav_upreg_plot", height = 300),
      p("In a real app, compression severity and Nav channel expression index would be linked and displayed with imaging biomarkers.")
    )
  ),

  ## ---- Tab 7: Safety ----
  tabPanel("Safety",
    sidebarLayout(
      sidebarPanel(selectInput("safety_scenario", "Scenario", choices = SCENARIOS)),
      mainPanel(
        plotOutput("sodium_plot", height = 300),
        plotOutput("sedation_plot", height = 300),
        helpText("Hyponatraemia warning threshold: Na+ < 135 mEq/L (higher risk with OXC than CBZ).")
      )
    )
  ),

  ## ---- Tab 8: References ----
  tabPanel("References",
    mainPanel(
      h4("Key References"),
      includeMarkdown("tn_references.md")
    )
  )
)

server <- function(input, output, session) {

  run_sim <- function(scenario_key, end_time = 4320) {
    mod2 <- mod %>% param(NVC_SEVERITY = input$nvc_severity)

    dosing <- switch(scenario_key,
      untreated      = ev(time = 0, amt = 0, cmt = "CBZ_GUT"),
      cbz_mono       = ev(time = seq(0, end_time, 12), amt = 200, cmt = "CBZ_GUT", ii = 12, addl = 359),
      oxc_mono       = ev(time = seq(0, end_time, 12), amt = 300, cmt = "OXC_GUT", ii = 12, addl = 359),
      cbz_baclofen   = seq(ev(time = seq(0, end_time, 12), amt = 200, cmt = "CBZ_GUT", ii = 12, addl = 359),
                            ev(time = seq(0, end_time, 8), amt = 10, cmt = "BAC_GUT", ii = 8, addl = 539)),
      mvd_post       = ev(time = seq(0, 312, 12), amt = 200, cmt = "CBZ_GUT", ii = 12, addl = 25),
      gbp_pgb_switch = seq(ev(time = seq(0, end_time, 8), amt = 300, cmt = "GBP_GUT", ii = 8, addl = 539),
                            ev(time = seq(0, end_time, 12), amt = 75, cmt = "PGB_GUT", ii = 12, addl = 359)),
      rf_rhizotomy   = ev(time = 0, amt = 0, cmt = "CBZ_GUT")
    )

    if (scenario_key == "mvd_post") mod2 <- mod2 %>% param(MVD_ON = 1, MVD_TIME = 336)
    if (scenario_key == "rf_rhizotomy") mod2 <- mod2 %>% param(RF_ON = 1, RF_TIME = 720)

    mod2 %>% ev(dosing) %>% mrgsim(end = end_time) %>% as_tibble()
  }

  output$patient_summary <- renderText({
    paste0("NVC severity: ", input$nvc_severity, "\nBranch: ", input$branch,
           "\nSecondary (MS): ", ifelse(input$secondary_ms, "Yes", "No"))
  })

  output$baseline_severity_plot <- renderPlot({
    df <- run_sim("untreated")
    ggplot(df, aes(time/24, PAROX)) + geom_line(color = "firebrick") +
      labs(x = "Day", y = "Attack Frequency (per day)", title = "Untreated Natural Course") + theme_minimal()
  })

  output$pk_plot <- renderPlot({
    df <- run_sim(input$pk_scenario)
    df_long <- df %>% select(time, all_of(input$pk_drugs)) %>%
      pivot_longer(-time, names_to = "analyte", values_to = "conc")
    ggplot(df_long, aes(time/24, conc, color = analyte)) + geom_line() +
      labs(x = "Day", y = "Concentration", color = "") + theme_minimal()
  })
  output$pk_table <- renderDT({ run_sim(input$pk_scenario) %>% select(time, all_of(input$pk_drugs)) %>% head(50) })

  output$navblock_plot <- renderPlot({
    df <- run_sim(input$pd_scenario)
    ggplot(df, aes(time/24, navblock_tot)) + geom_line(color = "steelblue") +
      labs(x = "Day", y = "Nav Channel Block Fraction", title = "Drug-Induced Nav Blockade") + theme_minimal()
  })
  output$centsens_plot <- renderPlot({
    df <- run_sim(input$pd_scenario)
    ggplot(df, aes(time/24, CENTSENS)) + geom_line(color = "purple") +
      labs(x = "Day", y = "Central Sensitization Index") + theme_minimal()
  })

  output$paroxysm_plot <- renderPlot({
    df <- run_sim(input$clin_scenario)
    ggplot(df, aes(time/24, PAROX)) + geom_line(color = "firebrick") +
      labs(x = "Day", y = "Attack Frequency (per day)") + theme_minimal()
  })
  output$pain_nrs_plot <- renderPlot({
    df <- run_sim(input$clin_scenario)
    ggplot(df, aes(time/24, PAIN)) + geom_line(color = "darkorange") +
      labs(x = "Day", y = "Pain NRS (0-10)") + theme_minimal()
  })
  output$bni_summary <- renderText({
    df <- run_sim(input$clin_scenario)
    last_pain <- tail(df$PAIN, 1)
    bni <- cut(last_pain, breaks = c(-1, 0.5, 3, 6, 8, 11), labels = c("I", "II", "III", "IV", "V"))
    paste0("Final Timepoint BNI Pain Scale Approximation: ", bni)
  })

  output$compare_paroxysm_plot <- renderPlot({
    dfs <- lapply(input$compare_scenarios, function(s) run_sim(s) %>% mutate(scenario = s))
    df_all <- bind_rows(dfs)
    ggplot(df_all, aes(time/24, PAROX, color = scenario)) + geom_line() +
      labs(x = "Day", y = "Attack Frequency (per day)", color = "Scenario") + theme_minimal()
  })
  output$compare_pain_plot <- renderPlot({
    dfs <- lapply(input$compare_scenarios, function(s) run_sim(s) %>% mutate(scenario = s))
    df_all <- bind_rows(dfs)
    ggplot(df_all, aes(time/24, PAIN, color = scenario)) + geom_line() +
      labs(x = "Day", y = "Pain NRS", color = "Scenario") + theme_minimal()
  })
  output$compare_summary_table <- renderDT({
    dfs <- lapply(input$compare_scenarios, function(s) {
      d <- run_sim(s); tibble(scenario = s, final_paroxysm = tail(d$PAROX, 1),
                               final_pain = tail(d$PAIN, 1), final_Na = tail(d$NA_PLASMA, 1))
    })
    bind_rows(dfs)
  })

  output$nav_upreg_plot <- renderPlot({
    df <- run_sim("untreated")
    ggplot(df, aes(time/24, NAV_UPREG)) + geom_line(color = "brown") +
      labs(x = "Day", y = "Nav Channel Expression/Ectopic Excitability Index") + theme_minimal()
  })

  output$sodium_plot <- renderPlot({
    df <- run_sim(input$safety_scenario)
    ggplot(df, aes(time/24, NA_PLASMA)) + geom_line(color = "darkgreen") +
      geom_hline(yintercept = 135, linetype = "dashed", color = "red") +
      labs(x = "Day", y = "Plasma Na+ (mEq/L)") + theme_minimal()
  })
  output$sedation_plot <- renderPlot({
    df <- run_sim(input$safety_scenario)
    ggplot(df, aes(time/24, SEDATION)) + geom_line(color = "gray30") +
      labs(x = "Day", y = "Drowsiness/Ataxia Score (0-10)") + theme_minimal()
  })
}

shinyApp(ui, server)
