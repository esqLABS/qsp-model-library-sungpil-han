## =============================================================================
## Retinitis Pigmentosa (RP) QSP Shiny App — skeleton
## Requires: shiny, mrgsolve, dplyr, tidyr, ggplot2, DT
## Companion model: rp_mrgsolve_model.R
## =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

mod <- mread("rp_mrgsolve_model.R")

SCENARIOS <- c(
  "1. Natural history - RHO-adRP"                       = "natural_rho",
  "2. Natural history - RPGR-XLRP"                      = "natural_xlrp",
  "3. Natural history - RPE65-LCA/EOSRD"                = "natural_rpe65",
  "4. Voretigene neparvovec (Luxturna, RPE65)"   = "voretigene",
  "5. Investigational RPGR gene therapy (XLRP)"             = "rpgr_gt",
  "6. MCO-010 optogenetic therapy (end-stage, genotype-independent)" = "mco010",
  "7. CNTF implantable sustained-release implant"                = "cntf",
  "8. N-acetylcysteine (oral antioxidant)"          = "nac",
  "9. Vitamin A palmitate + DHA (oral)"            = "vita",
  "10. Voretigene neparvovec + NAC combination"         = "combo"
)

GENOTYPES <- c(
  "USH2A (mild-to-moderate, slowly progressive)" = "0.75",
  "RHO-autosomal dominant RP (moderate)"     = "1.0",
  "RPGR-X-linked RP (rapid progression)"        = "1.6",
  "RPE65-LCA/early-onset severe retinal dystrophy" = "2.5"
)

ui <- navbarPage(
  title = "Retinitis Pigmentosa (RP) QSP Explorer",

  ## ---- Tab 1: Patient profile ----
  tabPanel("Patient profile",
    sidebarLayout(
      sidebarPanel(
        h4("Patient / genotype parameters"),
        selectInput("genotype", "Genotype (severity)", choices = GENOTYPES, selected = "1.0"),
        numericInput("age_yr", "Age (years)", value = 20, min = 0, max = 90),
        numericInput("horizon_yr", "Simulation horizon (years)", value = 15, min = 1, max = 40),
        hr(),
        actionButton("simulate", "Run simulation", class = "btn-primary")
      ),
      mainPanel(
        h4("Patient summary"),
        verbatimTextOutput("patient_summary"),
        plotOutput("baseline_rodcone_plot")
      )
    )
  ),

  ## ---- Tab 2: PK (gene vector expression · drug concentration) ----
  tabPanel("PK (pharmacokinetics)",
    sidebarLayout(
      sidebarPanel(
        selectInput("pk_scenario", "Scenario", choices = SCENARIOS),
        checkboxGroupInput("pk_analytes", "Displayed items",
          choices = c("RPE65 expression (GT65_EXPR)" = "GT65_EXPR", "RPGR expression (GTRPGR_EXPR)" = "GTRPGR_EXPR",
                      "MCO opsin expression (OP_EXPR)" = "OP_EXPR", "CNTF tissue concentration" = "CNTF_DEV",
                      "NAC plasma concentration" = "NAC_CENT", "Vitamin A store" = "VITA_STORE"),
          selected = c("GT65_EXPR", "CNTF_DEV"))
      ),
      mainPanel(plotOutput("pk_plot", height = 450), DTOutput("pk_table"))
    )
  ),

  ## ---- Tab 3: PD key metrics (photoreceptor survival / oxidative stress) ----
  tabPanel("PD key metrics",
    sidebarLayout(
      sidebarPanel(selectInput("pd_scenario", "Scenario", choices = SCENARIOS)),
      mainPanel(
        plotOutput("rodcone_plot", height = 300),
        plotOutput("ros_microglia_plot", height = 300)
      )
    )
  ),

  ## ---- Tab 4: Clinical endpoints ----
  tabPanel("Clinical endpoints",
    sidebarLayout(
      sidebarPanel(selectInput("clin_scenario", "Scenario", choices = SCENARIOS)),
      mainPanel(
        plotOutput("erg_plot", height = 280),
        plotOutput("vf_bcva_plot", height = 280),
        plotOutput("fst_mlmt_plot", height = 280),
        verbatimTextOutput("clin_summary")
      )
    )
  ),

  ## ---- Tab 5: Scenario comparison ----
  tabPanel("Scenario comparison",
    sidebarLayout(
      sidebarPanel(
        checkboxGroupInput("compare_scenarios", "Scenarios to compare",
          choices = SCENARIOS, selected = c("natural_rpe65", "voretigene", "mco010", "cntf"))
      ),
      mainPanel(
        plotOutput("compare_mlmt_plot", height = 350),
        plotOutput("compare_fst_plot", height = 350),
        DTOutput("compare_summary_table")
      )
    )
  ),

  ## ---- Tab 6: Biomarkers (RGC / macular oedema) ----
  tabPanel("Biomarkers",
    sidebarLayout(
      sidebarPanel(selectInput("bio_scenario", "Scenario", choices = SCENARIOS),
                   checkboxInput("cai_on", "Concomitant CAI therapy (dorzolamide / acetazolamide)", FALSE)),
      mainPanel(
        plotOutput("rgc_plot", height = 300),
        plotOutput("cme_plot", height = 300),
        helpText("RGC_FRAC: retinal ganglion cell survival rate (determines responsiveness to optogenetic/prosthetic therapy). CME_CST: cystoid macular oedema central retinal thickness.")
      )
    )
  ),

  ## ---- Tab 7: Gene therapy / optogenetic response ----
  tabPanel("Gene therapy · optogenetics",
    sidebarLayout(
      sidebarPanel(selectInput("gt_scenario", "Scenario", choices = SCENARIOS,
                                selected = "voretigene")),
      mainPanel(
        plotOutput("visualcycle_plot", height = 300),
        plotOutput("bypass_plot", height = 300),
        helpText("Voretigene neparvovec: visual cycle recovery indicator. MCO-010: bypass photosensitisation signal independent of photoreceptor survival (proportional to RGC survival rate).")
      )
    )
  ),

  ## ---- Tab 8: References ----
  tabPanel("References",
    mainPanel(
      h4("Key references"),
      includeMarkdown("rp_references_en.md")
    )
  )
)

server <- function(input, output, session) {

  run_sim <- function(scenario_key, end_time = NULL) {
    if (is.null(end_time)) end_time <- 24*365*input$horizon_yr

    mod2 <- mod %>% param(SEVERITY = as.numeric(input$genotype))

    is_rpe65_geno <- as.numeric(input$genotype) == 2.5
    is_xlrp_geno  <- as.numeric(input$genotype) == 1.6

    mod2 <- switch(scenario_key,
      natural_rho    = mod2 %>% param(IS_RPE65 = 0, IS_XLRP = 0),
      natural_xlrp   = mod2 %>% param(IS_RPE65 = 0, IS_XLRP = 1),
      natural_rpe65  = mod2 %>% param(IS_RPE65 = 1, IS_XLRP = 0),
      voretigene     = mod2 %>% param(IS_RPE65 = 1, IS_XLRP = 0),
      rpgr_gt        = mod2 %>% param(IS_RPE65 = 0, IS_XLRP = 1),
      mco010         = mod2 %>% param(IS_RPE65 = as.numeric(is_rpe65_geno), IS_XLRP = as.numeric(is_xlrp_geno)),
      cntf           = mod2 %>% param(CNTF_PROD = 0.050),
      nac            = mod2,
      vita           = mod2 %>% param(VITA_KIN = 0.0020),
      combo          = mod2 %>% param(IS_RPE65 = 1, IS_XLRP = 0)
    )
    mod2 <- mod2 %>% param(CAI_ON = as.numeric(input$cai_on))

    dosing <- switch(scenario_key,
      natural_rho    = ev(time = 0, amt = 0, cmt = "GT65_VG"),
      natural_xlrp   = ev(time = 0, amt = 0, cmt = "GT65_VG"),
      natural_rpe65  = ev(time = 0, amt = 0, cmt = "GT65_VG"),
      voretigene     = ev(time = 0, amt = 100, cmt = "GT65_VG"),
      rpgr_gt        = ev(time = 0, amt = 100, cmt = "GTRPGR_VG"),
      mco010         = ev(time = 0, amt = 100, cmt = "OP_VG"),
      cntf           = ev(time = 0, amt = 0, cmt = "GT65_VG"),
      nac            = ev(time = seq(0, end_time, 8), amt = 600, cmt = "NAC_GUT", ii = 8, addl = length(seq(0, end_time, 8))),
      vita           = ev(time = 0, amt = 0, cmt = "GT65_VG"),
      combo          = ev(time = 0, amt = 100, cmt = "GT65_VG")
    )

    out <- mod2 %>% ev(dosing) %>% mrgsim(end = end_time, delta = 24) %>% as_tibble()
    out
  }

  output$patient_summary <- renderText({
    paste0("Genotype/severity coefficient: ", input$genotype, "\nAge: ", input$age_yr, " years\n",
           "Simulation horizon: ", input$horizon_yr, " years")
  })

  output$baseline_rodcone_plot <- renderPlot({
    df <- run_sim("natural_rho")
    ggplot(df, aes(time/(24*365))) +
      geom_line(aes(y = ROD_FRAC, color = "Rod photoreceptor survival rate")) +
      geom_line(aes(y = CONE_FRAC, color = "Cone photoreceptor survival rate")) +
      labs(x = "Year", y = "Survival fraction", color = "", title = "Natural history (untreated)") + theme_minimal()
  })

  output$pk_plot <- renderPlot({
    df <- run_sim(input$pk_scenario)
    df_long <- df %>% select(time, all_of(input$pk_analytes)) %>%
      pivot_longer(-time, names_to = "analyte", values_to = "value")
    ggplot(df_long, aes(time/(24*365), value, color = analyte)) + geom_line() +
      labs(x = "Year", y = "Expression / concentration (relative units)", color = "") + theme_minimal()
  })
  output$pk_table <- renderDT({ run_sim(input$pk_scenario) %>% select(time, all_of(input$pk_analytes)) %>% head(50) })

  output$rodcone_plot <- renderPlot({
    df <- run_sim(input$pd_scenario)
    ggplot(df, aes(time/(24*365))) +
      geom_line(aes(y = ROD_FRAC, color = "Rod photoreceptor survival rate")) +
      geom_line(aes(y = CONE_FRAC, color = "Cone photoreceptor survival rate")) +
      labs(x = "Year", y = "Survival fraction", color = "") + theme_minimal()
  })
  output$ros_microglia_plot <- renderPlot({
    df <- run_sim(input$pd_scenario)
    ggplot(df, aes(time/(24*365))) +
      geom_line(aes(y = ROS, color = "Oxidative stress (ROS)")) +
      geom_line(aes(y = MICROGLIA, color = "Microglial activity")) +
      labs(x = "Year", y = "Relative index", color = "") + theme_minimal()
  })

  output$erg_plot <- renderPlot({
    df <- run_sim(input$clin_scenario)
    ggplot(df, aes(time/(24*365))) +
      geom_line(aes(y = ERG_ROD, color = "ERG rod b-wave")) +
      geom_line(aes(y = ERG_CONE, color = "ERG cone amplitude")) +
      labs(x = "Year", y = "Amplitude (uV)", color = "") + theme_minimal()
  })
  output$vf_bcva_plot <- renderPlot({
    df <- run_sim(input$clin_scenario)
    ggplot(df, aes(time/(24*365))) +
      geom_line(aes(y = VF_AREA/10, color = "Visual field area/10 (deg2)")) +
      geom_line(aes(y = BCVA*500, color = "BCVA x500 (logMAR)")) +
      labs(x = "Year", y = "Value (rescaled)", color = "") + theme_minimal()
  })
  output$fst_mlmt_plot <- renderPlot({
    df <- run_sim(input$clin_scenario)
    ggplot(df, aes(time/(24*365))) +
      geom_line(aes(y = FST, color = "FST (dB, lower is better)")) +
      geom_line(aes(y = MLMT*10, color = "MLMT x10")) +
      labs(x = "Year", y = "Value", color = "") + theme_minimal()
  })
  output$clin_summary <- renderText({
    df <- run_sim(input$clin_scenario)
    paste0(input$horizon_yr, "-year ERG rod amplitude: ", round(tail(df$ERG_ROD, 1), 1), " uV\n",
           input$horizon_yr, "-year visual field area: ", round(tail(df$VF_AREA, 1), 0), " deg2\n",
           input$horizon_yr, "-year BCVA: ", round(tail(df$BCVA, 1), 2), " logMAR\n",
           input$horizon_yr, "-year MLMT: ", round(tail(df$MLMT, 1), 2))
  })

  output$compare_mlmt_plot <- renderPlot({
    dfs <- lapply(input$compare_scenarios, function(s) run_sim(s) %>% mutate(scenario = s))
    df_all <- bind_rows(dfs)
    ggplot(df_all, aes(time/(24*365), MLMT, color = scenario)) + geom_line() +
      labs(x = "Year", y = "MLMT score", color = "Scenario") + theme_minimal()
  })
  output$compare_fst_plot <- renderPlot({
    dfs <- lapply(input$compare_scenarios, function(s) run_sim(s) %>% mutate(scenario = s))
    df_all <- bind_rows(dfs)
    ggplot(df_all, aes(time/(24*365), FST, color = scenario)) + geom_line() +
      labs(x = "Year", y = "FST (dB)", color = "Scenario") + theme_minimal()
  })
  output$compare_summary_table <- renderDT({
    dfs <- lapply(input$compare_scenarios, function(s) {
      d <- run_sim(s)
      tibble(scenario = s, final_ROD = tail(d$ROD_FRAC, 1), final_CONE = tail(d$CONE_FRAC, 1),
             final_MLMT = tail(d$MLMT, 1), final_FST = tail(d$FST, 1))
    })
    bind_rows(dfs)
  })

  output$rgc_plot <- renderPlot({
    df <- run_sim(input$bio_scenario)
    ggplot(df, aes(time/(24*365), RGC_FRAC)) + geom_line(color = "darkgreen") +
      labs(x = "Year", y = "RGC survival fraction") + theme_minimal()
  })
  output$cme_plot <- renderPlot({
    df <- run_sim(input$bio_scenario)
    ggplot(df, aes(time/(24*365), CME_CST)) + geom_line(color = "darkorange") +
      labs(x = "Year", y = "Central retinal thickness (um)") + theme_minimal()
  })

  output$visualcycle_plot <- renderPlot({
    df <- run_sim(input$gt_scenario)
    ggplot(df, aes(time/(24*365), visual_cycle)) + geom_line(color = "steelblue") +
      labs(x = "Year", y = "Visual cycle flux (relative)", title = "RPE65 visual cycle recovery") + theme_minimal()
  })
  output$bypass_plot <- renderPlot({
    df <- run_sim(input$gt_scenario)
    ggplot(df, aes(time/(24*365), op_bypass)) + geom_line(color = "purple") +
      labs(x = "Year", y = "Optogenetic bypass signal (relative)", title = "MCO-010 bypass pathway (proportional to RGC survival rate)") + theme_minimal()
  })
}

shinyApp(ui, server)
