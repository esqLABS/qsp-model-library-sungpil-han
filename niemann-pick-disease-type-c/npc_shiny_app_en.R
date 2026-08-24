## =============================================================================
##  Niemann-Pick disease type C (NPC) — QSP Shiny dashboard
##  Niemann-Pick disease type C · interactive dashboard
##
##  9 tabs:
##    1. Patient profile          genotype x developmental vulnerability
##    2. Drug PK                  for all four agents
##    3. Lysosome · lipids        lysosomal biology (the mechanism)
##    4. Cerebellum · reserve     the damage-reserve integral
##    5. Clinical endpoints       SARA, NPCCSS 17/5/4, saccades
##    6. Scenario comparison      18 regimens
##    7. Biomarkers               and the compartment claim
##    8. Trial-design experiment  the design claim
##    9. Calibration table        model vs published targets, including the misses
##
##  Run:  shiny::runApp("npc_shiny_app.R")
##  Requires: shiny, mrgsolve, dplyr, tidyr, ggplot2, DT
## =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

MODEL_FILE <- "npc_mrgsolve_model_en.R"
YR <- 365.25

mod_global <- mread_cache("npc", MODEL_FILE)

## ---------------------------------------------------------------------------
## genotype and archetype tables (kept in sync with the .R model's helpers)
## ---------------------------------------------------------------------------
npc_genotypes <- list(
  `WT (normal)`                = list(f_null = 0.00, theta0 = 1.000, f_npc2 = 1.00),
  `I1061T/I1061T (juvenile)`   = list(f_null = 0.00, theta0 = 0.055, f_npc2 = 1.00),
  `I1061T/null (late-infantile)` = list(f_null = 0.50, theta0 = 0.055, f_npc2 = 1.00),
  `null/null (perinatal)`      = list(f_null = 1.00, theta0 = 0.055, f_npc2 = 1.00),
  `mild/mild (adult)`          = list(f_null = 0.00, theta0 = 0.322, f_npc2 = 1.00),
  `NPC2 disease`               = list(f_null = 0.00, theta0 = 1.000, f_npc2 = 0.03)
)

DRUGS <- c("miglustat", "arimoclomol", "levacetylleucine", "IT adrabetadex")

make_events <- function(drugs, start_yr, dur_yr,
                        mig_mg = 200, ari_mg = 124, nal_g = 4,
                        cd_mg = 900, cd_weeks = 2) {
  t0 <- start_yr * YR
  dur <- dur_yr * YR
  es <- list()
  if ("miglustat" %in% drugs)
    es <- c(es, list(ev(time = t0, amt = mig_mg, cmt = "MIG_GUT",
                        ii = 1/3, addl = ceiling(dur * 3) - 1)))
  if ("arimoclomol" %in% drugs)
    es <- c(es, list(ev(time = t0, amt = ari_mg, cmt = "ARI_GUT",
                        ii = 1/3, addl = ceiling(dur * 3) - 1)))
  if ("levacetylleucine" %in% drugs)
    es <- c(es, list(ev(time = t0, amt = nal_g * 1000 / 3, cmt = "NAL_GUT",
                        ii = 1/3, addl = ceiling(dur * 3) - 1)))
  if ("IT adrabetadex" %in% drugs)
    es <- c(es, list(ev(time = t0, amt = cd_mg, cmt = "CD_CSF",
                        ii = 7 * cd_weeks,
                        addl = max(0, floor(dur / (7 * cd_weeks)) - 1))))
  if (!length(es)) return(ev(time = 0, amt = 0, cmt = "MIG_GUT"))
  Reduce(`+`, es)
}

run_sim <- function(geno, v_dev, liver_risk, years, events, delta = 7, ...) {
  g <- npc_genotypes[[geno]]
  m <- param(mod_global, f_null = g$f_null, theta0 = g$theta0,
             f_npc2 = g$f_npc2, v_dev = v_dev, liver_risk = liver_risk, ...)
  m %>% mrgsim(events = events, end = years * YR, delta = delta,
               hmax = 0.5, atol = 1e-8, rtol = 1e-8) %>% as_tibble()
}

at_age <- function(out, col, age) approx(out$AGE_YR, out[[col]], xout = age, rule = 2)$y

theme_npc <- function() {
  theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          legend.position = "bottom",
          plot.title = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(size = 10, colour = "grey30"))
}

longify <- function(out, cols) {
  out %>% select(AGE_YR, all_of(cols)) %>%
    pivot_longer(-AGE_YR, names_to = "variable", values_to = "value")
}

## =============================================================================
##  UI
## =============================================================================
ui <- fluidPage(
  titlePanel("Niemann-Pick disease type C (NPC) — QSP dashboard"),
  tags$p(style = "color:#555; margin-top:-8px;",
         HTML("The plasma markers read <b>visceral</b> storage while the disease is in the <b>cerebellum</b> · ",
              "function = 1 − D<sub>rev</sub> − D<sub>irr</sub> · ",
              "each drug was tested in a design that could only see its own mechanism")),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("① Patient"),
      selectInput("geno", "Genotype",
                  choices = names(npc_genotypes),
                  selected = "I1061T/I1061T (juvenile)"),
      sliderInput("v_dev", "Developmental vulnerability v_dev", 0, 3, 0.6, step = 0.1),
      helpText(HTML("<small>Residual NPC1 activity alone does not fix the age-of-onset form. ",
                    "Because the storage phenotype saturates, siblings with the same genotype ",
                    "can differ by ten years.</small>")),
      checkboxInput("liver", "Perinatal liver disease risk (liver_risk)", FALSE),
      numericInput("years", "Follow-up horizon (years)", 30, min = 2, max = 70),
      hr(),
      h4("② Therapy"),
      checkboxGroupInput("drugs", "Drugs", choices = DRUGS, selected = character(0)),
      sliderInput("start", "Age at start of therapy (years)", 0, 30, 13, step = 0.5),
      sliderInput("dur", "Duration of therapy (years)", 0.25, 30, 3, step = 0.25),
      hr(),
      h4("③ Dose"),
      numericInput("mig_mg", "Miglustat mg/dose (tid)", 200, min = 0, max = 400, step = 25),
      numericInput("ari_mg", "Arimoclomol mg/dose (tid)", 124, min = 0, max = 372, step = 31),
      numericInput("nal_g",  "Levacetylleucine g/day", 4, min = 0, max = 6, step = 0.5),
      numericInput("cd_mg",  "Adrabetadex mg (IT)", 900, min = 0, max = 1800, step = 100),
      numericInput("cd_wk",  "IT dosing interval (weeks)", 2, min = 1, max = 8, step = 1),
      hr(),
      h4("④ Mechanism switches"),
      checkboxInput("nal_sym", "Levacetylleucine symptomatic component ON", TRUE),
      checkboxInput("nal_dm",  "Levacetylleucine disease-modifying component ON", TRUE),
      helpText(HTML("<small>Turning the two switches off separately lets you reproduce ",
                    "the design experiment of tab 8 by hand.</small>"))
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("1. Patient profile",
                 br(), fluidRow(column(6, plotOutput("p_profile", height = 300)),
                                column(6, plotOutput("p_func", height = 300))),
                 hr(), h5("Summary measures"), tableOutput("t_profile"),
                 helpText(HTML("f_NPC1 = lysosomal NPC1 function relative to normal. ",
                               "NPC1 and NPC2 act <b>in series</b> on a single flux, so their ",
                               "functional fractions multiply."))),
        tabPanel("2. Drug PK",
                 br(), plotOutput("p_pk", height = 420),
                 hr(), plotOutput("p_pd_target", height = 300),
                 helpText(HTML("Miglustat's UGCG inhibition is <b>about twice as strong in the viscera as in the CNS</b> ",
                               "(brain:plasma partition 0.45). This is a prediction, not an assumption, and it is ",
                               "why miglustat's clearest effects are the visceral and bulbar ones."))),
        tabPanel("3. Lysosome · lipids",
                 br(), plotOutput("p_lipid", height = 380),
                 hr(), plotOutput("p_lyso", height = 320)),
        tabPanel("4. Cerebellum · reserve",
                 br(), plotOutput("p_dam", height = 320),
                 hr(), plotOutput("p_pc", height = 340),
                 helpText(HTML("<b>DAM is the integral of stress and has no repair term.</b> ",
                               "The delay to onset is <i>inversely proportional</i> to stress, while the ",
                               "death rate after the gate opens <i>saturates</i> — which is why the ",
                               "progression slope is almost independent of the age of onset (PMID 19415691)."))),
        tabPanel("5. Clinical endpoints",
                 br(), plotOutput("p_scales", height = 380),
                 hr(), fluidRow(column(6, plotOutput("p_sacc", height = 300)),
                                column(6, plotOutput("p_surv", height = 300)))),
        tabPanel("6. Scenario comparison",
                 br(),
                 helpText("Every scenario in the table below is applied to the patient set in the left panel."),
                 actionButton("run_scen", "Run the 18 scenarios", class = "btn-primary"),
                 br(), br(), DTOutput("t_scen"),
                 hr(), plotOutput("p_scen", height = 420)),
        tabPanel("7. Biomarkers",
                 br(), plotOutput("p_bio", height = 380),
                 hr(), plotOutput("p_compart", height = 320),
                 helpText(HTML("<b>The compartment claim.</b> C-triol production saturates. In patients it ",
                               "already sits at about 82% of its own ceiling, so it cannot reflect any ",
                               "further storage. That is why the published triol–NPCCSS5 correlation is ",
                               "only ρ = 0.265 (PMID 33228797) — not because of noise, but because it is ",
                               "measuring a <i>different compartment</i>."))),
        tabPanel("8. Trial-design experiment",
                 br(),
                 helpText(HTML("A <b>purely symptomatic drug</b> and a <b>purely disease-modifying drug</b>, ",
                               "matched to the same 12-month efficacy, are each run through two published designs.")),
                 numericInput("dsg_weeks", "Design A: observation window (weeks)", 12, min = 4, max = 104),
                 actionButton("run_dsg", "Run the design experiment", class = "btn-primary"),
                 br(), br(), tableOutput("t_design"),
                 hr(), plotOutput("p_design", height = 340),
                 helpText(HTML("A 12-week window sees only D<sub>rev</sub>. A 12-month change-from-baseline ",
                               "cannot distinguish a constant offset from a change in slope. ",
                               "The only design that separates the two components is <b>randomised withdrawal</b>."))),
        tabPanel("9. Calibration table",
                 br(), h4("Model vs published values"), DTOutput("t_calib"),
                 hr(),
                 h4("Items that do not match (reported without adjustment)"),
                 htmlOutput("misses"))
      )
    )
  )
)

## =============================================================================
##  SERVER
## =============================================================================
server <- function(input, output, session) {

  ev_now <- reactive({
    make_events(input$drugs, input$start, input$dur,
                input$mig_mg, input$ari_mg, input$nal_g, input$cd_mg, input$cd_wk)
  })

  sim <- reactive({
    run_sim(input$geno, input$v_dev, as.numeric(input$liver),
            input$years, ev_now(), delta = 7,
            nal_sym_on = as.numeric(input$nal_sym),
            nal_dm_on  = as.numeric(input$nal_dm))
  })

  sim_untr <- reactive({
    run_sim(input$geno, input$v_dev, as.numeric(input$liver),
            input$years, ev(time = 0, amt = 0, cmt = "MIG_GUT"), delta = 7)
  })

  ## ---- tab 1 --------------------------------------------------------------
  output$p_profile <- renderPlot({
    o <- sim()
    longify(o, c("f_NPC1_out", "f_eg", "theta")) %>%
      ggplot(aes(AGE_YR, value, colour = variable)) +
      geom_line(linewidth = 0.9) +
      labs(title = "NPC1 folding yield and egress capacity",
           subtitle = "theta = ER folding yield · f_eg = NPC1 x NPC2 (in series)",
           x = "Age (years)", y = "Fraction of normal", colour = NULL) +
      theme_npc()
  })

  output$p_func <- renderPlot({
    o <- sim()
    longify(o, c("D_REV", "D_IRR", "FUNC")) %>%
      ggplot(aes(AGE_YR, value, colour = variable)) +
      geom_line(linewidth = 0.9) +
      scale_colour_manual(values = c(D_REV = "#e67e22", D_IRR = "#c0392b",
                                     FUNC = "#1a7f37")) +
      labs(title = "function = 1 − D_rev − D_irr",
           subtitle = "D_rev comes back and D_irr does not",
           x = "Age (years)", y = NULL, colour = NULL) +
      theme_npc()
  })

  output$t_profile <- renderTable({
    o <- sim(); u <- sim_untr()
    ages <- c(5, 10, 15, 20)
    ages <- ages[ages <= input$years]
    data.frame(
      `Age` = ages,
      `f_NPC1` = round(sapply(ages, function(a) at_age(o, "f_NPC1_out", a)), 4),
      `CHOL_C (x normal)` = round(sapply(ages, function(a) at_age(o, "CHOL_C_fold", a)), 1),
      `DAM/reserve` = round(sapply(ages, function(a) at_age(o, "DAM", a)) / 5707.9, 2),
      `NPCCSS5` = round(sapply(ages, function(a) at_age(o, "NPCCSS5", a)), 2),
      `SARA` = round(sapply(ages, function(a) at_age(o, "SARA", a)), 2),
      `Untreated NPCCSS5` = round(sapply(ages, function(a) at_age(u, "NPCCSS5", a)), 2),
      check.names = FALSE)
  })

  ## ---- tab 2 --------------------------------------------------------------
  output$p_pk <- renderPlot({
    o <- sim() %>% filter(AGE_YR >= input$start, AGE_YR <= input$start + 0.12)
    if (!nrow(o)) o <- sim()
    longify(o, c("C_MIG", "C_ARI", "C_NAL", "CD_CSF_CONC")) %>%
      ggplot(aes(AGE_YR, value)) +
      geom_line(linewidth = 0.7, colour = "#8e44ad") +
      facet_wrap(~variable, scales = "free_y", ncol = 2) +
      labs(title = "Drug concentrations (the first 6 weeks after start of therapy)",
           subtitle = "C_MIG uM · C_ARI ng/mL · C_NAL mg/L · CD_CSF mg/L",
           x = "Age (years)", y = NULL) +
      theme_npc()
  })

  output$p_pd_target <- renderPlot({
    o <- sim()
    longify(o, c("I_mig_v", "I_mig_c", "HSP70", "E_nal_sym", "E_nal_dm")) %>%
      ggplot(aes(AGE_YR, value, colour = variable)) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~variable, scales = "free_y", ncol = 3) +
      labs(title = "Pharmacodynamics at the level of the target",
           subtitle = "I_mig_v = visceral UGCG inhibition · I_mig_c = CNS UGCG inhibition",
           x = "Age (years)", y = NULL) +
      theme_npc() + theme(legend.position = "none")
  })

  ## ---- tab 3 --------------------------------------------------------------
  output$p_lipid <- renderPlot({
    o <- sim()
    longify(o, c("CHOL_V_fold", "CHOL_C_fold", "GSL_V", "GSL_C", "SPH")) %>%
      ggplot(aes(AGE_YR, value, colour = variable)) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~variable, scales = "free_y", ncol = 3) +
      labs(title = "Lipid storage — a time constant of days in the viscera, years in the brain",
           subtitle = "tau_cns = 25 : the same steady state, reached 25 times more slowly",
           x = "Age (years)", y = NULL) +
      theme_npc() + theme(legend.position = "none")
  })

  output$p_lyso <- renderPlot({
    o <- sim()
    longify(o, c("CA_LY", "HYD", "AUTOPH", "MITO", "ROS", "GLUC_METAB")) %>%
      ggplot(aes(AGE_YR, value, colour = variable)) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~variable, scales = "free_y", ncol = 3) +
      labs(title = "Lysosomal and metabolic function",
           subtitle = "sphingosine → blockade of the acidic Ca²⁺ store → reduced fusion and hydrolysis (PMID 18953351)",
           x = "Age (years)", y = NULL) +
      theme_npc() + theme(legend.position = "none")
  })

  ## ---- tab 4 --------------------------------------------------------------
  output$p_dam <- renderPlot({
    o <- sim()
    o %>% mutate(reserve_used = DAM / 5707.9) %>%
      select(AGE_YR, reserve_used, gate, stress) %>%
      pivot_longer(-AGE_YR) %>%
      ggplot(aes(AGE_YR, value, colour = name)) +
      geom_line(linewidth = 0.9) +
      geom_hline(yintercept = 1, linetype = 2, colour = "grey40") +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(title = "The damage integral and the reserve gate",
           subtitle = "past reserve_used = 1 Purkinje death begins, and cannot be undone",
           x = "Age (years)", y = NULL) +
      theme_npc() + theme(legend.position = "none")
  })

  output$p_pc <- renderPlot({
    o <- sim()
    longify(o, c("PC", "PC_S", "PC_LOST", "INFL", "SYN", "CBL")) %>%
      ggplot(aes(AGE_YR, value, colour = variable)) +
      geom_line(linewidth = 0.9) +
      labs(title = "Purkinje cell pools · neuroinflammation · cerebellar volume",
           subtitle = "PC + PC_S + PC_LOST = 1 (a conserved quantity)",
           x = "Age (years)", y = "Fraction / index", colour = NULL) +
      theme_npc()
  })

  ## ---- tab 5 --------------------------------------------------------------
  output$p_scales <- renderPlot({
    o <- sim(); u <- sim_untr()
    bind_rows(mutate(longify(o, c("SARA","NPCCSS5","NPCCSS4","NPCCSS17")), arm = "Treated"),
              mutate(longify(u, c("SARA","NPCCSS5","NPCCSS4","NPCCSS17")), arm = "Untreated")) %>%
      ggplot(aes(AGE_YR, value, colour = arm)) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~variable, scales = "free_y", ncol = 2) +
      scale_colour_manual(values = c(`Treated` = "#8e44ad", `Untreated` = "grey40")) +
      labs(title = "Clinical scales", x = "Age (years)", y = NULL, colour = NULL) +
      theme_npc()
  })

  output$p_sacc <- renderPlot({
    o <- sim(); u <- sim_untr()
    bind_rows(mutate(longify(o, c("SACCADE","SWALLOW","HEARING_dB")), arm = "Treated"),
              mutate(longify(u, c("SACCADE","SWALLOW","HEARING_dB")), arm = "Untreated")) %>%
      ggplot(aes(AGE_YR, value, colour = arm)) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~variable, scales = "free_y", ncol = 1) +
      scale_colour_manual(values = c(`Treated` = "#8e44ad", `Untreated` = "grey40")) +
      labs(title = "Saccade velocity · swallowing · hearing threshold shift",
           subtitle = "the hearing threshold shift is cyclodextrin ototoxicity (the same mechanism as its efficacy)",
           x = "Age (years)", y = NULL, colour = NULL) +
      theme_npc()
  })

  output$p_surv <- renderPlot({
    o <- sim(); u <- sim_untr()
    bind_rows(mutate(select(o, AGE_YR, SURV), arm = "Treated"),
              mutate(select(u, AGE_YR, SURV), arm = "Untreated")) %>%
      ggplot(aes(AGE_YR, SURV, colour = arm)) +
      geom_line(linewidth = 0.9) +
      geom_hline(yintercept = 0.5, linetype = 2, colour = "grey50") +
      scale_colour_manual(values = c(`Treated` = "#8e44ad", `Untreated` = "grey40")) +
      labs(title = "Survival probability", subtitle = "the hazard is proportional to the square of swallowing function (PMID 23039766)",
           x = "Age (years)", y = "S(t)", colour = NULL) +
      theme_npc()
  })

  ## ---- tab 6 --------------------------------------------------------------
  scen <- eventReactive(input$run_scen, {
    defs <- list(
      list(id = "S01", lab = "Untreated",                  d = character(0), s = 13, u = 3),
      list(id = "S06", lab = "Miglustat",                  d = "miglustat",  s = 13, u = 3),
      list(id = "S07", lab = "Levacetylleucine",           d = "levacetylleucine", s = 13, u = 3),
      list(id = "S08", lab = "Arimoclomol",                d = "arimoclomol", s = 13, u = 3),
      list(id = "S09", lab = "Arimoclomol+miglustat",      d = c("arimoclomol","miglustat"), s = 13, u = 3),
      list(id = "S10", lab = "Levacetylleucine+miglustat", d = c("levacetylleucine","miglustat"), s = 13, u = 3),
      list(id = "S11", lab = "Triple therapy",             d = c("levacetylleucine","miglustat","arimoclomol"), s = 13, u = 3),
      list(id = "S12", lab = "IT cyclodextrin",            d = "IT adrabetadex", s = 13, u = 3),
      list(id = "S13", lab = "IT CD + miglustat",          d = c("IT adrabetadex","miglustat"), s = 13, u = 3),
      list(id = "S14", lab = "Triple, from age 2",         d = c("levacetylleucine","miglustat","arimoclomol"), s = 2,  u = 18),
      list(id = "S15", lab = "Triple, from age 8",         d = c("levacetylleucine","miglustat","arimoclomol"), s = 8,  u = 12),
      list(id = "S16", lab = "Triple, from age 12",        d = c("levacetylleucine","miglustat","arimoclomol"), s = 12, u = 8)
    )
    withProgress(message = "Running scenarios...", {
      purrr::map_dfr(seq_along(defs), function(i) {
        s <- defs[[i]]; incProgress(1 / length(defs))
        o <- run_sim(input$geno, input$v_dev, as.numeric(input$liver), 20,
                     make_events(s$d, s$s, s$u, input$mig_mg, input$ari_mg,
                                 input$nal_g, input$cd_mg, input$cd_wk), delta = 14)
        mutate(o, scenario = s$id, label = s$lab)
      })
    })
  })

  output$t_scen <- renderDT({
    scen() %>% group_by(scenario, label) %>%
      summarise(`NPCCSS5@20` = round(at_age(cur_data_all(), "NPCCSS5", 20), 2),
                `SARA@20` = round(at_age(cur_data_all(), "SARA", 20), 2),
                `D_rev@20` = round(at_age(cur_data_all(), "D_REV", 20), 3),
                `D_irr@20` = round(at_age(cur_data_all(), "D_IRR", 20), 3),
                `PC_lost@20` = round(at_age(cur_data_all(), "PC_LOST", 20), 3),
                `Hearing dB@20` = round(at_age(cur_data_all(), "HEARING_dB", 20), 1),
                .groups = "drop") %>%
      datatable(options = list(pageLength = 12, dom = "t"), rownames = FALSE)
  })

  output$p_scen <- renderPlot({
    ggplot(scen(), aes(AGE_YR, NPCCSS5, colour = label)) +
      geom_line(linewidth = 0.8) +
      labs(title = "5-domain NPCCSS by scenario",
           subtitle = "what an early start buys is D_rev; what is already lost is D_irr",
           x = "Age (years)", y = "5-domain NPCCSS", colour = NULL) +
      theme_npc()
  })

  ## ---- tab 8 --------------------------------------------------------------
  design <- eventReactive(input$run_dsg, {
    wk <- input$dsg_weeks
    horizon <- 13 * YR + wk * 7
    run1 <- function(sym, dm, drug_on) {
      e <- if (drug_on) make_events("levacetylleucine", 13, wk * 7 / YR,
                                    nal_g = input$nal_g)
           else ev(time = 0, amt = 0, cmt = "MIG_GUT")
      run_sim(input$geno, input$v_dev, as.numeric(input$liver),
              horizon / YR, e, delta = 3.5,
              nal_sym_on = sym, nal_dm_on = dm)
    }
    pbo <- run1(1, 1, FALSE); sym <- run1(1, 0, TRUE); dmo <- run1(0, 1, TRUE)
    tab <- tibble(
      arm = c("Placebo", "Pure symptomatic", "Pure disease-modifying"),
      dSARA = c(at_age(pbo,"SARA",horizon/YR) - at_age(pbo,"SARA",13),
                at_age(sym,"SARA",horizon/YR) - at_age(sym,"SARA",13),
                at_age(dmo,"SARA",horizon/YR) - at_age(dmo,"SARA",13)),
      dN5   = c(at_age(pbo,"NPCCSS5",horizon/YR) - at_age(pbo,"NPCCSS5",13),
                at_age(sym,"NPCCSS5",horizon/YR) - at_age(sym,"NPCCSS5",13),
                at_age(dmo,"NPCCSS5",horizon/YR) - at_age(dmo,"NPCCSS5",13))
    ) %>% mutate(`SARA vs placebo` = dSARA - dSARA[1],
                 `NPCCSS5 vs placebo` = dN5 - dN5[1])
    list(tab = tab, traj = bind_rows(mutate(pbo, arm = "Placebo"),
                                     mutate(sym, arm = "Pure symptomatic"),
                                     mutate(dmo, arm = "Pure disease-modifying")))
  })

  output$t_design <- renderTable({ design()$tab }, digits = 3)

  output$p_design <- renderPlot({
    design()$traj %>%
      ggplot(aes(AGE_YR, SARA, colour = arm)) +
      geom_line(linewidth = 0.9) +
      labs(title = "The same 12-month efficacy, entirely different short-term SARA trajectories",
           subtitle = "a short window sees only the symptomatic component",
           x = "Age (years)", y = "SARA", colour = NULL) +
      theme_npc()
  })

  ## ---- tab 9 --------------------------------------------------------------
  output$t_calib <- renderDT({
    oj <- run_sim("I1061T/I1061T (juvenile)", 0.6, 0,
                  35, ev(time = 0, amt = 0, cmt = "MIG_GUT"), delta = 7)
    ow <- run_sim("WT (normal)", 0.6, 0, 45,
                  ev(time = 0, amt = 0, cmt = "MIG_GUT"), delta = 30)
    tibble(
      `Target` = c("T1 plasma triol, patient (ng/mL)", "T2 plasma triol, control (ng/mL)",
                 "T3 5-domain NPCCSS slope (points/year)", "T4 17-domain slope (points/year)",
                 "T11 SARA at trial entry"),
      `Role` = c("calibration", "calibration", "calibration", "validation", "calibration"),
      `Published` = c(88.31, 5.97, 1.50, 2.80, 15.91),
      `Model` = round(c(at_age(oj, "TRIOL", 10), at_age(ow, "TRIOL", 10),
                       at_age(oj,"NPCCSS5",13) - at_age(oj,"NPCCSS5",12),
                       at_age(oj,"NPCCSS17",13) - at_age(oj,"NPCCSS17",12),
                       at_age(oj, "SARA", 13)), 3)
    ) %>% mutate(`Relative error %` = round(100 * (`Model` - `Published`) / `Published`, 1)) %>%
      datatable(options = list(dom = "t"), rownames = FALSE)
  })

  output$misses <- renderUI({
    HTML(paste0(
      "<ul>",
      "<li><b>Arimoclomol effect size.</b> The published 12-month 5-domain difference is −1.40 ",
      "(95% CI −2.76 ~ −0.03). Even with the folding yield restored <i>completely</i>, the model ",
      "removes only about 48% of the progression, and stays around −0.7 points in a patient entering at 13. ",
      "The model value is <i>inside</i> the published confidence interval but far smaller than the point estimate. ",
      "The model's reading: the reserve integral cannot be undone, so <b>repairing the primary defect after ",
      "the reserve has already been spent cannot return the rate of progression to normal.</b> ",
      "In patients entering earlier the effect is larger (see tab 4).</li>",
      "<li><b>Levacetylleucine long-term extension.</b> The published change from baseline in SARA is ",
      "−1.88 at 12 months and −1.64 at 18 months. The model reproduces the 12-week difference (−1.28) exactly, ",
      "but with a symptomatic effect of that size the change from baseline should already have crossed zero ",
      "by 12 months. The untreated SARA progression implied by the published values is 0.22 points/year, ",
      "whereas the value the model <i>derives</i> from the NPCCSS targets is 1.5 points/year — ",
      "a discrepancy of about sevenfold. There is no published untreated SARA slope for NPC, so ",
      "this discrepancy is left unresolved.</li>",
      "<li><b>Separating the late-infantile and juvenile forms.</b> Because the storage phenotype saturates, ",
      "a twofold difference in residual NPC1 activity barely changes the steady-state burden. ",
      "The two age forms are therefore specified <i>separately</i> through the developmental vulnerability ",
      "v_dev, and are not derived from the genotype.</li>",
      "</ul>"))
  })
}

shinyApp(ui, server)
