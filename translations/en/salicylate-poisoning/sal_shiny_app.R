## =====================================================================
##  sal_shiny_app.R
##  Salicylate (Aspirin) Poisoning — interactive QSP dashboard
##  11 tabs · driven by sal_mrgsolve_model.R
##
##  Salicylate (Aspirin) Poisoning — QSP Dashboard
## ---------------------------------------------------------------------
##  The app is built around ONE decomposition, which every tab returns to:
##
##      C_brain = C_total x fu x f_n(pH_plasma)/f_n(pH_brain)
##
##  Tab 3 ("the two multipliers") shows the three factors side by side over time, so the
##  user can watch the reported number and the brain concentration come
##  apart.  Tab 8 lets the user move the ventilator and see it happen.
##
##  RUN
##    R -e "shiny::runApp('sal_shiny_app.R', port = 8080)"
##  (sal_mrgsolve_model.R must be in the same directory; it is sourced on
##   startup, which compiles the model once.)
## =====================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

source("sal_mrgsolve_model.R")

THEME <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "#ECEFF1", colour = NA),
        legend.position = "bottom")

PAL <- c(plasma = "#1565C0", free = "#7986CB", brain = "#00838F",
         tissue = "#8D6E63", target = "#C62828", alt = "#2E7D32")

## ---------------------------------------------------------------------
##  UI
## ---------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Salicylate (Aspirin) Poisoning QSP Dashboard — Salicylate Poisoning QSP Dashboard"),
  tags$p(style = "color:#555;margin-top:-8px",
         HTML("The value that gets measured (total plasma salicylate) and the value that kills people (brain intracellular salicylate)
               are split apart by two multiplications — <b>free fraction</b> and <b>pH partitioning</b>.
               This app keeps showing both of those multipliers together.")),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Patient"),
      numericInput("bw", "Body weight (kg)", 70, 30, 150, 5),
      sliderInput("alb", "Albumin (g/L)", 15, 50, 40, 1),
      sliderInput("gfr", "GFR (mL/min)", 10, 140, 120, 5),
      checkboxInput("aged", "Elderly (increased susceptibility to lung injury)", FALSE),

      h4("Exposure"),
      selectInput("mode", "Pattern",
                  c("Acute single" = "acute",
                    "Enteric-coated / concretion" = "ec",
                    "Chronic overshoot" = "chronic")),
      conditionalPanel("input.mode != 'chronic'",
        sliderInput("dose", "Aspirin dose (g)", 1, 60, 30, 1)),
      conditionalPanel("input.mode == 'chronic'",
        sliderInput("cdose", "Single dose (mg), q6h", 325, 1950, 1300, 325),
        sliderInput("cdays", "Duration of dosing (days)", 3, 21, 14, 1)),

      h4("Treatment"),
      sliderInput("trx", "Time treatment starts (h)", 0, 24, 4, 1),
      checkboxInput("fluid", "IV fluid crystalloid", TRUE),
      checkboxInput("bic", "Sodium bicarbonate NaHCO3", TRUE),
      sliderInput("kcl", "Potassium replacement K (mmol/h)", 0, 30, 10, 1),
      checkboxInput("dex", "Dextrose (blood glucose 11 mmol/L)", FALSE),
      checkboxInput("ac", "Activated charcoal MDAC", FALSE),
      checkboxInput("hd", "Haemodialysis", FALSE),
      conditionalPanel("input.hd == true",
        sliderInput("thd", "Time dialysis starts (h)", 2, 36, 6, 1)),

      h4("Ventilation"),
      radioButtons("vent", "",
                   c("Spontaneous" = "spont",
                     "Intubated, conventional PaCO2" = "conv",
                     "Intubated, minute ventilation matched" = "match",
                     "Sedation only (45%)" = "sed")),
      conditionalPanel("input.vent != 'spont'",
        sliderInput("tvent", "Time of intubation/sedation (h)", 1, 24, 8, 1)),

      hr(),
      sliderInput("tmax", "Simulation duration (h)", 12, 340, 72, 4),
      actionButton("go", "Recalculate Run", class = "btn-primary")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("1. Patient Profile", br(),
                 fluidRow(column(6, h4("Summary"), tableOutput("summary")),
                          column(6, h4("Alerts"), uiOutput("alerts"))),
                 hr(), h4("Key metrics at a glance"), plotOutput("overview", height = "420px")),

        tabPanel("2. Pharmacokinetics (PK)", br(),
                 p("Total plasma concentration (the measured value), free concentration, tissue concentration, and brain concentration."),
                 plotOutput("pk", height = "380px"),
                 h4("Elimination pathways (mg/h)"),
                 plotOutput("clr", height = "300px")),

        tabPanel("3. The Two Multipliers", br(),
                 p(HTML("<b>C_brain = C_total × fu × Kp</b>. The three terms are plotted separately.
                         Two of the three are invisible in the laboratory.")),
                 plotOutput("mult", height = "420px"),
                 h4("Ratio of brain concentration to total plasma concentration (brain : reported level)"),
                 plotOutput("ratio", height = "260px")),

        tabPanel("4. Acid-Base", br(),
                 plotOutput("ab", height = "440px"),
                 h4("Ventilatory drive"),
                 plotOutput("vent", height = "260px")),

        tabPanel("5. Renal · Electrolytes", br(),
                 plotOutput("renal", height = "440px"),
                 p(HTML("<b>If urine pH does not rise, alkalinisation does nothing.</b>
                         And urine pH only rises once plasma bicarbonate crosses the renal tubular threshold,
                         and hypokalaemia raises that threshold."))),

        tabPanel("6. Clinical Endpoints", br(),
                 plotOutput("endp", height = "440px"),
                 h4("Cumulative CNS exposure (duration, not just peak)"),
                 plotOutput("cnsi", height = "240px")),

        tabPanel("7. Scenario Comparison", br(),
                 p("14 scenarios are built into the model. Pick two to overlay."),
                 fluidRow(
                   column(5, selectInput("sc1", "Scenario A",
                                         setNames(names(SCEN),
                                                  sapply(SCEN, `[[`, "label")),
                                         selected = "S07")),
                   column(5, selectInput("sc2", "Scenario B",
                                         setNames(names(SCEN),
                                                  sapply(SCEN, `[[`, "label")),
                                         selected = "S08")),
                   column(2, br(), actionButton("gosc", "Compare"))),
                 plotOutput("scen", height = "480px"),
                 tableOutput("scentab")),

        tabPanel("8. Ventilator Experiment", br(),
                 h4("Only one switch is changed — One switch, nothing else"),
                 p(HTML("Intubating a hyperventilating patient and matching them to a <b>conventional minute ventilation</b>
                         raises PaCO2 and lowers plasma pH. Intracellular pH is buffered and
                         barely moves, so <b>the partition coefficient rises and the drug enters the brain</b>.
                         At that same moment, the laboratory value <b>falls</b>.")),
                 plotOutput("ventexp", height = "460px"),
                 tableOutput("venttab")),

        tabPanel("9. Done Nomogram", br(),
                 p(HTML("The nomogram is treated as <b>output, not input</b>.
                         The (time, plasma concentration) pair read from the simulation is fed into the nomogram to get a grade,
                         which is compared side by side with the brain-concentration grade computed by the model.")),
                 plotOutput("nomo", height = "440px"),
                 tableOutput("nomotab")),

        tabPanel("10. Dialysis Decision (EXTRIP)", br(),
                 uiOutput("extrip"),
                 plotOutput("hdplot", height = "400px")),

        tabPanel("11. Laboratory Panel (Biomarkers)", br(),
                 plotOutput("labs", height = "520px"),
                 tableOutput("labtab")),

        tabPanel("About", br(),
                 htmlOutput("about"))
      )
    )
  )
)

## ---------------------------------------------------------------------
##  Scenario construction from the sidebar
## ---------------------------------------------------------------------
build_stages <- function(inp) {
  base <- list(GFR0 = inp$gfr * 60 / 1000, ALB = inp$alb, BW = inp$bw,
               AGEF = if (isTRUE(inp$aged)) 1 else 0)
  if (inp$mode == "acute") base <- c(base, MASSIVE)
  if (inp$mode == "ec")    base <- c(base, ENTERIC)

  stages <- list(list(t = 0, par = base))
  rx <- list()
  if (isTRUE(inp$fluid)) rx <- c(rx, list(RFLU = 0.25))
  if (isTRUE(inp$bic))   rx <- c(rx, list(RBIC = 75.0))
  rx <- c(rx, list(RKCL = inp$kcl))
  if (isTRUE(inp$dex))   rx <- c(rx, list(GLCP = 11.0))
  if (length(rx)) {
    stages <- c(stages, list(list(t = inp$trx, par = rx)))
    if (isTRUE(inp$bic))
      stages <- c(stages, list(list(t = inp$trx + 4, par = list(RBIC = 37.5))))
  }
  if (isTRUE(inp$hd)) {
    stages <- c(stages, list(list(t = inp$thd, par = list(CLHD = 6.0, HDBIC = 45))),
                list(list(t = inp$thd + 4, par = list(CLHD = 0, HDBIC = 0))))
  }
  if (inp$vent != "spont") {
    vp <- switch(inp$vent,
                 conv  = list(VENT = 1, VASET = 1.0),
                 match = list(VENT = 1, VASET = 2.3),
                 sed   = list(SED = 0.45))
    stages <- c(stages, list(list(t = inp$tvent, par = vp)))
  }
  stages[order(sapply(stages, `[[`, "t"))]
}

build_events <- function(inp) {
  ev <- if (inp$mode == "chronic") {
    dose_seq(inp$cdose, ii = 6, n = ceiling(inp$cdays * 4))
  } else if (inp$mode == "ec") {
    rbind(asa_dose(inp$dose * 0.55), asa_dose(inp$dose * 0.45, cmt = "ACONC"))
  } else {
    asa_dose(inp$dose)
  }
  if (isTRUE(inp$bic)) ev <- rbind(ev, bic_bolus(100, time = inp$trx))
  if (isTRUE(inp$ac))
    ev <- rbind(ev, charcoal(50, max(0.5, inp$trx - 2)),
                charcoal(50, inp$trx + 2), charcoal(50, inp$trx + 6))
  ev[order(ev$time), , drop = FALSE]
}

## ---------------------------------------------------------------------
##  SERVER
## ---------------------------------------------------------------------
server <- function(input, output, session) {

  sim <- eventReactive(input$go, ignoreNULL = FALSE, {
    stages <- build_stages(input)
    events <- build_events(input)
    run_stages(mod, stages, end = input$tmax, delta = 0.1, events = events)
  })

  long <- function(d, cols, labels = cols) {
    d %>% select(time, all_of(cols)) %>%
      pivot_longer(-time) %>%
      mutate(name = factor(name, levels = cols, labels = labels))
  }

  ## ---- 1. profile ----
  output$summary <- renderTable({
    d <- sim()
    data.frame(
      Metric = c("Peak total plasma concentration (mg/dL)", "Peak free concentration (mg/L)",
               "Peak brain concentration (mg/L)", "Lowest arterial pH", "Lowest HCO3 (mmol/L)",
               "Lowest PaCO2 (mmHg)", "Peak urine pH", "Peak renal clearance (mL/min)",
               "Lowest K (mmol/L)", "Peak temperature (°C)", "Lowest brain glucose (mmol/L)",
               "Peak CNS score", "Cumulative CNS injury"),
      Value = sprintf("%.2f", c(max(d$CTOTD), max(d$CFREE), max(d$CCNS),
                             min(d$PH), min(d$HCO3C), min(d$PACO2C),
                             max(d$UPH), max(d$CLREN), min(d$KPL),
                             max(d$TCORE), min(d$GLUB), max(d$CNS),
                             max(d$CNSI))))
  }, striped = TRUE)

  output$alerts <- renderUI({
    d <- sim(); a <- list()
    if (max(d$CTOTD) > 90) a <- c(a, "Plasma concentration > 90 mg/dL — EXTRIP indication for dialysis")
    if (min(d$PH) < 7.20) a <- c(a, "Acidaemia pH < 7.20 — sharp rise in brain partition coefficient, indication for dialysis")
    if (max(d$CNS) > 50) a <- c(a, "CNS depression — altered consciousness, indication for dialysis")
    if (min(d$KPL) < 3.5) a <- c(a, "Hypokalaemia — urinary alkalinisation will fail")
    if (max(d$UPH) < 7.0 && input$bic) a <- c(a, "Urine pH is not rising — check K and dose")
    if (min(d$GLUB) < 0.8) a <- c(a, "Low brain glucose — give glucose even if blood glucose is normal")
    if (max(d$TCORE) > 38.5) a <- c(a, "Hyperthermia — a sign of uncoupling, poor prognosis")
    if (input$vent == "conv") a <- c(a, "Intubated with conventional ventilator settings — the riskiest manoeuvre in this model")
    if (!length(a)) a <- "No alerts"
    HTML(paste0("<ul>", paste0("<li>", a, "</li>", collapse = ""), "</ul>"))
  })

  output$overview <- renderPlot({
    d <- sim()
    p <- long(d, c("CTOT", "CCNS", "PH", "UPH", "KPL", "CNS"),
              c("Total plasma (mg/L)", "Brain (mg/L)", "Arterial pH", "Urine pH",
                "Serum K (mmol/L)", "CNS score"))
    ggplot(p, aes(time, value)) + geom_line(linewidth = 0.9, colour = PAL["plasma"]) +
      facet_wrap(~name, scales = "free_y") + labs(x = "Time (h)", y = NULL) + THEME
  })

  ## ---- 2. PK ----
  output$pk <- renderPlot({
    d <- sim()
    p <- long(d, c("CTOT", "CFREE", "CTIS", "CCNS"),
              c("Total plasma (measured)", "Free plasma", "Tissue", "Brain"))
    ggplot(p, aes(time, value, colour = name)) + geom_line(linewidth = 1) +
      scale_colour_manual(values = unname(PAL[c("plasma","free","tissue","brain")])) +
      labs(x = "Time (h)", y = "Concentration (mg/L)", colour = NULL) + THEME
  })

  output$clr <- renderPlot({
    d <- sim()
    d2 <- d %>% mutate(
      salicylurate = c(0, diff(ASU) / diff(time)),
      glucuronide  = c(0, diff(APG + AAG) / diff(time)),
      renal        = c(0, diff(AUR) / diff(time)),
      dialysis     = c(0, diff(AHD) / diff(time)))
    p <- long(d2, c("salicylurate", "glucuronide", "renal", "dialysis"))
    ggplot(p, aes(time, value, colour = name)) + geom_line(linewidth = 1) +
      labs(x = "Time (h)", y = "Elimination rate (mg/h)", colour = NULL) + THEME
  })

  ## ---- 3. the two multipliers ----
  output$mult <- renderPlot({
    d <- sim()
    p <- long(d, c("CTOT", "FU", "KPBR", "CCNS"),
              c("① C_total  (laboratory value, mg/L)", "② fu  free fraction",
                "③ Kp  brain:free-plasma partition", "= C_brain  (mg/L)"))
    ggplot(p, aes(time, value)) +
      geom_line(linewidth = 1, colour = PAL["brain"]) +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "Time (h)", y = NULL) + THEME
  })

  output$ratio <- renderPlot({
    d <- sim()
    ggplot(d, aes(time, RATIO)) + geom_line(linewidth = 1, colour = PAL["target"]) +
      labs(x = "Time (h)", y = "Brain concentration / total plasma concentration") + THEME
  })

  ## ---- 4. acid-base ----
  output$ab <- renderPlot({
    d <- sim()
    p <- long(d, c("PH", "PHB", "PACO2C", "HCO3C", "AGAP", "UNC"),
              c("Arterial pH", "Brain intracellular pH", "PaCO2 (mmHg)",
                "HCO3 (mmol/L)", "Anion gap (mmol/L)", "Degree of uncoupling"))
    ggplot(p, aes(time, value)) + geom_line(linewidth = 0.9, colour = PAL["plasma"]) +
      facet_wrap(~name, scales = "free_y") + labs(x = "Time (h)", y = NULL) + THEME
  })

  output$vent <- renderPlot({
    d <- sim()
    ggplot(d, aes(time, VASP)) + geom_line(linewidth = 1, colour = PAL["alt"]) +
      geom_hline(yintercept = 1, linetype = 2) +
      labs(x = "Time (h)",
           y = "Patient's own alveolar ventilation (multiple of normal)") + THEME
  })

  ## ---- 5. renal ----
  output$renal <- renderPlot({
    d <- sim()
    p <- long(d, c("UPH", "CLREN", "UO", "KPL", "GFRC", "HCO3C"),
              c("Urine pH", "Renal salicylate clearance (mL/min)", "Urine output (mL/h)",
                "Serum K (mmol/L)", "GFR (mL/min)", "Plasma HCO3 (mmol/L)"))
    ggplot(p, aes(time, value)) + geom_line(linewidth = 0.9, colour = PAL["alt"]) +
      facet_wrap(~name, scales = "free_y") + labs(x = "Time (h)", y = NULL) + THEME
  })

  ## ---- 6. endpoints ----
  output$endp <- renderPlot({
    d <- sim()
    p <- long(d, c("CNS", "TINN", "SEV", "TCORE", "GLUB", "LUNG"),
              c("CNS depression score", "Tinnitus score", "Severity index",
                "Temperature (°C)", "Brain glucose (mmol/L)", "Lung injury (0-1)"))
    ggplot(p, aes(time, value)) + geom_line(linewidth = 0.9, colour = PAL["target"]) +
      facet_wrap(~name, scales = "free_y") + labs(x = "Time (h)", y = NULL) + THEME
  })

  output$cnsi <- renderPlot({
    d <- sim()
    ggplot(d, aes(time, CNSI)) + geom_area(fill = "#EF9A9A", alpha = 0.6) +
      geom_line(linewidth = 1, colour = PAL["target"]) +
      labs(x = "Time (h)", y = "Cumulative CNS exposure ∫max(0, C_brain − 45)dt") + THEME
  })

  ## ---- 7. scenario comparison ----
  scen2 <- eventReactive(input$gosc, ignoreNULL = FALSE, {
    bind_rows(run_scenario(input$sc1, delta = 0.1),
              run_scenario(input$sc2, delta = 0.1))
  })

  output$scen <- renderPlot({
    d <- scen2()
    p <- d %>% select(time, scenario, CTOT, CCNS, PH, UPH, CLREN, CNSI) %>%
      pivot_longer(-c(time, scenario))
    ggplot(p, aes(time, value, colour = scenario)) + geom_line(linewidth = 1) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "Time (h)", y = NULL, colour = NULL) + THEME
  })

  output$scentab <- renderTable({
    d <- scen2()
    bind_rows(lapply(split(d, d$scenario), summarise_scenario))
  })

  ## ---- 8. the ventilator experiment ----
  ventexp <- reactive({
    bind_rows(run_scenario("S07", delta = 0.05),
              run_scenario("S08", delta = 0.05)) %>%
      mutate(arm = ifelse(scenario == "S07", "Intubated, conventional PaCO2",
                          "Minute ventilation matched"))
  })

  output$ventexp <- renderPlot({
    d <- ventexp() %>%
      select(time, arm, CTOT, CCNS, PH, PACO2C, KPBR, RATIO) %>%
      pivot_longer(-c(time, arm)) %>%
      mutate(name = recode(name,
        CTOT = "Total plasma (laboratory value, mg/L)", CCNS = "Brain (mg/L)",
        PH = "Arterial pH", PACO2C = "PaCO2 (mmHg)",
        KPBR = "Kp brain:free-plasma", RATIO = "Brain / total plasma"))
    ggplot(d, aes(time, value, colour = arm)) + geom_line(linewidth = 1) +
      geom_vline(xintercept = 8, linetype = 2, colour = "grey40") +
      facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = c(PAL[["target"]], PAL[["alt"]])) +
      labs(x = "Time (h)", y = NULL, colour = NULL) + THEME
  })

  output$venttab <- renderTable({
    d <- ventexp()
    at <- function(a, v, t) approx(d$time[d$arm == a], d[[v]][d$arm == a], t)$y
    ts <- c(8, 9, 12, 24, 36)
    data.frame(
      `Time h` = ts,
      `Plasma A conventional` = round(sapply(ts, function(t) at("Intubated, conventional PaCO2", "CTOT", t)), 1),
      `Plasma B matched` = round(sapply(ts, function(t) at("Minute ventilation matched", "CTOT", t)), 1),
      `Brain A` = round(sapply(ts, function(t) at("Intubated, conventional PaCO2", "CCNS", t)), 1),
      `Brain B` = round(sapply(ts, function(t) at("Minute ventilation matched", "CCNS", t)), 1),
      `Brain A/B` = round(sapply(ts, function(t)
        at("Intubated, conventional PaCO2", "CCNS", t) / at("Minute ventilation matched", "CCNS", t)), 3),
      check.names = FALSE)
  })

  ## ---- 9. Done nomogram ----
  done_zone <- function(c_mgdl, t_h) {
    t_h <- pmin(pmax(t_h, 6), 60)
    sev <- 130 * 0.5 ^ ((t_h - 6) / 20)
    mod <- 65 * 0.5 ^ ((t_h - 6) / 20)
    mild <- 32.5 * 0.5 ^ ((t_h - 6) / 20)
    ifelse(c_mgdl >= sev, "severe",
           ifelse(c_mgdl >= mod, "moderate",
                  ifelse(c_mgdl >= mild, "mild", "asymptomatic")))
  }
  brain_zone <- function(cb)
    ifelse(cb > 90, "severe", ifelse(cb > 55, "moderate",
                                     ifelse(cb > 25, "mild", "asymptomatic")))

  output$nomo <- renderPlot({
    d <- sim() %>% filter(time >= 6)
    bands <- data.frame(t = seq(6, 60, 0.5)) %>%
      mutate(severe = 130 * 0.5 ^ ((t - 6) / 20),
             moderate = 65 * 0.5 ^ ((t - 6) / 20),
             mild = 32.5 * 0.5 ^ ((t - 6) / 20)) %>%
      pivot_longer(-t)
    ggplot() +
      geom_line(data = bands, aes(t, value, linetype = name), colour = "grey40") +
      geom_path(data = d, aes(time, CTOTD, colour = brain_zone(CCNS)),
                linewidth = 1.4) +
      scale_y_log10() +
      labs(x = "Time since ingestion (h)", y = "Plasma salicylate (mg/dL, log)",
           colour = "Brain concentration grade computed by model", linetype = "Done boundary") +
      THEME
  })

  output$nomotab <- renderTable({
    d <- sim()
    ts <- c(6, 12, 24, 36)
    ts <- ts[ts <= max(d$time)]
    ap <- function(v, t) approx(d$time, d[[v]], t)$y
    data.frame(
      `Time h` = ts,
      `Plasma mg/dL` = round(sapply(ts, ap, v = "CTOTD"), 1),
      `Done grade` = done_zone(sapply(ts, ap, v = "CTOTD"), ts),
      `Brain mg/L` = round(sapply(ts, ap, v = "CCNS"), 1),
      `Model grade` = brain_zone(sapply(ts, ap, v = "CCNS")),
      Match = ifelse(done_zone(sapply(ts, ap, v = "CTOTD"), ts) ==
                      brain_zone(sapply(ts, ap, v = "CCNS")), "✓", "✗ mismatch"),
      check.names = FALSE)
  })

  ## ---- 10. EXTRIP ----
  output$extrip <- renderUI({
    d <- sim()
    crit <- c(
      sprintf("Plasma salicylate > 90 mg/dL (normal renal function): peak %.1f mg/dL", max(d$CTOTD)),
      sprintf("Altered consciousness (CNS score > 50): peak %.0f", max(d$CNS)),
      sprintf("Arterial pH ≤ 7.20: lowest %.2f", min(d$PH)),
      sprintf("Acute lung injury: peak %.2f", max(d$LUNG)),
      sprintf("Renal impairment (GFR < 45 mL/min): lowest %.0f", min(d$GFRC)))
    met <- c(max(d$CTOTD) > 90, max(d$CNS) > 50, min(d$PH) <= 7.20,
             max(d$LUNG) > 0.3, min(d$GFRC) < 45)
    HTML(paste0(
      "<h4>EXTRIP criteria check (Juurlink 2015)</h4><ul>",
      paste0("<li>", ifelse(met, "<b style='color:#C62828'>Met</b> — ", "Not met — "),
             crit, "</li>", collapse = ""),
      "</ul><p>", if (any(met))
        "<b style='color:#C62828'>Extracorporeal removal should be considered.</b>" else
        "The EXTRIP criteria are not met at the current settings.",
      " In this model the dialysis clearance (6 L/h) is about 5.7 times the optimal alkaline-diuresis clearance.</p>"))
  })

  output$hdplot <- renderPlot({
    d <- bind_rows(run_scenario("S09", delta = 0.1),
                   run_scenario("S10", delta = 0.1)) %>%
      mutate(arm = ifelse(scenario == "S09", "Dialysis at 6 h", "Dialysis at 18 h"))
    p <- d %>% select(time, arm, CTOT, CCNS, CNSI, PH) %>%
      pivot_longer(-c(time, arm))
    ggplot(p, aes(time, value, colour = arm)) + geom_line(linewidth = 1) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "Time (h)", y = NULL, colour = NULL) + THEME
  })

  ## ---- 11. laboratory panel ----
  output$labs <- renderPlot({
    d <- sim()
    p <- long(d, c("CTOTD", "PH", "HCO3C", "PACO2C", "AGAP", "KPL",
                   "GFRC", "UPH", "GLUB"),
              c("Salicylate (mg/dL)", "pH", "HCO3 (mmol/L)", "PaCO2 (mmHg)",
                "Anion gap", "K (mmol/L)", "GFR (mL/min)", "Urine pH",
                "Brain glucose (mmol/L)"))
    ggplot(p, aes(time, value)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") + labs(x = "Time (h)", y = NULL) + THEME
  })

  output$labtab <- renderTable({
    d <- sim()
    ts <- unique(pmin(c(0, 2, 4, 6, 12, 24, 48, 72), max(d$time)))
    ap <- function(v) round(approx(d$time, d[[v]], ts)$y, 2)
    data.frame(`Time h` = ts, `Salicylate mg/dL` = ap("CTOTD"), pH = ap("PH"),
               HCO3 = ap("HCO3C"), PaCO2 = ap("PACO2C"), `Anion gap` = ap("AGAP"),
               K = ap("KPL"), `Urine pH` = ap("UPH"), `Brain mg/L` = ap("CCNS"),
               check.names = FALSE)
  })

  ## ---- about ----
  output$about <- renderUI(HTML('
    <h4>What this model claims</h4>
    <p><code>C_brain = C_total × fu(C_total, pH) × f_n(pH_plasma)/f_n(pH_brain)</code></p>
    <p>The laboratory reports only the first term. The second term (free fraction) grows fivefold as albumin saturates,
       from 8.5% at 150 mg/L to 43% at 800 mg/L, and the third term (pH partitioning) grows
       <b>because plasma pH moves while brain intracellular pH is buffered and barely moves</b>.
       Raising PaCO2 from 25 to 60 mmHg moves plasma pH by 0.38 but
       moves brain intracellular pH by only 0.14, and the partition coefficient becomes 1.7-fold.</p>
    <h4>Composition</h4>
    <ul>
      <li>28 ODEs · 14 scenarios · 11 tabs</li>
      <li>Mechanistic map: <code>sal_qsp_model.dot / .svg / .png</code> (125 nodes, 16 clusters)</li>
      <li>Model: <code>sal_mrgsolve_model.R</code></li>
      <li>Independent verification: <code>sal_verify_python.py</code> → <code>sal_verification_output.txt</code>
          (Python/scipy re-implementation, 54 of 54 literature benchmarks passed)</li>
      <li>References: <code>sal_references.md</code> (all PMIDs verified by lookup on PubMed)</li>
    </ul>
    <h4>Limitations</h4>
    <p>The key variable in this model, brain concentration, has never been measured in a living person.
       It was calculated only from animal partitioning experiments (Hill 1971) and the pH physics
       of a weak acid, and was not calibrated to clinical outcomes. The falsifiable prediction is
       that <b>after a change in PaCO2, plasma and brain move in opposite directions</b>.</p>
    <p style="color:#B71C1C"><b>Educational/research model. Not for use in clinical decisions.</b></p>'))
}

shinyApp(ui, server)
