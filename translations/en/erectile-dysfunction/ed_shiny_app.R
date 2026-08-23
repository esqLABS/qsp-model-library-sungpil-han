## =====================================================================
##  ed_shiny_app.R
##  Erectile Dysfunction QSP model — interactive dashboard
##  Erectile Dysfunction QSP Model — Interactive Dashboard
##
##  14 tabs:
##   1  Patient Profile          patient profile and derived baseline
##   2  Pharmacokinetics                 PK of the four oral PDE5 inhibitors
##   3  Erection Dynamics              intracavernosal pressure / rigidity time course
##   4  Signal Transduction Cascade    NO -> cGMP -> PKG -> Ca -> MLC -> relaxation
##   5  Haemodynamics                 inflow, outflow, veno-occlusion, P-V loop
##   6  Dose-Response              dose-response and the saturation of the gain
##   7  Time Window of Action            time window after a single dose (the 36 h claim)
##   8  Clinical Endpoints        IIEF-EF, SEP2/SEP3, EHS, MCID
##   9  Scenario Comparison          the 16 therapeutic scenarios side by side
##  10  Structural Remodelling          the slow arm: SMI, collagen, smooth muscle
##  11  Penile Rehabilitation              REACTT-style rehabilitation with washout
##  12  Safety                 MAP, nitrate interaction, PDE6/PDE11 occupancy
##  13  Virtual Population              virtual population and responder fraction
##  14  Biomarkers             ADMA, ROS, endocrine axis, Doppler proxies
##
##  Run:  shiny::runApp("ed_shiny_app.R")
##  Requires: shiny, mrgsolve, dplyr, tidyr, ggplot2, DT
## =====================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)

source("ed_mrgsolve_model.R")     # provides ed_mod, ED_DRUGS, ED_ARCH, helpers

THEME <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "#ECEFF1"),
        legend.position = "bottom")

PAL <- c("#1E88E5", "#D81B60", "#43A047", "#FB8C00", "#8E24AA", "#00897B",
         "#6D4C41", "#546E7A")

## ============================== UI ==================================
ui <- fluidPage(
  titlePanel("Erectile Dysfunction QSP Model — 45 ODE / 4 PDE5 Inhibitors"),
  tags$p(style = "color:#546E7A;margin-top:-8px",
         paste("An erection is a threshold readout sitting on a saturating amplifier:",
               "nitric-oxide pulse × cGMP gain × veno-occlusive ceiling → rigidity threshold.",
               "PDE5 inhibitors multiply only the gain.")),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectInput("arch", "Patient archetype",
                  choices = names(ED_ARCH), selected = "vasculo"),
      selectInput("drug", "Oral PDE5 inhibitor",
                  choices = names(ED_DRUGS), selected = "sildenafil"),
      sliderInput("dose", "Dose (mg)", 0, 200, 100, step = 5),
      sliderInput("lead", "Dose → attempt interval (h)", 0.25, 48, 1, step = 0.25),
      hr(),
      h5("Comorbidities / Risk Factors"),
      sliderInput("agey", "Age (years)", 20, 85, 60),
      sliderInput("hba1c", "HbA1c (%)", 4.5, 13, 5.6, step = 0.1),
      sliderInput("bmi", "BMI (kg/m²)", 18, 45, 29, step = 0.5),
      sliderInput("ldl", "LDL-C (mg/dL)", 50, 250, 145, step = 5),
      sliderInput("smoke", "Smoking (0 = non-smoker, 1 = current)", 0, 1, 0.5, step = 0.1),
      sliderInput("raas", "RAAS activity (-)", 1, 1.8, 1.25, step = 0.05),
      hr(),
      h5("Neural · Endocrine · Psychological"),
      sliderInput("nrv", "Cavernosal nerve integrity NRV (0–1)", 0.05, 1, 0.43, step = 0.01),
      sliderInput("nrvmax", "Nerve-sparing recovery ceiling NRVMAX", 0.10, 1, 1.0, step = 0.05),
      sliderInput("kstl", "Leydig cell function KSTL", 40, 260, 229.6, step = 5),
      sliderInput("stress", "Psychological stress STRESS", 0, 1.5, 0, step = 0.05),
      hr(),
      h5("Adjunct · Combination"),
      sliderInput("pge", "Alprostadil ICI (µg)", 0, 40, 0, step = 2),
      sliderInput("pap", "Papaverine (mg)", 0, 40, 0, step = 2),
      sliderInput("phen", "Phentolamine (mg)", 0, 2, 0, step = 0.1),
      sliderInput("gtn", "Nitroglycerin (ng/mL equivalent)", 0, 3, 0, step = 0.05),
      sliderInput("cact", "sGC activator (cinaciguat equivalent)", 0, 3, 0, step = 0.1),
      sliderInput("cstim", "sGC stimulator (riociguat equivalent)", 0, 3, 0, step = 0.1),
      sliderInput("rocki", "ROCK inhibitor", 0, 3, 0, step = 0.1),
      checkboxInput("exer", "Exercise programme", FALSE),
      checkboxInput("statin", "Statin", FALSE),
      actionButton("go", "Run simulation", class = "btn-primary")
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs", type = "tabs",
        tabPanel("1 Patient Profile",
                 h4("Derived Baseline State (all values computed from equations)"),
                 tableOutput("t_base"),
                 plotOutput("p_base", height = "330px"),
                 helpText(paste("SMI = smooth muscle/(smooth muscle+collagen); it sets the veno-occlusive",
                                "ceiling. STEN is the stenosis fraction computed from",
                                "the accumulated oxidative burden."))),
        tabPanel("2 Pharmacokinetics",
                 plotOutput("p_pk", height = "330px"),
                 plotOutput("p_pk2", height = "300px"),
                 tableOutput("t_pk"),
                 helpText("Cu/IC50 determines residual PDE5 activity, 1/(1+Cu/IC50).")),
        tabPanel("3 Erection Dynamics",
                 plotOutput("p_icp", height = "380px"),
                 fluidRow(column(6, tableOutput("t_att")),
                          column(6, plotOutput("p_ehs", height = "260px")))),
        tabPanel("4 Signal Transduction",
                 plotOutput("p_sig", height = "460px"),
                 helpText(paste("NO is generated only while the nerve is firing.",
                                "PDE5 inhibitors do not make NO; they only multiply the cGMP",
                                "gain."))),
        tabPanel("5 Haemodynamics",
                 fluidRow(column(7, plotOutput("p_flow", height = "340px")),
                          column(5, plotOutput("p_pv", height = "340px"))),
                 plotOutput("p_vocc", height = "280px")),
        tabPanel("6 Dose-Response",
                 plotOutput("p_dr", height = "360px"),
                 tableOutput("t_dr"),
                 helpText(paste("Because the amplifier saturates, there is almost no",
                                "change above 100 mg — the structural reason for the labelled maximum dose."))),
        tabPanel("7 Time Window of Action",
                 plotOutput("p_win", height = "380px"),
                 tableOutput("t_win")),
        tabPanel("8 Clinical Endpoints",
                 fluidRow(column(6, plotOutput("p_end", height = "330px")),
                          column(6, plotOutput("p_mcid", height = "330px"))),
                 tableOutput("t_end")),
        tabPanel("9 Scenario Comparison",
                 plotOutput("p_scn", height = "420px"),
                 tableOutput("t_scn")),
        tabPanel("10 Structural Remodelling",
                 sliderInput("months", "Observation period (months)", 1, 36, 12),
                 plotOutput("p_struct", height = "420px"),
                 tableOutput("t_struct")),
        tabPanel("11 Penile Rehabilitation",
                 helpText(paste("REACTT design: 9 months of tadalafil 5 mg daily /",
                                "20 mg on-demand / placebo → 6-week drug-free washout, followed by",
                                "unassisted erection assessment.")),
                 plotOutput("p_rehab", height = "400px"),
                 tableOutput("t_rehab")),
        tabPanel("12 Safety",
                 plotOutput("p_map", height = "330px"),
                 fluidRow(column(6, plotOutput("p_iso", height = "300px")),
                          column(6, tableOutput("t_safe")))),
        tabPanel("13 Virtual Population",
                 sliderInput("npop", "Number of virtual patients", 11, 81, 31, step = 10),
                 sliderInput("sigma", "Log-normal SD of NO-generating capacity", 0.2, 1.0, 0.55,
                             step = 0.05),
                 plotOutput("p_pop", height = "380px"),
                 tableOutput("t_pop")),
        tabPanel("14 Biomarkers",
                 plotOutput("p_bio", height = "420px"),
                 tableOutput("t_bio"))
      )
    )
  )
)

## ============================ server ================================
server <- function(input, output, session) {

  ## ---- assemble the patient ----------------------------------------
  patient <- eventReactive(input$go, {
    A <- ED_ARCH[[input$arch]]
    par <- modifyList(A$par, list(
      AGEY = input$agey, HBA1C = input$hba1c, BMI = input$bmi,
      LDL = input$ldl, SMOKE = input$smoke, RAAS = input$raas,
      KSTL = input$kstl, STRESS = input$stress, NRVMAX = input$nrvmax,
      EXER = as.numeric(input$exer), STATIN = as.numeric(input$statin),
      CACT = input$cact, CSTIM = input$cstim, CROCKI = input$rocki))
    arch2 <- A; arch2$par <- par; arch2$nrv <- input$nrv
    ED_ARCH[["_ui"]] <<- arch2
    ed_baseline(ed_mod, "_ui", input$drug)
  }, ignoreNULL = FALSE)

  attempt <- reactive({
    ed_attempt(patient(), dose = input$dose, lead = input$lead,
               pge = input$pge, pap = input$pap, phen = input$phen,
               gtn = input$gtn, window = input$lead + 4)
  })

  ## ---- 1 patient profile -------------------------------------------
  output$t_base <- renderTable({
    y <- patient()$ini
    d <- as.data.frame(mrgsim(patient()$mod, end = 24, delta = 0.02,
                              hmax = 0.02))
    data.frame(
      `Metric` = c("ROS (fold)", "eNOS-bound fraction", "ADMA (µmol/L)", "Oxidised sGC fraction",
               "Stenosis STEN", "Smooth muscle SM", "Collagen COL", "SMI",
               "Veno-occlusive ceiling", "Corporal pO₂ (mmHg)", "nNOS", "PDE5A expression",
               "Total testosterone (ng/dL)", "LH (IU/L)", "Performance anxiety PA",
               "Peak nocturnal erection ICP (mmHg)"),
      `Value` = sprintf("%.3g", c(y$ROS, y$ECPL, y$ADMA, y$SGCOX, y$STEN, y$SM,
                             y$COL, y$SM / (y$SM + y$COL) * 2, max(d$VOCMAXo),
                             y$PO2, y$NNOS, y$PDE5E, y$TT, y$LH,
                             tail(d$PA, 1), max(d$ICPo))))
  })
  output$p_base <- renderPlot({
    y <- patient()$ini
    ref <- c(ROS = 1, ECPL = 0.714, ADMA = 0.45, SM = 1, COL = 1,
             NNOS = 1, PDE5E = 1, PO2 = 95)
    v <- c(ROS = y$ROS, ECPL = y$ECPL, ADMA = y$ADMA, SM = y$SM, COL = y$COL,
           NNOS = y$NNOS, PDE5E = y$PDE5E, PO2 = y$PO2)
    ggplot(data.frame(k = names(v), rel = as.numeric(v / ref[names(v)])),
           aes(reorder(k, rel), rel, fill = rel > 1)) +
      geom_col(width = .65) + geom_hline(yintercept = 1, linetype = 2) +
      coord_flip() + scale_fill_manual(values = c("#1E88E5", "#D81B60")) +
      labs(x = NULL, y = "Ratio to healthy reference value", fill = "> Normal") + THEME
  })

  ## ---- 2 PK ---------------------------------------------------------
  pkall <- reactive({
    do.call(rbind, lapply(names(ED_DRUGS), function(dg) {
      ds <- c(sildenafil = 100, tadalafil = 20, vardenafil = 20,
              avanafil = 200)[[dg]]
      m <- ed_drug(ed_mod, dg) %>% param(NPTON = 0, ATTT = -1)
      d <- as.data.frame(mrgsim(m, events = ev(amt = ds, cmt = "AGUT"),
                                end = 60, delta = 0.05))
      data.frame(drug = dg, dose = ds, time = d$time, Cp = d$CP1 * 1000,
                 Cu = d$CUo, ratio = d$RATIO, res = d$P5RESo,
                 p6 = d$P6INH, p11 = d$P11INH)
    }))
  })
  output$p_pk <- renderPlot({
    ggplot(pkall(), aes(time, Cu, colour = drug)) + geom_line(linewidth = 1) +
      scale_y_log10() + scale_colour_manual(values = PAL) +
      labs(x = "Time (h)", y = "Unbound plasma concentration Cu (nM)",
           title = "Unbound Exposure at the Labelled Maximum Dose") + THEME
  })
  output$p_pk2 <- renderPlot({
    ggplot(pkall(), aes(time, res, colour = drug)) + geom_line(linewidth = 1) +
      scale_colour_manual(values = PAL) + ylim(0, 1) +
      labs(x = "Time (h)", y = "Residual PDE5 activity 1/(1+Cu/IC50)",
           title = "All Four Drugs Start in the Saturated Region") + THEME
  })
  output$t_pk <- renderTable({
    pkall() %>% group_by(drug, dose) %>%
      summarise(`Cmax (ng/mL)` = max(Cp), `Cu,max (nM)` = max(Cu),
                `Cu/IC50` = max(ratio), `Minimum residual activity` = min(res),
                `PDE6 occupancy (%)` = max(p6), `PDE11 occupancy (%)` = max(p11),
                .groups = "drop")
  }, digits = 3)

  ## ---- 3 erection dynamics -----------------------------------------
  output$p_icp <- renderPlot({
    d <- attempt()$sim
    d %>% select(time, ICPo, RIGo, VSIN) %>%
      pivot_longer(-time) %>%
      mutate(name = recode(name, ICPo = "Intracavernosal pressure ICP (mmHg)",
                           RIGo = "Axial rigidity (%)",
                           VSIN = "Corporal blood volume (mL)")) %>%
      ggplot(aes(time, value)) + geom_line(colour = "#1E88E5", linewidth = 1) +
      geom_hline(data = data.frame(name = "Axial rigidity (%)", y = 60),
                 aes(yintercept = y), linetype = 2, colour = "#D81B60") +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      labs(x = "Time (h)", y = NULL) + THEME
  })
  output$t_att <- renderTable({
    r <- attempt()
    data.frame(`Metric` = c("Peak ICP (mmHg)", "Peak rigidity (%)",
                        "Time above threshold (min)", "Peak cGMP", "Peak NO (nM)",
                        "Peak relaxation R", "Cu (nM)", "Cu/IC50",
                        "SEP2 probability", "SEP3 probability", "IIEF-EF",
                        "Lowest MAP (mmHg)"),
               `Value` = sprintf("%.3g", c(r$ICPpk, r$RIGpk, r$TAEmin, r$CGMPpk,
                                      r$NOpk, r$Rpk, r$CUpk, r$RATIO,
                                      r$PS2, r$PS3, r$IIEF, r$MAPmin)))
  })
  output$p_ehs <- renderPlot({
    d <- attempt()$sim
    ggplot(d, aes(time, EHS)) + geom_step(colour = "#43A047", linewidth = 1) +
      ylim(1, 4) + labs(x = "Time (h)", y = "Erection Hardness Score EHS") + THEME
  })

  ## ---- 4 signalling -------------------------------------------------
  output$p_sig <- renderPlot({
    attempt()$sim %>%
      select(time, NO, CGMP, CAMP, CAI, MLCP, Ro) %>%
      pivot_longer(-time) %>%
      mutate(name = factor(name, c("NO", "CGMP", "CAMP", "CAI", "MLCP", "Ro"),
                           c("NO (nM)", "cGMP (fold)", "cAMP (fold)",
                             "Intracellular Ca²⁺ (fold)", "Phosphorylated MLC20 fraction",
                             "Relaxation R"))) %>%
      ggplot(aes(time, value)) + geom_line(colour = "#8E24AA", linewidth = 1) +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "Time (h)", y = NULL) + THEME
  })

  ## ---- 5 haemodynamics ---------------------------------------------
  output$p_flow <- renderPlot({
    attempt()$sim %>% select(time, QINo, QOUTo) %>% pivot_longer(-time) %>%
      mutate(name = recode(name, QINo = "Arterial inflow Q_in",
                           QOUTo = "Venous outflow Q_out")) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 1) +
      scale_colour_manual(values = PAL) +
      labs(x = "Time (h)", y = "Flow (mL/h)", colour = NULL) + THEME
  })
  output$p_pv <- renderPlot({
    ggplot(attempt()$sim, aes(VSIN, ICPo, colour = time)) +
      geom_path(linewidth = 1) + scale_colour_viridis_c(option = "C") +
      labs(x = "Corporal blood volume (mL)", y = "ICP (mmHg)", colour = "Time (h)",
           title = "Pressure-Volume Loop (Exponential Stiffening of the Tunica)") + THEME
  })
  output$p_vocc <- renderPlot({
    attempt()$sim %>% select(time, VOCCo, GVENo) %>% pivot_longer(-time) %>%
      mutate(name = recode(name, VOCCo = "Veno-occlusive engagement VOCC",
                           GVENo = "Venous leak conductance G_ven")) %>%
      ggplot(aes(time, value)) + geom_line(colour = "#00897B", linewidth = 1) +
      facet_wrap(~name, scales = "free_y") + labs(x = "Time (h)", y = NULL) +
      THEME
  })

  ## ---- 6 dose-response ---------------------------------------------
  drtab <- reactive({
    b <- patient()
    do.call(rbind, lapply(c(0, 5, 10, 25, 50, 100, 150, 200), function(ds) {
      r <- ed_attempt(b, dose = ds, lead = input$lead)
      data.frame(dose = ds, Cu = r$CUpk, ratio = r$RATIO, cGMP = r$CGMPpk,
                 R = r$Rpk, ICP = r$ICPpk, RIG = r$RIGpk, IIEF = r$IIEF)
    }))
  })
  output$p_dr <- renderPlot({
    drtab() %>% select(dose, cGMP, R, ICP, IIEF) %>% pivot_longer(-dose) %>%
      ggplot(aes(dose, value)) +
      geom_line(colour = "#FB8C00", linewidth = 1) + geom_point() +
      facet_wrap(~name, scales = "free_y") +
      labs(x = sprintf("%s dose (mg)", input$drug), y = NULL) + THEME
  })
  output$t_dr <- renderTable(drtab(), digits = 3)

  ## ---- 7 duration window -------------------------------------------
  wintab <- reactive({
    b <- patient()
    lags <- c(0.5, 1, 2, 4, 6, 8, 12, 18, 24, 36, 48, 72)
    do.call(rbind, lapply(lags, function(L) {
      r <- ed_attempt(b, dose = input$dose, lead = L, window = L + 3)
      data.frame(lag = L, Cu = r$CUpk, ICP = r$ICPpk, RIG = r$RIGpk,
                 SEP3 = r$PS3, IIEF = r$IIEF)
    }))
  })
  output$p_win <- renderPlot({
    ggplot(wintab(), aes(lag, RIG)) +
      geom_line(colour = "#1E88E5", linewidth = 1) + geom_point() +
      geom_hline(yintercept = 60, linetype = 2, colour = "#D81B60") +
      scale_x_log10() +
      labs(x = "Time since dose (h, log)", y = "Peak rigidity (%)",
           title = "Effective Window After a Single Dose") + THEME
  })
  output$t_win <- renderTable(wintab(), digits = 3)

  ## ---- 8 endpoints --------------------------------------------------
  output$p_end <- renderPlot({
    b <- patient()
    r0 <- ed_attempt(b); r1 <- attempt()
    df <- data.frame(
      arm = rep(c("No treatment", "Treatment"), each = 4),
      k = rep(c("SEP2", "SEP3", "IIEF-EF/30", "Rigidity/100"), 2),
      v = c(r0$PS2, r0$PS3, r0$IIEF / 30, r0$RIGpk / 100,
            r1$PS2, r1$PS3, r1$IIEF / 30, r1$RIGpk / 100))
    ggplot(df, aes(k, v, fill = arm)) +
      geom_col(position = "dodge", width = .7) +
      scale_fill_manual(values = PAL) + ylim(0, 1) +
      labs(x = NULL, y = "Proportion / normalised score", fill = NULL) + THEME
  })
  output$p_mcid <- renderPlot({
    b <- patient(); r0 <- ed_attempt(b); r1 <- attempt()
    d <- r1$IIEF - r0$IIEF
    sev <- if (r0$IIEF < 11) "Severe" else if (r0$IIEF < 17) "Moderate" else "Mild"
    mcid <- c("Severe" = 7, "Moderate" = 5, "Mild" = 2)[[sev]]
    ggplot(data.frame(x = c("ΔIIEF-EF", "MCID"), y = c(d, mcid)),
           aes(x, y, fill = x)) + geom_col(width = .55) +
      scale_fill_manual(values = PAL) +
      labs(x = NULL, y = "IIEF-EF score",
           title = sprintf("Baseline severity: %s (MCID %d points)", sev, mcid)) + THEME
  })
  output$t_end <- renderTable({
    b <- patient(); r0 <- ed_attempt(b); r1 <- attempt()
    data.frame(`Category` = c("No treatment", "Treatment", "Difference"),
               `IIEF-EF` = c(r0$IIEF, r1$IIEF, r1$IIEF - r0$IIEF),
               SEP2 = c(r0$PS2, r1$PS2, r1$PS2 - r0$PS2),
               SEP3 = c(r0$PS3, r1$PS3, r1$PS3 - r0$PS3),
               `Peak ICP` = c(r0$ICPpk, r1$ICPpk, r1$ICPpk - r0$ICPpk),
               check.names = FALSE)
  }, digits = 3)

  ## ---- 9 scenarios --------------------------------------------------
  scn <- reactive({
    b <- patient()
    rows <- list(
      c("No treatment", 0, 0, 0, 0, 0),
      c("PDE5i low dose", input$dose / 2, 0, 0, 0, 0),
      c("PDE5i standard dose", input$dose, 0, 0, 0, 0),
      c("Alprostadil 10 µg ICI", 0, 10, 0, 0, 0),
      c("Alprostadil 20 µg ICI", 0, 20, 0, 0, 0),
      c("Trimix (triple mixture)", 0, 10, 15, 0.5, 0),
      c("PDE5i + nitroglycerin", input$dose, 0, 0, 0, 0.76))
    do.call(rbind, lapply(rows, function(z) {
      r <- ed_attempt(b, dose = as.numeric(z[2]), lead = input$lead,
                      pge = as.numeric(z[3]), pap = as.numeric(z[4]),
                      phen = as.numeric(z[5]), gtn = as.numeric(z[6]))
      data.frame(`Scenario` = z[1], ICP = r$ICPpk, `Rigidity` = r$RIGpk,
                 SEP3 = r$PS3, `IIEF-EF` = r$IIEF, `Lowest MAP` = r$MAPmin,
                 check.names = FALSE)
    }))
  })
  output$p_scn <- renderPlot({
    scn() %>% mutate(`Scenario` = factor(`Scenario`, `Scenario`)) %>%
      ggplot(aes(`Scenario`, `IIEF-EF`, fill = ICP)) +
      geom_col(width = .7) + coord_flip() +
      scale_fill_viridis_c(option = "D") +
      labs(x = NULL, y = "IIEF-EF", fill = "Peak ICP") + THEME
  })
  output$t_scn <- renderTable(scn(), digits = 3)

  ## ---- 10 structural arm -------------------------------------------
  strun <- reactive({
    m <- patient()$mod %>% param(ATTEVERY = 120)
    as.data.frame(mrgsim(m, end = 24 * 30 * input$months, delta = 12,
                         hmax = 0.05))
  })
  output$p_struct <- renderPlot({
    strun() %>% mutate(day = time / 24) %>%
      select(day, SM, COL, SMIo, TGFB, PO2, LEN, NRV, ERFR) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value)) + geom_line(colour = "#6D4C41", linewidth = 1) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "Day", y = NULL, title = "The Slow Arm: Week-to-Month Time Constants") + THEME
  })
  output$t_struct <- renderTable({
    d <- strun()
    data.frame(`Metric` = c("SMI start", "SMI end", "Collagen end", "Smooth muscle end",
                        "Penile length (cm)", "NRV end", "Corporal pO₂"),
               `Value` = sprintf("%.3g", c(d$SMIo[1], tail(d$SMIo, 1),
                                      tail(d$COL, 1), tail(d$SM, 1),
                                      tail(d$LEN, 1), tail(d$NRV, 1),
                                      tail(d$PO2, 1))))
  })

  ## ---- 11 rehabilitation -------------------------------------------
  rehab <- reactive({
    tab <- ed_s9_rehab(ed_mod)
    tab
  })
  output$p_rehab <- renderPlot({
    rehab() %>% select(arm, SMI_onTx, SMI_wash, ICP_unassisted,
                       IIEF_unassisted) %>%
      pivot_longer(-arm) %>%
      ggplot(aes(arm, value, fill = arm)) + geom_col(width = .6) +
      facet_wrap(~name, scales = "free_y") +
      scale_fill_manual(values = PAL) +
      labs(x = NULL, y = NULL, fill = NULL,
           title = "9-Month Treatment Followed by 6-Week Washout") + THEME
  })
  output$t_rehab <- renderTable(rehab(), digits = 3)

  ## ---- 12 safety ----------------------------------------------------
  output$p_map <- renderPlot({
    b <- patient()
    arms <- list(`No treatment` = list(0, 0), `Nitro only` = list(0, 0.76),
                 `PDE5i only` = list(input$dose, 0),
                 `PDE5i + nitro` = list(input$dose, 0.76))
    d <- do.call(rbind, lapply(names(arms), function(a) {
      r <- ed_attempt(b, dose = arms[[a]][[1]], lead = input$lead,
                      gtn = arms[[a]][[2]], window = input$lead + 3)
      data.frame(arm = a, time = r$sim$time, MAP = r$sim$MAP)
    }))
    ggplot(d, aes(time, MAP, colour = arm)) + geom_line(linewidth = 1) +
      scale_colour_manual(values = PAL) +
      labs(x = "Time (h)", y = "Mean arterial pressure (mmHg)", colour = NULL,
           title = "The Same Product of Two Factors Read Out in the Systemic Vasculature") + THEME
  })
  output$p_iso <- renderPlot({
    d <- attempt()$sim
    d %>% select(time, P6INH, P11INH) %>% pivot_longer(-time) %>%
      mutate(name = recode(name, P6INH = "PDE6 occupancy (%) — vision",
                           P11INH = "PDE11 occupancy (%) — myalgia")) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 1) +
      scale_colour_manual(values = PAL) +
      labs(x = "Time (h)", y = "Inhibition (%)", colour = NULL) + THEME
  })
  output$t_safe <- renderTable({
    r <- attempt()
    data.frame(`Metric` = c("Lowest MAP (mmHg)", "MAP change (mmHg)",
                        "Maximum PDE6 occupancy (%)", "Maximum PDE11 occupancy (%)",
                        "Maximum systemic cGMP (fold)"),
               `Value` = sprintf("%.3g", c(r$MAPmin, min(r$sim$DMAP),
                                      max(r$sim$P6INH), max(r$sim$P11INH),
                                      max(r$sim$SCG))))
  })

  ## ---- 13 virtual population ---------------------------------------
  poptab <- reactive({
    q <- (seq_len(input$npop) - 0.5) / input$npop
    fn <- exp(input$sigma * qnorm(q))
    do.call(rbind, lapply(fn, function(f) {
      b <- ed_baseline(ed_mod, "_ui", input$drug, fnoi = f)
      r0 <- ed_attempt(b); r1 <- ed_attempt(b, dose = input$dose,
                                            lead = input$lead)
      data.frame(FNOI = f, IIEF0 = r0$IIEF, IIEF1 = r1$IIEF,
                 ICP0 = r0$ICPpk, ICP1 = r1$ICPpk,
                 SEP30 = r0$PS3, SEP31 = r1$PS3)
    }))
  })
  output$p_pop <- renderPlot({
    poptab() %>% select(FNOI, IIEF0, IIEF1) %>% pivot_longer(-FNOI) %>%
      mutate(name = recode(name, IIEF0 = "No treatment", IIEF1 = "Treatment")) %>%
      ggplot(aes(FNOI, value, colour = name)) +
      geom_line(linewidth = 1) + geom_point(size = 1) + scale_x_log10() +
      scale_colour_manual(values = PAL) +
      labs(x = "Individual NO-generating-capacity multiple FNOI (log)", y = "IIEF-EF",
           colour = NULL,
           title = "Threshold Readout + Inter-Individual Variation = Population Average") + THEME
  })
  output$t_pop <- renderTable({
    d <- poptab()
    data.frame(`Category` = c("No treatment", "Treatment"),
               `Mean IIEF-EF` = c(mean(d$IIEF0), mean(d$IIEF1)),
               `Mean SEP3` = c(mean(d$SEP30), mean(d$SEP31)),
               `Responder fraction (SEP3≥0.5)` = c(mean(d$SEP30 >= .5),
                                            mean(d$SEP31 >= .5)),
               check.names = FALSE)
  }, digits = 3)

  ## ---- 14 biomarkers ------------------------------------------------
  output$p_bio <- renderPlot({
    strun() %>% mutate(day = time / 24) %>%
      select(day, ROS, ADMA, ECPL, OXD, SGCOX, TT, LH, HCT, PSV) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value)) + geom_line(colour = "#546E7A", linewidth = 1) +
      facet_wrap(~name, scales = "free_y") + labs(x = "Day", y = NULL) + THEME
  })
  output$t_bio <- renderTable({
    d <- strun()
    data.frame(`Biomarker` = c("ROS (fold)", "ADMA (µmol/L)", "eNOS binding",
                              "Oxidative damage", "Oxidised sGC", "Testosterone",
                              "LH", "Haematocrit (%)",
                              "Cavernosal artery peak velocity surrogate (cm/s)"),
               `Value` = sprintf("%.3g", c(tail(d$ROS, 1), tail(d$ADMA, 1),
                                      tail(d$ECPL, 1), tail(d$OXD, 1),
                                      tail(d$SGCOX, 1), tail(d$TT, 1),
                                      tail(d$LH, 1), tail(d$HCT, 1),
                                      max(d$PSV))))
  })
}

shinyApp(ui, server)
