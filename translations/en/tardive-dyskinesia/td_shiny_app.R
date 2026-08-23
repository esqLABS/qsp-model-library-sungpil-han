## =============================================================================
## Tardive Dyskinesia QSP — Shiny dashboard
## Tardive Dyskinesia QSP Model Interactive Dashboard (9 tabs)
##
##   shiny::runApp("td_shiny_app.R")
##
## The model itself lives in td_mrgsolve_model.R; this file only loads it
## (TD_SKIP_RUN prevents the batch analyses from running on source) and wires
## it to a UI.
##
## Tabs
##   1 Patient · Prescription Profile   patient, offending drug, adjuncts
##   2 PK · Receptor Occupancy     concentrations, striatal D2 and VMAT2 occupancy
##   3 Plasticity State          RUP / ROS / SDAM with the latch threshold drawn
##   4 Basal Ganglia Circuit          IND-GPe-STN-GPi-THAL and the EXC decomposition
##   5 AIMS · Clinical Endpoints AIMS, psychosis, parkinsonism, mood, QTc, adherence
##   6 Scenario Comparison        up to four strategies on one axis
##   7 Window of Reversibility            withdrawal-timing scan -> point of no return
##   8 Opposing-Lever Frontier   AIMS vs psychosis / parkinsonism trade-off
##   9 Calibration Criteria · Literature Anchors  calibration anchors table
## =============================================================================

library(shiny)
TD_SKIP_RUN <- TRUE
source("td_mrgsolve_model.R")

PAL <- c(aims = "#c0392b", rup = "#8e44ad", sdam = "#d35400",
         ros = "#16a085", occ = "#2471a3", occv = "#117864",
         psy = "#7f8c8d", park = "#b7950b", depr = "#5b2c6f",
         a = "#2471a3", b = "#c0392b", c = "#117864", d = "#b7950b")

lineplot <- function(x, ys, cols, labs, ylab, main, ylim = NULL, hlines = NULL,
                     vlines = NULL, lwd = 2.2) {
  if (is.null(ylim)) ylim <- range(0, unlist(ys), na.rm = TRUE)
  plot(x, ys[[1]], type = "n", ylim = ylim, xlab = "Day", ylab = ylab,
       main = main, las = 1, bty = "l")
  grid(col = "gray90")
  if (!is.null(hlines))
    for (h in seq_along(hlines))
      abline(h = hlines[[h]]$y, lty = 3, col = hlines[[h]]$col)
  if (!is.null(vlines)) for (v in vlines) abline(v = v, lty = 2, col = "gray55")
  for (i in seq_along(ys)) lines(x, ys[[i]], col = cols[i], lwd = lwd)
  legend("topleft", legend = labs, col = cols, lwd = lwd, bty = "n",
         cex = 0.85)
}

## -----------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Tardive Dyskinesia QSP Model Dashboard"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Patient (Host)"),
      sliderInput("age", "Age (year)", 20, 85, 58, step = 1),
      checkboxInput("fga", "First-generation antipsychotic-level oxidative burden (FGA)", TRUE),
      checkboxInput("dm", "Diabetes · metabolic syndrome", FALSE),
      checkboxInput("estro", "Estrogen replete (neuroprotective)", FALSE),
      sliderInput("gen", "Genetic risk multiplier (GEN_RISK)", 0.8, 1.5, 1.0,
                  step = 0.05),
      hr(),
      h4("Offending drug"),
      sliderInput("dose1", "Initial dose (risperidone-eq mg/day)", 0, 16, 8,
                  step = 0.5),
      sliderInput("dose2", "Dose after change (mg/day)", 0, 16, 8, step = 0.5),
      sliderInput("tsw", "Timing of dose change (day)", 0, 1825, 730, step = 5),
      sliderInput("flai", "Fraction administered as LAI", 0, 1, 0, step = 0.1),
      checkboxInput("esc", "Clinical up-titration policy (PSYCH>0.20 → increase dose)", FALSE),
      hr(),
      h4("TD-directed therapy"),
      selectInput("vmat", "VMAT2 inhibitor",
                  c("None" = "none", "Valbenazine" = "val",
                    "Deutetrabenazine" = "dtb")),
      sliderInput("vdose", "VMAT2 inhibitor dose (mg/day)", 0, 120, 80, step = 4),
      sliderInput("tvmat", "VMAT2 inhibitor start (day)", 0, 1825, 730, step = 5),
      selectInput("cyp", "CYP2D6 phenotype",
                  c("Normal (1.0x)" = "1.0", "UM (1.6x)" = "1.6",
                    "IM (0.7x)" = "0.7", "PM (0.5x)" = "0.5",
                    "PM + CYP3A4 inhibition (0.35x)" = "0.35")),
      sliderInput("clz", "Clozapine switch dose (mg/day, 0=none)", 0, 600, 0,
                  step = 25),
      hr(),
      h4("Adjuncts"),
      checkboxInput("gkb", "Ginkgo biloba EGb761 240 mg (antioxidant)", FALSE),
      sliderInput("ama", "Amantadine (mg/day)", 0, 400, 0, step = 50),
      sliderInput("ach", "Benztropine (mg/day)", 0, 4, 0, step = 0.5),
      checkboxInput("bont", "Botulinum toxin (local)", FALSE),
      hr(),
      sliderInput("days", "Simulation duration (day)", 365, 3650, 1825,
                  step = 65)
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        type = "tabs",
        tabPanel("1 Patient · Prescription", br(), verbatimTextOutput("profile"),
                 plotOutput("p_summary", height = "340px")),
        tabPanel("2 PK · Occupancy", br(), plotOutput("p_pk", height = "300px"),
                 plotOutput("p_occ", height = "300px")),
        tabPanel("3 Plasticity State", br(),
                 plotOutput("p_plast", height = "340px"),
                 verbatimTextOutput("latchtxt")),
        tabPanel("4 Basal Ganglia Circuit", br(), plotOutput("p_bg", height = "300px"),
                 plotOutput("p_exc", height = "300px")),
        tabPanel("5 AIMS · Endpoints", br(),
                 plotOutput("p_aims", height = "300px"),
                 plotOutput("p_ep", height = "300px")),
        tabPanel("6 Scenario Comparison", br(),
                 plotOutput("p_scen", height = "420px"),
                 tableOutput("t_scen")),
        tabPanel("7 Window of Reversibility", br(),
                 sliderInput("wgrid", "Withdrawal-timing grid (day)", 60, 1460,
                             c(60, 730), step = 10, width = "60%"),
                 plotOutput("p_window", height = "380px"),
                 tableOutput("t_window")),
        tabPanel("8 Opposing Lever", br(), plotOutput("p_front", height = "420px"),
                 tableOutput("t_front")),
        tabPanel("9 Calibration Criteria", br(), tableOutput("t_anchor"),
                 verbatimTextOutput("disclaim"))
      )
    )
  )
)

## -----------------------------------------------------------------------------
server <- function(input, output, session) {

  pset <- reactive({
    p <- list(
      AGE = input$age, RISK_FGA = if (input$fga) 1.6 else 1.0,
      DM_RISK = as.numeric(input$dm), ESTROGEN = as.numeric(input$estro),
      GEN_RISK = input$gen,
      DOSE_AP = input$dose1, DOSE_AP2 = input$dose2, TSW_AP = input$tsw,
      F_LAI = input$flai, ESC_ON = as.numeric(input$esc),
      CYP2D6 = as.numeric(input$cyp),
      DOSE_AMA = input$ama, TSTART_AMA = if (input$ama > 0) 0 else 1e6,
      DOSE_ACH = input$ach, TSTART_ACH = if (input$ach > 0) 0 else 1e6,
      GKB_ON = as.numeric(input$gkb),
      TSTART_GKB = if (input$gkb) 0 else 1e6,
      BONT_ON = as.numeric(input$bont),
      TSTART_BONT = if (input$bont) input$tvmat else 1e6)
    if (input$vmat == "val") {
      p$DOSE_VAL <- input$vdose; p$TSTART_VAL <- input$tvmat
    } else if (input$vmat == "dtb") {
      p$DOSE_DTB <- input$vdose; p$TSTART_DTB <- input$tvmat
    }
    if (input$clz > 0) {
      p$DOSE_CLZ <- input$clz; p$TSTART_CLZ <- input$tsw
    }
    p
  })

  run <- reactive(do.call(sim, c(list(days = input$days), pset())))

  ## ---- tab 1 -------------------------------------------------------------
  output$profile <- renderText({
    d <- run(); n <- nrow(d)
    lat <- d$time[which(d$SDAM > 0.49)[1]]
    paste0(
      sprintf("Age %.0f · %s · genetic risk %.2f\n", input$age,
              if (input$fga) "FGA-level oxidative burden" else "SGA level", input$gen),
      sprintf("D2 occupancy: start %.3f → end %.3f | VMAT2 occupancy end %.3f\n",
              at_day(d, "OCCo", min(60, input$days)),
              d$OCCo[n], d$OCCVo[n]),
      sprintf("AIMS: 1yr %.2f · 2yr %.2f · end %.2f (max %.2f)\n",
              at_day(d, "AIMS", 365), at_day(d, "AIMS", 730), d$AIMS[n],
              max(d$AIMS)),
      sprintf("RUP %.3f · ROS %.3f · SDAM %.3f (latch threshold 0.50-0.55)\n",
              d$RUP[n], d$ROS[n], d$SDAM[n]),
      sprintf("Time of latch crossing: %s\n",
              ifelse(is.na(lat), "Not crossed (remains reversible)",
                     paste0(round(lat), " days"))),
      sprintf("PSYCH %.3f · PARK %.3f · DEPR %.3f · QTc %.2f ms · ADHER %.3f\n",
              d$PSYCH[n], d$PARK[n], d$DEPR[n], d$QTC[n], d$ADHER[n]),
      sprintf("Cumulative occupancy-days %.0f · plasticity-drive-days %.0f\n",
              d$CUMOCC[n],
              sum(diff(d$time) * (head(d$PLAST_DRIVE, -1) +
                                    tail(d$PLAST_DRIVE, -1)) / 2)))
  })

  output$p_summary <- renderPlot({
    d <- run()
    lineplot(d$time, list(d$AIMS, 26 * d$SDAM, 10 * d$RUP),
             c(PAL["aims"], PAL["sdam"], PAL["rup"]),
             c("AIMS (0-26)", "SDAM x 26", "RUP x 10"),
             "Index", "Disease trajectory summary",
             vlines = c(input$tsw, input$tvmat))
  })

  ## ---- tab 2 -------------------------------------------------------------
  output$p_pk <- renderPlot({
    d <- run()
    lineplot(d$time, list(d$C_APo, d$C_NBIo, d$C_HTBo, d$C_CLZo / 10),
             c(PAL["occ"], PAL["c"], PAL["b"], PAL["psy"]),
             c("Antipsychotic active moiety (ng/mL)", "NBI-98782 (ng/mL)",
               "(a+b)-HTBZ (ng/mL)", "Clozapine /10 (ng/mL)"),
             "Concentration (ng/mL)", "Drug concentration (Pharmacokinetics)",
             vlines = c(input$tsw, input$tvmat))
  })
  output$p_occ <- renderPlot({
    d <- run()
    lineplot(d$time, list(d$OCCo, d$OCCVo, d$PLAST_DRIVE),
             c(PAL["occ"], PAL["occv"], PAL["rup"]),
             c("Striatal D2 occupancy", "VMAT2 occupancy",
               "Plasticity drive Hill(OCC;0.70,4)"),
             "Fraction (-)", "Receptor occupancy and plasticity drive", ylim = c(0, 1),
             hlines = list(list(y = 0.65, col = "gray40"),
                           list(y = 0.80, col = "gray40")),
             vlines = c(input$tsw, input$tvmat))
    text(input$days * 0.02, 0.67, "Therapeutic threshold 65%", cex = 0.75, adj = 0)
    text(input$days * 0.02, 0.82, "EPS threshold 80%", cex = 0.75, adj = 0)
  })

  ## ---- tab 3 -------------------------------------------------------------
  output$p_plast <- renderPlot({
    d <- run()
    lineplot(d$time, list(d$RUP, d$SDAM, d$ROS / 4),
             c(PAL["rup"], PAL["sdam"], PAL["ros"]),
             c("RUP (D2 supersensitivity)", "SDAM (structural damage)", "ROS / 4"),
             "State (-)", "Two memories: reversible RUP vs quasi-irreversible SDAM",
             hlines = list(list(y = 0.52, col = PAL["sdam"])),
             vlines = c(input$tsw, input$tvmat))
    text(input$days * 0.02, 0.55, "Latch threshold SDAM* ~ 0.52", cex = 0.8,
         adj = 0, col = PAL["sdam"])
  })
  output$latchtxt <- renderText({
    d <- run(); n <- nrow(d)
    paste0(
      "The SDAM equation is bistable: 0 (repair-dominant) · threshold 0.50-0.55 · ",
      "upper 0.85-0.90 (self-sustaining).\n",
      sprintf("Current SDAM = %.3f → %s\n", d$SDAM[n],
              ifelse(d$SDAM[n] > 0.55, "Latched (full recovery impossible even with offending-drug withdrawal)",
                     ifelse(d$SDAM[n] > 0.30,
                            "Approaching threshold (now is the time to intervene)",
                            "Reversible zone"))),
      sprintf("Once the drive disappears, RUP decays with τ = %.0f days (slower when SDAM is higher).\n",
              300 / (1 + 0.8 * d$SDAM[n])))
  })

  ## ---- tab 4 -------------------------------------------------------------
  output$p_bg <- renderPlot({
    d <- run()
    lineplot(d$time, list(d$IND, d$GPE, d$STN, d$GPI, d$THAL),
             c(PAL["a"], PAL["b"], PAL["c"], PAL["d"], PAL["aims"]),
             c("Indirect-pathway MSN (IND)", "GPe", "STN", "GPi/SNr", "Thalamus THAL"),
             "Relative activity (-)", "Basal ganglia-thalamic loop",
             hlines = list(list(y = 1, col = "gray40")),
             vlines = c(input$tsw, input$tvmat))
  })
  output$p_exc <- renderPlot({
    d <- run()
    da_n <- d$DA_No
    stim <- da_n * (1 + 2.2 * d$RUP) * (1 - 0.55 * d$OCCo)
    lineplot(d$time, list(da_n, stim, pmax(stim - 1, 0) + 1.2 * d$SDAM,
                          1.2 * d$SDAM),
             c(PAL["occv"], PAL["rup"], PAL["aims"], PAL["sdam"]),
             c("Synaptic dopamine (baseline=1)", "D2STIM", "EXC (total overstimulation)",
               "Contribution of structural damage to EXC"),
             "Index (-)", "Decomposition of overstimulation: masking effect and structural floor",
             vlines = c(input$tsw, input$tvmat))
  })

  ## ---- tab 5 -------------------------------------------------------------
  output$p_aims <- renderPlot({
    d <- run()
    lineplot(d$time, list(d$AIMS), PAL["aims"], "AIMS total motor score (0-26)",
             "AIMS", "AIMS trajectory (includes measurement delay τ = 7 d)",
             vlines = c(input$tsw, input$tvmat))
  })
  output$p_ep <- renderPlot({
    d <- run()
    lineplot(d$time, list(d$PSYCH, d$PARK, d$DEPR, d$ADHER, d$QTC / 10),
             c(PAL["psy"], PAL["park"], PAL["depr"], PAL["c"], PAL["b"]),
             c("Psychosis burden PSYCH", "Parkinsonism PARK", "Mood burden DEPR",
               "Adherence ADHER", "QTc/10 (ms)"),
             "Index (-)", "Concomitant endpoints (competing harms)",
             vlines = c(input$tsw, input$tvmat))
  })

  ## ---- tab 6 -------------------------------------------------------------
  scen_runs <- reactive({
    base <- pset(); base$DOSE_AP2 <- base$DOSE_AP; base$TSW_AP <- 1e6
    base$DOSE_VAL <- NULL; base$DOSE_DTB <- NULL; base$DOSE_CLZ <- NULL
    tsw <- input$tsw
    arms <- list(
      "Continue" = list(),
      "50% dose reduction" = list(DOSE_AP2 = input$dose1 / 2, TSW_AP = tsw),
      "Complete discontinuation" = list(DOSE_AP2 = 0, TSW_AP = tsw),
      "Switch to clozapine" = list(DOSE_AP2 = 0, TSW_AP = tsw, DOSE_CLZ = 350,
                        TSTART_CLZ = tsw),
      "Add valbenazine 80 mg" = list(DOSE_VAL = 80, TSTART_VAL = tsw),
      "Clozapine + valbenazine" = list(DOSE_AP2 = 0, TSW_AP = tsw,
                            DOSE_CLZ = 350, TSTART_CLZ = tsw,
                            DOSE_VAL = 80, TSTART_VAL = tsw))
    lapply(arms, function(a)
      do.call(sim, c(list(days = input$days), modifyList(base, a))))
  })

  output$p_scen <- renderPlot({
    rs <- scen_runs()
    cols <- c("#7f8c8d", "#2471a3", "#c0392b", "#117864", "#b7950b",
              "#8e44ad")
    lineplot(rs[[1]]$time, lapply(rs, function(d) d$AIMS), cols, names(rs),
             "AIMS (0-26)", "AIMS trajectory by strategy (same patient, same timing)",
             vlines = input$tsw)
  })
  output$t_scen <- renderTable({
    rs <- scen_runs()
    data.frame(
      Strategy = names(rs),
      `AIMS at switch` = sapply(rs, function(d) at_day(d, "AIMS", input$tsw)),
      `AIMS +6wk` = sapply(rs, function(d) at_day(d, "AIMS", input$tsw + 42)),
      `AIMS end` = sapply(rs, function(d) tail(d$AIMS, 1)),
      `PSYCH end` = sapply(rs, function(d) tail(d$PSYCH, 1)),
      `PARK end` = sapply(rs, function(d) tail(d$PARK, 1)),
      `RUP end` = sapply(rs, function(d) tail(d$RUP, 1)),
      `SDAM end` = sapply(rs, function(d) tail(d$SDAM, 1)),
      check.names = FALSE)
  }, digits = 3)

  ## ---- tab 7 -------------------------------------------------------------
  window_scan <- reactive({
    g <- round(seq(input$wgrid[1], input$wgrid[2], length.out = 12))
    p <- pset(); p$DOSE_AP2 <- 0
    res <- lapply(g, function(s) {
      p$TSW_AP <- s
      d <- do.call(sim, c(list(days = s + 2190), p))
      post <- d[d$time >= s, ]
      data.frame(stop = s, sdam = at_day(d, "SDAM", s),
                 aims_stop = at_day(d, "AIMS", s),
                 peak = max(post$AIMS),
                 peak_day = post$time[which.max(post$AIMS)] - s,
                 aims6 = at_day(d, "AIMS", s + 2190))
    })
    do.call(rbind, res)
  })

  output$p_window <- renderPlot({
    w <- window_scan()
    plot(w$stop, w$aims6, type = "b", pch = 19, col = PAL["aims"], lwd = 2,
         xlab = "Timing of offending-drug withdrawal (day)",
         ylab = "AIMS 6 years after withdrawal", las = 1, bty = "l",
         main = "Window of reversibility: how late can withdrawal still help",
         ylim = c(0, max(w$aims6, w$peak)))
    grid(col = "gray90")
    lines(w$stop, w$peak, type = "b", pch = 1, col = PAL["park"], lwd = 2)
    lines(w$stop, w$aims_stop, type = "b", pch = 2, col = PAL["occ"],
          lwd = 1.5)
    abline(h = 2, lty = 3)
    legend("topleft", c("AIMS 6 years after withdrawal (residual)", "Peak AIMS after withdrawal (withdrawal-emergent worsening)",
                        "AIMS at withdrawal"),
           col = c(PAL["aims"], PAL["park"], PAL["occ"]), pch = c(19, 1, 2),
           lwd = 2, bty = "n", cex = 0.85)
  })
  output$t_window <- renderTable({
    w <- window_scan()
    w$outcome <- ifelse(w$aims6 > 2, "Persistent TD", "Recovery")
    names(w) <- c("Stop day", "SDAM at stop", "AIMS at stop", "Peak AIMS",
                  "Day of peak", "+6yr AIMS", "Outcome")
    w
  }, digits = 3)

  ## ---- tab 8 -------------------------------------------------------------
  frontier <- reactive({
    base <- pset(); base$DOSE_AP2 <- base$DOSE_AP; base$TSW_AP <- 1e6
    base$DOSE_VAL <- NULL; base$DOSE_DTB <- NULL; base$DOSE_CLZ <- NULL
    tsw <- input$tsw; readout <- min(tsw + 365, input$days)
    arms <- c(
      lapply(c(0, 0.25, 0.5, 0.75, 1), function(f)
        list(nm = sprintf("D2 blockade %.0f%% retained", 100 * (1 - f)),
             lever = "Offending-drug dose reduction",
             p = list(DOSE_AP2 = input$dose1 * (1 - f), TSW_AP = tsw))),
      lapply(c(20, 40, 80, 120), function(x)
        list(nm = sprintf("Valbenazine %d mg", x), lever = "Dopamine-supply blockade",
             p = list(DOSE_VAL = x, TSTART_VAL = tsw))),
      list(list(nm = "Switch to clozapine 350", lever = "Switch to low-occupancy agent",
                p = list(DOSE_AP2 = 0, TSW_AP = tsw, DOSE_CLZ = 350,
                         TSTART_CLZ = tsw))))
    do.call(rbind, lapply(arms, function(a) {
      d <- do.call(sim, c(list(days = readout + 5), modifyList(base, a$p)))
      data.frame(strategy = a$nm, lever = a$lever,
                 AIMS = at_day(d, "AIMS", readout),
                 PSYCH = at_day(d, "PSYCH", readout),
                 PARK = at_day(d, "PARK", readout),
                 DEPR = at_day(d, "DEPR", readout),
                 ADHER = at_day(d, "ADHER", readout),
                 RUP = at_day(d, "RUP", readout))
    }))
  })

  output$p_front <- renderPlot({
    f <- frontier()
    cols <- c("Offending-drug dose reduction" = PAL[["occ"]],
              "Dopamine-supply blockade" = PAL[["occv"]],
              "Switch to low-occupancy agent" = PAL[["aims"]])
    plot(f$AIMS, f$PSYCH, pch = 19, cex = 1.6, col = cols[f$lever],
         xlab = "AIMS (lower is better)", ylab = "Psychosis burden PSYCH",
         main = "Opposing levers: the same AIMS reduction pays a different price",
         las = 1, bty = "l", xlim = range(f$AIMS) + c(-1, 1),
         ylim = c(0, max(f$PSYCH) + 0.1))
    grid(col = "gray90")
    text(f$AIMS, f$PSYCH, f$strategy, pos = 3, cex = 0.78)
    points(f$AIMS, f$PARK, pch = 1, cex = 1.4, col = cols[f$lever])
    legend("topright", c(names(cols), "Open circle = parkinsonism PARK"),
           col = c(cols, "black"), pch = c(19, 19, 19, 1), bty = "n",
           cex = 0.85)
  })
  output$t_front <- renderTable(frontier(), digits = 3)

  ## ---- tab 9 -------------------------------------------------------------
  output$t_anchor <- renderTable({
    data.frame(
      `Calibration target` = c(
        "Risperidone active-moiety PK", "Striatal D2 occupancy EC50",
        "Therapeutic/EPS occupancy thresholds", "Clozapine D2 occupancy",
        "Valbenazine 80 mg exposure", "Deutetrabenazine 36 mg exposure",
        "KINECT-3 (valbenazine 80 mg, 6 wk)",
        "ARM-TD/AIM-TD (deutetrabenazine 36 mg)",
        "Annual TD incidence", "Timing of withdrawal-emergent dyskinesia",
        "Persistence after discontinuation", "Ginkgo biloba EGb761 240 mg"),
      `Literature value` = c(
        "CL/F ~ 5 L/h, t1/2 ~ 20 h, 4 mg → ~33 ng/mL",
        "~12 ng/mL (4 mg ~ 73%, 8 mg ~ 86%)",
        "65% / 78-80%", "350 mg → ~390 ng/mL, D2 ~ 20-40%",
        "NBI-98782 Cave ~ 24 ng/mL (about 2x in PM)",
        "(a+b)-HTBZ ~ 15 ng/mL, QTc ~ 4 ms",
        "AIMS -3.2 (baseline ~10)", "AIMS ~ -3.0",
        "FGA 5-6%/yr, SGA 3-4%/yr, 25-30%/yr in age ≥55",
        "Peaks within 4-6 weeks after dose reduction",
        "About 1/3 recover, the rest persist for years", "AIMS improvement (RCT)"),
      `Model parameter / output` = c(
        "CL_AP 120 L/day, V2 100 L", "EC50_D2 = 12, HILL_D2 = 1.25",
        "OCC50_P = 0.58, OCC50_PK = 0.78",
        "CL_CLZ 900 L/day, EC50_D2_CLZ = 900 → OCC 0.28",
        "FM_VAL 0.30, CL_NBI 1000 L/day, CYP2D6 x0.5",
        "FM_DTB 0.50, CL_HTB 1200 L/day, QT_HTB 0.25",
        "Model: -3.51 (-33%), baseline 10.59",
        "Model: -3.09 (-29%) at 6 wk",
        "RISK_FGA 1.6, RISKMOD = 1+0.015(age-40)",
        "Model peak time = 30 days after withdrawal",
        "SDAM bistable (threshold 0.52), KOUT_S = 3e-4",
        "ANTIOX_MAX 0.85 → latch 400 days → 780 days"),
      check.names = FALSE)
  })
  output$disclaim <- renderText(paste(
    "All values in this dashboard are outputs of a QSP model for educational/research purposes.",
    "\nThey must not be used for clinical decision-making, prescribing, or regulatory submission.",
    "\nSee td_references.md for the literature list, and",
    "td_reference_check.py for independent numerical verification."))
}

shinyApp(ui, server)
