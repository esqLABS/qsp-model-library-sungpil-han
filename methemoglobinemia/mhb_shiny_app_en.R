# =============================================================================
#  mhb_shiny_app.R
#  Methaemoglobinaemia QSP model — interactive dashboard
#  Methaemoglobinaemia QSP model — interactive dashboard
# =============================================================================
#
#  WHAT THIS APP IS FOR
#  --------------------
#  Most dashboards in this library let you turn a knob and watch a curve.  This
#  one is built around a single disagreement, and almost every tab exists to
#  make that disagreement visible:
#
#        the number on the monitor,
#        the number in the blood gas,
#        and the number at the tissue
#        are three different numbers, and they move apart as the patient
#        gets worse.
#
#  Tab 3 (the two liars) puts SpO2, co-oximetry SaO2 and PvO2 on one figure.
#  Tab 4 converts %MetHb into the anaemia it is haemodynamically equivalent to.
#  Tab 7 shows the antidote's dose-response bending back on itself.
#  Tab 8 shows why it does that, in the currency the model actually spends.
#
#  RUN
#  ---
#    install.packages(c("shiny","mrgsolve","dplyr","tidyr","ggplot2","DT","gridExtra"))
#    shiny::runApp("mhb_shiny_app.R")
#
#  The model file `mhb_mrgsolve_model_en.R` must sit beside this one.
# =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

mod <- mread_cache("mhb", "mhb_mrgsolve_model_en.R")

MWMB <- 319.85
mb_ev <- function(mgkg, t = 0, wt = 70)
  ev(time = t, amt = mgkg * wt * 1000 / MWMB, cmt = "MBC")

PAL <- c(SpO2 = "#2b6cb0", `co-ox SaO2` = "#b3541e", PvO2 = "#2f6b34",
         MetHb = "#8b2f39", Hb = "#5b3a80")
thm <- theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        legend.position = "bottom")

# ---------------------------------------------------------------------------
#  build the event sequence and parameter set from the UI
# ---------------------------------------------------------------------------
build_sim <- function(input, tend = NULL) {
  p <- list(
    HB0    = input$hb,
    G6PD   = input$g6pd,
    EB5SET = input$eb5,
    FHBM   = input$fhbm,
    PAO2   = if (input$hbo) input$pao2 else 95,
    VO2    = input$vo2,
    ALPHAM = input$alpham
  )
  e <- NULL
  add <- function(a, b) if (is.null(a)) b else c(a, b)

  if (input$oxidant == "benzocaine")
    e <- add(e, ev(time = 0, amt = input$bzc, cmt = "BZCD"))
  if (input$oxidant == "dapsone_chronic")
    e <- add(e, ev(time = 0, amt = input$dapdose, cmt = "DAPG",
                   ii = 24, addl = max(0, ceiling(input$tend / 24) - 1)))
  if (input$oxidant == "dapsone_od")
    e <- add(e, ev(time = 0, amt = input$dapod, cmt = "DAPG"))
  if (input$oxidant == "nitrite")
    e <- add(e, ev(time = 0, amt = input$nit * 1000 * 1000 / 69, cmt = "NITG"))

  if (input$cimetidine)
    e <- add(e, ev(time = 0, amt = 400, cmt = "CIMG", ii = 8,
                   addl = max(0, ceiling(input$tend / 8) - 1)))
  if (input$ssri)
    e <- add(e, ev(time = 0, amt = 2.5, cmt = "SSRI"))
  if (input$ascorbate)
    e <- add(e, ev(time = input$mbtime, amt = 10 * 1e6 / 176.1, cmt = "ASCC"))

  if (input$nmb > 0)
    for (k in seq_len(input$nmb))
      e <- add(e, mb_ev(input$mbdose, input$mbtime + (k - 1) * input$mbint))

  list(p = p, e = e, tend = if (is.null(tend)) input$tend else tend)
}

run_sim <- function(s, delta = 0.05) {
  m <- param(mod, s$p)
  if (is.null(s$e)) mrgsim(m, end = s$tend, delta = delta) %>% as_tibble()
  else mrgsim(m, events = s$e, end = s$tend, delta = delta) %>% as_tibble()
}

# EQUIV Hb: invert the oxygen-transport block for a patient with normal Hb
equiv_hb <- function(pvo2, p50 = 26.8, n = 2.7, vo2 = 250, co = 5) {
  f <- function(hb) {
    sao2 <- 95^n / (p50^n + 95^n)
    cao2 <- 1.34 * hb * sao2 + 0.003 * 95
    cvo2 <- cao2 - vo2 / (co * 10)
    svo2 <- pmin(pmax(cvo2 / (1.34 * hb), 1e-6), 1 - 1e-6)
    p50 * (svo2 / (1 - svo2))^(1 / n) - pvo2
  }
  out <- numeric(length(pvo2))
  for (i in seq_along(pvo2)) {
    out[i] <- tryCatch(uniroot(f, c(0.2, 30))$root, error = function(e) NA_real_)
  }
  out
}

# static oxygen-transport block, for the sweep tabs (no integration needed)
transport <- function(fmet, hb = 15, p50 = 26.8, nh = 2.7, alpham = 0.5,
                      pao2 = 95, vo2 = 250, co = 5, bpg = 1, bbpg = 0.55,
                      e660m = 15, e940m = 15) {
  hbf <- hb * (1 - fmet)
  fred <- fmet
  neff <- 1 + (nh - 1) * (1 - fred)
  p50e <- pmax(p50 * (1 + bbpg * (bpg - 1)) * (1 - alpham * fred), 1)
  sao2 <- pao2^neff / (p50e^neff + pao2^neff)
  cao2 <- 1.34 * hbf * sao2 + 0.003 * pao2
  cvo2 <- cao2 - vo2 / (co * 10)
  caph <- 1.34 * hbf
  svo2 <- ifelse(cvo2 >= caph, 1, pmin(pmax(cvo2 / caph, 1e-6), 1 - 1e-6))
  pvo2 <- ifelse(cvo2 >= caph, (cvo2 - caph) / 0.003,
                 p50e * (svo2 / (1 - svo2))^(1 / neff))
  fO2 <- hbf * sao2 / hb; fdx <- hbf * (1 - sao2) / hb
  num <- fO2 * 0.32 + fdx * 3.23 + fmet * e660m
  den <- fO2 * 1.21 + fdx * 0.69 + fmet * e940m
  R <- num / den
  spo2 <- pmin(pmax(110 - 25 * R, 0), 100)
  tibble(fmet = fmet, MetPct = 100 * fmet, hbf = hbf, P50E = p50e, NEFF = neff,
         SAO2 = 100 * sao2, PVO2 = pmax(pvo2, 0), R = R, SPO2 = spo2,
         SAO2CO = 100 * fO2, GAP = spo2 - 100 * fO2,
         CYAN = (hbf * (1 - sao2)) / 5 + hb * fmet / 1.5)
}

# =============================================================================
#  UI
# =============================================================================
ui <- fluidPage(
  titlePanel("Methaemoglobinaemia QSP model"),
  tags$p(style = "color:#555;margin-top:-8px",
         HTML("<b>It is a disease of delivery, not of saturation.</b> ",
              "The number on the monitor · the number in the blood gas · the number in the tissue are different numbers, ",
              "and they move further apart as the patient worsens.")),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Patient"),
      sliderInput("hb", "Total haemoglobin Hb (g/dL)", 5, 19, 15, 0.5),
      sliderInput("g6pd", "G6PD activity (fraction)", 0.01, 1.2, 1.0, 0.01),
      helpText("0.02 = Mediterranean (class I–II) · 0.15 = A− (class III)"),
      sliderInput("eb5", "CYB5R3 activity (in vivo flux fraction)", 0.01, 1.2, 1.0, 0.005),
      helpText("0.035 = congenital type I phenotype"),
      sliderInput("fhbm", "HbM variant fraction (non-reducible)", 0, 0.4, 0, 0.05),
      sliderInput("vo2", "Oxygen consumption VO₂ (mL/min)", 150, 600, 250, 10),

      h4("Oxidant exposure"),
      selectInput("oxidant", NULL,
                  c("None" = "none",
                    "Benzocaine spray" = "benzocaine",
                    "Dapsone chronic administration" = "dapsone_chronic",
                    "Dapsone overdose" = "dapsone_od",
                    "Sodium nitrite ingestion" = "nitrite"),
                  selected = "benzocaine"),
      conditionalPanel("input.oxidant == 'benzocaine'",
        sliderInput("bzc", "Absorbed amount (mg)", 20, 600, 250, 10)),
      conditionalPanel("input.oxidant == 'dapsone_chronic'",
        sliderInput("dapdose", "Daily dose (mg)", 25, 400, 100, 25)),
      conditionalPanel("input.oxidant == 'dapsone_od'",
        sliderInput("dapod", "Ingested dose (mg)", 500, 6000, 2000, 250)),
      conditionalPanel("input.oxidant == 'nitrite'",
        sliderInput("nit", "NaNO₂ (g)", 0.5, 20, 5, 0.5)),

      h4("Treatment"),
      sliderInput("nmb", "Number of methylene blue doses", 0, 6, 1, 1),
      sliderInput("mbdose", "Dose per administration (mg/kg)", 0.25, 10, 1, 0.25),
      sliderInput("mbtime", "Time of first dose (h)", 0, 24, 0.75, 0.25),
      sliderInput("mbint", "Dosing interval (h)", 0.5, 12, 4, 0.5),
      checkboxInput("cimetidine", "Cimetidine 400 mg tid (formation-inhibiting term)", FALSE),
      checkboxInput("ascorbate", "Ascorbic acid 10 g IV", FALSE),
      checkboxInput("hbo", "Hyperbaric oxygen / high-concentration oxygen", FALSE),
      conditionalPanel("input.hbo == true",
        sliderInput("pao2", "PaO₂ (mmHg)", 95, 2000, 1800, 25)),
      checkboxInput("ssri", "Concomitant SSRI (serotonin tab)", FALSE),

      h4("Model assumption"),
      sliderInput("alpham", "Allosteric left shift α", 0, 1.0, 0.5, 0.05),
      helpText(HTML("<b>This slider is the model's Achilles' heel.</b> If α = 0, ",
                    "%MetHb becomes a sufficient bedside variable, and this model's central claim disappears.")),
      sliderInput("tend", "Simulation duration (h)", 6, 720, 48, 6)
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        tabPanel("1. Patient profile",
          br(), DTOutput("profile"), br(),
          htmlOutput("verdict")),

        tabPanel("2. Time course",
          br(), plotOutput("tc_main", height = "620px")),

        tabPanel("3. The two liars",
          br(),
          tags$p(HTML("Blue is the <b>monitor</b>, orange is <b>co-oximetry</b>, ",
                      "green is the <b>tissue</b>. The gap that opens between the three is this disease.")),
          plotOutput("liars", height = "560px")),

        tabPanel("4. Equivalent anaemia",
          br(),
          tags$p(HTML("\"MetHb 30%\" is translated into <b>the normal haemoglobin concentration that produces the same tissue PO₂</b>. ",
                      "The dashed line is the conventional arithmetic Hb×(1−f), the solid line is the actual value. ",
                      "The gap between the two is the haemoglobin neutralised by the allosteric left shift.")),
          plotOutput("equiv", height = "520px"), br(), DTOutput("equivtab")),

        tabPanel("5. Pulse oximeter",
          br(),
          tags$p(HTML("The 85% floor is not a property of the blood but <b>a value from the calibration line</b>. ",
                      "The slope in the bottom-right figure below is the point of this tab.")),
          plotOutput("oxi", height = "540px")),

        tabPanel("6. Clinical endpoints",
          br(), plotOutput("endpoints", height = "520px"), br(),
          DTOutput("endtab")),

        tabPanel("7. Methylene blue dose–response",
          br(),
          tags$p(HTML("The model has <b>no</b> upper-dose-limit parameter. ",
                      "And yet the curve bends back on itself.")),
          plotOutput("mbdr", height = "560px")),

        tabPanel("8. The NADPH economy",
          br(),
          tags$p(HTML("The currency the model actually spends. Methylene blue and glutathione draw electrons ",
                      "from <b>the same wallet</b>.")),
          plotOutput("nadph", height = "560px")),

        tabPanel("9. G6PD gradient",
          br(), plotOutput("g6pdplot", height = "540px"), br(),
          DTOutput("g6pdtab")),

        tabPanel("10. Scenario comparison",
          br(), plotOutput("scen", height = "620px")),

        tabPanel("11. Biomarkers · haemolysis",
          br(), plotOutput("biomark", height = "600px")),

        tabPanel("12. Serotonin",
          br(), plotOutput("sero", height = "460px")),

        tabPanel("13. Model notes",
          br(), htmlOutput("notes"))
      )
    )
  )
)

# =============================================================================
#  SERVER
# =============================================================================
server <- function(input, output, session) {

  sim <- reactive({ run_sim(build_sim(input)) })

  # ---- 1. patient profile --------------------------------------------------
  output$profile <- renderDT({
    d <- sim(); pk <- d[which.max(d$MetPct), ]; en <- d[nrow(d), ]
    tibble(
      Item = c("Peak MetHb (%)", "Time of peak (h)", "SpO₂ at that time (%)",
               "Co-oximetry SaO₂ at that time (%)", "Saturation gap (points)",
               "Lowest PvO₂ (mmHg)", "Lowest SpO₂ (%)",
               "Cumulative hypoxic burden (mmHg·h)", "Lowest total Hb (g/dL)",
               "Peak Heinz body burden", "Peak lactate (mmol/L)",
               "Cyanosis index (≥1 = visible cyanosis)"),
      Value = round(c(pk$MetPct, pk$time, pk$SPO2, pk$SAO2CO, pk$SPO2 - pk$SAO2CO,
                   min(d$PVO2), min(d$SPO2), en$INJ, min(d$HBTOT),
                   max(d$HEINZ), max(d$LAC), max(d$CYAN)), 2))
  }, options = list(dom = "t", pageLength = 20), rownames = FALSE)

  output$verdict <- renderUI({
    d <- sim()
    pv <- min(d$PVO2); mx <- max(d$MetPct); hb <- min(d$HBTOT)
    msg <- c()
    if (pv < 20) msg <- c(msg, sprintf(
      "<span style='color:#b3541e'><b>Tissue PO₂ falls below the critical threshold (lowest %.1f mmHg).</b>
       In this patient, %%MetHb %.0f%% already means anaerobic metabolism.</span>", pv, mx))
    if (mx > 15 && pv >= 20) msg <- c(msg, sprintf(
      "MetHb rises to %.0f%%, but tissue PO₂ remains above the critical threshold at %.1f mmHg.", mx, pv))
    if (hb < input$hb - 1.5) msg <- c(msg, sprintf(
      "<span style='color:#8b2f39'><b>Haemolysis is occurring</b> (Hb %.1f → %.1f g/dL).
       A falling %%MetHb should not be read as recovery — oxidised red cells are being cleared first.</span>",
      input$hb, hb))
    if (input$g6pd < 0.1 && input$nmb > 0) msg <- c(msg,
      "<span style='color:#8b2f39'><b>Methylene blue was administered with very low G6PD activity.</b>
       The model has no rule forbidding this combination — the result below is computed from a single shared NADPH ceiling.</span>")
    HTML(paste0("<div style='font-size:15px;line-height:1.7'>",
                paste(msg, collapse = "<br>"), "</div>"))
  })

  # ---- 2. time course ------------------------------------------------------
  output$tc_main <- renderPlot({
    d <- sim() %>% select(time, MetPct, SPO2, SAO2CO, PVO2, HBTOT, LAC, CO,
                          CMBRBC, GSH) %>%
      pivot_longer(-time)
    lab <- c(MetPct = "MetHb (%)", SPO2 = "SpO₂ (%)", SAO2CO = "co-ox SaO₂ (%)",
             PVO2 = "PvO₂ (mmHg)", HBTOT = "Total Hb (g/dL)", LAC = "Lactate (mmol/L)",
             CO = "Cardiac output (L/min)", CMBRBC = "Erythrocyte MB (µM)",
             GSH = "GSH (µM)")
    d$name <- factor(lab[d$name], levels = lab)
    ggplot(d, aes(time, value)) +
      geom_line(linewidth = 0.9, colour = "#2b6cb0") +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "Time (h)", y = NULL) + thm
  })

  # ---- 3. the two liars ----------------------------------------------------
  output$liars <- renderPlot({
    d <- sim()
    g1 <- d %>% select(time, SpO2 = SPO2, `co-ox SaO2` = SAO2CO) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
      geom_line(linewidth = 1.1) +
      scale_colour_manual(values = PAL, name = NULL) +
      labs(title = "What the monitor sees vs what the blood actually carries",
           x = "Time (h)", y = "%") + thm
    g2 <- ggplot(d, aes(time, PVO2)) +
      geom_hline(yintercept = 20, linetype = "dashed", colour = "#b3541e") +
      annotate("text", x = max(d$time) * 0.75, y = 21.5,
               label = "Critical capillary PO₂", colour = "#b3541e", size = 3.6) +
      geom_line(linewidth = 1.1, colour = PAL[["PvO2"]]) +
      labs(title = "What the tissue actually sees", x = "Time (h)", y = "PvO₂ (mmHg)") + thm
    g3 <- ggplot(d, aes(time, SPO2 - SAO2CO)) +
      geom_line(linewidth = 1.1, colour = "#5b3a80") +
      geom_hline(yintercept = 5, linetype = "dotted") +
      labs(title = "Saturation gap (>5 suggests an abnormal haemoglobin)",
           x = "Time (h)", y = "SpO₂ − SaO₂ (points)") + thm
    gridExtra::grid.arrange(g1, g2, g3, ncol = 1)
  })

  # ---- 4. equivalent anaemia ----------------------------------------------
  eqdata <- reactive({
    tr <- transport(seq(0, 0.75, by = 0.01), hb = input$hb,
                    alpham = input$alpham, vo2 = input$vo2)
    tr$naive <- input$hb * (1 - tr$fmet)
    tr$equiv <- equiv_hb(tr$PVO2, vo2 = input$vo2)
    tr$lost <- tr$naive - tr$equiv
    tr
  })

  output$equiv <- renderPlot({
    tr <- eqdata()
    ggplot(tr, aes(MetPct)) +
      geom_ribbon(aes(ymin = equiv, ymax = naive), fill = "#f9d6c8", alpha = .8) +
      geom_line(aes(y = naive), linetype = "dashed", linewidth = 1) +
      geom_line(aes(y = equiv), linewidth = 1.2, colour = "#8b2f39") +
      annotate("text", x = 45, y = max(tr$naive) * 0.75,
               label = "This band = haemoglobin neutralised by the left shift",
               colour = "#8b2f39", size = 4.4) +
      labs(title = "Translating %MetHb into an equivalent anaemia",
           subtitle = sprintf("Hb %.1f g/dL, VO₂ %d mL/min, α = %.2f",
                              input$hb, input$vo2, input$alpham),
           x = "Methaemoglobin (%)", y = "Haemoglobin (g/dL)") + thm
  })

  output$equivtab <- renderDT({
    tr <- eqdata() %>% filter(round(MetPct) %in% c(0, 10, 15, 20, 25, 30, 40, 50, 60))
    tr %>% transmute(`MetHb (%)` = round(MetPct),
                     `Conventional arithmetic Hb×(1−f)` = round(naive, 2),
                     `Equivalent Hb` = round(equiv, 2),
                     `Neutralised Hb` = round(lost, 2),
                     `% of residual capacity` = round(100 * lost / naive, 1),
                     `PvO₂` = round(PVO2, 1),
                     `SpO₂` = round(SPO2, 1))
  }, options = list(dom = "t"), rownames = FALSE)

  # ---- 5. the oximeter -----------------------------------------------------
  output$oxi <- renderPlot({
    tr <- transport(seq(0, 0.9, by = 0.005), hb = input$hb, alpham = input$alpham)
    tr$slope <- c(NA, diff(tr$SPO2) / diff(tr$MetPct))
    g1 <- tr %>% select(MetPct, SpO2 = SPO2, `co-ox SaO2` = SAO2CO) %>%
      pivot_longer(-MetPct) %>%
      ggplot(aes(MetPct, value, colour = name)) +
      geom_hline(yintercept = 85, linetype = "dotted") +
      geom_line(linewidth = 1.2) +
      scale_colour_manual(values = PAL, name = NULL) +
      annotate("text", x = 60, y = 87.5, label = "110 − 25×1 = 85", size = 3.8) +
      labs(title = "The floor is a property of the instrument", x = "MetHb (%)", y = "%") + thm
    g2 <- ggplot(tr, aes(MetPct, R)) +
      geom_hline(yintercept = 1, linetype = "dotted") +
      geom_line(linewidth = 1.2, colour = "#5b3a80") +
      labs(title = "The only thing the instrument calculates: R = A₆₆₀/A₉₄₀",
           x = "MetHb (%)", y = "R") + thm
    g3 <- ggplot(tr[-1, ], aes(MetPct, -slope)) +
      geom_line(linewidth = 1.2, colour = "#b3541e") +
      labs(title = "Sensitivity collapse: −dSpO₂ / d(%MetHb)",
           subtitle = "In the range where the patient is dying, the instrument has effectively stopped responding",
           x = "MetHb (%)", y = "points / %") + thm
    g4 <- ggplot(tr, aes(MetPct, GAP)) +
      geom_line(linewidth = 1.2, colour = "#2f6b34") +
      labs(title = "Saturation gap", x = "MetHb (%)", y = "points") + thm
    gridExtra::grid.arrange(g1, g2, g3, g4, ncol = 2)
  })

  # ---- 6. clinical endpoints ----------------------------------------------
  output$endpoints <- renderPlot({
    grid <- expand.grid(hb = seq(6, 18, by = 0.5), fmet = seq(0, 0.7, by = 0.01))
    res <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i)
      transport(grid$fmet[i], hb = grid$hb[i], alpham = input$alpham,
                vo2 = input$vo2)))
    res$hb <- grid$hb
    ggplot(res, aes(MetPct, hb, fill = PVO2)) +
      geom_raster() +
      geom_contour(aes(z = PVO2), breaks = 20, colour = "white", linewidth = 1.2) +
      geom_contour(aes(z = CYAN), breaks = 1, colour = "#333333",
                   linetype = "dashed") +
      geom_point(aes(x = 30, y = 15), colour = "red", size = 3) +
      scale_fill_viridis_c(name = "PvO₂ (mmHg)") +
      labs(title = "The same percentage is not the same disease",
           subtitle = paste("White line = anaerobic threshold (PvO₂ 20 mmHg) · dotted line = visible cyanosis threshold ·",
                            "red dot = the textbook '30%'"),
           x = "Methaemoglobin (%)", y = "Total haemoglobin (g/dL)") + thm
  })

  output$endtab <- renderDT({
    do.call(rbind, lapply(c(17, 15, 12, 9, 7), function(h)
      do.call(rbind, lapply(c(0.15, 0.20, 0.30, 0.40), function(f) {
        t <- transport(f, hb = h, alpham = input$alpham, vo2 = input$vo2)
        tibble(`Hb` = h, `MetHb (%)` = 100 * f,
               `Functional Hb` = round(t$hbf, 2), `PvO₂` = round(t$PVO2, 1),
               `SpO₂` = round(t$SPO2, 1),
               `Anaerobic?` = ifelse(t$PVO2 < 20, "Yes", "No"))
      }))))
  }, options = list(dom = "tp", pageLength = 10), rownames = FALSE)

  # ---- 7. methylene blue dose-response ------------------------------------
  output$mbdr <- renderPlot({
    doses <- c(0, 0.5, 1, 1.5, 2, 3, 4, 5, 7, 10)
    res <- lapply(doses, function(dz) {
      s <- build_sim(input, tend = 48)
      s$e <- NULL
      if (input$oxidant == "benzocaine")
        s$e <- ev(time = 0, amt = input$bzc, cmt = "BZCD")
      if (dz > 0) s$e <- if (is.null(s$e)) mb_ev(dz, 0.75) else c(s$e, mb_ev(dz, 0.75))
      d <- run_sim(s, delta = 0.1)
      tibble(dose = dz,
             `MetHb @2h (%)` = d$MetPct[which.min(abs(d$time - 2))],
             `Lowest MetHb (%)` = min(d$MetPct[d$time > 0.75]),
             `Hb @48h (g/dL)` = d$HBTOT[nrow(d)],
             `Heinz body peak` = max(d$HEINZ),
             `ROS exposure (AUC)` = sum(d$ROS) * 0.1)
    }) %>% bind_rows()
    res %>% pivot_longer(-dose) %>%
      ggplot(aes(dose, value)) +
      geom_vline(xintercept = 7, linetype = "dotted", colour = "#b3541e") +
      geom_line(linewidth = 1.1, colour = "#2b6cb0") + geom_point(size = 2) +
      facet_wrap(~name, scales = "free_y") +
      labs(title = "The curve bends back even with no upper-limit parameter",
           subtitle = "Dotted line = the clinically known cumulative 7 mg/kg ceiling (not entered into the model)",
           x = "Methylene blue (mg/kg)", y = NULL) + thm
  })

  # ---- 8. the NADPH economy -----------------------------------------------
  output$nadph <- renderPlot({
    d <- sim()
    d$leuco <- d$LEUCOF
    g1 <- ggplot(d, aes(time, CMBRBC)) + geom_line(linewidth = 1.1) +
      labs(title = "Total methylene blue in erythrocytes", x = "h", y = "µM") + thm
    g2 <- ggplot(d, aes(time, leuco)) + geom_line(linewidth = 1.1, colour = "#2f6b34") +
      labs(title = "Leuco fraction — is NADPH keeping up?",
           subtitle = "When this value falls, what remains is oxidant",
           x = "h", y = "LMB / (LMB+MB)") + thm
    g3 <- ggplot(d, aes(time, GSH)) + geom_line(linewidth = 1.1, colour = "#5b3a80") +
      geom_hline(yintercept = 1000, linetype = "dotted") +
      labs(title = "Reduced glutathione — the share the membrane spends",
           x = "h", y = "µM") + thm
    g4 <- ggplot(d, aes(time, HEINZ)) + geom_line(linewidth = 1.1, colour = "#8b2f39") +
      labs(title = "Heinz body burden → haemolysis", x = "h", y = "0–1") + thm
    gridExtra::grid.arrange(g1, g2, g3, g4, ncol = 2)
  })

  # ---- 9. G6PD gradient ---------------------------------------------------
  g6data <- reactive({
    gs <- c(1.0, 0.6, 0.4, 0.25, 0.15, 0.08, 0.02)
    lapply(gs, function(g) {
      s <- build_sim(input, tend = 72); s$p$G6PD <- g
      d <- run_sim(s, delta = 0.1)
      tibble(G6PD = g,
             `MetHb @2h (%)` = d$MetPct[which.min(abs(d$time - 2))],
             `Lowest GSH (µM)` = min(d$GSH),
             `Heinz body peak` = max(d$HEINZ),
             `Hb @72h (g/dL)` = d$HBTOT[nrow(d)])
    }) %>% bind_rows()
  })

  output$g6pdplot <- renderPlot({
    g6data() %>% pivot_longer(-G6PD) %>%
      ggplot(aes(G6PD, value)) +
      geom_line(linewidth = 1.1, colour = "#2b6cb0") + geom_point(size = 2) +
      scale_x_log10() + facet_wrap(~name, scales = "free_y") +
      labs(title = "The model has no rule 'contraindicated in G6PD deficiency'",
           subtitle = "This whole gradient comes from a single shared NADPH ceiling",
           x = "G6PD activity (fraction, log)", y = NULL) + thm
  })
  output$g6pdtab <- renderDT({ g6data() %>% mutate(across(where(is.numeric), ~round(.x, 3))) },
                             options = list(dom = "t"), rownames = FALSE)

  # ---- 10. scenario comparison --------------------------------------------
  output$scen <- renderPlot({
    base <- ev(time = 0, amt = input$bzc, cmt = "BZCD")
    scn <- list(
      "No treatment" = list(e = base, p = list()),
      "MB 1 mg/kg" = list(e = c(base, mb_ev(1, 0.75)), p = list()),
      "MB 5 mg/kg" = list(e = c(base, mb_ev(5, 0.75)), p = list()),
      "MB + G6PD 2%" = list(e = c(base, mb_ev(2, 0.75)), p = list(G6PD = 0.02)),
      "Ascorbate" = list(e = c(base, ev(time = 0.75, amt = 10 * 1e6 / 176.1,
                                          cmt = "ASCC")), p = list()),
      "Hyperbaric oxygen 2.8 ATA" = list(e = base, p = list(PAO2 = 1800)),
      "Anaemia Hb 8 · no treatment" = list(e = base, p = list(HB0 = 8))
    )
    d <- bind_rows(lapply(names(scn), function(n) {
      m <- param(mod, modifyList(list(HB0 = input$hb, ALPHAM = input$alpham),
                                 scn[[n]]$p))
      mrgsim(m, events = scn[[n]]$e, end = 24, delta = 0.1) %>% as_tibble() %>%
        mutate(scenario = n)
    }))
    d %>% select(time, scenario, MetPct, SPO2, PVO2, HBTOT) %>%
      pivot_longer(c(MetPct, SPO2, PVO2, HBTOT)) %>%
      ggplot(aes(time, value, colour = scenario)) +
      geom_line(linewidth = 1) +
      facet_wrap(~name, scales = "free_y") +
      labs(title = "Seven scenarios, the same exposure",
           x = "Time (h)", y = NULL, colour = NULL) + thm
  })

  # ---- 11. biomarkers ------------------------------------------------------
  output$biomark <- renderPlot({
    d <- sim() %>% select(time, HBTOT, MetPct, HEINZ, FHB, HAPT, BILI, LDH,
                          GSH, LAC) %>% pivot_longer(-time)
    lab <- c(HBTOT = "Total Hb (g/dL)", MetPct = "MetHb (%)", HEINZ = "Heinz body burden",
             FHB = "Plasma free Hb (µM haem)", HAPT = "Haptoglobin (µM)",
             BILI = "Indirect bilirubin (mg/dL)", LDH = "LDH (U/L)",
             GSH = "GSH (µM)", LAC = "Lactate (mmol/L)")
    d$name <- factor(lab[d$name], levels = lab)
    ggplot(d, aes(time, value)) + geom_line(linewidth = 0.95, colour = "#8b2f39") +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(title = "Haemolysis markers — what to look at when %MetHb looks reassuring",
           x = "Time (h)", y = NULL) + thm
  })

  # ---- 12. serotonin -------------------------------------------------------
  output$sero <- renderPlot({
    variants <- list(
      "MB alone" = list(ssri = FALSE, mb = TRUE),
      "SSRI alone" = list(ssri = TRUE,  mb = FALSE),
      "MB + SSRI" = list(ssri = TRUE,  mb = TRUE))
    d <- bind_rows(lapply(names(variants), function(n) {
      v <- variants[[n]]
      e <- ev(time = 0, amt = input$bzc, cmt = "BZCD")
      if (v$ssri) e <- c(e, ev(time = 0, amt = 2.5, cmt = "SSRI"))
      if (v$mb)   e <- c(e, mb_ev(input$mbdose, 0.75))
      mrgsim(param(mod, list(HB0 = input$hb)), events = e, end = 12, delta = 0.05) %>%
        as_tibble() %>% mutate(scenario = n)
    }))
    ggplot(d, aes(time, HT5, colour = scenario)) +
      geom_hline(yintercept = 1, linetype = "dotted") +
      geom_line(linewidth = 1.2) +
      labs(title = "The pharmacology of the antidote itself: methylene blue is a potent MAO-A inhibitor",
           x = "Time (h)", y = "Synaptic serotonin (normal = 1)", colour = NULL) + thm
  })

  # ---- 13. notes -----------------------------------------------------------
  output$notes <- renderUI(HTML("
<div style='max-width:980px;line-height:1.75;font-size:15px'>
<h3>What this app is trying to show</h3>
<p>This model's parameter list contains no <b>severity scale</b>, no <b>7 mg/kg ceiling</b>,
no <b>G6PD contraindication switch</b>, no <b>rebound term</b>, and no <b>85% floor constant</b>.
All five of those are <b>outputs</b>. What went in is only five structures:</p>
<ol>
<li>A <b>catalytic co-oxidation cycle</b> whose coenzyme is glutathione — without this, 250 mg of benzocaine
could stoichiometrically produce at most MetHb 0.5%, and this disease could not exist.</li>
<li><b>A single finite NADPH supply</b> with two consumers.</li>
<li>A <b>leucomethylene blue fork</b> with a saturable productive branch and a non-saturating futile branch.</li>
<li>An <b>allosteric term</b> linking the ferrous fraction to P50 and Hill n.</li>
<li>Two extinction coefficients per pigment and <b>one linear calibration equation</b>.</li>
</ol>

<h3>The most vulnerable assumption (not hidden)</h3>
<p>The sidebar's <b>α slider</b> is this model's Achilles' heel. If α is set to 0, the band in
tab 4 disappears, the conventional arithmetic \"MetHb 30% ≡ Hb 10.5\" becomes correct, and this model's
central claim vanishes. Quantitative human data on precisely how much methaemoglobin lowers
the P50 of the remaining ferrous subunits is scarce, and α = 0.5 is merely a value consistent with
the literature, not one uniquely determined by it.</p>
<p><b>The falsification experiment is simple:</b> in the same specimen, measure co-oximetry %MetHb and
the measured P50 together, and look at the slope. If the slope is 0, this model is wrong.</p>

<h3>Validation</h3>
<p>Because the container had no R, the mrgsolve file underwent <b>equation verification rather than
compilation verification</b>. Every ODE and parameter has been independently reimplemented in Python/scipy in
<code>mhb_reference_check.py</code>, 139 parameters were cross-checked with 0 mismatches on both sides,
and all 42 anchors pass:</p>
<pre>python3 mhb_reference_check.py --params mhb_mrgsolve_model_en.R
python3 mhb_reference_check.py --all</pre>

<h3>Warning</h3>
<p>This is an educational/research model. Do not use it for clinical decision-making. In particular, for methylene blue's
dose, cumulative ceiling, and whether to use it in G6PD deficiency, always follow current clinical guidelines
and toxicology consultation.</p>
</div>"))
}

shinyApp(ui, server)
