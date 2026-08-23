## =====================================================================
##  CONTROLLED OVARIAN STIMULATION — Shiny dashboard
##  Interactive dashboard for the controlled ovarian stimulation QSP model
##
##  Run:  R -e "shiny::runApp('cos_shiny_app.R', port = 8080)"
##  Needs: shiny, mrgsolve, dplyr, tidyr, ggplot2, DT
##         and cos_mrgsolve_model.R in the same directory.
##
##  The app is organised around the model's claim, not around its outputs:
##  tab 5 ("Trigger") is the point of the whole thing — it shows the SAME
##  ligand pool being read by three kernels, and lets the user watch the
##  VEGF area collapse by an order of magnitude while the maturation edge
##  stays intact.  Tab 7 shows what that does to the patient.
## =====================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

source("cos_mrgsolve_model.R")     # defines mod, cos_protocol, run_cos, ...

THEME <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(colour = "grey35", size = 10))

PAL <- c("#2f6fb5", "#a52a2a", "#3d7a2f", "#b5761f", "#6b3fa0",
         "#22707f", "#95552a", "#96231f")

## ---------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Controlled Ovarian Stimulation — QSP Dashboard"),
  p(style = "color:#555;",
    strong("The claim of this model: "),
    "yield and OHSS are two readings of the same state variable (the growing granulosa mass), so ",
    "they cannot be separated by the FSH dose. The separation comes only from the ",
    strong("shape"), " of the trigger — maturation needs only the ", strong("leading edge"),
    " of the signal, while OHSS is proportional to its ", strong("area"), "."),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Patient"),
      sliderInput("afc", "AFC — antral follicle count (2–9 mm)", 3, 36, 12, step = 1),
      sliderInput("t50", "Median FSH threshold T50 (IU/L)", 6, 14, 9, step = 0.5),
      sliderInput("tone", "GnRH pulse tone (PCOS 1.7–1.8)", 0.8, 2.0, 1.0,
                  step = 0.1),
      sliderInput("age", "Age (years)", 24, 45, 32, step = 1),
      hr(),
      h4("Stimulation"),
      selectInput("gonad", "Gonadotrophin",
                  c("Daily rFSH" = "rfsh",
                    "corifollitropin alfa 150 µg single dose" = "cori",
                    "hp-hMG (with LH activity)" = "hmg")),
      sliderInput("dose", "Daily dose (IU)", 50, 450, 150, step = 25),
      sliderInput("fshstart", "Stimulation start (cycle day)", 2, 5, 2, step = 1),
      checkboxInput("coast", "coasting (3 days after stopping FSH)", FALSE),
      conditionalPanel("input.coast == true",
        sliderInput("coastday", "Day FSH is stopped (cycle day)", 6, 12, 8, step = 1)),
      hr(),
      h4("Surge prevention"),
      selectInput("supp", "Method",
                  c("GnRH antagonist" = "ant",
                    "Progestin priming (PPOS)" = "ppos",
                    "None" = "none")),
      conditionalPanel("input.supp == 'ant'",
        sliderInput("antstart", "Antagonist start (cycle day)", 1, 9, 6, step = 1)),
      hr(),
      h4("Trigger"),
      selectInput("trig", "Type",
                  c("hCG 10 000 IU" = "hcg10", "hCG 5 000 IU" = "hcg5",
                    "hCG 2 500 IU" = "hcg25", "hCG 1 500 IU" = "hcg15",
                    "GnRH agonist 0.2 mg" = "ago",
                    "Dual trigger" = "dual")),
      sliderInput("opu", "Trigger→retrieval interval (hours)", 26, 46, 36, step = 1),
      hr(),
      h4("Transfer & adjuvants"),
      radioButtons("transfer", "Transfer",
                   c("Fresh transfer" = "fresh",
                     "Freeze-all" = "freeze"), inline = TRUE),
      sliderInput("p4sup", "Luteal progesterone (mg/day)", 0, 800, 600,
                  step = 100),
      checkboxInput("cab", "Cabergoline 0.5 mg/day", FALSE),
      checkboxInput("letro", "Luteal letrozole 5 mg/day", FALSE),
      sliderInput("iv", "Supportive fluid (L/day)", 0, 3, 0, step = 0.5),
      hr(),
      actionButton("go", "Run simulation", class = "btn-primary",
                   style = "width:100%")
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        ## ---------------------------------------------------------- 1
        tabPanel("1. Patient profile",
          fluidRow(
            column(6, plotOutput("p_thresh", height = 320)),
            column(6, plotOutput("p_reserve", height = 320))),
          hr(),
          h4("What the model computed for this patient (outputs, not inputs)"),
          DTOutput("t_patient"),
          helpText("Neither the follicle count nor the days of stimulation nor the number of oocytes ",
                   "retrieved is entered. Only AFC and the threshold distribution are given; the rest is integrated.")),
        ## ---------------------------------------------------------- 2
        tabPanel("2. Pharmacokinetics (PK)",
          plotOutput("p_pk", height = 620),
          helpText("rFSH: F 0.80 · V/F 26 L · multiple-dose t½ 42 h → 11.9 IU/L in blood ",
                   "at 150 IU/day. Antagonist: Cmax 11.4 ng/mL, t½ 13 h. ",
                   "hCG 10 000 IU: peak 220 IU/L, t½ 33 h. ",
                   "Agonist (triptorelin): t½ 3 h.")),
        ## ---------------------------------------------------------- 3
        tabPanel("3. Follicular growth (folliculometry)",
          plotOutput("p_foll", height = 400),
          fluidRow(
            column(6, plotOutput("p_gran", height = 300)),
            column(6, plotOutput("p_count", height = 300))),
          helpText("Each line is one threshold-quantile slot (= AFC/10 follicles). ",
                   "Low-threshold slots grow first, and a slot whose FSH falls below its ",
                   "threshold undergoes atresia.")),
        ## ---------------------------------------------------------- 4
        tabPanel("4. Endocrine",
          plotOutput("p_endo", height = 640),
          helpText("In the natural cycle E2 and inhibin pull FSH down from 7.8 to 2.6 IU/L ",
                   "so that only one follicle ovulates. Exogenous FSH removes that fall, ",
                   "so many follicles grow — the same equations.")),
        ## ---------------------------------------------------------- 5
        tabPanel("5. ★ Trigger: one ligand, three kernels",
          fluidRow(
            column(7, plotOutput("p_kern", height = 380)),
            column(5, plotOutput("p_area", height = 380))),
          hr(),
          DTOutput("t_trig"),
          helpText(HTML(paste0(
            "<b>Support kernel</b> K = 2 IU/L (basal pulses suffice) · ",
            "<b>Maturation kernel</b> K = 25 IU/L, n = 6 (only the <i>leading edge</i> of the surge is needed, ",
            "then a 34-hour autonomous clock) · <b>VEGF kernel</b> K = 40 IU/L, n = 3 ",
            "(proportional to the <i>area</i>). hCG gives both the leading edge and the area, ",
            "whereas the agonist trigger gives only the leading edge.")))),
        ## ---------------------------------------------------------- 6
        tabPanel("6. Oocytes · embryos · clinical outcome",
          fluidRow(
            column(7, plotOutput("p_chain", height = 380)),
            column(5, plotOutput("p_clbr", height = 380))),
          hr(), DTOutput("t_end"),
          helpText("2PN = MII × 0.72 · blastocysts = 2PN × (0.55 − 0.011·(age−33)) ",
                   "· euploid fraction = 0.85/(1+exp((age−38.2)/4)) ",
                   "· CLBR = 1 − (1 − 0.52)^number of euploid blastocysts")),
        ## ---------------------------------------------------------- 7
        tabPanel("7. OHSS",
          plotOutput("p_ohss", height = 620),
          fluidRow(column(12, DTOutput("t_ohss"))),
          helpText("The Golan grade is not a parameter, it is computed from the fluid state: ",
                   "Hct ≥ 45% or ascites ≥ 1.5 L = severe, Hct ≥ 55% = critical. ",
                   "Late OHSS has pregnancy hCG as its only cause, so freeze-all abolishes it.")),
        ## ---------------------------------------------------------- 8
        tabPanel("8. Scenario comparison",
          h4("Six arms in the current patient, changing only the trigger and transfer strategy"),
          actionButton("runcmp", "Run comparison", class = "btn-success"),
          br(), br(),
          plotOutput("p_cmp", height = 380),
          DTOutput("t_cmp")),
        ## ---------------------------------------------------------- 9
        tabPanel("9. Dose sweep",
          h4("Does the FSH dose separate yield from OHSS?"),
          actionButton("runsw", "Run sweep (75–450 IU)", class = "btn-success"),
          br(), br(),
          plotOutput("p_sweep", height = 420),
          DTOutput("t_sweep"),
          helpText("Model result (AFC 12): a 3-fold dose → oocytes +15%, ascites +125%. ",
                   "The oocyte yield saturates at the cohort size but the VEGF area does ",
                   "not saturate.")),
        ## --------------------------------------------------------- 10
        tabPanel("10. Monitoring chart",
          h4("Ultrasound and hormone monitoring (clinic chart format)"),
          DTOutput("t_chart"),
          helpText("Trigger criterion: three follicles of 17 mm or more. Because this criterion ",
                   "sets the days of stimulation, 'days of stimulation' is a prediction, not an input.")),
        ## --------------------------------------------------------- 11
        tabPanel("11. Assumptions and limitations",
          h4("What this model reproduces"),
          tags$ul(
            tags$li("Single ovulation in the unstimulated cycle — with no follicle-count parameter"),
            tags$li("11–12 days of stimulation, 8.0 oocytes (AFC 12) / 18.5 (AFC 32)"),
            tags$li("Trigger-day P4 behaving as a counter of follicle number"),
            tags$li("Agonist trigger: identical maturation rate, 43-fold less ascites"),
            tags$li("The 34–38 hour retrieval window (the gap between the two clocks)"),
            tags$li("Without an antagonist, a premature LH surge on stimulation day 6–8 → cycle cancellation")),
          h4("What the model is not explicitly told"),
          tags$ul(
            tags$li("How many follicles grow, and how many oocytes come out"),
            tags$li("That hCG is 'more dangerous' than the agonist"),
            tags$li("Any OHSS severity scale or risk score"),
            tags$li("That PCOS is a high-responder group")),
          h4("Known limitations (revealed by validation)"),
          tags$ul(
            tags$li("The natural-cycle dominant follicle stalls at 16.9 mm (observed 20–22 mm)"),
            tags$li("Luteal deficiency after an agonist trigger is under-predicted (−41%, clinically larger)"),
            tags$li("The cabergoline effect is more conservative than the meta-analysis (RR 0.38), at −17%"),
            tags$li("Luteal letrozole lowers E2 by 79% but changes ascites by less than 4% ",
                    "— a falsifiable prediction of the structure that reads E2 as a marker and VEGF as the mediator")),
          h4("Disclaimer"),
          p(style = "color:#a52a2a;",
            "This is a semi-quantitative model for teaching and research. It must not be used for ",
            "real clinical decisions, prescribing, or regulatory submission."))
      )
    )
  )
)

## ---------------------------------------------------------------------
build_protocol <- function(input) {
  trig <- switch(input$trig, hcg10 = "hcg", hcg5 = "hcg", hcg25 = "hcg",
                 hcg15 = "hcg", ago = "ago", dual = "dual")
  hd <- switch(input$trig, hcg10 = 10000, hcg5 = 5000, hcg25 = 2500,
               hcg15 = 1500, 10000)
  cos_protocol(
    fsh_dose  = if (input$gonad == "cori") input$dose else input$dose,
    fsh_start = input$fshstart,
    fsh_stop  = if (isTRUE(input$coast)) input$coastday else NA,
    hmg_lh    = if (input$gonad == "hmg") 75 else 0,
    cori      = if (input$gonad == "cori") 150 else 0,
    ant_start = if (input$supp == "ant") input$antstart else 0,
    trigger   = trig, hcg_dose = hd,
    opu_delay = input$opu / 24,
    cab = input$cab, letro = input$letro,
    luteal_p4 = input$p4sup,
    fresh = (input$transfer == "fresh"),
    ppos = (input$supp == "ppos"),
    iv_fluid = input$iv,
    name = "current")
}

server <- function(input, output, session) {

  patient <- reactive(list(AFC = input$afc, T50 = input$t50,
                           TONE = input$tone, AGE = input$age))

  sim <- eventReactive(input$go, {
    withProgress(message = "Integrating (65 ODEs)...", {
      run_cos(mod, patient(), build_protocol(input), tend = 32)
    })
  }, ignoreNULL = FALSE)

  ends <- reactive(cos_endpoints(sim()))

  ## ---- 1 patient ----------------------------------------------------
  output$p_thresh <- renderPlot({
    z <- c(-1.6449, -1.0364, -0.6745, -0.3853, -0.1257,
           0.1257, 0.3853, 0.6745, 1.0364, 1.6449)
    th <- input$t50 * exp(0.45 * z)
    amh <- 0.0838 * (0.55 * input$afc + input$afc * (1 - (5/8)^4/(1 + (5/8)^4))) / 0.6
    th <- th * (1 + 0.15 * (amh / 2.5 - 1))
    d <- data.frame(slot = factor(1:10), th = th, m = input$afc / 10)
    ggplot(d, aes(slot, th)) +
      geom_col(aes(fill = th), width = 0.75) +
      geom_hline(yintercept = 11.9 + 0.3, linetype = 2, colour = "#a52a2a") +
      annotate("text", x = 8.5, y = 12.6, label = "Total FSH at 150 IU/day",
               colour = "#a52a2a", size = 3.4) +
      scale_fill_gradient(low = "#c8e6b8", high = "#2f6fb5", guide = "none") +
      labs(title = "FSH threshold distribution (10 slots)",
           subtitle = sprintf("%.1f follicles per slot · AMH %.2f ng/mL (computed)",
                              input$afc / 10, amh),
           x = "Slot (sensitive → insensitive)", y = "Threshold (IU/L)") + THEME
  })

  output$p_reserve <- renderPlot({
    o <- sim()
    d <- o %>% select(time, AMH, INHB, NF11, NF17) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      scale_colour_manual(values = PAL, guide = "none") +
      labs(title = "Ovarian reserve markers and follicle counts",
           subtitle = "AMH falls during stimulation (a smaller small-follicle pool)",
           x = "Cycle day (d)", y = NULL) + THEME
  })

  output$t_patient <- renderDT({
    e <- ends()
    datatable(data.frame(
      Item = c("Days of stimulation (predicted)", "Trigger-day E2 (pg/mL)",
               "Trigger-day P4 (ng/mL)", "Follicles 11 mm or more", "Follicles 17 mm or more",
               "Oocytes retrieved", "MII oocytes", "Maturation rate", "Blastocysts", "Euploid",
               "Cumulative live birth rate (%)", "OHSS grade (0–4)"),
      Value = round(c(e$stim_days, e$E2_trig, e$P4_trig, e$foll11, e$foll17,
                   e$oocytes, e$MII, e$MII_frac, e$blast, e$euploid,
                   100 * e$CLBR, e$OHSS), 2)),
      options = list(dom = "t", pageLength = 12), rownames = FALSE)
  })

  ## ---- 2 PK ---------------------------------------------------------
  output$p_pk <- renderPlot({
    o <- sim()
    d <- o %>% transmute(time,
        `Exogenous FSH (IU/L)` = FSHTOTAL,
        `Antagonist (ng/mL)` = CANTOUT,
        `hCG (IU/L)` = CHCGOUT,
        `Agonist (ng/mL)` = CAGOOUT,
        `LH (IU/L)` = LH,
        `LH-equivalent ligand (IU/L)` = LHEQOUT) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      geom_vline(xintercept = o$ttrig[1], linetype = 3) +
      scale_colour_manual(values = PAL, guide = "none") +
      labs(title = "Pharmacokinetics and the ligand pool",
           subtitle = "Dotted = time of trigger", x = "Cycle day (d)", y = NULL) + THEME
  })

  ## ---- 3 follicles --------------------------------------------------
  output$p_foll <- renderPlot({
    o <- sim()
    d <- o %>% select(time, D1:D10) %>% pivot_longer(-time)
    d$name <- factor(d$name, levels = paste0("D", 1:10))
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 0.8) +
      geom_hline(yintercept = c(11, 17), linetype = 2, colour = "grey45") +
      geom_vline(xintercept = c(o$ttrig[1], o$topu[1]), linetype = 3) +
      scale_colour_viridis_d(option = "D", end = 0.92, name = "Slot") +
      labs(title = "Follicle diameter (corresponding to ultrasound measurement)",
           subtitle = "Dashed horizontal lines = 11 mm (retrievable) / 17 mm (trigger criterion)",
           x = "Cycle day (d)", y = "Diameter (mm)") + THEME
  })

  output$p_gran <- renderPlot({
    o <- sim()
    d <- o %>% select(time, G1:G10) %>% pivot_longer(-time)
    ggplot(d, aes(time, value, group = name, colour = name)) +
      geom_line(linewidth = 0.7, show.legend = FALSE) +
      scale_colour_viridis_d(option = "D", end = 0.92) +
      labs(title = "Granulosa mass G_i (mature follicle = 1)",
           subtitle = "This mass makes E2 · P4 · oocytes · VEGF all at once",
           x = "Cycle day (d)", y = NULL) + THEME
  })

  output$p_count <- renderPlot({
    o <- sim()
    d <- o %>% transmute(time, `MG (total granulosa mass)` = MGOUT,
                         `Ovarian volume (mL)` = OVOL) %>% pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = PAL, guide = "none") +
      labs(title = "MG — the central state variable of the model", x = "Cycle day (d)", y = NULL) +
      THEME
  })

  ## ---- 4 endocrine --------------------------------------------------
  output$p_endo <- renderPlot({
    o <- sim()
    d <- o %>% transmute(time,
        `E2 (pg/mL)` = E2, `P4 (ng/mL)` = P4, `LH (IU/L)` = LH,
        `Endogenous FSH (IU/L)` = FSHE, `Inhibin B (pg/mL)` = INHB,
        `AMH (ng/mL)` = AMH, `Luteal mass` = CL, `Surge readiness RS` = RS) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      geom_vline(xintercept = o$ttrig[1], linetype = 3) +
      scale_colour_manual(values = rep(PAL, 2), guide = "none") +
      labs(title = "Endocrine profile", x = "Cycle day (d)", y = NULL) + THEME
  })

  ## ---- 5 the three kernels -----------------------------------------
  kern <- reactive({
    o <- sim()
    o %>% transmute(time, LHEQ = LHEQOUT,
      `Support (K 2)` = LHEQOUT / (LHEQOUT + 2),
      `Maturation (K 25, n 6)` = (LHEQOUT / 25)^6 / (1 + (LHEQOUT / 25)^6),
      `VEGF (K 40, n 3)` = (LHEQOUT / 40)^3 / (1 + (LHEQOUT / 40)^3))
  })

  output$p_kern <- renderPlot({
    o <- sim(); d <- kern() %>% select(-LHEQ) %>% pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 1.0) +
      geom_vline(xintercept = c(o$ttrig[1], o$topu[1]), linetype = 3) +
      coord_cartesian(xlim = c(o$ttrig[1] - 1, o$ttrig[1] + 9)) +
      scale_colour_manual(values = c("#2f6fb5", "#a52a2a", "#96231f"),
                          name = NULL) +
      labs(title = "One ligand, three kernel occupancies",
           subtitle = "Maturation takes only the leading edge, VEGF takes the area",
           x = "Cycle day (d)", y = "Occupancy (0–1)") + THEME
  })

  output$p_area <- renderPlot({
    o <- sim(); tr <- o$ttrig[1]
    a <- data.frame(
      kernel = c("Maturation area (d)", "VEGF area (d)"),
      value = c(tail(o$EXPLH, 1) - approx(o$time, o$EXPLH, tr)$y,
                tail(o$EXPV, 1) - approx(o$time, o$EXPV, tr)$y))
    ggplot(a, aes(kernel, value, fill = kernel)) +
      geom_col(width = 0.55, show.legend = FALSE) +
      geom_text(aes(label = sprintf("%.2f d", value)), vjust = -0.4) +
      scale_fill_manual(values = c("#a52a2a", "#96231f")) +
      labs(title = "Cumulative exposure after the trigger",
           subtitle = "hCG 10 000 IU: 6.2 / 5.2 · agonist 0.2 mg: 0.5 / 0.4",
           x = NULL, y = "Integrated area (days)") + THEME
  })

  output$t_trig <- renderDT({
    o <- sim(); tr <- o$ttrig[1]; e <- ends()
    post <- o$time >= tr
    datatable(data.frame(
      Item = c("Peak ligand (IU/L)", "Maturation kernel area (d)",
               "VEGF kernel area (d)", "Maturation rate (MII/oocyte)",
               "Peak Hct (%)", "Peak ascites (L)", "Luteal day 7 P4 (ng/mL)"),
      Value = round(c(max(o$LHEQOUT[post]), e$AUC_meiosis, e$AUC_vegf,
                   e$MII_frac, e$HCT_max, e$ASC_max, e$P4_lut7), 2)),
      options = list(dom = "t"), rownames = FALSE)
  })

  ## ---- 6 embryology -------------------------------------------------
  output$p_chain <- renderPlot({
    e <- ends()
    d <- data.frame(
      step = factor(c("Oocytes retrieved", "MII", "2PN", "Blastocyst", "Euploid"),
                    levels = c("Oocytes retrieved", "MII", "2PN", "Blastocyst",
                               "Euploid")),
      n = c(e$oocytes, e$MII, e$MII * 0.72, e$blast, e$euploid))
    ggplot(d, aes(step, n, fill = step)) +
      geom_col(width = 0.68, show.legend = FALSE) +
      geom_text(aes(label = sprintf("%.1f", n)), vjust = -0.4) +
      scale_fill_manual(values = PAL) +
      labs(title = "The oocyte → embryo chain (every step is a multiplication)",
           subtitle = sprintf("Age %d · euploid fraction %.2f",
                              input$age, 0.85/(1+exp((input$age-38.2)/4))),
           x = NULL, y = "Number") + THEME
  })

  output$p_clbr <- renderPlot({
    e <- ends(); n <- seq(0, 8, 0.1)
    d <- data.frame(n = n, clbr = 100 * (1 - 0.48^n))
    ggplot(d, aes(n, clbr)) + geom_line(linewidth = 1) +
      geom_point(data = data.frame(n = e$euploid, clbr = 100 * e$CLBR),
                 colour = "#a52a2a", size = 3.5) +
      labs(title = "Saturation of the cumulative live birth rate",
           subtitle = sprintf("This patient: euploid %.2f → CLBR %.1f%%",
                              e$euploid, 100 * e$CLBR),
           x = "Number of euploid blastocysts", y = "CLBR (%)") + THEME
  })

  output$t_end <- renderDT({
    e <- ends()
    datatable(as.data.frame(lapply(e, function(x)
      if (is.numeric(x)) round(x, 2) else x)),
      options = list(dom = "t", scrollX = TRUE), rownames = FALSE)
  })

  ## ---- 7 OHSS -------------------------------------------------------
  output$p_ohss <- renderPlot({
    o <- sim()
    d <- o %>% transmute(time,
        `VEGF signal` = VEGF, `Vascular permeability PERM` = PERM,
        `Ascites (L)` = ASC, `Plasma volume (L)` = VP,
        `Haematocrit (%)` = HCT, `Pregnancy hCG (IU/L)` = PHCG,
        `Ovarian volume (mL)` = OVOL, `Golan grade` = OHSSG) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      geom_vline(xintercept = c(o$ttrig[1], o$topu[1]), linetype = 3) +
      scale_colour_manual(values = rep(c("#96231f", "#a52a2a", "#b5761f",
                                         "#2f6fb5"), 2), guide = "none") +
      labs(title = "OHSS — reading the 'area' of the same mass",
           subtitle = "Early OHSS is exogenous hCG, late OHSS is pregnancy hCG",
           x = "Cycle day (d)", y = NULL) + THEME
  })

  output$t_ohss <- renderDT({
    o <- sim(); e <- ends()
    g <- c("None", "Mild", "Moderate", "Severe", "Critical")[e$OHSS + 1]
    datatable(data.frame(
      Item = c("Golan grade", "Peak Hct (%)", "Peak ascites (L)",
               "Minimum plasma volume (L)", "Peak ovarian volume (mL)",
               "VEGF kernel area (d)", "Day of peak early OHSS (days after retrieval)",
               "Late OHSS occurs"),
      Value = c(g, round(e$HCT_max, 1), round(e$ASC_max, 2),
             round(min(o$VP), 2), round(max(o$OVOL), 0),
             round(e$AUC_vegf, 2),
             round(o$time[which.max(o$ASC)] - o$topu[1], 1),
             ifelse(max(o$PHCG) > 5, "Yes (fresh transfer + pregnancy)", "No"))),
      options = list(dom = "t"), rownames = FALSE)
  })

  ## ---- 8 scenario comparison ---------------------------------------
  cmp <- eventReactive(input$runcmp, {
    base <- build_protocol(input)
    arms <- list(
      "hCG 10000 + fresh" = modifyList(base, list(trigger = "hcg",
          hcg_dose = 10000, fresh = TRUE, name = "hCG10-fresh")),
      "hCG 10000 + freeze-all" = modifyList(base, list(trigger = "hcg",
          hcg_dose = 10000, fresh = FALSE, name = "hCG10-freeze")),
      "hCG 5000" = modifyList(base, list(trigger = "hcg", hcg_dose = 5000,
          fresh = FALSE, name = "hCG5")),
      "Agonist trigger + freeze-all" = modifyList(base, list(trigger = "ago",
          fresh = FALSE, name = "agonist")),
      "Dual trigger" = modifyList(base, list(trigger = "dual", fresh = FALSE,
          name = "dual")),
      "hCG + cabergoline" = modifyList(base, list(trigger = "hcg",
          cab = TRUE, fresh = FALSE, name = "hCG-cab")))
    withProgress(message = "Integrating 6 arms...", {
      bind_rows(lapply(names(arms), function(k) {
        e <- cos_endpoints(run_cos(mod, patient(), arms[[k]], tend = 32))
        e$arm <- k; e }))
    })
  })

  output$p_cmp <- renderPlot({
    d <- cmp() %>% select(arm, oocytes, MII, ASC_max, HCT_max, AUC_vegf,
                          CLBR) %>%
      mutate(CLBR = 100 * CLBR) %>% pivot_longer(-arm)
    ggplot(d, aes(reorder(arm, value), value, fill = name)) +
      geom_col(show.legend = FALSE) + coord_flip() +
      facet_wrap(~name, scales = "free_x", ncol = 3) +
      scale_fill_manual(values = PAL) +
      labs(title = "Trigger and transfer strategies compared (same patient, same stimulation)",
           x = NULL, y = NULL) + THEME
  })

  output$t_cmp <- renderDT({
    datatable(cmp() %>% select(arm, stim_days, E2_trig, P4_trig, oocytes, MII,
                               MII_frac, AUC_vegf, HCT_max, ASC_max, OHSS,
                               euploid, CLBR) %>%
              mutate(across(where(is.numeric), ~round(.x, 2))),
              options = list(dom = "t", scrollX = TRUE), rownames = FALSE)
  })

  ## ---- 9 dose sweep -------------------------------------------------
  sw <- eventReactive(input$runsw, {
    base <- build_protocol(input)
    withProgress(message = "Integrating the dose sweep...", {
      bind_rows(lapply(c(75, 112, 150, 225, 300, 450), function(dd) {
        e <- cos_endpoints(run_cos(mod, patient(),
               modifyList(base, list(fsh_dose = dd,
                                     name = paste0(dd, " IU"))), tend = 32))
        e$dose <- dd; e }))
    })
  })

  output$p_sweep <- renderPlot({
    d <- sw()
    s <- max(d$oocytes) / max(d$ASC_max + 1e-6)
    ggplot(d, aes(dose)) +
      geom_line(aes(y = oocytes, colour = "Oocytes retrieved"), linewidth = 1.1) +
      geom_point(aes(y = oocytes, colour = "Oocytes retrieved"), size = 2) +
      geom_line(aes(y = ASC_max * s, colour = "Peak ascites (L)"),
                linewidth = 1.1) +
      geom_point(aes(y = ASC_max * s, colour = "Peak ascites (L)"), size = 2) +
      scale_y_continuous(name = "Oocytes retrieved (n)",
        sec.axis = sec_axis(~./s, name = "Peak ascites (L)")) +
      scale_colour_manual(values = c("Oocytes retrieved" = "#3d7a2f",
                                     "Peak ascites (L)" = "#96231f"),
                          name = NULL) +
      labs(title = "The FSH dose does not separate yield from risk",
           subtitle = "Oocytes saturate at the cohort size, ascites does not",
           x = "Daily rFSH dose (IU)") + THEME
  })

  output$t_sweep <- renderDT({
    datatable(sw() %>% select(dose, stim_days, E2_trig, oocytes, MII, euploid,
                              AUC_vegf, HCT_max, ASC_max, OHSS, CLBR) %>%
              mutate(across(where(is.numeric), ~round(.x, 2))),
              options = list(dom = "t", scrollX = TRUE), rownames = FALSE)
  })

  ## ---- 10 clinic chart ---------------------------------------------
  output$t_chart <- renderDT({
    o <- sim()
    d <- o %>% filter(abs(time - round(time)) < 1e-6, time >= 1) %>%
      transmute(`Cycle day` = time,
                `E2 (pg/mL)` = round(E2),
                `P4 (ng/mL)` = round(P4, 2),
                `LH (IU/L)` = round(LH, 2),
                `FSH total (IU/L)` = round(FSHTOTAL, 1),
                `≥11 mm` = round(NF11, 1),
                `≥14 mm` = round(NF14, 1),
                `≥17 mm` = round(NF17, 1),
                `Max diameter (mm)` = round(pmax(D1, D2, D3, D4, D5), 1),
                `Hct (%)` = round(HCT, 1),
                `Ascites (L)` = round(ASC, 2))
    datatable(d, options = list(pageLength = 20, scrollX = TRUE),
              rownames = FALSE)
  })
}

shinyApp(ui, server)
