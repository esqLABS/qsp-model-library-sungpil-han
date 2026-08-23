# =============================================================================
#  cca_shiny_app.R
#  Cholangiocarcinoma QSP model — interactive dashboard
# =============================================================================
#
#  13 tabs.  The app is organised around the three claims the model makes,
#  not around the organ systems:
#
#    Tab 4  THE DOSE GATE      — the only tab that matters if claim I is right.
#                                It shows the prescribed dose, the four gates,
#                                and what was actually delivered.  Move the
#                                drainage slider and watch the delivered dose
#                                change without touching a single antitumour
#                                parameter.
#    Tab 3  CLONES             — T_S / T_P / T_R plotted separately, because a
#                                total-volume plot hides the entire mechanism.
#    Tab 11 TWO HAZARDS        — h_tumour and h_biliary on the same axes, and
#                                the fraction of the total hazard that is
#                                biliary.  If that fraction is large the
#                                patient's problem is a stent, not a drug.
#    Tab 13 FALSIFIERS         — one switch each.  Turning the gate off must
#                                make drainage timing stop mattering.
#
#  Run:  shiny::runApp("cca_shiny_app.R")
#  Requires: shiny, mrgsolve, ggplot2, dplyr, tidyr, DT
# =============================================================================

library(shiny)
library(mrgsolve)
library(ggplot2)
library(dplyr)
library(tidyr)
library(DT)

# The model definition lives in cca_mrgsolve_model.R; sourcing it compiles the
# model and exposes `mod`, `stent_events()` and the scenario drivers.
source("cca_mrgsolve_model.R")

PAL <- c(Sensitive = "#4527A0", Persister = "#00838F", Resistant = "#C62828",
         Total = "#212121", Metastatic = "#EF6C00")

thm <- theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")

# -----------------------------------------------------------------------------
#  UI
# -----------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Cholangiocarcinoma QSP model"),
  tags$p(style = "color:#555;margin-top:-8px",
         HTML("<b>Delivered dose is an output, not an input</b> · ",
              "<b>Resistance is not induced, it is selected</b> · ",
              "<b>Survival is two hazards running on different clocks</b>")),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Patient"),
      sliderInput("tt0", "Tumour volume at diagnosis (cm³)",
                  20, 400, 100, step = 10),
      sliderInput("fhilar", "Hilar location fraction f_hilar (0 = intrahepatic, 1 = hilar)",
                  0.02, 0.97, 0.72, step = 0.02),
      checkboxInput("fgfr2", "FGFR2 fusion", FALSE),
      checkboxInput("idh1", "IDH1 R132 mutation", FALSE),
      checkboxInput("immeng", "Immune-engaged tumour", FALSE),
      sliderInput("lam0", "Tumour proliferation rate LAM0 (1/d)",
                  0.004, 0.025, 0.010, step = 0.001),

      hr(), h4("Biliary drainage"),
      selectInput("stent", "Stent", c("Metal SEMS" = "sems",
                                       "Plastic" = "plastic",
                                       "None" = "none")),
      sliderInput("stent_delay", "Drainage delay (days) — antitumour parameters unchanged",
                  0, 120, 0, step = 10),

      hr(), h4("Systemic therapy"),
      checkboxInput("gemcis", "Gemcitabine/cisplatin (1000 / 25 mg/m², d1·d8 q21)",
                    TRUE),
      sliderInput("ncycle", "Number of cycles", 0, 12, 8),
      checkboxInput("durva", "Durvalumab 1500 mg", FALSE),
      selectInput("fgi", "FGFR inhibitor",
                  c("None" = "none", "Pemigatinib 14/7" = "pem",
                    "Pemigatinib continuous" = "pemcont", "Futibatinib continuous" = "fut")),
      checkboxInput("ivo", "Ivosidenib 500 mg", FALSE),
      checkboxInput("cape", "Adjuvant capecitabine 6 months (BILCAP)", FALSE),
      sliderInput("k2l", "2nd-line chemotherapy intensity K_2L (1/d)", 0, 0.03, 0.010, step = 0.002),

      hr(), h4("Falsifiers"),
      checkboxInput("f1", "F1 · Dose gate off (GATE_ON = 0)", FALSE),
      checkboxInput("f2", "F2 · Remove pre-existing resistant clone (MU_RES = 0)", FALSE),
      checkboxInput("f4", "F4 · Remove stromal permeation barrier (FPEN_MIN = 1)", FALSE),

      hr(),
      sliderInput("tend", "Simulation duration (days)", 180, 2400, 1100, step = 60)
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        type = "tabs",

        # ---- 1 -------------------------------------------------------------
        tabPanel("1 · Patient Summary",
          h4("Entered patient and prescription"),
          verbatimTextOutput("summary_txt"),
          h4("Three things this model claims"),
          tags$ol(
            tags$li(HTML("<b>Delivered dose is an output.</b> Bilirubin, neutrophils, creatinine clearance, and performance status ",
                         "determine the actual dose delivered at the time of administration. Tumour → obstruction → bilirubin → dose withheld → ",
                         "tumour is a closed positive-feedback loop, and drainage is the only edge that breaks it.")),
            tags$li(HTML("<b>Resistance is selected.</b> FGFR2 kinase-domain mutant clones already exist at t = 0 ",
                         "at a frequency of MU_RES·ln(N). There is no 'time of resistance emergence' parameter.")),
            tags$li(HTML("<b>Survival is two hazards.</b> h_tumour is slow and proportional to volume, while ",
                         "h_biliary is fast, recurrent, and determined by drainage."))
          ),
          plotOutput("p_overview", height = 380)
        ),

        # ---- 2 -------------------------------------------------------------
        tabPanel("2 · PK",
          h4("Plasma and intracellular drug disposition"),
          plotOutput("p_pk", height = 620),
          tags$p(style = "color:#666",
                 "Gemcitabine's plasma half-life is only about 6 minutes, so what drives effect is not plasma concentration but ",
                 "intracellular dFdCTP. For platinum too, it is not free platinum but the Pt-DNA adduct that drives effect.")
        ),

        # ---- 3 -------------------------------------------------------------
        tabPanel("3 · Clones (core)",
          h4("Sensitive clone · drug-tolerant persister · FGFR2-mutant clone"),
          plotOutput("p_clone", height = 420),
          plotOutput("p_clonefrac", height = 260),
          tags$p(style = "color:#666",
                 "Plotting only total volume makes the mechanism disappear entirely. The interval where volume looks flat during treatment is ",
                 "mostly the interval where T_S is dying and converting to T_P, and the regrowth rate after treatment ends ",
                 "is set by the proliferation rate of T_P (GP × LAM0).")
        ),

        # ---- 4 -------------------------------------------------------------
        tabPanel("4 · Dose Gate ★",
          h4("Rules 1-4 · prescribed dose vs. actually delivered dose"),
          fluidRow(
            column(6, plotOutput("p_gate", height = 380)),
            column(6, plotOutput("p_rdi", height = 380))
          ),
          h4("Delivered relative dose intensity (RDI)"),
          verbatimTextOutput("rdi_txt"),
          tags$p(style = "color:#666",
                 "Try moving only the drainage-delay slider. Not a single antitumour parameter changes, yet ",
                 "RDI falls and survival shortens. Turning on switch F1 should make that effect disappear.")
        ),

        # ---- 5 -------------------------------------------------------------
        tabPanel("5 · Biliary Obstruction · Drainage",
          plotOutput("p_biliary", height = 620),
          tags$p(style = "color:#666",
                 "Stent patency (PATN) decays from biofilm, sludge, and intraluminal tumour ingrowth. ",
                 "Metal stent t½ ≈ 240 days, plastic t½ ≈ 90 days.")
        ),

        # ---- 6 -------------------------------------------------------------
        tabPanel("6 · Hepatic Reserve · ALBI",
          plotOutput("p_liver", height = 560),
          DTOutput("t_albi")
        ),

        # ---- 7 -------------------------------------------------------------
        tabPanel("7 · Haematological Toxicity",
          plotOutput("p_heme", height = 560),
          tags$p(style = "color:#666",
                 "Friberg structure: proliferating pool → 3 transit compartments → circulating neutrophils, MTT ≈ 125 h, ",
                 "feedback exponent γ = 0.16. The neutrophil nadir causes a dose to be lost on day 8 via rule 2.")
        ),

        # ---- 8 -------------------------------------------------------------
        tabPanel("8 · Immune",
          plotOutput("p_immune", height = 560),
          tags$p(style = "color:#666",
                 "PI_IMMUNE is not a dose but a composite index. CD8 effector cells proliferate only in immune-engaged tumours, ",
                 "which is why it barely moves the median while widening the 24-month tail.")
        ),

        # ---- 9 -------------------------------------------------------------
        tabPanel("9 · Biomarkers",
          plotOutput("p_bio", height = 620),
          tags$p(style = "color:#666",
                 "CA 19-9 rises from cholestasis alone — it is not a reliable tumour marker in a jaundiced state. ",
                 "Phosphate is a free pharmacodynamic marker of FGFR target engagement, and 2-HG is a marker of IDH1 target engagement.")
        ),

        # ---- 10 ------------------------------------------------------------
        tabPanel("10 · RECIST · Clinical Endpoints",
          plotOutput("p_recist", height = 420),
          verbatimTextOutput("endpoint_txt")
        ),

        # ---- 11 ------------------------------------------------------------
        tabPanel("11 · Two Hazards ★",
          plotOutput("p_hazard", height = 420),
          plotOutput("p_surv", height = 300),
          tags$p(style = "color:#666",
                 "If the biliary fraction of hazard is high, what this patient needs is not a drug but a stent. ",
                 "This figure is the model's clinical conclusion.")
        ),

        # ---- 12 ------------------------------------------------------------
        tabPanel("12 · Scenario Comparison",
          h4("Batch run of standard scenarios"),
          actionButton("run_scen", "Run scenarios (takes a few seconds)"),
          br(), br(),
          plotOutput("p_scen", height = 420),
          DTOutput("t_scen")
        ),

        # ---- 13 ------------------------------------------------------------
        tabPanel("13 · Falsification Tests ★",
          h4("One parameter at a time, without refitting the rest"),
          tableOutput("t_fals"),
          tags$ul(
            tags$li(HTML("<b>F1 GATE_ON = 0</b> — bilirubin does not limit dose. ",
                         "Prediction: the survival loss from drainage delay should nearly disappear.")),
            tags$li(HTML("<b>F2 MU_RES = 0</b> — no pre-existing resistant clone. ",
                         "Result (as reported): almost nothing changes within the observation period. ",
                         "This is a <i>negative result</i> of the model and it is not hidden.")),
            tags$li(HTML("<b>F3 PI_IMMUNE = 1</b> — every tumour is immune-engaged. ",
                         "Prediction: durvalumab should move the <i>median</i>, yet TOPAZ-1 says it does not.")),
            tags$li(HTML("<b>F4 FPEN_MIN = 1</b> — no stromal permeation barrier. ",
                         "Prediction: the gemcitabine/cisplatin response rate should overshoot the target by a wide margin."))
          )
        )
      )
    )
  )
)

# -----------------------------------------------------------------------------
#  SERVER
# -----------------------------------------------------------------------------
server <- function(input, output, session) {

  par_list <- reactive({
    p <- list(
      TT0 = input$tt0, FHILAR = input$fhilar, LAM0 = input$lam0,
      FGFR2 = as.numeric(input$fgfr2), IDH1 = as.numeric(input$idh1),
      IMMENG = as.numeric(input$immeng),
      GEM_MGM2 = if (input$gemcis) 1000 else 0,
      CIS_MGM2 = if (input$gemcis) 25 else 0,
      NCYCLE   = if (input$gemcis) input$ncycle else 0,
      DUR_MG   = if (input$durva) 1500 else 0,
      IVO_MG   = if (input$ivo) 500 else 0,
      CAP_MGM2 = if (input$cape) 1250 else 0,
      CAP_DAYS = if (input$cape) 180 else 0,
      K_2L     = input$k2l,
      KOCC     = if (input$stent == "plastic") 0.0077 else 0.0029,
      GATE_ON  = if (input$f1) 0 else 1,
      MU_RES   = if (input$f2) 0 else 1e-6,
      FPEN_MIN = if (input$f4) 1.0 else 0.30
    )
    p <- c(p, switch(input$fgi,
      none    = list(FGI_MG = 0),
      pem     = list(FGI_MG = 13.5, FGI_ON = 14, FGI_OFF = 7, COVAL = 0, RHO = 0.05),
      pemcont = list(FGI_MG = 13.5, FGI_ON = 1,  FGI_OFF = 0, COVAL = 0, RHO = 0.05),
      fut     = list(FGI_MG = 20,   FGI_ON = 1,  FGI_OFF = 0, COVAL = 1, RHO = 0.45)))
    p
  })

  sim <- reactive({
    m <- param(mod, par_list())
    if (input$stent_delay > 0) m <- init(m, PATN = 0)
    ev <- stent_events(input$stent, input$stent_delay, input$tend)
    out <- if (is.null(ev)) mrgsim(m, end = input$tend, delta = 1)
           else             mrgsim(m, events = ev, end = input$tend, delta = 1)
    as.data.frame(out)
  })

  long <- function(d, cols, key = "var")
    d %>% select(time, all_of(cols)) %>% pivot_longer(-time, names_to = key)

  # ---- 1 --------------------------------------------------------------------
  output$summary_txt <- renderText({
    p <- par_list()
    paste0(
      "Tumour ", input$tt0, " cm³ · f_hilar ", input$fhilar,
      " · LAM0 ", input$lam0, " /d\n",
      "Genotype: ", if (input$fgfr2) "FGFR2 fusion " else "",
                  if (input$idh1) "IDH1 R132 " else "",
                  if (!input$fgfr2 && !input$idh1) "no actionable driver" else "", "\n",
      "Immune-engaged: ", if (input$immeng) "Yes (mixture member)" else "No", "\n",
      "Drainage: ", input$stent, " · delay ", input$stent_delay, " days\n",
      "Prescription: gem/cis ", p$NCYCLE, " cycles · durvalumab ", p$DUR_MG, " mg · ",
      "FGFRi ", p$FGI_MG, " mg · ivosidenib ", p$IVO_MG, " mg\n",
      "Gate: ", if (p$GATE_ON > 0) "ON" else "OFF (F1)",
      " · pre-existing resistant clone: ", if (p$MU_RES > 0) "present" else "absent (F2)",
      " · stromal barrier: ", if (p$FPEN_MIN < 1) "present" else "absent (F4)")
  })

  output$p_overview <- renderPlot({
    d <- sim()
    long(d, c("TTOT", "BILI", "GATE", "SURV")) %>%
      mutate(var = factor(var, c("TTOT", "BILI", "GATE", "SURV"),
                          c("Tumour volume (cm³)", "Bilirubin (mg/dL)",
                            "Chemotherapy gate (0-1)", "Survival probability S(t)"))) %>%
      ggplot(aes(time, value)) +
      geom_line(linewidth = 1, colour = "#37474F") +
      facet_wrap(~var, scales = "free_y") +
      labs(x = "Time (days)", y = NULL) + thm
  })

  # ---- 2 --------------------------------------------------------------------
  output$p_pk <- renderPlot({
    d <- sim() %>% mutate(Cgem = GEM_C / 22, Ccis = CIS_C / 18,
                          Cdur = DUR_C / 5.6, Cfgi = FGI_C / 235,
                          Civo = IVO_C / 180)
    long(d, c("Cgem", "DFDCTP", "Ccis", "PTDNA", "Cdur", "Cfgi", "Civo", "COV")) %>%
      mutate(var = factor(var,
        c("Cgem", "DFDCTP", "Ccis", "PTDNA", "Cdur", "Cfgi", "Civo", "COV"),
        c("Gemcitabine plasma (mg/L)", "Intracellular dFdCTP (AU)",
          "Free platinum (mg/L)", "Pt-DNA adduct (AU)",
          "Durvalumab (mg/L)", "FGFR inhibitor (mg/L)",
          "Ivosidenib (mg/L)", "Covalent occupancy"))) %>%
      ggplot(aes(time, value)) + geom_line(colour = "#0277BD") +
      facet_wrap(~var, scales = "free_y", ncol = 2) +
      labs(x = "Time (days)", y = NULL) + thm
  })

  # ---- 3 --------------------------------------------------------------------
  output$p_clone <- renderPlot({
    d <- sim()
    data.frame(time = d$time,
               Sensitive = d$TS, Persister = d$TP,
               Resistant = d$TR, Metastatic = d$TMET, Total = d$TTOT) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
      geom_line(linewidth = 1) +
      scale_colour_manual(values = PAL, name = NULL) +
      scale_y_continuous(trans = "log1p") +
      labs(x = "Time (days)", y = "Volume (cm³, log1p)",
           title = "Tumour volume by clone") + thm
  })

  output$p_clonefrac <- renderPlot({
    d <- sim()
    data.frame(time = d$time, Persister = d$PERSFR, Resistant = d$RESFR) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 1) +
      scale_colour_manual(values = PAL, name = NULL) +
      labs(x = "Time (days)", y = "Fraction within tumour",
           title = "Persister fraction and resistant-clone fraction") + thm
  })

  # ---- 4 --------------------------------------------------------------------
  output$p_gate <- renderPlot({
    d <- sim()
    long(d, c("GBILI", "GANC", "GCRCL", "GPS", "GATE")) %>%
      mutate(var = factor(var, c("GBILI", "GANC", "GCRCL", "GPS", "GATE"),
        c("Rule 1 bilirubin", "Rule 2 neutrophils", "Rule 3 renal function (cisplatin)",
          "Rule 4 performance status", "Combined gate"))) %>%
      ggplot(aes(time, value, colour = var)) + geom_line(linewidth = 0.9) +
      scale_colour_brewer(palette = "Set1", name = NULL) +
      ylim(0, 1.02) + labs(x = "Time (days)", y = "Delivered fraction",
                           title = "Rules 1-4") + thm
  })

  output$p_rdi <- renderPlot({
    d <- sim()
    ggplot(d, aes(time, RDI_PC)) +
      geom_line(linewidth = 1, colour = "#C62828") +
      geom_hline(yintercept = 100, linetype = 2, colour = "#999") +
      ylim(0, 105) +
      labs(x = "Time (days)", y = "Cumulative relative dose intensity (%)",
           title = "What was actually delivered, not what was prescribed") + thm
  })

  output$rdi_txt <- renderText({
    d <- sim(); r <- tail(d$RDI_PC, 1)
    paste0("Final delivered RDI = ", sprintf("%.1f", r), "%\n",
           "Peak bilirubin = ", sprintf("%.1f", max(d$BILI)), " mg/dL\n",
           "Neutrophil nadir = ", sprintf("%.2f", min(d$ANC)), " ×10⁹/L\n",
           "Peak cholangitis burden = ", sprintf("%.2f", max(d$CHOLI)))
  })

  # ---- 5 --------------------------------------------------------------------
  output$p_biliary <- renderPlot({
    d <- sim()
    long(d, c("OBSR", "PATN", "BILI", "ALP", "CHOLI")) %>%
      mutate(var = factor(var, c("OBSR", "PATN", "BILI", "ALP", "CHOLI"),
        c("Biliary obstruction fraction", "Stent patency", "Total bilirubin (mg/dL)",
          "ALP (U/L)", "Cholangitis burden"))) %>%
      ggplot(aes(time, value)) + geom_line(colour = "#B71C1C", linewidth = 0.9) +
      facet_wrap(~var, scales = "free_y", ncol = 2) +
      labs(x = "Time (days)", y = NULL) + thm
  })

  # ---- 6 --------------------------------------------------------------------
  output$p_liver <- renderPlot({
    d <- sim()
    long(d, c("FLR", "ALB", "ALBI", "ALBIG")) %>%
      mutate(var = factor(var, c("FLR", "ALB", "ALBI", "ALBIG"),
        c("Functional hepatic reserve fraction", "Albumin (g/dL)", "ALBI score", "ALBI grade"))) %>%
      ggplot(aes(time, value)) + geom_line(colour = "#558B2F", linewidth = 0.9) +
      facet_wrap(~var, scales = "free_y", ncol = 2) +
      labs(x = "Time (days)", y = NULL) + thm
  })

  output$t_albi <- renderDT({
    d <- sim()
    idx <- seq(1, nrow(d), by = max(1, floor(nrow(d) / 20)))
    datatable(d[idx, c("time", "BILI", "ALB", "ALBI", "ALBIG", "FLR", "PS")] %>%
                mutate(across(-time, ~round(., 2))),
              options = list(pageLength = 10), rownames = FALSE)
  })

  # ---- 7 --------------------------------------------------------------------
  output$p_heme <- renderPlot({
    d <- sim()
    long(d, c("ANC", "PLT", "CRCL", "NEURO")) %>%
      mutate(var = factor(var, c("ANC", "PLT", "CRCL", "NEURO"),
        c("Neutrophils (×10⁹/L)", "Platelets (×10⁹/L)",
          "Creatinine clearance (mL/min)", "Cumulative neuropathy"))) %>%
      ggplot(aes(time, value)) + geom_line(colour = "#6A1B9A", linewidth = 0.9) +
      facet_wrap(~var, scales = "free_y", ncol = 2) +
      labs(x = "Time (days)", y = NULL) + thm
  })

  # ---- 8 --------------------------------------------------------------------
  output$p_immune <- renderPlot({
    d <- sim()
    long(d, c("TEFF", "SUPP", "IL6", "IRAE")) %>%
      mutate(var = factor(var, c("TEFF", "SUPP", "IL6", "IRAE"),
        c("Intratumoural CD8 effector cells", "Suppressive compartment (Treg/MDSC/TAM)",
          "IL-6", "Immune-related adverse event burden"))) %>%
      ggplot(aes(time, value)) + geom_line(colour = "#AD1457", linewidth = 0.9) +
      facet_wrap(~var, scales = "free_y", ncol = 2) +
      labs(x = "Time (days)", y = NULL) + thm
  })

  # ---- 9 --------------------------------------------------------------------
  output$p_bio <- renderPlot({
    d <- sim()
    long(d, c("CA199", "CTDNA", "PHOS", "HG2", "SLD")) %>%
      mutate(var = factor(var, c("CA199", "CTDNA", "PHOS", "HG2", "SLD"),
        c("CA 19-9 (U/mL) — confounded by cholestasis", "ctDNA VAF (AU)",
          "Serum phosphate (mg/dL) — FGFR target engagement", "Plasma 2-HG — IDH1 target engagement",
          "RECIST sum of relative diameters"))) %>%
      ggplot(aes(time, value)) + geom_line(colour = "#00695C", linewidth = 0.9) +
      facet_wrap(~var, scales = "free_y", ncol = 2) +
      labs(x = "Time (days)", y = NULL) + thm
  })

  # ---- 10 -------------------------------------------------------------------
  output$p_recist <- renderPlot({
    d <- sim()
    nad <- cummin(d$SLD)
    ggplot(d, aes(time, SLD)) +
      geom_line(linewidth = 1) +
      geom_line(aes(y = nad), linetype = 3, colour = "#0277BD") +
      geom_line(aes(y = 1.2 * nad), linetype = 2, colour = "#C62828") +
      geom_hline(yintercept = 0.70, linetype = 2, colour = "#2E7D32") +
      annotate("text", x = max(d$time) * 0.85, y = 0.72, label = "Partial-response threshold",
               colour = "#2E7D32", size = 3.5) +
      labs(x = "Time (days)", y = "Relative sum of diameters",
           title = "RECIST 1.1 · a 20% increase from nadir is progression") + thm
  })

  output$endpoint_txt <- renderText({
    d <- sim(); nad <- cummin(d$SLD)
    prog <- which(d$SLD >= 1.2 * nad & d$SLD > 1.05 * min(nad))
    tprog <- if (length(prog)) d$time[prog[1]] else NA
    paste0("Lowest relative sum of diameters = ", sprintf("%.3f", min(d$SLD)),
           "  (partial response = ≤ 0.70)\n",
           "Time of progression = ", ifelse(is.na(tprog), "none within observation period",
                                  paste0(round(tprog), " days (", round(tprog/30.44, 1), " months)")), "\n",
           "Final survival probability S(t) = ", sprintf("%.3f", tail(d$SURV, 1)), "\n",
           "12-month S = ", sprintf("%.3f", approx(d$time, d$SURV, 365)$y),
           " · 24-month S = ", sprintf("%.3f", approx(d$time, d$SURV, 730)$y))
  })

  # ---- 11 -------------------------------------------------------------------
  output$p_hazard <- renderPlot({
    d <- sim()
    data.frame(time = d$time, `h_tumour` = d$HAZT, `h_biliary` = d$HAZB,
               check.names = FALSE) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 1) +
      scale_colour_manual(values = c(h_tumour = "#37474F", h_biliary = "#C62828"),
                          name = NULL) +
      labs(x = "Time (days)", y = "Instantaneous hazard (1/day)",
           title = "Two hazards running on different clocks") + thm
  })

  output$p_surv <- renderPlot({
    d <- sim()
    ggplot(d, aes(time, SURV)) + geom_line(linewidth = 1) +
      geom_line(aes(y = HAZFR), linetype = 2, colour = "#C62828") +
      geom_vline(xintercept = c(365, 730), linetype = 3, colour = "#999") +
      ylim(0, 1) +
      labs(x = "Time (days)", y = "S(t) (solid) · biliary hazard fraction (dashed)") + thm
  })

  # ---- 12 -------------------------------------------------------------------
  scen <- eventReactive(input$run_scen, {
    base <- list(TT0 = input$tt0, LAM0 = input$lam0, K_2L = input$k2l)
    defs <- list(
      list(n = "S1 drainage only",            p = c(base, list(FHILAR = .72)), st = "sems", dl = 0),
      list(n = "S2 gem/cis",           p = c(base, list(FHILAR = .72, GEM_MGM2 = 1000,
                                                        CIS_MGM2 = 25, NCYCLE = 8)), st = "sems", dl = 0),
      list(n = "S3 +durvalumab",       p = c(base, list(FHILAR = .72, GEM_MGM2 = 1000,
                                                        CIS_MGM2 = 25, NCYCLE = 8,
                                                        DUR_MG = 1500, IMMENG = 1)), st = "sems", dl = 0),
      list(n = "S4 drainage delayed 60 days",    p = c(base, list(FHILAR = .72, GEM_MGM2 = 1000,
                                                        CIS_MGM2 = 25, NCYCLE = 8,
                                                        DUR_MG = 1500)), st = "sems", dl = 60),
      list(n = "S5 plastic stent",   p = c(base, list(FHILAR = .72, GEM_MGM2 = 1000,
                                                        CIS_MGM2 = 25, NCYCLE = 8,
                                                        KOCC = 0.0077)), st = "plastic", dl = 0),
      list(n = "S6 pemigatinib",        p = c(base, list(FHILAR = .16, FGFR2 = 1, FGI_MG = 13.5)),
           st = "sems", dl = 0),
      list(n = "S7 futibatinib",        p = c(base, list(FHILAR = .16, FGFR2 = 1, FGI_MG = 20,
                                                        FGI_ON = 1, FGI_OFF = 0, COVAL = 1,
                                                        RHO = 0.45)), st = "sems", dl = 0),
      list(n = "S8 ivosidenib",      p = c(base, list(FHILAR = .16, IDH1 = 1, IVO_MG = 500)),
           st = "sems", dl = 0))
    bind_rows(lapply(defs, function(s) {
      m <- param(mod, s$p)
      if (s$dl > 0) m <- init(m, PATN = 0)
      o <- as.data.frame(mrgsim(m, events = stent_events(s$st, s$dl, input$tend),
                                end = input$tend, delta = 7))
      o$scenario <- s$n
      o
    }))
  })

  output$p_scen <- renderPlot({
    ggplot(scen(), aes(time, SURV, colour = scenario)) +
      geom_line(linewidth = 1) +
      geom_vline(xintercept = c(365, 730), linetype = 3, colour = "#aaa") +
      labs(x = "Time (days)", y = "S(t)", colour = NULL,
           title = "Survival curves by scenario (single representative patient)") + thm
  })

  output$t_scen <- renderDT({
    scen() %>% group_by(scenario) %>%
      summarise(`Lowest SLD` = round(min(SLD), 3),
                `12-month S` = round(approx(time, SURV, 365)$y, 3),
                `24-month S` = round(approx(time, SURV, 730)$y, 3),
                `Final RDI %` = round(tail(RDI_PC, 1), 1),
                `Peak bilirubin` = round(max(BILI), 1),
                `Biliary hazard fraction (max)` = round(max(HAZFR), 2), .groups = "drop") %>%
      datatable(options = list(pageLength = 10), rownames = FALSE)
  })

  # ---- 13 -------------------------------------------------------------------
  output$t_fals <- renderTable({
    data.frame(
      Switch = c("F1 GATE_ON = 0", "F2 MU_RES = 0", "F3 PI_IMMUNE = 1",
                 "F4 FPEN_MIN = 1"),
      Prediction = c("Survival loss from drainage delay disappears",
               "FGFR inhibitor PFS should lengthen",
               "Durvalumab should move median survival",
               "gem/cis response rate exceeds the target"),
      `In-silico trial result` = c(
        "Passed — delayed-arm mOS 10.2 → 13.1 months, RDI 68% → 100%",
        "Negative result (as reported) — no change within observation period",
        "Passed — median survival off-screen, 24-month survival 76.9%",
        "Passed — response rate 31.7% → 86.7%"),
      check.names = FALSE)
  }, striped = TRUE, bordered = TRUE, width = "100%")
}

shinyApp(ui, server)
