# =============================================================================
#  aad_shiny_app.R
#  Acute Type B Aortic Dissection — QSP simulator dashboard
#  Acute type B aortic dissection — QSP simulator dashboard
# =============================================================================
#
#  The app exists to make one comparison hard to avoid: set any two arms to the
#  SAME systolic pressure and look at what else differs.  Tab 3 does nothing but
#  that.  Every other tab is there to show which mechanism carried the
#  difference.
#
#  Tabs
#    1  Patient & lesion profile
#    2  Pharmacokinetics & receptor occupancy
#    3  Pressure versus impulse — THE dissociation
#    4  The two-lumen divider
#    5  Wall stress & rupture hazard
#    6  False-lumen thrombus & the paradox
#    7  Branch perfusion & the renin trap
#    8  Spinal cord & the collateral network
#    9  Scenario comparison
#   10  Clinical endpoints & biomarkers
#
#  Run with:
#    shiny::runApp("aortic-dissection/aad_shiny_app.R")
# =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)

# -----------------------------------------------------------------------------
#  model
# -----------------------------------------------------------------------------
MODEL_DIR <- if (dir.exists("aortic-dissection")) "aortic-dissection" else "."
mod <- mread_cache("aad_mrgsolve_model", MODEL_DIR)

THEME <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")

# The 22 arms of the reference model, so the app and the paper agree.
SCENARIOS <- list(
  "S01 Untreated"                                                = list(),
  "S02 Nicardipine alone (15 mg/h)"                              = list(VD_AGENT = 1, VD_MAX = 15),
  "S03 Esmolol alone (titrated to blood pressure)"               = list(BB_AGENT = 1, BB_MAX = 1350, BB_ON_SBP = 1),
  "S04 Guideline sequence: esmolol → nicardipine"                = list(BB_AGENT = 1, BB_MAX = 900, VD_AGENT = 1, VD_MAX = 15),
  "S05 Nitroprusside alone"                                      = list(VD_AGENT = 3, VD_MAX = 13.5),
  "S06 Nitroprusside + esmolol (Wheat)"                          = list(VD_AGENT = 3, VD_MAX = 13.5, BB_AGENT = 1, BB_MAX = 900),
  "S07 Labetalol alone"                                          = list(BB_AGENT = 2, BB_MAX = 120, BB_ON_SBP = 1),
  "S08 Clevidipine + esmolol"                                    = list(VD_AGENT = 2, VD_MAX = 21, BB_AGENT = 1, BB_MAX = 900),
  "S09 Analgesia only (fentanyl 75 µg/h)"                        = list(R_FEN = 0.075),
  "S10 Full protocol (analgesia + β + CCB)"                      = list(R_FEN = 0.075, BB_AGENT = 1, BB_MAX = 900,
                                                        VD_AGENT = 1, VD_MAX = 15),
  "S21 Nicardipine over-dose (pressure-matched)"                 = list(VD_AGENT = 1, VD_MAX = 30),
  "S11 Aggressive target SBP 95 + renal malperfusion"            = list(VD_AGENT = 1, VD_MAX = 32, BB_AGENT = 1, BB_MAX = 900,
                                                        SBP_SET = 95, OBST_REN = 0.45, INV_REN = 1, INV_MES = 1),
  "S12 Renal malperfusion, standard target"                      = list(VD_AGENT = 1, VD_MAX = 32, BB_AGENT = 1, BB_MAX = 900,
                                                        OBST_REN = 0.45, INV_REN = 1, INV_MES = 1),
  "S13 High-dose SNP + renal impairment (cyanide)"               = list(VD_AGENT = 3, VD_MAX = 45, BB_AGENT = 1, BB_MAX = 900,
                                                        OBST_REN = 0.45, INV_REN = 1, SBP_SET = 100),
  "S14 Chronic 5 years: patent false lumen"                      = list(GRE0 = 2.40, D_MET = 100, D_LOS = 50, .end = 43800),
  "S15 Chronic 5 years: partial thrombosis (blind sac)"          = list(GRE0 = 0.05, D_MET = 100, D_LOS = 50, .end = 43800),
  "S16 Chronic 5 years: complete thrombosis (small entry tear)"  = list(GRE0 = 0.05, AEN0 = 9, D_MET = 100, D_LOS = 50, .end = 43800),
  "S17 Chronic 5 years: TEVAR (day 14)"                          = list(GRE0 = 0.05, TEVAR_T = 336, D_MET = 100, D_LOS = 50, .end = 43800),
  "S22 Chronic 5 years: TEVAR + perfusion pressure augmentation" = list(GRE0 = 0.05, TEVAR_T = 336, D_MET = 50, D_LOS = 25,
                                                        SBP_SET = 135, .end = 43800),
  "S18 Marfan phenotype + losartan (5 years)"                    = list(GEN_MMP = 2.2, GRE0 = 0.05, D_MET = 100, D_LOS = 50, .end = 43800),
  "S19 Marfan phenotype, β-blocker only (5 years)"               = list(GEN_MMP = 2.2, GRE0 = 0.05, D_MET = 100, D_LOS = 0, .end = 43800),
  "S20 Chronic 5 years: adherence 25%"                           = list(GRE0 = 0.05, D_MET = 100, D_LOS = 50, ADHERE = 0.25, .end = 43800)
)

simulate_arm <- function(pars, end = NULL) {
  e <- if (!is.null(pars$.end)) pars$.end else (end %||% 72)
  pars$.end <- NULL
  d <- if (e > 800) 24 else 0.1
  mod %>% param(pars) %>% mrgsim(end = e, delta = d) %>% as_tibble()
}
`%||%` <- function(a, b) if (is.null(a)) b else a

# -----------------------------------------------------------------------------
#  UI
# -----------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Acute Type B Aortic Dissection — QSP Simulator"),
  p(tags$b("One product:"),
    "ASI = (false-lumen systolic stress / reference) × (dP/dt max / reference).",
    "Pressure enters only the first term, sympathetic drive only the second.",
    "A vasodilator lowers the first term and raises the second."),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Patient and lesion"),
      sliderInput("AEN0", "Entry tear area A_en (mm²)", 4, 120, 42, step = 1),
      sliderInput("GRE0", "Re-entry conductance G_re (0 = blind sac)", 0, 3, 1, step = 0.05),
      sliderInput("DFL0", "False-lumen diameter (mm)", 8, 40, 17, step = 1),
      sliderInput("DTL0", "True-lumen diameter (mm)", 10, 32, 20, step = 1),
      sliderInput("ELN", "Elastic lamellar integrity (0-1)", 0.3, 1, 0.82, step = 0.01),
      sliderInput("GEN_MMP", "Genetic MMP-9 multiplier (Marfan 2.2)", 1, 3, 1, step = 0.1),
      sliderInput("BARO_REF", "Baroreceptor set point (mmHg)", 85, 125, 106, step = 1),
      checkboxInput("INV_REN", "Renal artery arises from the dissected segment", FALSE),
      sliderInput("OBST_REN", "Static renal artery obstruction (0-1)", 0, 0.8, 0, step = 0.05),

      hr(), h4("Acute-phase drugs"),
      selectInput("BB", "β-blocker",
                  c("None" = "0", "Esmolol" = "1", "Labetalol" = "2"), "1"),
      sliderInput("BB_MAX", "β-blocker maximum infusion (mg/h)", 0, 1350, 900, step = 25),
      checkboxInput("BB_ON_SBP", "Titrate the β-blocker to blood pressure (heart rate by default)", FALSE),
      selectInput("VD", "Vasodilator",
                  c("None" = "0", "Nicardipine" = "1", "Clevidipine" = "2",
                    "Nitroprusside" = "3"), "1"),
      sliderInput("VD_MAX", "Vasodilator maximum infusion (mg/h)", 0, 45, 15, step = 1),
      sliderInput("R_FEN", "Fentanyl (mg/h)", 0, 0.2, 0, step = 0.005),
      sliderInput("SBP_SET", "Systolic target (mmHg)", 85, 140, 110, step = 1),

      hr(), h4("Chronic phase and intervention"),
      sliderInput("D_MET", "Metoprolol (mg / 12 h)", 0, 200, 0, step = 25),
      sliderInput("D_LOS", "Losartan (mg / 12 h)", 0, 100, 0, step = 25),
      sliderInput("ADHERE", "Medication adherence", 0, 1, 1, step = 0.05),
      numericInput("TEVAR_T", "TEVAR time (h, −1 = not performed)", -1, min = -1, max = 8760),
      sliderInput("TEVAR_COVER", "Segmental artery coverage fraction", 0, 0.7, 0.30, step = 0.05),
      sliderInput("END", "Simulation duration (h)", 24, 43800, 72, step = 24),
      actionButton("go", "Run simulation", class = "btn-primary")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",

        # ---- 1 -------------------------------------------------------------
        tabPanel(
          "1 · Patient & lesion profile",
          h4("The lesion geometry, and what that geometry has already settled"),
          tableOutput("profile"),
          p(tags$b("How to read it:"), "the outer wall of the false lumen is outer media + adventitia only, so its thickness is 40% of the whole wall",
            "(1.90 mm → 0.76 mm). At the same pressure the circumferential stress in the false lumen is about 2.5 times that in the true lumen,",
            "and this asymmetry is settled by anatomy, not by drugs."),
          plotOutput("p_geom", height = "330px"),
          plotOutput("p_stressbar", height = "300px")
        ),

        # ---- 2 -------------------------------------------------------------
        tabPanel(
          "2 · Pharmacokinetics & receptor occupancy",
          fluidRow(column(6, plotOutput("p_pk", height = "340px")),
                   column(6, plotOutput("p_occ", height = "340px"))),
          plotOutput("p_rate", height = "290px"),
          p(tags$b("Titration is an integral controller."), "The infusion rate is set by integrating the error against the target,",
            "and saturation shows up in the table just as it is — nicardipine, at its labelled ceiling of 15 mg/h,",
            "cannot bring the systolic pressure below 121 mmHg."),
          tableOutput("t_exposure")
        ),

        # ---- 3 -------------------------------------------------------------
        tabPanel(
          "3 · Pressure vs impulse (the core)",
          h4("What differs at the same pressure"),
          fluidRow(column(4, selectInput("cmpA", "Comparison A", names(SCENARIOS),
                                         "S21 Nicardipine over-dose (pressure-matched)")),
                   column(4, selectInput("cmpB", "Comparison B", names(SCENARIOS),
                                         "S04 Guideline sequence: esmolol → nicardipine")),
                   column(4, selectInput("cmpC", "Reference C", names(SCENARIOS),
                                         "S01 Untreated"))),
          plotOutput("p_dissoc", height = "430px"),
          tableOutput("t_dissoc"),
          p(tags$b("What the model says:"), "nicardipine at its labelled ceiling lowers the systolic pressure by 33 mmHg and",
            "changes dP/dt max by 0.3%. At the same pressure the flap stress-impulse differs 2.8-fold,",
            "the 72-hour tear growth 14-fold, and the 72-hour mortality 9-fold.",
            "Mortality in the vasodilator-alone arm is higher than in the untreated arm.")
        ),

        # ---- 4 -------------------------------------------------------------
        tabPanel(
          "4 · The two-lumen divider",
          h4("Two conductances make three phenotypes"),
          plotOutput("p_divider", height = "380px"),
          plotOutput("p_gre_sweep", height = "330px"),
          p(tags$b("The arithmetic of the blind sac:"), "with G_re → 0 there is no outflow, so there is no pressure drop either",
            "(false-lumen mean pressure = systemic mean pressure), and the pulse pressure is low-pass filtered away by τ_FL, so",
            "the false-lumen diastolic pressure rises to the systemic mean pressure — higher than the true-lumen diastolic pressure.",
            "While the distal stagnant zone clots and closes the re-entry tear, the entry tear stays open as it was.")
        ),

        # ---- 5 -------------------------------------------------------------
        tabPanel(
          "5 · Wall stress and rupture hazard",
          plotOutput("p_sigma", height = "350px"),
          fluidRow(column(6, plotOutput("p_hazard", height = "320px")),
                   column(6, plotOutput("p_wall", height = "320px"))),
          p(tags$b("Rupture hazard is an exponential in σ/S:"),
            "h = 1.10e−4 · exp((σ/S − 1)/0.155) /h.",
            "calibrated to 2% per year at σ/S = 0.40 and 11% per year at 0.67.")
        ),

        # ---- 6 -------------------------------------------------------------
        tabPanel(
          "6 · False-lumen thrombus and the paradox",
          h4("Thrombus supports the wall linearly and thins it exponentially"),
          plotOutput("p_paradox", height = "360px"),
          p(tags$b("Analytic solution:"), "for growth ∝ (1 + A·THR)^n · (1 − s·THR),",
            "dg/dTHR = 0 → THR* = (nA − s)/((1+n)sA) = 0.426.",
            "The numerical optimum is 0.425. The worst false lumen is the one that has clotted only half way."),
          plotOutput("p_thr_time", height = "330px")
        ),

        # ---- 7 -------------------------------------------------------------
        tabPanel(
          "7 · Branch perfusion and the renin trap",
          plotOutput("p_perf", height = "350px"),
          plotOutput("p_raas", height = "330px"),
          p(tags$b("The trap:"), "lower the blood pressure and the kidney becomes more ischaemic, and renin restores",
            "the resistance by way of angiotensin II. On the same prescription and the same target it becomes SBP 141 versus 113, and",
            "the vasodilator is already at its ceiling. The answer is not a fourth vasodilator but opening the kidney."),
          tableOutput("t_cn")
        ),

        # ---- 8 -------------------------------------------------------------
        tabPanel(
          "8 · Spinal cord and the collateral network",
          plotOutput("p_cord", height = "350px"),
          plotOutput("p_collat", height = "320px"),
          p(tags$b("The collateral network is a clock:"), "segmental artery loss happens at once, but collateral recruitment",
            "takes about 5 days. The paraplegia risk sits inside that window, and staged coverage and",
            "perfusion pressure augmentation are procedures aimed at that window — not discharge prescriptions.")
        ),

        # ---- 9 -------------------------------------------------------------
        tabPanel(
          "9 · Scenario comparison",
          checkboxGroupInput("multi", "Scenarios to compare",
                             choices = names(SCENARIOS),
                             selected = c("S01 Untreated",
                                          "S02 Nicardipine alone (15 mg/h)",
                                          "S04 Guideline sequence: esmolol → nicardipine",
                                          "S07 Labetalol alone",
                                          "S21 Nicardipine over-dose (pressure-matched)"),
                             inline = TRUE),
          plotOutput("p_multi", height = "460px"),
          tableOutput("t_multi")
        ),

        # ---- 10 ------------------------------------------------------------
        tabPanel(
          "10 · Clinical endpoints & biomarkers",
          plotOutput("p_end", height = "380px"),
          fluidRow(column(6, plotOutput("p_bio", height = "320px")),
                   column(6, plotOutput("p_growth", height = "320px"))),
          tableOutput("t_end"),
          p(tags$b("Caution:"), "this model is for teaching and hypothesis generation and has not been validated for use in patient care.",
            "The predicted null effect of losartan and the net effect of perfusion pressure augmentation are the model's most fragile predictions, and",
            "they are flagged as such in the CALIBRATION NOTES of aad_mrgsolve_model_en.R.")
        )
      )
    )
  )
)

# -----------------------------------------------------------------------------
#  server
# -----------------------------------------------------------------------------
server <- function(input, output, session) {

  user_pars <- reactive({
    list(AEN0 = input$AEN0, AEN_REF = 42, GRE0 = input$GRE0,
         DFL0 = input$DFL0, DTL0 = input$DTL0, GEN_MMP = input$GEN_MMP,
         BARO_REF = input$BARO_REF,
         INV_REN = as.numeric(input$INV_REN), OBST_REN = input$OBST_REN,
         BB_AGENT = as.numeric(input$BB), BB_MAX = input$BB_MAX,
         BB_ON_SBP = as.numeric(input$BB_ON_SBP),
         VD_AGENT = as.numeric(input$VD), VD_MAX = input$VD_MAX,
         R_FEN = input$R_FEN, SBP_SET = input$SBP_SET,
         D_MET = input$D_MET, D_LOS = input$D_LOS, ADHERE = input$ADHERE,
         TEVAR_T = input$TEVAR_T, TEVAR_COVER = input$TEVAR_COVER)
  })

  sim <- eventReactive(input$go, {
    p <- user_pars()
    d <- if (input$END > 800) 24 else 0.1
    m <- mod %>% param(p)
    m <- m %>% init(ELN = input$ELN, DFL = input$DFL0, DTL = input$DTL0,
                    AEN = input$AEN0)
    m %>% mrgsim(end = input$END, delta = d) %>% as_tibble()
  }, ignoreNULL = FALSE)

  tail_mean <- function(d, col, hours = 6) {
    tt <- max(d$time) - hours
    mean(d[[col]][d$time >= tt], na.rm = TRUE)
  }

  # ---- 1 -------------------------------------------------------------------
  output$profile <- renderTable({
    d <- sim()
    h_fl <- 1.90 * 0.40
    data.frame(
      Item = c("False-lumen outer wall thickness h_FL (mm)", "True-lumen wall thickness (mm)",
               "Thickness ratio (true/false)", "Total diameter at presentation D_ao (mm)",
               "False-lumen fraction", "False-lumen systolic stress σ_FL (kPa)",
               "True-lumen systolic stress σ_TL (kPa)", "Stress ratio (false/true)",
               "τ_FL (s)", "False-lumen diastolic pressure (mmHg)", "ASI (at presentation)"),
      Value = c(sprintf("%.2f", h_fl), sprintf("%.2f", 1.90 * 0.60),
             sprintf("%.2f", 0.60 / 0.40),
             sprintf("%.1f", d$D_AORTA[1]),
             sprintf("%.2f", d$FL_FRACTION[1]),
             sprintf("%.1f", d$SIGMA_FL_SYS[1]),
             sprintf("%.1f", d$SIGMA_TL[1]),
             sprintf("%.2f", d$SIGMA_FL_SYS[1] / d$SIGMA_TL[1]),
             sprintf("%.3f", d$TAU_FL[1]),
             sprintf("%.1f", d$PFL_DIASTOLIC[1]),
             sprintf("%.2f", d$ASI[1])))
  }, striped = TRUE)

  output$p_geom <- renderPlot({
    d <- sim()
    d %>% select(time, `Total diameter D_ao` = D_AORTA, `False lumen D_FL` = DFL,
                 `True lumen D_TL` = DTL) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
      geom_line(linewidth = 0.9) +
      geom_hline(yintercept = 55, linetype = "dashed", colour = "grey40") +
      annotate("text", x = 0, y = 56.5, hjust = 0, size = 3.2,
               label = "Intervention threshold 55 mm") +
      labs(x = "Time (h)", y = "Diameter (mm)", colour = NULL,
           title = "Time course of the aortic geometry") + THEME
  })

  output$p_stressbar <- renderPlot({
    d <- sim()
    tibble(Compartment = c("False lumen (h 0.76 mm)", "True lumen (h 1.14 mm)", "Intimal flap"),
           Stress = c(tail_mean(d, "SIGMA_FL_SYS"), tail_mean(d, "SIGMA_TL"),
                    tail_mean(d, "SIGMA_FLAP"))) %>%
      ggplot(aes(reorder(Compartment, Stress), Stress, fill = Compartment)) +
      geom_col(width = 0.6, show.legend = FALSE) +
      geom_hline(yintercept = 118, linetype = "dashed") +
      annotate("text", x = 0.6, y = 126, hjust = 0, size = 3.2,
               label = "Normal aorta 118 kPa") +
      coord_flip() +
      labs(x = NULL, y = "Circumferential stress (kPa)",
           title = "Same pressure, different thickness — the asymmetry Laplace creates") + THEME
  })

  # ---- 2 -------------------------------------------------------------------
  output$p_pk <- renderPlot({
    d <- sim()
    d %>% select(time, Esmolol = A_ESM, Nicardipine = A_NICC,
                 Nitroprusside = A_SNP, Labetalol = A_LABC,
                 Clevidipine = A_CLV) %>%
      pivot_longer(-time) %>% filter(value > 1e-9) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      scale_y_log10() +
      labs(x = "Time (h)", y = "Amount in the body (mg, log)", colour = NULL,
           title = "Pharmacokinetics") + THEME
  })

  output$p_occ <- renderPlot({
    d <- sim()
    d %>% select(time, `β1 occupancy` = OCC_BETA1, `DHP occupancy` = OCC_CCB,
                 `NO donation` = OCC_NOD, `AT1 blockade` = OCC_AT1) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      ylim(0, 1) +
      labs(x = "Time (h)", y = "Occupancy / blockade", colour = NULL,
           title = "Receptor occupancy (competitive, additive)") + THEME
  })

  output$p_rate <- renderPlot({
    d <- sim()
    d %>% select(time, `Vasodilator infusion (mg/h)` = RATE_VD,
                 `β-blocker infusion (mg/h)` = RATE_BB) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "Time (h)", y = NULL, colour = NULL,
           title = "Closed-loop titration — the nurse's titration is an integral controller") +
      THEME + theme(legend.position = "none")
  })

  output$t_exposure <- renderTable({
    d <- sim()
    data.frame(
      Drug = c("Esmolol", "Nicardipine", "Nitroprusside", "Labetalol",
               "Clevidipine", "Cyanide", "Thiocyanate"),
      `Steady_state_concentration_mg_per_L` = sprintf(
        "%.4f", c(tail_mean(d, "A_ESM") / 255, tail_mean(d, "A_NICC") / 45,
                  tail_mean(d, "A_SNP") / 15, tail_mean(d, "A_LABC") / 420,
                  tail_mean(d, "A_CLV") / 12, tail_mean(d, "CYANIDE"),
                  tail_mean(d, "THIOCYANATE"))))
  }, striped = TRUE)

  # ---- 3 -------------------------------------------------------------------
  cmp <- reactive({
    keys <- c(input$cmpA, input$cmpB, input$cmpC)
    bind_rows(lapply(keys, function(k)
      simulate_arm(SCENARIOS[[k]]) %>% mutate(arm = k)))
  })

  output$p_dissoc <- renderPlot({
    cmp() %>%
      select(time, arm, `Systolic pressure (mmHg)` = SBP, `dP/dt max (mmHg/s)` = DPDT,
             `ASI (product)` = ASI, `Flap stress-impulse` = ASI_FLAP) %>%
      pivot_longer(c(-time, -arm)) %>%
      ggplot(aes(time, value, colour = arm)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "Time (h)", y = NULL, colour = NULL,
           title = "The pressure comes down but the impulse stays where it was") + THEME
  })

  output$t_dissoc <- renderTable({
    cmp() %>% group_by(arm) %>%
      filter(time >= max(time) - 6) %>%
      summarise(`SBP` = mean(SBP), `HR` = mean(HR), `dP/dt` = mean(DPDT),
                `ASI` = mean(ASI), `ASI_flap` = mean(ASI_FLAP),
                `σ_FL` = mean(SIGMA_FL_SYS), .groups = "drop") %>%
      left_join(cmp() %>% group_by(arm) %>%
                  summarise(`Tear growth %` = 100 * (last(AEN) / first(AEN) - 1),
                            `Mortality %` = last(MORTALITY), .groups = "drop"),
                by = "arm")
  }, digits = 3, striped = TRUE)

  # ---- 4 -------------------------------------------------------------------
  output$p_divider <- renderPlot({
    d <- sim()
    d %>% select(time, `True-lumen systolic pressure` = SBP, `True-lumen diastolic pressure` = DBP,
                 `False-lumen mean pressure` = PFL, `False-lumen systolic pressure` = PFL_SYSTOLIC,
                 `False-lumen diastolic pressure` = PFL_DIASTOLIC) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      labs(x = "Time (h)", y = "Pressure (mmHg)", colour = NULL,
           title = "Pressures in the two lumens — the diastolic pressure is the discriminating variable") + THEME
  })

  output$p_gre_sweep <- renderPlot({
    base <- user_pars()
    gres <- c(0.02, 0.05, 0.1, 0.25, 0.5, 1, 2, 3)
    res <- bind_rows(lapply(gres, function(g) {
      p <- base; p$GRE0 <- g
      d <- mod %>% param(p) %>% mrgsim(end = 24, delta = 1) %>% as_tibble()
      tibble(G_re = g,
             `False-lumen diastolic pressure` = tail(d$PFL_DIASTOLIC, 1),
             `False-lumen mean pressure` = tail(d$PFL, 1),
             `False-lumen flow` = tail(d$FLOW_FALSE, 1) * 100,
             `True-lumen diastolic pressure` = tail(d$DBP, 1))
    }))
    res %>% pivot_longer(-G_re) %>%
      ggplot(aes(G_re, value, colour = name)) +
      geom_line(linewidth = 0.9) + geom_point() + scale_x_log10() +
      labs(x = "Re-entry conductance G_re (log)", y = "mmHg  /  flow ×100",
           colour = NULL,
           title = "Close the re-entry tear and the false-lumen diastolic pressure rises to the systemic mean pressure") +
      THEME
  })

  # ---- 5 -------------------------------------------------------------------
  output$p_sigma <- renderPlot({
    d <- sim()
    d %>% select(time, `False-lumen systolic stress` = SIGMA_FL_SYS,
                 `False-lumen mean stress` = SIGMA_FL, `True-lumen stress` = SIGMA_TL,
                 `Flap stress` = SIGMA_FLAP, `Wall strength S` = WALL_STRENGTH) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      labs(x = "Time (h)", y = "Stress / strength (kPa)", colour = NULL,
           title = "Rupture is where the wall stress meets the wall strength") + THEME
  })

  output$p_hazard <- renderPlot({
    d <- sim()
    d %>% mutate(`Annual rupture hazard (%)` =
                   100 * (1 - exp(-8760 * 1.10e-4 *
                                    exp(pmin((STRESS_RATIO - 1) / 0.155, 30))))) %>%
      select(time, `σ/S` = STRESS_RATIO, `Annual rupture hazard (%)`) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      labs(x = "Time (h)", y = NULL, colour = NULL, title = "Rupture hazard") +
      THEME + theme(legend.position = "none")
  })

  output$p_wall <- renderPlot({
    d <- sim()
    d %>% select(time, `False-lumen outer wall thickness (mm)` = HFL, `Elastin ELN` = ELN,
                 `Adventitial collagen COL` = COL) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      labs(x = "Time (h)", y = NULL, colour = NULL,
           title = "Wall remodelling — collagen is the only negative feedback") +
      THEME + theme(legend.position = "none")
  })

  # ---- 6 -------------------------------------------------------------------
  output$p_paradox <- renderPlot({
    n <- 1.9; s <- 0.85; A <- 1.0
    thr <- seq(0, 1, by = 0.005)
    rel <- (1 + A * thr)^n * (1 - s * thr)
    star <- (n * A - s) / ((1 + n) * s * A)
    tibble(THR = thr, `Relative growth rate` = rel / rel[1]) %>%
      ggplot(aes(THR, `Relative growth rate`)) +
      geom_line(linewidth = 1.1, colour = "#b8860b") +
      geom_vline(xintercept = star, linetype = "dashed") +
      geom_hline(yintercept = 1, linetype = "dotted") +
      annotate("text", x = star + 0.02, y = 0.75, hjust = 0, size = 3.6,
               label = sprintf("Analytic solution THR* = %.3f", star)) +
      labs(x = "False-lumen thrombus fraction THR", y = "Growth rate (relative to no thrombus)",
           title = "The partial thrombosis paradox: support is linear, thinning is a power of 1.9") + THEME
  })

  output$p_thr_time <- renderPlot({
    d <- sim()
    d %>% select(time, `False-lumen body thrombus` = THRFL, `Distal thrombus` = THRD,
                 `Entry conductance` = G_ENTRY, `Re-entry conductance` = G_REENTRY,
                 `False-lumen flow` = FLOW_FALSE) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      labs(x = "Time (h)", y = "Fraction / conductance / flow", colour = NULL,
           title = "Time course of thrombosis — sealing does not start until 0.70 is passed") + THEME
  })

  # ---- 7 -------------------------------------------------------------------
  output$p_perf <- renderPlot({
    d <- sim()
    d %>% select(time, `Renal perfusion` = QREN, `Mesenteric perfusion` = Q_MESENTERIC,
                 `Spinal cord perfusion` = Q_SPINAL, `True-lumen compression` = TL_COLLAPSE) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      labs(x = "Time (h)", y = "Relative perfusion / compression index", colour = NULL,
           title = "Branch perfusion: static and dynamic obstruction") + THEME
  })

  output$p_raas <- renderPlot({
    d <- sim()
    d %>% select(time, `Renin activity` = PRA, `Angiotensin II` = ANGII,
                 `Aldosterone` = ALDO, `SVR` = SVR, `SBP` = SBP) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "Time (h)", y = NULL, colour = NULL,
           title = "The renin trap — lower the blood pressure and the kidney puts it back") +
      THEME + theme(legend.position = "none")
  })

  output$t_cn <- renderTable({
    d <- sim()
    data.frame(
      Item = c("Peak cyanide (mg/L)", "Toxicity threshold (mg/L)",
               "Thiocyanate (mg/L)", "Lactate (mmol/L)",
               "Creatinine (mg/dL)", "eGFR"),
      Value = sprintf("%.2f", c(max(d$CYANIDE), 1.0, tail_mean(d, "THIOCYANATE"),
                            tail_mean(d, "LACT"), tail(d$CREAT, 1),
                            tail_mean(d, "EGFR"))))
  }, striped = TRUE)

  # ---- 8 -------------------------------------------------------------------
  output$p_cord <- renderPlot({
    d <- sim()
    d %>% select(time, `Spinal cord perfusion index` = Q_SPINAL,
                 `Ischaemic burden SCI` = SCI, `Paraplegia probability` = PARAPLEGIA_PROB) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      labs(x = "Time (h)", y = NULL, colour = NULL,
           title = "Spinal cord ischaemia is a recoverable burden, and it has a window") +
      THEME + theme(legend.position = "none")
  })

  output$p_collat <- renderPlot({
    d <- sim()
    d %>% select(time, `Collateral network COLLAT` = COLLAT, `Mean arterial pressure` = MAP) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      labs(x = "Time (h)", y = NULL, colour = NULL,
           title = "Collateral recruitment (baseline 0.55, τ ≈ 5 days)") +
      THEME + theme(legend.position = "none")
  })

  # ---- 9 -------------------------------------------------------------------
  multi <- reactive({
    req(length(input$multi) > 0)
    bind_rows(lapply(input$multi, function(k)
      simulate_arm(SCENARIOS[[k]]) %>% mutate(arm = k)))
  })

  output$p_multi <- renderPlot({
    multi() %>%
      select(time, arm, SBP, HR, DPDT, ASI, SIGMA_FL_SYS, D_AORTA,
             STRESS_RATIO, MORTALITY) %>%
      pivot_longer(c(-time, -arm)) %>%
      ggplot(aes(time, value, colour = arm)) + geom_line(linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "Time (h)", y = NULL, colour = NULL,
           title = "Scenario comparison") +
      THEME + theme(legend.text = element_text(size = 8))
  })

  output$t_multi <- renderTable({
    multi() %>% group_by(arm) %>%
      summarise(`Duration (h)` = max(time),
                SBP = mean(SBP[time >= max(time) - 6]),
                HR = mean(HR[time >= max(time) - 6]),
                `dP/dt` = mean(DPDT[time >= max(time) - 6]),
                ASI = mean(ASI[time >= max(time) - 6]),
                `D_ao (mm)` = last(D_AORTA),
                `σ/S` = mean(STRESS_RATIO[time >= max(time) - 6]),
                `In target %` = 100 * last(TIT) / max(time),
                `Mortality %` = last(MORTALITY), .groups = "drop")
  }, digits = 3, striped = TRUE)

  # ---- 10 ------------------------------------------------------------------
  output$p_end <- renderPlot({
    d <- sim()
    d %>% select(time, `Pain 0-10` = PAIN, `Cumulative time in target (h)` = TIT,
                 `Cumulative hazard` = CUMH, `Survival probability` = SURV) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "Time (h)", y = NULL, colour = NULL,
           title = "Clinical endpoints") +
      THEME + theme(legend.position = "none")
  })

  output$p_bio <- renderPlot({
    d <- sim()
    d %>% select(time, `IL-6` = IL6, `MMP-9` = MMP9, `TGF-β` = TGFB,
                 `D-dimer` = DDIM, `Lactate` = LACT) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      labs(x = "Time (h)", y = "Relative value / concentration", colour = NULL,
           title = "Biomarkers") + THEME
  })

  output$p_growth <- renderPlot({
    d <- sim()
    yrs <- max(d$time) / 8760
    d %>% mutate(`Growth (mm)` = D_AORTA - first(D_AORTA)) %>%
      ggplot(aes(time / 8760, `Growth (mm)`)) +
      geom_line(linewidth = 1, colour = "#b8860b") +
      labs(x = "Time (years)", y = "Diameter increase (mm)",
           title = sprintf("Mean growth rate %.2f mm/year",
                           (last(d$D_AORTA) - first(d$D_AORTA)) / max(yrs, 1e-9))) +
      THEME
  })

  output$t_end <- renderTable({
    d <- sim()
    hit55 <- d$time[which(d$D_AORTA >= 55)][1]
    data.frame(
      Endpoint = c("Final systolic pressure (mmHg)", "Final heart rate (bpm)",
                 "dP/dt max (mmHg/s)", "ASI", "Flap stress-impulse",
                 "Tear growth (%)", "Total diameter (mm)", "Growth rate (mm/year)",
                 "Time to 55 mm (years)", "σ/S", "Time in target (%)",
                 "Cumulative mortality (%)", "Paraplegia probability"),
      Value = c(sprintf("%.1f", tail_mean(d, "SBP")),
             sprintf("%.1f", tail_mean(d, "HR")),
             sprintf("%.0f", tail_mean(d, "DPDT")),
             sprintf("%.3f", tail_mean(d, "ASI")),
             sprintf("%.3f", tail_mean(d, "ASI_FLAP")),
             sprintf("%.2f", 100 * (tail(d$AEN, 1) / d$AEN[1] - 1)),
             sprintf("%.1f", tail(d$D_AORTA, 1)),
             sprintf("%.2f", (tail(d$D_AORTA, 1) - d$D_AORTA[1]) /
                       max(max(d$time) / 8760, 1e-9)),
             ifelse(is.na(hit55), "Not reached", sprintf("%.2f", hit55 / 8760)),
             sprintf("%.3f", tail_mean(d, "STRESS_RATIO")),
             sprintf("%.1f", 100 * tail(d$TIT, 1) / max(d$time)),
             sprintf("%.2f", tail(d$MORTALITY, 1)),
             sprintf("%.3f", max(d$PARAPLEGIA_PROB))))
  }, striped = TRUE)
}

shinyApp(ui, server)
