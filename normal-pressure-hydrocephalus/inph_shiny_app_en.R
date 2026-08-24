## =============================================================================
##  inph_shiny_app.R
##  Interactive dashboard for the iNPH QSP model
##  Idiopathic normal pressure hydrocephalus (iNPH) QSP model — interactive dashboard
## =============================================================================
##
##  13 tabs:
##    1  Patient profile         patient phenotype and the hydraulic baseline
##    2  CSF dynamics            Davson / Marmarou — pressure, compliance, amplitude
##    3  Shunt hydraulics        valve + gravitational unit + posture, resolved
##    4  Valve titration map     the titration map (the central clinical result)
##    5  Pressure-volume curve   P-V curve, pulse amplitude, RAP
##    6  Ventricular morphology  Evans index, callosal angle, DESH, two-clock creep
##    7  White matter injury     recoverable vs permanent white matter
##    8  Clinical triad          gait / cognition / continence, fast + slow arms
##    9  Diagnostic tests        tap test and external lumbar drain simulator
##   10  Drug PK/PD              acetazolamide (saturable RBC binding) and others
##   11  CSF biomarkers          the dilution/flux decomposition
##   12  Scenario comparison     all 21 scenarios side by side
##   13  Safety                  subdural collection, headache, acidosis
##
##  Run:  shiny::runApp("inph_shiny_app.R")
##  Requires inph_mrgsolve_model_en.R in the same directory.
## =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)

source("inph_mrgsolve_model_en.R", local = FALSE)

THEME <- theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "#eef3f8"),
        plot.title = element_text(face = "bold"))

## Colour scale used consistently across tabs
PAL <- c("#2c6fb5", "#c0504d", "#4e8b4e", "#c98a2e", "#7b4fa8", "#2f8177",
         "#b0518a", "#8a8a3f")

fmt <- function(x, d = 2) formatC(x, format = "f", digits = d)

# --------------------------------------------------------------------------- #
#  UI
# --------------------------------------------------------------------------- #
ui <- fluidPage(
  titlePanel(paste("Idiopathic Normal Pressure Hydrocephalus (iNPH)",
                   "QSP dashboard")),
  tags$p(style = "color:#555;margin-top:-8px",
         paste("Mean intracranial pressure is normal by definition. The real",
               "variables of the pathophysiology are outflow resistance (Rout) · compliance (C) · pulsatility (AMP),",
               "and the safety margin of treatment is set not by biology but by hydrostatics.")),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Patient phenotype"),
      sliderInput("Rout_init", "CSF outflow resistance Rout [mmHg/(mL/min)]",
                  6, 26, 19, 0.5),
      ## E1 is defined against the P-V reference pressure P0_marm = -2 mmHg;
      ## healthy ~0.156, iNPH ~0.222, post-decompression floor 0.162.
      sliderInput("E1_init", "Elastance coefficient E1 [1/mL] (inverse of compliance)",
                  0.13, 0.34, 0.222, 0.005),
      sliderInput("Vv0", "Ventricular volume Vv [mL]", 30, 160, 95, 5),
      sliderInput("WM_init", "White matter functional integrity WMint", 0.35, 0.95, 0.68, 0.01),
      sliderInput("WMperm_init", "Permanent axonal loss WMperm", 0.0, 0.45, 0.10, 0.01),
      sliderInput("APOE", "APOE coefficient (1.0 = e3/e3, 2.2 = e4)",
                  1.0, 2.5, 1.0, 0.1),
      sliderInput("atrophy", "Cortical atrophy (size of the subdural space)", 0, 1, 0.55, 0.05),
      sliderInput("f_up", "Upright fraction of the day f_up", 0.2, 0.9, 0.60, 0.05),
      sliderInput("comorb", "Non-hydrocephalic gait impairment burden", 0, 0.5, 0, 0.05),

      hr(),
      h4("Shunt hardware"),
      checkboxInput("shunt", "Shunt inserted (VP shunt)", FALSE),
      sliderInput("Popen_cm", "Valve opening pressure Popen [cmH2O]", 0, 30, 10, 1),
      sliderInput("Ggrav_cm", "Gravitational unit [cmH2O] (active only when upright)",
                  0, 40, 0, 5),
      sliderInput("asd_eff", "Membrane anti-siphon efficiency", 0, 1, 0, 0.05),
      sliderInput("Rsh", "Shunt resistance Rsh [mmHg/(mL/min)]", 1, 40, 2.5, 0.5),
      sliderInput("Hcol_cm", "Hydrostatic column height [cmH2O]", 20, 60, 45, 1),
      checkboxInput("etv", "ETV (endoscopic third ventriculostomy)", FALSE),

      hr(),
      h4("Drugs"),
      selectInput("az_dose", "Acetazolamide (BID)",
                  c("None" = "0", "250 mg" = "250", "500 mg" = "500",
                    "1000 mg" = "1000"), "0"),
      checkboxInput("melatonin", "Melatonin 2 mg at bedtime", FALSE),
      checkboxInput("solifenacin", "Solifenacin 5 mg", FALSE),
      checkboxInput("donepezil", "Donepezil 10 mg", FALSE),
      checkboxInput("antithrombotic", "On antithrombotic therapy", FALSE),

      hr(),
      sliderInput("tend", "Simulation duration [months]", 3, 60, 24, 3),
      actionButton("go", "Run simulation", class = "btn-primary"),
      tags$p(style = "font-size:11px;color:#777;margin-top:10px",
             paste("This is a model for education and research. It must not be",
                   "used for clinical decision-making."))
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs", type = "tabs",

        tabPanel("1 Patient profile",
                 h4("Baseline hydraulics (closed-form flow balance)"),
                 tableOutput("tbl_baseline"),
                 h4("Comparison with a healthy 75-year-old control"),
                 tableOutput("tbl_vs_control"),
                 plotOutput("plt_profile", height = "320px"),
                 tags$small(paste("Rout rises more than twofold, yet mean intracranial",
                                  "pressure never leaves the normal 7-15 mmHg range.",
                                  "The pulse amplitude (AMP) amplifies that loss."))),

        tabPanel("2 CSF dynamics",
                 plotOutput("plt_hydro", height = "560px"),
                 tags$small(paste("Davson: ICP = Pss + If·Rout.",
                                  "Marmarou: C = 1/(E1·(P−P0)).",
                                  "AMP = ΔVp/C."))),

        tabPanel("3 Shunt hydraulics",
                 h4("Resolved by posture (supine vs upright)"),
                 tableOutput("tbl_posture"),
                 plotOutput("plt_posture", height = "380px"),
                 tags$small(paste("When upright, the ventriculo-peritoneal hydrostatic",
                                  "column adds directly to the driving pressure. A gravitational unit",
                                  "raises the opening pressure only when upright, cancelling that term."))),

        tabPanel("4 Valve titration map",
                 h4("Opening-pressure sweep — with and without a gravitational unit"),
                 actionButton("run_map", "Compute the titration map (takes tens of seconds)"),
                 plotOutput("plt_map", height = "560px"),
                 tableOutput("tbl_map")),

        tabPanel("5 Pressure-volume curve",
                 plotOutput("plt_pv", height = "420px"),
                 plotOutput("plt_amp", height = "300px")),

        tabPanel("6 Ventricular morphology",
                 plotOutput("plt_morph", height = "520px"),
                 tags$small(paste("Even after a successful shunt the Evans index barely",
                                  "changes — because the plastic component of the",
                                  "dilatation does not reverse."))),

        tabPanel("7 White matter injury",
                 plotOutput("plt_wm", height = "520px"),
                 tags$small(paste("WMint is the recoverable pool, WMperm the",
                                  "irreversible one. A shunt reverses the former and",
                                  "only stops the clock on the latter."))),

        tabPanel("8 Clinical triad",
                 plotOutput("plt_triad", height = "560px"),
                 tableOutput("tbl_triad")),

        tabPanel("9 Diagnostic tests",
                 h4("Tap test / continuous drainage simulator"),
                 fluidRow(
                   column(4, sliderInput("tap_vol", "Volume removed at tap [mL]",
                                         0, 60, 40, 5)),
                   column(4, sliderInput("eld_h", "Continuous drainage duration [hours]",
                                         0, 120, 72, 12)),
                   column(4, sliderInput("thr", "Response threshold [m/s]",
                                         0.02, 0.20, 0.10, 0.01))),
                 plotOutput("plt_dx", height = "440px"),
                 tableOutput("tbl_dx"),
                 tags$small(paste("τ = Rout·C ≈ 6 min is the linearised value and does",
                                  "not apply to a 40 mL tap. Below sinus pressure",
                                  "absorption falls to zero, so recovery is limited by the production rate",
                                  "(0.35 mL/min) and takes about 114 min.",
                                  "The fast component is asymmetric — 29 min to charge,",
                                  "21.6 h to discharge."))),

        tabPanel("10 Drug PK/PD",
                 plotOutput("plt_pk", height = "440px"),
                 tableOutput("tbl_pk"),
                 tags$small(paste("Acetazolamide binds erythrocyte carbonic anhydrase",
                                  "saturably — which is why the apparent",
                                  "half-life lengthens with dose. The PD has a low EC50, so it is",
                                  "on a plateau above 250 mg BID."))),

        tabPanel("11 CSF biomarkers",
                 plotOutput("plt_bio", height = "480px"),
                 tableOutput("tbl_bio"),
                 tags$small(paste("c = J(brain→CSF) / [(If+Qsh)/Vcsf + kdeg].",
                                  "In iNPH the numerator (glymphatic efflux) and the",
                                  "turnover in the denominator fall together, so both Aβ42 and p-tau",
                                  "are measured low."))),

        tabPanel("12 Scenario comparison",
                 actionButton("run_scen", "Run all 21 scenarios (takes minutes)"),
                 tableOutput("tbl_scen"),
                 plotOutput("plt_scen", height = "460px")),

        tabPanel("13 Safety",
                 plotOutput("plt_safety", height = "520px"),
                 tableOutput("tbl_safety"),
                 tags$small(paste("The incidence of subdural collections is ordinal in",
                                  "this model — because k_haz was not",
                                  "calibrated — and it must not be read as an absolute",
                                  "incidence.")))
      )
    )
  )
)

# --------------------------------------------------------------------------- #
#  Server
# --------------------------------------------------------------------------- #
server <- function(input, output, session) {

  ## ---- assemble the parameter overrides from the sidebar ---------------- ##
  overrides <- reactive({
    list(Rout_init = input$Rout_init, atrophy = input$atrophy,
         f_up = input$f_up, comorb = input$comorb, APOE = input$APOE,
         shunt = as.numeric(input$shunt), Popen_cm = input$Popen_cm,
         Ggrav_cm = input$Ggrav_cm, asd_eff = input$asd_eff,
         Rsh = input$Rsh, Hcol_cm = input$Hcol_cm,
         etv = as.numeric(input$etv),
         antithrombotic = as.numeric(input$antithrombotic))
  })

  ## The phenotype sliders that are INITIAL CONDITIONS rather than parameters
  ## have to be pushed into the init vector after inph_init() has built it.
  patient_now <- reactive({
    p <- as.list(param(mod))
    p[names(overrides())] <- overrides()
    y <- inph_init(p, healthy = FALSE)
    y["Vv"] <- input$Vv0
    y["Vplast"] <- p$f_plastic * (input$Vv0 - p$Vv_norm)
    y["E1"] <- input$E1_init
    y["WMint"] <- input$WM_init
    y["WMperm"] <- input$WMperm_init
    ## re-seat the clinical scores and the references on the edited state
    yl <- as.list(y)
    hy <- inph_hydro(yl, p)
    wex <- max(0, yl$Wpv - p$Wpv_norm)
    CBF <- p$CBF0 * (yl$Autoreg / p$Autoreg_norm) /
      (1 + p$kCBF_W * wex + p$kCBF_P * hy$dPtmP)
    perf <- 0.55 + 0.45 * min(1, CBF / p$CBF0)
    y["Gslow"] <- max(0.05, p$G_max * yl$WMint^p$G_wm_pow * perf *
                        (1 - 0.25 * yl$Ab_plq) * (1 - p$comorb))
    taub <- min(0.4, 0.4 * max(0, yl$Tau_isf - 1) / 3)
    y["Cslow"] <- max(5, p$MMSE_max - p$kcog_wm * (1 - yl$WMint) -
                        p$kcog_ad * min(1, 0.6 * yl$Ab_plq + taub))
    y["Urin"] <- max(0, min(p$Ur_max, p$kur_wm * (1 - yl$WMint)))
    p$AMP_ref <- hy$AMPday
    p$Pday_ref <- hy$Pday
    p$E1_floor_eff <- min(p$E1_floor, yl$E1)
    list(p = p, init = y, hyd = hy, CBF = CBF)
  })

  ## ---- drug events ------------------------------------------------------ ##
  drug_events <- reactive({
    tend <- input$tend * 30.4
    e <- NULL
    add <- function(a, b) if (is.null(a)) b else c(a, b)
    azd <- as.numeric(input$az_dose)
    if (azd > 0)
      e <- add(e, dose_ev("Aaz_g", azd, bioav = 0.9, ii = 0.5,
                          addl = ceiling(tend / 0.5) - 1))
    if (input$melatonin)
      e <- add(e, dose_ev("Ame_g", 2, bioav = 0.15, ii = 1,
                          addl = ceiling(tend) - 1))
    if (input$solifenacin)
      e <- add(e, dose_ev("Aso_g", 5, bioav = 0.9, ii = 1,
                          addl = ceiling(tend) - 1))
    if (input$donepezil)
      e <- add(e, dose_ev("Ado_g", 10, bioav = 1.0, ii = 1,
                          addl = ceiling(tend) - 1))
    e
  })

  ## ---- the main simulation --------------------------------------------- ##
  sim <- eventReactive(input$go, {
    pt <- patient_now()
    tend <- input$tend * 30.4
    m <- mod %>% param(pt$p) %>% init(pt$init)
    ev0 <- drug_events()
    out <- if (is.null(ev0)) mrgsim(m, end = tend, delta = 1, hmax = 0.25)
           else mrgsim(m, data = ev0, end = tend, delta = 1, hmax = 0.25)
    as_tibble(out)
  }, ignoreNULL = FALSE)

  ## ---- reference: the untreated / healthy comparators ------------------ ##
  control <- reactive({
    p <- as.list(param(mod))
    r <- inph_refs(p, healthy = TRUE)
    p[names(r)] <- r
    y <- inph_init(p, healthy = TRUE)
    hy <- inph_hydro(as.list(y), p)
    list(p = p, init = y, hyd = hy)
  })

  ## ===================== TAB 1 — patient profile ======================== ##
  output$tbl_baseline <- renderTable({
    pt <- patient_now(); h <- pt$hyd
    tibble(
      Item = c("Rout [mmHg/(mL/min)]", "E1 [1/mL]",
               "Mean ICP supine [mmHg]", "Mean ICP upright [mmHg]",
               "Daily mean ICP [mmHg]", "Compliance C (supine) [mL/mmHg]",
               "ICP pulse amplitude AMP [mmHg]", "Pulsatile transmantle gradient [mmHg]",
               "Mean transmantle gradient [mmHg]",
               "Shunt flow [mL/day]", "CSF production [mL/day]",
               "Periventricular white matter perfusion [mL/100g/min]",
               "τ = Rout·C (linearised) [min]",
               "Minimum recovery time after a 40 mL tap [min]"),
      Value = c(fmt(pt$init[["Rout"]], 1), fmt(pt$init[["E1"]], 3),
             fmt(h$Ps), fmt(h$Pu), fmt(h$Pday), fmt(h$Csup, 3),
             fmt(h$AMP), fmt(h$dPtmP, 3), fmt(h$dPtmM, 3),
             fmt(h$Qsh, 0), fmt(h$Ifd, 0), fmt(pt$CBF, 1),
             ## Rout [mmHg/(mL/min)] x C [mL/mmHg] is already in minutes
             fmt(pt$init[["Rout"]] * h$Csup, 1),
             ## but the tap-test recovery is PRODUCTION-limited, not
             ## resistance-limited: below sinus pressure nothing is absorbed
             fmt(40 / pt$p$If0, 0))
    )
  }, striped = TRUE)

  output$tbl_vs_control <- renderTable({
    pt <- patient_now(); ct <- control()
    a <- ct$hyd; b <- pt$hyd
    rat <- function(x, y) if (abs(x) > 1e-9) paste0(fmt(y / x), "x") else "n/a"
    tibble(
      Item = c("Rout", "Mean ICP (supine)", "AMP (pulse amplitude)",
               "Pulsatile transmantle gradient", "Evans index"),
      `Healthy 75 y` = c(fmt(ct$init[["Rout"]], 1), fmt(a$Ps), fmt(a$AMP),
                     fmt(a$dPtmP, 3),
                     fmt(0.20 + 0.0035 * (ct$init[["Vv"]] - 25), 3)),
      `Patient` = c(fmt(pt$init[["Rout"]], 1), fmt(b$Ps), fmt(b$AMP),
                fmt(b$dPtmP, 3),
                fmt(0.20 + 0.0035 * (pt$init[["Vv"]] - 25), 3)),
      `Ratio` = c(rat(ct$init[["Rout"]], pt$init[["Rout"]]),
              rat(a$Ps, b$Ps), rat(a$AMP, b$AMP),
              rat(a$dPtmP, b$dPtmP),
              rat(0.20 + 0.0035 * (ct$init[["Vv"]] - 25),
                  0.20 + 0.0035 * (pt$init[["Vv"]] - 25)))
    )
  }, striped = TRUE)

  output$plt_profile <- renderPlot({
    ## mean ICP and pulse amplitude as functions of Rout, at this patient's E1
    pt <- patient_now(); p <- pt$p
    grid <- lapply(seq(6, 26, 0.5), function(R) {
      y <- as.list(pt$init); y$Rout <- R
      h <- inph_hydro(y, p)
      tibble(Rout = R, `Mean ICP (supine)` = h$Ps, `AMP (pulse amplitude)` = h$AMP)
    }) %>% bind_rows() %>% pivot_longer(-Rout)
    ggplot(grid, aes(Rout, value, colour = name)) +
      geom_line(linewidth = 1) +
      geom_vline(xintercept = pt$init[["Rout"]], linetype = 2,
                 colour = "#555555") +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = 7, ymax = 15,
               alpha = 0.08, fill = "#2c6fb5") +
      annotate("text", x = 7.5, y = 15.6, hjust = 0, size = 3.4,
               label = "Normal ICP range 7-15 mmHg") +
      geom_hline(yintercept = 4, linetype = 3, colour = "#c0504d") +
      annotate("text", x = 24, y = 4.4, hjust = 1, size = 3.4,
               colour = "#c0504d", label = "AMP 4 mmHg threshold") +
      scale_colour_manual(values = PAL) +
      labs(x = "CSF outflow resistance Rout [mmHg/(mL/min)]", y = "mmHg", colour = NULL,
           title = "The disease barely shows in the variable its own name points to") +
      THEME
  })

  ## ===================== TAB 2 — CSF hydrodynamics ====================== ##
  output$plt_hydro <- renderPlot({
    d <- sim()
    d %>% select(time, ICP_sup, ICP_up, ICP_day, AMP, Cspine, dPtm_pulse,
                 dPtm_mean, Qsh_day, If_day) %>%
      pivot_longer(-time) %>%
      mutate(name = factor(name, levels = c(
        "ICP_sup", "ICP_up", "ICP_day", "AMP", "Cspine", "dPtm_pulse",
        "dPtm_mean", "Qsh_day", "If_day"))) %>%
      ggplot(aes(time / 30.4, value)) +
      geom_line(colour = PAL[1], linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "Months", y = NULL, title = "Time course of CSF dynamics") + THEME
  })

  ## ===================== TAB 3 — shunt hydraulics ======================= ##
  posture_table <- reactive({
    pt <- patient_now(); p <- pt$p; y <- as.list(pt$init)
    cm <- function(x) x / 1.36
    rows <- lapply(c(FALSE, TRUE), function(up) {
      pp <- p; pp$f_up <- if (up) 1 else 0
      h <- inph_hydro(y, pp)
      tibble(Posture = if (up) "Upright" else "Supine",
             `ICP [mmHg]` = fmt(if (up) h$Pu else h$Ps),
             `Hydrostatic column [mmHg]` = fmt(if (up)
               cm(p$Hcol_cm) * (1 - p$asd_eff) else 0),
             `Effective opening pressure [mmHg]` = fmt(cm(p$Popen_cm) +
                                         if (up) cm(p$Ggrav_cm) else 0),
             `Distal pressure [mmHg]` = fmt(cm(if (up) p$Pd_up_cm else p$Pd_sup_cm)),
             `Shunt flow [mL/day]` = fmt(h$Qsh, 0),
             `Fraction of production` = fmt(h$Qsh / h$Ifd))
    })
    bind_rows(rows)
  })

  output$tbl_posture <- renderTable(posture_table(), striped = TRUE)

  output$plt_posture <- renderPlot({
    pt <- patient_now(); p <- pt$p; y <- as.list(pt$init)
    grid <- lapply(seq(0, 30, 1), function(Po) {
      bind_rows(lapply(c(0, 30), function(G) {
        pp <- p; pp$shunt <- 1; pp$Popen_cm <- Po; pp$Ggrav_cm <- G
        h <- inph_hydro(y, pp)
        tibble(Popen = Po,
               arm = if (G > 0) "With gravitational unit 30 cmH2O" else "No protection",
               `Daily mean ICP` = h$Pday, `Upright ICP` = h$Pu,
               `Drainage mL/day` = h$Qsh)
      }))
    }) %>% bind_rows() %>% pivot_longer(-c(Popen, arm))
    ggplot(grid, aes(Popen, value, colour = arm)) +
      geom_line(linewidth = 1) +
      geom_hline(yintercept = 0, linetype = 3) +
      facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = PAL[c(2, 3)]) +
      labs(x = "Valve opening pressure [cmH2O]", y = NULL, colour = NULL,
           title = "Hydrostatics sets the safety margin") + THEME
  })

  ## ===================== TAB 4 — titration map ========================== ##
  map_data <- eventReactive(input$run_map, {
    withProgress(message = "Computing the titration map...", {
      pt <- patient_now()
      tend <- min(input$tend * 30.4, 730)
      rows <- list()
      grid <- expand.grid(G = c(0, 30),
                          Po = c(2, 4, 6, 8, 10, 12, 14, 16, 20, 24, 28))
      for (i in seq_len(nrow(grid))) {
        incProgress(1 / nrow(grid))
        p <- pt$p; p$shunt <- 1
        p$Popen_cm <- grid$Po[i]; p$Ggrav_cm <- grid$G[i]
        m <- mod %>% param(p) %>% init(pt$init)
        o <- as_tibble(mrgsim(m, end = tend, delta = 30, hmax = 0.25))
        rows[[i]] <- tibble(
          Ggrav = grid$G[i], Popen = grid$Po[i],
          ICP_day = last(o$ICP_day), ICP_up = last(o$ICP_up),
          AMP = last(o$AMP), dGait = last(o$gait) - first(o$gait),
          Vsdh = last(o$Vsdh), SDH_pct = 100 * last(o$SDH_inc),
          headache = last(o$Headx))
      }
      bind_rows(rows) %>%
        mutate(utility = dGait - 1.1 * SDH_pct / 100 - 0.03 * headache,
               arm = ifelse(Ggrav > 0, "With gravitational unit 30 cmH2O",
                            "No protection"))
    })
  })

  output$plt_map <- renderPlot({
    tm <- map_data()
    tm %>% select(arm, Popen, dGait, SDH_pct, Vsdh, utility) %>%
      pivot_longer(-c(arm, Popen)) %>%
      ggplot(aes(Popen, value, colour = arm)) +
      geom_line(linewidth = 1) + geom_point(size = 1.8) +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      scale_colour_manual(values = PAL[c(3, 2)]) +
      labs(x = "Valve opening pressure [cmH2O]", y = NULL, colour = NULL,
           title = "Valve titration map — the gravitational unit moves the optimum itself") +
      THEME
  })

  output$tbl_map <- renderTable({
    tm <- map_data()
    tm %>% group_by(arm) %>%
      slice_max(utility, n = 1) %>%
      transmute(arm, `Optimal opening pressure [cmH2O]` = Popen,
                `Gait gain [m/s]` = fmt(dGait, 3),
                `Subdural incidence [%]` = fmt(SDH_pct, 1),
                `Collection [mL]` = fmt(Vsdh, 1),
                `Daily mean ICP [mmHg]` = fmt(ICP_day)) %>% ungroup()
  }, striped = TRUE)

  ## ===================== TAB 5 — P-V curve ============================== ##
  output$plt_pv <- renderPlot({
    pt <- patient_now(); ct <- control()
    mk <- function(E1, lab) {
      P <- seq(-4, 40, 0.25)
      tibble(P = P, C = pmin(1 / (E1 * pmax(P - pt$p$P0_marm, 0.5)),
                             pt$p$C_max), arm = lab)
    }
    d <- bind_rows(mk(pt$init[["E1"]], sprintf("Patient E1 = %.3f",
                                               pt$init[["E1"]])),
                   mk(ct$init[["E1"]], sprintf("Healthy control E1 = %.3f",
                                               ct$init[["E1"]])))
    ggplot(d, aes(P, C, colour = arm)) + geom_line(linewidth = 1) +
      geom_vline(xintercept = pt$hyd$Ps, linetype = 2) +
      scale_colour_manual(values = PAL[c(1, 3)]) +
      labs(x = "ICP [mmHg]", y = "Compliance C [mL/mmHg]", colour = NULL,
           title = "Marmarou exponential pressure-volume curve",
           subtitle = paste("Dashed line = the patient's current operating point.",
                            "Low compliance turns the same arterial pulse into a",
                            "larger pressure waveform.")) + THEME
  })

  output$plt_amp <- renderPlot({
    d <- sim()
    d %>% select(time, AMP, AMP_day, Cspine) %>% pivot_longer(-time) %>%
      ggplot(aes(time / 30.4, value, colour = name)) +
      geom_line(linewidth = 0.9) +
      scale_colour_manual(values = PAL) +
      labs(x = "Months", y = NULL, colour = NULL,
           title = "Time course of pulse amplitude and compliance") + THEME
  })

  ## ===================== TAB 6 — morphology ============================= ##
  output$plt_morph <- renderPlot({
    d <- sim()
    d %>% select(time, Vv, Vsas, Vplast, EvansIdx, CallAngle, dPtm_pulse) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time / 30.4, value)) +
      geom_line(colour = PAL[4], linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "Months", y = NULL,
           title = "Ventricular morphology — the plastic dilatation does not reverse") + THEME
  })

  ## ===================== TAB 7 — white matter =========================== ##
  output$plt_wm <- renderPlot({
    d <- sim()
    d %>% select(time, WMint, WMperm, Myel, Wpv, AQ, CBF_pv, Astro, Micro) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time / 30.4, value)) +
      geom_line(colour = PAL[3], linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 4) +
      labs(x = "Months", y = NULL,
           title = "Two clocks: the recoverable pool and the irreversible pool") + THEME
  })

  ## ===================== TAB 8 — clinical triad ========================= ##
  output$plt_triad <- renderPlot({
    d <- sim()
    d %>% select(time, gait, MMSE, urin, iNPHGS, Gfast, Gslow) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time / 30.4, value)) +
      geom_line(colour = PAL[1], linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "Months", y = NULL,
           title = "Hakim triad — a fast component (pressure coupling) and a slow one (tissue recovery)") +
      THEME
  })

  output$tbl_triad <- renderTable({
    d <- sim()
    tibble(Item = c("Gait speed [m/s]", "MMSE", "Urinary urgency index", "iNPHGS"),
           Start = c(fmt(first(d$gait), 3), fmt(first(d$MMSE), 1),
                   fmt(first(d$urin)), fmt(first(d$iNPHGS), 0)),
           End = c(fmt(last(d$gait), 3), fmt(last(d$MMSE), 1),
                   fmt(last(d$urin)), fmt(last(d$iNPHGS), 0)),
           Change = c(fmt(last(d$gait) - first(d$gait), 3),
                   fmt(last(d$MMSE) - first(d$MMSE), 1),
                   fmt(last(d$urin) - first(d$urin)),
                   fmt(last(d$iNPHGS) - first(d$iNPHGS), 0)))
  }, striped = TRUE)

  ## ===================== TAB 9 — diagnostics ============================ ##
  dx_sim <- reactive({
    pt <- patient_now()
    m <- mod %>% param(pt$p) %>% init(pt$init)
    tap <- as_tibble(mrgsim(m, data = tap_ev(0.5, input$tap_vol),
                            end = 8, delta = 0.01, hmax = 0.02)) %>%
      mutate(arm = sprintf("Tap test %g mL", input$tap_vol))
    ph <- pt$p; ph$eld_rate <- 10 / 60
    me <- mod %>% param(ph) %>% init(pt$init)
    ## drain for input$eld_h hours, then stop
    tsw <- input$eld_h / 24
    o1 <- as_tibble(mrgsim(me, end = max(tsw, 0.01), delta = 0.01, hmax = 0.02))
    y2 <- unlist(o1[nrow(o1), names(pt$init)])
    m2 <- mod %>% param(pt$p) %>% init(y2)
    o2 <- as_tibble(mrgsim(m2, end = 8 - tsw, delta = 0.01, hmax = 0.02)) %>%
      mutate(time = time + tsw)
    eld <- bind_rows(o1, o2) %>%
      mutate(arm = sprintf("ELD 10 mL/h x %g h", input$eld_h))
    bind_rows(tap, eld)
  })

  output$plt_dx <- renderPlot({
    dx_sim() %>% select(time, arm, ICP_sup, AMP, gait) %>%
      pivot_longer(-c(time, arm)) %>%
      ggplot(aes(time * 24, value, colour = arm)) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      scale_colour_manual(values = PAL[c(1, 2)]) +
      labs(x = "Time [h]", y = NULL, colour = NULL,
           title = "Diagnostic perturbations: a volume step (tap) vs a flow step (continuous drainage)") +
      THEME
  })

  output$tbl_dx <- renderTable({
    d <- dx_sim()
    d %>% group_by(arm) %>%
      summarise(`Baseline gait` = fmt(first(gait), 3),
                `Peak gait` = fmt(max(gait), 3),
                `Peak gain` = fmt(max(gait) - first(gait), 3),
                `Threshold passed` = ifelse(max(gait) - first(gait) >= input$thr,
                                    "Positive", "Negative"),
                `Minimum ICP` = fmt(min(ICP_sup)), .groups = "drop")
  }, striped = TRUE)

  ## ===================== TAB 10 — PK/PD ================================= ##
  output$plt_pk <- renderPlot({
    pt <- patient_now()
    bind_rows(lapply(c(125, 250, 500, 1000), function(dz) {
      m <- mod %>% param(pt$p) %>% init(pt$init)
      as_tibble(mrgsim(m, data = dose_ev("Aaz_g", dz, bioav = 0.9, ii = 0.5,
                                         addl = 39),
                       end = 20, delta = 0.05, hmax = 0.02)) %>%
        mutate(dose = sprintf("%g mg BID", dz))
    })) %>%
      select(time, dose, AZ_free, AZ_eff, AZ_rbcfrac, If_day, ICP_sup) %>%
      pivot_longer(-c(time, dose)) %>%
      ggplot(aes(time, value, colour = dose)) +
      geom_line(linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = PAL) +
      labs(x = "Days", y = NULL, colour = NULL,
           title = "Acetazolamide: saturable erythrocyte binding + PD plateau") + THEME
  })

  output$tbl_pk <- renderTable({
    pt <- patient_now()
    bind_rows(lapply(c(0, 125, 250, 500, 1000), function(dz) {
      m <- mod %>% param(pt$p) %>% init(pt$init)
      o <- if (dz == 0) as_tibble(mrgsim(m, end = 20, delta = 1, hmax = 0.25))
           else as_tibble(mrgsim(m, data = dose_ev("Aaz_g", dz, bioav = 0.9,
                                                   ii = 0.5, addl = 39),
                                 end = 20, delta = 1, hmax = 0.02))
      tibble(`Dose (BID)` = if (dz == 0) "None" else sprintf("%g mg", dz),
             `Free concentration [mg/L]` = fmt(last(o$AZ_free), 4),
             `CSF production inhibition [%]` = fmt(100 * last(o$AZ_eff), 1),
             `Erythrocyte-bound fraction` = fmt(last(o$AZ_rbcfrac), 3),
             `CSF production [mL/day]` = fmt(last(o$If_day), 0),
             `ICP supine [mmHg]` = fmt(last(o$ICP_sup)),
             `HCO3 [mmol/L]` = fmt(last(o$HCO3), 1))
    }))
  }, striped = TRUE)

  ## ===================== TAB 11 — biomarkers ============================ ##
  output$plt_bio <- renderPlot({
    d <- sim()
    d %>% select(time, CSF_Ab42, CSF_pTau, CSF_NfL, CSF_LRG, turnover, AQ) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time / 30.4, value)) +
      geom_line(colour = PAL[8], linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "Months", y = NULL,
           title = "CSF biomarkers — dilution combined with reduced efflux") + THEME
  })

  output$tbl_bio <- renderTable({
    d <- sim(); ct <- control()
    pc <- as.list(ct$p); yc <- as.list(ct$init)
    Vc <- yc$Vv + yc$Vsas
    koutc <- pc$If0 * 1440 / Vc + pc$kdeg_csf
    glyc <- yc$AQ / (1 + 0.35 * yc$Wpv)
    tibble(
      Item = c("CSF turnover [/day]", "CSF Aβ42 [pg/mL]", "CSF p-tau [pg/mL]",
               "CSF NfL [pg/mL]", "CSF LRG [a.u.]", "AQP4 polarisation"),
      `Healthy control` = c(fmt(pc$If0 * 1440 / Vc, 3),
                     fmt(pc$kflux_ab * glyc * yc$Ab_isf / koutc / Vc, 0),
                     fmt(pc$kflux_tau * glyc * yc$Tau_isf / koutc / Vc, 1),
                     "-", "-", fmt(yc$AQ, 3)),
      Start = c(fmt(first(d$turnover), 3), fmt(first(d$CSF_Ab42), 0),
              fmt(first(d$CSF_pTau), 1), fmt(first(d$CSF_NfL), 0),
              fmt(first(d$CSF_LRG)), "-"),
      End = c(fmt(last(d$turnover), 3), fmt(last(d$CSF_Ab42), 0),
              fmt(last(d$CSF_pTau), 1), fmt(last(d$CSF_NfL), 0),
              fmt(last(d$CSF_LRG)), "-"))
  }, striped = TRUE)

  ## ===================== TAB 12 — scenarios ============================= ##
  scen_data <- eventReactive(input$run_scen, {
    withProgress(message = "Running 21 scenarios...", {
      res <- lapply(names(scen), function(nm) {
        incProgress(1 / length(scen), detail = nm)
        scen[[nm]]()
      })
      bind_rows(res)
    })
  })

  output$tbl_scen <- renderTable({
    summarise_scenarios(scen_data()) %>%
      mutate(across(where(is.numeric), ~ round(.x, 3)))
  }, striped = TRUE)

  output$plt_scen <- renderPlot({
    scen_data() %>%
      filter(time <= 730) %>%
      ggplot(aes(time / 30.4, gait, colour = scenario)) +
      geom_line(linewidth = 0.7, show.legend = FALSE) +
      facet_wrap(~scenario, ncol = 5) +
      labs(x = "Months", y = "Gait speed [m/s]",
           title = "Gait trajectories of the 21 scenarios") + THEME
  })

  ## ===================== TAB 13 — safety ================================ ##
  output$plt_safety <- renderPlot({
    d <- sim()
    d %>% mutate(SDH_pct = 100 * SDH_inc) %>%
      select(time, Vsdh, SDH_pct, Headx, HCO3, Kser, Occl) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time / 30.4, value)) +
      geom_line(colour = PAL[2], linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "Months", y = NULL,
           title = "Safety: overdrainage · low-pressure headache · metabolic acidosis · occlusion") + THEME
  })

  output$tbl_safety <- renderTable({
    d <- sim()
    tibble(Item = c("Subdural collection [mL]", "Symptomatic subdural cumulative incidence [%]",
                    "Low-pressure headache index", "Serum HCO3 [mmol/L]",
                    "Serum K [mmol/L]", "Shunt occlusion fraction",
                    "Cumulative drainage volume [mL]"),
           `Final value` = c(fmt(last(d$Vsdh), 1), fmt(100 * last(d$SDH_inc), 1),
                     fmt(last(d$Headx)), fmt(last(d$HCO3), 1),
                     fmt(last(d$Kser)), fmt(last(d$Occl), 3),
                     fmt(last(d$Vdr), 0)))
  }, striped = TRUE)
}

shinyApp(ui, server)
