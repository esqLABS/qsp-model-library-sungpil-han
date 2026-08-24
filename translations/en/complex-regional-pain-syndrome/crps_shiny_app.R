## =============================================================================
## CRPS QSP Shiny dashboard
## Complex Regional Pain Syndrome QSP Model Interactive Dashboard
##
##   library(shiny); library(mrgsolve); library(ggplot2); library(dplyr)
##   shiny::runApp("crps_shiny_app.R")
##
## The app loads crps_mrgsolve_model.R from the same directory and exposes the
## model's structure rather than only its output: the left panel separates
## PATIENT covariates (trait, sympathetic tone, insult, cast) from TREATMENT
## choices (arms, start day, intensity), because in this model the interaction
## between those two blocks -- not the potency of any arm -- is what decides
## the outcome.
##
## Tabs (8):
##   1. Patient profile          patient profile, attractor read-out, ring gain
##   2. Drug PK                all drug concentrations on one time axis
##   3. Peripheral node (PD)         SP/CGRP, cytokines, NGF, ROS, oedema, alpha1
##   4. Central node · Glial latch spinal, glial latch, descending tone, cortex
##   5. Clinical endpoints        NRS, CSS, ROM, temperature asymmetry
##   6. Treatment window (window)       start-day scan -> t*, with the arm decomposition
##   7. Scenario comparison          9 prebuilt scenarios side by side
##   8. Biomarkers · Bone        CTX-I, BMD, blister-fluid-like cytokine panel
## =============================================================================

library(shiny)
library(mrgsolve)
library(ggplot2)
library(dplyr)
library(tidyr)

MODEL_FILE <- "crps_mrgsolve_model.R"
mod <- mrgsolve::mread_cache("crps", MODEL_FILE)
E   <- mod@envir            # event builders + analysis functions
DAY <- 24

STATE_LABELS <- c(
  NP = "SP/CGRP neuropeptides", CYT = "Cytokines (IL-6/TNF/IL-1b)",
  NGF = "NGF", EDEMA = "Oedema", ROS = "Oxidative stress (ROS)",
  AAB = "Autoantibody titre", ALPHA1 = "alpha1-receptor upregulation",
  PSENS = "Peripheral sensitisation", PERF = "Perfusion index", HYPOX = "Tissue hypoxia/acidosis",
  SSENS = "Spinal central sensitisation", GLIA = "Glial activation (latch variable)",
  DINH = "Descending inhibitory tone", CORTEX = "Cortical map degradation", DISUSE = "Disuse/avoidance",
  PAIN = "Pain NRS", ROM = "Joint range of motion", OC = "Osteoclast activity",
  BMD = "Local BMD", CTXI = "Serum CTX-I")

theme_crps <- function() {
  theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold"),
          legend.position = "bottom")
}

## ---------------------------------------------------------------------------
## UI
## ---------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Complex Regional Pain Syndrome (CRPS) QSP Simulator — Peripheral Node · Behaviour-Cortex Loop · Glial Latch"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Patient covariates"),
      sliderInput("kfear", "Fear-avoidance gain, KFEAR (trait)",
                  min = 0.2, max = 1.6, value = 1.0, step = 0.05),
      sliderInput("symp", "Sympathetic tone, SYMP_TONE",
                  min = 0.5, max = 2.0, value = 1.0, step = 0.1),
      sliderInput("inj", "Inciting-event intensity, INJ_AMP",
                  min = 0.1, max = 1.6, value = 1.0, step = 0.05),
      sliderInput("immob", "Cast immobilisation duration (days)",
                  min = 0, max = 90, value = 30, step = 5),
      hr(),
      h4("Treatment"),
      sliderInput("start", "Treatment start day (day)",
                  min = 3, max = 365, value = 7, step = 1),
      checkboxGroupInput("arms", "Treatment arm",
                         choices = c("Prednisolone (oral taper)" = "pred",
                                     "NAC / DMSO (antioxidant)" = "nac",
                                     "Rehabilitation · GMI · Graded exposure" = "rehab",
                                     "Neridronate 100mg IV x4" = "nerid",
                                     "Ketamine, 100-hour infusion" = "ket",
                                     "Gabapentin 600mg TID" = "gbp",
                                     "Amitriptyline 50mg" = "amt",
                                     "IVIG 0.5 g/kg x6" = "ivig",
                                     "Spinal cord stimulation (SCS)" = "scs",
                                     "Sympathetic block" = "sb",
                                     "Vasodilator (PDE5i-like)" = "vaso"),
                         selected = c("pred", "nac", "rehab")),
      sliderInput("rehab_int", "Rehabilitation intensity/adherence (0-1)",
                  min = 0, max = 1, value = 1, step = 0.05),
      sliderInput("ketrate", "Ketamine max infusion rate (mg/h)",
                  min = 5, max = 90, value = 22, step = 1),
      sliderInput("years", "Simulation period (years)",
                  min = 1, max = 5, value = 3, step = 1),
      hr(),
      helpText("Turn on treatment, then move only the 'treatment start day'. The same regimen ",
               "shows complete remission before day 90 and no effect after day 95 — ",
               "that is the central result of this model.")
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("1. Patient Profile",
                 fluidRow(column(4, wellPanel(h4("Status at 3 years"), tableOutput("attractor"))),
                          column(8, plotOutput("profile_pain", height = "300px"))),
                 plotOutput("profile_ring", height = "260px")),
        tabPanel("2. Drug PK", plotOutput("pk", height = "620px")),
        tabPanel("3. Peripheral Node (PD)", plotOutput("periph", height = "620px")),
        tabPanel("4. Central Node · Glial Latch",
                 plotOutput("central", height = "420px"),
                 plotOutput("latch", height = "250px")),
        tabPanel("5. Clinical Endpoints", plotOutput("endpoints", height = "620px")),
        tabPanel("6. Treatment Window",
                 fluidRow(column(6, plotOutput("window", height = "380px")),
                          column(6, plotOutput("armscan", height = "380px"))),
                 tableOutput("window_tab")),
        tabPanel("7. Scenario Comparison",
                 plotOutput("scen_plot", height = "420px"),
                 tableOutput("scen_tab")),
        tabPanel("8. Biomarkers · Bone", plotOutput("bio", height = "620px"))
      )
    )
  )
)

## ---------------------------------------------------------------------------
## server
## ---------------------------------------------------------------------------
server <- function(input, output, session) {

  build_events <- function(arms, start_day, ketrate) {
    st <- start_day * DAY
    evs <- list()
    if ("pred"  %in% arms) evs <- c(evs, list(E$ev_prednisolone(st)))
    if ("nac"   %in% arms) evs <- c(evs, list(E$ev_nac(st, days = 180)))
    if ("nerid" %in% arms) evs <- c(evs, list(E$ev_neridronate(st)))
    if ("ket"   %in% arms) evs <- c(evs, list(E$ev_ketamine(st, rate_max = ketrate)))
    if ("gbp"   %in% arms) evs <- c(evs, list(E$ev_gabapentin(st, days = 180)))
    if ("amt"   %in% arms) evs <- c(evs, list(E$ev_amitriptyline(st, days = 180)))
    if ("ivig"  %in% arms) evs <- c(evs, list(E$ev_ivig(st)))
    if (!length(evs)) return(NULL)
    do.call(E$comb_ev, evs)
  }

  build_pars <- function(input) {
    st <- input$start * DAY
    p <- list(KFEAR = input$kfear, SYMP_TONE = input$symp,
              INJ_AMP = input$inj, IMMOB_DUR = input$immob * DAY)
    if ("rehab" %in% input$arms)
      p <- c(p, list(REHAB = input$rehab_int, REHAB_T0 = st,
                     REHAB_DUR = 180 * DAY))
    if ("scs"  %in% input$arms) p <- c(p, list(SCS_ON = 1, SCS_T0 = st))
    if ("sb"   %in% input$arms) p <- c(p, list(SYMPBLOCK = 1, SB_T0 = st))
    if ("vaso" %in% input$arms) p <- c(p, list(VASODIL = 1))
    p
  }

  simulate <- reactive({
    end <- input$years * 8760
    E$sim(mod, build_pars(input),
          build_events(input$arms, input$start, input$ketrate),
          end = end, delta = 6)
  })

  untreated <- reactive({
    end <- input$years * 8760
    E$sim(mod, list(KFEAR = input$kfear, SYMP_TONE = input$symp,
                    INJ_AMP = input$inj, IMMOB_DUR = input$immob * DAY),
          NULL, end = end, delta = 6)
  })

  long <- function(d, cols) {
    d %>% select(all_of(c("time", cols))) %>%
      mutate(day = time / 24) %>% select(-time) %>%
      pivot_longer(-day, names_to = "state", values_to = "value") %>%
      mutate(label = ifelse(state %in% names(STATE_LABELS),
                            paste0(state, " — ", STATE_LABELS[state]), state))
  }

  ## ---- tab 1 -------------------------------------------------------------
  output$attractor <- renderTable({
    f <- E$summarise_run(simulate())
    u <- E$summarise_run(untreated())
    tab <- data.frame(
      metric   = c("NRS", "CSS (0-16)", "ROM", "BMD", "GLIA", "CORTEX",
                   "DISUSE", "TEMP_ASYM (C)", "latched", "remission"),
      treated  = c(f$NRS, f$CSS, f$ROM, f$BMD, f$GLIA, f$CORTEX,
                   f$DISUSE, f$TEMP, f$LATCHED, f$REMISSION),
      untreated = c(u$NRS, u$CSS, u$ROM, u$BMD, u$GLIA, u$CORTEX,
                    u$DISUSE, u$TEMP, u$LATCHED, u$REMISSION))
    names(tab) <- c("Metric", "Treated", "Untreated")
    tab
  }, digits = 3)

  output$profile_pain <- renderPlot({
    d <- bind_rows(mutate(simulate(), arm = "Treated"),
                   mutate(untreated(), arm = "Untreated"))
    ggplot(d, aes(time / 24, PAIN, colour = arm)) +
      geom_line(linewidth = 0.9) +
      geom_vline(xintercept = input$start, linetype = 2, colour = "grey40") +
      labs(x = "day", y = "Pain NRS (0-10)", colour = NULL,
           title = "Pain course — dashed line = treatment start") +
      scale_colour_manual(values = c("Treated" = "#1565c0", "Untreated" = "#b71c1c")) +
      theme_crps()
  })

  output$profile_ring <- renderPlot({
    d <- simulate()
    ggplot(d, aes(time / 24, RING_GAIN)) +
      geom_hline(yintercept = 1, linetype = 2, colour = "#b71c1c") +
      geom_line(colour = "#6a1b9a", linewidth = 0.9) +
      labs(x = "day", y = "Ring gain",
           title = "Local gain of the pain→disuse→cortex→pain loop (>1 means self-sustaining is possible)") +
      theme_crps()
  })

  ## ---- tab 2 -------------------------------------------------------------
  output$pk <- renderPlot({
    d <- simulate() %>% mutate(day = time / 24) %>%
      select(day, `Ketamine (ng/mL)` = KET_ng, `Norketamine (ng/mL)` = NORK_ng,
             `Prednisolone (ng/mL)` = PRED_ng, `Gabapentin (ng/mL)` = GBP_ng,
             `Amitriptyline (ng/mL)` = AMT_ng, `NAC (ng/mL)` = NAC_ng,
             `IVIG (g/L)` = IVIG_gL, `Neridronate plasma (mg/L)` = NER_mgL) %>%
      pivot_longer(-day)
    ggplot(d, aes(day, value)) + geom_line(colour = "#1565c0") +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "day", y = NULL, title = "Drug concentration-time profile") + theme_crps()
  })

  ## ---- tab 3 -------------------------------------------------------------
  output$periph <- renderPlot({
    ggplot(long(simulate(), c("NP", "CYT", "NGF", "EDEMA", "ROS", "AAB",
                              "ALPHA1", "PSENS", "HYPOX")),
           aes(day, value)) +
      geom_line(colour = "#c62828") +
      facet_wrap(~label, scales = "free_y", ncol = 3) +
      labs(x = "day", y = "Activity index (0-1)",
           title = "Fast peripheral node — driven by the inciting event and self-extinguishing") +
      theme_crps()
  })

  ## ---- tab 4 -------------------------------------------------------------
  output$central <- renderPlot({
    ggplot(long(simulate(), c("SSENS", "GLIA", "DINH", "CORTEX", "DISUSE")),
           aes(day, value)) +
      geom_line(colour = "#4527a0") +
      facet_wrap(~label, scales = "free_y", ncol = 3) +
      labs(x = "day", y = "Index", title = "Slow central loop (weeks-to-months time constant)") +
      theme_crps()
  })

  output$latch <- renderPlot({
    g50 <- as.numeric(mrgsolve::param(mod)$GLIA50)
    d <- bind_rows(mutate(simulate(), arm = "Treated"),
                   mutate(untreated(), arm = "Untreated"))
    ggplot(d, aes(time / 24, GLIA, colour = arm)) +
      geom_hline(yintercept = g50, linetype = 2, colour = "#b71c1c") +
      geom_line(linewidth = 0.9) +
      annotate("text", x = Inf, y = g50, label = " GLIA50 (latch threshold)",
               hjust = 1.02, vjust = -0.6, size = 3.4, colour = "#b71c1c") +
      scale_colour_manual(values = c("Treated" = "#1565c0", "Untreated" = "#b71c1c")) +
      labs(x = "day", y = "GLIA", colour = NULL,
           title = "Glial latch variable: once this line is crossed, spinal sensitisation is sustained without afferent input") +
      theme_crps()
  })

  ## ---- tab 5 -------------------------------------------------------------
  output$endpoints <- renderPlot({
    d <- simulate() %>% mutate(day = time / 24) %>%
      select(day, `Pain NRS` = PAIN, `CRPS severity score (0-16)` = CSS,
             `Joint range of motion (ratio to normal)` = ROM,
             `Temperature asymmetry (°C, +warm/-cold)` = TEMP_ASYM,
             `Oedema index` = EDEMA, `Active CRPS flag` = ACTIVE_CRPS) %>%
      pivot_longer(-day)
    ggplot(d, aes(day, value)) + geom_line(colour = "#00695c") +
      geom_vline(xintercept = input$start, linetype = 2, colour = "grey50") +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "day", y = NULL, title = "Clinical endpoints") + theme_crps()
  })

  ## ---- tab 6 -------------------------------------------------------------
  window_res <- reactive({
    days <- c(3, 7, 14, 21, 30, 45, 60, 75, 90, 95, 105, 120, 180, 240, 365)
    base <- build_pars(input)
    arms <- intersect(input$arms, c("pred", "nac", "rehab"))
    if (!length(arms)) arms <- c("pred", "nac", "rehab")
    out <- lapply(days, function(dd) {
      st <- dd * DAY
      p <- base
      p$REHAB_T0 <- st
      if ("rehab" %in% arms) p$REHAB <- input$rehab_int else p$REHAB <- 0
      evs <- list()
      if ("pred" %in% arms) evs <- c(evs, list(E$ev_prednisolone(st)))
      if ("nac"  %in% arms) evs <- c(evs, list(E$ev_nac(st, days = 180)))
      ev <- if (length(evs)) do.call(E$comb_ev, evs) else NULL
      d <- E$sim(mod, p, ev, end = 3 * 8760, delta = 24)
      f <- E$final_row(d)
      data.frame(start_day = dd, NRS = f$PAIN, CSS = f$CSS,
                 GLIA = f$GLIA, BMD = f$BMD, ROM = f$ROM)
    })
    bind_rows(out)
  })

  output$window <- renderPlot({
    w <- window_res()
    ggplot(w, aes(start_day, NRS)) +
      geom_line(colour = "#2e7d32", linewidth = 1) +
      geom_point(colour = "#2e7d32") +
      labs(x = "Treatment start day (day)", y = "NRS at 3 years",
           title = "Treatment window: same regimen, only the start day differs") +
      theme_crps()
  })

  output$armscan <- renderPlot({
    days <- c(7, 30, 60, 120, 240)
    combos <- list(pred = "pred", nac = "nac", rehab = "rehab",
                   all = c("pred", "nac", "rehab"))
    res <- lapply(names(combos), function(nm) {
      sapply(days, function(dd) {
        p <- E$package(dd, combos[[nm]])
        d <- E$sim(mod, utils::modifyList(build_pars(input), p$pars),
                   p$events, end = 3 * 8760, delta = 24)
        E$final_row(d)$PAIN
      })
    })
    df <- data.frame(day = rep(days, length(combos)),
                     arm = rep(names(combos), each = length(days)),
                     NRS = unlist(res))
    ggplot(df, aes(day, NRS, colour = arm)) +
      geom_line(linewidth = 0.9) + geom_point() +
      labs(x = "Treatment start day", y = "NRS at 3 years", colour = "arm",
           title = "Which arm creates the window") + theme_crps()
  })

  output$window_tab <- renderTable(window_res(), digits = 3)

  ## ---- tab 7 -------------------------------------------------------------
  scen <- reactive(E$run_scenarios(mod, end = 3 * 8760))

  output$scen_plot <- renderPlot({
    s <- scen()
    d <- bind_rows(lapply(names(s), function(n)
      data.frame(day = s[[n]]$time / 24, NRS = s[[n]]$PAIN, scenario = n)))
    ggplot(d, aes(day, NRS, colour = scenario)) + geom_line(linewidth = 0.8) +
      labs(x = "day", y = "NRS", colour = NULL, title = "9 scenarios") +
      theme_crps()
  })

  output$scen_tab <- renderTable(attr(scen(), "summary"), digits = 3)

  ## ---- tab 8 -------------------------------------------------------------
  output$bio <- renderPlot({
    d <- simulate() %>% mutate(day = time / 24) %>%
      select(day, `Serum CTX-I (ng/mL)` = CTXI, `Local BMD (ratio to contralateral)` = BMD,
             `Osteoclast activity` = OC, `Cytokine index (blister-fluid analogue)` = CYT,
             `Autoantibody titre` = AAB, `Oxidative stress` = ROS) %>%
      pivot_longer(-day)
    ggplot(d, aes(day, value)) + geom_line(colour = "#6d4c41") +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "day", y = NULL, title = "Biomarkers and bone axis") + theme_crps()
  })
}

shinyApp(ui, server)
