# =============================================================================
#  ami_shiny_app.R
#  Acute Myocardial Infarction (STEMI) — QSP interactive dashboard
# =============================================================================
#
#  A front end for ami_mrgsolve_model.R.  The app is organised around the
#  model's thesis rather than around its variable list, so each tab answers one
#  question:
#
#    1  Patient · lesion          Who is this patient, and what did the artery do?
#    2  Pharmacokinetics          What is in the blood, and when?
#    3  Necrosis wavefront        Race 1 — does the wavefront beat the clock?
#    4  Reperfusion injury        Race 2 — acid washout versus the oxidant burst
#    5  Microvascular obstruction  Is the tissue reperfused, or only the artery?
#    6  Inflammation · scar       Does a scar get built before the wall thins?
#    7  Remodelling bifurcation   The bifurcation — converge or diverge?
#    8  Clinical endpoints        What would we have measured?
#    9  Scenario comparison       Head-to-head strategy comparison
#   10  Biomarkers                Troponin, CK-MB, CRP, BNP, fibrinogen
#   11  Sensitivity · structure   Cut one edge; see what the model was standing on
#
#  RUN
#    install.packages(c("shiny","mrgsolve","dplyr","tidyr","ggplot2",
#                       "DT","bslib","bsicons","scales"))
#    shiny::runApp("ami_shiny_app.R")
#
#  The app never displays a number the model did not compute.  Where a quantity
#  is a reporting convenience rather than a simulated state (the haemorrhage
#  index on tab 8), it is labelled as such in the UI, not only in the code.
# =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)
library(bslib)
library(bsicons)
library(scales)

MODEL_FILE <- "ami_mrgsolve_model.R"
NOT_GIVEN  <- 1e6          # the model's "this intervention was withheld" value

mod <- mread_cache("ami", MODEL_FILE)

# -----------------------------------------------------------------------------
#  helpers
# -----------------------------------------------------------------------------
theme_ami <- function() {
  theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank(),
          legend.position = "bottom",
          plot.title = element_text(face = "bold", size = 14),
          plot.subtitle = element_text(colour = "grey35", size = 11))
}

# palette: cool for flow/treatment, warm for injury, green for protection
PAL <- c(L1 = "#8B0000", L2 = "#C0392B", L3 = "#E67E22",
         L4 = "#7D9F35", L5 = "#1E8449",
         isch = "#C0392B", rep = "#7B241C", mvo = "#B8860B",
         edv = "#1A5276", ef = "#0E6251", scar = "#145A32")

#' Turn UI inputs into a param() list.  Withheld interventions become 1e6
#' rather than NA, because the model reads them as times inside a tanh switch.
build_params <- function(input) {
  gv <- function(flag, value) if (isTRUE(flag)) value else NOT_GIVEN
  list(
    AAR      = input$aar,
    COLL     = input$coll,
    PRECOND  = as.numeric(input$precond),
    AGE      = input$age,
    SBP0     = input$sbp0,
    HR0      = input$hr0,
    T_OCC    = 0,
    T_PCI    = gv(input$do_pci,  input$t_pci  / 60),
    T_LYT    = gv(input$do_lyt,  input$t_lyt  / 60),
    TNK_DOSE = input$tnk_dose,
    T_ASA    = gv(input$do_asa,  input$t_anti / 60),
    T_TIC    = gv(input$do_tic,  input$t_anti / 60),
    T_HEP    = gv(input$do_hep,  input$t_anti / 60),
    T_MET_IV = gv(input$do_ivbb, input$t_ivbb / 60),
    T_CSA    = gv(input$do_csa,  input$t_csa  / 60),
    T_MET_PO = gv(input$do_bb,   input$t_gdmt),
    T_RAM    = gv(input$do_acei, input$t_gdmt),
    T_EPL    = gv(input$do_mra,  input$t_gdmt + 24),
    T_EMP    = gv(input$do_sglt, input$t_gdmt),
    T_COLC   = gv(input$do_colc, input$t_colc),
    T_CAN    = gv(input$do_can,  input$t_can),
    ARNI     = as.numeric(input$arni)
  )
}

#' Simulate.  `delta` is deliberately fine early: the mPTP transient on tab 4
#' lives inside the first thirty minutes of reflow and a coarse grid hides it.
simulate <- function(pars, end_h, over = list()) {
  p <- modifyList(pars, over)
  fine <- mod %>% param(p) %>% mrgsim(end = min(end_h, 12), delta = 0.005) %>% as_tibble()
  if (end_h <= 12) return(fine)
  coarse <- mod %>% param(p) %>% mrgsim(end = end_h, delta = 0.25) %>% as_tibble()
  bind_rows(fine, filter(coarse, time > 12)) %>% arrange(time)
}

at_time <- function(df, h) df %>% slice(which.min(abs(time - h)))

fmt <- function(x, d = 2) formatC(x, format = "f", digits = d)

# -----------------------------------------------------------------------------
#  UI
# -----------------------------------------------------------------------------
ui <- page_sidebar(
  title = "Acute myocardial infarction (STEMI) QSP — two races and one bifurcation",
  theme = bs_theme(version = 5, bootswatch = "flatly"),

  sidebar = sidebar(
    width = 360,
    accordion(
      open = c("Patient · lesion", "Reperfusion strategy"),

      accordion_panel(
        "Patient · lesion",
        sliderInput("aar", "Area at risk AAR (%LV)", 8, 55, 35, step = 1,
                    post = " %"),
        sliderInput("coll", "Collateral flow COLL (relative to normal flow)",
                    0.01, 0.35, 0.10, step = 0.01),
        helpText("The collateral flow sets the 'clock speed' of the wavefront. ",
                 "A layer that exceeds the basal metabolic demand (0.20) survives even when it becomes akinetic."),
        sliderInput("age", "Age (years)", 35, 92, 62, step = 1),
        sliderInput("sbp0", "Baseline systolic pressure (mmHg)", 90, 180, 120, step = 5),
        sliderInput("hr0", "Baseline heart rate (1/min)", 50, 110, 72, step = 2),
        checkboxInput("precond", "Pre-infarction angina (ischaemic preconditioning)", FALSE)
      ),

      accordion_panel(
        "Reperfusion strategy",
        checkboxInput("do_pci", "Primary/rescue PCI", TRUE),
        sliderInput("t_pci", "PCI time (minutes after symptom onset)", 20, 720, 90, step = 5),
        checkboxInput("do_lyt", "Fibrinolysis (tenecteplase)", FALSE),
        sliderInput("t_lyt", "Lytic bolus time (minutes)", 20, 360, 45, step = 5),
        sliderInput("tnk_dose", "Tenecteplase dose (mg)", 10, 50, 40, step = 5),
        helpText("Turn both strategies on at once and it becomes a pharmaco-invasive strategy (STREAM-type).")
      ),

      accordion_panel(
        "Antiplatelet · anticoagulation",
        sliderInput("t_anti", "Administration time (minutes)", 5, 360, 60, step = 5),
        checkboxInput("do_asa", "Aspirin", TRUE),
        checkboxInput("do_tic", "P2Y12 inhibitor (ticagrelor)", TRUE),
        checkboxInput("do_hep", "Anticoagulation (heparin)", TRUE)
      ),

      accordion_panel(
        "Cardioprotection (the window is minutes)",
        checkboxInput("do_ivbb", "IV metoprolol", FALSE),
        sliderInput("t_ivbb", "IV beta-blocker time (minutes)", 10, 360, 40, step = 5),
        checkboxInput("do_csa", "mPTP inhibitor (ciclosporin)", FALSE),
        sliderInput("t_csa", "CsA time (minutes)", 10, 360, 85, step = 5),
        helpText("CsA only means anything immediately before reperfusion — the model has no time-window ",
                 "parameter, because the pore opens within a few minutes.")
      ),

      accordion_panel(
        "Long-term drugs (GDMT)",
        sliderInput("t_gdmt", "GDMT start (hours)", 6, 168, 24, step = 6),
        checkboxInput("do_acei", "ACE inhibitor (ramipril)", TRUE),
        checkboxInput("arni", "Replace with ARNI (sacubitril/valsartan)", FALSE),
        checkboxInput("do_bb", "Oral beta-blocker", TRUE),
        checkboxInput("do_mra", "MR antagonist (eplerenone)", TRUE),
        checkboxInput("do_sglt", "SGLT2 inhibitor (empagliflozin)", TRUE),
        checkboxInput("do_colc", "Colchicine", FALSE),
        sliderInput("t_colc", "Colchicine start (hours)", 6, 336, 24, step = 6),
        checkboxInput("do_can", "Anti-IL-1β (canakinumab)", FALSE),
        sliderInput("t_can", "Anti-IL-1β time (hours)", 1, 336, 24, step = 1)
      ),

      accordion_panel(
        "Simulation",
        radioButtons("horizon", "Time horizon",
                     c("48 hours (acute)" = 48, "30 days" = 720,
                       "180 days" = 4320, "1 year" = 8760),
                     selected = 4320),
        actionButton("run", "Run", class = "btn-primary w-100")
      )
    )
  ),

  navset_card_tab(
    id = "tabs",

    # ---- 1 ------------------------------------------------------------------
    nav_panel(
      "1 · Patient · lesion",
      layout_columns(
        col_widths = c(6, 6),
        card(card_header("Per-layer collateral flow versus basal metabolic demand"),
             plotOutput("p_supply", height = 330),
             card_footer("The horizontal line is the basal demand. Only the layers below it die — ",
                         "this one panel is the whole survival verdict of the model.")),
        card(card_header("Epicardial patency (PAT)"), plotOutput("p_patency", height = 330))
      ),
      card(card_header("Summary of this run"), DTOutput("t_summary"))
    ),

    # ---- 2 ------------------------------------------------------------------
    nav_panel(
      "2 · Pharmacokinetics",
      layout_columns(
        col_widths = c(6, 6),
        card(card_header("Thrombolysis · antiplatelet"), plotOutput("p_pk_acute", height = 330)),
        card(card_header("Long-term drug effect occupancy"), plotOutput("p_pk_chronic", height = 330))
      ),
      card(card_header("The fibrinolytic system — plasmin · plasminogen · fibrinogen"),
           plotOutput("p_lysis", height = 300),
           card_footer("The fibrinogen nadir is the surrogate for the systemic bleeding cost."))
    ),

    # ---- 3 ------------------------------------------------------------------
    nav_panel(
      "3 · Necrosis wavefront (Race 1)",
      layout_columns(
        col_widths = c(7, 5),
        card(card_header("Per-layer necrosis — from endocardium to epicardium"),
             plotOutput("p_wavefront", height = 380),
             card_footer("The only difference between the layers is the collateral flow. ",
                         "Nowhere was the model told the order of progression.")),
        card(card_header("Energy charge E"), plotOutput("p_energy", height = 380))
      ),
      layout_columns(
        col_widths = c(6, 6),
        card(card_header("Infarct size and salvage index"), plotOutput("p_is", height = 300)),
        card(card_header("Reperfusion time sweep — the time–muscle curve"),
             plotOutput("p_salvage_curve", height = 300),
             card_footer("Everything else in the current settings held fixed, changing only the PCI time."))
      )
    ),

    # ---- 4 ------------------------------------------------------------------
    nav_panel(
      "4 · Reperfusion injury (Race 2)",
      card(card_header("The first minutes after reperfusion — acid (H⁺) · succinate · ROS · Ca²⁺ · mPTP"),
           plotOutput("p_reperf", height = 400),
           card_footer("As the acid is washed out the gate of the pore opens, and at that same moment ",
                       "succinate is burned and the oxidant bursts. ",
                       "This is the pH paradox, and there is no time-window parameter.")),
      layout_columns(
        col_widths = c(6, 6),
        card(card_header("Ischaemic versus reperfusion necrosis"), plotOutput("p_split", height = 320)),
        card(card_header("The acid gate and pore opening"), plotOutput("p_gate", height = 320))
      )
    ),

    # ---- 5 ------------------------------------------------------------------
    nav_panel(
      "5 · Microvascular obstruction",
      layout_columns(
        col_widths = c(7, 5),
        card(card_header("MVO and effective layer flow — an open artery, a closed tissue"),
             plotOutput("p_mvo", height = 380),
             card_footer("MVO lowers the flow, and the lowered flow grows the MVO. ",
                         "No-reflow is positive feedback.")),
        card(card_header("MVO versus infarct size"), plotOutput("p_mvo_is", height = 380))
      )
    ),

    # ---- 6 ------------------------------------------------------------------
    nav_panel(
      "6 · Inflammation · scar",
      layout_columns(
        col_widths = c(6, 6),
        card(card_header("Cells — neutrophils · M1 · M2"), plotOutput("p_cells", height = 330)),
        card(card_header("Cytokines — IL-1β · IL-6 · TGF-β · CRP"),
             plotOutput("p_cyto", height = 330))
      ),
      layout_columns(
        col_widths = c(6, 6),
        card(card_header("Scar strength versus MMP activity — the race"),
             plotOutput("p_scar", height = 330),
             card_footer("If the scar is late the wall thins. ",
                         "Suppress inflammation too hard and you lose this race.")),
        card(card_header("Infarct wall thinning and thickness"), plotOutput("p_thin", height = 330))
      )
    ),

    # ---- 7 ------------------------------------------------------------------
    nav_panel(
      "7 · Remodelling bifurcation",
      card(card_header("The Laplace loop — wall stress · volume · mass"),
           plotOutput("p_remodel", height = 380),
           card_footer("Dilatation is positive feedback, concentric hypertrophy negative feedback. ",
                       "Convergence or divergence is a bifurcation, not a dose.")),
      layout_columns(
        col_widths = c(6, 6),
        card(card_header("EDV · ESV · EF trajectories"), plotOutput("p_volumes", height = 330)),
        card(card_header("Area at risk sweep — where is the separatrix?"),
             plotOutput("p_bifurc", height = 330),
             card_footer("Changing only the AAR and computing the 1-year EDV growth rate. ",
                         "There is no critical-infarct-size parameter in the model."))
      )
    ),

    # ---- 8 ------------------------------------------------------------------
    nav_panel(
      "8 · Clinical endpoints",
      layout_columns(
        col_widths = c(4, 4, 4),
        value_box(title = "Infarct size (%LV)", value = textOutput("vb_is"),
                  showcase = bsicons::bs_icon("heart-pulse")),
        value_box(title = "Salvage index (%)", value = textOutput("vb_salv"),
                  showcase = bsicons::bs_icon("shield-check")),
        value_box(title = "MVO (%LV)", value = textOutput("vb_mvo"),
                  showcase = bsicons::bs_icon("droplet"))
      ),
      layout_columns(
        col_widths = c(4, 4, 4),
        value_box(title = "EF, final time point (%)", value = textOutput("vb_ef"),
                  showcase = bsicons::bs_icon("activity")),
        value_box(title = "EDV, final time point (mL)", value = textOutput("vb_edv"),
                  showcase = bsicons::bs_icon("arrows-expand")),
        value_box(title = "Scar strength (0–1)", value = textOutput("vb_scar"),
                  showcase = bsicons::bs_icon("bandaid"))
      ),
      card(card_header("Haemodynamics and neurohormones"), plotOutput("p_neuro", height = 330)),
      card(card_header("Reporting indices (not a validated risk model)"),
           DTOutput("t_report"),
           card_footer(strong("Caution: "), "the intracranial haemorrhage index is a reporting figure built from plasmin exposure and ",
                       "age alone and must never be read as an absolute risk. ",
                       "The model does not compute event rates."))
    ),

    # ---- 9 ------------------------------------------------------------------
    nav_panel(
      "9 · Scenario comparison",
      card(card_header("Strategy comparison — the same patient, a different decision"),
           checkboxGroupInput(
             "scen", NULL, inline = TRUE,
             choices = c("No reperfusion" = "none",
                         "PCI 90 min" = "pci90",
                         "PCI 240 min" = "pci240",
                         "Prehospital lysis 45 min" = "lyt45",
                         "Lysis 45 min + PCI 4 h" = "pharminv",
                         "PCI 90 min + full GDMT" = "pci90gdmt",
                         "PCI 90 min + GDMT + cardioprotection" = "full"),
             selected = c("none", "pci90", "pci90gdmt", "full")),
           plotOutput("p_scen", height = 420)),
      card(card_header("Scenario endpoint table"), DTOutput("t_scen"))
    ),

    # ---- 10 -----------------------------------------------------------------
    nav_panel(
      "10 · Biomarkers",
      card(card_header("cTnI and CK-MB — the washout phenomenon"),
           plotOutput("p_markers", height = 360),
           card_footer("A reperfused artery gives an early, high-peak curve; an occluded artery gives ",
                       "a later, lower curve from a larger real infarct. ",
                       "That is why measuring an infarct by enzyme needs the area under the curve.")),
      layout_columns(
        col_widths = c(6, 6),
        card(card_header("Inflammation · wall-stress markers"), plotOutput("p_bio2", height = 320)),
        card(card_header("Peak versus AUC versus real infarct size"), DTOutput("t_markers"))
      )
    ),

    # ---- 11 -----------------------------------------------------------------
    nav_panel(
      "11 · Sensitivity · structure",
      card(card_header("Cutting one mechanism at a time"),
           helpText("Each row sends a single parameter to zero or to infinity and re-runs the same ",
                    "scenario. It shows what the model was ",
                    "standing on."),
           DTOutput("t_sens")),
      card(card_header("The effect of the structural cuts"), plotOutput("p_sens", height = 380))
    )
  )
)

# -----------------------------------------------------------------------------
#  server
# -----------------------------------------------------------------------------
server <- function(input, output, session) {

  pars <- eventReactive(input$run, build_params(input), ignoreNULL = FALSE)
  horizon <- eventReactive(input$run, as.numeric(input$horizon), ignoreNULL = FALSE)

  sim <- reactive({
    p <- pars(); p$AAR <- p$AAR / 100
    simulate(p, horizon())
  })

  acute <- reactive(filter(sim(), time <= 48))
  last_row <- reactive(sim() %>% slice_tail(n = 1))

  # ---- 1 --------------------------------------------------------------------
  output$p_supply <- renderPlot({
    p <- pars()
    gc <- c(0.05, 0.35, 0.75, 1.35, 2.50)
    d <- tibble(layer = factor(paste0("L", 1:5), levels = paste0("L", 1:5)),
                collateral = p$COLL * gc)
    ggplot(d, aes(layer, collateral, fill = layer)) +
      geom_col(width = .65) +
      geom_hline(yintercept = 0.20, linetype = "dashed",
                 colour = "#8B0000", linewidth = 1) +
      annotate("text", x = 4.2, y = 0.215, label = "Basal metabolic demand 0.20",
               colour = "#8B0000", size = 4) +
      scale_fill_manual(values = PAL[1:5], guide = "none") +
      labs(x = NULL, y = "Collateral flow (relative to normal)",
           title = "The layers that die and the layers that live",
           subtitle = "Below the dashed line = basal demand not met = the path to death") +
      theme_ami()
  })

  output$p_patency <- renderPlot({
    ggplot(acute(), aes(time * 60, PATENCY)) +
      geom_line(colour = PAL[["edv"]], linewidth = 1.1) +
      scale_x_continuous(limits = c(0, 480)) +
      labs(x = "After symptom onset (min)", y = "Epicardial patency (%)",
           title = "When did the artery open") +
      theme_ami()
  })

  output$t_summary <- renderDT({
    a <- at_time(sim(), 48); z <- last_row()
    tibble(
      Item = c("Area at risk (%LV)", "Collateral flow", "Infarct size 48 h (%LV)",
               "Salvage index (%)", "Ischaemic component (%LV)", "Reperfusion component (%LV)",
               "Reperfusion injury fraction (%)", "MVO (%LV)",
               "EF 48 h (%)", "EF final (%)", "EDV final (mL)",
               "Scar strength", "cTnI peak"),
      Value = c(fmt(input$aar, 0), fmt(pars()$COLL, 3), fmt(a$IS), fmt(a$SALV, 1),
             fmt(a$IS_ISCH), fmt(a$IS_REP), fmt(a$RPFRAC, 1), fmt(a$MVO_LV),
             fmt(a$EF, 1), fmt(z$EF, 1), fmt(z$EDV, 1), fmt(z$SCAR, 3),
             fmt(max(sim()$CTNI), 1))
    ) %>% datatable(rownames = FALSE, options = list(dom = "t", pageLength = 20))
  })

  # ---- 2 --------------------------------------------------------------------
  output$p_pk_acute <- renderPlot({
    acute() %>%
      transmute(time, `Tenecteplase (mg/L)` = TNK1 / 4.2,
                `Ticagrelor (mg/L)` = TICc / 88,
                `Metoprolol (mg/L)` = METc / 250) %>%
      pivot_longer(-time) %>% filter(value > 1e-9) %>%
      ggplot(aes(time, value, colour = name)) +
      geom_line(linewidth = 1) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      scale_colour_manual(values = unname(PAL[c(1, 4, 5)]), guide = "none") +
      labs(x = "Time (h)", y = "Blood concentration") + theme_ami()
  })

  output$p_pk_chronic <- renderPlot({
    sim() %>%
      transmute(time = time / 24,
                ACEi = RAMc, MRA = EPL, SGLT2i = EMP,
                `Colchicine` = COLC, `Anti-IL-1β` = CAN) %>%
      pivot_longer(-time) %>% group_by(name) %>%
      filter(max(value) > 1e-9) %>% mutate(value = value / max(value)) %>%
      ungroup() %>%
      ggplot(aes(time, value, colour = name)) +
      geom_line(linewidth = 1) +
      labs(x = "Time (days)", y = "Effect compartment (relative to maximum)", colour = NULL,
           title = "Long-term drug exposure") + theme_ami()
  })

  output$p_lysis <- renderPlot({
    acute() %>% filter(time <= 12) %>%
      transmute(time, `Plasmin` = PLN, `Plasminogen` = PLG,
                `Fibrinogen` = FIBX, `PAI-1` = PAI) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
      geom_line(linewidth = 1) +
      labs(x = "Time (h)", y = "Relative to normal / activity", colour = NULL) + theme_ami()
  })

  # ---- 3 --------------------------------------------------------------------
  output$p_wavefront <- renderPlot({
    acute() %>%
      select(time, NECRO1:NECRO5) %>%
      pivot_longer(-time, names_to = "layer") %>%
      mutate(layer = factor(sub("NECRO", "L", layer),
                            levels = paste0("L", 1:5))) %>%
      ggplot(aes(time, value, colour = layer)) +
      geom_line(linewidth = 1.1) +
      scale_colour_manual(values = PAL[1:5],
                          labels = c("L1 subendocardial", "L2", "L3 mid-wall",
                                     "L4", "L5 subepicardial")) +
      labs(x = "Time (h)", y = "Necrosis within layer (%)", colour = NULL,
           title = "The necrosis wavefront", subtitle = "The only difference between layers is the collateral flow") +
      theme_ami()
  })

  output$p_energy <- renderPlot({
    acute() %>% filter(time <= 12) %>%
      select(time, E1, E2, E3, E4, E5) %>%
      pivot_longer(-time, names_to = "layer") %>%
      mutate(layer = factor(sub("E", "L", layer), levels = paste0("L", 1:5))) %>%
      ggplot(aes(time, value, colour = layer)) +
      geom_line(linewidth = 1.1) +
      scale_colour_manual(values = PAL[1:5], guide = "none") +
      labs(x = "Time (h)", y = "Energy charge E",
           title = "Energy falls only where the basal demand cannot be met") + theme_ami()
  })

  output$p_is <- renderPlot({
    acute() %>%
      ggplot(aes(time)) +
      geom_line(aes(y = IS, colour = "Infarct size (%LV)"), linewidth = 1.2) +
      geom_line(aes(y = SALV * input$aar / 100, colour = "Salvage (%LV)"),
                linewidth = 1) +
      scale_colour_manual(values = c(unname(PAL[["isch"]]), unname(PAL[["scar"]]))) +
      labs(x = "Time (h)", y = "%LV", colour = NULL) + theme_ami()
  })

  salvage_curve <- reactive({
    p <- pars(); p$AAR <- p$AAR / 100
    tp <- c(15, 30, 45, 60, 90, 120, 180, 240, 360, 480)
    lapply(tp, function(m) {
      s <- simulate(p, 48, over = list(T_PCI = m / 60, T_ASA = 0.05,
                                       T_TIC = 0.05, T_HEP = 0.05,
                                       T_LYT = NOT_GIVEN)) %>% slice_tail(n = 1)
      tibble(delay = m, IS = s$IS, SALV = s$SALV, RI = s$RPFRAC, MVO = s$MVO_LV)
    }) %>% bind_rows()
  })

  output$p_salvage_curve <- renderPlot({
    salvage_curve() %>%
      ggplot(aes(delay, IS)) +
      geom_line(colour = PAL[["isch"]], linewidth = 1.2) +
      geom_point(size = 2.4, colour = PAL[["isch"]]) +
      labs(x = "Reperfusion delay (min)", y = "Final infarct size (%LV)",
           title = "The time–muscle curve is not a straight line") + theme_ami()
  })

  # ---- 4 --------------------------------------------------------------------
  output$p_reperf <- renderPlot({
    trep <- if (isTRUE(input$do_pci)) input$t_pci / 60 else input$t_lyt / 60
    acute() %>%
      filter(time >= max(trep - 0.25, 0), time <= trep + 2) %>%
      transmute(min = (time - trep) * 60,
                `H⁺ (L1)` = H1, `Succinate (L1)` = SU1, `ROS` = ROS,
                `Ca²⁺ (L1)` = C1, `mPTP (L1)` = P1 * 10) %>%
      pivot_longer(-min) %>%
      ggplot(aes(min, value, colour = name)) +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
      geom_line(linewidth = 1.1) +
      annotate("text", x = 2, y = Inf, vjust = 1.5, hjust = 0,
               label = "Reperfusion", colour = "grey30") +
      labs(x = "Relative to reperfusion (min)", y = "Value (mPTP is ×10)", colour = NULL,
           title = "The first few minutes after reperfusion",
           subtitle = "The moment the acid is washed out, the gate of the pore opens") +
      theme_ami()
  })

  output$p_split <- renderPlot({
    acute() %>% select(time, `Ischaemic` = IS_ISCH, `Reperfusion` = IS_REP) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, fill = name)) +
      geom_area(alpha = .85) +
      scale_fill_manual(values = unname(PAL[c("isch", "rep")])) +
      labs(x = "Time (h)", y = "%LV", fill = NULL,
           title = "The decomposition into two mechanisms — a result, not an assumption") + theme_ami()
  })

  output$p_gate <- renderPlot({
    acute() %>% filter(time <= 12) %>%
      transmute(time, `Acid gate 1/(1+H/0.3)` = 1 / (1 + H1 / 0.30),
                `mPTP opening (L1)` = P1) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
      geom_line(linewidth = 1.1) +
      scale_colour_manual(values = c(unname(PAL[["scar"]]), unname(PAL[["rep"]]))) +
      labs(x = "Time (h)", y = "0–1", colour = NULL) + theme_ami()
  })

  # ---- 5 --------------------------------------------------------------------
  output$p_mvo <- renderPlot({
    acute() %>% transmute(time, `MVO (%LV)` = MVO_LV,
                          `Epicardial patency (%)` = PATENCY) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
      geom_line(linewidth = 1.1) +
      scale_colour_manual(values = c(unname(PAL[["edv"]]), unname(PAL[["mvo"]]))) +
      labs(x = "Time (h)", y = NULL, colour = NULL,
           title = "The artery is open. And the tissue?") + theme_ami()
  })

  output$p_mvo_is <- renderPlot({
    salvage_curve() %>%
      ggplot(aes(IS, MVO, colour = delay)) +
      geom_path(linewidth = 1) + geom_point(size = 3) +
      scale_colour_viridis_c(option = "plasma") +
      labs(x = "Infarct size (%LV)", y = "MVO (%LV)", colour = "Delay (min)",
           title = "MVO is not a function of infarct size") + theme_ami()
  })

  # ---- 6 --------------------------------------------------------------------
  output$p_cells <- renderPlot({
    sim() %>% filter(time <= 24 * 21) %>%
      transmute(day = time / 24, `Neutrophils` = NEU, `M1` = M1, `M2` = M2) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value, colour = name)) + geom_line(linewidth = 1.1) +
      labs(x = "Days", y = "Relative level", colour = NULL,
           title = "From inflammation to repair") + theme_ami()
  })

  output$p_cyto <- renderPlot({
    sim() %>% filter(time <= 24 * 21) %>%
      transmute(day = time / 24, `IL-1β` = IL1, `IL-6` = IL6,
                `TGF-β` = TGF, `CRP` = CRPX) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value, colour = name)) + geom_line(linewidth = 1.1) +
      labs(x = "Days", y = "Relative level", colour = NULL) + theme_ami()
  })

  output$p_scar <- renderPlot({
    sim() %>% filter(time <= 24 * 60) %>%
      transmute(day = time / 24, `Scar strength` = SCAR, `MMP activity` = MMP) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value, colour = name)) + geom_line(linewidth = 1.2) +
      scale_colour_manual(values = c(unname(PAL[["mvo"]]), unname(PAL[["scar"]]))) +
      labs(x = "Days", y = NULL, colour = NULL,
           title = "The rate at which the scar is built versus the rate at which it is torn down") + theme_ami()
  })

  output$p_thin <- renderPlot({
    sim() %>% transmute(day = time / 24, `Thinning (%)` = THINPCT,
                        `Infarct thickness (cm)` = H_INF * 10) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value, colour = name)) + geom_line(linewidth = 1.1) +
      labs(x = "Days", y = "Value (thickness in mm)", colour = NULL) + theme_ami()
  })

  # ---- 7 --------------------------------------------------------------------
  output$p_remodel <- renderPlot({
    sim() %>% transmute(day = time / 24,
                        `Wall stress (remote)` = WS / first(WS),
                        `Wall stress (infarct)` = WS_INF / first(WS),
                        `EDV / EDV₀` = EDV / first(EDV),
                        `Mass / mass₀` = MASS / first(MASS)) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value, colour = name)) +
      geom_hline(yintercept = 1, linetype = "dotted") +
      geom_line(linewidth = 1.1) +
      labs(x = "Days", y = "Relative to baseline", colour = NULL,
           title = "The Laplace loop",
           subtitle = "Dilatation raises the stress, hypertrophy lowers it") + theme_ami()
  })

  output$p_volumes <- renderPlot({
    sim() %>% transmute(day = time / 24, EDV, ESV, EF) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value, colour = name)) + geom_line(linewidth = 1.1) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      scale_colour_manual(values = unname(PAL[c("edv", "ef", "isch")]),
                          guide = "none") +
      labs(x = "Days", y = NULL) + theme_ami()
  })

  bifurc <- reactive({
    p <- pars()
    aars <- c(0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50)
    lapply(aars, function(a) {
      s <- simulate(p, 24 * 365, over = list(AAR = a))
      e30 <- at_time(s, 24 * 30)$EDV
      z <- slice_tail(s, n = 1)
      tibble(AAR = a * 100, IS = z$IS,
             growth = 100 * (z$EDV - e30) / e30, EF = z$EF)
    }) %>% bind_rows()
  })

  output$p_bifurc <- renderPlot({
    bifurc() %>%
      ggplot(aes(IS, growth)) +
      geom_hline(yintercept = c(4, 12), linetype = c("dotted", "dashed"),
                 colour = c("grey50", "#8B0000")) +
      geom_line(colour = PAL[["edv"]], linewidth = 1.2) +
      geom_point(size = 2.6, colour = PAL[["edv"]]) +
      annotate("text", x = Inf, y = 12, hjust = 1.1, vjust = -0.5,
               label = "Divergence line", colour = "#8B0000", size = 3.6) +
      labs(x = "Infarct size (%LV)", y = "30-day→1-year EDV growth (%)",
           title = "The separatrix is a computed result") + theme_ami()
  })

  # ---- 8 --------------------------------------------------------------------
  output$vb_is   <- renderText(fmt(at_time(sim(), 48)$IS))
  output$vb_salv <- renderText(fmt(at_time(sim(), 48)$SALV, 1))
  output$vb_mvo  <- renderText(fmt(max(sim()$MVO_LV)))
  output$vb_ef   <- renderText(fmt(last_row()$EF, 1))
  output$vb_edv  <- renderText(fmt(last_row()$EDV, 1))
  output$vb_scar <- renderText(fmt(last_row()$SCAR, 3))

  output$p_neuro <- renderPlot({
    sim() %>% transmute(day = time / 24, NE, `Ang II` = ANG, `Aldosterone` = ALD,
                        `NP (BNP)` = NTPROBNP, `Volume` = VOL) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value, colour = name)) + geom_line(linewidth = 1.1) +
      geom_hline(yintercept = 1, linetype = "dotted") +
      labs(x = "Days", y = "Relative to normal", colour = NULL) + theme_ami()
  })

  output$t_report <- renderDT({
    s <- sim()
    pln_auc <- sum(diff(s$time) * head(s$PLN, -1))
    ich <- 0.9 * pln_auc * (1 + 0.055 * (input$age - 50))
    tibble(
      Metric = c("Plasmin exposure (AUC)", "Fibrinogen nadir (relative to normal)",
               "Intracranial haemorrhage index (reporting only)", "cTnI AUC", "CRP peak"),
      Value = c(fmt(pln_auc, 4), fmt(min(s$FIBX), 3), fmt(ich, 3),
             fmt(sum(diff(s$time) * head(s$CTNI, -1)), 0), fmt(max(s$CRPX), 2))
    ) %>% datatable(rownames = FALSE, options = list(dom = "t"))
  })

  # ---- 9 --------------------------------------------------------------------
  SCEN <- list(
    none      = list(lab = "No reperfusion",
                     ov = list(T_PCI = NOT_GIVEN, T_LYT = NOT_GIVEN)),
    pci90     = list(lab = "PCI 90 min",
                     ov = list(T_PCI = 1.5, T_LYT = NOT_GIVEN, T_ASA = 1,
                               T_TIC = 1, T_HEP = 1)),
    pci240    = list(lab = "PCI 240 min",
                     ov = list(T_PCI = 4, T_LYT = NOT_GIVEN, T_ASA = 3.5,
                               T_TIC = 3.5, T_HEP = 3.5)),
    lyt45     = list(lab = "Prehospital lysis 45 min",
                     ov = list(T_PCI = NOT_GIVEN, T_LYT = 0.75, T_ASA = 0.65,
                               T_TIC = 0.65, T_HEP = 0.65)),
    pharminv  = list(lab = "Lysis 45 min + PCI 4 h",
                     ov = list(T_PCI = 4, T_LYT = 0.75, T_ASA = 0.65,
                               T_TIC = 0.65, T_HEP = 0.65)),
    pci90gdmt = list(lab = "PCI 90 min + full GDMT",
                     ov = list(T_PCI = 1.5, T_LYT = NOT_GIVEN, T_ASA = 1,
                               T_TIC = 1, T_HEP = 1, T_RAM = 24,
                               T_MET_PO = 24, T_EPL = 48, T_EMP = 24)),
    full      = list(lab = "PCI 90 min + GDMT + cardioprotection",
                     ov = list(T_PCI = 1.5, T_LYT = NOT_GIVEN, T_ASA = 1,
                               T_TIC = 1, T_HEP = 1, T_MET_IV = 0.6,
                               T_CSA = 1.45, T_COLC = 24, T_RAM = 24,
                               T_MET_PO = 24, T_EPL = 48, T_EMP = 24))
  )

  scen_runs <- reactive({
    req(length(input$scen) > 0)
    p <- pars(); p$AAR <- p$AAR / 100
    lapply(input$scen, function(k) {
      simulate(p, horizon(), over = SCEN[[k]]$ov) %>%
        mutate(scenario = SCEN[[k]]$lab)
    }) %>% bind_rows()
  })

  output$p_scen <- renderPlot({
    scen_runs() %>%
      transmute(day = time / 24, scenario,
                `Infarct size (%LV)` = IS, `MVO (%LV)` = MVO_LV,
                `EF (%)` = EF, `EDV (mL)` = EDV) %>%
      pivot_longer(-c(day, scenario)) %>%
      ggplot(aes(day, value, colour = scenario)) +
      geom_line(linewidth = 1) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "Days", y = NULL, colour = NULL) + theme_ami()
  })

  output$t_scen <- renderDT({
    scen_runs() %>% group_by(scenario) %>%
      summarise(`Infarct size (%LV)` = round(IS[which.min(abs(time - 48))], 2),
                `Salvage index (%)` = round(SALV[which.min(abs(time - 48))], 1),
                `MVO (%LV)` = round(max(MVO_LV), 2),
                `cTnI peak` = round(max(CTNI), 1),
                `EF final (%)` = round(last(EF), 1),
                `EDV final (mL)` = round(last(EDV), 1),
                `Scar strength` = round(last(SCAR), 3),
                .groups = "drop") %>%
      datatable(rownames = FALSE, options = list(dom = "t"))
  })

  # ---- 10 -----------------------------------------------------------------
  output$p_markers <- renderPlot({
    p <- pars(); p$AAR <- p$AAR / 100
    variants <- list(
      list(lab = "PCI 60 min", ov = list(T_PCI = 1, T_LYT = NOT_GIVEN)),
      list(lab = "PCI 240 min", ov = list(T_PCI = 4, T_LYT = NOT_GIVEN)),
      list(lab = "No reperfusion", ov = list(T_PCI = NOT_GIVEN, T_LYT = NOT_GIVEN))
    )
    lapply(variants, function(v)
      simulate(p, 120, over = v$ov) %>% mutate(scenario = v$lab)) %>%
      bind_rows() %>%
      transmute(time, scenario, `cTnI` = CTNI, `CK-MB` = CKMBP) %>%
      pivot_longer(-c(time, scenario)) %>%
      ggplot(aes(time, value, colour = scenario)) +
      geom_line(linewidth = 1.1) + facet_wrap(~name, scales = "free_y") +
      labs(x = "Time (h)", y = "Blood concentration (arbitrary units)", colour = NULL,
           title = "The washout phenomenon",
           subtitle = "Reperfused artery = early, high peak. Occluded artery = later and lower, but the infarct is larger") +
      theme_ami()
  })

  output$p_bio2 <- renderPlot({
    sim() %>% transmute(day = time / 24, `CRP` = CRPX,
                        `NP (BNP)` = NTPROBNP, `Fibrinogen` = FIBX) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value, colour = name)) + geom_line(linewidth = 1.1) +
      labs(x = "Days", y = "Relative to normal", colour = NULL) + theme_ami()
  })

  output$t_markers <- renderDT({
    p <- pars(); p$AAR <- p$AAR / 100
    rows <- list(
      list(lab = "PCI 60 min", ov = list(T_PCI = 1, T_LYT = NOT_GIVEN)),
      list(lab = "PCI 120 min", ov = list(T_PCI = 2, T_LYT = NOT_GIVEN)),
      list(lab = "PCI 240 min", ov = list(T_PCI = 4, T_LYT = NOT_GIVEN)),
      list(lab = "No reperfusion", ov = list(T_PCI = NOT_GIVEN, T_LYT = NOT_GIVEN))
    )
    lapply(rows, function(v) {
      s <- simulate(p, 120, over = v$ov)
      auc <- sum(diff(s$time) * head(s$CTNI, -1))
      tibble(Scenario = v$lab,
             `TnI peak` = round(max(s$CTNI), 1),
             `Peak time (h)` = round(s$time[which.max(s$CTNI)], 2),
             `TnI AUC` = round(auc, 0),
             `Real infarct (%LV)` = round(last(s$IS), 2),
             `AUC / infarct` = round(auc / max(last(s$IS), 1e-6), 0))
    }) %>% bind_rows() %>%
      datatable(rownames = FALSE, options = list(dom = "t"))
  })

  # ---- 11 -----------------------------------------------------------------
  CUTS <- list(
    list(lab = "(intact model)", ov = list()),
    list(lab = "No mPTP (KP_ON = 0)", ov = list(KP_ON = 0)),
    list(lab = "No reperfusion Ca surge (KCA_RP = 0)", ov = list(KCA_RP = 0)),
    list(lab = "Acid gate removed (HP50 = 1e6)", ov = list(HP50 = 1e6)),
    list(lab = "No succinate (KSUC = 0)", ov = list(KSUC = 0)),
    list(lab = "No oxidant burst (KROS = 0)", ov = list(KROS = 0)),
    list(lab = "Collateral gradient flattened", ov = list(GC1 = 1, GC2 = 1, GC3 = 1,
                                                  GC4 = 1, GC5 = 1)),
    list(lab = "No collateral flow (COLL = 0)", ov = list(COLL = 0)),
    list(lab = "No MVO feedback", ov = list(GM1 = 0, GM2 = 0, GM3 = 0,
                                          GM4 = 0, GM5 = 0)),
    list(lab = "No dilatation loop (KDIL = 0)", ov = list(KDIL = 0)),
    list(lab = "No hypertrophy brake (KHYP = 0)", ov = list(KHYP = 0)),
    list(lab = "No anaerobic reserve (F_ANAERO = 0)", ov = list(F_ANAERO = 0)),
    list(lab = "Preconditioned host (PRECOND = 1)", ov = list(PRECOND = 1))
  )

  sens <- reactive({
    p <- pars(); p$AAR <- p$AAR / 100
    lapply(CUTS, function(c) {
      s <- simulate(p, max(horizon(), 24 * 180), over = c$ov)
      a <- at_time(s, 48); z <- slice_tail(s, n = 1)
      tibble(perturbation = c$lab, IS = a$IS, IS_REP = a$IS_REP,
             MVO = max(s$MVO_LV), EDV = z$EDV, EF = z$EF)
    }) %>% bind_rows()
  })

  output$t_sens <- renderDT({
    sens() %>%
      transmute(`Cut` = perturbation,
                `Infarct size (%LV)` = round(IS, 2),
                `Reperfusion component (%LV)` = round(IS_REP, 2),
                `MVO (%LV)` = round(MVO, 2),
                `EDV 180 days (mL)` = round(EDV, 1),
                `EF 180 days (%)` = round(EF, 1)) %>%
      datatable(rownames = FALSE, options = list(dom = "t", pageLength = 15))
  })

  output$p_sens <- renderPlot({
    base <- sens() %>% filter(perturbation == "(intact model)")
    sens() %>% filter(perturbation != "(intact model)") %>%
      mutate(dIS = IS - base$IS, dEDV = EDV - base$EDV) %>%
      pivot_longer(c(dIS, dEDV)) %>%
      mutate(name = recode(name, dIS = "Change in infarct size (%LV)",
                           dEDV = "Change in 180-day EDV (mL)")) %>%
      ggplot(aes(reorder(perturbation, value), value, fill = value > 0)) +
      geom_col() + coord_flip() +
      facet_wrap(~name, scales = "free_x") +
      scale_fill_manual(values = c("TRUE" = "#C0392B", "FALSE" = "#1E8449"),
                        guide = "none") +
      labs(x = NULL, y = "Change relative to the intact model",
           title = "What was the model standing on") + theme_ami()
  })
}

shinyApp(ui, server)
