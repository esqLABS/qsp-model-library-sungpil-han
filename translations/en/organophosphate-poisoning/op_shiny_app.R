## =============================================================================
##  op_shiny_app.R
##  Acute organophosphorus insecticide self-poisoning — QSP dashboard
##  Acute organophosphorus insecticide poisoning — QSP dashboard
##
##  Run:
##      library(shiny); library(mrgsolve); library(dplyr); library(tidyr)
##      library(ggplot2); library(DT)
##      shiny::runApp("op_shiny_app.R")
##
##  The app is organised around ONE question: for this patient, at this moment,
##  is the oxime capable of doing anything at all?  That question has a number
##  attached to it —
##
##       OMEGA = k_r_max / (k_i * C_oxon),  free-enzyme ceiling = OMEGA/(1+OMEGA)
##
##  — and every tab is either an input to that number, a consequence of it, or
##  a reminder that the ventilator does not care about it.
## =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

MODEL_FILE <- "op_mrgsolve_model.R"
mod <- mread_cache(file = MODEL_FILE, model = "op")

## ---------------------------------------------------------------------------
## Compound and oxime libraries: they live in the model file's $ENV block, so
## the app and the model can never drift apart.
## ---------------------------------------------------------------------------
MODENV <- tryCatch(mrgsolve::env_get(mod),
                   error = function(e) as.list(mod@envir))
OP_LIBRARY    <- MODENV$OP_LIBRARY
OXIME_LIBRARY <- MODENV$OXIME_LIBRARY
stopifnot(!is.null(OP_LIBRARY), !is.null(OXIME_LIBRARY))

op_params <- function(compound = "chlorpyrifos", oxime = NULL, volume_ml = 50,
                      pon1 = "RR", wt = 60) {
  L <- OP_LIBRARY[[compound]]
  p <- L[c("MW_OP", "KI", "KI_NMJ", "KI_NTE", "T12AGE", "T12SPO", "T12AGNTE",
           "KA_GUT", "FBIO", "V1TH", "VFTH", "QFTH", "VMAXBIO", "KMBIO",
           "CLTHOTH", "CLOXON", "VOXON", "QOXON", "VTOXON", "SOLVFR", "SOLVCV")]
  p$CONC_GL <- L$CONC_GL
  p$DOSE_ML <- volume_ml
  p$PON1SC  <- if (identical(pon1, "QQ")) L$PON1_QQ else 1
  p$WT      <- wt
  if (is.null(oxime) || identical(oxime, "none")) {
    p$OXTYPE <- 0; p$OXRATE <- 0
  } else {
    O <- OXIME_LIBRARY[[oxime]]
    k <- if (identical(O$oxime, "obi")) L$obi else L$pam
    p$KRMAX <- k$KRMAX; p$KDOX <- k$KDOX
    p$OXV1 <- O$V1; p$OXV2 <- O$V2; p$OXCL <- O$CL; p$OXQ <- O$Q
    p$OXMW <- O$MW
    p$OXTYPE <- 1
    p$OXRATE <- if (!is.null(O$inf_mg_kg_h)) O$inf_mg_kg_h * wt * 1000 / O$MW
                else if (!is.null(O$inf_mg_24h)) O$inf_mg_24h / 24 * 1000 / O$MW
                else 0
  }
  p
}

op_events <- function(oxime = NULL, wt = 60, diazepam = TRUE, magnesium = FALSE,
                      scavenger_mg = 0, oxime_start = 1.5, oxime_dur = 168) {
  e <- list()
  if (!is.null(oxime) && !identical(oxime, "none")) {
    O <- OXIME_LIBRARY[[oxime]]
    if (!is.null(O$load_mg_kg))
      e <- c(e, list(ev(time = oxime_start, amt = O$load_mg_kg * wt * 1000 / O$MW,
                        cmt = "A_pam_c")))
    if (!is.null(O$load_mg))
      e <- c(e, list(ev(time = oxime_start, amt = O$load_mg * 1000 / O$MW,
                        cmt = "A_pam_c")))
    if (!is.null(O$bolus_mg))
      e <- c(e, list(ev(time = oxime_start, amt = O$bolus_mg * 1000 / O$MW,
                        cmt = "A_pam_c", ii = O$interval_h,
                        addl = floor(oxime_dur / O$interval_h))))
  }
  if (diazepam)
    e <- c(e, list(ev(time = 0.6, amt = 10 * 1000 / 284.7, cmt = "A_dz",
                      ii = 6, addl = 2)))
  if (magnesium)
    e <- c(e, list(ev(time = 1.0, amt = 16.2, cmt = "A_mg")))
  if (scavenger_mg > 0)
    e <- c(e, list(ev(time = 1.0, amt = scavenger_mg / 1000 / 340000 * 1e6 * 4,
                      cmt = "SCAV")))
  if (!length(e)) return(ev(time = 0, amt = 0, cmt = "A_dz"))
  Reduce(`+`, e)
}

run_case <- function(compound, oxime, volume_ml, pon1, wt, oxime_start,
                     oxime_dur, atr_mode, ventilator, oxygen, icu,
                     glyco, magnesium, scavenger_mg, charcoal_t, tmax = 336) {
  p <- op_params(compound, oxime, volume_ml, pon1, wt)
  p$OXSTART <- oxime_start
  p$OXDUR   <- oxime_dur
  p$ATRMODE <- atr_mode
  p$VENTAVL <- as.numeric(ventilator)
  p$O2AVL   <- as.numeric(oxygen)
  p$ICUSUCT <- as.numeric(icu)
  p$BBBF    <- if (isTRUE(glyco)) 0.02 else 0.35
  p$CHARC_T <- if (is.null(charcoal_t) || charcoal_t < 0) -1 else charcoal_t
  e <- op_events(oxime, wt, TRUE, magnesium, scavenger_mg, oxime_start, oxime_dur)
  mod %>% param(p) %>% ev(e) %>%
    mrgsim(end = tmax, delta = 0.25) %>% as_tibble()
}

## ---------------------------------------------------------------------------
## Closed-form helpers — no ODE needed, these ARE the model's arithmetic
## ---------------------------------------------------------------------------
omega_fn   <- function(krmax, ki, C) krmax / (ki * C)
ceiling_fn <- function(krmax, ki, C) { w <- omega_fn(krmax, ki, C); w / (1 + w) }
achieved_fn <- function(krmax, KD, ki, ks, C, X) {
  kr <- krmax * X / (KD + X); (kr + ks) / (kr + ks + ki * C)
}
Ccrit_fn   <- function(krmax, ki, target = 0.30) (1 - target) / target * krmax / ki
Xreq_fn    <- function(krmax, KD, ki, ks, C, target = 0.30) {
  need <- target / (1 - target) * ki * C - ks
  if (need <= 0) return(0)
  if (need >= krmax) return(NA_real_)      # arithmetically impossible
  KD * need / (krmax - need)
}

THEME <- theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 13),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")

## ===========================================================================
## UI
## ===========================================================================
ui <- fluidPage(
  titlePanel("Acute OP Insecticide Self-Poisoning — QSP Dashboard"),
  tags$p(style = "color:#666;margin-top:-8px;",
         "51-state mrgsolve model · the oxime sufficiency number Ω = k_r_max/(k_i·C_oxon) decides everything the oxime can do"),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("① Exposure"),
      selectInput("compound", "Compound",
                  choices = c("Chlorpyrifos (diethyl)"   = "chlorpyrifos",
                              "Parathion (diethyl)"      = "parathion",
                              "Dimethoate (dimethyl)"    = "dimethoate",
                              "Fenthion (dimethyl, lipophilic)" = "fenthion"),
                  selected = "chlorpyrifos"),
      sliderInput("volume", "Volume ingested (mL of concentrate)", 1, 250, 50, step = 1),
      sliderInput("wt", "Body weight (kg)", 35, 100, 60, step = 1),
      radioButtons("pon1", "PON1 Q192R genotype", inline = TRUE,
                   choices = c("R192R (fast)" = "RR", "Q192Q (slow)" = "QQ")),

      hr(), h4("② Antidotes"),
      selectInput("oxime", "Oxime",
                  choices = c("None"                              = "none",
                              "Pralidoxime, WHO infusion"         = "pam_who",
                              "Pralidoxime 1 g q6h bolus"         = "pam_bolus",
                              "Obidoxime 250 mg + 750 mg/24h"     = "obidoxime"),
                  selected = "pam_who"),
      sliderInput("oxstart", "Oxime start time (h)", 0.25, 24, 1.5, step = 0.25),
      sliderInput("oxdur", "Oxime treatment duration (h)", 12, 336, 168, step = 12),
      radioButtons("atr", "Atropine titration", inline = FALSE,
                   choices = c("Rapid doubling protocol"                 = 1,
                               "Conventional titration (slow, ad hoc)"   = 2,
                               "None"                                   = 0),
                   selected = 1),
      checkboxInput("glyco", "Substitute glycopyrrolate (does not cross the BBB)", FALSE),
      checkboxInput("mg", "Magnesium sulphate 4 g", FALSE),
      sliderInput("scav", "Plasma BChE bioscavenger (mg)", 0, 5000, 0, step = 250),
      sliderInput("charcoal", "Time of activated charcoal (h; -1 = not given)", -1, 12, -1, step = 0.5),

      hr(), h4("③ Supportive care"),
      checkboxInput("vent", "Ventilator available", TRUE),
      checkboxInput("o2", "Oxygen available", TRUE),
      checkboxInput("icu", "ICU airway care (suctioning)", TRUE),

      hr(),
      actionButton("go", "Simulate", class = "btn-primary",
                   width = "100%"),
      tags$p(style = "font-size:11px;color:#888;margin-top:10px;",
             "For education and research only. Do not use for clinical decision-making.")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",

        ## ------------------------------------------------------------------
        tabPanel("① Patient profile",
          br(),
          fluidRow(
            column(3, wellPanel(h5("Moles ingested"), textOutput("kpi_dose"))),
            column(3, wellPanel(h5("Steady-state oxon"), textOutput("kpi_oxon"))),
            column(3, wellPanel(h5("Ω (24 h)"), textOutput("kpi_omega"))),
            column(3, wellPanel(h5("Oxime ceiling"), textOutput("kpi_ceiling")))
          ),
          fluidRow(
            column(3, wellPanel(h5("AChE nadir"), textOutput("kpi_nadir"))),
            column(3, wellPanel(h5("Cumulative atropine"), textOutput("kpi_atr"))),
            column(3, wellPanel(h5("Hours of ventilation"), textOutput("kpi_vent"))),
            column(3, wellPanel(h5("14-day mortality"), textOutput("kpi_mort")))
          ),
          hr(),
          h4("What the oxime can do in this patient"),
          verbatimTextOutput("verdict"),
          hr(),
          plotOutput("plot_profile", height = "420px")
        ),

        ## ------------------------------------------------------------------
        tabPanel("② Toxicokinetics (PK/TK)",
          br(),
          plotOutput("plot_tk", height = "700px"),
          hr(),
          tags$p("The parent compound (the thion) does not inhibit AChE. What inhibits it is the oxon made by CYP, ",
                 "and its formation saturates. So after a large ingestion the oxon plateau concentration converges on ",
                 "Vmax/CL_oxon and becomes independent of dose — which is why the oxime ceiling is a constant after a large ingestion.")
        ),

        ## ------------------------------------------------------------------
        tabPanel("③ The esterase switch",
          br(),
          plotOutput("plot_esterase", height = "760px"),
          hr(),
          h5("Three-state mass balance"),
          tags$pre("E  --k_i*C_oxon-->  EP  --k_a-->  EP_aged   (irreversible)\nE  <--k_s+k_r[X]--  EP\n\nThe ageing loss is not an event but an integral:  aged(T) = 1 - exp(-k_a * ∫ f_EP dt)")
        ),

        ## ------------------------------------------------------------------
        tabPanel("④ The cholinergic syndrome",
          br(),
          plotOutput("plot_signs", height = "760px"),
          hr(),
          tags$p("The atropine requirement grows as 1/E, not as (1-E). That is because the steady state of synaptic ACh is ",
                 "the hyperbola (1+leak)/(E+leak), and because the antagonism is competitive the concentration required is ",
                 "linear in ACh. That a patient at 2% AChE needs several times more atropine than one at 20% is not a ",
                 "matter of 'degree' of severity but of arithmetic.")
        ),

        ## ------------------------------------------------------------------
        tabPanel("⑤ Respiration · clinical endpoints",
          br(),
          plotOutput("plot_resp", height = "760px"),
          hr(),
          tags$p("Ventilation is the minimum of what the brain demands and what the muscles can deliver. ",
                 "Resting ventilation needs only about 28% of maximum inspiratory force (the neuromuscular safety factor), ",
                 "so neuromuscular block has far more margin than central depression. This asymmetry is the model's ",
                 "expression of the clinical observation that 'early deaths are central'.")
        ),

        ## ------------------------------------------------------------------
        tabPanel("⑥ Scenario comparison",
          br(),
          fluidRow(
            column(6, checkboxGroupInput("scen", "Scenarios to compare",
              choices = c("Supportive care only"                 = "sup",
                          "2-PAM WHO infusion"                   = "pam_who",
                          "2-PAM 1 g q6h bolus"                  = "pam_bolus",
                          "Obidoxime"                            = "obidoxime",
                          "No ventilator (2-PAM)"                = "novent",
                          "Slow atropine titration (2-PAM)"      = "slowatr",
                          "Glycopyrrolate (2-PAM)"               = "glyco",
                          "Activated charcoal at 1 h (2-PAM)"    = "char1",
                          "PON1 Q192Q (2-PAM)"                   = "ponqq"),
              selected = c("sup", "pam_who", "novent", "slowatr"))),
            column(6, radioButtons("scen_y", "Variable to display",
              choices = c("RBC AChE (%)"        = "AChE_RBC",
                          "Cumulative mortality (%)" = "MORTALITY",
                          "Ventilatory capacity VCAP" = "VCAP",
                          "Cumulative atropine (mg)" = "ATRCUM",
                          "POP severity score"   = "POP",
                          "Oxon (nM)"            = "C_oxon"),
              selected = "AChE_RBC"))
          ),
          plotOutput("plot_scen", height = "440px"),
          hr(),
          DTOutput("tbl_scen")
        ),

        ## ------------------------------------------------------------------
        tabPanel("⑦ Oxime ceiling explorer (Ω)",
          br(),
          fluidRow(
            column(4, sliderInput("ex_oxon", "Plasma oxon (nM, log axis)",
                                  min = 0, max = 4, value = 2.1, step = 0.05,
                                  pre = "10^")),
            column(4, sliderInput("ex_target", "Target free-enzyme fraction", 0.05, 0.60,
                                  0.30, step = 0.05)),
            column(4, sliderInput("ex_X", "Oxime concentration (µM)", 0, 400, 130, step = 5))
          ),
          plotOutput("plot_omega", height = "420px"),
          hr(),
          DTOutput("tbl_omega"),
          hr(),
          tags$p(strong("How to read this: "),
            "Where the curve cannot rise above the target line, the problem is not that the oxime is 'underdosed' ",
            "but that the target is unreachable with an oxime as a class of drug. No randomised trial can show a ",
            "benefit in patients in that region.")
        ),

        ## ------------------------------------------------------------------
        tabPanel("⑧ Dose-response sweep",
          br(),
          actionButton("go_sweep", "Run the sweep (8 doses × 2 arms)", class = "btn-info"),
          br(), br(),
          plotOutput("plot_sweep", height = "440px"),
          DTOutput("tbl_sweep")
        ),

        ## ------------------------------------------------------------------
        tabPanel("⑨ Virtual clinical trial",
          br(),
          fluidRow(
            column(3, numericInput("n_trial", "Number of patients", 40, min = 10, max = 200, step = 10)),
            column(3, sliderInput("trial_med", "Median volume ingested (mL)", 5, 80, 28, step = 1)),
            column(3, sliderInput("trial_cpf", "Proportion of diethyl OP", 0, 1, 0.55, step = 0.05)),
            column(3, br(), actionButton("go_trial", "Run the trial", class = "btn-warning"))
          ),
          plotOutput("plot_trial", height = "400px"),
          hr(),
          verbatimTextOutput("txt_trial"),
          tags$p("A real-world randomised trial enrols the patients who arrive, as they arrive. If the patients above ",
                 "the ceiling dilute the real effect in those below it, a negative result follows without ",
                 "assuming that the oxime is powerless.")
        ),

        ## ------------------------------------------------------------------
        tabPanel("⑩ Biomarkers · delayed syndromes",
          br(),
          plotOutput("plot_biomarker", height = "700px"),
          hr(),
          tags$p("Erythrocyte AChE is not resynthesised (it is replaced only as the red cells turn over). Muscle and brain AChE ",
                 "are resynthesised with a half-life of about 5 days. That RBC AChE and clinical state diverge during recovery ",
                 "is not measurement error; the two compartments simply turn over at different rates.")
        ),

        ## ------------------------------------------------------------------
        tabPanel("⑪ Documentation",
          br(),
          uiOutput("docs")
        )
      )
    )
  )
)

## ---------------------------------------------------------------------------
## small helper so the app does not hard-depend on patchwork/gridExtra
## ---------------------------------------------------------------------------
gridExtra_arrange <- function(...) {
  ps <- list(...)
  if (requireNamespace("patchwork", quietly = TRUE)) {
    patchwork::wrap_plots(ps, ncol = 1)
  } else if (requireNamespace("gridExtra", quietly = TRUE)) {
    gridExtra::grid.arrange(grobs = ps, ncol = 1)
  } else {
    ps[[1]]
  }
}

## ===========================================================================
## SERVER
## ===========================================================================
server <- function(input, output, session) {

  sim <- eventReactive(input$go, ignoreNULL = FALSE, {
    withProgress(message = "Integrating the 51-state ODE...", value = 0.5, {
      run_case(input$compound,
               if (identical(input$oxime, "none")) NULL else input$oxime,
               input$volume, input$pon1, input$wt, input$oxstart, input$oxdur,
               as.numeric(input$atr), input$vent, input$o2, input$icu,
               input$glyco, input$mg, input$scav, input$charcoal)
    })
  })

  cmp <- reactive(OP_LIBRARY[[input$compound]])

  ## ---------------------------------------------------------------- KPIs
  output$kpi_dose <- renderText({
    L <- cmp()
    mmol <- input$volume * L$CONC_GL / 1000 / L$MW_OP * 1000
    sprintf("%.1f mmol (%.1f g)", mmol, input$volume * L$CONC_GL / 1000)
  })
  output$kpi_oxon <- renderText({
    d <- sim(); sprintf("%.0f nM (24 h)", approx(d$time, d$C_oxon, 24)$y)
  })
  output$kpi_omega <- renderText({
    d <- sim(); sprintf("%.2f", approx(d$time, d$OMEGA, 24)$y)
  })
  output$kpi_ceiling <- renderText({
    d <- sim(); sprintf("%.0f %% free enzyme", 100 * approx(d$time, d$E_CEIL, 24)$y)
  })
  output$kpi_nadir <- renderText({ sprintf("%.2f %%", min(sim()$AChE_RBC)) })
  output$kpi_atr   <- renderText({ sprintf("%.0f mg / 14 days", max(sim()$ATRCUM)) })
  output$kpi_vent  <- renderText({ sprintf("%.0f h", max(sim()$VTIME)) })
  output$kpi_mort  <- renderText({ sprintf("%.1f %%", max(sim()$MORTALITY)) })

  output$verdict <- renderText({
    d <- sim(); L <- cmp()
    C  <- approx(d$time, d$C_oxon, 24)$y
    krmax <- if (identical(input$oxime, "obidoxime")) L$obi$KRMAX else L$pam$KRMAX
    KD    <- if (identical(input$oxime, "obidoxime")) L$obi$KDOX  else L$pam$KDOX
    ks <- log(2) / L$T12SPO
    ceil <- ceiling_fn(krmax, L$KI, C)
    Xr30 <- Xreq_fn(krmax, KD, L$KI, ks, C, 0.30)
    Xclin <- if (identical(input$oxime, "obidoxime")) 20 else 130
    ach  <- achieved_fn(krmax, KD, L$KI, ks, C, Xclin)
    phi  <- (log(2) / L$T12AGE) /
      (log(2) / L$T12AGE + ks + krmax * Xclin / (KD + Xclin))
    paste0(
      sprintf("Oxon at 24 h %.0f nM,  k_i·C_oxon = %.1f /h,  k_r_max = %.0f /h\n",
              C, L$KI * C, krmax),
      sprintf("Ω = %.2f  →  maximum free enzyme reachable at any dose = %.0f%%\n",
              omega_fn(krmax, L$KI, C), 100 * ceil),
      sprintf("At the clinically reachable concentration (%.0f µM) the actual free enzyme = %.0f%%\n",
              Xclin, 100 * ach),
      if (is.na(Xr30))
        "→ 30% free enzyme is arithmetically impossible with this oxime (it is below the k_r_max ceiling).\n"
      else if (Xr30 > Xclin)
        sprintf("→ Reaching 30%% would need %.0f µM, but the clinically reachable concentration is %.0f µM.\n",
                Xr30, Xclin)
      else
        sprintf("→ The concentration needed to reach 30%% is %.0f µM, which the current regimen can reach.\n",
                Xr30),
      sprintf("Ageing ratchet φ = %.3f (probability that a single inhibition event ends irreversibly)\n", phi),
      sprintf("Predicted 14-day mortality %.1f%%,  mechanical ventilation %.0f h,  cumulative atropine %.0f mg",
              max(d$MORTALITY), max(d$VTIME), max(d$ATRCUM))
    )
  })

  ## -------------------------------------------------------------- ① profile
  output$plot_profile <- renderPlot({
    d <- sim()
    d %>% select(time, `RBC AChE (%)` = AChE_RBC, `Ω` = OMEGA,
                 `Ventilatory capacity VCAP` = VCAP, `Mortality (%)` = MORTALITY) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) +
      geom_line(linewidth = 1.0, colour = "#1f6fb2") +
      facet_wrap(~name, scales = "free_y") +
      scale_x_continuous(breaks = c(0, 24, 72, 168, 336)) +
      labs(x = "Time (h)", y = NULL,
           title = "Core trajectories") + THEME
  })

  ## ------------------------------------------------------------------ ② TK
  output$plot_tk <- renderPlot({
    d <- sim()
    bind_rows(
      tibble(time = d$time, value = d$C_thion,  panel = "Parent thion (µM)"),
      tibble(time = d$time, value = d$C_oxon,   panel = "Plasma oxon (nM)"),
      tibble(time = d$time, value = d$C_oxonT,  panel = "Tissue oxon (nM)"),
      tibble(time = d$time, value = d$C_oxime_, panel = "Oxime (µM)"),
      tibble(time = d$time, value = d$Ce_atr,   panel = "Atropine effect site (nM)"),
      tibble(time = d$time, value = d$kinh_,    panel = "Inhibition rate k_i·C_oxon (1/h)"),
      tibble(time = d$time, value = d$kr_,      panel = "Reactivation rate k_r[X] (1/h)"),
      tibble(time = d$time, value = d$OMEGA,    panel = "Ω = k_r_max/(k_i·C_oxon)")
    ) %>%
      ggplot(aes(time, value)) +
      geom_line(linewidth = 0.9, colour = "#12507b") +
      facet_wrap(~panel, scales = "free_y", ncol = 2) +
      scale_x_continuous(breaks = c(0, 24, 72, 168, 336)) +
      labs(x = "Time (h)", y = NULL,
           title = "Toxicokinetics and the two rate constants") +
      THEME
  })

  ## ----------------------------------------------------------- ③ esterase
  output$plot_esterase <- renderPlot({
    d <- sim()
    p1 <- d %>% select(time, `Free E` = E_rbc, `Phosphorylated EP` = EP_rbc,
                       `Aged EP-aged` = EA_rbc) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, 100 * value, fill = name)) +
      geom_area(alpha = 0.85) +
      scale_fill_manual(values = c("Free E" = "#2e9e5b", "Phosphorylated EP" = "#f0a202",
                                   "Aged EP-aged" = "#c0392b")) +
      scale_x_continuous(breaks = c(0, 24, 72, 168, 336)) +
      labs(x = "Time (h)", y = "Fraction of erythrocyte AChE (%)", fill = NULL,
           title = "The three states of erythrocyte AChE — the red area cannot be undone") + THEME
    p2 <- d %>% select(time, `RBC AChE` = E_rbc, `NMJ AChE` = E_nmj,
                       `Brain AChE` = E_cns, `Plasma BChE` = B_free) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, 100 * value, colour = name)) +
      geom_line(linewidth = 1.0) +
      geom_hline(yintercept = 30, linetype = 2, colour = "#888") +
      annotate("text", x = 250, y = 33, label = "Symptom threshold ~30%",
               size = 3.4, colour = "#666") +
      scale_x_continuous(breaks = c(0, 24, 72, 168, 336)) +
      labs(x = "Time (h)", y = "Relative to baseline (%)", colour = NULL,
           title = "Esterases by compartment — the oxime reaches the periphery and barely reaches the brain") +
      THEME
    p3 <- d %>% select(time, `Ω` = OMEGA, `Ceiling E_CEIL` = E_CEIL,
                       `Quasi-steady state E_QSS` = E_QSS, `Ageing ratchet φ` = PHI) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.9, colour = "#8e44ad") +
      facet_wrap(~name, scales = "free_y", nrow = 1) +
      scale_x_continuous(breaks = c(0, 72, 168, 336)) +
      labs(x = "Time (h)", y = NULL, title = "The two dimensionless numbers that govern the model") + THEME
    gridExtra_arrange(p1, p2, p3)
  })

  ## -------------------------------------------------------------- ④ signs
  output$plot_signs <- renderPlot({
    d <- sim()
    bind_rows(
      tibble(time = d$time, value = d$ACh_m,   panel = "Muscarinic synaptic ACh (× normal)"),
      tibble(time = d$time, value = d$MUSX,    panel = "Muscarinic signal excess (0-1)"),
      tibble(time = d$time, value = d$NICX,    panel = "Nicotinic signal excess (0-1)"),
      tibble(time = d$time, value = d$CNSX,    panel = "Central muscarinic excess (0-1)"),
      tibble(time = d$time, value = d$SEC,     panel = "Airway secretions (mL)"),
      tibble(time = d$time, value = d$HR,      panel = "Heart rate (bpm)"),
      tibble(time = d$time, value = d$MAP,     panel = "Mean arterial pressure (mmHg)"),
      tibble(time = d$time, value = d$ATRCUM,  panel = "Cumulative atropine (mg)"),
      tibble(time = d$time, value = d$Rm_down, panel = "Receptor downregulation (tolerance)"),
      tibble(time = d$time, value = d$POP,     panel = "POP severity score (0-11)")
    ) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.9, colour = "#b7791f") +
      facet_wrap(~panel, scales = "free_y", ncol = 2) +
      scale_x_continuous(breaks = c(0, 24, 72, 168, 336)) +
      labs(x = "Time (h)", y = NULL,
           title = "The cholinergic syndrome and closed-loop atropine") + THEME
  })

  ## --------------------------------------------------------------- ⑤ resp
  output$plot_resp <- renderPlot({
    d <- sim()
    bind_rows(
      tibble(time = d$time, value = d$RESPD,     panel = "Central respiratory drive RESPD"),
      tibble(time = d$time, value = d$MSTR,      panel = "Respiratory muscle strength MSTR"),
      tibble(time = d$time, value = d$VCAP,      panel = "Spontaneous ventilatory capacity VCAP"),
      tibble(time = d$time, value = d$VENT_ON,   panel = "Mechanical ventilation on (0/1)"),
      tibble(time = d$time, value = d$PACO2,     panel = "PaCO₂ (mmHg)"),
      tibble(time = d$time, value = d$PAO2_,     panel = "PaO₂ (mmHg)"),
      tibble(time = d$time, value = d$LUNG,      panel = "Aspiration lung injury (0-1)"),
      tibble(time = d$time, value = d$MORTALITY, panel = "Cumulative mortality (%)")
    ) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.9, colour = "#a03e6f") +
      facet_wrap(~panel, scales = "free_y", ncol = 2) +
      scale_x_continuous(breaks = c(0, 24, 72, 168, 336)) +
      labs(x = "Time (h)", y = NULL,
           title = "Respiratory failure — the final common path of this poisoning") + THEME
  })

  ## ----------------------------------------------------------- ⑥ scenarios
  scen_runs <- reactive({
    base <- list(compound = input$compound, volume_ml = input$volume,
                 pon1 = input$pon1, wt = input$wt, oxime_start = input$oxstart,
                 oxime_dur = input$oxdur, atr_mode = 1, ventilator = TRUE,
                 oxygen = TRUE, icu = TRUE, glyco = FALSE, magnesium = FALSE,
                 scavenger_mg = 0, charcoal_t = -1)
    defs <- list(
      sup       = modifyList(base, list(oxime = NULL)),
      pam_who   = modifyList(base, list(oxime = "pam_who")),
      pam_bolus = modifyList(base, list(oxime = "pam_bolus")),
      obidoxime = modifyList(base, list(oxime = "obidoxime")),
      novent    = modifyList(base, list(oxime = "pam_who", ventilator = FALSE,
                                        icu = FALSE)),
      slowatr   = modifyList(base, list(oxime = "pam_who", atr_mode = 2)),
      glyco     = modifyList(base, list(oxime = "pam_who", glyco = TRUE)),
      char1     = modifyList(base, list(oxime = "pam_who", charcoal_t = 1)),
      ponqq     = modifyList(base, list(oxime = "pam_who", pon1 = "QQ"))
    )
    sel <- intersect(input$scen, names(defs))
    withProgress(message = "Simulating the scenarios", value = 0, {
      out <- lapply(seq_along(sel), function(i) {
        incProgress(1 / length(sel), detail = sel[i])
        a <- defs[[sel[i]]]
        do.call(run_case, a) %>% mutate(scenario = sel[i])
      })
    })
    bind_rows(out)
  })

  output$plot_scen <- renderPlot({
    d <- scen_runs()
    req(nrow(d) > 0)
    d$y <- d[[input$scen_y]]
    ggplot(d, aes(time, y, colour = scenario)) +
      geom_line(linewidth = 1.0) +
      scale_x_continuous(breaks = c(0, 24, 72, 168, 336)) +
      labs(x = "Time (h)", y = input$scen_y, colour = NULL,
           title = "Scenario comparison") + THEME
  })

  output$tbl_scen <- renderDT({
    d <- scen_runs(); req(nrow(d) > 0)
    d %>% group_by(scenario) %>%
      summarise(`AChE nadir (%)` = round(min(AChE_RBC), 2),
                `AChE 72h (%)`  = round(approx(time, AChE_RBC, 72)$y, 1),
                `Aged 72h (%)`  = round(approx(time, AGED_RBC, 72)$y, 1),
                `Ω 24h`         = round(approx(time, OMEGA, 24)$y, 2),
                `Atropine (mg)` = round(max(ATRCUM)),
                `Ventilation (h)` = round(max(VTIME)),
                `IMS index (%)` = round(max(IMS_IDX)),
                `POP peak`      = round(max(POP), 1),
                `Mortality (%)` = round(max(MORTALITY), 1), .groups = "drop") %>%
      datatable(options = list(dom = "t", pageLength = 12), rownames = FALSE)
  })

  ## --------------------------------------------------------- ⑦ Ω explorer
  omega_tbl <- reactive({
    C <- 10^input$ex_oxon
    rows <- lapply(names(OP_LIBRARY), function(nm) {
      L <- OP_LIBRARY[[nm]]; ks <- log(2) / L$T12SPO
      lapply(c("pam", "obi"), function(ox) {
        k <- if (ox == "obi") L$obi else L$pam
        Xr <- Xreq_fn(k$KRMAX, k$KDOX, L$KI, ks, C, input$ex_target)
        data.frame(
          compound = nm, oxime = ifelse(ox == "obi", "obidoxime", "2-PAM"),
          k_i = L$KI, k_r_max = k$KRMAX,
          Omega = round(omega_fn(k$KRMAX, L$KI, C), 3),
          ceiling_pct = round(100 * ceiling_fn(k$KRMAX, L$KI, C), 1),
          achieved_pct = round(100 * achieved_fn(k$KRMAX, k$KDOX, L$KI, ks, C,
                                                 input$ex_X), 1),
          C_crit_nM = round(Ccrit_fn(k$KRMAX, L$KI, input$ex_target)),
          X_required_uM = ifelse(is.na(Xr), NA, round(Xr)),
          verdict = ifelse(is.na(Xr), "Arithmetically impossible",
                    ifelse(Xr <= input$ex_X, "Reached at the current concentration", "Higher concentration needed")))
      }) %>% bind_rows()
    }) %>% bind_rows()
    rows
  })

  output$plot_omega <- renderPlot({
    Cg <- 10^seq(0, 4, length.out = 240)
    dat <- lapply(names(OP_LIBRARY), function(nm) {
      L <- OP_LIBRARY[[nm]]; ks <- log(2) / L$T12SPO
      bind_rows(
        data.frame(compound = nm, oxime = "2-PAM", C = Cg,
                   ceiling = ceiling_fn(L$pam$KRMAX, L$KI, Cg),
                   achieved = achieved_fn(L$pam$KRMAX, L$pam$KDOX, L$KI, ks,
                                          Cg, input$ex_X)),
        data.frame(compound = nm, oxime = "obidoxime", C = Cg,
                   ceiling = ceiling_fn(L$obi$KRMAX, L$KI, Cg),
                   achieved = achieved_fn(L$obi$KRMAX, L$obi$KDOX, L$KI, ks,
                                          Cg, input$ex_X)))
    }) %>% bind_rows() %>%
      pivot_longer(c(ceiling, achieved), names_to = "which", values_to = "E")
    ggplot(dat, aes(C, 100 * E, colour = oxime, linetype = which)) +
      geom_line(linewidth = 0.95) +
      geom_hline(yintercept = 100 * input$ex_target, linetype = 3,
                 colour = "#c0392b", linewidth = 0.8) +
      geom_vline(xintercept = 10^input$ex_oxon, linetype = 2, colour = "#555") +
      facet_wrap(~compound, nrow = 1) +
      scale_x_log10() +
      scale_linetype_manual(values = c(ceiling = 1, achieved = 2),
                            labels = c(achieved = "Actual, at the current oxime concentration",
                                       ceiling = "Ceiling at infinite dose")) +
      labs(x = "Plasma oxon concentration (nM, log)", y = "Free AChE (%)",
           colour = NULL, linetype = NULL,
           title = "The oxime ceiling — where the target line passes below the solid curve, the oxime is not a drug for that patient") +
      THEME
  })

  output$tbl_omega <- renderDT({
    datatable(omega_tbl(), options = list(dom = "t", pageLength = 10),
              rownames = FALSE)
  })

  ## ------------------------------------------------------------- ⑧ sweep
  sweep_res <- eventReactive(input$go_sweep, {
    vols <- c(2, 5, 10, 20, 35, 50, 100, 200)
    withProgress(message = "Dose-response sweep", value = 0, {
      lapply(vols, function(v) {
        incProgress(1 / length(vols), detail = paste0(v, " mL"))
        a <- run_case(input$compound, NULL, v, input$pon1, input$wt, 1.5, 168,
                      1, TRUE, TRUE, TRUE, FALSE, FALSE, 0, -1)
        b <- run_case(input$compound, "pam_who", v, input$pon1, input$wt, 1.5,
                      168, 1, TRUE, TRUE, TRUE, FALSE, FALSE, 0, -1)
        tibble(volume_ml = v,
               oxon_24h  = approx(b$time, b$C_oxon, 24)$y,
               Omega     = approx(b$time, b$OMEGA, 24)$y,
               ceiling   = 100 * approx(b$time, b$E_CEIL, 24)$y,
               AChE72_supportive = approx(a$time, a$AChE_RBC, 72)$y,
               AChE72_oxime      = approx(b$time, b$AChE_RBC, 72)$y,
               mort_supportive   = max(a$MORTALITY),
               mort_oxime        = max(b$MORTALITY))
      }) %>% bind_rows()
    })
  })

  output$plot_sweep <- renderPlot({
    d <- sweep_res()
    d %>% select(volume_ml, `Supportive care` = AChE72_supportive,
                 `2-PAM` = AChE72_oxime, `Ceiling` = ceiling) %>%
      pivot_longer(-volume_ml) %>%
      ggplot(aes(volume_ml, value, colour = name)) +
      geom_line(linewidth = 1.0) + geom_point(size = 2) +
      geom_hline(yintercept = 30, linetype = 2, colour = "#888") +
      scale_x_log10() +
      labs(x = "Volume ingested (mL, log)", y = "RBC AChE at 72 h (%)", colour = NULL,
           title = "The oxime effect vanishes as the dose rises — not because of the drug but because of the oxon") +
      THEME
  })

  output$tbl_sweep <- renderDT({
    datatable(sweep_res() %>% mutate(across(where(is.numeric), ~round(.x, 2))),
              options = list(dom = "t", pageLength = 10), rownames = FALSE)
  })

  ## ------------------------------------------------------------- ⑨ trial
  trial_res <- eventReactive(input$go_trial, {
    set.seed(20260805)
    n <- input$n_trial
    withProgress(message = "Virtual trial", value = 0, {
      lapply(seq_len(n), function(i) {
        incProgress(1 / n, detail = paste0(i, "/", n))
        ml  <- min(250, max(1, exp(rnorm(1, log(input$trial_med), 0.85))))
        cmpn <- if (runif(1) < input$trial_cpf) "chlorpyrifos" else "dimethoate"
        pon <- if (runif(1) < 0.25) "QQ" else "RR"
        a <- run_case(cmpn, NULL, ml, pon, 60, 1.5, 168, 1, TRUE, TRUE, TRUE,
                      FALSE, FALSE, 0, -1)
        b <- run_case(cmpn, "pam_who", ml, pon, 60, 1.5, 168, 1, TRUE, TRUE,
                      TRUE, FALSE, FALSE, 0, -1)
        tibble(id = i, op = cmpn, ml = ml, pon1 = pon,
               mort_placebo = max(a$MORTALITY), mort_oxime = max(b$MORTALITY),
               vent_placebo = max(a$VTIME),     vent_oxime = max(b$VTIME))
      }) %>% bind_rows()
    })
  })

  output$plot_trial <- renderPlot({
    d <- trial_res()
    d %>% mutate(delta = mort_placebo - mort_oxime) %>%
      ggplot(aes(ml, delta, colour = op)) +
      geom_hline(yintercept = 0, colour = "#999") +
      geom_point(size = 2.4, alpha = 0.85) +
      geom_smooth(se = FALSE, method = "loess", formula = y ~ x,
                  linewidth = 0.8) +
      scale_x_log10() +
      labs(x = "Volume ingested (mL, log)", y = "Reduction in mortality by the oxime (%p)",
           colour = NULL,
           title = "The oxime effect in individual patients — above and below zero mixed inside one trial") +
      THEME
  })

  output$txt_trial <- renderText({
    d <- trial_res()
    strat <- d %>% mutate(s = ifelse(ml <= 15, "<=15 mL", ">15 mL")) %>%
      group_by(op, s) %>%
      summarise(n = n(), pl = mean(mort_placebo), ox = mean(mort_oxime),
                .groups = "drop")
    paste0(
      sprintf("n = %d  |  overall placebo mortality %.1f%%  oxime %.1f%%  RR = %.2f\n",
              nrow(d), mean(d$mort_placebo), mean(d$mort_oxime),
              mean(d$mort_oxime) / mean(d$mort_placebo)),
      sprintf("Hours of mechanical ventilation: placebo %.0f h  oxime %.0f h\n",
              mean(d$vent_placebo), mean(d$vent_oxime)),
      paste(apply(strat, 1, function(r)
        sprintf("  %-14s %-8s n=%-4s placebo %5.1f%%  oxime %5.1f%%",
                r[["op"]], r[["s"]], r[["n"]],
                as.numeric(r[["pl"]]), as.numeric(r[["ox"]]))), collapse = "\n")
    )
  })

  ## --------------------------------------------------------- ⑩ biomarkers
  output$plot_biomarker <- renderPlot({
    d <- sim()
    bind_rows(
      tibble(time = d$time, value = d$AChE_RBC, panel = "Erythrocyte AChE (%)"),
      tibble(time = d$time, value = d$BChE_PL,  panel = "Plasma BChE (%)"),
      tibble(time = d$time, value = d$AGED_RBC, panel = "Irreversibly aged fraction (%)"),
      tibble(time = d$time, value = 100 * d$E_nmj, panel = "NMJ AChE (%)"),
      tibble(time = d$time, value = d$IMS_IDX,  panel = "Intermediate syndrome index (nAChR block %)"),
      tibble(time = d$time, value = 100 * d$N_aged, panel = "Aged NTE (%) — precursor of OPIDN"),
      tibble(time = d$time, value = 100 * d$OPIDN, panel = "OPIDN score (%)"),
      tibble(time = d$time, value = d$POP, panel = "POP severity score")
    ) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.9, colour = "#2e7d5b") +
      facet_wrap(~panel, scales = "free_y", ncol = 2) +
      scale_x_continuous(breaks = c(0, 24, 72, 168, 336)) +
      labs(x = "Time (h)", y = NULL,
           title = "Biomarkers and the delayed syndromes") + THEME
  })

  ## ---------------------------------------------------------------- ⑪ docs
  output$docs <- renderUI({
    if (file.exists("README.md")) includeMarkdown("README.md")
    else tags$p("README.md could not be found.")
  })
}

shinyApp(ui, server)
