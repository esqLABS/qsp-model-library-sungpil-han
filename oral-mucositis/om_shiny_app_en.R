# =============================================================================
#  Oral Mucositis (OM) QSP — Shiny dashboard
#  Oral mucositis quantitative systems pharmacology dashboard
# =============================================================================
#
#  The app is organised around the model's ONE claim: onset belongs to the
#  INSULT term and duration to the REGENERATION term, and the two belong to
#  different drugs.  Tabs 8-11 are not plots of the model, they are the
#  EXPERIMENTS that test that claim — the same experiments that
#  om_analysis.py runs offline.
#
#  Run:
#    install.packages(c("shiny","mrgsolve","dplyr","tidyr","ggplot2",
#                       "DT","tibble"))
#    shiny::runApp("om_shiny_app.R")
#
#  NOTE: no R toolchain was available where this file was written, so it has
#  not been executed.  It mirrors om_python_reference.py, which was.
# =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(DT)

MODEL_FILE <- "om_mrgsolve_model_en.R"
mod <- mread_cache("om", MODEL_FILE)

# -----------------------------------------------------------------------------
# schedule builders (mirrors the block at the bottom of the mrgsolve model)
# -----------------------------------------------------------------------------
FX_DUR <- 10 / 1440          # a 2 Gy fraction is delivered over 10 minutes

blank_cov <- function(times) {
  tibble(ID = 1, time = times, evid = 0, cmt = 0, amt = 0, rate = 0)
}

rt_records <- function(nfx = 35, dpf = 2, per_week = 5, t0 = 0,
                       gap = NULL) {
  day <- t0; given <- 0; starts <- c()
  while (given < nfx && (day - t0) < 250) {
    in_gap <- !is.null(gap) && (day - t0) >= gap[1] && (day - t0) < gap[2]
    if ((round(day - t0) %% 7) < per_week && !in_gap) {
      starts <- c(starts, day); given <- given + 1
    }
    day <- day + 1
  }
  bind_rows(lapply(starts, function(s)
    tibble(ID = 1, time = c(s, s + FX_DUR), evid = 0, cmt = 0, amt = 0,
           rate = 0, RRATE = c(dpf / FX_DUR, 0), DPF = dpf)))
}

ev_inf <- function(times, amt, cmt, dur_d) {
  as_tibble(as.data.frame(
    ev(ID = 1, time = times, amt = amt, cmt = cmt, rate = amt / dur_d)))
}

build_schedule <- function(inp) {
  recs <- list()

  if (inp$regimen == "HDM") {
    amt <- inp$mel_dose * inp$bsa * 1000
    recs[[length(recs) + 1]] <- ev_inf(0, amt, "A_mel_c", 0.5 / 24)

  } else if (inp$regimen == "TBI") {
    dur <- FX_DUR
    st <- c(0, 0.35, 1, 1.35, 2, 2.35, 3, 3.35)
    recs[[length(recs) + 1]] <- bind_rows(lapply(st, function(s)
      tibble(ID = 1, time = c(s, s + dur), evid = 0, cmt = 0, amt = 0,
             rate = 0, RRATE = c(1.5 / dur, 0), DPF = 1.5)))
    amt <- inp$cy_equiv * inp$bsa * 1000
    recs[[length(recs) + 1]] <- ev_inf(c(4, 6), amt / 2, "A_mel_c", 0.5 / 24)

  } else if (inp$regimen == "chemoRT") {
    gap <- if (inp$gap_len > 0) c(inp$gap_start, inp$gap_start + inp$gap_len)
           else NULL
    recs[[length(recs) + 1]] <- rt_records(inp$nfx, inp$dpf, inp$fx_week,
                                           gap = gap)
    if (inp$cis_dose > 0) {
      amt <- inp$cis_dose * inp$bsa * 1000
      recs[[length(recs) + 1]] <- ev_inf(c(0, 21, 42), amt, "A_cis_c", 2 / 24)
    }

  } else if (inp$regimen == "5FU_bolus") {
    amt <- inp$fu_dose * inp$bsa * 1000
    recs[[length(recs) + 1]] <- ev_inf(0:4, amt, "A_5fu_c", 5 / 1440)

  } else if (inp$regimen == "5FU_CI") {
    amt <- inp$fu_total * inp$bsa * 1000
    recs[[length(recs) + 1]] <- ev_inf(0, amt, "A_5fu_c", inp$fu_days)
  }

  # ---- interventions --------------------------------------------------------
  if (inp$cryo_h > 0) {
    recs[[length(recs) + 1]] <- tibble(ID = 1, time = c(0, inp$cryo_h / 24),
                                       evid = 0, cmt = 0, amt = 0, rate = 0,
                                       CRYO = c(1, 0))
  }
  if (inp$palifermin) {
    amt <- 60 * inp$wt * 1000
    d0 <- inp$pal_gap
    days <- c(-2, -1, 0) - d0 + inp$pal_anchor
    days2 <- inp$pal_post + c(0, 1, 2)
    recs[[length(recs) + 1]] <-
      as_tibble(as.data.frame(ev(ID = 1, time = c(days, days2), amt = amt,
                                 cmt = "A_pal_c")))
  }
  if (inp$pbm) {
    d <- seq(0, inp$pbm_days, by = 1)
    recs[[length(recs) + 1]] <- bind_rows(lapply(d, function(x)
      tibble(ID = 1, time = c(x, x + 5 / 1440), evid = 0, cmt = 0, amt = 0,
             rate = 0, PBMR = c(inp$pbm_fluence / (5 / 1440), 0))))
  }
  if (inp$benzydamine) {
    recs[[length(recs) + 1]] <- tibble(ID = 1, time = c(0, inp$t_end),
                                       evid = 0, cmt = 0, amt = 0, rate = 0,
                                       BZDL = c(1, 0))
  }
  if (inp$glutamine) {
    recs[[length(recs) + 1]] <- tibble(ID = 1, time = c(0, inp$t_end),
                                       evid = 0, cmt = 0, amt = 0, rate = 0,
                                       GLNL = c(1, 0))
  }

  d <- bind_rows(recs)
  for (nm in c("RRATE", "CRYO", "PBMR", "BZDL", "GLNL")) {
    if (!nm %in% names(d)) d[[nm]] <- 0
    d[[nm]][is.na(d[[nm]])] <- 0
  }
  if (!"DPF" %in% names(d)) d$DPF <- 2
  d$DPF[is.na(d$DPF)] <- 2
  for (nm in c("evid", "cmt", "amt", "rate")) {
    if (!nm %in% names(d)) d[[nm]] <- 0
    d[[nm]][is.na(d[[nm]])] <- 0
  }
  # LOCF the covariates so they persist between records
  d <- d %>% arrange(time)
  d
}

run_sim <- function(inp, pset = list(), sched = NULL, t_end = NULL) {
  s <- if (is.null(sched)) build_schedule(inp) else sched
  te <- if (is.null(t_end)) inp$t_end else t_end
  m <- mod
  pset$sens <- (pset$sens %||% 1) * inp$sens
  if (length(pset)) m <- param(m, pset)
  m %>% data_set(s) %>%
    mrgsim(end = te, delta = 0.05, recsort = 3) %>% as_tibble()
}

`%||%` <- function(a, b) if (is.null(a)) b else a

summarise_run <- function(out, sev = 3) {
  dt <- 0.05
  sm <- out$WHO >= sev
  tibble(
    onset_sev = if (any(sm)) min(out$time[sm]) else NA_real_,
    dur_sev   = sum(sm) * dt,
    peak_area = max(out$ULC),
    peak_who  = max(out$WHO),
    area_auc  = sum(out$ULC) * dt,
    pain_auc  = sum(out$VAS) * dt,
    opioid_pk = max(out$MEDMG),
    anc_nadir = min(out$ANC))
}

# =============================================================================
# UI
# =============================================================================
ui <- fluidPage(
  titlePanel(HTML(paste0(
    "<b>Oral Mucositis QSP Dashboard</b><br>",
    "<span style='font-size:14px;color:#555'>",
    "One asymmetry: every cytotoxic agent strikes the <b>proliferative ",
    "compartments</b> of the epithelial belt and leaves the <b>differentiated ",
    "barrier</b> untouched → <b>onset is set by the insult term</b>, <b>duration by the regeneration term</b>",
    "</span>"))),
  tags$style(HTML("
    .well { background:#fafbfc; }
    .nav-tabs > li > a { padding:6px 10px; font-size:13px; }
    .clockbox { border-left:4px solid #b9770e; padding:6px 12px;
                background:#fef9e7; margin-bottom:8px; font-size:13px; }
  ")),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("① Regimen"),
      selectInput("regimen", NULL,
                  c("High-dose melphalan (HDM, autologous HSCT)" = "HDM",
                    "TBI + VP16/Cy conditioning" = "TBI",
                    "Head and neck chemoradiotherapy (70 Gy + cisplatin)" = "chemoRT",
                    "5-FU bolus d1-5" = "5FU_bolus",
                    "5-FU 96 h continuous infusion" = "5FU_CI"),
                  selected = "HDM"),
      conditionalPanel("input.regimen == 'HDM'",
        sliderInput("mel_dose", "Melphalan (mg/m²)", 100, 200, 200, 10)),
      conditionalPanel("input.regimen == 'TBI'",
        sliderInput("cy_equiv", "Alkylator equivalent dose (mg/m² melphalan equivalent)",
                    20, 300, 120, 10)),
      conditionalPanel("input.regimen == 'chemoRT'",
        sliderInput("nfx", "Number of fractions", 20, 70, 35, 1),
        sliderInput("dpf", "Dose per fraction (Gy)", 1.0, 2.5, 2.0, 0.1),
        sliderInput("fx_week", "Fractions per week", 5, 10, 5, 1),
        sliderInput("cis_dose", "Cisplatin (mg/m², 0 = RT alone)",
                    0, 100, 100, 10),
        sliderInput("gap_start", "Treatment gap start day", 0, 45, 28, 1),
        sliderInput("gap_len", "Treatment gap length (days)", 0, 21, 0, 1)),
      conditionalPanel("input.regimen == '5FU_bolus'",
        sliderInput("fu_dose", "5-FU bolus (mg/m²/day)", 200, 600, 425, 25)),
      conditionalPanel("input.regimen == '5FU_CI'",
        sliderInput("fu_total", "5-FU total dose (mg/m²)", 1000, 6000, 4000, 250),
        sliderInput("fu_days", "Infusion duration (days)", 1, 7, 4, 1)),

      hr(),
      h4("② Insult-arm intervention"),
      div(class = "clockbox",
          "Clock 1 — insult: moves onset time and incidence"),
      sliderInput("cryo_h", "Oral cryotherapy (hours)", 0, 8, 0, 0.5),

      hr(),
      h4("③ Regeneration-arm intervention"),
      div(class = "clockbox",
          "Clock 2 — regeneration: moves duration"),
      checkboxInput("palifermin", "Palifermin 60 µg/kg/day × 3+3", FALSE),
      conditionalPanel("input.palifermin",
        sliderInput("pal_gap",
                    "Interval between the last pre-conditioning dose and cytotoxic therapy (days)",
                    -1, 3, 1, 0.5),
        sliderInput("pal_anchor", "Time of cytotoxic dosing (days)", 0, 8, 0, 1),
        sliderInput("pal_post", "Post-transplant re-dosing start day", 1, 14, 3, 1)),
      checkboxInput("pbm", "Photobiomodulation (daily)", FALSE),
      conditionalPanel("input.pbm",
        sliderInput("pbm_fluence", "Fluence (J/cm²)", 1, 12, 6, 1),
        sliderInput("pbm_days", "Treatment duration (days)", 5, 40, 20, 1)),
      checkboxInput("glutamine", "Oral glutamine", FALSE),

      hr(),
      h4("④ Anti-inflammatory · host"),
      checkboxInput("benzydamine", "Benzydamine mouthwash", FALSE),
      sliderInput("sens", "Individual susceptibility multiplier (sens)", 0.5, 2.0, 1.0, 0.05),
      sliderInput("bsa", "Body surface area (m²)", 1.4, 2.2, 1.8, 0.05),
      sliderInput("wt", "Body weight (kg)", 45, 110, 75, 1),
      sliderInput("t_end", "Simulation duration (days)", 20, 130, 45, 5),
      helpText(HTML("For the provenance and calibration history of every parameter see ",
                    "<code>calib.log</code> and ",
                    "<code>om_reference_output.txt</code>."))
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        type = "tabs",

        tabPanel("1. Patient · regimen",
          h4("Currently configured dosing / irradiation schedule"),
          p(HTML("Radiation fractions are delivered over a <b>10-minute window</b>. If the ",
                 "integrator steps over that window, the dose actually delivered is governed ",
                 "not by the prescription but by where the integrator placed its steps, so ",
                 "every fraction boundary is stated explicitly as a data set record.")),
          DTOutput("tbl_sched"),
          h4("Summary metrics"), DTOutput("tbl_summary")),

        tabPanel("2. Pharmacokinetics (PK)",
          h4("Plasma concentrations"), plotOutput("p_pk", height = 300),
          h4("Mucosal compartment concentrations — the only point where cryotherapy acts"),
          p(HTML("The mucosa is a perfusion-limited compartment with an equilibration ",
                 "time constant of about 2 minutes. So <b>reducing blood flow alone</b> ",
                 "cannot lower exposure by much — tab 9 breaks this down.")),
          plotOutput("p_muc", height = 300)),

        tabPanel("3. Epithelial belt",
          h4("S → P1 → P2 → P3 → D — the spine of the model"),
          p(HTML("Cytotoxicity strikes <b>only the proliferative compartments (S, P)</b>. ",
                 "The differentiated barrier D carries no cytotoxic term. ",
                 "That is why a latent period exists, and its length is not a ",
                 "signalling mystery but <b>1/k_shed</b>.")),
          plotOutput("p_belt", height = 380),
          plotOutput("p_freg", height = 240)),

        tabPanel("4. Ulceration · severity",
          plotOutput("p_ulcer", height = 340),
          plotOutput("p_scales", height = 300),
          p(HTML("<b>The WHO grade is ordinal and saturates at grade 4</b>. ",
                 "OMAS and area do not saturate — see tab 12."))),

        tabPanel("5. Inflammatory cascade",
          h4("Sonis phases 2·3 — a single fast signalling loop"),
          p(HTML("The five 'stages' are not sequential events but ",
                 "<b>a single fast signalling loop</b> laid on top of a slow tissue variable (the belt). ",
                 "The loop gain of the TNF → NF-κB positive feedback must be less ",
                 "than 1, otherwise the model diverges.")),
          plotOutput("p_infl", height = 340),
          h4("Microbial colonisation and ceramide at the ulcer surface"),
          plotOutput("p_micro", height = 260)),

        tabPanel("6. Pain · analgesia",
          plotOutput("p_pain", height = 320),
          h4("Opioid analgesic requirement (titrated to a VAS ≤ 4 target)"),
          plotOutput("p_opioid", height = 260)),

        tabPanel("7. Myelosuppression · infection",
          plotOutput("p_anc", height = 300),
          h4("The ulcer is a portal of entry — risk ∝ area × bacterial load / ANC"),
          p(HTML("The mucosal curve and the ANC curve are driven by the <b>same</b> ",
                 "melphalan exposure but have completely different time constants. So the ",
                 "risk is the product of two curves that peak at different times.")),
          plotOutput("p_inf", height = 260)),

        tabPanel("8. Two clocks",
          h4("Testing the model's central claim as a 2×2"),
          p(HTML("Perturb the insult term (<code>sens</code>) and the regeneration term (<code>lamS</code>) ",
                 "by ±40% each and watch how <b>onset</b> and <b>duration</b> ",
                 "move. If the claim is right, each term should move only its own ",
                 "endpoint appreciably.")),
          actionButton("run_clocks", "Run", class = "btn-primary"),
          br(), br(),
          DTOutput("tbl_clocks"),
          plotOutput("p_clocks", height = 320)),

        tabPanel("9. Cryotherapy criterion",
          h4("Ice works only while it is in the mouth"),
          p(HTML("So the benefit = (1 − f_cryo) × (the fraction of the insult ",
                 "AUC that falls inside the ice window), and for a drug of half-life t½ cooled for T ",
                 "immediately after dosing that fraction is exactly <b>1 − 2^(−T/t½)</b>. ",
                 "The pattern of 'which regimens it works in' that the guidelines report is ",
                 "<b>derived</b> here.")),
          actionButton("run_cryo", "Run", class = "btn-primary"),
          br(), br(),
          plotOutput("p_cryo", height = 330),
          DTOutput("tbl_cryo")),

        tabPanel("10. Palifermin scheduling",
          h4("Turning regeneration on also enlarges the target"),
          p(HTML("KGF raises the proliferation rate — that is the point. ",
                 "But <b>the enlarged proliferative pool is a bigger target for ",
                 "cell-cycle-specific cytotoxicity</b>. Given concurrently the benefit shrinks and ",
                 "can reverse, and this is what gives rise to the label's 24-hour separation ",
                 "requirement. The sweep below <b>reproduces</b> that requirement.")),
          actionButton("run_pal", "Run", class = "btn-primary"),
          br(), br(),
          plotOutput("p_pal", height = 330),
          DTOutput("tbl_pal")),

        tabPanel("11. Fractionation · treatment gap",
          h4("Mucosal recovery and tumour repopulation compete over the same time"),
          actionButton("run_frac", "Run", class = "btn-primary"),
          br(), br(),
          DTOutput("tbl_frac"),
          plotOutput("p_frac", height = 300)),

        tabPanel("12. Scale saturation · test sensitivity",
          h4("The WHO scale loses power exactly where the disease is worst"),
          p(HTML("Measure the same intervention with the WHO grade and with ulcer area, ",
                 "changing baseline severity through <code>sens</code> alone. ",
                 "If the WHO column collapses to 0 while the area column does not, ",
                 "what changed is not the drug but the <b>instrument</b>.")),
          actionButton("run_sat", "Run", class = "btn-primary"),
          br(), br(),
          DTOutput("tbl_sat"),
          plotOutput("p_sat", height = 300)),

        tabPanel("13. Scenario comparison",
          h4("Run all 17 scenarios at once"),
          actionButton("run_scen", "Run", class = "btn-primary"),
          br(), br(),
          DTOutput("tbl_scen"),
          plotOutput("p_scen", height = 340)),

        tabPanel("14. Virtual population",
          h4("Between-subject variability sits only in the three places the structure allows"),
          p(HTML("Cytotoxic susceptibility · regenerative capacity · barrier threshold. ",
                 "Everything else is shared.")),
          sliderInput("npop", "Population size", 20, 400, 100, 20),
          actionButton("run_pop", "Run", class = "btn-primary"),
          br(), br(),
          DTOutput("tbl_pop"),
          plotOutput("p_pop", height = 320)),

        tabPanel("15. Documentation",
          h4("The numbers this model used"),
          verbatimTextOutput("txt_prov"))
      )
    )
  )
)

# =============================================================================
# server
# =============================================================================
server <- function(input, output, session) {

  sim <- reactive({
    req(input$regimen)
    run_sim(input)
  })

  output$tbl_sched <- renderDT({
    datatable(build_schedule(input), options = list(pageLength = 12,
                                                    scrollX = TRUE))
  })

  output$tbl_summary <- renderDT({
    datatable(summarise_run(sim()) %>%
                mutate(across(everything(), ~ round(.x, 3))),
              options = list(dom = "t"))
  })

  # ---- 2. PK ---------------------------------------------------------------
  output$p_pk <- renderPlot({
    d <- sim() %>% select(time, CMEL, C5FU, CCIS, CPAL, CMOR) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value)) + geom_line(linewidth = 0.7) +
      facet_wrap(~ name, scales = "free_y", nrow = 1) +
      labs(x = "Day", y = "Concentration") + theme_bw()
  })

  output$p_muc <- renderPlot({
    d <- sim() %>% select(time, Cm_mel, Cm_5fu, Cm_cis, Cm_mtx) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value)) + geom_line(colour = "#148f77",
                                            linewidth = 0.8) +
      facet_wrap(~ name, scales = "free_y", nrow = 1) +
      labs(x = "Day", y = "Mucosal concentration") + theme_bw()
  })

  # ---- 3. belt -------------------------------------------------------------
  output$p_belt <- renderPlot({
    d <- sim() %>%
      transmute(time, `S clonogenic basal cells` = S, `Sd lethally damaged cells` = Sd,
                `P1+P2+P3 amplifying compartments` = P1 + P2 + P3,
                `D differentiated barrier` = BARRIER) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) +
      geom_line(linewidth = 0.9) +
      geom_hline(yintercept = 0.34, linetype = "dashed",
                 colour = "#c0392b") +
      annotate("text", x = 0, y = 0.36, hjust = 0, size = 3,
               colour = "#c0392b", label = "D_crit — ulceration threshold") +
      labs(x = "Day", y = "Normalised value", colour = NULL) +
      theme_bw() + theme(legend.position = "top")
  })

  output$p_freg <- renderPlot({
    d <- sim() %>% transmute(time, KGFe, PBMe, GLNe) %>% pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 0.8) +
      labs(x = "Day", y = "Regeneration modulator pools", colour = NULL) + theme_bw()
  })

  # ---- 4. ulcer ------------------------------------------------------------
  output$p_ulcer <- renderPlot({
    ggplot(sim(), aes(time, ULC)) +
      geom_area(fill = "#e59866", alpha = 0.6) +
      geom_line(linewidth = 0.8) +
      labs(x = "Day", y = "Ulcerated area fraction A_ulc") + theme_bw()
  })

  output$p_scales <- renderPlot({
    d <- sim() %>% transmute(time, `WHO grade (0-4)` = WHO,
                             `OMAS (0-5)` = OMAS,
                             `Ulcer area ×5` = 5 * ULC) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) +
      geom_step(data = ~ filter(.x, name == "WHO grade (0-4)"),
                linewidth = 0.9) +
      geom_line(data = ~ filter(.x, name != "WHO grade (0-4)"),
                linewidth = 0.8) +
      labs(x = "Day", y = NULL, colour = NULL) +
      theme_bw() + theme(legend.position = "top")
  })

  # ---- 5. inflammation -----------------------------------------------------
  output$p_infl <- renderPlot({
    d <- sim() %>% select(time, ROS, NFkB, TNF, IL1b, IL6) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 0.8) +
      facet_wrap(~ name, scales = "free_y") +
      labs(x = "Day", y = NULL) + theme_bw() +
      theme(legend.position = "none")
  })

  output$p_micro <- renderPlot({
    d <- sim() %>% select(time, CER, MB) %>% pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 0.8) +
      labs(x = "Day", y = NULL, colour = NULL) + theme_bw()
  })

  # ---- 6. pain -------------------------------------------------------------
  output$p_pain <- renderPlot({
    d <- sim() %>% transmute(time, `VAS before analgesia` = VAS_RAW,
                             `VAS after analgesia` = VAS) %>% pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      geom_hline(yintercept = 4, linetype = "dashed") +
      labs(x = "Day", y = "VAS (0-10)", colour = NULL) +
      theme_bw() + theme(legend.position = "top")
  })

  output$p_opioid <- renderPlot({
    ggplot(sim(), aes(time, MEDMG)) + geom_area(fill = "#f0b27a",
                                                alpha = 0.7) +
      labs(x = "Day", y = "Intravenous morphine equivalent (mg/day)") + theme_bw()
  })

  # ---- 7. marrow -----------------------------------------------------------
  output$p_anc <- renderPlot({
    ggplot(sim(), aes(time, ANC)) + geom_line(linewidth = 0.9) +
      geom_hline(yintercept = 0.5, linetype = "dashed",
                 colour = "#c0392b") +
      labs(x = "Day", y = "ANC (×10⁹/L)") + theme_bw()
  })

  output$p_inf <- renderPlot({
    d <- sim() %>% transmute(time, `Ulcer area` = ULC,
                             `Bacterial load` = MB,
                             `Cumulative bacteraemia risk` = cInf) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 0.8) +
      labs(x = "Day", y = NULL, colour = NULL) +
      theme_bw() + theme(legend.position = "top")
  })

  # ---- 8. two clocks -------------------------------------------------------
  clocks <- eventReactive(input$run_clocks, {
    base <- summarise_run(run_sim(input))
    rows <- list(mutate(base, perturbation = "Baseline"))
    for (r in list(list("Insult ×0.6", "sens", 0.6),
                   list("Insult ×1.4", "sens", 1.4),
                   list("Regeneration ×0.6", "lamS", 0.6),
                   list("Regeneration ×1.4", "lamS", 1.4))) {
      pv <- as.numeric(param(mod)[[r[[2]]]]) * r[[3]]
      p <- setNames(list(pv), r[[2]])
      rows[[length(rows) + 1]] <- summarise_run(run_sim(input, pset = p)) %>%
        mutate(perturbation = r[[1]])
    }
    bind_rows(rows) %>% relocate(perturbation)
  })
  output$tbl_clocks <- renderDT(datatable(
    clocks() %>% mutate(across(where(is.numeric), ~ round(.x, 3))),
    options = list(dom = "t")))
  output$p_clocks <- renderPlot({
    d <- clocks() %>% select(perturbation, onset_sev, dur_sev) %>%
      pivot_longer(-perturbation)
    ggplot(d, aes(perturbation, value, fill = name)) +
      geom_col(position = "dodge") + coord_flip() +
      labs(x = NULL, y = "Days", fill = NULL) + theme_bw()
  })

  # ---- 9. cryotherapy ------------------------------------------------------
  cryo <- eventReactive(input$run_cryo, {
    th <- 1.2 / 24
    hrs <- c(0, 0.25, 0.5, 1, 2, 4, 6, 8)
    bind_rows(lapply(hrs, function(h) {
      i2 <- modifyList(reactiveValuesToList(input), list(cryo_h = h))
      o <- run_sim(i2)
      s <- summarise_run(o)
      tibble(ice_h = h, dur_sev = s$dur_sev, peak_area = s$peak_area,
             auc_muc = max(o$cAUCmuc),
             analytic_cover = 1 - 2 ^ (-(h / 24) / th))
    })) %>% mutate(auc_ratio = auc_muc / auc_muc[1])
  })
  output$p_cryo <- renderPlot({
    d <- cryo()
    ggplot(d, aes(ice_h)) +
      geom_line(aes(y = 1 - auc_ratio, colour = "Model (AUC reduction)"),
                linewidth = 1) +
      geom_line(aes(y = analytic_cover * (1 - 0.22),
                    colour = "Analytic (1−2^(−T/t½))×(1−f)"),
                linetype = "dashed", linewidth = 1) +
      labs(x = "Ice retention time (h)", y = "Fraction of insult averted",
           colour = NULL) + theme_bw() + theme(legend.position = "top")
  })
  output$tbl_cryo <- renderDT(datatable(
    cryo() %>% mutate(across(where(is.numeric), ~ round(.x, 4))),
    options = list(dom = "t")))

  # ---- 10. palifermin ------------------------------------------------------
  pal <- eventReactive(input$run_pal, {
    gaps <- c(-1, -0.5, 0, 0.5, 1, 2, 3)
    ref <- summarise_run(run_sim(modifyList(reactiveValuesToList(input),
                                            list(palifermin = FALSE))))
    bind_rows(lapply(gaps, function(g) {
      i2 <- modifyList(reactiveValuesToList(input),
                       list(palifermin = TRUE, pal_gap = g))
      s <- summarise_run(run_sim(i2))
      tibble(gap_d = g, dur_sev = s$dur_sev, peak_area = s$peak_area,
             delta_pct = 100 * (s$dur_sev - ref$dur_sev) / ref$dur_sev)
    })) %>% mutate(no_pal_dur = ref$dur_sev)
  })
  output$p_pal <- renderPlot({
    ggplot(pal(), aes(gap_d, dur_sev)) +
      geom_line(linewidth = 1, colour = "#117864") +
      geom_point(size = 2) +
      geom_hline(aes(yintercept = no_pal_dur), linetype = "dashed",
                 colour = "#c0392b") +
      annotate("text", x = -1, y = Inf, vjust = 1.4, hjust = 0, size = 3.4,
               colour = "#c0392b", label = "No palifermin") +
      labs(x = "Last pre-dose ↔ cytotoxic therapy interval (days)",
           y = "Duration of severe mucositis (days)") + theme_bw()
  })
  output$tbl_pal <- renderDT(datatable(
    pal() %>% mutate(across(where(is.numeric), ~ round(.x, 3))),
    options = list(dom = "t")))

  # ---- 11. fractionation ---------------------------------------------------
  frac <- eventReactive(input$run_frac, {
    arms <- list(
      list("70 Gy/35 fx/7 weeks (standard)", 35, 2.0, 5, 0),
      list("81.6 Gy/68 fx b.i.d.", 68, 1.2, 10, 0),
      list("70 Gy/35 fx/6 weeks (accelerated)", 35, 2.0, 6, 0),
      list("70 Gy + 7-day gap", 35, 2.0, 5, 7))
    bind_rows(lapply(arms, function(a) {
      i2 <- modifyList(reactiveValuesToList(input),
                       list(regimen = "chemoRT", nfx = a[[2]], dpf = a[[3]],
                            fx_week = a[[4]], gap_len = a[[5]],
                            t_end = 120))
      o <- run_sim(i2); s <- summarise_run(o)
      tibble(arm = a[[1]], onset = s$onset_sev, dur = s$dur_sev,
             peak_area = s$peak_area, BED_tumour = max(o$cBEDt),
             BED_mucosa = max(o$cBEDm))
    }))
  })
  output$tbl_frac <- renderDT(datatable(
    frac() %>% mutate(across(where(is.numeric), ~ round(.x, 2))),
    options = list(dom = "t")))
  output$p_frac <- renderPlot({
    ggplot(frac(), aes(BED_tumour, dur, label = arm)) +
      geom_point(size = 3, colour = "#1f618d") +
      geom_text(vjust = -0.9, size = 3.4) +
      labs(x = "Tumour BED (Gy10, repopulation-corrected)",
           y = "Duration of severe mucositis (days)") + theme_bw()
  })

  # ---- 12. WHO saturation --------------------------------------------------
  sat <- eventReactive(input$run_sat, {
    bind_rows(lapply(c(0.7, 1.0, 1.3, 1.7), function(sf) {
      i0 <- modifyList(reactiveValuesToList(input),
                       list(sens = sf, cryo_h = 0))
      i1 <- modifyList(reactiveValuesToList(input),
                       list(sens = sf, cryo_h = 6))
      a <- summarise_run(run_sim(i0)); b <- summarise_run(run_sim(i1))
      tibble(sens = sf,
             area_ctrl = a$area_auc, area_trt = b$area_auc,
             d_area_pct = 100 * (b$area_auc - a$area_auc) / a$area_auc,
             who_ctrl = a$dur_sev, who_trt = b$dur_sev,
             d_who_pct = 100 * (b$dur_sev - a$dur_sev) /
               pmax(a$dur_sev, 1e-9))
    }))
  })
  output$tbl_sat <- renderDT(datatable(
    sat() %>% mutate(across(where(is.numeric), ~ round(.x, 2))),
    options = list(dom = "t")))
  output$p_sat <- renderPlot({
    d <- sat() %>% select(sens, `Area-based effect %` = d_area_pct,
                          `WHO-based effect %` = d_who_pct) %>%
      pivot_longer(-sens)
    ggplot(d, aes(sens, value, colour = name)) +
      geom_line(linewidth = 1) + geom_point(size = 2) +
      labs(x = "Baseline severity (sens)", y = "Measured effect size (%)",
           colour = NULL) + theme_bw() + theme(legend.position = "top")
  })

  # ---- 13. scenarios -------------------------------------------------------
  scen <- eventReactive(input$run_scen, {
    S <- list(
      list("HDM 200, no prophylaxis", list(regimen = "HDM", mel_dose = 200)),
      list("HDM 140 (dose-reduced)", list(regimen = "HDM", mel_dose = 140)),
      list("HDM 200 + cryotherapy 30 min", list(regimen = "HDM", cryo_h = 0.5)),
      list("HDM 200 + cryotherapy 6 h", list(regimen = "HDM", cryo_h = 6)),
      list("HDM 200 + palifermin (separated)",
           list(regimen = "HDM", palifermin = TRUE, pal_gap = 1)),
      list("HDM 200 + palifermin (concurrent)",
           list(regimen = "HDM", palifermin = TRUE, pal_gap = -1)),
      list("HDM 200 + photobiomodulation", list(regimen = "HDM", pbm = TRUE)),
      list("HDM 200 + glutamine", list(regimen = "HDM", glutamine = TRUE)),
      list("HDM 200 + cryotherapy + palifermin",
           list(regimen = "HDM", cryo_h = 6, palifermin = TRUE)),
      list("TBI-VP16-Cy (placebo)", list(regimen = "TBI")),
      list("TBI-VP16-Cy + palifermin",
           list(regimen = "TBI", palifermin = TRUE)),
      list("Head and neck 70 Gy + cisplatin",
           list(regimen = "chemoRT", t_end = 120)),
      list("Head and neck 70 Gy alone",
           list(regimen = "chemoRT", cis_dose = 0, t_end = 120)),
      list("Hyperfractionated 81.6 Gy/68 fx",
           list(regimen = "chemoRT", nfx = 68, dpf = 1.2, fx_week = 10,
                cis_dose = 0, t_end = 120)),
      list("Head and neck chemoradiotherapy + benzydamine",
           list(regimen = "chemoRT", benzydamine = TRUE, t_end = 120)),
      list("5-FU bolus d1-5", list(regimen = "5FU_bolus", t_end = 40)),
      list("5-FU 96 h continuous infusion", list(regimen = "5FU_CI", t_end = 40)))
    bind_rows(lapply(seq_along(S), function(i) {
      i2 <- modifyList(reactiveValuesToList(input), S[[i]][[2]])
      summarise_run(run_sim(i2)) %>% mutate(no = i, scenario = S[[i]][[1]])
    })) %>% relocate(no, scenario)
  })
  output$tbl_scen <- renderDT(datatable(
    scen() %>% mutate(across(where(is.numeric), ~ round(.x, 3))),
    options = list(pageLength = 17, dom = "t")))
  output$p_scen <- renderPlot({
    ggplot(scen(), aes(onset_sev, dur_sev, label = scenario)) +
      geom_point(size = 3, colour = "#b9770e") +
      geom_text(vjust = -0.8, size = 3, check_overlap = TRUE) +
      labs(x = "Onset day of severe mucositis (clock 1: insult)",
           y = "Duration of severe mucositis (clock 2: regeneration)") + theme_bw()
  })

  # ---- 14. virtual population ----------------------------------------------
  pop <- eventReactive(input$run_pop, {
    n <- input$npop
    set.seed(20260806)
    idat <- tibble(
      ID    = 1:n,
      sens  = input$sens * exp(rnorm(n, 0, 0.38) - 0.5 * 0.38 ^ 2),
      lamS  = as.numeric(param(mod)[["lamS"]]) *
                exp(rnorm(n, 0, 0.30) - 0.5 * 0.30 ^ 2),
      Dcrit = 0.34 * exp(rnorm(n, 0, 0.07) - 0.5 * 0.07 ^ 2))
    s <- build_schedule(input)
    s <- bind_rows(lapply(1:n, function(i) mutate(s, ID = i)))
    o <- mod %>% idata_set(idat) %>% data_set(s) %>%
      mrgsim(end = input$t_end, delta = 0.05, recsort = 3) %>% as_tibble()
    o %>% group_by(ID) %>%
      summarise(dur_sev = sum(WHO >= 3) * 0.05,
                peak_who = max(WHO), peak_area = max(ULC),
                .groups = "drop")
  })
  output$tbl_pop <- renderDT({
    d <- pop()
    datatable(tibble(
      `Severe (≥3) incidence` = mean(d$dur_sev > 0),
      `Grade 4 incidence`    = mean(d$peak_who >= 4),
      `Median duration (days)` = median(d$dur_sev[d$dur_sev > 0]),
      `Mean peak area`   = mean(d$peak_area)) %>%
        mutate(across(everything(), ~ round(.x, 3))),
      options = list(dom = "t"))
  })
  output$p_pop <- renderPlot({
    ggplot(pop(), aes(dur_sev)) +
      geom_histogram(bins = 25, fill = "#5b5b77", colour = "white") +
      labs(x = "Duration of severe mucositis (days)", y = "Number of patients") + theme_bw()
  })

  # ---- 15. provenance ------------------------------------------------------
  output$txt_prov <- renderText({
    paste(readLines("om_reference_output.txt", warn = FALSE)[1:60],
          collapse = "\n")
  })
}

shinyApp(ui, server)
