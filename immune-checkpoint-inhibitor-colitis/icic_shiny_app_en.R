## ---------------------------------------------------------------------------
## icic_shiny_app.R
## ===========================================================================
## Interactive dashboard for the IMMUNE CHECKPOINT INHIBITOR COLITIS QSP model.
##
## The app is built around the model's central claim, so the controls are
## arranged to let a user break it:
##
##   Phi_col = (a_eff*Teff + a_trm*Trm) / (a_reg*Treg + kappa)
##
##   - move the ANTI-PD-1 dose slider across its whole range and watch the
##     occupancy trace refuse to move (tab 3).  That is the flat dose-toxicity
##     curve, drawn.
##   - move the IPILIMUMAB dose slider and watch the ADCC trace move
##     proportionally while its own occupancy trace stays flat (tab 3).
##   - watch S_eff fall for weeks with a flat stool chart (tab 4).  The dashed
##     line is S* = L/A_max: nothing is reported until the curve crosses it.
##   - choose a rescue agent and read the onset time off tab 5, and the
##     tumour cost off tab 7.  They are not the same ranking.
##
## Requires: shiny, mrgsolve, dplyr, tidyr, ggplot2, DT
## Run:  shiny::runApp("icic_shiny_app.R")
## ---------------------------------------------------------------------------

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

MODEL_FILE <- "icic_mrgsolve_model.R"
BW <- 70

mod <- mread_cache("icic", MODEL_FILE)

## Phenotype presets. The MEDIAN patient does not get colitis at any of these
## doses -- only ~7-12% of ipilimumab-treated patients do -- so a dashboard
## that opened on the median would look broken. It opens on the susceptible
## phenotype and says so.
PHENO <- list(
  "Median patient (does NOT develop colitis — this is correct)" =
    list(par = list(ag_scale = 1.00, kappa = 0.15, phi_FCGR = 1.00),
         ini = list(Dv = 1.00, Bu = 1.00, Trm = 0.02)),
  "Susceptible (FCGR3A V/V, thin regulatory floor, prior antibiotics)" =
    list(par = list(ag_scale = 1.45, kappa = 0.11, phi_FCGR = 1.60),
         ini = list(Dv = 0.70, Bu = 0.70^1.3, Trm = 0.06)),
  "Highly susceptible (steroid-refractory phenotype)" =
    list(par = list(ag_scale = 1.75, kappa = 0.095, phi_FCGR = 1.60),
         ini = list(Dv = 0.60, Bu = 0.60^1.3, Trm = 0.12)),
  "Pre-existing microscopic colitis (latent primed pool)" =
    list(par = list(ag_scale = 1.30, kappa = 0.13, phi_FCGR = 1.00),
         ini = list(Dv = 0.85, Bu = 0.85^1.3, Trm = 0.30, Ent = 0.85)),
  "FCGR3A F/F (low-affinity FcgammaRIIIa — inefficient ADCC)" =
    list(par = list(ag_scale = 1.45, kappa = 0.11, phi_FCGR = 0.55),
         ini = list(Dv = 0.70, Bu = 0.70^1.3, Trm = 0.06))
)

## Prednisolone as a continuous mg/day input: a 3 h half-life drug with a
## ~12 h genomic effect half-life is pharmacodynamically indistinguishable
## from a daily dose, and this keeps the taper a smooth covariate.
pred_rate <- function(tau, top, plateau, taper_wk) {
  ifelse(tau < 0, 0,
         ifelse(tau < plateau, top,
                {
                  wk <- (tau - plateau) / 7
                  ifelse(wk >= taper_wk, 0, top * (1 - wk / taper_wk))
                }))
}

build_data <- function(inp) {
  tt <- seq(0, inp$tend, by = 1)
  cov <- data.frame(ID = 1, time = tt, evid = 0, amt = 0, cmt = 0,
                    PRED_RATE = pred_rate(tt - inp$ster_day, inp$ster_mgkg * BW,
                                          inp$ster_plateau, inp$ster_taper),
                    JAK_RATE = ifelse(inp$rescue == "tofacitinib" &
                                        tt >= inp$resc_day &
                                        tt < inp$resc_day + 60, 20, 0),
                    ABX_ON = ifelse(tt < inp$abx_days, 1, 0))

  ev <- data.frame()
  add <- function(times, amt, cmt)
    data.frame(ID = 1, time = times, amt = amt, cmt = cmt, evid = 1)

  if (inp$ipi_mgkg > 0)
    ev <- bind_rows(ev, add(seq(0, by = 21, length.out = inp$ipi_n),
                            inp$ipi_mgkg * BW, 1))
  if (inp$pd1_mgkg > 0)
    ev <- bind_rows(ev, add(seq(0, by = 14, length.out = inp$pd1_n),
                            inp$pd1_mgkg * BW, 3))
  if (inp$rescue == "infliximab")
    ev <- bind_rows(ev, add(inp$resc_day + c(0, 14, 42), 5 * BW, 5))
  if (inp$rescue == "vedolizumab")
    ev <- bind_rows(ev, add(inp$resc_day + c(0, 14, 42), 300, 7))
  if (inp$rescue == "tocilizumab")
    ev <- bind_rows(ev, add(inp$resc_day + c(0, 28), 8 * BW, 9))

  if (nrow(ev)) {
    ev$PRED_RATE <- NA; ev$JAK_RATE <- NA; ev$ABX_ON <- NA
    out <- bind_rows(cov, ev) %>% arrange(time, desc(evid)) %>%
      fill(PRED_RATE, JAK_RATE, ABX_ON, .direction = "downup")
  } else out <- cov
  out
}

run_sim <- function(inp) {
  ph <- PHENO[[inp$pheno]]
  m <- mod %>% param(ph$par) %>% init(ph$ini)
  if (inp$fmt_day > 0) {
    ## FMT is a reset of the antigen/SCFA arm, so it enters as a step change
    ## in Dv rather than as a drug.
    d1 <- build_data(inp) %>% filter(time <= inp$fmt_day)
    r1 <- m %>% data_set(d1) %>% mrgsim(end = inp$fmt_day, delta = 0.5)
    st <- as.list(tail(as_tibble(r1), 1))
    st$Dv <- min(1, st$Dv + 0.60)
    keep <- names(init(m))
    m2 <- m %>% init(st[keep])
    d2 <- build_data(inp) %>% filter(time > inp$fmt_day) %>%
      mutate(time = time - inp$fmt_day)
    r2 <- m2 %>% data_set(d2) %>% mrgsim(end = inp$tend - inp$fmt_day,
                                         delta = 0.5)
    return(bind_rows(as_tibble(r1),
                     as_tibble(r2) %>% mutate(time = time + inp$fmt_day)))
  }
  m %>% data_set(build_data(inp)) %>% mrgsim(end = inp$tend, delta = 0.5) %>%
    as_tibble()
}

THEME <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom", legend.title = element_blank(),
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(colour = "grey35", size = 10))

## ===========================================================================
## UI — eight tabs
## ===========================================================================
ui <- fluidPage(
  titlePanel(HTML(paste0(
    "<b>Immune Checkpoint Inhibitor-Induced Colitis QSP Dashboard</b>",
    "<br><span style='font-size:14px;color:#555'>Immune Checkpoint Inhibitor",
    "-Induced Colitis — a ratio with a saturated numerator and a depletable",
    " denominator</span>"))),
  tags$hr(),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectInput("pheno", "Patient phenotype",
                  choices = names(PHENO), selected = names(PHENO)[2]),
      tags$hr(),
      tags$b("Checkpoint inhibitor"),
      sliderInput("ipi_mgkg", "ipilimumab (mg/kg q3w)", 0, 10, 3, step = 0.5),
      sliderInput("ipi_n", "Number of ipilimumab doses", 0, 8, 4, step = 1),
      sliderInput("pd1_mgkg", "anti-PD-1 (mg/kg q2w)", 0, 10, 0, step = 0.5),
      sliderInput("pd1_n", "Number of anti-PD-1 doses", 0, 20, 13, step = 1),
      helpText(HTML(paste0("<i>Even moving the anti-PD-1 slider all the way, the occupancy curve",
                           " on tab 3 barely changes. That is the flat",
                           " dose-toxicity relationship.</i>"))),
      tags$hr(),
      tags$b("Rescue therapy"),
      sliderInput("ster_day", "Steroid start day", 0, 180, 53, step = 1),
      sliderInput("ster_mgkg", "prednisolone (mg/kg/day)", 0, 2, 1, step = 0.1),
      sliderInput("ster_plateau", "Days at fixed dose", 0, 21, 7, step = 1),
      sliderInput("ster_taper", "Taper duration (weeks)", 1, 12, 4, step = 1),
      radioButtons("rescue", "Additional rescue agent",
                   c("none", "infliximab", "vedolizumab", "tocilizumab",
                     "tofacitinib"), selected = "none"),
      sliderInput("resc_day", "Rescue agent start day", 0, 200, 60, step = 1),
      sliderInput("fmt_day", "FMT day (0 = not performed)", 0, 200, 0,
                  step = 1),
      tags$hr(),
      sliderInput("abx_days", "Pre-treatment antibiotic exposure (days)", 0, 30, 0, step = 1),
      sliderInput("tend", "Observation period (days)", 90, 365, 210, step = 15),
      tags$hr(),
      helpText(HTML(paste0("<span style='font-size:11px'>This is a model for educational/research ",
                           "purposes. Do not use it for clinical decision-making.",
                           "</span>")))
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        type = "tabs",
        ## ---- 1 -----------------------------------------------------------
        tabPanel(
          "1. Patient Profile",
          br(),
          fluidRow(column(12, uiOutput("verdict"))),
          br(),
          fluidRow(column(6, plotOutput("p_ratio", height = 300)),
                   column(6, plotOutput("p_cells", height = 300))),
          br(),
          fluidRow(column(12, DTOutput("t_summary")))
        ),
        ## ---- 2 -----------------------------------------------------------
        tabPanel(
          "2. Pharmacokinetics (PK)",
          br(),
          plotOutput("p_pk", height = 360),
          br(),
          plotOutput("p_clifx", height = 300),
          helpText(HTML(paste0("<b>The disease eats its own antidote.</b> ",
                               "Protein-losing enteropathy from severe colitis drops albumin,",
                               " and albumin is the biggest single covariate of monoclonal antibody clearance",
                               " (CL ∝ ALB<sup>-0.9</sup>). The same 5 mg/kg delivers",
                               " the lowest exposure to the sickest patient.")))
        ),
        ## ---- 3 -----------------------------------------------------------
        tabPanel(
          "3. Target Engagement — Two Different Exposure-Response Shapes",
          br(),
          plotOutput("p_occ", height = 340),
          br(),
          fluidRow(column(6, plotOutput("p_doseresp", height = 320)),
                   column(6, plotOutput("p_G", height = 320))),
          helpText(HTML(paste0("PD-1 occupancy is 0.98-0.999 across the entire clinical dose",
                               " range — there is no dose-response to find. Because ADCC does not saturate,",
                               " ipilimumab's dose-response comes entirely from the <b>denominator</b>",
                               ".")))
        ),
        ## ---- 4 -----------------------------------------------------------
        tabPanel(
          "4. Absorptive Reserve and Censored Symptoms",
          br(),
          plotOutput("p_reserve", height = 380),
          br(),
          plotOutput("p_stool", height = 300),
          helpText(HTML(paste0("<b>S* = L/A_max = 0.389.</b> 61% of colonic absorptive function can ",
                               "be lost alongside a completely normal bowel-movement",
                               " record. Grade is a late, high-loss trigger, and symptom",
                               "-driven treatment is inherently late.")))
        ),
        ## ---- 5 -----------------------------------------------------------
        tabPanel(
          "5. Clinical Endpoints",
          br(),
          plotOutput("p_grade", height = 320),
          br(),
          fluidRow(column(6, plotOutput("p_epi", height = 300)),
                   column(6, plotOutput("p_ulcer", height = 300))),
          br(),
          DTOutput("t_endpoints")
        ),
        ## ---- 6 -----------------------------------------------------------
        tabPanel(
          "6. Biomarkers",
          br(),
          plotOutput("p_biomark", height = 360),
          br(),
          plotOutput("p_cyto", height = 320),
          helpText(HTML(paste0("Calprotectin and endoscopy read the <b>lesion</b>, so they are ",
                               "not censored. That is why they move before symptoms.")))
        ),
        ## ---- 7 -----------------------------------------------------------
        tabPanel(
          "7. Tumour — Selectivity Is the Currency",
          br(),
          plotOutput("p_tumor", height = 340),
          br(),
          fluidRow(column(6, plotOutput("p_imm", height = 300)),
                   column(6, plotOutput("p_haz", height = 300))),
          helpText(HTML(paste0("χ (non-gut-selectivity): prednisolone 1.00 · ",
                               "infliximab 0.60 · vedolizumab 0.10. ",
                               "α4β7:MAdCAM-1 is a <b>gut address</b> that the tumour does not use,",
                               " so the same colitis control can be bought much more cheaply.")))
        ),
        ## ---- 8 -----------------------------------------------------------
        tabPanel(
          "8. Scenario Comparison",
          br(),
          checkboxGroupInput(
            "cmp", "Scenarios to compare",
            choices = c("ipilimumab 1 mg/kg", "ipilimumab 3 mg/kg",
                        "ipilimumab 10 mg/kg", "anti-PD-1 3 mg/kg",
                        "anti-PD-1 10 mg/kg", "ipi 3 + nivo 1"),
            selected = c("ipilimumab 1 mg/kg", "ipilimumab 3 mg/kg",
                         "ipilimumab 10 mg/kg", "anti-PD-1 3 mg/kg"),
            inline = TRUE),
          plotOutput("p_cmp", height = 420),
          br(),
          DTOutput("t_cmp"),
          helpText(HTML(paste0("Check two things. (1) anti-PD-1 3 vs 10 ",
                               "mg/kg essentially overlap. (2) The combination is not the sum of the two monotherapies",
                               " — because grade is a threshold indicator, and ipilimumab has already",
                               " pushed the already-vulnerable patient past that threshold.")))
        )
      )
    )
  )
)

## ===========================================================================
## SERVER
## ===========================================================================
server <- function(input, output, session) {

  sim <- reactive({
    inp <- reactiveValuesToList(input)
    run_sim(inp)
  })

  ## ---- 1. verdict banner --------------------------------------------------
  output$verdict <- renderUI({
    s <- sim()
    pg <- max(s$grade_o); mn <- min(s$S_eff_o)
    d2 <- s$time[which(s$grade_o >= 2)][1]
    dc <- s$time[which(s$Calpro >= 200)][1]
    col <- c("#2e7d32", "#f9a825", "#ef6c00", "#c62828")[pg + 1]
    HTML(sprintf(paste0(
      "<div style='padding:14px;border-left:6px solid %s;background:#fafafa'>",
      "<span style='font-size:20px;font-weight:bold;color:%s'>Peak grade G%d",
      "</span><br><span style='color:#444'>Lowest absorptive function S_eff = %.3f ",
      "(threshold S* = %.3f) · remaining reserve %.0f%%<br>",
      "Calprotectin &gt;200 µg/g: day %s · grade 2 reached: day %s · <b>lead time %s days</b>",
      "<br>Lowest albumin %.2f g/dL · maximum ulcer index %.2f · cumulative steroid %.0f mg",
      "</span></div>"),
      col, col, pg, mn, s$S_star_o[1], 100 * max(0, mn - s$S_star_o[1]) /
        s$S_star_o[1],
      ifelse(is.na(dc), "—", format(dc)), ifelse(is.na(d2), "—", format(d2)),
      ifelse(is.na(d2) || is.na(dc), "—", format(d2 - dc)),
      min(s$Alb), max(s$Ulcer), max(s$SterCum)))
  })

  output$p_ratio <- renderPlot({
    s <- sim()
    s %>% select(time, `G = Reg0/Reg (regulatory release)` = G_o,
                 `Ppd1 (PD-1 release)` = Ppd1_o) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      geom_hline(yintercept = 1, linetype = 3) +
      labs(x = "Day", y = "Fold (baseline = 1)",
           title = "The two factors of the ratio",
           subtitle = "The denominator depletes (G rises hyperbolically), the numerator saturates") +
      THEME
  })

  output$p_cells <- renderPlot({
    s <- sim()
    s %>% select(time, Teff, Trm, Treg, Nclone) %>% pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      labs(x = "Day", y = "Relative amount", title = "Cell compartments",
           subtitle = "Trm outlasts the drug — it is a state, not a trigger") +
      THEME
  })

  output$t_summary <- renderDT({
    s <- sim()
    data.frame(
      Metric = c("Peak CTCAE diarrhoea grade", "Lowest S_eff", "Censoring threshold S*",
               "Peak calprotectin (µg/g)", "Peak CRP (mg/L)",
               "Lowest albumin (g/dL)", "Maximum ulcer index", "Lowest crypt (ISC) index",
               "d180 Trm", "Cumulative steroid (mg)", "Cumulative infection risk",
               "d180 tumour burden"),
      Value = c(max(s$grade_o), round(min(s$S_eff_o), 3), round(s$S_star_o[1], 3),
             round(max(s$Calpro)), round(max(s$CRP), 1), round(min(s$Alb), 2),
             round(max(s$Ulcer), 2), round(min(s$ISC), 3),
             round(s$Trm[nrow(s)], 3), round(max(s$SterCum)),
             round(max(s$InfHaz), 3), round(s$Tumor[nrow(s)], 1))) %>%
      datatable(rownames = FALSE, options = list(dom = "t", pageLength = 20))
  })

  ## ---- 2. PK --------------------------------------------------------------
  output$p_pk <- renderPlot({
    s <- sim()
    s %>% select(time, ipilimumab = C_ipi_o, `anti-PD-1` = C_pd1_o,
                 infliximab = C_ifx_o, vedolizumab = C_vdz_o) %>%
      pivot_longer(-time) %>% filter(value > 1e-6) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      labs(x = "Day", y = "Plasma concentration (µg/mL)", title = "Monoclonal antibody pharmacokinetics") +
      THEME
  })

  output$p_clifx <- renderPlot({
    s <- sim()
    ggplot(s, aes(time)) +
      geom_line(aes(y = CL_ifx_o, colour = "infliximab CL (L/day)"),
                linewidth = 0.9) +
      geom_line(aes(y = Alb / 10, colour = "Albumin / 10 (g/dL)"),
                linewidth = 0.9) +
      labs(x = "Day", y = NULL,
           title = "Disease severity raises the clearance of its own antidote",
           subtitle = "CL_ifx = 0.40 · (ALB/4.0)^-0.9 · (1 + 0.006·CRP)") +
      THEME
  })

  ## ---- 3. target engagement ----------------------------------------------
  output$p_occ <- renderPlot({
    s <- sim()
    s %>% select(time, `PD-1 occupancy (saturated)` = O_PD1_o,
                 `CTLA-4 occupancy (saturated)` = O_CTLA4_o,
                 `ADCC drive (unsaturated)` = A_adcc_o) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 1) +
      ylim(0, 1) +
      labs(x = "Day", y = "Fraction of maximum",
           title = "Why one drug has a dose-response and the other does not",
           subtitle = paste0("Occupancy is saturated at every clinical dose. ",
                             "Only ADCC does not saturate.")) +
      THEME
  })

  output$p_doseresp <- renderPlot({
    ## computed analytically from the same equations, not simulated
    p <- as.list(param(mod))
    d <- seq(0.1, 10, length.out = 120)
    pd1 <- d * BW / (p$CL_pd1 * 14); ipi <- d * BW / (p$CL_ipi * 21)
    tibble(
      dose = rep(d, 2),
      value = c((pd1 * p$BDC * 1e3 / p$MW_pd1) /
                  (pd1 * p$BDC * 1e3 / p$MW_pd1 + p$KD_PD1_nM),
                (ipi * p$BDC) / (ipi * p$BDC + p$EC50_ADCC)),
      what = rep(c("anti-PD-1 occupancy", "ipilimumab ADCC drive"), each = 120)
    ) %>% ggplot(aes(dose, value, colour = what)) +
      geom_line(linewidth = 1) + ylim(0, 1) +
      scale_x_log10(breaks = c(0.1, 0.3, 1, 3, 10)) +
      labs(x = "Dose (mg/kg, log)", y = "Fraction of maximum",
           title = "Two exposure-response shapes",
           subtitle = "One is flat across the entire 30-fold dose range") +
      THEME
  })

  output$p_G <- renderPlot({
    p <- as.list(param(mod))
    treg <- seq(0.02, 1, length.out = 200)
    tibble(treg, G = (p$a_reg * 1 + p$kappa) / (p$a_reg * treg + p$kappa)) %>%
      ggplot(aes(treg, G)) + geom_line(linewidth = 1, colour = "#c0392b") +
      geom_hline(yintercept = (1 + p$kappa) / p$kappa, linetype = 2) +
      annotate("text", x = 0.55, y = (1 + p$kappa) / p$kappa * 0.93, size = 3.4,
               label = sprintf("Ceiling G_max = Reg0/kappa = %.1f\n(set by the non-Treg floor kappa)",
                               (1 + p$kappa) / p$kappa)) +
      labs(x = "Residual colonic Treg (vs baseline)", y = "Regulatory-release coefficient G",
           title = "The depletable denominator is a hyperbola",
           subtitle = "G keeps rising even after every occupancy term has flattened") +
      THEME
  })

  ## ---- 4. reserve ---------------------------------------------------------
  output$p_reserve <- renderPlot({
    s <- sim()
    ggplot(s, aes(time)) +
      geom_ribbon(aes(ymin = S_star_o, ymax = pmax(S_eff_o, S_star_o)),
                  fill = "#d4f0ef", alpha = 0.65) +
      geom_line(aes(y = S_eff_o), linewidth = 1, colour = "#17807a") +
      geom_hline(aes(yintercept = S_star_o), linetype = 2) +
      annotate("text", x = max(s$time) * 0.02, y = s$S_star_o[1] + 0.04,
               hjust = 0, size = 3.6,
               label = "S* = L/A_max = 0.389 — only once it falls below this does the bowel-movement record move") +
      ylim(0, 1.02) +
      labs(x = "Day", y = "Effective absorptive function S_eff",
           title = "Reserve exhausted before the first symptom",
           subtitle = "The shaded region is remaining reserve. Symptoms are zero and the lesion is progressing") +
      THEME
  })

  output$p_stool <- renderPlot({
    s <- sim()
    ggplot(s, aes(time)) +
      geom_line(aes(y = Stool_o, colour = "Stool water content (mL/day)"), linewidth = 0.9) +
      geom_line(aes(y = Calpro, colour = "Calprotectin (µg/g)"), linewidth = 0.9) +
      geom_hline(yintercept = 200, linetype = 3) +
      labs(x = "Day", y = NULL,
           title = "A censored trigger and an uncensored trigger",
           subtitle = "Calprotectin moves first. That is the basis for early intervention") +
      THEME
  })

  ## ---- 5. endpoints -------------------------------------------------------
  output$p_grade <- renderPlot({
    s <- sim()
    ggplot(s, aes(time)) +
      geom_step(aes(y = grade_o), linewidth = 1, colour = "#7e5109") +
      geom_line(aes(y = dfreq_o), linewidth = 0.7, colour = "#b9770e",
                linetype = 2) +
      scale_y_continuous("CTCAE grade (solid) · bowel movements/day increase (dashed)",
                         breaks = 0:12) +
      labs(x = "Day", title = "Clinical endpoints",
           subtitle = "G1 <4 · G2 4-6 · G3 ≥7 (increase over baseline)") +
      THEME
  })

  output$p_epi <- renderPlot({
    s <- sim()
    s %>% select(time, `Colonocyte mass Ent` = Ent, `Tight junction TJ` = TJ,
                 `Crypt stem cells ISC` = ISC, `Mucus Muc` = Muc) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      ylim(0, 1.05) +
      labs(x = "Day", y = "vs baseline", title = "Epithelial compartments",
           subtitle = "ISC recovery has a 20-day time constant — the thing steroids cannot reverse") +
      THEME
  })

  output$p_ulcer <- renderPlot({
    s <- sim()
    ggplot(s, aes(time)) +
      geom_line(aes(y = Ulcer, colour = "Ulcer index"), linewidth = 0.9) +
      geom_line(aes(y = Alb, colour = "Albumin (g/dL)"), linewidth = 0.9) +
      labs(x = "Day", y = NULL, title = "Ulceration and protein loss",
           subtitle = "Deep ulceration predicts steroid-refractoriness") +
      THEME
  })

  output$t_endpoints <- renderDT({
    s <- sim()
    first_at <- function(v, thr) {
      i <- which(v >= thr)[1]; if (is.na(i)) NA else s$time[i]
    }
    data.frame(
      Trigger = c("Lesion (S_eff below 90% of baseline)", "Calprotectin >150 µg/g",
               "Calprotectin >200 µg/g", "CTCAE grade 1", "CTCAE grade 2",
               "CTCAE grade 3"),
      FirstDay = c(s$time[which(s$S_eff_o < 0.9 * s$S_eff_o[1])[1]],
                 first_at(s$Calpro, 150), first_at(s$Calpro, 200),
                 first_at(s$grade_o, 1), first_at(s$grade_o, 2),
                 first_at(s$grade_o, 3))) %>%
      datatable(rownames = FALSE, options = list(dom = "t"),
                caption = "The censoring ladder: read left to right, and this is reserve being exhausted")
  })

  ## ---- 6. biomarkers ------------------------------------------------------
  output$p_biomark <- renderPlot({
    s <- sim()
    s %>% select(time, `Calprotectin (µg/g)` = Calpro, `CRP (mg/L)` = CRP,
                 `Albumin ×100 (g/dL)` = Alb) %>%
      mutate(`Albumin ×100 (g/dL)` = `Albumin ×100 (g/dL)` * 100) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "Day", y = NULL, title = "Biomarkers") +
      THEME + theme(legend.position = "none")
  })

  output$p_cyto <- renderPlot({
    s <- sim()
    s %>% select(time, `IFN-γ` = IFNg, `TNF-α` = TNFa, `IL-6` = IL6,
                 `IL-15 (Trm fuel)` = IL15, `CXCL10` = CXCL10,
                 `MAdCAM-1` = MADCAM) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.85) +
      labs(x = "Day", y = "Fold vs baseline", title = "Cytokine network",
           subtitle = "IFN-γ sets the programme, TNF-α executes it") +
      THEME
  })

  ## ---- 7. tumour ----------------------------------------------------------
  output$p_tumor <- renderPlot({
    s <- sim()
    ggplot(s, aes(time)) +
      geom_line(aes(y = Tumor, colour = "Tumour burden"), linewidth = 1) +
      geom_line(aes(y = TumTeff * 200, colour = "Intratumoural effector ×200"),
                linewidth = 0.85) +
      labs(x = "Day", y = NULL, title = "The tumour cost of rescue therapy",
           subtitle = "The clone that damages the colon is the clone that kills the tumour") +
      THEME
  })

  output$p_imm <- renderPlot({
    s <- sim()
    s %>% select(time, `Systemic immunosuppression (χ-weighted)` = Imm_sys_o,
                 `E_pred` = E_pred_o, `E_ifx` = E_ifx_o,
                 `E_vdz (gut-selective)` = E_vdz_o) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      ylim(0, 1) +
      labs(x = "Day", y = "Effect fraction", title = "Composition of immunosuppression",
           subtitle = "vedolizumab draws on the systemic account at only χ=0.10") +
      THEME
  })

  output$p_haz <- renderPlot({
    s <- sim()
    ggplot(s, aes(time)) +
      geom_line(aes(y = InfHaz, colour = "Cumulative infection risk"), linewidth = 0.9) +
      geom_line(aes(y = SterCum / 20000, colour = "Cumulative steroid / 20000"),
                linewidth = 0.9) +
      labs(x = "Day", y = NULL, title = "The price of immunosuppression") + THEME
  })

  ## ---- 8. scenario comparison --------------------------------------------
  cmp_sim <- reactive({
    inp <- reactiveValuesToList(input)
    specs <- list(
      "ipilimumab 1 mg/kg"  = list(ipi_mgkg = 1,  pd1_mgkg = 0),
      "ipilimumab 3 mg/kg"  = list(ipi_mgkg = 3,  pd1_mgkg = 0),
      "ipilimumab 10 mg/kg" = list(ipi_mgkg = 10, pd1_mgkg = 0),
      "anti-PD-1 3 mg/kg"   = list(ipi_mgkg = 0,  pd1_mgkg = 3),
      "anti-PD-1 10 mg/kg"  = list(ipi_mgkg = 0,  pd1_mgkg = 10),
      "ipi 3 + nivo 1"      = list(ipi_mgkg = 3,  pd1_mgkg = 1)
    )
    bind_rows(lapply(input$cmp, function(nm) {
      i2 <- modifyList(inp, specs[[nm]])
      run_sim(i2) %>% mutate(arm = nm)
    }))
  })

  output$p_cmp <- renderPlot({
    s <- cmp_sim(); req(nrow(s) > 0)
    ggplot(s, aes(time, S_eff_o, colour = arm)) + geom_line(linewidth = 0.9) +
      geom_hline(aes(yintercept = S_star_o), linetype = 2) + ylim(0, 1.02) +
      labs(x = "Day", y = "Effective absorptive function S_eff",
           title = "Scenario comparison — does it cross the threshold, or not",
           subtitle = "anti-PD-1 3 vs 10 mg/kg overlap. That is the point") +
      THEME
  })

  output$t_cmp <- renderDT({
    s <- cmp_sim(); req(nrow(s) > 0)
    s %>% group_by(arm) %>%
      summarise(`Peak grade` = max(grade_o),
                `Lowest S_eff` = round(min(S_eff_o), 3),
                `Day grade 2 reached` = suppressWarnings(min(time[grade_o >= 2])),
                `Peak calprotectin` = round(max(Calpro)),
                `Lowest Treg` = round(min(Treg), 3),
                `Peak G` = round(max(G_o), 2),
                `Peak ADCC drive` = round(max(A_adcc_o), 3),
                `d180 tumour` = round(last(Tumor), 1), .groups = "drop") %>%
      mutate(`Day grade 2 reached` = ifelse(is.infinite(`Day grade 2 reached`), NA,
                                     `Day grade 2 reached`)) %>%
      datatable(rownames = FALSE, options = list(dom = "t"))
  })
}

shinyApp(ui, server)
