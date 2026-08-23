## =====================================================================
##  Cancer Anorexia-Cachexia Syndrome (CACS) — QSP Shiny dashboard
##  cacs_shiny_app.R
## ---------------------------------------------------------------------
##  Front end for cacs_mrgsolve_model.R.
##
##  The dashboard is organised around the model's central claim: that mass
##  and function are separate endpoints driven by separate arms, and that
##  a drug's position in the mechanism determines which one it can move.
##  Every tab is there to make one part of that argument inspectable.
##
##    1  Patient profile      Patient       tumour, host, staging, the driver
##    2  Pharmacokinetics     PK            all eleven agents, one panel
##    3  ARM A Appetite       Intake        GDF-15 -> GFRAL -> melanocortin -> kcal
##    4  ARM B Catabolism     Catabolism    STAT3/SMAD -> FoxO -> UPS, and REE
##    5  AXIS C Muscle quality Quality      the axis mass-only drugs never touch
##    6  Body composition     Composition   muscle, fat, water, and what DXA sees
##    7  Clinical endpoints   Endpoints     weight, grip, FAACT, ECOG, staging
##    8  Scenario comparison  Compare       run several arms side by side
##    9  Biomarkers           Biomarkers    CRP, albumin, mGPS, IGF-1, cytokines
##   10  Dissociation & hysteresis Dissociation the LBM-vs-function grid, and timing
##
##  Run:
##    setwd("cancer-cachexia"); shiny::runApp("cacs_shiny_app.R")
##  Requires: shiny, mrgsolve, ggplot2, dplyr, tidyr, DT
## =====================================================================

suppressPackageStartupMessages({
  library(shiny)
  library(mrgsolve)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(DT)
})

## ---- load the model without executing its demonstration block ---------
model_env <- new.env()
local({
  src <- readLines("cacs_mrgsolve_model.R")
  cut <- grep("^if \\(identical\\(environment\\(\\), globalenv\\(\\)\\)\\)", src)
  if (length(cut)) src <- src[seq_len(cut[1] - 1)]
  eval(parse(text = paste(src, collapse = "\n")), envir = model_env)
})
mod        <- model_env$mod
CMT        <- model_env$CMT
pat_nsclc  <- model_env$pat_nsclc
pat_panc   <- model_env$pat_panc
pat_trial  <- model_env$pat_trial
pat_gdfhi  <- model_env$pat_gdfhi
pat_early  <- model_env$pat_early
pat_refr   <- model_env$pat_refr

THEME <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")

PAL <- c("#2f5f9e", "#b04a4a", "#2f7a3f", "#e07b1f", "#7d3cb0",
         "#c2568b", "#3a7f86", "#8a6b2f", "#5c6670", "#a03c6a")

DRUGS <- list(
  anamorelin    = list(label = "Anamorelin 100 mg PO qd",        cmt = "ANAD", dose = 100,  F1 = 0.35, ii = 1),
  megestrol     = list(label = "Megestrol acetate 800 mg/d",     cmt = "MEGD", dose = 800,  F1 = 0.30, ii = 1),
  dexamethasone = list(label = "Dexamethasone 4 mg/d",           cmt = "DEXD", dose = 4,    F1 = 0.80, ii = 1),
  olanzapine    = list(label = "Olanzapine 2.5 mg qhs",          cmt = "OLZD", dose = 2.5,  F1 = 0.60, ii = 1),
  enobosarm     = list(label = "Enobosarm 3 mg PO qd",           cmt = "ENOD", dose = 3,    F1 = 0.60, ii = 1),
  espindolol    = list(label = "Espindolol 10 mg PO bid",        cmt = "ESPD", dose = 10,   F1 = 0.85, ii = 0.5),
  epa           = list(label = "EPA 2 g/d",                      cmt = "EPAD", dose = 2000, F1 = 0.90, ii = 1),
  celecoxib     = list(label = "Celecoxib 200 mg bid",           cmt = "CELD", dose = 200,  F1 = 0.40, ii = 0.5)
)

build_events <- function(sel, days, pons_mg, tcz_on, bim_on) {
  parts <- list()
  for (nm in sel) {
    d <- DRUGS[[nm]]
    if (is.null(d)) next
    parts[[length(parts) + 1]] <- ev(amt = d$dose * d$F1, cmt = CMT(d$cmt),
                                     ii = d$ii, addl = ceiling(days / d$ii) - 1, time = 0)
  }
  if (pons_mg > 0)
    parts[[length(parts) + 1]] <- ev(amt = pons_mg / 147000 * 1e6 * 0.65,
                                     cmt = CMT("POND"), ii = 28,
                                     addl = floor(days / 28), time = 0)
  if (isTRUE(tcz_on))
    parts[[length(parts) + 1]] <- ev(amt = 162 * 0.80, cmt = CMT("TCZD"),
                                     ii = 7, addl = floor(days / 7), time = 0)
  if (isTRUE(bim_on))
    parts[[length(parts) + 1]] <- ev(amt = 10 * 70, cmt = CMT("BIMC"),
                                     ii = 28, addl = floor(days / 28), time = 0)
  if (!length(parts)) return(NULL)
  Reduce(`+`, parts)
}

simulate <- function(patient, pars, events, days) {
  m <- param(mod, modifyList(patient, pars))
  tg <- c(seq(0, 14, 0.25), seq(15, days, 1))
  if (is.null(events)) mrgsim_df(m, end = days, add = tg, delta = 1)
  else mrgsim_df(m, events = events, end = days, add = tg, delta = 1)
}

at_day <- function(d, t) d[which.min(abs(d$time - t)), ]

long <- function(d, vars, labels = vars) {
  out <- d[, c("time", vars), drop = FALSE]
  out <- pivot_longer(out, -time, names_to = "var", values_to = "value")
  out$var <- factor(out$var, levels = vars, labels = labels)
  out
}

facet_plot <- function(d, vars, labels = vars, ncol = 3, ylab = NULL) {
  ggplot(long(d, vars, labels), aes(time, value)) +
    geom_line(linewidth = 0.8, colour = PAL[1]) +
    facet_wrap(~var, scales = "free_y", ncol = ncol) +
    labs(x = "Day", y = ylab) + THEME
}

## =====================================================================
##  UI
## =====================================================================
ui <- fluidPage(
  titlePanel("Cancer Anorexia-Cachexia Syndrome — QSP simulator (CACS)"),
  tags$p(style = "color:#555;margin-top:-10px",
         HTML("<b>ARM A</b> Appetite (GDF-15 &rarr; GFRAL &rarr; melanocortin) &middot; ",
              "<b>ARM B</b> Catabolism + hypermetabolism (STAT3/SMAD &rarr; FoxO) &middot; ",
              "<b>AXIS C</b> Muscle quality (force per kg) &nbsp;&mdash;&nbsp; ",
              "GRIP = MASS &times; QUALITY")),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectInput("patient", "Patient archetype",
                  c("Advanced non-small-cell lung cancer (NSCLC)"       = "nsclc",
                    "Pancreatic cancer, high cachexin drive (pancreatic)"     = "panc",
                    "Early precachexia (precachexia)"      = "early",
                    "Refractory cachexia (refractory)"       = "refr",
                    "Registration-trial population (trial-like)"         = "trial",
                    "GDF-15-high population (enriched)"      = "gdfhi",
                    "Healthy control (no tumour)"              = "healthy"),
                  selected = "nsclc"),
      sliderInput("days", "Simulation period (days)", 28, 540, 168, step = 14),
      hr(),
      h5("Tumour"),
      sliderInput("tvol", "Initial tumour burden (g)", 0, 400, 25, step = 5),
      sliderInput("kgrow", "Gompertz growth rate (1/d)", 0.001, 0.020, 0.006, step = 0.001),
      sliderInput("fcach", "Tissue-type cachexin drive", 0.3, 2.0, 1.0, step = 0.1),
      sliderInput("cxsens", "Host sensitivity, CXSENS", 0.5, 1.8, 1.0, step = 0.05),
      hr(),
      h5("Anticancer therapy"),
      sliderInput("onceff", "Tumour kill rate (1/d)", 0, 0.08, 0, step = 0.005),
      sliderInput("oncstart", "Treatment start day", 0, 400, 14, step = 7),
      sliderInput("myotox", "Direct myotoxicity (0-1)", 0, 1, 0, step = 0.05),
      hr(),
      h5("Drugs"),
      checkboxGroupInput("drugs", NULL,
                         choices = setNames(names(DRUGS), sapply(DRUGS, `[[`, "label"))),
      selectInput("pons", "Ponsegromab SC q4w",
                  c("None" = "0", "100 mg" = "100", "200 mg" = "200", "400 mg" = "400"),
                  selected = "0"),
      checkboxInput("tcz", "Tocilizumab 162 mg SC weekly", FALSE),
      checkboxInput("bim", "Bimagrumab 10 mg/kg IV q4w", FALSE),
      hr(),
      h5("Non-pharmacologic"),
      sliderInput("ons", "ONS supplemental calories (kcal/d)", 0, 900, 0, step = 100),
      sliderInput("prot", "Protein prescription (g/kg/d)", 0.6, 2.0, 1.0, step = 0.1),
      sliderInput("exres", "Resistance exercise (0-1)", 0, 1, 0, step = 0.1),
      sliderInput("exaer", "Aerobic exercise (0-1)", 0, 1, 0, step = 0.1),
      checkboxInput("leu", "Leucine/HMB-enriched bolus", FALSE),
      hr(),
      h5("Symptom burden"),
      sliderInput("gibar", "Mechanical/mucosal barrier, GIBAR", 0, 0.6, 0, step = 0.05)
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("1. Patient profile",
                 br(), htmlOutput("summary_box"), hr(),
                 plotOutput("p_driver", height = "420px"),
                 hr(), DTOutput("t_stage")),
        tabPanel("2. PK",
                 br(), plotOutput("p_pk", height = "560px"),
                 helpText("All concentrations are shown only for the selected drugs. ",
                          "Monoclonal-antibody clearance increases, via albumin-dependent FcRn ",
                          "recycling, in hypoalbuminaemia — the sickest patients get the least exposure.")),
        tabPanel("3. ARM A — Appetite/intake",
                 br(), plotOutput("p_armA", height = "560px"),
                 helpText("A fall in leptin should normally switch on AgRP, but central ",
                          "IL-1b/PGE2 blocks that compensation. This is the difference between starvation and cachexia.")),
        tabPanel("4. ARM B — Catabolism/hypermetabolism",
                 br(), plotOutput("p_armB", height = "560px"),
                 helpText("REE/pREE > 1.10 is hypermetabolism. That intake falls while ",
                          "expenditure rises is the defining paradox of this disease.")),
        tabPanel("5. AXIS C — Muscle quality",
                 br(), plotOutput("p_axisC", height = "520px"),
                 hr(), htmlOutput("qual_note")),
        tabPanel("6. Body composition",
                 br(), plotOutput("p_body", height = "520px"),
                 hr(), htmlOutput("dxa_note")),
        tabPanel("7. Clinical endpoints",
                 br(), plotOutput("p_end", height = "560px"),
                 hr(), DTOutput("t_end")),
        tabPanel("8. Scenario comparison",
                 br(),
                 checkboxGroupInput("cmp", "Scenarios to compare",
                                    choices = c("Untreated natural history" = "none",
                                                "Nutrition only" = "nut",
                                                "Nutrition + resistance exercise" = "nutex",
                                                "Anamorelin" = "ana",
                                                "Ponsegromab 400 mg" = "pons",
                                                "Megestrol" = "meg",
                                                "Enobosarm" = "eno",
                                                "Tocilizumab" = "tcz",
                                                "Espindolol" = "esp",
                                                "Multimodal" = "multi"),
                                    selected = c("none", "nut", "ana", "pons", "multi"),
                                    inline = TRUE),
                 plotOutput("p_cmp", height = "560px"),
                 hr(), DTOutput("t_cmp")),
        tabPanel("9. Biomarkers",
                 br(), plotOutput("p_bio", height = "560px"),
                 hr(), htmlOutput("tcz_note")),
        tabPanel("10. Dissociation & hysteresis",
                 br(), h4("Dissociation of lean body mass (LBM) from function"),
                 plotOutput("p_dissoc", height = "380px"),
                 DTOutput("t_dissoc"),
                 hr(), h4("The point of no return"),
                 helpText("Muscle mass recovered by day 540 when the same curative ",
                          "anticancer therapy is started progressively later."),
                 plotOutput("p_hyst", height = "340px"))
      )
    )
  )
)

## =====================================================================
##  SERVER
## =====================================================================
server <- function(input, output, session) {

  base_patient <- reactive({
    switch(input$patient,
           nsclc = pat_nsclc, panc = pat_panc, early = pat_early,
           refr = pat_refr, trial = pat_trial, gdfhi = pat_gdfhi,
           healthy = list(TVOL0 = 0, KGROW = 0.006, FCACHEX = 1, CXSENS = 1))
  })

  observeEvent(input$patient, {
    p <- base_patient()
    updateSliderInput(session, "tvol",   value = p$TVOL0)
    updateSliderInput(session, "kgrow",  value = p$KGROW)
    updateSliderInput(session, "fcach",  value = p$FCACHEX)
    updateSliderInput(session, "cxsens", value = p$CXSENS)
  })

  pars <- reactive({
    list(TVOL0 = input$tvol, KGROW = input$kgrow, FCACHEX = input$fcach,
         CXSENS = input$cxsens, ONCEFF = input$onceff,
         ONCSTART = if (input$onceff > 0) input$oncstart else 1e6,
         MYOTOX = input$myotox, ONSKCAL = input$ons, PROTTG = input$prot,
         EXRES = input$exres, EXAER = input$exaer,
         LEUBOL = as.numeric(isTRUE(input$leu)), GIBAR = input$gibar)
  })

  sim <- reactive({
    ev0 <- build_events(input$drugs, input$days,
                        as.numeric(input$pons), input$tcz, input$bim)
    simulate(list(), pars(), ev0, input$days)
  })

  ## ---- 1. patient profile --------------------------------------------
  output$summary_box <- renderUI({
    d <- sim(); a <- at_day(d, 0); b <- at_day(d, input$days)
    card <- function(lab, val, sub = "")
      sprintf("<div style='display:inline-block;min-width:150px;margin:4px 10px 4px 0;
               padding:8px 12px;border:1px solid #dcdcdc;border-radius:8px'>
               <div style='font-size:11px;color:#666'>%s</div>
               <div style='font-size:19px;font-weight:600'>%s</div>
               <div style='font-size:11px;color:#888'>%s</div></div>", lab, val, sub)
    stage_lab <- c("None", "Precachexia", "Cachexia", "Refractory")[b$STAGE + 1]
    HTML(paste0(
      card("Weight (kg)", sprintf("%.1f", b$BW), sprintf("Baseline %.1f", a$BW)),
      card("6-month weight loss", sprintf("%.1f %%", b$PCTWL), sprintf("BMI %.1f", b$BMI)),
      card("Lean body mass, LBM (kg)", sprintf("%.1f", b$LBM), sprintf("Skeletal muscle %.1f kg", b$MUSC)),
      card("Grip strength (kg)", sprintf("%.1f", b$GRIP), sprintf("Muscle quality %.2f", b$MYOQ)),
      card("Intake (kcal/d)", sprintf("%.0f", b$INTK), sprintf("TEE %.0f", b$TEE)),
      card("REE/pREE", sprintf("%.2f", b$RQREE),
           if (b$RQREE > 1.10) "Hypermetabolic" else "Normal range"),
      card("Fearon stage", stage_lab, sprintf("Martin grade %d", as.integer(b$MGRADE))),
      card("ECOG", sprintf("%.1f", b$ECOG), sprintf("FAACT %.0f/48", b$FAACT)),
      card("Survival probability", sprintf("%.0f %%", 100 * b$SURV), sprintf("Day %d", input$days))
    ))
  })

  output$p_driver <- renderPlot({
    facet_plot(sim(),
      c("TUMOR", "GDF15", "IL6", "ACTA", "CRP", "ALB"),
      c("Tumour burden (g)", "Total GDF-15 (ng/mL)", "IL-6 (pg/mL)",
        "Activin A (pg/mL)", "CRP (mg/L)", "Albumin (g/L)"))
  })

  output$t_stage <- renderDT({
    d <- sim()
    ts <- unique(pmin(c(0, 28, 56, 84, 112, 168, 252, 365, 540), input$days))
    tab <- do.call(rbind, lapply(ts, function(t) {
      b <- at_day(d, t)
      data.frame(Day = t, BW = round(b$BW, 1), `%WL` = round(b$PCTWL, 1),
                 BMI = round(b$BMI, 1), SMI = round(b$SMI, 1),
                 Sarcopenia = ifelse(b$SARCO > 0.5, "yes", "no"),
                 Fearon = c("none", "pre", "cachexia", "refractory")[b$STAGE + 1],
                 Martin = as.integer(b$MGRADE), mGPS = as.integer(b$MGPS),
                 ECOG = round(b$ECOG, 2), Survival = round(b$SURV, 3),
                 check.names = FALSE)
    }))
    datatable(tab, rownames = FALSE, options = list(dom = "t", pageLength = 20))
  })

  ## ---- 2. PK -----------------------------------------------------------
  output$p_pk <- renderPlot({
    d <- sim()
    v <- c(); l <- c()
    add <- function(cond, var, lab) { if (cond) { v <<- c(v, var); l <<- c(l, lab) } }
    add("anamorelin"    %in% input$drugs, "CANA", "Anamorelin (ng/mL)")
    add("megestrol"     %in% input$drugs, "CMEG", "Megestrol (ng/mL)")
    add("dexamethasone" %in% input$drugs, "CDEX", "Dexamethasone (ng/mL)")
    add("olanzapine"    %in% input$drugs, "COLZ", "Olanzapine (ng/mL)")
    add("enobosarm"     %in% input$drugs, "CENO", "Enobosarm (ng/mL)")
    add("espindolol"    %in% input$drugs, "CESP", "Espindolol (ng/mL)")
    add(as.numeric(input$pons) > 0,       "CPON", "Ponsegromab (nM)")
    add(isTRUE(input$tcz),                "CTCZ", "Tocilizumab (mg/L)")
    add(isTRUE(input$bim),                "CBIM", "Bimagrumab (mg/L)")
    if (!length(v))
      return(ggplot() + annotate("text", 0, 0, label = "Please select a drug", size = 6) +
               theme_void())
    facet_plot(d, v, l, ncol = 3, ylab = "Concentration")
  })

  ## ---- 3. ARM A ---------------------------------------------------------
  output$p_armA <- renderPlot({
    facet_plot(sim(),
      c("GDFFREE", "BSSIG", "NAUS", "AGRP", "POMC", "ANOR", "LEPT", "GHRL", "INTK"),
      c("Free GDF-15 (ng/mL)", "Brainstem aversion signal", "Nausea (0-10)",
        "AgRP/NPY tone", "POMC tone", "Integrated anorexia drive",
        "Leptin (ng/mL)", "Ghrelin (pg/mL)", "Energy intake (kcal/d)"))
  })

  ## ---- 4. ARM B ---------------------------------------------------------
  output$p_armB <- renderPlot({
    facet_plot(sim(),
      c("STAT3", "SMAD", "FOXO", "UPS", "AUTOP", "ARES",
        "IGF1", "CORT", "RQREE"),
      c("Muscle STAT3", "SMAD2/3", "FoxO activity",
        "Proteasome flux", "Autophagy flux", "Anabolic resistance index",
        "IGF-1 (ng/mL)", "Cortisol (ug/dL)", "REE / predicted REE"))
  })

  ## ---- 5. AXIS C --------------------------------------------------------
  output$p_axisC <- renderPlot({
    d <- sim()
    facet_plot(d, c("PGC1", "MITO", "ROS", "MYOQ", "SATC", "GRIP"),
               c("PGC-1alpha", "Mitochondrial content", "Oxidative stress (ROS)",
                 "Muscle quality, MYOQ (force/kg)", "Satellite-cell pool, SATC",
                 "Grip strength (kg) = quality x contractile mass"), ncol = 3)
  })

  output$qual_note <- renderUI({
    d <- sim(); b <- at_day(d, input$days)
    HTML(sprintf(
      "<div style='background:#f5f0fa;padding:12px;border-radius:8px'>
       <b>Mass DXA sees vs mass that pulls</b><br/>
       Of skeletal muscle <b>%.2f kg</b>, contractile mass is <b>%.2f kg</b>; GH-derived lean water is
       <b>%.2f kg</b>. Of what DXA counts as lean mass,
       <b>%.2f kg</b> generates no force.<br/>
       Muscle quality, MYOQ = <b>%.3f</b> — inflammation lowers this axis, and only loading (resistance exercise) raises it.
       Purely anabolic drugs cannot move this value, so they can hit the lean-mass endpoint yet
       miss the functional one.</div>",
      b$MUSC, b$MUSC - (b$NONCON - b$LNW), b$LNW, b$NONCON, b$MYOQ))
  })

  ## ---- 6. body composition ---------------------------------------------
  output$p_body <- renderPlot({
    d <- sim()
    ggplot(long(d, c("MUSC", "FATM", "OLBM", "LIVM", "LNW", "EDEM"),
                c("Skeletal muscle", "Fat", "Other lean mass", "Liver", "GH-derived lean water", "Oedema")),
           aes(time, value, fill = var)) +
      geom_area(alpha = 0.85) +
      scale_fill_manual(values = PAL) +
      labs(x = "Day", y = "kg", fill = NULL,
           title = "Body-composition breakdown — each stacked layer makes up the number on the scale") +
      THEME
  })

  output$dxa_note <- renderUI({
    d <- sim(); a <- at_day(d, 0); b <- at_day(d, input$days)
    HTML(sprintf(
      "<div style='background:#f0f4f8;padding:12px;border-radius:8px'>
       Weight change <b>%+.2f kg</b> = muscle <b>%+.2f</b> + fat <b>%+.2f</b> +
       other lean mass <b>%+.2f</b> + liver <b>%+.2f</b> + lean water <b>%+.2f</b> +
       oedema <b>%+.2f</b>.<br/>
       Hypoalbuminaemia and progestins expand body fluid, <b>hiding real tissue loss from the scale</b>.
       A stable weight does not mean stable tissue.</div>",
      b$BW - a$BW, b$MUSC - a$MUSC, b$FATM - a$FATM, b$OLBM - a$OLBM,
      b$LIVM - a$LIVM, b$LNW - a$LNW, b$EDEM - a$EDEM))
  })

  ## ---- 7. endpoints ------------------------------------------------------
  output$p_end <- renderPlot({
    facet_plot(sim(),
      c("BW", "PCTWL", "LBM", "GRIP", "SMI", "FAACT", "ECOG", "MGRADE", "SURV"),
      c("Weight (kg)", "6-month weight loss (%)", "Lean body mass (kg)", "Grip strength (kg)",
        "L3 SMI (cm2/m2)", "FAACT A/CS (0-48)", "ECOG PS", "Martin grade",
        "Survival probability"))
  })

  output$t_end <- renderDT({
    d <- sim(); a <- at_day(d, 0)
    ts <- unique(pmin(c(28, 56, 84, 112, 168, 252, 365), input$days))
    tab <- do.call(rbind, lapply(ts, function(t) {
      b <- at_day(d, t)
      data.frame(Day = t,
                 `dBW` = round(b$BW - a$BW, 2),
                 `dLBM` = round(b$LBM - a$LBM, 2),
                 `dMuscle` = round(b$MUSC - a$MUSC, 2),
                 `dFat` = round(b$FATM - a$FATM, 2),
                 `dGrip` = round(b$GRIP - a$GRIP, 2),
                 `dFAACT` = round(b$FAACT - a$FAACT, 1),
                 Intake = round(b$INTK), TEE = round(b$TEE),
                 `REE ratio` = round(b$RQREE, 3),
                 check.names = FALSE)
    }))
    datatable(tab, rownames = FALSE, options = list(dom = "t", pageLength = 20))
  })

  ## ---- 8. scenario comparison ---------------------------------------------
  cmp_defs <- list(
    none  = list(lab = "Untreated",             pars = list(), ev = function(dy) NULL),
    nut   = list(lab = "Nutrition only",          pars = list(ONSKCAL = 600, PROTTG = 1.5),
                 ev = function(dy) NULL),
    nutex = list(lab = "Nutrition + resistance exercise",    pars = list(ONSKCAL = 600, PROTTG = 1.5,
                                                        EXRES = 1, EXAER = 0.6, LEUBOL = 1),
                 ev = function(dy) NULL),
    ana   = list(lab = "Anamorelin",         pars = list(),
                 ev = function(dy) build_events("anamorelin", dy, 0, FALSE, FALSE)),
    pons  = list(lab = "Ponsegromab 400 mg", pars = list(),
                 ev = function(dy) build_events(character(0), dy, 400, FALSE, FALSE)),
    meg   = list(lab = "Megestrol",          pars = list(),
                 ev = function(dy) build_events("megestrol", dy, 0, FALSE, FALSE)),
    eno   = list(lab = "Enobosarm",          pars = list(),
                 ev = function(dy) build_events("enobosarm", dy, 0, FALSE, FALSE)),
    tcz   = list(lab = "Tocilizumab",        pars = list(),
                 ev = function(dy) build_events(character(0), dy, 0, TRUE, FALSE)),
    esp   = list(lab = "Espindolol",         pars = list(),
                 ev = function(dy) build_events("espindolol", dy, 0, FALSE, FALSE)),
    multi = list(lab = "Multimodal",           pars = list(ONSKCAL = 600, PROTTG = 1.5,
                                                        EXRES = 1, EXAER = 0.8, LEUBOL = 1),
                 ev = function(dy) build_events("epa", dy, 400, FALSE, FALSE))
  )

  cmp_runs <- reactive({
    req(length(input$cmp) > 0)
    base <- pars()
    do.call(rbind, lapply(input$cmp, function(k) {
      cd <- cmp_defs[[k]]
      o <- simulate(list(), modifyList(base, cd$pars), cd$ev(input$days), input$days)
      o$arm <- cd$lab
      o
    }))
  })

  output$p_cmp <- renderPlot({
    d <- cmp_runs()
    v <- c("BW", "LBM", "MUSC", "GRIP", "INTK", "FAACT", "MYOQ", "ECOG", "SURV")
    l <- c("Weight (kg)", "Lean body mass (kg)", "Skeletal muscle (kg)", "Grip strength (kg)",
           "Intake (kcal/d)", "FAACT", "Muscle quality, MYOQ", "ECOG", "Survival probability")
    dd <- pivot_longer(d[, c("time", "arm", v)], -c(time, arm),
                       names_to = "var", values_to = "value")
    dd$var <- factor(dd$var, levels = v, labels = l)
    ggplot(dd, aes(time, value, colour = arm)) +
      geom_line(linewidth = 0.85) +
      facet_wrap(~var, scales = "free_y", ncol = 3) +
      scale_colour_manual(values = PAL) +
      labs(x = "Day", y = NULL, colour = NULL) + THEME
  })

  output$t_cmp <- renderDT({
    d <- cmp_runs()
    tab <- do.call(rbind, lapply(split(d, d$arm), function(x) {
      a <- at_day(x, 0); b <- at_day(x, input$days)
      data.frame(Arm = x$arm[1],
                 dBW = round(b$BW - a$BW, 2), dLBM = round(b$LBM - a$LBM, 2),
                 dMuscle = round(b$MUSC - a$MUSC, 2), dFat = round(b$FATM - a$FATM, 2),
                 dGrip = round(b$GRIP - a$GRIP, 2), MYOQ = round(b$MYOQ, 3),
                 Intake = round(b$INTK), `REE ratio` = round(b$RQREE, 3),
                 CRP = round(b$CRP, 1), ECOG = round(b$ECOG, 2),
                 Survival = round(b$SURV, 3), check.names = FALSE)
    }))
    datatable(tab, rownames = FALSE, options = list(dom = "t", pageLength = 12))
  })

  ## ---- 9. biomarkers -------------------------------------------------------
  output$p_bio <- renderPlot({
    facet_plot(sim(),
      c("CRP", "ALB", "MGPS", "IL6", "TNFA", "MSTN", "IGF1", "TEST", "GDFTOT"),
      c("CRP (mg/L)", "Albumin (g/L)", "mGPS (0-2)", "IL-6 (pg/mL)",
        "TNF-alpha (pg/mL)", "Myostatin (ng/mL) - paradoxical decline",
        "IGF-1 (ng/mL)", "Testosterone (ng/dL)", "Total GDF-15 (ng/mL)"))
  })

  output$tcz_note <- renderUI({
    HTML("<div style='background:#fdf6ec;padding:12px;border-radius:8px'>
          <b>Two traps in this panel</b><br/>
          &bull; <b>Plasma IL-6 rises with tocilizumab.</b> Blocking the IL-6 receptor
          removes receptor-mediated clearance, so ligand accumulates. The rise in IL-6 is not
          treatment failure but evidence of target engagement; what should actually be judged is CRP.<br/>
          &bull; <b>Myostatin falls in cachexia.</b> Because muscle itself is the source
          of it, this is a consequence of muscle wasting, not its cause, and must not be misread as a biomarker.
          The ActRIIB ligand that drives catabolism is activin A.</div>")
  })

  ## ---- 10. dissociation and hysteresis --------------------------------------
  dissoc <- reactive({
    base <- pars()
    arms <- list(
      list(k = "Anamorelin",   p = list(), e = function(dy) build_events("anamorelin", dy, 0, FALSE, FALSE)),
      list(k = "Enobosarm",    p = list(), e = function(dy) build_events("enobosarm", dy, 0, FALSE, FALSE)),
      list(k = "Bimagrumab",   p = list(), e = function(dy) build_events(character(0), dy, 0, FALSE, TRUE)),
      list(k = "Ponsegromab",  p = list(), e = function(dy) build_events(character(0), dy, 400, FALSE, FALSE)),
      list(k = "Megestrol",    p = list(), e = function(dy) build_events("megestrol", dy, 0, FALSE, FALSE)),
      list(k = "Espindolol",   p = list(), e = function(dy) build_events("espindolol", dy, 0, FALSE, FALSE)),
      list(k = "Resistance exercise",     p = list(EXRES = 1, EXAER = 0.6), e = function(dy) NULL),
      list(k = "Multimodal",     p = list(ONSKCAL = 600, PROTTG = 1.5, EXRES = 1, EXAER = 0.8),
           e = function(dy) build_events("epa", dy, 400, FALSE, FALSE))
    )
    ref <- simulate(list(), base, NULL, 84)
    r0 <- at_day(ref, 0); r8 <- at_day(ref, 84)
    do.call(rbind, lapply(arms, function(a) {
      o <- simulate(list(), modifyList(base, a$p), a$e(84), 84)
      x0 <- at_day(o, 0); x8 <- at_day(o, 84)
      dl <- (x8$LBM - x0$LBM) - (r8$LBM - r0$LBM)
      dg <- (x8$GRIP - x0$GRIP) - (r8$GRIP - r0$GRIP)
      data.frame(Intervention = a$k, dLBM = round(dl, 2), dGrip = round(dg, 2),
                 `Grip per kg LBM` = ifelse(abs(dl) < 0.05, NA, round(dg / dl, 2)),
                 `% of pure-mass expectation` =
                   ifelse(abs(dl) < 0.05, NA, round(100 * (dg / dl) / (40 / 25))),
                 check.names = FALSE)
    }))
  })

  output$p_dissoc <- renderPlot({
    d <- dissoc()
    ggplot(d, aes(dLBM, dGrip, label = Intervention)) +
      geom_abline(slope = 40 / 25, intercept = 0, linetype = "dashed", colour = "#999") +
      annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.6, size = 3.4,
               colour = "#777", label = "Dashed line = expected line if all added mass were normal muscle") +
      geom_hline(yintercept = 0, colour = "#ccc") +
      geom_vline(xintercept = 0, colour = "#ccc") +
      geom_point(size = 3.4, colour = PAL[1]) +
      geom_text(vjust = -0.9, size = 3.6) +
      labs(x = "12-week change in lean body mass, vs untreated (kg)",
           y = "12-week change in grip strength, vs untreated (kg)") +
      expand_limits(x = c(-1.5, 4), y = c(-2, 9)) + THEME
  })

  output$t_dissoc <- renderDT(
    datatable(dissoc(), rownames = FALSE, options = list(dom = "t", pageLength = 10)))

  output$p_hyst <- renderPlot({
    base <- pars()
    delays <- c(0, 60, 120, 180, 240, 300)
    res <- do.call(rbind, lapply(delays, function(dl) {
      o <- simulate(list(), modifyList(base, list(ONCEFF = 0.045, ONCSTART = dl)),
                    NULL, 540)
      b <- o[nrow(o), ]
      data.frame(start = dl, muscle = b$MUSC, satc = b$SATC, grip = b$GRIP,
                 pct = 100 * b$MUSC / o$MUSC[1])
    }))
    ggplot(res, aes(start, pct)) +
      geom_line(linewidth = 1, colour = PAL[2]) +
      geom_point(size = 3, colour = PAL[2]) +
      geom_text(aes(label = sprintf("SATC %.2f", satc)), vjust = -1.1, size = 3.3,
                colour = "#666") +
      labs(x = "Curative anticancer therapy start day",
           y = "Muscle mass at day 540 (% of pre-disease)",
           title = "Same therapy, different timing — once the satellite-cell pool is exhausted, it cannot be reversed") +
      expand_limits(y = c(60, 105)) + THEME
  })
}

shinyApp(ui, server)
