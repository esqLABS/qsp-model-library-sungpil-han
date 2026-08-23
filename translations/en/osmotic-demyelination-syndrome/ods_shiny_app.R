## =====================================================================
##  ods_shiny_app.R
##  Osmotic Demyelination Syndrome (ODS) — interactive QSP dashboard
##  Osmotic Demyelination Syndrome — Interactive Dashboard (11 tabs)
##
##  The app is organised around one question that the model answers and a
##  bedside chart cannot:  WHO IS SETTING THE CORRECTION RATE?
##
##  Tab 1  Patient Profile          Patient & phenotype
##  Tab 2  Electrolytes · Fluid PK         Sodium, water, urine
##  Tab 3  The Kidney — What Sets the Rate The kidney: AVP, U_osm, EFWC
##  Tab 4  Brain Osmolytes (PD Core)     Brain osmolytes and OMEGA
##  Tab 5  Injury Cascade         Astrocyte -> microglia -> myelin
##  Tab 6  Clinical Endpoints            Deficit, MRI, the biphasic course
##  Tab 7  Central Experiment (AVP Switch)  The counterfactual
##  Tab 8  Correction-Rate Map         The FOSM x rate safety map
##  Tab 9  Relowering Deadline         The relowering deadline
##  Tab 10 Scenario Comparison          Scenario comparison
##  Tab 11 Biomarkers & Model Explanation  Biomarkers and what is fitted
##
##  Requires: shiny, mrgsolve, dplyr, tidyr, ggplot2, DT
##  Run:  shiny::runApp("ods_shiny_app.R")
## =====================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

source_model <- function() {
  e <- new.env()
  sys.source("ods_mrgsolve_model.R", envir = e)
  e
}
E   <- source_model()
mod <- E$mod

THEME <- theme_bw(base_size = 12) +
  theme(strip.background = element_rect(fill = "grey93"),
        legend.position = "bottom", panel.grid.minor = element_blank())

## ---------------------------------------------------------------------
##  Phenotype presets.  Each one is a different KIDNEY at the same sodium.
## ---------------------------------------------------------------------
PHENO <- list(
  "SIADH (euvolaemic)"          = list(kind = "siadh",     days = 21,
    note = "AVP is secreted autonomously. Urine osmolality is fixed high, and the urine (Na+K)/serum Na ratio exceeds 1, so fluid restriction alone will not correct it."),
  "Hypovolaemic (vomiting/diuretics)"      = list(kind = "hypovol",   days = 7,
    note = "★ The phenotype in which ODS most readily occurs. AVP is elevated by the volume stimulus, and once normal saline replenishes volume, that stimulus disappears and the kidney raises sodium on its own."),
  "Thiazide-induced"              = list(kind = "thiazide",  days = 21,
    note = "The diluting segment (NCC) is blocked, so minimum urine osmolality is elevated. Diluting capacity returns the moment the drug is stopped."),
  "Beer potomania / tea-and-toast" = list(kind = "potomania", days = 21,
    note = "Solute intake is <250 mOsm/d, so free-water excretion capacity itself is low. Restoring the diet alone can cause overcorrection."),
  "Adrenal insufficiency"                   = list(kind = "adrenal",   days = 21,
    note = "Cortisol deficiency fails to suppress AVP. Water diuresis opens up the moment hydrocortisone is given."),
  "Acute (<12 h, exercise/MDMA)"     = list(kind = "acute",     days = 8/24,
    note = "★ The opposite-direction disease. Organic osmolytes have not yet been extruded, so the brain is swollen. The danger is herniation, not demyelination.")
)

build_patient <- function(kind, days, target = 110) {
  if (kind == "acute") {
    w <- uniroot(function(x) {
      o <- mod %>% param(AVPSIADH = 6, ACUTE = 1, WLOAD = x, TWLEND = days) %>%
             mrgsim(end = days, delta = 1/96)
      tail(o$SODIUM, 1) - target
    }, c(1, 60), tol = 1e-4)$root
    p <- list(AVPSIADH = 6, ACUTE = 0)
    o <- mod %>% param(AVPSIADH = 6, ACUTE = 1, WLOAD = w, TWLEND = days) %>%
           mrgsim(end = days, delta = 1/96)
    return(list(param = p, state = as.numeric(tail(o, 1)[, names(init(mod))]), build = o))
  }
  w <- E$solve_win(mod, kind, target, days)
  r <- E$make_chronic(mod, kind, win = w, days = days)
  list(param = r$param, state = r$state, build = r$out)
}

start_at <- function(m, st)
  do.call(init, c(list(m), as.list(setNames(st, names(init(mod))))))

## =====================================================================
##  UI
## =====================================================================
ui <- fluidPage(
  titlePanel("Osmotic Demyelination Syndrome (ODS) — QSP Dashboard"),
  tags$p(style = "color:#555;margin-top:-8px",
    HTML("<b>Ω = ORG_set(Osm_eff) − ORG</b> — the organic-osmolyte deficit. ",
         "Zero in the normal brain, zero even in an <i>adapted</i> hyponatraemic brain, ",
         "and becomes positive only when plasma tonicity moves faster than transcription.")),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectInput("pheno", "Phenotype", names(PHENO)),
      sliderInput("na0", "Serum [Na] at presentation (mmol/L)", 100, 128, 110, step = 1),
      hr(),
      h5("Correction prescription"),
      radioButtons("mode", NULL,
        c("Rate-targeted 3% NaCl titration" = "ctrl",
          "0.9% normal saline only"       = "ns",
          "Fluid restriction only"             = "fr",
          "Tolvaptan 15 mg/day"         = "tlv",
          "Oral urea 30 g/day"       = "urea"), selected = "ctrl"),
      conditionalPanel("input.mode == 'ctrl'",
        sliderInput("rate", "Target correction rate (mmol/L/24 h)", 2, 24, 6, step = 1),
        sliderInput("cap",  "Correction endpoint [Na]", 120, 145, 135, step = 1)),
      conditionalPanel("input.mode == 'ns'",
        sliderInput("r09", "0.9% NaCl (L/day)", 0, 4, 2, step = 0.5)),
      conditionalPanel("input.mode == 'fr'",
        sliderInput("win", "Fluid intake (L/day)", 0.3, 2.5, 1.0, step = 0.1)),
      hr(),
      h5("Rate control"),
      checkboxInput("ddavp", "Pre-emptive desmopressin clamp 2 µg q8h", FALSE),
      checkboxInput("rescue", "Relowering on overcorrection (D5W + DDAVP)", FALSE),
      conditionalPanel("input.rescue",
        sliderInput("trescue", "Relowering start (hours)", 6, 96, 24, step = 6),
        sliderInput("nares",   "Relowering target [Na]",   110, 128, 118, step = 1)),
      hr(),
      h5("Risk modifiers"),
      sliderInput("fosm", "Organic-osmolyte transport capacity FOSM", 0.35, 1.0, 1.0, step = 0.05),
      helpText(HTML("<small>1.00 normal · 0.55 alcoholism/malnutrition · 0.70 liver disease</small>")),
      sliderInput("fnut", "Astrocyte energy reserve FNUT", 0.6, 1.0, 1.0, step = 0.05),
      sliderInput("kdef", "Total body potassium deficit (mmol)", 0, 700, 0, step = 50),
      sliderInput("rkcl", "KCl replacement (mmol/day)", 0, 200, 0, step = 20),
      hr(),
      h5("Adjunct treatment (preclinical evidence)"),
      checkboxInput("dex",  "Dexamethasone 16 mg/day", FALSE),
      checkboxInput("mino", "Minocycline 200 mg/day", FALSE),
      hr(),
      sliderInput("tend", "Simulation duration (days)", 7, 90, 45, step = 1),
      actionButton("go", "Run", class = "btn-primary btn-block")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",

        tabPanel("① Patient Profile",
          br(), htmlOutput("phenoNote"),
          h4("The Adapted Brain the Model Builds for Itself"),
          p(HTML("The values below are <b>an output, not an input</b>. When hyponatraemia is generated by simulation, the brain sheds osmolytes on its own.")),
          DTOutput("tblPatient"), br(),
          plotOutput("plotBuild", height = "340px")),

        tabPanel("② Electrolytes · Fluid",
          br(), plotOutput("plotNa", height = "300px"),
          plotOutput("plotFluid", height = "300px"),
          h5("Maximum 24-hour rise"), verbatimTextOutput("txtRate")),

        tabPanel("③ The Kidney — What Sets the Rate",
          br(),
          p(HTML("<b>EFWC = V<sub>urine</sub> × (1 − (U<sub>Na</sub>+U<sub>K</sub>)/S<sub>Na</sub>)</b>. ",
                 "It is this value, not the prescription, that sets the correction rate.")),
          plotOutput("plotKidney", height = "420px"),
          h5("Furst ratio (U_Na+U_K)/S_Na — above 1, fluid restriction fails"),
          plotOutput("plotFurst", height = "220px")),

        tabPanel("④ Brain Osmolytes (Ω)",
          br(),
          p(HTML("It leaves through channels and returns <b>via transcription</b>. That asymmetry is the disease.")),
          plotOutput("plotOsmo", height = "360px"),
          plotOutput("plotOmega", height = "300px")),

        tabPanel("⑤ Injury Cascade",
          br(),
          p(HTML("The order matters: <b>astrocyte → microglia·BBB → oligodendrocyte → myelin</b>.")),
          plotOutput("plotInjury", height = "480px")),

        tabPanel("⑥ Clinical Endpoints",
          br(),
          p(HTML("<b>Biphasic course</b>: sodium returns to normal and the patient improves, and only then does deterioration arrive. ",
                 "and MRI lags even further behind.")),
          plotOutput("plotClin", height = "380px"),
          DTOutput("tblTiming")),

        tabPanel("⑦ Central Experiment (AVP Switch)",
          br(),
          p(HTML("Same patient, same fluids — hypertonic saline is <b>0 mL in both arms</b>. ",
                 "Only one arm has AVP responding physiologically.")),
          actionButton("goCF", "Run counterfactual experiment", class = "btn-warning"),
          br(), br(), plotOutput("plotCF", height = "480px"),
          verbatimTextOutput("txtCF")),

        tabPanel("⑧ Correction-Rate Map",
          br(),
          p(HTML("The threshold Ω* = 8 mOsm/kg is <b>a single number, the same for every patient</b>. ",
                 "Risk factors act only through FOSM (transporter capacity).")),
          actionButton("goMap", "Calculate safety map", class = "btn-warning"),
          br(), br(), plotOutput("plotMap", height = "440px"),
          DTOutput("tblMap")),

        tabPanel("⑨ Relowering Deadline",
          br(),
          p(HTML("After overcorrection has occurred, how late can relowering (D5W + DDAVP + fluid restriction) be started?")),
          actionButton("goResc", "Calculate deadline", class = "btn-warning"),
          br(), br(), plotOutput("plotResc", height = "400px"),
          DTOutput("tblResc")),

        tabPanel("⑩ Scenario Comparison",
          br(),
          checkboxGroupInput("cmp", "Prescriptions to compare", inline = TRUE,
            choices = c("3% NaCl +6/day", "3% NaCl +12/day", "3% NaCl +20/day",
                        "0.9% NaCl 2 L/day", "Fluid restriction 1 L/day",
                        "Tolvaptan 15 mg", "DDAVP clamp + 3% NaCl"),
            selected = c("3% NaCl +6/day", "3% NaCl +12/day", "0.9% NaCl 2 L/day")),
          actionButton("goCmp", "Run comparison", class = "btn-warning"),
          br(), br(), plotOutput("plotCmp", height = "460px"),
          DTOutput("tblCmp")),

        tabPanel("⑪ Biomarkers & Model Explanation",
          br(),
          h4("What Can Be Measured and What Cannot"),
          plotOutput("plotBio", height = "340px"),
          hr(),
          h4("What Was Fitted and What Was Predicted"),
          htmlOutput("txtFit"))
      )
    )
  )
)

## =====================================================================
##  SERVER
## =====================================================================
server <- function(input, output, session) {

  patient <- eventReactive(input$go, {
    ph <- PHENO[[input$pheno]]
    withProgress(message = "Simulating the adapted brain...", {
      build_patient(ph$kind, ph$days, input$na0)
    })
  }, ignoreNULL = FALSE)

  rx_param <- reactive({
    p <- list(FOSM = input$fosm, FNUT = input$fnut, RKCL = input$rkcl,
              DEXON = as.numeric(input$dex), MINOON = as.numeric(input$mino))
    switch(input$mode,
      ctrl = c(p, list(CTRLON = 1, RATETGT = input$rate, NASTART = input$na0,
                       NACAP = input$cap, WIN = 1.0)),
      ns   = c(p, list(R09 = input$r09)),
      fr   = c(p, list(WIN = input$win)),
      tlv  = p,
      urea = c(p, list(UREADOSE = 30, WIN = 1.0)))
  })

  sim <- eventReactive(input$go, {
    pt <- patient()
    st <- pt$state
    st[which(names(init(mod)) == "KE")] <- st[which(names(init(mod)) == "KE")] - input$kdef
    pp <- c(pt$param, rx_param())
    if (PHENO[[input$pheno]]$kind == "hypovol")
      pp <- c(pp, list(NALOSS = 0, KLOSS = 0, WLOSS = 0, TLOSSEND = 0))
    ev <- NULL
    if (input$ddavp) ev <- E$ddavp_q8(0, 6)
    if (input$mode == "tlv") ev <- c(ev, E$tolvaptan(7))
    m <- start_at(mod, st) %>% param(pp)
    if (!input$rescue) return(as_tibble(m %>% mrgsim(events = ev, end = input$tend, delta = 0.25)))
    ## the rescue is a change of PRESCRIPTION, so it is run as a second phase
    th <- input$trescue / 24
    a  <- m %>% mrgsim(events = ev, end = th, delta = 0.05)
    st2 <- as.numeric(tail(a, 1)[, names(init(mod))])
    pr <- c(pp, list(RESCUE = 1, TRESCUE = th, NARES = input$nares,
                     DURRES = 1.0, WIN = 0.5, R09 = 0))
    b <- start_at(mod, st2) %>% param(pr) %>%
      mrgsim(events = E$ddavp_q8(th, th + 2.5), end = input$tend, delta = 0.25)
    bind_rows(as_tibble(a), as_tibble(b) %>% filter(time > 0) %>% mutate(time = time))
  }, ignoreNULL = FALSE)

  gg <- function(d, vars, labs_) {
    d %>% select(time, all_of(vars)) %>% pivot_longer(-time) %>%
      mutate(name = factor(name, levels = vars, labels = labs_)) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.85, colour = "#2b6cb0") +
      facet_wrap(~name, scales = "free_y") + labs(x = "Day", y = NULL) + THEME
  }

  ## --- Tab 1 -----------------------------------------------------------
  output$phenoNote <- renderUI({
    HTML(paste0("<div style='background:#f0f6ff;padding:10px;border-left:4px solid #2b6cb0'>",
                PHENO[[input$pheno]]$note, "</div>"))
  })
  output$tblPatient <- renderDT({
    b <- as_tibble(patient()$build); e <- tail(b, 1)
    datatable(tibble(
      Item = c("Serum [Na] (mmol/L)", "Effective osmolality (mOsm/kg)", "Brain water (mL/100 g, normal 80.00)",
               "Total organic osmolytes (mOsm/kg, normal 48.0)", "myo-inositol (normal 7.0)",
               "Inorganic ions (normal 227)", "Ω (organic-osmolyte deficit)", "Plasma AVP (pg/mL)",
               "Urine osmolality (mOsm/kg)", "Urine output (L/day)", "Furst ratio (U_Na+U_K)/S_Na"),
      Value = round(c(e$SODIUM, e$TONICITY, e$BRAINH2O, e$ORGTOT, e$MYOINOS,
                   e$ORGTOT * 0 + (e$TONICITY - 10 - e$ORGTOT), e$OMEGAc,
                   e$AVPpg, e$UOSMOL, e$UVOL, e$FURST), 2)),
      rownames = FALSE, options = list(dom = "t", pageLength = 20))
  })
  output$plotBuild <- renderPlot({
    gg(as_tibble(patient()$build),
       c("SODIUM","BRAINH2O","ORGTOT","MYOINOS"),
       c("Serum [Na]","Brain water (mL/100 g)","Total organic osmolytes","myo-inositol")) +
      ggtitle("The Brain Sheds Osmolytes While Hyponatraemia Develops")
  })

  ## --- Tab 2 -----------------------------------------------------------
  output$plotNa <- renderPlot(
    gg(sim(), c("SODIUM","TONICITY","POTASSIUM","BUNmgdL"),
       c("Serum [Na] (mmol/L)","Effective osmolality","Serum [K]","BUN (mg/dL)")))
  output$plotFluid <- renderPlot(
    gg(sim(), c("UVOL","UOSMOL","EFWCL","BRAINH2O"),
       c("Urine output (L/day)","Urine osmolality","Free-water clearance EFWC (L/day)","Brain water")))
  output$txtRate <- renderText({
    d <- sim(); r <- E$na_rise_24(d$time, d$SODIUM)
    lim <- if (input$fosm < 0.8) 8 else 12
    sprintf("Maximum rise in any 24-hour window = %.1f mmol/L   (this patient's model limit ≈ %d)\n%s",
            r, lim, if (r > lim) "*** limit exceeded ***" else "within limit")
  })

  ## --- Tab 3 -----------------------------------------------------------
  output$plotKidney <- renderPlot(
    gg(sim(), c("AVPpg","UOSMOL","UVOL","EFWCL","UNAKc","SODIUM"),
       c("Plasma AVP (pg/mL)","Urine osmolality (mOsm/kg)","Urine output (L/day)",
         "EFWC (L/day)","Urine (Na+K) (mmol/L)","Serum [Na]")))
  output$plotFurst <- renderPlot({
    ggplot(sim(), aes(time, FURST)) + geom_line(linewidth = 0.9, colour = "#b0416a") +
      geom_hline(yintercept = 1, linetype = 2) +
      labs(x = "Day", y = "(U_Na+U_K)/S_Na") + THEME
  })

  ## --- Tab 4 -----------------------------------------------------------
  output$plotOsmo <- renderPlot({
    d <- as_tibble(sim())
    d %>% select(time, ORGTOT, ORGSETC, MYOINOS) %>%
      pivot_longer(-time) %>%
      mutate(name = recode(name, ORGTOT = "Actual organic osmolytes (ORG)",
                           ORGSETC = "Value current tonicity demands (ORG_set)",
                           MYOINOS = "myo-inositol")) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.95) +
      labs(x = "Day", y = "mOsm/kg brain water", colour = NULL,
           title = "The Disease Is Failing to Keep Up") + THEME
  })
  output$plotOmega <- renderPlot({
    d <- as_tibble(sim())
    ggplot(d, aes(time)) +
      geom_ribbon(aes(ymin = THRESH, ymax = pmax(STRESSc, THRESH)),
                  fill = "#f4b8b8", alpha = 0.7) +
      geom_line(aes(y = STRESSc), linewidth = 1, colour = "#b0416a") +
      geom_line(aes(y = OMEGAc), linewidth = 0.8, colour = "#2f7a53") +
      geom_hline(aes(yintercept = THRESH), linetype = 2) +
      labs(x = "Day", y = "mOsm/kg",
           title = "Green = Ω (organic-osmolyte deficit) · Red = total osmotic stress · Dashed = threshold 8.0") + THEME
  })

  ## --- Tab 5 -----------------------------------------------------------
  output$plotInjury <- renderPlot(
    gg(sim(), c("STRESSc","INJRATE","BRAINH2O","SHRINKc","MRIPOS","DEFICIT"),
       c("Osmotic stress","Injury rate","Brain water","Brain shrinkage (%)","MRI signal","Clinical deficit")))

  ## --- Tab 6 -----------------------------------------------------------
  output$plotClin <- renderPlot({
    d <- as_tibble(sim())
    d %>% select(time, SODIUM, DEFICIT, MRIPOS) %>% pivot_longer(-time) %>%
      mutate(name = recode(name, SODIUM = "Serum [Na]",
                           DEFICIT = "Clinical deficit (0-100)", MRIPOS = "MRI signal")) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.95, colour = "#2b6cb0") +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      labs(x = "Day", y = NULL,
           title = "Sodium Normalises First, and the Patient Worsens Afterward") + THEME
  })
  output$tblTiming <- renderDT({
    d <- as_tibble(sim())
    f <- function(v, th) { i <- which(v > th)[1]; if (is.na(i)) NA else d$time[i] }
    datatable(tibble(
      Event = c("Serum [Na] normalised (>135)", "Symptom onset (deficit >10)",
               "MRI positive (signal >0.25)", "Maximum deficit"),
      `Day of occurrence` = c(f(d$SODIUM, 135), f(d$DEFICIT, 10), f(d$MRIPOS, 0.25),
                   d$time[which.max(d$DEFICIT)])),
      rownames = FALSE, options = list(dom = "t"))
  })

  ## --- Tab 7 : the counterfactual ---------------------------------------
  cf <- eventReactive(input$goCF, {
    pt <- build_patient("hypovol", 7, input$na0)
    base <- c(pt$param, list(NALOSS = 0, KLOSS = 0, WLOSS = 0, TLOSSEND = 0, R09 = 2))
    iavp <- which(names(init(mod)) == "AVP")
    A <- start_at(mod, pt$state) %>% param(base) %>% mrgsim(end = 30, delta = 0.1)
    B <- start_at(mod, pt$state) %>%
      param(c(base, list(AVPFREEZE = 1, TFREEZE = 0, AVPFRZV = pt$state[iavp]))) %>%
      mrgsim(end = 30, delta = 0.1)
    bind_rows(as_tibble(A) %>% mutate(arm = "AVP responds (physiological)"),
              as_tibble(B) %>% mutate(arm = "AVP fixed (counterfactual)"))
  })
  output$plotCF <- renderPlot({
    cf() %>% filter(time <= 10) %>%
      select(time, arm, SODIUM, AVPpg, UOSMOL, EFWCL, STRESSc, DEFICIT) %>%
      pivot_longer(-c(time, arm)) %>%
      mutate(name = recode(name, SODIUM = "Serum [Na]", AVPpg = "Plasma AVP",
                           UOSMOL = "Urine osmolality", EFWCL = "EFWC (L/day)",
                           STRESSc = "Osmotic stress", DEFICIT = "Clinical deficit")) %>%
      ggplot(aes(time, value, colour = arm)) + geom_line(linewidth = 1) +
      facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = c("#b0416a", "#2b6cb0")) +
      labs(x = "Day", y = NULL, colour = NULL,
           title = "Hypertonic Saline Is 0 mL in Both Arms. The Kidney Makes All the Difference.") + THEME
  })
  output$txtCF <- renderText({
    d <- cf()
    s <- d %>% group_by(arm) %>%
      summarise(rise24 = E$na_rise_24(time, SODIUM), str = max(STRESSc),
                def = max(DEFICIT), mye = min(MYE), .groups = "drop")
    paste(apply(s, 1, function(r)
      sprintf("%-24s 24h rise %5.1f   max stress %5.1f   deficit %5.1f   myelin nadir %5.3f",
              r[1], as.numeric(r[2]), as.numeric(r[3]), as.numeric(r[4]), as.numeric(r[5]))),
      collapse = "\n")
  })

  ## --- Tab 8 : the safety map -------------------------------------------
  smap <- eventReactive(input$goMap, {
    pt <- build_patient("siadh", 21, input$na0)
    grid <- expand.grid(fosm = c(1.00, 0.80, 0.65, 0.55, 0.45),
                        rate = c(4, 6, 8, 10, 12, 14, 16, 20))
    withProgress(message = "Calculating safety map...", {
      grid$stress <- mapply(function(f, r) {
        o <- start_at(mod, pt$state) %>%
          param(c(pt$param, list(CTRLON = 1, RATETGT = r, NASTART = input$na0,
                                 NACAP = 140, WIN = 1, FOSM = f))) %>%
          mrgsim(end = 6, delta = 0.1)
        max(o$STRESSc)
      }, grid$fosm, grid$rate)
    })
    grid
  })
  output$plotMap <- renderPlot({
    ggplot(smap(), aes(factor(rate), factor(fosm), fill = stress)) +
      geom_tile(colour = "white") +
      geom_text(aes(label = sprintf("%.1f", stress)), size = 3.4) +
      scale_fill_gradient2(low = "#cfe9da", mid = "#f7e7c2", high = "#e79a9a",
                           midpoint = 8, name = "Maximum stress") +
      labs(x = "Target correction rate (mmol/L/24 h)",
           y = "FOSM (organic-osmolyte transport capacity)",
           title = "Cells Above the Threshold of 8.0 Are Dangerous — the Guideline's Two Numbers Come From Here") + THEME
  })
  output$tblMap <- renderDT({
    smap() %>% mutate(Safety = ifelse(stress < 8, "Safe", "Threshold exceeded")) %>%
      datatable(rownames = FALSE, options = list(pageLength = 10))
  })

  ## --- Tab 9 : the relowering deadline ------------------------------------
  resc <- eventReactive(input$goResc, {
    pt <- build_patient("hypovol", 7, input$na0)
    pA <- c(pt$param, list(NALOSS = 0, KLOSS = 0, WLOSS = 0, TLOSSEND = 0, R09 = 2))
    Ts <- c(8, 12, 18, 24, 36, 48, 72)
    withProgress(message = "Calculating deadline...", {
      out <- lapply(Ts, function(T) {
        th <- T / 24
        a <- start_at(mod, pt$state) %>% param(pA) %>% mrgsim(end = th, delta = 0.05)
        st2 <- as.numeric(tail(a, 1)[, names(init(mod))])
        b <- start_at(mod, st2) %>%
          param(c(pA, list(RESCUE = 1, TRESCUE = th, NARES = 118,
                           DURRES = 1, WIN = 0.5, R09 = 0))) %>%
          mrgsim(events = E$ddavp_q8(th, th + 2.5), end = 90, delta = 0.25)
        tibble(T = T, deficit = max(b$DEFICIT), astro = min(b$AST))
      })
      none <- start_at(mod, pt$state) %>% param(pA) %>% mrgsim(end = 90, delta = 0.25)
      bind_rows(bind_rows(out),
                tibble(T = NA, deficit = max(none$DEFICIT), astro = min(none$AST)))
    })
  })
  output$plotResc <- renderPlot({
    d <- resc(); n <- d$deficit[is.na(d$T)]
    ggplot(filter(d, !is.na(T)), aes(T, deficit)) +
      geom_hline(yintercept = n, linetype = 2, colour = "#b0416a") +
      geom_line(linewidth = 1, colour = "#2b6cb0") + geom_point(size = 3) +
      annotate("text", x = 60, y = n, vjust = -0.6, colour = "#b0416a",
               label = "No relowering") +
      labs(x = "Relowering start time (hours)", y = "Maximum clinical deficit",
           title = "The Deadline Is Set by the Time Constant of Astrocyte Death") + THEME
  })
  output$tblResc <- renderDT(datatable(resc(), rownames = FALSE,
                                       options = list(dom = "t")))

  ## --- Tab 10 : scenario comparison -----------------------------------------
  cmp <- eventReactive(input$goCmp, {
    pt <- patient(); base <- pt$param
    if (PHENO[[input$pheno]]$kind == "hypovol")
      base <- c(base, list(NALOSS = 0, KLOSS = 0, WLOSS = 0, TLOSSEND = 0))
    spec <- list(
      "3% NaCl +6/day"        = list(p = list(CTRLON=1, RATETGT=6,  NASTART=input$na0, NACAP=135, WIN=1), e = NULL),
      "3% NaCl +12/day"       = list(p = list(CTRLON=1, RATETGT=12, NASTART=input$na0, NACAP=135, WIN=1), e = NULL),
      "3% NaCl +20/day"       = list(p = list(CTRLON=1, RATETGT=20, NASTART=input$na0, NACAP=135, WIN=1), e = NULL),
      "0.9% NaCl 2 L/day"     = list(p = list(R09 = 2), e = NULL),
      "Fluid restriction 1 L/day"       = list(p = list(WIN = 1.0), e = NULL),
      "Tolvaptan 15 mg"          = list(p = list(), e = E$tolvaptan(7)),
      "DDAVP clamp + 3% NaCl"= list(p = list(CTRLON=1, RATETGT=6, NASTART=input$na0, NACAP=135, WIN=1),
                                     e = E$ddavp_q8(0, 6)))
    withProgress(message = "Comparing scenarios...", {
      bind_rows(lapply(input$cmp, function(k) {
        s <- spec[[k]]
        as_tibble(start_at(mod, pt$state) %>%
          param(c(base, s$p, list(FOSM = input$fosm, FNUT = input$fnut))) %>%
          mrgsim(events = s$e, end = input$tend, delta = 0.25)) %>% mutate(arm = k)
      }))
    })
  })
  output$plotCmp <- renderPlot({
    cmp() %>% select(time, arm, SODIUM, EFWCL, OMEGAc, STRESSc, MYE, DEFICIT) %>%
      pivot_longer(-c(time, arm)) %>%
      mutate(name = recode(name, SODIUM = "Serum [Na]", EFWCL = "EFWC (L/day)",
                           OMEGAc = "Ω", STRESSc = "Osmotic stress",
                           MYE = "Myelin content", DEFICIT = "Clinical deficit")) %>%
      ggplot(aes(time, value, colour = arm)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "Day", y = NULL, colour = NULL) + THEME
  })
  output$tblCmp <- renderDT({
    cmp() %>% group_by(arm) %>%
      summarise(`24h rise` = round(E$na_rise_24(time, SODIUM), 1),
                `[Na] 48h` = round(approx(time, SODIUM, 2)$y, 1),
                `Maximum Ω` = round(max(OMEGAc), 1),
                `Maximum stress` = round(max(STRESSc), 1),
                `Myelin nadir` = round(min(MYE), 3),
                `Maximum deficit` = round(max(DEFICIT), 1), .groups = "drop") %>%
      datatable(rownames = FALSE, options = list(dom = "t"))
  })

  ## --- Tab 11 -----------------------------------------------------------------
  output$plotBio <- renderPlot(
    gg(sim(), c("MYOINOS","BRAINH2O","MRIPOS","DEFICIT"),
       c("¹H-MRS myo-inositol (measurable)", "Brain water (hard to measure)",
         "MRI signal (delayed)", "Clinical deficit")))
  output$txtFit <- renderUI(HTML("
<table class='table table-sm'>
<tr><th>Fitted to normal physiology</th><td>Water/solute steady state (WIN was left free so the healthy model reaches an exact steady state: max|dy/dt| = 2.5e-11), the Edelman regression, urine-concentrating range, plasma urea, AVP osmotic threshold and gain</td></tr>
<tr><th>Fitted to adaptation data</th><td>Osmotic-response coefficients β<sub>i</sub> (at [Na] 110 adaptation, total organic osmolytes −40%, myo-inositol −66%), efflux/influx time constants (myo-inositol recovery ~5–6 days)</td></tr>
<tr><th>Fitted to the injury dose-response</th><td>KINJ · KAST · KOLI · KDEM, four parameters. The threshold Ω* = 8 was calibrated once so that <b>the normal-risk limit falls at 10–12</b>.</td></tr>
<tr><th style='color:#2f7a53'>Not fitted (= predicted)</th><td>High-risk limit of 6–8 (a consequence of FOSM 0.55), acute/chronic asymmetry, the kidney's autonomous overcorrection, the size of the DDAVP clamp effect, the relowering deadline, the Furst-ratio rule, the adverse effect of normal saline when urine osmolality &gt; 308, MRI delay</td></tr>
<tr><th style='color:#b0416a'>The most exposed assumption</th><td>Urea does not reduce the osmotic injury itself (Ω) and acts only through the BBB/microglial arm. If, in a rate-matched experiment, urea reduces an early marker of astrocyte death, this model is wrong.</td></tr>
<tr><th>Known limitations</th><td>R is not available in the build container, so this model has been <b>verified by equation but not by compilation</b>. Every ODE was independently reimplemented and integrated in Python/scipy (ods_verify_python.py). Ω has never been measured in a living human brain.</td></tr>
</table>"))
}

shinyApp(ui, server)
