## =============================================================================
##  nec_shiny_app.R
##  Necrotising Enterocolitis (NEC) — interactive QSP dashboard
##  Neonatal necrotising enterocolitis QSP dashboard
##
##  Requires nec_mrgsolve_model.R in the same directory.
##      shiny::runApp("nec_shiny_app.R")
##
##  DESIGN INTENT
##  -------------
##  The app is organised around the ONE STRUCTURAL CLAIM of the model, not
##  around the list of state variables.  Every tab answers one question:
##
##    ① Patient     who is this infant, and how far is their own threshold?
##    ② Core loop   is the positive loop closing?  (E -> Pb -> Jtr -> INJ -> E)
##    ③ Gut ecology which of the two ecological basins did the gut land in?
##    ④ Crit. load  where are B_lo and B_hi for THIS infant, and where is B now?
##    ⑤ Lesion      has the lesion passed NECa_crit = INJth/wNEC?
##    ⑥ Endpoints   Bell stage, surgery, growth, cholestasis, NDI index
##    ⑦ Scenarios   arm-vs-arm comparison under identical provocation
##    ⑧ Biomarkers  CRP / platelets / lactate / I-FABP-like readouts and their
##                 lead time relative to Bell II — the point being that they
##                 mostly do NOT lead
##    ⑨ PK          plasma AND luminal exposure, because only the luminal
##                 concentration reshapes the microbiome
##
##  DISCLAIMER: teaching / research tool.  Not for clinical use.
## =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)

source("nec_mrgsolve_model.R", local = TRUE)

PAL <- c("#2f6fb0", "#c0392b", "#2e8b57", "#b8860b", "#6a51a3", "#0e7c86")

thin <- function(d, every = 4) d[seq(1, nrow(d), by = every), , drop = FALSE]

long_plot <- function(d, vars, labs_map, title, ylab = NULL, free = TRUE) {
  d %>%
    select(time, all_of(vars)) %>%
    pivot_longer(-time) %>%
    mutate(name = factor(labs_map[name], levels = unname(labs_map[vars]))) %>%
    ggplot(aes(time, value)) +
    geom_line(colour = PAL[1], linewidth = 0.65) +
    facet_wrap(~name, scales = if (free) "free_y" else "fixed") +
    labs(x = "Postnatal day", y = ylab, title = title) +
    theme_bw(base_size = 11) +
    theme(strip.background = element_rect(fill = "#eef2f7"),
          panel.grid.minor = element_blank())
}

## =============================================================================
##  UI
## =============================================================================
ui <- fluidPage(
  tags$head(tags$style(HTML("
    .well { background:#f7f9fc; border-color:#dbe3ee; }
    h4 { color:#1f456f; margin-top:0; }
    .keynum { font-size:22px; font-weight:700; color:#1f456f; }
    .keylab { font-size:11px; color:#667; }
    .claim { background:#fff6f5; border-left:4px solid #c0392b;
             padding:8px 12px; margin-bottom:10px; font-size:12px; }
  "))),

  titlePanel("Neonatal Necrotising Enterocolitis (NEC) — QSP Simulator"),

  div(class = "claim",
      HTML("<b>The structural claim of this model:</b> NEC is a closed positive feedback loop hanging on <i>E</i>,
      enterocyte integrity, alone, and the gain of that loop is set by the intraluminal pathogen load
      <i>B</i>. &nbsp;
      <code>E → Pb → Jtr → TLR4s → (cell loss↑, proliferation↓) → E</code> &nbsp;
      At <b>B_lo</b> NEC becomes <i>possible</i> and at <b>B_hi</b> it becomes <i>inevitable</i>;
      between the two it is not the load but the <b>precipitating insult</b> that decides fate.")),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("① Patient"),
      sliderInput("GA", "Gestational age GA (weeks)", 23, 34, 27, step = 0.5),
      sliderInput("BW", "Birth weight (kg)", 0.40, 2.20, 0.92, step = 0.02),
      radioButtons("binfant", "Carries an HMO-metabolising strain (B. infantis)",
                   c("Carrier (~35%)" = "1",
                     "Non-carrier (~65%)" = "0.05"),
                   selected = "1"),
      sliderInput("dysbio", "Dysbiotic pressure (environmental influx multiple)",
                  0.5, 3.0, 1.5, step = 0.1),

      h4("② Nutrition"),
      selectInput("milk", "Feed type",
                  c("Mother's own milk" = "MOM",
                    "Pasteurised donor milk" = "DONOR",
                    "Preterm formula" = "FORMULA",
                    "Mixed" = "MIXED"), selected = "MOM"),
      sliderInput("feedrate", "Advancement rate (mL/kg/d per day)", 0, 40, 20, step = 5),
      sliderInput("feedmax", "Target feed volume (mL/kg/d)", 60, 200, 160, step = 10),
      checkboxInput("npo", "Maintain complete nil by mouth (NPO + TPN)", FALSE),

      h4("③ Drugs"),
      sliderInput("abx", "Empirical ampicillin+gentamicin (days)", 0, 14, 3, step = 1),
      sliderInput("mtz", "Metronidazole (days)", 0, 14, 0, step = 1),
      sliderInput("ind", "Indomethacin for PDA (days)", 0, 3, 0, step = 1),
      sliderInput("ibu", "Ibuprofen for PDA (days)", 0, 3, 0, step = 1),
      sliderInput("prob", "Probiotic daily dose (1e9)", 0, 2, 0, step = 0.1),
      sliderInput("thr", "Threshold-raising agent (thr_boost)", 1.0, 2.5, 1.0, step = 0.1),

      h4("④ Precipitating insult"),
      sliderInput("insday", "Day of the deep hypotensive episode", 3, 35, 14, step = 1),
      sliderInput("insdep", "Depth of the episode (fall in mesenteric perfusion)", 0, 0.6, 0.40, step = 0.02),
      sliderInput("insdur", "Duration of the episode (days)", 0.1, 1.5, 0.6, step = 0.1),
      checkboxInput("tx", "Two transfusions (day 10 · day 20)", TRUE),

      h4("⑤ Response"),
      sliderInput("dxlag", "Diagnosis–treatment delay dx_lag (days)", 0, 3, 1.2, step = 0.1),
      checkboxInput("treat", "NPO + broad-spectrum antibiotics at Bell ≥ II", TRUE),
      hr(),
      sliderInput("tend", "Simulation duration (days)", 21, 70, 56, step = 7),
      helpText(HTML("<small>This app is for teaching and research. Do not use it
                     for clinical decisions.</small>"))
    ),

    mainPanel(
      width = 9,
      fluidRow(
        column(2, div(class = "keylab", "Peak Bell stage"),
                  div(class = "keynum", textOutput("kBell", inline = TRUE))),
        column(2, div(class = "keylab", "Day of diagnosis (Bell II)"),
                  div(class = "keynum", textOutput("kOnset", inline = TRUE))),
        column(2, div(class = "keylab", "B_lo (possible)"),
                  div(class = "keynum", textOutput("kBlo", inline = TRUE))),
        column(2, div(class = "keylab", "B_hi (inevitable)"),
                  div(class = "keynum", textOutput("kBhi", inline = TRUE))),
        column(2, div(class = "keylab", "Peak B"),
                  div(class = "keynum", textOutput("kBmax", inline = TRUE))),
        column(2, div(class = "keylab", "Point of no return crossed"),
                  div(class = "keynum", textOutput("kCrit", inline = TRUE)))
      ),
      hr(),

      tabsetPanel(
        id = "tabs", type = "tabs",

        tabPanel("① Patient profile",
          br(),
          fluidRow(column(6, plotOutput("pMat", height = 300)),
                   column(6, plotOutput("pFeed", height = 300))),
          h4("What maturation sets"),
          tableOutput("tMat"),
          helpText(HTML("As corrected age rises TLR4 over-expression falls, and PAF-AH·IL-10·
          defensins·tight junctions rise together. The values in this table set <b>how large a
          response this infant will mount to the same bacterial load</b>."))
        ),

        tabPanel("② The core feedback loop",
          br(),
          plotOutput("pLoop", height = 620),
          helpText(HTML("Read the loop once round, from top to bottom:
          <code>E</code> → permeability <code>Pb</code> → translocation flux <code>Jtr</code>
          → NF-κB → total injury <code>INJ</code> → back to <code>E</code>.
          The lesion grows only while the injury exceeds <code>INJth</code>."))
        ),

        tabPanel("③ Gut ecology",
          br(),
          plotOutput("pEco", height = 460),
          h4("Which basin did it fall into"),
          tableOutput("tEco"),
          helpText(HTML("HMOs are a <b>private substrate</b> that only the commensals can use.
          So where the strain that metabolises HMOs is absent (<i>non-carrier</i>), even a breast-fed infant
          can fall into an Enterobacteriaceae-dominant state."))
        ),

        tabPanel("④ Critical load · watershed",
          br(),
          fluidRow(column(7, plotOutput("pField", height = 420)),
                   column(5, plotOutput("pWindow", height = 420))),
          helpText(HTML("Left: the reduced one-dimensional field <code>g(E)</code>. The roots
          are the fates, and the unstable middle root is the <b>watershed E*</b>.
          Right: this infant's <code>B_lo</code>–<code>B_hi</code> window and the trajectory of the
          actual load. Inside the window it is the insult, not the load, that decides."))
        ),

        tabPanel("⑤ Lesion · stage",
          br(),
          plotOutput("pLesion", height = 480),
          helpText(HTML("Necrotic tissue is itself an inflammatory stimulus (DAMPs), so the lesion
          sustains the injury it created. A <b>point of no return
          NECa_crit = INJth / wNEC</b> therefore exists, and whether that line is crossed during the
          diagnosis–treatment delay <code>dx_lag</code> is what separates medical NEC from surgical
          NEC. Intramural gas (pneumatosis) appears only after the mucosa has already been breached
          (<code>Pb &gt; PbPNE</code>)."))
        ),

        tabPanel("⑥ Clinical endpoints",
          br(),
          fluidRow(column(6, plotOutput("pBell", height = 300)),
                   column(6, plotOutput("pOrgan", height = 300))),
          h4("Endpoint summary"),
          tableOutput("tEnd")
        ),

        tabPanel("⑦ Scenario comparison",
          br(),
          checkboxGroupInput("arms", "Arms to compare (identical precipitating insult)",
            choices = c("Mother's own milk" = "MOM",
                        "Donor milk" = "DONOR",
                        "Formula" = "FORMULA",
                        "Formula + probiotic" = "FORM_PROB",
                        "Formula + threshold-raising agent" = "FORM_THR",
                        "Formula + both" = "FORM_BOTH",
                        "Nil by mouth + TPN" = "NPO"),
            selected = c("MOM", "FORMULA", "FORM_PROB", "FORM_BOTH"),
            inline = TRUE),
          actionButton("go", "Run scenarios", class = "btn-primary"),
          br(), br(),
          plotOutput("pArms", height = 420),
          tableOutput("tArms"),
          helpText(HTML("<b>A drug that moves the threshold</b> and <b>a drug that moves the state</b>
          do not add, they multiply. Breast milk is both, so the fair comparator is not one
          single-target drug but <b>the two combined</b>."))
        ),

        tabPanel("⑧ Biomarkers · lead time",
          br(),
          plotOutput("pBio", height = 420),
          h4("Lead time relative to Bell II (negative = moves later)"),
          tableOutput("tLead"),
          helpText(HTML("What to see here is not that the markers <b>rise well</b>, but that
          most of them rise <b>only after the switch has flipped</b>. What moves earlier
          are the state variables (<code>E</code>, <code>Jtr</code>), and that is the structural reason
          why prevention comes before early diagnosis."))
        ),

        tabPanel("⑨ Pharmacokinetics (plasma vs lumen)",
          br(),
          plotOutput("pPK", height = 460),
          h4("Luminal exposure (the concentration the microbes actually see)"),
          tableOutput("tPK"),
          helpText(HTML("What changes the microbial community is not the plasma concentration but the
          <b>luminal concentration</b>. Ampicillin (f<sub>lum</sub>≈0.45) and
          metronidazole (≈0.80) reach the lumen, whereas an intravenous aminoglycoside
          (≈0.05) barely does. The asymmetry by which empirical ampicillin+gentamicin
          devastates only the anaerobes while leaving Enterobacteriaceae relatively intact
          comes from exactly this."))
        )
      )
    )
  )
)

## =============================================================================
##  SERVER
## =============================================================================
server <- function(input, output, session) {

  ## ---- insult schedule from the sidebar ----------------------------------
  insults <- reactive({
    z <- list(INS1A = input$insday, INS1B = input$insday + input$insdur,
              INS1DEP = input$insdep, INS1NFK = 0.010)
    if (isTRUE(input$tx)) {
      z <- c(z, list(INS2A = 10.0, INS2B = 10.5, INS2DEP = 0.16, INS2NFK = 0.030,
                     INS3A = 20.0, INS3B = 20.5, INS3DEP = 0.16, INS3NFK = 0.030))
    } else {
      z <- c(z, list(INS2A = 1e6, INS2B = 1e6, INS2DEP = 0, INS2NFK = 0,
                     INS3A = 1e6, INS3B = 1e6, INS3DEP = 0, INS3NFK = 0))
    }
    z
  })

  ## ---- the index simulation ---------------------------------------------
  sim <- reactive({
    run_nec(GA = input$GA, BW = input$BW, milk = input$milk,
            feedrate = input$feedrate, feedmax = input$feedmax,
            abx_days = input$abx, mtz_days = input$mtz,
            ind_days = input$ind, ibu_days = input$ibu,
            prob_dose = input$prob, thrboost = input$thr,
            dysbio = input$dysbio,
            binfant = as.numeric(input$binfant),
            NPO_ON = as.numeric(isTRUE(input$npo)),
            dx_lag = input$dxlag,
            npo_on_nec = as.numeric(isTRUE(input$treat)),
            insults = insults(), end = input$tend)
  })

  onset  <- reactive({ d <- sim(); if (any(d$BellStage >= 2))
                       min(d$time[d$BellStage >= 2]) else NA_real_ })
  critNa <- reactive({ first(sim()$NECa_crit) })

  loads <- reactive({
    pp  <- par_list()
    pp$Ktlr0 <- pp$Ktlr0 * input$thr
    ctx <- nec_context(pp, milk = input$milk, GA = input$GA,
                       feed = input$feedmax, probiotic = input$prob > 0)
    bifurcation_loads(pp, ctx)
  })

  ## ---- key numbers -------------------------------------------------------
  output$kBell  <- renderText(as.character(max(sim()$BellStage)))
  output$kOnset <- renderText(ifelse(is.na(onset()), "—",
                                     sprintf("d%.1f", onset())))
  output$kBlo   <- renderText(sprintf("%.1f", loads()[["B_lo"]]))
  output$kBhi   <- renderText(sprintf("%.1f", loads()[["B_hi"]]))
  output$kBmax  <- renderText(sprintf("%.1f", max(sim()$B)))
  output$kCrit  <- renderText(ifelse(any(sim()$NECA > critNa()), "Yes", "No"))

  ## ---- ① patient ---------------------------------------------------------
  output$pMat <- renderPlot({
    d <- thin(sim())
    pp <- par_list()
    d %>% mutate(TLR4expr = pp$TLR4max*(1 - pp$phiTLR*Maturation),
                 PAF_AH   = pp$dPAF0 + pp$dPAFm*Maturation) %>%
      select(time, Maturation, TLR4expr, PAF_AH) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
      geom_line(linewidth = 0.7) +
      scale_colour_manual(values = PAL[1:3], name = NULL) +
      labs(x = "Postnatal day", y = NULL, title = "Maturation and the coefficients it sets") +
      theme_bw(base_size = 11)
  })

  output$pFeed <- renderPlot({
    d <- thin(sim())
    d %>% select(time, FeedVol, SUB, Ischaemia) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) + geom_line(colour = PAL[4], linewidth = 0.7) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      labs(x = "Postnatal day", y = NULL,
           title = "Feeding cuts both ways: nutrition · substrate · oxygen demand") +
      theme_bw(base_size = 11)
  })

  output$tMat <- renderTable({
    d <- sim(); pp <- par_list()
    m0 <- first(d$Maturation); m1 <- last(d$Maturation)
    data.frame(
      Item = c("Maturation MAT", "TLR4 over-expression coefficient", "PAF-AH capacity",
               "Signal threshold Ktlr", "Trophic stimulus Ftroph (max)"),
      Start = c(sprintf("%.3f", m0),
               sprintf("%.3f", pp$TLR4max*(1 - pp$phiTLR*m0)),
               sprintf("%.2f", pp$dPAF0 + pp$dPAFm*m0),
               sprintf("%.1f", first(d$Threshold)), "—"),
      End = c(sprintf("%.3f", m1),
               sprintf("%.3f", pp$TLR4max*(1 - pp$phiTLR*m1)),
               sprintf("%.2f", pp$dPAF0 + pp$dPAFm*m1),
               sprintf("%.1f", max(d$Threshold)), "—"))
  })

  ## ---- ② the loop --------------------------------------------------------
  output$pLoop <- renderPlot({
    d <- thin(sim())
    lm <- c(E = "① Enterocyte integrity E", Permeability = "② Permeability Pb",
            Translocat = "③ Translocation flux Jtr", TLRS = "④ NF-κB activity",
            Injury = "⑤ Total injury INJ", NECA = "⑥ Full-thickness necrosis NECa")
    pp <- par_list()
    p <- d %>% select(time, E, Permeability, Translocat, TLRS, Injury, NECA) %>%
      pivot_longer(-time) %>%
      mutate(name = factor(lm[name], levels = unname(lm))) %>%
      ggplot(aes(time, value)) +
      geom_line(colour = PAL[2], linewidth = 0.7) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      labs(x = "Postnatal day", y = NULL,
           title = "Once round the closed positive feedback loop") +
      theme_bw(base_size = 11)
    if (!is.na(onset())) p <- p + geom_vline(xintercept = onset(),
                                            linetype = 2, colour = "#888")
    p
  })

  ## ---- ③ ecology --------------------------------------------------------
  output$pEco <- renderPlot({
    d <- thin(sim())
    lm <- c(B = "Pathogen B (Enterobacteriaceae)", C = "Commensal anaerobes C",
            P = "Probiotic P", SCFA = "SCFA (mM)",
            HMO = "HMO pool (g/kg)", GAS = "Luminal gas")
    long_plot(d, c("B", "C", "P", "SCFA", "HMO", "GAS"), lm,
              "Gut ecology — which of the two basins did it fall into")
  })

  output$tEco <- renderTable({
    d <- sim(); L <- loads()
    data.frame(
      Metric = c("Final B", "Final C", "Final SCFA", "B_lo", "B_hi",
               "Position of peak B within the window"),
      Value = c(sprintf("%.2f", last(d$B)), sprintf("%.2f", last(d$C)),
             sprintf("%.1f mM", last(d$SCFA)),
             sprintf("%.1f", L[["B_lo"]]), sprintf("%.1f", L[["B_hi"]]),
             sprintf("%.0f %%", 100*(max(d$B) - L[["B_lo"]]) /
                                (L[["B_hi"]] - L[["B_lo"]]))))
  })

  ## ---- ④ field and window ------------------------------------------------
  output$pField <- renderPlot({
    pp <- par_list(); pp$Ktlr0 <- pp$Ktlr0 * input$thr
    ctx <- nec_context(pp, milk = input$milk, GA = input$GA,
                       feed = input$feedmax, probiotic = input$prob > 0)
    L <- loads()
    Bs <- round(c(0.5*L[["B_lo"]], L[["B_lo"]],
                  0.5*(L[["B_lo"]] + L[["B_hi"]]), L[["B_hi"]]), 1)
    grid <- expand.grid(E = seq(0.02, 0.99, by = 0.005), B = Bs)
    grid$g <- mapply(reduced_field, grid$E, grid$B,
                     MoreArgs = list(pp = pp, ctx = ctx))
    ggplot(grid, aes(E, g, colour = factor(B))) +
      geom_hline(yintercept = 0, linetype = 2, colour = "#777") +
      geom_line(linewidth = 0.75) +
      scale_colour_manual(values = PAL, name = "B (1e9/g)") +
      labs(x = "Enterocyte integrity E", y = "g(E) = dE/dt",
           title = "The reduced one-dimensional field: the root is the fate") +
      theme_bw(base_size = 11)
  })

  output$pWindow <- renderPlot({
    d <- thin(sim()); L <- loads()
    ggplot(d, aes(time, B)) +
      annotate("rect", xmin = -Inf, xmax = Inf,
               ymin = L[["B_lo"]], ymax = L[["B_hi"]],
               fill = "#f6d9d5", alpha = 0.7) +
      annotate("text", x = max(d$time)*0.5,
               y = (L[["B_lo"]] + L[["B_hi"]])/2,
               label = "Bistable window: the insult decides fate",
               size = 3.4, colour = "#8a2418") +
      geom_hline(yintercept = L[["B_lo"]], linetype = 2, colour = "#c0392b") +
      geom_hline(yintercept = L[["B_hi"]], linetype = 2, colour = "#8a2418") +
      geom_line(colour = PAL[1], linewidth = 0.8) +
      labs(x = "Postnatal day", y = "Pathogen load B (1e9 CFU/g)",
           title = "This infant's critical loads and actual load") +
      theme_bw(base_size = 11)
  })

  ## ---- ⑤ lesion ---------------------------------------------------------
  output$pLesion <- renderPlot({
    d <- thin(sim()); cr <- critNa()
    lm <- c(Injury = "Total injury INJ", NECA = "Full-thickness necrosis NECa",
            PNEU = "Intramural gas PNEU", PERF = "Mesenteric perfusion PERF")
    p <- d %>% select(time, Injury, NECA, PNEU, PERF) %>%
      pivot_longer(-time) %>%
      mutate(name = factor(lm[name], levels = unname(lm))) %>%
      ggplot(aes(time, value)) +
      geom_line(colour = PAL[2], linewidth = 0.7) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "Postnatal day", y = NULL,
           title = sprintf("The lesion and the point of no return (NECa_crit = %.3f)", cr)) +
      theme_bw(base_size = 11)
    p
  })

  ## ---- ⑥ endpoints ------------------------------------------------------
  output$pBell <- renderPlot({
    d <- thin(sim())
    ggplot(d, aes(time, BellStage)) +
      geom_step(colour = PAL[2], linewidth = 0.8) +
      scale_y_continuous(breaks = 0:3, limits = c(0, 3)) +
      labs(x = "Postnatal day", y = "Bell stage", title = "Clinical staging over time") +
      theme_bw(base_size = 11)
  })

  output$pOrgan <- renderPlot({
    d <- thin(sim())
    lm <- c(BILI = "Cholestasis (mg/dL)", NIN = "Neuroinflammation index",
            KIN = "Tubular injury index", WtGain_g = "Weight gain (g)")
    long_plot(d, c("BILI", "NIN", "KIN", "WtGain_g"), lm, "Extra-intestinal outcomes")
  })

  output$tEnd <- renderTable({
    d <- sim()
    data.frame(
      Endpoint = c("Peak Bell stage", "Day of diagnosis (Bell II)", "Reached Bell III",
               "Point of no return crossed", "Minimum E", "Final cholestatic bilirubin",
               "Weight gain", "Neuroinflammation index", "Tubular injury index"),
      Value = c(as.character(max(d$BellStage)),
             ifelse(is.na(onset()), "—", sprintf("d%.1f", onset())),
             ifelse(any(d$BellStage >= 3),
                    sprintf("d%.1f", min(d$time[d$BellStage >= 3])), "—"),
             ifelse(any(d$NECA > critNa()), "Yes", "No"),
             sprintf("%.3f", min(d$E)),
             sprintf("%.2f mg/dL", last(d$BILI)),
             sprintf("%.0f g", last(d$WtGain_g)),
             sprintf("%.3f", last(d$NIN)),
             sprintf("%.3f", last(d$KIN))))
  })

  ## ---- ⑦ scenario comparison --------------------------------------------
  armdefs <- list(
    MOM       = list(label = "Mother's own milk", args = list(milk = "MOM")),
    DONOR     = list(label = "Donor milk", args = list(milk = "DONOR")),
    FORMULA   = list(label = "Formula", args = list(milk = "FORMULA")),
    FORM_PROB = list(label = "Formula + probiotic",
                     args = list(milk = "FORMULA", prob_dose = 1.10)),
    FORM_THR  = list(label = "Formula + threshold-raising agent",
                     args = list(milk = "FORMULA", thrboost = 1.6)),
    FORM_BOTH = list(label = "Formula + both",
                     args = list(milk = "FORMULA", prob_dose = 1.10,
                                 thrboost = 1.6)),
    NPO       = list(label = "Nil by mouth + TPN",
                     args = list(milk = "MOM", NPO_ON = 1))
  )

  arms_sim <- eventReactive(input$go, {
    req(length(input$arms) > 0)
    withProgress(message = "Running scenarios...", {
      dplyr::bind_rows(lapply(input$arms, function(k) {
        incProgress(1/length(input$arms), detail = armdefs[[k]]$label)
        a <- c(list(GA = input$GA, BW = input$BW, feedrate = input$feedrate,
                    feedmax = input$feedmax, abx_days = input$abx,
                    mtz_days = input$mtz, ind_days = input$ind,
                    dysbio = input$dysbio,
                    binfant = as.numeric(input$binfant),
                    dx_lag = input$dxlag,
                    insults = insults(), end = input$tend),
               armdefs[[k]]$args)
        do.call(run_nec, a) %>% mutate(arm = armdefs[[k]]$label)
      }))
    })
  }, ignoreNULL = FALSE)

  output$pArms <- renderPlot({
    d <- arms_sim(); req(nrow(d) > 0)
    d <- d[seq(1, nrow(d), by = 4), ]
    lm <- c(E = "Enterocyte integrity E", B = "Pathogen load B",
            Injury = "Total injury INJ", NECA = "Full-thickness necrosis NECa")
    d %>% select(time, arm, E, B, Injury, NECA) %>%
      pivot_longer(-c(time, arm)) %>%
      mutate(name = factor(lm[name], levels = unname(lm))) %>%
      ggplot(aes(time, value, colour = arm)) +
      geom_line(linewidth = 0.7) +
      facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = PAL, name = NULL) +
      labs(x = "Postnatal day", y = NULL,
           title = "Arms compared under an identical precipitating insult") +
      theme_bw(base_size = 11) + theme(legend.position = "bottom")
  })

  output$tArms <- renderTable({
    d <- arms_sim(); req(nrow(d) > 0)
    d %>% group_by(arm) %>%
      summarise(`Peak Bell` = max(BellStage),
                `Day of diagnosis` = ifelse(any(BellStage >= 2),
                                  sprintf("d%.1f", min(time[BellStage >= 2])), "—"),
                `Point of no return crossed` = ifelse(any(NECA > first(NECa_crit)), "Yes", "No"),
                `Minimum E` = sprintf("%.3f", min(E)),
                `Peak B` = sprintf("%.1f", max(B)),
                `Final SCFA` = sprintf("%.1f", last(SCFA)),
                `Weight gain g` = sprintf("%.0f", last(WtGain_g)),
                .groups = "drop")
  })

  ## ---- ⑧ biomarkers -----------------------------------------------------
  output$pBio <- renderPlot({
    d <- thin(sim())
    lm <- c(CRP = "CRP (mg/L)", PLT = "Platelets (1e9/L)", LAC = "Lactate (mmol/L)",
            LPSP = "Endotoxin index", PNEU = "Intramural gas", Translocat = "Jtr (mechanistic)")
    p <- long_plot(d, c("CRP", "PLT", "LAC", "LPSP", "PNEU", "Translocat"), lm,
                   "Clinical markers and mechanistic quantities")
    if (!is.na(onset())) p <- p + geom_vline(xintercept = onset(),
                                            linetype = 2, colour = "#c0392b")
    p
  })

  output$tLead <- renderTable({
    d <- sim(); o <- onset()
    cross <- function(v, thr, below = FALSE) {
      hit <- if (below) which(v < thr) else which(v > thr)
      if (!length(hit)) NA_real_ else d$time[hit[1]]
    }
    mk <- list("CRP > 12"      = cross(d$CRP, 12),
               "Platelets < 100" = cross(d$PLT, 100, TRUE),
               "Lactate > 4"   = cross(d$LAC, 4),
               "Intramural gas" = cross(d$PNEU, 0.30),
               "Jtr > 3 (mechanistic)" = cross(d$Translocat, 3))
    data.frame(
      Marker = names(mk),
      Crossing_day = vapply(mk, function(x) ifelse(is.na(x), "—",
                                             sprintf("d%.2f", x)), ""),
      `Lead_time_h` = vapply(mk, function(x)
        ifelse(is.na(x) || is.na(o), "—", sprintf("%+.1f", (o - x)*24)), ""),
      row.names = NULL)
  })

  ## ---- ⑨ PK -------------------------------------------------------------
  output$pPK <- renderPlot({
    d <- thin(sim(), 2); pp <- par_list()
    d %>% transmute(time,
                    `Ampicillin plasma` = AMPC,
                    `Ampicillin lumen` = AMPC*pp$fl_amp,
                    `Gentamicin plasma` = GENC,
                    `Gentamicin lumen` = GENC*pp$fl_gen,
                    `Metronidazole plasma` = MTZC,
                    `Metronidazole lumen` = MTZC*pp$fl_mtz) %>%
      pivot_longer(-time) %>%
      mutate(drug = sub(" .*", "", name),
             site = ifelse(grepl("lumen", name), "Lumen", "Plasma")) %>%
      ggplot(aes(time, value, colour = site)) +
      geom_line(linewidth = 0.65) +
      facet_wrap(~drug, scales = "free_y") +
      scale_colour_manual(values = c(`Lumen` = PAL[2],
                                     `Plasma` = PAL[1]), name = NULL) +
      labs(x = "Postnatal day", y = "Concentration (mg/L)",
           title = "Plasma and lumen are different exposures") +
      theme_bw(base_size = 11) + theme(legend.position = "bottom")
  })

  output$tPK <- renderTable({
    d <- sim(); pp <- par_list()
    data.frame(
      Drug = c("Ampicillin", "Gentamicin", "Metronidazole",
               "Indomethacin", "Ibuprofen"),
      `f_lum` = c(sprintf("%.2f", pp$fl_amp), sprintf("%.2f", pp$fl_gen),
                  sprintf("%.2f", pp$fl_mtz), "—", "—"),
      `Plasma_Cmax` = c(sprintf("%.1f", max(d$AMPC)), sprintf("%.2f", max(d$GENC)),
                      sprintf("%.1f", max(d$MTZC)), sprintf("%.3f", max(d$INDC)),
                      sprintf("%.2f", max(d$IBUC))),
      `Lumen_Cmax` = c(sprintf("%.1f", max(d$AMPC)*pp$fl_amp),
                      sprintf("%.2f", max(d$GENC)*pp$fl_gen),
                      sprintf("%.1f", max(d$MTZC)*pp$fl_mtz), "—", "—"),
      `Main_target` = c("Commensal anaerobes C (the harmful direction)", "Pathogen B (the helpful direction)",
                    "Commensal anaerobes C (the harmful direction)", "Prostanoids → mesenteric perfusion",
                    "Prostanoids (smaller effect)"))
  })
}

shinyApp(ui, server)
