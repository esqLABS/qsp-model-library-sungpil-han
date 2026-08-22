## =====================================================================
##  CVID QSP MODEL — SHINY DASHBOARD
##  Common Variable Immunodeficiency
##  ---------------------------------------------------------------------
##  10 tabs, organised around the two arms of the disease:
##
##    ARM 1  antibody deficiency  — tabs 2, 3, 4  (replaceable)
##    ARM 2  immune dysregulation — tabs 6, 7     (not replaceable)
##    the layer neither arm reverses — tab 5      (irreversible)
##
##  The controls are deliberately arranged so that the two questions a
##  CVID clinician actually asks are one slider apart:
##      "how much immunoglobulin, given how?"        (ARM 1 panel)
##      "and what do I do about the dysregulation?"  (ARM 2 panel)
##
##  Run:  shiny::runApp("cvid_shiny_app.R")
##  Requires: shiny, mrgsolve, dplyr, tidyr, ggplot2, DT
## =====================================================================

suppressPackageStartupMessages({
  library(shiny); library(mrgsolve); library(dplyr)
  library(tidyr); library(ggplot2); library(DT)
})

## ---------------------------------------------------------------------
##  Load the model. The model file ends in an invisible list, and its
##  pre-runs solve the healthy and CVID baselines, so sourcing it gives
##  us the compiled model plus both sets of initial conditions.
## ---------------------------------------------------------------------
MODEL_FILE <- "cvid_mrgsolve_model.R"
message("Sourcing ", MODEL_FILE, " (compiles the model and solves baselines)...")
BUILD <- source(MODEL_FILE, local = new.env())$value
mod        <- BUILD$mod
CVID_SS    <- BUILD$cvid
HEALTHY_SS <- BUILD$healthy
CVID_GENO  <- list(HEALTHY = 0, DAMAGEON = 1, FCSR = 0.06, FPCD = 0.60,
                   FTACI = 0.35, FTFH = 0.55)

THEME <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom", legend.title = element_blank(),
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(colour = "grey35", size = 10))

PAL <- c("#1565C0", "#C62828", "#2E7D32", "#6A1B9A", "#EF6C00", "#00838F")

## =====================================================================
##  UI
## =====================================================================
ui <- fluidPage(
  titlePanel(HTML(paste0(
    "<b>Common Variable Immunodeficiency (CVID) — QSP simulator</b>",
    "<br><span style='font-size:13px;color:#555'>",
    "ARM 1 antibody deficiency (replaceable) &nbsp;∥&nbsp; ARM 2 immune dysregulation (not replaceable) ",
    "&nbsp;→&nbsp; irreversible structural damage</span>"))),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      h4("Patient"),
      sliderInput("wt", "Body weight (kg)", 40, 120, 70, 5),
      sliderInput("delay", "Diagnostic delay (years)", 0, 15, 2, 1),
      sliderInput("years", "Simulation duration (years)", 5, 30, 20, 5),
      helpText(HTML("<small>During the delay the ratchet turns <b>with no brake</b>.
                     This slider is the most powerful variable in the model.</small>")),

      hr(), h4("ARM 1 · Immunoglobulin replacement"),
      radioButtons("route", "Route",
                   c("IVIG (q3-4 weeks)" = "IV",
                     "SCIG (1-2 per week)" = "SC",
                     "fSCIG (hyaluronidase, q4 weeks)" = "FSC"),
                   selected = "IV"),
      sliderInput("dose", "Dose (mg/kg per 4 weeks, IV-equivalent)",
                  200, 1500, 500, 50),
      conditionalPanel("input.route == 'IV'",
        radioButtons("tau_iv", "Dosing interval", c("q3 weeks" = 21, "q4 weeks" = 28), 28)),
      conditionalPanel("input.route == 'SC'",
        radioButtons("tau_sc", "Dosing frequency",
                     c("Once weekly" = 7, "Twice weekly" = 3.5, "Daily" = 1), 7)),
      checkboxInput("scadj", "Apply the SC dose adjustment factor 1.37 (AUC matching)", TRUE),

      hr(), h4("ARM 2 · Immune dysregulation"),
      sliderInput("dysgeno", "Dysregulation load", 0, 1, 0, 0.05),
      selectInput("pheno", "Dominant phenotype (Chapel phenotype)",
                  c("None (CVID without complications)" = "none",
                    "GLILD / granuloma" = "glild",
                    "Autoimmune cytopenia (ITP/AIHA)" = "cytopenia",
                    "Enteropathy / protein loss" = "enteropathy",
                    "Polyclonal lymphoproliferation" = "lproc")),
      selectInput("geno", "Genotype",
                  c("Unidentified (polygenic)" = "poly",
                    "CTLA4 haploinsufficiency (CHAI)" = "ctla4",
                    "LRBA deficiency" = "lrba",
                    "APDS (PIK3CD GOF)" = "apds")),

      hr(), h4("Therapies"),
      checkboxGroupInput("tx", NULL,
        c("Rituximab" = "rtx",
          "Azathioprine/MMF" = "aza",
          "Prednisone" = "gc",
          "Abatacept" = "aba",
          "Sirolimus" = "siro",
          "Leniolisib (APDS)" = "lenio",
          "JAK inhibitor" = "jak",
          "Eltrombopag" = "elt",
          "Splenectomy" = "splenec",
          "Azithromycin prophylaxis" = "azm",
          "Airway clearance physiotherapy" = "physio")),
      sliderInput("tstart", "ARM 2 therapy start (years)", 0, 20, 5, 1),

      hr(),
      actionButton("go", "Run simulation", class = "btn-primary btn-block"),
      checkboxInput("cmp", "Compare with the standard regimen (IVIG 500 q4w)", TRUE)
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs", type = "tabs",

        ## ---- 1 -----------------------------------------------------
        tabPanel("1 · Patient profile",
          br(), fluidRow(
            column(4, wellPanel(h5("Immune phenotype"),
                                tableOutput("t_pheno"))),
            column(4, wellPanel(h5("Outcomes at 20 years"),
                                tableOutput("t_outcome"))),
            column(4, wellPanel(h5("Model structure"),
                                tableOutput("t_struct")))),
          hr(), h5("EUROclass / Freiburg position"),
          plotOutput("p_euro", height = "300px"),
          hr(),
          HTML("<div style='background:#FFFDE7;border-left:4px solid #F9A825;
                padding:10px'><b>The premise of this model</b><br>
                CVID is two diseases that share one laboratory definition,
                and the therapy we give treats <b>only one</b> of them.
                <b>ARM 1</b> (antibody deficiency → infection) is quantitatively
                predicted by serum IgG and replacement very nearly solves it.
                <b>ARM 2</b> (immune dysregulation → autoimmunity · GLILD · enteropathy · lymphoma)
                is a T-cell and tolerance disorder on which replacement has
                essentially no effect. And it is the latter that kills people
                (Resnick 2012: mortality RR ≈ 11-fold for non-infectious complications).</div>")),

        ## ---- 2 -----------------------------------------------------
        tabPanel("2 · IgG PK",
          br(), plotOutput("p_pk", height = "330px"),
          hr(), fluidRow(
            column(6, plotOutput("p_pk_zoom", height = "300px")),
            column(6, br(), tableOutput("t_pk"),
                   HTML("<small><b>Same monthly dose, different trough.</b>
                    The trough is not set by the monthly amount alone but by <i>how
                    it was divided up</i>. IVIG q4w produces a 600-900 mg/dL
                    sawtooth, while the same AUC as weekly SCIG is nearly flat, so
                    <b>the same milligrams buy a higher trough</b>.
                    </small>")))),

        ## ---- 3 -----------------------------------------------------
        tabPanel("3 · Exposure-Response",
          br(), fluidRow(
            column(6, plotOutput("p_er", height = "340px")),
            column(6, plotOutput("p_er2", height = "340px"))),
          hr(),
          HTML("<div style='background:#E3F2FD;border-left:4px solid #1565C0;
                padding:10px'><b>The slope in Orange 2010 is right; the corollary
                is wrong.</b><br>
                Meta-regression of 17 IVIG studies: the pneumonia rate falls about
                <b>27%</b> for every 100 mg/dL rise in IgG trough. That slope is robust.
                The 'zero pneumonia at 1400 mg/dL' quoted alongside it is an artefact
                of <i>fitting a straight line to a rate</i>. This model keeps the slope
                but implements it as a Hill form, so protection improves <b>without
                ever reaching sterility</b> and a knee of diminishing returns appears
                around 800-1000 mg/dL. The linear form says 'push to 1400 and forget
                it'; the Hill form says 'above that there is almost nothing left to buy'.</div>"),
          hr(), h5("Convexity (Jensen) penalty — the excess risk that comes from the profile shape alone"),
          plotOutput("p_jensen", height = "300px")),

        ## ---- 4 -----------------------------------------------------
        tabPanel("4 · Infection",
          br(), plotOutput("p_inf", height = "330px"),
          hr(), fluidRow(
            column(7, plotOutput("p_divergence", height = "300px")),
            column(5, br(), HTML("<b>The dose-response of pneumonia and of sinusitis are different.</b>
              <br><small>Serum IgG penetrates the airway surface poorly, and
              secretory IgA is <b>replaced by no product at all</b>.
              A <b>floor independent of IgG</b> therefore exists in the mucosal infection rate,
              and patients on replacement still have 2-3 episodes of sinusitis or bronchitis a year.
              Consequence: escalating the dose is worth it for recurrent <b>pneumonia</b>,
              and largely pointless for persistent <b>sinusitis</b>.
              What acts on that floor is azithromycin.</small>")))),

        ## ---- 5 -----------------------------------------------------
        tabPanel("5 · Irreversible damage (The Ratchet)",
          br(), plotOutput("p_ratchet", height = "340px"),
          hr(), fluidRow(
            column(6, plotOutput("p_delay", height = "320px")),
            column(6, plotOutput("p_lung", height = "320px"))),
          hr(),
          HTML("<div style='background:#FFEBEE;border-left:4px solid #B71C1C;
                padding:10px'><b>The ratchet.</b><br>
                Colonisation → IL-8 → neutrophil elastase → airway wall destruction →
                <b>bronchiectasis (irreversible)</b> → reduced clearance → mucus stasis →
                <i>more colonisation</i>. A positive feedback loop that carries an
                irreversible state variable inside it.<br>
                Two patients with the same steady-state trough and everything else identical,
                one diagnosed in year 1 and the other in year 7, end up for that reason alone with a
                <b>permanently different FEV1</b>. Replacement stops the ratchet
                but does not reverse it.<br>
                So the most valuable intervention in CVID is not a better immunoglobulin
                but an <b>earlier serum IgG measurement</b>.</div>")),

        ## ---- 6 -----------------------------------------------------
        tabPanel("6 · ARM 2 · Immune dysregulation",
          br(), plotOutput("p_arm2", height = "340px"),
          hr(), fluidRow(
            column(6, plotOutput("p_cyto", height = "300px")),
            column(6, plotOutput("p_glild", height = "300px"))),
          hr(),
          HTML("<div style='background:#F3E5F5;border-left:4px solid #4A148C;
                padding:10px'><b>A structural null.</b> The GLILD equation contains
                <b>no IgG term directly</b>. Across the whole
                IVIG 300 → 1000 mg/kg range the model moves GLILD activity barely 8%,
                and even that is not a direct effect but the
                <b>indirect</b> route infection → chronic activation → CD21low B → lymphoid aggregates.
                In the same model rituximab + azathioprine moves it 80%
                (about a 10-fold difference). That means GLILD stays essentially
                unchanged even at maximal replacement, and this follows not from a parameter
                but from the <b>structure</b>.</div>")),

        ## ---- 7 -----------------------------------------------------
        tabPanel("7 · The cost of immunosuppression",
          br(), plotOutput("p_isb", height = "330px"),
          hr(), fluidRow(
            column(6, plotOutput("p_tradeoff", height = "310px")),
            column(6, br(), tableOutput("t_tradeoff"),
              HTML("<small><b>Rituximab is relatively cheap in CVID.</b>
              In other diseases the main humoral cost of rituximab is secondary
              hypogammaglobulinaemia. In CVID that cost has <b>already been
              paid and is already being replaced</b>. In the model the same
              B-cell depletion (-98%) lowers IgG by -11% in a normal host but
              by only -3% in a CVID patient on replacement. Rituximab's real
              cost is not humoral but <b>cellular</b>.<br><br>
              There is a cost in the opposite direction too. Every ARM 2 therapy raises
              net immunosuppression, and that comes straight back as ARM 1
              infection risk.</small>")))),

        ## ---- 8 -----------------------------------------------------
        tabPanel("8 · Scenario comparison",
          br(), fluidRow(
            column(3, selectInput("cmp_set", "Comparison set",
              c("Dose (300-1000 mg/kg)" = "dose",
                "Route · interval (AUC matched)" = "route",
                "Diagnostic delay (1/4/7/15 years)" = "delay",
                "4 strategies for refractory ITP" = "itp",
                "4 strategies for GLILD" = "glild"))),
            column(3, selectInput("cmp_y", "Endpoint",
              c("Serum IgG" = "IGG", "Pneumonia rate" = "PNEU_YR",
                "Sinusitis rate" = "SINO_YR", "Bronchiectasis (Reiff)" = "BE",
                "FEV1 %pred" = "FEV1", "DLCO %pred" = "DLCO",
                "GLILD activity" = "GLILD", "Platelets" = "PLT",
                "ARM 2 burden" = "ARM2", "Quality of life" = "QOL"))),
            column(6, br(), helpText("Each set is simulated afresh when you run."))),
          plotOutput("p_cmp", height = "380px"),
          hr(), DTOutput("t_cmp")),

        ## ---- 9 -----------------------------------------------------
        tabPanel("9 · Biomarkers",
          br(), plotOutput("p_bio", height = "360px"),
          hr(), fluidRow(
            column(6, plotOutput("p_bcell", height = "300px")),
            column(6, plotOutput("p_baff", height = "300px"))),
          hr(), HTML("<small><b>The markers replacement masks and the ones it cannot
            mask.</b> Serum IgG becomes uninterpretable once dosing starts. But
            <b>IgA · IgM</b>, <b>switched memory B %</b>, <b>CD21low %</b> and
            <b>free light chains</b> are unaffected by exogenous IgG, so residual
            endogenous B-cell function and the dysregulation trajectory can still be read.
            Soluble <b>BAFF</b> is elevated in CVID because the B-cell sink shrinks,
            and rises further after rituximab.</small>")),

        ## ---- 10 ----------------------------------------------------
        tabPanel("10 · Validation & diagnostics",
          br(), h5("Comparison against published literature anchors (Validation anchors)"),
          DTOutput("t_anchor"),
          hr(), h5("Model diagnostics (34 tests, every one of which can fail)"),
          DTOutput("t_diag"),
          hr(), HTML("<small>Diagnostic D7 was originally written expecting an
            <i>exact zero</i>, and the model <b>refuted</b> that expectation. The GLILD
            equation has no IgG term, but a residual indirect sensitivity of about 8%
            survives through the infection → chronic activation route. It is a
            biologically real pathway, so it is <b>reported</b> rather than removed. D15 (leniolisib),
            with the arbitrary lymph-node shrinkage term removed and only PI3Kδ
            inhibition left, <b>under-predicts</b>: -21% against an observed -39%.
            That too is reported as it stands rather than fitted.</small>"))
      )
    )
  )
)

## =====================================================================
##  SERVER
## =====================================================================
server <- function(input, output, session) {

  ## ---- assemble parameters from the UI ------------------------------
  pars_of <- function() {
    p <- CVID_GENO
    p$WT <- input$wt
    p$VC <- 0.055 * input$wt * 10      # dL
    p$VP <- 0.045 * input$wt * 10
    p$DYSGENO <- input$dysgeno

    ## Chapel phenotypes are largely mutually exclusive: emphasise one arm
    p$KINAA <- 0.006; p$KINLA <- 0.0020; p$KINGUT <- 0.0080
    if (input$pheno == "glild")      p$KINLA  <- 0.0030
    if (input$pheno == "cytopenia")  p$KINAA  <- 0.0300
    if (input$pheno == "enteropathy")p$KINGUT <- 0.0120
    if (input$pheno == "lproc")    { p$KINSPL <- 0.0045; p$KINLAD <- 0.0100 }

    if (input$geno == "ctla4") p$CTLA4G  <- 0.5
    if (input$geno == "lrba")  p$CTLA4G  <- 0.3
    if (input$geno == "apds")  p$PI3KGOF <- 0.8

    t0 <- input$tstart * 365
    if ("aza"     %in% input$tx) { p$AZAON <- 1; p$AZAT0 <- t0 }
    if ("jak"     %in% input$tx) { p$JAKON <- 1; p$JAKT0 <- t0 }
    if ("azm"     %in% input$tx) { p$AZMON <- 1; p$AZMT0 <- input$delay*365 }
    if ("physio"  %in% input$tx)   p$PHYSIOON <- 1
    if ("splenec" %in% input$tx) { p$SPLENEC <- 1; p$TSPLENEC <- t0 }
    if (input$route == "FSC")      p$HYAL <- 1
    p
  }

  ## ---- assemble the event table -------------------------------------
  ev_of <- function(dose = input$dose, route = input$route) {
    yrs <- input$years; wt <- input$wt; d0 <- input$delay * 365
    tau <- switch(route, IV = as.numeric(input$tau_iv),
                         SC = as.numeric(input$tau_sc), FSC = 28)
    ## dose is specified as IV-equivalent mg/kg per 4 weeks
    per_dose <- dose * tau / 28
    if (route == "SC"  && isTRUE(input$scadj)) per_dose <- per_dose / 0.73
    if (route == "FSC" && isTRUE(input$scadj)) per_dose <- per_dose / 0.93
    cmtn <- if (route == "IV") "IGG_C" else "IGG_SC"
    nn <- max(1, floor((yrs*365 - d0)/tau) + 1)
    e <- as.data.frame(ev(amt = per_dose*wt, cmt = cmtn, time = d0,
                          ii = tau, addl = nn - 1, evid = 1))
    t0 <- input$tstart * 365; add <- list()
    if ("rtx" %in% input$tx)
      add[[length(add)+1]] <- as.data.frame(ev(
        amt = 700, cmt = "RTX_C", time = t0, ii = 7, addl = 3))
      if ("rtx" %in% input$tx)
        add[[length(add)+1]] <- as.data.frame(ev(
          amt = 700, cmt = "RTX_C", time = t0 + 182, ii = 182,
          addl = max(0, floor((yrs*365 - t0)/182) - 1)))
    if ("gc" %in% input$tx)
      add[[length(add)+1]] <- as.data.frame(ev(
        amt = 20, cmt = "PRED_C", time = t0, ii = 1,
        addl = max(0, yrs*365 - t0)))
    if ("aba" %in% input$tx)
      add[[length(add)+1]] <- as.data.frame(ev(
        amt = 10*wt, cmt = "ABA_C", time = t0, ii = 28,
        addl = max(0, ceiling((yrs*365 - t0)/28))))
    if ("siro" %in% input$tx)
      add[[length(add)+1]] <- as.data.frame(ev(
        amt = 2000, cmt = "SIRO_C", time = t0, ii = 1,
        addl = max(0, yrs*365 - t0)))
    if ("lenio" %in% input$tx)
      add[[length(add)+1]] <- as.data.frame(ev(
        amt = 70000, cmt = "LENIO_C", time = t0, ii = 0.5,
        addl = max(0, (yrs*365 - t0)*2)))
    if ("elt" %in% input$tx)
      add[[length(add)+1]] <- as.data.frame(ev(
        amt = 50, cmt = "ELT_C", time = t0, ii = 1,
        addl = max(0, yrs*365 - t0)))
    out <- bind_rows(c(list(e), add)); out$ID <- 1
    out %>% arrange(time)
  }

  sim1 <- function(pars, events, yrs, delta = 7) {
    mod %>% param(pars) %>% init(CVID_SS) %>%
      mrgsim(data = events, end = yrs*365, delta = delta,
             atol = 1e-8, rtol = 1e-6) %>% as_tibble()
  }

  R <- eventReactive(input$go, ignoreNULL = FALSE, {
    withProgress(message = "Simulating...", value = 0.3, {
      main <- sim1(pars_of(), ev_of(), input$years) %>% mutate(arm = "Selected regimen")
      incProgress(0.4)
      ref <- if (isTRUE(input$cmp)) {
        p <- CVID_GENO; p$WT <- input$wt
        p$VC <- 0.055*input$wt*10; p$VP <- 0.045*input$wt*10
        p$DYSGENO <- input$dysgeno
        sim1(p, ev_of(500, "IV"), input$years) %>%
          mutate(arm = "Standard IVIG 500 q4w")
      } else NULL
      list(main = main, both = bind_rows(main, ref))
    })
  })

  yr <- function(d) d$time / 365

  ## ================= TAB 1 ==========================================
  output$t_pheno <- renderTable({
    d <- R()$main; l <- tail(d, 1)
    tibble(Metric = c("Serum IgG (mg/dL)", "Switched memory B (% of B)",
                    "CD21low B (% of B)", "Naive B (% of B)",
                    "Total B cells (/uL)", "Soluble BAFF (x normal)"),
           Value = sprintf("%.1f", c(l$IGG, l$SMBPCT, l$CD21PCT,
                                 l$NAIVEPCT, l$BTOTAL, l$BAFF)))
  }, colnames = TRUE)

  output$t_outcome <- renderTable({
    d <- R()$main; l <- tail(d, 1)
    w <- d %>% filter(time >= max(time) - 28)
    tibble(Metric = c("IgG trough (mg/dL)", "Pneumonia (/yr)", "Sinusitis (/yr)",
                    "Bronchiectasis (Reiff 0-18)", "FEV1 (%pred)",
                    "DLCO (%pred)", "Cumulative pneumonia", "Quality of life"),
           Value = sprintf("%.2f", c(min(w$IGG), mean(w$PNEU_YR),
                                 mean(w$SINO_YR), l$BE, l$FEV1, l$DLCO,
                                 l$PNEU_CUM, l$QOL)))
  })

  output$t_struct <- renderTable({
    tibble(Item = c("ODE compartments", "Parameters", "Scenarios", "Diagnostic tests",
                    "Map nodes", "References"),
           Value = c("74", "224", "30", "34 (34 PASS)", "234", "118"))
  })

  output$p_euro <- renderPlot({
    l <- tail(R()$main, 1)
    grid <- expand.grid(smb = seq(0, 12, 0.25), c21 = seq(0, 30, 0.5)) %>%
      mutate(cls = case_when(
        smb < 2 & c21 > 10 ~ "smB− CD21low↑ (Freiburg Ia)\nsplenomegaly · granuloma · autoimmunity",
        smb < 2            ~ "smB− CD21low normal (Freiburg Ib)",
        c21 > 10           ~ "smB+ CD21low↑",
        TRUE               ~ "smB+ CD21low normal (mild phenotype)"))
    ggplot(grid, aes(smb, c21, fill = cls)) + geom_raster(alpha = 0.5) +
      geom_vline(xintercept = 2, linetype = 2) +
      geom_hline(yintercept = 10, linetype = 2) +
      annotate("point", x = min(l$SMBPCT, 12), y = min(l$CD21PCT, 30),
               size = 6, shape = 21, fill = "#C62828", colour = "black",
               stroke = 1.2) +
      annotate("text", x = min(l$SMBPCT, 12) + 0.4, y = min(l$CD21PCT, 30) + 1.4,
               label = "this patient", fontface = "bold", size = 4) +
      scale_fill_manual(values = c("#FFCDD2", "#FFE0B2", "#C5CAE9", "#C8E6C9")) +
      labs(x = "Switched memory B cells (% of B)", y = "CD21low B cells (% of B)",
           title = "EUROclass / Freiburg classification position",
           subtitle = "smB− together with CD21low>10% predicts the dysregulation phenotype") +
      THEME
  })

  ## ================= TAB 2 ==========================================
  output$p_pk <- renderPlot({
    d <- R()$both
    ggplot(d, aes(yr(d), IGG, colour = arm)) + geom_line(linewidth = 0.5) +
      geom_hline(yintercept = c(500, 700, 1000), linetype = 3,
                 colour = "grey45") +
      annotate("text", x = 0.3, y = c(520, 720, 1020), hjust = 0, size = 3,
               colour = "grey35", label = c("500", "700", "1000 mg/dL")) +
      scale_colour_manual(values = PAL) +
      labs(x = "Time (years)", y = "Serum IgG (mg/dL)",
           title = "Immunoglobulin replacement PK — whole course",
           subtitle = "About 5 half-lives (3-6 months) to steady state. During the delay only endogenous IgG is present") +
      THEME
  })

  output$p_pk_zoom <- renderPlot({
    d <- R()$both; te <- max(d$time)
    dz <- d %>% filter(time >= te - 84)
    ggplot(dz, aes((time - (te - 84)), IGG, colour = arm)) +
      geom_line(linewidth = 0.8) +
      scale_colour_manual(values = PAL) +
      labs(x = "Last 12 weeks (days)", y = "Serum IgG (mg/dL)",
           title = "Steady-state profile (last 12 weeks)",
           subtitle = "The amplitude of the sawtooth is the size of the convexity penalty") + THEME
  })

  output$t_pk <- renderTable({
    d <- R()$both; te <- max(d$time)
    d %>% filter(time >= te - 28) %>% group_by(arm) %>%
      summarise(`Mean` = round(mean(IGG)), `trough` = round(min(IGG)),
                `peak` = round(max(IGG)),
                `Swing` = round(max(IGG) - min(IGG)),
                `CL (dL/day)` = round(mean(CLTOT), 2),
                `Pneumonia/yr` = round(mean(PNEU_YR), 4), .groups = "drop")
  })

  ## ================= TAB 3 ==========================================
  ER <- reactive({
    pp <- as.list(param(mod))
    cc <- seq(50, 1600, 10)
    brd <- pp$BRDPOOL*cc/(cc + pp$BRDC50) + pp$BRDENDOG*(3.2/pp$MEMREF)
    ops <- cc*(0.5 + 0.5*brd/(pp$BRDPOOL + pp$BRDENDOG))
    iga <- 0.24
    mp  <- ops*(0.35 + 0.65*iga)
    tibble(IgG = cc,
           pneu = pp$RMAXPNEU/(1 + (ops/pp$C50PNEU)^pp$HILLPNEU),
           sino = pp$RSINOFLR*(1 + 0.9*(1 - iga)) +
                  pp$RMAXSINO/(1 + (mp/pp$C50SINO)^pp$HILLSINO),
           lin  = pmax(0, 0.113*(1 - 0.0027*(cc - 500))))
  })

  output$p_er <- renderPlot({
    e <- ER(); l <- tail(R()$main, 1)
    ggplot(e, aes(IgG)) +
      geom_line(aes(y = pneu, colour = "Hill form (this model)"), linewidth = 1) +
      geom_line(aes(y = lin, colour = "Linear form (Orange extrapolation)"),
                linewidth = 0.8, linetype = 2) +
      geom_point(data = tibble(IgG = 500, y = 0.113), aes(y = y),
                 size = 3, colour = "black") +
      annotate("text", x = 560, y = 0.125, hjust = 0, size = 3.2,
               label = "anchor: 0.113/yr @ 500 mg/dL") +
      geom_vline(xintercept = l$IGG, colour = "#C62828", linetype = 3) +
      geom_vline(xintercept = 1400, colour = "grey60", linetype = 3) +
      annotate("text", x = 1395, y = 0.28, hjust = 1, size = 3.2,
               colour = "grey40", label = "where the linear form reaches 0\n(an artefact)") +
      scale_colour_manual(values = c("#1565C0", "#C62828")) +
      labs(x = "Serum IgG (mg/dL)", y = "Pneumonia (events/yr)",
           title = "ARM 1 exposure-response",
           subtitle = "The two forms almost coincide over 500-900 and say entirely different things outside it") +
      THEME
  })

  output$p_er2 <- renderPlot({
    e <- ER() %>% mutate(slope = c(NA, diff(log(pneu))/diff(IgG))) %>%
      mutate(pct100 = (exp(slope*100) - 1)*100)
    ggplot(e, aes(IgG, pct100)) + geom_line(linewidth = 1, colour = "#1565C0") +
      geom_hline(yintercept = -27, linetype = 2, colour = "#C62828") +
      annotate("text", x = 1450, y = -25, hjust = 1, size = 3.4,
               colour = "#C62828", label = "Orange 2010: −27% / 100 mg/dL") +
      labs(x = "Serum IgG (mg/dL)", y = "Change in pneumonia rate per 100 mg/dL increment (%)",
           title = "The slope varies with concentration",
           subtitle = "The knee of diminishing returns is around 800-1000 mg/dL") + THEME
  })

  output$p_jensen <- renderPlot({
    e <- ER(); m <- 700; a <- 350
    seg <- tibble(x = c(m - a, m + a),
                  y = approx(e$IgG, e$pneu, c(m - a, m + a))$y)
    flat <- approx(e$IgG, e$pneu, m)$y
    ggplot(e %>% filter(IgG > 250, IgG < 1200), aes(IgG, pneu)) +
      geom_line(linewidth = 1, colour = "#1565C0") +
      geom_line(data = seg, aes(x, y), linetype = 2, colour = "#C62828") +
      geom_point(data = seg, aes(x, y), colour = "#C62828", size = 2.5) +
      annotate("point", x = m, y = flat, colour = "#2E7D32", size = 3.5) +
      annotate("point", x = m, y = mean(seg$y), colour = "#C62828", size = 3.5) +
      annotate("segment", x = m, xend = m, y = flat, yend = mean(seg$y),
               arrow = arrow(length = unit(2, "mm"), ends = "both"),
               colour = "black") +
      annotate("text", x = m + 25, y = (flat + mean(seg$y))/2, hjust = 0,
               size = 3.6, fontface = "bold",
               label = sprintf("convexity penalty\n%+.0f%%",
                               100*(mean(seg$y)/flat - 1))) +
      labs(x = "Serum IgG (mg/dL)", y = "Pneumonia (events/yr)",
           title = "Why a flat profile is better at the same AUC — Jensen's inequality",
           subtitle = "Because the curve is convex, a profile that fluctuates about the same mean carries a higher mean risk") +
      THEME
  })

  ## ================= TAB 4 ==========================================
  output$p_inf <- renderPlot({
    d <- R()$both %>%
      select(time, arm, `Pneumonia (/yr)` = PNEU_YR,
             `Sinus and bronchial infection (/yr)` = SINO_YR,
             `Cumulative pneumonia` = PNEU_CUM, `Cumulative invasive infection` = INVAS_CUM) %>%
      pivot_longer(-c(time, arm))
    ggplot(d, aes(time/365, value, colour = arm)) + geom_line(linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = PAL) +
      labs(x = "Time (years)", y = NULL, title = "Infection burden (ARM 1)",
           subtitle = "The step during the delay marks the start of replacement") + THEME
  })

  output$p_divergence <- renderPlot({
    e <- ER() %>% mutate(pneu_rel = pneu/pneu[which.min(abs(IgG - 400))],
                         sino_rel = sino/sino[which.min(abs(IgG - 400))]) %>%
      filter(IgG >= 300, IgG <= 1500) %>%
      select(IgG, `Pneumonia (relative)` = pneu_rel, `Sinusitis (relative)` = sino_rel) %>%
      pivot_longer(-IgG)
    ggplot(e, aes(IgG, value, colour = name)) + geom_line(linewidth = 1.1) +
      scale_colour_manual(values = c("#C62828", "#00838F")) +
      scale_y_continuous(limits = c(0, 1.35)) +
      labs(x = "Serum IgG (mg/dL)", y = "Relative incidence versus IgG 400 mg/dL",
           title = "Escalating the dose abolishes pneumonia but barely touches sinusitis",
           subtitle = "The mucosal floor comes from the absence of secretory IgA and cannot be reached at any dose") +
      THEME
  })

  ## ================= TAB 5 ==========================================
  output$p_ratchet <- renderPlot({
    d <- R()$main %>%
      select(time, `Airway colonisation` = COLON, `Airway inflammation` = AIRINF,
             `Elastase activity` = NEACT,
             `Bronchiectasis (Reiff)` = BE,
             `Annual progression (Reiff/yr)` = BEFLUX, `Mucus load` = MUCUS) %>%
      pivot_longer(-time)
    ggplot(d, aes(time/365, value)) +
      geom_line(linewidth = 0.8, colour = "#B71C1C") +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "Time (years)", y = NULL,
           title = "The stages of the ratchet",
           subtitle = "Once an adequate trough is reached the progression rate falls to 0, and the accumulated damage remains") +
      THEME
  })

  output$p_delay <- renderPlot({
    withProgress(message = "Delay scenarios...", {
      p <- pars_of()
      dd <- lapply(c(1, 4, 7, 15), function(k) {
        d0 <- k*365; tau <- 28
        nn <- max(1, floor((input$years*365 - d0)/tau) + 1)
        e <- as.data.frame(ev(amt = 500*input$wt, cmt = "IGG_C", time = d0,
                              ii = tau, addl = nn - 1, evid = 1)); e$ID <- 1
        sim1(p, e, input$years, delta = 28) %>%
          mutate(lab = sprintf("delay %2d yr", k))
      }) %>% bind_rows()
    })
    ggplot(dd, aes(time/365, FEV1, colour = lab)) + geom_line(linewidth = 0.9) +
      scale_colour_manual(values = PAL) +
      labs(x = "Time (years)", y = "FEV1 (% predicted)",
           title = "The same therapy, a different starting point",
           subtitle = "The troughs are all identical. The whole difference is irreversible") + THEME
  })

  output$p_lung <- renderPlot({
    d <- R()$both %>%
      select(time, arm, `FEV1 %pred` = FEV1, `FVC %pred` = FVC,
             `DLCO %pred` = DLCO) %>% pivot_longer(-c(time, arm))
    ggplot(d, aes(time/365, value, colour = arm, linetype = name)) +
      geom_line(linewidth = 0.8) + scale_colour_manual(values = PAL) +
      labs(x = "Time (years)", y = "% predicted", title = "Lung function",
           subtitle = "DLCO is more sensitive to GLILD/fibrosis than to bronchiectasis") +
      THEME + theme(legend.box = "vertical")
  })

  ## ================= TAB 6 ==========================================
  output$p_arm2 <- renderPlot({
    d <- R()$both %>%
      select(time, arm, `GLILD activity` = GLILD, `Interstitial fibrosis` = FIB,
             `Splenomegaly` = SPLEEN, `Lymphadenopathy` = LAD,
             `ARM 2 total burden` = ARM2, `Survival probability` = SURV) %>%
      pivot_longer(-c(time, arm))
    ggplot(d, aes(time/365, value, colour = arm)) + geom_line(linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = PAL) +
      labs(x = "Time (years)", y = NULL,
           title = "ARM 2 — the disease replacement cannot touch",
           subtitle = "If the two curves overlap, that is precisely the structural null") +
      THEME
  })

  output$p_cyto <- renderPlot({
    d <- R()$main %>%
      select(time, `Platelets (10^9/L)` = PLT, `Haemoglobin (g/dL)` = HGB,
             `Albumin (g/dL)` = ALB) %>% pivot_longer(-time)
    ggplot(d, aes(time/365, value)) + geom_line(linewidth = 0.9,
                                                colour = "#AD1457") +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      labs(x = "Time (years)", y = NULL, title = "Cytopenias and protein loss") + THEME
  })

  output$p_glild <- renderPlot({
    d <- R()$main %>%
      select(time, `Lymphoid aggregates` = LYMPHAGG, `Granuloma` = GRAN,
             `Enteropathy activity` = GUT, `Protein loss (PLE)` = PLE) %>%
      pivot_longer(-time)
    ggplot(d, aes(time/365, value, colour = name)) + geom_line(linewidth = 0.9) +
      scale_colour_manual(values = PAL) +
      labs(x = "Time (years)", y = NULL,
           title = "GLILD components and enteropathy",
           subtitle = "Once PLE switches on, IgG clearance rises and the trough at the same dose collapses") +
      THEME
  })

  ## ================= TAB 7 ==========================================
  output$p_isb <- renderPlot({
    d <- R()$main %>%
      select(time, `Net immunosuppression burden` = ISBURDEN,
             `Susceptibility multiplier` = SUSC, `Pneumonia (/yr)` = PNEU_YR,
             `IgG clearance (dL/day)` = CLTOT) %>% pivot_longer(-time)
    ggplot(d, aes(time/365, value)) + geom_line(linewidth = 0.9,
                                                colour = "#4527A0") +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "Time (years)", y = NULL,
           title = "The route by which ARM 2 therapy comes back as ARM 1 disease",
           subtitle = "Net immunosuppression raises susceptibility, and susceptibility raises the infection rate") +
      THEME
  })

  output$p_tradeoff <- renderPlot({
    withProgress(message = "Computing the rituximab asymmetry...", {
      rtxev <- as.data.frame(ev(amt = 700, cmt = "RTX_C", time = 0,
                                ii = 7, addl = 3)); rtxev$ID <- 1
      nrm <- mod %>% param(HEALTHY = 1, DAMAGEON = 0) %>% init(HEALTHY_SS) %>%
        mrgsim(data = rtxev, end = 3*365, delta = 7) %>% as_tibble() %>%
        mutate(arm = "Normal host")
      sc <- as.data.frame(ev(amt = 500/4/0.73*input$wt, cmt = "IGG_SC",
                             time = 0, ii = 7, addl = 155))
      cv <- mod %>% param(CVID_GENO) %>% param(DAMAGEON = 0) %>%
        init(CVID_SS) %>%
        mrgsim(data = (bind_rows(rtxev, sc) %>% mutate(ID = 1) %>%
                       arrange(time)),
               end = 3*365, delta = 7) %>% as_tibble() %>%
        mutate(arm = "CVID (on replacement)")
      dd <- bind_rows(nrm, cv) %>% group_by(arm) %>%
        mutate(IgG_rel = 100*IGG/IGG[which.min(abs(time - 180))],
               B_rel = 100*BTOTAL/BTOTAL[1]) %>% ungroup()
    })
    ggplot(dd, aes(time/365)) +
      geom_line(aes(y = B_rel, colour = arm, linetype = "Total B cells"),
                linewidth = 0.9) +
      geom_line(aes(y = IgG_rel, colour = arm, linetype = "Serum IgG"),
                linewidth = 0.9) +
      scale_colour_manual(values = c("#C62828", "#1565C0")) +
      labs(x = "Time after rituximab (years)", y = "Relative to baseline (%)",
           title = "Rituximab: the same depletion, a different cost",
           subtitle = "In CVID the humoral cost has already been paid and is already being replaced") +
      THEME + theme(legend.box = "vertical")
  })

  output$t_tradeoff <- renderTable({
    d <- R()$main; w <- d %>% filter(time >= max(time) - 28)
    tibble(Item = c("Net immunosuppression burden", "Susceptibility multiplier",
                    "Cumulative immunosuppression exposure", "Cumulative invasive infection",
                    "ARM 2 burden", "Model survival probability"),
           Value = sprintf("%.3f", c(mean(w$ISBURDEN), mean(w$SUSC),
                                 tail(d$ISCUM,1), tail(d$INVAS_CUM,1),
                                 tail(d$ARM2,1), tail(d$SURV,1))))
  })

  ## ================= TAB 8 ==========================================
  CMP <- reactive({
    p <- pars_of(); yrs <- input$years; wt <- input$wt
    mkiv <- function(dose, tau = 28, d0 = input$delay*365) {
      nn <- max(1, floor((yrs*365 - d0)/tau) + 1)
      e <- as.data.frame(ev(amt = dose*tau/28*wt, cmt = "IGG_C", time = d0,
                            ii = tau, addl = nn - 1, evid = 1)); e$ID <- 1; e
    }
    mksc <- function(dose, tau) {
      d0 <- input$delay*365
      nn <- max(1, floor((yrs*365 - d0)/tau) + 1)
      e <- as.data.frame(ev(amt = dose*tau/28/0.73*wt, cmt = "IGG_SC",
                            time = d0, ii = tau, addl = nn - 1, evid = 1))
      e$ID <- 1; e
    }
    withProgress(message = "Simulating the comparison set...", {
      switch(input$cmp_set,
        dose = lapply(c(300, 400, 500, 600, 800, 1000), function(k)
          sim1(p, mkiv(k), yrs, 14) %>% mutate(lab = paste0(k, " mg/kg q4w"))),
        route = list(
          sim1(p, mkiv(500, 28), yrs, 2) %>% mutate(lab = "IVIG q4w"),
          sim1(p, mkiv(500, 21), yrs, 2) %>% mutate(lab = "IVIG q3w"),
          sim1(p, mksc(500, 7),  yrs, 2) %>% mutate(lab = "SCIG weekly"),
          sim1(p, mksc(500, 3.5),yrs, 2) %>% mutate(lab = "SCIG twice weekly"),
          sim1(modifyList(p, list(HYAL = 1)), mksc(500, 28), yrs, 2) %>%
            mutate(lab = "fSCIG q4w")),
        delay = lapply(c(1, 4, 7, 15), function(k)
          sim1(p, mkiv(500, 28, k*365), yrs, 14) %>%
            mutate(lab = sprintf("delay %2d yr", k))),
        itp = { pi <- modifyList(p, list(DYSGENO = 0.6, KINAA = 0.030,
                                         KINLA = 0.004))
          t0 <- input$tstart*365; base <- mkiv(500)
          add <- function(x) (bind_rows(base, as.data.frame(x)) %>%
                              mutate(ID = 1) %>% arrange(time))
          list(
            sim1(pi, base, yrs, 14) %>% mutate(lab = "IgG replacement only"),
            sim1(pi, add(ev(amt = 20, cmt = "PRED_C", time = t0, ii = 1,
                            addl = yrs*365 - t0)), yrs, 14) %>%
              mutate(lab = "Prednisone"),
            sim1(pi, add(ev(amt = 700, cmt = "RTX_C", time = t0, ii = 365,
                            addl = 14)), yrs, 14) %>%
              mutate(lab = "Rituximab"),
            sim1(modifyList(pi, list(SPLENEC = 1, TSPLENEC = t0)), base,
                 yrs, 14) %>% mutate(lab = "Splenectomy"),
            sim1(pi, add(ev(amt = 50, cmt = "ELT_C", time = t0, ii = 1,
                            addl = yrs*365 - t0)), yrs, 14) %>%
              mutate(lab = "Eltrombopag")) },
        glild = { pg <- modifyList(p, list(DYSGENO = 0.7, KINLA = 0.003))
          t0 <- input$tstart*365; base <- mkiv(500)
          add <- function(x) (bind_rows(base, as.data.frame(x)) %>%
                              mutate(ID = 1) %>% arrange(time))
          list(
            sim1(pg, base, yrs, 14) %>% mutate(lab = "IgG replacement only"),
            sim1(pg, add(ev(amt = 20, cmt = "PRED_C", time = t0, ii = 1,
                            addl = yrs*365 - t0)), yrs, 14) %>%
              mutate(lab = "Prednisone"),
            sim1(modifyList(pg, list(AZAON = 1, AZAT0 = t0)),
                 add(ev(amt = 700, cmt = "RTX_C", time = t0, ii = 182,
                        addl = 30)), yrs, 14) %>%
              mutate(lab = "Rituximab+AZA"),
            sim1(modifyList(pg, list(JAKON = 1, JAKT0 = t0)), base,
                 yrs, 14) %>% mutate(lab = "JAK inhibitor")) }
      ) %>% bind_rows()
    })
  })

  output$p_cmp <- renderPlot({
    d <- CMP(); yv <- input$cmp_y
    ggplot(d, aes(time/365, .data[[yv]], colour = lab)) +
      geom_line(linewidth = 0.85) +
      scale_colour_manual(values = rep(PAL, 3)) +
      labs(x = "Time (years)", y = yv,
           title = paste("Scenario comparison —", yv)) + THEME
  })

  output$t_cmp <- renderDT({
    d <- CMP(); te <- max(d$time)
    d %>% group_by(lab) %>%
      summarise(`IgG trough` = round(min(IGG[time >= te - 28])),
                `Swing` = round(diff(range(IGG[time >= te - 28]))),
                `Pneumonia/yr` = round(mean(PNEU_YR[time >= te - 28]), 4),
                `Sinusitis/yr` = round(mean(SINO_YR[time >= te - 28]), 2),
                `Reiff` = round(last(BE), 2), `FEV1` = round(last(FEV1), 1),
                `DLCO` = round(last(DLCO), 1),
                `GLILD` = round(last(GLILD), 3),
                `Platelets` = round(last(PLT)),
                `ARM2` = round(last(ARM2), 2),
                `QoL` = round(last(QOL), 3),
                `Survival` = round(last(SURV), 3), .groups = "drop") %>%
      datatable(options = list(dom = "t", pageLength = 12), rownames = FALSE)
  })

  ## ================= TAB 9 ==========================================
  output$p_bio <- renderPlot({
    d <- R()$main %>%
      select(time, `Serum IgG (mg/dL)` = IGG,
             `Effective opsonin (mg/dL equivalent)` = OPSONIN,
             `Mucosal protection (mg/dL equivalent)` = MUCPROT,
             `Soluble BAFF (x normal)` = BAFF,
             `Switched memory B (%)` = SMBPCT, `CD21low B (%)` = CD21PCT) %>%
      pivot_longer(-time)
    ggplot(d, aes(time/365, value)) + geom_line(linewidth = 0.85,
                                                colour = "#00695C") +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "Time (years)", y = NULL, title = "Biomarker trajectories",
           subtitle = "The gap between opsonin and mucosal protection is the IgA void") + THEME
  })

  output$p_bcell <- renderPlot({
    d <- R()$main %>%
      select(time, `Naive B` = NAIVEPCT, `Switched memory B` = SMBPCT,
             `CD21low B` = CD21PCT) %>% pivot_longer(-time)
    ggplot(d, aes(time/365, value, fill = name)) +
      geom_area(position = "stack", alpha = 0.85) +
      scale_fill_manual(values = c("#90CAF9", "#EF9A9A", "#A5D6A7")) +
      labs(x = "Time (years)", y = "B-cell subset (% of B)",
           title = "B-cell compartment composition") + THEME
  })

  output$p_baff <- renderPlot({
    d <- R()$main %>%
      select(time, `BAFF` = BAFF, `Total B cells (/uL)` = BTOTAL) %>%
      pivot_longer(-time)
    ggplot(d, aes(time/365, value)) + geom_line(linewidth = 0.9,
                                                colour = "#00838F") +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      labs(x = "Time (years)", y = NULL,
           title = "BAFF and the B-cell sink",
           subtitle = "As the sink shrinks BAFF rises — 2-4 fold after rituximab") + THEME
  })

  ## ================= TAB 10 =========================================
  output$t_anchor <- renderDT({
    datatable(BUILD$anchors, options = list(dom = "t", pageLength = 20),
              rownames = FALSE)
  })
  output$t_diag <- renderDT({
    BUILD$diagnostics %>%
      mutate(pass = ifelse(pass, "PASS", "FAIL")) %>%
      datatable(options = list(dom = "tp", pageLength = 15), rownames = FALSE)
  })
}

shinyApp(ui, server)
