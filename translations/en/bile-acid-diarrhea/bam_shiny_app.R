# =====================================================================
#  Bile Acid Diarrhoea / Malabsorption (BAM) — QSP Shiny dashboard
# ---------------------------------------------------------------------
#  The whole app is organised around ONE control panel idea: the disease
#  is written as two numbers, and everything else is a consequence.
#
#      phi    surviving ileal ASBT absorptive capacity  (slider)
#      kappa  ileal FGF19 sensor gain                   (slider)
#
#  Tab 2 ("The central experiment") lets the user FREEZE the FGF19 ->
#  CYP7A1 loop and watch the diarrhoea disappear while the malabsorption
#  stays exactly where it was.  That is the point of the whole model.
#
#  Run:  shiny::runApp("bam_shiny_app.R")
#  Requires: shiny, mrgsolve, dplyr, tidyr, ggplot2, DT
#  The model is defined in bam_mrgsolve_model.R and sourced here.
# =====================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

source("bam_mrgsolve_model.R")   # defines mod, ev_*, bam_readout, sim_bam

THEME <- theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey93"),
        legend.position = "bottom")

PAL <- c(healthy = "#4f9d69", patient = "#b3272d", treated = "#1f6fb2",
         frozen  = "#8e6fc4", live = "#cc8400")

# ---------------------------------------------------------------------
ui <- fluidPage(
  titlePanel(HTML(paste0(
    "<b>Bile Acid Diarrhoea / Bile Acid Malabsorption QSP Dashboard</b>",
    "<br><span style='font-size:14px;color:#555'>",
    "Bile Acid Diarrhoea — the loop measures what is ABSORBED, ",
    "the symptom is made by what is SPILLED</span>"))),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("① The two axes of the lesion"),
      sliderInput("phi", HTML("<b>&phi;</b> — surviving ileal ASBT absorptive capacity"),
                  min = 0.02, max = 1.0, value = 0.20, step = 0.01),
      sliderInput("kappa", HTML("<b>&kappa;</b> — ileal FGF19 sensor gain"),
                  min = 0.05, max = 1.2, value = 1.00, step = 0.05),
      helpText(HTML(paste0(
        "<small>Type 1 = &phi;&darr; · Type 2 = &kappa;&darr; · ",
        "Type 3 = cholecystectomy · Type 4 = drug-induced</small>"))),
      hr(),
      h4("② Structural modifiers"),
      checkboxInput("chole", "Cholecystectomy", FALSE),
      sliderInput("entmass", "Ileal epithelial functional mass (ENTMASS)",
                  min = 0.3, max = 1.0, value = 1.0, step = 0.05),
      sliderInput("bai", HTML("Microbial 7&alpha;-dehydroxylating capacity (BAI, relative)"),
                  min = 0.1, max = 2.5, value = 1.0, step = 0.1),
      sliderInput("aTR", "Transit-absorption positive feedback gain (aTR)",
                  min = 0.4, max = 3.6, value = 1.55, step = 0.05),
      hr(),
      h4("③ Treatment"),
      selectInput("seq", "Bile acid sequestrant (class A)",
                  c("None" = "0", "colesevelam 1.875 g BID" = "1875_2",
                    "colesevelam 3.75 g BID" = "3750_2",
                    "colestyramine 4 g TID" = "4000_3")),
      selectInput("fxr", "FXR agonist (class B)",
                  c("None" = "0", "obeticholic acid 25 mg QD" = "oca25",
                    "obeticholic acid 10 mg QD" = "oca10",
                    "tropifexor 90 µg QD" = "tro90")),
      selectInput("elo", "ASBT inhibitor (class C)",
                  c("None" = "0", "elobixibat 5 mg QD" = "5",
                    "elobixibat 10 mg QD" = "10",
                    "elobixibat 15 mg QD" = "15")),
      selectInput("tra", "Transit modifier (class D)",
                  c("None" = "0", "loperamide 2 mg BID" = "lop2",
                    "loperamide 4 mg BID" = "lop4",
                    "ondansetron 4 mg BID" = "ond4")),
      checkboxInput("rif", "rifaximin 550 mg TID (class E)", FALSE),
      hr(),
      numericInput("days", "Simulation duration (days)", 56, min = 21, max = 120),
      actionButton("go", "Run simulation", class = "btn-primary btn-block"),
      hr(),
      helpText(HTML(paste0(
        "<small>The provenance and calibration basis of every parameter is in the ",
        "PROVENANCE block at the foot of <code>bam_mrgsolve_model.R</code>, ",
        "and the verification figures are in <code>bam_verification_output.txt</code>.",
        "</small>")))
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs", type = "tabs",

        # ---------------------------------------------------------
        tabPanel(
          "1 · Patient profile",
          br(),
          fluidRow(
            column(4, wellPanel(h4("Diagnostic panel"), tableOutput("dxpanel"))),
            column(4, wellPanel(h4("Symptoms"), tableOutput("sxpanel"))),
            column(4, wellPanel(h4("Systemic effects"), tableOutput("orgpanel")))
          ),
          h4("Subtype classification — two tests read the two lesions separately"),
          plotOutput("quadrant", height = "420px"),
          helpText(HTML(paste0(
            "The horizontal axis, SeHCAT, reads <b>&phi;</b> alone; the vertical axis, C4, reads ",
            "<b>S = f(&phi;,&kappa;)</b>. Top left = classic type 1, top right = <b>SeHCAT-negative bile acid diarrhoea</b> (type 2), ",
            "bottom left = compensated malabsorption (asymptomatic), bottom right = normal.")))
        ),

        # ---------------------------------------------------------
        tabPanel(
          "2 · The central experiment (freezing the feedback)",
          br(),
          h4("The same malabsorption, with the feedback switched on and off"),
          helpText(HTML(paste0(
            "On the left is the physiological patient with the FGF19&rarr;CYP7A1 feedback alive; ",
            "on the right is the same patient with the feedback <b>frozen</b> at the normal operating point. ",
            "The malabsorption (SeHCAT) collapses identically in both cases, but ",
            "the diarrhoea appears only on the side where the feedback is alive."))),
          plotOutput("frozen", height = "460px"),
          br(),
          DT::dataTableOutput("frozentab")
        ),

        # ---------------------------------------------------------
        tabPanel(
          "3 · Bile acid PK — the enterohepatic circulation",
          br(),
          fluidRow(
            column(6, plotOutput("pkpool", height = "330px")),
            column(6, plotOutput("pkduo", height = "330px"))
          ),
          fluidRow(
            column(6, plotOutput("pkile", height = "330px")),
            column(6, plotOutput("pkserum", height = "330px"))
          )
        ),

        # ---------------------------------------------------------
        tabPanel(
          "4 · The sensor axis (FXR–FGF19–CYP7A1)",
          br(),
          fluidRow(
            column(6, plotOutput("fgf", height = "330px")),
            column(6, plotOutput("cyp", height = "330px"))
          ),
          fluidRow(
            column(6, plotOutput("c4plot", height = "330px")),
            column(6, plotOutput("loopgain", height = "330px"))
          ),
          helpText(HTML(paste0(
            "The bottom right panel is the closed-form steady-state solution ",
            "<code>S = Smax / (1 + (a&middot;&kappa;&middot;x&middot;S)^h)</code>, ",
            "x = f/(1-f). Because the feedback is logarithmically stiff, a 65-fold fall in the recycling ",
            "ratio raises synthesis by only about 4-fold.")))
        ),

        # ---------------------------------------------------------
        tabPanel(
          "5 · The colon — the 3 mM threshold",
          br(),
          fluidRow(
            column(7, plotOutput("colspecies", height = "360px")),
            column(5, plotOutput("threshold", height = "360px"))
          ),
          helpText(HTML(paste0(
            "The secretory driving force D<sub>s</sub> is the potency-weighted sum of the <b>free</b> bile acids ",
            "(CA 0.15 / CDCA 1.00 / DCA 1.15 / LCA 0.30). ",
            "The 3 mM threshold comes directly from the human colonic perfusion experiments of Mekjian 1971 ",
            "(PMID 4938344)."))),
          plotOutput("transit", height = "320px")
        ),

        # ---------------------------------------------------------
        tabPanel(
          "6 · Clinical endpoints",
          br(),
          fluidRow(
            column(6, plotOutput("stoolplot", height = "330px")),
            column(6, plotOutput("bmplot", height = "330px"))
          ),
          fluidRow(
            column(6, plotOutput("fatplot", height = "330px")),
            column(6, plotOutput("ldlplot", height = "330px"))
          )
        ),

        # ---------------------------------------------------------
        tabPanel(
          "7 · Comparison of the drug classes",
          br(),
          h4("The same patient, the entry points of five classes"),
          plotOutput("drugclass", height = "440px"),
          br(),
          DT::dataTableOutput("drugtab"),
          helpText(HTML(paste0(
            "<b>A sequestrant raises by itself the very driving force it is meant to treat</b> ",
            "— bind the luminal bile acids and the sensor ligand disappears too, so C4 rises. ",
            "An FXR agonist gives that signal back, which is why the combination is synergistic.")))
        ),

        # ---------------------------------------------------------
        tabPanel(
          "8 · The 100 cm rule",
          br(),
          h4("Resection length → phenotype switch is a derived result"),
          helpText(HTML(paste0(
            "There are only three inputs: the distal skew of ASBT density (decay length 60 cm), ",
            "the CYP7A1 compensation ceiling (6-fold) and the effective critical micellar concentration (1.5 mM). ",
            "No resection outcome data were used at all."))),
          plotOutput("rule100", height = "480px"),
          br(),
          DT::dataTableOutput("ruletab")
        ),

        # ---------------------------------------------------------
        tabPanel(
          "9 · Biomarkers & the microbiota",
          br(),
          fluidRow(
            column(6, plotOutput("baiplot", height = "340px")),
            column(6, plotOutput("biomk", height = "340px"))
          ),
          helpText(HTML(paste0(
            "The microbiota barely changes the <b>quantity</b> of bile acid reaching the colon; it changes its ",
            "<b>potency</b>. So the same faecal bile acid output can go with symptoms differing ",
            "up to 3-fold — the second reason the test and the symptom diverge."))),
          plotOutput("sehcatcurve", height = "330px")
        ),

        # ---------------------------------------------------------
        tabPanel(
          "10 · Scenario library",
          br(),
          h4("The 14 predefined scenarios"),
          actionButton("runlib", "Run the whole library (takes minutes)",
                       class = "btn-warning"),
          br(), br(),
          DT::dataTableOutput("libtab"),
          br(),
          helpText(HTML(paste0(
            "The figures in this table are reproduced independently by <code>bam_verify_python.py</code> ",
            "(<code>bam_verification_output.txt</code>).")))
        )
      )
    )
  )
)

# ---------------------------------------------------------------------
server <- function(input, output, session) {

  build_ev <- function() {
    e <- NULL
    add <- function(a, b) if (is.null(a)) b else a + b
    d <- input$days
    if (input$seq != "0") {
      pr <- strsplit(input$seq, "_")[[1]]
      e <- add(e, ev_seq(as.numeric(pr[1]), as.numeric(pr[2]), d))
    }
    if (input$fxr == "oca25") e <- add(e, ev_oca(25, d))
    if (input$fxr == "oca10") e <- add(e, ev_oca(10, d))
    if (input$fxr == "tro90") e <- add(e, ev_tro(0.090, d))
    if (input$elo != "0")     e <- add(e, ev_elo(as.numeric(input$elo), d))
    if (input$tra == "lop2")  e <- add(e, ev_lop(2, d))
    if (input$tra == "lop4")  e <- add(e, ev_lop(4, d))
    if (input$tra == "ond4")  e <- add(e, ev_ond(4, d))
    if (isTRUE(input$rif))    e <- add(e, ev_rif(550, d))
    e
  }

  pars <- reactive({
    p <- list(phi = input$phi, kappa = input$kappa,
              ENTMASS = input$entmass, aTR = input$aTR,
              kbai = 0.16 * input$bai)
    if (isTRUE(input$chole)) p$GBfun <- 0
    p
  })

  # ---- main simulations, recomputed only on button press -----------
  simset <- eventReactive(input$go, {
    withProgress(message = "Simulating...", value = 0, {
      p <- pars(); e <- build_ev()
      incProgress(0.15, detail = "healthy control")
      h <- do.call(sim_bam, c(list(days = input$days, events = NULL),
                              list(phi = 1, kappa = 1)))
      incProgress(0.3, detail = "patient, untreated")
      u <- do.call(sim_bam, c(list(days = input$days, events = NULL), p))
      incProgress(0.3, detail = "patient, treated")
      tr <- if (is.null(e)) u else
        do.call(sim_bam, c(list(days = input$days, events = e), p))
      incProgress(0.25, detail = "feedback frozen")
      fz <- do.call(sim_bam, c(list(days = input$days, events = NULL),
                               c(p, list(FEEDBACK = 0))))
      list(healthy = h, untreated = u, treated = tr, frozen = fz)
    })
  }, ignoreNULL = FALSE)

  ro <- reactive({
    s <- simset()
    list(healthy = bam_readout(s$healthy), untreated = bam_readout(s$untreated),
         treated = bam_readout(s$treated), frozen = bam_readout(s$frozen))
  })

  last_day <- function(out) {
    d <- as.data.frame(out); te <- max(d$time)
    dplyr::filter(d, time >= te - 48) %>% mutate(tod = time - (te - 48))
  }

  # ================= TAB 1 =========================================
  output$dxpanel <- renderTable({
    r <- ro()
    data.frame(
      Item = c("SeHCAT 7-day retention (%)", "Fasting serum C4 (ng/mL)",
               "Fasting FGF19 (pg/mL)", "Faecal bile acids (µmol/day)",
               "Bile acid pool (g)", "Single-pass ileal conservation f"),
      Patient = c(sprintf("%.1f", r$treated$SeHCAT_pct),
               sprintf("%.1f", r$treated$C4),
               sprintf("%.0f", r$treated$FGF19),
               sprintf("%.0f", r$treated$fecBA_umold),
               sprintf("%.2f", r$treated$pool_g),
               sprintf("%.4f", 1 - r$treated$kturn_day / 4.5)),
      Normal = c(sprintf("%.1f", r$healthy$SeHCAT_pct),
               sprintf("%.1f", r$healthy$C4),
               sprintf("%.0f", r$healthy$FGF19),
               sprintf("%.0f", r$healthy$fecBA_umold),
               sprintf("%.2f", r$healthy$pool_g), "~0.97"))
  }, striped = TRUE, width = "100%")

  output$sxpanel <- renderTable({
    r <- ro()
    data.frame(
      Item = c("Stool water (mL/day)", "Bowel movements (per day)",
               "Bristol score", "Colonic secretory driving force (mM-eq)",
               "Colonic transit (relative)"),
      Patient = c(sprintf("%.0f", r$treated$stool_mL),
               sprintf("%.2f", r$treated$BM_per_day),
               sprintf("%.1f", r$treated$bristol),
               sprintf("%.2f", r$treated$Ds),
               sprintf("%.2f", r$treated$TRANS)),
      Normal = c(sprintf("%.0f", r$healthy$stool_mL),
               sprintf("%.2f", r$healthy$BM_per_day),
               sprintf("%.1f", r$healthy$bristol),
               sprintf("%.2f", r$healthy$Ds), "1.00"))
  }, striped = TRUE, width = "100%")

  output$orgpanel <- renderTable({
    r <- ro()
    data.frame(
      Item = c("Faecal fat (g/day)", "Peak duodenal [BA] (mM)",
               "LDL-C (mg/dL)", "Urinary oxalate (mg/day)"),
      Patient = c(sprintf("%.1f", r$treated$fecfat_g),
               sprintf("%.1f", r$treated$Cduo_peak),
               sprintf("%.0f", r$treated$LDLC),
               sprintf("%.0f", r$treated$UOX)),
      Normal = c(sprintf("%.1f", r$healthy$fecfat_g),
               sprintf("%.1f", r$healthy$Cduo_peak),
               sprintf("%.0f", r$healthy$LDLC),
               sprintf("%.0f", r$healthy$UOX)))
  }, striped = TRUE, width = "100%")

  output$quadrant <- renderPlot({
    r <- ro()
    pt <- data.frame(x = r$treated$SeHCAT_pct, y = r$treated$C4,
                     lab = "Current patient")
    hc <- data.frame(x = r$healthy$SeHCAT_pct, y = r$healthy$C4,
                     lab = "Healthy control")
    ggplot() +
      annotate("rect", xmin = 0, xmax = 15, ymin = 48, ymax = Inf,
               fill = "#ffd9b3", alpha = .5) +
      annotate("rect", xmin = 15, xmax = Inf, ymin = 48, ymax = Inf,
               fill = "#d9c9f0", alpha = .5) +
      annotate("rect", xmin = 0, xmax = 15, ymin = 0, ymax = 48,
               fill = "#e0f0d8", alpha = .5) +
      annotate("text", x = 7, y = 90, label = "Type 1\n(ileal lesion)", size = 5) +
      annotate("text", x = 40, y = 90,
               label = "Type 2 — SeHCAT-negative\nbile acid diarrhoea", size = 5) +
      annotate("text", x = 7, y = 20, label = "Compensated malabsorption\n(few symptoms)", size = 4.5) +
      annotate("text", x = 40, y = 20, label = "Normal", size = 5) +
      geom_vline(xintercept = 15, linetype = 2) +
      geom_hline(yintercept = 48, linetype = 2) +
      geom_point(data = hc, aes(x, y), size = 5, colour = PAL["healthy"]) +
      geom_point(data = pt, aes(x, y), size = 6, colour = PAL["patient"]) +
      geom_text(data = pt, aes(x, y, label = lab), vjust = -1.4, size = 4.5) +
      scale_x_continuous("SeHCAT 7-day retention (%)  →  reads φ",
                         limits = c(0, 60)) +
      scale_y_continuous("Fasting serum C4 (ng/mL)  →  reads S = f(φ, κ)",
                         limits = c(0, 120)) + THEME
  })

  # ================= TAB 2 =========================================
  output$frozen <- renderPlot({
    r <- ro()
    d <- bind_rows(
      data.frame(cond = "Feedback LIVE", metric = "SeHCAT (%)",
                 v = r$untreated$SeHCAT_pct),
      data.frame(cond = "Feedback FROZEN", metric = "SeHCAT (%)",
                 v = r$frozen$SeHCAT_pct),
      data.frame(cond = "Feedback LIVE", metric = "Colonic bile acid load (µmol/d)",
                 v = r$untreated$fecBA_umold),
      data.frame(cond = "Feedback FROZEN", metric = "Colonic bile acid load (µmol/d)",
                 v = r$frozen$fecBA_umold),
      data.frame(cond = "Feedback LIVE", metric = "Serum C4 (ng/mL)",
                 v = r$untreated$C4),
      data.frame(cond = "Feedback FROZEN", metric = "Serum C4 (ng/mL)",
                 v = r$frozen$C4),
      data.frame(cond = "Feedback LIVE", metric = "Stool water (mL/d)",
                 v = r$untreated$stool_mL),
      data.frame(cond = "Feedback FROZEN", metric = "Stool water (mL/d)",
                 v = r$frozen$stool_mL))
    ggplot(d, aes(cond, v, fill = cond)) +
      geom_col(width = .6) +
      geom_text(aes(label = sprintf("%.1f", v)), vjust = -.4, size = 4.2) +
      facet_wrap(~metric, scales = "free_y", nrow = 1) +
      scale_fill_manual(values = c("Feedback LIVE" = unname(PAL["live"]),
                                   "Feedback FROZEN" = unname(PAL["frozen"]))) +
      labs(x = NULL, y = NULL, fill = NULL) +
      expand_limits(y = 0) + THEME
  })

  output$frozentab <- DT::renderDataTable({
    r <- ro()
    tab <- bind_rows(
      cbind(Condition = "Feedback LIVE", r$untreated),
      cbind(Condition = "Feedback FROZEN", r$frozen),
      cbind(Condition = "Healthy control", r$healthy))
    DT::datatable(tab, options = list(dom = "t", scrollX = TRUE),
                  rownames = FALSE) %>%
      DT::formatRound(2:ncol(tab), 2)
  })

  # ================= TAB 3 =========================================
  plot_ts <- function(var, ylab, title) {
    renderPlot({
      s <- simset()
      d <- bind_rows(
        cbind(grp = "Healthy", last_day(s$healthy)),
        cbind(grp = "Patient (untreated)", last_day(s$untreated)),
        cbind(grp = "Patient (treated)", last_day(s$treated)))
      ggplot(d, aes(tod, .data[[var]], colour = grp)) +
        geom_line(linewidth = .9) +
        scale_colour_manual(values = c("Healthy" = unname(PAL["healthy"]),
                                       "Patient (untreated)" = unname(PAL["patient"]),
                                       "Patient (treated)" = unname(PAL["treated"]))) +
        scale_x_continuous("Last 48 hours (h)", breaks = seq(0, 48, 8)) +
        labs(y = ylab, title = title, colour = NULL) + THEME
    })
  }
  output$pkpool   <- plot_ts("POOLtot", "µmol", "Circulating bile acid pool")
  output$pkduo    <- plot_ts("CDUO", "mM", "Duodenal bile acid concentration (CMC 1.5 mM)")
  output$pkile    <- plot_ts("ILE1", "µmol", "Terminal ileal lumen (ASBT compartment 1)")
  output$pkserum  <- plot_ts("SERUMBA", "µM", "Systemic serum bile acids")
  output$fgf      <- plot_ts("FGF19", "pg/mL", "Plasma FGF19 (postprandial 90-120 min peak)")
  output$cyp      <- plot_ts("CYP7A1", "Relative", "Hepatic CYP7A1 (ceiling = 6-fold)")
  output$c4plot   <- plot_ts("C4", "ng/mL", "Serum C4 (cut-off 48)")
  output$stoolplot<- plot_ts("WCOL", "mL", "Colonic luminal water")
  output$fatplot  <- plot_ts("FATEFF", "Fraction", "Fat absorption efficiency")
  output$ldlplot  <- plot_ts("LDLC", "mg/dL", "LDL-C")
  output$transit  <- plot_ts("TRANS", "Relative", "Colonic transit rate (positive feedback)")

  output$bmplot <- renderPlot({
    r <- ro()
    d <- data.frame(
      grp = factor(c("Healthy", "Patient (untreated)", "Patient (treated)"),
                   levels = c("Healthy", "Patient (untreated)", "Patient (treated)")),
      bm  = c(r$healthy$BM_per_day, r$untreated$BM_per_day, r$treated$BM_per_day),
      br  = c(r$healthy$bristol, r$untreated$bristol, r$treated$bristol))
    ggplot(d, aes(grp, bm, fill = grp)) + geom_col(width = .6) +
      geom_text(aes(label = sprintf("%.2f per day\nBristol %.1f", bm, br)),
                vjust = -.3, size = 4) +
      scale_fill_manual(values = unname(PAL[c("healthy", "patient", "treated")])) +
      labs(x = NULL, y = "Bowel movements (per day)", title = "Clinical endpoints") +
      expand_limits(y = 0) + guides(fill = "none") + THEME
  })

  # ================= TAB 4: analytic loop ==========================
  output$loopgain <- renderPlot({
    S0 <- 700; Smax <- 6 * S0; h <- 1.5; f0 <- 0.975
    x0 <- f0 / (1 - f0)
    a <- uniroot(function(a) S0 * (1 + (a * x0 * S0)^h) - Smax,
                 c(1e-14, 1))$root
    solveS <- function(f, k)
      uniroot(function(S) S * (1 + (a * k * f / (1 - f) * S)^h) - Smax,
              c(1e-6, Smax))$root
    ff <- seq(0.55, 0.99, length.out = 120)
    d <- bind_rows(lapply(c(1, 0.45, 0.25), function(k)
      data.frame(f = ff, kappa = factor(k),
                 S = sapply(ff, solveS, k = k))))
    ggplot(d, aes(f, S, colour = kappa)) + geom_line(linewidth = 1) +
      geom_hline(yintercept = Smax, linetype = 2) +
      annotate("text", x = 0.6, y = Smax * .96,
               label = "Compensation ceiling Smax = 6 × baseline", size = 4) +
      scale_x_continuous("Single-pass ileal conservation f") +
      scale_y_continuous("Steady-state synthesis = colonic load (µmol/day)") +
      labs(colour = "κ", title = "Closed-form steady-state solution") + THEME
  })

  # ================= TAB 5 =========================================
  output$colspecies <- renderPlot({
    s <- simset()
    d <- last_day(s$treated) %>%
      select(tod, CCA, CCDCA, CDCA, CLCA) %>%
      pivot_longer(-tod, names_to = "Species", values_to = "µmol") %>%
      mutate(Species = recode(Species, CCA = "cholate (w 0.15)",
                         CCDCA = "chenodeoxycholate (w 1.00)",
                         CDCA = "deoxycholate (w 1.15)",
                         CLCA = "lithocholate (w 0.30)"))
    ggplot(d, aes(tod, µmol, fill = Species)) + geom_area(alpha = .85) +
      scale_x_continuous("Last 48 hours (h)", breaks = seq(0, 48, 8)) +
      labs(title = "Colonic bile acid composition — potency differs 8-fold between species", fill = NULL) +
      THEME
  })

  output$threshold <- renderPlot({
    s <- simset()
    d <- bind_rows(
      cbind(grp = "Healthy", last_day(s$healthy)),
      cbind(grp = "Patient (treated)", last_day(s$treated)))
    ggplot(d, aes(tod, DS, colour = grp)) + geom_line(linewidth = 1) +
      geom_hline(yintercept = 3, linetype = 2, colour = "#b3272d") +
      annotate("text", x = 24, y = 3.25,
               label = "Secretory threshold 3 mM (Mekjian 1971)", size = 4,
               colour = "#b3272d") +
      scale_colour_manual(values = unname(PAL[c("healthy", "treated")])) +
      scale_x_continuous("Last 48 hours (h)", breaks = seq(0, 48, 8)) +
      labs(y = "Ds (mM CDCA equivalents)", colour = NULL,
           title = "Colonic secretory driving force") + THEME
  })

  # ================= TAB 7: drug classes ===========================
  drugcmp <- eventReactive(input$go, {
    withProgress(message = "Comparing the drug classes...", value = 0, {
      base <- list(phi = input$phi, kappa = input$kappa)
      sc <- list(
        "Untreated"                        = NULL,
        "A sequestrant (colesevelam)"      = ev_seq(1875, 2, 56),
        "B FXR agonist (OCA 25 mg)"        = ev_oca(25, 56),
        "B FXR agonist (tropifexor)"       = ev_tro(0.090, 56),
        "C ASBT inhibitor (elobixibat)"    = ev_elo(10, 56),
        "D transit (loperamide)"           = ev_lop(2, 56),
        "D transit (ondansetron)"          = ev_ond(4, 56),
        "E microbiota (rifaximin)"         = ev_rif(550, 56),
        "A + B combination"                = ev_seq(1875, 2, 56) + ev_oca(25, 56))
      out <- lapply(names(sc), function(nm) {
        incProgress(1 / length(sc), detail = nm)
        o <- do.call(sim_bam, c(list(days = 56, events = sc[[nm]]), base))
        cbind(Treatment = nm, bam_readout(o))
      })
      bind_rows(out)
    })
  }, ignoreNULL = FALSE)

  output$drugclass <- renderPlot({
    d <- drugcmp()
    base <- d$stool_mL[d$Treatment == "Untreated"]
    d <- d %>% mutate(Change = 100 * (stool_mL - base) / base,
                      C4change = 100 * (C4 - C4[Treatment == "Untreated"]) /
                        C4[Treatment == "Untreated"])
    dd <- d %>% select(Treatment, `Stool water change (%)` = Change,
                       `Serum C4 change (%)` = C4change) %>%
      pivot_longer(-Treatment)
    ggplot(dd, aes(reorder(Treatment, value), value, fill = value > 0)) +
      geom_col() + coord_flip() + facet_wrap(~name, scales = "free_x") +
      geom_hline(yintercept = 0) +
      scale_fill_manual(values = c("TRUE" = "#b3272d", "FALSE" = "#1f6fb2")) +
      labs(x = NULL, y = "Change vs baseline (%)") +
      guides(fill = "none") + THEME
  })

  output$drugtab <- DT::renderDataTable({
    d <- drugcmp()
    DT::datatable(d, options = list(dom = "t", scrollX = TRUE, pageLength = 12),
                  rownames = FALSE) %>% DT::formatRound(2:ncol(d), 2)
  })

  # ================= TAB 8: the 100 cm rule ========================
  rule <- eventReactive(input$go, {
    withProgress(message = "Resection length sweep...", value = 0, {
      Ls <- c(0, 20, 40, 60, 80, 100, 120, 150, 180)
      bind_rows(lapply(Ls, function(L) {
        incProgress(1 / length(Ls), detail = paste0(L, " cm"))
        o <- sim_bam(days = 56, events = NULL, phi = phi_from_resection(L))
        cbind(resection_cm = L, phi = round(phi_from_resection(L), 3), bam_readout(o))
      }))
    })
  }, ignoreNULL = FALSE)

  output$rule100 <- renderPlot({
    d <- rule()
    dd <- d %>% select(resection_cm, `Faecal fat (g/day)` = fecfat_g,
                       `Stool water (mL/day)` = stool_mL,
                       `Bile acid pool (g)` = pool_g,
                       `Peak duodenal [BA] (mM)` = Cduo_peak) %>%
      pivot_longer(-resection_cm)
    hl <- data.frame(name = c("Faecal fat (g/day)", "Peak duodenal [BA] (mM)"),
                     y = c(7, 1.5))
    ggplot(dd, aes(resection_cm, value)) +
      geom_line(linewidth = 1, colour = PAL["patient"]) +
      geom_point(size = 2.2) +
      geom_hline(data = hl, aes(yintercept = y), linetype = 2,
                 colour = "#b3272d") +
      geom_vline(xintercept = 100, linetype = 3, linewidth = .9) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "Length of terminal ileum resected (cm)", y = NULL,
           title = "Dotted line = the 100 cm criterion of the clinical textbooks") + THEME
  })

  output$ruletab <- DT::renderDataTable({
    d <- rule()
    DT::datatable(d, options = list(dom = "t", scrollX = TRUE, pageLength = 10),
                  rownames = FALSE) %>% DT::formatRound(3:ncol(d), 2)
  })

  # ================= TAB 9 =========================================
  baisweep <- eventReactive(input$go, {
    withProgress(message = "Microbiota sweep...", value = 0, {
      bs <- c(0.2, 0.5, 1.0, 1.5, 2.0)
      bind_rows(lapply(bs, function(b) {
        incProgress(1 / length(bs), detail = paste0("BAI ", b))
        o <- sim_bam(days = 40, events = NULL, phi = input$phi,
                     kappa = input$kappa, kbai = 0.16 * b)
        cbind(BAI = b, bam_readout(o))
      }))
    })
  }, ignoreNULL = FALSE)

  output$baiplot <- renderPlot({
    d <- baisweep() %>%
      mutate(`Faecal bile acids (relative)` = fecBA_umold / fecBA_umold[BAI == 1],
             `Stool water (relative)` = stool_mL / stool_mL[BAI == 1]) %>%
      select(BAI, `Faecal bile acids (relative)`, `Stool water (relative)`) %>%
      pivot_longer(-BAI)
    ggplot(d, aes(BAI, value, colour = name)) +
      geom_line(linewidth = 1) + geom_point(size = 2.5) +
      geom_hline(yintercept = 1, linetype = 2) +
      labs(x = "7α-dehydroxylating capacity (BAI, relative)", y = "Fold vs BAI = 1",
           colour = NULL, title = "The load stays the same; only the potency changes") + THEME
  })

  output$biomk <- renderPlot({
    d <- baisweep()
    ggplot(d, aes(Ds, stool_mL)) +
      geom_path(linewidth = 1, colour = PAL["treated"]) +
      geom_point(aes(size = BAI), colour = PAL["patient"]) +
      geom_vline(xintercept = 3, linetype = 2, colour = "#b3272d") +
      labs(x = "Colonic secretory driving force Ds (mM-eq)", y = "Stool water (mL/day)",
           title = "Crossing the threshold is what makes the symptom") + THEME
  })

  output$sehcatcurve <- renderPlot({
    d <- data.frame(k = seq(0.02, 1.2, length.out = 300)) %>%
      mutate(R7 = 100 * exp(-7 * k))
    ggplot(d, aes(k, R7)) + geom_line(linewidth = 1) +
      geom_hline(yintercept = c(5, 10, 15), linetype = 2,
                 colour = c("#b3272d", "#cc8400", "#4f9d69")) +
      annotate("text", x = 1.0, y = c(8, 13, 18),
               label = c("5 % severe", "10 % moderate", "15 % lower limit of normal"),
               size = 4) +
      scale_x_continuous("Pool turnover k = faecal loss / pool (/day)") +
      scale_y_continuous("SeHCAT 7-day retention (%)") +
      labs(title = "SeHCAT is a compressed scale — the whole clinical grading crowds into a narrow band of k") +
      THEME
  })

  # ================= TAB 10 ========================================
  lib <- eventReactive(input$runlib, {
    withProgress(message = "Running the scenario library...", value = 0.1, {
      bam_scenarios()
    })
  })
  output$libtab <- DT::renderDataTable({
    d <- lib()
    DT::datatable(d, options = list(scrollX = TRUE, pageLength = 15),
                  rownames = FALSE) %>% DT::formatRound(2:ncol(d), 2)
  })
}

shinyApp(ui, server)
