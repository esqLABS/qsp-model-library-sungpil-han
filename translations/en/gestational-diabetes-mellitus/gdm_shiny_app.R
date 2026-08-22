# =============================================================================
#  gdm_shiny_app.R
#  Interactive dashboard for the Gestational Diabetes Mellitus QSP model
#  Gestational diabetes mellitus QSP model · Shiny dashboard (9 tabs)
#  ---------------------------------------------------------------------------
#  Run:  shiny::runApp("gdm_shiny_app.R")
#  Needs gdm_mrgsolve_model.R in the same directory (it is sourced for the
#  compiled model object `gdm` and the event/helper functions).
#
#  DESIGN INTENT
#  -------------
#  The app is built around the model's central claim: GDM is a race between the
#  placental contra-insulin clock and the beta-cell compensation clock. So the
#  primary control is NOT "how high is her glucose" — it is BCAP (beta-cell
#  adaptive capacity) and pre-gestational BMI. Glucose is an OUTPUT. Every tab
#  is arranged to make that causal direction visible, and the fetal-exposure tab
#  exists because that is the only axis on which the three drug options
#  genuinely differ.
# =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)

# ---- load the model ---------------------------------------------------------
# Sourcing the model file compiles `gdm` and defines meal_events(), met_events(),
# insulin_events(), glyb_events(), combine_ev(), fasting_series(),
# titrate_insulin(), run_scenario(), at_delivery().
source("gdm_mrgsolve_model.R", local = FALSE, echo = FALSE)

THEME <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.position  = "bottom",
        plot.title       = element_text(face = "bold"))

TARGET_FPG <- 95    # mg/dL, ADA/ACOG fasting target in pregnancy
TARGET_1H  <- 140   # mg/dL, 1-hour post-prandial target

# =============================================================================
#  UI
# =============================================================================
ui <- fluidPage(
  titlePanel("Gestational Diabetes Mellitus (GDM) QSP model"),
  tags$p(style = "color:#555;",
         "Placental contra-insulin drive vs beta-cell compensation — whether it is GDM is not an input but an output.",
         tags$em("Whether GDM occurs is an output of beta-cell adaptive capacity, not an input.")),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      h4("1. Patient"),
      sliderInput("bcap", "Beta-cell adaptive capacity BCAP",
                  min = 0.25, max = 1.20, value = 0.55, step = 0.05),
      helpText(tags$small("1.00 normal · 0.85 secretion-deficient · 0.55 mixed · 0.30 severe")),
      sliderInput("bmi", "Pre-pregnancy BMI (kg/m²)",
                  min = 18, max = 42, value = 31, step = 0.5),
      sliderInput("age", "Maternal age (years) — for risk display",
                  min = 18, max = 48, value = 34, step = 1),
      checkboxInput("lact", "Lactation ≥3 months postpartum", FALSE),

      h4("2. Lifestyle / MNT"),
      sliderInput("cho", "Daily carbohydrate (g/day)",
                  min = 140, max = 280, value = 200, step = 5),
      checkboxInput("lowgi", "Low-GI diet (slowed absorption)", FALSE),
      sliderInput("ex", "Exercise effect EXEFF",
                  min = 0, max = 0.30, value = 0, step = 0.02),
      helpText(tags$small("0.20 ≈ at least 150 min/week of moderate exercise (AMPK pathway, insulin-independent)")),

      h4("3. Pharmacotherapy"),
      radioButtons("drug", NULL,
                   choices = c("None (MNT only)"         = "none",
                               "Metformin"               = "met",
                               "Insulin (basal-bolus)"   = "ins",
                               "Glyburide"               = "glyb",
                               "Metformin + add insulin"  = "both"),
                   selected = "none"),
      sliderInput("tstart", "Treatment start week (GW)",
                  min = 20, max = 36, value = 28, step = 1),
      conditionalPanel("input.drug == 'met' || input.drug == 'both'",
        sliderInput("metdose", "Metformin dose per administration (mg, BID)",
                    min = 250, max = 1250, value = 1000, step = 250),
        sliderInput("metadh", "Adherence (reflecting GI side effects)",
                    min = 0.4, max = 1.0, value = 1.0, step = 0.05)),
      conditionalPanel("input.drug == 'ins' || input.drug == 'both'",
        checkboxInput("titrate", "Auto-titrate to a fasting glucose target of 95 mg/dL", TRUE),
        conditionalPanel("!input.titrate",
          sliderInput("insdose", "Total daily insulin (U/day)",
                      min = 4, max = 140, value = 40, step = 2))),
      conditionalPanel("input.drug == 'glyb'",
        sliderInput("glybdose", "Glyburide dose per administration (mg, BID)",
                    min = 1.25, max = 10, value = 5, step = 1.25)),

      h4("4. Delivery"),
      sliderInput("tdel", "Delivery week (GW)",
                  min = 36, max = 41, value = 39, step = 0.5),

      hr(),
      actionButton("run", "Run simulation", class = "btn-primary",
                   width = "100%"),
      br(), br(),
      checkboxInput("addref", "Also show a normal-pregnancy control", TRUE),
      helpText(tags$small("The first run may take 20-40 s while the model compiles."))
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs", type = "tabs",

        # ------------------------------------------------------------- TAB 1 --
        tabPanel(
          "1. Patient profile",
          br(),
          fluidRow(
            column(6, h4("Diagnostic summary"),
                   tableOutput("tbl_profile")),
            column(6, h4("Position among the Powe physiological subtypes"),
                   plotOutput("plt_subtype", height = "300px"),
                   helpText(tags$small(
                     "Powe (Diabetes Care 2016;39:1052) divided GDM into insulin-resistance-predominant, ",
                     "secretion-deficient and mixed subtypes. In this model the subtype emerges automatically ",
                     "from BCAP and BMI.")))
          ),
          hr(),
          h4("Metabolic adaptation across gestation"),
          plotOutput("plt_adapt", height = "330px")
        ),

        # ------------------------------------------------------------- TAB 2 --
        tabPanel(
          "2. Placental endocrine axis",
          br(),
          h4("Placental growth and the contra-insulin hormones (placental clock)"),
          plotOutput("plt_placenta", height = "420px"),
          helpText("hPL · progesterone · oestradiol · TNF-α follow the placental mass, ",
                   "and the contra-insulin index (CID) that sums them determines insulin sensitivity. ",
                   "This curve is nearly identical in every pregnancy — what differs is the beta-cell side."),
          hr(),
          h4("Adipokines"),
          plotOutput("plt_adipo", height = "300px")
        ),

        # ------------------------------------------------------------- TAB 3 --
        tabPanel(
          "3. Maternal glucose · insulin",
          br(),
          fluidRow(
            column(6, h4("Fasting glucose trajectory"),
                   plotOutput("plt_fpg", height = "310px")),
            column(6, h4("Fasting insulin"),
                   plotOutput("plt_ins", height = "310px"))
          ),
          hr(),
          h4("24-hour glucose profile (selected week) — CGM-like output"),
          sliderInput("cgmweek", "Week displayed (GW)", min = 12, max = 40,
                      value = 34, step = 1, width = "60%"),
          plotOutput("plt_cgm", height = "320px"),
          tableOutput("tbl_cgm")
        ),

        # ------------------------------------------------------------- TAB 4 --
        tabPanel(
          "4. Beta-cell compensation",
          br(),
          h4("Compensation vs failure"),
          plotOutput("plt_beta", height = "400px"),
          helpText("BCM = beta-cell mass x function, KG50E = EC50 of the glucose-sensing curve ",
                   "(leftward shift = increased sensitivity). DI = SI x BCM (disposition index)."),
          hr(),
          fluidRow(
            column(6, h4("Disposition index trajectory"),
                   plotOutput("plt_di", height = "300px")),
            column(6, h4("Glucose-stimulated insulin secretion curve"),
                   plotOutput("plt_isr", height = "300px"))
          )
        ),

        # ------------------------------------------------------------- TAB 5 --
        tabPanel(
          "5. Drug PK · fetal exposure",
          br(),
          h4("Maternal versus fetal drug concentration"),
          plotOutput("plt_pk", height = "360px"),
          tableOutput("tbl_pk"),
          hr(),
          tags$div(
            style = "background:#fff6f6; border-left:4px solid #c04040; padding:10px;",
            tags$b("Why this tab exists:"), br(),
            "The three drugs are similar to one another in their maternal glucose-lowering effect. What actually differs is ",
            tags$b("fetal exposure"), " — insulin does not cross the placenta (structurally zero), ",
            "metformin crosses freely to give umbilical:maternal ≈ 1.0, while glyburide ",
            "reaches cord:maternal ≈ 0.7 despite BCRP/ABCG2 efflux."
          ),
          br(),
          h4("Insulin requirement (the result of titration)"),
          verbatimTextOutput("txt_titrate")
        ),

        # ------------------------------------------------------------- TAB 6 --
        tabPanel(
          "6. Fetal growth",
          br(),
          h4("Fetal weight and composition"),
          plotOutput("plt_fetal", height = "360px"),
          helpText("Grey dashed line = reference (50th percentile) curve. The overgrowth is ", tags$b("asymmetric"),
                   " — fat is far more insulin-elastic than lean mass."),
          hr(),
          fluidRow(
            column(6, h4("Fat vs lean"),
                   plotOutput("plt_comp", height = "300px")),
            column(6, h4("Fetal glucose · insulin (the Pedersen chain)"),
                   plotOutput("plt_fetal_gi", height = "300px"))
          )
        ),

        # ------------------------------------------------------------- TAB 7 --
        tabPanel(
          "7. Clinical endpoints",
          br(),
          h4("Endpoints at delivery"),
          tableOutput("tbl_endpoint"),
          plotOutput("plt_endpoint", height = "330px"),
          hr(),
          h4("Against the calibration targets"),
          tableOutput("tbl_calib"),
          helpText("Calibration basis: HAPO (NEJM 2008;358:1991) continuous glucose–LGA gradient, ",
                   "the Landon/MFMU treatment effect in mild GDM (NEJM 2009;361:1339), ",
                   "ACHOIS (NEJM 2005;352:2477).")
        ),

        # ------------------------------------------------------------- TAB 8 --
        tabPanel(
          "8. Scenario comparison",
          br(),
          h4("Treatment strategy comparison (the same patient physiology, a different strategy)"),
          checkboxGroupInput(
            "cmp", NULL, inline = TRUE,
            choices = c("No treatment" = "none", "MNT+exercise" = "mnt",
                        "Metformin" = "met", "Insulin" = "ins",
                        "Glyburide" = "glyb"),
            selected = c("none", "mnt", "ins")),
          actionButton("runcmp", "Run comparison", class = "btn-default"),
          br(), br(),
          plotOutput("plt_cmp", height = "330px"),
          tableOutput("tbl_cmp")
        ),

        # ------------------------------------------------------------- TAB 9 --
        tabPanel(
          "9. Biomarkers · postpartum",
          br(),
          fluidRow(
            column(6, h4("Glycation indices (HbA1c, mean glucose)"),
                   plotOutput("plt_a1c", height = "300px")),
            column(6, h4("Free fatty acids · maternal fat mass"),
                   plotOutput("plt_bio", height = "300px"))
          ),
          hr(),
          h4("Postpartum trajectory and 5-year type 2 diabetes risk"),
          plotOutput("plt_pp", height = "330px"),
          tableOutput("tbl_pp"),
          helpText("When the placenta involutes at delivery the CID collapses and insulin sensitivity recovers, ",
                   "but the beta-cell deficit does not — the postpartum disposition index ",
                   "reveals that deficit (Bellamy Lancet 2009;373:1773, RR 7.43).")
        ),

        # ---------------------------------------------------------- TAB 10 --
        tabPanel(
          "10. HAPO validation",
          br(),
          h4("Reproducing the HAPO continuous gradient (no threshold)"),
          actionButton("runhapo", "Run sweep (BCAP x BMI grid — about 1-3 min)",
                       class = "btn-default"),
          br(), br(),
          plotOutput("plt_hapo", height = "380px"),
          tableOutput("tbl_hapo"),
          helpText("Observed: HAPO category 1 (FPG <75 mg/dL) LGA 5.3% → category 7 (≥100) 26.3%. ",
                   "What the model has to reproduce is not a particular cut-point but ", tags$b("continuity without a threshold"),
                   ".")
        )
      )
    )
  )
)

# =============================================================================
#  SERVER
# =============================================================================
server <- function(input, output, session) {

  # ---- assemble parameters and events from the UI ---------------------------
  build <- function(drug = input$drug, bcap = input$bcap, bmi = input$bmi,
                    ex = input$ex, cho = input$cho, lowgi = input$lowgi) {
    tdel <- input$tdel*7
    pars <- list(BCAP = bcap, BMI = bmi, EXEFF = ex, CHOD = cho,
                 TDEL = tdel, LACT = as.numeric(input$lact),
                 KTR_G = if (lowgi) 16 else 26)
    end  <- max(tdel + 84, 364)
    ev   <- meal_events(cho, 56, end)
    titr <- NA_real_

    if (drug %in% c("met", "both")) {
      pars$METADH <- input$metadh
      ev <- combine_ev(ev, met_events(input$metdose, input$tstart, input$tdel))
    }
    if (drug %in% c("ins", "both")) {
      if (isTRUE(input$titrate)) {
        tt <- titrate_insulin(pars, TARGET_FPG, input$tstart, input$tdel, cho,
                              u_start = if (drug == "both") 12 else 20)
        titr <- tt$total_u
      } else {
        titr <- input$insdose
      }
      ev <- combine_ev(ev, insulin_events(titr, input$tstart, input$tdel))
    }
    if (drug == "glyb") {
      ev <- combine_ev(ev, glyb_events(input$glybdose, TRUE,
                                       input$tstart, input$tdel))
    }
    list(pars = pars, ev = ev, end = end, tdel = tdel, titr = titr)
  }

  simulate <- function(spec, delta = 0.05) {
    gdm |> param(spec$pars) |> data_set(spec$ev) |>
      mrgsim(start = 56, end = spec$end, delta = delta, hmax = 0.05) |>
      as.data.frame()
  }

  # ---- main reactive --------------------------------------------------------
  sim <- eventReactive(input$run, {
    withProgress(message = "Simulating...", value = 0.3, {
      spec <- build()
      incProgress(0.4)
      d <- simulate(spec)
      d$arm <- "Patient"
      out <- list(spec = spec, main = d, titr = spec$titr)

      if (isTRUE(input$addref)) {
        rspec <- build(drug = "none", bcap = 1.00, bmi = 22, ex = 0,
                       cho = input$cho, lowgi = FALSE)
        r <- simulate(rspec)
        r$arm <- "Normal control (NGT)"
        out$ref <- r
      }
      incProgress(0.3)
      out
    })
  }, ignoreNULL = FALSE)

  allsim <- reactive({
    s <- sim()
    if (is.null(s$ref)) s$main else rbind(s$main, s$ref)
  })

  fast <- reactive({
    allsim() |>
      mutate(day = floor(time), frac = time - day) |>
      filter(frac > 0.20, frac < 0.29) |>
      group_by(arm, day) |>
      summarise(across(c(GLU, INS, MBG, HBA1C, FFA, BCM, SIP, DI, DIREF,
                         CID, KG50E, FATM, ADIPO, TNFA, HPL, PROG, E2,
                         PLAC, GLUF, INSF, BCF, FWT, FWREF, BWZ, FATPCT,
                         CMET, CMETF, CGLY, CGLYF, NHR, PLGA, PNHYP,
                         PSD, PCS, PPE, CIT2D, HYPOAUC),
                       mean), .groups = "drop") |>
      mutate(GW = day/7)
  })

  del <- reactive({
    s <- sim(); td <- s$spec$tdel
    allsim() |> group_by(arm) |> slice(which.min(abs(time - td))) |> ungroup()
  })

  # ========================================================== TAB 1 ==========
  output$tbl_profile <- renderTable({
    d <- del(); f <- fast()
    ft <- f |> filter(GW <= input$tdel) |> group_by(arm) |>
      summarise(FPG28 = mean(GLU[GW >= 27 & GW <= 29]),
                FPGterm = mean(tail(GLU, 7)),
                INSterm = mean(tail(INS, 7)), .groups = "drop")
    data.frame(
      Item = c("BCAP (beta-cell adaptive capacity)", "Pre-pregnancy BMI", "Maternal age",
               "GW 28 fasting glucose (mg/dL)", "Fasting glucose at delivery (mg/dL)",
               "Fasting insulin at delivery (µU/mL)", "Insulin sensitivity (% of pre-pregnancy)",
               "Beta-cell mass x function BCM", "Disposition index (relative)",
               "HbA1c (%)", "Above the IADPSG fasting threshold (5.1 mmol/L = 92)"),
      Value = c(sprintf("%.2f", input$bcap), sprintf("%.1f", input$bmi),
             sprintf("%d", input$age),
             sprintf("%.1f", ft$FPG28[ft$arm == "Patient"]),
             sprintf("%.1f", ft$FPGterm[ft$arm == "Patient"]),
             sprintf("%.1f", ft$INSterm[ft$arm == "Patient"]),
             sprintf("%.0f%%", 100*d$MATSUDA[d$arm == "Patient"]),
             sprintf("%.2f", d$BCM[d$arm == "Patient"]),
             sprintf("%.2f", (d$DI/d$DIREF)[d$arm == "Patient"]),
             sprintf("%.2f", d$HBA1C[d$arm == "Patient"]),
             ifelse(ft$FPG28[ft$arm == "Patient"] > 92, "Yes", "No")),
      check.names = FALSE)
  }, striped = TRUE, spacing = "xs")

  output$plt_subtype <- renderPlot({
    d <- del() |> filter(arm == "Patient")
    grid <- expand.grid(BCAP = seq(0.3, 1.15, 0.05), BMI = seq(20, 40, 1))
    grid$SI  <- 1/(1 + 1.40)*exp(-0.055*(grid$BMI - 22))
    grid$SEC <- grid$BCAP
    ggplot(grid, aes(SEC, SI)) +
      geom_raster(aes(fill = SEC*SI), alpha = 0.55) +
      scale_fill_gradient(low = "#c94040", high = "#3a7d44", guide = "none") +
      annotate("point", x = input$bcap,
               y = 1/(1 + 1.40)*exp(-0.055*(input$bmi - 22)),
               size = 5, shape = 21, fill = "white", stroke = 1.4) +
      annotate("text", x = 0.42, y = 0.40, label = "secretion-deficient",
               size = 3.4) +
      annotate("text", x = 1.05, y = 0.13, label = "insulin-resistant",
               size = 3.4) +
      annotate("text", x = 0.45, y = 0.13, label = "mixed", size = 3.4) +
      labs(x = "Secretory capacity (BCAP)", y = "Late-gestation insulin sensitivity (relative)",
           title = "Physiological subtype coordinates") + THEME
  })

  output$plt_adapt <- renderPlot({
    fast() |>
      select(arm, GW, `Insulin sensitivity SI` = SIP, `Beta-cell BCM` = BCM,
             `Contra-insulin CID` = CID, `Fasting glucose FPG` = GLU) |>
      pivot_longer(-c(arm, GW)) |>
      ggplot(aes(GW, value, colour = arm)) +
      geom_line(linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "Gestational week (GW)", y = NULL, colour = NULL,
           title = "Two clocks: placental drive vs beta-cell compensation") + THEME
  })

  # ========================================================== TAB 2 ==========
  output$plt_placenta <- renderPlot({
    fast() |>
      transmute(arm, GW,
                `Placental mass (norm.)` = PLAC,
                `hPL (mg/L)` = HPL,
                `Progesterone (ng/mL)` = PROG,
                `Oestradiol (ng/mL)` = E2/1000,
                `TNF-α (pg/mL)` = TNFA,
                `Contra-insulin index CID` = CID) |>
      pivot_longer(-c(arm, GW)) |>
      ggplot(aes(GW, value, colour = arm)) +
      geom_line(linewidth = 0.8) +
      geom_vline(xintercept = input$tdel, linetype = 3) +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "Gestational week (GW)", y = NULL, colour = NULL,
           title = "The placental clock — nearly identical in every pregnancy") + THEME
  })

  output$plt_adipo <- renderPlot({
    fast() |>
      transmute(arm, GW, `Adiponectin (µg/mL)` = ADIPO,
                `Maternal fat mass (kg)` = FATM) |>
      pivot_longer(-c(arm, GW)) |>
      ggplot(aes(GW, value, colour = arm)) +
      geom_line(linewidth = 0.8) + facet_wrap(~name, scales = "free_y") +
      labs(x = "Gestational week (GW)", y = NULL, colour = NULL) + THEME
  })

  # ========================================================== TAB 3 ==========
  output$plt_fpg <- renderPlot({
    ggplot(fast(), aes(GW, GLU, colour = arm)) +
      geom_line(linewidth = 0.9) +
      geom_hline(yintercept = TARGET_FPG, linetype = 2, colour = "#b03030") +
      geom_hline(yintercept = 92, linetype = 3, colour = "#3060b0") +
      geom_vline(xintercept = input$tdel, linetype = 3) +
      annotate("text", x = 14, y = TARGET_FPG + 3, size = 3,
               label = "Treatment target 95", colour = "#b03030") +
      annotate("text", x = 14, y = 89, size = 3,
               label = "IADPSG 92", colour = "#3060b0") +
      labs(x = "Gestational week (GW)", y = "Fasting glucose (mg/dL)", colour = NULL) + THEME
  })

  output$plt_ins <- renderPlot({
    ggplot(fast(), aes(GW, INS, colour = arm)) +
      geom_line(linewidth = 0.9) +
      geom_vline(xintercept = input$tdel, linetype = 3) +
      labs(x = "Gestational week (GW)", y = "Fasting insulin (µU/mL)", colour = NULL) + THEME
  })

  output$plt_cgm <- renderPlot({
    d0 <- floor(input$cgmweek*7)
    allsim() |> filter(time >= d0, time < d0 + 2) |>
      mutate(hour = (time - d0)*24) |>
      ggplot(aes(hour, GLU, colour = arm)) +
      geom_line(linewidth = 0.8) +
      geom_hline(yintercept = c(63, TARGET_1H), linetype = 2, colour = "grey40") +
      scale_x_continuous(breaks = seq(0, 48, 6)) +
      labs(x = sprintf("Hours after GW %d (h)", input$cgmweek),
           y = "Glucose (mg/dL)", colour = NULL,
           title = "48-hour glucose profile (CGM-like)",
           subtitle = "Dashed lines = 63 and 140 mg/dL (TIR range)") + THEME
  })

  output$tbl_cgm <- renderTable({
    d0 <- floor(input$cgmweek*7)
    allsim() |> filter(time >= d0, time < d0 + 7) |>
      group_by(arm) |>
      summarise(`Mean glucose (mg/dL)` = round(mean(GLU), 1),
                `Peak glucose` = round(max(GLU), 1),
                `Nadir glucose` = round(min(GLU), 1),
                `TIR 63-140 (%)` = round(100*mean(GLU >= 63 & GLU <= 140), 1),
                `Time >140 (%)` = round(100*mean(GLU > 140), 1),
                `Time <70 (%)` = round(100*mean(GLU < 70), 1),
                .groups = "drop")
  }, striped = TRUE, spacing = "xs")

  # ========================================================== TAB 4 ==========
  output$plt_beta <- renderPlot({
    fast() |>
      transmute(arm, GW, `Beta-cell BCM` = BCM,
                `Glucose EC50 KG50E (mg/dL)` = KG50E,
                `Disposition index (rel.)` = DI/DIREF,
                `Fasting glucose (mg/dL)` = GLU) |>
      pivot_longer(-c(arm, GW)) |>
      ggplot(aes(GW, value, colour = arm)) +
      geom_line(linewidth = 0.8) + facet_wrap(~name, scales = "free_y") +
      geom_vline(xintercept = input$tdel, linetype = 3) +
      labs(x = "Gestational week (GW)", y = NULL, colour = NULL,
           title = "Did the compensation succeed, or fail") + THEME
  })

  output$plt_di <- renderPlot({
    ggplot(fast(), aes(SIP, BCM, colour = arm)) +
      geom_path(linewidth = 0.9, arrow = arrow(length = unit(0.15, "cm"))) +
      labs(x = "Insulin sensitivity SI", y = "Beta-cell BCM", colour = NULL,
           title = "The trajectory on the hyperbolic relationship",
           subtitle = "In the normal case BCM rises enough as SI falls that the point travels along the curve") +
      THEME
  })

  output$plt_isr <- renderPlot({
    d <- del() |> filter(arm == "Patient")
    g <- seq(50, 250, 2)
    curves <- rbind(
      data.frame(GLU = g, arm = "At delivery (adapted)",
                 ISR = d$BCM*g^2/(d$KG50E^2 + g^2)),
      data.frame(GLU = g, arm = "Pre-pregnancy (non-adapted)",
                 ISR = 1.0*g^2/(130^2 + g^2)))
    ggplot(curves, aes(GLU, ISR, colour = arm)) +
      geom_line(linewidth = 0.9) +
      geom_vline(xintercept = d$GLU, linetype = 2) +
      labs(x = "Glucose (mg/dL)", y = "Relative insulin secretion rate", colour = NULL,
           title = "Glucose-stimulated insulin secretion curve",
           subtitle = "Pregnancy moves the curve up (mass) and to the left (sensitivity)") +
      THEME
  })

  # ========================================================== TAB 5 ==========
  output$plt_pk <- renderPlot({
    d <- allsim() |> filter(arm == "Patient",
                            time > input$tstart*7 - 1, time < input$tdel*7 + 1)
    long <- rbind(
      data.frame(GW = d$GWk, conc = d$CMET,  who = "Maternal", drug = "Metformin (mg/L)"),
      data.frame(GW = d$GWk, conc = d$CMETF, who = "Fetal",    drug = "Metformin (mg/L)"),
      data.frame(GW = d$GWk, conc = d$CGLY,  who = "Maternal", drug = "Glyburide (mg/L)"),
      data.frame(GW = d$GWk, conc = d$CGLYF, who = "Fetal",    drug = "Glyburide (mg/L)"))
    long <- long |> group_by(drug) |> filter(max(conc) > 1e-6) |> ungroup()
    if (nrow(long) == 0) {
      return(ggplot() + annotate("text", 0, 0, size = 5,
        label = "No oral drug has been selected.\nInsulin does not cross the placenta, so fetal exposure is structurally zero.") +
        theme_void())
    }
    ggplot(long, aes(GW, conc, colour = who)) +
      geom_line(linewidth = 0.7) + facet_wrap(~drug, scales = "free_y") +
      labs(x = "Gestational week (GW)", y = "Blood concentration", colour = NULL,
           title = "The fetus is exposed to whatever the mother is exposed to") + THEME
  })

  output$tbl_pk <- renderTable({
    d <- allsim() |> filter(arm == "Patient",
                            time > input$tdel*7 - 7, time <= input$tdel*7)
    data.frame(
      Drug = c("Metformin", "Glyburide", "Insulin"),
      `Maternal Cmax` = c(round(max(d$CMET), 3), round(max(d$CGLY), 4), NA),
      `Maternal Cavg` = c(round(mean(d$CMET), 3), round(mean(d$CGLY), 4), NA),
      `Fetal Cavg` = c(round(mean(d$CMETF), 3), round(mean(d$CGLYF), 4), 0),
      `Fetal:maternal ratio` = c(
        ifelse(mean(d$CMET) > 1e-6, round(mean(d$CMETF)/mean(d$CMET), 2), NA),
        ifelse(mean(d$CGLY) > 1e-8, round(mean(d$CGLYF)/mean(d$CGLY), 2), NA),
        0),
      `Reported ratio (literature)` = c("~1.0 (Vanky 2005)", "~0.7 (Hemauer 2010)", "0 (does not cross)"),
      check.names = FALSE)
  }, striped = TRUE, spacing = "xs")

  output$txt_titrate <- renderPrint({
    s <- sim()
    if (is.na(s$titr)) {
      cat("This is a scenario that does not use insulin.\n")
    } else {
      wt <- input$bmi*1.62^2
      cat(sprintf("Total daily insulin: %.0f U/day  (%.2f U/kg/day, assuming a body weight of %.0f kg)\n",
                  s$titr, s$titr/wt, wt))
      cat(sprintf("basal:bolus = 50:50, the bolus split across 3 meals\n"))
      cat("Note: the third-trimester requirement usually rises to 0.9-1.0 U/kg/day.\n")
      cat("This value was not prescribed; it is the result of titration to reach the\n")
      cat("fasting glucose target of 95 mg/dL — that is, a product of the patient's physiology.\n")
    }
  })

  # ========================================================== TAB 6 ==========
  output$plt_fetal <- renderPlot({
    f <- fast() |> filter(GW <= input$tdel)
    ggplot(f, aes(GW, FWT, colour = arm)) +
      geom_line(linewidth = 0.9) +
      geom_line(aes(y = FWREF), colour = "grey45", linetype = 2) +
      labs(x = "Gestational week (GW)", y = "Fetal weight (g)", colour = NULL,
           title = "Fetal growth trajectory",
           subtitle = "Grey dashed line = reference 50th percentile") + THEME
  })

  output$plt_comp <- renderPlot({
    allsim() |> filter(time <= input$tdel*7) |>
      transmute(arm, GW = GWk, `Fat (g)` = FATF, `Lean (g)` = LEANF) |>
      pivot_longer(-c(arm, GW)) |>
      ggplot(aes(GW, value, colour = arm, linetype = name)) +
      geom_line(linewidth = 0.8) +
      labs(x = "Gestational week (GW)", y = "Mass (g)", colour = NULL, linetype = NULL,
           title = "Asymmetric overgrowth: fat is the elastic compartment") + THEME
  })

  output$plt_fetal_gi <- renderPlot({
    fast() |> filter(GW <= input$tdel) |>
      transmute(arm, GW, `Fetal glucose (mg/dL)` = GLUF,
                `Fetal insulin (µU/mL)` = INSF,
                `Fetal beta-cell BCF` = BCF,
                `Fetal fat fraction (%)` = FATPCT) |>
      pivot_longer(-c(arm, GW)) |>
      ggplot(aes(GW, value, colour = arm)) +
      geom_line(linewidth = 0.8) + facet_wrap(~name, scales = "free_y") +
      labs(x = "Gestational week (GW)", y = NULL, colour = NULL,
           title = "The Pedersen chain") + THEME
  })

  # ========================================================== TAB 7 ==========
  ep_tbl <- reactive({
    del() |> transmute(
      arm,
      `Birth weight (g)` = round(FWT),
      `Weight z-score` = round(BWZ, 2),
      `Fetal fat (%)` = round(FATPCT, 1),
      `Cord insulin ratio (NHR)` = round(NHR, 2),
      `LGA (%)` = round(100*PLGA, 1),
      `Macrosomia >4000g (%)` = round(100*PMACR, 1),
      `Neonatal hypoglycaemia (%)` = round(100*PNHYP, 1),
      `Shoulder dystocia (%)` = round(100*PSD, 1),
      `Pre-eclampsia (%)` = round(100*PPE, 1),
      `Caesarean section (%)` = round(100*PCS, 1),
      `Maternal hypoglycaemia exposure (mg/dL·d)` = round(HYPOAUC, 0),
      `5-year T2DM (%)` = round(100*CIT2D, 1))
  })

  output$tbl_endpoint <- renderTable(ep_tbl(), striped = TRUE, spacing = "xs")

  output$plt_endpoint <- renderPlot({
    del() |>
      transmute(arm, LGA = 100*PLGA, `Macrosomia` = 100*PMACR,
                `Neonatal hypoglycaemia` = 100*PNHYP, `Shoulder dystocia` = 100*PSD,
                `Pre-eclampsia` = 100*PPE, `Caesarean section` = 100*PCS) |>
      pivot_longer(-arm) |>
      ggplot(aes(reorder(name, -value), value, fill = arm)) +
      geom_col(position = "dodge") +
      labs(x = NULL, y = "Probability (%)", fill = NULL,
           title = "Clinical endpoints at delivery") + THEME
  })

  output$tbl_calib <- renderTable({
    d <- del() |> filter(arm == "Patient")
    data.frame(
      Endpoint = c("LGA", "Macrosomia >4000 g", "Shoulder dystocia", "Pre-eclampsia/HDP", "Caesarean section"),
      `Model (%)` = c(round(100*d$PLGA, 1), round(100*d$PMACR, 1),
                     round(100*d$PSD, 1), round(100*d$PPE, 1),
                     round(100*d$PCS, 1)),
      `MFMU untreated (%)` = c(14.5, 14.3, 4.0, 13.6, 33.8),
      `MFMU treated (%)`   = c(7.1, 5.9, 1.5, 8.6, 26.9),
      Source = rep("Landon NEJM 2009;361:1339", 5),
      check.names = FALSE)
  }, striped = TRUE, spacing = "xs")

  # ========================================================== TAB 8 ==========
  cmp <- eventReactive(input$runcmp, {
    withProgress(message = "Comparing strategies...", value = 0.1, {
      arms <- input$cmp
      out <- list()
      n <- max(1, length(arms))
      for (a in arms) {
        spec <- switch(
          a,
          none = build(drug = "none", ex = 0, cho = input$cho, lowgi = FALSE),
          mnt  = build(drug = "none", ex = 0.20, cho = 175, lowgi = TRUE),
          met  = build(drug = "met",  ex = 0.20, cho = 175, lowgi = TRUE),
          ins  = build(drug = "ins",  ex = 0.20, cho = 175, lowgi = TRUE),
          glyb = build(drug = "glyb", ex = 0.20, cho = 175, lowgi = TRUE))
        d <- simulate(spec, delta = 0.05)
        d$arm <- c(none = "No treatment", mnt = "MNT+exercise", met = "Metformin",
                   ins = "Insulin", glyb = "Glyburide")[[a]]
        d$titr <- spec$titr
        out[[a]] <- d
        incProgress(1/n)
      }
      do.call(rbind, out)
    })
  })

  output$plt_cmp <- renderPlot({
    cmp() |>
      mutate(day = floor(time), frac = time - day) |>
      filter(frac > 0.20, frac < 0.29) |>
      group_by(arm, day) |> summarise(FPG = mean(GLU), .groups = "drop") |>
      ggplot(aes(day/7, FPG, colour = arm)) +
      geom_line(linewidth = 0.9) +
      geom_hline(yintercept = TARGET_FPG, linetype = 2) +
      labs(x = "Gestational week (GW)", y = "Fasting glucose (mg/dL)", colour = NULL,
           title = "The same patient, a different strategy") + THEME
  })

  output$tbl_cmp <- renderTable({
    td <- input$tdel*7
    cmp() |> group_by(arm) |> slice(which.min(abs(time - td))) |>
      transmute(Strategy = arm,
                `Fasting glucose` = round(GLU, 1),
                `Mean glucose` = round(MBG, 1),
                `HbA1c` = round(HBA1C, 2),
                `Birth weight (g)` = round(FWT),
                `z-score` = round(BWZ, 2),
                `Fetal fat %` = round(FATPCT, 1),
                `LGA %` = round(100*PLGA, 1),
                `Neonatal hypoglycaemia %` = round(100*PNHYP, 1),
                `Shoulder dystocia %` = round(100*PSD, 1),
                `Fetal metformin` = round(CMETF, 3),
                `Fetal glyburide` = round(CGLYF, 4),
                `5-year T2DM %` = round(100*CIT2D, 1)) |>
      ungroup()
  }, striped = TRUE, spacing = "xs")

  # ========================================================== TAB 9 ==========
  output$plt_a1c <- renderPlot({
    fast() |> transmute(arm, GW, `Mean glucose MBG (mg/dL)` = MBG,
                        `HbA1c (%)` = HBA1C) |>
      pivot_longer(-c(arm, GW)) |>
      ggplot(aes(GW, value, colour = arm)) +
      geom_line(linewidth = 0.8) + facet_wrap(~name, scales = "free_y") +
      geom_vline(xintercept = input$tdel, linetype = 3) +
      labs(x = "Gestational week (GW)", y = NULL, colour = NULL) + THEME
  })

  output$plt_bio <- renderPlot({
    fast() |> transmute(arm, GW, `Free fatty acids (mmol/L)` = FFA,
                        `Maternal fat mass (kg)` = FATM) |>
      pivot_longer(-c(arm, GW)) |>
      ggplot(aes(GW, value, colour = arm)) +
      geom_line(linewidth = 0.8) + facet_wrap(~name, scales = "free_y") +
      labs(x = "Gestational week (GW)", y = NULL, colour = NULL) + THEME
  })

  output$plt_pp <- renderPlot({
    fast() |>
      transmute(arm, GW, `Insulin sensitivity SI` = SIP, `Beta-cell BCM` = BCM,
                `Disposition index (rel.)` = DI/DIREF,
                `Fasting glucose (mg/dL)` = GLU) |>
      pivot_longer(-c(arm, GW)) |>
      ggplot(aes(GW, value, colour = arm)) +
      geom_line(linewidth = 0.8) +
      geom_vline(xintercept = input$tdel, linetype = 2, colour = "#b03030") +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "Gestational week (GW) — after delivery these are postpartum weeks", y = NULL, colour = NULL,
           title = "What remains once the placenta is gone at delivery (red line)",
           subtitle = "SI recovers but the beta-cell deficit does not") + THEME
  })

  output$tbl_pp <- renderTable({
    td <- input$tdel*7
    allsim() |> group_by(arm) |>
      slice(which.min(abs(time - (td + 56)))) |>   # ~8 weeks post partum
      transmute(arm,
                `Postpartum fasting glucose (mg/dL)` = round(GLU, 1),
                `Postpartum insulin (µU/mL)` = round(INS, 1),
                `Postpartum SI (% of pre-pregnancy)` = round(100*MATSUDA, 0),
                `Postpartum BCM` = round(BCM, 2),
                `Postpartum DI (relative)` = round(DI/DIREF, 2),
                `5-year T2DM cumulative incidence (%)` = round(100*CIT2D, 1),
                `Lactation included` = ifelse(input$lact, "Yes", "No")) |>
      ungroup()
  }, striped = TRUE, spacing = "xs")

  # ========================================================= TAB 10 ==========
  hapo <- eventReactive(input$runhapo, {
    withProgress(message = "Running the HAPO sweep (BCAP x BMI)...", value = 0.1, {
      grid <- expand.grid(BCAP = seq(0.45, 1.15, 0.10), BMI = c(22, 27, 32, 37))
      res <- lapply(seq_len(nrow(grid)), function(i) {
        p <- list(BCAP = grid$BCAP[i], BMI = grid$BMI[i], TDEL = 273)
        s <- gdm |> param(p) |> data_set(meal_events(200, 56, 280)) |>
          mrgsim(start = 56, end = 280, delta = 0.25, hmax = 0.05) |>
          as.data.frame()
        f <- fasting_series(s, "GLU")
        d <- s[which.min(abs(s$time - 273)), ]
        incProgress(1/nrow(grid))
        data.frame(BCAP = grid$BCAP[i], BMI = grid$BMI[i],
                   FPG = mean(tail(f$GLU, 7)), MBG = d$MBG,
                   BW = d$FWT, BWZ = d$BWZ, FATPCT = d$FATPCT,
                   LGA = 100*d$PLGA, NEOHYPO = 100*d$PNHYP,
                   SD = 100*d$PSD, T2DM = 100*d$CIT2D)
      })
      do.call(rbind, res)
    })
  })

  output$plt_hapo <- renderPlot({
    h <- hapo()
    obs <- data.frame(FPG = c(72, 77, 81, 85, 89, 94, 102),
                      LGA = c(5.3, 7.9, 9.4, 11.5, 14.4, 19.3, 26.3))
    ggplot(h, aes(FPG, LGA)) +
      geom_point(aes(colour = factor(BMI)), size = 2.4) +
      geom_smooth(se = FALSE, method = "loess", formula = y ~ x,
                  colour = "grey30", linewidth = 0.8) +
      geom_point(data = obs, shape = 4, size = 3.6, stroke = 1.2,
                 colour = "#b03030") +
      geom_line(data = obs, colour = "#b03030", linetype = 2) +
      labs(x = "Late-gestation fasting glucose (mg/dL)", y = "LGA probability (%)",
           colour = "Pre-pregnancy BMI",
           title = "The HAPO continuous gradient — there is no threshold",
           subtitle = "Red X / dashed line = HAPO observed (categories 1-7); points = model") + THEME
  })

  output$tbl_hapo <- renderTable({
    hapo() |> arrange(FPG) |>
      transmute(BCAP, BMI, `Fasting glucose` = round(FPG, 1),
                `Mean glucose` = round(MBG, 1), `Birth weight` = round(BW),
                `z-score` = round(BWZ, 2), `Fetal fat %` = round(FATPCT, 1),
                `LGA %` = round(LGA, 1), `Neonatal hypoglycaemia %` = round(NEOHYPO, 1),
                `Shoulder dystocia %` = round(SD, 1), `5-year T2DM %` = round(T2DM, 1))
  }, striped = TRUE, spacing = "xs")
}

shinyApp(ui, server)

# =============================================================================
#  NOTE ON HONEST PRESENTATION
#  ---------------------------------------------------------------------------
#  Tab 7 shows the model's endpoint probabilities SIDE BY SIDE with the observed
#  MFMU trial rates rather than in isolation, and tab 10 overlays the actual
#  HAPO category means on the model sweep. That is deliberate: a QSP dashboard
#  that only shows its own output invites the reader to mistake calibration for
#  validation. Where the model disagrees with the trial numbers, the disagreement
#  should be visible on the same axes.
# =============================================================================
