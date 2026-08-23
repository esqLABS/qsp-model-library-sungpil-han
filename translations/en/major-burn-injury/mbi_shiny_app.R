## =====================================================================
##  mbi_shiny_app.R
##  MAJOR THERMAL BURN INJURY — QSP interactive dashboard
##  =====================================================================
##
##  The app is organised around the model's three statements rather than
##  around organ systems, because the statements are the model:
##
##   (1) THE PRODUCT      dVP = V_infused x f_ret(Pi_p), and Pi_p is
##                        CONVEX in protein, so f_ret decays as you work.
##   (2) THE CONTROLLER   UO-titrated resuscitation is a closed loop whose
##                        gain IS f_ret.  Fluid creep is its fixed point.
##   (3) THE TWO CLOCKS   leak tau ~ 9.5 h vs hypermetabolism tau ~ 104 d,
##                        and the slow clock is driven by A_open(t).
##
##  Tabs
##   1  Patient & protocol      — set the case and the resuscitation plan
##   2  The convex function     — COP vs protein, and where you are on it
##   3  Fluid & Starling        — compartments, fluxes, retention fraction
##   4  The controller          — rate, urine output, loop gain, creep
##   5  Resuscitation morbidity — IAP, ACS, lung water, weight gain
##   6  Hypermetabolism         — REE, catecholamines, core temperature
##   7  Body composition        — lean mass, the S-B difference, bone, fat
##   8  Wound & closure         — A_open, grafts, donor sites, scar
##   9  Inflammation & infection— cytokines, HLA-DR, colony count, sepsis
##  10  Pharmacology & PK       — propranolol, oxandrolone, vancomycin/ARC
##  11  Endpoints & validation  — the calibration table and rBaux
##  12  Scenario comparison     — any subset of the 22 scenarios overlaid
##
##  Run:
##      source("mbi_mrgsolve_model.R")   # defines mod, SCEN, run_scenario
##      shiny::runApp("mbi_shiny_app.R")
## =====================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)

if (!exists("mod") || !exists("SCEN")) {
  source("mbi_mrgsolve_model.R")
}

THEME <- theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 13),
        legend.position = "bottom", legend.title = element_blank())

LP <- function(C) 2.1*C + 0.16*C^2 + 0.009*C^3

TARGETS <- tibble::tribble(
  ~metric,                                   ~lit,                     ~lo,   ~hi,
  "24-h volume (mL/kg/%TBSA)",               "5.2-6.7",                5.2,   6.7,
  "in:out ratio vs Parkland",                "1.2-1.6",                1.2,   1.6,
  "plasma volume nadir (%)",                 "60-80",                 60,    80,
  "plasma COP at 24 h (mmHg)",               "10-16",                 10,    16,
  "peak weight gain (%)",                    "+15 to +30",            15,    30,
  "peak REE (% predicted)",                  "120-180",              120,   180,
  "lean mass at 14 d, no drug (%)",          "about -9",             -13,    -5,
  "24-h volume (mL/kg)",                     "IAH above ~250",       150,   250
)

## =====================================================================
##  UI
## =====================================================================
ui <- fluidPage(
  titlePanel("Major Thermal Burn Injury — QSP model"),
  tags$p(style = "color:#555;margin-top:-8px",
         HTML("One injury · two clocks · and <b>a controller whose gain decays as it works</b>")),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      h4("Patient"),
      sliderInput("WT",   "Weight (kg)",      30, 150, 80,  step = 5),
      sliderInput("AGE",  "Age (years)",       1,  95, 35,  step = 1),
      sliderInput("TBSA", "Burn extent %TBSA",        5,  95, 45,  step = 5),
      checkboxInput("INH", "Inhalation injury", FALSE),

      hr(),
      h4("Resuscitation"),
      selectInput("formula", "Starting formula",
                  c("Parkland 4 mL/kg/%"        = 4,
                    "Modified Brooke 2 mL/kg/%" = 2,
                    "ISBI start-low 2 mL/kg/%"  = 2.0001,
                    "No resuscitation (historical)"   = 0), selected = 4),
      checkboxInput("titrate", "Titrated to urine output (closed loop)", TRUE),
      sliderInput("uotgt", "Urine output target (mL/kg/h)", 0.3, 1.5, 0.5, step = 0.1),
      sliderInput("gctrl", "Controller gain",      0.05, 0.9, 0.32, step = 0.01),
      sliderInput("rmax",  "Rate ceiling (× formula)", 1.5, 6.0, 3.0, step = 0.5),
      sliderInput("kopi",  "Opioid-induced fall in urine output", 0, 0.6, 0, step = 0.05),

      hr(),
      h4("Colloid · antioxidant (Starling terms)"),
      sliderInput("fcol",  "Fraction as 5% albumin", 0, 0.6, 0, step = 0.05),
      sliderInput("albstart", "Albumin start (h)", 0, 24, 8, step = 1),
      checkboxInput("vitc", "High-dose ascorbic acid 66 mg/kg/h × 24 h", FALSE),

      hr(),
      h4("Surgery"),
      sliderInput("texc", "Time of first excision (day)", 1, 21, 5, step = 1),
      sliderInput("excint", "Interval between operations (h)", 24, 168, 72, step = 24),
      sliderInput("topical", "Topical antimicrobial effect (/h)", 0.000, 0.045, 0.022, step = 0.002),

      hr(),
      h4("Drugs"),
      sliderInput("prop", "Propranolol (mg/kg/day)", 0, 6, 0, step = 0.5),
      sliderInput("oxa",  "Oxandrolone (mg BID)",       0, 20, 0, step = 5),
      sliderInput("tamb", "Ambient temperature (°C)",           20, 34, 31, step = 1),

      hr(),
      sliderInput("end", "Observation period (days)", 3, 90, 60, step = 1),
      actionButton("go", "Run", class = "btn-primary btn-block")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",

        tabPanel("1 · Patient · protocol",
                 br(), htmlOutput("summary_box"), br(),
                 h4("Summary on admission"), tableOutput("patient_tbl"),
                 h4("What these settings change"),
                 tags$ul(
                   tags$li(HTML("<b>The formula and the urine output target</b> change only the <i>set point</i> of the controller. What is actually given is decided by the leak.")),
                   tags$li(HTML("<b>The albumin fraction</b> does not add volume; it <b>restores f_ret</b>.")),
                   tags$li(HTML("<b>Ascorbic acid</b> acts only on the oxidative component of the Kf lesion.")),
                   tags$li(HTML("<b>The timing of excision</b> removes <b>the driver itself</b> of the slow clock.")),
                   tags$li(HTML("<b>Propranolol</b> blocks not the driver but <b>the transmitter</b> — the two are not additive."))
                 )),

        tabPanel("2 · The convex function",
                 br(), plotOutput("p_convex", height = 420),
                 HTML("<p>Landis–Pappenheimer: <code>Pi = 2.1C + 0.16C² + 0.009C³</code>.
                       Dilute the protein by half and <b>62 %</b> of the colloid osmotic pressure is gone.
                       The red points show where the current simulation travels along this curve
                       over 24 hours. Because the curve is convex, <b>every litre
                       costs more than the one before it</b> — this is the mathematical engine of fluid creep.</p>"),
                 plotOutput("p_fret", height = 300)),

        tabPanel("3 · Fluid · Starling",
                 br(), plotOutput("p_compart", height = 340),
                 plotOutput("p_flux", height = 320),
                 plotOutput("p_pressures", height = 300)),

        tabPanel("4 · The controller",
                 br(), plotOutput("p_ctrl", height = 400),
                 htmlOutput("creep_box"),
                 plotOutput("p_ctrl2", height = 320)),

        tabPanel("5 · Resuscitation morbidity",
                 br(), plotOutput("p_iap", height = 340),
                 htmlOutput("acs_box"),
                 plotOutput("p_wt", height = 320)),

        tabPanel("6 · Hypermetabolism",
                 br(), plotOutput("p_ree", height = 360),
                 htmlOutput("ree_box"),
                 plotOutput("p_neuro", height = 320)),

        tabPanel("7 · Body composition",
                 br(), plotOutput("p_lbm", height = 360),
                 HTML("<p>The net balance is <b>the difference of two large numbers</b>. Move either one by 20 %
                       and the sign changes. That is why propranolol (breakdown ↓, synthesis ↑) and
                       oxandrolone (synthesis ↑) move things in the same direction while touching
                       different terms.</p>"),
                 plotOutput("p_body", height = 320)),

        tabPanel("8 · Wound · closure",
                 br(), plotOutput("p_wound", height = 380),
                 htmlOutput("closure_box"),
                 plotOutput("p_scar", height = 300)),

        tabPanel("9 · Inflammation · infection",
                 br(), plotOutput("p_cyto", height = 340),
                 plotOutput("p_infect", height = 340)),

        tabPanel("10 · Pharmacology · PK",
                 br(), plotOutput("p_pk", height = 340),
                 htmlOutput("pk_box"),
                 plotOutput("p_arc", height = 300)),

        tabPanel("11 · Endpoints · validation",
                 br(), h4("The current simulation against the calibration targets"),
                 tableOutput("valid_tbl"),
                 h4("Revised Baux score (an external validator)"),
                 htmlOutput("rbaux_box"),
                 plotOutput("p_rbaux", height = 340)),

        tabPanel("12 · Scenario comparison",
                 br(),
                 checkboxGroupInput("scens", "Scenarios to compare",
                                    choices = names(SCEN),
                                    selected = c("parkland_fixed", "parkland_titrated",
                                                 "creep_uncapped", "albumin_0h", "vitc"),
                                    inline = TRUE),
                 actionButton("gocmp", "Run the comparison", class = "btn-primary"),
                 br(), br(),
                 tableOutput("cmp_tbl"),
                 plotOutput("p_cmp1", height = 340),
                 plotOutput("p_cmp2", height = 340))
      )
    )
  )
)

## =====================================================================
##  SERVER
## =====================================================================
server <- function(input, output, session) {

  sim <- eventReactive(input$go, {
    fm <- as.numeric(input$formula)
    par <- list(
      WT = input$WT, AGE = input$AGE, TBSA = input$TBSA,
      INH = as.numeric(input$INH),
      FORMULA = if (abs(fm - 2.0001) < 1e-6) 2 else fm,
      TITRATE = as.numeric(input$titrate),
      UOTGT = input$uotgt, GCTRL = input$gctrl, RMAX = input$rmax,
      KOPI = input$kopi,
      FCOL = input$fcol,
      ALBSTART = if (input$fcol > 0) input$albstart else 1e9,
      TEXC = input$texc*24, EXCINT = input$excint,
      KTOP = input$topical, TAMB = input$tamb)

    ev <- surg_events(WT = input$WT, TBSA = input$TBSA,
                      texc = input$texc*24, excint = input$excint)
    dv <- drug_events(WT = input$WT,
                      prop_mgkgday = input$prop, prop_start = 120,
                      prop_n = ceiling(input$end*24/6),
                      oxa_mg = input$oxa, oxa_start = 120,
                      oxa_n = ceiling(input$end*24/12),
                      vitc_mgkgh = if (input$vitc) 66 else 0)
    ev <- if (is.null(ev)) dv else if (is.null(dv)) ev else rbind(ev, dv)

    m <- mod %>% param(par)
    out <- if (is.null(ev)) {
      m %>% mrgsim(end = input$end*24, delta = 1, hmax = 0.05)
    } else {
      m %>% ev(as.ev(ev)) %>% mrgsim(end = input$end*24, delta = 1, hmax = 0.05)
    }
    as_tibble(out)
  }, ignoreNULL = FALSE)

  at24 <- function(d) d[which.min(abs(d$time - 24)), ]

  ## ---------------- 1. patient ----------------
  output$patient_tbl <- renderTable({
    d <- sim(); a <- at24(d)
    tibble::tibble(
      Item = c("Weight", "Age", "%TBSA", "Inhalation injury", "Revised Baux score",
               "rBaux predicted mortality", "24-h volume given", "in:out ratio",
               "Plasma volume nadir", "COP at 24 h", "Peak REE", "Wound closure", "Model mortality"),
      Value = c(sprintf("%.0f kg", d$WT[1]), sprintf("%.0f y", d$AGE[1]),
             sprintf("%.0f %%", d$TBSA[1]), ifelse(d$INH[1] > 0, "Yes", "No"),
             sprintf("%.0f", a$RBAUX), sprintf("%.1f %%", a$RBMORT),
             sprintf("%.2f mL/kg/%%  (%.0f mL/kg)", a$MLKGP, a$MLKG),
             sprintf("%.2f", a$INOUT),
             sprintf("%.0f %% of baseline", min(d$VPPCT[d$time <= 24])),
             sprintf("%.1f mmHg", a$COPt),
             sprintf("%.0f %% of predicted", max(d$REEPCT)),
             { z <- which(d$AOPEN <= 0.05*d$TBSA[1])
               if (length(z)) sprintf("%.0f days", d$time[z[1]]/24) else "Not closed" },
             sprintf("%.1f %%", tail(d$MORT, 1))))
  })

  output$summary_box <- renderUI({
    d <- sim(); a <- at24(d)
    col <- if (a$INOUT > 1.6) "#c0392b" else if (a$INOUT > 1.3) "#b9770e" else "#1e8449"
    HTML(sprintf(
      "<div style='padding:12px;border-left:6px solid %s;background:#fbfbfb'>
       <b>24-h volume %.2f mL/kg/%%TBSA</b> (<b>%.2f times</b> the Parkland prescription) ·
       plasma volume nadir <b>%.0f %%</b> · colloid osmotic pressure at 24 h <b>%.1f mmHg</b> ·
       peak intra-abdominal pressure <b>%.1f mmHg</b> · peak REE <b>%.0f %%</b> ·
       60-day mortality <b>%.1f %%</b> (rBaux prediction %.1f %%)</div>",
      col, a$MLKGP, a$INOUT, min(d$VPPCT[d$time <= 24]), a$COPt,
      max(d$IAPt), max(d$REEPCT), tail(d$MORT, 1), a$RBMORT))
  })

  ## ---------------- 2. the convex function ----------------
  output$p_convex <- renderPlot({
    d <- sim()
    cc <- seq(0.5, 8, by = 0.05)
    traj <- d %>% filter(time <= 24)
    ggplot() +
      geom_line(aes(cc, LP(cc)), linewidth = 1.1, colour = "#2471a3") +
      geom_path(data = traj, aes(CPt, COPt), colour = "#c0392b", linewidth = 1.2) +
      geom_point(data = traj[c(1, nrow(traj)), ], aes(CPt, COPt),
                 colour = "#c0392b", size = 3) +
      geom_segment(aes(x = 7, xend = 3.5, y = LP(7), yend = LP(7)),
                   linetype = 3, colour = "grey40") +
      geom_segment(aes(x = 3.5, xend = 3.5, y = LP(7), yend = LP(3.5)),
                   linetype = 3, colour = "grey40") +
      annotate("text", x = 5.2, y = LP(7) + 1.5,
               label = "protein -50 %  ->  COP -62 %", colour = "grey30") +
      labs(x = "Total protein (g/dL)",
           y = "Colloid osmotic pressure COP (mmHg)",
           title = "The cost convexity creates: every litre is dearer than the one before") + THEME
  })

  output$p_fret <- renderPlot({
    d <- sim() %>% filter(time <= 48)
    ggplot(d, aes(time)) +
      geom_line(aes(y = COPt, colour = "Plasma COP (mmHg)"), linewidth = 1) +
      geom_line(aes(y = CPt*4, colour = "Total protein ×4 (g/dL)"), linewidth = 1) +
      labs(x = "Time (h)", y = "", title = "The collapse of COP — which is the collapse of f_ret") + THEME
  })

  ## ---------------- 3. fluid ----------------
  output$p_compart <- renderPlot({
    d <- sim() %>% filter(time <= 168)
    d %>% select(time, VP, VIB, VIU, VASC, EVLW) %>%
      pivot_longer(-time) %>%
      mutate(name = recode(name, VP = "Plasma", VIB = "Interstitium (burned)",
                           VIU = "Interstitium (unburned)", VASC = "Abdominal cavity", EVLW = "Lung water")) %>%
      ggplot(aes(time, value/1000, fill = name)) +
      geom_area(alpha = 0.85) +
      labs(x = "Time (h)", y = "Volume (L)", title = "Where did the litres go") + THEME
  })

  output$p_flux <- renderPlot({
    d <- sim() %>% filter(time <= 72)
    ggplot(d, aes(time)) +
      geom_line(aes(y = VPPCT, colour = "Plasma volume (% of baseline)"), linewidth = 1) +
      geom_hline(yintercept = 100, linetype = 3) +
      labs(x = "Time (h)", y = "%", title = "Plasma volume — not what was given but what was retained") + THEME
  })

  output$p_pressures <- renderPlot({
    d <- sim() %>% filter(time <= 48)
    ggplot(d, aes(time)) +
      geom_line(aes(y = COPt, colour = "Plasma COP"), linewidth = 1) +
      geom_line(aes(y = IAPt, colour = "Intra-abdominal pressure"), linewidth = 1) +
      geom_hline(yintercept = c(12, 20), linetype = c(3, 2), colour = "#c0392b") +
      annotate("text", x = 2, y = 21, label = "ACS", colour = "#c0392b", hjust = 0) +
      labs(x = "Time (h)", y = "mmHg", title = "The pressures") + THEME
  })

  ## ---------------- 4. controller ----------------
  output$p_ctrl <- renderPlot({
    d <- sim() %>% filter(time <= 48)
    tgt <- input$uotgt*input$WT
    ggplot(d, aes(time)) +
      geom_line(aes(y = RSTATE*100, colour = "Rate multiplier ×100"), linewidth = 1.1) +
      geom_line(aes(y = UOWIN, colour = "Measured urine output (mL/h)"), linewidth = 1) +
      geom_line(aes(y = COPt*8, colour = "Plasma COP ×8"), linewidth = 1) +
      geom_hline(yintercept = tgt, linetype = 2, colour = "grey30") +
      annotate("text", x = 44, y = tgt + 8, label = "Set point", colour = "grey30") +
      labs(x = "Time (h)", y = "",
           title = "The controller pushes harder on the very gain it is destroying") + THEME
  })

  output$creep_box <- renderUI({
    d <- sim(); a <- at24(d)
    msg <- if (a$INOUT > 1.6)
      "<b style='color:#c0392b'>Fluid creep</b> — the controller has hit its ceiling."
    else if (a$INOUT > 1.3)
      "<b style='color:#b9770e'>Borderline range</b> — within the in:out 1.2–1.6 of the literature."
    else "<b style='color:#1e8449'>Close to the prescribed volume</b>."
    HTML(sprintf("<div style='padding:10px;background:#f7f9fa'>%s Final in:out <b>%.2f</b>,
                  total <b>%.0f mL/kg</b>. Above 250 mL/kg the risk of intra-abdominal hypertension rises steeply
                  (model peak intra-abdominal pressure %.1f mmHg).</div>",
                 msg, a$INOUT, a$MLKG, max(d$IAPt)))
  })

  output$p_ctrl2 <- renderPlot({
    d <- sim() %>% filter(time <= 72)
    ggplot(d, aes(time)) +
      geom_line(aes(y = UOC/1000, colour = "Cumulative urine (L)"), linewidth = 1) +
      geom_line(aes(y = FIN/1000, colour = "Cumulative infused (L)"), linewidth = 1) +
      labs(x = "Time (h)", y = "L", title = "What went in and what came out") + THEME
  })

  ## ---------------- 5. morbidity ----------------
  output$p_iap <- renderPlot({
    d <- sim() %>% filter(time <= 168)
    ggplot(d, aes(time, IAPt)) + geom_line(linewidth = 1.1, colour = "#c0392b") +
      geom_hline(yintercept = c(12, 20), linetype = c(3, 2)) +
      annotate("text", x = 5, y = 12.8, label = "Intra-abdominal hypertension (IAH)", hjust = 0, size = 3.5) +
      annotate("text", x = 5, y = 20.8, label = "Abdominal compartment syndrome (ACS)", hjust = 0, size = 3.5) +
      labs(x = "Time (h)", y = "Intra-abdominal pressure (mmHg)",
           title = "Everything in this graph was made by the treatment") + THEME
  })

  output$acs_box <- renderUI({
    d <- sim()
    HTML(sprintf("<div style='padding:10px;background:#fdf2f0'>Peak intra-abdominal pressure
      <b>%.1f mmHg</b> · peak lung water <b>%.0f mL</b> · peak weight gain <b>%.0f %%</b>.
      Intra-abdominal pressure cuts renal perfusion and so lowers urine output, and the controller reads that as hypovolaemia and
      <b>gives more</b> — this is the second, positive feedback loop.</div>",
      max(d$IAPt), max(d$EVLW), max(d$WTGAIN)))
  })

  output$p_wt <- renderPlot({
    d <- sim() %>% filter(time <= 336)
    ggplot(d, aes(time/24)) +
      geom_line(aes(y = WTGAIN, colour = "Fluid weight gain (%)"), linewidth = 1.1) +
      geom_line(aes(y = OEDEMA, colour = "Extracellular oedema (L)"), linewidth = 1) +
      labs(x = "Days", y = "", title = "The formation and mobilisation of oedema") + THEME
  })

  ## ---------------- 6. hypermetabolism ----------------
  output$p_ree <- renderPlot({
    d <- sim()
    ggplot(d, aes(time/24)) +
      geom_line(aes(y = REEPCT, colour = "REE (% of predicted)"), linewidth = 1.2) +
      geom_line(aes(y = AOPENT + 100, colour = "Open wound area + 100 (%TBSA)"),
                linewidth = 1, linetype = 2) +
      geom_hline(yintercept = 100, linetype = 3) +
      labs(x = "Days", y = "",
           title = "The slow clock: the driver is not the %TBSA on admission but today's open area") + THEME
  })

  output$ree_box <- renderUI({
    d <- sim()
    z <- which(d$AOPEN <= 0.05*d$TBSA[1])
    cl <- if (length(z)) d$time[z[1]]/24 else NA
    HTML(sprintf("<div style='padding:10px;background:#fefaf0'>Peak REE
      <b>%.0f %%</b> (on day %s) · at end of observation <b>%.0f %%</b> · wound closure <b>%s</b>.
      The hypermetabolism persists <b>after</b> the wound has closed (model decay time constant about 104 days) —
      this is the best-documented feature of the syndrome.</div>",
      max(d$REEPCT), round(d$time[which.max(d$REEPCT)]/24),
      tail(d$REEPCT, 1), ifelse(is.na(cl), "not closed", sprintf("%.0f days", cl))))
  })

  output$p_neuro <- renderPlot({
    d <- sim() %>% filter(time <= 720)
    ggplot(d, aes(time/24)) +
      geom_line(aes(y = CAT, colour = "Catecholamines (× normal)"), linewidth = 1) +
      geom_line(aes(y = COR/4, colour = "Cortisol /4 (µg/dL)"), linewidth = 1) +
      geom_line(aes(y = TCORE - 33, colour = "Core temperature −33 (°C)"), linewidth = 1) +
      labs(x = "Days", y = "", title = "The transmitter layer") + THEME
  })

  ## ---------------- 7. body composition ----------------
  output$p_lbm <- renderPlot({
    d <- sim()
    ggplot(d, aes(time/24, LBMPCT)) +
      geom_line(linewidth = 1.2, colour = "#7d3c98") +
      geom_hline(yintercept = 0, linetype = 3) +
      geom_hline(yintercept = -10, linetype = 2, colour = "#c0392b") +
      annotate("text", x = 2, y = -11.5, hjust = 0, size = 3.5, colour = "#c0392b",
               label = "−10 %: the point at which wound healing starts to fail") +
      labs(x = "Days", y = "Change in lean body mass (%)",
           title = "The net balance is the difference of two large numbers — hence extremely parameter-sensitive") + THEME
  })

  output$p_body <- renderPlot({
    d <- sim() %>% filter(time <= 1440)
    ggplot(d, aes(time/24)) +
      geom_line(aes(y = LBM, colour = "Lean body mass (kg)"), linewidth = 1) +
      geom_line(aes(y = FATM, colour = "Fat mass (kg)"), linewidth = 1) +
      geom_line(aes(y = BMC*20, colour = "Bone mineral content ×20"), linewidth = 1) +
      geom_line(aes(y = HFAT/40, colour = "Liver fat /40 (g)"), linewidth = 1) +
      labs(x = "Days", y = "", title = "Body composition and the liver") + THEME
  })

  ## ---------------- 8. wound ----------------
  output$p_wound <- renderPlot({
    d <- sim()
    d %>% select(time, AOPEN, AGRF, AHEAL, ADON) %>%
      pivot_longer(-time) %>%
      mutate(name = recode(name, AOPEN = "Open wound", AGRF = "Grafted, awaiting take",
                           AHEAL = "Healed", ADON = "Donor site (open)")) %>%
      ggplot(aes(time/24, value, fill = name)) + geom_area(alpha = 0.85) +
      labs(x = "Days", y = "%TBSA",
           title = "Wound closure — the donor site is the rate-limiting resource in a large burn") + THEME
  })

  output$closure_box <- renderUI({
    d <- sim()
    z <- which(d$AOPEN <= 0.05*d$TBSA[1])
    pool <- (100 - d$TBSA[1])*0.60*3
    HTML(sprintf("<div style='padding:10px;background:#f1f9f4'>Available donor site
      (60 %% of unburned skin, meshed 1:3) = equivalent of <b>%.0f %%TBSA</b> · area required
      <b>%.0f %%TBSA</b> · time to closure <b>%s</b>. When the area required exceeds the donor site, one must wait
      for the donor site to heal again (about 7 days), and that is why time to closure is
      severely non-linear in %%TBSA.</div>",
      pool, d$TBSA[1],
      ifelse(length(z), sprintf("%.0f days", d$time[z[1]]/24), "not closed within the observation period")))
  })

  output$p_scar <- renderPlot({
    d <- sim()
    ggplot(d, aes(time/24, SCARC)) + geom_line(linewidth = 1.1, colour = "#c0392b") +
      labs(x = "Days", y = "Scar collagen (a.u.)",
           title = "Hypertrophic scar — the risk follows the time taken to closure") + THEME
  })

  ## ---------------- 9. inflammation / infection ----------------
  output$p_cyto <- renderPlot({
    d <- sim() %>% filter(time <= 720)
    ggplot(d, aes(time/24)) +
      geom_line(aes(y = IL6, colour = "IL-6 (pg/mL)"), linewidth = 1) +
      geom_line(aes(y = TNFA*4, colour = "TNF-α ×4"), linewidth = 1) +
      geom_line(aes(y = CRP/2, colour = "CRP /2 (mg/L)"), linewidth = 1) +
      geom_line(aes(y = HLADR*2, colour = "Monocyte HLA-DR ×2 (%)"), linewidth = 1) +
      labs(x = "Days", y = "",
           title = "IL-6 suppresses albumin synthesis — albumin is a negative acute-phase protein") + THEME
  })

  output$p_infect <- renderPlot({
    d <- sim() %>% filter(time <= 1440)
    ggplot(d, aes(time/24)) +
      geom_line(aes(y = BWD, colour = "Wound bacterial count (log₁₀ CFU/g)"), linewidth = 1.1) +
      geom_line(aes(y = BSYS, colour = "Systemic bacterial load (a.u.)"), linewidth = 1.1) +
      geom_hline(yintercept = 5, linetype = 2, colour = "#c0392b") +
      annotate("text", x = 1, y = 5.3, hjust = 0, size = 3.5, colour = "#c0392b",
               label = "10⁵ CFU/g — the classical quantitative threshold for invasive infection") +
      labs(x = "Days", y = "", title = "Density is density; invasion is area") + THEME
  })

  ## ---------------- 10. pharmacology ----------------
  output$p_pk <- renderPlot({
    d <- sim() %>% filter(time <= 720)
    ggplot(d, aes(time/24)) +
      geom_line(aes(y = CPRt, colour = "Propranolol (ng/mL)"), linewidth = 1) +
      geom_line(aes(y = CVANt, colour = "Vancomycin (mg/L)"), linewidth = 1) +
      labs(x = "Days", y = "Concentration", title = "Drug concentrations") + THEME
  })

  output$pk_box <- renderUI({
    d <- sim()
    i <- which.min(abs(d$time - 288))
    cp <- d$CPt[i]
    HTML(sprintf("<div style='padding:10px;background:#eefaf8'>Total protein on day 12
      <b>%.2f g/dL</b>. For a 90 %%-bound drug the unbound fraction rises from 0.10 to <b>%.3f</b>.
      For a low-extraction drug <code>CL = fu × CL_int</code>, so at steady state the <b>unbound</b>
      concentration does not change and only the <b>total</b> concentration falls — TDM on total
      concentration therefore <b>over-diagnoses</b> underexposure. Conversely augmented renal clearance (×%.2f)
      lowers the unbound concentration itself, so TDM <b>under-diagnoses</b> it. The two
      characteristic PK changes of burn patients act with <b>opposite signs</b> on the quantity that matters.</div>",
      cp, min(0.10*(7/max(cp, 0.5)), 1), 1 + 0.55))
  })

  output$p_arc <- renderPlot({
    d <- sim() %>% filter(time <= 720)
    ggplot(d, aes(time/24)) +
      geom_line(aes(y = SCR, colour = "Serum creatinine (mg/dL)"), linewidth = 1) +
      geom_line(aes(y = LAC, colour = "Lactate (mmol/L)"), linewidth = 1) +
      geom_line(aes(y = GLC/50, colour = "Blood glucose /50 (mg/dL)"), linewidth = 1) +
      labs(x = "Days", y = "", title = "Renal function · perfusion · glycaemia") + THEME
  })

  ## ---------------- 11. validation ----------------
  output$valid_tbl <- renderTable({
    d <- sim(); a <- at24(d)
    vals <- c(a$MLKGP, a$INOUT, min(d$VPPCT[d$time <= 24]), a$COPt,
              max(d$WTGAIN), max(d$REEPCT),
              d$LBMPCT[which.min(abs(d$time - 336))], a$MLKG)
    TARGETS %>% mutate(
      Model = sprintf("%.2f", vals),
      Verdict = ifelse(vals >= lo & vals <= hi, "○ within range", "△ outside range")) %>%
      select(Metric = metric, Literature = lit, Model, Verdict)
  })

  output$rbaux_box <- renderUI({
    d <- sim(); a <- at24(d)
    HTML(sprintf("<div style='padding:10px;background:#f7f9fa'>Revised Baux score
      = age + %%TBSA + 17×inhalation = <b>%.0f</b> → logistic predicted mortality
      <b>%.1f %%</b>. The mechanistic mortality of the model is <b>%.1f %%</b>.
      <b>rBaux is not an input to the model</b> — the risk is computed from the area of open wound and
      the time it stays open, age, inhalation injury, sepsis, compartment syndrome, ARDS and loss of lean body mass,
      and rBaux is used only to score the result.</div>",
      a$RBAUX, a$RBMORT, tail(d$MORT, 1)))
  })

  output$p_rbaux <- renderPlot({
    d <- sim()
    ggplot(d, aes(time/24, MORT)) + geom_line(linewidth = 1.2, colour = "#34495e") +
      geom_hline(aes(yintercept = RBMORT), linetype = 2, colour = "#c0392b") +
      annotate("text", x = 1, y = d$RBMORT[1] + 2, hjust = 0, size = 3.5,
               colour = "#c0392b", label = "rBaux prediction") +
      labs(x = "Days", y = "Cumulative mortality (%)",
           title = "Mechanistically computed mortality vs rBaux") + THEME
  })

  ## ---------------- 12. scenario comparison ----------------
  cmp <- eventReactive(input$gocmp, {
    req(length(input$scens) > 0)
    bind_rows(lapply(input$scens, function(s)
      run_scenario(s, WT = input$WT, AGE = input$AGE, TBSA = input$TBSA,
                   INH = as.numeric(input$INH), end = input$end*24)))
  })

  output$cmp_tbl <- renderTable({
    cmp() %>% group_by(scenario) %>% group_modify(~summarise_run(.x)) %>%
      ungroup() %>% select(-scenario1) %>%
      mutate(across(where(is.numeric), ~round(.x, 2)))
  })

  output$p_cmp1 <- renderPlot({
    cmp() %>% filter(time <= 48) %>%
      ggplot(aes(time, MLKGP, colour = scenario)) + geom_line(linewidth = 1) +
      geom_hline(yintercept = c(4, 5.2, 6.7), linetype = c(2, 3, 3)) +
      labs(x = "Time (h)", y = "Cumulative mL/kg/%TBSA",
           title = "24-h fluid volume — dotted lines are the Parkland prescription (4) and the observed range (5.2–6.7)") + THEME
  })

  output$p_cmp2 <- renderPlot({
    cmp() %>% ggplot(aes(time/24, REEPCT, colour = scenario)) +
      geom_line(linewidth = 1) +
      labs(x = "Days", y = "REE (% of predicted)",
           title = "The slow clock") + THEME
  })
}

shinyApp(ui, server)
