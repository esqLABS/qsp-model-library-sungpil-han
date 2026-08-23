##############################################################################
## csdh_shiny_app.R
## Chronic Subdural Haematoma (cSDH) QSP Interactive Dashboard
## Chronic Subdural Haematoma — QSP interactive dashboard
##
## The app is organised around the model's one structural claim: the haematoma
## volume is the FIXED POINT of
##        dV/dt = J_exudation + J_rebleed - J_absorption
## so the first thing every tab shows is which of those three terms the user's
## chosen therapy is actually moving.  Tab 3 (Flux Balance) is the tab to read
## first; everything else is downstream of it.
##
## Run:
##   shiny::runApp("csdh_shiny_app.R")
## Requires: shiny, mrgsolve, dplyr, tidyr, ggplot2, DT
##############################################################################

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

source("csdh_mrgsolve_model.R")   # defines mod, csdh_sim, ev_* , CSDH_* helpers

THEME <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "#EEF2F7", colour = NA),
        legend.position = "bottom")

PAL <- c("#2C5F8D", "#8B2635", "#276749", "#9C5B1E", "#553C7B",
         "#B8860B", "#C9548A", "#333333")

##############################################################################
## UI
##############################################################################

ui <- fluidPage(
  titlePanel("Chronic Subdural Haematoma (cSDH) QSP Simulator"),
  tags$p(style = "color:#8B0000;font-weight:600;margin-top:-8px;",
         "dV/dt = J_exudation + J_rebleed − J_absorption  ·  ",
         "The haematoma is not a clot that failed to be absorbed, but a 'secreting organ'"),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      h4("① Patient"),
      sliderInput("age", "Age (yr)", 45, 95, 76, 1),
      sliderInput("vres", "Intracranial reserve space V_RES (mL) — atrophy",
                  0, 60, 30, 1),
      helpText(tags$small(
        "V_RES enters the model twice, with opposite signs: it buffers pressure, letting a larger ",
        "haematoma be tolerated (beneficial), while also slowing brain re-expansion, leaving a space ",
        "behind after drainage (harmful).")),
      sliderInput("vpres", "Volume at presentation V_pres (mL)", 40, 140, 78, 2),

      h4("② Procedure"),
      selectInput("surg", "Surgery",
                  c("None (conservative)" = "none",
                    "Burr-hole drainage (no drain)" = "bh",
                    "Burr-hole drainage + drain 48h" = "bhd",
                    "Craniotomy + membrane excision" = "crani")),
      sliderInput("wash", "Wash degree wash (lower = cleaner)", 0.05, 1.0, 0.15, 0.05),
      helpText(tags$small("wash lowers the 'concentration' inside the cavity — the manipulation that switches off the haem switch.")),
      checkboxInput("mmae", "MMA embolisation (MMAE)", FALSE),
      sliderInput("mmae_t", "MMAE timing (days)", 0, 60, 0, 1),

      h4("③ Drugs"),
      checkboxInput("dex", "Dexamethasone 16 mg/d × 2-week taper", FALSE),
      checkboxInput("atv", "Atorvastatin 20 mg × 8 weeks", FALSE),
      checkboxInput("txa", "Tranexamic acid 750 mg/d", FALSE),
      checkboxInput("doac", "Resume apixaban", FALSE),
      sliderInput("doac_t", "DOAC resumption timing (days)", 0, 90, 30, 1),

      h4("④ Simulation"),
      sliderInput("tend", "Observation period (days)", 60, 365, 180, 5),
      actionButton("go", "Run", class = "btn-primary btn-block"),
      hr(),
      helpText(tags$small(
        "Every curve is output from the 50-ODE model; parameters were calibrated to Santarius 2009 · ",
        "ATOCH 2018 · Dex-CSDH 2020 · EMBOLISE 2024. ",
        "It is for education and research and must not be used for clinical decision-making."))
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs", type = "tabs",

        ## ------------------------------------------------------------
        tabPanel("1. Patient Profile",
          br(),
          fluidRow(
            column(6, h4("Natural history up to presentation (run-in)"),
                   plotOutput("p_runin", height = "300px")),
            column(6, h4("Summary"), tableOutput("t_summary"))),
          hr(),
          h4("The Two Signs of V_RES in This Patient"),
          plotOutput("p_tworole", height = "300px"),
          helpText("The left axis (midline shift at presentation) and the right axis (reoperation probability) ",
                   "move in opposite directions as atrophy increases.")
        ),

        ## ------------------------------------------------------------
        tabPanel("2. Drug PK",
          br(), h4("Plasma and Intra-Haematoma Concentrations"),
          plotOutput("p_pk", height = "420px"),
          hr(),
          h4("TXA: Plasma t½ 3 Hours vs Site-of-Action t½ ≈ 6 Days"),
          plotOutput("p_txa", height = "300px"),
          helpText("TXA enters the cavity not by diffusion but carried in the exudate (J_ex), ",
                   "and is washed out by absorption (J_abs). The site half-life is therefore set by V/J_abs, ",
                   "and takes 1–2 weeks regardless of plasma pharmacokinetics.")
        ),

        ## ------------------------------------------------------------
        tabPanel("3. Flux Balance ★",
          br(),
          h4("dV/dt = J_ex + J_rb − J_abs — Time Course of the Three Terms"),
          plotOutput("p_flux", height = "360px"),
          hr(),
          fluidRow(
            column(6, h4("Net Rate of Change dV/dt"), plotOutput("p_net", height="280px")),
            column(6, h4("Starling Factors"), plotOutput("p_starling", height="280px"))),
          helpText("Which term a treatment touches is the whole point of this model. Surgery acts on STOCK, ",
                   "embolisation on SOURCE, and statins move Lp and σ simultaneously.")
        ),

        ## ------------------------------------------------------------
        tabPanel("4. Haematoma Volume · Imaging",
          br(),
          fluidRow(
            column(6, h4("Volume (mL)"), plotOutput("p_vol", height="290px")),
            column(6, h4("Maximum Thickness (mm)"), plotOutput("p_thick", height="290px"))),
          hr(),
          fluidRow(
            column(6, h4("CT Density (HU)"), plotOutput("p_hu", height="290px")),
            column(6, h4("Membrane Area · Patency"), plotOutput("p_area", height="290px")))
        ),

        ## ------------------------------------------------------------
        tabPanel("5. Neomembrane · Neovascularisation",
          br(),
          fluidRow(
            column(6, h4("Vessel Populations"), plotOutput("p_vessel", height="290px")),
            column(6, h4("Maturity and Leakiness Factor"), plotOutput("p_mature", height="290px"))),
          hr(),
          h4("Cytokines · Proteases"),
          plotOutput("p_cyto", height = "300px")
        ),

        ## ------------------------------------------------------------
        tabPanel("6. Haem Switch (Bistability) ★",
          br(),
          h4("Haem Drive haem_drive = hb²/(1+hb²)"),
          plotOutput("p_switch", height = "300px"),
          hr(),
          h4("Changing the Wash Degree Splits the Same Surgery into Healing or Recurrence"),
          plotOutput("p_bistable", height = "320px"),
          helpText("Lowering the Hill exponent to 1 collapses every curve in this figure onto one, ",
                   "and the model then has only a single fixed point that no surgery can cure.")
        ),

        ## ------------------------------------------------------------
        tabPanel("7. Fibrinolysis Loop",
          br(),
          fluidRow(
            column(6, h4("t-PA · Plasmin"), plotOutput("p_lysis", height="290px")),
            column(6, h4("FDP and Haemostatic Capacity H"), plotOutput("p_fdp", height="290px"))),
          hr(),
          h4("Amplification Factor 1/(1−(1−H)) — How Many-Fold the Loop Multiplies Rebleeding"),
          plotOutput("p_gain", height = "280px")
        ),

        ## ------------------------------------------------------------
        tabPanel("8. Brain Mechanics",
          br(),
          fluidRow(
            column(6, h4("ICP (Bidirectional)"), plotOutput("p_icp", height="290px")),
            column(6, h4("Midline Shift (mm)"), plotOutput("p_mls", height="290px"))),
          hr(),
          fluidRow(
            column(6, h4("Compression Deficit and Effective Reserve"),
                   plotOutput("p_comp", height="290px")),
            column(6, h4("Residual Space = Substrate for Recurrence"),
                   plotOutput("p_space", height="290px")))
        ),

        ## ------------------------------------------------------------
        tabPanel("9. Clinical Endpoints",
          br(),
          fluidRow(
            column(6, h4("Reoperation Probability"), plotOutput("p_reop", height="290px")),
            column(6, h4("Probability of Favourable Outcome (mRS 0–3)"),
                   plotOutput("p_fav", height="290px"))),
          hr(),
          fluidRow(
            column(6, h4("Symptoms · Cognition"), plotOutput("p_symp", height="290px")),
            column(6, h4("Safety: Glucose · Infection · Myopathy · Thrombosis"),
                   plotOutput("p_safety", height="290px")))
        ),

        ## ------------------------------------------------------------
        tabPanel("10. Scenario Comparison",
          br(),
          h4("16 Standard Scenarios"),
          checkboxGroupInput("scn_pick", NULL,
            choices = c("Conservative"="S01","Burr-hole drainage"="S02","Burr-hole + drain"="S03",
                        "MMAE alone"="S05","Surgery+MMAE"="S06","Dexamethasone"="S08",
                        "Surgery+dex"="S09","Atorvastatin"="S10",
                        "TXA"="S12","Surgery+TXA"="S13","All"="S16"),
            selected = c("S01","S03","S06","S10"), inline = TRUE),
          plotOutput("p_scn", height = "380px"),
          hr(), h4("Endpoint Table"), DTOutput("t_scn")
        ),

        ## ------------------------------------------------------------
        tabPanel("11. MMAE Time Course ★",
          br(),
          h4("Why the 30-Day Endpoint Is Negative and the 180-Day Endpoint Is Positive"),
          plotOutput("p_mmae", height = "340px"),
          hr(),
          h4("Time Course of the Absolute Risk Reduction"),
          plotOutput("p_mmae_ard", height = "300px"),
          helpText("Embolisation eliminates SOURCE, not STOCK. Existing neovessels ",
                   "decay with τ ≈ 9–18 days, so the reduction in J_ex only shows up through integration, ",
                   "and the difference at an early endpoint is therefore inevitably small. ",
                   "The discrepancy between EMBOLISE (180 days, positive) and MAGIC-MT (90 days, negative) ",
                   "follows from this structure and was not built in as an assumption.")
        ),

        ## ------------------------------------------------------------
        tabPanel("12. Literature Comparison",
          br(), h4("Model vs. Published Trial Figures"),
          DTOutput("t_ledger"),
          hr(),
          helpText("This table is the model's report card, and rows that disagree are not hidden. ",
                   "See also the 'most exposed prediction' section of README.md.")
        )
      )
    )
  )
)

##############################################################################
## SERVER
##############################################################################

server <- function(input, output, session) {

  ## ---- build the requested scenario --------------------------------
  sim <- eventReactive(input$go, {
    doses <- NULL
    if (input$dex)  doses <- c(doses, ev_dex())
    if (input$atv)  doses <- c(doses, ev_atv())
    if (input$txa)  doses <- c(doses, ev_txa(days = input$tend))
    if (input$doac) doses <- c(doses, ev_doac(t0 = input$doac_t,
                                              days = input$tend - input$doac_t))
    surg <- switch(input$surg, none = NA, bh = 0, bhd = 0, crani = 0)
    dd   <- if (identical(input$surg, "bhd")) 2 else 0
    pars <- list(V_RES = input$vres, AGE = input$age)
    if (identical(input$surg, "crani")) {
      pars$K_FUSEX <- 0.20     # membranectomy removes extent
      pars$K_SEPT  <- 0.002
    }
    withProgress(message = "Integrating the 50-ODE system...", value = 0.4, {
      csdh_sim("user", tend = input$tend, surgery = surg,
               drain_days = dd, wash = input$wash,
               embolise_at = if (input$mmae) input$mmae_t else NA,
               doses = doses, pars = pars, V_pres = input$vpres)
    })
  }, ignoreNULL = FALSE)

  runin <- eventReactive(input$go, {
    m <- param(mod, list(V_RES = input$vres))
    m %>% init(csdh_fresh()) %>% mrgsim(end = 260, delta = 0.5) %>% as_tibble()
  }, ignoreNULL = FALSE)

  lineplot <- function(d, cols, labs, ylab, hline = NULL) {
    dd <- d %>% select(time, all_of(cols)) %>%
      pivot_longer(-time, names_to = "v", values_to = "y") %>%
      mutate(v = factor(v, levels = cols, labels = labs))
    g <- ggplot(dd, aes(time, y, colour = v)) +
      geom_line(linewidth = 0.85) +
      scale_colour_manual(values = PAL[seq_along(cols)], name = NULL) +
      labs(x = "Time (days)", y = ylab) + THEME
    if (!is.null(hline))
      g <- g + geom_hline(yintercept = hline, linetype = 2, colour = "#8B0000")
    g
  }

  ## ---- 1. patient profile ------------------------------------------
  output$p_runin <- renderPlot({
    d <- runin()
    lineplot(d, c("VHEM"), c("Haematoma volume"), "mL",
             hline = input$vpres) +
      annotate("text", x = 5, y = input$vpres * 1.06,
               label = "Presentation threshold", hjust = 0, size = 3.2, colour = "#8B0000")
  })

  output$t_summary <- renderTable({
    d <- sim(); r <- runin()
    i <- which(r$VHEM >= input$vpres)[1]
    data.frame(
      Item = c("Latency to presentation (days)", "Thickness at presentation (mm)",
               "Midline shift at presentation (mm)", "Volume at end of observation (mL)",
               "Thickness at end of observation (mm)", "Reoperation probability",
               "Favourable outcome probability", "Maximum glucose (mmol/L)"),
      Value = c(sprintf("%.1f", if (is.na(i)) NA else r$time[i]),
             sprintf("%.1f", d$dmax_[1]), sprintf("%.2f", d$MLS_[1]),
             sprintf("%.1f", tail(d$VHEM, 1)), sprintf("%.1f", tail(d$dmax_, 1)),
             sprintf("%.3f", tail(d$P_REOP, 1)),
             sprintf("%.3f", tail(d$P_FAV, 1)), sprintf("%.1f", max(d$GLU))))
  }, striped = TRUE)

  output$p_tworole <- renderPlot({
    vs <- c(5, 12, 20, 30, 40, 50, 60)
    rows <- lapply(vs, function(vr) {
      p <- list(V_RES = vr)
      d <- csdh_sim("s", 180, surgery = 0, drain_days = 2, pars = p)
      data.frame(V_RES = vr, MLS = d$MLS_[1],
                 reop = tail(d$P_REOP, 1))
    })
    dd <- bind_rows(rows)
    ggplot(dd, aes(V_RES)) +
      geom_line(aes(y = MLS, colour = "MLS at presentation (mm) — atrophy protects"),
                linewidth = 1) +
      geom_point(aes(y = MLS, colour = "MLS at presentation (mm) — atrophy protects")) +
      geom_line(aes(y = reop * 25, colour = "Reoperation probability ×25 — atrophy harms"),
                linewidth = 1) +
      geom_point(aes(y = reop * 25, colour = "Reoperation probability ×25 — atrophy harms")) +
      scale_y_continuous("Midline shift at presentation (mm)",
                         sec.axis = sec_axis(~./25, name = "Reoperation probability")) +
      scale_colour_manual(values = c("#2C5F8D", "#8B2635"), name = NULL) +
      labs(x = "Intracranial reserve space V_RES (mL)") + THEME
  })

  ## ---- 2. PK --------------------------------------------------------
  output$p_pk <- renderPlot({
    d <- sim()
    p <- as.list(param(mod))
    dd <- d %>% transmute(time,
      `Dexamethasone plasma (µg/L)` = DEXC / p$DEX_V1 * 1000,
      `Atorvastatin plasma (µg/L)` = ATVC / p$ATV_V1 * 1000,
      `TXA plasma (mg/L)` = TXAC / p$TXA_V1,
      `TXA cavity (mg/L)` = TXAH / pmax(VHEM, 1e-4) * 1000,
      `Apixaban plasma (ng/mL)` = DOACC / p$DOAC_V1 * 1000,
      `Apixaban cavity (ng/mL)` = DOACH / pmax(VHEM, 1e-4) * 1e6) %>%
      pivot_longer(-time, names_to = "v", values_to = "y")
    ggplot(dd, aes(time, y, colour = v)) + geom_line(linewidth = 0.8) +
      facet_wrap(~v, scales = "free_y", ncol = 3) +
      scale_colour_manual(values = rep(PAL, 3), guide = "none") +
      labs(x = "Time (days)", y = NULL) + THEME
  })

  output$p_txa <- renderPlot({
    d <- csdh_sim("txa", 120, surgery = 0, drain_days = 2,
                  doses = ev_txa(days = 120),
                  pars = list(V_RES = input$vres))
    p <- as.list(param(mod))
    dd <- d %>% transmute(time, Plasma = TXAC / p$TXA_V1,
                          Cavity = TXAH / pmax(VHEM, 1e-4) * 1000) %>%
      pivot_longer(-time, names_to = "v", values_to = "y")
    ggplot(dd, aes(time, y, colour = v)) + geom_line(linewidth = 0.9) +
      scale_colour_manual(values = c("#2C5F8D", "#8B2635"), name = NULL) +
      labs(x = "Time (days)", y = "TXA concentration (mg/L)") + THEME
  })

  ## ---- 3. flux ledger ----------------------------------------------
  output$p_flux <- renderPlot({
    lineplot(sim(), c("JEX", "JRB", "JABS"),
             c("J_ex exudation", "J_rb rebleeding", "J_abs absorption"), "mL/day")
  })
  output$p_net <- renderPlot({
    d <- sim() %>% mutate(net = JEX + JRB - JABS)
    ggplot(d, aes(time, net)) +
      geom_hline(yintercept = 0, linetype = 2, colour = "#8B0000") +
      geom_line(linewidth = 0.9, colour = "#2C5F8D") +
      labs(x = "Time (days)", y = "dV/dt (mL/day)") + THEME
  })
  output$p_starling <- renderPlot({
    lineplot(sim(), c("DPNET", "LPEFF", "SIGEFF"),
             c("ΔP net filtration pressure (mmHg)", "Lp effective conductance", "σ reflection coefficient"),
             "Value")
  })

  ## ---- 4. volume / imaging -----------------------------------------
  output$p_vol   <- renderPlot(lineplot(sim(), "VHEM", "Volume", "mL"))
  output$p_thick <- renderPlot(lineplot(sim(), "dmax_", "Maximum thickness",
                                        "mm", hline = 10))
  output$p_hu    <- renderPlot(lineplot(sim(), "HU_", "CT density", "HU"))
  output$p_area  <- renderPlot(lineplot(sim(), c("A_mem_", "f_pat_"),
                                        c("Membrane area (cm²)", "Patency f_pat"), "Value"))

  ## ---- 5. membrane / vessels ---------------------------------------
  output$p_vessel <- renderPlot(lineplot(sim(), c("NCAP", "NMAT", "NMAC"),
      c("Immature vessels N_CAP", "Mature vessels N_MAT", "Macrophages N_MAC"), "Density index"))
  output$p_mature <- renderPlot(lineplot(sim(), c("SIGEFF", "LPEFF", "NPC"),
      c("σ reflection coefficient", "Lp conductance", "Pericyte coverage"), "Value"))
  output$p_cyto <- renderPlot({
    d <- sim() %>% transmute(time, `VEGF (pg/mL)` = CVEGF,
      `Ang-2 (ng/mL)` = CANG2, `MMP-9 (ng/mL)` = CMMP9,
      `IL-6 (pg/mL)` = CIL6, `IL-8 (pg/mL)` = CIL8) %>%
      pivot_longer(-time, names_to = "v", values_to = "y")
    ggplot(d, aes(time, y, colour = v)) + geom_line(linewidth = 0.8) +
      facet_wrap(~v, scales = "free_y", ncol = 3) +
      scale_colour_manual(values = rep(PAL, 2), guide = "none") +
      labs(x = "Time (days)", y = NULL) + THEME
  })

  ## ---- 6. the switch ------------------------------------------------
  output$p_switch <- renderPlot({
    lineplot(sim(), c("HDRIVE", "CHB_"),
             c("Haem drive (0–1)", "Cavity haemoglobin (g/dL)"), "Value")
  })
  output$p_bistable <- renderPlot({
    ws <- c(0.05, 0.15, 0.30, 0.55, 1.00)
    dd <- bind_rows(lapply(ws, function(w) {
      d <- csdh_sim("w", 180, surgery = 0, drain_days = 2, wash = w,
                    pars = list(V_RES = input$vres))
      data.frame(time = d$time, V = d$VHEM, wash = sprintf("wash %.2f", w))
    }))
    ggplot(dd, aes(time, V, colour = wash)) + geom_line(linewidth = 0.9) +
      geom_hline(yintercept = 0, linetype = 3) +
      scale_colour_manual(values = PAL[seq_along(ws)], name = NULL) +
      labs(x = "Time (days)", y = "Haematoma volume (mL)") + THEME
  })

  ## ---- 7. fibrinolysis ----------------------------------------------
  output$p_lysis <- renderPlot(lineplot(sim(), c("CTPA_", "CPLS", "CPAI"),
      c("t-PA (ng/mL)", "Plasmin activity", "PAI-1 (ng/mL)"), "Value"))
  output$p_fdp <- renderPlot(lineplot(sim(), c("CFDP_", "HCOMP"),
      c("FDP (µg/mL)", "Haemostatic capacity H"), "Value"))
  output$p_gain <- renderPlot({
    d <- sim() %>% mutate(gain = 1 / pmax(HCOMP, 1e-3))
    ggplot(d, aes(time, gain)) + geom_line(linewidth = 0.9, colour = "#553C7B") +
      geom_hline(yintercept = 1, linetype = 2, colour = "#8B0000") +
      labs(x = "Time (days)", y = "1/H  (rebleeding amplification factor)") + THEME
  })

  ## ---- 8. mechanics -------------------------------------------------
  output$p_icp <- renderPlot({
    p <- as.list(param(mod))
    lineplot(sim(), "ICP_", "ICP", "mmHg", hline = p$ICP0)
  })
  output$p_mls  <- renderPlot(lineplot(sim(), "MLS_", "Midline shift", "mm"))
  output$p_comp <- renderPlot(lineplot(sim(), c("VCOMP", "R_"),
      c("Compression deficit V_comp (mL)", "Effective reserve R (mL)"), "mL"))
  output$p_space <- renderPlot({
    d <- sim() %>% mutate(space = pmax(R_ - VHEM, 0))
    ggplot(d, aes(time, space)) +
      geom_area(fill = "#FFD8A8", colour = "#B8860B") +
      labs(x = "Time (days)", y = "Residual subdural space (mL)") + THEME
  })

  ## ---- 9. endpoints --------------------------------------------------
  output$p_reop <- renderPlot(lineplot(sim(), "P_REOP", "Reoperation probability", "Probability"))
  output$p_fav  <- renderPlot(lineplot(sim(), "P_FAV", "Favourable outcome", "Probability"))
  output$p_symp <- renderPlot(lineplot(sim(), c("SSYMP", "SCOG", "NNEUR"),
      c("Symptom burden", "Cognitive deficit", "Neurological integrity"), "Value"))
  output$p_safety <- renderPlot(lineplot(sim(), c("GLU", "XINF", "XMYO", "P_THR"),
      c("Glucose (mmol/L)", "Infection burden", "Myopathy burden", "Thrombosis probability"), "Value"))

  ## ---- 10. scenarios --------------------------------------------------
  scn_data <- reactive({
    keys <- names(CSDH_SCENARIOS)
    sel  <- keys[substr(keys, 1, 3) %in% input$scn_pick]
    if (!length(sel)) return(NULL)
    bind_rows(lapply(sel, function(k) {
      d <- CSDH_SCENARIOS[[k]]()
      data.frame(time = d$time, V = d$VHEM, dmax = d$dmax_,
                 reop = d$P_REOP, fav = d$P_FAV,
                 scenario = attr(d, "name"))
    }))
  })

  output$p_scn <- renderPlot({
    dd <- scn_data(); req(dd)
    ggplot(dd, aes(time, V, colour = scenario)) + geom_line(linewidth = 0.9) +
      scale_colour_manual(values = rep(PAL, 3), name = NULL) +
      labs(x = "Time (days)", y = "Haematoma volume (mL)") + THEME +
      guides(colour = guide_legend(ncol = 2))
  })

  output$t_scn <- renderDT({
    dd <- scn_data(); req(dd)
    dd %>% group_by(scenario) %>%
      summarise(`V at 90 days` = round(approx(time, V, 90, rule=2)$y, 1),
                `V at 180 days` = round(approx(time, V, 180, rule=2)$y, 1),
                `Thickness at 180 days (mm)` = round(approx(time, dmax, 180, rule=2)$y, 1),
                `Reoperation probability` = round(max(reop), 3),
                `Favourable outcome` = round(min(fav), 3), .groups = "drop") %>%
      datatable(options = list(dom = "t", pageLength = 20), rownames = FALSE)
  })

  ## ---- 11. MMAE time course --------------------------------------------
  mmae_pair <- reactive({
    p <- list(V_RES = input$vres)
    list(a = csdh_sim("Surgery only", 180, surgery = 0, drain_days = 2, pars = p),
         b = csdh_sim("Surgery+MMAE", 180, surgery = 0, drain_days = 2,
                      embolise_at = 0, pars = p))
  })

  output$p_mmae <- renderPlot({
    z <- mmae_pair()
    dd <- bind_rows(
      data.frame(time = z$a$time, V = z$a$VHEM, NCAP = z$a$NCAP, arm = "Surgery only"),
      data.frame(time = z$b$time, V = z$b$VHEM, NCAP = z$b$NCAP, arm = "Surgery+MMAE")) %>%
      pivot_longer(c(V, NCAP), names_to = "v", values_to = "y") %>%
      mutate(v = factor(v, c("V", "NCAP"),
                        c("Haematoma volume (mL)", "Immature vessels N_CAP")))
    ggplot(dd, aes(time, y, colour = arm)) + geom_line(linewidth = 0.9) +
      facet_wrap(~v, scales = "free_y") +
      geom_vline(xintercept = c(30, 90, 180), linetype = 3, colour = "#666666") +
      scale_colour_manual(values = c("#8B2635", "#276749"), name = NULL) +
      labs(x = "Time (days)", y = NULL) + THEME
  })

  output$p_mmae_ard <- renderPlot({
    z <- mmae_pair()
    tt <- seq(7, 180, by = 1)
    ard <- approx(z$a$time, z$a$P_REOP, tt, rule=2)$y -
           approx(z$b$time, z$b$P_REOP, tt, rule=2)$y
    ggplot(data.frame(time = tt, ard = ard), aes(time, ard)) +
      geom_hline(yintercept = 0, linetype = 2) +
      geom_line(linewidth = 1, colour = "#276749") +
      geom_vline(xintercept = c(30, 90, 180), linetype = 3, colour = "#666666") +
      annotate("text", x = c(30, 90, 180), y = max(ard) * 0.15,
               label = c("30 days", "90 days\n(MAGIC-MT)", "180 days\n(EMBOLISE)"),
               size = 3.1, hjust = -0.05) +
      labs(x = "Time (days)", y = "Absolute risk reduction (reoperation)") + THEME
  })

  ## ---- 12. ledger -------------------------------------------------------
  output$t_ledger <- renderDT({
    withProgress(message = "Calculating literature comparison...", value = 0.5, {
      rows <- CSDH_trial_ledger()
    })
    data.frame(Endpoint = sapply(rows, `[[`, 1),
               Model = round(as.numeric(sapply(rows, `[[`, 2)), 3),
               Published = sapply(rows, `[[`, 3)) %>%
      datatable(options = list(dom = "t", pageLength = 20), rownames = FALSE)
  })
}

shinyApp(ui, server)
