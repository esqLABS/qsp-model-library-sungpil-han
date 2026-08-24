## =============================================================================
##  myp_shiny_app.R
##  Interactive dashboard for the progressive-myopia QSP model
##
##  10 tabs:
##    1  Patient profile        patient profile, baseline optics, risk stratum
##    2  Growth trajectory      axial length and refraction over time
##    3  Optical decomposition  the three optical elements, and why SER misleads
##    4  Choroid vs sclera      the two actuators that share one measurement
##    5  Atropine PK            ocular and systemic exposure
##    6  Dose-response          efficacy vs side effect, and the ratio
##    7  Scenario comparison    up to six arms side by side
##    8  Rebound                stopping, tapering, receptor up-regulation
##    9  Scleral biology        the matrix cascade behind the creep term
##   10  Clinical endpoints     lifetime risk, CARE, NNT
##
##  Run:  shiny::runApp("myp_shiny_app.R")
##  Needs: shiny, mrgsolve, dplyr, tidyr, ggplot2, DT
## =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

`%||%` <- function(a, b) if (is.null(a)) b else a

MOD <- mread_cache("myp", "myp_mrgsolve_model.R")

## --- the baseline axial length that reproduces a prescribed refraction -----
## The optics block is inverted by bisection, exactly as solve_AL0() does in
## the Python reference model, so that t = 0 reproduces the requested SER.
ser_of <- function(AL, ACD = 3.60, LT = 3.45, PL = 22.60, CR = 7.80,
                   naq = 1.336, nvit = 1.336, vtx = 0.012) {
  PCORN <- (naq - 1) / (CR / 1000)
  d1 <- ACD + LT / 2
  d2 <- pmax(AL - d1, 5)
  L2 <- 1000 * nvit / d2 - PL
  L1 <- L2 / (1 + (d1 / 1000 / naq) * L2)
  FC <- L1 - PCORN
  FC / (1 + vtx * FC)
}
solve_AL0 <- function(SER0, ...) {
  lo <- 18; hi <- 34
  for (i in 1:200) {
    mid <- (lo + hi) / 2
    if (ser_of(mid, ...) > SER0) lo <- mid else hi <- mid
  }
  (lo + hi) / 2
}

DEVICE <- c("Single-vision spectacles"           = 0.15,
            "Progressive addition lens (PAL)"    = -0.30,
            "MiSight dual-focus contact lens"    = -0.90,
            "DIMS spectacles"                    = -1.20,
            "HAL lenslet spectacles"             = -2.20,
            "Under-correction +0.75 D (harmful)" = 0.90)

## simulate one arm; `days` in days
run_arm <- function(age0, ser0, npar, grs, ethn, outd, neard,
                    atro, atrostop, atrotap, trtpd, okon, rlrl,
                    fuiris, adhfix, trtstart, days) {
  al0 <- solve_AL0(ser0)
  pars <- list(AGE0 = age0, SER0 = ser0, AL0 = al0, NPAR = npar, GRS = grs,
               ETHN = ethn, OUTD = outd, NEARD = neard,
               ATROPCT = atro, ATROSTOP = atrostop, ATROTAP = atrotap,
               TRTPD = trtpd, OKON = okon, RLRLON = rlrl,
               FUIRISR = fuiris, ADHFIX = adhfix, TRTSTART = trtstart,
               OPTSTOP = if (atrostop < 1e8 && atro == 0) atrostop else 1e9)
  ev_drops <- if (atro > 0)
    ev(time = 0, amt = atro * 10 * 30 / 289.37 * 1000, cmt = "ATEAR",
       ii = 1, addl = max(1, floor(min(days, atrostop))) - 1) else ev()
  MOD %>% param(pars) %>%
    mrgsim(events = ev_drops, end = days, delta = 7, hmax = 1) %>%
    as_tibble() %>% mutate(yr = time / 365)
}

lifetime_vi <- function(AL) 1 / (1 + exp(-(-2.060 + 0.950 * (AL - 26))))

THEME <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom", legend.title = element_blank(),
        strip.text = element_text(face = "bold"))

## =============================================================================
##  UI
## =============================================================================
ui <- fluidPage(
  titlePanel("Paediatric Myopia Progression — QSP Dashboard"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Patient"),
      sliderInput("age0", "Starting age (yr)", 6, 15, 9.7, 0.1),
      sliderInput("ser0", "Baseline refraction SER (D)", -8, 1, -3, 0.25),
      selectInput("npar", "Number of myopic parents", c(0, 1, 2), selected = 2),
      sliderInput("grs", "Polygenic risk score", 0, 1, 0.70, 0.05),
      radioButtons("ethn", "Ancestry", c("East Asian" = 1, "European" = 0), inline = TRUE),
      sliderInput("outd", "Time outdoors (h/day)", 0.2, 4, 1.0, 0.1),
      sliderInput("neard", "Near work (h/day)", 0.5, 8, 3.0, 0.5),
      radioButtons("iris", "Iris colour", c("Dark" = 1, "Medium" = 2, "Light" = 4),
                   inline = TRUE),
      hr(),
      h4("Treatment"),
      selectInput("atro", "Atropine concentration (% w/v)",
                  c(0, 0.01, 0.025, 0.05, 0.1, 0.5, 1.0), selected = 0),
      selectInput("dev", "Optical prescription", names(DEVICE), selected = names(DEVICE)[1]),
      checkboxInput("okon", "Orthokeratology (ortho-K)", FALSE),
      checkboxInput("rlrl", "Repeated low-level red light (RLRL)", FALSE),
      sliderInput("trtstart", "Treatment start (yr)", 0, 5, 0, 0.5),
      sliderInput("stopyr", "Treatment stop (yr; 10 = never stopped)", 0.5, 10, 10, 0.5),
      sliderInput("tapyr", "Taper duration (yr)", 0, 2, 0, 0.25),
      sliderInput("adh", "Adherence (0 = computed by the model)", 0, 1, 0, 0.05),
      hr(),
      sliderInput("years", "Follow-up (yr)", 1, 11, 8, 1),
      helpText("Every curve is referenced to the measured axial length (cornea→RPE), and the ",
               "refraction is not a state variable but the output of the optical equations.")
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        # ------------------------------------------------------------------
        tabPanel("1. Patient profile",
          br(),
          fluidRow(
            column(6, h4("Baseline optics"), tableOutput("t_base")),
            column(6, h4("Risk stratum"), tableOutput("t_risk"))),
          hr(),
          plotOutput("p_profile", height = "330px"),
          helpText("Once the AL/CR ratio exceeds 3.0 the eye is classified as myopic even ",
                   "while the refraction is still normal. Axial length moves before refraction does.")),
        # ------------------------------------------------------------------
        tabPanel("2. Growth trajectory",
          br(), plotOutput("p_traj", height = "560px"),
          hr(), h4("Annual progression"),
          DTOutput("t_traj")),
        # ------------------------------------------------------------------
        tabPanel("3. Optical decomposition",
          br(),
          h4("The refraction is nothing but the difference between three optical elements"),
          plotOutput("p_optics", height = "400px"),
          hr(),
          h4("dSER/dAL — a derived quantity, not an assumed constant"),
          plotOutput("p_dsda", height = "300px"),
          helpText("At the corneal plane the sensitivity falls by 31% as the eye lengthens, but ",
                   "the vertex conversion to the spectacle plane almost cancels that, leaving only 7%. ",
                   "That is why the '2.7 D/mm' rule of thumb fits only spectacle refraction.")),
        # ------------------------------------------------------------------
        tabPanel("4. Choroid vs sclera",
          br(),
          h4("Two actuators sharing a single measurement"),
          plotOutput("p_chor", height = "400px"),
          hr(),
          h4("The reversible (choroidal) share of the measured one-year effect"),
          tableOutput("t_decomp"),
          helpText("A biometer measures from the cornea to the RPE, so a thickening ",
                   "choroid makes the axial length look as though it had shortened. ",
                   "This component reverses when treatment stops.")),
        # ------------------------------------------------------------------
        tabPanel("5. Atropine PK",
          br(), plotOutput("p_pk", height = "420px"),
          hr(), h4("Systemic exposure"),
          tableOutput("t_pk"),
          helpText("Plasma C_max is about 8 pg/mL at 0.01% and about 765 pg/mL at 1%, ",
                   "consistent with the 300–900 pg/mL range measured after 1% instillation.")),
        # ------------------------------------------------------------------
        tabPanel("6. Dose-response",
          br(),
          h4("Efficacy and side effect have different dose-response curves"),
          plotOutput("p_dr", height = "420px"),
          hr(), DTOutput("t_dr"),
          helpText("The Hill slope of the axial response is 1.25 and the ED50 about 0.05%, ",
                   "effectively the same as the concentration LAMP recommended. ",
                   "The efficacy/mydriasis ratio is worst at 0.01%.")),
        # ------------------------------------------------------------------
        tabPanel("7. Scenario comparison",
          br(),
          checkboxGroupInput("arms", "Treatments to compare (up to 6)",
            choices = c("No treatment", "Atropine 0.01%", "Atropine 0.05%",
                        "Atropine 1%", "DIMS spectacles", "HAL spectacles",
                        "MiSight", "Orthokeratology", "RLRL",
                        "Ortho-K + atropine 0.01%"),
            selected = c("No treatment", "Atropine 0.05%", "DIMS spectacles",
                         "Orthokeratology"), inline = TRUE),
          plotOutput("p_cmp", height = "480px"),
          hr(), DTOutput("t_cmp")),
        # ------------------------------------------------------------------
        tabPanel("8. Rebound",
          br(),
          h4("After stopping: the fast collapse of the choroid plus the slow hypersensitivity of receptor up-regulation"),
          plotOutput("p_reb", height = "440px"),
          hr(), h4("Washout-year progression by mode of stopping"),
          tableOutput("t_reb"),
          helpText("Rebound is the sum of two components. Receptor up-regulation ",
                   "decays with τ = 90 days, so a taper reduces the rebound.")),
        # ------------------------------------------------------------------
        tabPanel("9. Scleral biology",
          br(), plotOutput("p_scl", height = "520px"),
          helpText("Every coupling is a power law, so no variable can ever become ",
                   "negative. The MMP-2:TIMP-2 balance is the single scalar that sets the direction ",
                   "of matrix turnover, and scleral creep is the final common pathway.")),
        # ------------------------------------------------------------------
        tabPanel("10. Clinical endpoints",
          br(),
          fluidRow(
            column(6, h4("Cumulative absolute reduction (CARE)"),
                   plotOutput("p_care", height = "320px")),
            column(6, h4("Lifetime visual impairment risk"),
                   plotOutput("p_vi", height = "320px"))),
          hr(),
          h4("The price of 1 mm is not constant"),
          tableOutput("t_mm"),
          helpText("At an axial length of 23 mm, 1 mm means 1.1%p of lifetime visual ",
                   "impairment risk, but at 28 mm it means 22.8%p — a 20-fold difference. ",
                   "So '% reduction' is the wrong currency."))
      )
    )
  )
)

## =============================================================================
##  SERVER
## =============================================================================
server <- function(input, output, session) {

  days <- reactive(input$years * 365)
  stopd <- reactive(if (input$stopyr >= 10) 1e9 else input$stopyr * 365)

  base_args <- reactive(list(
    age0 = input$age0, ser0 = input$ser0, npar = as.numeric(input$npar),
    grs = input$grs, ethn = as.numeric(input$ethn), outd = input$outd,
    neard = input$neard, fuiris = as.numeric(input$iris),
    adhfix = if (input$adh > 0) input$adh else -1,
    trtstart = input$trtstart * 365, days = days()))

  sim_one <- function(atro = 0, trtpd = 0.15, okon = 0, rlrl = 0,
                      atrostop = 1e9, atrotap = 0) {
    a <- base_args()
    do.call(run_arm, c(a, list(atro = atro, atrostop = atrostop,
                               atrotap = atrotap, trtpd = trtpd,
                               okon = okon, rlrl = rlrl)))
  }

  cur <- reactive(sim_one(atro = as.numeric(input$atro),
                          trtpd = DEVICE[[input$dev]],
                          okon = as.numeric(input$okon),
                          rlrl = as.numeric(input$rlrl),
                          atrostop = stopd(),
                          atrotap = input$tapyr * 365))
  ctl <- reactive(sim_one())

  ## ---- tab 1 ------------------------------------------------------------
  output$t_base <- renderTable({
    al0 <- solve_AL0(input$ser0)
    data.frame(
      Item = c("Baseline refraction SER (D)", "Baseline axial length (mm)", "Corneal radius (mm)",
               "Corneal power (D)", "Lens power (D)", "AL/CR ratio",
               "dSER/dAL (spectacle plane, D/mm)"),
      Value = c(sprintf("%.2f", input$ser0), sprintf("%.3f", al0), "7.80",
             sprintf("%.2f", 0.336 / 0.00780), "22.60",
             sprintf("%.3f", al0 / 7.80),
             sprintf("%.3f", (ser_of(al0 + 0.01) - ser_of(al0)) / 0.01)))
  }, striped = TRUE)

  output$t_risk <- renderTable({
    d <- cur(); e <- tail(d, 1)
    data.frame(
      Item = c("Final axial length (mm)", "Final SER (D)", "Reached high myopia (≤ −6 D)",
               "Reached 26 mm axial length", "Lifetime visual impairment risk",
               "Choroidal thickness (µm)", "Posterior scleral thickness (µm)"),
      Value = c(sprintf("%.3f", e$AL), sprintf("%.2f", e$SERTR),
             ifelse(any(d$HIGHMYO > 0), "Yes", "No"),
             ifelse(any(d$AL26 > 0), "Yes", "No"),
             sprintf("%.1f%%", 100 * lifetime_vi(e$AL)),
             sprintf("%.0f", e$CHT), sprintf("%.0f", e$SCT)))
  }, striped = TRUE)

  output$p_profile <- renderPlot({
    d <- cur()
    d %>% select(yr, AL, ALCR) %>%
      pivot_longer(-yr) %>%
      mutate(name = recode(name, AL = "Measured axial length (mm)",
                           ALCR = "AL/CR ratio")) %>%
      ggplot(aes(yr, value)) +
      geom_line(linewidth = 1.1, colour = "#2f6fb5") +
      geom_hline(data = data.frame(name = "AL/CR ratio", y = 3.0),
                 aes(yintercept = y), linetype = 2, colour = "#b8564a") +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "Follow-up (yr)", y = NULL) + THEME
  })

  ## ---- tab 2 ------------------------------------------------------------
  output$p_traj <- renderPlot({
    bind_rows(mutate(ctl(), arm = "No treatment"),
              mutate(cur(), arm = "This setting")) %>%
      select(yr, arm, AL, SERTR, CHT, CREEPS) %>%
      pivot_longer(c(AL, SERTR, CHT, CREEPS)) %>%
      mutate(name = recode(name,
        AL = "Measured axial length (mm)", SERTR = "Refraction SER (D)",
        CHT = "Choroidal thickness (µm)", CREEPS = "Scleral creep (relative)")) %>%
      ggplot(aes(yr, value, colour = arm)) +
      geom_line(linewidth = 1.1) +
      facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = c("#8a8f98", "#2f6fb5")) +
      labs(x = "Follow-up (yr)", y = NULL) + THEME
  })

  output$t_traj <- renderDT({
    d <- cur(); c0 <- ctl()
    yrs <- seq_len(input$years)
    pick <- function(df, y) df[which.min(abs(df$yr - y)), ]
    tibble(
      `Year` = yrs,
      `Age` = sapply(yrs, function(y) round(pick(d, y)$AGEy, 1)),
      `Axial length (mm)` = sapply(yrs, function(y) round(pick(d, y)$AL, 3)),
      `Annual axial growth (mm/yr)` = sapply(yrs, function(y)
        round(pick(d, y)$AL - pick(d, y - 1)$AL, 3)),
      `SER (D)` = sapply(yrs, function(y) round(pick(d, y)$SERTR, 2)),
      `Untreated axial length (mm)` = sapply(yrs, function(y) round(pick(c0, y)$AL, 3)),
      `CARE (mm)` = sapply(yrs, function(y)
        round(pick(c0, y)$AL - pick(d, y)$AL, 3))
    ) %>% datatable(rownames = FALSE, options = list(dom = "t", pageLength = 12))
  })

  ## ---- tab 3 ------------------------------------------------------------
  output$p_optics <- renderPlot({
    d <- cur()
    al0 <- solve_AL0(input$ser0)
    ## contribution of each element to the refractive change, in dioptres
    d %>% transmute(
      yr,
      `Axial length contribution (element 1)` = ser_of(AL, 3.60, 3.45, 22.60, 7.80) -
        ser_of(al0, 3.60, 3.45, 22.60, 7.80),
      `Lens contribution (element 3)` = ser_of(al0, 3.60, 3.45, PLENS, 7.80) -
        ser_of(al0, 3.60, 3.45, 22.60, 7.80),
      `Net SER change` = SERTR - input$ser0) %>%
      pivot_longer(-yr) %>%
      ggplot(aes(yr, value, colour = name)) +
      geom_hline(yintercept = 0, colour = "grey70") +
      geom_line(linewidth = 1.1) +
      scale_colour_manual(values = c("#2f6fb5", "#b8564a", "#1f1f1f")) +
      labs(x = "Follow-up (yr)", y = "Refractive change (D)") + THEME
  })

  output$p_dsda <- renderPlot({
    al <- seq(22.5, 31, 0.1)
    spec <- (ser_of(al + 0.01) - ser_of(al)) / 0.01
    ## the corneal-plane sensitivity, obtained by undoing the vertex conversion
    fc <- function(s) s / (1 - 0.012 * s)
    corn <- (fc(ser_of(al + 0.01)) - fc(ser_of(al))) / 0.01
    tibble(al, `Spectacle plane` = spec,
           `Corneal plane` = corn) %>%
      pivot_longer(-al) %>%
      ggplot(aes(al, value, colour = name)) +
      geom_line(linewidth = 1.1) +
      scale_colour_manual(values = c("#b8564a", "#2f6fb5")) +
      labs(x = "Axial length (mm)", y = "dSER/dAL (D/mm)") + THEME
  })

  ## ---- tab 4 ------------------------------------------------------------
  output$p_chor <- renderPlot({
    d <- cur()
    al0 <- solve_AL0(input$ser0)
    d %>% transmute(yr,
      `Axial length, sclera only (permanent)` = al0 + AXCUM,
      `Measured axial length (what a biometer reads)` = AL,
      `Choroidal thickness (µm, not on its own axis)` = CHT / 1000 + al0) %>%
      pivot_longer(-yr) %>%
      ggplot(aes(yr, value, colour = name)) +
      geom_line(linewidth = 1.1) +
      scale_colour_manual(values = c("#c9955f", "#1f1f1f", "#2f6fb5")) +
      labs(x = "Follow-up (yr)", y = "mm") + THEME
  })

  output$t_decomp <- renderTable({
    a <- base_args(); a$days <- 365
    c0 <- do.call(run_arm, c(a, list(atro = 0, atrostop = 1e9, atrotap = 0,
                                     trtpd = 0.15, okon = 0, rlrl = 0)))
    al0 <- solve_AL0(input$ser0)
    arms <- list(`Atropine 0.01%` = list(atro = 0.01),
                 `Atropine 0.05%` = list(atro = 0.05),
                 `Atropine 1%`    = list(atro = 1.00),
                 `DIMS spectacles` = list(trtpd = -1.20),
                 `Orthokeratology` = list(okon = 1),
                 `RLRL`           = list(rlrl = 1))
    rows <- lapply(names(arms), function(nm) {
      z <- arms[[nm]]
      d <- do.call(run_arm, c(a, list(
        atro = z$atro %||% 0, atrostop = 1e9, atrotap = 0,
        trtpd = z$trtpd %||% 0.15, okon = z$okon %||% 0, rlrl = z$rlrl %||% 0)))
      dm <- (tail(c0, 1)$AL - al0) - (tail(d, 1)$AL - al0)
      ds <- tail(c0, 1)$AXCUM - tail(d, 1)$AXCUM
      data.frame(Treatment = nm,
                 `Measured reduction (mm)` = round(dm, 4),
                 `Scleral reduction (mm)` = round(ds, 4),
                 `Choroidal component (mm)` = round(dm - ds, 4),
                 `Reversible fraction` = sprintf("%.1f%%", 100 * (dm - ds) / dm),
                 check.names = FALSE)
    })
    bind_rows(rows)
  }, striped = TRUE)

  ## ---- tab 5 ------------------------------------------------------------
  output$p_pk <- renderPlot({
    if (as.numeric(input$atro) == 0)
      return(ggplot() + annotate("text", 0, 0, label = "Select an atropine concentration") +
               theme_void())
    a <- base_args(); a$days <- 14
    d <- do.call(run_arm, c(a, list(atro = as.numeric(input$atro),
                                    atrostop = 1e9, atrotap = 0, trtpd = 0.15,
                                    okon = 0, rlrl = 0)))
    d %>% select(time, CPLA_PG, CCHOR_NM, CSCL_NM) %>%
      pivot_longer(-time) %>%
      mutate(name = recode(name, CPLA_PG = "Plasma (pg/mL)",
                           CCHOR_NM = "Choroid (nM)", CSCL_NM = "Sclera (nM)")) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 1) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      labs(x = "Time (day)", y = NULL) + THEME
  })

  output$t_pk <- renderTable({
    a <- base_args(); a$days <- 30
    rows <- lapply(c(0.01, 0.025, 0.05, 0.1, 0.5, 1.0), function(pc) {
      d <- do.call(run_arm, c(a, list(atro = pc, atrostop = 1e9, atrotap = 0,
                                      trtpd = 0.15, okon = 0, rlrl = 0)))
      data.frame(`Concentration (%)` = pc,
                 `Plasma C_max (pg/mL)` = round(max(d$CPLA_PG), 1),
                 `Pupil (mm)` = round(tail(d, 1)$PUPD, 2),
                 `Accommodation (D)` = round(tail(d, 1)$AAMP, 2),
                 `Drive inhibition EATR` = round(tail(d, 1)$EATRo, 3),
                 check.names = FALSE)
    })
    bind_rows(rows)
  }, striped = TRUE)

  ## ---- tab 6 ------------------------------------------------------------
  dr_tab <- reactive({
    a <- base_args(); a$days <- 365
    al0 <- solve_AL0(input$ser0)
    c0 <- do.call(run_arm, c(a, list(atro = 0, atrostop = 1e9, atrotap = 0,
                                     trtpd = 0.15, okon = 0, rlrl = 0)))
    base <- tail(c0, 1)$AL - al0
    doses <- c(0.005, 0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 1.0)
    bind_rows(lapply(doses, function(pc) {
      d <- tail(do.call(run_arm, c(a, list(atro = pc, atrostop = 1e9,
        atrotap = 0, trtpd = 0.15, okon = 0, rlrl = 0))), 1)
      eff <- 1 - (d$AL - al0) / base
      myd <- (d$PUPD - 4.60) / 3.80
      tibble(dose = pc, efficacy = eff, mydriasis = myd,
             accom_loss = 1 - d$AAMP / 13.40,
             ratio = eff / pmax(myd, 1e-6), adherence = d$ADH)
    }))
  })

  output$p_dr <- renderPlot({
    dr_tab() %>%
      select(dose, `Axial efficacy` = efficacy,
             `Mydriasis` = mydriasis,
             `Accommodation loss` = accom_loss,
             `Efficacy/mydriasis ratio` = ratio) %>%
      pivot_longer(-dose) %>%
      ggplot(aes(dose, value, colour = name)) +
      geom_line(linewidth = 1.1) + geom_point(size = 1.8) +
      scale_x_log10(breaks = c(0.005, 0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 1)) +
      geom_vline(xintercept = 0.05, linetype = 2, colour = "grey60") +
      labs(x = "Atropine concentration (% w/v, log)", y = NULL) + THEME
  })

  output$t_dr <- renderDT({
    dr_tab() %>%
      transmute(`Concentration (%)` = dose,
                `Axial efficacy` = round(efficacy, 3),
                `Mydriasis fraction` = round(mydriasis, 3),
                `Accommodation loss fraction` = round(accom_loss, 3),
                `Efficacy/mydriasis ratio` = round(ratio, 2),
                `Adherence` = round(adherence, 2)) %>%
      datatable(rownames = FALSE, options = list(dom = "t"))
  })

  ## ---- tab 7 ------------------------------------------------------------
  ARMS <- list(
    `No treatment`             = list(),
    `Atropine 0.01%`           = list(atro = 0.01),
    `Atropine 0.05%`           = list(atro = 0.05),
    `Atropine 1%`              = list(atro = 1.00),
    `DIMS spectacles`          = list(trtpd = -1.20),
    `HAL spectacles`           = list(trtpd = -2.20),
    `MiSight`                  = list(trtpd = -0.90),
    `Orthokeratology`          = list(okon = 1, trtpd = 0),
    `RLRL`                     = list(rlrl = 1),
    `Ortho-K + atropine 0.01%` = list(okon = 1, trtpd = 0, atro = 0.01))

  cmp <- reactive({
    sel <- head(input$arms, 6)
    bind_rows(lapply(sel, function(nm) {
      z <- ARMS[[nm]]
      mutate(sim_one(atro = z$atro %||% 0, trtpd = z$trtpd %||% 0.15,
                     okon = z$okon %||% 0, rlrl = z$rlrl %||% 0), arm = nm)
    }))
  })

  output$p_cmp <- renderPlot({
    cmp() %>% select(yr, arm, AL, SERTR, CHT, VIRISK) %>%
      pivot_longer(c(AL, SERTR, CHT, VIRISK)) %>%
      mutate(name = recode(name, AL = "Measured axial length (mm)",
                           SERTR = "Refraction SER (D)",
                           CHT = "Choroidal thickness (µm)",
                           VIRISK = "Lifetime visual impairment risk")) %>%
      ggplot(aes(yr, value, colour = arm)) + geom_line(linewidth = 1) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "Follow-up (yr)", y = NULL) + THEME
  })

  output$t_cmp <- renderDT({
    al0 <- solve_AL0(input$ser0)
    cmp() %>% group_by(arm) %>% slice_tail(n = 1) %>% ungroup() %>%
      transmute(`Treatment` = arm,
                `Axial length (mm)` = round(AL, 3),
                `Total elongation (mm)` = round(AL - al0, 3),
                `SER (D)` = round(SERTR, 2),
                `Lifetime visual impairment risk` = sprintf("%.1f%%", 100 * VIRISK),
                `Pupil (mm)` = round(PUPD, 2),
                `Accommodation (D)` = round(AAMP, 2)) %>%
      arrange(`Axial length (mm)`) %>%
      datatable(rownames = FALSE, options = list(dom = "t"))
  })

  ## ---- tab 8 ------------------------------------------------------------
  output$p_reb <- renderPlot({
    a <- base_args(); a$days <- 3 * 365
    mk <- function(nm, pc, tap) mutate(do.call(run_arm, c(a, list(
      atro = pc, atrostop = 2 * 365, atrotap = tap, trtpd = 0.15,
      okon = 0, rlrl = 0))), arm = nm)
    bind_rows(
      mutate(do.call(run_arm, c(a, list(atro = 0, atrostop = 1e9, atrotap = 0,
        trtpd = 0.15, okon = 0, rlrl = 0))), arm = "No treatment"),
      mk("0.01% → stop", 0.01, 0),
      mk("0.5% → abrupt stop", 0.50, 0),
      mk("0.5% → 1-year taper", 0.50, 365)) %>%
      select(yr, arm, AL, SERTR, CHT) %>%
      pivot_longer(c(AL, SERTR, CHT)) %>%
      mutate(name = recode(name, AL = "Measured axial length (mm)",
                           SERTR = "Refraction SER (D)",
                           CHT = "Choroidal thickness (µm)")) %>%
      ggplot(aes(yr, value, colour = arm)) +
      geom_vline(xintercept = 2, linetype = 2, colour = "grey60") +
      geom_line(linewidth = 1) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      labs(x = "Follow-up (yr)", y = NULL) + THEME
  })

  output$t_reb <- renderTable({
    a <- base_args(); a$days <- 3 * 365
    pick <- function(df, y) df[which.min(abs(df$yr - y)), ]
    rows <- lapply(list(c(0.01, 0), c(0.10, 0), c(0.50, 0),
                        c(0.50, 90), c(0.50, 180), c(0.50, 365)),
      function(z) {
        d <- do.call(run_arm, c(a, list(atro = z[1], atrostop = 2 * 365,
          atrotap = z[2], trtpd = 0.15, okon = 0, rlrl = 0)))
        data.frame(`Concentration (%)` = z[1], `Taper (day)` = z[2],
                   `Washout SER change (D)` =
                     round(pick(d, 3)$SERTR - pick(d, 2)$SERTR, 3),
                   `Washout axial length (mm)` =
                     round(pick(d, 3)$AL - pick(d, 2)$AL, 3),
                   check.names = FALSE)
      })
    bind_rows(rows)
  }, striped = TRUE)

  ## ---- tab 9 ------------------------------------------------------------
  output$p_scl <- renderPlot({
    d <- cur()
    d %>% select(yr, MTRo, CREEPS, SCT, CHOX, PSTAPH, LACQ) %>%
      pivot_longer(-yr) %>%
      mutate(name = recode(name,
        MTRo = "MMP-2 : TIMP-2 balance", CREEPS = "Scleral creep (relative)",
        SCT = "Posterior scleral thickness (µm)", CHOX = "Choroidal hypoxia index",
        PSTAPH = "Posterior staphyloma index", LACQ = "Lacquer crack / atrophy index")) %>%
      ggplot(aes(yr, value)) + geom_line(linewidth = 1, colour = "#a05a95") +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "Follow-up (yr)", y = NULL) + THEME
  })

  ## ---- tab 10 -----------------------------------------------------------
  output$p_care <- renderPlot({
    c0 <- ctl(); d <- cur()
    pick <- function(df, y) df[which.min(abs(df$yr - y)), ]
    yrs <- seq_len(input$years)
    tibble(yr = yrs,
           CARE = sapply(yrs, function(y) pick(c0, y)$AL - pick(d, y)$AL)) %>%
      mutate(increment = CARE - lag(CARE, default = 0)) %>%
      pivot_longer(-yr) %>%
      mutate(name = recode(name, CARE = "CARE, cumulative (mm)",
                           increment = "Yearly increment (mm)")) %>%
      ggplot(aes(yr, value, fill = name)) +
      geom_col(position = "dodge") +
      scale_fill_manual(values = c("#2f6fb5", "#c9955f")) +
      labs(x = "Year", y = "mm") + THEME
  })

  output$p_vi <- renderPlot({
    bind_rows(mutate(ctl(), arm = "No treatment"), mutate(cur(), arm = "This setting")) %>%
      ggplot(aes(yr, 100 * VIRISK, colour = arm)) +
      geom_line(linewidth = 1.1) +
      scale_colour_manual(values = c("#8a8f98", "#2f6fb5")) +
      labs(x = "Follow-up (yr)", y = "Lifetime visual impairment risk (%)") + THEME
  })

  output$t_mm <- renderTable({
    al <- 23:30
    data.frame(`Axial length (mm)` = al,
               `Lifetime risk` = sprintf("%.1f%%", 100 * lifetime_vi(al)),
               `After +1 mm` = sprintf("%.1f%%", 100 * lifetime_vi(al + 1)),
               `Increase (%p)` = round(100 * (lifetime_vi(al + 1) -
                                            lifetime_vi(al)), 2),
               check.names = FALSE)
  }, striped = TRUE)
}

shinyApp(ui, server)
