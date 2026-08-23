## =====================================================================
##  Hereditary Spherocytosis QSP — Shiny dashboard
##  Hereditary Spherocytosis QSP Dashboard (14 tabs)
##
##  Front end for hsph_mrgsolve_model.R.  Like the model file, this has NOT
##  been executed (no R toolchain in the build environment); it mirrors the
##  Python reference implementation, which has.
##
##  run:  shiny::runApp("hsph_shiny_app.R")
## =====================================================================
library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)

MOD <- mread_cache("hsph_mrgsolve_model", ".")

## ---------------------------------------------------------------- helpers
## the geometry kernel, in R, so the geometry tab does not need the ODEs
vsph  <- function(A) A^1.5 / (6 * sqrt(pi))
dcrit <- function(A, V) {
  s  <- pmin(pmax(V / vsph(A), 1e-9), 1)
  2 * sqrt(A / pi) * cos(acos(-s) / 3 - 2 * pi / 3)
}
ofpoint <- function(A, V, mchc) {          # 50% lysis, %NaCl
  vs <- 0.75 * mchc * V / 100 + 0.04 * V
  vm <- vsph(A)
  ifelse(vm > vs, 290 * (V - vs) / (vm - vs) * (0.9 / 308), 0.9)
}

GENO <- list(
  "Normal"                              = c(fdef = 0.00, f_b3ves = 1.00),
  "Carrier (HS trait)"                          = c(fdef = 0.12, f_b3ves = 1.00),
  "Mild HS"                             = c(fdef = 0.20, f_b3ves = 1.00),
  "Moderate ANK1/SPTB"                = c(fdef = 0.30, f_b3ves = 1.00),
  "Moderate SLC4A1 band 3"            = c(fdef = 0.30, f_b3ves = 0.15),
  "Severe SPTA1 (recessive)"             = c(fdef = 0.45, f_b3ves = 1.00)
)

run_one <- function(par, end = 1400, delta = 2, burn = 1400) {
  m <- MOD %>% param(par)
  m %>% mrgsim(end = end, delta = delta, hmax = 2) %>% as_tibble()
}

## run to steady state, then apply an intervention at t = 0 of a second run
run_intervention <- function(base_par, iv_par, pre = 1400, post = 1460,
                             delta = 2) {
  s0 <- MOD %>% param(base_par) %>% mrgsim(end = pre, delta = pre) %>%
    as_tibble() %>% tail(1)
  ini <- MOD %>% param(base_par) %>%
    mrgsim(end = pre, delta = pre, obsonly = FALSE)
  MOD %>% param(c(base_par, iv_par)) %>%
    init(as.list(ini@data[nrow(ini@data), names(init(MOD))])) %>%
    mrgsim(end = post, delta = delta, hmax = 2) %>% as_tibble()
}

thm <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92"),
        legend.position = "bottom")

## ------------------------------------------------------------------- UI
ui <- fluidPage(
  titlePanel("Hereditary Spherocytosis QSP Model"),
  tags$p(tags$em(paste(
    "The red blood cell is two numbers: membrane area S and volume V. These two numbers determine the minimum cylinder diameter",
    "D_c, and D_c is the only thing the splenic sinusoidal wall can see.",
    "Area loss → D_c rises → prolonged splenic cord retention → area loss —",
    "this feedback loop is the disease."))),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectInput("geno", "Genotype / Severity",
                  names(GENO), selected = "Moderate ANK1/SPTB"),
      sliderInput("fdef", "Vertical linkage defect fdef = 1 − spectrin content",
                  0, 0.6, 0.30, 0.01),
      sliderInput("fb3", "Vesicle band 3 content f_b3ves",
                  0, 1, 1.00, 0.05,
                  post = "  (1 = spectrin/ankyrin-type, 0.15 = band 3-type)"),
      hr(),
      radioButtons("spleen", "Spleen",
                   c("Normal (intact)" = "1",
                     "Partial splenectomy, 20% remaining" = "0.2",
                     "Partial splenectomy, 10% remaining" = "0.1",
                     "Total splenectomy" = "0"), selected = "1"),
      checkboxInput("regrow", "Regrowth of residual spleen after partial splenectomy", TRUE),
      hr(),
      sliderInput("mita", "Mitapivat (mg BID)", 0, 100, 0, 5),
      selectInput("ugt", "UGT1A1",
                  c("Normal (wild type)" = "1",
                    "Gilbert syndrome UGT1A1*28/*28" = "0.287")),
      checkboxInput("folate", "Folate supplementation (1–5 mg/day)", TRUE),
      hr(),
      checkboxInput("parvo", "Parvovirus B19 aplastic crisis", FALSE),
      numericInput("parvo_t", "  B19 onset timing (day)", 30, 0, 300),
      hr(),
      checkboxInput("tx", "Scheduled transfusion (2 units q4w)", FALSE),
      checkboxInput("chel", "Iron chelation (deferasirox)", FALSE),
      hr(),
      sliderInput("tend", "Simulation duration (days)", 200, 3650, 730, 30),
      actionButton("go", "Run", class = "btn-primary")
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("1 · Patient Profile",
                 h4("Steady State Summary"),
                 tableOutput("profile"),
                 h4("Severity Spectrum — determined by fdef alone"),
                 plotOutput("spectrum", height = 340)),
        tabPanel("2 · Red Cell Geometry",
                 h4("D_c = 2√(S/π)·cos(arccos(−V/V_sph)/3 − 2π/3)"),
                 tags$p("A closed-form equation with nothing fitted."),
                 plotOutput("geomap", height = 380),
                 h4("The (S, V) Trajectory with Cell Age"),
                 plotOutput("geotraj", height = 300)),
        tabPanel("3 · Splenic Filter",
                 h4("A Filter That Sees Only D_c"),
                 plotOutput("filter", height = 380),
                 tableOutput("filtertab")),
        tabPanel("4 · Blood Tests (CBC)",
                 plotOutput("cbc", height = 620)),
        tabPanel("5 · Haemolysis Markers",
                 plotOutput("hemol", height = 560)),
        tabPanel("6 · Clearance Mechanism Decomposition",
                 tags$p(paste("Normal red cells die by being labelled, HS",
                              "red cells die because of their shape. Same organ,",
                              "different mechanism.")),
                 plotOutput("clearance", height = 420),
                 h4("Where Membrane Loss Happens"),
                 tableOutput("stripping")),
        tabPanel("7 · Splenectomy",
                 h4("It Fixes the Anaemia, Not the Cell"),
                 plotOutput("splx", height = 430),
                 tableOutput("splxtab")),
        tabPanel("8 · Genotype × Splenectomy",
                 tags$p(paste("In this model, only one parameter separates the two genotypes:",
                              "whether band 3 is carried away in the vesicle.")),
                 plotOutput("genoplot", height = 400),
                 tableOutput("genotab")),
        tabPanel("9 · Diagnostic Tests",
                 h4("Osmotic Fragility Curve"),
                 plotOutput("offcurve", height = 340),
                 h4("EMA Binding · Osmotic Gradient Ektacytometry"),
                 plotOutput("emaplot", height = 300)),
        tabPanel("10 · Crises",
                 h4("Parvovirus B19: ΔHb ≈ Hb × arrest days / lifespan"),
                 plotOutput("crisis", height = 420),
                 tableOutput("crisistab")),
        tabPanel("11 · Drug (mitapivat)",
                 plotOutput("drug", height = 560),
                 tableOutput("drugtab")),
        tabPanel("12 · Gallstones · Iron",
                 plotOutput("stone", height = 480),
                 tableOutput("stonetab")),
        tabPanel("13 · Scenario Comparison",
                 plotOutput("scen", height = 620),
                 tableOutput("scentab")),
        tabPanel("14 · Virtual Patient Cohort",
                 numericInput("npop", "Number of patients", 100, 20, 500, 20),
                 plotOutput("pop", height = 480),
                 tableOutput("poptab"))
      )
    )
  )
)

## --------------------------------------------------------------- server
server <- function(input, output, session) {

  observeEvent(input$geno, {
    g <- GENO[[input$geno]]
    updateSliderInput(session, "fdef", value = unname(g["fdef"]))
    updateSliderInput(session, "fb3",  value = unname(g["f_b3ves"]))
  })

  par_now <- reactive({
    sf <- as.numeric(input$spleen)
    c(fdef      = input$fdef,
      f_b3ves   = input$fb3,
      spl_frac  = sf,
      k_regrow  = if (input$regrow && sf > 0 && sf < 1) 1 / 700 else 0,
      dose_m    = input$mita,
      ugt_f     = as.numeric(input$ugt),
      fol_ok    = if (input$folate) 1.0 else 0.35,
      parvo_t   = if (input$parvo) input$parvo_t else -1,
      tx_start  = if (input$tx) 0 else -1,
      k_chel    = if (input$chel) 0.0016 else 0)
  })

  sim <- eventReactive(input$go, {
    run_one(par_now(), end = input$tend)
  }, ignoreNULL = FALSE)

  ss <- reactive(tail(sim(), 1))

  ## ---- 1 profile -------------------------------------------------------
  output$profile <- renderTable({
    s <- ss()
    tibble(
      Metric = c("Haemoglobin Hb (g/dL)", "RBC (10¹²/L)", "MCV (fL)",
               "MCHC (g/dL)", "Reticulocytes (%)", "RBC lifespan (days)",
               "Total bilirubin (mg/dL)", "Spleen volume (mL)",
               "Membrane area S (µm²)", "Minimum cylinder diameter D_c (µm)",
               "EMA binding ratio", "50% lysis (%NaCl)",
               "RBC-bound IgG (molecules/cell)"),
      Value = c(s$Hb, s$RBC, s$MCV, s$MCHC, s$RETpct, s$LIFE, s$TBIL,
             s$SPLEEN, s$AREA, s$DCRIT, s$EMA, s$MCF, s$IGGC)
    )
  }, digits = 3)

  output$spectrum <- renderPlot({
    fs <- seq(0, 0.5, by = 0.05)
    d <- lapply(fs, function(f) {
      p <- par_now(); p["fdef"] <- f
      tail(run_one(p, end = 1200, delta = 1200), 1) %>% mutate(fdef = f)
    }) %>% bind_rows()
    d %>% select(fdef, Hb, RETpct, LIFE, TBIL, AREA, DCRIT) %>%
      pivot_longer(-fdef) %>%
      ggplot(aes(fdef, value)) + geom_line(linewidth = 1, colour = "#c62828") +
      geom_point() + facet_wrap(~name, scales = "free_y") +
      labs(x = "fdef = 1 − spectrin content", y = NULL) + thm
  })

  ## ---- 2 geometry ------------------------------------------------------
  output$geomap <- renderPlot({
    g <- expand.grid(S = seq(90, 150, 1), V = seq(70, 105, 1)) %>%
      mutate(Dc = dcrit(S, V)) %>% filter(V < vsph(S))
    pts <- tibble(
      lab = c("Normal young", "Normal old", "Mild HS", "Moderate HS", "Severe HS"),
      S = c(140, 125.3, 128, 118, 106),
      V = c(94, 86.1, 89, 86, 80))
    ggplot(g, aes(S, V)) +
      geom_raster(aes(fill = Dc)) +
      geom_contour(aes(z = Dc), colour = "white", breaks = seq(2.6, 4.2, .2)) +
      scale_fill_viridis_c(option = "magma", direction = -1) +
      geom_point(data = pts, size = 3, colour = "cyan") +
      geom_text(data = pts, aes(label = lab), nudge_y = 2.4,
                colour = "cyan", size = 3.4) +
      labs(x = "Membrane area S (µm²)", y = "Volume V (µm³)",
           fill = "D_c (µm)",
           title = "Losing area without losing volume makes D_c rise sharply") + thm
  })

  output$geotraj <- renderPlot({
    s <- sim()
    s %>% select(time, AREA, VOLc, DCRIT, SPHER) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "Time (days)", y = NULL) + thm
  })

  ## ---- 3 splenic filter ------------------------------------------------
  output$filter <- renderPlot({
    p  <- as.list(param(MOD))
    Dc <- seq(2.5, 4.2, 0.01)
    ps <- 1 / (1 + exp(-(Dc - p$D50) / p$wD))
    tc <- p$tau0 * exp((Dc - p$Dc_ref) / p$w_esc)
    Rc <- pmin(p$R_MAX, p$f_pass0 * ps * tc)
    tibble(Dc,
           `p_slow (fraction entering the cord)` = ps,
           `τ_cord (days)` = tc,
           `R_cord (fraction of time in cord)` = Rc,
           `Membrane loss multiplier (1+41·R)` = 1 + (p$cordamp - 1) * Rc) %>%
      pivot_longer(-Dc) %>%
      ggplot(aes(Dc, value)) + geom_line(linewidth = 1) +
      geom_vline(xintercept = c(2.88, 3.27), linetype = 2,
                 colour = c("#2e7d32", "#c62828")) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "Minimum cylinder diameter D_c (µm)", y = NULL,
           caption = "Green = normal red cell, red = moderate HS") + thm
  })

  output$filtertab <- renderTable({
    s <- ss()
    tibble(`Cord retention fraction R_cord` = s$RCORD,
           `Membrane loss occurring in the spleen (%)` =
             100 * 41 * s$RCORD / (1 + 41 * s$RCORD),
           `Spleen volume (mL)` = s$SPLEEN)
  }, digits = 4)

  ## ---- 4 CBC -----------------------------------------------------------
  output$cbc <- renderPlot({
    sim() %>% select(time, Hb, Hct, RBC, MCV, MCHC, RETpct, RETabs, LIFE) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.9, colour = "#1565c0") +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "Time (days)", y = NULL) + thm
  })

  ## ---- 5 haemolysis ----------------------------------------------------
  output$hemol <- renderPlot({
    sim() %>% select(time, TBIL, BILU, BILC, LDH, HPT, BRPROD, SPLEEN) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.9, colour = "#6a1b9a") +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "Time (days)", y = NULL) + thm
  })

  ## ---- 6 clearance decomposition ---------------------------------------
  output$clearance <- renderPlot({
    sim() %>% select(time, FGEOM, FOPS, FSEN, FLIV, FLYS) %>%
      pivot_longer(-time, names_to = "mech") %>%
      mutate(mech = recode(mech,
                           FGEOM = "Geometric (splenic cord)", FOPS = "Opsonic (IgG)",
                           FSEN = "Senescence marking", FLIV = "Liver (geometric)",
                           FLYS = "Intravascular haemolysis")) %>%
      ggplot(aes(time, value, fill = mech)) +
      geom_area() + scale_fill_brewer(palette = "Set2") +
      labs(x = "Time (days)", y = "Fraction destroyed", fill = NULL) + thm
  })

  output$stripping <- renderTable({
    lapply(names(GENO), function(g) {
      p <- par_now(); p["fdef"] <- GENO[[g]]["fdef"]
      p["f_b3ves"] <- GENO[[g]]["f_b3ves"]
      s <- tail(run_one(p, end = 1200, delta = 1200), 1)
      tibble(Genotype = g, R_cord = s$RCORD,
             `Membrane loss in the spleen (%)` = 100 * 41 * s$RCORD /
               (1 + 41 * s$RCORD),
             `Geometric clearance` = s$FGEOM, `Opsonic clearance` = s$FOPS,
             `Senescence clearance` = s$FSEN)
    }) %>% bind_rows()
  }, digits = 3)

  ## ---- 7 splenectomy ---------------------------------------------------
  output$splx <- renderPlot({
    base <- par_now(); base["spl_frac"] <- 1
    arms <- list(`No surgery` = c(spl_frac = 1),
                 `Partial splenectomy 20%` = c(spl_frac = 0.2, k_regrow = 1 / 700),
                 `Partial splenectomy 10%` = c(spl_frac = 0.1, k_regrow = 1 / 700),
                 `Total splenectomy` = c(spl_frac = 0))
    lapply(names(arms), function(a)
      run_intervention(base, arms[[a]], post = 1825, delta = 5) %>%
        mutate(arm = a)) %>% bind_rows() %>%
      select(time, arm, Hb, RETpct, TBIL, SPLEEN, AREA, LIFE) %>%
      pivot_longer(c(-time, -arm)) %>%
      ggplot(aes(time / 365, value, colour = arm)) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "Time since surgery (years)", y = NULL, colour = NULL) + thm
  })

  output$splxtab <- renderTable({
    base <- par_now(); base["spl_frac"] <- 1
    pre  <- tail(run_one(base, end = 1400, delta = 1400), 1)
    post <- tail(run_intervention(base, c(spl_frac = 0), post = 1400,
                                  delta = 1400), 1)
    tibble(Metric = c("Hb", "Reticulocytes %", "Lifespan (days)", "Total bilirubin",
                    "Membrane area S", "EMA binding ratio", "50% lysis %NaCl", "MCHC"),
           PreOp = c(pre$Hb, pre$RETpct, pre$LIFE, pre$TBIL, pre$AREA,
                      pre$EMA, pre$MCF, pre$MCHC),
           PostOp = c(post$Hb, post$RETpct, post$LIFE, post$TBIL, post$AREA,
                      post$EMA, post$MCF, post$MCHC))
  }, digits = 3)

  ## ---- 8 genotype x splenectomy ----------------------------------------
  output$genoplot <- renderPlot({
    lapply(c(1.0, 0.15), function(fb) {
      p <- par_now(); p["f_b3ves"] <- fb; p["spl_frac"] <- 1
      lapply(c(1, 0), function(sf) {
        q <- p; q["spl_frac"] <- sf
        tail(run_one(q, end = 1400, delta = 1400), 1) %>%
          mutate(genotype = ifelse(fb > 0.5, "Spectrin/ankyrin", "Band 3"),
                 spleen = ifelse(sf > 0.5, "Spleen present", "Splenectomy"))
      }) %>% bind_rows()
    }) %>% bind_rows() %>%
      select(genotype, spleen, Hb, LIFE, IGGC, FOPS, FGEOM) %>%
      pivot_longer(c(-genotype, -spleen)) %>%
      ggplot(aes(spleen, value, fill = genotype)) +
      geom_col(position = "dodge") +
      facet_wrap(~name, scales = "free_y") +
      scale_fill_manual(values = c("#1565c0", "#c62828")) +
      labs(x = NULL, y = NULL, fill = NULL) + thm
  })

  output$genotab <- renderTable({
    lapply(c(1.0, 0.15), function(fb) {
      p <- par_now(); p["f_b3ves"] <- fb
      p["spl_frac"] <- 1; a <- tail(run_one(p, end = 1400, delta = 1400), 1)
      p["spl_frac"] <- 0; b <- tail(run_one(p, end = 1400, delta = 1400), 1)
      tibble(Genotype = ifelse(fb > .5, "Spectrin/ankyrin", "Band 3"),
             `IgG/cell (spleen present)` = a$IGGC,
             `IgG/cell (post-splenectomy)` = b$IGGC,
             `Hb before` = a$Hb, `Hb after` = b$Hb, `ΔHb` = b$Hb - a$Hb)
    }) %>% bind_rows()
  }, digits = 2)

  ## ---- 9 diagnostics ---------------------------------------------------
  output$offcurve <- renderPlot({
    s <- ss()
    nacl <- seq(0.1, 0.85, 0.005)
    cells <- tibble(
      grp = c("Normal", "Patient"),
      S = c(132, s$AREA), V = c(89, s$VOLc),
      MCHC = c(33.7, s$MCHC))
    d <- lapply(seq_len(nrow(cells)), function(i) {
      m <- ofpoint(cells$S[i], cells$V[i], cells$MCHC[i])
      tibble(grp = cells$grp[i], nacl,
             lysis = 100 / (1 + exp((nacl - m) / 0.045)))
    }) %>% bind_rows()
    ggplot(d, aes(nacl, lysis, colour = grp)) + geom_line(linewidth = 1.1) +
      scale_x_reverse() +
      geom_vline(xintercept = c(0.43, 0.60), linetype = 3) +
      labs(x = "NaCl (%)", y = "Lysis (%)", colour = NULL,
           title = "The curve's position is determined by V_sph = S^1.5/(6√π) alone") + thm
  })

  output$emaplot <- renderPlot({
    lapply(names(GENO), function(g) {
      p <- par_now(); p["fdef"] <- GENO[[g]]["fdef"]
      p["f_b3ves"] <- GENO[[g]]["f_b3ves"]
      tail(run_one(p, end = 1200, delta = 1200), 1) %>% mutate(grp = g)
    }) %>% bind_rows() %>%
      mutate(`EMA reduction (%)` = 100 * (1 - EMA)) %>%
      select(grp, `EMA reduction (%)`, MCF, MCHC) %>%
      pivot_longer(-grp) %>%
      ggplot(aes(grp, value)) + geom_col(fill = "#37474f") +
      facet_wrap(~name, scales = "free_y") +
      geom_hline(data = tibble(name = "EMA reduction (%)", y = 16),
                 aes(yintercept = y), linetype = 2, colour = "red") +
      coord_flip() + labs(x = NULL, y = NULL) + thm
  })

  ## ---- 10 crises -------------------------------------------------------
  output$crisis <- renderPlot({
    lapply(names(GENO)[c(1, 3, 4, 6)], function(g) {
      p <- par_now(); p["fdef"] <- GENO[[g]]["fdef"]
      p["f_b3ves"] <- GENO[[g]]["f_b3ves"]; p["parvo_t"] <- -1
      run_intervention(p, c(parvo_t = 10, parvo_dur = 8),
                       post = 90, delta = 0.5) %>% mutate(grp = g)
    }) %>% bind_rows() %>%
      select(time, grp, Hb, RETpct, TBIL) %>% pivot_longer(c(-time, -grp)) %>%
      ggplot(aes(time, value, colour = grp)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      labs(x = "Days after infection", y = NULL, colour = NULL) + thm
  })

  output$crisistab <- renderTable({
    lapply(names(GENO)[c(1, 3, 4, 6)], function(g) {
      p <- par_now(); p["fdef"] <- GENO[[g]]["fdef"]
      p["f_b3ves"] <- GENO[[g]]["f_b3ves"]; p["parvo_t"] <- -1
      b <- tail(run_one(p, end = 1400, delta = 1400), 1)
      d <- run_intervention(p, c(parvo_t = 10, parvo_dur = 8),
                            post = 90, delta = 0.5)
      tibble(Genotype = g, `Baseline Hb` = b$Hb, `Lifespan (days)` = b$LIFE,
             `Predicted ΔHb = Hb·8/lifespan` = b$Hb * 8 / b$LIFE,
             `Model lowest Hb` = min(d$Hb),
             `Actual ΔHb` = b$Hb - min(d$Hb))
    }) %>% bind_rows()
  }, digits = 2)

  ## ---- 11 drug ---------------------------------------------------------
  output$drug <- renderPlot({
    lapply(c(0, 5, 20, 50, 100), function(d) {
      p <- par_now(); p["dose_m"] <- d
      tail(run_one(p, end = 900, delta = 900), 1) %>% mutate(dose = d)
    }) %>% bind_rows() %>%
      select(dose, MITA, Hb, RETpct, TBIL, AREA, LIFE, DCRIT) %>%
      pivot_longer(-dose) %>%
      ggplot(aes(dose, value)) + geom_line(linewidth = 1, colour = "#01579b") +
      geom_point() + facet_wrap(~name, scales = "free_y") +
      labs(x = "Mitapivat (mg BID)", y = NULL) + thm
  })

  output$drugtab <- renderTable({
    base <- par_now()
    arms <- list(`Baseline` = c(),
                 `Mitapivat 100 mg` = c(dose_m = 100),
                 `Splenectomy` = c(spl_frac = 0),
                 `Combination` = c(dose_m = 100, spl_frac = 0))
    lapply(names(arms), function(a) {
      p <- base; if (length(arms[[a]])) p[names(arms[[a]])] <- arms[[a]]
      s <- tail(run_one(p, end = 1200, delta = 1200), 1)
      tibble(Regimen = a, Hb = s$Hb, `Reticulocytes %` = s$RETpct, `Lifespan` = s$LIFE,
             `ATP` = NA_real_, `Total bilirubin` = s$TBIL)
    }) %>% bind_rows()
  }, digits = 2)

  ## ---- 12 stones and iron ----------------------------------------------
  output$stone <- renderPlot({
    lapply(c("1", "0.287"), function(u) {
      p <- par_now(); p["ugt_f"] <- as.numeric(u)
      run_one(p, end = 40 * 365, delta = 90) %>%
        mutate(grp = ifelse(u == "1", "UGT1A1 normal", "Gilbert syndrome"))
    }) %>% bind_rows() %>%
      select(time, grp, STONEPCT, TBIL, FELIV = TBIL) %>%
      pivot_longer(c(-time, -grp)) %>%
      ggplot(aes(time / 365, value, colour = grp)) +
      geom_line(linewidth = 1) + facet_wrap(~name, scales = "free_y") +
      labs(x = "Age (years)", y = NULL, colour = NULL) + thm
  })

  output$stonetab <- renderTable({
    s <- ss()
    tibble(`Gallstone index (%)` = s$STONEPCT, `Total bilirubin` = s$TBIL,
           `Bilirubin production (mg/day)` = s$BRPROD)
  }, digits = 2)

  ## ---- 13 scenarios ----------------------------------------------------
  output$scen <- renderPlot({
    base <- par_now()
    arms <- list(
      `1 Natural history` = c(),
      `2 Total splenectomy` = c(spl_frac = 0),
      `3 Partial splenectomy 20%` = c(spl_frac = 0.2, k_regrow = 1 / 700),
      `4 Mitapivat 100 mg` = c(dose_m = 100),
      `5 Splenectomy + mitapivat` = c(spl_frac = 0, dose_m = 100),
      `6 Scheduled transfusion` = c(tx_start = 0),
      `7 Folate deficiency` = c(fol_ok = 0.30),
      `8 Comorbid Gilbert syndrome` = c(ugt_f = 0.287))
    lapply(names(arms), function(a) {
      p <- base; if (length(arms[[a]])) p[names(arms[[a]])] <- arms[[a]]
      run_one(p, end = input$tend, delta = 5) %>% mutate(arm = a)
    }) %>% bind_rows() %>%
      select(time, arm, Hb, RETpct, TBIL, LIFE, AREA, SPLEEN) %>%
      pivot_longer(c(-time, -arm)) %>%
      ggplot(aes(time, value, colour = arm)) + geom_line(linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "Time (days)", y = NULL, colour = NULL) + thm
  })

  output$scentab <- renderTable({
    base <- par_now()
    arms <- list(`Natural history` = c(), `Total splenectomy` = c(spl_frac = 0),
                 `Mitapivat` = c(dose_m = 100),
                 `Splenectomy + mitapivat` = c(spl_frac = 0, dose_m = 100))
    lapply(names(arms), function(a) {
      p <- base; if (length(arms[[a]])) p[names(arms[[a]])] <- arms[[a]]
      s <- tail(run_one(p, end = 1200, delta = 1200), 1)
      tibble(Scenario = a, Hb = s$Hb, `Reticulocytes %` = s$RETpct,
             `Lifespan` = s$LIFE, `Bilirubin` = s$TBIL, `Spleen mL` = s$SPLEEN)
    }) %>% bind_rows()
  }, digits = 2)

  ## ---- 14 virtual population -------------------------------------------
  pop <- eventReactive(input$go, {
    set.seed(20260806)
    n <- input$npop
    tibble(id = seq_len(n),
           fdef = pmin(pmax(rnorm(n, 0.28, 0.09), 0.05), 0.55),
           fb3  = ifelse(runif(n) < 0.22, 0.15, 1.00),
           ugt  = ifelse(runif(n) < 0.12, 0.287, 1.00)) %>%
      rowwise() %>%
      mutate(res = list(tail(run_one(c(par_now(),
                                       fdef = fdef, f_b3ves = fb3,
                                       ugt_f = ugt),
                                     end = 1100, delta = 1100), 1))) %>%
      tidyr::unnest(res) %>% ungroup()
  })

  output$pop <- renderPlot({
    pop() %>% select(Hb, RETpct, TBIL, LIFE, MCHC, EMA) %>%
      mutate(`EMA reduction %` = 100 * (1 - EMA)) %>% select(-EMA) %>%
      pivot_longer(everything()) %>%
      ggplot(aes(value)) + geom_histogram(bins = 25, fill = "#1565c0") +
      facet_wrap(~name, scales = "free") + labs(x = NULL, y = "Number of patients") + thm
  })

  output$poptab <- renderTable({
    d <- pop()
    tibble(Metric = c("Hb", "Reticulocytes %", "Bilirubin", "Lifespan", "EMA reduction %"),
           p10 = c(quantile(d$Hb, .1), quantile(d$RETpct, .1),
                   quantile(d$TBIL, .1), quantile(d$LIFE, .1),
                   quantile(100 * (1 - d$EMA), .1)),
           Median = c(median(d$Hb), median(d$RETpct), median(d$TBIL),
                      median(d$LIFE), median(100 * (1 - d$EMA))),
           p90 = c(quantile(d$Hb, .9), quantile(d$RETpct, .9),
                   quantile(d$TBIL, .9), quantile(d$LIFE, .9),
                   quantile(100 * (1 - d$EMA), .9)))
  }, digits = 2)
}

shinyApp(ui, server)
