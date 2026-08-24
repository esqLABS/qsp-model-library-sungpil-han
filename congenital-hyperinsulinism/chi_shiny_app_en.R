## =====================================================================
##  Congenital Hyperinsulinism (CHI) — QSP explorer (Shiny)
##  ---------------------------------------------------------------------
##  The app is built around ONE control: the slider for g, the fraction of
##  beta-cell K_ATP conductance that survives.  Every tab is a different
##  way of reading the consequences of that one number.
##
##  Run:  shiny::runApp("chi_shiny_app.R")
##  Requires: shiny, mrgsolve, dplyr, tidyr, ggplot2, DT
##  The model code is sourced from chi_mrgsolve_model.R (same directory).
## =====================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

## ---------------------------------------------------------------------
## Load the compiled model + helpers from the model file.  We source it in
## a throw-away environment and lift out what we need, so that running the
## model file's own scenarios does not happen here.
## ---------------------------------------------------------------------
.src <- new.env()
local({
  code <- readLines("chi_mrgsolve_model.R")
  ## keep everything up to (but not including) the scenario section
  cut <- grep("^##  SCENARIO 1", code)
  if (length(cut)) code <- code[seq_len(cut[1] - 3)]
  eval(parse(text = paste(code, collapse = "\n")), envir = .src)
})
mod  <- .src$mod
geno <- .src$geno

FEED_MGKGMIN <- 845 * 8 / 1440    # 4.69 mg/kg/min enteral baseline

GENO_LABEL <- c(
  "Normal newborn (normal)"                              = "normal",
  "K_ATP dominant (ABCC8 dominant, g≈0.60)"               = "katp_dominant",
  "K_ATP recessive diffuse (biallelic diffuse, g≈0.02)"      = "katp_recess_diff",
  "K_ATP focal type (focal, paternal UPD 11p15)"           = "katp_focal",
  "GDH-HI / HI-HA (GLUD1 activation)"                     = "gdh_hi",
  "GCK activation (set-point shift)"                        = "gck_activating",
  "HNF4A (transient, late MODY)"                          = "hnf4a",
  "SCHAD/HADH deficiency"                                   = "schad",
  "SLC16A1 (exercise-induced)"                               = "slc16a1",
  "After subtotal pancreatectomy (remnant 5%)"                      = "post_pancreatect"
)

feeds     <- function(days = 4)  ev(amt = 845, cmt = "GGUT", time = 0,
                                   ii = 3, addl = days * 8 - 1)
diazoxide <- function(d, days = 6) ev(amt = d/3, cmt = "DZXg", time = 0,
                                   ii = 8, addl = days * 3 - 1)
octreotide<- function(d, days = 6) ev(amt = d/4, cmt = "OCTs", time = 0,
                                   ii = 6, addl = days * 4 - 1)
ersodetug <- function(d)          ev(amt = d/0.055, cmt = "ERS", time = 0)
nifedipine<- function(d, days=6)  ev(amt = d/1.2, cmt = "NIF", time = 0,
                                   ii = 8, addl = days * 3 - 1)
sirolimus <- function(d, days=10) ev(amt = d, cmt = "SIRg", time = 0,
                                   ii = 24, addl = days - 1)

## ---------------------------------------------------------------------
## one simulation, assembled from the UI state
## ---------------------------------------------------------------------
simulate <- function(gt, g_override = NA, dzx = 0, oct = 0, ers = 0, nif = 0,
                     sir = 0, gcg = 0, remnant = 1, fed = TRUE,
                     closed_loop = TRUE, target = 70, end = 72,
                     fast_from = NA) {
  p <- geno[[gt]]
  if (!is.na(g_override) && gt != "normal") p$g_ab <- g_override
  p$BMASS0 <- remnant
  p$GCG_inf <- gcg
  if (closed_loop) { p$loop <- 1; p$Gtarget <- target }

  e <- NULL
  nfeed <- if (is.na(fast_from)) ceiling(end/24) + 1 else max(1, floor(fast_from/24) + 1)
  if (fed) {
    nf <- if (is.na(fast_from)) nfeed * 8 else max(1, floor(fast_from/3))
    e <- ev(amt = 845, cmt = "GGUT", time = 0, ii = 3, addl = nf - 1)
  }
  add <- function(x) if (is.null(e)) x else c(e, x)
  if (dzx > 0) e <- add(diazoxide(dzx, ceiling(end/24) + 1))
  if (oct > 0) e <- add(octreotide(oct, ceiling(end/24) + 1))
  if (ers > 0) e <- add(ersodetug(ers))
  if (nif > 0) e <- add(nifedipine(nif, ceiling(end/24) + 1))
  if (sir > 0) e <- add(sirolimus(sir, ceiling(end/24) + 1))

  out <- mod %>% param(p)
  if (!is.null(e)) out <- out %>% ev(e)
  out %>% mrgsim(end = end, delta = 0.05) %>% as_tibble()
}

tail_stats <- function(d, frac = 0.25) {
  tmax <- max(d$time); d <- filter(d, time > tmax * (1 - frac))
  tibble(
    glucose_mean = mean(d$GLU), glucose_min = min(d$GLU),
    insulin      = mean(d$INS), cpeptide = mean(d$CPEP),
    GIR          = mean(d$GIRout), total_glc = mean(d$TOTGLC),
    BOHB         = mean(d$BOHB), FFA = mean(d$FFA),
    ammonia      = mean(d$NH3), glycogen = mean(d$GLY),
    fuel_ratio   = mean(d$fuelrat), Vm = mean(d$Vm_ab),
    pct_below70  = 100 * mean(d$GLU < 70), pct_below54 = 100 * mean(d$GLU < 54)
  )
}

THEME <- theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom", plot.title = element_text(face = "bold"))

## =====================================================================
## UI
## =====================================================================
ui <- fluidPage(
  titlePanel("Congenital Hyperinsulinism — QSP Explorer"),
  tags$p(style = "color:#555;margin-top:-8px",
    HTML("Every result in this app is the arithmetic consequence of <b>a single number</b>: the residual beta-cell K<sub>ATP</sub> conductance <b>g</b>. ",
         "Diazoxide <b>multiplies</b> g, octreotide <b>adds to</b> the voltage divider, ",
         "and glucagon · glucose · ersodetug · surgery are <b>independent of</b> g.")),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectInput("gt", "Genotype", choices = GENO_LABEL,
                  selected = "katp_recess_diff"),
      sliderInput("g", HTML("Residual K<sub>ATP</sub> conductance <b>g</b>"),
                  min = 0, max = 1, value = 0.02, step = 0.01),
      helpText(HTML("<small>g≈0–0.05 recessive diffuse · 0.5–0.8 dominant · 1.0 normal<br>",
                    "at g&lt;0.1, <b>severity saturates</b> (the membrane potential is already pinned at the leak potential).</small>")),
      hr(),
      h5("Therapy"),
      sliderInput("dzx", "Diazoxide (mg/kg/day)", 0, 20, 0, step = 2.5),
      sliderInput("oct", "Octreotide (µg/kg/day)", 0, 40, 0, step = 5),
      sliderInput("gcg", "Glucagon continuous infusion (µg/kg/h)", 0, 20, 0, step = 2.5),
      sliderInput("ers", "Ersodetug (mg/kg IV)", 0, 12, 0, step = 3),
      sliderInput("nif", "Nifedipine (mg/kg q8h)", 0, 1, 0, step = 0.25),
      sliderInput("sir", "Sirolimus (mg/kg/day)", 0, 0.12, 0, step = 0.02),
      hr(),
      h5("Surgery / Management"),
      sliderInput("remnant", "Residual beta-cell mass (post-surgery)", 0.02, 1.0, 1.0, step = 0.02),
      checkboxInput("loop", "IV glucose autotitration (closed loop)", TRUE),
      sliderInput("target", "Target glucose (mg/dL)", 50, 100, 70, step = 5),
      checkboxInput("fed", "Enteral feeding q3h", TRUE),
      sliderInput("end", "Simulation duration (h)", 24, 168, 72, step = 12)
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        type = "tabs",
        tabPanel("① Patient profile",
          br(), fluidRow(column(12, DTOutput("profile"))),
          br(), plotOutput("betacell", height = "330px"),
          helpText(HTML("<small>Left: the voltage divider. When K<sub>ATP</sub> conductance dominates the leak, the cell is silent; ",
                        "when it disappears, the cell is pinned at −20 mV and secretes insulin <b>at every glucose level</b>.</small>"))),
        tabPanel("② PK (drug concentration)",
          br(), plotOutput("pk", height = "560px"),
          helpText(HTML("<small>Diazoxide t½ ≈ 20 h (up to 30 h in newborns) — steady state reached in 3–5 days. ",
                        "Octreotide t½ ≈ 1.7 h, so q6h dosing shows a large amplitude.</small>"))),
        tabPanel("③ PD Key Metrics",
          br(), plotOutput("pd", height = "620px")),
        tabPanel("④ Clinical Endpoints",
          br(), fluidRow(column(12, DTOutput("endpoints"))),
          br(), plotOutput("gir", height = "330px")),
        tabPanel("⑤ Scenario Comparison",
          br(),
          helpText("Compares IV glucose requirement when each agent is applied alone in the same genotype."),
          plotOutput("compare", height = "380px"),
          br(), DTOutput("compare_tbl")),
        tabPanel("⑥ Biomarkers",
          br(), plotOutput("biomarker", height = "620px"),
          helpText(HTML("<small>The absence of ketones is CHI's <b>diagnostic core</b>. Because lipolysis IC₅₀(12) &lt; ketogenesis(15) ",
                        "&lt; glycolysis(30) &lt; gluconeogenesis(45) µU/mL, an insulin level capable of causing hypoglycaemia has ",
                        "already switched off ketones. Ammonia rises only in GDH-HI and does not respond to diazoxide.</small>"))),
        tabPanel("⑦ Cerebral Fuel (key)",
          br(), plotOutput("fuel", height = "400px"),
          br(), plotOutput("fuelthresh", height = "330px"),
          helpText(HTML("<small>What the model computes: with zero ketones, the brain runs out of fuel at <b>43.7 mg/dL</b>, ",
                        "but at BOHB 2 mM it can tolerate down to <b>25.3 mg/dL</b>. Hyperinsulinism erases ",
                        "an <b>18.4 mg/dL safety margin</b> before glucose is even taken away — this is why ketotic hypoglycaemia ",
                        "is tolerated into the 40s, yet the CHI target is 70 mg/dL.</small>"))),
        tabPanel("⑧ g-Axis Exploration",
          br(),
          helpText(HTML("The model's central claim: <b>same endpoint, same dose, opposite conclusions.</b>")),
          plotOutput("gsweep", height = "420px"),
          br(), DTOutput("gsweep_tbl")),
        tabPanel("⑨ Surgical Planning",
          br(), plotOutput("surgery", height = "380px"),
          br(), DTOutput("surgery_tbl"),
          helpText(HTML("<small>Surgery reduces <b>B (mass)</b> and leaves <b>g</b> unchanged. So a single equation ",
                        "produces both persistent hypoglycaemia and post-surgical diabetes. A remnant of 0.20–0.30 is the only window, ",
                        "and ≤0.10 is immediate diabetes. Moreover, growth re-reads the same remnant: a ",
                        "2% remnant that is sufficient at 3.5 kg is only 0.35% of normal at 20 kg.</small>"))),
        tabPanel("⑩ Documentation",
          br(),
          h4("What the model calibrated and what it predicted"),
          HTML("<p>The numbers used for calibration are <b>9</b>, and 8 of them are <b>normal newborn physiology</b>:
                insulin 5 µU/mL at a glucose of 75 mg/dL, the <i>shape</i> of the normal glucose→insulin dose-response
                (relative to 75 mg/dL: 0.25-fold at 45, 1.8-fold at 90, 12-fold at 150, 24-fold at 250),
                newborn cerebral glucose consumption of 4.2 mg/kg/min, term-infant hepatic glycogen of about 2000 mg/kg,
                and the ranking of insulin IC₅₀ values. Only one number was taken from
                CHI itself (<b>gGIRK=1.2</b>, the observation that octreotide reduces the glucose requirement of
                severe diffuse CHI by about half).</p>
                <p>So everything below is a <b>prediction</b>: the total glucose requirement of severe diffuse CHI at
                15.9 mg/kg/min · the genotype separation of diazoxide response rate (0.1% at g=0.02, 100% at g≥0.3) ·
                severity saturation below g&lt;0.1 · the absence of ketones · a positive glucagon stimulation test (+62 mg/dL) ·
                the narrow window for resection extent · GDH-HI's leucine sensitivity and ammonia that does not respond to diazoxide ·
                nifedipine's failure (the mechanism is right, but EC₅₀ cannot be reached).</p>
                <h4>A tension honestly disclosed</h4>
                <p>The predicted plasma insulin (≈83 µU/mL) is
                <b>higher</b> than reported critical-sample values (10–50 µU/mL). This is not a choice but forced arithmetic — given
                newborn insulin sensitivity, a glucose requirement of 16 mg/kg/min cannot be explained by 20 µU/mL. This is the most
                exposed parameter in this model, and it is documented in full as T1–T4 in the README.</p>
                <p style='color:#a33'><b>Disclaimer:</b> this is a QSP model for educational/research purposes.
                Do not use it for clinical decision-making or dosing.</p>")
        )
      )
    )
  )
)

## =====================================================================
## SERVER
## =====================================================================
server <- function(input, output, session) {

  observeEvent(input$gt, {
    gdef <- geno[[input$gt]]$g_ab
    if (!is.null(gdef)) updateSliderInput(session, "g", value = gdef)
    rem <- geno[[input$gt]]$BMASS0
    updateSliderInput(session, "remnant", value = if (is.null(rem)) 1 else rem)
  })

  sim <- reactive({
    simulate(input$gt, g_override = input$g, dzx = input$dzx, oct = input$oct,
             ers = input$ers, nif = input$nif, sir = input$sir, gcg = input$gcg,
             remnant = input$remnant, fed = input$fed,
             closed_loop = input$loop, target = input$target, end = input$end)
  })

  ## ---------------- ① profile ----------------
  output$profile <- renderDT({
    s <- tail_stats(sim())
    datatable(tibble(
      Item = c("Residual K_ATP conductance g", "Beta-cell membrane potential (mV)",
               "Mean glucose (mg/dL)", "Lowest glucose (mg/dL)",
               "Plasma insulin (µU/mL)", "C-peptide (ng/mL)",
               "IV glucose (mg/kg/min)", "Total glucose supply (mg/kg/min)",
               "Hepatic glycogen (mg/kg)", "Cerebral fuel supply/demand ratio"),
      Value = c(sprintf("%.2f", input$g), sprintf("%.1f", s$Vm),
             sprintf("%.1f", s$glucose_mean), sprintf("%.1f", s$glucose_min),
             sprintf("%.1f", s$insulin), sprintf("%.2f", s$cpeptide),
             sprintf("%.2f", s$GIR), sprintf("%.2f", s$total_glc),
             sprintf("%.0f", s$glycogen), sprintf("%.3f", s$fuel_ratio))
    ), options = list(dom = "t", pageLength = 12), rownames = FALSE)
  })

  output$betacell <- renderPlot({
    gg <- seq(0, 1, 0.01)
    p  <- as.list(param(mod))
    Gs <- c(45, 60, 90, 150)
    d <- expand.grid(g = gg, G = Gs) %>% as_tibble() %>% rowwise() %>%
      mutate({
        Gb <- G/18.016
        Rt <- p$Rbas + (p$Rmax - p$Rbas)*Gb^p$hG/(p$KG^p$hG + Gb^p$hG) +
              p$kAA*(p$AA0/(p$KAA + p$AA0))
        Po <- 1/(1 + (Rt/p$R50)^p$nR)
        GK <- p$gKmax*g*Po
        tibble(Vm = ((GK)*p$EK + p$gleak*p$Eleak)/(GK + p$gleak))
      }) %>% ungroup()
    ggplot(d, aes(g, Vm, colour = factor(G))) +
      geom_hline(yintercept = -50, linetype = "dashed", colour = "grey40") +
      annotate("text", x = 0.8, y = -47, label = "Ca_v activation threshold (−50 mV)",
               size = 3.4, colour = "grey30") +
      geom_vline(xintercept = input$g, colour = "#d1495b", linewidth = 1) +
      geom_line(linewidth = 1.1) + scale_x_reverse() +
      labs(title = "g alone sets the membrane potential",
           x = "Residual K_ATP conductance g (right = normal)", y = "Resting membrane potential (mV)",
           colour = "Glucose (mg/dL)") + THEME
  })

  ## ---------------- ② PK ----------------
  output$pk <- renderPlot({
    d <- sim() %>%
      select(time, `Diazoxide (µg/mL)` = DZX, `Octreotide (ng/mL)` = OCT,
             `Ersodetug (µg/mL)` = ERS, `Sirolimus (mg/L)` = SIR,
             `Nifedipine (mg/L)` = NIF, `Glucagon (pg/mL)` = GCG) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value)) + geom_line(linewidth = 0.9, colour = "#2f6f9f") +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "Time (h)", y = NULL, title = "Pharmacokinetics") + THEME
  })

  ## ---------------- ③ PD ----------------
  output$pd <- renderPlot({
    d <- sim() %>%
      select(time, `Glucose (mg/dL)` = GLU, `CGM interstitial fluid (mg/dL)` = GLUi,
             `Insulin (µU/mL)` = INS, `Receptor-bound insulin X (µU/mL)` = insact,
             `C-peptide (ng/mL)` = CPEP,
             `Beta-cell membrane potential (mV)` = Vm_ab) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value)) + geom_line(linewidth = 0.9, colour = "#c0504d") +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "Time (h)", y = NULL,
           title = "Pharmacodynamics — glucose · insulin · membrane potential") + THEME
  })

  ## ---------------- ④ endpoints ----------------
  output$endpoints <- renderDT({
    s <- tail_stats(sim()); d <- sim()
    datatable(tibble(
      Endpoint = c("Mean glucose", "Lowest glucose", "% time below 70 mg/dL",
                     "% time below 54 mg/dL", "Hypoglycaemia area (mg/dL·h)",
                     "Cumulative cerebral fuel deficit (h)", "Neurodevelopmental deficit index",
                     "IV glucose requirement", "Total glucose supply"),
      Value = c(sprintf("%.1f mg/dL", s$glucose_mean),
             sprintf("%.1f mg/dL", s$glucose_min),
             sprintf("%.1f %%", s$pct_below70), sprintf("%.1f %%", s$pct_below54),
             sprintf("%.0f", max(d$AUCHYPO)), sprintf("%.2f", max(d$TFUEL)),
             sprintf("%.2f", max(d$DEV)),
             sprintf("%.2f mg/kg/min", s$GIR),
             sprintf("%.2f mg/kg/min", s$total_glc))
    ), options = list(dom = "t", pageLength = 10), rownames = FALSE)
  })

  output$gir <- renderPlot({
    d <- sim()
    ggplot(d, aes(time)) +
      geom_ribbon(aes(ymin = 0, ymax = GIRout), fill = "#93cdf7", alpha = 0.6) +
      geom_line(aes(y = GLU/10), colour = "#c0504d", linewidth = 1) +
      scale_y_continuous(name = "IV glucose (mg/kg/min)",
        sec.axis = sec_axis(~.*10, name = "Glucose (mg/dL)")) +
      labs(x = "Time (h)", title = "Glucose requirement and glucose") + THEME
  })

  ## ---------------- ⑤ compare ----------------
  cmp <- reactive({
    arms <- list("No treatment" = list(), "Diazoxide 15" = list(dzx = 15),
                 "Octreotide 30" = list(oct = 30),
                 "Glucagon 15" = list(gcg = 15),
                 "Ersodetug 9" = list(ers = 9),
                 "Nifedipine 0.5" = list(nif = 0.5),
                 "Octreotide+glucagon" = list(oct = 30, gcg = 15))
    lapply(names(arms), function(n) {
      a <- arms[[n]]
      d <- do.call(simulate, c(list(gt = input$gt, g_override = input$g,
                                    remnant = input$remnant, fed = input$fed,
                                    closed_loop = TRUE, target = input$target,
                                    end = 72), a))
      tail_stats(d) %>% mutate(arm = n)
    }) %>% bind_rows()
  })

  output$compare <- renderPlot({
    d <- cmp() %>% mutate(arm = factor(arm, levels = arm))
    ggplot(d, aes(arm, GIR, fill = GIR)) +
      geom_col(width = 0.65) + coord_flip() +
      scale_fill_gradient(low = "#4c9f70", high = "#d1495b", guide = "none") +
      labs(x = NULL, y = "IV glucose requirement (mg/kg/min)",
           title = sprintf("Treatment response at g = %.2f", input$g)) + THEME
  })

  output$compare_tbl <- renderDT({
    d <- cmp()
    base <- d$GIR[1]
    datatable(d %>% transmute(Treatment = arm,
        `IV glucose` = round(GIR, 2), `% reduction` = round(100*(base - GIR)/max(base, 1e-9), 1),
        `Mean glucose` = round(glucose_mean, 1), `Insulin` = round(insulin, 1),
        `BOHB` = round(BOHB, 3)),
      options = list(dom = "t", pageLength = 10), rownames = FALSE)
  })

  ## ---------------- ⑥ biomarkers ----------------
  output$biomarker <- renderPlot({
    d <- sim() %>%
      select(time, `BOHB (mmol/L)` = BOHB, `Free fatty acids (mmol/L)` = FFA,
             `Hepatic glycogen (mg/kg)` = GLY, `Ammonia (µmol/L)` = NH3,
             `Lactate (mmol/L)` = LAC, `Glucagon (pg/mL)` = GCG) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value)) + geom_line(linewidth = 0.9, colour = "#8a7355") +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "Time (h)", y = NULL, title = "Biomarkers") + THEME
  })

  ## ---------------- ⑦ cerebral fuel ----------------
  output$fuel <- renderPlot({
    p <- as.list(param(mod)); d <- sim()
    dd <- d %>% mutate(
      glc  = p$UbrMax*GLU/(p$Kbr + GLU)*(1 + p$acbf*pmax(0, p$Gcbf - GLU)/p$Gcbf),
      ket  = p$VketMax*BOHB/(p$Kket + BOHB),
      lact = p$VlacMax*LAC/(p$Klac + LAC)) %>%
      select(time, `Glucose` = glc, `Ketone` = ket, `Lactate` = lact) %>%
      pivot_longer(-time)
    ggplot(dd, aes(time, value, fill = name)) +
      geom_area() +
      geom_hline(yintercept = p$CMRreq, linetype = "dashed", linewidth = 1) +
      annotate("text", x = max(dd$time)*0.1, y = p$CMRreq*1.06,
               label = "Cerebral energy requirement", size = 3.6) +
      scale_fill_manual(values = c(`Glucose` = "#f39c9c", `Ketone` = "#cbb99a",
                                   `Lactate` = "#d9ecf7")) +
      labs(x = "Time (h)", y = "Cerebral fuel supply (mg/kg/min, glucose-equivalent)",
           fill = NULL, title = "The brain feeds on the sum") + THEME
  })

  output$fuelthresh <- renderPlot({
    p <- as.list(param(mod))
    f <- function(G, B) {
      p$UbrMax*G/(p$Kbr + G)*(1 + p$acbf*max(0, p$Gcbf - G)/p$Gcbf) +
        p$VketMax*B/(p$Kket + B) + p$VlacMax*1.3/(p$Klac + 1.3)
    }
    d <- tibble(BOHB = seq(0, 5, 0.05)) %>% rowwise() %>%
      mutate(thr = {
        lo <- 1; hi <- 200
        for (i in 1:60) { m <- (lo + hi)/2; if (f(m, BOHB) < p$CMRreq) lo <- m else hi <- m }
        (lo + hi)/2 }) %>% ungroup()
    ggplot(d, aes(BOHB, thr)) +
      geom_line(linewidth = 1.2, colour = "#b5423f") +
      geom_point(data = tibble(BOHB = 0, thr = 43.7), size = 3, colour = "#d1495b") +
      annotate("text", x = 0.9, y = 46,
               label = "Zero ketones (CHI): 43.7 mg/dL", size = 3.8, colour = "#d1495b") +
      annotate("text", x = 2.6, y = 27,
               label = "BOHB 2 mM (ketotic): 25.3 mg/dL", size = 3.8, colour = "#8a7355") +
      labs(x = "Plasma BOHB (mmol/L)",
           y = "Glucose at which cerebral fuel is exhausted (mg/dL)",
           title = "Hyperinsulinism erases the safety margin before glucose is even taken away") + THEME
  })

  ## ---------------- ⑧ g sweep ----------------
  gsw <- reactive({
    gs <- c(0.6, 0.4, 0.3, 0.2, 0.1, 0.05, 0.02)
    lapply(gs, function(g) {
      b <- simulate(input$gt, g_override = g, closed_loop = TRUE,
                    target = input$target, end = 72)
      dz <- simulate(input$gt, g_override = g, dzx = 15, closed_loop = TRUE,
                     target = input$target, end = 72)
      oc <- simulate(input$gt, g_override = g, oct = 30, closed_loop = TRUE,
                     target = input$target, end = 72)
      b0 <- tail_stats(b)$GIR
      tibble(g = g, GIR0 = b0,
             `Diazoxide` = if (b0 > 0.05) 100*(b0 - tail_stats(dz)$GIR)/b0 else NA,
             `Octreotide` = if (b0 > 0.05) 100*(b0 - tail_stats(oc)$GIR)/b0 else NA)
    }) %>% bind_rows()
  })

  output$gsweep <- renderPlot({
    d <- gsw() %>% pivot_longer(c(`Diazoxide`, `Octreotide`))
    ggplot(d, aes(g, value, colour = name)) +
      geom_line(linewidth = 1.3) + geom_point(size = 2.6) + scale_x_reverse() +
      scale_colour_manual(values = c(`Diazoxide` = "#2f6f9f",
                                     `Octreotide` = "#c0504d")) +
      labs(x = "Residual K_ATP conductance g (right = closer to normal)",
           y = "Reduction in IV glucose requirement (%)", colour = NULL,
           title = "One equation, two conclusions",
           subtitle = "Diazoxide multiplies g (vanishes as g→0); octreotide adds to the divider") +
      THEME
  })

  output$gsweep_tbl <- renderDT({
    datatable(gsw() %>% transmute(g, `Untreated glucose` = round(GIR0, 2),
        `Diazoxide % reduction` = round(`Diazoxide`, 1),
        `Octreotide % reduction` = round(`Octreotide`, 1)),
      options = list(dom = "t"), rownames = FALSE)
  })

  ## ---------------- ⑨ surgery ----------------
  surg <- reactive({
    lapply(c(1.0, 0.5, 0.3, 0.2, 0.1, 0.05, 0.02), function(b) {
      d <- simulate(input$gt, g_override = input$g, remnant = b,
                    closed_loop = TRUE, target = input$target, end = 48)
      tail_stats(d) %>% mutate(remnant = b)
    }) %>% bind_rows()
  })

  output$surgery <- renderPlot({
    d <- surg()
    ggplot(d, aes(remnant)) +
      geom_hline(yintercept = c(70, 126), linetype = "dashed", colour = "grey45") +
      geom_line(aes(y = glucose_mean), linewidth = 1.3, colour = "#c0504d") +
      geom_point(aes(y = glucose_mean), size = 2.6, colour = "#c0504d") +
      geom_col(aes(y = GIR*10), fill = "#93cdf7", alpha = 0.55, width = 0.02) +
      scale_x_log10() +
      scale_y_continuous(name = "Mean glucose (mg/dL)",
        sec.axis = sec_axis(~./10, name = "IV glucose (mg/kg/min)")) +
      labs(x = "Residual beta-cell mass (log)",
           title = "Resection reduces B and leaves g behind — so the window is narrow") + THEME
  })

  output$surgery_tbl <- renderDT({
    datatable(surg() %>% transmute(`Remnant` = remnant,
        `IV glucose` = round(GIR, 2), `Mean glucose` = round(glucose_mean, 1),
        `Insulin` = round(insulin, 1),
        Interpretation = ifelse(GIR > 0.3, "Still glucose-dependent",
                 ifelse(glucose_mean < 126, "Normal glucose", "Diabetes"))),
      options = list(dom = "t"), rownames = FALSE)
  })
}

shinyApp(ui, server)
