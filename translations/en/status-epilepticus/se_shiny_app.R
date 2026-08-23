## =====================================================================
##  se_shiny_app.R
##  Status Epilepticus QSP model — interactive dashboard
##  Status Epilepticus QSP Model — Shiny Dashboard
##
##  10 tabs, each answering one question the model exists to answer:
##
##   1  Patient & Aetiology   — set the DRIVE, see where the
##                                          thresholds land for this patient
##   2  Treatment Timeline — build a drug schedule minute by
##                                          minute and watch the clock run
##   3  The Receptor Clock   — R_SYN out, NR_SYN in, R_EXTRA flat
##   4  Drug PK (Pharmacokinetics)         — 11 drugs, plasma and effect site
##   5  The Price of a Minute — CREQ_BZD(t) and its asymptote
##   6  The Crossover             — benzodiazepine vs ketamine target
##   7  Motor-EEG Dissociation (Motor vs EEG)      — what the bedside sees vs the brain
##   8  Systemic Physiology    — the phase I / phase II switch
##   9  Brain Injury (Injury & outcome)         — the second integrator
##  10  Virtual Population     — response RATES, ESETT / RAMPART
##
##  run:  shiny::runApp("se_shiny_app.R")
##  (loads se_mrgsolve_model.R from the same directory)
## =====================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

# ---------------------------------------------------------------------
#  Load the model without executing the script's own reporting section
# ---------------------------------------------------------------------
.src <- readLines("se_mrgsolve_model.R")
.cut <- grep("^#  MAIN$", .src)
if (length(.cut)) .src <- .src[1:(.cut[1] - 2)]
eval(parse(text = paste(.src, collapse = "\n")))

theme_se <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold"),
        legend.position = "bottom")

PALETTE <- c("#1E88E5", "#E65100", "#1B5E20", "#880E4F", "#4527A0",
             "#00838F", "#B71C1C", "#546E7A", "#F9A825", "#6A1B9A")

# =====================================================================
#  UI
# =====================================================================
ui <- fluidPage(
  titlePanel("Status Epilepticus QSP Model — a receptor-trafficking clock"),
  tags$p(style = "color:#555;margin-top:-8px",
         HTML("<b>59 ODEs · 11 drugs · 22 scenarios.</b> One clock drives two receptor pools in opposite directions and leaves a third alone. Every curve below is that arithmetic.")),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      h4("① Patient & Aetiology"),
      selectInput("aetiology", "Aetiology",
                  choices = c("AED withdrawal"          = 1.15,
                              "Metabolic"             = 1.35,
                              "Typical adult SE"      = 1.55,
                              "Acute stroke"     = 1.80,
                              "CNS infection"       = 2.10,
                              "Autoimmune encephalitis (anti-NMDAR)"     = 2.40,
                              "FIRES / NORSE"                  = 2.60,
                              "Hypoxic (post-anoxic)"         = 3.20),
                  selected = 1.55),
      checkboxInput("reversible", "Aetiology is treatable (drive removed at 30 min)", FALSE),
      sliderInput("alb", "Serum albumin (g/dL)", 1.8, 5.0, 4.0, 0.1),
      sliderInput("crcl", "Creatinine clearance (mL/min)", 10, 140, 90, 5),
      checkboxInput("anak", "IL-1 blockade (anakinra)", FALSE),
      checkboxInput("gluciv", "IV glucose + thiamine", FALSE),
      checkboxInput("cool", "Targeted temperature management (TTM)", FALSE),

      hr(),
      h4("② First-line agent"),
      selectInput("bzd", "benzodiazepine",
                  c("None" = "none",
                    "lorazepam 4 mg IV"        = "lzp4",
                    "lorazepam 2 mg IV (low dose)" = "lzp2",
                    "midazolam 10 mg IM"       = "mdzim",
                    "midazolam 0.2 mg/kg IV"   = "mdziv",
                    "diazepam 10 mg IV"        = "dzp"),
                  selected = "lzp4"),
      sliderInput("t_bzd", "Administration time (min)", 0, 120, 8, 1),
      checkboxInput("bzd_rep", "Second dose 4 minutes later", TRUE),

      hr(),
      h4("③ Second-line agent"),
      selectInput("second", "agent",
                  c("None" = "none",
                    "levetiracetam 60 mg/kg" = "lev",
                    "fosphenytoin 20 mg PE/kg" = "fos",
                    "valproate 40 mg/kg" = "vpa",
                    "phenobarbital 20 mg/kg" = "pb"),
                  selected = "lev"),
      sliderInput("t_second", "Administration time (min)", 5, 240, 30, 5),

      hr(),
      h4("④ Third-line agent (anaesthesia)"),
      checkboxGroupInput("third", NULL,
                         c("midazolam infusion 0.4 mg/kg/h" = "mdzinf",
                           "propofol 2 mg/kg + 6 mg/kg/h"   = "prop",
                           "ketamine 2 mg/kg + 3 mg/kg/h"   = "ket",
                           "allopregnanolone 86 ug/kg/h"    = "allo")),
      sliderInput("t_third", "Start time (min)", 30, 480, 120, 10),

      hr(),
      sliderInput("tend", "Simulation length (min)", 120, 1440, 480, 60),
      helpText(HTML("<small>All tabs read the same single simulation. Changing a parameter updates all 10 tabs at once.</small>"))
    ),

    mainPanel(
      width = 9,
      fluidRow(
        column(12,
          div(style = "background:#F1F8E9;border-left:5px solid #33691E;padding:8px 12px;margin-bottom:10px",
              htmlOutput("verdict")))
      ),
      tabsetPanel(
        id = "tabs",
        tabPanel("① Patient & Aetiology",
                 br(), plotOutput("p_gain", height = 320),
                 br(), tableOutput("t_thresholds"),
                 helpText("G = excitatory gain / inhibitory tone. If G > 1 the seizure grows, if G < 1 it stops. Aetiology moves only the numerator; drugs mostly move only the denominator.")),
        tabPanel("② Treatment Timeline",
                 br(), plotOutput("p_timeline", height = 430),
                 br(), tableOutput("t_summary")),
        tabPanel("③ The Receptor Clock",
                 br(), plotOutput("p_clock", height = 400),
                 br(), plotOutput("p_pools", height = 250),
                 helpText("While the red curve (NR_SYN) rises, the blue curve (R_SYN × F_BZS) falls. Green (R_EXTRA) barely moves — the drugs that still work late are exactly the ones that use this pool.")),
        tabPanel("④ Drug PK",
                 br(), plotOutput("p_pk", height = 430),
                 br(), plotOutput("p_occ", height = 260)),
        tabPanel("⑤ The Price of a Minute",
                 br(), plotOutput("p_creq", height = 380),
                 br(), tableOutput("t_creq"),
                 helpText("CREQ_BZD(t) = the benzodiazepine effect-site concentration (as a multiple of EC50) needed right now to bring G back to 1. Vertical asymptote = the time point beyond which no dose can do it.")),
        tabPanel("⑥ The Crossover",
                 br(), plotOutput("p_cross", height = 380),
                 br(), tableOutput("t_cross"),
                 helpText("How many % each drug cuts G by, at the same occupancy (0.75). Dose and potency were held equal so only the target is compared.")),
        tabPanel("⑦ Motor-EEG Dissociation",
                 br(), plotOutput("p_diss", height = 400),
                 br(), tableOutput("t_diss"),
                 helpText("VA Cooperative: the same lorazepam is 64.9% in overt SE, 7.7-24.2% in subtle SE. What changed is not the drug but the window being read.")),
        tabPanel("⑧ Systemic Physiology",
                 br(), plotOutput("p_sys", height = 470),
                 helpText("Around 30 minutes, the compensated phase (rising catecholamines · increased CBF) crosses over into the decompensated phase (autoregulation failure · hypotension · hypoglycaemia · hyperthermia).")),
        tabPanel("⑨ Brain Injury",
                 br(), plotOutput("p_injury", height = 400),
                 br(), tableOutput("t_injury"),
                 helpText("SEIZ and INJURY are different state variables. 'The seizure stopped' and 'the brain was protected' are not the same statement.")),
        tabPanel("⑩ Virtual Population",
                 br(),
                 fluidRow(
                   column(4, numericInput("npop", "Population size", 200, 50, 600, 50)),
                   column(4, sliderInput("meddrive", "EDRIVE median", 1.2, 2.2, 1.55, 0.05)),
                   column(4, sliderInput("sddrive", "EDRIVE log SD", 0.10, 0.60, 0.30, 0.02))),
                 actionButton("runpop", "Run population simulation", class = "btn-primary"),
                 br(), br(),
                 tableOutput("t_pop"),
                 plotOutput("p_pop", height = 320),
                 helpText("Calibration target: first-line benzodiazepine 59-73% (PHTSE / RAMPART), second-line agent in benzodiazepine failures 45-47% (ESETT)."))
      )
    )
  )
)

# =====================================================================
#  SERVER
# =====================================================================
server <- function(input, output, session) {

  build_events <- reactive({
    e <- NULL
    add <- function(x) if (is.null(e)) x else e + x

    tb <- input$t_bzd
    if (input$bzd != "none") {
      e <- switch(input$bzd,
                  lzp4  = d_lzp(tb, 4),
                  lzp2  = d_lzp(tb, 2),
                  mdzim = d_mdz_im(tb, 10),
                  mdziv = d_mdz_iv(tb, 0.2),
                  dzp   = d_dzp(tb, 10))
      if (input$bzd_rep) {
        e2 <- switch(input$bzd,
                     lzp4  = d_lzp(tb + 4, 4),
                     lzp2  = d_lzp(tb + 4, 2),
                     mdzim = d_mdz_im(tb + 4, 10),
                     mdziv = d_mdz_iv(tb + 4, 0.2),
                     dzp   = d_dzp(tb + 4, 10))
        e <- e + e2
      }
    }
    if (input$second != "none") {
      s2 <- switch(input$second,
                   lev = d_lev(input$t_second, 60),
                   fos = d_fos(input$t_second, 20),
                   vpa = d_vpa(input$t_second, 40),
                   pb  = d_pb(input$t_second, 20))
      e <- add(s2)
    }
    tt <- input$t_third
    dur <- max(30, input$tend - tt)
    if ("mdzinf" %in% input$third)
      e <- add(d_mdz_iv(tt, 0.2) + d_mdz_inf(tt, dur, 0.4))
    if ("prop" %in% input$third)
      e <- add(d_pro(tt, 2) + d_pro_inf(tt, dur, 6))
    if ("ket" %in% input$third)
      e <- add(d_ket(tt, 2) + d_ket_inf(tt, dur, 3))
    if ("allo" %in% input$third)
      e <- add(d_allo_inf(tt, dur, 86))
    if (is.null(e)) e <- ev(time = 0, amt = 0, cmt = 1)
    e
  })

  sim <- reactive({
    p <- list(EDRIVE = as.numeric(input$aetiology),
              ALB = input$alb, CRCL = input$crcl,
              ANAK = as.numeric(input$anak),
              GLUCIV = as.numeric(input$gluciv),
              COOL = as.numeric(input$cool),
              TDRIVE = if (input$reversible) 30 else 1e6,
              DRESID = 0.10)
    param(mod, p) %>%
      mrgsim_df(events = build_events(), end = input$tend, delta = 0.5,
                maxsteps = 2000000)
  })

  # untreated counterfactual with the same aetiology — the reference clock
  sim0 <- reactive({
    param(mod, EDRIVE = as.numeric(input$aetiology)) %>%
      mrgsim_df(events = ev(time = 0, amt = 0, cmt = 1),
                end = max(input$tend, 720), delta = 1, maxsteps = 2000000)
  })

  tstop <- reactive({
    o <- sim(); n5 <- 10
    r <- rle(o$SEIZ < 0.05)
    i <- which(r$values & r$lengths >= n5)
    if (!length(i)) return(NA_real_)
    pos <- if (i[1] == 1) 1 else sum(r$lengths[1:(i[1] - 1)]) + 1
    o$time[pos]
  })

  # ---------------- verdict banner ----------------
  output$verdict <- renderUI({
    o <- sim(); ts <- tstop()
    loss <- max(o$NEURLOSS); burden <- max(o$TSEIZ)
    if (is.na(ts)) {
      HTML(sprintf("<b style='color:#B71C1C'>Seizure continues (not terminated within %.0f min)</b> &nbsp;·&nbsp; cumulative seizure time %.0f min &nbsp;·&nbsp; hippocampal neuron loss %.1f%% &nbsp;·&nbsp; lowest MAP %.0f mmHg",
                   input$tend, burden, loss, min(o$MAP)))
    } else {
      HTML(sprintf("<b style='color:#1B5E20'>Seizure terminated @ %.1f min</b> &nbsp;·&nbsp; cumulative seizure time %.1f min &nbsp;·&nbsp; hippocampal neuron loss %.1f%% &nbsp;·&nbsp; lowest MAP %.0f mmHg &nbsp;·&nbsp; intubation %s",
                   ts, burden, loss, min(o$MAP),
                   ifelse(any(o$INTUB > 0.5), "required", "not required")))
    }
  })

  # ---------------- ① gain ----------------
  output$p_gain <- renderPlot({
    o <- sim()
    d <- o %>% select(time, oEGAIN, oINH_TOT, GNET) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) +
      geom_hline(yintercept = 1, linetype = 2, colour = "grey40") +
      geom_line(linewidth = 1) +
      facet_wrap(~ name, scales = "free_y", ncol = 3,
                 labeller = as_labeller(c(oEGAIN = "Excitatory gain E",
                                          oINH_TOT = "Inhibitory tone I",
                                          GNET = "G = E / I"))) +
      scale_colour_manual(values = PALETTE) +
      labs(x = "Time (min)", y = NULL) + theme_se + theme(legend.position = "none")
  })

  output$t_thresholds <- renderTable({
    o0 <- sim0()
    i <- function(tt) which.min(abs(o0$time - tt))
    tv <- c(5, 15, 30, 60, 120, 240)
    data.frame(`Time (min)` = tv,
               `Untreated G` = round(o0$GNET[sapply(tv, i)], 2),
               `Benzo target pool` = round(o0$oTGT_BZD[sapply(tv, i)], 3),
               `NMDA conductance` = round(o0$oNMDA[sapply(tv, i)], 3),
               `Required benzo concentration (x EC50)` =
                 ifelse(o0$CREQ_BZD[sapply(tv, i)] < 0, "Not possible",
                        sprintf("%.2f", o0$CREQ_BZD[sapply(tv, i)])),
               check.names = FALSE)
  }, digits = 3)

  # ---------------- ② timeline ----------------
  output$p_timeline <- renderPlot({
    o <- sim()
    d <- o %>% select(time, EEG = EEGOUT, Motor = MOTOROUT,
                      `G (gain)` = GNET, `Benzo pool` = oTGT_BZD,
                      `NMDA cond.` = oNMDA, `Injury (%)` = NEURLOSS) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) +
      geom_line(linewidth = 1) +
      facet_wrap(~ name, scales = "free_y", ncol = 2) +
      scale_colour_manual(values = PALETTE) +
      labs(x = "Time (min)", y = NULL, title = "One simulation, six reading windows") +
      theme_se + theme(legend.position = "none")
  })

  output$t_summary <- renderTable({
    o <- sim()
    data.frame(
      Item = c("Seizure termination time (min)", "Cumulative seizure time (min)",
               "Maximum intracellular Cl- (mM)", "60-min R_SYN", "60-min NMDA conductance",
               "Hippocampal neuron loss (%)", "Peak CK (U/L)", "Lowest MAP (mmHg)",
               "Lowest glucose (mmol/L)", "Peak temperature (°C)", "Intubation", "Epileptogenesis burden"),
      Value = c(ifelse(is.na(tstop()), "Not terminated", sprintf("%.1f", tstop())),
             sprintf("%.1f", max(o$TSEIZ)),
             sprintf("%.1f", max(o$CLI)),
             sprintf("%.3f", o$RSYN[which.min(abs(o$time - 60))]),
             sprintf("%.3f", o$oNMDA[which.min(abs(o$time - 60))]),
             sprintf("%.1f", max(o$NEURLOSS)),
             sprintf("%.0f", max(o$CK)),
             sprintf("%.0f", min(o$MAP)),
             sprintf("%.2f", min(o$GLUCP)),
             sprintf("%.2f", max(o$TEMP)),
             ifelse(any(o$INTUB > 0.5), "Required", "Not required"),
             sprintf("%.2f", max(o$EPG))))
  })

  # ---------------- ③ receptor clock ----------------
  output$p_clock <- renderPlot({
    o <- sim()
    d <- o %>% transmute(time,
                         `R_SYN × F_BZS (benzo target)` = oTGT_BZD,
                         `R_SYN (synaptic GABA-A)` = RSYN,
                         `R_EXTRA (extrasynaptic δ)` = REXTRA,
                         `NR_SYN × Mg release (NMDA)` = oNMDA) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) +
      geom_hline(yintercept = 1, linetype = 3, colour = "grey60") +
      geom_line(linewidth = 1.2) +
      scale_colour_manual(values = c("#1E88E5", "#64B5F6", "#1B5E20", "#E65100")) +
      labs(x = "Time (min)", y = "Fraction of baseline",
           colour = NULL, title = "One clock, two receptors moving in opposite directions") +
      theme_se
  })

  output$p_pools <- renderPlot({
    o <- sim()
    d <- o %>% transmute(time, `KCC2` = KCC2, `[Cl-]i / 20 mM` = CLI / 20,
                         `Cl coefficient (CLFAC)` = oCLFAC,
                         `Neuropeptide balance` = NETPEP,
                         `Adenosine` = ADO) %>% pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 1) +
      scale_colour_manual(values = PALETTE[c(6, 4, 1, 5, 3)]) +
      labs(x = "Time (min)", y = NULL, colour = NULL,
           title = "The remaining three pathways that erode inhibition") + theme_se
  })

  # ---------------- ④ PK ----------------
  output$p_pk <- renderPlot({
    o <- sim()
    d <- o %>% transmute(time,
                         `lorazepam (ng/mL)` = CP_LZP,
                         `midazolam (ng/mL)` = CP_MDZ,
                         `diazepam (ng/mL)` = CP_DZP,
                         `levetiracetam (mg/L)` = CP_LEV,
                         `phenytoin total (mg/L)` = CP_PHT,
                         `valproate total (mg/L)` = CP_VPA,
                         `phenobarbital (mg/L)` = CP_PB,
                         `ketamine (ng/mL)` = CP_KET,
                         `propofol (mg/L)` = CP_PRO,
                         `allopregnanolone (ng/mL)` = CP_ALLO) %>%
      pivot_longer(-time) %>% group_by(name) %>% filter(max(value) > 1e-6) %>% ungroup()
    if (!nrow(d)) return(NULL)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 1) +
      facet_wrap(~ name, scales = "free_y") +
      scale_colour_manual(values = PALETTE) +
      labs(x = "Time (min)", y = "Plasma concentration") + theme_se +
      theme(legend.position = "none")
  })

  output$p_occ <- renderPlot({
    o <- sim()
    d <- o %>% transmute(time, benzodiazepine = oOCC_BZD, ketamine = oOCC_KET,
                         levetiracetam = oOCC_LEV, phenytoin = oOCC_PHT,
                         valproate = oOCC_VPA, phenobarbital = oOCC_PB,
                         propofol = oOCC_PRO) %>%
      pivot_longer(-time) %>% group_by(name) %>% filter(max(value) > 1e-4) %>% ungroup()
    if (!nrow(d)) return(NULL)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 1) +
      scale_colour_manual(values = PALETTE) + ylim(0, 1) +
      labs(x = "Time (min)", y = "Site occupancy", colour = NULL,
           title = "The drug determines occupancy; the effect is occupancy × remaining target") +
      theme_se
  })

  # ---------------- ⑤ price of a minute ----------------
  output$p_creq <- renderPlot({
    o0 <- sim0()
    d <- o0 %>% filter(time <= 240) %>%
      transmute(time,
                benzodiazepine = ifelse(CREQ_BZD < 0, NA, CREQ_BZD),
                phenobarbital  = ifelse(CREQ_PB  < 0, NA, CREQ_PB)) %>%
      pivot_longer(-time)
    tinf <- o0$time[which(o0$CREQ_BZD < 0)[1]]
    ggplot(d, aes(time, value, colour = name)) +
      geom_line(linewidth = 1.2, na.rm = TRUE) +
      { if (!is.na(tinf)) geom_vline(xintercept = tinf, linetype = 2, colour = "#B71C1C") } +
      { if (!is.na(tinf)) annotate("text", x = tinf, y = Inf, vjust = 1.4, hjust = -0.05,
                                   label = sprintf("No benzo dose is possible past: %.0f min", tinf),
                                   colour = "#B71C1C", size = 4) } +
      scale_y_log10() +
      scale_colour_manual(values = c("#1E88E5", "#E65100")) +
      labs(x = "Time (min)", y = "Effect-site concentration needed to make G=1 (x EC50, log)",
           colour = NULL, title = "The price of a minute is not constant") + theme_se
  })

  output$t_creq <- renderTable({
    o0 <- sim0(); tv <- c(5, 10, 20, 30, 45, 60, 120, 240)
    i <- sapply(tv, function(tt) which.min(abs(o0$time - tt)))
    data.frame(`Time (min)` = tv,
               `R_SYN` = round(o0$RSYN[i], 3),
               `F_BZS` = round(o0$FBZS[i], 3),
               `Benzo target` = round(o0$oTGT_BZD[i], 3),
               `Required benzo (x EC50)` = ifelse(o0$CREQ_BZD[i] < 0, "Not possible", sprintf("%.2f", o0$CREQ_BZD[i])),
               `Required phenobarbital (x EC50)` = ifelse(o0$CREQ_PB[i] < 0, "Not possible", sprintf("%.2f", o0$CREQ_PB[i])),
               check.names = FALSE)
  })

  # ---------------- ⑥ crossover ----------------
  cross <- reactive({
    o <- sim0(); p <- as.list(param(mod)); occ <- 0.75
    tv <- seq(5, min(720, max(o$time)), by = 5)
    i <- sapply(tv, function(tt) which.min(abs(o$time - tt)))
    inh0 <- o$oINH_TOT[i]; eg0 <- o$oEGAIN[i]; g0 <- eg0 / inh0
    tgt_b <- o$oTGT_BZD[i]
    tgt_p <- p$FSYNPB * o$RSYN[i] + (1 - p$FSYNPB) * o$REXTRA[i]
    blk <- p$EMAXKET * occ * (p$UDMIN + (1 - p$UDMIN) * o$SEIZ[i])
    data.frame(time = tv,
               benzodiazepine = 1 - (eg0 / (inh0 + p$EMAXBZD * occ * tgt_b)) / g0,
               phenobarbital  = 1 - (eg0 / (inh0 + p$EMAXPB  * occ * tgt_p)) / g0,
               ketamine       = blk * o$KETSHARE[i])
  })

  output$p_cross <- renderPlot({
    d <- cross(); dl <- pivot_longer(d, -time)
    xc <- d$time[which(d$ketamine > d$benzodiazepine)[1]]
    ggplot(dl, aes(time, value, colour = name)) + geom_line(linewidth = 1.2) +
      { if (!is.na(xc)) geom_vline(xintercept = xc, linetype = 2, colour = "grey30") } +
      { if (!is.na(xc)) annotate("text", x = xc, y = 0.05, hjust = -0.05,
                                 label = sprintf("Crossover: %.0f min", xc), size = 4) } +
      scale_colour_manual(values = c("#1E88E5", "#E65100", "#1B5E20")) +
      labs(x = "Time (min)", y = "Reduction in G at the same occupancy (0.75)", colour = NULL,
           title = "Comparing target alone — one departs while one arrives") + theme_se
  })

  output$t_cross <- renderTable({
    d <- cross(); tv <- c(5, 30, 60, 120, 240, 480)
    d <- d[sapply(tv, function(tt) which.min(abs(d$time - tt))), ]
    data.frame(`Time (min)` = d$time,
               `Benzo` = round(d$benzodiazepine, 3),
               `Phenobarbital` = round(d$phenobarbital, 3),
               `Ketamine` = round(d$ketamine, 3), check.names = FALSE)
  })

  # ---------------- ⑦ dissociation ----------------
  output$p_diss <- renderPlot({
    o <- sim()
    d <- o %>% transmute(time, `EEG burden (electrographic)` = EEGOUT,
                         `Motor expression (motor)` = MOTOROUT,
                         `Motor gain MOTG` = MOTG,
                         `Benzodiazepine occupancy` = oOCC_BZD) %>% pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 1.2) +
      geom_hline(yintercept = 0.15, linetype = 3, colour = "grey50") +
      scale_colour_manual(values = c("#880E4F", "#B71C1C", "#546E7A", "#1E88E5")) +
      labs(x = "Time (min)", y = NULL, colour = NULL,
           title = "What the bedside sees and what the brain is doing") + theme_se
  })

  output$t_diss <- renderTable({
    o <- sim(); tv <- c(15, 30, 60, 90, 120, 180, 240)
    tv <- tv[tv <= max(o$time)]
    i <- sapply(tv, function(tt) which.min(abs(o$time - tt)))
    data.frame(`Time (min)` = tv,
               `Motor output` = round(o$MOTOROUT[i], 3),
               `EEG burden` = round(o$EEGOUT[i], 3),
               `Bedside verdict` = ifelse(o$MOTOROUT[i] < 0.15, "Appears controlled", "Seizure continues"),
               `Actual` = ifelse(o$EEGOUT[i] < 0.05, "Controlled", "Seizure continues"),
               check.names = FALSE)
  })

  # ---------------- ⑧ systemic ----------------
  output$p_sys <- renderPlot({
    o <- sim()
    d <- o %>% transmute(time, `MAP (mmHg)` = MAP, `CBF (relative)` = oCBF,
                         `Supply-demand mismatch` = oMISMATCH,
                         `Glucose (mmol/L)` = GLUCP, `Lactate (mmol/L)` = LAC,
                         `Temperature (°C)` = TEMP, `CK (U/L)` = CK,
                         `Autoregulation AUTO` = AUTO, `Respiratory burden` = RESPD) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 1) +
      facet_wrap(~ name, scales = "free_y", ncol = 3) +
      scale_colour_manual(values = PALETTE) +
      labs(x = "Time (min)", y = NULL) + theme_se + theme(legend.position = "none")
  })

  # ---------------- ⑨ injury ----------------
  output$p_injury <- renderPlot({
    o <- sim()
    d <- o %>% transmute(time, `Glutamate (uM)` = GLU, `[Ca2+]i` = CAI,
                         `Cumulative injury` = INJURY, `Hippocampal neuron loss (%)` = NEURLOSS,
                         `IL-1beta` = IL1B, `Epileptogenesis` = EPG) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 1) +
      facet_wrap(~ name, scales = "free_y", ncol = 3) +
      scale_colour_manual(values = PALETTE) +
      labs(x = "Time (min)", y = NULL) + theme_se + theme(legend.position = "none")
  })

  output$t_injury <- renderTable({
    o <- sim()
    data.frame(
      Item = c("Cumulative seizure time (min)", "Peak [Ca2+]i", "Cumulative injury (a.u.)",
               "Hippocampal neuron loss (%)", "Peak IL-1beta", "BBB permeability fold",
               "P-gp induction fold", "Brain:plasma partition index", "Epileptogenesis burden"),
      Value = round(c(max(o$TSEIZ), max(o$CAI), max(o$INJURY), max(o$NEURLOSS),
                   max(o$IL1B), max(o$BBBP), max(o$PGP), min(o$BPRATIO),
                   max(o$EPG)), 3))
  })

  # ---------------- ⑩ virtual population ----------------
  popres <- eventReactive(input$runpop, {
    n <- input$npop
    set.seed(20260804)
    id <- data.frame(ID = 1:n, EDRIVE = input$meddrive * exp(rnorm(n, 0, input$sddrive)))
    o1 <- mod %>% mrgsim_df(events = d_lzp(20, 4), idata = id, end = 100, delta = 1,
                            maxsteps = 2000000, obsonly = TRUE)
    f1 <- o1 %>% group_by(ID) %>%
      summarise(ok = any(SEIZ[time >= 20 & time <= 80] < 0.05), .groups = "drop")
    id2 <- id[id$ID %in% f1$ID[!f1$ok], ]
    rate2 <- function(e) {
      if (!nrow(id2)) return(NA_real_)
      o <- mod %>% mrgsim_df(events = e, idata = id2, end = 150, delta = 1,
                             maxsteps = 2000000, obsonly = TRUE)
      100 * mean((o %>% group_by(ID) %>%
                    summarise(ok = any(SEIZ[time >= 65 & time <= 125] < 0.05),
                              .groups = "drop"))$ok)
    }
    im <- mod %>% mrgsim_df(events = d_mdz_im(8, 10), idata = id, end = 100, delta = 1,
                            maxsteps = 2000000, obsonly = TRUE)
    iv <- mod %>% mrgsim_df(events = d_lzp(12.5, 4), idata = id, end = 100, delta = 1,
                            maxsteps = 2000000, obsonly = TRUE)
    rr <- function(o, a, b) 100 * mean((o %>% group_by(ID) %>%
      summarise(ok = any(SEIZ[time >= a & time <= b] < 0.05), .groups = "drop"))$ok)
    list(
      tab = data.frame(
        endpoint = c("First-line benzodiazepine (20 min)", "ESETT levetiracetam",
                     "ESETT fosphenytoin", "ESETT valproate",
                     "RAMPART IM midazolam", "RAMPART IV lorazepam"),
        model_pct = round(c(100 * mean(f1$ok),
                            rate2(d_lzp(20, 4) + d_lev(65, 60)),
                            rate2(d_lzp(20, 4) + d_fos(65, 20)),
                            rate2(d_lzp(20, 4) + d_vpa(65, 40)),
                            rr(im, 8, 68), rr(iv, 12.5, 72.5)), 1),
        observed_pct = c("59-73", "47", "45", "46", "73.4", "63.4")),
      first = do.call(rbind, lapply(c(5, 10, 20, 30, 45, 60), function(tt) {
        o <- mod %>% mrgsim_df(events = d_lzp(tt, 4), idata = id, end = tt + 70,
                               delta = 1, maxsteps = 2000000, obsonly = TRUE)
        data.frame(line = "1st line (benzodiazepine)", t = tt, pct = rr(o, tt, tt + 60))
      })),
      second = if (nrow(id2)) do.call(rbind, lapply(c(30, 45, 60, 90, 120), function(tt) {
        o <- mod %>% mrgsim_df(events = d_lzp(20, 4) + d_lev(tt, 60), idata = id2,
                               end = tt + 70, delta = 1, maxsteps = 2000000, obsonly = TRUE)
        data.frame(line = "2nd line (levetiracetam, benzo failures)", t = tt, pct = rr(o, tt, tt + 60))
      })) else NULL)
  })

  output$t_pop <- renderTable({ popres()$tab })

  output$p_pop <- renderPlot({
    r <- popres(); d <- rbind(r$first, r$second)
    ggplot(d, aes(t, pct, colour = line)) +
      geom_line(linewidth = 1.2) + geom_point(size = 2.5) +
      scale_colour_manual(values = c("#1E88E5", "#1B5E20")) + ylim(0, 100) +
      labs(x = "Administration time (min)", y = "Response rate (%)", colour = NULL,
           title = "Two cliffs — and the 40 minutes between them") + theme_se
  })
}

shinyApp(ui, server)
