## =============================================================================
##  dlb_shiny_app.R
##  Dementia with Lewy Bodies — interactive QSP dashboard
##
##  Run with:
##      setwd("dementia-with-lewy-bodies")
##      shiny::runApp("dlb_shiny_app.R")
##
##  Requires: shiny, ggplot2, dplyr, tidyr, mrgsolve  (and dlb_mrgsolve_model.R
##  in the same directory).
##
##  ---------------------------------------------------------------------------
##  WHAT THIS APP IS FOR
##  ---------------------------------------------------------------------------
##  Not "here are ten line plots of a simulation".  Each tab is built around one
##  claim the model makes, and is laid out so that the claim can be CHECKED --
##  usually by putting the thing that is supposed to explain the effect directly
##  next to the effect.
##
##    Tab 3  puts the three transducers side by side, with the receptor-density
##           equation for each, so the sign difference is visible rather than
##           asserted.
##    Tab 4  plots the attention cubic's potential landscape next to the
##           measured fluctuation score, so "fluctuation is a variance" is a
##           picture and not a sentence.
##    Tab 5  decomposes the hallucination drive into its three factors and shows
##           what happens to the product when you knock out one at a time.
##    Tab 6  overlays D2 occupancy, postsynaptic reserve, and the resulting
##           deficit for whichever phenotype you have selected, so neuroleptic
##           sensitivity can be traced from cause to consequence in one screen.
##    Tab 7  lets you switch ambroxol on at any day and watch the GBA1 loop.
##
##  Licence: see repository LICENSE.  EDUCATIONAL / RESEARCH USE ONLY.
## =============================================================================

library(shiny)
library(mrgsolve)
library(ggplot2)
library(dplyr)
library(tidyr)

source("dlb_mrgsolve_model.R")

## ----------------------------------------------------------------- theming --
qsp_theme <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.text  = element_text(face = "bold"),
        plot.title  = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(colour = "grey35", size = 10),
        legend.position = "bottom",
        legend.title = element_blank())

PAL <- c(DLB = "#1f4a7a", `PD-dementia` = "#b8860b", AD = "#7a2f4a",
         treated = "#1f7a4a", untreated = "#8a8a99")

yr <- function(d) d$time / 365.25

## Every plot in the app goes through this, so that the x axis, the theme and
## the "years since diagnosis" convention are identical everywhere.
trace_plot <- function(d, vars, labels = vars, title = "", subtitle = "",
                       ylab = "", hlines = NULL) {
  long <- d %>%
    select(time, scenario, all_of(vars)) %>%
    pivot_longer(all_of(vars), names_to = "var", values_to = "y") %>%
    mutate(var = factor(var, levels = vars, labels = labels),
           yrs = time / 365.25)
  p <- ggplot(long, aes(yrs, y, colour = scenario, linetype = scenario)) +
    geom_line(linewidth = 0.75) +
    facet_wrap(~var, scales = "free_y") +
    labs(title = title, subtitle = subtitle,
         x = "Years since diagnosis", y = ylab) +
    qsp_theme
  if (!is.null(hlines))
    p <- p + geom_hline(yintercept = hlines, linetype = "dotted",
                        colour = "grey50")
  p
}

## =============================================================================
##  UI
## =============================================================================
ui <- fluidPage(
  titlePanel("Dementia with Lewy Bodies QSP Dashboard"),
  tags$p(style = "color:#555;margin-top:-8px;",
         "66-compartment mrgsolve model · three-transmitter presynaptic/postsynaptic dissociation · ",
         tags$b("educational / research use only")),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      h4("① Patient"),
      selectInput("pheno", "Phenotype",
                  choices = c("DLB" = 0, "PD-dementia" = 1, "AD" = 2),
                  selected = 0),
      helpText(tags$small(
        "The phenotype is not a label but the ", tags$b("initial distribution of pathology"),
        ". Only the starting values of brainstem/limbic/neocortical α-synuclein and of tau·amyloid change, and the ",
        "66 equations are completely identical.")),
      sliderInput("gbaf", "GBA1 functional allele dosage (1 = WT)",
                  min = 0.30, max = 1.00, value = 1.00, step = 0.05),
      sliderInput("apoe4", "Number of APOE ε4 alleles", min = 0, max = 2, value = 0, step = 1),
      sliderInput("snca", "SNCA gene dosage", min = 0.8, max = 1.6, value = 1.0, step = 0.05),
      sliderInput("antich", "Anticholinergic burden",
                  min = 0, max = 3, value = 0, step = 1),

      hr(), h4("② Treatment"),
      selectInput("chei", "Cholinesterase inhibitor",
                  choices = c("None" = "none",
                              "Rivastigmine 6 mg BID oral" = "riv_oral",
                              "Rivastigmine 9.5 mg/24h patch" = "riv_patch",
                              "Donepezil 10 mg qd" = "don")),
      numericInput("chei_start", "ChEI start day (day)", value = 180, min = 0, max = 2000, step = 30),
      checkboxInput("mem",  "Memantine 20 mg qd", FALSE),
      checkboxInput("pim",  "Pimavanserin 34 mg qd (day 365~)", FALSE),
      selectInput("apd", "Antipsychotic (day 730~, 180 days)",
                  choices = c("None" = "none",
                              "Risperidone 1 mg qd" = "risp",
                              "Quetiapine 50 mg qd" = "quet")),
      checkboxInput("ld",   "Levodopa/carbidopa 150 mg TID (day 365~)", FALSE),
      checkboxInput("zon",  "Zonisamide 25 mg qd (day 365~)", FALSE),
      checkboxInput("amb",  "Ambroxol 1.26 g/day", FALSE),
      numericInput("amb_start", "Ambroxol start day (day)", value = 0, min = 0, max = 2500, step = 30),
      checkboxInput("mab",  "Anti-α-synuclein antibody 4500 mg IV q4w", FALSE),
      checkboxInput("mel",  "Melatonin/clonazepam (RBD)", FALSE),
      checkboxInput("drox", "Droxidopa/midodrine (orthostatic hypotension)", FALSE),
      checkboxInput("mod",  "Modafinil (daytime sleepiness)", FALSE),

      hr(), h4("③ Simulation"),
      sliderInput("years", "Simulation duration (years)", min = 2, max = 10, value = 8, step = 1),
      checkboxInput("show_ref", "Show the untreated control alongside", TRUE),
      actionButton("go", "Run", class = "btn-primary btn-block"),
      hr(),
      tags$small(tags$b("Note."), " This is a semi-quantitative teaching model assembled from the public literature, and ",
                 "it must not be used for clinical decisions, prescribing, or regulatory submission.")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",

        ## ------------------------------------------------------- 1. profile
        tabPanel("1 · Patient profile",
                 br(), h4("State at diagnosis"),
                 tableOutput("tbl_profile"),
                 h4("Core & indicative features"),
                 plotOutput("p_profile", height = 420),
                 helpText("MIBG H/M < 2.0, a reduced DaTSCAN binding ratio and an EEG dominant frequency of 5.6–7.9 Hz (pre-alpha) are ",
                          "the indicative biomarkers of the fourth consensus criteria. All three values are computed by the model, ",
                          "not constants specified separately.")),

        ## ------------------------------------------------------------ 2. PK
        tabPanel("2 · Pharmacokinetics (PK)",
                 br(), plotOutput("p_pk", height = 480),
                 h4("Rivastigmine: the plasma half-life and the pharmacodynamic half-life come apart"),
                 plotOutput("p_carb", height = 320),
                 helpText("The plasma half-life of rivastigmine is about 1.5 hours, but recovery of the carbamylated enzyme is ",
                          "set by the rate of enzyme resynthesis (about 9 hours). The model holds the carbamylated enzyme as ",
                          "a separate state variable, and that is why twelve-hourly dosing works. ",
                          "The patch achieves the same AChE inhibition at a much lower Cmax, and because ",
                          "the gastrointestinal adverse-effect index is driven by Cmax, it comes out lower on the patch.")),

        ## ------------------------------------------- 3. the three transducers
        tabPanel("3 · The three transducers ★",
                 br(),
                 div(style = "background:#fffbe6;border-left:4px solid #b8860b;padding:10px;",
                     tags$b("This tab is the central claim of the model."), br(),
                     "All three ascending transmitter systems pass through the same three-factor product — ",
                     tags$code("presynaptic capacity × postsynaptic density × (1 − occupancy)"), ". ",
                     "The only difference is the differential equation of the middle factor, and its sign differs from system to system."),
                 br(),
                 plotOutput("p_transducer", height = 300),
                 h4("Postsynaptic receptor density — three equations with different signs"),
                 plotOutput("p_receptor", height = 300),
                 tableOutput("tbl_transducer")),

        ## --------------------------------------------- 4. cognition & CAF
        tabPanel("4 · Cognition and fluctuation ★",
                 br(), plotOutput("p_cog", height = 340),
                 h4("The attention potential landscape"),
                 plotOutput("p_potential", height = 340),
                 helpText("Left: the cubic potential at the selected time point. Two wells means the system is inside the bistable region, ",
                          "where state switching occurs and is therefore scored in clinic as 'cognitive fluctuation'. ",
                          "Right: the bistability depth (BISTAB) and the actual CAF score. ",
                          "A ChEI raises the drive and pushes the patient out of this band, ",
                          "which is why it moves the variance (CAF) sooner, and further, than the mean (MMSE).")),

        ## -------------------------------------------------- 5. hallucinations
        tabPanel("5 · Visual hallucinations (the PAD model)",
                 br(), plotOutput("p_vh", height = 340),
                 h4("Decomposition of the three-factor product"),
                 plotOutput("p_pad", height = 340),
                 helpText("The hallucination drive is not a sum but a ", tags$b("product"), ": ",
                          tags$code("(1−bottom-up evidence fidelity) × (1−top-down attentional binding) × 5-HT2A gain"), ". ",
                          "Because it is a product, no single factor can either explain the hallucinations or abolish them — ",
                          "which is why pimavanserin does not abolish hallucinations but only reduces them by about a third.")),

        ## ---------------------------------- 6. motor & neuroleptic sensitivity
        tabPanel("6 · Motor features and neuroleptic sensitivity ★",
                 br(), plotOutput("p_motor", height = 340),
                 h4("The causal chain of the sensitivity on one screen (cause → consequence)"),
                 plotOutput("p_nsens", height = 360),
                 helpText("D2 occupancy is identical in all three phenotypes. The only thing that differs is ",
                          tags$b("the postsynaptic reserve (RESERVE)"), ", and ",
                          "the drug-induced deficit (DEFICIT) is amplified the smaller that reserve is. ",
                          "Off drug the occupancy is zero, so the deficit is zero too — ",
                          "that is, the sensitivity reaction is not a milestone of disease severity but a ", tags$b("drug event"), ".")),

        ## ---------------------------------------------- 7. pathology & GBA1
        tabPanel("7 · Pathological progression and GBA1",
                 br(), plotOutput("p_path", height = 340),
                 h4("The GBA1 positive feedback loop"),
                 plotOutput("p_gba", height = 340),
                 helpText("GCase↓ → GlcCer↑ → delayed oligomer clearance → oligomers↑ → blocked GCase trafficking → GCase↓. ",
                          "At the default parameters the gain of this loop is less than one, so it is ",
                          tags$b("monostable"), " — stop the ambroxol and it comes back. ",
                          "What this model actually supports is therefore not 'miss the deadline and it is over' but ",
                          "'the earlier you pull, the longer the lever'. ",
                          "Try varying the ambroxol start day and see.")),

        ## ------------------------------------------------------ 8. biomarkers
        tabPanel("8 · Biomarkers",
                 br(), plotOutput("p_bio", height = 480),
                 tableOutput("tbl_bio")),

        ## ------------------------------------------------ 9. scenario compare
        tabPanel("9 · Scenario comparison",
                 br(),
                 helpText("Runs all 23 scenarios defined in the model file and compares them. ",
                          "The first run takes around 30 seconds."),
                 actionButton("go_all", "Run all scenarios", class = "btn-warning"),
                 br(), br(),
                 tableOutput("tbl_all"),
                 plotOutput("p_all", height = 460)),

        ## --------------------------------------------------------- 10. survival
        tabPanel("10 · Survival and hazard",
                 br(), plotOutput("p_surv", height = 400),
                 tableOutput("tbl_surv"),
                 helpText("The hazard function is log-linear in cognitive loss, motor score and fall index, and ",
                          "neuroleptic sensitivity enters ", tags$b("multiplicatively"), ". ",
                          "It is exactly this multiplicative structure that reproduces the signal of a 2–3 fold ",
                          "mortality after antipsychotic exposure.")),

        ## ------------------------------------------------------ 11. model card
        tabPanel("11 · Model card",
                 br(), htmlOutput("model_card"))
      )
    )
  )
)

## =============================================================================
##  SERVER
## =============================================================================
server <- function(input, output, session) {

  ## ------------------------------------------------------------ dose builder
  build_dose <- function() {
    end <- input$years * 365
    ev_list <- list()
    add <- function(e) ev_list[[length(ev_list) + 1]] <<- e

    st <- input$chei_start
    if (input$chei == "riv_oral")  add(riv_oral(6,   st, max(end - st, 1)))
    if (input$chei == "riv_patch") add(riv_patch(9.5, st, max(end - st, 1)))
    if (input$chei == "don")       add(don_oral(10,  st, max(end - st, 1)))
    if (input$mem)  add(mem_oral(20, st, max(end - st, 1)))
    if (input$pim)  add(pim_oral(34, 365, max(end - 365, 1)))
    if (input$ld)   add(ld_oral(150, 365, max(end - 365, 1)))
    if (input$zon)  add(zon_oral(25, 365, max(end - 365, 1)))
    if (input$amb)  add(amb_oral(420, input$amb_start,
                                 max(end - input$amb_start, 1)))
    if (input$mab)  add(mab_iv(4500, 180, max(floor((end - 180)/28), 1)))
    if (input$apd == "risp") add(apd_oral(1,  730, 180))
    if (input$apd == "quet") add(apd_oral(50, 730, 180))

    if (!length(ev_list)) return(NULL)
    Reduce(c, ev_list)
  }

  build_param <- function() {
    p <- list(PHENO  = as.numeric(input$pheno),
              GBAF   = input$gbaf,
              APOE4  = input$apoe4,
              SNCADOSE = input$snca,
              ANTICH = input$antich,
              MELON  = as.numeric(input$mel),
              DROXON = as.numeric(input$drox),
              MODON  = as.numeric(input$mod))
    if (input$apd == "quet")
      p <- c(p, list(EC50D2 = 900, CLAPD = 1400, VAPD = 700,
                     APDA1 = 1.0, APDSED = 1.0))
    p
  }

  sim <- eventReactive(input$go, {
    end <- input$years * 365
    m   <- param(mod, build_param())
    d   <- build_dose()
    trt <- if (is.null(d)) mrgsim_df(m, end = end, delta = 1)
           else            mrgsim_df(m, events = d, end = end, delta = 1)
    trt$scenario <- "treated"
    if (isTRUE(input$show_ref)) {
      ref <- mrgsim_df(param(mod, build_param()), end = end, delta = 1)
      ref$scenario <- "untreated"
      rbind(ref, trt)
    } else trt
  }, ignoreNULL = FALSE)

  ph_label <- reactive(c("0" = "DLB", "1" = "PD-dementia", "2" = "AD")[input$pheno])

  ## ----------------------------------------------------------- 1. profile --
  output$tbl_profile <- renderTable({
    d <- sim(); d0 <- d[d$time == 0 & d$scenario == d$scenario[1], ][1, ]
    data.frame(
      `Item` = c("Phenotype", "MMSE", "MDS-UPDRS III", "CAF (cognitive fluctuation)",
               "NPI total", "Hallucination burden", "RBD severity", "Orthostatic SBP fall (mmHg)",
               "Epworth sleepiness", "MIBG H/M", "DaTSCAN binding ratio", "EEG dominant frequency (Hz)",
               "Neocortical fibril load", "Limbic fibril load", "Brainstem fibril load"),
      `Value` = c(ph_label(),
             sprintf("%.1f", d0$MMSE), sprintf("%.1f", d0$MOT), sprintf("%.1f", d0$CAF),
             sprintf("%.1f", d0$NPI), sprintf("%.2f", d0$VHB), sprintf("%.1f", d0$RBDS),
             sprintf("%.1f", d0$AUTS), sprintf("%.1f", d0$EDSS), sprintf("%.2f", d0$MIBG),
             sprintf("%.2f", d0$DATSBR), sprintf("%.2f", d0$EEGF),
             sprintf("%.3f", d0$FIBN), sprintf("%.3f", d0$FIBL), sprintf("%.3f", d0$FIBB)),
      check.names = FALSE)
  })

  output$p_profile <- renderPlot({
    trace_plot(sim(), c("MIBG", "DATSBR", "EEGF", "RBDS", "AUTS", "EDSS"),
               c("MIBG H/M", "DaTSCAN binding ratio", "EEG dominant frequency (Hz)",
                 "RBD severity", "Orthostatic SBP fall (mmHg)", "Epworth"),
               title = "Time course of the indicative biomarkers and the supportive features",
               subtitle = paste0("Phenotype: ", ph_label())) +
      scale_colour_manual(values = PAL)
  })

  ## ---------------------------------------------------------------- 2. PK --
  output$p_pk <- renderPlot({
    d <- sim(); d <- d[d$scenario == "treated", ]
    trace_plot(d, c("CRIV", "CDON", "CPIM", "CLD", "CAPD", "CAMBB"),
               c("Rivastigmine (ng/mL)", "Donepezil (ng/mL)", "Pimavanserin (ng/mL)",
                 "Levodopa (ng/mL)", "Antipsychotic (ng/mL)", "Brain ambroxol"),
               title = "Drug concentration time course",
               subtitle = "The output grid is one day, so the peaks of short half-life drugs may not be displayed") +
      scale_colour_manual(values = PAL)
  })

  output$p_carb <- renderPlot({
    d <- sim(); d <- d[d$scenario == "treated", ]
    trace_plot(d, c("INHACHE", "INHBCHE", "CHDRIVE"),
               c("AChE inhibition", "BuChE inhibition", "Cholinergic transducer output"),
               title = "Carbamylated enzyme — the pharmacodynamic half-life is set by enzyme resynthesis",
               hlines = 0.6) + scale_colour_manual(values = PAL)
  })

  ## ------------------------------------------------------- 3. transducers --
  output$p_transducer <- renderPlot({
    trace_plot(sim(), c("CHDRIVE", "DADRIVE", "HT2ASIG"),
               c("Cholinergic (ACh × M1)", "Dopaminergic (DA × D2 × postsynaptic)",
                 "Serotonergic (5-HT × 5-HT2A)"),
               title = "The output of the three transducers — same product, different trajectories",
               subtitle = "All three systems are 'presynaptic × postsynaptic × (1 − occupancy)'") +
      scale_colour_manual(values = PAL)
  })

  output$p_receptor <- renderPlot({
    trace_plot(sim(), c("M1R", "D2R", "HT2A"),
               c("M1/M4 (responds to tau only → preserved in DLB)",
                 "D2 (limbic α-syn prevents up-regulation)",
                 "5-HT2A (denervation + neocortical fibrils → rises)"),
               title = "Postsynaptic receptor density — three equations whose signs differ",
               subtitle = "ChEI superiority · levodopa failure · neuroleptic sensitivity · pimavanserin efficacy all follow from this figure") +
      scale_colour_manual(values = PAL)
  })

  output$tbl_transducer <- renderTable({
    d <- sim(); d <- d[d$scenario == tail(unique(d$scenario), 1), ]
    idx <- sapply(c(0, 1, 3, 5) * 365.25, function(t) which.min(abs(d$time - t)))
    data.frame(
      `Year` = c(0, 1, 3, 5),
      `ACh × M1` = round(d$CHDRIVE[idx], 3),
      `M1 density`  = round(d$M1R[idx], 3),
      `DA signal`  = round(d$DADRIVE[idx], 3),
      `D2 density`  = round(d$D2R[idx], 3),
      `Postsynaptic reserve` = round(d$RESERVE[idx], 3),
      `5-HT2A signal` = round(d$HT2ASIG[idx], 3),
      `5-HT2A density` = round(d$HT2A[idx], 3),
      check.names = FALSE)
  })

  ## --------------------------------------------------------- 4. cognition --
  output$p_cog <- renderPlot({
    trace_plot(sim(), c("MMSE", "CAF", "ATTM", "BISTAB"),
               c("MMSE", "CAF (cognitive fluctuation, 0–16)", "Attention state ATTM", "Bistability depth BISTAB"),
               title = "Cognition: the mean and the variance measure different things") +
      scale_colour_manual(values = PAL)
  })

  output$p_potential <- renderPlot({
    d <- sim()
    dt <- d[d$scenario == tail(unique(d$scenario), 1), ]
    pick <- c(0.5, 1.5, 3, 5)
    A <- seq(-1.6, 1.6, by = 0.01)
    alpha <- as.numeric(param(mod)$ALPHAB)
    off   <- as.numeric(param(mod)$ATTOFF)
    grid <- do.call(rbind, lapply(pick, function(t) {
      i <- which.min(abs(dt$time - t * 365.25))
      D <- dt$DRIVEA[i] - off
      data.frame(A = A,
                 V = A^4/4 - alpha*A^2/2 - D*A,
                 lab = sprintf("%.1f yr  (drive %.3f)", t, dt$DRIVEA[i]))
    }))
    grid <- grid %>% group_by(lab) %>% mutate(V = V - min(V)) %>% ungroup()
    ggplot(grid, aes(A, V, colour = lab)) +
      geom_line(linewidth = 0.8) +
      labs(title = "The attention-state potential V(A) = A⁴/4 − αA²/2 − D·A",
           subtitle = "While there are two wells a state switch is possible, and that is the clinical 'cognitive fluctuation'",
           x = "Attention state A", y = "Potential (relative)") +
      qsp_theme
  })

  ## ---------------------------------------------------- 5. hallucinations --
  output$p_vh <- renderPlot({
    trace_plot(sim(), c("VHB", "VHDRIVE", "HT2ASIG", "OCC2A"),
               c("Hallucination burden (NPI-hall)", "Hallucination drive (product of the three factors)",
                 "5-HT2A signal", "Pimavanserin 5-HT2A occupancy"),
               title = "Visual hallucinations") + scale_colour_manual(values = PAL)
  })

  output$p_pad <- renderPlot({
    d <- sim()
    long <- d %>%
      mutate(`Bottom-up deficit (1−BU)` = 1 - BOTTOMUP,
             `Top-down deficit (1−TD)` = 1 - TOPDOWN,
             `5-HT2A gain`      = HT2ASIG,
             yrs = time/365.25) %>%
      select(yrs, scenario, `Bottom-up deficit (1−BU)`, `Top-down deficit (1−TD)`, `5-HT2A gain`) %>%
      pivot_longer(-c(yrs, scenario))
    ggplot(long, aes(yrs, value, colour = scenario)) +
      geom_line(linewidth = 0.75) + facet_wrap(~name, scales = "free_y") +
      labs(title = "The three factors of the PAD product",
           subtitle = "All three factors must grow together for hallucinations to appear — which is why single-mechanism explanations and single-mechanism treatments both fail",
           x = "Years since diagnosis", y = "") +
      scale_colour_manual(values = PAL) + qsp_theme
  })

  ## ------------------------------------------------------------- 6. motor --
  output$p_motor <- renderPlot({
    trace_plot(sim(), c("MOT", "DADRIVE", "OCCD2", "NSENS"),
               c("MDS-UPDRS III", "Dopaminergic signal", "Striatal D2 occupancy",
                 "Neuroleptic sensitivity index"),
               title = "Motor symptoms") + scale_colour_manual(values = PAL)
  })

  output$p_nsens <- renderPlot({
    end <- input$years * 365
    d <- do.call(rbind, lapply(0:2, function(ph) {
      m <- param(mod, c(build_param(), list(PHENO = ph)))
      x <- mrgsim_df(m, events = if (input$apd == "quet") apd_oral(50, 730, 180)
                                 else apd_oral(1, 730, 180),
                     end = min(end, 1460), delta = 1)
      x$scenario <- c("DLB", "PD-dementia", "AD")[ph + 1]; x
    }))
    trace_plot(d, c("OCCD2", "RESERVE", "DEFICIT", "NSENS", "MOT", "SURV"),
               c("D2 occupancy (identical in all three)", "Postsynaptic reserve (differs)",
                 "Drug-induced deficit = occupancy × reserve amplification", "Sensitivity index",
                 "MDS-UPDRS III", "Survival probability"),
               title = "The same risperidone exposure, three receptor landscapes, three outcomes",
               subtitle = "Nothing but PHENO was changed, and PHENO sets only the initial distribution of pathology") +
      scale_colour_manual(values = PAL)
  })

  ## --------------------------------------------------------- 7. pathology --
  output$p_path <- renderPlot({
    trace_plot(sim(), c("FIBB", "FIBL", "FIBN", "SEED"),
               c("Brainstem fibrils", "Limbic fibrils", "Neocortical fibrils", "Interstitial seed pool"),
               title = "Regional stage progression",
               subtitle = "Uptake in a downstream region is gated by the load in the upstream region, which yields sequential staging") +
      scale_colour_manual(values = PAL)
  })

  output$p_gba <- renderPlot({
    trace_plot(sim(), c("GCASE", "GLCCER", "OLIGTOT", "GCTGT"),
               c("GCase activity", "Glucosylceramide", "Total oligomers", "GCase trafficking target"),
               title = "The GBA1 positive feedback loop",
               subtitle = "Stop the ambroxol and it comes back = monostable. A question of 'lever length', not of a 'deadline'") +
      scale_colour_manual(values = PAL)
  })

  ## -------------------------------------------------------- 8. biomarkers --
  output$p_bio <- renderPlot({
    trace_plot(sim(), c("MIBG", "DATSBR", "EEGF", "PTAU", "ABETA", "QTC"),
               c("MIBG H/M", "DaTSCAN binding ratio", "EEG dominant frequency (Hz)",
                 "tau pathology (au)", "Amyloid (Centiloid)", "QTc change (ms)"),
               title = "Biomarkers") + scale_colour_manual(values = PAL)
  })

  output$tbl_bio <- renderTable({
    d <- sim(); d <- d[d$scenario == tail(unique(d$scenario), 1), ]
    idx <- sapply(c(0, 2, 4, 6) * 365.25, function(t) which.min(abs(d$time - t)))
    data.frame(`Year` = c(0, 2, 4, 6),
               `MIBG H/M` = round(d$MIBG[idx], 2),
               `DaTSCAN`  = round(d$DATSBR[idx], 2),
               `EEG Hz`   = round(d$EEGF[idx], 2),
               `p-tau`    = round(d$PTAU[idx], 3),
               `Centiloid`= round(d$ABETA[idx], 1),
               `QTc (ms)` = round(d$QTC[idx], 1),
               check.names = FALSE)
  })

  ## -------------------------------------------------- 9. all the scenarios --
  all_sc <- eventReactive(input$go_all, {
    withProgress(message = "Running 23 scenarios...", value = 0, {
      run_all()
    })
  })

  output$tbl_all <- renderTable({ summarise_scn(all_sc()) })

  output$p_all <- renderPlot({
    d <- all_sc()
    ggplot(d, aes(time/365.25, MMSE, colour = scenario)) +
      geom_line(linewidth = 0.6) +
      labs(title = "MMSE trajectories of every scenario",
           x = "Years since diagnosis", y = "MMSE") +
      qsp_theme + theme(legend.position = "right",
                        legend.text = element_text(size = 7)) +
      guides(colour = guide_legend(ncol = 1))
  })

  ## --------------------------------------------------------- 10. survival --
  output$p_surv <- renderPlot({
    trace_plot(sim(), c("SURV", "HAZ", "FALLS", "CUMH"),
               c("Survival probability", "Instantaneous hazard", "Fall index", "Cumulative hazard"),
               title = "Survival", hlines = NULL) + scale_colour_manual(values = PAL)
  })

  output$tbl_surv <- renderTable({
    d <- sim()
    do.call(rbind, lapply(split(d, d$scenario), function(x) {
      i <- which(x$SURV <= 0.5)[1]
      data.frame(`Scenario` = x$scenario[1],
                 `Median survival (years)` = if (is.na(i)) NA_real_ else round(x$time[i]/365.25, 2),
                 `5-year survival` = round(x$SURV[which.min(abs(x$time - 1825))], 3),
                 check.names = FALSE)
    }))
  })

  ## ------------------------------------------------------- 11. model card --
  output$model_card <- renderUI({
    HTML('
<h3>Model card</h3>
<table border="0" cellpadding="6">
<tr><td><b>Compartments</b></td><td>66 ODE</td></tr>
<tr><td><b>Time unit</b></td><td>days</td></tr>
<tr><td><b>Integrator</b></td><td>LSODA (mrgsolve), rtol 1e-8 / atol 1e-8</td></tr>
<tr><td><b>Drugs</b></td><td>rivastigmine (oral·patch), donepezil, galantamine, memantine,
pimavanserin (+AC-279), levodopa/carbidopa, a generic antipsychotic compartment, zonisamide,
ambroxol, anti-α-synuclein antibody, droxidopa, melatonin/clonazepam, modafinil</td></tr>
<tr><td><b>Phenotype</b></td><td>DLB · PD-dementia · AD — only the <i>initial distribution of pathology</i> differs; the equations are identical</td></tr>
<tr><td><b>Inter-individual variation</b></td><td>residual nbM/PPN·SNc/LC at diagnosis, limbic·neocortical load, amyloid</td></tr>
</table>

<h4>What the model <b>claims</b></h4>
<ol>
<li>The three transmitter systems share one transducer, and only the <b>sign</b> of postsynaptic plasticity differs.
ChEI superiority·levodopa failure·neuroleptic sensitivity·pimavanserin efficacy are not separate rules but
<b>consequences</b> of this structure.</li>
<li>Cognitive fluctuation is not severity but <b>variance</b> — the switching frequency of a bistable attention state.
A ChEI therefore improves the fluctuation measure proportionally more than MMSE (model calculation: about 3.5-fold).</li>
<li>Hallucinations are the <b>product</b> of three deficits. So they are neither explained nor abolished by any single mechanism.</li>
<li>Neuroleptic sensitivity is a <b>drug event</b>, not a severity milestone —
off drug the D2 occupancy is zero, so the deficit is zero as well.</li>
</ol>

<h4>What the model does <b>not</b> claim</h4>
<ul>
<li><b>The GBA1 loop is not bistable.</b> The draft tried to put a saddle-node bifurcation in, but
with parameters in the range the literature supports the loop gain did not reach one. Stop the ambroxol and
GCase comes back (check it directly in tab 7). So what this model supports is not "miss the deadline and it is over" but
"the earlier you start, the longer the lever".</li>
<li>The rate constants for neuronal loss were not measured individually; they were fitted simultaneously to
natural-history reference points and are correlated with each other. No biological meaning should be read into the individual values.</li>
<li>None of the predictions has been validated prospectively.</li>
</ul>

<p style="color:#a00;"><b>Not for clinical use.</b> This is a semi-quantitative model for teaching and research purposes.</p>
')
  })
}

shinyApp(ui, server)
