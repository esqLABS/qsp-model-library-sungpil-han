# =============================================================================
#  X-ALD QSP — Shiny Dashboard
#  X-linked Adrenoleukodystrophy · interactive QSP explorer
#
#  Run:
#     setwd("x-linked-adrenoleukodystrophy")
#     shiny::runApp("xald_shiny_app.R")
#
#  Required packages: shiny, mrgsolve, ggplot2, dplyr, tidyr, DT
#
#  Design intent of this app
#  -----------------------------------------------------------------------------
#  Tab 2 ("the two variables") is this app's reason for existing. Whatever combination the user runs,
#  they will see the left side of the screen (measured plasma C26:0) and the right side (unmeasured
#  cerebral inflammation switch) **move independently of each other**. Lorenzo's oil moves only the left,
#  gene therapy moves only the right. That mismatch is the clinical reality of this disease.
# =============================================================================

library(shiny)
library(mrgsolve)
library(ggplot2)
library(dplyr)
library(tidyr)
library(DT)

source("xald_mrgsolve_model_en.R")

MOD <- xald_mod()
assign(".XALD_MOD", MOD, envir = globalenv())

THEME <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 13),
        legend.position = "bottom")

# -----------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("X-linked Adrenoleukodystrophy (X-ALD) — QSP Exploration Tool"),
  tags$p(tags$em(paste(
    "One lesion (ABCD1 loss), two independent variables —",
    "measurable VLCFA accumulation (reachable from the periphery), and",
    "an unmeasurable bistable cerebral inflammation switch (reachable only from inside the CNS)."))),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("① Patient (Genotype)"),
      sliderInput("mutres", "Residual ALDP function (MUTRES)", 0, 1, 0.02, step = 0.01),
      sliderInput("fmos", "Fraction of ALDP-deficient cells (female carrier FMOS)", 0, 1, 0, step = 0.05),
      sliderInput("susc", "Cerebral susceptibility (SUSCTV)", 0.05, 1.6, 0.46, step = 0.01),
      helpText(HTML(paste0(
        "Threshold <b>SUSC* &asymp; 0.41</b> (an output of the model, found by bisection).<br>",
        "Below it = pure AMN, above it = cerebral. The closer to the threshold, the <b>later</b> the onset",
        " (critical slowing down)."))),
      hr(),
      h4("② Treatment"),
      checkboxInput("lo",   "Lorenzo's oil (oral, peripheral competitive ELOVL1 inhibition)", FALSE),
      sliderInput("lo_start", "Lorenzo's oil start (years)", 0.5, 20, 1.5, step = 0.5),
      checkboxInput("diet", "VLCFA-restricted diet", FALSE),
      checkboxInput("gt",   "Gene therapy / HSCT", FALSE),
      sliderInput("gt_loes", "Transplant timing (Loes score)", 1, 20, 2, step = 1),
      sliderInput("vcn", "Vector copy number VCN (0 = allogeneic transplant)", 0, 2, 0.8, step = 0.1),
      sliderInput("bu", "Busulfan exposure multiplier", 0.5, 1.5, 1.0, step = 0.02),
      checkboxInput("hc",   "Hydrocortisone supplementation", FALSE),
      checkboxInput("leri", "Leriglitazone (CNS-penetrant PPARγ)", FALSE),
      checkboxInput("sob",  "Sobetirome class (CNS-penetrant metabolic target)", FALSE),
      checkboxInput("ster", "High-dose steroids", FALSE),
      checkboxInput("anti", "Triple antioxidant + DMF", FALSE),
      hr(),
      sliderInput("years", "Observation period (years)", 5, 60, 20, step = 5),
      actionButton("run", "Run simulation", class = "btn-primary btn-block")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("1. Patient Summary",
                 br(), verbatimTextOutput("summary"),
                 h5("Phenotype determination"), verbatimTextOutput("phenotype")),
        tabPanel("2. Two Variables (Core)",
                 br(), plotOutput("twovar", height = "460px"),
                 helpText(paste("The left is what is measured, the right is what is not measured.",
                                "Turn on treatment and check whether the two panels move together."))),
        tabPanel("3. VLCFA Metabolism · Drug Exposure",
                 br(), plotOutput("metab", height = "440px")),
        tabPanel("4. Cerebral Switch · Gain",
                 br(), plotOutput("switchplot", height = "440px"),
                 helpText("If LOOPGAIN > 1, inflammation sustains itself without any ignition input.")),
        tabPanel("5. MRI · Neurological Endpoints",
                 br(), plotOutput("neuro", height = "440px")),
        tabPanel("6. Adrenal Axis",
                 br(), plotOutput("adrenal", height = "440px")),
        tabPanel("7. Spinal Form (AMN)",
                 br(), plotOutput("cord", height = "440px")),
        tabPanel("8. Transplant · Gene Therapy",
                 br(), plotOutput("transplant", height = "440px"),
                 helpText(paste("Microglial replacement is slow (KMGREP).",
                                "That slowness is the cause of the 12-18 month lag window after transplant."))),
        tabPanel("9. Scenario Comparison",
                 br(),
                 selectInput("scens", "Scenarios to compare", multiple = TRUE,
                             choices = names(xald_scenarios()),
                             selected = c("s10_elicel_early",
                                          "s11_elicel_early_control",
                                          "s07_lo_ccald")),
                 selectInput("scvar", "Variable to display",
                             choices = c("LOES","NFS","C26","C26LYSO","EDSS",
                                         "MGLCORR","CORTISOL","CHIMER"),
                             selected = "LOES"),
                 actionButton("runsc", "Run scenario"),
                 plotOutput("scenplot", height = "420px")),
        tabPanel("10. Population · Phenotype Distribution",
                 br(),
                 numericInput("npop", "Number of individuals", 150, min = 20, max = 600, step = 10),
                 actionButton("runpop", "Run population simulation"),
                 plotOutput("popplot", height = "380px"),
                 verbatimTextOutput("poptext")),
        tabPanel("11. Validation Anchors",
                 br(), actionButton("runval", "Run validation suite (A1-A23)"),
                 DTOutput("valtab"))
      )
    )
  )
)

# -----------------------------------------------------------------------------
server <- function(input, output, session) {

  build <- reactive({
    p <- list(MUTRES = input$mutres, FMOS = input$fmos, SUSCTV = input$susc,
              VCNIN = input$vcn,
              DIETSC = if (input$diet) 0.35 else 1.0,
              STEROIDX = if (input$ster) 1.0 else 0.0,
              NACD = if (input$anti) 1.0 else 0.0,
              DMFD = if (input$anti) 1.0 else 0.0)
    end <- input$years * YR
    ev_list <- list(); txday <- NULL
    if (input$lo)   ev_list <- c(ev_list, list(.lo_regimen(input$lo_start*YR, end)))
    if (input$hc)   ev_list <- c(ev_list, list(.hc_regimen(7.5*YR, end)))
    if (input$leri) ev_list <- c(ev_list, list(.leri_regimen(2*YR, end)))
    if (input$sob)  ev_list <- c(ev_list, list(.sob_regimen(2*YR, end)))
    if (input$gt) {
      # The transplant timing is the day the user-chosen Loes value is reached — i.e., an output of the matched natural history
      p0  <- list(MUTRES = input$mutres, FMOS = input$fmos, SUSCTV = input$susc)
      nat <- MOD %>% param(p0) %>% init(xald_birth(MOD, p0)) %>% zero_re() %>%
        mrgsim(end = end, delta = 10, maxsteps = 1e6) %>% as.data.frame()
      i <- which(nat$LOES >= input$gt_loes)
      if (length(i)) {
        ev_list <- c(ev_list, list(.hsct_regimen(nat$time[i[1]],
                                                 vcn = input$vcn,
                                                 bu_scale = input$bu)))
        txday <- nat$time[i[1]]
      }
    }
    e <- if (length(ev_list)) do.call(c, ev_list) else NULL
    list(param = p, ev = e, end = end, txday = txday)
  })

  simdata <- eventReactive(input$run, {
    b <- build()
    init0 <- xald_birth(MOD, b$param)
    m <- MOD %>% param(b$param) %>% init(init0) %>% zero_re()
    d <- if (is.null(b$ev)) {
      m %>% mrgsim(end = b$end, delta = 10, maxsteps = 1e6)
    } else {
      m %>% mrgsim(events = b$ev, end = b$end, delta = 10, maxsteps = 1e6)
    }
    d <- as.data.frame(d)
    d$age <- d$time / YR
    attr(d, "txday") <- b$txday
    d
  }, ignoreNULL = FALSE)

  # ---- 1. summary ----
  output$summary <- renderPrint({
    d <- simdata()
    last <- d[nrow(d), ]
    cat(sprintf("Observation period    : %.0f years\n", max(d$age)))
    cat(sprintf("Plasma C26:0          : %.3f umol/L  (normal 0.35, elevated %.2fx)\n",
                last$C26, last$C26FOLD))
    cat(sprintf("C26:0-lysoPC         : %.3f umol/L\n", last$C26LYSO))
    cat(sprintf("Brain VLCFA permeation gate : %.3f  (saturated in patients -> cannot discriminate phenotype)\n",
                last$VLCGATE))
    cat(sprintf("Peak Loes / NFS       : %.1f / %.1f\n", max(d$LOES), max(d$NFS)))
    cat(sprintf("Peak EDSS             : %.2f\n", max(d$EDSS)))
    cat(sprintf("Cortisol (final)      : %.0f nmol/L,  ACTH %.0f pg/mL\n",
                last$CORTISOL, last$ACTHOUT))
    cat(sprintf("Adrenal reserve       : %.3f  (adrenal insufficiency threshold 0.36)\n", last$ADRRESV))
    cat(sprintf("Monocyte chimerism    : %.3f,  corrected microglia %.3f\n",
                last$CHIMER, last$MGLCORR))
    if (max(d$PMDS) > 0)
      cat(sprintf("MDS/AML cumulative probability : %.3f\n", max(d$PMDS)))
    if (!is.null(attr(d, "txday")))
      cat(sprintf("Transplant timing     : %.2f years\n", attr(d, "txday")/YR))
  })

  output$phenotype <- renderPrint({
    d <- simdata()
    fired <- max(d$MGP) > 0.15
    cat(if (input$mutres >= 0.9) "Non-carrier / normal\n"
        else if (fired) "Cerebral form (cerebral ALD) — switch fired\n"
        else if (input$fmos > 0) "Female carrier phenotype (spinal-form predominant)\n"
        else "Pure AMN (spinal form) — switch not fired\n")
    if (fired) {
      i <- which(d$LOES >= 1)
      if (length(i)) cat(sprintf("Cerebral onset (Loes >= 1) : %.2f years\n", d$age[i[1]]))
      j <- which(d$LOES >= 15)
      if (length(i) && length(j))
        cat(sprintf("Loes 1 -> 15 elapsed  : %.2f years\n", (d$age[j[1]] - d$age[i[1]])))
    }
    k <- which(d$ADRINSUF >= 1)
    if (length(k)) cat(sprintf("Adrenal insufficiency onset : %.2f years\n", d$age[k[1]]))
    cat(sprintf("Major functional disability (MFD) : %s\n", if (max(d$MFD) > 0) "Present" else "Absent"))
  })

  # ---- 2. the two variables ----
  output$twovar <- renderPlot({
    d <- simdata()
    a <- ggplot(d, aes(age)) +
      geom_line(aes(y = C26, colour = "Plasma C26:0"), linewidth = 1.1) +
      geom_line(aes(y = C26LYSO*2, colour = "C26:0-lysoPC (x2)"), linewidth = 1) +
      geom_hline(yintercept = 0.35, linetype = 2, colour = "grey40") +
      annotate("text", x = max(d$age)*0.6, y = 0.40, label = "Normal C26:0",
               size = 3.2, colour = "grey30") +
      scale_colour_manual(values = c("Plasma C26:0" = "#c0392b",
                                     "C26:0-lysoPC (x2)" = "#e67e22")) +
      labs(title = "Variable (1) The measured one — reachable from the periphery",
           x = "Age (years)", y = "umol/L", colour = NULL) + THEME
    b <- ggplot(d, aes(age)) +
      geom_line(aes(y = MGP, colour = "Active microglia (switch)"), linewidth = 1.2) +
      geom_line(aes(y = DEMYELIN, colour = "Demyelination fraction"), linewidth = 1) +
      geom_line(aes(y = GADENH, colour = "Gadolinium enhancement"), linewidth = 0.9,
                linetype = 3) +
      scale_colour_manual(values = c("Active microglia (switch)" = "#8e44ad",
                                     "Demyelination fraction" = "#2471a3",
                                     "Gadolinium enhancement" = "#16a085")) +
      labs(title = "Variable (2) The unmeasured one — reachable only from inside the CNS",
           x = "Age (years)", y = "Fraction", colour = NULL) + THEME
    gridExtra_arrange(a, b)
  })

  # ---- 3. metabolism ----
  output$metab <- renderPlot({
    d <- simdata()
    long <- d %>%
      select(age, `Plasma C26:0` = C26, `Brain white matter VLCFA` = CBR,
             `Spinal cord VLCFA` = CSC, `Adrenal cortex VLCFA` = CADR,
             `Erucic acid (umol/L /100)` = ERUCIC) %>%
      mutate(`Erucic acid (umol/L /100)` = `Erucic acid (umol/L /100)`/100) %>%
      pivot_longer(-age)
    ggplot(long, aes(age, value, colour = name)) +
      geom_line(linewidth = 1.05) +
      labs(title = "VLCFA and drug exposure by tissue (tissue pools normalised to normal = 1.0)",
           subtitle = paste("About 88% of brain/spinal cord VLCFA is local synthesis —",
                            "erucic acid cannot cross the BBB, so it cannot reach here"),
           x = "Age (years)", y = "Relative value", colour = NULL) + THEME
  })

  # ---- 4. switch ----
  output$switchplot <- renderPlot({
    d <- simdata()
    a <- ggplot(d, aes(age)) +
      geom_line(aes(y = MGP), colour = "#8e44ad", linewidth = 1.2) +
      geom_hline(yintercept = 0.15, linetype = 2, colour = "grey40") +
      labs(title = "Active microglia fraction (bistable state variable)",
           x = "Age (years)", y = "MGP") + THEME
    b <- ggplot(d, aes(age)) +
      geom_line(aes(y = LOOPGAIN), colour = "#c0392b", linewidth = 1.2) +
      geom_hline(yintercept = 1, linetype = 2) +
      labs(title = "Self-amplifying loop gain (1 = the point of no return)",
           x = "Age (years)", y = "LOOPGAIN") + THEME
    gridExtra_arrange(a, b)
  })

  # ---- 5. neuro endpoints ----
  output$neuro <- renderPlot({
    d <- simdata()
    long <- d %>% select(age, `Loes (0-34)` = LOES, `NFS (0-25)` = NFS) %>%
      pivot_longer(-age)
    ggplot(long, aes(age, value, colour = name)) +
      geom_line(linewidth = 1.15) +
      geom_hline(yintercept = 9, linetype = 3, colour = "#e67e22") +
      annotate("text", x = min(d$age) + 1, y = 10.2, hjust = 0,
               label = "Loes 9 = upper bound of the transplant treatment window", size = 3.2,
               colour = "#e67e22") +
      labs(title = "Imaging · neurological function scores (accumulate in one direction only)",
           x = "Age (years)", y = "Score", colour = NULL) + THEME
  })

  # ---- 6. adrenal ----
  output$adrenal <- renderPlot({
    d <- simdata()
    a <- ggplot(d, aes(age)) +
      geom_line(aes(y = CORTISOL, colour = "Cortisol (nmol/L)"), linewidth = 1.1) +
      geom_line(aes(y = ACTHOUT, colour = "ACTH (pg/mL)"), linewidth = 1.1) +
      scale_colour_manual(values = c("Cortisol (nmol/L)" = "#2471a3",
                                     "ACTH (pg/mL)" = "#c0392b")) +
      labs(title = "Adrenal cortex axis", x = "Age (years)", y = NULL, colour = NULL) + THEME
    b <- ggplot(d, aes(age)) +
      geom_line(aes(y = ADRRESV, colour = "Steroidogenic reserve"), linewidth = 1.1) +
      geom_line(aes(y = POTASSIUM/10, colour = "Serum K+ / 10"), linewidth = 1.1) +
      geom_hline(yintercept = 0.36, linetype = 2, colour = "grey40") +
      scale_colour_manual(values = c("Steroidogenic reserve" = "#8e44ad",
                                     "Serum K+ / 10" = "#e67e22")) +
      labs(title = "Reserve and electrolytes (dashed line = adrenal insufficiency threshold)",
           x = "Age (years)", y = NULL, colour = NULL) + THEME
    gridExtra_arrange(a, b)
  })

  # ---- 7. cord ----
  output$cord <- renderPlot({
    d <- simdata()
    long <- d %>% select(age, `EDSS` = EDSS, `Spinal cord injury fraction` = CORDLOSS,
                         `6-minute walk (m/100)` = WALK6MIN) %>%
      mutate(`6-minute walk (m/100)` = `6-minute walk (m/100)`/100) %>%
      pivot_longer(-age)
    ggplot(long, aes(age, value, colour = name)) +
      geom_line(linewidth = 1.1) +
      labs(title = "Spinal axonopathy (AMN) — progresses independently of the switch",
           subtitle = "The longest axons go first (LENCST > LENDC), CNS axon regeneration is 0",
           x = "Age (years)", y = NULL, colour = NULL) + THEME
  })

  # ---- 8. transplant ----
  output$transplant <- renderPlot({
    d <- simdata()
    long <- d %>% select(age, `Monocyte chimerism` = CHIMER,
                         `Corrected microglia` = MGLCORR,
                         `Total microglia` = MGLTOT,
                         `Neutrophils` = NEUTRO,
                         `MDS/AML cumulative probability` = PMDS) %>%
      pivot_longer(-age)
    p <- ggplot(long, aes(age, value, colour = name)) +
      geom_line(linewidth = 1.05) +
      labs(title = "Post-transplant reconstitution — fast monocytes, slow microglia",
           subtitle = "The gap between the two curves is exactly the treatment-lag window",
           x = "Age (years)", y = "Fraction / probability", colour = NULL) + THEME
    if (!is.null(attr(d, "txday")))
      p <- p + geom_vline(xintercept = attr(d, "txday")/YR, linetype = 2)
    p
  })

  # ---- 9. scenarios ----
  scendata <- eventReactive(input$runsc, {
    req(input$scens)
    bind_rows(lapply(input$scens, function(nm) {
      d <- xald_run(nm, MOD, delta = 20)
      data.frame(age = d$time/YR, value = d[[input$scvar]], scenario = nm)
    }))
  })
  output$scenplot <- renderPlot({
    ggplot(scendata(), aes(age, value, colour = scenario)) +
      geom_line(linewidth = 1.1) +
      labs(title = paste0("Scenario comparison — ", input$scvar),
           x = "Age (years)", y = input$scvar, colour = NULL) + THEME
  })

  # ---- 10. population ----
  popdata <- eventReactive(input$runpop, {
    xald_population(MOD, n = input$npop)
  })
  output$popplot <- renderPlot({
    p <- popdata()
    ggplot(p, aes(SUSC, maxLOES, colour = factor(cerebral))) +
      geom_point(alpha = 0.75, size = 2) +
      geom_vline(xintercept = 0.41, linetype = 2) +
      scale_colour_manual(values = c("0" = "#2471a3", "1" = "#c0392b"),
                          labels = c("Spinal form (AMN)", "Cerebral form"), name = NULL) +
      labs(title = "Inter-individual variation in a single susceptibility parameter creates two diseases",
           subtitle = "Dashed line = threshold found by bisection. Plasma C26:0 is identical at every point.",
           x = "Susceptibility SUSC", y = "Peak Loes score") + THEME
  })
  output$poptext <- renderPrint({
    p <- popdata()
    cat(sprintf("Cerebral conversion rate : %.1f%%   (literature 35-40%%)\n", 100*mean(p$cerebral)))
    cat(sprintf("Adrenal insufficiency rate : %.1f%%   (literature about 80%% — the model over-predicts)\n",
                100*mean(p$adrinsuf)))
    o <- p$onset[!is.na(p$onset)]
    if (length(o))
      cat(sprintf("Cerebral onset age : median %.1f years (range %.1f-%.1f)  (literature 4-8 years)\n",
                  median(o), min(o), max(o)))
  })

  # ---- 11. validation ----
  valdata <- eventReactive(input$runval, { xald_validate(MOD, verbose = FALSE) })
  output$valtab <- renderDT({
    datatable(valdata(), rownames = FALSE,
              options = list(pageLength = 25, dom = "t")) %>%
      formatStyle("pct", color = styleInterval(c(-20, 20),
                                              c("#c0392b", "#1a7f37", "#c0392b")))
  })
}

# Two-panel helper without taking a hard dependency on gridExtra
gridExtra_arrange <- function(p1, p2) {
  if (requireNamespace("gridExtra", quietly = TRUE)) {
    gridExtra::grid.arrange(p1, p2, ncol = 1)
  } else {
    p1   # if gridExtra is unavailable, show only the first panel
  }
}

shinyApp(ui, server)
