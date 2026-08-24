###############################################################################
##  Kidney Transplant Rejection & Immunosuppression — QSP interactive dashboard
##  ---------------------------------------------------------------------------
##  ktx_shiny_app.R      companion to ktx_mrgsolve_model.R / ktx_qsp_model.dot
##
##  The point of this app is to let you sit where the transplant clinician sits:
##  you cannot see the alloimmune response, you can only see a trough level, a
##  creatinine, a DSA titre and a viral load — and you have exactly one dial,
##  net immunosuppression.  Turn it down and the graft rejects; turn it up and
##  the virus wins.  Tab 8 makes that trade-off explicit by sweeping the dial.
##
##  Ten tabs:
##    1  Donor & recipient      the fixed risk set, and the resulting net IS
##    2  Drug exposure & TDM    troughs vs protocol targets, MPA AUC, biologics
##    3  Cellular alloimmunity  naive/blast/effector/memory/Treg, graft infiltrate
##    4  DSA, complement, NK    the humoral arm and why rituximab under-performs
##    5  Histology (Banff)      t, i, g, ptc, cg, ci, ah, C4d over 5 years
##    6  Graft function         eGFR, creatinine, proteinuria, dd-cfDNA
##    7  Infection & safety     BK, CMV, leucocytes, glucose, cumulative IS
##    8  The U-curve            sweep net immunosuppression, see both failure modes
##    9  Scenario comparison    stack any number of saved runs side by side
##   10  Parameters & notes     the full annotated parameter set
##
##  Run:  Rscript -e 'shiny::runApp("ktx_shiny_app.R")'
###############################################################################

suppressPackageStartupMessages({
  library(shiny)
  library(mrgsolve)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

## --------------------------------------------------------------------------
## Load the model definition straight out of ktx_mrgsolve_model.R so that the
## app and the batch script can never drift apart.
## --------------------------------------------------------------------------
MODEL_FILE <- "ktx_mrgsolve_model.R"
if (!file.exists(MODEL_FILE)) stop("ktx_mrgsolve_model.R must sit next to this app")
src <- readLines(MODEL_FILE, warn = FALSE)
i1  <- grep("^code <- '", src)[1]
i2  <- grep("^mod <- mcode", src)[1]
eval(parse(text = paste(src[i1:(i2 - 1)], collapse = "\n")))
mod <- mcode("ktx_shiny", code, soloc = tempdir())

THEME <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")

PAL <- c("#2d5b8c", "#9c2d57", "#1f6b52", "#8a5a20", "#553f85",
         "#7d3269", "#25566f", "#6a5f1e", "#9b3535", "#46632c")

yrs <- function(d) d / 365.25

## a horizontal reference band helper
band <- function(lo, hi, fill = "#4caf50", alpha = 0.10) {
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = lo, ymax = hi,
           fill = fill, alpha = alpha)
}

###############################################################################
## UI
###############################################################################
ui <- fluidPage(
  titlePanel("Kidney Transplant Rejection QSP Model"),
  tags$p(style = "color:#666;margin-top:-8px",
         "59-ODE mechanistic model · alloimmunity vs infection, balanced on one dial: net immunosuppression"),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      tags$h4("Donor"),
      sliderInput("KDPI", "KDPI (%)", 0, 100, 40, 5),
      sliderInput("CIT", "Cold ischaemia (h)", 1, 30, 14, 1),
      checkboxInput("DCD", "DCD donor (donation after circulatory death)", FALSE),

      tags$h4("Immunological risk"),
      sliderInput("HLAMM", "HLA A/B/DR mismatches (0-6)", 0, 6, 3, 1),
      sliderInput("EPLET", "HLA-DR/DQ eplet mismatch load", 0, 35, 12, 1),
      numericInput("PREDSA", "Preformed DSA (MFI)", 0, min = 0, max = 20000, step = 500),
      sliderInput("FCD28N", "CD28-null memory T-cell fraction (belatacept resistance)",
                  0, 1, 0.45, 0.05),
      sliderInput("TMEM0", "Heterologous cross-reactive memory T-cell pool (heterologous memory)",
                  0, 1.5, 0.70, 0.05),

      tags$h4("Maintenance"),
      selectInput("regimen", "Backbone regimen",
                  c("TAC + MPA + steroid"      = "tac",
                    "CsA + MPA + steroid"      = "csa",
                    "Belatacept + MPA + steroid (CNI-free)" = "bela",
                    "TAC (reduced) + everolimus"   = "evr")),
      selectInput("induction", "Induction", c("Basiliximab" = "bas", "rATG" = "atg", "None" = "none")),
      sliderInput("TACTGT2", "Tacrolimus target trough, months 3-12 (ng/mL)", 3, 12, 7, 0.5),
      checkboxInput("TDM", "Perform TDM (concentration-based dose adjustment)", TRUE),
      checkboxInput("CYP3A5", "CYP3A5 expresser (*1 carrier, 1.9-fold clearance)", FALSE),
      sliderInput("MMFDOSE", "MMF dose (mg/day)", 0, 3000, 2000, 250),
      numericInput("PREDSTOP", "Steroid discontinuation day (large value if none)", 100000, min = 5, step = 1),

      tags$h4("Adherence"),
      sliderInput("ADHFINAL", "Fraction of dose taken during non-adherence period", 0, 1, 1, 0.05),
      numericInput("ADHSTART", "Non-adherence start day", 100000, min = 0, step = 30),
      numericInput("ADHEND", "Adherence recovery day", 100000, min = 0, step = 30),
      sliderInput("IPVAMP", "Tacrolimus exposure variability, IPV (amplitude)", 0, 0.8, 0, 0.05),

      tags$h4("Infection"),
      checkboxInput("BKSUSC", "BK virus reactivation predisposition", FALSE),
      checkboxInput("BKSCREEN", "BK screening + pre-emptive immunosuppression reduction (closed loop)", TRUE),
      checkboxInput("CMVDR", "CMV D+/R- mismatch", FALSE),
      numericInput("VGCSTOP", "Valganciclovir prophylaxis end day", 90, min = 0, step = 30),

      tags$h4("Rejection therapy"),
      checkboxInput("MPON", "Methylprednisolone pulse therapy", FALSE),
      numericInput("MPSTART", "  Pulse therapy start day", 95, min = 0, step = 5),
      checkboxInput("ATG2ON", "ATG rescue therapy (steroid-resistant)", FALSE),
      numericInput("ATG2START", "  ATG rescue start day", 105, min = 0, step = 5),
      checkboxInput("PLEXON", "Plasma exchange (PLEX)", FALSE),
      checkboxInput("IVIGON", "IVIG 2 g/kg", FALSE),
      checkboxInput("RTXON", "Rituximab (anti-CD20)", FALSE),
      checkboxInput("TCZON", "Tocilizumab/Clazakimab (anti-IL-6)", FALSE),
      checkboxInput("FZBON", "Felzartamab (anti-CD38)", FALSE),
      numericInput("ABMRDAY", "ABMR treatment start day", 730, min = 0, step = 30),

      tags$hr(),
      actionButton("save", "Save this run to the comparison tab", class = "btn-primary"),
      actionButton("clear", "Clear comparison list")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("1 Donor · Recipient", plotOutput("p_profile", height = "300px"),
                 tableOutput("t_profile"), verbatimTextOutput("txt_profile")),
        tabPanel("2 Drug Exposure · TDM", plotOutput("p_pk", height = "620px")),
        tabPanel("3 Cellular Alloimmunity", plotOutput("p_cell", height = "620px")),
        tabPanel("4 DSA · Complement · NK", plotOutput("p_dsa", height = "620px")),
        tabPanel("5 Histology (Banff)", plotOutput("p_banff", height = "620px")),
        tabPanel("6 Graft Function", plotOutput("p_func", height = "620px")),
        tabPanel("7 Infection · Safety", plotOutput("p_inf", height = "620px")),
        tabPanel("8 The U-curve (immunosuppression balance)",
                 tags$p("Sweeps immunosuppression intensity and computes rejection burden and infection burden simultaneously. Takes a few seconds."),
                 actionButton("sweep", "Compute U-curve", class = "btn-warning"),
                 plotOutput("p_ucurve", height = "540px")),
        tabPanel("9 Scenario Comparison", plotOutput("p_cmp", height = "620px"),
                 tableOutput("t_cmp")),
        tabPanel("10 Parameters · Notes",
                 tags$p("The model's full parameter set ($PARAM block of ktx_mrgsolve_model.R)."),
                 tableOutput("t_par"))
      )
    )
  )
)

###############################################################################
## SERVER
###############################################################################
server <- function(input, output, session) {

  saved <- reactiveVal(list())

  ## a tiny stacker so the app has no gridExtra dependency
  gridExtra_like <- function(plots, heights) {
    n <- length(plots)
    grid::grid.newpage()
    tot <- sum(heights)
    y0 <- 1
    for (i in seq_len(n)) {
      h <- heights[i] / tot
      print(plots[[i]], vp = grid::viewport(x = 0.5, y = y0 - h / 2,
                                            width = 1, height = h))
      y0 <- y0 - h
    }
  }


  pars <- reactive({
    p <- list(
      KDPI = input$KDPI, CIT = input$CIT, DCD = as.numeric(input$DCD),
      HLAMM = input$HLAMM, EPLET = input$EPLET, PREDSA = input$PREDSA,
      FCD28N = input$FCD28N, TMEM0 = input$TMEM0,
      TACTGT2 = input$TACTGT2, TACTGT3 = max(4, input$TACTGT2 - 1),
      TDM = as.numeric(input$TDM), CYP3A5 = as.numeric(input$CYP3A5),
      MMFDOSE = input$MMFDOSE, MMFON = as.numeric(input$MMFDOSE > 0),
      PREDSTOP = input$PREDSTOP,
      ADHFINAL = input$ADHFINAL, ADHSTART = input$ADHSTART, ADHEND = input$ADHEND,
      IPVAMP = input$IPVAMP, PIPV = 21,
      BKSUSC = as.numeric(input$BKSUSC), BKSCREEN = as.numeric(input$BKSCREEN),
      CMVDR = as.numeric(input$CMVDR), CMVSUSC = as.numeric(input$CMVDR),
      VGCSTOP = input$VGCSTOP,
      MPON = as.numeric(input$MPON), MPSTART = input$MPSTART,
      ATG2ON = as.numeric(input$ATG2ON), ATG2START = input$ATG2START,
      PLEXON = as.numeric(input$PLEXON), PLEXSTART = input$ABMRDAY,
      IVIGON = as.numeric(input$IVIGON), IVIGSTART = input$ABMRDAY + 10,
      RTXON  = as.numeric(input$RTXON),  RTXDAY   = input$ABMRDAY + 12,
      TCZON  = as.numeric(input$TCZON),  TCZSTART = input$ABMRDAY, TCZN = 18,
      FZBON  = as.numeric(input$FZBON),  FZBSTART = input$ABMRDAY
    )
    # maintenance backbone
    p$TACON <- 0; p$CSAON <- 0; p$BELAON <- 0; p$EVRON <- 0
    if (input$regimen == "tac")  p$TACON <- 1
    if (input$regimen == "csa")  p$CSAON <- 1
    if (input$regimen == "bela") { p$BELAON <- 1; p$TDM <- 0 }
    if (input$regimen == "evr") {
      p$TACON <- 1; p$EVRON <- 1; p$MMFON <- 0
      p$TACTGT1 <- 5; p$TACTGT2 <- 4; p$TACTGT3 <- 4
    }
    # induction
    p$BASON <- as.numeric(input$induction == "bas")
    p$ATGON <- as.numeric(input$induction == "atg")
    p
  })

  sim <- reactive({
    mod %>% param(pars()) %>%
      mrgsim(end = 1825, delta = 2, hmax = 3, atol = 1e-8, rtol = 1e-5) %>%
      as_tibble()
  })

  long <- function(d, vars, labs) {
    d %>% select(time, all_of(vars)) %>%
      rename_with(~labs, all_of(vars)) %>%
      pivot_longer(-time) %>%
      mutate(name = factor(name, levels = labs))
  }

  gg <- function(d, ncol = 2, ylab = NULL) {
    ggplot(d, aes(yrs(time), value)) +
      geom_line(linewidth = 0.7, colour = PAL[1]) +
      facet_wrap(~name, scales = "free_y", ncol = ncol) +
      labs(x = "Time since transplant (years)", y = ylab) + THEME
  }

  ## ---- 1  profile ---------------------------------------------------------
  output$p_profile <- renderPlot({
    d <- sim()
    long(d, c("NISr", "CNINHr", "MPAINHr", "STINHr", "COSBLKr", "ATGDr"),
         c("Net immunosuppression (NIS)", "Calcineurin inhibition",
           "IMPDH inhibition", "Glucocorticoid effect",
           "Costimulation blockade", "Lymphodepletion")) %>%
      gg(ncol = 3, ylab = "effect (0-1)")
  })

  output$t_profile <- renderTable({
    d <- sim()
    idx <- function(day) which.min(abs(d$time - day))
    data.frame(
      Time_point = c("3 months", "12 months", "3 years", "5 years"),
      eGFR = round(d$EGFRr[c(idx(90), idx(365), idx(1095), idx(1825))], 1),
      `Serum Cr (mg/dL)` = round(d$SCRDL[c(idx(90), idx(365), idx(1095), idx(1825))], 2),
      `TAC C0` = round(d$TACC0r[c(idx(90), idx(365), idx(1095), idx(1825))], 1),
      `DSA (MFI)` = round(d$DSA[c(idx(90), idx(365), idx(1095), idx(1825))], 0),
      `Banff t` = round(d$BANFF_T[c(idx(90), idx(365), idx(1095), idx(1825))], 2),
      `g+ptc` = round(d$MVISUM[c(idx(90), idx(365), idx(1095), idx(1825))], 2),
      `Graft survival probability` = round(d$GS[c(idx(90), idx(365), idx(1095), idx(1825))], 3),
      check.names = FALSE
    )
  })

  output$txt_profile <- renderPrint({
    d <- sim()
    cat("Expected delayed graft function (DGF):", ifelse(d$DGF[1] > 0.5, "yes", "no"), "\n")
    cat("Lowest eGFR:", round(min(d$EGFRr), 1), "mL/min/1.73 m2  (day",
        round(d$time[which.min(d$EGFRr)]), ")\n")
    cat("Peak Banff t:", round(max(d$BANFF_T), 2),
        " / peak g+ptc:", round(max(d$MVISUM), 2), "\n")
    cat("Peak DSA:", round(max(d$DSA)), "MFI\n")
    cat("Peak BK viral load:", round(max(d$BKLOG), 2), "log10 copies/mL\n")
    cat("5-year eGFR slope:",
        round((d$EGFRr[which.min(abs(d$time - 1825))] -
               d$EGFRr[which.min(abs(d$time - 365))]) / 4, 2),
        "mL/min/1.73 m2/year\n")
  })

  ## ---- 2  PK --------------------------------------------------------------
  output$p_pk <- renderPlot({
    d <- sim()
    p1 <- ggplot(d, aes(yrs(time), TACC0r)) + band(5, 10) +
      geom_line(colour = PAL[1], linewidth = 0.7) +
      labs(x = NULL, y = "Tacrolimus trough (ng/mL)",
           title = "TAC C0 (green = 5-10 ng/mL target range)") + THEME
    p2 <- ggplot(d, aes(yrs(time), MPAAUC)) + band(30, 60) +
      geom_line(colour = PAL[3], linewidth = 0.7) +
      labs(x = NULL, y = "MPA AUC0-12 (mg*h/L)", title = "MPA exposure (green = 30-60)") + THEME
    p3 <- long(d, c("CPREDr", "EVRC0", "CSAC0r"),
               c("Prednisolone-eq (ng/mL)", "Everolimus C0 (ng/mL)",
                 "Cyclosporine C0 (ng/mL)")) %>% gg(ncol = 3)
    p4 <- long(d, c("CATGr", "CBELAr", "CRTXr", "CTCZr", "CFZBr"),
               c("rATG (ug/mL)", "Belatacept (ug/mL)", "Rituximab (ug/mL)",
                 "anti-IL-6 (ug/mL)", "anti-CD38 (ug/mL)")) %>% gg(ncol = 5)
    gridExtra_like(list(p1, p2, p3, p4), c(1, 1, 1, 1))
  })

  ## ---- 3  cellular --------------------------------------------------------
  output$p_cell <- renderPlot({
    d <- sim()
    long(d, c("TN", "TACT", "TEFF", "TMEM", "TREG", "TINF", "IL2", "IFNG", "AG"),
         c("Naive alloreactive T (TN)", "Activated blasts (TACT)",
           "Effector T (TEFF)", "Memory T (TMEM)", "Treg", "Graft-infiltrating T (TINF)",
           "IL-2", "IFN-gamma", "Antigen presentation (AG)")) %>%
      gg(ncol = 3, ylab = "relative units")
  })

  ## ---- 4  humoral ---------------------------------------------------------
  output$p_dsa <- renderPlot({
    d <- sim()
    p1 <- ggplot(d, aes(yrs(time), DSA)) +
      geom_hline(yintercept = 1000, linetype = 2, colour = "#9c2d57") +
      geom_line(colour = PAL[2], linewidth = 0.8) +
      labs(x = NULL, y = "DSA (MFI)",
           title = "Donor-specific antibody (dashed line = conventional detection threshold, 1000 MFI)") + THEME
    p2 <- long(d, c("BN", "BGC", "PB", "LLPC", "C4D", "NKACT", "ENDO", "MVI"),
               c("Naive/memory B (BN)", "Germinal-centre B (BGC)",
                 "Plasmablasts (PB)", "Long-lived plasma cells (LLPC)",
                 "C4d deposition", "Activated NK", "Endothelial activation",
                 "Microvascular inflammation")) %>% gg(ncol = 4)
    gridExtra_like(list(p1, p2), c(1, 2))
  })

  ## ---- 5  Banff -----------------------------------------------------------
  output$p_banff <- renderPlot({
    d <- sim()
    dd <- long(d, c("BANFF_T", "BANFF_I", "BANFF_G", "BANFF_PTC",
                    "BANFF_CG", "BANFF_CI", "BANFF_AH", "BANFF_C4D"),
               c("t (tubulitis)", "i (interstitial)", "g (glomerulitis)",
                 "ptc (capillaritis)", "cg (transplant glomerulopathy)",
                 "ci (fibrosis)", "ah (arteriolar hyalinosis)", "C4d"))
    ggplot(dd, aes(yrs(time), value)) +
      geom_hline(yintercept = 1, linetype = 3, colour = "grey60") +
      geom_line(linewidth = 0.7, colour = PAL[4]) +
      facet_wrap(~name, ncol = 4) +
      coord_cartesian(ylim = c(0, 3)) +
      labs(x = "Time since transplant (years)", y = "Banff score (0-3)",
           title = "Banff lesion score over time") + THEME
  })

  ## ---- 6  function --------------------------------------------------------
  output$p_func <- renderPlot({
    d <- sim()
    p1 <- ggplot(d, aes(yrs(time), EGFRr)) +
      geom_hline(yintercept = 15, linetype = 2, colour = "#9b3535") +
      geom_line(colour = PAL[1], linewidth = 0.8) +
      labs(x = NULL, y = "eGFR (mL/min/1.73 m2)",
           title = "Graft function (dashed line = dialysis-resumption threshold, 15)") + THEME
    p2 <- long(d, c("SCRDL", "PROT", "CFDNA", "GS"),
               c("Serum creatinine (mg/dL)", "Proteinuria, UPCR (g/g)",
                 "dd-cfDNA (%)", "Death-censored graft survival probability")) %>% gg(ncol = 4)
    gridExtra_like(list(p1, p2), c(1, 1))
  })

  ## ---- 7  infection -------------------------------------------------------
  output$p_inf <- renderPlot({
    d <- sim()
    p1 <- ggplot(d, aes(yrs(time), BKLOG)) +
      geom_hline(yintercept = 4, linetype = 2, colour = "#9b3535") +
      geom_line(colour = PAL[3], linewidth = 0.8) +
      labs(x = NULL, y = "BK (log10 copies/mL)",
           title = "BK polyomavirus (dashed line = 10^4, pre-emptive-reduction threshold)") + THEME
    p2 <- long(d, c("CMVLOG", "WBC", "GLU", "NISAUC", "BKVAN", "NISr"),
               c("CMV (log10 IU/mL)", "Leucocytes (10^9/L)", "Fasting glucose index (mmol/L)",
                 "Cumulative immunosuppression exposure (NIS-days)", "BK nephropathy burden",
                 "Net immunosuppression (NIS)")) %>% gg(ncol = 3)
    gridExtra_like(list(p1, p2), c(1, 1.4))
  })

  ## ---- 8  the U-curve -----------------------------------------------------
  ucurve <- eventReactive(input$sweep, {
    base <- pars()
    base$BKSUSC <- 1        # a patient who CAN reactivate BK - otherwise there
    base$BKSCREEN <- 0      # is no other arm to the trade-off
    tg <- seq(2, 14, by = 1)
    withProgress(message = "Sweeping immunosuppression intensity...", value = 0, {
      do.call(rbind, lapply(seq_along(tg), function(i) {
        incProgress(1 / length(tg))
        p <- base; p$TACTGT1 <- tg[i] + 2; p$TACTGT2 <- tg[i]; p$TACTGT3 <- tg[i]
        s <- mod %>% param(p) %>%
          mrgsim(end = 1825, delta = 5, hmax = 3, atol = 1e-8, rtol = 1e-5) %>%
          as_tibble()
        data.frame(
          target   = tg[i],
          NIS      = mean(s$NISr),
          rejection = max(s$BANFF_T) + max(s$MVISUM) / 2,
          infection = max(s$BKLOG),
          eGFR5     = s$EGFRr[which.min(abs(s$time - 1825))],
          survival  = s$GS[which.min(abs(s$time - 1825))]
        )
      }))
    })
  })

  output$p_ucurve <- renderPlot({
    u <- ucurve()
    dd <- u %>%
      select(NIS, `Rejection burden (Banff t + g+ptc/2)` = rejection,
             `Infection burden (peak BK log10)` = infection,
             `5-year eGFR` = eGFR5, `5-year graft survival probability` = survival) %>%
      pivot_longer(-NIS)
    ggplot(dd, aes(NIS, value)) +
      geom_line(linewidth = 0.9, colour = PAL[2]) +
      geom_point(size = 2, colour = PAL[2]) +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "Mean net immunosuppression (NIS)", y = NULL,
           title = "The U-shaped trade-off in immunosuppression intensity",
           subtitle = "Go left and you get rejection; go right and you get the virus — the optimum lies between the two") +
      THEME
  })

  ## ---- 9  comparison ------------------------------------------------------
  observeEvent(input$save, {
    s <- sim()
    s$run <- paste0("run ", length(saved()) + 1, ": ",
                    input$regimen, "/", input$induction,
                    " TGT", input$TACTGT2,
                    if (input$ADHFINAL < 1) paste0(" ADH", input$ADHFINAL) else "",
                    if (input$BKSUSC) " BK" else "")
    saved(c(saved(), list(s)))
  })
  observeEvent(input$clear, saved(list()))

  output$p_cmp <- renderPlot({
    if (!length(saved())) return(NULL)
    d <- bind_rows(saved())
    d %>%
      select(time, run, eGFR = EGFRr, DSA, `Banff t` = BANFF_T,
             `g+ptc` = MVISUM, `BK log10` = BKLOG, NIS = NISr) %>%
      pivot_longer(-c(time, run)) %>%
      ggplot(aes(yrs(time), value, colour = run)) +
      geom_line(linewidth = 0.7) +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      scale_colour_manual(values = rep(PAL, 3)) +
      labs(x = "Time since transplant (years)", y = NULL) + THEME +
      guides(colour = guide_legend(ncol = 2))
  })

  output$t_cmp <- renderTable({
    if (!length(saved())) return(NULL)
    do.call(rbind, lapply(saved(), function(s) {
      i <- function(day) which.min(abs(s$time - day))
      data.frame(
        run = s$run[1],
        `eGFR 1 year` = round(s$EGFRr[i(365)], 1),
        `eGFR 5 years` = round(s$EGFRr[i(1825)], 1),
        `Slope/year` = round((s$EGFRr[i(1825)] - s$EGFRr[i(365)]) / 4, 2),
        `Peak Banff t` = round(max(s$BANFF_T), 2),
        `DSA 5 years` = round(s$DSA[i(1825)]),
        `cg 5 years` = round(s$BANFF_CG[i(1825)], 2),
        `BK peak` = round(max(s$BKLOG), 2),
        `Survival probability 5 years` = round(s$GS[i(1825)], 3),
        check.names = FALSE
      )
    }))
  })

  ## ---- 10  parameters -----------------------------------------------------
  output$t_par <- renderTable({
    b <- src[(grep("^\\$PARAM", src)[1] + 1):(grep("^\\$CMT", src)[1] - 1)]
    b <- b[grepl(":", b)]
    parts <- strsplit(b, "\\s*:\\s*")
    data.frame(
      parameter   = vapply(parts, `[`, "", 1),
      value       = vapply(parts, `[`, "", 2),
      description = vapply(parts, function(x) paste(x[-c(1, 2)], collapse = " : "), "")
    )
  })
}

shinyApp(ui, server)
