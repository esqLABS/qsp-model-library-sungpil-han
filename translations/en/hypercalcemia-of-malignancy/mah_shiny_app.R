## =====================================================================
##  mah_shiny_app.R
##  Hypercalcaemia of Malignancy — interactive QSP dashboard
##  10 tabs · driven by mah_mrgsolve_model.R
##
##  Hypercalcaemia of Malignancy — QSP Dashboard
## ---------------------------------------------------------------------
##  The app is organised around the model's three theses rather than
##  around its compartments:
##
##    THESIS 1  three arms, three clocks   -> tabs 3, 4, 5
##    THESIS 2  the ceiling and the fold   -> tabs 6, 7
##    THESIS 3  same number, different     -> tabs 2, 8
##              mechanism, different drug
##
##  Run with:  shiny::runApp("mah_shiny_app.R")
##  Requires:  shiny, mrgsolve, dplyr, tidyr, ggplot2, DT
## =====================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)
library(gridExtra)

## The model file defines the mrgsolve object `mah`, the dose constants, the
## helpers run_mah()/ev_bolus()/run_scenario() and the SCEN scenario list.
## Look for it next to this file first, then in the working directory.
.here <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) ".")
.model <- file.path(.here, "mah_mrgsolve_model.R")
if (!file.exists(.model)) .model <- "mah_mrgsolve_model.R"
source(.model)

PAL <- c(saline = "#0D47A1", calcitonin = "#E65100", zoledronate = "#004D40",
         denosumab = "#827717", prednisolone = "#880E4F",
         cinacalcet = "#4A148C", dialysis = "#01579B", none = "#616161")

theme_mah <- function() {
  theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(colour = "grey35", size = 10),
          strip.text = element_text(face = "bold"),
          legend.position = "bottom")
}

band <- function(lo, hi, fill = "#C8E6C9") {
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = lo, ymax = hi,
           fill = fill, alpha = 0.35)
}

## =====================================================================
##  UI
## =====================================================================
ui <- fluidPage(
  titlePanel(HTML(paste0(
    "<b>Hypercalcaemia of Malignancy QSP Dashboard</b>",
    "<br><span style='font-size:14px;color:#555'>Hypercalcaemia of Malignancy — ",
    "three arms, three clocks, and a ceiling that the disease pushes down</span>"))),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      h4("① Patient & tumour"),
      selectInput("aeti", "Mechanism",
                  choices = c("Humoral / PTHrP (≈80%)"        = "HHM",
                              "Local osteolytic (≈20%)"       = "LOH",
                              "Calcitriol / lymphoma (<1%)"   = "CTD",
                              "Ectopic PTH / parathyroid ca"  = "EPT",
                              "Healthy control"               = "NONE"),
                  selected = "HHM"),
      sliderInput("alb", "Serum albumin (g/L)", 15, 48, 40, step = 1),
      sliderInput("nfunc", "Residual nephron mass (fraction)",
                  0.15, 1.0, 1.0, step = 0.05),
      checkboxInput("immob", "Immobilisation", FALSE),
      checkboxInput("thia", "Thiazide", FALSE),
      sliderInput("tumg", "Tumour growth rate (1/day)",
                  0, 0.08, 0.030, step = 0.005),

      hr(),
      h4("② Therapy — starts day 12"),
      sliderInput("iv", "Saline (L/day, first 72 h)",
                  0, 8, 4.8, step = 0.2),
      sliderInput("ivdays", "High-rate infusion period (days)", 0, 7, 3, step = 1),
      checkboxInput("ctn", "Calcitonin 4 IU/kg q12h × 48 h", FALSE),
      selectInput("anti", "Antiresorptive",
                  choices = c("none", "zoledronate", "denosumab"),
                  selected = "zoledronate"),
      conditionalPanel("input.anti == 'zoledronate'",
                       sliderInput("zadose", "Zoledronate dose (mg)",
                                   0.5, 16, 4, step = 0.5)),
      checkboxInput("prd", "Prednisolone 60 mg/day × 20 d", FALSE),
      checkboxInput("cin", "Cinacalcet 90 mg bd", FALSE),
      checkboxInput("fur", "Furosemide 40 mg q6h × 72 h", FALSE),
      checkboxInput("dial", "Haemodialysis (Ca 1.0, 4 h × 2)", FALSE),
      sliderInput("chemo", "Antitumour treatment effect (1/day)",
                  0, 0.25, 0, step = 0.01),

      hr(),
      sliderInput("end", "Observation horizon (days)", 30, 120, 60, step = 5),
      actionButton("go", "Run simulation", class = "btn-primary",
                   width = "100%"),
      br(), br(),
      helpText(HTML(paste0(
        "<small>In every scenario the tumour appears on day 0 and <b>therapy begins on day 12</b>",
        ". The mechanism-specific secretion constants are calibrated so that total calcium at day 12 is 3.55 mmol/L, ",
        "so that the treatment comparisons are made at <b>the same calcium value</b>.<br>",
        "All mechanisms are calibrated to the SAME calcium at the moment ",
        "therapy starts, so the comparisons are like-for-like.</small>")))
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs", type = "tabs",

        ## ---- 1 ---------------------------------------------------------
        tabPanel("① Patient profile",
                 br(),
                 fluidRow(column(12, uiOutput("kpi"))),
                 plotOutput("p_profile", height = "560px"),
                 helpText(HTML(paste0(
                   "Shows total calcium, ionised calcium, symptom score, and ECF volume together. ",
                   "<b>Only ionised calcium is the value the CaSR reads</b>, and lowering albumin ",
                   "creates a state where total calcium is low but ionised calcium is high.")))),

        ## ---- 2 ---------------------------------------------------------
        tabPanel("② Mechanism discrimination panel",
                 br(),
                 plotOutput("p_mech", height = "480px"),
                 br(),
                 DTOutput("t_mech"),
                 helpText(HTML(paste0(
                   "Comparison of the four mechanisms at day 12 (just before therapy). <b>Even at the same total calcium, ",
                   "bone resorption and urinary calcium are completely different</b> — because PTHrP ",
                   "opens the tap on bone while simultaneously closing the drain in the kidney.")))),

        ## ---- 3 ---------------------------------------------------------
        tabPanel("③ Three arms",
                 br(),
                 plotOutput("p_arms", height = "520px"),
                 helpText(HTML(paste0(
                   "dCa/dt = J_res − J_form + J_gut − U_Ca. In the panel below, ",
                   "calcium keeps rising for as long as the <b>net input</b> stays above zero. ",
                   "Therapy touches <b>only one</b> of these four terms at a time.")))),

        ## ---- 4 ---------------------------------------------------------
        tabPanel("④ Drug PK",
                 br(),
                 plotOutput("p_pk", height = "560px"),
                 helpText(HTML(paste0(
                   "Zoledronic acid reaches its site of action <b>not via plasma but via the bone surface</b>, ",
                   "and the rate of delivery is proportional to the rate of bone resorption — the drug ",
                   "stopping itself is the very process that is also its delivery route.")))),

        ## ---- 5 ---------------------------------------------------------
        tabPanel("⑤ Bone / RANKL",
                 br(),
                 plotOutput("p_bone", height = "560px"),
                 helpText(HTML(paste0(
                   "PTH1R signalling raises RANKL while <b>simultaneously lowering</b> OPG. ",
                   "The two errors adding in the same direction is what uncoupling is, ",
                   "and is why P1NP fails to keep pace with CTX.")))),

        ## ---- 6 ---------------------------------------------------------
        tabPanel("⑥ Renal excretion ceiling",
                 br(),
                 fluidRow(
                   column(6, sliderInput("ceil_nf", "Nephron mass",
                                         0.2, 1.0, 1.0, step = 0.1)),
                   column(6, sliderInput("ceil_iv", "Saline (L/day)",
                                         0, 6, 0, step = 0.6))),
                 plotOutput("p_ceiling", height = "460px"),
                 helpText(HTML(paste0(
                   "Read after holding total calcium fixed and letting everything else equilibrate: ",
                   "<b>the steady-state urinary calcium</b>. This curve does not simply saturate — ",
                   "it <b>turns over and comes back down</b>. The peak is the excretion ceiling, and once net input ",
                   "exceeds that ceiling, no steady state exists.")))),

        ## ---- 7 ---------------------------------------------------------
        tabPanel("⑦ Saddle-node bifurcation (fold)",
                 br(),
                 plotOutput("p_fold", height = "440px"),
                 br(),
                 DTOutput("t_fold"),
                 helpText(HTML(paste0(
                   "By increasing the unregulated calcium input (J_EXO), we look for ",
                   "the point where the steady state disappears. <b>A hypercalcaemic crisis is not a large number ",
                   "but a bifurcation</b>, and what saline does is not to remove calcium ",
                   "but to <b>raise the ceiling back up</b>.")))),

        ## ---- 8 ---------------------------------------------------------
        tabPanel("⑧ Scenario comparison",
                 br(),
                 checkboxGroupInput(
                   "cmp", "Scenarios to overlay",
                   choices = c("saline alone" = "s03",
                               "saline + calcitonin" = "s04",
                               "saline + zoledronate" = "s05",
                               "guideline triple" = "s06",
                               "zoledronate, no saline" = "s07",
                               "saline + denosumab" = "s08",
                               "osteolytic + ZA" = "s11",
                               "lymphoma + ZA" = "s12",
                               "lymphoma + prednisolone" = "s13",
                               "inadequate saline + furosemide" = "s16",
                               "inadequate saline alone" = "s23",
                               "ZA + antitumour therapy" = "s21"),
                   selected = c("s03", "s05", "s06", "s08"), inline = TRUE),
                 plotOutput("p_cmp", height = "440px"),
                 br(), DTOutput("t_cmp")),

        ## ---- 9 ---------------------------------------------------------
        tabPanel("⑨ Biomarkers",
                 br(),
                 plotOutput("p_bio", height = "560px"),
                 helpText(HTML(paste0(
                   "CTX moves within 24–48 hours, and P1NP only ",
                   "follows days to weeks later. Using an antiresorptive drug drops <b>both</b>, ",
                   "which is the price of coupling.")))),

        ## ---- 10 --------------------------------------------------------
        tabPanel("⑩ Clinical endpoints · safety",
                 br(),
                 plotOutput("p_clin", height = "520px"),
                 br(),
                 htmlOutput("safety"))
      )
    )
  )
)

## =====================================================================
##  SERVER
## =====================================================================
server <- function(input, output, session) {

  build <- reactive({
    input$go
    isolate({
      tx <- 12
      pmod <- list(ALB = input$alb, TUM_G = input$tumg,
                   IMMOB = as.numeric(input$immob),
                   THIAZIDE = as.numeric(input$thia),
                   CHEMO = input$chemo)
      imod <- list(N_func = input$nfunc)

      iv <- NULL
      if (input$iv > 0 && input$ivdays > 0) {
        iv <- data.frame(
          time = c(tx, tx + input$ivdays, tx + input$ivdays + 4),
          IV_RATE = c(input$iv, input$iv / 2, 0))
      }
      if (input$dial) {
        iv <- merge_segs(iv, data.frame(time = tx + c(0, 1/6, 1, 1 + 1/6),
                                        DIAL = c(1, 0, 1, 0)))
      }

      ev <- NULL
      add <- function(e) if (is.null(ev)) e else rbind(ev, e)
      if (input$ctn) ev <- add(ev_bolus(tx + c(0, .5, 1, 1.5), "CTN_SC", CTN_4IU))
      if (input$anti == "zoledronate")
        ev <- add(ev_bolus(tx, "ZOL_C", 14.70 * input$zadose / 4))
      if (input$anti == "denosumab")  ev <- add(ev_bolus(tx, "DMB_SC", DMB_120))
      if (input$prd) ev <- add(ev_bolus(tx + 0:19, "PRD_C", PRD_60))
      if (input$cin) ev <- add(ev_bolus(tx + seq(0, 29.5, .5), "CIN_C", CIN_90))
      if (input$fur) ev <- add(ev_bolus(tx + seq(0, 2.75, .25), "FUR_C", FUR_40))

      list(out = run_mah(aeti = input$aeti, events = ev, iv = iv,
                         pmod = pmod, imod = imod, end = input$end),
           tx = tx)
    })
  })

  ## ---- KPI strip ------------------------------------------------------
  output$kpi <- renderUI({
    d <- build()$out; tx <- build()$tx
    post <- dplyr::filter(d, time >= tx)
    nad <- min(post$Ca_tot); tnad <- post$time[which.min(post$Ca_tot)] - tx
    tnorm <- suppressWarnings(min(post$time[post$Ca_tot < 2.70]) - tx)
    box <- function(lab, val, col) sprintf(
      paste0("<div style='display:inline-block;min-width:150px;margin:4px;",
             "padding:9px 13px;border-left:5px solid %s;background:#F5F6F8'>",
             "<div style='font-size:11px;color:#666'>%s</div>",
             "<div style='font-size:20px;font-weight:700'>%s</div></div>"),
      col, lab, val)
    HTML(paste0(
      box("Ca at therapy start", sprintf("%.2f", post$Ca_tot[1]), "#B71C1C"),
      box("ionised at start", sprintf("%.2f", post$iCa[1]), "#B71C1C"),
      box("nadir (total)", sprintf("%.2f", nad), "#004D40"),
      box("day of nadir", sprintf("%.1f", tnad), "#004D40"),
      box("days to Ca < 2.70",
          if (is.finite(tnorm)) sprintf("%.1f", tnorm) else "not reached",
          "#0D47A1"),
      box("min GFR (mL/min)", sprintf("%.0f", min(post$eGFR_ml)), "#E65100")))
  })

  ## ---- 1. profile -----------------------------------------------------
  output$p_profile <- renderPlot({
    d <- build()$out; tx <- build()$tx
    long <- d %>%
      select(time, `Total calcium (mmol/L)` = Ca_tot,
             `Ionised calcium (mmol/L)` = iCa,
             `Payne-corrected (mmol/L)` = Ca_corr,
             `Symptom score (0-10)` = SYM,
             `ECF volume (fraction)` = vol,
             `eGFR (mL/min)` = eGFR_ml) %>%
      pivot_longer(-time)
    ggplot(long, aes(time, value)) +
      geom_vline(xintercept = tx, linetype = 2, colour = "#B71C1C") +
      geom_line(linewidth = 0.9, colour = "#1A237E") +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "day", y = NULL,
           title = "Patient profile",
           subtitle = "red dashed line = therapy start (day 12)") +
      theme_mah()
  })

  ## ---- 2. mechanism panel --------------------------------------------
  mech_tbl <- reactive({
    rows <- lapply(c("HHM", "LOH", "CTD", "EPT"), function(k) {
      o <- run_mah(aeti = k, end = 12.05, delta = 0.05)
      r <- o[nrow(o), ]
      data.frame(mechanism = k, Ca_tot = r$Ca_tot, iCa = r$iCa,
                 PTH = r$PTH, PTHrP = r$PTHRP, D125 = r$D125, D25 = r$D25,
                 phosphate = r$Pi_c, J_res = r$J_res, J_gut = r$J_gut,
                 U_Ca = r$U_Ca, FE_Ca = 100 * r$FE_Ca, Tm = r$Tm_dct,
                 CTX = r$CTX, P1NP = r$P1NP)
    })
    do.call(rbind, rows)
  })

  output$p_mech <- renderPlot({
    m <- mech_tbl() %>%
      select(mechanism, `bone resorption\n(mmol/day)` = J_res,
             `net gut absorption\n(mmol/day)` = J_gut,
             `urinary calcium\n(mmol/day)` = U_Ca,
             `FE_Ca (%)` = FE_Ca,
             `distal Ca Tm\n(mmol/day)` = Tm,
             `phosphate (mmol/L)` = phosphate) %>%
      pivot_longer(-mechanism)
    ggplot(m, aes(mechanism, value, fill = mechanism)) +
      geom_col(width = 0.68) +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      scale_fill_manual(values = c(HHM = "#B71C1C", LOH = "#3E2723",
                                   CTD = "#FF6F00", EPT = "#01579B")) +
      labs(x = NULL, y = NULL, title = "Four mechanisms, same time point (day 12, before treatment)",
           subtitle = paste("HHM = PTHrP humoral · LOH = local osteolytic ·",
                            "CTD = calcitriol · EPT = ectopic PTH")) +
      theme_mah() + theme(legend.position = "none")
  })

  output$t_mech <- renderDT({
    datatable(mech_tbl() %>% mutate(across(where(is.numeric), ~round(.x, 3))),
              options = list(dom = "t", scrollX = TRUE), rownames = FALSE)
  })

  ## ---- 3. three arms --------------------------------------------------
  output$p_arms <- renderPlot({
    d <- build()$out; tx <- build()$tx
    d2 <- d %>% mutate(net = J_res - J_form + J_gut - U_Ca)
    a <- d2 %>% select(time, `bone resorption` = J_res,
                       `bone formation` = J_form,
                       `net gut absorption` = J_gut,
                       `urinary excretion` = U_Ca) %>% pivot_longer(-time)
    p1 <- ggplot(a, aes(time, value, colour = name)) +
      geom_vline(xintercept = tx, linetype = 2, colour = "#B71C1C") +
      geom_line(linewidth = 0.95) +
      scale_colour_manual(values = c("bone resorption" = "#6D4C41",
                                     "bone formation" = "#A1887F",
                                     "net gut absorption" = "#33691E",
                                     "urinary excretion" = "#006064")) +
      labs(x = NULL, y = "mmol/day", colour = NULL,
           title = "The three arms of the balance") + theme_mah()
    p2 <- ggplot(d2, aes(time, net)) +
      geom_hline(yintercept = 0, colour = "grey40") +
      geom_vline(xintercept = tx, linetype = 2, colour = "#B71C1C") +
      geom_area(fill = "#9FA8DA", alpha = 0.6) + geom_line(linewidth = 0.9) +
      labs(x = "day", y = "mmol/day",
           subtitle = "net input into the ECF pool; calcium rises while this is above zero") +
      theme_mah()
    gridExtra::grid.arrange(p1, p2, heights = c(3, 2))
  })

  ## ---- 4. PK ----------------------------------------------------------
  output$p_pk <- renderPlot({
    d <- build()$out; tx <- build()$tx
    long <- d %>%
      select(time, `Zoledronate plasma (umol)` = ZOL_C,
             `Zoledronate on bone surface (umol)` = ZOL_B,
             `Intra-osteoclast zoledronate (umol)` = ZOL_OC,
             `Denosumab free (pmol/L)` = DMB_C,
             `Calcitonin biophase (IU/L)` = CTN_E,
             `Calcitonin receptor R_CT` = R_CT,
             `Prednisolone biophase (mg/L)` = PRD_E,
             `Cinacalcet (ng/mL)` = CIN_C) %>%
      pivot_longer(-time)
    ggplot(long, aes(time, value)) +
      geom_vline(xintercept = tx, linetype = 2, colour = "#B71C1C") +
      geom_line(linewidth = 0.9, colour = "#004D40") +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "day", y = NULL, title = "Drug PK and site-of-action exposure") + theme_mah()
  })

  ## ---- 5. bone --------------------------------------------------------
  output$p_bone <- renderPlot({
    d <- build()$out; tx <- build()$tx
    long <- d %>% mutate(ratio = RKL / pmax(OPG, 1e-6)) %>%
      select(time, `free RANKL (pmol/L)` = RKL, `OPG (pmol/L)` = OPG,
             `RANKL / OPG` = ratio, `osteoclasts (relative)` = OC,
             `osteoblasts (relative)` = OB, `PTH1R signal` = PS,
             `CTX (ng/mL)` = CTX, `P1NP (ug/L)` = P1NP) %>%
      pivot_longer(-time)
    ggplot(long, aes(time, value)) +
      geom_vline(xintercept = tx, linetype = 2, colour = "#B71C1C") +
      geom_line(linewidth = 0.9, colour = "#3E2723") +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "day", y = NULL, title = "Bone remodelling unit") +
      theme_mah()
  })

  ## ---- 6. the ceiling -------------------------------------------------
  ceiling_curve <- reactive({
    cas <- seq(2.4, 4.2, by = 0.1)
    res <- vapply(cas, function(ca) {
      o <- mah %>%
        param(CLAMP_CA = ca, N_func = input$ceil_nf, IV_RATE = input$ceil_iv) %>%
        init(A_Ca = ca * 15, N_func = input$ceil_nf) %>%
        mrgsim_df(end = 250, delta = 5)
      r <- o[nrow(o), ]
      c(U_Ca = r$U_Ca, GFR = r$eGFR_ml, vol = r$vol)
    }, numeric(3))
    data.frame(Ca = cas, U_Ca = res["U_Ca", ], GFR = res["GFR", ],
               vol = res["vol", ])
  })

  output$p_ceiling <- renderPlot({
    cc <- ceiling_curve()
    top <- cc[which.max(cc$U_Ca), ]
    ggplot(cc, aes(Ca, U_Ca)) +
      geom_line(linewidth = 1.2, colour = "#4527A0") +
      geom_point(size = 2, colour = "#4527A0") +
      geom_point(data = top, size = 5, shape = 21, fill = "#FFC107",
                 colour = "#311B92", stroke = 1.3) +
      geom_text(data = top, aes(label = sprintf("ceiling %.1f mmol/day at Ca %.1f",
                                                U_Ca, Ca)),
                vjust = -1.4, size = 4.2, colour = "#311B92") +
      labs(x = "total plasma calcium held at (mmol/L)",
           y = "steady-state urinary calcium (mmol/day)",
           title = "The renal escape-valve ceiling",
           subtitle = paste0("nephron mass ", input$ceil_nf,
                             " · saline ", input$ceil_iv, " L/day",
                             " — the curve TURNS OVER; it does not merely saturate")) +
      expand_limits(y = 0) + theme_mah()
  })

  ## ---- 7. the fold ----------------------------------------------------
  fold_tbl <- reactive({
    conds <- list(
      list("baseline (no support)", list()),
      list("+ saline 1.2 L/day", list(IV_RATE = 1.2)),
      list("+ saline 2.4 L/day", list(IV_RATE = 2.4)),
      list("distal Tm +60% (PTHrP renal effect)", list(TM_DCT = 54.4)),
      list("thiazide", list(THIAZIDE = 1)),
      list("nephron mass 60%", list()),
      list("nephron mass 30%", list()))
    nf <- c(1, 1, 1, 1, 1, 0.6, 0.3)
    collapses <- function(je, pm, n) {
      o <- mah %>% param(c(list(J_EXO = je), pm)) %>% init(N_func = n) %>%
        mrgsim_df(end = 150, delta = 1)
      any(o$vol <= 0.62 | o$Ca_tot >= 5 | !is.finite(o$Ca_tot))
    }
    out <- lapply(seq_along(conds), function(i) {
      lo <- 0; hi <- 80
      for (k in 1:11) {
        m <- (lo + hi) / 2
        if (collapses(m, conds[[i]][[2]], nf[i])) hi <- m else lo <- m
      }
      data.frame(condition = conds[[i]][[1]], fold_mmol_per_day = round(lo, 2))
    })
    do.call(rbind, out)
  })

  output$p_fold <- renderPlot({
    ft <- fold_tbl()
    ft$condition <- factor(ft$condition, levels = rev(ft$condition))
    ggplot(ft, aes(fold_mmol_per_day, condition)) +
      geom_col(fill = "#7E57C2", width = 0.65) +
      geom_vline(xintercept = ft$fold_mmol_per_day[1], linetype = 2,
                 colour = "#B71C1C") +
      geom_text(aes(label = sprintf("%.1f", fold_mmol_per_day)),
                hjust = -0.2, size = 4) +
      expand_limits(x = max(ft$fold_mmol_per_day) * 1.15) +
      labs(x = "unregulated calcium input at which no steady state exists (mmol/day)",
           y = NULL, title = "The saddle-node bifurcation point (the fold)",
           subtitle = "red line = baseline; anything that moves the bar LEFT is a precipitant") +
      theme_mah()
  })
  output$t_fold <- renderDT(datatable(fold_tbl(), options = list(dom = "t"),
                                      rownames = FALSE))

  ## ---- 8. scenario comparison ----------------------------------------
  cmp_runs <- reactive({
    ids <- input$cmp
    if (!length(ids)) return(NULL)
    do.call(rbind, lapply(ids, function(i) {
      o <- run_scenario(i)
      o$scenario <- SCEN[[i]][[1]]
      o
    }))
  })

  output$p_cmp <- renderPlot({
    r <- cmp_runs(); req(r)
    ggplot(r, aes(time, Ca_tot, colour = scenario)) +
      band(2.15, 2.55) +
      geom_hline(yintercept = 2.70, linetype = 3) +
      geom_vline(xintercept = 12, linetype = 2, colour = "#B71C1C") +
      geom_line(linewidth = 1.0) +
      coord_cartesian(ylim = c(1.8, 5.0)) +
      labs(x = "day", y = "total calcium (mmol/L)", colour = NULL,
           title = "Scenario comparison",
           subtitle = "green band = reference interval · dotted line = 2.70 mmol/L") +
      theme_mah() + guides(colour = guide_legend(ncol = 2))
  })

  output$t_cmp <- renderDT({
    r <- cmp_runs(); req(r)
    s <- r %>% filter(time >= 12) %>% group_by(scenario) %>%
      summarise(`Ca at start` = round(first(Ca_tot), 3),
                nadir = round(min(Ca_tot), 3),
                `day of nadir` = round(time[which.min(Ca_tot)] - 12, 2),
                `days < 2.70` = round(sum(Ca_tot < 2.70) * 0.05, 2),
                `CTX nadir` = round(min(CTX), 3),
                `min eGFR` = round(min(eGFR_ml), 1), .groups = "drop")
    datatable(s, options = list(dom = "t", scrollX = TRUE), rownames = FALSE)
  })

  ## ---- 9. biomarkers --------------------------------------------------
  output$p_bio <- renderPlot({
    d <- build()$out; tx <- build()$tx
    long <- d %>%
      select(time, `CTX-I (ng/mL)` = CTX, `P1NP (ug/L)` = P1NP,
             `PTH (pmol/L)` = PTH, `PTHrP (pmol/L)` = PTHRP,
             `1,25(OH)2D (pmol/L)` = D125, `25(OH)D (nmol/L)` = D25,
             `FGF23 (pg/mL)` = FGF23, `phosphate (mmol/L)` = Pi_c,
             `magnesium (mmol/L)` = Mg_c) %>%
      pivot_longer(-time)
    ggplot(long, aes(time, value)) +
      geom_vline(xintercept = tx, linetype = 2, colour = "#B71C1C") +
      geom_line(linewidth = 0.9, colour = "#F57F17") +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "day", y = NULL, title = "Biomarkers") + theme_mah()
  })

  ## ---- 10. clinical endpoints ----------------------------------------
  output$p_clin <- renderPlot({
    d <- build()$out; tx <- build()$tx
    long <- d %>%
      select(time, `Symptom score (0-10)` = SYM,
             `Adapted set point (mmol/L)` = ADAPT,
             `QTc (ms)` = QTc, `KDIGO AKI stage` = KDIGO,
             `Urine output (L/day)` = Q_u, `AQP2 (fraction)` = AQP2,
             `Nephron mass (fraction)` = N_func,
             `Cumulative iCa exposure > 1.40` = AUC_CA) %>%
      pivot_longer(-time)
    ggplot(long, aes(time, value)) +
      geom_vline(xintercept = tx, linetype = 2, colour = "#B71C1C") +
      geom_line(linewidth = 0.9, colour = "#263238") +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "day", y = NULL, title = "Clinical endpoints") + theme_mah()
  })

  output$safety <- renderUI({
    d <- build()$out; tx <- build()$tx
    post <- dplyr::filter(d, time >= tx)
    msgs <- character(0)
    if (min(post$Ca_tot) < 2.05)
      msgs <- c(msgs, sprintf(paste0("<b style='color:#B71C1C'>Hypocalcaemia alert</b> ",
        "— total calcium fell to %.2f mmol/L. Antiresorptive over-suppression ",
        "is the model's own predicted iatrogenic harm."), min(post$Ca_tot)))
    if (min(post$vol) <= 0.62)
      msgs <- c(msgs, sprintf(paste0("<b style='color:#B71C1C'>Model valid-range ",
        "departure (collapse boundary)</b> — the ECF deficit exceeded 38%% on day ",
        "%.1f. Beyond this the trajectory is outside the calibrated region; ",
        "it is NOT a prediction of death."),
        post$time[which(post$vol <= 0.62)[1]]))
    if (max(post$Ca_tot) >= 5)
      msgs <- c(msgs, "<b style='color:#B71C1C'>Model valid-range departure</b> — total calcium exceeded 5.0 mmol/L.")
    if (min(post$K_c) < 3.0)
      msgs <- c(msgs, sprintf("<b style='color:#E65100'>Hypokalaemia</b> — potassium %.2f mmol/L.", min(post$K_c)))
    if (min(post$Mg_c) < 0.60)
      msgs <- c(msgs, sprintf("<b style='color:#E65100'>Hypomagnesaemia</b> — magnesium %.2f mmol/L.", min(post$Mg_c)))
    if (max(post$pHa) > 7.50)
      msgs <- c(msgs, sprintf(paste0("<b style='color:#33691E'>Metabolic alkalosis</b> — pH %.2f. ",
        "Note this LOWERS the ionised calcium: a protective arm inside the vicious cycle."),
        max(post$pHa)))
    if (min(post$N_func) < 0.9)
      msgs <- c(msgs, sprintf("<b style='color:#E65100'>Permanent nephron loss</b> — nephron mass %.2f.", min(post$N_func)))
    if (!length(msgs)) msgs <- "<span style='color:#1B5E20'>No safety signal triggered in this run.</span>"
    HTML(paste0("<div style='padding:10px;background:#F5F6F8;border-left:5px solid #263238'>",
                paste(msgs, collapse = "<br><br>"), "</div>"))
  })
}

shinyApp(ui, server)
