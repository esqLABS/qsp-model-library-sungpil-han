## =====================================================================
##  icp_shiny_app.R
##  Intrahepatic Cholestasis of Pregnancy (ICP) — QSP dashboard
##  10 tabs over icp_mrgsolve_model.R
##
##  DESIGN NOTE — why this app is laid out the way it is
##  ---------------------------------------------------
##  The app is built to make ONE thing hard to miss: the bile acid axis
##  and the itch axis are different systems, and the assay the clinic
##  uses reports neither of them cleanly.  So:
##    * Tab 2 always plots the assay value AND the endogenous fraction on
##      the same axes, because the gap between them is the drug effect
##      nobody attributes to the drug.
##    * Tab 4 (fetal) and Tab 6 (itch) can be driven from the same run
##      and deliberately show opposite responses to the same regimen.
##    * Tab 7 refuses to report a single "risk": it decomposes the
##      PITCHES composite so you can see which component moved.
##    * Tab 8 plots the delivery decision as two curves crossing, not as
##      a threshold, because that is what it is.
##
##  USAGE
##      install.packages(c("shiny","mrgsolve","dplyr","tidyr","ggplot2",
##                         "DT","patchwork"))
##      shiny::runApp("icp_shiny_app.R")
## =====================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

MODFILE <- "icp_mrgsolve_model"
mod <- mread_cache(MODFILE, ".")

MW_UDCA <- 392.57
ga2d <- function(ga) (ga - 20) * 7

## ---- susceptibility vectors -----------------------------------------
GENOTYPES <- list(
  "Normal pregnancy (wild type)"     = list(GBSEP = 1.000, GMDR3 = 1.000, GSULT = 1.00),
  "ICP · TBA <40 (mild, ~16)"        = list(GBSEP = 0.795, GMDR3 = 0.695, GSULT = 1.72),
  "ICP · TBA <40 (upper end, ~26)"   = list(GBSEP = 0.735, GMDR3 = 0.620, GSULT = 1.95),
  "ICP · TBA 40-99 (severe, ~58)"    = list(GBSEP = 0.650, GMDR3 = 0.500, GSULT = 2.40),
  "ICP · TBA >=100 (very severe, ~138)" = list(GBSEP = 0.520, GMDR3 = 0.350, GSULT = 3.10),
  "ABCB11 V444A homozygous"          = list(GBSEP = 0.550, GMDR3 = 0.950, GSULT = 2.20),
  "ABCB4 severe variant (MDR3 30%)"  = list(GBSEP = 0.920, GMDR3 = 0.300, GSULT = 2.20),
  "ATP8B1 variant"                   = list(GFIC1 = 0.680, GMDR3 = 0.800, GSULT = 2.40),
  "NR1H4(FXR) loss of function"      = list(GFXR = 0.550, GBSEP = 0.650, GMDR3 = 0.500, GSULT = 2.40),
  "Reduced placental transport (GP 0.6)" = list(GBSEP = 0.650, GMDR3 = 0.500, GSULT = 2.40, GP = 0.60)
)

## ---- dosing regimen builders ----------------------------------------
build_events <- function(udca_mg, udca_from, rif_mg, rif_from, chol_g, chol_from,
                         same_mg, ibat_umol, ntx_mg, ah_mg, bet_ga, del_ga) {
  e <- NULL
  add <- function(x) if (is.null(e)) x else c(e, x)
  if (udca_mg > 0 && udca_from < del_ga) {
    n <- max(1, round((del_ga - udca_from) * 14))
    e <- add(ev(amt = udca_mg / 2 / MW_UDCA * 1000, cmt = "UDDEP",
                ii = 0.5, addl = n - 1, time = ga2d(udca_from)))
  }
  if (rif_mg > 0 && rif_from < del_ga) {
    n <- max(1, round((del_ga - rif_from) * 14))
    e <- add(ev(amt = rif_mg / 2, cmt = "RIFDEP", ii = 0.5, addl = n - 1,
                time = ga2d(rif_from)))
  }
  if (chol_g > 0 && chol_from < del_ga) {
    n <- max(1, round((del_ga - chol_from) * 28))
    e <- add(ev(amt = chol_g / 4, cmt = "CHOLL", ii = 0.25, addl = n - 1,
                time = ga2d(chol_from)))
  }
  if (same_mg > 0) {
    n <- max(1, round(del_ga - 32))
    e <- add(ev(amt = same_mg, cmt = "SAMC", ii = 1, addl = n - 1,
                time = ga2d(32)))
  }
  if (ibat_umol > 0) {
    n <- max(1, round(del_ga - 32))
    e <- add(ev(amt = ibat_umol, cmt = "ODEV", ii = 1, addl = n - 1,
                time = ga2d(32)))
  }
  if (ntx_mg > 0) {
    n <- max(1, round(del_ga - 32))
    e <- add(ev(amt = ntx_mg, cmt = "NTXC", ii = 1, addl = n - 1,
                time = ga2d(32)))
  }
  if (ah_mg > 0) {
    n <- max(1, round((del_ga - 32) * 3))
    e <- add(ev(amt = ah_mg / 3, cmt = "AHC", ii = 1/3, addl = n - 1,
                time = ga2d(32)))
  }
  if (!is.na(bet_ga) && bet_ga > 0) {
    e <- add(ev(amt = 12, cmt = "BETC", ii = 1, addl = 1, time = ga2d(bet_ga)))
  }
  e
}

run_case <- function(gen, ev_obj, del_ga, twin = 1, extra = list()) {
  p <- c(gen, extra, list(TDEL = ga2d(del_ga), TWIN = twin))
  m <- mod %>% param(p)
  out <- if (is.null(ev_obj)) {
    m %>% mrgsim(end = ga2d(del_ga) + 21, delta = 0.5)
  } else {
    m %>% mrgsim(events = ev_obj, end = ga2d(del_ga) + 21, delta = 0.5)
  }
  as_tibble(out) %>% mutate(GA = 20 + time / 7,
                            phase = ifelse(time < ga2d(del_ga),
                                           "during pregnancy", "after delivery"))
}

THEME <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 10, colour = "grey35"))

## =====================================================================
##  UI
## =====================================================================
ui <- fluidPage(
  titlePanel("Intrahepatic cholestasis of pregnancy (ICP) — QSP dashboard"),
  p(style = "color:#555;",
    strong("This is a model for education and research."),
    " Do not use it for clinical decisions. The core claim of this model is that the two causal chains (the bile acid axis · the pruritus axis)",
    " respond to different drugs, and putting tab 4 and tab 6 side by side under",
    " the same prescription is what reveals that separation."),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Patient profile"),
      selectInput("gen", "Genetic susceptibility (susceptibility vector)",
                  choices = names(GENOTYPES), selected = names(GENOTYPES)[4]),
      checkboxInput("twin", "Twin pregnancy (sex hormone load ×1.55)", FALSE),
      sliderInput("del_ga", "Gestational age at delivery (GA, weeks)", 34, 41, 39, step = 0.5),
      numericInput("bet_ga", "Betamethasone administration week (0 = none)", 0, 0, 40, 0.5),
      hr(),
      h4("Bile acid axis therapy"),
      sliderInput("udca_mg", "UDCA (mg/day, split in 2)", 0, 2000, 0, step = 250),
      sliderInput("udca_from", "UDCA start week", 26, 38, 30, step = 0.5),
      sliderInput("chol_g", "Cholestyramine (g/day, split in 4)", 0, 24, 0, step = 4),
      sliderInput("chol_from", "Cholestyramine start week", 26, 38, 32, step = 0.5),
      sliderInput("ibat", "IBAT inhibitor (hypothetical, intestinal µmol/day)", 0, 8, 0, step = 1),
      sliderInput("same_mg", "SAMe (mg/day IV)", 0, 1600, 0, step = 200),
      checkboxInput("vitk", "Vitamin K 10 mg/day supplementation (from week 34)", FALSE),
      hr(),
      h4("Pruritus axis therapy"),
      sliderInput("rif_mg", "Rifampicin (mg/day, split in 2)", 0, 600, 0, step = 150),
      sliderInput("rif_from", "Rifampicin start week", 26, 38, 32, step = 0.5),
      sliderInput("ntx_mg", "Naltrexone (mg/day)", 0, 50, 0, step = 25),
      sliderInput("ah_mg", "Antihistamine (mg/day)", 0, 24, 0, step = 4),
      hr(),
      h4("Comparator (tab 9)"),
      selectInput("cmp_gen", "Comparator genetic susceptibility", choices = names(GENOTYPES),
                  selected = names(GENOTYPES)[4]),
      sliderInput("cmp_udca", "Comparator UDCA (mg/day)", 0, 2000, 1000, step = 250),
      sliderInput("cmp_rif", "Comparator rifampicin (mg/day)", 0, 600, 0, step = 150),
      sliderInput("cmp_del", "Comparator gestational age at delivery", 34, 41, 39, step = 0.5)
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs", type = "tabs",

        ## ---- 1 -------------------------------------------------------
        tabPanel(
          "1. Patient profile",
          br(),
          fluidRow(column(6, plotOutput("p_steroid", height = 300)),
                   column(6, plotOutput("p_inhib", height = 300))),
          hr(),
          h4("Summary at the current settings"),
          DTOutput("t_summary"),
          br(),
          div(style = "background:#f7f7fb;padding:10px;border-radius:6px;",
              strong("How to read this."),
              " The two graphs on the left hold the whole temporal structure of this model.",
              " Placental sex hormones rise exponentially towards term, and their metabolites (E2-17G and",
              " progesterone sulphate) inhibit BSEP. That is why the same genotype is",
              " silent at week 20, crosses the diagnostic threshold in the third trimester, and disappears",
              " within days of delivery — without a parameter being changed at either boundary.")
        ),

        ## ---- 2 -------------------------------------------------------
        tabPanel(
          "2. The bile acid assay vs reality",
          br(),
          plotOutput("p_tba", height = 340),
          hr(),
          plotOutput("p_species", height = 300),
          br(),
          div(style = "background:#fff4e5;padding:10px;border-radius:6px;",
              strong("Why this tab exists."),
              " The total bile acid assay used clinically is a 3α-hydroxysteroid dehydrogenase",
              " method, so it also counts UDCA and hyocholic acid. The gap that opens between the solid line",
              "(the assay value) and the dashed line (the endogenous fraction) in the graph above is precisely the drug effect",
              " that nobody attributes to the drug. The lower graph is the species composition, showing how",
              " lithocholic acid (LCA, cytotoxicity weight 1.0) changes when UDCA is given",
              " — it is a 7β-dehydroxylation product of the gut flora.")
        ),

        ## ---- 3 -------------------------------------------------------
        tabPanel(
          "3. Liver function · drug PK",
          br(),
          fluidRow(column(6, plotOutput("p_alt", height = 300)),
                   column(6, plotOutput("p_pk", height = 300))),
          hr(),
          fluidRow(column(6, plotOutput("p_transport", height = 300)),
                   column(6, plotOutput("p_cyp", height = 300))),
          br(),
          div(style = "background:#f7f7fb;padding:10px;border-radius:6px;",
              strong("Note."),
              " The CYP3A4 induction at the lower right is rifampicin's therapeutic mechanism (6α-hydroxylation →",
              " hyocholic acid, SULT2A1 induction → renal excretion) while at the same time co-inducing UGT1A1",
              " and so raising E2-17G, a BSEP inhibitor. The model does not hide this conflict,",
              " it computes it as it is.")
        ),

        ## ---- 4 -------------------------------------------------------
        tabPanel(
          "4. Fetal exposure and myocardium",
          br(),
          fluidRow(column(6, plotOutput("p_fetal", height = 300)),
                   column(6, plotOutput("p_wbar", height = 300))),
          hr(),
          fluidRow(column(6, plotOutput("p_heart", height = 300)),
                   column(6, plotOutput("p_arri", height = 300))),
          br(),
          div(style = "background:#ffecec;padding:10px;border-radius:6px;",
              strong("Where the threshold actually is."),
              " Across three stratified pregnancies maternal total bile acids differ 8.7-fold while the arrhythmia",
              " index differs 37-fold. Writing stillbirth risk on maternal bile acids needs an exponent of 2.95,",
              " on fetal hydrophobic load 2.55, and on the arrhythmia index only 1.55 is needed.",
              " In other words the steepness of the 100 µmol/L threshold lies not in placental transport or in the dose-response",
              " but in the threshold character of connexin-43 uncoupling. The map of this model was drawn",
              " expecting the opposite, and the fit corrected it.")
        ),

        ## ---- 5 -------------------------------------------------------
        tabPanel(
          "5. Placental transport saturation",
          br(),
          plotOutput("p_placenta", height = 340),
          hr(),
          fluidRow(column(6, plotOutput("p_ratio", height = 300)),
                   column(6, plotOutput("p_hyp", height = 300))),
          br(),
          div(style = "background:#f7f7fb;padding:10px;border-radius:6px;",
              strong("What placental transport explains and what it does not."),
              " When diffusive maternal→fetal influx and fetal synthesis itself exhaust the",
              " transporter capacity V_P, the cord:maternal ratio falls and then rises again (lower left).",
              " That sets the absolute level of fetal exposure. The ablation experiment",
              "(section I of icp_calibration.py), however, shows that removing this saturation does not remove the curvature",
              " of the risk curve — the curvature is in the myocardium of tab 4.")
        ),

        ## ---- 6 -------------------------------------------------------
        tabPanel(
          "6. The pruritus axis (a separate system)",
          br(),
          fluidRow(column(6, plotOutput("p_vas", height = 300)),
                   column(6, plotOutput("p_atx", height = 300))),
          hr(),
          plotOutput("p_axes", height = 320),
          br(),
          div(style = "background:#ecf0ff;padding:10px;border-radius:6px;",
              strong("Why UDCA barely works on pruritus."),
              " Autotaxin is driven mainly by progesterone sulphate and oestradiol and responds",
              " to bile acids only weakly, and saturably. This structure exists in order to satisfy four",
              " observations at once — circulating autotaxin does not fall with UDCA",
              " but does with rifampicin, the degree of pruritus correlates poorly with total bile acids,",
              " in a substantial fraction of patients pruritus precedes the biochemical",
              " abnormality by weeks, and asymptomatic hypercholanaemia exists. A model that drives pruritus from bile acids",
              " reproduces none of these.")
        ),

        ## ---- 7 -------------------------------------------------------
        tabPanel(
          "7. Clinical endpoints · the PITCHES decomposition",
          br(),
          fluidRow(column(6, plotOutput("p_haz", height = 300)),
                   column(6, plotOutput("p_comp", height = 300))),
          hr(),
          h4("Composition of the composite endpoint (at delivery)"),
          DTOutput("t_endpoints"),
          br(),
          div(style = "background:#fff4e5;padding:10px;border-radius:6px;",
              strong("The arithmetic of why PITCHES was negative."),
              " Two of the three components of the composite (delivery before 37 weeks, neonatal intensive care",
              " admission) are functions of gestational age at delivery and of the decision to deliver, not functions of bile acids.",
              " Stillbirth accounts for less than 1% of the placebo-arm composite, so abolishing stillbirth entirely",
              " moves the composite by only 0.2 percentage points. With 604 women it could not",
              " detect its own mechanism, and so this result does not contradict the stillbirth",
              " signal from the individual-patient-data meta-analysis — the two studies are measuring different things.")
        ),

        ## ---- 8 -------------------------------------------------------
        tabPanel(
          "8. Optimisation of delivery timing",
          br(),
          sliderInput("wmorb", "Neonatal morbidity : stillbirth utility ratio (WMORB)",
                      0.005, 0.15, 0.055, step = 0.005, width = "60%"),
          plotOutput("p_deliv", height = 380),
          hr(),
          DTOutput("t_deliv"),
          br(),
          div(style = "background:#ffecec;padding:10px;border-radius:6px;",
              strong("This tab agrees with the guidelines only in part."),
              " In the ≥100 µmol/L band the model computes 36–37 weeks as optimal, pointing the same way as",
              " the guidelines, and with treatment the optimum is pushed to 39 weeks. In the",
              " 40–99 band, however, it gives 39–40 weeks and does not reproduce the guidelines'",
              " 37–38 weeks. The reason is arithmetic: the excess stillbirth risk in that band is 0.28% against",
              " 0.13%, that is 0.15 percentage points, and under no reasonable weighting does it exceed the neonatal",
              " cost of delivering 2–3 weeks early. Take the slider above below 0.01 and",
              " 37 weeks becomes optimal, but then the <40 band is pulled forward with it and the conclusion is one",
              " no guideline recommends. The likeliest thing the model has missed is",
              " within-individual variability of bile acids — someone measured at 70 µmol/L can spend days above",
              " 100 between visits, and this model runs a smooth trajectory.")
        ),

        ## ---- 9 -------------------------------------------------------
        tabPanel(
          "9. Scenario comparison",
          br(),
          plotOutput("p_cmp", height = 480),
          hr(),
          DTOutput("t_cmp"),
          br(),
          div(style = "background:#f7f7fb;padding:10px;border-radius:6px;",
              strong("The comparison to try."),
              " Set the main setting to 'rifampicin 300 mg twice daily' and the comparator to 'UDCA 1000 mg',",
              " then look at the TBA panel and the VAS panel side by side. The ranking of the two drugs is inverted between the two",
              " panels. That is the whole of what this model has to say.")
        ),

        ## ---- 10 ------------------------------------------------------
        tabPanel(
          "10. Biomarkers · safety",
          br(),
          fluidRow(column(6, plotOutput("p_vitk", height = 300)),
                   column(6, plotOutput("p_mec", height = 300))),
          hr(),
          fluidRow(column(6, plotOutput("p_uterus", height = 300)),
                   column(6, plotOutput("p_surf", height = 300))),
          br(),
          div(style = "background:#eefaf6;padding:10px;border-radius:6px;",
              strong("A trade-off with a price on it."),
              " Cholestyramine lowers bile acids and at the same time lowers the intraluminal bile acid concentration,",
              " reducing micelle formation and so worsening vitamin K status (INR, upper left).",
              " Rifampicin does the same thing by shrinking the pool itself. Both drugs trade postpartum haemorrhage",
              " risk against bile acids, and the model computes that exchange rate. UDCA and",
              " the IBAT inhibitor do not pay this price.")
        )
      )
    )
  )
)

## =====================================================================
##  SERVER
## =====================================================================
server <- function(input, output, session) {

  gen_main <- reactive(GENOTYPES[[input$gen]])

  ev_main <- reactive({
    build_events(input$udca_mg, input$udca_from, input$rif_mg, input$rif_from,
                 input$chol_g, input$chol_from, input$same_mg, input$ibat,
                 input$ntx_mg, input$ah_mg,
                 if (input$bet_ga > 0) input$bet_ga else NA, input$del_ga)
  })

  sim <- reactive({
    extra <- if (input$vitk) list(VKDOSE = 1.6, VKSTART = ga2d(34)) else list()
    run_case(gen_main(), ev_main(), input$del_ga,
             twin = if (input$twin) 1.55 else 1, extra = extra)
  })

  sim_cmp <- reactive({
    e <- build_events(input$cmp_udca, 30, input$cmp_rif, 32, 0, 32, 0, 0, 0, 0,
                      NA, input$cmp_del)
    run_case(GENOTYPES[[input$cmp_gen]], e, input$cmp_del,
             twin = if (input$twin) 1.55 else 1)
  })

  at_delivery <- function(d, del_ga) {
    d %>% filter(GA <= del_ga + 1e-6) %>% slice_tail(n = 1)
  }

  vline_del <- function(del_ga)
    geom_vline(xintercept = del_ga, linetype = "dotted", colour = "grey40")

  ## ---- tab 1 --------------------------------------------------------
  output$p_steroid <- renderPlot({
    d <- sim()
    d %>% select(GA, `Oestradiol (nmol/L)` = E2,
                 `Progesterone sulphate (nmol/L)` = P4S,
                 `E2-17G (nmol/L ×10)` = E2G) %>%
      mutate(`E2-17G (nmol/L ×10)` = `E2-17G (nmol/L ×10)` * 10) %>%
      pivot_longer(-GA) %>%
      ggplot(aes(GA, value, colour = name)) +
      geom_line(linewidth = 0.9) + vline_del(input$del_ga) +
      labs(title = "Placental sex hormones and their BSEP-inhibiting metabolites",
           subtitle = "Exponential rise to term, rapid clearance after delivery (dotted = delivery)",
           x = "Gestational week", y = "Concentration", colour = NULL) + THEME
  })

  output$p_inhib <- renderPlot({
    d <- sim()
    d %>% mutate(`Total bile acid assay` = TBA,
                 `Diagnostic threshold 10` = 10, `Severe 40` = 40, `Very severe 100` = 100) %>%
      select(GA, `Total bile acid assay`, `Diagnostic threshold 10`, `Severe 40`, `Very severe 100`) %>%
      pivot_longer(-GA) %>%
      ggplot(aes(GA, value, colour = name,
                 linetype = name == "Total bile acid assay")) +
      geom_line(linewidth = 0.9) + vline_del(input$del_ga) +
      scale_linetype_manual(values = c("TRUE" = "solid", "FALSE" = "dashed"),
                            guide = "none") +
      labs(title = "When the same genotype crosses the diagnostic criterion",
           subtitle = "Both the timing of onset and the postpartum resolution emerge without any parameter change",
           x = "Gestational week", y = "µmol/L", colour = NULL) + THEME
  })

  output$t_summary <- renderDT({
    d <- at_delivery(sim(), input$del_ga)
    post <- sim() %>% filter(GA > input$del_ga, TBA < 10) %>% slice_head(n = 1)
    onset <- sim() %>% filter(GA <= input$del_ga, TBA >= 10) %>% slice_head(n = 1)
    tibble(
      Item = c("Total bile acid assay at delivery (µmol/L)",
               "Endogenous bile acids at delivery (µmol/L)",
               "Plasma UDCA (µmol/L)", "Hyocholic acid (µmol/L)",
               "Cord blood total bile acids (µmol/L)", "Fetal hydrophobic load FCL",
               "Fetal pool hydrophobicity w̄", "Arrhythmia index ARRI",
               "ALT (U/L)", "Pruritus VAS (0-10)", "Autotaxin (relative)", "INR",
               "Week the diagnostic threshold (10) is crossed", "Falls below 10 µmol/L after delivery"),
      Value = c(sprintf("%.1f", d$TBA), sprintf("%.1f", d$TBA_ENDO),
             sprintf("%.1f", d$UDCA_P), sprintf("%.1f", d$HCA_P),
             sprintf("%.1f", d$CORD), sprintf("%.2f", d$FCLo),
             sprintf("%.3f", d$WBAR), sprintf("%.4f", d$ARRIo),
             sprintf("%.0f", d$ALT), sprintf("%.1f", d$VAS),
             sprintf("%.2f", d$ATX), sprintf("%.2f", d$INR),
             if (nrow(onset)) sprintf("%.1f wk", onset$GA) else "not crossed",
             if (nrow(post)) sprintf("%.0f d", (post$GA - input$del_ga) * 7)
             else "not reached within 21 days")
    )
  }, options = list(dom = "t", pageLength = 20), rownames = FALSE)

  ## ---- tab 2 --------------------------------------------------------
  output$p_tba <- renderPlot({
    d <- sim()
    d %>% select(GA, `Assay value (includes UDCA · HCA)` = TBA,
                 `Endogenous fraction` = TBA_ENDO, `Plasma UDCA` = UDCA_P,
                 `Hyocholic acid` = HCA_P) %>%
      pivot_longer(-GA) %>%
      ggplot(aes(GA, value, colour = name, linetype = name)) +
      geom_line(linewidth = 1) +
      geom_hline(yintercept = c(10, 40, 100), linetype = "dotted",
                 colour = "grey60") +
      vline_del(input$del_ga) +
      scale_linetype_manual(values = c("Assay value (includes UDCA · HCA)" = "solid",
                                       "Endogenous fraction" = "longdash",
                                       "Plasma UDCA" = "dotted",
                                       "Hyocholic acid" = "dotdash")) +
      labs(title = "What the clinic measures and what the clinic believes it is measuring",
           subtitle = "The gap between the two lines = the drug effect nobody attributes to the drug",
           x = "Gestational week", y = "µmol/L", colour = NULL, linetype = NULL) + THEME
  })

  output$p_species <- renderPlot({
    d <- sim()
    d %>% mutate(CA = SPCA / 18, CDCA = SPCD / 18, DCA = SPDC / 18,
                 LCA = SPLC / 18, UDCA = SPUD / 18) %>%
      select(GA, CA, CDCA, DCA, LCA, UDCA) %>%
      pivot_longer(-GA) %>%
      mutate(name = factor(name, levels = c("CA", "CDCA", "DCA", "LCA", "UDCA"))) %>%
      ggplot(aes(GA, value, fill = name)) +
      geom_area(position = "stack", alpha = 0.85) + vline_del(input$del_ga) +
      scale_fill_brewer(palette = "RdYlBu", direction = -1) +
      labs(title = "Species composition of maternal circulating bile acids",
           subtitle = "Cytotoxicity weights CA 0.20 · CDCA 0.60 · DCA 0.72 · LCA 1.00 · UDCA 0.02",
           x = "Gestational week", y = "µmol/L", fill = NULL) + THEME
  })

  ## ---- tab 3 --------------------------------------------------------
  output$p_alt <- renderPlot({
    d <- sim()
    d %>% select(GA, `ALT (U/L)` = ALT,
                 `Oxidative stress ×100` = ROS) %>%
      mutate(`Oxidative stress ×100` = `Oxidative stress ×100` * 100) %>%
      pivot_longer(-GA) %>%
      ggplot(aes(GA, value, colour = name)) + geom_line(linewidth = 0.9) +
      vline_del(input$del_ga) +
      labs(title = "Hepatocellular injury", subtitle = "A rising biliary bile acid:PC ratio is the mechanism of the ABCB4 variants",
           x = "Gestational week", y = NULL, colour = NULL) + THEME
  })

  output$p_pk <- renderPlot({
    d <- sim()
    d %>% mutate(`Rifampicin (mg/L)` = RIFC / 45,
                 `Plasma UDCA (µmol/L)` = UDCA_P,
                 `SAMe (mg/L)` = SAMC / 22) %>%
      select(GA, `Rifampicin (mg/L)`, `Plasma UDCA (µmol/L)`, `SAMe (mg/L)`) %>%
      pivot_longer(-GA) %>%
      ggplot(aes(GA, value, colour = name)) + geom_line(linewidth = 0.8) +
      vline_del(input$del_ga) +
      labs(title = "Drug exposure", subtitle = "Dosing stops at delivery",
           x = "Gestational week", y = NULL, colour = NULL) + THEME
  })

  output$p_transport <- renderPlot({
    d <- sim()
    d %>% select(GA, `Basolateral efflux MRP3/4/OSTαβ` = MRP4,
                 `SULT2A1 induction` = SULT, `CYP7A1` = CYP7A1,
                 `FGF19` = FGF19) %>%
      pivot_longer(-GA) %>%
      ggplot(aes(GA, value, colour = name)) + geom_line(linewidth = 0.9) +
      vline_del(input$del_ga) +
      labs(title = "Adaptive response", subtitle = "The escape valves are quiet at baseline and strongly induced in cholestasis",
           x = "Gestational week", y = "Relative", colour = NULL) + THEME
  })

  output$p_cyp <- renderPlot({
    d <- sim()
    d %>% select(GA, `CYP3A4 induction` = CYP3A, `E2-17G (nmol/L)` = E2G) %>%
      pivot_longer(-GA) %>%
      ggplot(aes(GA, value, colour = name)) + geom_line(linewidth = 0.9) +
      vline_del(input$del_ga) +
      labs(title = "The two conflicting effects of rifampicin",
           subtitle = "CYP3A4 induction (therapeutic) and UGT1A1 co-induction → rising E2-17G (counterproductive)",
           x = "Gestational week", y = NULL, colour = NULL) + THEME
  })

  ## ---- tab 4 --------------------------------------------------------
  output$p_fetal <- renderPlot({
    d <- sim()
    d %>% select(GA, `Maternal assay value` = TBA, `Cord blood total bile acids` = CORD,
                 `Fetal hydrophobic load FCL` = FCLo) %>%
      pivot_longer(-GA) %>%
      ggplot(aes(GA, value, colour = name)) + geom_line(linewidth = 1) +
      vline_del(input$del_ga) +
      labs(title = "Maternal · fetal exposure", subtitle = "The fetal load is a hydrophobicity-weighted sum, not a total",
           x = "Gestational week", y = "µmol/L", colour = NULL) + THEME
  })

  output$p_wbar <- renderPlot({
    d <- sim() %>% filter(GA <= input$del_ga)
    ggplot(d, aes(GA, WBAR)) + geom_line(linewidth = 1, colour = "#6a1b9a") +
      labs(title = "Mean hydrophobicity of the fetal pool w̄",
           subtitle = "The hydrophobic species have higher diffusive permeability, so the pool does not merely grow, it gets worse",
           x = "Gestational week", y = "FCL / cord blood total bile acids") + THEME
  })

  output$p_heart <- renderPlot({
    d <- sim() %>% filter(GA <= input$del_ga)
    d %>% select(GA, `Connexin-43 coupling GJ` = GJ, `Calcium overload index` = FCA) %>%
      pivot_longer(-GA) %>%
      ggplot(aes(GA, value, colour = name)) + geom_line(linewidth = 1) +
      labs(title = "Two states of the fetal myocardium", subtitle = "GJ is threshold-like (Hill 1.6), Ca is a gently saturating form",
           x = "Gestational week", y = "Relative (normal = 1)", colour = NULL) + THEME
  })

  output$p_arri <- renderPlot({
    d <- sim() %>% filter(GA <= input$del_ga)
    ggplot(d, aes(GA, ARRIo)) +
      geom_line(linewidth = 1.1, colour = "#b71c1c") +
      labs(title = "Arrhythmia index ARRI = (1 − GJ) × Ca",
           subtitle = "Stillbirth hazard h(t) = HSB0 + HSBSC · ARRI^1.55 · (1 + 0.55·hypoxia)",
           x = "Gestational week", y = "ARRI") + THEME
  })

  ## ---- tab 5 --------------------------------------------------------
  output$p_placenta <- renderPlot({
    d <- sim() %>% filter(GA <= input$del_ga)
    d %>% mutate(fCA = FPCA, fCD = FPCD, fDC = FPDC, fLC = FPLC) %>%
      select(GA, `Fetal CA` = fCA, `Fetal CDCA` = fCD, `Fetal DCA` = fDC,
             `Fetal LCA` = fLC) %>%
      pivot_longer(-GA) %>%
      ggplot(aes(GA, value, colour = name)) + geom_line(linewidth = 0.9) +
      labs(title = "Bile acid species in the fetal compartment (µmol)",
           subtitle = "Diffusion is bidirectional, so the fetal concentration saturates slightly above the maternal one",
           x = "Gestational week", y = "µmol", colour = NULL) + THEME
  })

  output$p_ratio <- renderPlot({
    d <- sim() %>% filter(GA <= input$del_ga, TBA > 0.1)
    ggplot(d, aes(TBA, CORD / TBA)) +
      geom_path(linewidth = 1, colour = "#d84315") +
      geom_vline(xintercept = 100, linetype = "dotted") +
      labs(title = "Cord : maternal ratio", subtitle = "The point where it stops falling and rises again is transporter saturation",
           x = "Maternal total bile acids (µmol/L)", y = "Cord / maternal") + THEME
  })

  output$p_hyp <- renderPlot({
    d <- sim() %>% filter(GA <= input$del_ga)
    d %>% select(GA, `Chorionic plate vasoconstriction` = VASO, `Trophoblast injury` = TROPH,
                 `Fetal hypoxia index` = HYP) %>%
      pivot_longer(-GA) %>%
      ggplot(aes(GA, value, colour = name)) + geom_line(linewidth = 0.9) +
      labs(title = "Placental vasculature and hypoxia", x = "Gestational week", y = "Relative",
           colour = NULL) + THEME
  })

  ## ---- tab 6 --------------------------------------------------------
  output$p_vas <- renderPlot({
    d <- sim()
    ggplot(d, aes(GA, VAS)) + geom_line(linewidth = 1.1, colour = "#283593") +
      vline_del(input$del_ga) + ylim(0, 10) +
      labs(title = "Pruritus VAS (0-10 cm)",
           subtitle = "Normal pregnancy ~1.3 · mild ICP ~4 · severe ~5.6 · very severe ~7",
           x = "Gestational week", y = "VAS") + THEME
  })

  output$p_atx <- renderPlot({
    d <- sim()
    d %>% select(GA, `Autotaxin ATX` = ATX, `LPA` = LPA,
                 `Central itch state` = ITCHC) %>%
      pivot_longer(-GA) %>%
      ggplot(aes(GA, value, colour = name)) + geom_line(linewidth = 0.9) +
      vline_del(input$del_ga) +
      labs(title = "Inside the pruritus axis", subtitle = "Only rifampicin lowers ATX (PXR-dependent)",
           x = "Gestational week", y = "Relative", colour = NULL) + THEME
  })

  output$p_axes <- renderPlot({
    d <- sim() %>% filter(GA <= input$del_ga)
    rng <- function(x) if (diff(range(x)) < 1e-9) rep(0.5, length(x)) else
      (x - min(x)) / diff(range(x))
    d %>% mutate(`Bile acid axis (normalised)` = rng(TBA_ENDO),
                 `Fetal load FCL (normalised)` = rng(FCLo),
                 `Pruritus axis VAS (normalised)` = rng(VAS)) %>%
      select(GA, `Bile acid axis (normalised)`, `Fetal load FCL (normalised)`,
             `Pruritus axis VAS (normalised)`) %>%
      pivot_longer(-GA) %>%
      ggplot(aes(GA, value, colour = name)) + geom_line(linewidth = 1.1) +
      labs(title = "The two axes overlaid on the same scale",
           subtitle = "When you change the prescription, watch whether the two lines move together or separately",
           x = "Gestational week", y = "Normalised to each one's own range", colour = NULL) + THEME
  })

  ## ---- tab 7 --------------------------------------------------------
  output$p_haz <- renderPlot({
    d <- sim() %>% filter(GA <= input$del_ga)
    d %>% select(GA, `Cumulative stillbirth risk %` = SBRISK,
                 `Cumulative spontaneous preterm birth risk %` = PTBRISK,
                 `Cumulative meconium risk %` = MECRISK) %>%
      mutate(across(-GA, ~ .x * 100)) %>%
      pivot_longer(-GA) %>%
      ggplot(aes(GA, value, colour = name)) + geom_line(linewidth = 1) +
      labs(title = "Cumulative risk", x = "Gestational week", y = "%", colour = NULL) + THEME
  })

  output$p_comp <- renderPlot({
    d <- at_delivery(sim(), input$del_ga)
    ptb_any <- if (input$del_ga < 37) 1 else d$PTBRISK
    tibble(Component = c("Stillbirth", "Delivery before 37 weeks", "Neonatal intensive care ≥4h"),
           Probability = c(d$SBRISK, ptb_any, d$NICURISK) * 100) %>%
      ggplot(aes(reorder(Component, Probability), Probability, fill = Component)) +
      geom_col(width = 0.6, show.legend = FALSE) +
      geom_text(aes(label = sprintf("%.2f%%", Probability)), hjust = -0.1, size = 4) +
      coord_flip(clip = "off") +
      scale_fill_manual(values = c("Stillbirth" = "#b71c1c",
                                   "Delivery before 37 weeks" = "#f9a825",
                                   "Neonatal intensive care ≥4h" = "#1565c0")) +
      labs(title = "Composition of the composite endpoint",
           subtitle = "Stillbirth is under 1% of the composite — abolishing it barely moves the composite",
           x = NULL, y = "%") + THEME +
      theme(plot.margin = margin(5, 40, 5, 5))
  })

  output$t_endpoints <- renderDT({
    d <- at_delivery(sim(), input$del_ga)
    ptb_any <- if (input$del_ga < 37) 1 else d$PTBRISK
    comp <- 1 - (1 - d$SBRISK) * (1 - ptb_any) * (1 - d$NICURISK)
    tibble(
      Endpoint = c("Stillbirth (cumulative)", "Spontaneous preterm birth (cumulative)", "Iatrogenic preterm birth",
                     "Meconium-stained amniotic fluid", "Neonatal respiratory distress", "Neonatal intensive care ≥4h",
                     "Composite endpoint (PITCHES definition)",
                     "Share of the composite taken by stillbirth"),
      `Probability %` = sprintf("%.3f", 100 * c(d$SBRISK, d$PTBRISK,
                                          if (input$del_ga < 37) 1 else 0,
                                          d$MECRISK, d$RDSRISK, d$NICURISK,
                                          comp, d$SBRISK / comp))
    )
  }, options = list(dom = "t", pageLength = 10), rownames = FALSE)

  ## ---- tab 8 --------------------------------------------------------
  deliv_curve <- reactive({
    gen <- gen_main()
    gas <- seq(35, 40, by = 0.5)
    res <- lapply(gas, function(g) {
      e <- build_events(input$udca_mg, input$udca_from, input$rif_mg,
                        input$rif_from, input$chol_g, input$chol_from,
                        input$same_mg, input$ibat, input$ntx_mg, input$ah_mg,
                        if (input$bet_ga > 0) input$bet_ga else NA, g)
      d <- run_case(gen, e, g, twin = if (input$twin) 1.55 else 1)
      a <- at_delivery(d, g)
      tibble(GA = g, Stillbirth = a$SBRISK,
             Morbidity = a$NICURISK + a$RDSRISK)
    })
    bind_rows(res) %>% mutate(TotalLoss = Stillbirth + input$wmorb * Morbidity)
  })

  output$p_deliv <- renderPlot({
    d <- deliv_curve()
    opt <- d$GA[which.min(d$TotalLoss)]
    d %>% mutate(`Stillbirth risk` = Stillbirth * 1000,
                 `Neonatal morbidity cost (weighted)` = input$wmorb * Morbidity * 1000,
                 `Total` = TotalLoss * 1000) %>%
      select(GA, `Stillbirth risk`, `Neonatal morbidity cost (weighted)`, `Total`) %>%
      pivot_longer(-GA) %>%
      ggplot(aes(GA, value, colour = name, linewidth = name == "Total")) +
      geom_line() +
      geom_vline(xintercept = opt, linetype = "dashed", colour = "#b71c1c") +
      annotate("text", x = opt, y = Inf, vjust = 1.5, hjust = -0.1,
               label = sprintf("optimum %.1f wk", opt), colour = "#b71c1c") +
      scale_linewidth_manual(values = c("TRUE" = 1.4, "FALSE" = 0.9),
                             guide = "none") +
      labs(title = "Delivery timing is not a threshold but the crossing of two curves",
           subtitle = "The further left you go the smaller stillbirth becomes and the more steeply the neonatal cost rises",
           x = "Gestational week at delivery", y = "milli-expected-loss", colour = NULL) + THEME
  })

  output$t_deliv <- renderDT({
    deliv_curve() %>%
      transmute(`Delivery week` = GA,
                `Stillbirth %` = sprintf("%.3f", 100 * Stillbirth),
                `Morbidity %` = sprintf("%.1f", 100 * Morbidity),
                `Total loss (×1000)` = sprintf("%.2f", 1000 * TotalLoss))
  }, options = list(dom = "t", pageLength = 12), rownames = FALSE)

  ## ---- tab 9 --------------------------------------------------------
  output$p_cmp <- renderPlot({
    a <- sim() %>% mutate(arm = "Main setting")
    b <- sim_cmp() %>% mutate(arm = "Comparator")
    bind_rows(a, b) %>%
      select(GA, arm, `Total bile acid assay` = TBA, `Endogenous bile acids` = TBA_ENDO,
             `Fetal hydrophobic load FCL` = FCLo, `Pruritus VAS` = VAS,
             `ALT` = ALT, `INR` = INR) %>%
      pivot_longer(c(-GA, -arm)) %>%
      mutate(name = factor(name, levels = c("Total bile acid assay", "Endogenous bile acids",
                                            "Fetal hydrophobic load FCL",
                                            "Pruritus VAS", "ALT", "INR"))) %>%
      ggplot(aes(GA, value, colour = arm)) + geom_line(linewidth = 1) +
      facet_wrap(~ name, scales = "free_y", ncol = 3) +
      scale_colour_manual(values = c("Main setting" = "#1565c0",
                                     "Comparator" = "#c62828")) +
      labs(title = "Comparison of the two scenarios",
           subtitle = "Check whether the ranking is inverted between the bile acid panel and the pruritus panel",
           x = "Gestational week", y = NULL, colour = NULL) + THEME
  })

  output$t_cmp <- renderDT({
    a <- at_delivery(sim(), input$del_ga)
    b <- at_delivery(sim_cmp(), input$cmp_del)
    f <- function(x) sprintf("%.2f", x)
    tibble(
      Item = c("Delivery week", "Total bile acid assay", "Endogenous bile acids",
               "Cord blood total bile acids", "Fetal hydrophobic load FCL", "Arrhythmia index ARRI",
               "Pruritus VAS", "Autotaxin", "ALT", "INR",
               "Cumulative stillbirth %", "Neonatal intensive care %"),
      `Main setting` = c(f(input$del_ga), f(a$TBA), f(a$TBA_ENDO), f(a$CORD),
                    f(a$FCLo), sprintf("%.4f", a$ARRIo), f(a$VAS), f(a$ATX),
                    sprintf("%.0f", a$ALT), f(a$INR),
                    sprintf("%.3f", 100 * a$SBRISK),
                    sprintf("%.1f", 100 * a$NICURISK)),
      `Comparator` = c(f(input$cmp_del), f(b$TBA), f(b$TBA_ENDO), f(b$CORD),
                   f(b$FCLo), sprintf("%.4f", b$ARRIo), f(b$VAS), f(b$ATX),
                   sprintf("%.0f", b$ALT), f(b$INR),
                   sprintf("%.3f", 100 * b$SBRISK),
                   sprintf("%.1f", 100 * b$NICURISK))
    )
  }, options = list(dom = "t", pageLength = 15), rownames = FALSE)

  ## ---- tab 10 -------------------------------------------------------
  output$p_vitk <- renderPlot({
    d <- sim()
    d %>% select(GA, `Vitamin K status` = VITK, `INR` = INR,
                 `Micelle formation surrogate (PIVKA)` = PIVKA) %>%
      pivot_longer(-GA) %>%
      ggplot(aes(GA, value, colour = name)) + geom_line(linewidth = 0.9) +
      vline_del(input$del_ga) +
      labs(title = "Fat-soluble vitamins and coagulation",
           subtitle = "Cholestyramine and rifampicin trade bile acids against vitamin K",
           x = "Gestational week", y = NULL, colour = NULL) + THEME
  })

  output$p_mec <- renderPlot({
    d <- sim() %>% filter(GA <= input$del_ga)
    d %>% select(GA, `Meconium bile acid accumulation (µmol/100)` = MEC,
                 `Cumulative meconium risk %` = MECRISK) %>%
      mutate(`Meconium bile acid accumulation (µmol/100)` = `Meconium bile acid accumulation (µmol/100)` / 100,
             `Cumulative meconium risk %` = `Cumulative meconium risk %` * 100) %>%
      pivot_longer(-GA) %>%
      ggplot(aes(GA, value, colour = name)) + geom_line(linewidth = 0.9) +
      labs(title = "Meconium", subtitle = "The risk accumulates only near term (maturation of gut motility)",
           x = "Gestational week", y = NULL, colour = NULL) + THEME
  })

  output$p_uterus <- renderPlot({
    d <- sim() %>% filter(GA <= input$del_ga)
    d %>% select(GA, `Oxytocin receptor density` = OTR, `Uterine activity index` = UTA) %>%
      pivot_longer(-GA) %>%
      ggplot(aes(GA, value, colour = name)) + geom_line(linewidth = 0.9) +
      labs(title = "Uterine muscle", subtitle = "Bile acids raise OTR and so create the spontaneous preterm birth risk",
           x = "Gestational week", y = "Relative", colour = NULL) + THEME
  })

  output$p_surf <- renderPlot({
    gas <- seq(34, 41, by = 0.25)
    d <- at_delivery(sim(), input$del_ga)
    tibble(GA = gas,
           `Respiratory distress (no steroid)` =
             (0.0028 + 1 / (1 + exp(1.15 * (gas - 33.6)))) * 100,
           `Respiratory distress (current SURF)` =
             (0.0028 + 1 / (1 + exp(1.15 * (gas - 33.6)))) /
             (1 + 1.55 * max(0, d$SURF - 1)) * 100) %>%
      pivot_longer(-GA) %>%
      ggplot(aes(GA, value, colour = name)) + geom_line(linewidth = 1) +
      geom_vline(xintercept = input$del_ga, linetype = "dotted") +
      labs(title = "The cost of iatrogenic preterm birth",
           subtitle = "The steepness of this curve is half of the delivery timing optimisation",
           x = "Gestational week at delivery", y = "%", colour = NULL) + THEME
  })
}

shinyApp(ui, server)
