## =============================================================================
##  ntm_shiny_app.R — Interactive dashboard for the MAC-PD QSP model
##  Nontuberculous Mycobacterial Lung Disease (Mycobacterium avium complex)
##
##  Run with:   shiny::runApp("ntm_shiny_app.R")
##  Requires:   shiny, mrgsolve, dplyr, tidyr, ggplot2, DT, bslib
##              plus ntm_mrgsolve_model.R in the same directory
##
##  11 tabs:
##    1  Patient profile       Patient profile & host phenotype
##    2  Drug PK               Drug PK across plasma / ELF / intracellular
##    3  Bacterial niches      Where the bacteria actually are (4 niches)
##    4  Clinical endpoints    Culture conversion, QOL-B, BMI, cavity
##    5  Resistance evolution  rrl mutant selection and the per-niche gate
##    6  Scenario comparison   12 regimens side by side
##    7  Biomarkers            Immune / inflammatory / airway readouts
##    8  Toxicity · safety     Ototoxicity, optic, renal, hepatic, QTc
##    9  pH paradox explorer   The one input, two consequences, explored
##   10  Sensitivity analysis  Local sensitivity on the conversion endpoint
##   11  Model information     Structure, assumptions, references
## =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

## ---------------------------------------------------------------------------
## Load the model. ntm_mrgsolve_model.R defines `ntm_code` and compiles `mod`;
## we source it in a sandbox environment so the scenario simulations at the
## bottom of that file do not run every time the app starts.
## ---------------------------------------------------------------------------
build_model <- function() {
  src  <- readLines("ntm_mrgsolve_model.R")
  stop_at <- grep("^##  COMPARTMENT NUMBERS FOR DOSING", src)[1]
  if (is.na(stop_at)) stop_at <- length(src)
  e <- new.env()
  eval(parse(text = paste(src[1:(stop_at - 2)], collapse = "\n")), envir = e)
  e$mod
}
mod <- build_model()

CMT <- setNames(seq_along(mrgsolve::cmt(mod)), mrgsolve::cmt(mod))

THEME <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "#EEF2F7", colour = NA),
        legend.position  = "bottom")

PAL <- c("#0277BD", "#C62828", "#2E7D32", "#EF6C00",
         "#5E35B1", "#00838F", "#AD1457", "#455A64")

## ---------------------------------------------------------------------------
## Regimen builder — turns UI switches into an mrgsolve event object
## ---------------------------------------------------------------------------
build_events <- function(tx_days, macro, macro_freq, emb, rif, cfz,
                         amk_iv, amk_iv_wks, alis, wt) {
  evs <- list()
  add <- function(cmt, amt, ii, until, rate = 0) {
    if (amt <= 0 || until <= 0) return(NULL)
    ev(amt = amt, cmt = cmt, ii = ii, addl = max(0, floor(until / ii)),
       time = 0, rate = rate)
  }
  ii_m <- if (macro_freq == "TIW") 7/3 else if (macro_freq == "BID") 0.5 else 1
  if (macro != "none")
    evs <- c(evs, list(add(CMT["MGUT"],
                           if (macro == "CLR") 500 else if (ii_m > 1) 500 else 250,
                           ii_m, tx_days)))
  if (emb)  evs <- c(evs, list(add(CMT["EGUT"],
                                   if (ii_m > 1) 25 * wt else 15 * wt, ii_m, tx_days)))
  if (rif)  evs <- c(evs, list(add(CMT["FGUT"], 600, ii_m, tx_days)))
  if (cfz)  evs <- c(evs, list(add(CMT["CGUT"], 100, 1, tx_days)))
  if (amk_iv) evs <- c(evs, list(add(CMT["KCEN"], 15 * wt, 7/3,
                                     amk_iv_wks * 7, rate = 15 * wt * 24)))
  if (alis) evs <- c(evs, list(add(CMT["KLIP"], 590, 1, tx_days)))
  evs <- Filter(Negate(is.null), evs)
  if (!length(evs)) return(NULL)
  Reduce(c, evs)
}

run_sim <- function(pars, events, end = 540) {
  s <- mod %>% param(pars)
  out <- if (is.null(events)) mrgsim(s, end = end, delta = 1)
         else mrgsim(s, events = events, end = end, delta = 1)
  as_tibble(out)
}

## ---------------------------------------------------------------------------
## Pre-defined comparison scenarios (tab 6)
## ---------------------------------------------------------------------------
SCENARIOS <- list(
  "1. Watchful waiting" =
    list(p = list(CAVFLAG = 0), e = NULL),
  "2. AZM+EMB+RIF three times weekly (nodular-bronchiectatic)" =
    list(p = list(CAVFLAG = 0, MACTYPE = 0),
         e = list(tx = 365, macro = "AZM", freq = "TIW", emb = TRUE, rif = TRUE,
                  cfz = FALSE, iv = FALSE, wks = 0, alis = FALSE)),
  "3. AZM+EMB+RIF daily (fibrocavitary)" =
    list(p = list(CAVFLAG = 1, MACTYPE = 0),
         e = list(tx = 365, macro = "AZM", freq = "QD", emb = TRUE, rif = TRUE,
                  cfz = FALSE, iv = FALSE, wks = 0, alis = FALSE)),
  "4. + IV amikacin 3 months (fibrocavitary)" =
    list(p = list(CAVFLAG = 1),
         e = list(tx = 365, macro = "AZM", freq = "QD", emb = TRUE, rif = TRUE,
                  cfz = FALSE, iv = TRUE, wks = 13, alis = FALSE)),
  "5. + ALIS inhaled 590 mg/day (CONVERT)" =
    list(p = list(CAVFLAG = 1),
         e = list(tx = 365, macro = "AZM", freq = "QD", emb = TRUE, rif = TRUE,
                  cfz = FALSE, iv = FALSE, wks = 0, alis = TRUE)),
  "6. Azithromycin monotherapy (contraindicated)" =
    list(p = list(CAVFLAG = 0),
         e = list(tx = 365, macro = "AZM", freq = "QD", emb = FALSE, rif = FALSE,
                  cfz = FALSE, iv = FALSE, wks = 0, alis = FALSE)),
  "7. CLR+EMB+RIF (CYP3A interaction)" =
    list(p = list(CAVFLAG = 1, MACTYPE = 1),
         e = list(tx = 365, macro = "CLR", freq = "BID", emb = TRUE, rif = TRUE,
                  cfz = FALSE, iv = FALSE, wks = 0, alis = FALSE)),
  "8. AZM+EMB+CFZ+ALIS (rifampicin excluded)" =
    list(p = list(CAVFLAG = 1),
         e = list(tx = 365, macro = "AZM", freq = "QD", emb = TRUE, rif = FALSE,
                  cfz = TRUE, iv = FALSE, wks = 0, alis = TRUE)),
  "9. No. 8 + airway clearance + nutritional support" =
    list(p = list(CAVFLAG = 1, ACT = 1, NUTR = 1),
         e = list(tx = 365, macro = "AZM", freq = "QD", emb = TRUE, rif = FALSE,
                  cfz = TRUE, iv = FALSE, wks = 0, alis = TRUE)),
  "10. Anti-IFN-γ autoantibody host + ALIS" =
    list(p = list(CAVFLAG = 1, IFNCAP = 0.15),
         e = list(tx = 365, macro = "AZM", freq = "QD", emb = TRUE, rif = TRUE,
                  cfz = FALSE, iv = FALSE, wks = 0, alis = TRUE)),
  "11. Fibrocavitary AZM monotherapy (normal host)" =
    list(p = list(CAVFLAG = 1),
         e = list(tx = 365, macro = "AZM", freq = "QD", emb = FALSE, rif = FALSE,
                  cfz = FALSE, iv = FALSE, wks = 0, alis = FALSE)),
  "12. Fibrocavitary monotherapy + anti-IFN-γ (resistance sweep)" =
    list(p = list(CAVFLAG = 1, IFNCAP = 0.15),
         e = list(tx = 365, macro = "AZM", freq = "QD", emb = FALSE, rif = FALSE,
                  cfz = FALSE, iv = FALSE, wks = 0, alis = FALSE))
)

scenario_events <- function(s, wt = 52) {
  if (is.null(s$e)) return(NULL)
  with(s$e, build_events(tx, macro, freq, emb, rif, cfz, iv, wks, alis, wt))
}

## =============================================================================
##  UI
## =============================================================================
ui <- fluidPage(
  tags$head(tags$style(HTML("
    body { font-family: -apple-system, 'Segoe UI', 'Noto Sans KR', sans-serif; }
    .banner { background:linear-gradient(90deg,#0277BD,#5E35B1);
              color:white; padding:14px 18px; border-radius:8px;
              margin-bottom:14px; }
    .banner h3 { margin:0 0 4px 0; font-weight:700; }
    .banner p  { margin:0; font-size:12.5px; opacity:.93; }
    .keybox { background:#FFF8E1; border-left:5px solid #F9A825;
              padding:10px 14px; border-radius:5px; margin:10px 0;
              font-size:13px; }
    .well { background:#F7F9FB; }
  "))),

  div(class = "banner",
      h3("Nontuberculous mycobacterial lung disease (MAC-PD) QSP model dashboard"),
      p("Nontuberculous Mycobacterial Lung Disease · Mycobacterium avium complex — ",
        "47-ODE niche-resolved PK/PD model  |  ",
        "Core structure: the bacteria are split across four physical niches, and phagosomal pH 5.2 alone ",
        "produces both macrolide accumulation (151-fold) and loss of activity (12.6-fold) at once.")),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      h4("Patient phenotype"),
      radioButtons("cav", "Disease phenotype",
                   c("Nodular-bronchiectatic" = 0,
                     "Fibrocavitary" = 1),
                   selected = 1),
      sliderInput("wt", "Body weight (kg)", 38, 80, 52, 1),
      sliderInput("ifncap", "Th1/IFN-γ axis function (1 = normal, ↓ = anti-IFN-γ antibody · MSMD)",
                  0.05, 1, 1, 0.05),
      sliderInput("mcc0", "Baseline mucociliary clearance MCC₀", 0.20, 0.95, 0.75, 0.05),

      hr(), h4("Drug regimen"),
      selectInput("macro", "Macrolide",
                  c("Azithromycin (AZM)" = "AZM",
                    "Clarithromycin (CLR)" = "CLR",
                    "None" = "none"), selected = "AZM"),
      selectInput("freq", "Dosing frequency",
                  c("Daily (QD)" = "QD", "Three times weekly (TIW)" = "TIW",
                    "Twice daily (BID, CLR)" = "BID"), selected = "QD"),
      checkboxInput("emb", "Ethambutol (EMB)", TRUE),
      checkboxInput("rif", "Rifampicin (RIF)", TRUE),
      checkboxInput("cfz", "Clofazimine (CFZ) 100 mg/day", FALSE),
      checkboxInput("alis", "ALIS inhaled 590 mg/day (liposomal amikacin)", FALSE),
      checkboxInput("ivamk", "IV amikacin 15 mg/kg three times weekly", FALSE),
      conditionalPanel("input.ivamk == true",
                       sliderInput("ivwks", "IV amikacin duration (weeks)", 2, 26, 12, 1)),
      sliderInput("tx", "Total treatment duration (days)", 90, 540, 365, 30),

      hr(), h4("Non-pharmacological interventions"),
      checkboxInput("act", "Airway clearance therapy (HFCWO + hypertonic saline)", FALSE),
      checkboxInput("nutr", "Nutritional support", FALSE),

      hr(),
      sliderInput("phag", "Phagosomal pH (the structural key parameter)", 4.5, 7.4, 5.2, 0.1),
      helpText("This one slider changes both the intracellular accumulation ratio of the macrolide ",
               "and the MIC rise at the same time. See the result directly in tab 9."),

      hr(),
      actionButton("go", "Run simulation", class = "btn-primary btn-block")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs", type = "tabs",

        ## ---- TAB 1 ---------------------------------------------------------
        tabPanel("1. Patient profile",
          br(),
          div(class = "keybox",
              strong("Why the host phenotype matters — "),
              "MAC-PD is defined less by the virulence of the pathogen than by the airway clearance and body habitus of the host. ",
              "A slender habitus, thoracic dysmorphism, CFTR heterozygosity and ciliary functional variants lower MCC, ",
              "and a lowered MCC becomes the entry point of the Cole vicious circle that closes as bacterial stasis → ",
              "biofilm → airway damage → MCC falling again."),
          fluidRow(
            column(6, plotOutput("p_profile", height = 320)),
            column(6, plotOutput("p_mcc", height = 320))),
          hr(),
          h4("Settings summary"), DTOutput("t_profile")),

        ## ---- TAB 2 ---------------------------------------------------------
        tabPanel("2. Drug PK",
          br(),
          div(class = "keybox",
              strong("One drug, four concentrations. "),
              "Plasma concentration predicts almost nothing in this disease. ",
              "In the figures below, note that the plasma / ELF / intracellular concentrations of the same drug ",
              "differ by orders of magnitude, and that the intracellular concentration of IV amikacin is ",
              "exactly zero."),
          fluidRow(
            column(6, plotOutput("p_pk_macro", height = 300)),
            column(6, plotOutput("p_pk_amk",  height = 300))),
          fluidRow(
            column(6, plotOutput("p_pk_other", height = 300)),
            column(6, plotOutput("p_enz",      height = 300))),
          hr(), h4("Steady-state exposure summary (day 180–365)"), DTOutput("t_pk")),

        ## ---- TAB 3 ---------------------------------------------------------
        tabPanel("3. Bacterial niches",
          br(),
          div(class = "keybox",
              strong("Fibrocavitary versus nodular-bronchiectatic is not a different model but a ",
                     "different initial condition of the same model. "),
              "CAVFLAG changes exactly one thing — the initial value of B_C (the caseum compartment). ",
              "Because oral therapy does not reach B_C (macrolide penetration 0.15, ",
              "amikacin 0.02), the gap in culture conversion rate between the two phenotypes arises arithmetically."),
          plotOutput("p_niche", height = 380),
          fluidRow(
            column(6, plotOutput("p_sputum", height = 300)),
            column(6, plotOutput("p_nichebar", height = 300)))),

        ## ---- TAB 4 ---------------------------------------------------------
        tabPanel("4. Clinical endpoints",
          br(),
          fluidRow(
            column(6, plotOutput("p_conv",  height = 300)),
            column(6, plotOutput("p_qolb",  height = 300))),
          fluidRow(
            column(6, plotOutput("p_cav",   height = 300)),
            column(6, plotOutput("p_bmi",   height = 300))),
          hr(), h4("Key endpoints"), DTOutput("t_end")),

        ## ---- TAB 5 ---------------------------------------------------------
        tabPanel("5. Resistance evolution",
          br(),
          div(class = "keybox",
              strong("Resistance is not a switch but a division. "),
              "Φ = (macrolide kill rate) / (total kill rate). Where Φ approaches 1 in a niche, the regimen ",
              "in that niche is 'functional macrolide monotherapy', and because MAC carries a single copy ",
              "of the rrl gene, one point mutation confers complete ",
              "resistance. Inside the phagosome the companion drugs do not arrive, so Φ_I is structurally ",
              "high — adding ALIS fills that denominator for the first time. ",
              br(), br(),
              "On verification, however, Φ approaching 1 is not by itself enough for resistance to take hold. ",
              "For resistant organisms actually to predominate, ", strong("two conditions simultaneously"), " are required — ",
              "① selection pressure (Φ→1) and ② host conditions under which the net growth rate of the resistant population is positive",
              "(high bacterial burden + impaired IFN-γ / mucociliary clearance). Only in scenario 12 does ",
              "resistance sweep to 100%."),
          fluidRow(
            column(6, plotOutput("p_res",  height = 320)),
            column(6, plotOutput("p_phi",  height = 320))),
          plotOutput("p_respop", height = 300)),

        ## ---- TAB 6 ---------------------------------------------------------
        tabPanel("6. Scenario comparison",
          br(),
          checkboxGroupInput("scn", "Select the regimens to compare",
                             choices = names(SCENARIOS),
                             selected = names(SCENARIOS)[c(1, 2, 3, 5, 6)],
                             inline = FALSE),
          actionButton("go_scn", "Run scenario comparison", class = "btn-primary"),
          br(), br(),
          plotOutput("p_scn_sput", height = 340),
          fluidRow(
            column(6, plotOutput("p_scn_res", height = 300)),
            column(6, plotOutput("p_scn_oto", height = 300))),
          hr(), h4("Scenario summary table"), DTOutput("t_scn")),

        ## ---- TAB 7 ---------------------------------------------------------
        tabPanel("7. Biomarkers",
          br(),
          fluidRow(
            column(6, plotOutput("p_immune", height = 300)),
            column(6, plotOutput("p_protease", height = 300))),
          fluidRow(
            column(6, plotOutput("p_airway", height = 300)),
            column(6, plotOutput("p_symp",   height = 300)))),

        ## ---- TAB 8 ---------------------------------------------------------
        tabPanel("8. Toxicity · safety",
          br(),
          div(class = "keybox",
              strong("Ototoxicity follows the plasma, not the lung. "),
              "Perilymph amikacin enters from KCEN (plasma), while efficacy comes from KELF/KMAC (lung). ",
              "Because the two routes are separated, 'high in the lung, low in the ear' for ALIS is ",
              "a result rather than an assumption. Check the lung:ear ratio below."),
          fluidRow(
            column(6, plotOutput("p_oto",  height = 300)),
            column(6, plotOutput("p_lungear", height = 300))),
          fluidRow(
            column(6, plotOutput("p_tox_other", height = 300)),
            column(6, plotOutput("p_qtc", height = 300))),
          hr(), h4("Peak toxicity indices"), DTOutput("t_tox")),

        ## ---- TAB 9 ---------------------------------------------------------
        tabPanel("9. pH paradox explorer",
          br(),
          div(class = "keybox",
              strong("One input, two consequences. "),
              "Changing the phagosomal pH moves (i) the Henderson-Hasselbalch ion-trapping ratio R_trap and ",
              "(ii) the macrolide MIC-rise factor at the same time. Below the point where the two curves ",
              "cross lies the region in which 'the more that enters the cell, the less it works'."),
          fluidRow(
            column(7, plotOutput("p_paradox", height = 380)),
            column(5, br(), DTOutput("t_paradox"))),
          hr(),
          sliderInput("pka", "Macrolide pKa", 7.5, 9.8, 8.7, 0.1, width = "60%"),
          sliderInput("gm2", "MIC rise per pH unit (log10, macrolide GM)",
                      0.1, 1.2, 0.5, 0.05, width = "60%"),
          plotOutput("p_paradox_net", height = 300)),

        ## ---- TAB 10 --------------------------------------------------------
        tabPanel("10. Sensitivity analysis",
          br(),
          helpText("For the currently configured regimen, how much the 12-month sputum burden ",
                   "(log10 CFU/mL) moves when each parameter is varied by ±30%."),
          actionButton("go_sens", "Run sensitivity analysis (takes a few seconds)",
                       class = "btn-warning"),
          br(), br(),
          plotOutput("p_sens", height = 460),
          DTOutput("t_sens")),

        ## ---- TAB 11 --------------------------------------------------------
        tabPanel("11. Model information",
          br(),
          h4("Structure"),
          tags$ul(
            tags$li(strong("46 ODEs"), " — 6 drug PK (oral macrolide · EMB · RIF · CFZ, ",
                    "IV amikacin, inhaled liposomal amikacin) + CYP3A induction + 6 bacterial compartments ",
                    "+ host immune/airway/structure + 5 toxicity compartments"),
            tags$li(strong("4 physical niches"), " — B_E (free organisms in the airway lumen) · B_B (biofilm) · ",
                    "B_I (macrophage phagosome) · B_C (caseum/cavity wall). Each drug has a different ",
                    "effective concentration in each niche."),
            tags$li(strong("3 resistance compartments"), " — rrl A2058G mutants are generated in proportion to replication in the extracellular / biofilm / ",
                    "intracellular compartments respectively, and are selected by the per-niche kill ratio. ",
                    "The initial value is MU × the burden of that compartment, so where the burden is low ",
                    "the mutant does not exist in the first place."),
            tags$li(strong("12 scenarios"), " — watchful waiting, guideline three-drug (TIW/daily), IV amikacin, ",
                    "ALIS, monotherapy, CYP3A interaction, rifampicin exclusion, non-pharmacological intervention, ",
                    "immunocompromised host, resistance sweep")),
          h4("Derived values (not parameters)"),
          tags$ul(
            tags$li("R_trap = (1+10^(pKa−pH_in))/(1+10^(pKa−pH_out)) ≈ 151× — ",
                    "macrophage accumulation of azithromycin"),
            tags$li("MIC rise = 10^(GM×(7.4−pH_in)) ≈ 12.6× — loss of activity at the same pH"),
            tags$li("Net potency gain = 151 / 12.6 ≈ 12× — the real value of '100–1000-fold lung penetration'"),
            tags$li("Φ_I = macrolide kill fraction inside the phagosome — the size of the resistance selection pressure")),
          h4("Limitations"),
          tags$ul(
            tags$li("This is a deterministic single-patient model. Inter-individual variability is not included, so ",
                    "culture conversion is expressed as threshold crossing rather than as a probability. Stochastic extinction ",
                    "is approximated by an extinction floor (EXFLOOR) that stops replication below 1 CFU."),
            tags$li("Reinfection and relapse are not distinguished (a distinction that requires genotyping)."),
            tags$li("Parameters are literature-based approximations and have not been fitted to ",
                    "patient data. They cannot be used for clinical decision-making.")),
          h4("References"),
          p("For the full list see ", code("ntm_references.md"), " (76 entries, 12 sections)."),
          hr(),
          p(em("Claude Code Routine · QSP Disease Model Library · 2026-08-01")))
      )
    )
  )
)

## =============================================================================
##  SERVER
## =============================================================================
server <- function(input, output, session) {

  cur_params <- reactive({
    list(CAVFLAG = as.numeric(input$cav),
         MACTYPE = ifelse(input$macro == "CLR", 1, 0),
         WTBL    = input$wt,
         IFNCAP  = input$ifncap,
         MCC0    = input$mcc0,
         PHPHAG  = input$phag,
         ACT     = as.numeric(input$act),
         NUTR    = as.numeric(input$nutr))
  })

  sim <- eventReactive(input$go, {
    ev_obj <- build_events(input$tx, input$macro, input$freq, input$emb,
                           input$rif, input$cfz, input$ivamk,
                           if (is.null(input$ivwks)) 0 else input$ivwks,
                           input$alis, input$wt)
    run_sim(cur_params(), ev_obj)
  }, ignoreNULL = FALSE)

  ## ---- helpers -------------------------------------------------------------
  lp <- function(d, ycol, ylab, title, sub = NULL, logy = FALSE, col = PAL[1]) {
    g <- ggplot(d, aes(time, .data[[ycol]])) +
      geom_line(linewidth = 0.9, colour = col) +
      labs(x = "Time (days)", y = ylab, title = title, subtitle = sub) + THEME
    if (logy) g <- g + scale_y_log10()
    g
  }

  ## ================= TAB 1 =================================================
  output$p_profile <- renderPlot({
    d <- sim()
    d %>% select(time, WT, BMI) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      scale_colour_manual(values = PAL) +
      labs(x = "Time (days)", y = NULL,
           title = "Body weight / BMI — the TNF-driven wasting loop",
           subtitle = "Low body weight is both a risk factor and a consequence (a bidirectional arrow)") +
      THEME + theme(legend.position = "none")
  })

  output$p_mcc <- renderPlot({
    d <- sim()
    d %>% select(time, MCC, MUC, BRO) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      scale_colour_manual(values = PAL[c(3, 4, 2)]) +
      labs(x = "Time (days)", y = NULL,
           title = "Airway defence and structural damage (the Cole vicious circle)",
           subtitle = "MCC ↓ → mucus ↑ → bronchiectasis ↑ → MCC ↓ again") +
      THEME + theme(legend.position = "none")
  })

  output$t_profile <- renderDT({
    p <- cur_params()
    tibble(Item = c("Disease phenotype", "Body weight", "Th1/IFN-γ function", "Baseline MCC",
                    "Phagosomal pH", "Macrolide", "Dosing frequency", "EMB", "RIF",
                    "CFZ", "ALIS", "IV amikacin", "Airway clearance", "Nutritional support",
                    "Treatment duration"),
           Value = c(ifelse(input$cav == 1, "Fibrocavitary", "Nodular-bronchiectatic"),
                  paste(input$wt, "kg"), input$ifncap, input$mcc0, input$phag,
                  input$macro, input$freq,
                  ifelse(input$emb, "Yes", "No"),
                  ifelse(input$rif, "Yes", "No"),
                  ifelse(input$cfz, "Yes", "No"),
                  ifelse(input$alis, "Yes", "No"),
                  ifelse(input$ivamk, paste0("Yes (", input$ivwks, " weeks)"), "No"),
                  ifelse(input$act, "Yes", "No"),
                  ifelse(input$nutr, "Yes", "No"),
                  paste(input$tx, "days"))) %>%
      datatable(rownames = FALSE, options = list(dom = "t", pageLength = 20))
  })

  ## ================= TAB 2 =================================================
  output$p_pk_macro <- renderPlot({
    sim() %>% select(time, `Plasma` = CMPo, `ELF` = CMEo, `Intracellular (macrophage)` = CMIo) %>%
      pivot_longer(-time) %>% filter(value > 0) %>%
      ggplot(aes(time, value, colour = name)) +
      geom_line(linewidth = 0.9) + scale_y_log10() +
      scale_colour_manual(values = PAL[c(1, 3, 5)]) +
      labs(x = "Time (days)", y = "Concentration (mg/L, log)",
           title = "Macrolide: the same drug, three orders of magnitude apart",
           subtitle = "The intracellular concentration is derived from the pKa and the phagosomal pH", colour = NULL) +
      THEME
  })

  output$p_pk_amk <- renderPlot({
    sim() %>% select(time, `Plasma` = CKPo, `ELF (free)` = CKEo,
                     `Intracellular (via liposome)` = CKIo, `Perilymph` = KPERI) %>%
      pivot_longer(-time) %>%
      mutate(value = pmax(value, 1e-4)) %>%
      ggplot(aes(time, value, colour = name)) +
      geom_line(linewidth = 0.9) + scale_y_log10() +
      scale_colour_manual(values = PAL[c(1, 3, 5, 2)]) +
      labs(x = "Time (days)", y = "Concentration (mg/L, log)",
           title = "Amikacin: the same molecule, a different address",
           subtitle = "Give IV only and the 'intracellular' curve sits on the floor",
           colour = NULL) + THEME
  })

  output$p_pk_other <- renderPlot({
    sim() %>% select(time, `Ethambutol` = CEPo, `Rifampicin` = CFPo,
                     `Clofazimine` = CCPo) %>%
      pivot_longer(-time) %>% mutate(value = pmax(value, 1e-4)) %>%
      ggplot(aes(time, value, colour = name)) +
      geom_line(linewidth = 0.9) + scale_y_log10() +
      scale_colour_manual(values = PAL[c(3, 4, 8)]) +
      labs(x = "Time (days)", y = "Plasma concentration (mg/L, log)",
           title = "Companion drug plasma concentrations", colour = NULL) + THEME
  })

  output$p_enz <- renderPlot({
    lp(sim(), "ENZ", "Relative CYP3A4 enzyme amount (baseline = 1)",
       "CYP3A induction by rifampicin",
       "Choose clarithromycin and this curve shaves that drug's exposure directly",
       col = PAL[4])
  })

  output$t_pk <- renderDT({
    d <- sim() %>% filter(time >= 180, time <= 365)
    tibble(
      Metric = c("Macrolide plasma (mg/L)", "Macrolide ELF (mg/L)",
               "Macrolide intracellular (mg/L)", "Intracellular:plasma ratio",
               "Amikacin ELF (mg/L)", "Amikacin intracellular (mg/L)",
               "Amikacin plasma (mg/L)", "Perilymph amikacin (mg/L)",
               "Extracellular C/MIC", "Intracellular C/MIC", "CYP3A induction fold"),
      Mean = round(c(mean(d$CMPo), mean(d$CMEo), mean(d$CMIo),
                     mean(d$CMIo) / max(mean(d$CMPo), 1e-9),
                     mean(d$CKEo), mean(d$CKIo), mean(d$CKPo), mean(d$KPERI),
                     mean(d$CME_MIC), mean(d$CMI_MIC), mean(d$ENZ)), 3)) %>%
      datatable(rownames = FALSE, options = list(dom = "t", pageLength = 15))
  })

  ## ================= TAB 3 =================================================
  output$p_niche <- renderPlot({
    sim() %>% select(time, `B_E airway lumen` = BE, `B_B biofilm` = BB,
                     `B_I intracellular` = BI, `B_C caseum/cavity` = BC) %>%
      pivot_longer(-time, names_to = "niche", values_to = "cfu") %>%
      mutate(cfu = pmax(cfu, 1)) %>%
      ggplot(aes(time, log10(cfu), colour = niche)) +
      geom_line(linewidth = 1) +
      scale_colour_manual(values = c("#FFB74D", "#FF8A65", "#FF7043", "#D84315")) +
      labs(x = "Time (days)", y = "log10 CFU", colour = NULL,
           title = "Where the bacteria are — the time course of the four niches",
           subtitle = "The curve that rises again after oral therapy ends is the reservoir of relapse") +
      THEME
  })

  output$p_sputum <- renderPlot({
    d <- sim()
    ggplot(d, aes(time, LOGSPUT)) +
      geom_line(linewidth = 1, colour = PAL[1]) +
      geom_hline(yintercept = 1, linetype = "dashed", colour = PAL[2]) +
      annotate("text", x = max(d$time) * 0.7, y = 1.35,
               label = "Culture-negative threshold", size = 3.4, colour = PAL[2]) +
      labs(x = "Time (days)", y = "Sputum log10 CFU/mL",
           title = "Sputum burden and culture conversion") + THEME
  })

  output$p_nichebar <- renderPlot({
    d <- sim()
    tp <- c(0, 90, 180, 365, max(d$time))
    d %>% filter(time %in% tp) %>%
      select(time, BE, BB, BI, BC) %>%
      pivot_longer(-time, names_to = "niche", values_to = "cfu") %>%
      mutate(cfu = pmax(cfu, 1), time = factor(time)) %>%
      ggplot(aes(time, log10(cfu), fill = niche)) +
      geom_col(position = "dodge") +
      scale_fill_manual(values = c("#FFB74D", "#FF8A65", "#FF7043", "#D84315")) +
      labs(x = "Time point (days)", y = "log10 CFU", fill = NULL,
           title = "Niche composition by time point") + THEME
  })

  ## ================= TAB 4 =================================================
  output$p_conv <- renderPlot({
    d <- sim()
    ggplot(d, aes(time, TNEG)) +
      geom_line(linewidth = 1, colour = PAL[3]) +
      labs(x = "Time (days)", y = "Cumulative culture-negative days",
           title = "Duration of maintained culture conversion",
           subtitle = "Guideline target: maintained for 12 months (365 days) after conversion") + THEME
  })
  output$p_qolb <- renderPlot(lp(sim(), "QOLB", "QOL-B respiratory domain (0–100)",
                                 "Quality of life", "MCID ≈ 10 points", col = PAL[5]))
  output$p_cav  <- renderPlot(lp(sim(), "CAV", "Cavity volume (cm³)",
                                 "Cavity size", col = PAL[2]))
  output$p_bmi  <- renderPlot(lp(sim(), "BMI", "BMI (kg/m²)",
                                 "Body mass index", col = PAL[6]))

  output$t_end <- renderDT({
    d <- sim()
    at <- function(t, col) d[[col]][which.min(abs(d$time - t))]
    tibble(
      Endpoint = c("Sputum log10 CFU/mL (baseline)", "Sputum (6 months)", "Sputum (12 months)",
                     "Sputum (18 months, after end of treatment)",
                     "Culture conversion within 6 months", "Culture conversion within 12 months",
                     "Cumulative negative days", "Resistant fraction (12 months)",
                     "Cavity volume (12 months, cm³)", "QOL-B (12 months)", "BMI (12 months)"),
      Value = c(round(at(0, "LOGSPUT"), 2), round(at(180, "LOGSPUT"), 2),
             round(at(365, "LOGSPUT"), 2), round(at(540, "LOGSPUT"), 2),
             ifelse(any(d$CULTNEG[d$time <= 180] == 1), "Yes", "No"),
             ifelse(any(d$CULTNEG[d$time <= 365] == 1), "Yes", "No"),
             round(max(d$TNEG), 0), signif(at(365, "RESFRAC"), 3),
             round(at(365, "CAV"), 1), round(at(365, "QOLB"), 1),
             round(at(365, "BMI"), 2))) %>%
      datatable(rownames = FALSE, options = list(dom = "t", pageLength = 15))
  })

  ## ================= TAB 5 =================================================
  output$p_res <- renderPlot({
    sim() %>% mutate(RESFRAC = pmax(RESFRAC, 1e-12)) %>%
      ggplot(aes(time, RESFRAC)) +
      geom_line(linewidth = 1, colour = PAL[2]) + scale_y_log10() +
      labs(x = "Time (days)", y = "Resistant fraction (log)",
           title = "rrl A2058G resistant fraction",
           subtitle = "MAC carries a single rrl copy — one point mutation gives complete resistance") + THEME
  })

  output$p_phi <- renderPlot({
    sim() %>% select(time, `Φ_E extracellular` = PHI_E, `Φ_I intracellular` = PHI_I) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
      geom_line(linewidth = 1) +
      geom_hline(yintercept = 0.9, linetype = "dashed", colour = "grey40") +
      scale_colour_manual(values = PAL[c(1, 2)]) + ylim(0, 1) +
      labs(x = "Time (days)", y = "Φ = macrolide kill fraction", colour = NULL,
           title = "The per-niche resistance selection gate",
           subtitle = "Above the dashed line (0.9) = effectively macrolide monotherapy") + THEME
  })

  output$p_respop <- renderPlot({
    sim() %>% select(time, `Susceptible extracellular` = BE, `Resistant extracellular` = RE,
                     `Susceptible biofilm` = BB, `Resistant biofilm` = RB,
                     `Susceptible intracellular` = BI, `Resistant intracellular` = RI) %>%
      pivot_longer(-time) %>% mutate(value = pmax(value, 1e-3)) %>%
      ggplot(aes(time, log10(value), colour = name)) +
      geom_line(linewidth = 0.9) +
      scale_colour_manual(values = PAL[c(1, 2, 3, 4, 6, 7)]) +
      labs(x = "Time (days)", y = "log10 CFU", colour = NULL,
           title = "Susceptible versus resistant populations (extracellular / biofilm / intracellular)") + THEME
  })

  ## ================= TAB 6 =================================================
  scn_sim <- eventReactive(input$go_scn, {
    req(length(input$scn) > 0)
    withProgress(message = "Running scenario simulations...", {
      lapply(input$scn, function(nm) {
        s <- SCENARIOS[[nm]]
        p <- modifyList(list(WTBL = input$wt, PHPHAG = input$phag), s$p)
        incProgress(1 / length(input$scn), detail = nm)
        run_sim(p, scenario_events(s, input$wt)) %>% mutate(scenario = nm)
      }) %>% bind_rows()
    })
  })

  output$p_scn_sput <- renderPlot({
    ggplot(scn_sim(), aes(time, LOGSPUT, colour = scenario)) +
      geom_line(linewidth = 0.9) +
      geom_hline(yintercept = 1, linetype = "dashed") +
      labs(x = "Time (days)", y = "Sputum log10 CFU/mL", colour = NULL,
           title = "Sputum burden by regimen") +
      THEME + theme(legend.text = element_text(size = 8))
  })

  output$p_scn_res <- renderPlot({
    scn_sim() %>% mutate(RESFRAC = pmax(RESFRAC, 1e-12)) %>%
      ggplot(aes(time, RESFRAC, colour = scenario)) +
      geom_line(linewidth = 0.9) + scale_y_log10() +
      labs(x = "Time (days)", y = "Resistant fraction", colour = NULL,
           title = "Emergence of resistance by regimen") +
      THEME + theme(legend.position = "none")
  })

  output$p_scn_oto <- renderPlot({
    ggplot(scn_sim(), aes(time, HEARDB, colour = scenario)) +
      geom_line(linewidth = 0.9) +
      labs(x = "Time (days)", y = "Hearing threshold shift (dB)", colour = NULL,
           title = "Ototoxicity by regimen") +
      THEME + theme(legend.position = "none")
  })

  output$t_scn <- renderDT({
    scn_sim() %>% group_by(scenario) %>%
      summarise(`Sputum 6 mo` = round(LOGSPUT[which.min(abs(time - 180))], 2),
                `Sputum 12 mo` = round(LOGSPUT[which.min(abs(time - 365))], 2),
                `Sputum 18 mo` = round(LOGSPUT[which.min(abs(time - 540))], 2),
                `6 mo conversion` = ifelse(any(CULTNEG[time <= 180] == 1), "○", "×"),
                `12 mo conversion` = ifelse(any(CULTNEG[time <= 365] == 1), "○", "×"),
                `Negative days` = round(max(TNEG)),
                `Resistant fraction (12M)` = signif(RESFRAC[which.min(abs(time - 365))], 3),
                `Cavity (12M)` = round(CAV[which.min(abs(time - 365))], 1),
                `QOL-B(12M)` = round(QOLB[which.min(abs(time - 365))], 1),
                `Hearing (dB)` = round(max(HEARDB), 1),
                `QTc(ms)` = round(max(QTCMS), 1),
                .groups = "drop") %>%
      datatable(rownames = FALSE, options = list(dom = "t", scrollX = TRUE))
  })

  ## ================= TAB 7 =================================================
  output$p_immune <- renderPlot({
    sim() %>% select(time, `Macrophages` = MPH, `IFN-γ` = IFNG, `TNF-α` = TNF) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      scale_colour_manual(values = PAL[c(1, 3, 2)]) +
      labs(x = "Time (days)", y = "Relative units", colour = NULL,
           title = "Immune axis") + THEME
  })
  output$p_protease <- renderPlot({
    sim() %>% select(time, `Neutrophils` = NEU, `MMP-1/9` = MMP) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      scale_colour_manual(values = PAL[c(4, 7)]) +
      labs(x = "Time (days)", y = "Relative units", colour = NULL,
           title = "Neutrophil–protease axis",
           subtitle = "The anti-inflammatory action of the macrolide lowers this before the bacteria die") +
      THEME
  })
  output$p_airway <- renderPlot({
    sim() %>% select(time, `Mucus burden` = MUC, `MCC` = MCC,
                     `Bronchiectasis score` = BRO) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      scale_colour_manual(values = PAL[c(6, 3, 2)]) +
      labs(x = "Time (days)", y = NULL, title = "Airway state") +
      THEME + theme(legend.position = "none")
  })
  output$p_symp <- renderPlot(lp(sim(), "SYM", "Symptom score (0–10)",
                                 "Symptoms (cough · sputum · fatigue · fever)", col = PAL[7]))

  ## ================= TAB 8 =================================================
  output$p_oto <- renderPlot(lp(sim(), "HEARDB", "Hearing threshold shift (dB)",
                                "Ototoxicity — driven by cumulative perilymph exposure", col = PAL[2]))
  output$p_lungear <- renderPlot({
    sim() %>% mutate(LUNGEAR = pmax(LUNGEAR, 1e-3)) %>%
      ggplot(aes(time, LUNGEAR)) +
      geom_line(linewidth = 1, colour = PAL[5]) + scale_y_log10() +
      labs(x = "Time (days)", y = "ELF : perilymph amikacin ratio (log)",
           title = "Lung : ear ratio — the reason ALIS exists",
           subtitle = "With IV this ratio is low; with ALIS it rises by orders of magnitude") +
      THEME
  })
  output$p_tox_other <- renderPlot({
    sim() %>% select(time, `Optic nerve (EMB)` = OPT, `Renal tubule (amikacin)` = NEP,
                     `Liver (RIF)` = HEP) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      scale_colour_manual(values = PAL[c(3, 5, 4)]) +
      labs(x = "Time (days)", y = "Injury index", colour = NULL,
           title = "Other organ toxicities") + THEME
  })
  output$p_qtc <- renderPlot(lp(sim(), "QTCMS", "QTc prolongation (ms)",
                                "ECG — macrolide hERG block",
                                "Clarithromycin > azithromycin", col = PAL[8]))

  output$t_tox <- renderDT({
    d <- sim()
    tibble(Toxicity = c("Hearing threshold shift (dB)", "Optic nerve injury index (0–1)",
                    "Renal tubular injury index (0–1)", "Hepatic injury index", "Peak QTc prolongation (ms)",
                    "Peak perilymph amikacin (mg/L)", "Peak plasma amikacin (mg/L)"),
           Peak = round(c(max(d$HEARDB), max(d$OPT), max(d$NEP), max(d$HEP),
                            max(d$QTCMS), max(d$KPERI), max(d$CKPo)), 3)) %>%
      datatable(rownames = FALSE, options = list(dom = "t"))
  })

  ## ================= TAB 9 =================================================
  paradox_df <- reactive({
    pH <- seq(4.5, 7.4, by = 0.05)
    tibble(pH = pH,
           R_trap = (1 + 10^(input$pka - pH)) / (1 + 10^(input$pka - 7.4)),
           MIC_fold = 10^(input$gm2 * (7.4 - pH))) %>%
      mutate(net = R_trap / MIC_fold)
  })

  output$p_paradox <- renderPlot({
    paradox_df() %>%
      select(pH, `R_trap (accumulation)` = R_trap, `MIC rise factor (loss of activity)` = MIC_fold) %>%
      pivot_longer(-pH) %>%
      ggplot(aes(pH, value, colour = name)) +
      geom_line(linewidth = 1.1) + scale_y_log10() +
      geom_vline(xintercept = input$phag, linetype = "dashed", colour = "grey30") +
      annotate("text", x = input$phag, y = 1.2, hjust = -0.1, size = 3.5,
               label = paste0("Current phagosomal pH = ", input$phag)) +
      scale_colour_manual(values = c(PAL[1], PAL[2])) +
      labs(x = "Phagosomal pH", y = "Fold (log)", colour = NULL,
           title = "One input, two consequences",
           subtitle = "Lower the pH and more enters (blue) while at the same time it works less well (red)") +
      THEME
  })

  output$p_paradox_net <- renderPlot({
    paradox_df() %>%
      ggplot(aes(pH, net)) +
      geom_line(linewidth = 1.2, colour = PAL[5]) +
      geom_hline(yintercept = 1, linetype = "dashed", colour = PAL[2]) +
      geom_vline(xintercept = input$phag, linetype = "dashed", colour = "grey30") +
      scale_y_log10() +
      labs(x = "Phagosomal pH", y = "Net potency gain = R_trap / MIC rise factor (log)",
           title = "What actually remains",
           subtitle = "Below the red line (=1) intracellular penetration is a net loss") +
      THEME
  })

  output$t_paradox <- renderDT({
    d <- paradox_df()
    row <- d[which.min(abs(d$pH - input$phag)), ]
    tibble(Item = c("Phagosomal pH", "Macrolide pKa", "Ion-trapping accumulation ratio R_trap",
                    "MIC rise factor", "Net potency gain",
                    "Real gain versus extracellular"),
           Value = c(input$phag, input$pka, round(row$R_trap, 1),
                  round(row$MIC_fold, 2), round(row$net, 1),
                  paste0(round(row$net, 1), " fold"))) %>%
      datatable(rownames = FALSE, options = list(dom = "t"))
  })

  ## ================= TAB 10 ================================================
  sens_res <- eventReactive(input$go_sens, {
    ps <- c("PHPHAG", "GM", "GK", "PCSM", "PCSK", "PBFK", "PBFL", "KUPT",
            "TOLCS", "TOLBF", "MU", "IFNCAP", "MCC0", "EMAXPERM", "RMELF")
    ev_obj <- build_events(input$tx, input$macro, input$freq, input$emb,
                           input$rif, input$cfz, input$ivamk,
                           if (is.null(input$ivwks)) 0 else input$ivwks,
                           input$alis, input$wt)
    base_pars <- cur_params()
    base_out <- run_sim(base_pars, ev_obj, end = 365)
    base_val <- base_out$LOGSPUT[nrow(base_out)]
    withProgress(message = "Running sensitivity analysis...", {
      lapply(ps, function(p) {
        incProgress(1 / length(ps), detail = p)
        b <- as.numeric(param(mod)[[p]])
        if (p %in% names(base_pars)) b <- as.numeric(base_pars[[p]])
        vapply(c(0.7, 1.3), function(f) {
          pp <- base_pars; pp[[p]] <- b * f
          o <- run_sim(pp, ev_obj, end = 365)
          o$LOGSPUT[nrow(o)]
        }, numeric(1)) -> v
        tibble(parameter = p, low = v[1], high = v[2], base = base_val,
               swing = abs(v[2] - v[1]))
      }) %>% bind_rows() %>% arrange(desc(swing))
    })
  })

  output$p_sens <- renderPlot({
    d <- sens_res()
    d %>% mutate(parameter = factor(parameter, levels = rev(d$parameter))) %>%
      ggplot() +
      geom_segment(aes(y = parameter, yend = parameter, x = low, xend = high),
                   linewidth = 5, colour = PAL[1], alpha = 0.55) +
      geom_point(aes(y = parameter, x = base), colour = PAL[2], size = 2.6) +
      labs(x = "12-month sputum log10 CFU/mL", y = NULL,
           title = "Local sensitivity — parameters ±30%",
           subtitle = "Red point = reference value. The longer the bar, the more the parameter moves the outcome") +
      THEME
  })

  output$t_sens <- renderDT({
    sens_res() %>% mutate(across(where(is.numeric), ~round(.x, 3))) %>%
      datatable(rownames = FALSE, options = list(dom = "tp", pageLength = 15))
  })
}

shinyApp(ui, server)
