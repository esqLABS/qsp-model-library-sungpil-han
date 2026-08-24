# =============================================================================
# Drug-Induced Liver Injury (DILI) — QSP simulator (Shiny)
# =============================================================================
#
# A dashboard over dili_mrgsolve_model.R. The app is deliberately organised
# around the model's thesis rather than around its compartments:
#
#   Tab 2  the drug arrives at some RATE
#   Tab 3  the liver disposes of it at some RATE, and that rate has a ceiling
#   Tab 4  when production outruns disposal, a POSITIVE FEEDBACK latches
#   Tab 5  the latch converts stressed hepatocytes into lost mass
#   Tab 6  ALT reports the rate, bilirubin reports the reserve — Hy's Law
#   Tab 7  the bile-acid and immune arms reach the same death node by other routes
#   Tab 8  the antidote window is where the trajectory crosses the separatrix
#   Tab 9  side-by-side scenarios
#
# Run:   R -e "shiny::runApp('dili_shiny_app.R', port=8080)"
# Needs: shiny, mrgsolve, ggplot2, dplyr, tidyr, DT
# =============================================================================

library(shiny)
library(mrgsolve)
library(ggplot2)
library(dplyr)
library(tidyr)
library(DT)

source("dili_mrgsolve_model.R", local = TRUE)   # defines mod, sim_dili, ...

THEME <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey93", colour = NA),
        legend.position = "bottom")

PAL <- c("#1F77B4", "#D62728", "#2CA02C", "#FF7F0E", "#9467BD",
         "#8C564B", "#17BECF", "#7F7F7F")

# -----------------------------------------------------------------------------
# UI
# -----------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Drug-Induced Liver Injury (DILI) QSP Simulator"),
  tags$p(style = "color:#8B0000; font-weight:600; margin-top:-8px;",
         "Core thesis: liver injury is a matter of rate rather than dose, and ",
         "because of the JNK–Sab positive feedback loop the system is bistable. ",
         "So the same total exposure can be asymptomatic or it can be fatal."),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      h4("① Drug"),
      selectInput("archetype", "Drug archetype",
                  choices = c("A · bioactivation type (APAP)"     = "A",
                              "B · BSEP inhibitor (cholestasis)" = "B",
                              "C · idiosyncratic immune (HLA)"   = "C"),
                  selected = "A"),
      numericInput("dose", "Dose per administration (mg/kg)", 350, min = 0, max = 1000, step = 10),
      selectInput("regimen", "Dosing pattern",
                  choices = c("Single acute ingestion" = "single",
                              "Every 6 h x 7 days"     = "q6h7d",
                              "Every 12 h x 28 days"   = "q12h28d",
                              "Every 24 h x 56 days"   = "q24h56d",
                              "Spread over a set time" = "spread"),
                  selected = "single"),
      conditionalPanel("input.regimen == 'spread'",
        sliderInput("spread_h", "Total ingestion time (h)", 1, 72, 12, step = 1)),

      hr(),
      h4("② Host"),
      sliderInput("fcyp", "CYP2E1 induction fold (chronic alcohol)", 1, 3, 1, step = 0.1),
      sliderInput("fgsh", "GSH pool fold (fasting · malnutrition)", 0.4, 1.2, 1, step = 0.05),
      sliderInput("cysbase", "Cysteine supply set-point", 1.0, 2.0, 1.65, step = 0.05),
      sliderInput("agef", "Regenerative capacity (1 = young, 0.5 = elderly)", 0.4, 1.2, 1, step = 0.05),

      hr(),
      h4("③ Drug liabilities"),
      sliderInput("logki", "BSEP inhibition Ki (log10 µM; 9 = none)",
                  -1, 9, 9, step = 0.25),
      sliderInput("func", "Uncoupling / β-oxidation block burden", 0, 0.6, 0, step = 0.05),
      checkboxInput("hla", "Carries an HLA risk allele", FALSE),
      checkboxInput("ici", "Concomitant checkpoint inhibitor (tolerance removed)", FALSE),

      hr(),
      h4("④ Treatment"),
      checkboxInput("use_nac", "N-acetylcysteine (IV, Prescott)", FALSE),
      conditionalPanel("input.use_nac",
        sliderInput("nac_t", "NAC start time (h)", 0, 48, 8, step = 1)),
      checkboxInput("use_ster", "Corticosteroid", FALSE),
      conditionalPanel("input.use_ster",
        sliderInput("ster_t", "Steroid start (day)", 1, 60, 42, step = 1)),

      hr(),
      sliderInput("tend", "Simulation duration (h)", 168, 2400, 336, step = 24),
      actionButton("go", "Run simulation", class = "btn-primary btn-block")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs", type = "tabs",

        # ---------------------------------------------------------------- 1
        tabPanel("① Patient profile / verdict",
          br(),
          fluidRow(
            column(4, wellPanel(h4("Outcome verdict"), htmlOutput("verdict"))),
            column(4, wellPanel(h4("Key metrics"), tableOutput("keytab"))),
            column(4, wellPanel(h4("What actually happened in this run"),
                                htmlOutput("narrative")))
          ),
          hr(),
          h4("Production flux vs disposal flux"),
          p("This figure is the whole model. The area by which the red line (reactive metabolite ",
            "production) exceeds the green line (GSH resynthesis capacity) is the injury."),
          plotOutput("fluxplot", height = 320)
        ),

        # ---------------------------------------------------------------- 2
        tabPanel("② Pharmacokinetics (PK)",
          br(),
          fluidRow(
            column(6, plotOutput("pk_conc", height = 300)),
            column(6, plotOutput("pk_nomogram", height = 300))
          ),
          hr(),
          h4("Metabolic pathway partitioning — the bioactivated fraction rises by itself with dose"),
          p("Sulfation (PAPS) and glucuronidation (UDPGA) saturate as their cofactors are depleted, ",
            "but the CYP route has a high Km and is almost linear. So at high doses the fraction ",
            "going to the reactive metabolite grows automatically — there is no separate switch."),
          plotOutput("pk_partition", height = 300)
        ),

        # ---------------------------------------------------------------- 3
        tabPanel("③ Detoxification capacity (thiol defence)",
          br(),
          fluidRow(
            column(6, plotOutput("gsh_plot", height = 300)),
            column(6, plotOutput("rm_plot", height = 300))
          ),
          hr(),
          h4("GSH is the gate"),
          p("The reactive metabolite is divided between GSH conjugation and covalent protein ",
            "binding purely by competing rates. At GSH 6.5 mM only about 0.6% binds protein, ",
            "whereas at 0.1 mM about 25% does. There is no threshold parameter."),
          plotOutput("gate_plot", height = 300)
        ),

        # ---------------------------------------------------------------- 4
        tabPanel("④ Mitochondria · JNK loop",
          br(),
          fluidRow(
            column(6, plotOutput("mito_plot", height = 300)),
            column(6, plotOutput("jnk_plot", height = 300))
          ),
          hr(),
          h4("Bistability and the separatrix"),
          p("Below is the JNK–Sab loop taken in isolation. When the redox reserve ",
            "g = (0.30 + 0.70·GSH/GSH₀)·NRF2 falls, the OFF state disappears. ",
            "The red dashed line is the separatrix, and once the simulated trajectory rises above it ",
            "the injury is irreversible."),
          plotOutput("bistab_plot", height = 340)
        ),

        # ---------------------------------------------------------------- 5
        tabPanel("⑤ Hepatocyte mass · regeneration",
          br(),
          plotOutput("mass_plot", height = 320),
          hr(),
          fluidRow(
            column(6, plotOutput("damp_plot", height = 280)),
            column(6, plotOutput("regen_plot", height = 280))
          )
        ),

        # ---------------------------------------------------------------- 6
        tabPanel("⑥ Clinical laboratory · Hy's Law",
          br(),
          fluidRow(
            column(6, plotOutput("lab_enz", height = 300)),
            column(6, plotOutput("lab_fn", height = 300))
          ),
          hr(),
          h4("Hy's Law = a joint test of rate × reserve"),
          p("ALT is a low-pass filter on the 'rate of dying', while bilirubin is set by the ",
            "'remaining reserve' (clearance is proportional to surviving mass). Demanding both arms ",
            "at once is demanding that the liver be dying fast and out of reserve at the same time, ",
            "which is why only the conjunction of the two arms predicts death."),
          plotOutput("hylaw_plot", height = 340)
        ),

        # ---------------------------------------------------------------- 7
        tabPanel("⑦ Bile acids · immunity",
          br(),
          fluidRow(
            column(6, plotOutput("bile_plot", height = 300)),
            column(6, plotOutput("rratio_plot", height = 300))
          ),
          hr(),
          fluidRow(
            column(6, plotOutput("innate_plot", height = 280)),
            column(6, plotOutput("adaptive_plot", height = 280))
          )
        ),

        # ---------------------------------------------------------------- 8
        tabPanel("⑧ Antidote time window",
          br(),
          p("This sweeps the NAC start time. The '8-hour rule' is not written into the model as a ",
            "constant — it comes out as the time at which the trajectory crosses the separatrix, ",
            "and so it moves with dose and with host state."),
          actionButton("run_nac", "Compute the NAC time window (about 40 simulations)",
                       class = "btn-warning"),
          br(), br(),
          plotOutput("nacwin_plot", height = 380),
          DTOutput("nacwin_tab")
        ),

        # ---------------------------------------------------------------- 9
        tabPanel("⑨ Scenario comparison",
          br(),
          checkboxGroupInput("scen_pick", "Scenarios to compare",
            choices = c("S2 150 mg/kg untreated"             = "S2",
                        "S3 250 mg/kg untreated"             = "S3",
                        "S4 350 mg/kg untreated"             = "S4",
                        "S5 350 mg/kg + NAC 8h"          = "S5",
                        "S6 350 mg/kg + NAC 16h"         = "S6",
                        "S8 150 mg/kg alcohol · fasting host" = "S8",
                        "S10 BSEP inhibitor 28 days"         = "S10",
                        "S11 HLA+ idiosyncratic 56 days"     = "S11",
                        "S12 + checkpoint inhibitor"         = "S12",
                        "S13 + corticosteroid"               = "S13"),
            selected = c("S3", "S4", "S5", "S8"), inline = TRUE),
          actionButton("run_scen", "Run the selected scenarios", class = "btn-warning"),
          br(), br(),
          plotOutput("scen_plot", height = 420),
          DTOutput("scen_tab")
        ),

        # --------------------------------------------------------------- 10
        tabPanel("⑩ Rate vs dose",
          br(),
          p("The total dose is held fixed and only the rate of arrival is changed. If liver injury ",
            "were a matter of dose, the bars below would all be the same height."),
          actionButton("run_rate", "Run the rate experiment", class = "btn-warning"),
          br(), br(),
          plotOutput("rate_plot", height = 360),
          DTOutput("rate_tab")
        )
      )
    )
  )
)

# -----------------------------------------------------------------------------
# Server
# -----------------------------------------------------------------------------
server <- function(input, output, session) {

  arche_params <- reactive({
    switch(input$archetype,
      A = list(),
      B = DRUG_B,
      C = DRUG_C)
  })

  dose_times <- reactive({
    switch(input$regimen,
      single   = 0,
      q6h7d    = seq(0, 162, by = 6),
      q12h28d  = seq(0, 660, by = 12),
      q24h56d  = seq(0, 1320, by = 24),
      spread   = seq(0, input$spread_h, length.out = max(2, input$spread_h)))
  })

  run <- eventReactive(input$go, {
    n  <- length(dose_times())
    ki <- 10^input$logki
    args <- c(
      list(dose_mgkg  = input$dose * n,
           dose_times = dose_times(),
           nac_start  = if (input$use_nac) input$nac_t else NULL,
           tend       = input$tend,
           FCYP    = input$fcyp,
           FGSH    = input$fgsh,
           CYSBASE = input$cysbase,
           AGEF    = input$agef,
           FUNC    = input$func,
           HLA     = as.numeric(input$hla),
           ICI     = as.numeric(input$ici),
           STER    = as.numeric(input$use_ster),
           STER_T0 = if (input$use_ster) input$ster_t * 24 else 1e9),
      arche_params())
    # explicit BSEP slider overrides the archetype default
    if (input$logki < 9) args$KI_BSEP <- ki
    as.data.frame(do.call(sim_dili, args))
  }, ignoreNULL = FALSE)

  smry <- reactive({
    d <- run()
    list(alt = max(d$ALT), ast = max(d$AST), alp = max(d$ALP),
         tbil = max(d$TBIL), inr = max(d$INR),
         R = max(d$ALT)/40 / (max(d$ALP)/120),
         hy = any(d$HYLAW > 0),
         gsh = min(d$GSHPCT), jnk = max(d$JNK), atp = min(d$ATP),
         lost = max(d$LOST), ba = max(d$BAP), mir = max(d$MIR),
         gmin = min(d$REDOXG), fb = max(d$FBIOACT))
  })

  # ----------------------------------------------------------------- tab 1
  output$verdict <- renderUI({
    s <- smry()
    patt <- if (s$R >= 5) "Hepatocellular" else
            if (s$R <= 2) "Cholestatic" else "Mixed"
    sev <- if (s$lost < 0.01) list("Adaptation — no clinical injury", "#1E7B34")
      else if (s$lost < 0.10) list("Mild injury — recovery expected", "#B8860B")
      else if (s$lost < 0.40) list("Moderate injury", "#E07B00")
      else if (s$lost < 0.70) list("Severe injury", "#B03A2E")
      else list("Fulminant hepatic failure range", "#7D1128")
    HTML(sprintf(
      "<div style='font-size:16px;font-weight:700;color:%s'>%s</div>
       <hr style='margin:8px 0'>
       <b>Injury pattern</b>: %s (R = %.2f)<br>
       <b>Hy's Law</b>: <span style='color:%s;font-weight:700'>%s</span><br>
       <b>Hepatocyte mass lost</b>: %.1f%%<br>
       <b>JNK loop</b>: %s (peak %.3f)<br>
       <b>Minimum redox reserve g</b>: %.2f",
      sev[[2]], sev[[1]], patt, s$R,
      if (s$hy) "#B71C1C" else "#1E7B34", if (s$hy) "Met (YES)" else "Not met",
      100*s$lost,
      if (s$jnk > 0.3) "Latched ON" else "Not latched", s$jnk,
      s$gmin))
  })

  output$keytab <- renderTable({
    s <- smry()
    data.frame(
      Metric = c("Peak ALT (U/L)", "Peak AST (U/L)", "Peak ALP (U/L)",
               "Peak total bilirubin (mg/dL)", "Peak INR", "Minimum GSH (% of baseline)",
               "Minimum ATP (fold)", "Peak serum bile acids (µM)", "Peak miR-122 (fold)",
               "Maximum bioactivated fraction"),
      Value = c(sprintf("%.0f", s$alt), sprintf("%.0f", s$ast),
             sprintf("%.0f", s$alp), sprintf("%.2f", s$tbil),
             sprintf("%.2f", s$inr), sprintf("%.1f", s$gsh),
             sprintf("%.3f", s$atp), sprintf("%.1f", s$ba),
             sprintf("%.0f", s$mir), sprintf("%.1f%%", 100*s$fb)))
  }, striped = TRUE, width = "100%")

  output$narrative <- renderUI({
    s <- smry(); d <- run()
    steps <- c()
    steps <- c(steps, sprintf("At most %.0f%% of hepatic metabolism went to the reactive metabolite.",
                              100*s$fb))
    steps <- c(steps, sprintf("GSH fell to %.0f%% of baseline.", s$gsh))
    steps <- c(steps, if (s$jnk > 0.3)
      sprintf("The redox reserve fell to g = %.2f and the JNK loop latched on.", s$gmin)
      else "The redox reserve was maintained and the JNK loop stayed OFF.")
    steps <- c(steps, if (s$ba > 20)
      sprintf("Bile acids rose to %.0f µM in serum and bile duct injury occurred.", s$ba)
      else "Bile acid homeostasis was maintained.")
    steps <- c(steps, sprintf("In the end %.1f%% of the hepatocyte mass was lost.",
                              100*s$lost))
    HTML(paste0("<ol style='padding-left:18px;margin-bottom:0'>",
                paste0("<li>", steps, "</li>", collapse = ""), "</ol>"))
  })

  output$fluxplot <- renderPlot({
    d <- run()
    p <- as.list(param(mod))
    # reconstruct the two fluxes from captured states
    prod <- p$VMAX_CYP * d$CU/(p$KM_CYP + d$CU) * input$fcyp / p$VLIV / 1000
    disp <- rep(p$VMAX_GSH * input$fgsh, nrow(d)) * (d$NRF2)
    df <- data.frame(time = d$time,
                     `Reactive metabolite production (mM/h)` = prod,
                     `GSH resynthesis capacity (mM/h)`       = disp,
                     check.names = FALSE) %>%
      pivot_longer(-time)
    ggplot(df, aes(time, value, colour = name)) +
      geom_line(linewidth = 1) +
      scale_colour_manual(values = c("#1E7B34", "#B03A2E")) +
      labs(x = "Time (h)", y = "flux (mM/h)", colour = NULL) + THEME
  })

  # ----------------------------------------------------------------- tab 2
  output$pk_conc <- renderPlot({
    d <- run()
    ggplot(d, aes(time)) +
      geom_line(aes(y = CP, colour = "Plasma"), linewidth = 1) +
      geom_line(aes(y = CU, colour = "Free in liver"),
                linewidth = 1) +
      scale_colour_manual(values = PAL[1:2]) +
      labs(title = "Drug concentrations", x = "Time (h)", y = "µM", colour = NULL) + THEME
  })

  output$pk_nomogram <- renderPlot({
    d <- run()
    t <- seq(4, 24, by = 0.5)
    nomo <- data.frame(time = t, y = 150 * exp(-log(2)/4 * (t - 4)))
    ggplot() +
      geom_line(data = nomo, aes(time, y), colour = "#B03A2E",
                linetype = 2, linewidth = 1) +
      geom_line(data = subset(d, time <= 24), aes(time, CP_UGML),
                colour = PAL[1], linewidth = 1) +
      annotate("text", x = 16, y = 60, colour = "#B03A2E",
               label = "Rumack-Matthew treatment line\n(4 h, 150 µg/mL)") +
      scale_y_log10() +
      labs(title = "Rumack-Matthew nomogram (valid only for a single acute ingestion)",
           x = "Time after ingestion (h)", y = "Plasma concentration (µg/mL, log)") + THEME
  })

  output$pk_partition <- renderPlot({
    d <- run(); p <- as.list(param(mod))
    ap <- arche_params()
    vs <- if (!is.null(ap$VMAX_SULT)) ap$VMAX_SULT else p$VMAX_SULT
    vu <- if (!is.null(ap$VMAX_UGT))  ap$VMAX_UGT  else p$VMAX_UGT
    vc <- if (!is.null(ap$VMAX_CYP))  ap$VMAX_CYP  else p$VMAX_CYP
    df <- data.frame(
      time = d$time,
      `Sulfation (SULT)`      = vs*d$CU/(p$KM_SULT+d$CU)*d$PAPS,
      `Glucuronidation (UGT)` = vu*d$CU/(p$KM_UGT +d$CU)*d$UDPGA,
      `CYP bioactivation`     = vc*d$CU/(p$KM_CYP +d$CU)*input$fcyp,
      check.names = FALSE) %>% pivot_longer(-time)
    ggplot(df, aes(time, value, fill = name)) +
      geom_area(position = "fill") +
      scale_fill_manual(values = c("#B03A2E", PAL[3], PAL[4])) +
      scale_y_continuous(labels = scales::percent) +
      labs(x = "Time (h)", y = "Fraction of metabolic flux", fill = NULL) + THEME
  })

  # ----------------------------------------------------------------- tab 3
  output$gsh_plot <- renderPlot({
    d <- run()
    ggplot(d, aes(time)) +
      geom_hline(yintercept = 30, linetype = 3, colour = "grey40") +
      geom_line(aes(y = GSHPCT, colour = "GSH (% of baseline)"), linewidth = 1) +
      geom_line(aes(y = 100*CYS, colour = "Cysteine availability (%)"), linewidth = 1) +
      scale_colour_manual(values = PAL[c(3,4)]) +
      labs(title = "Glutathione and its rate-limiting substrate", x = "Time (h)",
           y = "% of baseline", colour = NULL) + THEME
  })

  output$rm_plot <- renderPlot({
    d <- run()
    ggplot(d, aes(time)) +
      geom_line(aes(y = RM*1000, colour = "Reactive metabolite (nM)"), linewidth = 1) +
      geom_line(aes(y = ADD, colour = "Mitochondrial protein adducts (au)"), linewidth = 1) +
      scale_colour_manual(values = c("#7D1128", "#D62728")) +
      labs(title = "The reactive metabolite and its covalent adducts",
           x = "Time (h)", y = NULL, colour = NULL) + THEME
  })

  output$gate_plot <- renderPlot({
    p <- as.list(param(mod))
    g <- 10^seq(-2, 1, length.out = 300)
    frac <- p$KBIND / (p$KBIND + p$KGST*g + p$KNQO)
    d <- run()
    ggplot(data.frame(g, frac), aes(g, 100*frac)) +
      geom_line(linewidth = 1.1, colour = "#B03A2E") +
      geom_vline(xintercept = min(d$GSH), linetype = 2, colour = PAL[1]) +
      annotate("text", x = min(d$GSH), y = 40, hjust = -0.1, colour = PAL[1],
               label = sprintf("Minimum GSH in this run = %.2f mM", min(d$GSH))) +
      scale_x_log10() +
      labs(title = "The GSH gate: the fraction of reactive metabolite that binds protein covalently",
           x = "Hepatic GSH (mM, log)", y = "Protein-bound fraction (%)") + THEME
  })

  # ----------------------------------------------------------------- tab 4
  output$mito_plot <- renderPlot({
    d <- run()
    ggplot(d, aes(time)) +
      geom_hline(yintercept = 0.45, linetype = 3, colour = "grey40") +
      geom_line(aes(y = MITO, colour = "Mitochondrial function"), linewidth = 1) +
      geom_line(aes(y = ATP,  colour = "ATP"), linewidth = 1) +
      geom_line(aes(y = ROS/5, colour = "ROS / 5"), linewidth = 1) +
      scale_colour_manual(values = PAL[c(1,3,2)]) +
      labs(title = "Bioenergetics (dotted = the half-opening point of the MPT ATP gate)",
           x = "Time (h)", y = "Fold", colour = NULL) + THEME
  })

  output$jnk_plot <- renderPlot({
    d <- run()
    ggplot(d, aes(time)) +
      geom_line(aes(y = JNK,    colour = "Active JNK"), linewidth = 1.1) +
      geom_line(aes(y = REDOXG, colour = "Redox reserve g"), linewidth = 1) +
      geom_line(aes(y = NRF2/3, colour = "NRF2 / 3"), linewidth = 1) +
      scale_colour_manual(values = c("#B03A2E", PAL[1], PAL[3])) +
      labs(title = "The JNK loop and the reserve that restrains it",
           x = "Time (h)", y = NULL, colour = NULL) + THEME
  })

  output$bistab_plot <- renderPlot({
    bb <- analysis_bistability(seq(0.6, 1.6, by = 0.02))
    d  <- run()
    ggplot(bb, aes(g)) +
      geom_line(aes(y = on_state, colour = "Stable ON state"), linewidth = 1.2,
                na.rm = TRUE) +
      geom_line(aes(y = separatrix, colour = "Separatrix (unstable)"),
                linewidth = 1.2, linetype = 2, na.rm = TRUE) +
      geom_point(data = data.frame(g = d$REDOXG, J = d$JNK),
                 aes(g, J), colour = "grey35", size = 0.35, alpha = 0.5) +
      scale_colour_manual(values = c("#B03A2E", "#7D1128")) +
      labs(title = "Bifurcation diagram of the JNK–Sab loop with this run's trajectory (grey)",
           x = "Redox reserve g", y = "Active JNK", colour = NULL) + THEME
  })

  # ----------------------------------------------------------------- tab 5
  output$mass_plot <- renderPlot({
    d <- run()
    df <- data.frame(time = d$time,
                     `Viable hepatocytes (HEP)`    = d$HEP,
                     `Stressed hepatocytes (HEPS)` = d$HEPS,
                     `Cumulative necrosis (NECR)`  = d$NECR,
                     check.names = FALSE) %>% pivot_longer(-time)
    ggplot(df, aes(time, value, colour = name)) +
      geom_line(linewidth = 1.1) +
      scale_colour_manual(values = c(PAL[3], PAL[4], "#B03A2E")) +
      labs(title = "The hepatocyte life cycle", x = "Time (h)",
           y = "Fraction of the normal liver", colour = NULL) + THEME
  })

  output$damp_plot <- renderPlot({
    d <- run()
    ggplot(d, aes(time)) +
      geom_line(aes(y = DAMP, colour = "DAMP"), linewidth = 1) +
      geom_line(aes(y = KC,   colour = "Kupffer cells (activated)"), linewidth = 1) +
      scale_colour_manual(values = PAL[c(2,5)]) +
      labs(title = "Initiation of sterile inflammation", x = "Time (h)", y = NULL,
           colour = NULL) + THEME
  })

  output$regen_plot <- renderPlot({
    d <- run()
    ggplot(d, aes(time)) +
      geom_line(aes(y = TNF,  colour = "TNF-α (regenerative priming + injury)"),
                linewidth = 1) +
      geom_line(aes(y = IL10, colour = "IL-10 (anti-inflammatory)"), linewidth = 1) +
      scale_colour_manual(values = c("#D62728", PAL[3])) +
      labs(title = "The double edge of TNF: it both amplifies injury and signals regeneration",
           x = "Time (h)", y = "au", colour = NULL) + THEME
  })

  # ----------------------------------------------------------------- tab 6
  output$lab_enz <- renderPlot({
    d <- run()
    ggplot(d, aes(time)) +
      geom_hline(yintercept = 120, linetype = 3, colour = "#B71C1C") +
      geom_line(aes(y = ALT, colour = "ALT"), linewidth = 1.1) +
      geom_line(aes(y = AST, colour = "AST"), linewidth = 1) +
      geom_line(aes(y = ALP, colour = "ALP"), linewidth = 1) +
      geom_line(aes(y = MIR*10, colour = "miR-122 × 10"), linewidth = 1) +
      scale_y_log10() +
      scale_colour_manual(values = PAL[c(2,4,5,1)]) +
      labs(title = "Enzyme and mechanistic biomarkers (dotted = ALT 3×ULN)",
           x = "Time (h)", y = "U/L or fold (log)", colour = NULL) + THEME
  })

  output$lab_fn <- renderPlot({
    d <- run()
    ggplot(d, aes(time)) +
      geom_hline(yintercept = 2.4, linetype = 3, colour = "#B71C1C") +
      geom_line(aes(y = TBIL, colour = "Total bilirubin (mg/dL)"), linewidth = 1.1) +
      geom_line(aes(y = INR,  colour = "INR"), linewidth = 1.1) +
      scale_colour_manual(values = c("#B8860B", "#7D1128")) +
      labs(title = "Synthetic and excretory function (dotted = TBIL 2×ULN)",
           x = "Time (h)", y = NULL, colour = NULL) + THEME
  })

  output$hylaw_plot <- renderPlot({
    d <- run()
    ggplot(d, aes(ALT, TBIL)) +
      annotate("rect", xmin = 120, xmax = Inf, ymin = 2.4, ymax = Inf,
               fill = "#B71C1C", alpha = 0.10) +
      annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.5,
               colour = "#B71C1C", fontface = 2, label = "Hy's Law zone") +
      geom_path(aes(colour = time), linewidth = 1.1) +
      geom_vline(xintercept = 120, linetype = 3) +
      geom_hline(yintercept = 2.4, linetype = 3) +
      scale_x_log10() +
      scale_colour_viridis_c(name = "Time (h)") +
      labs(title = "ALT–bilirubin trajectory (rate axis × reserve axis)",
           x = "ALT (U/L, log)", y = "Total bilirubin (mg/dL)") + THEME
  })

  # ----------------------------------------------------------------- tab 7
  output$bile_plot <- renderPlot({
    d <- run()
    ggplot(d, aes(time)) +
      geom_hline(yintercept = 100, linetype = 3, colour = "grey40") +
      geom_line(aes(y = BAH, colour = "Hepatocellular bile acids (µM)"), linewidth = 1.1) +
      geom_line(aes(y = BAP, colour = "Serum total bile acids (µM)"), linewidth = 1.1) +
      scale_colour_manual(values = c("#B8860B", PAL[1])) +
      labs(title = "Bile acid homeostasis (dotted = the 100 µM cytotoxic threshold)",
           x = "Time (h)", y = "µM", colour = NULL) + THEME
  })

  output$rratio_plot <- renderPlot({
    d <- subset(run(), ALT > 0 & ALP > 0)
    ggplot(d, aes(time, RRATIO)) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = 5, ymax = Inf,
               fill = PAL[2], alpha = 0.08) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = 2,
               fill = PAL[4], alpha = 0.10) +
      geom_line(linewidth = 1.1, colour = "grey20") +
      scale_y_log10() +
      labs(title = "R ratio (above = hepatocellular, below = cholestatic) — a computed result",
           x = "Time (h)", y = "R = (ALT/ULN)/(ALP/ULN), log") + THEME
  })

  output$innate_plot <- renderPlot({
    d <- run()
    ggplot(d, aes(time)) +
      geom_line(aes(y = KC,   colour = "Kupffer cells"), linewidth = 1) +
      geom_line(aes(y = TNF,  colour = "TNF-α"), linewidth = 1) +
      geom_line(aes(y = IL10, colour = "IL-10"), linewidth = 1) +
      scale_colour_manual(values = PAL[c(5,2,3)]) +
      labs(title = "Innate immunity", x = "Time (h)", y = "au", colour = NULL) + THEME
  })

  output$adaptive_plot <- renderPlot({
    d <- run()
    ggplot(d, aes(time)) +
      geom_line(aes(y = TCELL, colour = "Drug-specific effector T cells"), linewidth = 1.1) +
      geom_line(aes(y = TREG,  colour = "Tolerance (Treg/PD-1)"), linewidth = 1.1) +
      scale_colour_manual(values = c("#B03A2E", PAL[3])) +
      labs(title = "Adaptive immunity — what decides idiosyncratic DILI is tolerance",
           x = "Time (h)", y = "Fraction", colour = NULL) + THEME
  })

  # ----------------------------------------------------------------- tab 8
  nacwin <- eventReactive(input$run_nac, {
    withProgress(message = "Computing the NAC time window...", {
      analysis_nac_window()
    })
  })

  output$nacwin_plot <- renderPlot({
    d <- nacwin()
    d$x <- ifelse(is.na(d$nac_start), 60, d$nac_start)
    ggplot(d, aes(x, 100*lost_mass, colour = host)) +
      geom_line(linewidth = 1.1) + geom_point(size = 2) +
      scale_colour_manual(values = PAL[1:3]) +
      scale_x_continuous(breaks = c(0,8,16,24,32,48,60),
                         labels = c("0","8","16","24","32","48","Untreated")) +
      labs(title = "Hepatocyte loss as a function of NAC start time — the cliff moves",
           x = "NAC start time (h)", y = "Hepatocyte mass lost (%)",
           colour = NULL) + THEME
  })

  output$nacwin_tab <- renderDT({
    datatable(nacwin(), options = list(pageLength = 15, scrollX = TRUE)) %>%
      formatRound(columns = which(sapply(nacwin(), is.numeric)), digits = 3)
  })

  # ----------------------------------------------------------------- tab 9
  scen_res <- eventReactive(input$run_scen, {
    sc <- scenarios()
    pick <- input$scen_pick
    withProgress(message = "Running scenarios...", {
      do.call(rbind, lapply(pick, function(k) {
        r <- summarise_run(do.call(sim_dili, sc[[k]]$args))
        cbind(id = k, label = sc[[k]]$label, r)
      }))
    })
  })

  output$scen_plot <- renderPlot({
    d <- scen_res()
    df <- d %>%
      select(label, `Peak ALT` = peak_ALT, `Peak TBIL x1000` = peak_TBIL,
             `Lost mass %` = lost_mass) %>%
      mutate(`Peak TBIL x1000` = `Peak TBIL x1000`*1000,
             `Lost mass %` = `Lost mass %`*100) %>%
      pivot_longer(-label)
    ggplot(df, aes(reorder(label, value), value, fill = name)) +
      geom_col(position = "dodge") + coord_flip() +
      scale_fill_manual(values = PAL[c(2,4,1)]) +
      labs(x = NULL, y = NULL, fill = NULL) + THEME
  })

  output$scen_tab <- renderDT({
    datatable(scen_res(), options = list(pageLength = 12, scrollX = TRUE)) %>%
      formatRound(columns = which(sapply(scen_res(), is.numeric)), digits = 3)
  })

  # ---------------------------------------------------------------- tab 10
  rate_res <- eventReactive(input$run_rate, {
    withProgress(message = "Running the rate experiment...", {
      analysis_rate_vs_dose(dose = input$dose)
    })
  })

  output$rate_plot <- renderPlot({
    d <- rate_res()
    ggplot(d, aes(factor(spread_h), 100*lost_mass)) +
      geom_col(fill = "#B03A2E") +
      geom_text(aes(label = sprintf("%.1f%%", 100*lost_mass)), vjust = -0.4) +
      labs(title = sprintf("A total of %g mg/kg given at different rates",
                           input$dose),
           x = "Total ingestion time (h; 0 = single bolus)",
           y = "Hepatocyte mass lost (%)") + THEME
  })

  output$rate_tab <- renderDT({
    datatable(rate_res(), options = list(pageLength = 10, scrollX = TRUE)) %>%
      formatRound(columns = which(sapply(rate_res(), is.numeric)), digits = 3)
  })
}

shinyApp(ui, server)
