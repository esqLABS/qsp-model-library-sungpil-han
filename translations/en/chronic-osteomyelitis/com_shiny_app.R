## ============================================================================
##  com_shiny_app.R
##  Chronic osteomyelitis · implant-associated bone infection QSP model — interactive dashboard
##
##  To run:
##      shiny::runApp("com_shiny_app.R")
##  or
##      Rscript -e 'shiny::runApp("com_shiny_app.R", port = 8080)'
##
##  This app is designed to answer one question:
##
##      "In this patient, can this regimen sterilise, or only suppress?"
##
##  The answer is not a figure obtained by moving the sliders; it is the single number
##  RSTER kept permanently at the top left of the screen. If RSTER <= 1, no amount of
##  added duration brings sterilisation — tab ⑩ shows that experimentally.
##
##  Tab layout (10 tabs)
##    ①  Patient · intervention summary   — what changes what, on one screen
##    ②  The two penalties                — the thesis of this model (spatial vs phenotypic penalty)
##    ③  Target-site pharmacokinetics     — plasma → bone interstitium → sequestrum → intracellular
##    ④  Bacterial subpopulations         — log10 trajectories of the 5 pools
##    ⑤  Killing kinetics · RSTER         — drug · immunity · growth in competition
##    ⑥  Bone · perfusion (vicious cycle) — the PK barrier the disease builds for itself
##    ⑦  Clinical endpoints · biomarkers  — CRP vs ESR, pain, reservoir, relapse probability
##    ⑧  Toxicity · safety                — efficacy that needs duration vs toxicity that duration creates
##    ⑨  Scenario comparison              — the 11 standard scenarios side by side
##    ⑩  Duration sweep                   — when duration helps and when it is meaningless
##
##  For education and research. Do not use for clinical decision-making.
## ============================================================================

library(shiny)
source("com_mrgsolve_model.R")

PAL <- c(plank = "#2563eb", bfa = "#f97316", bfp = "#b91c1c",
         ic   = "#7c3aed", res = "#0f766e", tot = "#111827",
         drugA = "#2563eb", drugB = "#16a34a", loc = "#0891b2",
         bone = "#78716c", seq = "#57534e", perf = "#db2777",
         crp  = "#ca8a04", esr = "#a16207", pain = "#9333ea",
         good = "#16a34a", bad = "#dc2626", grey = "#94a3b8")

panelplot <- function(...) {
  par(mar = c(4.0, 4.4, 2.4, 1.0), mgp = c(2.5, 0.7, 0), las = 1,
      cex.axis = 0.85, cex.lab = 0.95, cex.main = 1.0, bty = "l")
}
vlines <- function(tsurg, tstop) {
  if (is.finite(tsurg) && tsurg < 1e5)
    abline(v = tsurg, col = PAL[["good"]], lty = 2, lwd = 1.4)
  if (is.finite(tstop) && tstop < 1e5)
    abline(v = tstop, col = PAL[["grey"]], lty = 3, lwd = 1.4)
}

## ============================================================================
##  UI
## ============================================================================

ui <- fluidPage(
  titlePanel("Chronic osteomyelitis · implant-associated bone infection — QSP simulator"),
  tags$p(style = "color:#475569;margin-top:-8px;",
    tags$b("Thesis: "),
    "Treatment failure arises from the product of two independent penalties — ",
    tags$b("a spatial penalty"), " (diffusion into unperfused sequestrum: eta = tanh(phi)/phi, which dose can push through) and ",
    tags$b("a phenotypic penalty"), " (growth-rate-dependent killing: as mu→0, kill→0, which no dose can push through). ",
    "So ", tags$b("debridement sets the initial condition"), ", ", tags$b("rifampicin sets the killing ceiling"), ", and ",
    tags$b("duration sets only the time to arrive"), "."),
  hr(),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("① Focus · host"),
      sliderInput("seq0", "Initial sequestrum volume SEQ (cm3)", 0.2, 20, 2.0, 0.2),
      sliderInput("aimp0", "Implant surface area (cm2) — 0 = no foreign body", 0, 80, 30, 5),
      sliderInput("perf0", "Initial perfused fraction PERF", 0.05, 1.0, 0.37, 0.01),
      sliderInput("padsev", "PAD severity (diabetic foot)", 0, 0.9, 0, 0.05),
      sliderInput("hostf", "Host immune competence (1 = normal)", 0.3, 1.2, 1.0, 0.05),
      sliderInput("egfr0", "eGFR (mL/min/1.73m2)", 20, 120, 90, 5),

      h4("② Antibiotics"),
      selectInput("drugA", "Backbone (growth-dependent)", names(drug_A_library), "vancomycin"),
      numericInput("doseA", "Backbone dose per administration (mg)", 1250, 0, 6000, 50),
      selectInput("iiA", "Backbone dosing interval", c("q8h" = 1/3, "q12h" = 1/2, "q24h" = 1), 1/2),
      checkboxInput("ivA", "Intravenous (uncheck for oral)", TRUE),
      selectInput("drugB", "Partner (growth-independent)", names(drug_B_library), "rifampicin"),
      numericInput("doseB", "Partner dose per administration (mg)", 450, 0, 1200, 50),
      selectInput("iiB", "Partner dosing interval", c("q12h" = 1/2, "q24h" = 1), 1/2),
      sliderInput("weeks", "Treatment duration (weeks)", 1, 52, 6, 1),

      h4("③ Surgery · local"),
      checkboxInput("surg", "Debridement performed", TRUE),
      sliderInput("tsurg", "Time of debridement (days)", 1, 180, 2, 1),
      sliderInput("deblog", "log10 bacterial reduction achieved by debridement (P0)", 0, 7, 5.5, 0.5),
      sliderInput("seqf", "Fraction of sequestrum resected", 0, 0.99, 0.95, 0.01),
      checkboxInput("imprem", "Implant removal", TRUE),
      sliderInput("flap", "Perfusion recovery from flap / muscle graft", 0, 0.35, 0.20, 0.05),
      checkboxInput("revasc", "Revascularisation (PAD corrected)", FALSE),
      numericInput("aloc", "Local antibiotic depot (mg, 0 = none)", 0, 0, 6000, 250),

      h4("④ Simulation"),
      sliderInput("tend", "Observation period (days)", 100, 900, 400, 50),
      actionButton("go", "Recalculate", class = "btn-primary")
    ),

    mainPanel(
      width = 9,
      htmlOutput("verdict"),
      tabsetPanel(
        id = "tabs", type = "tabs",
        tabPanel("① Patient · intervention summary",
                 br(), tableOutput("summ"),
                 plotOutput("ovw", height = "430px"),
                 htmlOutput("interp")),
        tabPanel("② The two penalties",
                 br(), htmlOutput("penexp"),
                 plotOutput("pen", height = "560px")),
        tabPanel("③ Target-site pharmacokinetics", br(), plotOutput("pk", height = "620px"),
                 htmlOutput("pkexp")),
        tabPanel("④ Bacterial subpopulations", br(), plotOutput("bact", height = "560px"),
                 htmlOutput("bactexp")),
        tabPanel("⑤ Killing kinetics · RSTER", br(), plotOutput("kill", height = "600px"),
                 htmlOutput("killexp")),
        tabPanel("⑥ Bone · perfusion (vicious cycle)", br(), plotOutput("bone", height = "620px"),
                 htmlOutput("boneexp")),
        tabPanel("⑦ Clinical endpoints", br(), plotOutput("clin", height = "600px"),
                 htmlOutput("clinexp")),
        tabPanel("⑧ Toxicity", br(), plotOutput("tox", height = "560px"),
                 htmlOutput("toxexp")),
        tabPanel("⑨ Scenario comparison", br(),
                 actionButton("runscn", "Run the 11 standard scenarios (tens of seconds)",
                              class = "btn-warning"),
                 br(), br(), tableOutput("scntab"), plotOutput("scnplot", height = "440px")),
        tabPanel("⑩ Duration sweep", br(),
                 actionButton("runsweep", "Run the duration sweep (2·4·6·12·26 weeks × with/without rifampicin)",
                              class = "btn-warning"),
                 br(), br(), tableOutput("swtab"), plotOutput("swplot", height = "420px"),
                 htmlOutput("swexp"))
      )
    )
  )
)

## ============================================================================
##  SERVER
## ============================================================================

server <- function(input, output, session) {

  built <- eventReactive(input$go, ignoreNULL = FALSE, {
    days <- input$weeks * 7
    m <- set_drugs(mod, input$drugA, input$drugB)
    tsurg <- if (input$surg) input$tsurg else 1e6

    m <- param(m,
      HOSTF   = input$hostf,
      PADSEV  = input$padsev,
      EGFR0   = input$egfr0,
      TSURG   = tsurg,
      SURGDUR = 0.25,
      DEBLOG  = if (input$surg) input$deblog else 0,
      DEBSEQF = if (input$surg) input$seqf   else 0,
      DEBPUSF = if (input$surg) 0.95 else 0,
      DEBICF  = if (input$surg) min(0.95, input$deblog / 6) else 0,
      IMPREM  = if (input$surg && input$imprem) 1 else 0,
      FLAPB   = if (input$surg) input$flap else 0,
      TREVASC = if (input$revasc) max(1, tsurg + 1) else 1e6,
      REVASCB = if (input$revasc) 0.75 else 0,
      TABX0   = 0,
      TABX1   = days
    )
    m <- update(m, init = list(SEQ = input$seq0, AIMP = input$aimp0,
                               PERF = input$perf0, EGFRC = input$egfr0))

    evs <- list()
    if (input$doseA > 0) {
      iiA <- as.numeric(input$iiA)
      evs[[length(evs) + 1]] <- if (input$ivA)
        ev_iv(input$doseA, iiA, days, "CENA", inf_h = 1.5)
      else ev_po(input$doseA, iiA, days, "GUTA")
    }
    if (input$doseB > 0 && input$drugB != "none")
      evs[[length(evs) + 1]] <- ev_po(input$doseB, as.numeric(input$iiB), days, "GUTB")
    ## The local depot is inserted at debridement, but in a scenario without
    ## debridement it goes in on day 0. tsurg is 1e6 when there is no debridement, so
    ## using that value directly makes the solver try to integrate to day one million and stall.
    if (input$aloc > 0) {
      tloc <- if (input$surg) min(input$tsurg, input$tend) else 0
      evs[[length(evs) + 1]] <- ev_bolus(input$aloc, "ALOC", time = tloc)
    }

    e <- if (length(evs)) do.call(ev_bind, evs) else ev_none()
    out <- as.data.frame(run_scn(m, e, end = input$tend, delta = 0.25))
    list(sim = out, tsurg = tsurg, tstop = days, mod = m)
  })

  idx_at <- function(d, t) which.min(abs(d$time - t))

  ## ---- always-visible verdict -------------------------------------------
  output$verdict <- renderUI({
    b <- built(); d <- b$sim
    i7 <- idx_at(d, min(7, b$tstop))
    itx <- idx_at(d, b$tstop)
    rs <- d$RSTER[i7]; res <- d$RESV[itx]; pr <- d$PRLPS[itx]
    col <- if (rs > 1) PAL[["good"]] else PAL[["bad"]]
    msg <- if (rs > 1)
      "Sterilisation possible (RSTER &gt; 1) — the only question left is duration."
    else
      "Sterilisation impossible (RSTER &le; 1) — this regimen is suppressive. Extending the duration will not bring cure."
    HTML(sprintf(
      "<div style='border-left:6px solid %s;background:#f8fafc;padding:10px 14px;margin-bottom:10px;'>
       <span style='font-size:1.25em;font-weight:700;color:%s'>RSTER (day 7) = %.3f</span>
       &nbsp;&nbsp;<span style='color:%s;font-weight:600'>%s</span><br>
       <span style='color:#475569'>Reservoir at end of treatment = <b>%.3g CFU</b> ·
       relapse probability = <b>%.1f%%</b> · sequestrum penetration efficiency eta (backbone) = <b>%.3f</b> ·
       diffusion distance L = <b>%.3f cm</b> · persister killing = <b>%.3f/d</b> (drug) + <b>%.3f/d</b> (immunity)</span></div>",
      col, col, rs, col, msg, res, 100 * pr,
      d$ETAA[i7], d$LEFFO[i7], d$KILLPS[i7], d$KIMMPS[i7]))
  })

  ## ---- ① summary ---------------------------------------------------------
  output$summ <- renderTable({
    b <- built(); d <- b$sim
    i0 <- 1; i7 <- idx_at(d, min(7, b$tstop)); itx <- idx_at(d, b$tstop); ie <- nrow(d)
    pick <- function(i) c(d$LBTOT[i], d$LBFP[i], d$LBIC[i], d$LBRES[i], d$SEQ[i],
                          d$PERF[i], d$BMLOSS[i], d$CRP[i], d$ESR[i], d$PAIN[i])
    tb <- data.frame(
      metric = c("Metric — total bacteria", "", "", "", "", "", "", "", "", ""),
      d0 = pick(i0), d7 = pick(i7), dtx = pick(itx), dend = pick(ie),
      stringsAsFactors = FALSE)
    tb$metric <- c("Total bacteria log10", "Persisters log10",
                   "Intracellular log10", "rpoB resistance log10",
                   "Sequestrum SEQ (cm3)", "Perfusion PERF",
                   "Cortical bone loss (%)",
                   "CRP (mg/L)", "ESR (mm/h)", "Pain NRS")
    names(tb) <- c("Metric", "0d", "7d",
                   "End of therapy", "End of follow-up")
    tb
  }, digits = 3, striped = TRUE, hover = TRUE, width = "100%")

  output$ovw <- renderPlot({
    b <- built(); d <- b$sim
    par(mfrow = c(1, 3)); panelplot()
    plot(d$time, d$LBTOT, type = "l", lwd = 2.4, col = PAL[["tot"]], ylim = c(-6.5, 10.5),
         xlab = "Days", ylab = "log10 CFU", main = "Total bacterial burden")
    lines(d$time, d$LBFP, lwd = 2, col = PAL[["bfp"]])
    abline(h = 0, col = PAL[["grey"]], lty = 3)
    vlines(b$tsurg, b$tstop)
    legend("topright", c("Total", "Persisters", "1 CFU"), col = c(PAL[["tot"]], PAL[["bfp"]], PAL[["grey"]]),
           lty = c(1, 1, 3), lwd = 2, bty = "n", cex = 0.85)

    plot(d$time, d$RSTER, type = "l", lwd = 2.4, col = PAL[["good"]],
         xlab = "Days", ylab = "RSTER", main = "Sterilisation ratio RSTER")
    abline(h = 1, col = PAL[["bad"]], lty = 2, lwd = 1.6)
    vlines(b$tsurg, b$tstop)

    plot(d$time, d$PERF, type = "l", lwd = 2.4, col = PAL[["perf"]], ylim = c(0, 1),
         xlab = "Days", ylab = "PERF / eta", main = "Perfusion and penetration efficiency")
    lines(d$time, d$ETAA, lwd = 2, col = PAL[["drugA"]])
    lines(d$time, d$ETAB, lwd = 2, col = PAL[["drugB"]], lty = 2)
    vlines(b$tsurg, b$tstop)
    legend("bottomright", c("PERF", "eta backbone", "eta partner"),
           col = c(PAL[["perf"]], PAL[["drugA"]], PAL[["drugB"]]),
           lty = c(1, 1, 2), lwd = 2, bty = "n", cex = 0.85)
  })

  output$interp <- renderUI(HTML(
    "<p style='color:#475569'>Vertical dotted lines: <span style='color:#16a34a'>green = debridement</span>,
     <span style='color:#94a3b8'>grey = antibiotics stopped</span>.
     If the total keeps falling after they are stopped, the reservoir really was emptied;
     if it climbs again after they are stopped, it regrew from a surviving reservoir.
     The eta in the right-hand panel shows that <b>when the geometry of the focus changes, the penetration of the same drug changes</b>
     — part of the effect of debridement is not that it reduced the bacteria but that it shortened the diffusion distance.</p>"))

  ## ---- ② the two penalties ----------------------------------------------
  output$penexp <- renderUI(HTML(
    "<p style='color:#475569'><b>On the left (the spatial penalty)</b> the problem is concentration — penetration gets worse
     the further right you go along the curve, but raise the concentration at the <b>rim</b> 100-fold and the deep concentration rises 100-fold too.
     <b>On the right (the phenotypic penalty)</b> the problem is not concentration — even at saturating concentration (setting h=1)
     the height of the persister bar is decided by Emax_gi alone. Only rifampicin has that bar, and
     that is the one reason rifampicin is special in this disease.</p>"))

  output$pen <- renderPlot({
    b <- built(); d <- b$sim
    par(mfrow = c(2, 2)); panelplot()

    ## (a) Thiele penalty curves
    Ls <- seq(0.005, 0.6, length.out = 200)
    pA <- drug_A_library[[input$drugA]]
    eta <- function(L, k, D) { p <- L * sqrt(k / D); ifelse(p < 1e-6, 1, tanh(p) / p) }
    plot(Ls, eta(Ls, pA$AKLS, pA$ADEFF), type = "l", lwd = 2.6, col = PAL[["drugA"]],
         ylim = c(0, 1), xlab = "Diffusion distance L (cm)", ylab = "Penetration efficiency eta",
         main = "(a) Spatial penalty: eta = tanh(phi)/phi")
    pB <- drug_B_library[[input$drugB]]
    lines(Ls, eta(Ls, pB$BKLS, pB$BDEFF), lwd = 2.6, col = PAL[["drugB"]], lty = 2)
    Lnow <- d$LEFFO[idx_at(d, min(7, b$tstop))]
    abline(v = Lnow, col = PAL[["bad"]], lty = 3, lwd = 1.6)
    text(Lnow, 0.95, sprintf(" current L = %.3f", Lnow), col = PAL[["bad"]], cex = 0.85, adj = 0)
    legend("topright", c(input$drugA, input$drugB), col = c(PAL[["drugA"]], PAL[["drugB"]]),
           lty = c(1, 2), lwd = 2.4, bty = "n", cex = 0.85)

    ## (b) Emax collapse across pools, at saturating concentration
    pools <- c("Planktonic\n(mu/mumax=0.35)", "Biofilm rim\n(0.12)", "Persisters\n(0)")
    kA <- c(pA$AEMGD * 0.35 + pA$AEMGI, pA$AEMGD * 0.12 + pA$AEMGI, pA$AEMGI)
    kB <- if (input$drugB == "none") c(0, 0, 0) else
      c(pB$BEMGD * 0.35 + pB$BEMGI, pB$BEMGD * 0.12 + pB$BEMGI, pB$BEMGI)
    mm <- rbind(kA, kB)
    rownames(mm) <- c("A", "B")
    barplot(mm, beside = TRUE, names.arg = pools, col = c(PAL[["drugA"]], PAL[["drugB"]]),
            ylab = "Killing rate at saturating concentration (1/d)",
            main = "(b) Phenotypic penalty: Emax collapse (assuming saturating concentration)")
    legend("topright", c("Backbone", "Partner"), fill = c(PAL[["drugA"]], PAL[["drugB"]]),
           bty = "n", cex = 0.85)
    abline(h = 6 * 0.12, col = PAL[["bad"]], lty = 2)
    text(1, 6 * 0.12, " biofilm growth rate", col = PAL[["bad"]], adj = 0, cex = 0.8)

    ## (c) the two penalties multiplied, over time
    plot(d$time, d$ETAA, type = "l", lwd = 2.4, col = PAL[["drugA"]], ylim = c(0, 1),
         xlab = "Days", ylab = "Dimensionless", main = "(c) The two penalties over time")
    lines(d$time, d$MUBF / 6, lwd = 2.4, col = PAL[["bfa"]])
    lines(d$time, d$ETAA * d$MUBF / 6, lwd = 3, col = PAL[["bad"]])
    vlines(b$tsurg, b$tstop)
    legend("topright", c("eta (spatial)", "mu/mumax (phenotypic)", "product (effective killing capacity)"),
           col = c(PAL[["drugA"]], PAL[["bfa"]], PAL[["bad"]]), lwd = 2.4, bty = "n", cex = 0.8)

    ## (d) T_ster map
    kwake <- as.numeric(param(b$mod)$KWAKE)
    kills <- seq(0, 1.5, length.out = 200)
    plot(kills, log(1e8) / (kills + kwake), type = "l", lwd = 2.4, col = PAL[["bfp"]],
         log = "y", ylim = c(1, 400), xlab = "Persister killing rate KILLPS (1/d)",
         ylab = "T_ster (days, log axis)", main = "(d) T_ster = ln(P0)/(KILLPS + KWAKE)")
    lines(kills, log(1e5) / (kills + kwake), lwd = 2.4, col = PAL[["bfa"]])
    lines(kills, log(1e2) / (kills + kwake), lwd = 2.4, col = PAL[["good"]])
    abline(h = c(42, 84), col = PAL[["grey"]], lty = 3)
    text(1.2, 45, "6 weeks", col = PAL[["grey"]], cex = 0.8)
    text(1.2, 90, "12 weeks", col = PAL[["grey"]], cex = 0.8)
    kp <- d$KILLPS[idx_at(d, min(7, b$tstop))] + d$KIMMPS[idx_at(d, min(7, b$tstop))]
    abline(v = kp, col = "#111827", lty = 2, lwd = 1.6)
    legend("topright", c("P0 = 1e8 (no debridement)", "P0 = 1e5", "P0 = 1e2 (radical debridement)"),
           col = c(PAL[["bfp"]], PAL[["bfa"]], PAL[["good"]]), lwd = 2.4, bty = "n", cex = 0.8)
  })

  ## ---- ③ target-site PK --------------------------------------------------
  output$pk <- renderPlot({
    b <- built(); d <- b$sim
    par(mfrow = c(2, 2)); panelplot()
    w <- d$time <= min(max(d$time), b$tstop + 14)

    plot(d$time[w], d$CAOBS[w], type = "l", lwd = 1.6, col = PAL[["drugA"]],
         xlab = "Days", ylab = "mg/L", main = "(a) Plasma concentrations")
    lines(d$time[w], d$CBOBS[w], lwd = 1.6, col = PAL[["drugB"]])
    legend("topright", c("Backbone", "Partner"), col = c(PAL[["drugA"]], PAL[["drugB"]]),
           lwd = 2, bty = "n", cex = 0.85)

    ymax <- max(1e-3, d$CBA[w], d$CSA[w], d$CIA[w], na.rm = TRUE)
    plot(d$time[w], pmax(d$CBA[w], 1e-4), type = "l", lwd = 2, col = PAL[["drugA"]],
         log = "y", ylim = c(max(1e-4, ymax / 1e4), ymax * 1.3),
         xlab = "Days", ylab = "mg/L (log axis)", main = "(b) Backbone: concentration by compartment")
    lines(d$time[w], pmax(d$CSA[w], 1e-4), lwd = 2, col = PAL[["seq"]])
    lines(d$time[w], pmax(d$CIA[w], 1e-4), lwd = 2, col = PAL[["ic"]])
    if (max(d$CLOC) > 0) lines(d$time[w], pmax(d$CLOC[w], 1e-4), lwd = 2, col = PAL[["loc"]], lty = 2)
    legend("topright", c("Bone interstitium", "Sequestrum/biofilm", "Intracellular", "Local wound fluid"),
           col = c(PAL[["drugA"]], PAL[["seq"]], PAL[["ic"]], PAL[["loc"]]),
           lty = c(1, 1, 1, 2), lwd = 2, bty = "n", cex = 0.8)

    ymax2 <- max(1e-4, d$CBB[w], d$CSB[w], d$CIB[w], na.rm = TRUE)
    plot(d$time[w], pmax(d$CBB[w], 1e-5), type = "l", lwd = 2, col = PAL[["drugB"]],
         log = "y", ylim = c(max(1e-5, ymax2 / 1e4), ymax2 * 1.3),
         xlab = "Days", ylab = "mg/L (log axis)", main = "(c) Partner: concentration by compartment")
    lines(d$time[w], pmax(d$CSB[w], 1e-5), lwd = 2, col = PAL[["seq"]])
    lines(d$time[w], pmax(d$CIB[w], 1e-5), lwd = 2, col = PAL[["ic"]])
    legend("topright", c("Bone interstitium", "Sequestrum/biofilm", "Intracellular"),
           col = c(PAL[["drugB"]], PAL[["seq"]], PAL[["ic"]]), lwd = 2, bty = "n", cex = 0.8)

    plot(d$time[w], d$BSRA[w], type = "l", lwd = 2.2, col = PAL[["bone"]],
         ylim = c(0, max(0.5, max(d$BSRA[w], na.rm = TRUE))),
         xlab = "Days", ylab = "Ratio", main = "(d) Apparent bone:serum ratio · induction state")
    lines(d$time[w], d$ETAA[w], lwd = 2, col = PAL[["drugA"]], lty = 2)
    vlines(b$tsurg, b$tstop)
    legend("topright", c("Bone:serum ratio (perfusion-dependent)", "eta (geometry-dependent)"),
           col = c(PAL[["bone"]], PAL[["drugA"]]), lty = c(1, 2), lwd = 2, bty = "n", cex = 0.8)
  })

  output$pkexp <- renderUI(HTML(
    "<p style='color:#475569'>The gaps between the three lines in (b)/(c) are the pharmacokinetics of this disease.
     The concentration is cut multiplicatively once from plasma into the bone interstitium, and again from the bone interstitium into the depth of the sequestrum.
     A local depot (the blue dotted line) bypasses <b>only the first multiplication</b> — the second (eta) remains untouched.
     That is why in S7 of tab ⑨ the biofilm survives even when the wound fluid concentration is tens of times the MIC.</p>"))

  ## ---- ④ bacterial pools --------------------------------------------------
  output$bact <- renderPlot({
    b <- built(); d <- b$sim
    par(mfrow = c(1, 2)); panelplot()
    plot(d$time, d$LBTOT, type = "l", lwd = 3, col = PAL[["tot"]], ylim = c(-6.5, 10.5),
         xlab = "Days", ylab = "log10 CFU", main = "(a) Trajectories by subpopulation")
    lines(d$time, d$LBPL,  lwd = 2, col = PAL[["plank"]])
    lines(d$time, d$LBFA,  lwd = 2, col = PAL[["bfa"]])
    lines(d$time, d$LBFP,  lwd = 2.4, col = PAL[["bfp"]])
    lines(d$time, d$LBIC,  lwd = 2, col = PAL[["ic"]])
    lines(d$time, d$LBRES, lwd = 2, col = PAL[["res"]], lty = 2)
    abline(h = 0, col = PAL[["grey"]], lty = 3)
    vlines(b$tsurg, b$tstop)
    legend("bottomleft", c("Total", "Planktonic", "Biofilm rim", "Persisters", "Intracellular", "rpoB-resistant"),
           col = c(PAL[["tot"]], PAL[["plank"]], PAL[["bfa"]], PAL[["bfp"]], PAL[["ic"]], PAL[["res"]]),
           lty = c(1, 1, 1, 1, 1, 2), lwd = 2, bty = "n", cex = 0.8, ncol = 2)

    plot(d$time, d$RESV, type = "l", lwd = 2.6, col = PAL[["bfp"]], log = "y",
         ylim = c(max(1e-6, min(d$RESV[d$RESV > 0])), max(1e10, max(d$RESV))),
         xlab = "Days", ylab = "Reservoir (CFU, log axis)",
         main = "(b) The raw material of relapse: reservoir and relapse probability")
    abline(h = 1, col = PAL[["bad"]], lty = 2, lwd = 1.6)
    par(new = TRUE)
    plot(d$time, d$PRLPS, type = "l", lwd = 2, col = PAL[["grey"]], axes = FALSE,
         xlab = "", ylab = "", ylim = c(0, 1))
    axis(4, col = PAL[["grey"]], col.axis = PAL[["grey"]])
    mtext("Relapse probability", side = 4, line = 2.2, col = PAL[["grey"]], cex = 0.85, las = 0)
    vlines(b$tsurg, b$tstop)
  })

  output$bactexp <- renderUI(HTML(
    "<p style='color:#475569'>If the persisters (red) come down later than the other lines, that regimen has no
     Emax_gi. In that case the persisters disappear <b>only through awakening (KWAKE)</b>, so the
     duration required is on the scale of hundreds of days — clinically the same statement as 'it does not cure'.
     The red horizontal line (1 CFU) in the right-hand panel is the practical boundary of cure: a deterministic ODE cannot
     represent anything below 1 CFU, so below it the relapse probability 1-exp(-p*N) takes over.</p>"))

  ## ---- ⑤ killing ----------------------------------------------------------
  output$kill <- renderPlot({
    b <- built(); d <- b$sim
    par(mfrow = c(2, 2)); panelplot()

    plot(d$time, d$KILLBF, type = "l", lwd = 2.4, col = PAL[["drugA"]],
         ylim = c(0, max(0.5, max(d$KILLBF, d$MUBF, d$KIMMBF, na.rm = TRUE))),
         xlab = "Days", ylab = "1/d", main = "(a) Biofilm rim: drug vs immunity vs growth")
    lines(d$time, d$KIMMBF, lwd = 2.4, col = PAL[["ic"]])
    lines(d$time, d$MUBF, lwd = 2.4, col = PAL[["bad"]], lty = 2)
    vlines(b$tsurg, b$tstop)
    legend("topright", c("Drug killing", "Immune killing", "Growth (the line to beat)"),
           col = c(PAL[["drugA"]], PAL[["ic"]], PAL[["bad"]]),
           lty = c(1, 1, 2), lwd = 2.2, bty = "n", cex = 0.8)

    plot(d$time, d$KILLPS, type = "l", lwd = 2.4, col = PAL[["drugA"]],
         ylim = c(0, max(0.2, max(d$KILLPS + d$KIMMPS, na.rm = TRUE))),
         xlab = "Days", ylab = "1/d", main = "(b) Persisters: what actually kills them")
    lines(d$time, d$KIMMPS, lwd = 2.4, col = PAL[["ic"]])
    abline(h = as.numeric(param(b$mod)$KWAKE), col = PAL[["grey"]], lty = 3, lwd = 1.6)
    vlines(b$tsurg, b$tstop)
    legend("topright", c("Drug", "Immunity (incl. clearance of infected cells)", "Awakening rate KWAKE"),
           col = c(PAL[["drugA"]], PAL[["ic"]], PAL[["grey"]]),
           lty = c(1, 1, 3), lwd = 2.2, bty = "n", cex = 0.8)

    plot(d$time, d$RSTER, type = "l", lwd = 2.8, col = PAL[["good"]],
         xlab = "Days", ylab = "RSTER", main = "(c) Sterilisation ratio — whether it exceeds 1 is everything")
    abline(h = 1, col = PAL[["bad"]], lty = 2, lwd = 1.8)
    vlines(b$tsurg, b$tstop)

    ts <- pmin(d$TSTER, 600)
    plot(d$time, ts, type = "l", lwd = 2.4, col = PAL[["bfp"]],
         xlab = "Days", ylab = "T_ster (days)", main = "(d) Time still needed for sterilisation")
    abline(h = c(42, 84), col = PAL[["grey"]], lty = 3)
    vlines(b$tsurg, b$tstop)
  })

  output$killexp <- renderUI(HTML(
    "<p style='color:#475569'>In (a) the blue line plus the purple line must rise above the red dotted line for the biofilm to shrink.
     (b) is the most clinical figure in this model: it separates whether what kills the persisters is <b>the drug</b>,
     <b>immunity</b>, or <b>merely waiting (KWAKE)</b>. In a regimen without rifampicin the
     blue line lies below the grey dotted line, and cure is then decided in effect by the single P0 that debridement left behind.</p>"))

  ## ---- ⑥ bone & perfusion -------------------------------------------------
  output$bone <- renderPlot({
    b <- built(); d <- b$sim
    par(mfrow = c(2, 2)); panelplot()

    plot(d$time, d$BM, type = "l", lwd = 2.4, col = PAL[["bone"]], ylim = c(0, 105),
         xlab = "Days", ylab = "Relative value", main = "(a) Viable cortical bone and sequestrum")
    par(new = TRUE)
    plot(d$time, d$SEQ, type = "l", lwd = 2.4, col = PAL[["seq"]], axes = FALSE,
         xlab = "", ylab = "")
    axis(4); mtext("SEQ (cm3)", side = 4, line = 2.2, cex = 0.85, las = 0)
    vlines(b$tsurg, b$tstop)
    legend("right", c("Viable bone (left axis)", "Sequestrum (right axis)"),
           col = c(PAL[["bone"]], PAL[["seq"]]), lwd = 2.2, bty = "n", cex = 0.8)

    plot(d$time, d$PERF, type = "l", lwd = 2.6, col = PAL[["perf"]], ylim = c(0, 1),
         xlab = "Days", ylab = "PERF", main = "(b) Perfusion — the common gate for drug and immunity")
    abline(h = as.numeric(param(b$mod)$PERFCRIT), col = PAL[["bad"]], lty = 2)
    text(max(d$time) * 0.6, as.numeric(param(b$mod)$PERFCRIT) + 0.03,
         "below this the bone dies", col = PAL[["bad"]], cex = 0.8)
    vlines(b$tsurg, b$tstop)

    plot(d$time, d$OC, type = "l", lwd = 2.2, col = PAL[["bad"]],
         xlab = "Days", ylab = "Relative value", main = "(c) Osteoclasts · osteoblasts · abscess",
         ylim = c(0, max(d$OC, d$OB, d$PUS, na.rm = TRUE)))
    lines(d$time, d$OB, lwd = 2.2, col = PAL[["good"]])
    lines(d$time, d$PUS, lwd = 2.2, col = PAL[["ic"]])
    vlines(b$tsurg, b$tstop)
    legend("topright", c("Osteoclasts", "Osteoblasts", "Abscess (mL)"),
           col = c(PAL[["bad"]], PAL[["good"]], PAL[["ic"]]), lwd = 2.2, bty = "n", cex = 0.8)

    plot(d$SEQ, d$ETAA, type = "l", lwd = 2.6, col = PAL[["drugA"]],
         xlab = "Sequestrum SEQ (cm3)", ylab = "Penetration efficiency eta",
         main = "(d) A section through the vicious cycle: as SEQ grows the drug cannot get in")
    points(d$SEQ[1], d$ETAA[1], pch = 19, col = PAL[["bad"]], cex = 1.3)
    points(d$SEQ[nrow(d)], d$ETAA[nrow(d)], pch = 17, col = PAL[["good"]], cex = 1.3)
    legend("topright", c("Start", "End"), pch = c(19, 17),
           col = c(PAL[["bad"]], PAL[["good"]]), bty = "n", cex = 0.85)
  })

  output$boneexp <- renderUI(HTML(
    "<p style='color:#475569'>(d) puts the central loop of this model on a single panel.
     The infection kills bone (the sequestrum grows), the dead bone lengthens the diffusion distance and lowers the drug's penetration,
     and the lowered penetration sustains the infection. Debridement is the only intervention that moves the point
     <b>up and to the left</b> along this curve, and where that point sits changes the effect of the same drug several-fold.</p>"))

  ## ---- ⑦ clinical --------------------------------------------------------
  output$clin <- renderPlot({
    b <- built(); d <- b$sim
    par(mfrow = c(2, 2)); panelplot()

    plot(d$time, d$CRP, type = "l", lwd = 2.6, col = PAL[["crp"]],
         ylim = c(0, max(d$CRP, d$ESR, na.rm = TRUE) * 1.05),
         xlab = "Days", ylab = "CRP (mg/L) / ESR (mm/h)",
         main = "(a) Same event, different clocks: CRP 19h vs ESR 10d")
    lines(d$time, d$ESR, lwd = 2.6, col = PAL[["esr"]], lty = 2)
    vlines(b$tsurg, b$tstop)
    legend("topright", c("CRP", "ESR"), col = c(PAL[["crp"]], PAL[["esr"]]),
           lty = c(1, 2), lwd = 2.4, bty = "n", cex = 0.85)

    plot(d$time, d$CRP / max(d$CRP), type = "l", lwd = 2.4, col = PAL[["crp"]],
         ylim = c(0, 1.05), xlab = "Days", ylab = "Normalised (0-1)",
         main = "(b) The biomarkers cannot see the reservoir")
    lines(d$time, d$RESV / max(d$RESV), lwd = 2.6, col = PAL[["bfp"]])
    vlines(b$tsurg, b$tstop)
    legend("topright", c("CRP (normalised)", "Reservoir (normalised)"),
           col = c(PAL[["crp"]], PAL[["bfp"]]), lwd = 2.4, bty = "n", cex = 0.85)

    plot(d$time, d$PAIN, type = "l", lwd = 2.6, col = PAL[["pain"]], ylim = c(0, 10),
         xlab = "Days", ylab = "NRS", main = "(c) Pain")
    vlines(b$tsurg, b$tstop)

    plot(d$time, d$PRLPS, type = "l", lwd = 2.8, col = PAL[["bad"]], ylim = c(0, 1),
         xlab = "Days", ylab = "Probability", main = "(d) Relapse probability from the reservoir")
    vlines(b$tsurg, b$tstop)
  })

  output$clinexp <- renderUI(HTML(
    "<p style='color:#475569'>(b) is where clinical practice is most often deceived. CRP follows the planktonic bacteria
     and the cytokines, so <b>it can normalise while the biofilm is still entirely intact</b>.
     'The inflammatory markers have come down, so we can stop' is dangerous exactly because these two curves diverge.
     ESR errs in the opposite direction — it stays high for weeks even after the infection has actually resolved.</p>"))

  ## ---- ⑧ toxicity --------------------------------------------------------
  output$tox <- renderPlot({
    b <- built(); d <- b$sim
    par(mfrow = c(2, 2)); panelplot()
    plot(d$time, d$EGFRC, type = "l", lwd = 2.4, col = PAL[["bad"]],
         xlab = "Days", ylab = "mL/min/1.73m2", main = "(a) eGFR (vancomycin AKI)")
    vlines(b$tsurg, b$tstop)
    plot(d$time, d$AUC24A, type = "l", lwd = 2.4, col = PAL[["drugA"]],
         xlab = "Days", ylab = "mg*h/L", main = "(b) Backbone AUC24 (target band 400-600)")
    abline(h = c(400, 600), col = PAL[["grey"]], lty = c(3, 2))
    vlines(b$tsurg, b$tstop)
    plot(d$time, d$PLT, type = "l", lwd = 2.4, col = PAL[["ic"]],
         xlab = "Days", ylab = "10^9/L", main = "(c) Platelets (linezolid, duration-dependent)")
    abline(h = 100, col = PAL[["bad"]], lty = 2); abline(v = 28, col = PAL[["grey"]], lty = 3)
    vlines(b$tsurg, b$tstop)
    plot(d$time, d$ALT, type = "l", lwd = 2.4, col = PAL[["esr"]],
         xlab = "Days", ylab = "U/L", main = "(d) ALT (rifampicin)")
    abline(h = 3 * 25, col = PAL[["bad"]], lty = 2)
    vlines(b$tsurg, b$tstop)
  })

  output$toxexp <- renderUI(HTML(
    "<p style='color:#475569'>Efficacy <b>demands</b> duration and toxicity is <b>created by</b> duration.
     So the optimal duration is T* = argmax [P_cure(T) - sum w_k Tox_k(T)], and
     debridement pushes the P_cure(T) curve to the left and so shortens T* — surgery is itself a reduction in toxicity.
     Vancomycin erodes renal function, and that falling renal function raises exposure again, a feedback ((a)+(b)).</p>"))

  ## ---- ⑨ scenarios -------------------------------------------------------
  scn_res <- eventReactive(input$runscn, {
    scns <- list(scn01(), scn02(), scn03(), scn04(), scn04b(), scn05(),
                 scn06(FALSE), scn06(TRUE), scn07(), scn08(FALSE), scn08(TRUE),
                 scn10(), scn11(2), scn11(120))
    lapply(scns, summarise_scn, verbose = FALSE)
  })

  output$scntab <- renderTable({
    r <- scn_res()
    tb <- do.call(rbind, lapply(r, function(x) x$summary))
    tb[, c("scenario", "RSTER_d7", "KILLBF_d7", "KILLPS_d7", "ETAA_d7",
           "logB_tx", "logPersist", "logRes_tx", "reservoir", "P_relapse",
           "logB_end", "CRP_tx", "SEQ_end", "PERF_end")]
  }, digits = 3, striped = TRUE, hover = TRUE, width = "100%")

  output$scnplot <- renderPlot({
    r <- scn_res()
    panelplot(); par(mfrow = c(1, 2))
    tb <- do.call(rbind, lapply(r, function(x) x$summary))
    cols <- ifelse(tb$P_relapse < 0.5, PAL[["good"]], PAL[["bad"]])
    barplot(tb$RSTER_d7, horiz = TRUE, col = cols, border = NA,
            names.arg = sub(" .*", "", tb$scenario), las = 1, cex.names = 0.7,
            xlab = "RSTER (day 7)", main = "Sterilisation potential by regimen")
    abline(v = 1, col = "#111827", lty = 2, lwd = 1.8)
    plot(tb$RSTER_d7, pmax(tb$reservoir, 1e-8), log = "y", pch = 19, col = cols, cex = 1.5,
         xlab = "RSTER (day 7)", ylab = "Reservoir at end of treatment (CFU, log axis)",
         main = "Whether RSTER exceeds 1 is what separates the outcomes")
    abline(v = 1, col = "#111827", lty = 2); abline(h = 1, col = PAL[["bad"]], lty = 2)
    text(tb$RSTER_d7, pmax(tb$reservoir, 1e-8), sub(" .*", "", tb$scenario),
         pos = 4, cex = 0.65, col = "#475569")
  })

  ## ---- ⑩ duration sweep --------------------------------------------------
  sw_res <- eventReactive(input$runsweep, {
    out <- list()
    for (w in c(2, 4, 6, 12, 26)) for (rif in c(TRUE, FALSE)) {
      s <- scn09(w, with_rif = rif, deblog = 3.0, impkeep = TRUE)
      x <- summarise_scn(s, verbose = FALSE)$summary
      x$weeks <- w; x$rif <- rif
      out[[length(out) + 1]] <- x
    }
    do.call(rbind, out)
  })

  output$swtab <- renderTable({
    sw <- sw_res()
    sw$regimen <- ifelse(sw$rif, "Levofloxacin + rifampicin", "Levofloxacin alone")
    sw[, c("regimen", "weeks", "RSTER_d7", "KILLPS_d7", "logB_tx", "logPersist",
           "reservoir", "P_relapse", "logB_end")]
  }, digits = 3, striped = TRUE, hover = TRUE, width = "100%")

  output$swplot <- renderPlot({
    sw <- sw_res(); panelplot(); par(mfrow = c(1, 2))
    a <- sw[sw$rif, ]; bq <- sw[!sw$rif, ]
    plot(a$weeks, pmax(a$reservoir, 1e-8), type = "b", log = "y", pch = 19, lwd = 2.4,
         col = PAL[["good"]], ylim = c(1e-8, 1e10),
         xlab = "Treatment duration (weeks)", ylab = "Reservoir at end of treatment (CFU, log axis)",
         main = "When does duration help")
    lines(bq$weeks, pmax(bq$reservoir, 1e-8), type = "b", pch = 17, lwd = 2.4, col = PAL[["bad"]])
    abline(h = 1, col = "#111827", lty = 2, lwd = 1.6)
    legend("topright", c("With rifampicin (RSTER > 1)", "Alone (RSTER < 1)"),
           col = c(PAL[["good"]], PAL[["bad"]]), pch = c(19, 17), lwd = 2.2, bty = "n", cex = 0.85)
    plot(a$weeks, a$P_relapse, type = "b", pch = 19, lwd = 2.4, col = PAL[["good"]],
         ylim = c(0, 1), xlab = "Treatment duration (weeks)", ylab = "Relapse probability",
         main = "Relapse probability vs duration")
    lines(bq$weeks, bq$P_relapse, type = "b", pch = 17, lwd = 2.4, col = PAL[["bad"]])
  })

  output$swexp <- renderUI(HTML(
    "<p style='color:#475569'>This tab is the conclusion of the whole model. The arm given rifampicin (green)
     has a <b>threshold</b> — 2 weeks is not enough, from 4 weeks on the reservoir falls below 1 CFU, and beyond that
     extending it changes nothing. The arm without rifampicin (red) <b>does not respond monotonically</b> to duration:
     from 2 weeks to 26 weeks the reservoir actually grows. When RSTER &le; 1, duration is not insufficient but
     <b>irrelevant</b>, and what is needed then is not longer antibiotics but the operating theatre again.</p>"))
}

shinyApp(ui, server)
