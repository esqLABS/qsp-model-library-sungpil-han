# =====================================================================
# Short Bowel Syndrome with Chronic Intestinal Failure (SBS-IF)
# Short Bowel Syndrome · Chronic Intestinal Failure — Shiny dashboard (10 tabs)
# ---------------------------------------------------------------------
# The dashboard is organised around the model's central claim: PARENTERAL
# SUPPORT VOLUME IS AN OUTPUT, NOT AN INPUT. So the first thing the user
# does is set an ANATOMY, and the app then SOLVES for the prescription
# that closes the fluid and energy balances. Nothing here lets you type in
# a PN volume, because in this model you cannot.
#
#   Tab 1   Patient profile     anatomy, diet, drinking behaviour -> phenotype
#   Tab 2   Two conservation equations the fluid and energy budgets, term by term
#   Tab 3   The water paradox    per-stream jejunal driving force and ORS sweep
#   Tab 4   PK                  GLP-2 analogue concentration and occupancy
#   Tab 5   Mucosal adaptation   villus, transporters, contact time, citrulline
#   Tab 6   Clinical endpoints   PN volume, infusion nights, 20% responder
#   Tab 7   Scenario comparison  drugs and adjuncts head to head
#   Tab 8   Biomarkers           citrulline, Mg, bicarbonate, micronutrients
#   Tab 9   Long-term complications IFALD, catheter/access, kidney, bone over 2 y
#   Tab 10  Weaning protocol     the weaning algorithm, and the placebo artefact
#
# RUN:  Rscript -e 'shiny::runApp("sbs_shiny_app.R", port = 8080)'
#       (the mrgsolve model is compiled once at start-up; allow ~30 s)
#
# NOTE: educational / research tool. Not validated for clinical
# decision-making, prescribing or regulatory submission.
# =====================================================================

library(shiny)

# ---------------------------------------------------------------------
# Load the model without triggering its scenario/diagnostic suite
# ---------------------------------------------------------------------
Sys.setenv(SBS_NORUN = "1")
MODEL_FILE <- "sbs_mrgsolve_model.R"
if (!file.exists(MODEL_FILE)) {
  # when launched from another working directory, look next to this file
  here <- tryCatch(dirname(normalizePath(sys.frames()[[1]]$ofile)),
                   error = function(e) ".")
  MODEL_FILE <- file.path(here, "sbs_mrgsolve_model.R")
}
suppressMessages(source(MODEL_FILE))   # also defines %||%, mp(), burnin(), sim()

PAL <- c(base = "#3f6d99", drug = "#a01050", alt = "#4a7a4a",
         warn = "#a04040", gold = "#a08000", grey = "#777777",
         purple = "#6a3a8a", orange = "#a04a10")

# ---------------------------------------------------------------------
# Small plotting helpers (base graphics only, no external dependencies)
# ---------------------------------------------------------------------
tsplot <- function(dl, col, ylab, main, cols = NULL, lty = NULL, hline = NULL,
                   legpos = "topright", xlab = "days") {
  dl <- dl[!vapply(dl, is.null, logical(1))]
  if (!length(dl)) { plot.new(); return(invisible()) }
  rng <- range(unlist(lapply(dl, function(d) d[[col]])), na.rm = TRUE)
  if (!is.finite(rng[1])) rng <- c(0, 1)
  if (diff(rng) < 1e-9) rng <- rng + c(-1, 1) * max(1e-3, abs(rng[1]) * 0.1)
  if (is.null(cols)) cols <- unname(PAL)[seq_along(dl)]
  if (is.null(lty))  lty  <- rep(1, length(dl))
  xr <- range(unlist(lapply(dl, function(d) d$time)))
  par(mar = c(4.2, 4.6, 3.0, 1.2))
  plot(NA, xlim = xr, ylim = rng, xlab = xlab, ylab = ylab, main = main,
       las = 1, cex.main = 1.02, font.main = 1)
  grid(col = "#e6e6e6", lty = 1)
  if (!is.null(hline)) abline(h = hline, col = "#bbbbbb", lty = 2)
  for (i in seq_along(dl))
    lines(dl[[i]]$time, dl[[i]][[col]], col = cols[i], lwd = 2.1, lty = lty[i])
  if (length(dl) > 1 || !is.null(names(dl)))
    legend(legpos, legend = names(dl), col = cols, lty = lty, lwd = 2.1,
           bty = "n", cex = 0.86)
}

barplot2 <- function(v, main, ylab, cols = NULL, horiz = FALSE, digits = 2) {
  par(mar = if (horiz) c(4.2, 11, 3, 1.2) else c(7.5, 4.6, 3, 1.2))
  if (is.null(cols)) cols <- rep(PAL["base"], length(v))
  b <- barplot(v, main = main, ylab = if (horiz) "" else ylab, col = cols,
               border = NA, las = if (horiz) 1 else 2, horiz = horiz,
               cex.names = 0.82, cex.main = 1.02, font.main = 1)
  if (!horiz) abline(h = 0, col = "#888888")
  text(if (horiz) v else b, if (horiz) b else v,
       labels = formatC(v, format = "f", digits = digits),
       pos = if (horiz) 4 else ifelse(v >= 0, 3, 1), cex = 0.76, xpd = TRUE)
  invisible(b)
}

waterfall <- function(labels, values, main, ylab) {
  par(mar = c(8.5, 4.8, 3, 1.2))
  cum <- c(0, cumsum(values))
  ylim <- range(c(cum, 0)) * 1.12
  plot(NA, xlim = c(0.4, length(values) + 0.6), ylim = ylim, xaxt = "n",
       xlab = "", ylab = ylab, main = main, las = 1, cex.main = 1.02,
       font.main = 1)
  grid(col = "#ececec", lty = 1); abline(h = 0, col = "#666666")
  for (i in seq_along(values)) {
    rect(i - 0.36, cum[i], i + 0.36, cum[i + 1], border = NA,
         col = if (values[i] >= 0) "#5b9bd5" else "#d9776a")
    text(i, max(cum[i], cum[i + 1]),
         sprintf("%+.2f", values[i]), pos = 3, cex = 0.72, xpd = TRUE)
  }
  lines(seq_along(values) + 0.36, cum[-1], type = "s", col = "#888888", lty = 3)
  axis(1, at = seq_along(values), labels = labels, las = 2, cex.axis = 0.78)
  invisible(cum)
}

kv_table <- function(m) {
  renderTable({ m }, rownames = FALSE, striped = TRUE, spacing = "xs",
              digits = 3)
}

# ---------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------
ui <- fluidPage(
  tags$head(tags$style(HTML("
    body { font-family: 'Helvetica Neue', Helvetica, Arial, sans-serif; }
    h4 { margin-top: 4px; }
    .thesis { background:#fffbe8; border-left:5px solid #a08000;
              padding:9px 13px; margin-bottom:11px; font-size:13.5px; }
    .caution { background:#fdf0f0; border-left:5px solid #a04040;
               padding:8px 12px; margin:9px 0; font-size:12.8px; }
    .note { color:#555; font-size:12.6px; margin:6px 0 12px 0; }
    .well { background:#f7f9fb; }
  "))),

  titlePanel("Short Bowel Syndrome with Chronic Intestinal Failure (SBS-IF) QSP dashboard"),

  div(class = "thesis",
      strong("The premise of this model: parenteral support volume is an output, not an input."),
      " What you set on the left is the ",
      em("anatomy, the diet and the drinking behaviour"), "; the app ",
      "then ", strong("computes"), " the prescription that closes the water, sodium and energy conservation equations. There is no box in which to enter a PN volume directly — ",
      "because in this model that is not possible.",
      br(),
      tags$span(style = "color:#555",
        "The thing you set is the anatomy; the prescription is SOLVED from the ",
        "fluid and energy budgets. There is deliberately no input box for PN volume.")),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("① Residual anatomy (Anatomy = parameter)"),
      selectInput("anat", "Messing anatomical type",
        c("Type 2 jejunocolic"                       = "jejunocolic",
          "Type 1 end-jejunostomy"                   = "end_jejunostomy",
          "Type 3 jejuno-ileo-colic"                 = "jejunoileocolic",
          "Custom"                                   = "custom"),
        selected = "jejunocolic"),
      sliderInput("sbl", "Residual small bowel length SBL (cm)", 20, 400, 80, 5),
      sliderInput("ileum", "Residual ileum (cm) — B12 · bile acids · L cells", 0, 200, 0, 5),
      sliderInput("colon", "Fraction of colon in continuity COLONFRAC", 0, 1, 0.5, 0.05),
      checkboxInput("icv", "Ileocaecal valve preserved", FALSE),
      sliderInput("mucqual", "Quality of the residual mucosa (1 = normal)", 0.5, 1, 1, 0.05),

      h4("② Diet and drinking"),
      sliderInput("ors", "ORS fraction (between-meal fluid)", 0, 0.95, 0.40, 0.05),
      sliderInput("drink", "Between-meal fluid intake (L/d)", 0.5, 5, 2.5, 0.1),
      sliderInput("kcal", "Oral energy target (kcal/d)", 1200, 4000, 2400, 100),
      sliderInput("entfrac", "Enteral intake multiplier (0 = complete fasting)", 0, 1.2, 1.0, 0.02),
      sliderInput("dietna", "Dietary sodium (mmol/d)", 60, 300, 160, 10),

      h4("③ Drug"),
      selectInput("drug", "GLP-2 therapy",
        c("None" = "none",
          "Teduglutide 0.05 mg/kg/day SC" = "ted",
          "Apraglutide 5 mg once weekly"       = "apra",
          "Glepaglutide 10 mg once weekly"     = "glep",
          "Glepaglutide 10 mg twice weekly"    = "glep2",
          "Native GLP-2, equimolar once daily" = "nat1",
          "Native GLP-2, equimolar twice daily" = "nat2"),
        selected = "ted"),
      checkboxGroupInput("adj", "Adjunctive therapy (adjuncts)",
        c("Loperamide 16 mg/day" = "lop",
          "High-dose PPI"        = "ppi",
          "Octreotide 100 µg tid" = "oct",
          "Somatropin 0.1 mg/kg/day x 4 weeks" = "gh",
          "Colestyramine 4 g tid" = "chol",
          "Rifaximin 550 mg bid x 14 days" = "rifax",
          "Potassium citrate · low-oxalate diet" = "cit",
          "Taurolidine catheter lock" = "lock")),
      selectInput("ile", "PN lipid emulsion",
        c("Soybean oil (phytosterol ~350)" = "soy",
          "SMOF composite"                 = "smof",
          "Fish-oil based (phytosterol 0)" = "fish",
          "Lipid minimisation"             = "min"),
        selected = "soy"),

      h4("④ Weaning protocol"),
      checkboxInput("weanon", "Urine-output-driven weaning algorithm active", TRUE),
      sliderInput("urtrig", "PN taper trigger (urine output / baseline urine output)", 1.02, 1.40, 1.10, 0.01),
      sliderInput("ktaper", "Maximum taper rate (1/day)", 0.001, 0.020, 0.006, 0.001),
      sliderInput("days", "Simulation duration (days)", 60, 1095, 168, 7),
      actionButton("go", "Run simulation", class = "btn-primary",
                   style = "width:100%; margin-top:8px;"),
      div(class = "note",
          "Change the anatomy and the model first re-finds the steady state (burn-in), ",
          "then runs the scenario. Every comparison therefore starts from a self-consistent patient.")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs", type = "tabs",

        # ---------------- TAB 1 -------------------------------------
        tabPanel("1 · Patient profile",
          h4("The phenotype the anatomy produces"),
          div(class = "note",
            "None of the values below is an input. Only the residual bowel length, colon, ileocaecal valve, diet and drinking ",
            "behaviour were given; everything else is what the model computed at steady state."),
          fluidRow(column(6, tableOutput("t1_pheno")),
                   column(6, tableOutput("t1_absorp"))),
          hr(),
          fluidRow(column(6, plotOutput("t1_anat", height = 300)),
                   column(6, plotOutput("t1_thresh", height = 300))),
          div(class = "caution",
            strong("Anatomy is a parameter, not a state."),
            " Messing's thresholds for permanent intestinal failure (end-jejunostomy ~115 cm, jejunocolic ~60 cm, ",
            "jejuno-ileo-colic ~35 cm) differ by nearly three-fold because the colon recovers 3-4 L of ",
            "water and up to ~1000 kcal a day.")
        ),

        # ---------------- TAB 2 -------------------------------------
        tabPanel("2 · Two conservation equations",
          h4("The water, sodium and energy budgets (the budgets whose residual the central line fills)"),
          fluidRow(column(6, plotOutput("t2_fluid", height = 340)),
                   column(6, plotOutput("t2_energy", height = 340))),
          hr(),
          fluidRow(column(7, tableOutput("t2_terms")),
                   column(5, plotOutput("t2_na", height = 300))),
          div(class = "note",
            "However far the last bar of the left-hand waterfall departs from 0 is exactly the deficit that has to be ",
            "filled intravenously. What the weaning algorithm does is not to recompute that deficit ",
            "but to track it through an observable signal, the urine output.")
        ),

        # ---------------- TAB 3 -------------------------------------
        tabPanel("3 · The water paradox",
          h4("The sign change in the jejunal Na⁺ driving force (where the sign changes)"),
          div(class = "note",
            "The fluid that arrives with a meal (stream M) and the fluid drunk between meals (stream D) ",
            "differ completely in sodium concentration and glucose content. The driving force is computed for each ",
            "stream separately and only then averaged by volume fraction. Average first and the paradox disappears."),
          fluidRow(column(6, plotOutput("t3_drive", height = 340)),
                   column(6, plotOutput("t3_output", height = 340))),
          hr(),
          fluidRow(column(7, tableOutput("t3_tab")),
                   column(5, plotOutput("t3_pn", height = 320))),
          div(class = "caution",
            strong("A thirsty end-jejunostomy patient who drinks 2 L of plain water LOSES sodium ",
            "and volume."), " Give the same 2 L as a glucose-saline of Na⁺ 90 mmol/L and the patient gains ",
            "both. This is not a warning notice but a sign change in the flux equation.")
        ),

        # ---------------- TAB 4 -------------------------------------
        tabPanel("4 · PK",
          h4("GLP-2 analogue pharmacokinetics and receptor occupancy"),
          fluidRow(column(6, plotOutput("t4_conc", height = 330)),
                   column(6, plotOutput("t4_occ", height = 330))),
          hr(),
          fluidRow(column(6, plotOutput("t4_zoom", height = 320)),
                   column(6, tableOutput("t4_tab"))),
          div(class = "caution",
            strong("The point is that the time constant of the mucosal response (weeks) is much longer ",
            "than the half-life of the drug (hours)."), " That is why a drug with a t½ of 2 hours works ",
            "given once daily. Native GLP-2, by contrast, is cleaved by DPP-4 at Ala2-Gly3 and has ",
            "a t½ of about 7 minutes, so an equimolar once-daily dose cannot sustain occupancy.")
        ),

        # ---------------- TAB 5 -------------------------------------
        tabPanel("5 · Mucosal adaptation",
          h4("Structural and functional adaptation — a multiplicative structure gated by nutrient"),
          fluidRow(column(6, plotOutput("t5_vill", height = 330)),
                   column(6, plotOutput("t5_trans", height = 330))),
          hr(),
          fluidRow(column(6, plotOutput("t5_gate", height = 330)),
                   column(6, plotOutput("t5_contact", height = 330))),
          div(class = "note",
            "TROPHIC = (luminal nutrient) x (1 + GLP2R occupancy term) x (1 + IGF-1 term). Being a product, ",
            "an enteral intake of 0 makes the nutrient signal 0 as well and the trophic arm stops ",
            "completely. The motility (ENS) arm of GLP-2, however, does not require luminal ",
            "nutrient, so the drug effect does not disappear entirely — a design intention the model ",
            "refuted by itself, and a testable structural prediction.")
        ),

        # ---------------- TAB 6 -------------------------------------
        tabPanel("6 · Clinical endpoints",
          h4("PN volume · infusion nights per week · ≥20% responder"),
          fluidRow(column(6, plotOutput("t6_pn", height = 330)),
                   column(6, plotOutput("t6_nights", height = 330))),
          hr(),
          fluidRow(column(5, tableOutput("t6_tab")),
                   column(7, plotOutput("t6_out", height = 320))),
          div(class = "note",
            "The STEPS primary endpoint is a patient whose PN volume fell by at least 20% at both week 20 ",
            "and week 24. In this app that value comes out of a simulation of the trial protocol, ",
            "not out of a regression. Note too that it is not the litres but the ",
            strong("nights"), " that have to fall before catheter-days fall.")
        ),

        # ---------------- TAB 7 -------------------------------------
        tabPanel("7 · Scenario comparison",
          h4("Drugs and adjuncts head to head (matched anatomy, matched protocol)"),
          checkboxGroupInput("cmp", "Scenarios to compare",
            c("No drug" = "none",
              "Teduglutide" = "ted",
              "Apraglutide once weekly" = "apra",
              "Glepaglutide once weekly" = "glep",
              "Glepaglutide twice weekly" = "glep2",
              "Native GLP-2, equimolar once daily" = "nat1",
              "ORS optimised (0.9)" = "ors",
              "Loperamide" = "lop",
              "High-dose PPI" = "ppi",
              "Octreotide" = "oct",
              "Teduglutide + ORS + loperamide" = "combo"),
            selected = c("none", "ted", "ors", "combo"), inline = TRUE),
          fluidRow(column(7, plotOutput("t7_pn", height = 350)),
                   column(5, plotOutput("t7_bar", height = 350))),
          hr(),
          tableOutput("t7_tab"),
          div(class = "caution",
            strong("Octreotide buys today's output and sells tomorrow's adaptation."),
            " In exchange for reducing the secretory volume it suppresses IGF-1 and suppresses ",
            "pancreatic exocrine secretion, lowering fat absorption. See for yourself whether the ",
            "sign of the benefit changes with time.")
        ),

        # ---------------- TAB 8 -------------------------------------
        tabPanel("8 · Biomarkers",
          h4("Citrulline · magnesium · acid-base · micronutrients"),
          fluidRow(column(6, plotOutput("t8_citr", height = 320)),
                   column(6, plotOutput("t8_mg", height = 320))),
          hr(),
          fluidRow(column(6, plotOutput("t8_acid", height = 320)),
                   column(6, plotOutput("t8_micro", height = 320))),
          div(class = "caution",
            strong("Mg deficiency impairs PTH secretion itself and at the same time produces ",
            "target-organ resistance."), " That is why hypocalcaemia does not respond to calcium ",
            "administration until the magnesium has been replaced. This is not asserted in the model but ",
            "generated by two Mg-dependent terms. Citrulline <20 µmol/L is treated as a direct ",
            "index of enterocyte mass.")
        ),

        # ---------------- TAB 9 -------------------------------------
        tabPanel("9 · Long-term complications",
          h4("Two-year course: IFALD · catheter and venous access · kidney · bone"),
          div(class = "note",
            "In modern SBS-IF what determines death is not the bowel that is left but the plumbing of the treatment. ",
            "The competing-risk structure in which the bowel improves while the liver and the veins deteriorate is the real decision problem."),
          fluidRow(column(6, plotOutput("t9_liver", height = 320)),
                   column(6, plotOutput("t9_cvc", height = 320))),
          hr(),
          fluidRow(column(6, plotOutput("t9_renal", height = 320)),
                   column(6, plotOutput("t9_bone", height = 320))),
          tableOutput("t9_tab")
        ),

        # ---------------- TAB 10 ------------------------------------
        tabPanel("10 · Weaning protocol",
          h4("The urine-output-driven taper algorithm, and the placebo response as an artefact"),
          div(class = "note",
            "In STEPS what decided the PN taper was not the investigator but the protocol: ",
            "if the 24-hour urine output rises more than 10% above baseline, reduce the PN; if the weight or the electrolytes ",
            "fall, hold. Implement that algorithm as an actual closed-loop controller and ",
            "two things follow — the taper stops itself (self-limiting), and a patient with absorptive reserve ",
            "left is weaned by the protocol regardless of the drug."),
          fluidRow(column(6, plotOutput("t10_loop", height = 330)),
                   column(6, plotOutput("t10_placebo", height = 330))),
          hr(),
          fluidRow(column(6, plotOutput("t10_trig", height = 330)),
                   column(6, tableOutput("t10_tab"))),
          div(class = "caution",
            strong("Much of the placebo-arm response in SBS trials is not a placebo effect but ",
            "the effect of well-organised care."), " On trial entry both arms have their oral rehydration and ",
            "diet standardised, and the urine-output trigger translates that improvement into a taper.")
        )
      ),
      hr(),
      div(class = "note",
        strong("Disclaimer:"), " This is a QSP model for education and research. It must not be used for ",
        "clinical decision-making, prescribing or regulatory submission. The absorption and fluid coefficients are ",
        "population parameters calibrated to reproduce published phenotypes and cannot be applied to an individual patient. ",
        "For the evidence see ", code("sbs_references.md"), " (236 papers, every PMID verified); for the structure see ",
        code("sbs_qsp_model.svg"), " (23 clusters · 276 nodes).")
    )
  )
)

# ---------------------------------------------------------------------
# SERVER
# ---------------------------------------------------------------------
server <- function(input, output, session) {

  observeEvent(input$anat, {
    a <- input$anat
    if (a == "end_jejunostomy") {
      updateSliderInput(session, "sbl", value = 90)
      updateSliderInput(session, "ileum", value = 0)
      updateSliderInput(session, "colon", value = 0)
      updateCheckboxInput(session, "icv", value = FALSE)
    } else if (a == "jejunocolic") {
      updateSliderInput(session, "sbl", value = 80)
      updateSliderInput(session, "ileum", value = 0)
      updateSliderInput(session, "colon", value = 0.5)
      updateCheckboxInput(session, "icv", value = FALSE)
    } else if (a == "jejunoileocolic") {
      updateSliderInput(session, "sbl", value = 70)
      updateSliderInput(session, "ileum", value = 25)
      updateSliderInput(session, "colon", value = 1)
      updateCheckboxInput(session, "icv", value = TRUE)
    }
  }, ignoreInit = TRUE)

  # ---- parameter set from the UI -------------------------------------
  basepars <- reactive({
    ilepar <- switch(input$ile,
      soy  = list(PHYTOCONC = 350, OMEGA3 = 0, FLIPIDPN = 1,   ILEDOSE = 0.90),
      smof = list(PHYTOCONC = 100, OMEGA3 = 1, FLIPIDPN = 1,   ILEDOSE = 0.90),
      fish = list(PHYTOCONC = 0,   OMEGA3 = 1, FLIPIDPN = 1,   ILEDOSE = 0.90),
      min  = list(PHYTOCONC = 350, OMEGA3 = 0, FLIPIDPN = 1,   ILEDOSE = 0.14))
    adj <- input$adj %||% character(0)
    c(list(SBL = input$sbl, ILEUMLEN = min(input$ileum, input$sbl),
           COLONFRAC = input$colon, ICV = as.numeric(isTRUE(input$icv)),
           MUCQUAL = input$mucqual, ORSFRAC = input$ors,
           DRINKVOL = input$drink, KCALORAL = input$kcal,
           ENTFRAC = input$entfrac, DIETNA = input$dietna,
           PPIDOSE = if ("ppi" %in% adj) 2 else 0,
           OXCITRATE = as.numeric("cit" %in% adj),
           LOCKFAC = if ("lock" %in% adj) 0.25 else 1,
           WEANON = as.numeric(isTRUE(input$weanon)),
           URTRIG = input$urtrig, KTAPER = input$ktaper),
      ilepar)
  })

  # burn-in depends only on the physiology, not on the drug
  bstate <- reactive({
    p <- basepars()
    p$WEANON <- 1   # the run-in must be allowed to find the prescription
    withProgress(message = "Finding the steady state (burn-in)...", value = 0.4,
                 burnin(p))
  })

  runpars <- reactive({
    st <- bstate(); p <- basepars()
    modifyList(p, list(URBASE = st[["URINESM"]],
                       WTTARGET = st[["LEAN"]] + st[["FATM"]]))
  })

  wt0 <- reactive({ st <- bstate(); st[["LEAN"]] + st[["FATM"]] })

  drug_ev <- function(which, days, wt) {
    switch(which,
      none  = NULL,
      ted   = ev_ted(wt, days),
      apra  = ev_apra(days),
      glep  = ev_glep(days, ii = 7),
      glep2 = ev_glep(days, ii = 3.5),
      nat1  = ev_native(wt, days, 1),
      nat2  = ev_native(wt, days, 2),
      NULL)
  }
  adj_ev <- function(adj, days, wt) {
    e <- NULL
    add <- function(x) if (is.null(e)) x else c(e, x)
    if ("lop"   %in% adj) e <- add(ev_lop(days))
    if ("oct"   %in% adj) e <- add(ev_oct(days))
    if ("gh"    %in% adj) e <- add(ev_gh(wt, min(28, days)))
    if ("chol"  %in% adj) e <- add(ev_chol(days))
    if ("rifax" %in% adj) e <- add(ev_rifax(min(14, days)))
    e
  }

  simrun <- eventReactive(input$go, {
    st <- bstate(); p <- runpars(); d <- input$days; w <- wt0()
    adj <- input$adj %||% character(0)
    ev1 <- drug_ev(input$drug, d, w); ev2 <- adj_ev(adj, d, w)
    ee <- if (is.null(ev1)) ev2 else if (is.null(ev2)) ev1 else c(ev1, ev2)
    withProgress(message = "Running the simulation...", value = 0.6, {
      list(drug = sim(st, p, ee, end = d),
           ref  = sim(st, p, NULL, end = d),
           st = st, p = p, wt = w, days = d)
    })
  }, ignoreNULL = FALSE)

  # =================== TAB 1 ==========================================
  output$t1_pheno <- kv_table({
    r <- simrun(); d <- r$ref
    data.frame(
      `Metric` = c("Prescribed PN volume (L/week)", "Infusion nights per week",
        "Prescribed PN energy (kcal/day)", "Faecal/stomal output (L/day)",
        "Stomal effluent [Na⁺] (mmol/L)", "24-hour urine output (L/day)",
        "Urine sodium (mmol/day)", "Body weight (kg)", "Plasma citrulline (µmol/L)",
        "Hyperphagia multiplier", "Mucosal contact time (relative)"),
      `Value` = round(c(d$PNVOLWK[1], d$PNNIGHTS[1], d$PNKCALD[1],
        d$OUT_LD[1], d$STOMANA[1], d$URINE_LD[1], d$UNA[1], d$WT[1],
        d$CITR[1], d$HYPERPH[1], d$CONTC[1]), 2),
      check.names = FALSE)
  })

  output$t1_absorp <- kv_table({
    r <- simrun(); d <- r$ref
    data.frame(
      `Absorption` = c("Carbohydrate absorption (%)", "Protein absorption (%)",
        "Fat absorption (%)", "Absorbed energy (kcal/day)", "Colonic SCFA salvage (kcal/day)",
        "Total energy expenditure (kcal/day)", "Unabsorbed fat (g/day)",
        "Net jejunal water flux (L/day)", "Water recovered by the colon (L/day)",
        "Bile acid pool (g)", "Micellar function (relative)"),
      `Value` = round(c(100*d$FCHOC[1], 100*d$FPROC[1], 100*d$FFATC[1],
        d$ABSKC[1], d$SCFAKC[1], d$TEEC[1], d$FATMALC[1], d$VJEJC[1],
        d$VCOLC[1], d$BAPOOLC[1], d$MICELC[1]), 2),
      check.names = FALSE)
  })

  output$t1_anat <- renderPlot({
    r <- simrun(); p <- r$p
    v <- c(`residual small bowel\nSBL (cm)` = p$SBL,
           `of which ileum\nileum (cm)` = p$ILEUMLEN,
           `colon %\ncolon` = 100*p$COLONFRAC,
           `ileocaecal valve\nICV` = 100*p$ICV)
    barplot2(v, "The anatomy you set", "cm or %",
             cols = c(PAL["base"], PAL["alt"], PAL["purple"], PAL["gold"]),
             digits = 0)
  })

  output$t1_thresh <- renderPlot({
    st <- bstate(); p <- basepars()
    ls <- seq(30, 250, by = 20)
    pn <- vapply(ls, function(L) {
      q <- modifyList(p, list(SBL = L, ILEUMLEN = min(p$ILEUMLEN, L)))
      s <- burnin(q, days = 500); s[["PNVOL"]] * 7
    }, numeric(1))
    par(mar = c(4.4, 4.8, 3, 1.2))
    plot(ls, pn, type = "b", pch = 19, col = PAL["base"], lwd = 2.1, las = 1,
         xlab = "Residual small bowel length (cm)", ylab = "Steady-state PN volume (L/week)",
         main = "Length versus PN dependence (at the current colon/ICV setting)",
         cex.main = 1.0, font.main = 1)
    grid(col = "#ececec", lty = 1)
    abline(v = p$SBL, col = PAL["drug"], lty = 2)
    abline(h = 0.35 * 7, col = "#999999", lty = 3)
    legend("topright", c("Current patient", "Approximate enteral autonomy threshold"),
           col = c(PAL["drug"], "#999999"), lty = c(2, 3), bty = "n", cex = 0.85)
  })

  # =================== TAB 2 ==========================================
  output$t2_fluid <- renderPlot({
    r <- simrun(); d <- r$ref; i <- 1
    lab <- c("Oral intake", "PN", "Output", "Insensible loss", "Urine")
    val <- c(d$VORALC[i], d$PNVOLWK[i]/7, -d$OUT_LD[i], -0.70, -d$URINE_LD[i])
    waterfall(lab, val, "Water conservation equation (L/day) — in balance when the final residual is 0",
              "L/day")
  })

  output$t2_energy <- renderPlot({
    r <- simrun(); d <- r$ref; i <- 1
    lab <- c("Absorbed energy", "Colonic SCFA salvage", "PN energy",
             "Total expenditure\n-TEE")
    val <- c(d$ABSKC[i], d$SCFAKC[i], d$PNKCALD[i], -d$TEEC[i])
    waterfall(lab, val, "Energy conservation equation (kcal/day)", "kcal/day")
  })

  output$t2_terms <- kv_table({
    r <- simrun(); d <- r$ref; i <- 1
    data.frame(
      `Term` = c("Volume delivered to the jejunum VDEL (L/day)", "Luminal [Na⁺] (mmol/L)",
        "Jejunal absorbed fraction FRACJ", "Net jejunal flux (L/day)", "Colonic inflow (L/day)",
        "Colonic recovery (L/day)", "Output (L/day)", "PN volume (L/day)",
        "Urine output (L/day)", "Body water deviation (L)"),
      `Value` = round(c(d$VORALC[i] + 5.2, d$CNALUM[i], d$FRACJC[i], d$VJEJC[i],
        d$VORALC[i] + 5.2 - d$VJEJC[i], d$VCOLC[i], d$OUT_LD[i],
        d$PNVOLWK[i]/7, d$URINE_LD[i], d$TBWC[i] - 36), 3),
      check.names = FALSE)
  })

  output$t2_na <- renderPlot({
    r <- simrun()
    tsplot(list(`Urine sodium (mmol/day)` = r$ref, `On drug` = r$drug),
           "UNA", "Urine sodium (mmol/day)",
           "The earliest index of sodium depletion: urine Na⁺ < 10",
           cols = c(PAL["grey"], PAL["drug"]), hline = 10)
  })

  # =================== TAB 3 ==========================================
  ors_sweep <- reactive({
    st <- bstate(); p <- runpars()
    of <- seq(0, 0.95, by = 0.0950)
    do.call(rbind, lapply(of, function(x) {
      d <- sim(st, modifyList(p, list(ORSFRAC = x, WEANON = 0)), end = 45)
      j <- nrow(d)
      data.frame(ORS = x, CNAD = d$CNAD[j], DRIVEM = d$DRIVEM[j],
                 DRIVED = d$DRIVED[j], DRIVE = d$DRIVEC[j],
                 FRACJ = d$FRACJC[j], OUT = d$OUT_LD[j], VJEJ = d$VJEJC[j])
    }))
  })

  output$t3_drive <- renderPlot({
    s <- ors_sweep()
    par(mar = c(4.4, 4.8, 3.2, 1.2))
    plot(s$ORS, s$DRIVED, type = "b", pch = 19, lwd = 2.2, col = PAL["warn"],
         ylim = range(c(s$DRIVED, s$DRIVEM, 0)), las = 1,
         xlab = "ORS fraction of the between-meal fluid",
         ylab = "Driving force [Na⁺]lum − CEQ (mmol/L)",
         main = "Driving force by stream: the sign changes only in stream D",
         cex.main = 1.0, font.main = 1)
    grid(col = "#ececec", lty = 1); abline(h = 0, lwd = 2, col = "#333333")
    lines(s$ORS, s$DRIVEM, type = "b", pch = 17, lwd = 2.2, col = PAL["alt"])
    lines(s$ORS, s$DRIVE,  type = "b", pch = 15, lwd = 2.2, col = PAL["base"])
    text(0.06, 6, "absorbs", cex = 0.8, col = "#333333", adj = 0)
    text(0.06, -8, "SECRETES", cex = 0.8, col = PAL["warn"], adj = 0)
    legend("bottomright", c("stream D (between-meal drinking)", "stream M (meal + secretion)",
                            "volume-weighted mean"),
           col = c(PAL["warn"], PAL["alt"], PAL["base"]), pch = c(19, 17, 15),
           lwd = 2.2, bty = "n", cex = 0.84)
  })

  output$t3_output <- renderPlot({
    s <- ors_sweep()
    par(mar = c(4.4, 4.8, 3.2, 4.6))
    plot(s$ORS, s$OUT, type = "b", pch = 19, lwd = 2.2, col = PAL["orange"],
         las = 1, xlab = "ORS fraction", ylab = "Output (L/day)",
         main = "The same drinking volume, changing only its composition: output and jejunal flux",
         cex.main = 1.0, font.main = 1)
    grid(col = "#ececec", lty = 1)
    par(new = TRUE)
    plot(s$ORS, s$VJEJ, type = "b", pch = 17, lwd = 2.2, col = PAL["base"],
         axes = FALSE, xlab = "", ylab = "")
    axis(4, las = 1); mtext("Net jejunal water flux (L/day)", side = 4, line = 3.2,
                            cex = 0.9)
    legend("topright", c("Output (left)", "Jejunal flux (right)"),
           col = c(PAL["orange"], PAL["base"]), pch = c(19, 17), lwd = 2.2,
           bty = "n", cex = 0.84)
  })

  output$t3_tab <- kv_table({
    s <- ors_sweep()
    data.frame(
      `ORS fraction` = round(s$ORS, 2),
      `Drink [Na⁺]` = round(s$CNAD, 1),
      `stream D driving force` = round(s$DRIVED, 1),
      `Direction` = ifelse(s$DRIVED < 0, "secretes", "absorbs"),
      `FRACJ` = round(s$FRACJ, 3),
      `Output L/day` = round(s$OUT, 2),
      check.names = FALSE)
  })

  output$t3_pn <- renderPlot({
    st <- bstate(); p <- runpars(); d <- max(84, input$days)
    a <- sim(st, modifyList(p, list(ORSFRAC = 0)),    end = d)
    b <- sim(st, modifyList(p, list(ORSFRAC = 0.45)), end = d)
    cc<- sim(st, modifyList(p, list(ORSFRAC = 0.90)), end = d)
    tsplot(list(`Plain water only (ORS 0)` = a, `ORS 0.45` = b, `ORS 0.90` = cc),
           "PNVOLWK", "PN volume (L/week)",
           "The PN trajectory when only the oral rehydration is changed",
           cols = c(PAL["warn"], PAL["gold"], PAL["alt"]))
  })

  # =================== TAB 4 ==========================================
  pk_all <- reactive({
    st <- bstate(); p <- runpars(); w <- wt0(); d <- min(60, input$days)
    m <- init(param(mod, p), as.list(st))
    f <- function(ee) as.data.frame(mrgsim(m, ee, end = d, delta = 0.02))
    list(`Teduglutide od` = f(ev_ted(w, d)),
         `Apraglutide qw` = f(ev_apra(d)),
         `Glepaglutide qw` = f(ev_glep(d, ii = 7)),
         `Native GLP-2 equimolar od` = f(ev_native(w, d, 1)))
  })

  output$t4_conc <- renderPlot({
    l <- pk_all()
    dl <- list(`Teduglutide` = l[[1]], `Apraglutide` = l[[2]],
               `Glepaglutide` = l[[3]])
    tsplot(dl, "CPTED", "Teduglutide concentration (mg/L)",
           "Concentration after SC dosing (teduglutide axis)",
           cols = c(PAL["drug"], PAL["base"], PAL["alt"]))
  })

  output$t4_occ <- renderPlot({
    l <- pk_all()
    tsplot(l, "OCCGC", "GLP2R occupancy (0-1)",
           "Receptor occupancy — the effect-site compartment integrates the time scale",
           cols = c(PAL["drug"], PAL["base"], PAL["alt"], PAL["warn"]))
  })

  output$t4_zoom <- renderPlot({
    l <- pk_all()
    dl <- lapply(l, function(d) d[d$time <= 3, ])
    tsplot(dl, "OCCGC", "GLP2R occupancy",
           "First 3 days magnified: with a t½ of 7 minutes, native GLP-2 does not sustain occupancy",
           cols = c(PAL["drug"], PAL["base"], PAL["alt"], PAL["warn"]))
  })

  output$t4_tab <- kv_table({
    l <- pk_all()
    do.call(rbind, lapply(names(l), function(n) {
      d <- l[[n]]; tail30 <- d[d$time >= max(d$time) - 7, ]
      data.frame(`Regimen` = n,
        `Mean occupancy` = round(mean(d$OCCGC), 3),
        `Peak occupancy` = round(max(tail30$OCCGC), 3),
        `Trough occupancy` = round(min(tail30$OCCGC), 3),
        `Trough/peak` = round(min(tail30$OCCGC)/max(max(tail30$OCCGC), 1e-9), 3),
        check.names = FALSE)
    }))
  })

  # =================== TAB 5 ==========================================
  output$t5_vill <- renderPlot({
    r <- simrun()
    tsplot(list(`No drug` = r$ref, `On drug` = r$drug), "VILLC",
           "Effective villus height (relative)", "Structural adaptation",
           cols = c(PAL["grey"], PAL["drug"]))
  })
  output$t5_trans <- renderPlot({
    r <- simrun()
    tsplot(list(`No drug` = r$ref, `On drug` = r$drug), "MUCRELC",
           "Mucosal mass index MUCREL", "Mucosal mass (the determinant of citrulline)",
           cols = c(PAL["grey"], PAL["drug"]))
  })
  output$t5_gate <- renderPlot({
    st <- bstate(); p <- runpars(); w <- wt0(); d <- input$days
    out <- lapply(c(1.0, 0.5, 0.15, 0.02), function(ef) {
      s <- burnin(modifyList(p, list(ENTFRAC = ef)), days = 600)
      q <- modifyList(p, list(ENTFRAC = ef, URBASE = s[["URINESM"]],
                              WTTARGET = s[["LEAN"]] + s[["FATM"]]))
      sim(s, q, ev_ted(s[["LEAN"]] + s[["FATM"]], d), end = d)
    })
    names(out) <- c("Enteral 100%", "Enteral 50%", "Enteral 15%", "Complete fasting 2%")
    tsplot(out, "TROPHC", "Nutrient signal TROPHIC",
           "Nutrient gating: being a product, fasting stops the trophic arm",
           cols = c(PAL["alt"], PAL["gold"], PAL["orange"], PAL["warn"]))
  })
  output$t5_contact <- renderPlot({
    r <- simrun()
    tsplot(list(`No drug` = r$ref, `On drug` = r$drug), "CONTC",
           "Mucosal contact time (relative)",
           "Contact time through the ileal brake — a multiplicative factor on every absorbed fraction",
           cols = c(PAL["grey"], PAL["drug"]))
  })

  # =================== TAB 6 ==========================================
  output$t6_pn <- renderPlot({
    r <- simrun()
    tsplot(list(`No drug (protocol only)` = r$ref, `Selected regimen` = r$drug),
           "PNVOLWK", "PN volume (L/week)",
           "PN volume — the trajectory the algorithm produced",
           cols = c(PAL["grey"], PAL["drug"]),
           hline = 0.8 * r$ref$PNVOLWK[1])
  })
  output$t6_nights <- renderPlot({
    r <- simrun()
    tsplot(list(`No drug` = r$ref, `Selected regimen` = r$drug), "PNNIGHTS",
           "Infusion nights per week", "Catheter-days fall only when the nights fall, not the litres",
           cols = c(PAL["grey"], PAL["drug"]))
  })
  output$t6_tab <- kv_table({
    r <- simrun(); b0 <- r$ref$PNVOLWK[1]; n <- nrow(r$drug)
    wk20 <- which.min(abs(r$drug$time - 140)); wk24 <- which.min(abs(r$drug$time - 168))
    red <- function(d, i) 100 * (b0 - d$PNVOLWK[i]) / b0
    data.frame(
      `Endpoint` = c("Baseline PN volume (L/week)", "Final PN volume (L/week)",
        "Change in PN volume (L/week)", "PN reduction (%)",
        "Reduction at week 20 (%)", "Reduction at week 24 (%)",
        "≥20% responder (week 20 AND week 24)", "Change in nights per week",
        "Cumulative catheter-days", "CRBSI per year"),
      `On drug` = c(round(b0, 2), round(r$drug$PNVOLWK[n], 2),
        round(r$drug$PNVOLWK[n] - b0, 2), round(red(r$drug, n), 1),
        round(red(r$drug, wk20), 1), round(red(r$drug, wk24), 1),
        ifelse(resp20(r$drug, b0) == 1, "YES", "no"),
        round(r$drug$PNNIGHTS[n] - r$drug$PNNIGHTS[1], 2),
        round(r$drug$CRBSICUM[n] * 0 + r$drug$time[n] * r$drug$PNNIGHTS[n]/7, 0),
        round(r$drug$CRBSIYR[n], 3)),
      `No drug` = c(round(b0, 2), round(r$ref$PNVOLWK[n], 2),
        round(r$ref$PNVOLWK[n] - b0, 2), round(red(r$ref, n), 1),
        round(red(r$ref, wk20), 1), round(red(r$ref, wk24), 1),
        ifelse(resp20(r$ref, b0) == 1, "YES", "no"),
        round(r$ref$PNNIGHTS[n] - r$ref$PNNIGHTS[1], 2),
        round(r$ref$time[n] * r$ref$PNNIGHTS[n]/7, 0),
        round(r$ref$CRBSIYR[n], 3)),
      check.names = FALSE)
  })
  output$t6_out <- renderPlot({
    r <- simrun()
    tsplot(list(`No drug` = r$ref, `Selected regimen` = r$drug), "OUT_LD",
           "Faecal/stomal output (L/day)", "Output — the quantity that can be measured",
           cols = c(PAL["grey"], PAL["drug"]))
  })

  # =================== TAB 7 ==========================================
  cmp_runs <- eventReactive(list(input$go, input$cmp), {
    st <- bstate(); p <- runpars(); w <- wt0(); d <- input$days
    sel <- input$cmp %||% "none"
    withProgress(message = "Computing the scenario comparison...", value = 0.5, {
      out <- lapply(sel, function(k) {
        switch(k,
          none  = sim(st, p, NULL, end = d),
          ted   = sim(st, p, ev_ted(w, d), end = d),
          apra  = sim(st, p, ev_apra(d), end = d),
          glep  = sim(st, p, ev_glep(d, ii = 7), end = d),
          glep2 = sim(st, p, ev_glep(d, ii = 3.5), end = d),
          nat1  = sim(st, p, ev_native(w, d, 1), end = d),
          ors   = sim(st, modifyList(p, list(ORSFRAC = 0.9)), NULL, end = d),
          lop   = sim(st, p, ev_lop(d), end = d),
          ppi   = sim(st, modifyList(p, list(PPIDOSE = 2)), NULL, end = d),
          oct   = sim(st, p, ev_oct(d), end = d),
          combo = sim(st, modifyList(p, list(ORSFRAC = 0.9)),
                      c(ev_ted(w, d), ev_lop(d)), end = d))
      })
      names(out) <- sel; out
    })
  }, ignoreNULL = FALSE)

  LBL <- c(none = "No drug", ted = "Teduglutide", apra = "Apraglutide qw",
           glep = "Glepaglutide qw", glep2 = "Glepaglutide 2x/wk",
           nat1 = "Native GLP-2 od", ors = "ORS 0.9", lop = "Loperamide",
           ppi = "High-dose PPI", oct = "Octreotide",
           combo = "Ted + ORS + loperamide")

  output$t7_pn <- renderPlot({
    l <- cmp_runs(); if (!length(l)) return(invisible())
    names(l) <- LBL[names(l)]
    tsplot(l, "PNVOLWK", "PN volume (L/week)", "PN trajectory by scenario")
  })
  output$t7_bar <- renderPlot({
    l <- cmp_runs(); if (!length(l)) return(invisible())
    v <- vapply(l, function(d) d$PNVOLWK[nrow(d)] - d$PNVOLWK[1], numeric(1))
    names(v) <- LBL[names(l)]
    barplot2(v, "Change in final PN volume (L/week)", "L/week",
             cols = ifelse(v < 0, PAL["alt"], PAL["warn"]), horiz = TRUE)
  })
  output$t7_tab <- kv_table({
    l <- cmp_runs(); if (!length(l)) return(data.frame())
    do.call(rbind, lapply(names(l), function(k) {
      d <- l[[k]]; n <- nrow(d)
      data.frame(`Scenario` = LBL[[k]],
        `PN L/week` = round(d$PNVOLWK[n], 2),
        `Change` = round(d$PNVOLWK[n] - d$PNVOLWK[1], 2),
        `Nights/week` = round(d$PNNIGHTS[n], 2),
        `Output L/day` = round(d$OUT_LD[n], 2),
        `Fat absorption %` = round(100*d$FFATC[n], 1),
        `Villus` = round(d$VILLC[n], 3),
        `Citrulline` = round(d$CITR[n], 1),
        `Weight kg` = round(d$WT[n], 1),
        `≥20% response` = ifelse(resp20(d, d$PNVOLWK[1]) == 1, "yes", "-"),
        check.names = FALSE)
    }))
  })

  # =================== TAB 8 ==========================================
  output$t8_citr <- renderPlot({
    r <- simrun()
    tsplot(list(`No drug` = r$ref, `On drug` = r$drug), "CITR",
           "Plasma citrulline (µmol/L)",
           "Citrulline — the index of enterocyte mass (<20 = suggests permanent intestinal failure)",
           cols = c(PAL["grey"], PAL["drug"]), hline = 20)
  })
  output$t8_mg <- renderPlot({
    r <- simrun()
    tsplot(list(`Serum Mg (no drug)` = r$ref, `Serum Mg (on drug)` = r$drug),
           "SMGC", "Serum magnesium (mmol/L)",
           "Magnesium — lost in proportion to the output and buffered from bone",
           cols = c(PAL["grey"], PAL["drug"]), hline = 0.70)
  })
  output$t8_acid <- renderPlot({
    r <- simrun()
    tsplot(list(`Bicarbonate (no drug)` = r$ref, `Bicarbonate (on drug)` = r$drug),
           "HCO3C", "Serum bicarbonate (mmol/L)",
           "Chronic metabolic acidosis — the common path to bone loss and stones",
           cols = c(PAL["grey"], PAL["drug"]), hline = 22)
  })
  output$t8_micro <- renderPlot({
    r <- simrun(); d <- r$ref; n <- nrow(d)
    v <- c(`Vitamin D` = d$VITDC[n], `Zinc (÷1400)` = d$ZNC[n]/1400,
           `B12 (÷3000)` = d$B12C[n]/3000, `Essential fatty acids (÷450)` = d$EFAC[n]/450,
           `PTH (÷45)` = d$PTHC[n]/45, `Bone mineral density` = d$BMDC[n])
    barplot2(v, "Micronutrient and bone status (1.0 = reference)", "Ratio to reference",
             cols = ifelse(v < 0.7, PAL["warn"], PAL["base"]))
  })

  # =================== TAB 9 ==========================================
  long_runs <- eventReactive(input$go, {
    st <- bstate(); p <- runpars(); w <- wt0()
    withProgress(message = "Computing the 2-year trajectories...", value = 0.5, {
      list(
        `Soybean-oil ILE` = sim(st, modifyList(p, list(PHYTOCONC = 350, OMEGA3 = 0)),
                          NULL, end = 730),
        `SMOF ILE` = sim(st, modifyList(p, list(PHYTOCONC = 100, OMEGA3 = 1)),
                          NULL, end = 730),
        `Fish-oil ILE` = sim(st, modifyList(p, list(PHYTOCONC = 0, OMEGA3 = 1)),
                          NULL, end = 730),
        `Soybean oil + teduglutide` = sim(st, modifyList(p, list(PHYTOCONC = 350,
                          OMEGA3 = 0)), ev_ted(w, 730), end = 730),
        `Soybean oil + taurolidine lock` = sim(st, modifyList(p, list(PHYTOCONC = 350,
                          OMEGA3 = 0, LOCKFAC = 0.25, TECHFAC = 0.7)),
                          NULL, end = 730))
    })
  }, ignoreNULL = FALSE)

  output$t9_liver <- renderPlot({
    tsplot(long_runs(), "BILIC", "Total bilirubin (mg/dL)",
           "IFALD: the composition of the lipid emulsion IS the dose", hline = 2)
  })
  output$t9_cvc <- renderPlot({
    tsplot(long_runs(), "VEINSC", "Usable central veins",
           "Access is a consumable resource")
  })
  output$t9_renal <- renderPlot({
    tsplot(long_runs(), "GFRC", "eGFR (mL/min/1.73m²)",
           "Renal function — the sum of chronic dehydration and oxalate load")
  })
  output$t9_bone <- renderPlot({
    tsplot(long_runs(), "UOXC", "Urine oxalate (mg/day)",
           "Enteric hyperoxaluria — it takes a colon to occur", hline = 40)
  })
  output$t9_tab <- kv_table({
    l <- long_runs()
    do.call(rbind, lapply(names(l), function(k) {
      d <- l[[k]]; n <- nrow(d)
      data.frame(`At 2 years` = k,
        `PN L/week` = round(d$PNVOLWK[n], 2),
        `Bilirubin` = round(d$BILIC[n], 2),
        `Fibrosis stage` = round(d$FIBC[n], 2),
        `Cumulative CRBSI` = round(d$CRBSICUM[n], 2),
        `Veins remaining` = round(d$VEINSC[n], 2),
        `eGFR` = round(d$GFRC[n], 1),
        `Stone burden` = round(d$STONEC[n], 3),
        `Bone mineral density` = round(d$BMDC[n], 3),
        `Survival %` = round(d$SURVP[n], 1),
        check.names = FALSE)
    }))
  })

  # =================== TAB 10 =========================================
  output$t10_loop <- renderPlot({
    r <- simrun(); d <- r$drug
    par(mar = c(4.4, 4.8, 3.2, 4.8))
    plot(d$time, d$URSM / r$p$URBASE, type = "l", lwd = 2.2, col = PAL["base"],
         las = 1, xlab = "days", ylab = "Urine output / baseline urine output",
         main = "Closed loop: taper when the urine output exceeds the trigger, stop when it does not",
         cex.main = 1.0, font.main = 1)
    grid(col = "#ececec", lty = 1)
    abline(h = input$urtrig, col = PAL["warn"], lty = 2, lwd = 2)
    abline(h = 1, col = "#999999", lty = 3)
    par(new = TRUE)
    plot(d$time, d$PNVOLWK, type = "l", lwd = 2.2, col = PAL["drug"],
         axes = FALSE, xlab = "", ylab = "")
    axis(4, las = 1); mtext("PN volume (L/week)", side = 4, line = 3.3, cex = 0.9)
    legend("topright", c("Urine output ratio (left)", "Taper trigger", "PN volume (right)"),
           col = c(PAL["base"], PAL["warn"], PAL["drug"]),
           lty = c(1, 2, 1), lwd = 2.2, bty = "n", cex = 0.84)
  })

  output$t10_placebo <- renderPlot({
    p <- basepars(); w <- NULL
    pre <- modifyList(p, list(ORSFRAC = 0.25))
    s <- burnin(pre, days = 700); wi <- s[["LEAN"]] + s[["FATM"]]
    q <- modifyList(p, list(URBASE = s[["URINESM"]], WTTARGET = wi))
    d <- max(168, input$days)
    a <- sim(s, modifyList(q, list(ORSFRAC = 0.60)), NULL, end = d)
    b <- sim(s, modifyList(q, list(ORSFRAC = 0.60)), ev_ted(wi, d), end = d)
    cc<- sim(s, modifyList(q, list(ORSFRAC = 0.25, WEANON = 0)), NULL, end = d)
    tsplot(list(`Placebo + protocolised care` = a, `Teduglutide + protocol` = b,
                `No protocol, no drug` = cc),
           "PNVOLWK", "PN volume (L/week)",
           "Where the placebo-arm response comes from: the protocol, not the placebo",
           cols = c(PAL["gold"], PAL["drug"], PAL["grey"]))
  })

  output$t10_trig <- renderPlot({
    st <- bstate(); p <- runpars(); w <- wt0(); d <- max(168, input$days)
    tr <- c(1.03, 1.05, 1.10, 1.20, 1.30)
    v <- vapply(tr, function(x) {
      dd <- sim(st, modifyList(p, list(URTRIG = x)), ev_ted(w, d), end = d)
      dd$PNVOLWK[nrow(dd)] - dd$PNVOLWK[1]
    }, numeric(1))
    names(v) <- sprintf("Trigger %.2f", tr)
    barplot2(v, "The taper-trigger threshold changes the measured drug effect", "Change in PN (L/week)",
             cols = PAL["purple"])
  })

  output$t10_tab <- kv_table({
    r <- simrun(); d <- r$drug; n <- nrow(d)
    data.frame(
      `Controller state` = c("Baseline urine output URBASE (L/day)", "Current smoothed urine output (L/day)",
        "Urine output ratio", "Taper trigger", "Re-escalation threshold", "Maximum taper rate (1/day)",
        "Weight guardrail (kg)", "Current weight (kg)", "Dehydration index", "Thirst response"),
      `Value` = round(c(r$p$URBASE, d$URSM[n], d$URSM[n]/r$p$URBASE,
        r$p$URTRIG, 0.85, r$p$KTAPER, 0.95 * r$p$WTTARGET, d$WT[n],
        d$DEHYDC[n], d$THIRSTC[n]), 3),
      check.names = FALSE)
  })
}

shinyApp(ui, server)
