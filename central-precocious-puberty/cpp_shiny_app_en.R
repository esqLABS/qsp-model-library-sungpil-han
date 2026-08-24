# =============================================================================
#  cpp_shiny_app.R
#  Central Precocious Puberty (CPP) — interactive QSP dashboard
#  ---------------------------------------------------------------------------
#  Central Precocious Puberty QSP Dashboard · 10 tabs
#
#  Run with:
#      shiny::runApp("cpp_shiny_app.R")
#  (the file sources cpp_mrgsolve_model_en.R from the same directory, which
#   compiles the 44-ODE mrgsolve model)
#
#  DESIGN PRINCIPLE OF THIS DASHBOARD
#  ----------------------------------
#  The app is built around the one thing the model exists to show: the SIGN of
#  the treatment effect on adult height is not a property of the drug, it is a
#  property of how much growth-plate reserve remains.  Every tab therefore
#  displays the treated trajectory AGAINST ITS OWN UNTREATED COUNTERFACTUAL
#  (same phenotype, no drug), because a treated trajectory on its own is
#  uninterpretable — and against the NORMAL reference girl/boy, because that is
#  the height the family is actually hoping for.
#
#  A second deliberate choice: the "Monitoring" tab plots growth velocity and
#  dBA/dCA side by side and labels which one is trustworthy.  Growth velocity
#  FALLS on effective therapy; a dashboard that shows it without that context
#  invites exactly the wrong clinical conclusion.
# =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

source("cpp_mrgsolve_model_en.R", local = TRUE, chdir = TRUE)

THEME <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92"),
        legend.position = "bottom")

PAL <- c(untreated = "#c0392b", treated = "#2471a3", normal = "#7f8c8d",
         alt = "#e67e22", alt2 = "#16a085")

# =============================================================================
#  UI
# =============================================================================
ui <- fluidPage(
  titlePanel("Central Precocious Puberty QSP Dashboard"),
  p(strong("The model's core proposition:"),
    "Estradiol enters the adult-height integral equation ", strong("twice, with opposite signs"),
    " — (+) it raises growth velocity via GH/IGF-1 amplification,",
    "(−) it consumes growth-plate reserve, ending the growth period.",
    "Because a GnRH agonist removes both arms simultaneously, whether treatment is beneficial is determined",
    strong("not by the drug but by the remaining growth-plate reserve"), "."),
  hr(),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("① Patient (Phenotype)"),
      radioButtons("sex", "Sex", c("Girl" = "F", "Boy" = "M"),
                   inline = TRUE),
      sliderInput("mk0", "MKRN3 brake remaining (MK0) — puberty onset is an outcome",
                  min = 0.20, max = 1.00, value = 0.40, step = 0.02),
      helpText("MK0 0.40 → thelarche at 6.5 yr · 1.00 → 10.3 yr · 0.66 → 8.6 yr (slowly progressive)"),
      sliderInput("ht0", "Height at age 5 (cm)", 96, 120, 108, step = 0.5),
      sliderInput("ba0", "Bone age at age 5 (yr)", 4.0, 7.0, 5.0, step = 0.1),
      checkboxInput("mas", "GnRH-independent (McCune-Albright, GNAS)", FALSE),
      checkboxInput("lesion", "Hypothalamic lesion/hamartoma (ectopic drive)", FALSE),

      h4("② Therapy"),
      selectInput("drug", "Formulation", c(
        "None (untreated)"                  = "none",
        "Leuprolide 3.75 mg IM q28d"        = "leup1m",
        "Leuprolide 7.5 mg IM q28d"         = "leup1mH",
        "Leuprolide 11.25 mg IM q12wk"      = "leup3m",
        "Leuprolide 30 mg IM q12wk"         = "leup3mH",
        "Triptorelin 11.25 mg IM q12wk"      = "trip3m",
        "Triptorelin 22.5 mg IM q24wk"       = "trip6m",
        "Histrelin 50 mg implant q12mo"   = "hist",
        "Nafarelin nasal spray 1800 µg/day"     = "naf",
        "GnRH antagonist 18 mg q28d"            = "antag"),
        selected = "leup3m"),
      sliderInput("start", "Age at treatment start (yr)", 5.5, 12.0, 7.4, step = 0.1),
      sliderInput("dur", "Treatment duration (yr)", 0.5, 8.0, 5.0, step = 0.5),
      sliderInput("late", "Actual dosing-interval delay multiplier (adherence proxy)",
                  1.0, 2.0, 1.0, step = 0.05),
      helpText("1.0 = exactly on schedule, 1.5 = a 28-day formulation given every 42 days"),

      h4("③ Add-on"),
      checkboxInput("ai", "Aromatase inhibitor (anastrozole 1 mg/day)", FALSE),
      checkboxInput("ai_potent", "  └ High-potency AI (letrozole-type)", FALSE),
      checkboxInput("gh", "rhGH 0.043 mg/kg/day", FALSE),
      checkboxInput("tam", "Tamoxifen 20 mg/day", FALSE),
      checkboxInput("cavd", "Calcium + vitamin D", FALSE),

      h4("④ Simulation"),
      sliderInput("years", "Simulation duration (yr from age 5)", 10, 18, 16.5,
                  step = 0.5),
      actionButton("go", "Run simulation", class = "btn-primary btn-block")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        # ---------------------------------------------------------------
        tabPanel(
          "1 · Patient profile",
          h4("The timing of puberty onset is an output, not an input"),
          p("Changing just the one MKRN3 brake (MK0) moves thelarche · menarche · bone age · final height",
            "all together. The table below shows the currently configured patient, that same patient's",
            "untreated control, and a normal reference child, side by side."),
          tableOutput("tbl_profile"),
          hr(),
          h4("MK0 sweep — phenotype spectrum"),
          plotOutput("plt_mk_sweep", height = "330px"),
          helpText("The slowly progressive variant (MK0 0.60-0.76) presents with early thelarche, but",
                   "the treatment benefit is very small — the quantitative form of the overtreatment problem.")
        ),
        # ---------------------------------------------------------------
        tabPanel(
          "2 · PK · Receptor",
          h4("The agonist does not block the receptor — it destroys the pulsatile code"),
          plotOutput("plt_pk", height = "620px"),
          helpText("fa = agonist occupancy, fx = antagonist occupancy, RS = fraction of sensitised receptors.",
                   "Because the agonist has an intrinsic activity of AINT = 1.6, stimulation rises above",
                   "baseline at the first dose (flare) — suppression only comes after RS collapses.",
                   "The antagonist has no flare, but also no desensitisation, so occupancy must be",
                   "maintained continuously.")
        ),
        # ---------------------------------------------------------------
        tabPanel(
          "3 · HPG axis",
          h4("KNDy pulse generator → pulse frequency → LH/FSH → gonadal steroids"),
          plotOutput("plt_hpg", height = "620px"),
          helpText("PULS (pulses/24h) is an explicit state variable. Continuous agonist exposure",
                   "drives the frequency information to zero — this is the principle of treatment.")
        ),
        # ---------------------------------------------------------------
        tabPanel(
          "4 · Growth · Bone age",
          h4("Two clocks: the height clock and the bone clock"),
          plotOutput("plt_growth", height = "640px"),
          helpText("Bone age is defined as 'maturation rate relative to a normal child of the same age'",
                   "(exactly the Greulich-Pyle atlas definition). So a normal child",
                   "by definition has ΔBA/ΔCA = 1.0, and a suppressed child appears to 'halt' at CA 10-12,",
                   "showing 0.6-0.7 — there is no halt term in the model.")
        ),
        # ---------------------------------------------------------------
        tabPanel(
          "5 · Clinical endpoints",
          h4("Final adult height, and the gain versus the untreated control"),
          tableOutput("tbl_endpoints"),
          hr(),
          plotOutput("plt_height", height = "380px"),
          helpText("BP-PAH = Bayley-Pinneau predicted adult height (the value computed clinically),",
                   "true final = the final height the model actually integrates. The difference between",
                   "the two values is exactly the bias of the clinical prediction formula.")
        ),
        # ---------------------------------------------------------------
        tabPanel(
          "6 · Sign flip (key result)",
          h4("The sign of the treatment benefit is set by the bone age at treatment start"),
          p("Nowhere in the curve below is a 'treat before age 8' rule coded in.",
            "The rule is an ", strong("output"), " that emerges from the competition between the (+) arm and the (−) arm."),
          plotOutput("plt_signflip", height = "420px"),
          tableOutput("tbl_signflip")
        ),
        # ---------------------------------------------------------------
        tabPanel(
          "7 · Monitoring (the trap)",
          h4("Growth velocity is an unreliable marker"),
          p("Effective suppression drives growth velocity even below the normal prepubertal value",
            "(because sex steroid and its GH/IGF-1 amplification are lost at the same time).",
            strong("A falling growth velocity is a sign of treatment success, not treatment failure.")),
          plotOutput("plt_monitor", height = "560px"),
          hr(),
          tableOutput("tbl_monitor")
        ),
        # ---------------------------------------------------------------
        tabPanel(
          "8 · Formulation comparison",
          h4("The determining variable is not potency but trough-interval coverage"),
          plotOutput("plt_forms", height = "420px"),
          tableOutput("tbl_forms"),
          helpText("Correctly administered depot formulations are almost indistinguishable from each other.",
                   "Differences arise with (a) nasal spray, (b) delayed injections, (c) low dose.")
        ),
        # ---------------------------------------------------------------
        tabPanel(
          "9 · Biomarkers",
          h4("Target tissue response — uterine volume · Tanner breast · endometrium"),
          plotOutput("plt_biomarker", height = "620px"),
          helpText("Pubic hair (adrenal androgen) is not suppressed by the GnRH agonist —",
                   "it continues to progress during treatment, so this must be mentioned when counselling caregivers.")
        ),
        # ---------------------------------------------------------------
        tabPanel(
          "10 · Safety · Quality of life",
          h4("The cost of treatment: bone density · body composition · hot flashes · psychosocial indices"),
          plotOutput("plt_safety", height = "620px"),
          helpText("The 'fall' in bone density Z-score is mostly an advanced baseline regressing",
                   "toward normal, and peak bone mass recovers.")
        )
      )
    )
  )
)

# =============================================================================
#  SERVER
# =============================================================================
server <- function(input, output, session) {

  base_pars <- reactive({
    male <- input$sex == "M"
    p <- if (male) pheno_boy(input$mk0) else pheno_girl(input$mk0)
    p$HT00 <- input$ht0
    p$BA00 <- input$ba0
    p$CAVD <- as.numeric(input$cavd)
    if (input$mas) { p$AUTSET <- 0.42; p$KSEC <- 0.22 }
    if (input$lesion) p$KNDLES <- 0.25
    if (input$ai_potent) { p$IMAXAI <- 0.992; p$IC50AI <- 2.2 }
    p
  })

  build_events <- reactive({
    a0 <- input$start; a1 <- input$start + input$dur
    L <- input$late
    ev <- switch(input$drug,
      none    = NULL,
      leup1m  = ev_leup(3.75, 28 * L, a0, a1),
      leup1mH = ev_leup(7.5, 28 * L, a0, a1),
      leup3m  = ev_leup(11.25, 84 * L, a0, a1),
      leup3mH = ev_leup(30, 84 * L, a0, a1),
      trip3m  = ev_trip(11.25, 84 * L, a0, a1),
      trip6m  = ev_trip(22.5, 168 * L, a0, a1),
      hist    = ev_hist(a0, a1),
      naf     = ev_naf(1800, a0, a1),
      antag   = ev_antag(18, 28 * L, a0, a1))
    if (input$ai)  ev <- c(ev, ev_ai(if (input$ai_potent) 2.5 else 1.0, a0, a1))
    if (input$gh)  ev <- c(ev, ev_gh(0.043 * 30, a0, a1))
    if (input$tam) ev <- c(ev, ev_tam(20, a0, a1))
    ev
  })

  sims <- eventReactive(input$go, {
    p <- base_pars()
    yrs <- input$years
    dT <- sim_cpp(p, build_events(), years = yrs, delta = 0.5)
    dU <- sim_cpp(p, NULL, years = yrs, delta = 0.5)
    pN <- p; pN$MK0 <- MK_NORMAL
    pN$AUTSET <- 0; pN$KSEC <- 0; pN$KNDLES <- 0
    dN <- sim_cpp(pN, NULL, years = yrs, delta = 0.5)
    list(treated = dT, untreated = dU, normal = dN, male = input$sex == "M")
  }, ignoreNULL = FALSE)

  long3 <- function(cols) {
    s <- sims()
    bind_rows(
      mutate(s$treated[, c("CA_out", cols)], arm = "treated"),
      mutate(s$untreated[, c("CA_out", cols)], arm = "untreated"),
      mutate(s$normal[, c("CA_out", cols)], arm = "normal")) %>%
      pivot_longer(all_of(cols), names_to = "var", values_to = "value")
  }

  panel3 <- function(cols, labs, logy = FALSE) {
    d <- long3(cols)
    d$var <- factor(d$var, levels = cols, labels = labs)
    g <- ggplot(d, aes(CA_out, value, colour = arm)) +
      geom_line(linewidth = 0.7) +
      facet_wrap(~var, scales = "free_y") +
      scale_colour_manual(values = PAL, name = NULL) +
      labs(x = "Chronological age (yr)", y = NULL) + THEME
    if (logy) g <- g + scale_y_log10()
    g
  }

  # ---- Tab 1 -----------------------------------------------------------
  output$tbl_profile <- renderTable({
    s <- sims()
    f <- function(d) c(
      `Thelarche/puberty onset (yr)` = thelarche(d),
      `Menarche (yr)` = menarche(d),
      `Bone maturation complete (yr)` = fusion_age(d),
      `Final adult height (cm)` = adult_height(d),
      `Peak growth velocity (cm/yr)` = max(d$GVo),
      `Peak E2 (pg/mL)` = max(d$E2),
      `Bone age at CA 8 (yr)` = at_age(d, 8, "BA"),
      `ΔBA/ΔCA at 7-9 yr` = ba_ca_ratio(d, 7, 9),
      `Cumulative E2 (pg/mL·yr)` = tail(d$CUME2, 1),
      `Adult BMD Z` = tail(d$BMDZ, 1))
    m <- cbind(`Current setting (treated)` = f(s$treated),
               `Same patient · untreated` = f(s$untreated),
               `Normal reference child` = f(s$normal))
    as.data.frame(round(m, 2))
  }, rownames = TRUE, digits = 2)

  output$plt_mk_sweep <- renderPlot({
    male <- input$sex == "M"
    mks <- seq(0.28, 0.92, by = 0.08)
    res <- do.call(rbind, lapply(mks, function(mk) {
      p <- if (male) pheno_boy(mk) else pheno_girl(mk)
      p$HT00 <- input$ht0; p$BA00 <- input$ba0
      dU <- sim_cpp(p, NULL, years = 17, delta = 2)
      a <- min(max(ifelse(is.na(thelarche(dU)), 12, thelarche(dU)) + 0.5, 6.3), 11.5)
      dT <- sim_cpp(p, ev_leup(11.25, 84, a, a + 5), years = 17, delta = 2)
      data.frame(MK0 = mk, onset = thelarche(dU),
                 noRx = adult_height(dU), Rx = adult_height(dT))
    }))
    res$gain <- res$Rx - res$noRx
    ggplot(res, aes(onset)) +
      geom_hline(yintercept = 0, linetype = 3) +
      geom_line(aes(y = gain), colour = PAL["treated"], linewidth = 1) +
      geom_point(aes(y = gain), colour = PAL["treated"], size = 2.4) +
      geom_line(aes(y = (noRx - min(res$noRx))), colour = PAL["untreated"],
                linetype = 2) +
      annotate("text", x = max(res$onset, na.rm = TRUE), y = 0.4,
               label = "Treatment benefit = 0", hjust = 1, size = 3.4) +
      labs(x = "Age at puberty onset if untreated (yr)",
           y = "Adult height gain (cm, solid line) / relative untreated final height (dashed line)") + THEME
  })

  # ---- Tab 2 -----------------------------------------------------------
  output$plt_pk <- renderPlot({
    s <- sims()
    cols <- c("CPTOT", "CANTo", "FAo", "FXo", "RS", "So")
    labs <- c("Agonist concentration (ng/mL)", "Antagonist concentration (ng/mL)",
              "Agonist occupancy fa", "Antagonist occupancy fx",
              "Sensitised receptor RS", "Gonadotrope stimulation S")
    d <- bind_rows(mutate(s$treated[, c("CA_out", cols)], arm = "treated"),
                   mutate(s$untreated[, c("CA_out", cols)], arm = "untreated")) %>%
      pivot_longer(all_of(cols), names_to = "var", values_to = "value")
    d$var <- factor(d$var, levels = cols, labels = labs)
    ggplot(d, aes(CA_out, value, colour = arm)) +
      geom_line(linewidth = 0.6) +
      facet_wrap(~var, scales = "free_y", ncol = 2) +
      scale_colour_manual(values = PAL, name = NULL) +
      labs(x = "Chronological age (yr)", y = NULL) + THEME
  })

  # ---- Tab 3 -----------------------------------------------------------
  output$plt_hpg <- renderPlot({
    panel3(c("KND", "PULS", "LH", "FSH", "E2", "INHB"),
           c("KNDy drive (-)", "GnRH pulse frequency (pulses/24h)",
             "LH (IU/L)", "FSH (IU/L)", "Estradiol (pg/mL)",
             "Inhibin B (pg/mL)"))
  })

  # ---- Tab 4 -----------------------------------------------------------
  output$plt_growth <- renderPlot({
    s <- sims()
    cols <- c("HT", "GVo", "IGF1", "GPRES", "BA", "DBADCA")
    labs <- c("Height (cm)", "Growth velocity (cm/yr)", "IGF-1 (ng/mL)",
              "Growth-plate reserve GPRES (-)", "Bone age (yr)", "ΔBA/ΔCA (-)")
    d <- long3(cols); d$var <- factor(d$var, levels = cols, labels = labs)
    ref <- data.frame(var = factor("Bone age (yr)", levels = labs),
                      x = range(s$treated$CA_out))
    ggplot(d, aes(CA_out, value, colour = arm)) +
      geom_line(linewidth = 0.7) +
      geom_abline(data = ref, aes(slope = 1, intercept = 0), linetype = 3,
                  inherit.aes = FALSE) +
      facet_wrap(~var, scales = "free_y", ncol = 2) +
      scale_colour_manual(values = PAL, name = NULL) +
      labs(x = "Chronological age (yr)", y = NULL,
           caption = "The dotted line in the bone-age panel is BA = CA (the normal reference line)") + THEME
  })

  # ---- Tab 5 -----------------------------------------------------------
  output$tbl_endpoints <- renderTable({
    s <- sims()
    male <- s$male
    f <- function(d) c(
      `Final adult height (cm)` = adult_height(d),
      `Age at bone maturation complete (yr)` = fusion_age(d),
      `Age at menarche (yr)` = menarche(d),
      `Final Tanner breast stage` = tail(d$BST, 1),
      `BP-PAH at treatment end (cm)` = pah_bp(at_age(d, min(input$start + input$dur,
                                                       max(d$CA_out)), "HT"),
                                          at_age(d, min(input$start + input$dur,
                                                        max(d$CA_out)), "BA"),
                                          male),
      `Peak psychosocial index (0-10)` = max(d$QOL),
      `Peak hot-flash index (0-10)` = max(d$HF),
      `BMD Z nadir` = min(d$BMDZ))
    m <- cbind(`Treated` = f(s$treated), `Untreated` = f(s$untreated),
               `Normal reference` = f(s$normal))
    out <- as.data.frame(round(m, 2))
    out["Adult height gain (cm)", ] <- c(
      round(adult_height(s$treated) - adult_height(s$untreated), 2), 0,
      round(adult_height(s$normal) - adult_height(s$untreated), 2))
    out
  }, rownames = TRUE, digits = 2)

  output$plt_height <- renderPlot({
    s <- sims()
    d <- long3(c("HT"))
    bp <- bind_rows(
      data.frame(CA_out = s$treated$CA_out, arm = "treated",
                 value = pah_bp(s$treated$HT, s$treated$BA, s$male)),
      data.frame(CA_out = s$untreated$CA_out, arm = "untreated",
                 value = pah_bp(s$untreated$HT, s$untreated$BA, s$male)))
    ggplot() +
      geom_line(data = d, aes(CA_out, value, colour = arm), linewidth = 0.8) +
      geom_line(data = bp, aes(CA_out, value, colour = arm), linetype = 2) +
      scale_colour_manual(values = PAL, name = NULL) +
      labs(x = "Chronological age (yr)", y = "Height (solid line) / BP-predicted adult height (dashed line), cm") +
      THEME
  })

  # ---- Tab 6 -----------------------------------------------------------
  signflip <- reactive({
    p <- base_pars()
    dU <- sim_cpp(p, NULL, years = 17, delta = 2)
    hU <- adult_height(dU)
    starts <- seq(6.4, 11.2, by = 0.4)
    do.call(rbind, lapply(starts, function(a) {
      dT <- sim_cpp(p, ev_leup(11.25, 84, a, a + 6), years = 17, delta = 2)
      data.frame(start_CA = a, start_BA = at_age(dU, a, "BA"),
                 GPRES = at_age(dU, a, "GPRES"),
                 gain = adult_height(dT) - hU,
                 GV_on_Rx = gv_mean(dT, a + 0.3, min(a + 2, 16.4)))
    }))
  })

  output$plt_signflip <- renderPlot({
    r <- signflip()
    ggplot(r, aes(start_BA, gain)) +
      geom_hline(yintercept = 0, linetype = 3) +
      geom_line(linewidth = 1, colour = PAL["treated"]) +
      geom_point(size = 2.4, colour = PAL["treated"]) +
      geom_vline(xintercept = 12, linetype = 2, colour = PAL["untreated"]) +
      annotate("text", x = 12.05, y = max(r$gain), hjust = 0, size = 3.6,
               label = "Clinical consensus: no benefit after bone age 12") +
      labs(x = "Bone age at treatment start (yr)",
           y = "Adult height gain vs untreated (cm)") + THEME
  })

  output$tbl_signflip <- renderTable({ signflip() }, digits = 2)

  # ---- Tab 7 -----------------------------------------------------------
  output$plt_monitor <- renderPlot({
    s <- sims()
    cols <- c("GVo", "DBADCA", "LH", "UTV")
    labs <- c("Growth velocity (cm/yr) — unreliable",
              "ΔBA/ΔCA — reliable", "Basal LH (IU/L)",
              "Uterine volume (mL)")
    d <- long3(cols); d$var <- factor(d$var, levels = cols, labels = labs)
    ggplot(d, aes(CA_out, value, colour = arm)) +
      geom_line(linewidth = 0.7) +
      facet_wrap(~var, scales = "free_y", ncol = 2) +
      scale_colour_manual(values = PAL, name = NULL) +
      labs(x = "Chronological age (yr)", y = NULL) + THEME
  })

  output$tbl_monitor <- renderTable({
    s <- sims()
    lo <- input$start + 0.35; hi <- input$start + input$dur
    f <- function(d, col) mean(d[[col]][d$CA_out >= lo & d$CA_out <= hi])
    data.frame(
      surrogate = c("Basal LH (IU/L)", "E2 (pg/mL)", "Uterine volume (mL)",
                    "Growth velocity (cm/yr)", "ΔBA/ΔCA", "Tanner breast", "IGF-1"),
      treated = c(f(s$treated, "LH"), f(s$treated, "E2"), f(s$treated, "UTV"),
                  f(s$treated, "GVo"), f(s$treated, "DBADCA"),
                  f(s$treated, "BST"), f(s$treated, "IGF1")),
      untreated = c(f(s$untreated, "LH"), f(s$untreated, "E2"),
                    f(s$untreated, "UTV"), f(s$untreated, "GVo"),
                    f(s$untreated, "DBADCA"), f(s$untreated, "BST"),
                    f(s$untreated, "IGF1")),
      check.names = FALSE)
  }, digits = 3)

  # ---- Tab 8 -----------------------------------------------------------
  forms <- reactive({
    p <- base_pars()
    a0 <- input$start; a1 <- input$start + input$dur
    regs <- list(
      "Leuprolide 3.75 q28d"   = ev_leup(3.75, 28, a0, a1),
      "Leuprolide 11.25 q12wk" = ev_leup(11.25, 84, a0, a1),
      "Triptorelin 22.5 q24wk"  = ev_trip(22.5, 168, a0, a1),
      "Histrelin implant"    = ev_hist(a0, a1),
      "Nafarelin nasal spray"      = ev_naf(1800, a0, a1),
      "3.75 mg every 42 days (delayed)" = ev_leup(3.75, 42, a0, a1),
      "1.875 mg q28d (low dose)"  = ev_leup(1.875, 28, a0, a1),
      "GnRH antagonist 18 mg"      = ev_antag(18, 28, a0, a1))
    hU <- adult_height(sim_cpp(p, NULL, years = 16.5, delta = 2))
    list(hU = hU, tab = do.call(rbind, lapply(names(regs), function(nm) {
      d <- sim_cpp(p, regs[[nm]], years = 16.5, delta = 1)
      data.frame(formulation = nm,
                 mean_Cp = mean(d$CPTOT[d$CA_out >= a0 + 0.25 & d$CA_out <= a1]),
                 mean_fa = mean(d$FAo[d$CA_out >= a0 + 0.25 & d$CA_out <= a1]),
                 pct_LH_esc = pct_time(d, a0 + 0.25, a1, "LH", 0.5),
                 pct_E2_esc = pct_time(d, a0 + 0.25, a1, "E2", 10),
                 adult_height = adult_height(d),
                 gain = adult_height(d) - hU)
    })))
  })

  output$plt_forms <- renderPlot({
    r <- forms()$tab
    ggplot(r, aes(reorder(formulation, gain), gain, fill = pct_E2_esc)) +
      geom_col() + coord_flip() +
      scale_fill_gradient(low = "#2471a3", high = "#c0392b",
                          name = "Time with E2 > 10 pg/mL (%)") +
      labs(x = NULL, y = "Adult height gain vs untreated (cm)") + THEME
  })

  output$tbl_forms <- renderTable({ forms()$tab }, digits = 3)

  # ---- Tab 9 -----------------------------------------------------------
  output$plt_biomarker <- renderPlot({
    panel3(c("UTV", "BST", "ENDO", "DHEAS", "TESTO", "FOL"),
           c("Uterine volume (mL)", "Tanner breast stage (1-5)",
             "Endometrial status (-)", "DHEAS (µg/dL) — not suppressed",
             "Testosterone (ng/dL)", "Follicular functional mass (-)"))
  })

  # ---- Tab 10 ----------------------------------------------------------
  output$plt_safety <- renderPlot({
    panel3(c("BMDZ", "BMIZ", "HF", "QOL", "GH", "CUME2"),
           c("Lumbar BMD Z-score", "BMI Z-score", "Hot-flash index (0-10)",
             "Psychosocial distress index (0-10)", "GH pulse amplitude (units)",
             "Cumulative E2 exposure (pg/mL·yr)"))
  })
}

shinyApp(ui, server)
