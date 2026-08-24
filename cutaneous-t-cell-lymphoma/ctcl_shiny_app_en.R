## =============================================================================
##  ctcl_shiny_app.R
##  Cutaneous T-cell lymphoma (mycosis fungoides / Sézary syndrome)
##  Interactive dashboard for the QSP model in ctcl_mrgsolve_model.R
##
##  The app is organised around the model's one claim: the clone lives in three
##  boxes, each therapy reaches a different subset of them, and each response
##  criterion reads only one of them.  Tab 3 (compartments) and tab 6 (global
##  response) are therefore the two tabs that matter; the rest support them.
##
##  Run with:  shiny::runApp("ctcl_shiny_app.R")
## =============================================================================

suppressPackageStartupMessages({
  library(shiny)
  library(mrgsolve)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(DT)
})

source("ctcl_mrgsolve_model.R", local = TRUE)

THEME <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(colour = "grey35", size = 10))

PAL <- c(skin = "#c0392b", blood = "#1f6fb2", node = "#7f8c8d",
         resistant = "#e07a80", transformed = "#7b241c",
         surveillance = "#1f6fb2", drug = "#3d8b40")

REGIMENS <- c("None (natural history)", "Narrowband UVB + topical steroid",
              "Chlormethine gel", "Low-dose TSEB 12 Gy", "Bexarotene 300 mg/m²",
              "Interferon alfa 3 MU tiw", "Bexarotene + interferon",
              "Vorinostat 400 mg qd", "Romidepsin 14 mg/m²",
              "Mogamulizumab 1 mg/kg", "Mogamulizumab 3 mg/kg",
              "Brentuximab vedotin ×16", "Romidepsin ×2 → brentuximab",
              "Gemcitabine ×6", "Low-dose methotrexate",
              "Low-dose alemtuzumab SC", "ECP q4w", "ECP + interferon",
              "Mogamulizumab + ECP", "Anti-staphylococcal antibiotic (d90–100)")

build_ev <- function(name, wt, bsa, dur) {
  switch(name,
    "None (natural history)"          = NULL,
    "Narrowband UVB + topical steroid" = c(ev_nbuvb(min(dur, 180)), ev_steroid(min(dur, 180))),
    "Chlormethine gel"                 = ev_chlormethine(dur),
    "Low-dose TSEB 12 Gy"              = ev_tseb(12, 8),
    "Bexarotene 300 mg/m²"             = ev_bexarotene(bsa, dur),
    "Interferon alfa 3 MU tiw"         = ev_interferon(dur),
    "Bexarotene + interferon"          = c(ev_bexarotene(bsa, dur), ev_interferon(dur)),
    "Vorinostat 400 mg qd"             = ev_vorinostat(dur),
    "Romidepsin 14 mg/m²"              = ev_romidepsin(bsa, ceiling(dur / 28)),
    "Mogamulizumab 1 mg/kg"            = ev_mogamulizumab(wt, dur),
    "Mogamulizumab 3 mg/kg"            = c(ev(time = 0, amt = 3 * wt, cmt = "MOG1", ii = 7, addl = 3),
                                        ev(time = 28, amt = 3 * wt, cmt = "MOG1",
                                           ii = 14, addl = max(0, floor((dur - 28) / 14)))),
    "Brentuximab vedotin ×16"          = ev_brentuximab(wt, 16),
    "Romidepsin ×2 → brentuximab"      = c(ev_romidepsin(bsa, 2),
                                        ev(time = 56, amt = 1.8 * wt, cmt = "BV1",
                                           ii = 21, addl = 12)),
    "Gemcitabine ×6"                   = ev_gemcitabine(bsa, 6),
    "Low-dose methotrexate"            = ev_methotrexate(dur),
    "Low-dose alemtuzumab SC"          = ev_alemtuzumab(min(dur, 180)),
    "ECP q4w"                       = ev_ecp(dur),
    "ECP + interferon"                 = c(ev_ecp(dur), ev_interferon(dur)),
    "Mogamulizumab + ECP"              = c(ev_mogamulizumab(wt, dur), ev_ecp(dur)),
    "Anti-staphylococcal antibiotic (d90–100)" = ev_antibiotic(90, 10))
}

## =============================================================================
##  UI
## =============================================================================
ui <- fluidPage(
  titlePanel("Cutaneous T-cell lymphoma (MF / Sézary) QSP dashboard — one clone, three compartments, two host variables"),
  tags$p(style = "color:#555;margin-top:-8px;",
         "Cutaneous T-cell lymphoma QSP dashboard · mSWAT reads only the skin, B-score reads only the blood. ",
         "The GLOBAL response is the logical AND of the two."),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Patient"),
      selectInput("stage", "Stage (TNMB stage)",
                  c("IA", "IB", "IIA", "IIB", "IIB_LCT", "IIIA", "IIIB",
                    "IVA1", "IVA2", "IVB"), selected = "IIB"),
      radioButtons("kind", "Phenotype (residency phenotype)",
                   c("MF (CD69 high · skin-resident)" = "MF",
                     "Sézary (CD69 low · recirculating)" = "SS"), selected = "MF"),
      sliderInput("cd69", "CD69 (tissue residency)", 0.05, 0.97, 0.85, 0.01),
      sliderInput("ccr7", "CCR7 (recirculation programme)", 0.02, 0.98, 0.10, 0.01),
      sliderInput("ccr4", "CCR4-positive fraction", 0.10, 1.0, 0.85, 0.05),
      sliderInput("cd30", "CD30 expression", 0.0, 0.95, 0.10, 0.05),
      sliderInput("fres", "Drug-resistant clone fraction FRES", 0.01, 0.90, 0.30, 0.01),
      sliderInput("sensf", "Drug sensitivity multiplier", 0.2, 3.0, 1.0, 0.1),
      sliderInput("nkf", "NK density multiplier", 0.2, 2.5, 1.0, 0.1),
      sliderInput("deficit", "Surveillance deficit (disease aggressiveness)", 0.0, 0.45, 0.25, 0.01),
      numericInput("wt", "Weight (kg)", 75, 40, 140, 1),
      numericInput("bsa", "Body surface area (m²)", 1.85, 1.3, 2.5, 0.05),
      hr(),
      h4("Treatment (Regimen)"),
      selectInput("reg", "Regimen A", REGIMENS, selected = "Mogamulizumab 1 mg/kg"),
      selectInput("reg2", "Regimen B (comparison)", REGIMENS, selected = "Vorinostat 400 mg qd"),
      sliderInput("dur", "Treatment duration (days)", 84, 730, 546, 7),
      sliderInput("end", "Simulation duration (days)", 180, 1095, 730, 30),
      actionButton("go", "Run simulation", class = "btn-primary btn-block")
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("① Patient profile",
                 br(), htmlOutput("profile"), br(),
                 plotOutput("p_profile", height = "330px"),
                 DT::dataTableOutput("t_profile")),
        tabPanel("② Pharmacokinetics (PK)",
                 br(), plotOutput("p_pk", height = "330px"),
                 h5("Tissue exposure and receptor occupancy — what's lacking in the skin is not the drug"),
                 plotOutput("p_occ", height = "280px")),
        tabPanel("③ Tumour burden by compartment",
                 br(), plotOutput("p_clone", height = "380px"),
                 h5("Sensitive · resistant subclones"),
                 plotOutput("p_subclone", height = "280px")),
        tabPanel("④ Skin endpoint (mSWAT)",
                 br(), plotOutput("p_mswat", height = "340px"),
                 h5("Differing clearance time constants for patch · plaque · tumour"),
                 plotOutput("p_morph", height = "280px")),
        tabPanel("⑤ Blood · lymph node",
                 br(), plotOutput("p_blood", height = "340px"),
                 plotOutput("p_node", height = "250px")),
        tabPanel("⑥ GLOBAL response (Olsen logical AND)",
                 br(), htmlOutput("respbox"),
                 plotOutput("p_resp", height = "300px"),
                 DT::dataTableOutput("t_resp")),
        tabPanel("⑦ Surveillance variables · durability",
                 br(),
                 tags$p("ORR reads 'how much was killed', while TTNT is set by 'who remains to keep watch'."),
                 plotOutput("p_surv", height = "340px"),
                 plotOutput("p_ttnt", height = "260px")),
        tabPanel("⑧ Pruritus · quality of life",
                 br(), plotOutput("p_itch", height = "330px"),
                 h5("Decoupling of tumour burden and pruritus"),
                 plotOutput("p_itchdecouple", height = "280px")),
        tabPanel("⑨ Serum biomarkers",
                 br(), plotOutput("p_biomk", height = "420px")),
        tabPanel("⑩ Barrier · staphylococcus · infection",
                 br(), plotOutput("p_barrier", height = "340px"),
                 plotOutput("p_infect", height = "260px")),
        tabPanel("⑪ Toxicity · safety",
                 br(), plotOutput("p_tox", height = "420px"),
                 DT::dataTableOutput("t_tox")),
        tabPanel("⑫ Scenario comparison",
                 br(), plotOutput("p_cmp", height = "380px"),
                 DT::dataTableOutput("t_cmp")),
        tabPanel("⑬ Structural experiments",
                 br(),
                 tags$p(strong("E1"), " — Mogamulizumab's compartment gap is due to effector cells, not exposure."),
                 DT::dataTableOutput("t_e1"), br(),
                 tags$p(strong("E6"), " — Moving CD69 alone, from MF to Sézary."),
                 plotOutput("p_e6", height = "280px"))
      )
    )
  )
)

## =============================================================================
##  SERVER
## =============================================================================
server <- function(input, output, session) {

  observeEvent(input$kind, {
    if (input$kind == "MF") {
      updateSliderInput(session, "cd69", value = 0.85)
      updateSliderInput(session, "ccr7", value = 0.10)
    } else {
      updateSliderInput(session, "cd69", value = 0.20)
      updateSliderInput(session, "ccr7", value = 0.85)
    }
  })

  covs <- reactive(list(CD69 = input$cd69, CCR7 = input$ccr7, CCR4 = input$ccr4,
                        CD30 = input$cd30, FRES = input$fres, SENSF = input$sensf,
                        NKF = input$nkf))

  patient <- eventReactive(input$go, {
    withProgress(message = "Equilibrating patient (stage = the surveillance switch's fixed point)", value = 0.3, {
      ctcl_patient(input$stage, input$kind, input$deficit, covs())
    })
  }, ignoreNULL = FALSE)

  simA <- eventReactive(input$go, {
    pt <- patient(); e <- build_ev(input$reg, input$wt, input$bsa, input$dur)
    d <- if (is.null(e)) mrgsim(pt$mod, end = input$end, delta = 1)
         else mrgsim(pt$mod, events = e, end = input$end, delta = 1)
    as.data.frame(d)
  }, ignoreNULL = FALSE)

  simB <- eventReactive(input$go, {
    pt <- patient(); e <- build_ev(input$reg2, input$wt, input$bsa, input$dur)
    d <- if (is.null(e)) mrgsim(pt$mod, end = input$end, delta = 1)
         else mrgsim(pt$mod, events = e, end = input$end, delta = 1)
    as.data.frame(d)
  }, ignoreNULL = FALSE)

  both <- reactive(bind_rows(mutate(simA(), arm = input$reg),
                             mutate(simB(), arm = input$reg2)))

  ## ---- ① profile ----------------------------------------------------------
  output$profile <- renderUI({
    pt <- patient(); x <- simA()
    HTML(sprintf(
      "<div style='background:#f6f6f9;padding:10px;border-radius:6px'>
       <b>Stage</b> %s · <b>Phenotype</b> %s ·
       <b>Equilibrium immune fitness IMMF<sub>eq</sub></b> %.2f → <b>applied value</b> %.2f
       (surveillance deficit %.0f%%)<br>
       <b>Baseline</b> mSWAT %.1f · Sézary %.0f/µL (B%d) · N-score %d ·
       Pruritus NRS %.1f · TARC %.0f pg/mL<br>
       <span style='color:#666'>Stage is not an initial value but is computed as the surveillance switch's fixed point,
       and the skin·blood split is <i>derived</i> from CD69/CCR7.</span></div>",
      pt$stage, pt$kind, pt$immf_eq, pt$immf, 100 * input$deficit,
      x$MSWAT[1], x$SEZCT[1], as.integer(x$BSCORE[1]), as.integer(x$NSCORE[1]),
      x$PRUR[1], x$TARC[1]))
  })

  output$p_profile <- renderPlot({
    x <- simA()[1, ]
    d <- data.frame(
      compartment = factor(c("Skin sensitive", "Skin resistant", "Transformed subclone",
                             "Blood sensitive", "Blood resistant", "Lymph node", "Viscera"),
                           levels = c("Skin sensitive", "Skin resistant", "Transformed subclone",
                                      "Blood sensitive", "Blood resistant", "Lymph node", "Viscera")),
      cells = c(x$NSK, x$NSKR, x$NTR, x$NBL, x$NBLR, x$NLN, x$NVS))
    ggplot(d, aes(compartment, cells, fill = compartment)) +
      geom_col(width = 0.65, show.legend = FALSE) +
      labs(title = "Baseline clone distribution (10⁹ cells)",
           subtitle = "The MF/SS difference comes from the CD69·CCR7 values, not the equations",
           x = NULL, y = "10⁹ cells") + THEME
  })

  output$t_profile <- DT::renderDataTable({
    x <- simA()[1, ]
    DT::datatable(data.frame(
      Item = c("mSWAT", "Patch %BSA", "Plaque %BSA", "Tumour %BSA", "Sézary /µL",
               "B-score", "N-score", "M-score", "T-stage", "Pruritus NRS",
               "TARC pg/mL", "sIL-2R U/mL", "LDH U/L", "Barrier integrity", "CD4 /µL"),
      Value = round(c(x$MSWAT, x$APAT, x$APLQ, x$ATUM, x$SEZCT, x$BSCORE, x$NSCORE,
                   x$MSCORE, x$TSKIN, x$PRUR, x$TARC, x$SIL2R, x$LDH, x$BARR,
                   x$CD4N), 2)),
      options = list(dom = "t", pageLength = 20), rownames = FALSE)
  })

  ## ---- ② PK ---------------------------------------------------------------
  output$p_pk <- renderPlot({
    x <- simA()
    d <- x %>% select(time, `Mogamulizumab plasma (µg/mL)` = CMOGP,
                      `Mogamulizumab skin (µg/mL)` = CMOGSKO,
                      `BV ADC plasma (µg/mL)` = CBVP,
                      `Free MMAE (ng/mL)` = CMMP) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value)) + geom_line(colour = PAL[["drug"]], linewidth = 0.7) +
      facet_wrap(~name, scales = "free_y") +
      labs(title = "Pharmacokinetics", subtitle = "Skin concentration = 0.157 × plasma (Shah & Betts 2013, an ill-fitting constant)",
           x = "Day", y = NULL) + THEME
  })

  output$p_occ <- renderPlot({
    x <- simA()
    d <- x %>% select(time, `Blood CCR4 occupancy` = OCCBL, `Skin CCR4 occupancy` = OCCSK) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 0.8) +
      scale_y_continuous(limits = c(0, 1)) +
      labs(title = "CCR4 receptor occupancy", x = "Day", y = "Occupancy", colour = NULL,
           subtitle = "Skin occupancy is also near 1 — why raising the dose does not move the skin response") + THEME
  })

  ## ---- ③ compartments -----------------------------------------------------
  output$p_clone <- renderPlot({
    d <- both() %>%
      transmute(time, arm, Skin = NSK + NSKR + NTR, Blood = NBL + NBLR,
                Node = NLN) %>% pivot_longer(c(Skin, Blood, Node))
    ggplot(d, aes(time, value, colour = name, linetype = arm)) +
      geom_line(linewidth = 0.8) + scale_y_log10() +
      scale_colour_manual(values = c(Skin = PAL[["skin"]], Blood = PAL[["blood"]],
                                     Node = PAL[["node"]])) +
      labs(title = "Clone burden by compartment", subtitle = "The same drug acts at a different rate in each compartment",
           x = "Day", y = "10⁹ cells (log)", colour = NULL, linetype = NULL) + THEME
  })

  output$p_subclone <- renderPlot({
    d <- simA() %>% transmute(time, Sensitive = NSK, Resistant = NSKR, Transformed = NTR) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value, fill = name)) + geom_area(alpha = 0.85) +
      scale_fill_manual(values = c(Sensitive = PAL[["skin"]], Resistant = PAL[["resistant"]],
                                   Transformed = PAL[["transformed"]])) +
      labs(title = "Subpopulations of the skin clone", x = "Day", y = "10⁹ cells", fill = NULL,
           subtitle = "The resistant subclone resists the drug but not CD8 surveillance") + THEME
  })

  ## ---- ④ mSWAT ------------------------------------------------------------
  output$p_mswat <- renderPlot({
    d <- both()
    b <- d %>% group_by(arm) %>% slice(1) %>% transmute(arm, base = MSWAT)
    d <- left_join(d, b, by = "arm")
    ggplot(d, aes(time, MSWAT, colour = arm)) + geom_line(linewidth = 0.9) +
      geom_hline(aes(yintercept = 0.5 * base, colour = arm), linetype = "dashed",
                 alpha = 0.5) +
      labs(title = "mSWAT (skin compartment score)", x = "Day", y = "mSWAT", colour = NULL,
           subtitle = "Dashed line = each arm's 50% reduction line (skin partial-response criterion)") + THEME
  })

  output$p_morph <- renderPlot({
    d <- simA() %>% transmute(time, Patch = APAT, Plaque = APLQ, Tumour = ATUM) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 0.8) +
      labs(title = "Area by lesion morphology", x = "Day", y = "%BSA", colour = NULL,
           subtitle = "τ: patch 8 days · plaque 25 days · tumour 12 days — why response confirmation needs 4 weeks") + THEME
  })

  ## ---- ⑤ blood / node -----------------------------------------------------
  output$p_blood <- renderPlot({
    ggplot(both(), aes(time, SEZCT, colour = arm)) + geom_line(linewidth = 0.9) +
      geom_hline(yintercept = c(250, 1000), linetype = "dotted", colour = "grey40") +
      scale_y_log10() +
      labs(title = "Circulating Sézary cell count", x = "Day", y = "cells/µL (log)", colour = NULL,
           subtitle = "Dashed lines = B1 (250) and B2 (1000) boundaries") + THEME
  })
  output$p_node <- renderPlot({
    ggplot(both(), aes(time, NLN, colour = arm)) + geom_line(linewidth = 0.9) +
      labs(title = "Lymph node clone burden", x = "Day", y = "10⁹ cells", colour = NULL) + THEME
  })

  ## ---- ⑥ global response --------------------------------------------------
  respA <- reactive(ctcl_endpoints(simA()))
  respB <- reactive(ctcl_endpoints(simB()))

  output$respbox <- renderUI({
    a <- respA()
    HTML(sprintf(
      "<div style='background:#f3eefa;padding:10px;border-radius:6px'>
       <b>%s</b> — best mSWAT reduction <b>%.1f%%</b> ·
       Skin response day %s · blood response day %s · <b>GLOBAL response day %s</b> ·
       PFS %s days · TTNT %s days<br>
       <span style='color:#666'>GLOBAL requires the response to hold simultaneously in every involved compartment.
       GLOBAL does not move if only one compartment improves.</span></div>",
      input$reg, a$best_skin_pct,
      ifelse(is.na(a$skin_response_day), "—", a$skin_response_day),
      ifelse(is.na(a$blood_response_day), "—", a$blood_response_day),
      ifelse(is.na(a$global_response_day), "—", a$global_response_day),
      ifelse(is.na(a$PFS_day), "Not reached", a$PFS_day),
      ifelse(is.na(a$TTNT_day), "Not reached", a$TTNT_day)))
  })

  output$p_resp <- renderPlot({
    x <- simA()
    b <- c(x$MSWAT[1], x$SEZCT[1], x$NLN[1])
    d <- data.frame(time = x$time,
                    Skin = 100 * (1 - x$MSWAT / b[1]),
                    Blood = 100 * (1 - x$SEZCT / max(b[2], 1e-9)),
                    Node = 100 * (1 - x$NLN / max(b[3], 1e-9))) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      geom_hline(yintercept = 50, linetype = "dashed", colour = "grey40") +
      scale_colour_manual(values = c(Skin = PAL[["skin"]], Blood = PAL[["blood"]],
                                     Node = PAL[["node"]])) +
      labs(title = "Reduction from baseline by compartment", x = "Day", y = "% reduction", colour = NULL,
           subtitle = "Dashed line = 50% (partial response). GLOBAL holds only if all three lines clear it") + THEME
  })

  output$t_resp <- DT::renderDataTable({
    DT::datatable(bind_rows(cbind(arm = input$reg, respA()),
                            cbind(arm = input$reg2, respB())),
                  options = list(dom = "t", scrollX = TRUE), rownames = FALSE)
  })

  ## ---- ⑦ surveillance -----------------------------------------------------
  output$p_surv <- renderPlot({
    d <- both() %>% group_by(arm) %>%
      transmute(time, `Skin CD8 surveillance` = SURVSK / first(SURVSK),
                `Regulatory T cells` = TREG / first(TREG),
                `NK cells` = NKB / first(NKB),
                `DC activity` = DCA / first(DCA)) %>% ungroup() %>%
      pivot_longer(-c(time, arm))
    ggplot(d, aes(time, value, colour = name, linetype = arm)) +
      geom_line(linewidth = 0.8) + geom_hline(yintercept = 1, colour = "grey70") +
      labs(title = "Second state variable — surveillance (relative to baseline)", x = "Day", y = "Fold",
           colour = NULL, linetype = NULL,
           subtitle = "Chemotherapy erases both the clone and surveillance together. ECP/interferon raise only surveillance.") + THEME
  })

  output$p_ttnt <- renderPlot({
    d <- both() %>% group_by(arm) %>%
      transmute(time, mSWAT = MSWAT / first(MSWAT) * 100) %>% ungroup()
    ggplot(d, aes(time, mSWAT, colour = arm)) + geom_line(linewidth = 0.9) +
      geom_hline(yintercept = 75, linetype = "dashed", colour = "grey40") +
      labs(title = "Time-to-next-treatment (TTNT) baseline", x = "Day",
           y = "mSWAT relative to baseline (%)", colour = NULL,
           subtitle = "Dashed line: time to re-reach 75% = TTNT") + THEME
  })

  ## ---- ⑧ pruritus ---------------------------------------------------------
  output$p_itch <- renderPlot({
    ggplot(both(), aes(time, PRUR, colour = arm)) + geom_line(linewidth = 0.9) +
      scale_y_continuous(limits = c(0, 10)) +
      labs(title = "Pruritus NRS", x = "Day", y = "NRS 0–10", colour = NULL,
           subtitle = "The sum of IL-31 · Th2 · barrier · superantigen. Central sensitisation (SENS) creates the slow tail.") + THEME
  })
  output$p_itchdecouple <- renderPlot({
    x <- simA()
    d <- data.frame(time = x$time,
                    `Tumour burden (baseline=1)` = (x$NSK + x$NSKR + x$NTR) / (x$NSK[1] + x$NSKR[1] + x$NTR[1]),
                    `Pruritus (baseline=1)` = x$PRUR / x$PRUR[1],
                    `IL-31 (baseline=1)` = x$IL31 / x$IL31[1],
                    check.names = FALSE) %>% pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 0.8) +
      labs(title = "Pruritus does not simply track tumour burden", x = "Day",
           y = "Relative to baseline", colour = NULL) + THEME
  })

  ## ---- ⑨ biomarkers -------------------------------------------------------
  output$p_biomk <- renderPlot({
    d <- both() %>% select(time, arm, `TARC (pg/mL)` = TARC, `sIL-2R (U/mL)` = SIL2R,
                           `LDH (U/L)` = LDH, `Th2 tone` = TH2, `IL-31` = IL31,
                           `IFN-γ` = IFNG) %>% pivot_longer(-c(time, arm))
    ggplot(d, aes(time, value, colour = arm)) + geom_line(linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y") +
      labs(title = "Serum biomarkers and cytokines", x = "Day", y = NULL, colour = NULL) + THEME
  })

  ## ---- ⑩ barrier ----------------------------------------------------------
  output$p_barrier <- renderPlot({
    d <- both() %>% select(time, arm, `Barrier integrity` = BARR,
                           `S. aureus` = SAUR, `Superantigen` = SAG) %>%
      pivot_longer(-c(time, arm))
    ggplot(d, aes(time, value, colour = arm)) + geom_line(linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y") +
      labs(title = "Barrier – staphylococcus – superantigen feedback", x = "Day", y = NULL, colour = NULL,
           subtitle = "Antibiotics break this loop without killing the clone") + THEME
  })
  output$p_infect <- renderPlot({
    d <- both() %>% select(time, arm, `Cumulative infection hazard` = HZI,
                           `Cumulative disease hazard` = HZD, `Survival probability` = SURV) %>%
      pivot_longer(-c(time, arm))
    ggplot(d, aes(time, value, colour = arm)) + geom_line(linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y") +
      labs(title = "Hazard integral and survival", x = "Day", y = NULL, colour = NULL) + THEME
  })

  ## ---- ⑪ toxicity ---------------------------------------------------------
  output$p_tox <- renderPlot({
    d <- both() %>% select(time, arm, `Neutrophils (10⁹/L)` = CIRC, `Platelets (10⁹/L)` = PLT,
                           `Peripheral neuropathy` = PN, `CD4 (/µL)` = CD4N,
                           `Free T4 (ng/dL)` = FT4, `Triglycerides (mg/dL)` = TG,
                           `Moga-associated rash` = MAR) %>% pivot_longer(-c(time, arm))
    ggplot(d, aes(time, value, colour = arm)) + geom_line(linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y") +
      labs(title = "Toxicity indicators", x = "Day", y = NULL, colour = NULL) + THEME
  })
  output$t_tox <- DT::renderDataTable({
    f <- function(x, nm) data.frame(
      arm = nm, ANCMin = round(min(x$CIRC), 2), PltMin = round(min(x$PLT), 0),
      NeuropathyMax = round(max(x$PN), 2), CD4Min = round(min(x$CD4N), 0),
      FT4Min = round(min(x$FT4), 2), TGMax = round(max(x$TG), 0),
      RashMax = round(max(x$MAR), 2))
    DT::datatable(bind_rows(f(simA(), input$reg), f(simB(), input$reg2)),
                  options = list(dom = "t"), rownames = FALSE)
  })

  ## ---- ⑫ comparison -------------------------------------------------------
  output$p_cmp <- renderPlot({
    d <- both() %>% group_by(arm) %>%
      transmute(time, mSWAT = MSWAT / first(MSWAT),
                Sezary = SEZCT / max(first(SEZCT), 1e-9),
                Surveillance = SURVSK / first(SURVSK), NRS = PRUR / first(PRUR)) %>%
      ungroup() %>% pivot_longer(-c(time, arm))
    ggplot(d, aes(time, value, colour = arm)) + geom_line(linewidth = 0.85) +
      facet_wrap(~name, scales = "free_y") + geom_hline(yintercept = 1, colour = "grey80") +
      labs(title = "Regimen A vs B (all baseline = 1)", x = "Day", y = "Relative to baseline", colour = NULL) + THEME
  })
  output$t_cmp <- DT::renderDataTable({
    DT::datatable(bind_rows(cbind(arm = input$reg, respA()),
                            cbind(arm = input$reg2, respB())),
                  options = list(dom = "t", scrollX = TRUE), rownames = FALSE)
  })

  ## ---- ⑬ structural experiments -------------------------------------------
  output$t_e1 <- DT::renderDataTable({
    DT::datatable(E1_effector_not_exposure(), options = list(dom = "t"), rownames = FALSE)
  })
  output$p_e6 <- renderPlot({
    d <- E6_residency_sweep()
    ggplot(d, aes(CD69, Sezary_per_uL)) +
      geom_line(colour = PAL[["blood"]], linewidth = 1) + geom_point(size = 2) +
      geom_hline(yintercept = c(250, 1000), linetype = "dotted") +
      scale_y_log10() + scale_x_reverse() +
      labs(title = "From MF to Sézary by CD69 alone",
           subtitle = "Same stage, same equations — only tissue residency changes",
           x = "CD69 (tissue residency, lower to the right)", y = "Sézary /µL (log)") + THEME
  })
}

shinyApp(ui, server)
