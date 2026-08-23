## =====================================================================
##  Multiple System Atrophy (MSA) — QSP Shiny dashboard
##  ---------------------------------------------------------------------
##  10 tabs, driven by the 61-ODE model in msa_mrgsolve_model.R.
##
##  DESIGN NOTE — why the compute functions are top-level.
##  Every panel is produced by a plain function (`tab_*`) that takes a
##  settings list and returns data / draws a plot.  Nothing that computes
##  lives inside renderPlot().  That makes the whole dashboard testable
##  head-less: msa_shiny_selftest() below walks all 10 tabs on a null PNG
##  device and asserts that each one produces finite numbers.  Run it with
##      Rscript -e 'source("msa_shiny_app.R"); msa_shiny_selftest()'
##  which needs mrgsolve but NOT a browser and NOT shiny.
##
##  THE ORGANISING IDEA THE DASHBOARD IS BUILT AROUND.
##  Tab 1 does not show "disease severity".  It shows the SEVEN MULTIPLICANDS
##  of upright blood pressure separately, because the whole clinical logic of
##  MSA is that different diseases break different factors of the same
##  product, and that is what decides which drug works.  Every other tab is
##  a consequence of that decomposition.
##
##  Launch:  Rscript -e 'shiny::runApp("msa_shiny_app.R", port = 8080)'
##  Educational / research use only.  NOT for clinical decision-making.
## =====================================================================

suppressPackageStartupMessages({
  library(mrgsolve)
})

## The model file defines: mod, pheno, pp, msa_history, msa_state_at,
## msa_tilt, tilt_metrics, msa_day, day_metrics, rx, rx_chronic,
## msa_history_rx.  Source it without triggering its own __main__ block.
.msa_src <- local({
  f <- "msa_mrgsolve_model.R"
  if (!file.exists(f)) f <- file.path(dirname(sys.frame(1)$ofile %||% "."), f)
  f
})
`%||%` <- function(a, b) if (is.null(a)) b else a
local({
  src <- readLines("msa_mrgsolve_model.R", warn = FALSE)
  cut <- grep('^if \\(!interactive\\(\\)', src)
  if (length(cut)) src <- src[seq_len(cut[1] - 1L)]
  eval(parse(text = paste(src, collapse = "\n")), envir = globalenv())
})

PH_LABEL <- c(MSA_P = "MSA-P (striatonigral predominant)",
              MSA_C = "MSA-C (olivopontocerebellar predominant)",
              PAF   = "PAF (pure autonomic failure, postganglionic lesion)",
              PD    = "PD (Parkinson disease, postsynaptic intact)")

PAL <- c("#2b6cb0", "#c53030", "#2f855a", "#b7791f", "#6b46c1",
         "#0987a0", "#97266d", "#4a5568")

## ---------------------------------------------------------------------
##  settings -> simulations (cached per settings signature)
## ---------------------------------------------------------------------
.cache <- new.env(parent = emptyenv())

sim_key <- function(s) paste(unlist(s), collapse = "|")

drug_events <- function(s, base_time = 1) {
  e <- NULL
  add <- function(x) if (is.null(e)) x else c(e, x)
  if (isTRUE(s$mido))  e <- add(ev(amt = s$mido_mg, cmt = "MIDG",
                                   time = base_time, ii = 5, addl = 2))
  if (isTRUE(s$mido_hs)) e <- add(ev(amt = s$mido_mg, cmt = "MIDG",
                                     time = 22, ii = 24, addl = 6))
  if (isTRUE(s$drox))  e <- add(ev(amt = s$drox_mg, cmt = "DRXG",
                                   time = base_time, ii = 5, addl = 2))
  if (isTRUE(s$carbi)) e <- add(ev(amt = 50, cmt = "CBDG",
                                   time = base_time, ii = 5, addl = 2))
  if (isTRUE(s$atx))   e <- add(ev(amt = s$atx_mg, cmt = "ATXG",
                                   time = base_time, ii = 12, addl = 3))
  if (isTRUE(s$fludro))e <- add(ev(amt = 0.1, cmt = "FLUG",
                                   time = base_time, ii = 24, addl = 6))
  if (isTRUE(s$pyrido))e <- add(ev(amt = 60, cmt = "PYRG",
                                   time = base_time, ii = 5, addl = 2))
  if (isTRUE(s$ddavp)) e <- add(ev(amt = 2, cmt = "DDA",
                                   time = 22, ii = 24, addl = 6))
  if (isTRUE(s$octreo))e <- add(ev(amt = 50, cmt = "OCT", time = base_time + 6.5))
  if (isTRUE(s$water)) e <- add(ev(amt = 500, cmt = "WB", time = base_time + 6))
  if (isTRUE(s$ldopa)) e <- add(ev(amt = s$ldopa_mg, cmt = "LDG",
                                   time = base_time, ii = 5, addl = 2))
  e
}

extra_params <- function(s) {
  p <- list(COQ2F = s$coq2, FTON = s$fton, MEALAMP = if (isTRUE(s$fasting)) 0 else 1)
  if (isTRUE(s$binder))  p$BINDER <- 1
  if (isTRUE(s$counter)) p$COUNTERM <- 1
  if (isTRUE(s$headup))  p$HUTF <- 0.22
  p
}

get_history <- function(s) {
  k <- paste0("hist|", sim_key(s))
  if (!is.null(.cache[[k]])) return(.cache[[k]])
  d <- msa_history(14, s$pheno, p = extra_params(s), delta = 72)
  .cache[[k]] <- d
  d
}

get_state <- function(s) {
  k <- paste0("state|", sim_key(s))
  if (!is.null(.cache[[k]])) return(.cache[[k]])
  d <- msa_state_at(s$year, s$pheno, p = extra_params(s))
  .cache[[k]] <- d
  d
}

get_tilt <- function(s, with_drug = TRUE) {
  k <- paste0("tilt|", with_drug, "|", sim_key(s))
  if (!is.null(.cache[[k]])) return(.cache[[k]])
  st <- get_state(s)
  ev0 <- if (with_drug) drug_events(s, base_time = 0) else NULL
  d <- msa_tilt(st, s$pheno, p = extra_params(s), dose = ev0)
  .cache[[k]] <- d
  d
}

get_day <- function(s, with_drug = TRUE) {
  k <- paste0("day|", with_drug, "|", sim_key(s))
  if (!is.null(.cache[[k]])) return(.cache[[k]])
  st <- get_state(s)
  ev0 <- if (with_drug) drug_events(s, base_time = 1) else NULL
  d <- msa_day(st, s$pheno, p = extra_params(s), dose = ev0, days = 2, delta = 0.1)
  .cache[[k]] <- d
  d
}

## ---------------------------------------------------------------------
##  small plotting helpers (base graphics only — no extra dependencies)
## ---------------------------------------------------------------------
.par <- function(...) par(mar = c(4.2, 4.4, 2.6, 1.2), mgp = c(2.5, 0.7, 0),
                          las = 1, cex.axis = 0.9, cex.lab = 0.95, ...)

lineplot <- function(x, ys, labels, xlab, ylab, main, cols = PAL,
                     lty = 1, ylim = NULL, hlines = NULL, vshade = NULL) {
  ys <- if (is.list(ys)) ys else list(ys)
  if (is.null(ylim)) {
    r <- range(unlist(ys)[is.finite(unlist(ys))])
    if (!all(is.finite(r)) || diff(r) == 0) r <- r[1] + c(-1, 1)
    ylim <- r + c(-0.06, 0.10)*diff(r)
  }
  .par()
  plot(NA, xlim = range(x), ylim = ylim, xlab = xlab, ylab = ylab, main = main)
  if (!is.null(vshade)) for (i in seq_len(nrow(vshade)))
    rect(vshade[i, 1], ylim[1], vshade[i, 2], ylim[2],
         col = "#f0f3f8", border = NA)
  grid(col = "#e2e8f0", lty = 1)
  if (!is.null(hlines)) abline(h = hlines, col = "#a0aec0", lty = 3)
  for (i in seq_along(ys))
    lines(x, ys[[i]], col = cols[(i - 1) %% length(cols) + 1], lwd = 2.1,
          lty = if (length(lty) > 1) lty[i] else lty)
  if (!is.null(labels))
    legend("topleft", legend = labels, col = cols[seq_along(ys)], lwd = 2.1,
           lty = if (length(lty) > 1) lty[seq_along(ys)] else lty,
           bty = "n", cex = 0.82)
  box()
  invisible(TRUE)
}

## Dosing records add extra output rows, so a drug arm and its drug-free
## comparator do NOT share a time grid.  Interpolate the comparator onto the
## drug arm's grid instead of assuming equal lengths.
onto <- function(xt, d, col) approx(d$time, d[[col]], xout = xt, rule = 2)$y

night_bands <- function(tmax) {
  b <- NULL
  d <- 0
  while (d*24 < tmax) {
    b <- rbind(b, c(d*24, d*24 + 7), c(d*24 + 22, min(d*24 + 24, tmax)))
    d <- d + 1
  }
  b[b[, 1] < tmax, , drop = FALSE]
}

## =====================================================================
##  TAB 1 — Patient profile: the SEVEN MULTIPLICANDS, shown separately
## =====================================================================
tab1_gains <- function(s) {
  st <- get_state(s)
  gc <- st$NMED^2.20 * st$NIML^1.60
  data.frame(
    term = c("Baroreceptor afferent", "Central integration\nNTS/RVLM",
             "Preganglionic output\nIML", "Postganglionic NE storage\nPOSTG",
             "α1 receptor density\nA1R", "Effective circulating volume\nECF",
             "Venous capacitance regulation"),
    value = c(1.0, st$NMED^2.20, st$NIML^1.60, st$POSTG,
              min(st$A1R/1.6, 1), min(st$ECF/14, 1.05), 1.0),
    stringsAsFactors = FALSE
  ) -> d
  d$product_gain <- c(1, gc, 1, st$POSTG, st$A1R/1.6, st$ECF/14, 1)
  attr(d, "G_CENT") <- gc
  attr(d, "state") <- st
  d
}

tab1_plot <- function(s) {
  d <- tab1_gains(s)
  .par(mar = c(7.5, 4.4, 3.0, 1.2))
  bp <- barplot(d$value, names.arg = d$term, las = 2, ylim = c(0, 1.15),
                col = ifelse(d$value < 0.5, "#c53030",
                      ifelse(d$value < 0.8, "#b7791f", "#2f855a")),
                border = NA, ylab = "Residual gain (1 = normal)",
                main = sprintf("%s · %.1f years post-onset — multiplicative decomposition of standing blood pressure",
                               PH_LABEL[[s$pheno]], s$year), cex.names = 0.72)
  abline(h = c(0.5, 0.8), col = "#a0aec0", lty = 3)
  text(bp, pmin(d$value + 0.06, 1.12), sprintf("%.2f", d$value), cex = 0.8)
  mtext(sprintf("Central gain G_CENT = %.3f   |   Postganglionic integrity POSTG = %.3f   →   %s",
                attr(d, "G_CENT"), attr(d, "state")$POSTG,
                if (attr(d, "state")$POSTG > 2.5*attr(d, "G_CENT"))
                  "Central-lesion predominant: conditions favour a NET inhibitor (atomoxetine)"
                else "Postganglionic-lesion predominant: only a direct α1 agonist (midodrine) can be expected to work"),
        side = 1, line = 6.1, cex = 0.78)
  invisible(d)
}

## =====================================================================
##  TAB 2 — Orthostatic challenge (head-up tilt / active stand)
## =====================================================================
tab2_plot <- function(s) {
  tt <- get_tilt(s, with_drug = TRUE)
  t0 <- get_tilt(s, with_drug = FALSE)
  x <- (tt$time - 1.0)*60          # minutes relative to tilt
  b <- function(col) onto(tt$time, t0, col)
  layout(matrix(1:4, 2, 2, byrow = TRUE))
  lineplot(x, list(b("SBP"), tt$SBP, b("DBP"), tt$DBP),
           c("SBP no drug", "SBP drug", "DBP no drug", "DBP drug"),
           "Time after standing (min)", "Blood pressure (mmHg)", "Tilt-table blood pressure response",
           cols = c("#2b6cb0", "#c53030", "#90a4c8", "#e08a8a"),
           lty = c(2, 1, 2, 1), hlines = 0)
  abline(v = 0, col = "#4a5568", lty = 2)
  lineplot(x, list(b("HR"), tt$HR), c("HR no drug", "HR drug"),
           "Time after standing (min)", "Heart rate (bpm)", "Heart-rate response (neurogenic index)",
           lty = c(2, 1))
  abline(v = 0, col = "#4a5568", lty = 2)
  lineplot(x, list(b("NEPL"), tt$NEPL), c("Plasma NE no drug", "Plasma NE drug"),
           "Time after standing (min)", "Plasma NE (pg/mL)",
           "Plasma norepinephrine — normal supine / blunted rise on standing", lty = c(2, 1))
  abline(v = 0, col = "#4a5568", lty = 2)
  lineplot(x, list(b("CBFMARG"), tt$CBFMARG), c("No drug", "Drug"),
           "Time after standing (min)", "MAP − lower limit of autoregulation (mmHg)",
           "Cerebral blood flow autoregulation margin", lty = c(2, 1), hlines = 0)
  abline(v = 0, col = "#4a5568", lty = 2)
  layout(1)
  invisible(TRUE)
}

tab2_table <- function(s) {
  m0 <- tilt_metrics(get_tilt(s, FALSE))
  m1 <- tilt_metrics(get_tilt(s, TRUE))
  ## NOTE: column NAMES are assigned as strings, never used as identifiers.
  ## Non-ASCII identifiers (e.g. a literal Korean column name) fail to parse outside a UTF-8
  ## locale, and Rscript defaults to LC_CTYPE=C in many containers.
  out <- data.frame(
    metric = c("Supine SBP (mmHg)", "3-min standing SBP (mmHg)", "ΔSBP (mmHg)",
               "ΔHR (bpm)", "ΔHR/ΔSBP (neurogenic: <0.5)",
               "nOH criteria met (ΔSBP ≤ −20)", "Supine plasma NE (pg/mL)",
               "Standing NE rise (%)"),
    no_drug = c(round(m0$SBP_supine, 1), round(m0$SBP_3min, 1),
                round(m0$dSBP, 1), round(m0$dHR, 1),
                round(m0$HR_SBP_ratio, 2), ifelse(m0$nOH, "Yes", "No"),
                round(m0$NE_supine, 0), round(m0$NE_rise_pct, 1)),
    on_drug = c(round(m1$SBP_supine, 1), round(m1$SBP_3min, 1),
                round(m1$dSBP, 1), round(m1$dHR, 1),
                round(m1$HR_SBP_ratio, 2), ifelse(m1$nOH, "Yes", "No"),
                round(m1$NE_supine, 0), round(m1$NE_rise_pct, 1)),
    stringsAsFactors = FALSE)
  names(out) <- c("Metric", "No drug", "Drug")
  out
}

## =====================================================================
##  TAB 3 — 24-h BP and THE NOCTURNAL LOOP
## =====================================================================
tab3_plot <- function(s) {
  dd <- get_day(s, TRUE)
  d0 <- get_day(s, FALSE)
  nb <- night_bands(max(dd$time))
  b <- function(col) onto(dd$time, d0, col)
  layout(matrix(1:4, 2, 2, byrow = TRUE))
  lineplot(dd$time, list(b("SBP"), dd$SBP), c("No drug", "Drug"),
           "Time (h, 0 = midnight)", "SBP (mmHg)",
           "48-hour SBP — grey = nocturnal supine",
           lty = c(2, 1), hlines = c(140, 90), vshade = nb)
  lineplot(dd$time, list(b("UNAV"), dd$UNAV), c("No drug", "Drug"),
           "Time (h)", "Urinary sodium excretion (mmol/h)",
           "Pressure natriuresis — nocturnal hypertension spends volume overnight",
           lty = c(2, 1), vshade = nb)
  lineplot(dd$time, list(b("UVOL"), dd$UVOL), c("No drug", "Drug"),
           "Time (h)", "Urine output (L/h)", "Nocturnal polyuria (desmopressin target)",
           lty = c(2, 1), vshade = nb)
  lineplot(dd$time, list(b("ECF"), dd$ECF), c("No drug", "Drug"),
           "Time (h)", "Extracellular fluid volume (L)",
           "Overnight volume loss → next-morning orthostatic hypotension",
           lty = c(2, 1), vshade = nb)
  layout(1)
  invisible(TRUE)
}

tab3_table <- function(s) {
  a <- day_metrics(get_day(s, FALSE))
  b <- day_metrics(get_day(s, TRUE))
  out <- data.frame(
    metric = c("Peak nocturnal supine SBP (mmHg)", "Morning 07:30 standing SBP (mmHg)",
               "Afternoon standing SBP (mmHg)", "Nocturnal sodium excretion (mmol)",
               "Nocturnal urine output (L)", "Early-morning extracellular fluid volume (L)", "OHSA symptom score"),
    no_drug = round(unlist(a[c(1, 2, 3, 4, 5, 6, 7)]), 2),
    on_drug = round(unlist(b[c(1, 2, 3, 4, 5, 6, 7)]), 2),
    stringsAsFactors = FALSE, row.names = NULL)
  names(out) <- c("Metric", "No drug", "Drug")
  out
}

## =====================================================================
##  TAB 4 — Drug pharmacokinetics
## =====================================================================
tab4_plot <- function(s) {
  dd <- get_day(s, TRUE)
  series <- list(); labs <- character(0)
  push <- function(v, l) { if (max(v, na.rm = TRUE) > 1e-9) {
    series[[length(series) + 1]] <<- v; labs <<- c(labs, l) } }
  push(dd$CDMM, "desglymidodrine (ng/mL)")
  push(dd$CDRX/10, "droxidopa (×10 ng/mL)")
  push(dd$CATX, "atomoxetine (ng/mL)")
  push(dd$CFLU, "fludrocortisone (ng/mL)")
  push(dd$CPYR, "pyridostigmine (ng/mL)")
  push(dd$CLDC/10, "levodopa (×10 ng/mL)")
  if (!length(series)) { series <- list(rep(0, nrow(dd))); labs <- "No drug" }
  layout(matrix(1:2, 2, 1))
  lineplot(dd$time, series, labs, "Time (h)", "Concentration",
           "Drug plasma concentration (48 hours)", vshade = night_bands(max(dd$time)))
  lineplot(dd$time, list(dd$OCC, dd$NETI, dd$CARBI),
           c("α1 occupancy OCC", "NET blockade (atomoxetine)",
             "Peripheral AADC blockade (carbidopa)"),
           "Time (h)", "Fraction (0–1)",
           "Target occupancy — carbidopa and droxidopa compete for the same enzyme",
           ylim = c(0, 1), vshade = night_bands(max(dd$time)))
  layout(1)
  invisible(TRUE)
}

## =====================================================================
##  TAB 5 — Autonomic axis (SNA, NE, renin, aldosterone, AVP, AQP2, A1R)
## =====================================================================
tab5_plot <- function(s) {
  dd <- get_day(s, TRUE)
  h  <- get_history(s)
  nb <- night_bands(max(dd$time))
  layout(matrix(1:4, 2, 2, byrow = TRUE))
  lineplot(dd$time, list(dd$SNA, dd$NES, dd$NEVES),
           c("Sympathetic activity SNA", "Synaptic NE", "Vesicular NE storage"),
           "Time (h)", "Normalised value", "Sympathetic drive and NE storage", vshade = nb)
  lineplot(dd$time, list(dd$PRA, dd$ALD, dd$MR),
           c("Plasma renin activity", "Aldosterone", "MR activation"),
           "Time (h)", "Normalised value",
           "RAAS — renin fails to rise on standing in MSA", vshade = nb)
  lineplot(dd$time, list(dd$AVP, dd$AQP2*10),
           c("Plasma AVP (pg/mL)", "AQP2 traffic (×10)"),
           "Time (h)", "Value", "AVP–AQP2 axis (baroreceptor arm lost)", vshade = nb)
  lineplot(h$YEARS, list(h$A1R, h$G_CENT + 1, h$POSTG + 1),
           c("α1 receptor density A1R", "Central gain + 1", "Postganglionic integrity + 1"),
           "Years post-onset", "Normalised value",
           "Denervation supersensitivity: receptors increase as drive dies")
  abline(v = s$year, col = "#c53030", lty = 2)
  layout(1)
  invisible(TRUE)
}

## =====================================================================
##  TAB 6 — Neuropathology progression (the glial cascade)
## =====================================================================
tab6_plot <- function(s) {
  h <- get_history(s)
  layout(matrix(1:4, 2, 2, byrow = TRUE))
  lineplot(h$YEARS, list(h$GCI, h$ASYNO, h$SEED*5, h$P25A),
           c("GCI inclusion burden", "α-syn oligomer", "Extracellular seed ×5",
             "p25α somal translocation"),
           "Years post-onset", "Arbitrary units", "Oligodendroglial α-synuclein cascade")
  abline(v = s$year, col = "#c53030", lty = 2)
  lineplot(h$YEARS, list(h$MYE, h$CQ, h$OXS/2, h$MGL/2),
           c("Myelin integrity MYE", "CoQ10 pool", "Oxidative stress ÷2",
             "Activated microglia ÷2"),
           "Years post-onset", "Normalised value", "Myelin · mitochondria · neuroinflammation")
  abline(v = s$year, col = "#c53030", lty = 2)
  lineplot(h$YEARS, list(h$NIML, h$NMED, h$POSTG),
           c("IML preganglionic neurons", "Brainstem cardiovascular neurons", "Postganglionic terminals (preserved)"),
           "Years post-onset", "Surviving fraction",
           "Autonomic axis: the centre collapses while the postganglionic arm remains", ylim = c(0, 1.05))
  abline(v = s$year, col = "#c53030", lty = 2)
  lineplot(h$YEARS, list(h$NMSN, h$NSN, h$NOPC, h$NONUF, h$NRESP),
           c("Striatal MSN (postsynaptic)", "Nigral DA (presynaptic)",
             "Olivopontocerebellar", "Onuf nucleus", "Respiratory neurons"),
           "Years post-onset", "Surviving fraction", "Regional vulnerability", ylim = c(0, 1.05))
  abline(v = s$year, col = "#c53030", lty = 2)
  layout(1)
  invisible(TRUE)
}

## =====================================================================
##  TAB 7 — Motor symptoms and the levodopa question
## =====================================================================
tab7_curve <- function(s, years = c(2, 4, 6, 8, 10)) {
  do.call(rbind, lapply(years, function(y) {
    s2 <- s; s2$year <- y
    st <- msa_state_at(y, s$pheno, p = extra_params(s))
    stp<- msa_state_at(y, "PD",    p = extra_params(s))
    g <- function(state, ph, mg) {
      e <- if (mg > 0) ev(amt = mg, cmt = "LDG", time = 1, ii = 5, addl = 2) else NULL
      d <- msa_day(state, ph, p = extra_params(s), dose = e, days = 1, delta = 0.5)
      mean(d$PARK)
    }
    off_m <- g(st, s$pheno, 0); on_m <- g(st, s$pheno, s$ldopa_mg)
    off_p <- g(stp, "PD", 0);   on_p <- g(stp, "PD", s$ldopa_mg)
    data.frame(year = y,
               G_post = st$NMSN^1.15,
               park_off = off_m, park_on = on_m,
               gain_pct = 100*(off_m - on_m)/max(off_m, 1e-6),
               pd_gain_pct = 100*(off_p - on_p)/max(off_p, 1e-6))
  }))
}

tab7_plot <- function(s) {
  h <- get_history(s)
  cv <- tab7_curve(s)
  layout(matrix(1:4, 2, 2, byrow = TRUE))
  lineplot(h$YEARS, list(h$PARK, h$ATAX, h$STRID*20),
           c("Parkinsonism (UMSARS-II score)", "Cerebellar ataxia (UMSARS-II score)",
             "Stridor index ×20"),
           "Years post-onset", "Score", "Composition of the motor phenotype")
  abline(v = s$year, col = "#c53030", lty = 2)
  lineplot(cv$year, list(cv$gain_pct, cv$pd_gain_pct),
           c(paste0(PH_LABEL[[s$pheno]], " response"), "PD comparison-group response"),
           "Years post-onset", "Reduction in parkinsonism from levodopa (%)",
           sprintf("Levodopa %d mg t.i.d. — same exposure, different target", s$ldopa_mg))
  lineplot(cv$year, list(cv$G_post), "Postsynaptic gain G_post",
           "Years post-onset", "G_post = NMSN^1.15",
           "Why the response disappears: the target disappears", ylim = c(0, 1.05))
  lineplot(h$YEARS, list(h$STRIAT, h$DA),
           c("Striatal output STRIAT", "Striatal dopamine DA"),
           "Years post-onset", "Normalised value",
           "Dopamine remains but output does not recover")
  abline(v = s$year, col = "#c53030", lty = 2)
  layout(1)
  invisible(cv)
}

## =====================================================================
##  TAB 8 — Clinical endpoints and survival
## =====================================================================
tab8_plot <- function(s) {
  h <- get_history(s)
  layout(matrix(1:4, 2, 2, byrow = TRUE))
  lineplot(h$YEARS, list(h$U1, h$U2, h$U1 + h$U2),
           c("UMSARS-I (ADL)", "UMSARS-II (motor exam)", "UMSARS I+II"),
           "Years post-onset", "Score", "UMSARS progression")
  abline(v = s$year, col = "#c53030", lty = 2)
  lineplot(h$YEARS, list(h$SCOPA, h$OHSA*4, h$PVRES/10),
           c("SCOPA-AUT", "OHSA symptoms ×4", "Residual urine ÷10 (mL)"),
           "Years post-onset", "Score", "Autonomic symptom burden")
  abline(v = s$year, col = "#c53030", lty = 2)
  lineplot(h$YEARS, list(h$U4), "UMSARS-IV overall disability",
           "Years post-onset", "1–5", "Overall disability", ylim = c(1, 5.2))
  abline(v = s$year, col = "#c53030", lty = 2)
  lineplot(h$YEARS, list(h$SURV), "Survival probability",
           "Years post-onset", "S(t)", "Model survival curve", ylim = c(0, 1.02),
           hlines = 0.5)
  abline(v = s$year, col = "#c53030", lty = 2)
  layout(1)
  invisible(TRUE)
}

tab8_table <- function(s) {
  h <- get_history(s)
  at <- function(y, col) h[[col]][which.min(abs(h$YEARS - y))]
  yrs <- c(2, 4, 6, 8, 10, 12)
  out <- data.frame(year = yrs,
             UMSARS_I = round(sapply(yrs, at, "U1"), 1),
             UMSARS_II = round(sapply(yrs, at, "U2"), 1),
             SCOPA_AUT = round(sapply(yrs, at, "SCOPA"), 1),
             OHSA = round(sapply(yrs, at, "OHSA"), 2),
             PVR_mL = round(sapply(yrs, at, "PVRES"), 0),
             survival = round(sapply(yrs, at, "SURV"), 3),
             stringsAsFactors = FALSE)
  names(out) <- c("Years", "UMSARS-I", "UMSARS-II", "SCOPA-AUT", "OHSA",
                  "Residual urine (mL)", "Survival probability")
  out
}

## =====================================================================
##  TAB 9 — Scenario comparison (the therapeutic trade-off table)
## =====================================================================
SCENARIOS <- list(
  "No treatment"                        = list(),
  "Non-pharmacological (binder+counter-manoeuvre+water)"   = list(binder = TRUE, counter = TRUE, water = TRUE),
  "Head-up sleeping"                 = list(headup = TRUE),
  "Bedtime desmopressin"             = list(ddavp = TRUE),
  "Head-up + desmopressin"       = list(headup = TRUE, ddavp = TRUE),
  "midodrine 10 mg t.i.d."        = list(mido = TRUE, mido_mg = 10),
  "midodrine at bedtime (error)"    = list(mido_hs = TRUE, mido_mg = 10),
  "droxidopa 300 mg t.i.d."       = list(drox = TRUE, drox_mg = 300),
  "droxidopa + carbidopa"         = list(drox = TRUE, drox_mg = 300, carbi = TRUE),
  "atomoxetine 18 mg b.i.d."      = list(atx = TRUE, atx_mg = 18),
  "fludrocortisone 0.1 mg"        = list(fludro = TRUE),
  "pyridostigmine 60 mg t.i.d."   = list(pyrido = TRUE),
  "midodrine + fludro + non-pharmacological"   = list(mido = TRUE, mido_mg = 10, fludro = TRUE,
                                         binder = TRUE, counter = TRUE)
)

tab9_table <- function(s, which = names(SCENARIOS)) {
  do.call(rbind, lapply(which, function(nm) {
    s2 <- modifyList(s, SCENARIOS[[nm]])
    tm <- tilt_metrics(get_tilt(s2, TRUE))
    dm <- day_metrics(get_day(s2, TRUE))
    data.frame(scenario = nm,
               stand_SBP_3min = round(tm$SBP_3min, 1),
               dSBP = round(tm$dSBP, 1),
               night_supine_peak = round(dm$SBP_supine_peak_night, 1),
               morning_stand_SBP = round(dm$SBP_morning, 1),
               night_Na_mmol = round(dm$UNa_night_mmol, 1),
               OHSA = round(dm$OHSA, 2),
               stringsAsFactors = FALSE)
  }))
}

tab9_plot <- function(s, which = names(SCENARIOS)) {
  tb <- tab9_table(s, which)
  .par(mar = c(4.2, 14, 3.0, 1.2))
  tb <- tb[order(tb$morning_stand_SBP), ]
  xr <- range(c(tb$stand_SBP_3min, tb$night_supine_peak, tb$morning_stand_SBP))
  plot(NA, xlim = xr + c(-4, 8), ylim = c(0.5, nrow(tb) + 0.5),
       yaxt = "n", xlab = "SBP (mmHg)", ylab = "",
       main = "Therapeutic trade-off: standing BP rise vs. supine hypertension")
  grid(col = "#e2e8f0", lty = 1)
  axis(2, at = seq_len(nrow(tb)), labels = tb$scenario, las = 1, cex.axis = 0.8)
  abline(v = 140, col = "#c53030", lty = 2)
  for (i in seq_len(nrow(tb))) {
    segments(tb$stand_SBP_3min[i], i, tb$night_supine_peak[i], i,
             col = "#cbd5e0", lwd = 5)
    points(tb$stand_SBP_3min[i], i, pch = 19, col = "#2b6cb0", cex = 1.2)
    points(tb$morning_stand_SBP[i], i, pch = 17, col = "#2f855a", cex = 1.1)
    points(tb$night_supine_peak[i], i, pch = 19, col = "#c53030", cex = 1.2)
  }
  legend("bottomright", c("3-min standing SBP", "Morning 07:30 standing SBP",
                          "Peak nocturnal supine SBP", "Supine hypertension threshold 140"),
         pch = c(19, 17, 19, NA), lty = c(NA, NA, NA, 2),
         col = c("#2b6cb0", "#2f855a", "#c53030", "#c53030"),
         bty = "n", cex = 0.78)
  box()
  invisible(tb)
}

## =====================================================================
##  TAB 10 — Biomarkers and the differential diagnosis panel
## =====================================================================
tab10_plot <- function(s) {
  h <- get_history(s)
  layout(matrix(1:4, 2, 2, byrow = TRUE))
  lineplot(h$YEARS, list(h$NFL), "Plasma NfL (pg/mL)",
           "Years post-onset", "pg/mL",
           "Neurofilament light chain — a marker of loss rate")
  abline(v = s$year, col = "#c53030", lty = 2)
  lineplot(h$YEARS, list(h$SEED, h$GCI/10),
           c("Extracellular seed (RT-QuIC correlate)", "GCI burden ÷10"),
           "Years post-onset", "Arbitrary units", "Model correlate of the seed amplification assay")
  abline(v = s$year, col = "#c53030", lty = 2)
  lineplot(h$YEARS, list(h$POSTG, h$G_CENT),
           c("Postganglionic integrity (¹²³I-MIBG correlate)", "Central gain G_CENT"),
           "Years post-onset", "Fraction",
           "MIBG is normal yet orthostatic hypotension is severe = MSA", ylim = c(0, 1.05))
  abline(v = s$year, col = "#c53030", lty = 2)
  ## differential-diagnosis panel across all four phenotypes at this year
  ph <- names(pheno)
  vals <- t(sapply(ph, function(p) {
    st <- msa_state_at(s$year, p, p = extra_params(s))
    c(G_CENT = st$NMED^2.20*st$NIML^1.60, POSTG = st$POSTG, NMSN = st$NMSN)
  }))
  ## vals is phenotypes x metrics; barplot() groups by COLUMN, so pass it as
  ## is to get one group per metric with one bar per phenotype.
  .par(mar = c(5.2, 4.4, 3.0, 1.2))
  cols <- c("#2b6cb0", "#c53030", "#2f855a", "#b7791f")
  barplot(vals, beside = TRUE, ylim = c(0, 1.25), border = NA,
          col = cols[seq_len(nrow(vals))],
          names.arg = c("Central gain G_CENT", "Postganglionic integrity POSTG", "Striatal MSN"),
          ylab = "Residual fraction",
          main = sprintf("Phenotype differentiation (%.0f years post-onset)", s$year))
  legend("topright", rownames(vals), fill = cols[seq_len(nrow(vals))],
         bty = "n", cex = 0.78)
  layout(1)
  invisible(vals)
}

## =====================================================================
##  HEAD-LESS SELF-TEST — walks all 10 tabs with no browser and no shiny
## =====================================================================
default_settings <- function() list(
  pheno = "MSA_P", year = 7, coq2 = 1.0, fton = 0.22, fasting = FALSE,
  mido = FALSE, mido_hs = FALSE, mido_mg = 10,
  drox = FALSE, drox_mg = 300, carbi = FALSE,
  atx = FALSE, atx_mg = 18, fludro = FALSE, pyrido = FALSE,
  ddavp = FALSE, octreo = FALSE, water = FALSE,
  binder = FALSE, counter = FALSE, headup = FALSE,
  ldopa = FALSE, ldopa_mg = 200
)

msa_shiny_selftest <- function(quick = FALSE) {
  ok <- TRUE
  png(tempfile(fileext = ".png"), width = 1100, height = 800)
  on.exit({ dev.off(); layout(1) }, add = TRUE)
  s <- default_settings()
  s$year <- 7
  s$mido <- TRUE; s$fludro <- TRUE; s$ldopa <- TRUE
  chk <- function(name, expr) {
    r <- tryCatch({ v <- expr; TRUE }, error = function(e) { message("   ", conditionMessage(e)); FALSE })
    cat(sprintf("  [%s] %s\n", if (r) "PASS" else "FAIL", name))
    if (!r) ok <<- FALSE
    invisible(r)
  }
  cat("\n--- Shiny panel self-test (head-less) ---\n")
  chk("tab 1  patient profile / multiplicative-term decomposition", tab1_plot(s))
  chk("tab 2  tilt-table test (plot)",        tab2_plot(s))
  chk("tab 2  tilt-table test (table)",       stopifnot(nrow(tab2_table(s)) == 8))
  chk("tab 3  24-hour BP · nocturnal loop",     tab3_plot(s))
  chk("tab 3  nocturnal loop metric table",            stopifnot(nrow(tab3_table(s)) == 7))
  chk("tab 4  drug PK",                     tab4_plot(s))
  chk("tab 5  autonomic axis",                 tab5_plot(s))
  chk("tab 6  neuropathological progression",               tab6_plot(s))
  chk("tab 7  motor symptoms and levodopa",         tab7_plot(s))
  chk("tab 8  clinical endpoints (plot)",      tab8_plot(s))
  chk("tab 8  clinical endpoints (table)",     stopifnot(nrow(tab8_table(s)) == 6))
  sub <- if (quick) names(SCENARIOS)[c(1, 6, 7, 8, 9, 10)] else names(SCENARIOS)
  chk("tab 9  scenario comparison",               tab9_plot(s, sub))
  chk("tab 10 biomarkers · differential diagnosis",           tab10_plot(s))
  cat(sprintf("--- Shiny panel self-test %s ---\n\n", if (ok) "PASSED" else "FAILED"))
  invisible(ok)
}

## =====================================================================
##  UI
## =====================================================================
if (requireNamespace("shiny", quietly = TRUE)) {

library(shiny)

ui <- fluidPage(
  titlePanel("Multiple System Atrophy — QSP Dashboard"),
  tags$p(style = "color:#4a5568;font-size:13px;margin-top:-8px",
    "A 61-ODE model that describes standing blood pressure as the product of 'baroreceptor afferent × central integration × preganglionic output × postganglionic NE storage × α1 receptor density × effective circulating volume × venous capacitance'. It is for education and research and must not be used for clinical decision-making."),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectInput("pheno", "Phenotype (only changes vulnerability)",
                  choices = setNames(names(PH_LABEL), PH_LABEL),
                  selected = "MSA_P"),
      sliderInput("year", "Years post-onset", min = 1, max = 12, value = 7, step = 1),
      sliderInput("coq2", "CoQ10 biosynthetic capacity (<1 with COQ2 variant)",
                  min = 0.4, max = 1.0, value = 1.0, step = 0.05),
      sliderInput("fton", "Residual tonic sympathetic drive FTON (high = supine-hypertension phenotype)",
                  min = 0.10, max = 0.45, value = 0.22, step = 0.01),
      checkboxInput("fasting", "Fasting protocol (removes postprandial hypotension)", FALSE),
      tags$hr(),
      tags$b("Pressor drugs"),
      checkboxInput("mido", "midodrine t.i.d.", FALSE),
      sliderInput("mido_mg", "midodrine dose (mg)", 2.5, 10, 10, 2.5),
      checkboxInput("mido_hs", "midodrine dosed at bedtime (wrong timing)", FALSE),
      checkboxInput("drox", "droxidopa t.i.d.", FALSE),
      sliderInput("drox_mg", "droxidopa dose (mg)", 100, 600, 300, 100),
      checkboxInput("carbi", "Concomitant carbidopa (AADC inhibition)", FALSE),
      checkboxInput("atx", "atomoxetine b.i.d.", FALSE),
      sliderInput("atx_mg", "atomoxetine dose (mg)", 10, 18, 18, 2),
      checkboxInput("fludro", "fludrocortisone 0.1 mg", FALSE),
      checkboxInput("pyrido", "pyridostigmine 60 mg t.i.d.", FALSE),
      checkboxInput("ddavp", "Bedtime desmopressin", FALSE),
      checkboxInput("octreo", "octreotide (before meals)", FALSE),
      tags$hr(),
      tags$b("Non-pharmacological · Dopaminergic"),
      checkboxInput("binder", "Abdominal binder / compression stockings", FALSE),
      checkboxInput("counter", "Physical counter-manoeuvres", FALSE),
      checkboxInput("headup", "Head-up sleeping (10–30°)", FALSE),
      checkboxInput("water", "500 mL water bolus", FALSE),
      checkboxInput("ldopa", "levodopa t.i.d.", FALSE),
      sliderInput("ldopa_mg", "levodopa dose (mg)", 100, 400, 200, 50)
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        type = "tabs",
        tabPanel("1 Patient Profile", plotOutput("p1", height = "560px"),
                 tags$p(style = "font-size:12px;color:#4a5568",
                   "Each bar is an independent multiplicative term that produces standing blood pressure. Which term has collapsed determines which drug will work.")),
        tabPanel("2 Tilt-Table Test", plotOutput("p2", height = "620px"),
                 tableOutput("t2")),
        tabPanel("3 24-Hour BP · Nocturnal Loop", plotOutput("p3", height = "620px"),
                 tableOutput("t3"),
                 tags$p(style = "font-size:12px;color:#4a5568",
                   "Supine hypertension → nocturnal pressure natriuresis → overnight volume loss → worsening morning orthostatic hypotension. This cycle sets the ceiling on pressor therapy.")),
        tabPanel("4 Drug PK", plotOutput("p4", height = "620px")),
        tabPanel("5 Autonomic Axis", plotOutput("p5", height = "620px")),
        tabPanel("6 Neuropathological Progression", plotOutput("p6", height = "620px")),
        tabPanel("7 Motor Symptoms · Levodopa", plotOutput("p7", height = "620px"),
                 tableOutput("t7")),
        tabPanel("8 Clinical Endpoints", plotOutput("p8", height = "620px"),
                 tableOutput("t8")),
        tabPanel("9 Scenario Comparison", plotOutput("p9", height = "560px"),
                 tableOutput("t9")),
        tabPanel("10 Biomarkers · Differential Diagnosis", plotOutput("p10", height = "620px"))
      )
    )
  )
)

server <- function(input, output, session) {
  S <- reactive({
    s <- default_settings()
    for (nm in names(s)) if (!is.null(input[[nm]])) s[[nm]] <- input[[nm]]
    s
  })
  output$p1  <- renderPlot(tab1_plot(S()))
  output$p2  <- renderPlot(tab2_plot(S()))
  output$t2  <- renderTable(tab2_table(S()))
  output$p3  <- renderPlot(tab3_plot(S()))
  output$t3  <- renderTable(tab3_table(S()))
  output$p4  <- renderPlot(tab4_plot(S()))
  output$p5  <- renderPlot(tab5_plot(S()))
  output$p6  <- renderPlot(tab6_plot(S()))
  output$p7  <- renderPlot(tab7_plot(S()))
  output$t7  <- renderTable(tab7_curve(S()))
  output$p8  <- renderPlot(tab8_plot(S()))
  output$t8  <- renderTable(tab8_table(S()))
  output$p9  <- renderPlot(tab9_plot(S()))
  output$t9  <- renderTable(tab9_table(S()))
  output$p10 <- renderPlot(tab10_plot(S()))
}

if (identical(Sys.getenv("MSA_SHINY_RUN"), "1")) shinyApp(ui, server)

} else {
  message("shiny is not installed — the compute functions and ",
          "msa_shiny_selftest() still work without it.")
}
