# X-linked Adrenoleukodystrophy (X-ALD) — QSP Model
### X-linked Adrenoleukodystrophy · Quantitative Systems Pharmacology

| Deliverable | File |
|--------|------|
| 🗺️ Mechanistic map | [`xald_qsp_model_en.dot`](xald_qsp_model_en.dot) · [SVG](xald_qsp_model_en.svg) · [PNG](xald_qsp_model_en.png) — 169 nodes / 265 edges / 18 clusters |
| ⚙️ mrgsolve model | [`xald_mrgsolve_model_en.R`](xald_mrgsolve_model_en.R) — 63 ODEs, time unit = day, 24 scenarios |
| 📊 Shiny dashboard | [`xald_shiny_app_en.R`](xald_shiny_app_en.R) — 11 tabs |
| 📚 References | [`xald_references_en.md`](xald_references_en.md) — 52 PubMed citations (every PMID verified) |

---

## The Organising Thesis

**One lesion creates two independent variables.**

```
                     ABCD1 loss (one lesion)
                              │
              ┌───────────────┴───────────────┐
              ▼                               ▼
      Variable (1) VLCFA accumulation    Variable (2) Cerebral inflammation switch
      ─────────────────────          ─────────────────────
      Deterministic · dose-dependent      Bistable (a switch, not a rate)
      Measured (plasma C26:0)             Not measured (until MRI)
      Reachable peripherally              Reachable only within the CNS
```

The connection from (1) to (2) is written as **a threshold, not a rate
constant**, and whether that threshold is crossed is decided by the
susceptibility coefficient `SUSC`. Cerebral VLCFA enters through a
permeation gate (`VLCGATE`), and **this gate is already saturated in every
X-ALD patient** — verification gives patients 0.994 / non-carriers 0.081. A
variable behind a saturated gate cannot distinguish phenotype. So all of the
following are **outputs of the model, not coded rules**:

- **Why plasma C26:0 fails to predict phenotype.** This is the field's
  biggest negative result. In population simulations, plasma C26:0 in
  cerebral-form patients and pure-AMN patients is **the same to the
  decimal, 1.299 μmol/L**, and the only thing that diverges is `SUSC` alone.
- **Why Lorenzo's oil normalises the marker but fails to stop the
  disease.** Erucic acid does not cross the BBB (`FLOBBB` = 0), and about
  88% of CNS VLCFA is synthesised locally.
- **Why HSCT / eli-cel halts the disease without lowering the marker.**
  Corrected microglia reverse variable (2), but make up only 6% of the
  white-matter lipid mass, so they cannot move `CBR`.
- **Why a 12–18-month lag window exists after transplant.** This is a
  consequence of `KMGREP` (the microglial replacement rate), not something
  written in as a rule.

> **The drug that fixes the marker does not fix the disease, and the drug
> that fixes the disease does not fix the marker.**

Validation anchors **A17 = 0.0%** (Loes reduction from Lorenzo's oil) and
**A20 = 0.0%** (plasma C26:0 change after gene therapy) are the quantitative
form of that sentence. Neither is coded in the equations; both are computed
from comparison against a matched control.

---

## The Bistable Switch Really Is Bistable

`xald_bifurcation()` finds the saddle-node bifurcation point by bisection.
This value is not written into the equations; it is searched for.

```
SUSC* = 0.4086
   SUSC 0.40 → switch does not fire → pure AMN (peak Loes 0.3)
   SUSC 0.44 → fires → cerebral form (peak Loes 34)
```

**Age of onset is not a parameter.** It is determined by how far above the
critical value the system sits (critical slowing down). Same equations,
only `SUSC` varied:

| SUSC | Cerebral-form onset | Loes 1→15 | Phenotype |
|------|-----------|-----------|--------|
| 0.28 | — | — | Pure AMN |
| 0.40 | — | — | Pure AMN (just below critical) |
| 0.44 | 7.7 years | 2.66 years | Cerebral form (near critical → late and slow) |
| 0.46 | 6.7 years | 1.70 years | Cerebral form |
| 0.52 | 5.1 years | 0.85 years | Cerebral form |
| 1.55 | 4.1 years | 1.45 years | Cerebral form (far from critical → fast) |

The clinical observation (onset 3–10 years, mode 7 years) emerges from the
distribution of a single parameter.

### Population Simulation — the 35–40% Conversion Rate Is Also an Output

Result of running 300 individuals with between-subject variability
(`$OMEGA` = 0.36, log-normal) applied only to `SUSC`:

| Item | Model | Literature |
|------|------|------|
| Cerebral-form conversion rate | **37.0%** | 35–40% |
| Adrenal insufficiency rate | 100% | about 80% (**overpredicted** — see failure items below) |
| Cerebral-form onset age | median 4.0 years, range 2.2–8.9 years | 4–8 years, mode 7 years |

---

## Validation — Actual Run Values, Including Misses

`xald_validate()` runs A1–A23. **None of the values in the "Model" column
below are coded into the equations.** All are read out of mrgsolve 2.0.1.

| # | Anchor | Target | Model | Error |
|---|------|------|------|------|
| A1 | Normal plasma C26:0 (μmol/L) | 0.35 | **0.348** | −0.6% |
| A2 | Untreated X-ALD plasma C26:0 | 1.30 | **1.299** | −0.1% |
| A3 | C26:0 fold-elevation | 3.70 | **3.71** | +0.3% |
| A4 | Normal C26:0-lysoPC | 0.060 | **0.0594** | −1.0% |
| A5 | X-ALD C26:0-lysoPC | 0.550 | **0.557** | +1.2% |
| A6 | Fold-elevation of lysoPC at birth | 8.0 | **9.37** | +17.1% |
| A7 | Normal cortisol (nmol/L) | 300 | **295** | −1.7% |
| A8 | Age at onset of adrenal insufficiency (years) | 7.5 | **6.65** | −11.3% |
| A9 | ACTH at adrenal insufficiency (pg/mL) | 150 | **154** | +2.7% |
| A10 | Serum K⁺ at adrenal insufficiency | 5.8 | **6.15** | +6.0% |
| A11 | Cerebral-form onset age (years) | 7.0 | **6.71** | −4.1% |
| A12 | Loes 1 → 15 duration (years) | 1.50 | **1.70** | +13.3% |
| A13 | Untreated final Loes | 32 | **34** | +6.2% |
| A14 | Pure AMN peak Loes | 0.5 | **0.32** | −36% |
| A15 | AMN EDSS (age 45) | 4.0 | **3.99** | −0.2% |
| A16 | Lorenzo's oil C26:0 reduction | 50% | **45%** | −10% |
| **A17** | **Lorenzo's oil Loes reduction** | **0%** | **0.0%** | — |
| A18 | Loes rise 18 months after eli-cel | 2.5 points | **2.1 points** | −16% |
| A19 | eli-cel 24-month MFD-free | Achieved | **Achieved** | 0% |
| **A20** | **Plasma C26:0 change after eli-cel** | **0%** | **0.0%** | — |
| A21 | Neutrophil engraftment (days) | 13 | **10** | −23.1% |
| A22 | Monocyte chimerism (1 year) | 0.95 | **0.982** | +3.4% |
| A23 | Busulfan cumulative AUC (mg·h/L) | 90 | **90** | 0% |

**Within 20%: 19 / 21** (A17 and A20 have no percentage defined, since the
target is 0).

### What Was Left Unvalidated (Held-Out Failures — Not Hidden)

1. **Adrenal insufficiency prevalence 100% (literature about 80%).** No
   between-subject variability was given to the adrenal axis, so every
   simulated hemizygote has identical adrenal dynamics. This is the direct
   cost of putting `$OMEGA` on `SUSC` alone; adding a second `$OMEGA` to the
   adrenal axis would fix it, but that would put anchor A8 back up for
   recalibration, so in this release it is **reported, not fixed**.
2. **Median cerebral-form onset age 4.0 years (literature mode 7 years).**
   Individuals in the right tail of the log-normal `SUSC` distribution sit
   far from the critical value and so present early. The range (2.2–8.9
   years) overlaps the clinical range, but the centre is shifted left.
3. **A14: pure-AMN peak Loes 0.32 (target 0.5).** This misses in the
   direction of being *less* severe than the target.
4. **A21: neutrophil engraftment 10 days (observed median 13 days).**
   Haematopoietic reconstitution after conditioning is somewhat fast.
5. **The sobetirome class fails to prevent disease — this may be a
   prediction rather than a failure.** A CNS-penetrant thyromimetic lowers
   cerebral VLCFA from 3.53 → 2.47 (−30%), but because the permeation gate
   is saturated the phenotype does not change (peak Loes stays at 34). The
   model's prediction is **"metabolic correction must be nearly complete to
   close the gate"**, which is a testable claim. However, the gate's Hill
   exponent (`HVBR` = 6) is a STRUCT parameter, so the strength of this
   prediction depends on that choice. That dependence is stated explicitly
   here.
6. **Leriglitazone produces only a partial delay.** This agrees in
   direction with the ADVANCE trial's missed primary endpoint, but the
   trial's count of cerebral lesion events was not used as a numeric anchor
   (see references file §14).

### ACTH Stimulation Test (Scenario 24) — the Only Segment Solved with Real Hormone Dynamics

`HSPEED` is set back to 1 and `CIRCAMP` to 0.55, and a cosyntropin
stimulation is given at age 10. The state at the test point is obtained by
**actually ageing 10 years forward from the birth state** (skipping that
and simply setting `start = t0` tests an undamaged adrenal, giving a
reserve of 0.854 instead of 0.289 — this actually happened during
integration).

| | Baseline cortisol | Peak after stimulation | Reserve | Verdict |
|---|---|---|---|---|
| Non-carrier | 295 nmol/L | **536** | 0.965 | Normal (>500) |
| X-ALD (age 10) | 147 nmol/L | **168** | 0.289 | Adrenal insufficiency |

**A stated limitation.** Because `KMA` = 25 pg/mL, the adrenal is already
driven to about 50% of Vmax at baseline ACTH levels. Its stimulation
headroom is narrower than the real adrenal's, so a stimulation longer than
the clinical 60 minutes (3.6 hours) is needed to clear the 500 nmol/L
cutoff. Re-fitting `KMA` and `VST` to widen the headroom would break
anchors A8–A10, so **the headroom is left as is and the limitation is
recorded here.** The discriminating power itself, 536 vs 168, is large
enough.

---

## Nine Defects Found and Fixed During Integration — Including Why Each Was Wrong

This model did not work on the first run. Everything below is **a failure
that actually occurred**, and the same explanation is left in the
corresponding code comment.

1. **LSODA died in a healthy individual** (`t = day 2621`, corrector
   convergence failed, `fabs(h_) = hmin`). Cause: a `max(x, 0)`-style
   rectifier was placed on **a term that is exactly 0 at the healthy
   equilibrium**, so the derivative discontinuity (kink) sat right on the
   stable point. `HSCN`, `MYEL`, `MYESCV`, `LOESV`, `NFSV`, and potassium's
   `pos(1-ALDACT)` all had this. → rewrote the capacity-limiting/repair
   terms as **smooth logistics**.
2. **Plasma C26:0 in a healthy person converged to 0.232 instead of
   0.348.** `PXC = ALDPF + 0.45*ABCD2F` made the healthy peroxisomal
   capacity 1.45, and at the same time **gave an untreated null patient a
   free 32% of normal capacity** — ABCD2 was rescuing the disease at
   baseline, the exact opposite of why this disease exists. → ABCD2 now
   contributes **only its induced increment**, via `rect(ABCD2F - 1)`.
3. **Serum potassium in a healthy person fell to 0.46 mmol/L.** The
   aldosterone steady state overshot its target (0.30) and converged to
   0.84, and `KKSENS*(1-ALDACT)` blew up in one direction. → recalculated
   `VALD` + a two-sided bounded Hill.
4. **Cytokines in a healthy person blew up to 1e18.** The TNF → astrocyte →
   CCL2 → monocyte → TNF loop was written as an **unsaturated linear
   amplifier**, with linear gain `1.355 × BBBP`. Once BBB opening exceeds
   0.738, **no bounded steady state exists.** → made every recruitment step
   capacity-limited (logistic) plus saturating with respect to stimulus.
5. **A healthy 40-year-old scored Loes 34.** The baseline ROS value fed
   straight into the damage term, even though the healthy equilibrium was
   `MYEL` 0.97, `OLGP` 0.98. → redefined damage drive as **the excess above
   the healthy reference** (`ROSREF`, `ATPREF`), so health scores exactly
   0.
6. **The switch could not fire itself.** The ignition gate was placed on
   myelin debris, but debris production (8.8e-4/day) was three orders of
   magnitude below clearance capacity (1.76/day), so debris never reached
   the threshold. → moved the gate to **microglial deficit** (Eichler 2008:
   microglial apoptosis precedes demyelination). Debris was moved to an
   **amplification, not ignition,** term.
7. **Every patient became the cerebral form.** Microglial death was 6×
   self-renewal, so the resident pool went to 0 in every patient, taking
   clearance capacity with it. The switch became independent of
   susceptibility. → recalibrated so the pool settles on **a reduced
   plateau**. For the switch to remain a switch, some of the pool must
   survive.
8. **eli-cel came out identical to its own control.** Anti-inflammatory
   suppression alone only moved the high fixed point from MGP 0.79 → 0.61,
   and amplification still saturated. → capped it next with the
   **corrected-microglia fraction** `MGN/(MGN+MGC)`, but that made
   **untreated cerebral disease switch itself off** (the ratio collapses
   once the resident pool dies). → settled on an uncorrected myeloid
   fraction `UNCF`, with **infiltrating monocytes** included in the
   denominator. Untreated, `UNCF ≡ 1` regardless of how depleted the
   resident pool is; after gene therapy, monocyte chimerism (over weeks) is
   what first brings it down.
9. **A transplant patient whose inflammation had shut off still scored Loes
   30.** Loes was written as a "ratchet that slowly chases the current
   damage target", which **permanently preserved the highest target ever
   reached**. Even a lesion that had arrested months earlier kept climbing
   to the score implied by that transient peak. The next version
   (integrating a destruction flux) instead kept accumulating flux for
   years, because remyelinated myelin was repeatedly destroyed again,
   stretching **Loes 1→15 to 5–7 years**. → settled on a unidirectional
   tracker that chases current damage **quickly** (reflected within weeks,
   stops immediately once worsening stops, never decreases). The
   reversible contrast-enhancement score is kept out of the state variable
   and split off as a `$TABLE` readout.

It also emerged during integration that lesion spread must be
**autocatalytic**. A first-order loss term alone cannot satisfy both "Loes
2→15 in about 18 months" and "near-complete demyelination within 3–4 years"
simultaneously (a first-order term decelerates, while the clinical course
accelerates). A term (`KSPREAD`) was added in which the length of the
advancing lesion front grows with the lesion itself, making the late-stage
rate about 5× the initial rate.

---

## The Input/Output Distinction

**Inputs** (what the user sets): residual ALDP function `MUTRES`, the
deficient-cell fraction `FMOS` and its complement `CRSC` in female carriers,
susceptibility `SUSCTV`, diet `DIETSC`, vector copy number `VCNIN`, busulfan
exposure multiple, each drug's dose and start time, transplant timing
(specified by Loes threshold).

**Outputs** (what the model must produce, never told to it): plasma C26:0
and C26:0-lysoPC, tissue-specific VLCFA, permeation-gate saturation, switch
state and loop gain, the bifurcation point `SUSC*`, cerebral-form onset age,
Loes/NFS/MFD, EDSS and 6-minute walk, cortisol/ACTH/K⁺ and age at adrenal
insufficiency onset, monocyte chimerism and the corrected-microglia
trajectory, the length of the treatment lag window, cumulative MDS/AML
probability, and the population's phenotype distribution.

**Even the transplant date is not an input.** `.find_loes_day()` **first
runs the matched untreated natural history**, finds the day it reaches the
specified Loes value, and transplants on that day.

---

## Matched Controls

So that a treatment effect never becomes "a comparison of different
patients", every major treatment arm is run alongside an untreated control
that starts from **the same birth state, the same genotype, and the same
susceptibility**.

| Treatment arm | Matched control |
|--------|---------------|
| 07 Lorenzo's oil (cerebral-form susceptibility) | 08 same patient, untreated |
| 10 eli-cel (transplant at Loes 2) | 11 same patient, untreated |
| 16 hydrocortisone replacement | 17 same patient, no replacement |
| 12 eli-cel (transplant at Loes 15) | same patient as 10, only transplant timing differs |
| 14 / 15 busulfan ×0.72 / ×1.28 | same patient as 10, only exposure differs |

---

## Birth State Is Also an Output

Newborn screening is possible **because X-ALD infants are already born with
elevated C26:0-lysoPC**. This model does not input that value.
`xald_birth()` starts from the normal initial values, **simulates the
280-day in-utero period with the patient's genotype**, and uses the end
state as the birth state. Anchor A6 (9.37-fold elevation of lysoPC at
birth) is the value that emerges from that.

---

## Two Honest Confessions About Time Scale

1. **`HSPEED` (default 0.02).** Solving the endogenous hormone subsystem
   (ACTH half-life 10 minutes, cortisol 80 minutes) as-is over a 14,600-day
   horizon would require LSODA to integrate 40 years at minute-scale steps.
   Production and loss are scaled by **the same factor**, so every steady
   state is mathematically unchanged and only the approach speed slows
   down. Only the ACTH stimulation test scenario (24) restores `HSPEED = 1`
   and circadian drive (`CIRCAMP = 0.55`) to solve real dynamics.
2. **Chronic dosing is expressed as a constant-rate infusion, not a daily
   bolus.** Mean exposure and steady state are correct; within-day
   peak–trough variation is discarded. This is a choice made to avoid
   restarting the solver at each of 3,650 boluses across 40 years; acute
   scenarios where peak–trough matters use real boluses.

---

## Reproduce

```r
setwd("x-linked-adrenoleukodystrophy")
source("xald_mrgsolve_model_en.R")

xald_report()                      # validation anchors + bifurcation + population phenotype distribution
xald_validate()                    # A1-A23 table only
xald_bifurcation()                 # saddle-node bifurcation search (bisection)
d <- xald_run("s10_elicel_early")  # one scenario
xald_population(n = 300)           # population simulation

# dashboard
shiny::runApp("xald_shiny_app_en.R")
```

Render the map:

```bash
dot -Tsvg xald_qsp_model_en.dot -o xald_qsp_model_en.svg
dot -Tpng -Gdpi=150 xald_qsp_model_en.dot -o xald_qsp_model_en.png
```

Validation environment: R 4.3.3 · mrgsolve 2.0.1 · Graphviz 2.43.0.

---

## ⚠️ Disclaimer

This is a semi-quantitative QSP model for education and research purposes.
It has not been independently validated or certified and must not be used
for clinical decision-making, prescribing, or regulatory submission. The
classification of parameter sources (LIT / FIT / STRUCT) and citation
evidence are in [`xald_references_en.md`](xald_references_en.md).
