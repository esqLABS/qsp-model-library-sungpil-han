# Alagille Syndrome (ALGS) — QSP Model

**JAG1/NOTCH2 haploinsufficiency → biliary duct paucity → chronic cholestasis → pruritus · xanthoma · growth failure · biliary cirrhosis**

| Deliverable | File | Scale |
|---|---|---|
| Mechanistic map | [`algs_qsp_model.dot`](algs_qsp_model.dot) · [SVG](algs_qsp_model.svg) · [PNG](algs_qsp_model.png) | 184 nodes · 251 edges · 22 clusters |
| mrgsolve ODE model | [`algs_mrgsolve_model.R`](algs_mrgsolve_model.R) | 44 compartments · 16 scenarios |
| Shiny dashboard | [`algs_shiny_app.R`](algs_shiny_app.R) | 12 tabs · 23 outputs (all pass `testServer`) |
| References | [`algs_references.md`](algs_references.md) | 154 PMIDs (all verified against the live PubMed API) |

---

## The one structural fact everything follows from

ALGS cholestasis is not one quantity. It is the collision of two fluxes that
are set **independently** and that every clinical measurement mixes together:

- **`J_DUCT`** — how much bile the liver can physically push into the gut.
  Set developmentally by Notch dose, slowly repaired postnatally by ductular
  reaction, and moved by **no drug in this model**.
- **`R`** — how much of that bile comes straight back through ASBT. The
  **only** thing an IBAT inhibitor touches.

Write those two separately and seven consequences fall out of published
numbers rather than being asserted. Every figure below is printed by
`run_all_algs()`.

---

## Axis 1 — the drug's ceiling is a duct property, not a bile-acid property

Define **Φ_EHC = R / (S + R)**: the share of the hepatocyte's bile-acid input
that arrived back from the gut rather than being freshly synthesised. An ideal
IBAT inhibitor removes Φ_EHC of that input and not one percent more.

In an unobstructed liver Φ_EHC = **0.98**, because each molecule is recycled
~50 times per molecule synthesised. Duct paucity caps biliary output, which
caps `R` with it. Calibrated **only** to the ASSERT bile-acid ratio, the model
puts the average trial patient at

> **Φ_EHC = 0.42** (duct capacity `JCAP` = J_DUCT/S = 1.99)

Duct paucity has destroyed more than half the drug's target before the first
dose. Sweeping duct capacity moves the achievable 24-week response
monotonically from **−92%** (JCAP 6.1) to **−9%** (JCAP 0.38). The drug does
not fail because bile acids are too high. It fails because there are no ducts.

## Axis 2 — a hypothesis the model rejected, and the collision that replaced it

The design intent was that below some duct capacity an IBAT inhibitor would
**invert**: blocking ASBT removes ileal FGF19, derepresses CYP7A1, and the
newly synthesised bile acid would have no duct to leave by, so serum bile acid
would *rise* on treatment.

**The model says that does not happen** — for a reason worth more than the
hypothesis was. In ALGS the FGF19 signal is *already floored before treatment*
(baseline FGF19 = 0.04 of normal; synthesis already derepressed 2.14-fold), so
the synthesis reserve that would drive an inversion has already been spent.
Blockade in severe paucity is **futile, not harmful**.

What survives is a graded threshold, and it lands somewhere unexpected:

| model boundary | duct capacity | total bilirubin | GALA's independent cut-point |
|---|---|---|---|
| response falls below 30% | JCAP 1.25 | **5.42 mg/dL** | **5.0 mg/dL** |
| response falls below 15% (futile) | JCAP 0.60 | **10.70 mg/dL** | **10.0 mg/dL** |

GALA is a natural-history cohort keyed on native liver survival with no drug
in it. This model never saw those thresholds. The **ratio** of the two
boundaries is **1.97** here against GALA's **2.00** — and unlike the
boundaries themselves, the ratio is completely independent of the single
assumed calibration anchor (see `FAILURE 2`). Two datasets with nothing in
common appear to be measuring the same duct-capacity threshold from opposite
directions.

## Axis 3 — itch is not a bile-acid measurement, and ASSERT's placebo arm proves it arithmetically

In ASSERT the **placebo** arm's itch *fell* by 0.8 of the 1.7 points the drug
arm fell — 47% of the on-drug improvement is not drug — while the same placebo
arm's bile acids **rose** by 22 µmol/L. One axis has a large favourable
placebo response; the other has a negative one. They cannot be the same
variable. The consequence is quantitative:

| slope (itch points per µmol/L) | value | ×controlled |
|---|---|---|
| placebo-controlled (ASSERT) | 0.9/113 = 0.0080 | 1.00 |
| single-arm (ASSERT) | 1.7/90 = 0.0189 | 2.37 |
| single-arm (ICONIC) | 1.6/96 = 0.0167 | 2.09 |

The two **single-arm** slopes agree with each other to 13% and are both ~2.1×
the placebo-controlled truth. Any model calibrated on single-arm cholestasis
data over-attributes itch to bile acids about two-fold. This model carries the
components separately and reproduces the controlled slope (0.0083 vs 0.0080)
and the drug-attributable itch difference (**−0.90 vs −0.90**) exactly.

A psychophysical note that changed the pharmacology: fitting ASSERT **rejected
a logarithmic (Weber–Fechner) itch law.** Forced through a saturating 0–4 map
with a baseline at 2.8, the trial demands a Stevens exponent of **1.66** —
supralinear, which no itch psychophysics supports. Treating the 0–4 scale as a
ceiling on the *instrument* rather than on the *sensation* puts the operating
point in the near-linear range and the solved exponent falls to **0.95**. The
reported scale is not the latent variable, and pretending otherwise distorts
the drug effect.

## Axis 4 — the survival benefit is larger than its own bile-acid effect can explain

GALA supplies its own exposure–hazard gradient: relative to TB < 5.0 mg/dL the
transplant hazard is 4.8× at 5–10 and 15.6× at ≥10, i.e. a power law of
exponent **n = 1.5–1.85** (n = 1.85 reproduces the severe stratum best, 17.3
predicted vs 15.6 observed).

Running six years of maralixibat against this model's own natural history:

| exponent n | model HR | published HR | GALA TB≥10 back-predicted | observed |
|---|---|---|---|---|
| 1.60 | 0.445 | 0.305 | 11.8 | 15.6 |
| **1.85** | **0.414** | **0.305** | **17.3** | **15.6** |
| 2.40 | 0.356 | 0.305 | 40.3 | 15.6 |
| 2.90 | 0.314 | 0.305 | **87.1** | 15.6 |

**No single exponent fits both.** At the exponent that fits GALA's own strata,
the model predicts HR 0.414 where 0.305 was published; at the exponent that
reproduces 0.305, GALA's severe stratum is over-predicted **5.6-fold**. A
residual hazard ratio of 0.305/0.414 = **0.74 is not explained by bile acids.**

It gets worse, not better: gamma frailty (Axis 7) *attenuates* population
hazard ratios toward 1 over time, so the individual-level effect implied by an
observed population HR of 0.305 is larger still. The candidates — six-year
drug-persistence selection against a comparator matched on baseline labs only,
a bile-acid-independent benefit, or curvature the cross-sectional strata cannot
see — are **not separable by any published data**. The model ships with n = 1.6
so that it *under*-predicts the trial, and scenario 12 prints both errors side
by side.

## Axis 5 — where in the gut a drug acts decides whether it starves the child

ASBT sits in the **terminal ileum**, downstream of the duodenal micellar window
where fat and vitamins A/D/E/K are absorbed. At steady state a duct-limited
liver delivers `J_DUCT` to the duodenum whatever ASBT is doing. Two years, one
patient, severe paucity:

| | serum bile acid | fat absorption | vitamin D | INR | height z |
|---|---|---|---|---|---|
| untreated | 236 | 0.662 | 19.1 | 1.47 | −1.51 |
| **odevixibat** | −46% | **0.652** | **18.8** | 1.48 | **−1.42** (improves) |
| **cholestyramine** | −37% | **0.294** | **8.5** | **2.12** | **−2.38** |

Both drugs remove bile acids; only the proximal-acting one causes
steatorrhoea, and the contrast is derived from anatomy rather than asserted.
The model also finds **its own boundary case**: in *mild* paucity the duct is
not saturated, duodenal delivery does fall with the drug (6.2 → 4.3 mM), fat
absorption drops 0.844 → 0.752, and the contrast narrows. The safety argument
for IBAT inhibitors depends on the liver being duct-limited.

## Axis 6 — the population ceiling

Cardiac disease dominates ALGS mortality in year one and
vasculopathy/intracranial haemorrhage later; both are Notch-dose diseases no
cholestasis drug touches. To age 18 the liver carries **77.9%** of the fatal
hazard, cardiac 10.4%, vascular 11.7% — so a therapy achieving a
liver-specific HR of **0.386** moves all-cause hazard only to **0.522**.

## Axis 7 — the GALA curve is heterogeneity, not biology

No single trajectory fits GALA native liver survival: its hazard is strongly
front-loaded and every mechanistic trajectory here is not (best
single-trajectory fit SSE 3.3×10⁻², visibly wrong at age 5). Admitting a
**gamma frailty of variance 2.94** fits all three time points **65-fold
better** (SSE 5.0×10⁻⁴):

| age | model (population) | GALA |
|---|---|---|
| 5 yr | 67.7% | 66.8% |
| 10 yr | 52.6% | 54.4% |
| 18 yr | 41.3% | 40.3% |

A variance of 2.94 means the **standard deviation of individual hazard is 1.7×
its mean** — the notorious ALGS phenotype spread (same variant → infant
transplant or asymptomatic parent), measured rather than described.
Consequence: population curves must use the `_POP` outputs. Comparing the
individual-level `exp(-H)` outputs to a cohort study is a category error that
this file made once (reporting 1.5% 18-year native liver survival against a
true 40.3%) before it was caught.

---

## Calibration summary

| target | published | model |
|---|---|---|
| ASSERT sBA, drug ÷ placebo at wk 21–24 | 0.550 | **0.543** |
| ASSERT itch difference (drug − placebo) | −0.90 | **−0.896** |
| ASSERT sBA difference | −113 µmol/L | −108 µmol/L |
| ASSERT placebo-controlled itch slope | 0.0080 | **0.0083** |
| ICONIC week-48 sBA change | −96 µmol/L | −107 µmol/L |
| ICONIC week-48 ItchRO(Obs) reported | −1.6 | −1.49 |
| GALA native liver survival 5/10/18 yr | 66.8/54.4/40.3% | **67.7/52.6/41.3%** |
| GALA event-free 10/18 yr | 48.5/34.0% | 46.6/36.8% |
| daily-bolus vs continuous-input dosing | — | 4.0% on sBA, 5.7% on blockade |

---

## Scenarios (16)

| # | Scenario | Purpose |
|---|---------|------|
| 0 | Equilibrium verification | Confirms the 24-week trial reproduction is not an equilibration transient |
| 1 | Natural course (age 1→18) | Calibration against GALA |
| 2 | ICONIC reproduction (maralixibat 380 µg/kg/d) | Single-arm slope |
| 2b | Daily bolus vs continuous-infusion equivalence | Validates the continuous-input approximation |
| 3 | ASSERT reproduction — **both arms** (odevixibat 120 µg/kg/d) | Arithmetic proof of Axis 3 |
| 4 | Itch slope per bile acid | The 2.1-fold difference with vs without a control arm |
| 5 | Duct capacity sweep | Axis 2 — collision with GALA's cut-point |
| 6 | Φ_EHC identifiability | Quantification of FAILURE 1 |
| 6b | Bilirubin-anchor sensitivity | Quantification of FAILURE 2 |
| 7 | PEBD (surgical biliary diversion) vs IBAT inhibitor | The ceiling of complete blockade |
| 8 | Drug panel (24 weeks) | 9 regimens compared in the same patient |
| 9 | Micellar window / site of action | Axis 5 and its boundary case |
| 10 | Early vs delayed treatment start | Integral argument (exposure avoidance 39.7% → 17.4%) |
| 11 | Event-free survival (6 years) | The core test of Axis 4 |
| 12 | Exponent-contradiction sweep | No single n satisfies both |
| 13 | Competing-risk ceiling | Axis 6 |
| 14 | Non-responder phenotype | DPR 0.32 vs 0.15, full trajectory |
| 15 | Postnatal ductular regeneration | Decomposition of FAILURE 4 |

```r
source("algs_mrgsolve_model.R")
run_all_algs()          # run every scenario
shiny::runApp("algs_shiny_app.R")
```

Rendering:
```bash
dot -Tsvg algs_qsp_model.dot -o algs_qsp_model.svg
dot -Tpng -Gdpi=150 algs_qsp_model.dot -o algs_qsp_model.png
```

---

## What this model cannot do (reported, not patched)

1. **Φ_EHC is not individually identifiable from a single trial.** Duct
   capacity, the achieved ASBT blockade fraction, and synthetic reserve enter
   the steady-state bile-acid ratio only as one combined quantity. Refitting
   duct capacity while sweeping the assumed maximal blockade fraction from
   0.50→0.92, so that ASSERT is still reproduced exactly, moves JCAP by
   3.58→1.95 (1.8-fold) and Φ_EHC by 0.57→0.42 (1.36-fold). Φ_EHC is the more
   robust of the two, but neither is pinned down.
2. **The bilirubin boundary depends on a single assumed anchor.** The ASSERT
   cohort's mean baseline total bilirubin could not be found in the original
   abstract, so it was set to 3.5 mg/dL. Anchors of 3.0/3.5/4.0 move the 30%
   boundary to 4.58/5.42/6.26 and the 15% boundary to 8.93/10.70/12.48. This
   brackets the GALA cut-points but a point estimate is not evidence. **Only
   the ratio of the two boundaries (1.97) is anchor-independent**, which is
   why Axis 2 leads with the ratio.
3. **The bile-acid-independent pruritus fraction is assumed, not solved
   for.** The ASSERT placebo arm only quantifies the *non-drug* component
   (47%); it says nothing about what pruritus remains once bile acids reach
   zero. One mitigating factor: across every combination of Stevens exponent
   tried (0.75–1.05) and assumed fraction (8–18%) the residual comes out at
   12–17%, so the conclusion is less sensitive than the parameter. It is
   still an assumption.
4. **Spontaneous improvement is reproduced only halfway.** With the fibrosis
   term switched off, ductular regeneration alone drops bilirubin from
   3.50→2.20 between ages 1 and 18. In the full model fibrosis outpaces duct
   recovery and it rises instead, 3.50→4.48. Which path a real child follows
   is the central clinical question, and this model cannot settle it — no
   serial-biopsy DPR time series exists to fit the regeneration rate against.
   Serum bile acid barely moves either way (236.9→233.4), which is itself a
   prediction: partial ductular recovery improves **bilirubin long before
   bile acids**.
5. **Risk is calibrated on a cohort and applied to an individual.** The
   frailty variance needed to fit GALA is 2.94, which is very large, so the
   cohort curves are largely a statement about between-patient variance.
   Scenario 1 is a calibration check, not a prognosis.
6. **There is no cirrhotic pharmacology and no placebo-arm progression.**
   The same duct/synthesis structure is applied at FIB ≥ 3 without capturing
   the hepatocellular-failure physiology of decompensated ALGS liver. Also,
   the ASSERT placebo arm's bile acids **rose** by 22 µmol/L over 24 weeks,
   while this model's untreated arm is flat (−0.9), so it is missing
   whatever drives that rise, and the absolute change is overestimated by
   that much (−109 vs −90).

### 9 defects found by running the model (ones that would have changed the conclusions)

1. `$MAIN` reset the initial conditions on every run, silently overwriting
   the external `init()` call, so **every duct-capacity scenario reported
   the same baseline bilirubin of 3.5**, making Axis 2's threshold itself
   unmeasurable.
2. Scenario 1 compared individual-level `exp(-H)` survival against a cohort
   curve, reporting 18-year native liver survival as **1.5%** against the
   true 40.3%.
3. Φ_EHC originally included a hepatocyte→plasma→hepatocyte futile loop,
   pinning it at ~0.98 for every patient and destroying the quantity's
   meaning.
4. `EPS` and `THETA` are mrgsolve reserved words, so compilation failed.
5. An apostrophe inside a comment terminated the single-quoted R model
   string early.
6. `$TABLE` did not recompute at the output timepoint and instead captured
   local variables left over from the last `$ODE` derivative call.
7. The first fibrosis rate drove FIB to 2 within 400 days; the corrected
   version then froze it at 1.0 for 17 years.
8. The bolus-vs-continuous-dosing validation sampled at exactly integer-day
   multiples of the dosing interval, reporting a **39% discrepancy that was
   pure aliasing** (5.7% on a fine grid).
9. The cholestyramine arm bound bile acids only for micelle formation and
   not for ASBT uptake, making it look like a pure malabsorption drug with
   no bile-acid effect at all.

---

## Disclaimer

This model is intended for research and educational purposes. It has not
been clinically validated and is not a clinical decision-making tool. No
parameter should be used to guide patient management.
