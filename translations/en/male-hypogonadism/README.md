# Male Hypogonadism (MHG) QSP Model

**Male Hypogonadism · MHG** — a quantitative systems pharmacology (QSP) model
that unites the hypothalamic–pituitary–gonadal (HPG) axis, Vermeulen binding
equilibrium, intratesticular testosterone, erythropoiesis, bone, and body
composition into a single mechanistic model.

> ⚠️ **For educational and research purposes only.** Do not use for clinical decision-making, prescribing, or regulatory submission.

---

## The Claim

**What is measured and what acts are not the same value, and the map between them is nonlinear.**

What is measured in the clinic is total testosterone (TT), but what actually
reaches the receptor is free testosterone (FT). Between the two lies a
mass-action equilibrium governed by SHBG, and SHBG binding **saturates**. So
FT is a **convex function** of TT (SHBG 35 nmol/L, albumin 4.3 g/dL):

| Total T (ng/dL) | 100 | 300 | 600 | 1000 | 1500 | 2000 |
|---|---|---|---|---|---|---|
| Free T (pg/mL) | 17.5 | 56.0 | 122.9 | 228.2 | 381.1 | 550.8 |
| Free fraction (%) | 1.75 | 1.87 | 2.05 | 2.28 | 2.54 | **2.75** |

The free fraction **rises** together with total T. Three consequences follow
from this single line, and the model's analysis functions print each one as a
**number**, not as prose.

### ① The diagnostic threshold is frame-dependent

A single total-T threshold (300 ng/dL) is **not** a single free-T threshold
(65 pg/mL). Output of `MHG_diagnostic_frame()`:

| SHBG (nmol/L) | FT (pg/mL) at TT=300 | TT (ng/dL) at which FT=65 | Verdict |
|---|---|---|---|
| 15 | 83.5 | **237** | TT threshold overdiagnoses |
| 25 | 67.3 | 290 | Nearly matches |
| 35 | 56.0 | 344 | TT threshold underdiagnoses |
| 55 | 41.5 | 450 | TT threshold underdiagnoses |
| 90 | 28.2 | **635** | TT threshold underdiagnoses |

The same free-T threshold moves **237 → 635 ng/dL, a 2.7-fold range**, when
expressed as total T. This means there is a band in which what decides
whether a man is hypogonadal is not the gonad but the **binding protein**.
The direction cuts both ways — obese/insulin-resistant men (low SHBG) are
overdiagnosed, and older men (high SHBG) are underdiagnosed.
`MHG_shbg_paradox()` simulates the two patients side by side.

### ② The waveform matters far less than expected — the model falsifies its own hypothesis

By Jensen's inequality, E[FT(TT)] ≥ FT(E[TT]), so at the same mean total T, a
regimen with a larger peak-to-trough amplitude delivers a higher
time-averaged free T. This model was originally built on the hypothesis that
this convexity, **composed with a convex erythropoietic response**, would
explain the roughly threefold gap in polycythaemia between intramuscular and
transdermal therapy. **Running the model showed it does not.**

The arithmetic that `MHG_convexity_decomposition()` actually prints:

| Component | Model output | Interpretation |
|---|---|---|
| (a) Convex **binding** | gel +0.3% → IM q2wk **+4.4%** | Real, and grows with amplitude, exactly as the mathematics requires |
| (b) Composition with the convex **erythropoietic response** | **≈ 0 or negative** | The peak crosses the EC50 and enters the **concave** upper part of the Hill curve. A response that is convex from below is concave from above |
| (c) **Waveform effect** at matched dose | Hct **~0.1 point** | The actual difference at equal weekly dose with 2× amplitude (100 mg weekly vs 200 mg every other week) |
| (d) Applying a fixed threshold (Hct>54%) to a **population distribution** | Incidence **~2.5-fold** | A ~1-point shift in the mean is amplified into a several-fold difference in incidence |

In other words, the IM–transdermal gap in the literature is **mostly a dose
difference, and the rest is manufactured by the decision rule.** (b) is the
point where the model falsifies the very hypothesis it was built to express,
and it is printed as is, not tuned away.

> **Exposed assumption**: `EC50_HEPC = 300 pg/mL` (set above the normal
> free-T range of 50-210) is **a calibrated value, not a measured one.**
> `MHG_hepcidin_sensitivity()` tabulates how much of the above conclusion
> survives when this value is moved from 150 to 600.

### ③ Serum T and intratesticular T are different variables, and only one is measured

Intratesticular testosterone (ITT) is about 35-100 times serum (about 700
nmol/L), and spermatogenesis stops once ITT falls below about 20-30% of
normal. Exogenous testosterone normalises serum T while switching off LH and
collapsing ITT to a few percent.

`MHG_ITT_collapse()` puts both values in the same table to show the
dissociation. **The rows where serum T is normal but ITT is single-digit are
exactly the regimens that "replace the hormone while rendering the patient
infertile."** Only regimens that leave LHCGR signalling intact (hCG
co-therapy, clomiphene) keep ITT above threshold. Replacement and restoration
are different interventions.

### ④ Bonus: much of "testosterone's effect" is actually oestradiol's effect

`MHG_finkelstein()` reproduces the NEJM 2013 design (GnRH agonist background
+ graded T gel ± anastrozole). Lean mass and strength track T, but **fat mass
tracks E2.** So suppressing aromatase to prevent gynaecomastia is not free —
`MHG_aromatase_cost()` bills that cost in terms of lumbar spine trabecular
vBMD.

---

## Files

| File | Contents |
|------|------|
| [`mhg_qsp_model.dot`](mhg_qsp_model.dot) | Mechanistic map source — **224 nodes · 268 edges · 22 clusters** |
| [`mhg_qsp_model.svg`](mhg_qsp_model.svg) | Vector map (scalable) |
| [`mhg_qsp_model.png`](mhg_qsp_model.png) | Raster map (150 dpi) |
| [`mhg_mrgsolve_model.R`](mhg_mrgsolve_model.R) | **49 ODE compartments** · 8 T formulations + 5 non-androgen drugs · 14 scenarios · 11 analysis functions |
| [`mhg_shiny_app.R`](mhg_shiny_app.R) | 10-tab interactive dashboard |
| [`mhg_references.md`](mhg_references.md) | **78 references** (PubMed links) + "What the model does not reproduce" appendix |

---

## Mechanistic Map

[![MHG QSP map](mhg_qsp_model.png)](mhg_qsp_model.svg)

22 clusters:

| # | Cluster | # | Cluster |
|---|---|---|---|
| ① | Hypothalamus — GnRH pulse generator (KNDy) | ⑫ | Erythropoiesis (hepcidin · EPO · iron) |
| ② | Anterior pituitary — gonadotrophs | ⑬ | Bone (E2-dominant pathway) |
| ③ | Leydig cell steroidogenesis | ⑭ | Muscle · body composition |
| ④ | Sertoli cell · spermatogenesis | ⑮ | Central nervous system — libido · erection · mood |
| ⑤ | **Circulating binding equilibrium (mathematical centre)** | ⑯ | Prostate · LUTS (saturation model) |
| ⑥ | SHBG regulation | ⑰ | Cardiovascular · safety |
| ⑦ | Peripheral metabolism (5α-reduction · aromatisation) | ⑱ | Testosterone formulation PK |
| ⑧ | AR / ER receptor signalling | ⑲ | Non-androgen strategies |
| ⑨ | Primary (testicular) causes | ⑳ | Axis suppression by exogenous androgen |
| ⑩ | Secondary (central) causes | ㉑ | Diagnosis · monitoring |
| ⑪ | Functional — obesity · inflammation vicious cycle | ㉒ | Clinical endpoints · key trials |

The pink notes (⚑) on the map are the model's core quantitative claims, each
corresponding to one analysis function in the R model.

---

## mrgsolve Model (49 ODEs · 132 parameters)

```r
source("mhg_mrgsolve_model.R")
MHG_run_all()                      # run all 11 analysis functions

d <- MHG_scenario_im_q2wk()        # IM cypionate 200 mg q2wk
MHG_plot_overview(d)
MHG_plot_waveforms()               # compare waveforms across formulations
```

### Compartment Structure

| Module | Compartments |
|------|------|
| Testosterone PK | `DEP_IM` `DEP_TU` `DEP_SC` `DEP_GEL` `DEP_ORAL` `DEP_PEL` `DEP_NAS` `CENT` `PERIPH` |
| Non-androgen drug PK | `HCG_D/C` `FSHD_D/C` `CLO_D/C` `ANA_D/C` |
| HPG axis | `GNRHD` `LH` `FSH` `LEYCAP` `ITT` `INHB` |
| Hormones · binding | `SHBG` `E2` `DHT` |
| Spermatogenesis (74+14-day transit chain) | `SG1` `SG2` `SG3` `SG4` `EPID` |
| Erythropoiesis | `HEPC` `EPO` `PROG_E` `RETIC` `RBC` |
| Bone | `SCLERO` `OB` `OC` `BMD_TR` `BMD_CO` `RSP` (remodelling space) |
| Body composition · other | `LEAN` `FAT` `PSA` `LIBIDO` `VITAL` |
| Exposure integrals | `CUMFT` `CUMDRIVE` |

### Patient Archetypes — untreated baseline values

These are values obtained by actually running the model, and reproduce the
laboratory pattern expected in the clinic for each aetiology.

| Archetype | TT (ng/dL) | FT (pg/mL) | SHBG | LH | FSH | ITT (%) | Pattern |
|---|---|---|---|---|---|---|---|
| `PT_ORGANIC` primary testicular failure | 232 | 38.8 | 40 | **14.8** | 13.6 | 29 | Hypergonadotrophic |
| `PT_FUNCTIONAL` obesity-related functional | 300 | 58.9 | **32** | 1.5 | 1.9 | 47 | Hypogonadotrophic · low SHBG |
| `PT_ELDERLY` ageing-related | 272 | 40.4 | **49** | 6.1 | 6.2 | 31 | Mixed · high SHBG |
| `PT_SECONDARY` pituitary | 276 | 51.1 | 35 | 1.3 | 1.7 | 40 | Hypogonadotrophic |
| `PT_KLINEFELTER` 47,XXY | 172 | 30.7 | 35 | 15.7 | **31.4** | 22 | FSH ≫ LH |
| `PT_OPIOID` opioid | 261 | 44.6 | 40 | 1.1 | 1.5 | 34 | Hypogonadotrophic · reversible |

### Treatment Scenarios (14)

| # | Scenario | Function |
|---|---|---|
| 1 | Untreated natural history of functional hypogonadism (3 years) | `MHG_scenario_natural()` |
| 2 | IM cypionate 200 mg q2wk | `MHG_scenario_im_q2wk()` |
| 3 | IM cypionate 100 mg weekly (same weekly dose) | `MHG_scenario_im_weekly()` |
| 4 | Subcutaneous auto-injector 75 mg weekly | `MHG_scenario_sc()` |
| 5 | Transdermal gel 1.62% 81 mg/day | `MHG_scenario_gel()` |
| 6 | Oral undecanoate 237 mg BID | `MHG_scenario_oral()` |
| 7 | IM undecanoate 1000 mg q12wk | `MHG_scenario_tu_im()` |
| 8 | Subcutaneous pellet 750 mg q4 months | `MHG_scenario_pellet()` |
| 9 | hCG monotherapy 1500 IU three times weekly (fertility preservation) | `MHG_scenario_hcg_mono()` |
| 10 | Gel + hCG 500 IU EOD (ITT rescue) | `MHG_scenario_t_plus_hcg()` |
| 11 | Clomiphene 25 mg daily | `MHG_scenario_clomiphene()` |
| 12 | Weight loss alone (no drug) | `MHG_scenario_weight_loss()` |
| 13 | Opioid-induced → discontinued at day 180 | `MHG_scenario_opioid()` |
| 14 | Discontinuation after 3 years of TRT — axis and sperm recovery | `MHG_scenario_cessation()` |

### Analysis Functions (11)

| Function | Output |
|------|------|
| `MHG_free_T_nomogram()` | Free T, free fraction, and convexity verification over a TT × SHBG grid |
| `MHG_diagnostic_frame()` | Table showing the total-T threshold moving 2.7-fold with SHBG |
| `MHG_shbg_paradox()` | Obese vs elderly patients — TT and FT verdicts diverge |
| `MHG_convexity_decomposition()` | **Core ledger**: four-way decomposition into (a) binding · (b) response · (c) dose-matched waveform · (d) threshold |
| `MHG_hepcidin_sensitivity()` | Sensitivity of the conclusion to `EC50_HEPC` |
| `MHG_ITT_collapse()` | Serum T and ITT in the same table — replacement vs restoration |
| `MHG_recovery_curve()` | Sperm recovery after discontinuation (validated against Liu 2006) |
| `MHG_finkelstein()` | Dissociating T vs E2 contributions (reproducing the NEJM 2013 design) |
| `MHG_aromatase_cost()` | Converts the cost of aromatase inhibition into vBMD |
| `MHG_formulation_ledger()` | One-year comparison table across 7 formulations × 12 endpoints |
| `MHG_trial_ledger()` | Model vs published endpoints — **including reproduction failures** |

---

## Shiny Dashboard (10 tabs)

```r
setwd("male-hypogonadism"); shiny::runApp("mhg_shiny_app.R")
```

① Patient profile · ② Binding equilibrium (SHBG nomogram) · ③ Frame-dependence
of the diagnostic threshold · ④ PK waveforms · ⑤ Polycythaemia convexity
ledger · ⑥ Intratesticular T · fertility · ⑦ Bone · body composition · ⑧
Biomarkers · ⑨ Scenario comparison · ⑩ Clinical trial comparison table

The sidebar lets you directly move 6 patient types, 12 regimens (up to 3
compared simultaneously), age · SHBG · insulin resistance · fat mass, and
**the `EC50_HEPC` assumption itself**.

---

## What the Model Does **Not** Reproduce

An honest QSP model must state where it is wrong. The most important item:

> **TRAVERSE fracture substudy (NEJM 2024): clinical fracture 3.50% vs
> 2.46%, HR 1.43.** **More fractures occurred** on the testosterone arm. The
> model is calibrated to reproduce the T-Trials rise in lumbar spine
> trabecular vBMD (published +7.5%, model +5.7~+7.7%), so it structurally
> **cannot** predict an increase in fractures. This is not something a
> parameter adjustment can fix — it is **a failure of the assumption that
> treats BMD as a surrogate for fracture itself**, and `MHG_trial_ledger()`
> prints it without hiding it.
>
> Furthermore, the very mechanism by which the model generates that BMD rise
> reinforces this warning. Much of the rise is not new bone but **recovery
> of the remodelling space** — the as-yet-unfilled resorption cavities
> created by the high bone turnover of the hypogonadal state are filled in
> as androgen restoration slows bone turnover (the same mechanism as the
> early BMD rise seen with antiresorptive agents). There is no guarantee
> that the gain the densitometer records is a gain in strength.

Beyond this: TRAVERSE MACE (no event process in a deterministic model), LH
pulsatility (simplified to mean concentration), SHBG allosteric binding (the
simple Vermeulen model is used instead), and inter-individual variability (AR
CAG repeat length etc. not implemented). See the
[appendix in `mhg_references.md`](mhg_references.md#appendix-what-this-model-does-not-reproduce)
for the full list.

---

## Calibration Targets

| Endpoint | Source | Role in the model |
|------|------|-----------------|
| Production rate 6 mg/day, MCR ~1000 L/day | Southren 1965 | Anchors the `KSPILL × ITT0` identity |
| K_SHBG 1×10⁹, K_Alb 3.6×10⁴ M⁻¹ | Vermeulen 1999 | The binding equilibrium itself |
| ITT ≈ 35-100 times serum | Jarow 2001 · Roth 2010 | `ITT0 = 700 nmol/L` |
| T 200 mg weekly dosing → ITT −94% | Coviello 2005 | `HCG_POT` · `ITT50_S` |
| Spermatogenic cycle 74 days | Heller & Clermont 1964 | `TAU_SPG` |
| Recovery: 67% at 6 months · 90% at 12 months | Liu 2006 (Lancet) | `MHG_recovery_curve()` |
| Lumbar spine trabecular vBMD +7.5% (1 year) | T-Trials Bone (JAMA IM 2017) | `KFORM`/`KRES` |
| Fat-mass increase attributed to E2 deficiency | Finkelstein (NEJM 2013) | `SFAT_E` |
| IM vs transdermal polycythaemia ~3-fold | Ohlander 2018 | `EC50_HEPC` (**calibrated**) |
| MACE HR 0.96 | TRAVERSE (NEJM 2023) | Not attempted to reproduce (stated) |
| Fracture HR 1.43 | TRAVERSE fracture (NEJM 2024) | **Reproduction failure (stated)** |

---

## Reproduce

```bash
# render the map
dot -Tsvg mhg_qsp_model.dot -o mhg_qsp_model.svg
dot -Tpng -Gdpi=150 mhg_qsp_model.dot -o mhg_qsp_model.png

# run the model
Rscript -e 'source("mhg_mrgsolve_model.R"); MHG_run_all()'
```

Required packages: `mrgsolve`, `shiny`. Graphviz for map rendering.

---

> **Disclaimer** — This model is an educational and research QSP model built
> from public literature and clinical trial data. It has not been
> independently validated or certified, and must not be used for actual
> clinical decision-making.
