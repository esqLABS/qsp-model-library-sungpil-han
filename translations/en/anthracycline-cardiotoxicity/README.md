# Anthracycline-Induced Cardiotoxicity (AIC) — QSP Model

> **One-line summary**: cumulative dose in mg/m² is not the exposure metric.
> Within the myocardium there is separately an axis that **reads peak
> concentration** (Top2b–DNA cleavage) and an axis that **reads AUC**
> (iron·ROS), and the resulting deficit likewise splits into one that is
> **reversible** (functional deficit) and one that is **irreversible**
> (cardiomyocyte loss). LVEF moves late because it is masked by compensatory
> hypertrophy, and the reversible share can only be recovered while fibrosis
> is still growing.

| File | Contents |
|------|------|
| [`aic_qsp_model.dot`](aic_qsp_model.dot) · [`.svg`](aic_qsp_model.svg) · [`.png`](aic_qsp_model.png) | Mechanistic map — 18 clusters · 174 nodes · 241 edges |
| [`aic_mrgsolve_model.R`](aic_mrgsolve_model.R) | mrgsolve QSP model — 32 ODEs · 14 scenarios · 11 analysis functions |
| [`aic_shiny_app.R`](aic_shiny_app.R) | Shiny dashboard — 8 tabs |
| [`aic_references.md`](aic_references.md) | 62 references (PMIDs verified) |
| [`aic_reference_check.py`](../../../anthracycline-cardiotoxicity/aic_reference_check.py) | Independent numpy/RK4 transcription + virtual population — **source of every figure below** |
| [`aic_reference_check_output.txt`](../../../anthracycline-cardiotoxicity/aic_reference_check_output.txt) | Full run output of the script above (A0–A11) |

---

## 0. The Four Things This Model Needs to Explain

Anthracycline-induced cardiotoxicity is usually described by **a single
number**: cumulative dose in mg/m², and a threshold near 450–550. That
description fails to explain four things actually encountered in
cardio-oncology practice.

1. **LVEF stays flat for months, then falls rapidly.** So the dose-incidence
   curve looks like a threshold, but there is no biological switch that flips
   on at 450 mg/m².
2. **Giving 480 mg/m² over 72 hours, or as a liposomal formulation, is far
   less toxic.** Both the cumulative dose and the antitumour AUC are
   unchanged.
3. **Troponin and GLS (global longitudinal strain) move months before LVEF
   does.**
4. **The same heart-failure drugs recover a substantial share of patients if
   started within 2 months of detection, but recover almost no one if
   started 6 months later.**

---

## 1. The Model's Two Structural Claims

### Claim 1 — The Two Damage Axes Read Different Exposure Metrics

Drug within the myocardium is split into **two compartments, not one**.

| Pool | Half-life | What it reads | What it drives |
|---|---|---|---|
| `CHF` fast nuclear/free pool | **0.03 d** | plasma **peak** | Top2b–DNA double-strand break (near-saturating Hill, n=2 → near-squared dependence) |
| `CH` slow residual bound pool | 9 d | plasma **AUC** | iron·ROS oxidative damage |
| `CHM` myocardial doxorubicinol | 28 d | plasma AUC (**accumulates cycle to cycle**) | iron·ROS + calcium disturbance |

Because ∫C²dt is dominated by the peak, **giving the same AUC slowly barely
engages the first axis.** This is why schedule and formulation matter.

The nuclear damage axis is not transient. DSBs create **`P53`, a
genotoxic-stress memory with a 28-day half-life**, which accumulates across
cycles and enters the death rate **multiplicatively** as `(1 + 1.2·P53)`.
Because the two axes combine multiplicatively, **removing just one still buys
more benefit than its share alone would suggest.**

### Claim 2 — Two Deficits, and the Mask That Hides Them

```
CONT = MYO · (1 − FUNC) · (1 + HYP) / (1 + 1.2·FIB)
LVEF = 62 · CONT^0.9
```

- `MYO` **the irreversible cardiomyocyte pool**. Adult cardiomyocyte turnover
  is 0.5–1%/year, so `KREG = 2×10⁻⁵/day` — what dies does not come back.
- `FUNC` **the reversible per-cardiomyocyte contractile failure** (SERCA2a
  downregulation, RyR2 leak, myofibrillar disarray). Recovery half-life 58
  days.
- `(1 + HYP)` is **the mask**. Surviving cardiomyocytes hypertrophy up to
  `HYP_target = min(0.32, 1.6 × deficit)`, so while the deficit stays below
  roughly 20% the compensation nearly offsets the loss and **LVEF does not
  move.** Once compensation saturates, LVEF falls sharply. **The apparent
  cumulative-dose threshold is not a property of the drug but the exhaustion
  of this reserve.**
- `GLS` has no compensation term (`GLS = 20.5·(1 − 0.9·FUNC − 1.6·cardiomyocyte
  loss − 0.35·FIB)`), and troponin is proportional not to cumulative damage
  but to the **death rate**. That is why both move first.
- Finally, **fibrosis closes the door**: `FIB` slows the rate of functional
  recovery to `KFOUT/(1 + 3·FIB)`. This single term is **the reversibility
  window**.

Taken together, the two claims mean therapeutic levers are not
interchangeable. **Upstream levers** (dexrazoxane, continuous infusion,
liposomal formulation, statins) act only on the axis they reach, and only
while the drug is in the body. **Downstream levers** (ACEi/ARB, beta-blockers,
ARNI) keep acting afterwards, but only on the reversible deficit and the
fibrosis lock. **There is no lever at all for the cardiomyocytes already
lost.**

---

## 2. Computed Results

Everything below is **output of these equations**, not a literature value.
It was generated with `aic_reference_check.py` (an independent numpy/RK4
transcription, virtual population 500–800 patients). Reference patient = 55
years old, BSA 1.8 m², baseline LVEF 62%, doxorubicin 60 mg/m² q3w. CTRCD =
LVEF fall of ≥10 points & <50% (ASE/ESC definition).

### Result 1 — The Dose-Incidence Curve Is Emergent, Not a Threshold

| Cumulative dose (mg/m²) | Cycles | CTRCD | LVEF<40% | 12-month LVEF fall | Mean nadir LVEF |
|---|---|---|---|---|---|
| 120 | 2 | 5.4% | 0.2% | 2.24 | 58.2 |
| 240 | 4 | 15.1% | 0.6% | 4.94 | 55.2 |
| 300 | 5 | 20.9% | 1.1% | 6.03 | 54.1 |
| 360 | 6 | 26.0% | 2.8% | 7.04 | 53.1 |
| 420 | 7 | 31.2% | 4.2% | 8.06 | 52.0 |
| 480 | 8 | 37.0% | 7.0% | 9.13 | 50.9 |
| 600 | 10 | 48.2% | 14.9% | 11.52 | 48.7 |

**Nothing** changes at 450 mg/m² inside the model. The curve steepens because
of the `Hill(ROS; n=3)` term in the death rate and because reserve runs out.
(Literature anchor — Von Hoff 1979: 300 mg/m² ~2%, 400 ~3–5%, 550 ~7–18%.)

### Result 2 — Cumulative Dose Is Not the Exposure Metric

480 mg/m² **fixed**, only schedule/formulation varied (reference patient):

| Schedule | Peak nuclear pool CHF | Residual pool AUC | Peak p53 | 24-month cardiomyocytes | 12-month LVEF | CTRCD (population) |
|---|---|---|---|---|---|---|
| IV push q3w | **10.19** | 26.2 | **0.883** | 0.927 | 55.1 | **37.0%** |
| Weekly 20 mg/m² ×24 | 3.40 | 26.2 | 0.302 | 0.962 | 58.3 | 20.1% |
| 72-hour continuous infusion | **0.60** | **26.7** | **0.105** | 0.972 | 59.4 | **14.2%** |
| Liposomal (PLD) | 0.09 | 13.3 | 0.004 | 0.998 | 62.2 | 2.5% |
| Push + dexrazoxane | 10.19 | 26.2 | 0.244 | 0.968 | 58.9 | 16.2% |

**Residual-pool AUC is essentially the same at 26.2 vs 26.7 vs 26.2, yet peak
nuclear exposure differs 17-fold and p53 memory 8-fold.** In other words,
**giving 480 mg/m² as a 72-hour infusion is less toxic than giving 240 mg/m²
as a push** (CTRCD 14.2% vs 15.1%). PLD also lowers the residual AUC, because
of its free-drug release fraction, so both axes fall together.

### Result 3 — What Moves First (480 mg/m², n=600)

| Marker | Patients reaching it | Median day reached |
|---|---|---|
| hs-cTnI > 14 ng/L | 100% | **24 days** |
| GLS relative fall ≥15% | 88.0% | **112 days** |
| CTRCD (LVEF-based) | 36.7% | **165 days** |

- GLS → CTRCD median lead time **79 days**
- cTnI → CTRCD median lead time **146 days**
- Among patients whose cTnI exceeded 40 ng/L during chemotherapy, CTRCD was
  51.5% vs 17.3% in those who did not → **negative predictive value 82.7%**.
  The value of a troponin-guided strategy lies not in positive prediction but
  in **the reassurance a negative result gives**.

### Result 4 — Dexrazoxane's Mechanism and Its Dose-Sparing Equivalence

Modelling Top2b as a protein state variable with **drug-induced proteasomal
degradation (25/day) + resynthesis with a 2-day half-life**, dexrazoxane at
a 10:1 ratio gives:

- p53 memory **0.883 → 0.244 (−72%)**
- CTRCD **37.0% → 16.2%** (RR **0.44**), LVEF<40% **7.0% → 0.4%** (RR **0.06**)
- **Equivalent cardiotoxic dose: 540 mg/m² (2.2× the unprotected 240 mg/m²)**

The contribution of the ADR-925 iron-chelation term is small — exposure
lasts hours while the residual pool's assault continues for weeks. This is
consistent with preclinical observations that a Top2-inactive analogue
(ICRF-161) gives weak protection.

### Result 5 — The Reversibility Window Is a Fibrosis Clock

In patients who developed CTRCD after 480 mg/m², only the start time of
enalapril + carvedilol is varied:

| Start day | LVEF at start | LVEF at 12 months | Recovery fraction | FIB at start |
|---|---|---|---|---|
| 60 | 57.73 | 49.33 | **55.9%** | 0.024 |
| 90 | 55.01 | 48.92 | 53.8% | 0.066 |
| 120 | 51.61 | 48.53 | 49.0% | 0.126 |
| 180 | 44.55 | 47.84 | 43.4% | 0.266 |
| 270 | 42.71 | 47.46 | 41.3% | 0.379 |
| 365 | 43.10 | 47.21 | 38.5% | 0.403 |
| 540 | 43.40 | 46.83 | **30.8%** | 0.406 |

Recovery fraction tracks **not time itself but the fibrosis level at that
time.** The direction agrees with Cardinale 2010 (64% within 2 months → ~0%
after 6 months), but **the slope is shallower than in the literature** —
noted here as a model limitation.

### Result 6 — Upstream and Downstream Protection Are Not Interchangeable

480 mg/m², same drugs, **only the start time** differs (n=600):

| Strategy | CTRCD | 12-month LVEF | 24-month LVEF |
|---|---|---|---|
| No protection | 36.7% | 52.70 | 51.31 |
| Statin, **during** chemotherapy (day 0) | **17.5%** | 56.77 | 54.58 |
| Statin, **after** chemotherapy (day 170) | **35.5%** | 53.01 | 51.54 |
| ACEi+BB, during chemotherapy | 10.2% | 57.08 | 55.96 |
| ACEi+BB, after chemotherapy | 20.7% | 55.98 | 55.02 |
| ARNI, after chemotherapy | 19.5% | 57.20 | **56.42** |
| SGLT2i, during chemotherapy | 28.0% | 54.58 | 52.68 |

**Statins lose almost all their benefit when delayed** (RR 0.48 → 0.97).
Their target is drug-driven ROS generation, so once the drug is gone, so is
the target. **ACEi+BB retain half their benefit even when delayed** (RR 0.28
→ 0.56), because their targets, the reversible deficit and the fibrosis
lock, persist after chemotherapy ends. ARNI started after chemotherapy
performs best by 24-month LVEF. Even under the same label of
"cardioprotectant", **when a drug is given matters as much as what is
given.**

### Result 7 — Trastuzumab Interaction Is Additive When Sequential, Supra-Additive Only When Concurrent

| Strategy | CTRCD | LVEF<40% | Maximum LVEF fall | Recovery by 24 months | 24-month cardiomyocytes |
|---|---|---|---|---|---|
| 240 mg/m² alone | 14.3% | 0.8% | 6.08 | **−0.90** | 0.961 |
| Trastuzumab alone | 4.3% | 0.2% | 3.68 | **+0.66** | 0.992 |
| Sequential (from day 84, 1 year) | 28.0% | 6.2% | **9.46** | +2.48 | 0.951 |
| Concurrent (from day 0, 1 year) | 54.5% | 23.3% | **15.25** | +4.25 | 0.898 |

Additive prediction = 6.08 + 3.68 = **9.76 points**.

- **Sequential: observed 9.46 points → additive** (excess −0.30). There is
  no interaction.
- **Concurrent: observed 15.25 points → supra-additive** (excess **+5.49
  points**).

The interaction lies not in the drug combination itself but in **the
temporal overlap of exposures**. ErbB2 blockade slows DSB repair by
`÷(1+0.9·ETR)` and raises the death rate by `×(1+0.65·ETR)`, but this only
matters if it happens **at the very moment nuclear damage arrives**. This
structure is why sequential dosing is the standard, and clinical practice
avoids concurrent administration.

The "recovery by 24 months" column replaces the traditional Type I/II
dichotomy. Anthracycline-only damage does not recover (−0.90), while damage
with a large trastuzumab contribution recovers after washout (+0.66 ~
+4.25). **Reversibility is not two diseases but the ratio of the two
deficits.**

### Result 8 — Decomposing the LVEF Deficit into Components (480 mg/m², reference patient)

| Day | LVEF | Cardiomyocytes | FUNC | FIB | HYP | Size of the mask | Share still recoverable |
|---|---|---|---|---|---|---|---|
| 60 | 60.37 | 0.991 | 0.039 | 0.009 | 0.031 | +1.64 | +2.23 |
| 120 | 58.46 | 0.966 | 0.077 | 0.045 | 0.108 | +5.17 | +4.40 |
| 170 | 56.56 | 0.942 | 0.099 | 0.092 | 0.181 | **+7.85** | **+5.56** |
| 240 | 56.39 | 0.936 | 0.077 | 0.144 | 0.221 | +9.29 | +4.23 |
| 365 | 55.15 | 0.934 | 0.032 | 0.173 | 0.173 | +7.37 | +1.65 |
| 730 | 52.84 | 0.927 | 0.002 | 0.198 | 0.119 | +5.09 | +0.10 |

("Size of the mask" = actual LVEF − LVEF with compensatory hypertrophy
removed. "Share still recoverable" = LVEF with FUNC reset to 0 − actual
LVEF.)

- **The mask is worth 7.85 LVEF points at the end of chemotherapy (day
  170).** Without compensatory hypertrophy, this patient's LVEF would be
  **48.7**, not 56.6.
- **The recoverable share disappears over time**: 5.56 points at day 170 →
  1.65 points at day 365 → 0.10 points at day 730. This is Result 5's
  reversibility window seen at the component level.
- The point where LVEF has fallen by 5 points is **day 157**, by which time
  **5.2%** of cardiomyocytes are already gone and **|GLS| has fallen by
  19.7%**. This is the quantitative form of the statement that "LVEF is a
  lagging indicator".

### Result 9 — Metabolic Phenotype (CBR1/AKR1C3) Shifts the Whole Curve

| FM (doxorubicinol formation fraction) | Peak myocardial metabolite | Peak ROS | 12-month LVEF (reference patient) | CTRCD (population) |
|---|---|---|---|---|
| 0.12 | 0.260 | 0.859 | 57.94 | 22.6% |
| 0.18 | 0.389 | 1.007 | 56.56 | 28.8% |
| 0.25 (reference) | 0.541 | 1.182 | 55.15 | 36.8% |
| 0.34 | 0.736 | 1.415 | 53.69 | 45.4% |
| 0.45 | 0.974 | 1.709 | 52.03 | **60.0%** |

A formation fraction of 0.12 → 0.45 shifts CTRCD **22.6% → 60.0%**.
Because the metabolite is a pool that **accumulates** with a 28-day
half-life, a small difference in formation fraction is amplified across
cycles. This is the axis along which outcomes diverge between patients on
the same regimen and the same cumulative dose.

### Result 10 — Phenotype Overrides Dose (High Risk: Age ≥70 + Hypertension + Prior Thoracic Radiation)

| Risk group | Regimen | CTRCD | LVEF<40% | 12-month LVEF |
|---|---|---|---|---|
| Standard | 240 mg/m² | 16.0% | 0.3% | 57.00 |
| Standard | 480 mg/m² | 38.7% | 6.0% | 52.86 |
| Standard | 240 mg/m² + dexrazoxane | 7.9% | 0.1% | 59.20 |
| **High risk** | **240 mg/m²** | **69.1%** | **23.6%** | **46.42** |
| High risk | 480 mg/m² | 90.3% | 60.0% | 39.17 |
| High risk | 240 mg/m² + dexrazoxane | 53.6% | 14.6% | 49.42 |

**240 mg/m² in a high-risk patient (69.1%) is far worse than 480 mg/m² in a
standard-risk patient (38.7%).** In other words, knowing the phenotype
matters more than respecting a cumulative-dose ceiling. Moreover, even
adding dexrazoxane in the high-risk group (53.6%) does not reach the level
of unprotected 480 mg/m² in the standard-risk group (38.7%) — protective
agents reduce risk but do not change the phenotype. In the model, the
high-risk phenotype is not a separate "risk score" but a combination of
1.9× death sensitivity, 1.5× functional-damage sensitivity, 1.6× fibrotic
propensity, 0.4× regenerative capacity, and a −4-point baseline LVEF.

---

## 3. Parameter Calibration Anchors

| Reference | Target value | Model value |
|---|---|---|
| Von Hoff 1979 | Clinical heart failure 300 mg/m² ~2%, 400 ~3–5%, 550 ~7–18% | 300 → 1.1%, 420 → 4.2%, 480 → 7.0%, 600 → 14.9% |
| Cardinale 2015 | CTRCD 9%, 98% within 1 year, median 3.5 months | 15.1% at 240 mg/m², median 165 days |
| Cochrane (dexrazoxane) | Heart failure RR 0.29 | Heart failure RR 0.06 · CTRCD RR 0.44 |
| Cochrane (continuous infusion) | Heart failure RR 0.27 | Heart failure RR 0.03 · CTRCD RR 0.39 |
| Cochrane (liposomal) | Heart failure RR 0.20 | CTRCD RR 0.07 — **model is optimistic** |
| STOP-CA 2023 | Atorvastatin ≥10% decline 9% vs 22% | CTRCD RR 0.51 |
| Cardinale 2010 | Treatment within 2 months recovers 64%, after 6 months ~0% | 55.9% → 30.8% (direction matches, slope shallower) |
| PK | AUC 2.0 µg·h/mL @ 60 mg/m², terminal t½ 20–48 h | 1.99 µg·h/mL, matched by design |

---

## 4. Model Structure Summary

**32 ODEs** = PK/exposure 13 + disease·cardiac·endpoint 19. Time unit = day.

| Block | State variables |
|---|---|
| Doxorubicin PK | `A1` `A2` `A3` `ALIP`(liposomal carrier) |
| Metabolite | `AM` `AMP` (doxorubicinol) |
| Myocardial exposure | `CHF`(peak) `CH`(AUC) `CHM`(accumulated metabolite) |
| Protectant·antibody PK | `ADEX` `T2B`(Top2b protein) `TR1` `TR2` |
| Oral cardioprotectants | `EACE` `EBB` `ESTA` `EARNI` `ESGLT` |
| Nuclear damage axis | `DSB` `P53` |
| Redox axis | `ROS` `MITOD` `LIP` |
| The two deficits | `MYO`(irreversible) `FUNC`(reversible) |
| Remodelling | `FIB` `HYP`(mask) `NH` |
| Biomarkers·integrals | `TNI` `BNP` `AUCH` `CUMKILL` |

**14 scenarios**: 240/360/480 mg/m² · 72-hour continuous infusion · weekly
split-dosing · liposomal · dexrazoxane · statin (early/delayed) · ACEi+BB ·
full combination · trastuzumab (sequential/concurrent/alone).

**8-tab Shiny app**: patient profile · PK and myocardial exposure · damage
axes · the two deficits and the mask · clinical endpoints · scenario
comparison · biomarkers · reversibility window.

### Usage

```r
library(mrgsolve)
mod <- mread("aic", "aic_mrgsolve_model.R")
out <- sim_scenario(mod, "DOX480")        # one scenario
cmp <- compare_schedules(mod)             # Result 2
win <- reversibility_window(mod)          # Result 5
shiny::runApp("aic_shiny_app.R")          # dashboard
```

```bash
python3 aic_reference_check.py            # regenerates every figure above (about 15 minutes)
python3 aic_reference_check.py --quick    # reference-patient scenario only
```

---

## 5. Where the Model Overstates Benefit

An honest record. Two items should be treated as an **upper bound**.

- **Liposomal doxorubicin**: model CTRCD RR ≈ 0.07 vs Cochrane clinical
  heart-failure RR 0.20. Because free-drug peak is the dominant driver of
  the nuclear damage axis, the benefit of encapsulation is overestimated.
- **Combined protectants** (dexrazoxane + statin + ACEi/BB): the model
  lowers CTRCD at 480 mg/m² to **0.0%** (0 of 800) and the 12-month LVEF
  fall to 0.42 points. No randomised trial has tested this combination, so
  it should be read as an upper bound, not a prediction.

Also, the reversibility window's slope is shallower than the literature
(Result 5), and troponin elevation appears in almost all patients, giving
lower specificity than seen clinically.

---

## ⚠️ Disclaimer

This is a **qualitative/semi-quantitative QSP model** for education and
research purposes. It was built from published literature and clinical
trial data but has not been independently validated or certified, and
**must not be used directly for clinical decision-making, prescribing, or
regulatory submission.** Parameters are illustrative approximations;
separate fitting and validation against real patient data is required.
