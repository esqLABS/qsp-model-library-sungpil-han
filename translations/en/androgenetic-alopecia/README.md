# Androgenetic Alopecia (Male-Pattern · Female-Pattern Hair Loss) — QSP Model

**Hair loss is not a disease of losing hairs. It is a disease of follicles
shrinking, and there is exactly one moment in the whole cycle when size can
change — the moment of escaping telogen.**

This model's reference patient starts untreated at age 18 with an intact
scalp and runs until reaching the baseline value of the finasteride phase-3
trial (**876 hairs inside a 5.07 cm² circle**, Kaufman 1998, PMID 9777765).
At that point he still has **1006 of his original 1268 follicles.** Of the
392 "missing" follicles, only 262 are truly gone; **130 are living vellus
hairs.** Every drug appearing in this file is competing for those 130.

| File | Contents |
|------|------|
| [`aga_qsp_model.dot`](../../../androgenetic-alopecia/aga_qsp_model.dot) · [SVG](../../../androgenetic-alopecia/aga_qsp_model.svg) · [PNG](../../../androgenetic-alopecia/aga_qsp_model.png) | Mechanistic map — 162 nodes / 198 edges / 17 clusters |
| [`aga_mrgsolve_model.R`](../../../androgenetic-alopecia/aga_mrgsolve_model.R) | mrgsolve ODE model — 50 compartments (drug PK 12 + enzymes 8 + endocrine 4 + signalling 4 + systemic markers 5 + follicle 17), 18 scenarios |
| [`aga_shiny_app.R`](../../../androgenetic-alopecia/aga_shiny_app.R) | Shiny dashboard — 10 tabs |
| [`aga_references.md`](aga_references.md) | 73 references (every PMID individually verified via PubMed E-utilities) |

---

## What This Model Claims (Thirteen Points, All Numerically Derived)

### 1. The disease runs on **cycle time**, not the calendar

Size change happens only at the moment of exiting telogen. So "how much
worse in a year" is really a question of "how many cycles turn in a year."
The model's cycle lengths (androgen-driven ARD at 0.95, with anagen
shortened by 30%):

| Grade | Anagen T_A | Cycle length | Cycles in 5 years |
|---|---|---|---|
| Terminal | 1000 d | **832 d** | 2.19 |
| Intermediate | 560 d | 517 d | 3.53 |
| Small intermediate | 250 d | 296 d | 6.17 |
| Vellus | 70 d | **167 d** | **10.92** |

A terminal hair only gets **2** chances to shrink in five years, so hair
loss takes decades; a vellus hair gets **11** chances to come back, so
treatment response shows up **within a year**. Both facts come from the
same number.

### 2. 100 hairs a day is not an observation — it is an **identity**

> Daily shed = N_total × f_telogen / T_telogen = 100 000 × 0.090 / 100 d =
> **90/day**

The model's normal scalp, at a telogen fraction of 9.0%, produces exactly
this value. The rise in telogen fraction to **23.6%** in AGA patients is
not a separately added term either — it falls out automatically because
smaller grades have shorter anagen phases. **AGA's "increased shedding" is
a consequence of miniaturisation, not a separate disease process, and the
telogen % on a trichogram is effectively a measurement of anagen length.**

### 3. Why finasteride's dose-response curve is flat: **the target integrates, not the plasma**

Finasteride's plasma half-life is 6 hours. But the half-life of the
enzyme–NADP–dihydrofinasteride adduct it forms is 30 days, and the turnover
half-life of the 5AR protein itself is 2 days. Because capture is
effectively irreversible, the steady-state free enzyme is

> E_ss = k_deg / (k_deg + ⟨k_on × C⟩)

i.e. it is a function of **average concentration**, not peak concentration.
As a result, within a single day:

| State variable | Trough | Peak | Peak/trough |
|---|---|---|---|
| Plasma finasteride FINC | 0.70 ng/mL | 6.28 | **8.98×** |
| Free scalp 5AR-II | 0.0164 | 0.0376 | 2.29× |
| **Scalp DHT** | 0.3589 | 0.3673 | **1.02×** |
| AR nuclear signal ARN | 0.3818 | 0.3845 | 1.01× |

While plasma swings 9-fold, the signal driving the disease moves **2%.**

This reproduces the 5-point dose-response curve of Drake 1999
(PMID 10495374) (day 42, model vs observed):

| Dose | Scalp DHT model | Observed | Serum DHT model | Observed |
|---|---|---|---|---|
| 0.01 mg | −18.3% | −14.9% | −20.1% | — |
| 0.05 mg | −43.1% | −61.6% ✗ | −47.3% | −49.5% |
| 0.2 mg | −57.5% | −56.5% | −63.1% | −68.6% |
| **1 mg** | **−63.6%** | **−64.1%** | **−70.3%** | **−71.4%** |
| 5 mg | −67.4% | −69.4% | −74.9% | −72.2% |

Down at the clinical endpoint: **5 mg buys only 4 more hairs a year than
1 mg** (+105.2 vs +101.3 over placebo). Conversely, one-fifth the dose at
0.2 mg costs only 7 hairs (+94.2).

### 4. 5 mg only beats 1 mg in **scalp** DHT — and the reason can be derived

In this model finasteride is 2000-fold weaker against type 1 than type 2.
Free type 1 is 0.979 at 1 mg and 0.903 at 5 mg. Unlike blood, **scalp skin
carries the sebaceous gland's type 1 enzyme, so only the scalp compartment
can see that 7.6%.** The asymmetry reported in the literature (scalp
−64.1 → −69.4, while serum barely moves, −71.4 → −72.2) is reproduced
without any scalp-specific fitting parameter.

### 5. Only **recently** miniaturised follicles can be recovered, and that resource ages

The vellus grade is split into two — A4a/C4a/T4a (recently vellus,
re-enlargeable, weight 0.75) and A4b/C4b/T4b (fixed vellus, weight 0.05) —
with an **irreversible** one-way flux `p_age` between them. This single
structure simultaneously produces three things: the finasteride response
plateau, late-start failure, and the irreversibility of Norwood VI.

| Treatment start (hair count) | Apparent age | Vellus pool | Deficit recovered after 3 years |
|---|---|---|---|
| 1150 | 30.8 y | 77 | **58.8%** |
| 1000 | 39.9 y | 121 | 40.9% |
| **876** | 47.8 y | 130 | **30.1%** |

### 6. Hair **mass** moves at twice the rate of hair count — because mass scales as d²

Same simulation, 12 months, two endpoints:

| Treatment | Count change | Mass change | Ratio |
|---|---|---|---|
| Finasteride 1 mg | +9.9% | +18.6% | **1.88** |
| Minoxidil 5% (responders) | +2.9% | +5.9% | **2.05** |
| Dutasteride 0.5 mg | +12.3% | +20.3% | 1.65 |

**A trial that counts hairs is reading the smaller half of the effect.**

### 7. Minoxidil's initial shed ("dread shed") is not a side effect — it is **the mechanism audibly working**

Minoxidil sulfate shortens telogen (T_T from 100 to 49 d at maximal
effect). But exiting telogen requires first releasing the club hair. The
model contains no term that treats shedding as toxicity, and yet:

| | Baseline shed | Peak | Fold | Timing | Returns below placebo | 12 months |
|---|---|---|---|---|---|---|
| Responder (SULT 2.0) | 2.38/day | **3.88** | ×1.63 | **day 10** | day 140 | +25.1 |
| Non-responder (SULT 0.25) | 2.38 | 3.10 | ×1.31 | day 10 | day 210 | +3.9 |

Converted to the full scalp, that is 188 → 306 hairs/day. **Telling a
patient that "shedding is normal" is not reassurance — it is arithmetic.**

### 8. Because minoxidil is a prodrug, **a trial's average is a mixture of two populations**

Same 5% solution, 12 months, versus placebo: **responders +39.8 hairs /
non-responders +18.7 hairs.** A single value — follicular SULT1A1
activity — creates a 2.1-fold difference. A trial reporting an average is
reporting neither group.

### 9. Systemic exposure from topical minoxidil is **a quarter of oral 5 mg**

| | Absorbed | Plasma Cavg | Cmax | ΔMAP | Hypertrichosis index |
|---|---|---|---|---|---|
| Topical 5% BID | 1.13 mg/d | 1.54 ng/mL | 1.65 | −1.0 mmHg | 11.1 |
| Oral 1 mg | 0.90 | 1.11 | 4.15 | −0.7 | 8.8 |
| Oral 5 mg | 4.50 | 5.53 | 20.7 | **−3.6** | **43.9** |
| Oral 10 mg | 9.00 | 11.06 | 41.5 | −7.2 | 87.8 |

The hair effect of oral 5 mg is almost identical to topical 5%
(+36.7 vs +39.8), yet systemic exposure is 3.6-fold higher. This is why
hypertrichosis and oedema from low-dose oral minoxidil appear in a
dose-dependent stepwise pattern.

### 10. On discontinuation, DHT returns in **ten days**, while hair is lost over **years**

After 2 years of treatment, following discontinuation:

| Days after stopping | Scalp DHT | Free scalp 5AR-II | Hair count vs baseline |
|---|---|---|---|
| 0 | 0.364 | 0.036 | +109.9 |
| **10** | **0.980** | 0.967 | +110.2 |
| 45 | 1.000 | 1.000 | +109.8 |
| 1 year | 1.000 | 1.000 | +101.0 |
| 3 years | 1.000 | 1.000 | **+79.4** |

The pharmacological clock (enzyme resynthesis t½ 2 days) and the clinical
clock (follicle cycle, hundreds of days) differ by two orders of magnitude.
This is why a patient says "it was fine for a while after I stopped" — and
why that statement is dangerous.

### 11. The same arithmetic **forgives imperfect adherence**

Taking the pill only 4 days out of 7 raises scalp DHT to 0.389 instead of
0.364, and the 1-year hair count becomes **+91.4** instead of +101.3 over
placebo. **57% of the pills buys 90% of the effect,** because what is
being titrated is the enzyme, not the plasma.

### 12. Topical finasteride **separates** the two compartments

Giving the scalp enzyme pool and the systemic enzyme pool separate values
(the model's `E1S`/`E2S` versus `E1`/`E2`), 0.25% topical yields:

| | Scalp DHT | Serum DHT | Testosterone | PSA | Sexual side-effect index |
|---|---|---|---|---|---|
| Oral 1 mg | −63.6% | −70.3% | +13.1% | −42.2% | 0.617 |
| **Topical 0.25%** | **−47.6%** | **−27.4%** | +4.9% | −17.4% | **0.259** |
| Dutasteride 0.5 mg | −90.3% | −96.1% | +18.1% | −56.9% | 0.979 |

The observed serum DHT suppression from the topical formulation is
about −25% (Piraccini 2022, PMID 34634163). **78%** of the hair effect is
bought with **42%** of the systemic androgen signal.

### 13. Female-pattern hair loss fails to respond not because of the drug but because of **what is driving it**

The model has an androgen-independent inhibitory drive, `ARIND`. In the
postmenopausal female arm (follicular 5AR at 0.30×, aromatase at 1.6×),
changing only the ratio of `GS` to `ARIND` — **and leaving the drug model
entirely untouched** — while matching the same severity (apparent age
33.5 years, hair count ~899, vellus 216):

| Androgenic-drive fraction | Finasteride 1 mg | Spironolactone 100 mg | **Minoxidil 5%** |
|---|---|---|---|
| 34% | +53.2 | +136.6 | +67.2 |
| 22% | +35.2 | +99.9 | +67.1 |
| 11% | **+15.7** | +53.4 | **+67.1** |

**An antiandrogen's effect is proportional to the androgenic-drive
fraction, while minoxidil's effect is completely independent of it (equal
even to the decimal place).** Price 2000's negative postmenopausal
finasteride trial (PMID 11050579) reads not as the drug failing but as the
result of enrolling a low-androgen-fraction population — the same model
produces +101.3 in men.

---

## Clinical Trial Reproduction (Calibration)

| Target | Source | Observed | Model |
|---|---|---|---|
| Baseline hair count in target area | Kaufman 1998 [9777765] | 876 | 876 (burn-in stopping condition) |
| Finasteride 1 mg vs placebo, **1 year** | Kaufman 1998 | **+107** | **+101.3** |
| Finasteride 1 mg vs placebo, **2 years** | Kaufman 1998 | **+138** | **+139.2** |
| Scalp DHT, 1 mg, day 42 | Drake 1999 [10495374] | −64.1% | −63.6% |
| Serum DHT, 1 mg, day 42 | Drake 1999 | −71.4% | −70.3% |
| Serum DHT, dutasteride | Clark 2004 [15126539] | ~−94% | −96.1% |
| Dutasteride > finasteride (24 weeks) | Gubelin Harcha 2014 [24411083] | p=0.003 | +80.8 vs +63.4 |
| Testosterone rise | Product label | +9–15% | +13.1% |
| Normal daily shed count | Textbook | ~100 | 90 |
| Topical finasteride serum DHT | Piraccini 2022 [34634163] | ~−25% | −27.4% |

### Where the model disagrees with the literature (reported as-is, not tuned away)

- **Dutasteride's scalp DHT**: model −90.3%, literature −51 to −79%. This is
  because the scalp compartment's 5AR-independent floor value (`W_ALT` 6%)
  is small. Raising it fixes dutasteride but breaks the much better-measured
  finasteride 1 mg point. **Reported explicitly as an over-prediction.**
- **Drake's 0.05 mg scalp point** (−61.6%): reproduction fails (model
  −43.1%). The raw data are themselves non-monotonic — 0.05 mg is reported
  as suppressing more than 0.2 mg (−56.5%), which no mass-action model can
  produce. The **serum** point at the same dose (−49.5%) matches well at
  −47.3%.
- **Setipiprant**: at `PGD_EFF 0.18` the model predicts +46 hairs over
  placebo at 12 months, but the actual phase 2a trial failed. Rather than
  deleting the arm, it is used **in reverse** — to be consistent with a
  negative trial, the PGD2 axis must account for **5% or less** of the
  total inhibitory drive:

  | PGD_EFF | 0.18 | 0.10 | 0.05 | 0.02 |
  |---|---|---|---|---|
  | Setipiprant, 12 months (vs placebo) | +46.0 | +27.6 | +14.5 | +6.0 |

- **Placebo-arm rate of decline**: model −14.7 (1 year) / −29.3 (2 years),
  literature −21 / −55. The model's untreated arm is linear, while the real
  trial accelerates.

---

## Implementation Verification

Every equation was first written and run as an **independent Python/SciPy
implementation** (fixed-step RK4, 50 states) before being ported to
mrgsolve. This cross-implementation caught **four genuine defects**, all
four of which are already fixed in the current files:

1. **Missing grade-3-to-grade-2 upward flux.** Follicles vanished silently,
   breaking total conservation, and as a result finasteride came out
   **worse** than placebo (12 months −10.2). After the fix: +101.3.
2. **A 1000-fold concentration-unit error.** Volume of distribution was set
   in L while the binding constant was expressed per ng/mL, driving free
   type 2 enzyme to zero at every dose (scalp DHT −86%, actual −64%).
3. **DHT weighting of 0.65 in the testosterone feedback.** On finasteride,
   testosterone rose **+47%** (observed ~+10%). Changed to 0.30, giving
   +13.1%.
4. **An RK4 step of 0.05 days was unstable at 5 mg.** k_on × C_max =
   224/day exceeded the stability limit, turning the entire arm to NaN,
   which silently disappeared from the table. Fixed by reducing the step to
   0.005 days (not an issue for mrgsolve's LSODA).

---

## How to Run

```r
# Model
source("aga_mrgsolve_model.R")
P0  <- burn_in(mod)                       # age 47.8, hair count 876
fin <- run_arm(mod, P0, ev_fin(1))        # finasteride 1 mg, 5 years
mxt <- run_arm(mod, P0, ev_mxt(50), SULT = 2.0)

# Dashboard
shiny::runApp("aga_shiny_app.R")
```

```bash
# Re-render the mechanistic map
dot -Tsvg aga_qsp_model.dot -o aga_qsp_model.svg
dot -Tpng -Gdpi=150 aga_qsp_model.dot -o aga_qsp_model.png
```

---

*This model is a quantitative systems pharmacology model for research and
educational purposes; it is not a clinical guideline. Finasteride and
dutasteride are contraindicated in women of childbearing potential and are
teratogenic to male fetuses.*
