# Adrenocortical Carcinoma (ACC) — Quantitative Systems Pharmacology Model

<p align="center">
  <a href="../../../adrenocortical-carcinoma/acc_qsp_model.svg">
    <img src="../../../adrenocortical-carcinoma/acc_qsp_model.png" width="880" alt="ACC QSP mechanistic map">
  </a><br>
  <sub>Click to open the vector (SVG) original · 258 nodes · 20 clusters</sub>
</p>

---

## The organising thesis: plasma mitotane is not the dose but **the overflow of a reservoir**

The apparent volume of distribution of mitotane (o,p'-DDD, logP ≈ 5.9) is about **6,000 L** — some
75 times total body water, and most of it **fat**. For a constant intake the plasma trough follows

```
        Cp(t)  =  (D · F / CL) · [ 1 − exp(−t · CL / Vss) ]
                   ╰─────────╯    ╰──────────────────────╯
                    where it gets to          when it gets there
                   (WHERE: dose/clearance)   (WHEN: τ = Vss/CL ≈ 50–110 days)
```

And **efficacy is a threshold on Cp** (≥ 14 mg/L), while **neurotoxicity is a second threshold on the
same Cp** (> 20 mg/L).

Very nearly every clinical peculiarity this drug and the pharmacotherapy of this disease display
follows from one fact: **a single state variable — mitotane sitting in fat — is read in three ways
that pull in opposite directions.**

| | Reading | Mechanism | Time scale |
|---|---|---|---|
| **①** | **Benefit** | SOAT1 inhibition → accumulation of free cholesterol and oxysterols → ER stress → apoptosis (**adrenolysis**) | slow · **irreversible** |
| | | direct inhibition of CYP11A1 / CYP11B1 / 3β-HSD (**synthesis blockade**) | fast · **reversible** |
| **②** | **Harm #1** | PXR → ~4-fold induction of CYP3A4 → etoposide exposure halved, hydrocortisone clearance doubled, **and its own clearance raised as well** | intermediate (t½ ≈ 4–5 days) |
| **③** | **Harm #2** | the same molecule gives cerebellar and cognitive toxicity above 20 mg/L | recovery goes at the rate the reservoir empties |

The therapeutic window is a **slot 6 mg/L wide** sitting on a variable with **between-subject
variability of CV 81.5%**.

---

## What the model demonstrates by computation (8 results)

Every number is output from `acc_reference_run.R` and is stored verbatim in
`acc_reference_output.txt`.

### 0. The reservoir is real

| | |
|---|---|
| Apparent Vss | **5,966 L** (reported ~6,086 L) |
| The fat-dependent part of it | 4,466 L (**75%**) |
| Terminal half-life (before induction → fully induced) | **86 days → 42 days** (reported range 18–160 days) |
| The reservoir's share of total body content at day 180 | **93.6%** |

The plasma concentration is a thin overflow film floating on a lipid reservoir. Which is why it
responds not only to the dose but to **the size of the reservoir**.

### A. Fat determines not *where* but *when*

The same 6/4/3 g regimen, differing only in fat mass:

| Fat mass | Vss | Time to 14 mg/L | Peak concentration | **d360** | **d480** | Days above 20 mg/L |
|---|---|---|---|---|---|---|
| 10 kg | 3,530 L | **20 days** | 21.6 | 15.8 | 15.8 | **51 days** |
| 22 kg | 5,966 L | **34 days** | 18.1 | 15.8 | 15.8 | 0 |
| 45 kg | 10,635 L | **110 days** | 15.7 | 15.6 | 15.7 | 0 |

Change the fat mass 4.5-fold and **the destination is the same to within 0.3 mg/L**; only the time to
arrive spreads, by **5.5-fold**. And the patient who goes above 20 mg/L into the neurotoxicity risk
zone is the **lean** one. Checking against the closed form, Css is **identical** at 16.7 mg/L in all
three cases and only τ moves, 35.7 → 60.4 → 107.6 days.

> Clinical translation: the same prescription is **a drug needing dose reduction** in a lean patient
> and **a drug that has to be pushed for longer** in an obese one. The dose is the same and the clock
> is different.

### B. What the dose-initiation RCT was up against

| | Target reached (overall) | **Reached by day 90** | Median time to target | Median Cp90 |
|---|---|---|---|---|
| High-dose start | 81% | **76%** | 35 days | **16.0 mg/L** |
| Low-dose ramp | 62% | **15%** | 121 days | **10.5 mg/L** |

Variance decomposition of **the day-90 plasma concentration** (a continuous variable with no
censoring) in a virtual population (n = 240/arm):

| Term | Variance explained |
|---|---|
| Start regimen | **37.7%** |
| Fat mass | 24.7% |
| Clearance | 26.9% |
| Residual | 10.7% |

Two readings are both true and must not be confused.

- **Variance share:** the regimen is the largest single term (37.7%), but **the two unmeasured
  patient factors (fat 24.7% + clearance 26.9% = 51.6%) explain more than the dose does**.
- **Contrast size:** compared on the same scale, switching regimen moves the median Cp90 by
  **5.5 mg/L**, while moving fat mass from the 10th to the 90th percentile moves it by **9.1 mg/L**.
  The **larger contrast** belongs to the variable nobody randomised and nobody measured
  (the correlation between fat mass and Cp90 is r = −0.59).

On top of this comes the censoring problem. "Time to target", the endpoint that gets reported,
removes from the sample precisely the patients the covariate hurts most — **the high-fat patients who
never reach target at all.** In the one analysis that would reveal the covariate, the covariate
becomes invisible.

> This explains without contradiction the situation in which **a popPK model of the same drug
> recommends a high-dose start while a small RCT making the same comparison reports no significant
> difference.** Not because the dose does not matter, but because **the dose contrast that was tested
> is smaller than the between-patient contrast left in the residual.** The actionable conclusion:
> stratify by **body composition** and not by dose alone, and read the trough as **continuous
> exposure** rather than as pass/fail against 14 mg/L.

### C. Hypercortisolism digs the pit it will sit in beforehand — **before treatment even starts**

An untreated secretory tumour:

| Day | Free cortisol | Fat mass | Vss |
|---|---|---|---|
| 0 | 0.23 | 22.0 kg | 5,966 L |
| 90 | 2.28 | 25.2 kg | 6,623 L |
| 270 | 2.63 | **26.9 kg** | **6,959 L** |

A patient who spent roughly nine months in a Cushingoid state before diagnosis arrives with
**26.9 kg** of fat instead of 22.0 kg — that is, **a reservoir about 1,000 L larger**. On an identical
regimen, the time to 14 mg/L is pushed from **35 days to 41 days** and the day-90 concentration from
17.6 to 16.6 mg/L.

To be honest about it, this feedback **barely operates during treatment** (a difference of 0.1 kg in
fat mass at day 180). Because once mitotane brings cortisol under control, fat accumulation stops.
The whole of the effect is **concentrated in the period before diagnosis**, and that period is
precisely when nobody measures any of this. The disease has already lengthened its own treatment
delay before the first tablet is swallowed.

### D. EDP-M is weakest exactly where mitotane works best

The counterfactual experiment (`INDVIC = 0`) **leaves mitotane auto-induction in place** and removes
only the victim-drug interaction, so it isolates a single mechanism with mitotane exposure unchanged.

| | Etoposide CL fold | AUC per cycle | Change at nadir | d180 | Call |
|---|---|---|---|---|---|
| EDP-M, induction **ON** (reality) | **2.07** | **4.58** | −38.2% | −18.9% | PR |
| EDP-M, induction OFF (counterfactual) | 1.00 | **8.86** | −51.2% | −16.8% | PR |
| EDP alone (no mitotane) | 1.00 | 8.86 | −62.3% | −34.6% | PR |
| Sz-M (FIRM-ACT arm B) | — | — | −4.7% | +52.4% | PD |

Induction takes **half of the etoposide exposure (1.93-fold)** and a substantial part of the depth of
response with it. The FIRM-ACT ordering in which EDP-M beats Sz-M is reproduced, but the two arms of
the combination **cancel one another through a single shared state variable**. Because induction
tracks the reservoir rather than the dose, **the moment at which etoposide exposure is lowest is the
moment mitotane has finally reached therapeutic concentrations**.

> ⚠️ A caution on RECIST interpretation. Progression is judged from the **nadir**, so in a patient who
> responded more deeply the absolute reference for +20%/+5 mm is lower and nadir-based PFS can
> actually look shorter. That is why the table above also carries the time at which the baseline
> diameter is regained.

### E. Total cortisol lies — two errors of different kinds, in the same direction

A patient with a destroyed adrenal, mean over days 300–360:

| | CBG | Total | **Free** | vs normal | **Total/free ratio** |
|---|---|---|---|---|---|
| HC 20 mg/day, no mitotane | 1.00 | 8.47 | 0.501 | 100% | **16.9** |
| HC 20 mg/day + mitotane | 1.81 | 6.53 | **0.175** | **35%** | **37.2** |
| HC 50 mg/day + mitotane | 1.81 | 11.64 | 0.422 | 84% | 27.6 |

Decomposed, the two errors are **of different kinds**.

- **Induction = an error inside the patient.** It doubles hydrocortisone clearance and drags free
  cortisol down to **35%** of normal (34% with induction alone switched on). This is a real
  deficiency and the one that does physiological harm.
- **The rise in CBG = an error of the assay.** At steady state it barely changes the free
  concentration itself (83% with CBG alone switched on), but it raises **the total/free ratio from
  16.9 to 37.2**. **The total cortisol the clinician actually prescribes and reads off is reassuring**
  while free cortisol has fallen by half.

One error masks the other, and both push the clinician wrong in the same direction. Doubling the dose
to 50 mg/day recovers free cortisol to 84%. Titration must be on **free or salivary cortisol, not
total.**

### F. Cortisol is an immunosuppressant the tumour administers to itself

| | Free cortisol | Effector T cells | Nadir | d180 | Call |
|---|---|---|---|---|---|
| Pembrolizumab, **non-secretory** | 0.23 | **3.20** | −43.8% | −28.0% | **PR** |
| Pembrolizumab, **secretory** | 2.45 | 1.28 | 0.0% | **+93.6%** | PD |
| Pembrolizumab + mitotane | 0.15 | 2.58 | −8.4% | −0.3% | SD |
| No treatment, secretory | 2.59 | 0.39 | 0.0% | +156.0% | PD |

**PD-1 receptor occupancy is above 99% in all four cases.** Target engagement is never the problem.
The effector cells the antibody would release **are already suppressed by the tumour's own steroid
secretion**. Secretory ACC has, in effect, given itself a steroid pre-medication before it ever meets
a checkpoint inhibitor. Which is why **controlling the steroid is not an alternative to immunotherapy
but a precondition for it**.

### G. A target present in 90% of tumours, and yet a drug that failed (linsitinib)

| | Linsitinib | No drug |
|---|---|---|
| Steady-state concentration | 0.84 mg/L | — |
| Unblocked fraction of IGF1R | **0.23** (77% blocked) | 1.00 |
| Unblocked metabolic IR-B | 0.29 | 1.00 |
| **Unblocked tumour IR-A** | **0.64** (the escape route) | 1.00 |
| Blood glucose | **108.5** mg/dL | 92.0 |
| Insulin | **2.38** | 1.00 |
| **IGF/AKT/mTOR signal** | **0.93** | 1.39 |
| Tumour at day 300 | 1,080 mL | 1,460 mL |

The drug blocks 77% of IGF1R (target engagement achieved). But blocking the metabolic insulin
receptor along with it means **blood glucose rises → insulin rises 2.4-fold → and that insulin drives
the same PI3K–AKT node again through the comparatively less blocked tumour IR-A**. The target was
hit, the pathway was restored, and the trial failed. In the counterfactual with the escape route
closed (`IC50RA` set to the IGF1R level): signal 0.93 → **0.44**, tumour 1,080 → **695 mL**. That is,
**more than half of the potential effect is lost to insulin rescue**. The reason for the other half is
that ACC growth is not IGF-limited to begin with (`APROLIF` = 0.30).

### H. There is no washout — the long tail

Mitotane stopped at day 240:

| After stopping | Plasma | CYP3A4 | Etoposide CL fold | Adrenal mass |
|---|---|---|---|---|
| day 0 | 15.9 | 4.05 | 2.07 | 0.167 |
| day 60 | 6.6 | 3.14 | 1.75 | 0.260 |
| **day 120** | 3.3 | 2.38 | **1.48** | 0.414 |
| day 240 | 1.0 | 1.53 | 1.19 | 0.646 |

- Plasma below 2 mg/L: **167 days**
- CYP3A4 back to within ±10% of baseline: **453 days**
- Adrenal cortex recovering 50% of its mass: **158 days**
- **Etoposide clearance still 48% higher 120 days after stopping**

The reservoir goes on dosing the patient for months. **There is no such thing as washing mitotane out
before giving a CYP3A4 substrate.**

### I. Safety readout (reference fat mass, high-dose start, 8 cycles of EDP-M)

| Item | Value |
|---|---|
| Peak plasma mitotane | 19.6 mg/L |
| Days above 20 mg/L | 0 |
| Within the 14–20 mg/L window | **370 days / 400 days (93%)** |
| Peak CNS injury score | 0.104 |
| Peak ALT | 56 U/L |
| Free T4 nadir | 0.88 ng/dL (baseline 1.20) |
| Neutrophil nadir | 1.08 ×10⁹/L |
| Cumulative anthracycline | 320 mg/m² |
| eGFR (day 400) | 81 (baseline 92) |
| Residual adrenal cortex | **4.6%** |
| CBG / SHBG rise | 1.81 / 1.75-fold |

---

## Deliverables

| File | Contents |
|---|---|
| [`acc_qsp_model.dot`](../../../adrenocortical-carcinoma/acc_qsp_model.dot) · [`.svg`](../../../adrenocortical-carcinoma/acc_qsp_model.svg) · [`.png`](../../../adrenocortical-carcinoma/acc_qsp_model.png) | Mechanistic map — **258 nodes · 20 clusters** |
| [`acc_mrgsolve_model.R`](acc_mrgsolve_model.R) | **52 ODE compartments · 148 parameters · 18 scenarios** + virtual population and closed-form verification |
| [`acc_reference_run.R`](../../../adrenocortical-carcinoma/acc_reference_run.R) | The script that regenerates every number above |
| [`acc_reference_output.txt`](../../../adrenocortical-carcinoma/acc_reference_output.txt) | The actual run output (the source of the tables above) |
| [`acc_shiny_app.R`](../../../adrenocortical-carcinoma/acc_shiny_app.R) | **11-tab** interactive dashboard |
| [`acc_references.md`](acc_references.md) | **46 verified PubMed citations** + the complete list of estimated parameters |

### Running it

```bash
# render the map
dot -Tsvg acc_qsp_model.dot -o acc_qsp_model.svg
dot -Tpng -Gdpi=150 acc_qsp_model.dot -o acc_qsp_model.png

# regenerate all results (about 3-6 min)
Rscript acc_reference_run.R > acc_reference_output.txt

# dashboard
Rscript -e 'shiny::runApp("acc_shiny_app.R")'
```

Dependencies: `mrgsolve`, `dplyr` (model) / `shiny`, `ggplot2`, `tidyr`, `gridExtra`
(dashboard) / `graphviz` (map).

---

## Model structure (52 ODE compartments)

| Module | Compartments |
|---|---|
| Mitotane PK | `MITG` `MITC` **`MITA` (fat reservoir)** `MITE` `MITN` |
| Induction | `ENZ` (CYP3A4, gated through `ENZV` for victim drugs) |
| Binding proteins | `CBG` `SHBG` |
| Cytotoxic PK | `ETOC` `ETOP` `DOXC` `DOXP` `DOXCUM` `CISC` `CISP` `SZC` |
| Other drugs | `HCG` `PEMC` `LING` `LINC` |
| Tumour | `TUMS` `TUMR` |
| Steroidogenesis | `STCAP` `ADRN` `CORT` `ACTH` `GRO` `DHEAS` `OHP17` `ALDO` |
| Body composition · target organs | **`FATKG`** `MUSC` `BMD` |
| IGF · metabolism | `IGF2` `IGFSIG` `GLU` `INS` |
| Immunity | `TEFF` `TREG` |
| Myelosuppression (Friberg) | `PROLN` `TR1N` `TR2N` `TR3N` `CIRCN` |
| Toxicity · organs | `NTOX` `ALT` `FT4` `LVEF` `GFR` `AIHAZ` |
| Exposure accounting | `AUCMIT` `TIW` |

### The 18 scenarios

01–04 body composition · start regimen (lean/reference/obese × high dose/low dose) ·
05–06 secretory, monotherapy vs natural history ·
07–10 EDP-M / induction counterfactual / EDP alone / Sz-M ·
11–13 hydrocortisone 20 vs 50 mg vs a mitotane-free control ·
14–16 pembrolizumab secretory / non-secretory / with mitotane ·
17 linsitinib · 18 washout after stopping

---

## Verification

The following were checked before trusting the model, and several genuine bugs were found and fixed
in the process.

- **Is the healthy baseline state genuinely a fixed point.** Running 400 days with no drug confirmed
  that total cortisol 4.99 µg/dL, free 0.225, ACTH 1.000 and adrenal mass 1.000 **do not drift.** The
  initial ACTH feedback equation drifted to ACTH = 0.5 at normal cortisol and so was not a fixed
  point; it was fixed by normalising it. Glucose and insulin were handled the same way, with `GUPT`
  and `KINS` back-calculated so that 92 mg/dL / INS 1 is a fixed point.
- **Comparison against the closed form.** The AUC of a single etoposide dose agrees with the analytic
  `Dose/CL` (**2.763 vs 2.761**). The Css and τ of result A agree with the algebraic computation in
  `analytic_depot()`.
- **Boundedness of the feedback.** Setting tumour → cortisol and tumour → IGF2 linear in burden at
  first let free cortisol run away to 119 µg/dL. Both paths were made physiologically saturating
  (`TFRAC`).
- **Baseline value of the immune term.** In the initial setup, drug-free baseline immune killing
  (5.5%/day) exceeded tumour growth and **the untreated tumour was eradicated**. The
  immunosuppression term was normalised to equal 1 at normal cortisol and `SLPIMM` was recalibrated.
- **Distinguishing a numerical artefact from a real effect.** The etoposide AUC coming out 2.2 times
  the expected value was not the model but an error in **my own trapezoidal integration of a sharp
  peak on a 0.25-day grid** (solved with a 0.02-day grid). Conversely, the low trough obtained when
  the central compartment was set to 45 L was **a real structural problem** (the within-day amplitude
  too large), so `V1MIT` was changed to 400 L.
- **Cleanliness of the counterfactual.** Switching off mitotane auto-induction as well, with a single
  induction switch, changed the counterfactual arm's mitotane concentration from 15.8 to 28.2 mg/L and
  contaminated the comparison. A victim-drug-only switch (`INDVIC`) was separated out so that the
  comparison is made at **identical** mitotane exposure.
- **Does the Shiny app actually render.** All **32 outputs** were run headless through
  `shiny::testServer` and confirmed to pass (not merely parsing).

### Comparison with the literature

| Target | Literature | Model |
|---|---|---|
| Apparent Vss | ~6,086 L (CV 81.5%) | **5,966 L** |
| Terminal half-life | 18–160 days | **42–86 days** |
| Therapeutic window | 14–20 mg/L | implemented as a readout |
| Increase in etoposide CL | ~2-fold (estimated) | **2.07-fold** |
| CBG / SHBG rise | ~2-fold each | 1.81 / 1.75-fold |
| Hydrocortisone requirement | roughly 2-fold | 20 → 50 mg/day |
| Free cortisol fraction | ~4–5% | **4.5%** |
| Loss of adrenal cortex | effectively complete | 4.6% residual |
| FIRM-ACT ordering | EDP-M > Sz-M | reproduced |
| Pembrolizumab ORR | ~14–23% | subtype dependent (PR only in the non-secretory case) |

---

## Limitations — honestly

- **The etoposide induction coefficient is an estimate.** `FM3A4ETO = 0.35` is not a value measured in
  a mitotane–etoposide interaction study but was derived from the magnitude of CYP3A4 induction and
  the CYP3A4-metabolised fraction of etoposide. The **magnitude** of result D depends directly on
  this coefficient (the mechanism itself is established), and it is the first parameter that should be
  replaced by a measurement.
- **`KPFAT = 203 L/kg` is not a measured tissue partition value** but was back-calculated to reproduce
  the reported apparent Vss and the BMI covariate effect. The mechanism of result A (fat changes only
  τ) is structurally robust, but the exact fold changes depend on this value.
- **There is no circadian rhythm.** Cortisol is handled as a daily mean, and hydrocortisone is given
  three times daily so that the non-linearity of saturable binding actually operates. Instantaneous
  values in the scenario summary tables may be inter-dose troughs, so they **must be read as daily
  means**.
- **Tumour growth is a single Gompertz** and does not carry the real heterogeneity of ACC (doubling
  times of weeks to years). Response calls come from a single deterministic trajectory, so population
  statistics such as ORR should be read only as reproducing the direction of the subtype effect.
- **Linsitinib's IR-B < IR-A sensitivity difference is a mechanistically motivated assumption**
  (from the observation that hyperglycaemia arrives before the antitumour effect). It sets the
  magnitude of the rescue in result G.
- Radiotherapy, the extent of surgical resection (R0/R1/R2), locoregional therapy, and later lines
  such as cabozantinib and temozolomide are only partly on the map and are not in the ODEs.

**The complete list of estimated parameters and the nature of the evidence behind each** is tabulated
in [`acc_references.md` §9](acc_references.md).

---

## ⚠️ Disclaimer

This is a **QSP model for education and research**. It was assembled from the public literature and
clinical trial data but has not been independently validated or certified, and
**must not be used directly for real clinical decision-making, prescribing, or regulatory
submission.** The parameters and assumptions are illustrative approximations, and fitting and
validation against real patient data would be required separately.
