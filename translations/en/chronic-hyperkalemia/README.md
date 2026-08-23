# Chronic Hyperkalaemia — QSP Model

**Chronic hyperkalaemia in CKD and heart failure — the RAAS-inhibitor potassium dilemma**

<a href="../../../chronic-hyperkalemia/hk_qsp_model.svg"><img src="../../../chronic-hyperkalemia/hk_qsp_model.png" width="720" alt="chronic hyperkalaemia QSP map"></a>

---

## 0. The one idea of this model

**Serum potassium is not a reservoir, it is a ratio.**

```
K_total  =  Ce·V_ECF  +  Ci0 · LAMrel · (Ce/Ce0)^α · V_ICF
```

- `LAMrel` — the **partition** between the inside and the outside of the cell.
  Insulin, β2, pH and osmolality change it on a **timescale of minutes**.
- `K_total` — the **total** body potassium. Diet, kidney, colon, binders and dialysis
  change it on a **timescale of days**.
- `α = 0.25` — the exponent that makes the intracellular space a **buffer** rather
  than a *mirror*.

Every treatment for hyperkalaemia acts on **exactly one term** of this equation, and
which term it acts on — rather than the size of its effect — is what determines when
that drug is the right drug.

| Treatment | The term it acts on | Potassium actually removed |
|---|---|---|
| Insulin/glucose, salbutamol | `LAMrel` | **0 mmol** |
| Calcium gluconate | **neither term** (it only shifts the threshold potential) | **0 mmol** |
| Potassium binders, diuretics, dietary restriction, dialysis | `K_total` | real mmol |
| Alkali therapy | **renal secretion**, not a shift into the cell | real mmol |

---

## 1. Deliverables

| File | Contents |
|---|---|
| [`hk_qsp_model.dot`](../../../chronic-hyperkalemia/hk_qsp_model.dot) · [`.svg`](../../../chronic-hyperkalemia/hk_qsp_model.svg) · [`.png`](../../../chronic-hyperkalemia/hk_qsp_model.png) | Mechanistic map — **134 nodes / 15 clusters / 229 edges** |
| [`hk_mrgsolve_model.R`](../../../chronic-hyperkalemia/hk_mrgsolve_model.R) | mrgsolve ODE model — **30 compartments**, 12 scenario functions |
| [`hk_reference_model.py`](../../../chronic-hyperkalemia/hk_reference_model.py) | Pure-Python reference implementation (zero dependencies) — it **computes** every number below |
| [`hk_model_report.txt`](../../../chronic-hyperkalemia/hk_model_report.txt) | The full output of that script (346 lines) — calibration, validation, discrepancies |
| [`hk_shiny_app.R`](hk_shiny_app.R) | Shiny dashboard — **9 tabs** |
| [`hk_references.md`](hk_references.md) | **84 references**, every one cross-checked against PubMed |

```bash
python3 hk_reference_model.py      # regenerate the report (~2 min 30 s, no dependencies)
dot -Tsvg hk_qsp_model.dot -o hk_qsp_model.svg
R -e 'library(mrgsolve); mod <- mread("hk_mrgsolve_model.R"); scenario_all(mod)'
R -e 'shiny::runApp("hk_shiny_app.R")'
```

> **A note on honesty.** R was not installed in this container.
> `hk_mrgsolve_model.R` and `hk_shiny_app.R` **have not been verified by running
> them.** Every number below was obtained by actually running
> `hk_reference_model.py`, which holds the same equations, and the mrgsolve file is a
> 1:1 transcription of that Python implementation. That the two implementations agree
> is **the author's claim, not a verified fact**.

---

## 2. Calibration and validation

### 2.1 Only five parameters were fitted

| Parameter | Value | Fitting anchor |
|---|---|---|
| `S0` (ASDN secretory capacity) | 75.18 mmol/day | eGFR 100 → K 4.20 |
| `ADAPT_P` (per-nephron up-regulation exponent) | 0.879 | eGFR 20 → K 5.00 |
| `ADAPT_MX` (ceiling on up-regulation) | 6.36 | eGFR 12 → K 5.30 |
| `KI_MRA` (MR competitive inhibition constant) | 0.1029 mg/L | RALES ΔK **+0.30** |
| `N_HCO3` (suppression of distal secretion by acidosis) | 0.435 | ΔK **−0.30** for HCO3 18→24 |

### 2.2 Everything else is prediction (held-out)

| eGFR | Model | Cohort target | Error |
|---|---|---|---|
| 60 | 4.47 | 4.40 | +0.07 |
| 45 | 4.61 | 4.55 | +0.06 |
| 30 | 4.80 | 4.75 | +0.05 |

Mean absolute error **0.059 mmol/L**.

### 2.3 External validation — not one of these was used in the fitting

| Trial/observation | Observed | Model | Verdict |
|---|---|---|---|
| The classical deficit nomogram (serum K 3.0) | −200 to −400 mmol | **−302 mmol** | derived from α |
| OPAL-HK patiromer 16.8 g/day, 4 weeks | −1.01 mmol/L | −0.95 | agrees |
| HARMONIZE SZC 5/10/15 g, 28 days | 4.8 / 4.5 / 4.4 | 4.89 / 4.63 / 4.49 | agrees (consistently +0.1) |
| HARMONIZE SZC 10 g TID, 48 hours | about −1.1 mmol/L | −0.95 | agrees |
| Insulin 10 U + glucose | −0.6 to −1.0 (30–60 min) | **−0.85 (50 min)** | agrees |
| Salbutamol 20 mg nebulised | −0.6 to −1.0 | −0.74 (2–4 hours) | agrees |
| ACE inhibitor alone | +0.1 to +0.4 | +0.21 to +0.23 | agrees |
| **FIDELIO finerenone 20 mg** | **+0.23 mmol/L** | **+0.10** | **disagrees** |
| **Hypoglycaemia after insulin/glucose** | **15–20%** | **0%** | **not reproduced** |

---

## 3. What the model actually says

### 3.1 The potassium reserve — serum K is a lagging indicator

The **fractional excretion of potassium (FE_K)** rises from **10.3% to 51.2%** as
eGFR goes from 100 to 12. Over the same range serum K moves only from 4.20 to 5.30.

> The kidney does not fail gradually with respect to potassium. It **spends the whole
> of its reserve keeping serum K flat**, and only once that reserve is exhausted does
> serum K move. Serum potassium is a **lagging indicator** of a reserve that has been
> depleting for years.

### 3.2 MRAs do not raise K, they *raise the threshold*

The eGFR threshold at which steady-state K exceeds 5.5:

| Prescription | Threshold eGFR |
|---|---|
| No RAAS blockade | **11.3** |
| + ACE inhibitor | 12.0 |
| + ACEi and MRA | **20.6** |
| + ACEi, MRA, patiromer 16.8 g | **6.7** |

A patient at eGFR 21 is not "a bit high on potassium because of the
spironolactone" — that patient has been **converted into an eGFR of 11 in terms of
potassium-handling capacity**. A binder pushes the threshold back to 6.7, buying back
about **14 mL/min/1.73 of apparent potassium-handling capacity** without changing a
single nephron.

### 3.3 I was wrong and the model was right — the potassium price of an MRA is *flat*

When building this model I expected the potassium price of an MRA to grow as eGFR
fell. **It does not.** From eGFR 90 down to 15 the increments the model produces are
**+0.30, +0.30, +0.30, +0.29, +0.28** — essentially constant.

What grows is not the increment but **the baseline it is laid on top of**. That is an
entirely different clinical claim — what is dangerous is not that *the drug effect
steepens* but that *the distance to the threshold shortens*. This report adopts what
the model supports rather than what I expected beforehand.

### 3.4 Binders have a computable ceiling

A binder can only capture **the potassium that is in the intestinal lumen**, so the
maximum negative balance is structurally limited to
`φ_max × f_abs × dietary intake`.

| Agent | φ | Maximum removal on a diet of 80 mmol/day |
|---|---|---|
| SZC 5 g/day | 0.173 | 12.5 mmol/day |
| SZC 10 g/day | 0.250 | 18.0 mmol/day |
| Patiromer 8.4 g/day | 0.173 | 12.5 mmol/day |
| Patiromer 16.8 g/day | 0.245 | 17.6 mmol/day |
| Patiromer 25.2 g/day | 0.285 | 20.5 mmol/day |

> No binder at any dose can remove more than about **43 mmol/day** on a standard
> diet. For a patient in whom renal plus colonic excretion falls short of intake by
> more than that, **dialysis is not a preference but a requirement of mass
> conservation**.

### 3.5 The same number, the opposite reservoir

The exchange rate is derived from α — **224 mmol/(mmol/L) chronically, 66 acutely**.
(The acute figure of 66 explains directly why 40 mmol of intravenous KCl is dangerous
and eating 40 mmol is not.)

| | Patient A — the mass-balance type | Patient B — the partition type |
|---|---|---|
| Situation | eGFR 12, ACEi + spironolactone | diabetic ketoacidosis (pH 7.09, glucose 40) |
| Serum K | **6.80** | **4.68** — reassuringly normal |
| `LAMrel` | 0.999 (normal) | **0.861** |
| Total body potassium | **+482 mmol excess** | **−400 mmol deficit** |
| After `LAMrel` is normalised with insulin and alkali | — | **K 2.67 — and that is a cardiac arrest** |
| The right treatment | binder/dialysis | insulin + **potassium replacement** |

The laboratory number on its own carries no information. It becomes information
**only when the number and `LAMrel` are known together**.

### 3.6 The acute phase — read the upper picture and the lower one together

Each treatment given alone to patient A (eGFR 12, serum K 6.80, still eating):

| Intervention | 0.5h | 1h | 2h | 4h | 8h | 24h | 48h | K removed |
|---|---|---|---|---|---|---|---|---|
| No treatment | 6.74 | 6.73 | 6.71 | 6.70 | 6.69 | 6.63 | 6.55 | **0.0 mmol** |
| Calcium gluconate 1 g | 6.74 | 6.73 | 6.71 | 6.70 | 6.69 | 6.63 | 6.55 | **0.0 mmol** |
| Insulin 10 U + D50 25 g | 5.97 | 5.83 | 6.34 | 6.69 | 6.68 | 6.63 | **6.55** | **0.0 mmol** |
| Salbutamol 20 mg nebulised | 6.54 | 6.31 | 6.07 | 6.00 | 6.13 | 6.58 | 6.56 | **0.0 mmol** |
| SZC 10 g TID | 6.72 | 6.70 | 6.65 | 6.59 | 6.49 | 6.16 | **5.85** | **68.8 mmol** |

Insulin looks dramatically good, with a nadir of 5.82 at 50 minutes, but **at 48 hours
it is indistinguishable from no treatment at all (6.55 vs 6.55).** Because it removed
0 mmol.

**The membrane-potential arithmetic (why calcium comes first, and why it is absent
from the table above):**

| | Em | Threshold | Excitability margin | QRS |
|---|---|---|---|---|
| K 4.2 | −93.7 mV | −70.0 | **23.7 mV** | 90 ms |
| K 6.8 | −80.8 mV | −70.0 | **10.8 mV** (−54%) | 111 ms |
| + calcium 1 g | −80.8 mV (**unchanged**) | −66.2 | **14.5 mV** | 111 ms |

Calcium restores **29%** of the lost margin within a minute while removing
**0 mmol** of potassium. It buys time; it does not treat.

### 3.7 The RAASi–potassium dilemma — five years with the clinician in the loop

The clinician was modelled as a controller:
`dRD/dt = −k_down·max(0, K−5.5)·RD + k_up·max(0, 5.0−K)·(1−RD)`

A diabetic CKD 4 patient (eGFR 25 at entry, diet 80 mmol/day), over five years:

| Prescription | K (5 yr) | **RAASi dose** | eGFR (5 yr) | Days with K>5.5 | Mean hazard ratio |
|---|---|---|---|---|---|
| ACEi + MRA, no binder | 5.50 | **51%** | 14.1 | **1009** | 0.861 |
| + patiromer 16.8 g/day | 4.66 | **100%** | 15.0 | 0 | **0.639** |
| + SZC 10 g/day | 4.65 | **100%** | 15.0 | 0 | 0.638 |
| + low-potassium diet 50 mmol/day | 4.47 | **100%** | 15.0 | 0 | **0.638** |
| + furosemide 40 mg/day | 5.50 | 98% | 15.0 | 44 | 0.788 |
| + SGLT2 inhibitor | 5.48 | 100% | 15.0 | 0 | 0.777 |
| + oral alkali 70 mmol/day | 5.50 | 98% | 15.0 | 24 | 0.788 |
| ACEi alone, MRA never started | 5.48 | 100% | **12.0** | 0 | 0.891 |

> **Read the RAASi column, not the K column.** Every prescription arrives at an
> acceptable potassium in the end — that is what the controller is doing. The
> difference is **how much RAASi was lost as the price of getting there**. Most of the
> benefit in the binder arms comes not from preventing arrhythmia but from *being able
> to carry on with renoprotective therapy*.
>
> The lowest mean hazard is the **low-potassium diet (0.638)**, essentially identical
> to the binders. The model does not say that an expensive drug is needed — it says
> only **lower the potassium load**, and a binder is one of several ways of achieving
> that.

### 3.8 Alkali acts on the kidney and not on the cell (an unexpected result)

This result was not put in; it came out. **At steady state serum K is fixed by mass
balance, and the partition only sets the size of the reservoir behind it.**

| HCO3 | pH | LAMrel | Urinary K | Steady-state K |
|---|---|---|---|---|
| 16 | 7.344 | 0.994 | 48.1 | 5.44 |
| 20 | 7.377 | 0.999 | 49.3 | 5.20 |
| 24 | 7.401 | 1.003 | 50.2 | 5.03 |

Correcting HCO3 from 18 to 24 is worth **−0.28 mmol/L**, and **100%** of that comes
not from a shift into the cell but from **additional kaliuresis (+1.5 mmol/day)**.
The well-known negative result follows from that — **bicarbonate is not a treatment
for acute hyperkalaemia.** Because the mechanism it actually has needs a kidney, and
needs weeks.

### 3.9 A virtual population of 600 — who is it that a binder rescues?

A grid of eGFR 15–60 × dietary K 40–140 mmol/day × HCO3 16–26, all of them on
ACEi+MRA:

- without a binder, **315 (52.5%)** have K > 5.5
- adding patiromer 16.8 g/day returns **194 of those (61.6%)** to below 5.5
- **121 (20.2%)** remain above 5.5 despite the binder

| eGFR band | n | K>5.5 without a binder | With a binder |
|---|---|---|---|
| 15–29 | 180 | 57.2% | 26.7% |
| 30–44 | 180 | 53.3% | 20.6% |
| 45–59 | 180 | 48.9% | 15.6% |

The benefit of a binder is not uniform. It is concentrated in the patients who
**exceed the threshold by less than the binder's ceiling**, and that band is a
**computable quantity** rather than a clinical impression.

---

## 4. Model structure

### 4.1 The mechanistic map — 15 clusters

| # | Cluster | Key content |
|---|---|---|
| 1 | Diet · gastrointestinal tract | gut potassium sensor, absorption, colonic BK secretion, faeces |
| 2 | Whole-body distribution | ECF · fast ICF · slow ICF (muscle), α, buffer capacity |
| 3 | Movement across the cell membrane | insulin-PI3K-Na/K-ATPase, β2-cAMP, pH, osmolality, cell lysis |
| 4 | Glomerulus · proximal | filtered load, PCT, TAL, distal delivery, distal flow |
| 5 | **ASDN — the potassium reserve** | WNK-SPAK, NCC, ENaC, ROMK, BK, H-K-ATPase, S_cap |
| 6 | RAAS · MR | renin-AngII-aldosterone, direct stimulation by K, MR, SGK1 |
| 7 | Acid-base | NEAP, NAE, bone buffering, Winter's correction, alkali therapy |
| 8 | RAAS inhibitors | ACEi/ARB/ARNI/DRI, steroidal and non-steroidal MRAs, **the down-titration loop** |
| 9 | Potassium binders | SPS, patiromer (distal colon), SZC (stomach · proximal), the φ ceiling, RAASi enablement |
| 10 | Acute measures | calcium, insulin/glucose, salbutamol, loop diuretics, dialysis, **rebound** |
| 11 | Concomitant drugs | NSAIDs, trimethoprim, heparin, CNIs, β-blockers, SGLT2i, salt substitutes |
| 12 | Cardiac electrophysiology | Em, threshold, Na channel availability, QRS, IKr, T wave, sine-wave |
| 13 | Precipitants | CKD, heart failure, diabetes, AKI on CKD, DKA, pseudohyperkalaemia |
| 14 | Clinical endpoints | serum K, arrhythmia, hospitalisation, **RAASi discontinuation**, CKD progression, death |
| 15 | Observables | FE_K, TTKG, urinary K, aldosterone, threshold eGFR |

### 4.2 The mrgsolve compartments (30)

`KE` `KIF` `KIS` `HCO3A` `ALDO` `RASDN` `GUTK` `COLK` `AACE` `CACE` `AMRA`
`CMRA` `PATP` `PATC` `SZCP` `SZCC` `INS` `INSE` `GLU` `B2C` `B2E` `AFUR`
`CFUR` `CAE` `GFR` `RD` `BICG` `T55` `CHAZ` `KREM`

That `RD` (the prescribed RAASi dose) is a **state variable** is the structural
distinguishing feature of this model. The clinician's down-titration decision enters
not as an exogenous input but as part of a feedback loop coupled to potassium
physiology, and the results of §3.7 *came out of* that coupling rather than being
imposed on it.

### 4.3 The scenario functions (12)

`scenario_reserve()` `scenario_binder()` `scenario_acute()` `scenario_dilemma()`
`scenario_threshold()` `scenario_diet()` `scenario_alkali()`
`scenario_two_patients()` `scenario_population()` `scenario_ecg()`
`sim_steady()` `scenario_all()`

---

## 5. Where not to trust this model

1. **The outcome layer (hazard) is association, not causation.** `hazard_K()` was
   fitted to the U-shaped curve of the observational studies, and in that curve a high
   potassium is in part *a marker of the disease that caused it*. So the model almost
   certainly **overestimates the benefit of "lowering the number" and underestimates
   the benefit of "treating the cause of the rise"**. Every hazard ratio in §3.7 has
   to be read with that bias in mind. This is the strongest assumption in the whole
   model and also the most weakly grounded.
   (The STOP-ACEi trial is included as a counterexample at reference #65.)
2. **The FIDELIO discrepancy.** Assuming finerenone's MR loading to be 35% of
   spironolactone 25 mg gives +0.10, short of the observed +0.23. Working backwards,
   **72%** is required. The model cannot distinguish "the occupancy is higher" from
   "the occupancy is the same but the tissue distribution differs", and the latter is
   precisely what is claimed for the non-steroidal MRAs. So this is reported as **an
   unresolved non-identifiability rather than as a fitted parameter**.
3. **It fails to reproduce hypoglycaemia.** The insulin effect site decays together
   with the plasma concentration (real tissue action lasts several hours longer), and
   there is no between-patient variability in glycogen reserve. Even giving no glucose
   at all, the nadir stays above 4 mmol/L.
4. **The clinician controller is a caricature.** Real down-titration is discrete, is
   delayed by the interval between blood tests, and is often permanent. A model in
   which the dose climbs back up of its own accord is optimistic about exactly the
   behaviour the binder trials were trying to change.
5. **The intracellular space is two compartments and is still not muscle.** It is
   sufficient to reproduce the chronic (224) and acute (66) buffer capacities at the
   same time, but it cannot represent **exercise**, which raises venous K by more than
   1 mmol/L within a minute through accumulation in the interstitium of the working
   muscle, and equally cannot represent **release from the reservoir itself** in
   rhabdomyolysis and tumour lysis.

### What is not in the model at all

**AKI on CKD** (the commonest cause of real severe hyperkalaemia) · digoxin
toxicity · the ASDN actions of trimethoprim/heparin/calcineurin inhibitors (on the
map but not in the ODEs) · tumour lysis syndrome · rhabdomyolysis ·
**pseudohyperkalaemia** · **dialysis** (an intermittent, very-high-clearance removal
term) · paediatrics and pregnancy (the parameters are for a 70 kg adult).

---

## 6. References

[`hk_references.md`](hk_references.md) — **84 references**, in 13 sections.
Every PMID was cross-checked against the PubMed esummary API, and the **19 that were
inaccurate in the first draft were replaced with the actual papers**. Each entry is
marked with the role that reference played in the model: **[F]** fitting /
**[V]** validation / **[S]** structural / **[C]** clinical context. Appendix A holds
a table of sources per parameter and Appendix B a list of out-of-scope entries.

---

## ⚠️ Disclaimer

This is a **quantitative systems pharmacology model for educational and research
purposes**. It was built from the public literature and clinical-trial data but has
not been independently validated or certified, and **must not be used directly for
real clinical decision-making, prescribing, or regulatory submission.**
In particular, the hazard ratios of §3.7 must not be quoted without reading the
limitations listed in §5 — above all the fact that the outcome layer is not causal.
