# Rhabdomyolysis-Induced Acute Kidney Injury (AKI) — QSP Model

> **One-sentence summary** — tubular toxic exposure is a **product**:
> *filtered myoglobin load ÷ urine flow × f<sub>diss</sub>(urine pH)*, and each of the three
> standard treatments moves only **one different term** in that product. Because it is a
> product, the more successful one term becomes, the smaller the marginal gain from the others
> — and this single piece of arithmetic explains the long-standing discrepancy that "the
> mechanism for alkalinisation is solid, but the trials are negative."

<p align="center">
  <a href="rhab_qsp_model.svg"><img src="rhab_qsp_model.png" width="820" alt="Rhabdomyolysis QSP mechanistic map"></a><br>
  <sub>158 nodes · 18 subgraph clusters — click the <a href="rhab_qsp_model.svg">SVG</a> to view full size</sub>
</p>

---

## 1. What This Model Argues (Three structural theses)

### THESIS 1 — Toxicity is the **product** of three terms, and treatment touches only one term at a time

```
     TOX  ~  (MB_plasma × GFR × σ)      ← extracorporeal therapy (numerator)
              ÷ Q_urine                 ← crystalloid          (denominator)
              × f_diss(pH_urine)        ← alkalinisation        (chemistry)
```

The **absolute amount** by which the same bicarbonate dose (12.5 mmol/h) reduces cumulative
toxic exposure:

| Urine flow | No alkali | +12.5 mmol/h | +30 mmol/h | Absolute gain from 12.5 |
|---|---|---|---|---|
| 100 mL/h | 100.1 | 36.3 | 40.8 ⚠ | **−63.8** |
| 250 mL/h | 63.9 | 30.4 | 31.5 | −33.5 |
| 500 mL/h | 55.1 | 29.2 | 25.9 | −25.9 |
| 1000 mL/h | 53.6 | 32.0 | 26.5 | **−21.6** |

The same table read as peak creatinine (mg/dL) goes from 3.64 → 1.89 → 1.65 (100 mL/h) to
2.07 → 1.65 → 1.52 (1000 mL/h), **the absolute gain from 12.5 mmol/h shrinking from 1.75 to
0.42 mg/dL.** Scenario 8 shows the same story from the other direction — **100 mL/h combined
with alkali (creatinine 1.88) lands almost exactly where 500 mL/h without alkali does (2.15).**
This model's answer to why Brown 2004 and Homsi 1997 failed to detect a bicarbonate benefit is
that those trials were run on top of a flow term that was already near its ceiling.

⚠ Note the cell at the top right of the table: raising the dose from 12.5 to 30 mmol/h at
100 mL/h **increases** cumulative exposure (36.3 → 40.8). Urine pH reaches 7.17, and because the
urine is not dilute enough to keep calcium phosphate in solution, crystallisation begins. This
carries forward into the pH-optimum discussion in (4) below.

### THESIS 2 — **A single pressure variable** opens **two opposing gates**

```
     ΔP = MAP − P_compartment
     gate 1 (ischaemia, opens when ΔP < 30)  →  produces necrosis
     gate 2 = 1 − gate 1                     →  produces washout
```

While a limb is trapped, gate 1 is open and gate 2 is shut — the muscle dies, but **quietly**.
Extrication flips both gates at once. So the surge in K⁺ and myoglobin immediately after
extrication (the "smiling death") is an **output of the model, not an assumption fed into it**.
The same arithmetic produces the fasciotomy trade-off (scenarios 14 vs 15, 10 kg, a tight
compartment):

| | Necrotic muscle | Salvaged muscle | Peak CK | Peak creatinine | Filtered Mb | eGFR d90 |
|---|---|---|---|---|---|---|
| Fasciotomy at 10 h | 6.19 kg | **3.81 kg** | 105,782 | 2.72 | 13.5 g | 96.9 |
| No fasciotomy | 10.00 kg | 0 kg | 41,257 | **1.18** | 3.9 g | 120.4 |

A limb that is never reperfused does not poison the kidney — because gate 2 stays shut. The
model reproduces the real, uncomfortable tension between limb salvage and renal load.
**Caution: the model has no wound infection, sepsis, or amputation, so the "no fasciotomy" arm
must not be read as a clinical recommendation.**

### THESIS 3 — **One efflux flux**, **two clocks**

| | Molecular weight | Glomerular filtration | Half-life | Role |
|---|---|---|---|---|
| **Creatine kinase** | 86 kDa | No | ~36 h | The **integrator** of injury — peaks at 37.5 h |
| **Myoglobin** | 17.8 kDa | Freely | ~2.4 h | The **speedometer** of instantaneous rate — peaks at 11.8 h |

The **25.7-hour lag** between the two peaks is an output of the model. Toxicity is set by
clock 2, while the test a clinician orders sits on clock 1. That is why a CK threshold
discriminates poorly — in this model, **even with CK fixed at 63,000 U/L**, the outcome can be
KDIGO 1 or KDIGO 3 depending purely on when fluids were started.

---

## 2. Four products and sums that explain the evidence

**(1) One saturable transporter splits a single filtered load into two toxic arms.**
`Reab = FL / (1 + FL/Tm)`, Tm ≈ 40 mg/h. Below Tm, almost everything is endocytosed (arm A:
haem), so the urine dipstick can be negative — why 20-50% of real rhabdomyolysis is
myoglobinuria dipstick-negative — and above Tm, the excess spills over to the distal nephron
(arm B: casts). Moreover, because per-cell load is `Reab ÷ nephron mass`, **the surviving
tubules inherit the share of the dead ones**: per-nephron injury accelerates without any
additional assumption.

**(2) Cast formation scales with the square of the distal concentration.** Doubling urine flow
cuts cast formation **fourfold**, but cast **clearance** is first-order in flow and slow (t½ of
days). Prevention structurally beats reversal.

**(3) The three injury arms add, so blocking one arm is capped by that arm's *share*** — and the
model computes that share rather than leaving it to intuition. Against a background of
200 mL/h, the haem arm carries about two-thirds of cumulative exposure (69.2 → 22.8 with
complete blockade). So **complete** blockade is worth a great deal (peak creatinine
2.54 → 1.81). The problem is not the size of the ceiling but the shape of the dose-response:

| Haem-arm blockade | 0% | 50% | 80% | 95% |
|---|---|---|---|---|
| Peak creatinine (mg/dL) | 2.54 | 2.08 | 1.87 | 1.81 |
| Cumulative exposure | 69.2 | 44.0 | 28.9 | 22.8 |
| eGFR d90 (mL/min) | 98.6 | 108.5 | 116.3 | 120.2 |

**A 50% target occupancy** — close to what antioxidants and iron chelators actually achieve in
vivo — **only buys 2.54 → 2.08**, and once flow is adequate, this arm's own share shrinks
further. This is a more specific, more falsifiable account of the negative antioxidant trials
than simply "the effect is small."

**(4) Urine pH has an internal optimum.** Acidic urine liberates ferrihaemate and precipitates
uric acid, while alkaline urine precipitates calcium phosphate (the HPO₄²⁻ fraction rises
sharply above pH 6.8). The combined index of the three chemistries is minimised at **urine pH
6.65** — the clinical guideline target is **derived**, not asserted.

| Urine pH | Uric acid solubility (mg/dL) | HPO₄²⁻ fraction | f_diss(haem) | Combined risk index |
|---|---|---|---|---|
| 5.0 | 7.7 | 0.016 | 0.969 | 4.92 |
| 5.6 | 11.1 | 0.059 | 0.760 | 3.59 |
| 6.2 | 24.8 | 0.201 | 0.240 | 1.89 |
| **6.65** | **58.1** | **0.414** | **0.053** | **1.48 ← minimum** |
| 7.1 | 152 | 0.666 | 0.010 | 1.67 |
| 7.7 | 586 | 0.888 | 0.001 | 2.01 |

The safe upper bound **depends on flow**. In the product matrix above, an additional
17.5 mmol/h of bicarbonate lowers cumulative exposure at 500 and 1000 mL/h (29.2 → 25.9,
32.0 → 26.5) but **raises it at 100 mL/h** (36.3 → 40.8), because high flow washes out
calcium-phosphate supersaturation and low flow does not. In that cell the two metrics diverge —
peak creatinine still falls a little (1.89 → 1.65) while cumulative exposure rises. Not a clean
win but **the honest signature of a trade-off**, and the reason the target is a **band** around
6.5 rather than "as alkaline as possible."

---

## 3. File layout

| File | Contents |
|---|---|
| [`rhab_qsp_model.dot`](rhab_qsp_model.dot) | Mechanistic map source — 158 nodes · 18 clusters |
| [`rhab_qsp_model.svg`](rhab_qsp_model.svg) / [`.png`](rhab_qsp_model.png) | Rendered map (150 dpi) |
| [`rhab_mrgsolve_model_en.R`](rhab_mrgsolve_model_en.R) | 47 ODEs · 187 parameters · 18 scenarios + calibration notes |
| [`rhab_shiny_app_en.R`](rhab_shiny_app_en.R) | 10-tab interactive dashboard |
| [`rhab_references_en.md`](rhab_references_en.md) | 80 references (every PMID looked up and verified) |

```r
# Run the model
source("rhab_mrgsolve_model_en.R")
sims <- run_all()                      # 18 scenarios, 90 days
print(summarise_scenario(sims), width = Inf)

# Dashboard
shiny::runApp("rhab_shiny_app_en.R")
```

Regenerate the mechanistic map:

```bash
dot -Tsvg rhab_qsp_model.dot -o rhab_qsp_model.svg
dot -Tpng -Gdpi=150 rhab_qsp_model.dot -o rhab_qsp_model.png
```

---

## 4. Model structure (47 ODEs)

| Block | Compartments |
|---|---|
| Muscle injury and efflux | `MISCH` `MNEC` `MLYS` `MDEB` |
| The two clocks | `CK` `MBP` (+ `MBFILT` `CKAUC`) |
| Kidney | `NEPH` `TUBI` `ROS` `CAST` `FIB` `ET1` |
| Volume and pressure | `ECFV` `ICFV` `MSEQ` `UOUT` |
| Electrolytes and acid-base | `NAE` `CLE` `HCOE` `KE` `POE` `CAE` `CADEP` `URT` `KREM` |
| Nitrogenous metabolites | `CRE` `CRN` `URN` |
| Inflammation | `IL6` `MPH` |
| Drugs | `MANC` `MANCUM` `FURC` `FURP` `ACZ` `DANT` `INS` `AOX` `CUMBIC` |
| Cumulative metrics | `TOXAUC` `ANURH` `HKH` `ARRH` `PHUAUC` `CUMFL` |

**Computed as outputs** (not inputs): KDIGO stage, the McMahon score, eGFR, the CK:myoglobin
ratio, urine pH, the osmolar gap.

---

## 5. 18 treatment scenarios — actual output

The figures below were obtained by extracting the `$GLOBAL`/`$MAIN`/`$ODE`/`$TABLE` blocks from
`rhab_mrgsolve_model_en.R` verbatim and independently compiling and integrating them (RK4,
dt = 0.005 h, 90 days). That is, they are **the output of the code as written**.

| # | Scenario | Peak CK | Mb | Peak SCr | K⁺ | iCa | Cl⁻ | Urine pH | 24h urine volume | Anuria | KDIGO | McM | eGFR d90 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | Healthy control (no injury) | 120 | 0.0 | 0.90 | 4.00 | 1.20 | 103.0 | 5.99 | 1441 | 0 | 0 | 0 | 125.0 |
| 2 | Crush injury 6 kg, **no treatment** | 63,869 | 84.1 | 9.06 | 8.07 | 0.95 | 108.5 | 5.36 | 584 | 46 | 3 | 11.0 | 84.8 |
| 3 | NS 200 mL/h, started at 24 h | 63,761 | 84.1 | 5.78 | 8.07 | 0.91 | 112.2 | 5.36 | 584 | 14 | 3 | 11.0 | 92.7 |
| 4 | NS 500 mL/h, immediately after extrication | 63,372 | 59.0 | 2.15 | 5.78 | 0.88 | 110.3 | 5.46 | 4,878 | 0 | 2 | 7.0 | 100.7 |
| 5 | NS 500 + NaHCO₃ 12.5 | 63,362 | 58.1 | 1.65 | 5.59 | 0.86 | 107.4 | 6.38 | 5,648 | 0 | 1 | 7.0 | 112.5 |
| 6 | **Better protocol** (1 L/h on scene + alkali) | 63,360 | 53.9 | 1.59 | 4.95 | 0.84 | 107.5 | 6.52 | 9,796 | 0 | 1 | 7.0 | 113.8 |
| 7 | Hartmann's solution 500 mL/h | 63,383 | 59.3 | 1.69 | 6.04 | 1.20 | 106.1 | 6.50 | 4,513 | 0 | 1 | 5.0 | 112.9 |
| 8 | Mass-casualty — 100 mL/h + alkali | 63,470 | 63.4 | 1.88 | 6.62 | 0.98 | 105.3 | 6.56 | 1,813 | 0 | 2 | 8.0 | 112.5 |
| 9 | Add mannitol 4 g/h | 63,377 | 58.7 | 2.17 | 5.61 | 0.87 | 110.1 | 5.47 | 5,504 | 0 | 2 | 7.0 | 98.4 |
| 10 | Mannitol 12 g/h — **osmotic nephrosis** | 63,387 | 58.4 | 5.18 | 6.68 | 0.86 | 110.3 | 5.49 | 6,702 | 0 | 3 | 7.0 | **22.7** |
| 11 | Add furosemide 20 mg/h | 63,837 | 63.4 | 2.98 | 5.29 | 0.83 | **143.8** | 5.37 | 9,477 | 0 | 2 | 8.5 | 95.9 |
| 12 | High-cutoff CVVH 40 mL/min (9 h) | 63,852 | 61.2 | 2.31 | 6.35 | 1.01 | 106.4 | 5.74 | 611 | 12 | 2 | 8.0 | 103.1 |
| 13 | Insulin/glucose (K⁺ shift) | 63,761 | 84.1 | 5.78 | **7.51** | 0.91 | 112.2 | 5.36 | 584 | 14 | 3 | 11.0 | 92.7 |
| 14 | Tight compartment + fasciotomy at 10 h | 105,782 | 86.7 | 2.72 | 8.11 | 0.74 | 110.1 | 5.20 | 13,230 | 0 | 2 | 10.5 | 96.9 |
| 15 | Tight compartment, no fasciotomy | 41,257 | 80.7 | 1.18 | 6.13 | 0.72 | 110.1 | 5.14 | 12,973 | 0 | 0 | 7.0 | 120.4 |
| 16 | Statin myotoxicity, 0.45 kg | 6,352 | 2.4 | 1.08 | 4.05 | 1.00 | 104.8 | 5.84 | 2,540 | 0 | 0 | 0 | 118.8 |
| 17 | Exertional rhabdomyolysis | 41,404 | 28.5 | 1.89 | 4.42 | 0.90 | 105.8 | 5.65 | 7,103 | 0 | 2 | 2.0 | 103.6 |
| 18 | Exertional + dantrolene | 21,224 | 13.5 | 1.44 | 4.15 | 0.91 | 105.8 | 5.73 | 7,320 | 0 | 1 | 2.0 | 110.1 |

Units: CK U/L · Mb mg/L · SCr mg/dL · K⁺, iCa, Cl⁻ mmol/L · urine volume mL · anuria h ·
eGFR mL/min

### What to read from the table

- **Timing beats dose.** Starting the same 500 mL/h at 0, 4, 8, 12, 24, or 48 hours gives peak
  creatinine of 2.13 → 2.13 → 2.15 → 2.44 → 4.89 → 7.50 mg/dL. **Note where the cliff sits:
  almost nothing is lost in the first 4-8 hours, and almost everything is lost between 12 and
  24 hours.** The reason is the squared term in cast formation plus the slowness of cast
  clearance — exposure integrated during a low-flow window cannot later be undone.
- **Crystalloid choice is a trade between two terms.** At the same 500 mL/h, normal saline gives
  urine pH 5.46, Cl⁻ 110.3, creatinine 2.15; Hartmann's solution gives urine pH 6.50, Cl⁻ 106.1,
  creatinine 1.69. Tracking bicarbonate as an **amount** while letting a chloride-rich fluid
  expand the distribution volume makes hyperchloraemic acidosis **emerge from the arithmetic** —
  no "dilutional acidosis" term was written anywhere.
- **What extracorporeal therapy actually does.** Decomposing it by setting the membrane's
  myoglobin sieving coefficient to zero: no CRRT (Mb peak 64.7, creatinine 2.54, cumulative
  exposure 69.2) → CRRT with zero sieving coefficient (63.9, 1.35, 38.7) → full CRRT
  (55.4, 1.34, 35.8). Myoglobin removal accounts for **2.9 (about 9%)** of the 33.4 exposure
  units CRRT removes, and for **0.01 (essentially none)** of the 1.20 mg/dL creatinine
  improvement. The rest is bicarbonate and potassium correction. The ceiling, too, is
  arithmetic: 40 mL/min = 2.4 L/h, but endogenous non-renal clearance is already
  `KELMB × VMB` = 2.4 L/h, so even complete anuria caps the benefit at a twofold increase — and
  even that only starts after the myoglobin peak has passed. Starting at 24 hours moves the peak
  not at all.
- **Risk score and outcome move independently.** The McMahon score is identically 7.0 across
  every well-treated crush scenario, yet the 90-day eGFR of those same scenarios spans
  98-114 mL/min. Five of the seven points for a crush patient (3 for cause + 2 for CK >
  40,000) are fixed before any treatment decision is made. **It is a baseline-risk tool, not a
  treatment-response tool** — nailed down by the fact that scenario 15 scores 7.0 while
  achieving KDIGO 0.
- **Recovery is biphasic.** Ionised calcium falls to 0.88 mmol/L over the first two days, then
  rebounds to 1.41 mmol/L as macrophages clear the debris scaffold and mobilise up to 331 mmol
  of deposited calcium phosphate. This is why giving calcium early, in the absence of
  arrhythmia, magnifies the later rebound — the administered calcium enters the same deposition
  compartment.

---

## 6. Validation — eight defects found by an independent reimplementation

All 47 ODEs were **independently reimplemented** in Python/scipy and run against clinical anchor
points; the process surfaced and led to the correction of the defects below. This list is kept
because, left silent, defects of this kind turn part of a "treatment effect" into baseline
drift.

1. **The crush trigger never fired at all.** Burial ischaemia had been written as a *consequence*
   of compartment pressure, but compartment pressure is near zero at presentation. Fixed by
   making external compression an explicit pressure variable — the two-gate structure and the
   extrication surge then followed on their own.
2. **Renal solute excretion was a set-point power law independent of urine flow.** As a result,
   chloride reached 257 mmol/L and anuria failed to stop potassium excretion. Rewritten as
   `urine flow × a regulated urine concentration`.
3. **Urine was already crystallising at baseline.** Uric acid supersaturation was referenced to
   thermodynamic solubility rather than the metastable limit, and the calcium-phosphate
   threshold was set **below the normal Ca × PO₄ product**, so a healthy patient kept
   depositing calcium.
4. **The injury driver responded to the absolute value rather than the increment.** An uninjured
   patient lost 0.6% of nephrons per day. All three arms were redefined as **increments** above
   the healthy baseline.
5. **Serum creatinine carried a permanent increase in production proportional to cumulative
   lysed mass.** It never normalised at all. Replaced with an explicit, consumed creatine
   compartment.
6. **Capillary leak was driven by cumulative lysed mass, so third-spacing and compartment
   pressure never resolved.** This produced a **spurious late relapse**, with nephron mass
   falling from 0.75 to 0.33 between days 7 and 30. Leak and macrophage influx now track an
   **active debris pool**, and macrophages are recruited by the very debris they clear.
7. **Osmotic nephrosis was linear in cumulative mannitol.** A standard dose of 100 g/day was
   enough to require dialysis. Replaced with a threshold at roughly 200 g cumulative.
8. **Extracorporeal therapy removed myoglobin, potassium, and urea but not creatinine, and
   loop-diuretic free-water loss did not carry chloride along with it.** Both were simple
   omissions, and the second **flipped the sign** of the predicted acid-base effect.

Not a bug but a structural omission: there was no **intracellular water compartment**, so
nothing constrained plasma tonicity. ICF water is now tracked explicitly, and ICF ⇄ ECF osmotic
shifts buffer sodium the way they do in vivo.

**Cross-validation.** After incorporating the fixes above, the C++ block was extracted from
`rhab_mrgsolve_model_en.R`, compiled with `g++ -Wall` (no warnings, all 47 derivatives confirmed
assigned), and independently integrated with RK4; the result agreed with the Python
implementation to 3-4 significant figures (e.g. scenario 4 — CK 63,372 vs 63,368 U/L, creatinine
2.15 vs 2.15, K⁺ 5.78 vs 5.78, eGFR d90 100.7 vs 100.7).

---

## 7. Limitations

- **A single-compartment muscle** with one pooled compartment pressure. In reality each limb
  has several independent compartments.
- **No thermoregulation** — heat stroke and malignant hyperthermia enter only as inputs.
- **No sepsis, wound infection, or amputation.** This is the most important reason the "no
  fasciotomy" scenario cannot be read as a clinical recommendation.
- **No coagulopathy/DIC** (which genuinely accompanies crush syndrome).
- No **drug-specific PBPK** for statins, only a myotoxicity-driving term.
- **Death** is reported as a cumulative arrhythmia risk rather than a calibrated survival model.
- Scenario 10 (mannitol 12 g/h) reaches an HCO₃⁻ of 3.7 mmol/L, a **deliberately exaggerated
  counterfactual** outside the range of survival. In clinical practice, osmolar-gap monitoring
  would stop it well before that point.
- Parameters are fitted to **published central-tendency values**, not to any specific patient
  dataset.

**This model is intended for education and hypothesis generation, not as a clinical
decision-making tool.**

---

## 8. References

[`rhab_references_en.md`](rhab_references_en.md) — 80 references across 14
sections. Every PMID was looked up via the NCBI E-utilities to confirm title, journal, and year,
and none was fabricated. The same file also holds a table mapping the clinical anchor points the
model is calibrated to reproduce to their supporting sections.

Key references:
- Bosch X et al. *Rhabdomyolysis and acute kidney injury.* N Engl J Med 2009 — [PMID 19571284](https://pubmed.ncbi.nlm.nih.gov/19571284/)
- Better OS, Stein JH. *Early management of shock and prophylaxis of acute renal failure in traumatic rhabdomyolysis.* N Engl J Med 1990 — [PMID 2407958](https://pubmed.ncbi.nlm.nih.gov/2407958/)
- Brown CV et al. *Preventing renal failure in patients with rhabdomyolysis: do bicarbonate and mannitol make a difference?* J Trauma 2004 — [PMID 15211124](https://pubmed.ncbi.nlm.nih.gov/15211124/)
- Sanders PW et al. *Mechanisms of intranephronal proteinaceous cast formation by low molecular weight proteins.* J Clin Invest 1990 — [PMID 2298921](https://pubmed.ncbi.nlm.nih.gov/2298921/)
- Gburek J et al. *Renal uptake of myoglobin is mediated by the endocytic receptors megalin and cubilin.* Am J Physiol Renal Physiol 2003 — [PMID 12724130](https://pubmed.ncbi.nlm.nih.gov/12724130/)
- McMahon GM et al. *A risk prediction score for kidney failure or mortality in rhabdomyolysis.* JAMA Intern Med 2013 — [PMID 24000014](https://pubmed.ncbi.nlm.nih.gov/24000014/)
- Bywaters EG, Beall D. *Crush injuries with impairment of renal function (1941).* J Am Soc Nephrol 1998 reprint — [PMID 9527411](https://pubmed.ncbi.nlm.nih.gov/9527411/)
