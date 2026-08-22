# Tumour Lysis Syndrome (TLS) — QSP Model

> **A race, not a syndrome.** TLS is usually taught as "four laboratory
> abnormalities (K↑ · PO₄↑ · urate↑ · Ca↓) plus acute kidney injury". This model
> instead poses it as **a race between two rates of the same dimension and of
> comparable magnitude**, and closes that race into **a loop** through the kidney.

```
        release flux        J_rel = Q_i · k_lys · N_lys              [mmol/h]
        clearance capacity  C_i(GFR, urine pH, urine flow rate)      [mmol/h]

  J_rel > C_i  →  concentration↑  →  supersaturation↑  →  precipitation  →  GFR↓  →  C_i↓  →  ...
```

In this loop every therapy is classified by **which term it touches**, and the
central claim of this model is that **those classes do not substitute for one
another**.

| Operator class | What it does | Examples |
|---|---|---|
| **POOL** | removes solute that already exists | rasburicase, dialysis |
| **FLUX** | blocks new generation | allopurinol, febuxostat |
| **FLUX-SHAPING** | same total amount, lowers the *rate* of release | venetoclax ramp, steroid prephase |
| **DILUTION** | lowers only the tubular concentration | hydration, (conditionally) furosemide |
| **SPECIATION** | shifts the solubility curve | bicarbonate → urine pH |
| **REDISTRIBUTION** | movement between compartments (not removal) | insulin/glucose, β₂ agonists |
| **SEQUESTRATION** | excretion out of the body | sevelamer, SZC/patiromer |

---

## Deliverables

| File | Contents |
|---|---|
| [`tls_qsp_model.dot`](../../../tumour-lysis-syndrome/tls_qsp_model.dot) · [SVG](../../../tumour-lysis-syndrome/tls_qsp_model.svg) · [PNG](../../../tumour-lysis-syndrome/tls_qsp_model.png) | Mechanistic map — 22 clusters / 211 nodes / 320 edges. Therapy nodes are coloured by the operator classes above |
| [`tls_mrgsolve_model.R`](../../../tumour-lysis-syndrome/tls_mrgsolve_model.R) | 48-ODE mrgsolve model (PK 19 · tumour/release 3 · solutes 11 · kidney/crystal 7 · safety 7 · accounting 1), 12 scenarios, 18 analysis functions |
| [`tls_shiny_app.R`](../../../tumour-lysis-syndrome/tls_shiny_app.R) | 10-tab interactive dashboard |
| [`tls_references.md`](../../../tumour-lysis-syndrome/tls_references.md) | 107 references, sources classified by parameter |
| [`tls_reference_check.py`](../../../tumour-lysis-syndrome/tls_reference_check.py) | **An independent numpy/scipy port of the same equation system.** Every number in this README was computed here, and can be reproduced without R |

The build environment of this repository has no R toolchain. Rather than
publishing an ODE model that has never been integrated, the same system was
ported a second time and integrated with scipy LSODA, and every number below is
output from `python3 tls_reference_check.py`. If the two ports disagree, one of
them is wrong.

```bash
python3 tls_reference_check.py            # everything (A0–A18)
python3 tls_reference_check.py --only A1  # core result only
```

---

## 0. Calibration check (baseline)

The result of integrating 600 hours in the normal state, with neither tumour nor
drug (`A0`):

| Quantity | Model | Target |
|---|---|---|
| Uric acid | 5.17 mg/dL | 4–6 |
| Potassium | 4.10 mmol/L | 3.5–5.0 |
| Phosphate | 1.15 mmol/L | 0.8–1.45 |
| Total calcium / ionised calcium | 2.35 / 1.17 mmol/L | 2.2–2.6 / 1.12–1.30 |
| Creatinine | 0.86 mg/dL | 0.7–1.2 |
| eGFR | 120 mL/min | ~120 |
| Urine volume | 1.99 L/day | 1.5–2.5 |
| Uric acid excretion | 715 mg/day | ~700 |
| Potassium excretion | 85 mmol/day | 70–100 |
| Phosphate excretion | 23 mmol/day | 20–30 |
| Calcium excretion | 4.9 mmol/day | 2.5–7.5 |
| Urine uric acid supersaturation | 1.89 | 1.5–2.5 (metastable) |
| Tubular injury / nephron mass | 0.000 / 1.000 | stationary |

It matters that normal urine is *already* supersaturated with respect to uric
acid (SS 1.89). Crystallisation begins not because of **supersaturation** but
**when the metastable limit (SS ≈ 2.5) is exceeded**.

---

## 1. Converting the race into two concentrations — the central result of this model

For uric acid, two threshold concentrations can be solved in closed form.

- **UA_req** — the plasma uric acid required to offset the maximum purine release
  flux by *excretion*. (It reflects the fractional excretion of uric acid rising
  from 8% to ~30% as URAT1/GLUT9 reabsorption saturates. This is the only defence
  the kidney has.)
- **UA_crit** — the plasma uric acid at which the tubular fluid reaches the
  metastable limit and uric acid begins to precipitate.

**The disease is the region where UA_req > UA_crit: the state in which the kidney
cannot excrete the load without making crystals.**

### UA_crit is a property of the prescription, not of the tumour

| Prescription | UA_crit |
|---|---|
| 2 L/day, pH 5.9 | **6.8 mg/dL** |
| 3 L/m²/day, pH 5.9 | **9.0 mg/dL** |
| 4 L/m²/day, pH 5.9 | **10.5 mg/dL** |
| 3 L/m²/day, pH 7.0 | **33.2 mg/dL** |

The Cairo–Bishop uric acid criterion for laboratory TLS is **8.0 mg/dL**. This
model reproduces that value as a **crystallisation threshold** computed from the
solubility law, the urine flow rate, and saturable uric acid handling alone. It
was not tuned to the criterion — 8.0 is not an input to this model but a value
that falls inside the range of its output.

### The race as a function of tumour load (3 L/m²/day, pH 5.9)

| Load (10¹² cells) | Urate potential | J_rel (mmol/h) | UA_req | UA_crit | req/crit | Actual peak urate | Peak K | Cr fold | TLS |
|---|---|---|---|---|---|---|---|---|---|
| 0.05 | 0.2 g | 0.04 | 1.2 | 9.0 | 0.13 | 5.8 | 4.12 | 1.00 | – |
| 0.30 | 1.2 g | 0.24 | 6.8 | 9.0 | 0.76 | 7.8 | 4.24 | 1.01 | – |
| 0.50 | 2.0 g | 0.40 | 8.1 | 9.0 | 0.90 | 8.9 | 4.33 | 1.02 | – |
| 1.00 | 4.0 g | 0.81 | 11.1 | 9.0 | 1.23 | 11.3 | 4.55 | 1.05 | – |
| 2.00 | 8.1 g | 1.62 | 16.7 | 9.0 | 1.85 | 15.3 | 4.97 | 1.14 | L |
| 3.00 | 12.1 g | 2.42 | 22.2 | 9.0 | 2.46 | 20.2 | 5.48 | 1.55 | **C** |
| 5.00 | 20.2 g | 4.04 | 32.9 | 9.0 | 3.64 | 39.9 | 7.08 | 8.26 | **C** |
| 8.00 | 32.3 g | 6.46 | 48.8 | 9.0 | 5.40 | 86.1 | 9.71 | 17.63 | **C** |

(L = laboratory TLS, C = clinical TLS)

**Independent verification.** UA_req is obtained algebraically and the "actual
peak urate" by integrating the 48-ODE system, so the two computational paths are
entirely different. And yet:

| Load | UA_req (closed form) | Peak urate (ODE) | Ratio |
|---|---|---|---|
| 0.5 | 8.1 | 8.9 | 1.09 |
| 1.0 | 11.1 | 11.3 | 1.01 |
| 2.0 | 16.7 | 15.3 | 0.92 |
| 3.0 | 22.2 | 20.2 | 0.91 |

Above 3×10¹² cells the two values diverge, and they diverge because the ODE loses
GFR while the closed form assumes GFR is intact. **That discrepancy *is* the
loop.**

---

## 2. The ranking of the three solutes — only uric acid loses the race

Computing the same race for the three solutes (load 3×10¹² cells, GFR intact):

| Solute | Intracellular content per 10¹² cells | Max release flux | Renal capacity | Flux/capacity |
|---|---|---|---|---|
| Uric acid | 24.0 mmol | 2.42 mmol/h | 0.39 | **6.2** |
| Potassium | 38.5 mmol | 3.89 mmol/h | 9.25 | **0.4** |
| Phosphate | 60.0 mmol | 6.06 mmol/h | 12.0 | **0.5** |

**While GFR is intact, uric acid is the only one of them that loses the race.**
Potassium and phosphate are excretable loads. So in this model hyperkalaemia and
hyperphosphataemia are not release phenomena but **secondary phenomena that arise
after urate crystals have taken the GFR away**.

This is a falsifiable prediction. A drug with no pharmacological action
whatsoever on potassium or phosphate must lower potassium and phosphate.
Rasburicase is that test (`A5c`):

| Prescription | Peak urate | Peak K | Peak PO₄ | Ionised Ca | eGFR nadir | P(arrhythmia) |
|---|---|---|---|---|---|---|
| No prophylaxis, 2 L/day (3×10¹²) | 25.0 | 6.07 | 2.49 | 1.10 | 37 | 2.86% |
| **+ rasburicase alone** | 13.9 | **5.38** | **1.92** | 1.17 | 85 | 1.85% |
| + febuxostat alone | 11.9 | **5.41** | **1.95** | 1.17 | 98 | 2.00% |
| No prophylaxis (6×10¹²) | 72.1 | 9.06 | 7.40 | 0.44 | 5 | 99.07% |
| **+ rasb 0.4 mg/kg ×5 alone** | 20.2 | **7.00** | **3.14** | 0.97 | 66 | 27.00% |

In an arm given nothing but the urate enzyme, potassium and phosphate move. That
is, in this model the hyperkalaemia of TLS **does not track LDH; it tracks
creatinine.**

---

## 3. Urine pH — the trade the alkalinisation era never priced

Uric acid solubility rises with pH (pKa 5.75). But so does the driving force for
calcium phosphate precipitation — because the precipitating species is HPO₄²⁻ and
its pKa₂ is 6.80 — and systemic alkalosis simultaneously lowers ionised calcium.
Pricing both sides:

| Urine pH | Urate solubility (mmol/L) | UA_crit (mg/dL) | HPO₄²⁻ fraction | vs pH 5.9 |
|---|---|---|---|---|
| 5.5 | 1.02 | 7.6 | 0.048 | urate ×0.6 / CaP ×0.4 |
| 5.9 | 1.58 | 9.0 | 0.112 | reference |
| 6.5 | 4.34 | 15.5 | 0.334 | urate ×2.7 / CaP ×3.0 |
| 7.0 | 12.30 | 33.2 | 0.613 | urate ×7.8 / CaP ×5.5 |
| 7.5 | 37.49 | 87.9 | 0.834 | urate ×23.7 / CaP ×7.5 |
| 7.8 | 74.15 | 167.0 | 0.909 | urate ×46.9 / CaP ×8.1 |

Which side wins is not a matter of opinion; it turns on **which solute is
rate-limiting**. So the optimal pH differs from patient to patient (`A3`):

| Situation | Optimal urine pH | Cr fold there | Cr fold at pH 5.90 | at pH 7.46 |
|---|---|---|---|---|
| Load 3×10¹², no rasburicase | **6.60** | 1.15 | 1.55 | 1.31 |
| Load 3×10¹², with rasburicase | **5.90 (= no alkalinisation)** | 1.15 | 1.15 | 1.30 |
| Load 6×10¹², no rasburicase | **7.00** | 2.57 | 16.50 | 2.88 |
| Load 6×10¹², with rasburicase | **6.60** | 2.24 | 5.00 | 2.88 |

**Erase the benefit on the urate side and only the cost remains in the pH term.**
This is where the abandonment of routine alkalinisation in a patient who has
urate oxidase is *derived* — rather than asserted. The cost appears along two
branches: renal calcium phosphate deposition rises 0.00 → 1.97 mmol, and ionised
calcium falls 1.17 → 0.97 mmol/L.

The same logic accounts for the benefit of alkalinisation being dramatic at a
load of 6×10¹² (16.50 → 2.57). At that load uric acid is overwhelmingly the
rate-limiting step, so it more than covers the CaP price.

---

## 4. Rasburicase is a zero-order POOL operator **of fixed capacity**

The Km of urate oxidase for uric acid (~25 µmol/L) is two orders of magnitude
below the uric acid concentrations of TLS. The enzyme **runs saturated.** So one
dose buys not a *fixed fractional reduction* but a **fixed mmol/h** — and a fixed
mmol/h is a value a tumour can exceed (`A5b`).

| Dose | C_rasb | Processing capacity | Load it covers |
|---|---|---|---|
| 0.05 mg/kg | 0.44 mg/L | 0.66 mmol/h | 0.8 ×10¹² cells |
| 0.10 | 0.88 | 1.31 | 1.6 |
| 0.15 | 1.31 | 1.97 | 2.4 |
| **0.20 (licensed dose)** | 1.75 | **2.62** | **3.2** |
| 0.30 | 2.62 | 3.94 | 4.9 |
| 0.40 | 3.50 | 5.25 | 6.5 |

Running the licensed dose of 0.20 mg/kg across loads in the ODE:

| Load | 4-hour urate | 4-hour change | **Peak urate** | Cr fold |
|---|---|---|---|---|
| 0.5 ×10¹² | 0.13 | **−98%** | 5.2 | 1.02 |
| 1.0 | 0.19 | **−96%** | 5.2 | 1.05 |
| 3.0 | 0.58 | **−89%** | 14.8 | 1.15 |
| 6.0 | 1.46 | **−72%** | 37.1 | 6.25 |
| 10.0 | 2.94 | **−43%** | 94.7 | 18.64 |

**The 4-hour fall — the primary endpoint the registration trial used — remains
impressive at every load.** Four hours is time enough for the enzyme to empty the
uric acid pool that *already existed*. But what sets the renal outcome is what
the release flux does after that, and that is the difference between a peak urate
of 5.2 and one of 94.7 mg/dL. **That is, the registration endpoint measures
whether the drug can empty an already existing pool, and barely measures whether
the drug keeps up with the release flux.** This is the model's explanation for
rasburicase lowering uric acid dramatically while its benefit on hard endpoints
such as dialysis and death remains unestablished, and equally for why repeated
and escalated dosing is needed in the highest-risk patients.

Escalation is more effective than repetition (load 6×10¹², 2 L/day):

| Prescription | Cr fold | P(dialysis) |
|---|---|---|
| rasb 0.2 ×1 | 10.91 | 93.2% |
| rasb 0.2 ×5 | 5.00 | 41.8% |
| rasb 0.4 ×5 | **2.04** | **8.7%** |

---

## 5. Operator decomposition — which term of the loop each therapy touches

Load 3×10¹² cells (`A2`):

| Arm | Peak urate | Peak K | Peak PO₄ | Ionised Ca | **Cr fold** | Urate crystal | CaP crystal | P(dialysis) |
|---|---|---|---|---|---|---|---|---|
| Nothing at all (2 L/day) | 25.0 | 6.07 | 2.49 | 1.10 | **3.06** | 2.08 | 0.06 | 8.4% |
| DILUTION 3 L/m²/day | 20.2 | 5.48 | 1.99 | 1.17 | 1.55 | 0.06 | 0.00 | 0.9% |
| DILUTION 4 L/m²/day | 19.8 | 5.38 | 1.93 | 1.17 | **1.27** | 0.00 | 0.00 | 0.9% |
| DILUTION + furosemide | 20.7 | 5.17 | 2.02 | 1.17 | 1.61 | 0.06 | 0.00 | 0.9% |
| FLUX allopurinol t=0 | 17.2 | 5.39 | 1.93 | 1.17 | 1.26 | 0.00 | 0.00 | 0.9% |
| FLUX febuxostat t=0 | 13.5 | 5.36 | 1.92 | 1.17 | **1.16** | 0.00 | 0.00 | 0.9% |
| POOL rasburicase | 14.8 | 5.35 | 1.91 | 1.17 | **1.15** | 0.00 | 0.00 | 0.9% |
| POOL dialysis d2–5 | 20.2 | 5.48 | 1.99 | 1.17 | 1.52 | 0.00 | 0.00 | 0.9% |
| SPECIATION pH 7.5 | 21.0 | 5.37 | 1.66 | **1.03** | 1.18 | 0.00 | 0.35 | 0.9% |
| SEQUESTRATION phosphate binder | 20.2 | 5.48 | 1.89 | 1.17 | 1.55 | 0.06 | 0.00 | 0.9% |
| FLUX-SHAPE steroid prephase | 11.4 | 4.55 | 1.40 | 1.17 | **1.04** | 0.00 | 0.00 | 1.1% |

The two crystal columns have to be read side by side. **Not one arm lowers both
at the same time.**

It is also striking that the cheapest intervention is the most powerful. 2 L/day
→ 4 L/m²/day alone takes the Cr fold from 3.06 to 1.27. That is because
`UA_crit` rises 6.8 → 10.5 mg/dL, and it is obtained without touching any pool
and any flux.

**Furosemide has no net benefit in this model** (1.55 → 1.61). It raises the urine
flow rate, but the unreplaced output shrinks the extracellular fluid and shaves
the GFR haemodynamically, so the dilution benefit is cancelled by the volume
loss (potassium does come down 5.48 → 5.17, so it remains useful for potassium).
The clinical principle that furosemide is a dilution operator only in the
volume-replete patient is reproduced exactly by the model's volume state, and the
price is paid, inside the model, in GFR.

---

## 6. The calcium phosphate branch — and why phosphate binders are useless

The tubular injury driving force is the sum of four terms, and that decomposition
is output for every arm (`A6`). So the claim "there is residual injury" is
computed, not asserted.

| Arm | Peak driving force | Urate crystal | Xanthine | **CaP** | Soluble urate | Cr fold | P(dialysis) |
|---|---|---|---|---|---|---|---|
| Hydration only 3 L/m² (3×10¹²) | 1.65 | 91% | 0% | **0%** | 9% | 1.55 | 0.9% |
| Hydration only 2 L/day (3×10¹²) | 2.76 | 95% | 0% | **0%** | 5% | 3.06 | 8.4% |
| + rasburicase | 0.70 | 83% | 0% | 0% | 17% | 1.15 | 0.9% |
| + rasb + phosphate binder | 0.70 | 83% | 0% | 0% | 17% | **1.15** | 0.9% |
| High load 2 L/day (6×10¹²) | 6.88 | 74% | 0% | **22%** | 3% | 16.50 | 98.6% |
| High load + rasb ×1 | 5.69 | 73% | 0% | 23% | 4% | 10.91 | 93.2% |
| High load + rasb 0.4 ×5 | 4.09 | 79% | 0% | **17%** | 4% | 2.04 | 8.7% |
| … + maximal hydration | 1.94 | 75% | 0% | 16% | 9% | 1.36 | 2.1% |
| … + early dialysis | 1.94 | 75% | 0% | 16% | 9% | **1.31** | 1.2% |

**(1) The calcium phosphate branch exists, but it is gated by load.** 0% of the
peak driving force at 3×10¹², 22% at 6×10¹², and still 17% even after the enzyme
has suppressed uric acid as far as it can be suppressed. No urate-directed drug
touches that share.

**(2) Phosphate binders do not touch it either, and the reason is arithmetic, not
pharmacology.**

```
phosphate entering the ECF from lysis   12.12 mmol/h   (load 6×10¹²)
phosphate entering the ECF from diet     0.95 mmol/h
share an enteric binder can reach        7.3 %
```

And indeed, adding sevelamer leaves the Cr fold unmoved at 1.15 → 1.15.
**Sequestration is an operator acting on the intake term, and in acute TLS
phosphate does not come from the intake term.** This is where the reason
phosphate binders are effective in chronic kidney disease and useless in acute
TLS comes from — the only route that takes phosphate out of the body is dialysis.

**(3) The xanthine branch does not open in this parameterisation.** It is 0% in
every arm. At the level at which allopurinol inhibits XO by 84% and febuxostat by
88%, the residual XO capacity (5.4 and 4.0 mmol/h respectively) still exceeds the
purine flux (2.42 mmol/h), and because hypoxanthine is excreted well, at a
fractional excretion of 0.60, it carries the bypassed load instead. Peak xanthine
stays at 0.5 mg/dL. Xanthine nephropathy is a real, reported complication, but
this model positions it as **a branch that opens only when XO inhibition is
nearly complete.**

---

## 7. Allopurinol's lead time — a result opposite to the hypothesis posed

This analysis was written to express the hypothesis that "a flux operator cannot
empty a pool, so lead time matters". **The model does not support that
hypothesis** (`A4`).

| Allopurinol start | XO inhibition at t=0 | Peak urate | Peak xanthine | **Cr fold** |
|---|---|---|---|---|
| None | 0.000 | 20.2 | 0.1 | 1.55 |
| t = 0 | 0.000 | 17.2 | 0.4 | **1.26** |
| −12 h | 0.795 | 16.4 | 0.5 | **1.23** |
| −24 h | 0.731 | 16.4 | 0.5 | 1.21 |
| −48 h | 0.802 | 16.1 | 0.5 | 1.19 |
| −72 h | 0.824 | 15.9 | 0.5 | 1.19 |
| −120 h | 0.837 | 15.7 | 0.5 | 1.18 |
| −168 h | 0.840 | 15.7 | 0.5 | 1.18 |

**Twelve hours captures most of the benefit, and moving the start a full week
earlier gains only 0.05 more.** The arithmetic reason is printed directly by the
function:

```
pre-existing miscible urate pool     5.4 mmol  (0.90 g)
purine released by lysis            72.0 mmol  (12.1 g)
ratio                               13.4 : 1
```

The "pool it failed to empty" when the flux operator starts late is a mere 7% of
the amount released. So the value of lead time is set **not by the pool size but
by oxypurinol accumulation (t½ 23 h)** — dosing at t=0 starts from an XO
inhibition of 0.000, and starting 12 hours earlier starts from 0.795. The model's
reading is that the part of the customary "start 2–3 days beforehand"
recommendation that actually works is **the time for the active metabolite to
build up.**

---

## 8. The loop comes back around to the drug as well — oxypurinol in AKI

Oxypurinol is renally excreted and febuxostat is not. So **the renal injury
allopurinol is trying to prevent changes allopurinol's own exposure** (`A14`):

| Arm | d1 | d3 | **d7** | d7 XO inhibition | d7 eGFR |
|---|---|---|---|---|---|
| Allopurinol, low load | 5.50 | 5.75 | 5.82 mg/L | 0.841 | 120 |
| Allopurinol, high load | 6.10 | 14.80 | **30.33 mg/L** | 0.965 | **16** |
| Febuxostat, high load | 0.22 | 0.22 | 0.22 mg/L | 0.881 | 103 |

By day 7 oxypurinol has **accumulated 5-fold**. This is the model's explanation
for why dose reduction of allopurinol becomes necessary in precisely the patients
at the highest TLS risk, and for why the risk of allopurinol hypersensitivity
syndrome (HLA-B\*58:01, DRESS/SJS-TEN) rises in those patients. Febuxostat is
hepatically metabolised and so is not caught by the same feedback.

---

## 9. FLUX-SHAPING — why the venetoclax ramp exists

The same drug, the same target, and in the end the same kill. **Only the release
rate differs.** That is because the clearance system is a low-pass filter on the
release flux. CLL kinetics (doubling time 2000 h), load 5×10¹² cells, 2 L/day
(`A8`):

| Venetoclax schedule | **Max J_rel** | Peak K | Peak PO₄ | Peak urate | **Cr fold** | P(dialysis) | TLS | Total lysed | Residual at 6 weeks |
|---|---|---|---|---|---|---|---|---|---|
| 400 mg from day 1 | 2.85 | 6.62 | 2.90 | 32.5 | **6.85** | 74.2% | **C** | 5.04 | 0.0000 |
| 200 mg from day 1 | 2.37 | 6.33 | 2.66 | 29.2 | 4.91 | 49.6% | C | 5.05 | 0.0000 |
| 100 mg from day 1 | 1.78 | 5.97 | 2.36 | 24.6 | 3.38 | 19.2% | C | 5.07 | 0.0000 |
| 50 mg from day 1 | 1.40 | 5.51 | 1.94 | 18.8 | 2.57 | 3.3% | C | 5.10 | 0.0000 |
| 2-step 20/50 | 0.84 | 4.88 | 1.51 | 12.8 | 1.63 | 1.6% | C | 5.18 | 0.0000 |
| **5-week label ramp** | 0.84 | 4.88 | 1.51 | 12.8 | **1.63** | 1.6% | C | 5.18 | 0.0000 |
| 8-week slow ramp | 0.55 | 4.53 | 1.33 | 10.0 | **1.24** | 1.6% | – | 5.28 | 0.0671 |

Every arm lyses 5.0–5.3 ×10¹² cells and (all but the 8-week arm) has zero
residual tumour at week 6. **The ramp does not buy safety by giving up kill.**

And what stands out: **the 2-step ramp and the 5-week label ramp are completely
identical in every safety column.** That is because the maximum release flux is
set by the **first dose step** — by the time the later steps arrive, the tumour
that could have been lysed is already gone. This is the model's explanation for
why the label attaches its monitoring requirement to the 20 mg dose.

---

## 10. Redistribution is not removal — and it is worse than expected

Each rescue therapy is started at the potassium peak of the untreated trajectory
(t = 42 h, K 6.07, eGFR 44). Alongside serum potassium, the **total body
exchangeable potassium (extracellular plus what has been parked inside the
cells)** is output, so the difference between "shifting" and "removing" is
visible as a number rather than as an assertion (`A9`):

| Rescue therapy | K 2h | K 6h | **K 12h** | K 24h | 24h total body K | **ΔTBK** |
|---|---|---|---|---|---|---|
| Nothing at all | 6.07 | 6.02 | 5.90 | 5.61 | 94 | −8 |
| insulin 10 U + glucose | **5.47** | 5.45 | **6.06** | 5.78 | 97 | **−5** |
| insulin q4h ×6 | 5.51 | 5.32 | 5.49 | 5.52 | 113 | **+11** |
| salbutamol 20 mg nebulised | 5.75 | 5.68 | 5.66 | 5.56 | 105 | **+3** |
| insulin + salbutamol | 5.19 | 5.16 | 5.78 | 5.71 | 109 | **+7** |
| SZC 10 g q8h | 5.65 | 5.49 | 5.01 | 4.74 | 80 | **−22** |
| furosemide 40 mg q8h | 5.91 | 5.78 | 5.54 | 5.26 | 87 | **−15** |
| haemodialysis 4 h | 6.07 | **1.89** | 3.08 | 4.52 | 76 | **−26** |

Insulin's 2-hour fall is −0.60 mmol/L, which matches the reported −0.6 to −1.0.
But **at 12 hours it is 6.06, higher than the 5.90 of the arm that did nothing.**
And the ΔTBK column runs in the opposite direction to expectation: total body
potassium in the redistribution arms **rises more** than in the untreated arm
(−5, +11, +3, +7 versus −8).

This is not a result of the parameterisation but a result that comes out of the
structure. Distal tubular potassium secretion depends on the **concentration**
delivered to it. Park potassium inside the cells and serum potassium comes down;
and once serum potassium comes down, **the only route that was actually taking
potassium out switches off.** That is, redistribution lowers the number on the
chart, and it leaves the patient's potassium not merely where it was but a little
more accumulated than that. That is why the effect disappears, and it is the
quantitative statement of why shifting is **a measure that buys time** and not a
definitive one.

The only things that actually change total body potassium are dialysis, enteric
binders, and (in the volume-replete patient) diuretics.

---

## 11. The calcium reflex — correcting the hypocalcaemia feeds the crystal

The hypocalcaemia of TLS is a **consequence** of calcium phosphate precipitation.
So the calcium given in order to correct it is a **reactant** of that
precipitation reaction (`A10`, load 3×10¹²):

| Calcium infusion | Ionised Ca | Ca×PO₄ (mg²/dL²) | Renal CaP | **Cr fold** | eGFR nadir | **P(dialysis)** | P(seizure) |
|---|---|---|---|---|---|---|---|
| 0.0 mmol/h | 1.17 | 59 | 0.00 | 1.55 | 79 | 0.9% | 0.5% |
| 0.5 | 1.17 | 66 | 0.01 | 1.59 | 77 | 0.9% | 0.5% |
| 1.5 | 1.17 | 73 | 2.32 | 1.88 | 64 | 0.9% | 0.5% |
| 3.0 | 1.17 | 80 | 11.03 | **3.19** | 36 | **16.0%** | 0.5% |
| 6.0 | 1.17 | 89 | 10.75 | **7.89** | 10 | **82.0%** | 0.5% |

The seizure column is why that reflex exists, and the dialysis column is its
price. In this model the guideline recommendation not to give calcium for
asymptomatic hypocalcaemia is a direct consequence of the fact that **the Ca×PO₄
product is the extent of reaction.**

---

## 12. The loop is a knee, not a switch — a negative result

A positive feedback loop invites the expectation of bistability. **This loop is
not** (`A11`).

The nephron reserve axis (load fixed at 3×10¹²):

| Starting nephron mass | Starting eGFR | Cr fold | eGFR nadir | Nephrons at day 14 | P(dialysis) |
|---|---|---|---|---|---|
| 1.00 | 120 | 1.55 | 79.0 | 0.989 | 0.9% |
| 0.80 | 96 | 1.85 | 65.7 | 0.912 | 0.9% |
| 0.60 | 72 | 2.35 | 51.2 | 0.830 | 1.2% |
| 0.40 | 48 | 3.39 | 33.8 | 0.726 | 20.9% |
| 0.30 | 36 | 4.98 | 21.5 | 0.626 | 56.8% |

The insult axis (nephrons intact, 2 L/day):

| Load | Cr fold | eGFR nadir | Max crystal | Nephrons at day 14 | P(dialysis) |
|---|---|---|---|---|---|
| 1.0 | 1.08 | 112.2 | 3.05 | 0.997 | 0.9% |
| 2.0 | 1.71 | 69.7 | 9.93 | 0.987 | 0.9% |
| 3.0 | 3.06 | 37.1 | 17.48 | 0.975 | 8.4% |
| **3.5** | **4.86** | 21.4 | 26.58 | 0.956 | **45.9%** |
| 4.0 | 8.48 | 11.1 | 38.92 | 0.904 | 85.3% |
| 6.0 | 16.50 | 5.2 | 57.40 | 0.765 | 98.6% |

A hysteresis test — the load is fixed below threshold (2×10¹²) and crystal is
**seeded** into the kidney at t=0. In a bistable system a large enough seed ought
to push an otherwise safe patient onto the second branch and hold them there:

| Seed crystal | Cr fold | eGFR nadir | Crystal at day 14 | **eGFR at day 14** |
|---|---|---|---|---|
| 0 mmol | 1.14 | 109.7 | 0.02 | **119** |
| 5 | 1.36 | 91.8 | 0.04 | **119** |
| 10 | 1.63 | 70.8 | 0.06 | **119** |
| 20 | 2.45 | 31.8 | 0.13 | **119** |
| 40 | 10.25 | 9.0 | 5.21 | 106 |

Every seed comes back to the same eGFR by day 14. **The loop gain does not exceed
1, so there is only one attractor, and the steep rise around 3.5×10¹² cells is a
knee, not a switch.** The clinical implication is, if anything, optimistic —
crystal nephropathy has a narrow and steep risk band, but because it is not
self-sustaining it recovers once release ends. Even the untreated high-load
trajectory in this model comes back to an eGFR of 113 by day 14.

---

## 13. What the Cairo–Bishop definition actually detects

**A negative result is reported as a negative result.** The "two or more of four
within a 7-day window" rule is **robust** to sampling frequency in this model —
three of the four analytes have time constants on the order of days, and even the
fastest of them, potassium, has a peak broad enough that q24h sampling misses a
mere 0.06 mmol/L (`A12`).

What is meaningful is **which criterion fires**:

| Scenario | Urate | K | PO₄ | Ca | n | Peak urate | UA_crit |
|---|---|---|---|---|---|---|---|
| High load, no prophylaxis | ✓ | ✓ | ✓ | – | 3 | 25.0 | 9.0 |
| High load + rasburicase | ✓ | ✓ | ✓ | – | 3 | 14.8 | 9.0 |
| High load + full bundle | ✓ | ✓ | ✓ | – | 3 | 15.1 | 9.0 |
| Moderate + allo −72h | ✓ | – | – | – | 1 | 9.6 | 9.0 |
| Low load | ✓ | – | – | – | 1 | 5.8 | 9.0 |

The urate criterion (8 mg/dL) is effectively the same number as the UA_crit (9.0)
computed independently in §1. That is, in this model **the Cairo–Bishop urate
criterion is a crystallisation detector and the remaining three criteria are
outcome detectors.** Rasburicase turns clinical TLS into laboratory TLS but does
not change the *number* of criteria that fire — meaning that what the definition
counts and what the definition is trying to predict are not the same thing.

---

## 14. The safety branch — the oxidative load is proportional to the indication

Urate oxidase makes 1 mol of H₂O₂ per mol of uric acid it destroys. So **the
patient with the largest tumour takes the largest peroxide load — the toxicity is
proportional to the reason for giving the drug** (`A13`):

| Load | G6PD | Urate oxidised | **Peak MetHb** | Hb nadir |
|---|---|---|---|---|
| 0.3 ×10¹² | normal | 33.3 mmol | 0.0% | 13.5 |
| 0.3 | deficient | 33.3 | **8.2%** | 12.2 |
| 1.0 | deficient | 47.5 | **12.4%** | 11.6 |
| 3.0 | deficient | 56.8 | **16.4%** | 11.1 |
| 6.0 | deficient | 60.9 | **18.0%** | 11.0 |

The reason G6PD deficiency is an absolute contraindication comes out
quantitatively. There is also an implication on the sampling side: in blood drawn
after rasburicase, **uric acid keeps being degraded inside the tube as well**, so
ice, pre-chilled tubes, and immediate analysis are required (see cluster 22 of
the map).

---

## 15. The 12 scenarios

| # | Scenario | Peak urate | Peak K | Peak PO₄ | Ionised Ca | Cr fold | eGFR nadir | Peak LDH | Total lysed | Lysis at t<0 | P(dialysis) | TLS |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| S1 | Low-risk solid tumour (0.05×10¹²) | 5.8 | 4.12 | 1.16 | 1.17 | 1.00 | 120 | 341 | 0.06 | 0.00 | 0.9% | – |
| S2 | High load, hydration only 2 L/day | 25.0 | 6.07 | 2.49 | 1.10 | 3.06 | 37 | 8679 | 3.73 | 0.00 | 8.4% | **C** |
| S3 | + aggressive hydration 3 L/m² | 20.2 | 5.48 | 1.99 | 1.17 | 1.55 | 79 | 8679 | 3.73 | 0.00 | 0.9% | **C** |
| S4 | + allopurinol t=0 | 17.2 | 5.39 | 1.93 | 1.17 | 1.26 | 100 | 8679 | 3.73 | 0.00 | 0.9% | L |
| S5 | + allopurinol −72 h | 15.9 | 5.36 | 1.92 | 1.17 | 1.19 | 108 | 8679 | 3.73 | 0.00 | 1.0% | L |
| S6 | + rasburicase 0.2 mg/kg | 14.8 | 5.35 | 1.91 | 1.17 | 1.15 | 115 | 8679 | 3.73 | 0.00 | 0.9% | L |
| S7 | + rasb, maximal hydration, no alkali | 15.1 | 5.35 | 1.91 | 1.17 | 1.15 | 119 | 8679 | 3.73 | 0.00 | 0.9% | L |
| S8 | + alkalinisation pH 7.5 | 21.0 | 5.37 | 1.66 | 1.03 | 1.18 | 111 | 8679 | 3.73 | 0.00 | 0.9% | L |
| S9 | Steroid prephase 5 days | 11.4 | 4.55 | 1.40 | 1.17 | **1.04** | 120 | 3643 | 1.73 | **0.44** | 1.1% | – |
| S10 | Full bundle | 10.9 | 5.34 | 1.84 | 1.17 | 1.16 | 115 | 8910 | 3.73 | 0.00 | 1.0% | L |
| S11 | Stage 3 CKD host (nephrons 0.45) | 11.7 | 5.98 | 2.42 | 1.11 | 1.06 | 71 | 8679 | 3.73 | 0.00 | 1.0% | **C** |
| S12 | G6PD deficient + rasb | 14.8 | 5.35 | 1.91 | 1.17 | 1.15 | 115 | 8679 | 3.73 | 0.00 | 0.9% | L |

S9 is the only row whose `lysis at t<0` column is non-zero. That is the point of
a prephase — it moves lysis earlier and slower. This arm starts the same patient
at t = −120 h (the tumour having been wound back along its own growth curve), and
the benefit of the prephase lies in preventing the tumour from reaching 3×10¹² at
t = 0 in the first place.

S11 is interesting. The CKD host's Cr fold is only 1.06, yet it **is classified
as clinical TLS** — because the starting creatinine already exceeds 1.5×ULN. That
is, the CTLS definition fires structurally in the CKD patient, and this is a
property of the criterion rather than something the patient goes through.

---

## 16. Trial ledger

| Item | Model | Reported |
|---|---|---|
| Rasburicase 4-hour urate change | **−96%** | −86% (Goldman 2001, *Blood*) |
| Allopurinol 4-hour urate change | **+5%** | +2% (Goldman 2001) |
| Urate AUC₀₋₉₆ rasb/allo | **0.06×** | 0.39× — **model over-separates** |
| Urate held <8 mg/dL on rasburicase | yes | 87% response (Cortes 2010, *JCO*) |
| Urate held <8 mg/dL on allopurinol | no | 66% response (Cortes 2010) |
| Dialysis rate, high risk + full bundle | 1% | 1.5–5% overall (Coiffier 2008 panel) |
| Dialysis rate, high risk without prophylaxis | 8% (3×10¹²) / 99% (6×10¹²) | up to ~30% in historical Burkitt cohorts |
| Peak LDH, high-load Burkitt | 8679 U/L | commonly >2×ULN, often >5000 |
| Potassium shift with insulin 10 U | −0.60 mmol/L (2 h) | −0.6 to −1.0 (Allon 1990) |
| Urine pH under an NaHCO₃ load | 7.04 | alkalinisation-era target 6.5–7.5 |

**What the model fails to match is left in.** The urate AUC ratio over-separates
badly: 0.06× against a reported 0.39×. The most plausible cause is that the trial
population was of mixed risk, whereas this model's parameterisation is built
around the extreme high-load patient. It is recorded as a discrepancy rather than
erased by tuning.

---

## 17. Numerical hygiene

`A17` checks the mass balance and whether any state goes negative. Load 3×10¹² +
rasburicase:

```
purine released by lysis                  72.00 mmol
endogenous purine input                   62.50 mmol   (336 h)
terminal purine pool + crystal             5.13 mmol
cumulative urate oxidised (= H2O2)        56.79 mmol
potassium released by lysis              115.50 mmol
phosphate released by lysis              180.00 mmol
negative states: none
```

The purine balance closes through renal excretion and the HGPRT salvage pathway,
both of which are fluxes that do not accumulate as states.

**Two physical constraints** are written explicitly into the model.

1. **Crystal growth is capped by delivery.** The deposition rate cannot exceed
   0.95 times the tubular delivery rate of that solute. Without this cap,
   `C_urine = E/Q_urine` diverges as the urine flow rate goes to zero, and an
   obstructed kidney keeps depositing crystal it is no longer being supplied
   with.
2. **Tubular obstruction limits the urine flow rate independently of GFR.** This
   is what makes crystal nephropathy **oliguric** rather than merely azotaemic.

---

## 18. What is derived and what was calibrated

**From measured / literature values** — intracellular content (potassium
110 mmol/L of cell volume × 0.35 pL; nucleic acid 16 pg/cell at MW 330 → purine
24 mmol and nucleic acid phosphorus ~48 mmol per 10¹² cells), the uric acid
solubility law (pKa 5.75), saturable uric acid handling (basal FE 8%, excretion
700 mg/day), TmP/GFR, HPO₄²⁻ pKa₂ 6.80, all of the drug PK
(allopurinol/oxypurinol/febuxostat/rasburicase/venetoclax), the Km of urate
oxidase of 25 µmol/L, and the Ca×PO₄ threshold of 60 mg²/dL².

**Calibrated (not measured)** — the crystal nucleation rate constant, the
half-maximal crystal mass for obstruction (12 mmol), the interstitial calcium
phosphate deposition rate constant, the rates of nephron loss and recovery, and
all three hazard functions for arrhythmia, seizure, and dialysis. These were set
so that **an unprophylaxed high-load Burkitt patient reaches a peak Cr fold of
about 3 and recovers over two weeks.** **The absolute probabilities coming out of
the hazard functions must not be quoted** — they are only a shape that is flat at
baseline and steep in the pathological range.

**Not in the model** — spatial structure (crystal deposition is a single renal
lump, so the model cannot speak to a medullary/cortical distribution, and the
medullary concentrating factor of 1.4 is one lumped number standing in for the
entire axial gradient), coagulation, sepsis, and the nephrotoxicity of the drugs
themselves (methotrexate, aminoglycosides, contrast) — all of which occur
alongside TLS in real patients and would add to the injury driving force.

---

## 19. Summary — the ten claims this model makes

1. **TLS is a problem of two threshold concentrations.** UA_req (the
   concentration required to excrete the load) and UA_crit (the concentration at
   which precipitation begins). The disease is the region where the former
   exceeds the latter.
2. **UA_crit is a property of the prescription** — 6.8 at 2 L/day, 10.5 at
   4 L/m²/day, 33.2 mg/dL at pH 7.0. The Cairo–Bishop criterion of 8.0 mg/dL
   falls inside the range of this calculation.
3. **While GFR is intact, uric acid is the only solute that loses the race**
   (6.2 versus potassium 0.4 and phosphate 0.5). The potassium and phosphate
   abnormalities are secondary phenomena, and both come down in an arm given
   nothing but the urate enzyme.
4. **The optimal urine pH is not a constant** — 6.6 at an unprophylaxed 3×10¹²,
   7.0 at 6×10¹², and 5.9 (= no alkalinisation) if rasburicase is present. Erase
   the urate benefit and only the cost remains in the pH term.
5. **Rasburicase is a zero-order operator of fixed capacity.** 0.2 mg/kg =
   2.62 mmol/h = 3.2×10¹² cells. Once the load exceeds that, the 4-hour fall is
   still −72% while the peak urate becomes 37 mg/dL — what the registration
   endpoint measures and what sets the renal outcome are not the same thing.
6. **Phosphate binders are arithmetically useless in acute TLS.** Lysis sends
   12.1 mmol/h and diet sends 0.95 mmol/h, so the share an enteric binder reaches
   is 7%. The only route that takes phosphate out of the body is dialysis.
7. **Allopurinol's lead time matters far less than the hypothesis posed.** Twelve
   hours captures most of the benefit, and the reason is not the pool size (7% of
   the amount released) but oxypurinol accumulation (t½ 23 h).
8. **The loop comes back around to the drug as well.** Oxypurinol accumulates
   5-fold in AKI — the reason dose reduction becomes necessary in the patients at
   highest risk.
9. **Redistribution is not removal, and it is worse than expected.** Insulin
   gives −0.60 at 2 hours but is higher than untreated at 12 hours, and total
   body potassium rises **more** than untreated (+11 versus −8). Park potassium in
   the cells and the only route that was taking it out switches off.
10. **The loop is a knee, not a switch.** Even seeding 40 mmol of crystal comes
    back to the same eGFR by day 14. The risk band is narrow and steep, but it is
    not self-sustaining.

---

## ⚠️ Disclaimer

This is a qualitative / semi-quantitative QSP model for educational and research
purposes. It was built on the published literature and clinical trial data, but
it has not been independently verified or certified, and **must not be used
directly for actual clinical decision-making, prescribing, or regulatory
submission.** In particular, the event probabilities coming out of the hazard
functions and the calibrated parameters listed in §18 are illustrative shapes
only.
