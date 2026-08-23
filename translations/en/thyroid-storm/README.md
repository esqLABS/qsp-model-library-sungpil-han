# Thyroid storm (Thyroid Storm / Thyrotoxic Crisis)

> **The storm is a disease of loop gain, not of hormone concentration.**
> Total T4 and total T3 do not distinguish storm from plain thyrotoxicosis. What
> distinguishes it in this model is **whether the fast positive feedback loop is
> closed**, and that is written as a single heat-balance ratio — and that ratio is
> **a product of factors**.
>
> ```
>              Q_prod        80 W · M_thy · M_temp · M_sns · M_unc
>     Λ  =    ────────  =  ──────────────────────────────────────────────
>              Q_loss       (8 + 52·E) · Vol^1.5 · (1 + cool) · (Tc − 33)
> ```
>
> **Every drug used in thyroid storm multiplies one or two of these factors.**
> So the clinical order of treatment is derived as **arithmetic** rather than as a list.

<p align="center">
  <a href="ts_qsp_model.svg"><img src="ts_qsp_model.png" width="900" alt="Thyroid storm QSP mechanistic map"></a>
</p>

| Deliverable | File | Scale |
|---|---|---|
| Mechanistic map | [`ts_qsp_model.dot`](ts_qsp_model.dot) · [SVG](ts_qsp_model.svg) · [PNG](ts_qsp_model.png) | 142 nodes · 21 clusters · 216 edges |
| mrgsolve ODE model | [`ts_mrgsolve_model.R`](ts_mrgsolve_model.R) | **38 ODEs** · 18 scenarios · 5 diagnostic sweeps |
| Shiny app | [`ts_shiny_app.R`](../../../thyroid-storm/ts_shiny_app.R) | 11 tabs |
| References | [`ts_references.md`](ts_references.md) | 145 papers (PubMed links throughout) |
| Numerical verification | [`ts_verify_python.py`](ts_verify_python.py) → [`ts_verification_output.txt`](ts_verification_output.txt) | **The source of every number** in this README |

---

## 1. Why hormone concentration will not do

Observations from 1975–1980 are the starting point of this model: **total and free
T4/T3 overlap between storm patients and plain thyrotoxicosis patients**
(Brooks &amp; Waldstein 1980; Brooks 1975; Jacobs 1973). Severity then cannot be written
as a function of hormone concentration.

So this model writes severity as the state of **three positive feedback loops**.

| Loop | Route | Which drug cuts it |
|---|---|---|
| **A** | fT3 → β1 receptor density · lipolytic machinery → NEFA → displacement of T3 from TBG → free fraction → fT3 | **β blockade** (because lipolysis is β-mediated) |
| **B** | Tc → Q10 metabolic amplification → heat production → Tc | cooling, antipyresis (paracetamol) |
| **C** | Tc → hypothalamic damage → loss of heat-loss effectors → Tc | cooling, fluids, early temperature control |

When the loops close, **the stable operating point disappears.** That is the storm.

### The control experiment (the core evidence of the model)

**Left: raise the hormone only. No precipitant at all.**

| TRAb | Total T4 | Total T3 | Free T3 | Fd | Store | Tc | HR | BWPS | Λ | Verdict at 72 h |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|:--|
| 1.00 | 99.2 | 1.79 | 5.36 | 1.000 | 66.0 d | 37.00 | 70 | 0 | 0.928 | no storm |
| 1.50 | 116.3 | 2.11 | 6.38 | 1.006 | 55.6 d | 37.07 | 80 | 0 | 0.928 | no storm |
| 2.00 | 142.2 | 2.62 | 8.18 | 1.040 | 44.2 d | 37.16 | 92 | 0 | 0.929 | no storm |
| 3.00 | 204.9 | 3.94 | 14.17 | 1.198 | 27.7 d | 37.34 | 111 | 25 | 0.937 | no storm |
| **4.19** | **265.7** | **5.40** | **21.74** | **1.342** | **18.6 d** | **37.46** | **120** | **35** | **0.943** | **no storm** |
| 6.00 | 328.7 | 7.26 | 31.52 | 1.446 | 12.4 d | 37.55 | 126 | 40 | 0.949 | no storm |
| 9.00 | 403.0 | 9.16 | 41.37 | 1.506 | 7.9 d | 37.60 | 129 | 40 | 0.953 | no storm |
| 14.00 | 479.9 | 10.90 | 50.43 | 1.542 | 5.0 d | 37.64 | 130 | 45 | 0.955 | no storm |
| 22.00 | 544.4 | 12.37 | 58.00 | 1.563 | 3.1 d | 37.66 | 132 | 45 | 0.957 | no storm |
| 35.00 | 594.6 | 13.51 | 63.89 | 1.576 | 1.9 d | 37.68 | 132 | 45 | 0.958 | no storm |

**No hormone concentration produces a storm on its own.** There is a reason Λ saturates
at 0.958 and does not cross 1: the nuclear signal `Sig` saturates at 2, whereas the
heat-loss effector `h` can vary **7-fold**, from 8 to 60 W/K. That is, **the thyroid axis
alone is not strong enough to overwhelm human thermoregulation.** A second hit is
arithmetically **necessary.** The bold row is the index patient used throughout below
(total T3 = 5.40 nmol/L). That BWPS flattens out at 40–45 also means something — it sits
exactly in the range the clinic calls "impending storm" and never crosses over.

**Right: fix the hormone at total T3 = 5.40 nmol/L and vary only the precipitant.**

| Precipitant | Tmax | Tc@24h | HRmax | BWPS@24h | BWPSmax | Bilirubin | 7-day mortality | Verdict |
|---:|---:|---:|---:|---:|---:|---:|---:|:--|
| 0.00 | 37.46 | 37.46 | 120 | 45 | 45 | 0.97 | 0.04 % | survives |
| 0.30 | 37.84 | 37.78 | 132 | 55 | 60 | 0.98 | 0.09 % | survives |
| 0.60 | 38.21 | 38.10 | 141 | 60 | 65 | 1.05 | 0.11 % | survives |
| 0.90 | 38.59 | 38.43 | 148 | 75 | 75 | 1.38 | 0.19 % | survives |
| 1.10 | 38.84 | 38.77 | 152 | 85 | 90 | 1.71 | 0.30 % | survives |
| **1.20** | **39.00** | 38.99 | 154 | 90 | **100** | 2.24 | **0.47 %** | **survives** |
| **1.25** | **43.00** | 39.13 | 166 | 90 | 140 | 6.59 | **82.4 %** | **storm (runaway)** |
| **1.30** | **43.00** | 39.29 | 167 | 95 | 140 | 6.64 | **85.0 %** | **storm (runaway)** |
| 1.50 | 43.00 | 40.68 | 169 | 115 | 140 | 6.70 | 88.4 % | storm (runaway) |
| 2.00 | 43.00 | 43.00 | 171 | 120 | 140 | 6.74 | 90.9 % | storm (runaway) |

**There is a threshold between precipitant 1.20 and 1.25** — not a gradient. Below it the
patient is hot and tachycardic and **already in storm by BWPS (≥45)** but thermally stable,
and survives. Above it there is no reachable operating point and the temperature runs away.

> Total T4 and total T3 are **never once different** between the two tables above. One
> variable changed, and it is not a hormone.

---

## 2. The central treatment table — what each drug actually moves

Change at 24 hours starting from precipitant 1.30 (fulminant) and total T3 5.40 nmol/L
(against the untreated trajectory):

| Monotherapy (the axis it attacks) | Δtotal T3 | Δfree T3 | ΔHR | ΔTc | Runaway? | 7-day mortality |
|---|---:|---:|---:|---:|:--:|---:|
| PTU (synthesis + D1 conversion) | −13.8 % | −14.5 % | −1 | −0.14 | **yes** | **75.3 %** |
| Methimazole (synthesis only) | −1.8 % | −1.9 % | −0 | −0.01 | **yes** | **84.0 %** |
| Hydrocortisone (D1/D2 + steroid) | −8.6 % | −9.0 % | −1 | −0.09 | **yes** | **69.4 %** |
| **Potassium iodide** (release blockade) | **−24.1 %** | −25.2 % | −2 | −0.24 | no | **0.29 %** |
| **Propranolol** (delivery only) | **−3.5 %** | **−48.6 %** | **−58** | −0.89 | no | **0.05 %** |
| Esmolol (delivery only) | +0.0 % | −18.2 % | −20 | −0.57 | no | 0.11 % |
| Cooling + fluids (heat balance only) | +0.0 % | −2.9 % | −6 | −1.20 | no | 0.18 % |
| *(untreated reference values)* | *5.40* | *32.54* | *155* | *39.29* | *yes* | *85.0 %* |

**This table has to be read as follows.** The drugs that lower the measured hormone the
most (PTU, methimazole) do not stop the runaway and the patient dies. **The drug that
removes no hormone at all (a β blocker) ends the storm** — because it lowers total T3 by
3.5 % while **halving free T3 (−48.6 %).**

And **of the drugs that target the hormone, only iodine works.** The reason is structural:
iodine blocks **the exit from the store (release)**, whereas the thionamides block the
store's **entrance (synthesis)**. The store is already full.

### Two actions — β blockade cuts loop A twice

| t | Cpro (ng/mL) | β occupancy | NEFA | Fd | Total T3 | **Free T3** | HR | ‖ | Untreated NEFA | Fd | **Free T3** |
|---:|---:|---:|---:|---:|---:|---:|---:|:--|---:|---:|---:|
| 0.5 h | 24.7 | 39 % | 1.08 | 1.430 | 5.40 | 23.2 | 126 | ‖ | 1.37 | 1.624 | 26.3 |
| 2 h | 63.8 | 60 % | 0.74 | 1.160 | 5.39 | 18.8 | 116 | ‖ | 2.12 | 1.898 | 30.8 |
| 12 h | 89.0 | 66 % | 0.56 | 1.072 | 5.30 | 17.0 | 104 | ‖ | 2.16 | 1.960 | 31.8 |
| **24 h** | 100.2 | 68 % | **0.47** | **1.071** | **5.21** | **16.7** | **97** | ‖ | **2.13** | **2.009** | **32.5** |

Propranolol touches neither thyroid, liver nor kidney. And yet NEFA falls from 2.13 to
0.47 mmol/L, the displacement factor Fd falls from 2.009 to 1.071, and **free T3 halves
from 32.5 to 16.7 pmol/L.** Because lipolysis is β-mediated. This is the arithmetic of
loop A, and it is the central claim of this model.

### Factorisation of Λ (at 6 hours)

| Treatment | M_thy | M_temp | M_sns | M_unc | E | Vol^1.5 | 1+cool | h (W/K) | 7-day outcome |
|---|---:|---:|---:|---:|---:|---:|---:|---:|:--|
| Untreated | 1.706 | 1.161 | 1.495 | 1.000 | 0.705 | 0.894 | 1.00 | 39.9 | **runaway → death** |
| PTU alone | 1.693 | 1.160 | 1.494 | 1.000 | 0.700 | 0.894 | 1.00 | 39.7 | **runaway → death** |
| Methimazole alone | 1.706 | 1.161 | 1.495 | 1.000 | 0.705 | 0.894 | 1.00 | 39.9 | **runaway → death** |
| **Propranolol alone** | **1.516** | 1.127 | **1.121** | 1.000 | 0.447 | 0.914 | 1.00 | 28.6 | controlled |
| Cooling + fluids | 1.702 | 1.103 | 1.465 | 1.000 | 0.258 | 0.928 | **2.20** | 43.7 | controlled |
| Full bundle | **1.475** | 1.078 | **1.113** | 1.000 | 0.119 | 0.948 | **2.20** | 29.6 | controlled |
| Bundle + aspirin | **1.578**↑ | 1.088 | 1.115 | **1.099**↑ | 0.167 | 0.941 | 2.20 | 34.5 | controlled |

**Only β blockade lowers two factors at once** — `M_sns` directly, and `M_thy` through the
NEFA → displacement → free T3 route. **Aspirin raises two factors at once** — `M_unc`
(uncoupling of oxidative phosphorylation) and `M_thy` (displacement).

---

## 3. The store arithmetic — why a thionamide cannot be an acute drug

Methimazole alone. Organification is **89 %** blocked within an hour. Everything after
that is the arithmetic of **hormone already made**.

| t | TPO blockade | Synthesis | Store S | ΔS | Secretion | Total T4 | ΔT4 | Total T3 | ΔT3 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 h | 91.4 % | 1.44 | 6242 | −0.19 % | 13.93 | 265.7 | **−0.00 %** | 5.40 | **−0.00 %** |
| 6 h | 88.6 % | 2.65 | 6184 | −1.13 % | 13.44 | 265.5 | −0.05 % | 5.40 | −0.09 % |
| **24 h** | **89.3 %** | 2.52 | 6004 | **−4.00 %** | 10.72 | 262.7 | **−1.10 %** | 5.30 | **−1.83 %** |
| 72 h | 89.3 % | 1.32 | 5702 | −8.83 % | 5.78 | 239.8 | −9.73 % | 4.65 | −13.9 % |
| 168 h | 89.3 % | 0.73 | 5353 | −14.4 % | 3.64 | 183.9 | −30.8 % | 3.42 | −36.6 % |

At t = 0 the store is 6254 nmol and secretion 13.99 nmol/h → **18.6 days' worth of
hormone is already made.** Even a perfect synthesis blocker cannot touch it. At 24 hours
total T4 has fallen by only **1.1 %**.

> The only source of this number is calibration target N7, that **a normal person's
> thyroid store is 60–90 days' worth**. That it shrinks to 18.6 days' worth in florid
> Graves disease is **a result, not a calibration** (prediction P5). And even in a patient
> whose store has already been emptied 4-fold, a thionamide does nothing within 24 hours.

### The only reason PTU is better than methimazole

| t | Cptu | D1 blockade | Total T3 (PTU) | Δ | Total T3 (MMI) | Δ | **rT3 (PTU)** | **rT3 (MMI)** |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 3 h | 7.36 | 58.2 % | 5.27 | −2.5 % | 5.40 | −0.0 % | **0.928** | 0.745 |
| 12 h | 2.37 | 38.0 % | 4.97 | −7.9 % | 5.38 | −0.4 % | **1.168** | 0.745 |
| **24 h** | 2.31 | **37.5 %** | **4.65** | **−13.9 %** | 5.30 | **−1.8 %** | **1.267** | 0.740 |
| 72 h | 2.31 | 37.5 % | 3.67 | −32.0 % | 4.65 | −13.9 % | 1.255 | 0.682 |

Besides inhibiting TPO, PTU **also inhibits D1.** D1 makes 53 % of circulating T3, and
T3's half-life is 1 day (not T4's 6.3 days). So PTU reaches **the fast axis**. **rT3 is
the laboratory evidence for this** — because rT3 clearance depends chiefly on D1
(fD1_rT3 = 0.85), it rises from 0.74 to 1.27 nmol/L on PTU and does not move on
methimazole. This is a prediction testable at the bedside.

---

## 4. The second axis — the free fraction

**Holding** total T3 fixed at 5.40 nmol/L and varying NEFA alone:

| NEFA (mmol/L) | Fd | **Measured total T3** | Free T3 | Sig | **Equivalent total T3** |
|---:|---:|---:|---:|---:|---:|
| 0.40 | 1.000 | **5.40** | 16.2 | 1.500 | 5.40 |
| 0.80 | 1.197 | **5.40** | 19.4 | 1.564 | 6.46 |
| 1.50 | 1.688 | **5.40** | 27.3 | 1.670 | 9.11 |
| 2.00 | 1.860 | **5.40** | 30.1 | 1.696 | 10.05 |
| 3.00 | 2.000 | **5.40** | 32.4 | 1.714 | **10.80** |

**The "measured" column never moves once.** Two patients with the same test result can
differ 2-fold in the hormone the nucleus sees. This is the arithmetic reason total hormone
assays cannot diagnose storm, and it is also the reason **aspirin is contraindicated** in
storm (salicylate displaces T3 from TBG).

The displacement term was written as **a steep function with a Hill exponent of 2** —
fatty-acid displacement is negligible at everyday concentrations and becomes steep only at
the concentrations reached in critical illness, with heparin, and in storm (Lim 1988).

---

## 5. Aspirin — two independent harms, quantified

| Treatment | Tmax | Fd@24h | Free T3@24h | Q_prod@24h | BWPS@24h | BWPS@72h |
|---|---:|---:|---:|---:|---:|---:|
| Bundle | 38.20 | **1.063** | **11.9** | **122.8 W** | **25** | 20 |
| Bundle + **aspirin** | 38.28 | **1.468** | **16.5** | **153.3 W** | **40** | 25 |
| Bundle + paracetamol | 38.00 | 1.063 | 11.9 | 121.2 W | 25 | 20 |

Aspirin (1) **raises free T3 by 38 %**, from 11.9 to 16.5 pmol/L, through displacement,
and (2) **raises heat production by 25 %**, from 122.8 to 153.2 W, through uncoupling of
oxidative phosphorylation. Paracetamol does neither. Because the hyperthermia of storm is
not PGE₂-mediated fever but **an imbalance of heat production against heat loss**, the role
of an antipyretic is limited to begin with, and aspirin is wrong in two directions.

---

## 6. The trap of β blockade — a heart that depends on its rate

The same bundle, differing only in cardiac reserve at the moment the storm begins:

| Reserve CR(0) | Propranolol (t½ 4 h) | | | | Esmolol (t½ 9 min) | | | |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| | Perfusion nadir | Shock | BWPS@24h | Mortality | Perfusion nadir | Shock | BWPS@24h | Mortality |
| 1.00 | 0.606 | 0.135 | **25** | 0.02 % | 0.898 | 0.000 | 55 | 0.03 % |
| 0.60 | 0.517 | 0.262 | 35 | 0.07 % | 0.831 | 0.000 | 65 | 0.05 % |
| 0.45 | 0.439 | 0.373 | 50 | 0.13 % | 0.626 | 0.105 | 70 | 0.06 % |
| **0.35** | **0.374** | **0.466** | 50 | **0.20 %** | **0.487** | **0.304** | 70 | **0.08 %** |

**The trade-off is revealed quantitatively.** With sufficient reserve, propranolol is
overwhelmingly better (BWPS 25 against 55) — because being non-selective it also blocks β2
lipolysis and so cuts loop A harder. As the reserve falls, propranolol drops perfusion to
0.374 and the shock index reaches 0.466. Esmolol, with a half-life of 9 minutes, is safe
because it is **reversible**, and the price is that it cuts loop A less (BWPS 70).

---

## 7. Iodine — the sign of the same drug inverts with the thyroid

### An unexpected result (reported honestly)

**In the model the order of administration hardly mattered.** Omit the thionamide entirely
and the 24-hour total T3 worsens only from 3.74 to 4.13 nmol/L, about **10 %**. The reason
is arithmetic: **the very iodine that supplies the substrate simultaneously switches off
organification** (Wolff–Chaikoff acts on **both arms**, release and organification). A
thyroid whose release is inhibited also has its synthesis inhibited, so it cannot run away
on the substrate it has just received.

### When, then, is Jod-Basedow real?

**Only in a thyroid whose autoregulation has failed.** This is not a question of the order
of prescribing but of the **parameter region**. Iodine alone, no thionamide, 21 days:

| Thyroid | Total T4 d0 | d3 | d7 | **d21** | Store d21 | Synthesis d21 |
|---|---:|---:|---:|---:|---:|---:|
| Normal autoregulation (WCrel .78 / WCorg .90) | 265.7 | 207.8 | 155.3 | **85.9** ↓ | 6491 | 4.26 |
| Partial failure (WCrel .40 / WCorg .45) | 265.7 | 239.0 | 222.6 | **239.9** | 10394 | 19.74 |
| **Autonomous nodule / iodine deficiency** (WCrel .15 / WCorg .20) | 265.8 | 260.8 | 272.1 | **362.3** ↑↑ | **11681** | 28.66 |

In a thyroid with normal autoregulation, iodine brings T4 down from 265.7 to 85.9. In an
autonomous nodule it **raises it from 265.8 to 362.3.** The store swells from 6491 to 11681.
**The sign of iodine is set by the thyroid's autoregulation** — a falsifiable claim, and one
consistent with who actually experiences Jod-Basedow in the clinic.

---

## 8. What was calibrated, and therefore what is prediction

### Matched to normal physiology (20 items, section [N])

| | Calibration target | Target | Achieved |
|---|---|---|---|
| N1–N2 | Total T4 · total T3 | 100 nmol/L · 1.8 | 99.20 · 1.787 |
| N3–N4 | Free T4 · free T3 | 12–22 pmol/L · 3.5–6.5 | 19.84 · 5.361 |
| N5 | Reverse T3 | 0.15–0.45 nmol/L | 0.278 |
| N6 | Peripheral conversion's contribution to T3 | ~80 % | 80.6 % |
| **N7** | **Thyroid store** | **60–90 days' worth** | **66.0 days** |
| N8 | TSH | 1.5 mIU/L | 1.544 |
| N9–N10 | Core temperature · heart rate | 37.0 °C · 70 bpm | 37.00 · 69.9 |
| N11–N12 | NEFA · cortisol | 0.4 mmol/L · 400 nmol/L | 0.398 · 402 |
| N13 | Basal heat production = heat loss | 80 W | 79.65 = 79.65 |
| N14–N15 | T4 · T3 half-life | 6–7 d · 1.0 d | 6.28 · 1.006 |
| N16–N17 | Dietary iodine · T4 secretion | 150 µg/d · 80–110 µg/d | 145.9 · 85.1 |
| N18–N20 | Plasma inorganic iodide · h · cardiac load | 0.006–0.05 µmol/L · — · 1.0 | 0.0198 · 19.93 W/K · 1.09 |

**The part that is Euclidean-consistent:** normal heat balance (h = 20 W/K at 37.0 °C)
**follows algebraically** from BMR = 80 W, and dietary iodine is **determined** once the
iodine mass balance is closed. Neither is a calibration.

### Matched to non-storm drug data (7 items)

Thionamide TPO potency (methimazole roughly 10–20 times more potent per mg) · PTU's D1
potency (a 20–30 % fall in serum T3 at 24 hours in ordinary thyrotoxicosis) · propranolol's
β potency (80 mg giving a 25–30 % fall in heart rate) · iodine's release-inhibition Emax
(Lugol giving a 30–50 % fall in T4 at 24–48 hours in preoperative preparation) · the NIS
downregulation time constant · the D1/D2 potencies of glucocorticoids and iopanoic acid.

### Matched to storm data — **exactly one**

`h0_haz`. It was set so that the 7-day mortality of untreated fulminant disease is **85 %**
(the historical figure from before the 1970s). **Nothing else at all was adjusted against
storm data.**

### Therefore all of the following are predictions

| | Prediction | Model | Literature / clinic |
|---|---|---|---|
| **P1** | No hormone concentration produces a storm on its own | Λ saturates at 0.956 | a precipitant is present in 70–90 % of storms |
| **P2** | The boundary is a **threshold** in the precipitant | 0.47 % → 82.4 % between 1.20 and 1.25 | storm is a state, not a grade |
| **P3** | Thionamides do not stop a storm | total T4 −1.1 % at 24 h · mortality 84.0 % | 1–2 weeks to normalise the hormone |
| **P4** | β blockade gives total T3 −3.5 % and **free T3 −48.6 %** | 32.5 → 16.7 pmol/L | β blockade changes acute mortality the most |
| **P5** | The store in florid Graves disease is 18.6 days' worth | 66.0 d → 18.6 d | reduced thyroidal iodine in surgical specimens |
| **P6** | With adrenal output unchanged, cortisol ≈ 200 nmol/L | clearance 2.01-fold → 199 | relative adrenal insufficiency |
| **P7** | The sign of iodine is set by **thyroid autoregulation** | T4 85.9 against 362.3 (21 days) | Jod-Basedow in autonomous nodules |
| **P8** | At SSKI doses, iodine alone keeps working for 21 days | T4 265.7 → 85.9 | **disagrees** with the textbook escape at 10–14 days |
| **P9** | The propranolol trap occurs only below a reserve threshold | shock 0.466 at CR 0.35 | reports of β-blocker-induced circulatory collapse in storm |
| **P10** | Aspirin worsens things through **two factors** | free T3 +38 % · heat production +25 % | aspirin contraindicated in storm |

---

## 9. Honest limitations

| | Limitation |
|---|---|
| **L1** | **It underpredicts mortality in treated patients.** Modern reports give 8–25 % with full treatment, whereas this model's bundle is under 1 %. Because there are **no** comorbidities, multi-organ failure, thromboembolism, ventilator complications or deaths from drug toxicity (agranulocytosis, PTU liver failure) at all. Since **only the untreated figure** was calibrated, the treated numbers must be read as "the storm physiology was controlled" and not as survival estimates. |
| **L2** | The runaway stops at the numerical ceiling of Tc = 43.0 °C. The ceiling means "lethal hyperthermia" and is not a predicted measurement. |
| **L3** | The fast/slow separation is a modelling choice. β1 receptor density (36 h) and cardiac reserve were placed in the slow states, and the 72-hour verdict freezes them. |
| **L4** | The precipitant axis is **a single scalar**. Real precipitants differ qualitatively — an iodine load acts through the thyroid, DKA through volume, surgery through catecholamines. Those differences exist only on the map. |
| **L5** | Amiodarone-induced thyrotoxicosis (in which iodine is contraindicated) is only on the map and is not a simulated scenario. |
| **L6** | Because the Burch–Wartofsky score is computed from simulated signs, it inherits every simplification above. It is not an independent measurement. |
| **L7** | With no separation of plasma and extravascular compartments, plasmapheresis is represented only as **net systemic removal** (about 20 % of T4, about 5 % of T3). The large fall in plasma concentration immediately after exchange, and the subsequent rebound, are not visible. |
| **L8** | **Clinical Wolff–Chaikoff escape is not reproduced within 21 days.** The mechanism (NIS downregulation) is present and quantitatively large (a 7-fold fall in NIS, an 18-fold fall in intrathyroidal iodide). But the iodide pool does not fall below the Wolff–Chaikoff IC50, so the release blockade is not released. Either IC_WC is too low, or a second escape mechanism is missing. |
| **L9** | `Ithy` is a lumped, undifferentiated pool. Its **absolute value** under a pharmacological iodine load is not physiological. Only a saturating function of it enters the model, so the conclusions do not depend on the absolute value. |
| **L10** | The cortisol term has no ACTH feedback compensation. That basal cortisol is normal in real thyrotoxicosis is because ACTH compensates, and the model's P6 is a prediction under the condition "if adrenal output is unchanged". |

---

## 10. Using the files

```r
# the mrgsolve model
source("ts_mrgsolve_model.R")
res <- ts_run_all()        # the 18 scenarios
ts_table(res)              # the table in section 2 above
ts_plot(res)               # a quick plot
ts_sweep_hormone()         # the left-hand table of section 1 (hormone sweep)
ts_sweep_precipitant()     # the right-hand table of section 1 (precipitant sweep)
ts_sweep_reserve()         # section 6 (the β-blockade trap)
ts_wolff_chaikoff()        # section 7
ts_free_fraction()         # section 4
res <- ts_demo()           # all of it at once

# Shiny dashboard (11 tabs)
shiny::runApp("ts_shiny_app.R")
```

```bash
# numerical verification — regenerates every number in this README
python3 ts_verify_python.py      # → ts_verification_output.txt

# rendering the map
dot -Tsvg ts_qsp_model.dot -o ts_qsp_model.svg
dot -Tpng -Gdpi=150 ts_qsp_model.dot -o ts_qsp_model.png
```

---

## 11. Clinical summary — the treatment order the model derives

| Rank | Intervention | Factor it multiplies | Time to act | Why |
|---|---|---|---|---|
| 1 | **β blockade** (propranolol; esmolol if the reserve is low) | `M_sns` **and** `M_thy` | **hours** | Cuts loop A twice. Halves free T3. |
| 2 | **Cooling + fluids** | `1+cool`, `Vol^1.5`, `E` | **hours** | Lowers the gain of loops B and C directly. |
| 3 | **Iodine** (≥1 h after the thionamide) | `M_thy` (the store's **exit**) | **hours to a day** | The only hormone-targeting drug fast enough. |
| 4 | **Glucocorticoid** | `M_thy` (D1/D2) + steroid deficiency | hours to a day | Closes the cortisol scissors the thyroid hormone opened. |
| 5 | **Thionamide** (PTU in storm, because of D1) | `M_thy` (the store's **entrance**) | **weeks** | Prevents relapse. Does not end the storm. |
| 6 | Iopanoic acid · cholestyramine · plasmapheresis | `M_thy` | days to weeks | Adjunctive. Only plasmapheresis reaches the bound pool. |
| — | **Aspirin** | `M_unc`↑ **and** `M_thy`↑ | — | **Contraindicated.** Moves two factors in the wrong direction. |
| — | Definitive treatment (¹³¹I · thyroidectomy) | the store itself | after recovery | Removes the store. |

> **In one sentence:** drugs aimed at the store work on a timescale of **weeks**, and drugs
> aimed at the gain work on a timescale of **hours**. The storm is a disease of gain, so
> the gain has to be attacked first.

---

## ⚠️ Disclaimer

This model is a **quantitative systems pharmacology model for educational and research
purposes**. It was constructed from the public literature and normal physiology data but
has not been independently validated or certified, and
**must not be used for real clinical decision-making, prescribing, or regulatory
submission.** For clinical care, consult the ATA 2016 and JTA/JES 2016 guidelines directly.

Parent library → [../README.md](../README.md) · category: endocrine and metabolic
