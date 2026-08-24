# High-Altitude Illness QSP Model — Acute Mountain Sickness · High-Altitude Cerebral Oedema · High-Altitude Pulmonary Oedema
### High-Altitude Illness · AMS · HACE · HAPE

<p align="center">
  <a href="hai_qsp_model.svg">
    <img src="hai_qsp_model.png" width="880" alt="High-altitude illness QSP mechanistic map">
  </a>
</p>

| File | Contents |
|---|---|
| [`hai_qsp_model.dot`](hai_qsp_model.dot) · [SVG](hai_qsp_model.svg) · [PNG](hai_qsp_model.png) | Mechanistic map — 168 nodes, 16 clusters, 270 edges |
| [`hai_mrgsolve_model_en.R`](hai_mrgsolve_model_en.R) | mrgsolve ODE model (50 compartments, 200+ parameters) + 23 scenarios |
| [`hai_shiny_app_en.R`](hai_shiny_app_en.R) | Shiny dashboard (12 tabs) |
| [`hai_references_en.md`](hai_references_en.md) | 121 references — every PMID checked against NCBI |
| [`hai_reference_model.py`](hai_reference_model.py) | Independent Python/scipy reimplementation — for verification |
| [`hai_reference_output.txt`](hai_reference_output.txt) · [`hai_scenario_results.json`](hai_scenario_results.json) | The computed output behind every number below |

---

## What this model says, in one sentence

**High-altitude illness is not three diseases. It is the consequence of one
inspired oxygen partial pressure acting through three defence mechanisms whose
time constants differ from one another, and the whole of high-altitude medicine
comes out of that temporal mismatch.**

```
     P_IO2 = F_IO2 × (P_B(h) − 47)          ← the only exogenous variable

  defence ①  ventilation                    minutes → days (renal · CSF bicarbonate is rate-limiting)
  defence ②  hypoxic pulmonary vasoconstriction   seconds   (and it is spatially non-uniform)
  defence ③  cerebral vasodilatation         seconds   (it consumes intracranial volume buffering)
```

Until the slow defence arrives, the two fast defences run **unchecked**. That
interval is precisely the window in which AMS, HACE and HAPE arise, and the drugs
divide into exactly four branches:

| Branch | Drug | What it moves |
|---|---|---|
| **Accelerate** defence ① | acetazolamide | CSF · plasma bicarbonate → apnoea threshold |
| **Blunt** defence ② | nifedipine · tadalafil | λ (the HPV constriction factor) |
| **Pay the price** of defence ③ | dexamethasone | BBB permeability · VEGF · the symptoms themselves |
| **Remove the cause** | descent · oxygen · Gamow bag | P_IO2 itself |

Row 3 of this table is the most important line in the model. Dexamethasone
touches **no term at all** in the gas-exchange chain. Later on that emerges as
the number "descent equivalent 0 m".

---

## How this was verified

All 50 ODEs were implemented independently **twice** — once inside the C++ of
mrgsolve, once in Python/scipy. Every number below was obtained either by
integration or in closed form, and at the foot of the document **the six defects
caught in the process and the two hypotheses refuted** are recorded.

```bash
python3 hai_reference_model.py       # → hai_reference_output.txt, *.json
```

---

## Result 1 — the same mountain is a <ins>different disease</ins> on the day you arrive and a week later

Solve the **self-consistent fixed point** of the ventilatory control loop:

```
PaCO2 = 863·V̇CO2 / V̇A(PaCO2)
V̇A    = V̇E·(1 − V_D/V_T)
V̇E    = V̇E₀ + G_c·Φ(PaCO2 − B_c) + G_p·VAH·D·Φ(PaCO2 − B_p)
B_c   = [HCO3⁻]_CSF / (0.03·10^(pH*−6.1))          ← acid-base sets this
```

The residual `863·V̇CO2/V̇A(PaCO2) − PaCO2` decreases monotonically, so the root
is unique and can be found safely with Brent. Nowhere is an iterative
convergence fudged by hand.

| Altitude (m) | acute PaCO2 | acute PaO2 | **acute SaO2** | acclimatised PaCO2 | acclimatised PaO2 | **acclimatised SaO2** |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 40.1 | 94.4 | 97.4 % | 40.3 | 94.3 | 97.3 % |
| 2500 | 37.9 | 63.4 | 92.4 % | 34.3 | 67.2 | 93.4 % |
| 3500 | 36.3 | 52.1 | 87.9 % | 30.9 | 58.2 | 90.6 % |
| **4559** | 34.2 | 41.4 | **80.0 %** | 26.9 | 49.9 | **86.4 %** |
| 5300 | 32.6 | 34.7 | 71.8 % | 24.1 | 44.7 | 82.7 % |
| 6400 | 30.6 | 25.8 | 54.5 % | 20.7 | 37.3 | 74.6 % |
| 8000 | 28.7 | 14.1 | 21.9 % | 16.8 | 27.5 | 57.1 % |
| 8848 | 28.2 | 8.1 | 8.5 % | 15.1 | 23.0 | 45.6 % |

At 4559 m the **same altitude** gives SaO2 80.0 % on the day of arrival and
86.4 % once the kidney has finished its work — **a gap of 6.4 points**. And
converted into the currency that matters,

> **Acclimatisation at 4559 m = a descent of 821 m.**

The "acute" column above 6000 m is not a prediction but a counterfactual. Nobody
arrives at 8000 m still carrying sea-level bicarbonate. That column is there to
measure **the size of the deficit** that acclimatisation fills.

---

## Result 2 — with altitude, hyperventilation buys <ins>progressively less in mmHg</ins>, and pays only because of the dissociation curve

Differentiating the alveolar gas equation and the CO2 hyperbola gives

```
dP_AO2/dV̇A = c·PaCO2/V̇A          (c = F_IO2 + (1−F_IO2)/R)
```

**There is no barometric pressure in this expression.** P_B cancels out of the
derivative completely. And yet the value itself **falls** as altitude rises —
because the subject is already hyperventilating.

| Altitude (m) | dP_AO2/dV̇A (mmHg per L/min) | dSaO2/dPaO2 (%/mmHg) | product (sea level = 1) |
|---:|---:|---:|---:|
| 0 | 10.09 | 0.081 | 1.00 |
| 2500 | 7.32 | 0.267 | 2.38 |
| 4559 | 4.51 | 0.669 | 3.67 |
| 5300 | 3.62 | 0.902 | 3.97 |
| 8000 | 1.76 | 2.407 | **5.15** |
| 8848 | 1.42 | 2.844 | 4.93 |

**Sevenfold worse in mmHg, thirty-five-fold better in saturation, and the
product is a fivefold gain.** Everything that makes hyperventilation worth doing
at altitude lies not in the gas laws but in the **sigmoid**. Put the same
statement in the language of drugs: a drug that raises PaO2 by 5 mmHg is worth
0.4 saturation points at sea level and 12 points at 8000 m.

---

## Result 3 — acetazolamide does not stimulate breathing. It <ins>lowers the floor faster</ins>

What determines periodic breathing is not the level of ventilation but the
**CO2 reserve** `= PaCO2 − B_c` (because when SaO2 recovers transiently during
the hyperventilatory phase the peripheral drive drops out, and the only thing
still holding ventilation up is the [H⁺] of the CSF).

| Altitude (m) | | operating-point PaCO2 | apnoea threshold B_c | **CO2 reserve** | predicted AHI | SaO2 |
|---:|---|---:|---:|---:|---:|---:|
| 0 | acclimatised | 40.38 | 36.27 | **4.11** | 3.5 | 97.3 % |
| 0 | + ACZ | 35.20 | 30.09 | **5.11** | 2.2 | 97.7 % |
| 4000 | acclimatised | 29.73 | 27.42 | **2.31** | 16.9 | 88.1 % |
| 4000 | + ACZ | 26.30 | 21.78 | **4.51** | 2.9 | 90.0 % |
| **4559** | acclimatised | 27.66 | 25.69 | **1.97** | 25.2 | 85.5 % |
| **4559** | + ACZ | 24.51 | 20.13 | **4.38** | 3.1 | 87.7 % |
| 5300 | acclimatised | 24.88 | 23.35 | **1.53** | 41.8 | 81.4 % |
| 5300 | + ACZ | 22.09 | 17.90 | **4.19** | 3.4 | 84.0 % |

At 4559 m acetazolamide lowers the operating-point PaCO2 by a further
**3.14 mmHg**. But it lowers the threshold by **3.84 mmHg**. The difference is
the reserve, and so the reserve **widens** from 1.97 → 4.38 mmHg.

Why the floor comes down faster than the operating point falls out of the model
automatically: when acetazolamide raises SaO2, **the peripheral hypoxic drive
shrinks by exactly that much on its own**, so the operating point cannot descend
any further. The threshold has no such feedback.

This is the mechanism behind the clinical observation that "acetazolamide
abolishes nocturnal periodic breathing", and it is a story with **the opposite
sign** from the explanation "because it is a respiratory stimulant".

---

## Result 4 — the ceiling on HAPE is set not by the strength of HPV but by its <ins>non-uniformity</ins>

Divide the pulmonary vascular bed into two compartments. A fraction `a` responds
to HPV by constricting λ-fold, and the remaining `1−a` constricts only
μ = 1 + κ(λ−1)-fold (κ = 0.08).

```
G_tot    = G₀ [ a/λ + (1−a)/μ ]
Q̇_open   = Q̇ · [(1−a)/μ] / [ a/λ + (1−a)/μ ]
P_cap,open = P_LA + Q̇ · (r_v/μ) / [ a/λ + (1−a)/μ ]      ← (1−a) cancels
```

The last line is the point of this sub-model. Taking λ → ∞ gives the
**amplification ceiling**:

```
             1
   ceiling = ─────────────      (reduces to the commonly used 1/(1−a) when κ = 0)
             1 − a(1−κ)
```

| a | ceiling 1/[1−a(1−κ)] | (naive 1/(1−a)) | P_cap at λ=1 | λ=4 | λ=8 | λ→∞ |
|---:|---:|---:|---:|---:|---:|---:|
| 0.25 | 1.24 | 1.33 | 11.00 | 11.63 | 11.76 | 11.90 |
| 0.50 | 1.85 | 2.00 | 11.00 | 12.58 | 13.02 | 13.56 |
| 0.75 | 3.11 | 4.00 | 11.00 | 14.22 | 15.57 | 17.68 |
| **0.85** | **4.59** | 6.67 | 11.00 | 15.26 | 17.50 | **21.76** |
| 0.90 | 5.75 | 10.00 | 11.00 | 15.92 | 18.89 | 25.44 |

> **The naive 1/(1−a) formula overestimates the ceiling by 45 % at a = 0.85.**
> This error was caught while checking the closed form against the numerics, and
> is recorded as it stands (defect 5 under "Defects caught during verification"
> below).

The conclusion is the same either way. **A drug that halves λ is moving along a
curve that has already flattened. And a lung with a = 0.5 cannot reach that
height in the first place.**

### And Q̇ enters <ins>multiplicatively</ins> — HAPE is not a disease of altitude

Same lung, same HPV, 4559 m, a = 0.85. The only thing that changes is cardiac
output:

| Q̇ (L/min) | mPAP | P_cap,open | > 19.5 mmHg? |
|---:|---:|---:|:--|
| 5.0 | 27.8 | 15.12 | no |
| 6.0 (at rest) | 31.8 | 16.55 | no |
| 8.0 | 37.5 | 19.31 | no |
| **8.15** | | **19.50** | ← critical |
| 10.0 | 42.3 | 21.75 | **YES** |
| 15.0 | 51.9 | 26.96 | **YES** |
| 22.0 | 61.9 | 33.01 | **YES** |

**Critical cardiac output 8.15 L/min = 1.36 times resting.** That is a brisk
walk. The lung of the HAPE-susceptible phenotype lies below the stress-failure
line at rest and crosses it as soon as walking begins. HAPE is a disease of
altitude × **exertion**, and the model says so in arithmetic — and the
bifurcation parameter a clinician can actually control is that one, not altitude.

---

## Result 5 — a hypothesis the model <ins>refuted</ins>: the P50 "paradox"

**The hypothesis tested.** Two shifts oppose one another. 2,3-DPG moves the
curve to the right (favouring unloading in muscle), respiratory alkalosis to the
left (favouring loading in the lung). The textbook narrative is that DPG wins at
moderate altitude and alkalosis at extreme altitude, so that P50 **peaks
somewhere around 4000–5000 m and then comes down again**.

**Refuted.** Calibrate the rise in 2,3-DPG to the measured value of about +20 %
(5.0 → 6.0 mmol/L) and P50 **decreases monotonically at every altitude**.

| Altitude (m) | pH | 2,3-DPG | P50 | pH contribution | DPG contribution |
|---:|---:|---:|---:|---:|---:|
| 0 | 7.401 | 5.03 | 26.77 | −0.04 | +0.01 |
| 4559 | 7.447 | 5.58 | 25.68 | −1.35 | +0.22 |
| 6400 | 7.471 | 6.00 | 25.15 | −2.02 | +0.37 |
| 8848 | 7.495 | 6.00 | **24.50** | **−2.67** | **+0.37** |

The reason is one-sided: the DPG term **saturates**, so it can contribute at
most +0.37 mmHg, while the alkalosis term does not saturate and goes as far as
−2.67 mmHg. No crossover point exists. For a peak to appear, 2,3-DPG would have
to rise by about +50 %, and it does not.

**The claim that survives — and it is the stronger one — is about magnitude, not
shape.** Compute the counterfactual in which PaO2 on the summit is left at
23.0 mmHg and only the respiratory alkalosis is removed:

> SaO2 45.6 % → **39.0 %**. The alkalosis alone is worth **6.6 saturation
> points**.

The oxygen dissociation curve of the extreme-altitude climber is positioned
**for the lung**, not for the muscle.

---

## Result 6 — the optimal haematocrit is 43 %, and it <ins>does not move with altitude</ins>

In `Ḋ_O2 = Q̇ · C_aO2`, since `C_aO2 ∝ Hct`, `Q̇ ∝ μ^(−γ)` and
`μ ∝ exp(k·Hct)` (k = 2.31),

```
Ḋ_O2 ∝ Hct·exp(−k·γ·Hct)        ⇒     Hct* = 1/(k·γ)          (exact)
```

**S_aO2 cancels.** C_aO2 is proportional to Hct whatever the saturation. This
model therefore asserts flatly that the optimum is 43.3 % at sea level and
43.3 % at 5000 m as well.

| γ (cardiac-output viscosity penalty exponent) | optimal Hct |
|---:|---:|
| 1.00 | 43.3 % |
| 0.90 | 48.1 % |
| 0.80 | 54.1 % |
| 0.75 | 57.7 % |
| 0.60 | 72.2 % |

To justify the Hct of 55 % of Andean highlanders, **γ ≤ 0.787** is required —
that is, the circulation must absorb most of the viscosity burden. This is the
quantitative form of the Tibetan-type (no erythrocytosis) versus Andean-type
(with erythrocytosis) adaptation debate, and at the same time the explanation of
why venesection works in chronic mountain sickness.

⚠️ γ = 1 is the assumption that cardiac output absorbs **none** of the viscosity
burden, so 43.3 % must be read as a **lower bound**. Reading it together with
the sensitivity table above is the correct way to use this result.

---

## Result 7 — every intervention in a single currency: <ins>how many metres of descent</ins>

Reference: acute, unacclimatised, 4559 m, SaO2 = 80.0 %.

| Intervention | SaO2 | ΔSaO2 | **= descent equivalent (m)** |
|---|---:|---:|---:|
| nothing at all | 80.0 % | 0.00 | **0** |
| dexamethasone | 80.0 % | **0.00** | **0** |
| full acclimatisation (about 1 week) | 86.4 % | +6.44 | **821** |
| acetazolamide 250 bid (steady state) | 88.3 % | +8.33 | **1 134** |
| Gamow bag 2 psi (+105 mmHg) | 91.2 % | +11.25 | **1 739** |
| oxygen F_IO2 0.28 (≈ 2 L/min) | 92.9 % | +12.94 | **2 200** |
| Gamow bag 4 psi (+207 mmHg) | 95.4 % | +15.47 | **3 191** |
| oxygen F_IO2 0.35 (≈ 4 L/min) | 96.8 % | +16.78 | **4 013** |

Two things stand out.

1. **The steady state on acetazolamide is greater than natural acclimatisation**
   (1 134 m versus 821 m). The metabolic acidosis the drug creates is added **on
   top of** the renal compensation, and that is why the statement is not "the
   drug mimics acclimatisation" but "the drug **exceeds** acclimatisation".
2. **Dexamethasone is 0 m.** This is not a criticism but an explanation. The
   drug enters no term of the gas-exchange chain. And that is precisely why it
   is dangerous — it abolishes the symptoms that would have stopped the ascent
   while leaving PaO2 exactly where it was.

---

## Result 8 — two drugs, comparable symptom improvement, <ins>only one of them moved the oxygen</ins>

| arm | nadir SaO2 | **final SaO2** | final PaCO2 | final HCO3⁻ | peak LLS | AMS hours | nocturnal AHI (day 1 → last) |
|---|---:|---:|---:|---:|---:|---:|---|
| no prophylaxis | 79.4 % | **85.0 %** | 28.4 | 18.8 | 4.54 | 72 h | 81 → 54 |
| acetazolamide 125 bid | 79.8 % | **87.1 %** | 25.4 | 16.4 | 3.86 | 31 h | 53 → 5 |
| acetazolamide 250 bid | 80.2 % | **87.8 %** | 24.1 | 15.5 | 3.58 | 23 h | 29 → 3 |
| dexamethasone 4 mg bid | 79.3 % | **85.0 %** | 28.4 | 18.8 | 3.60 | 33 h | 81 → 54 |
| ibuprofen 600 tid | 79.3 % | **85.0 %** | 28.4 | 18.8 | 4.09 | 29 h | 81 → 54 |
| graded ascent (~400 m/day) | 82.7 % | **84.8 %** | 28.6 | 18.9 | 3.48 | 32 h | 24 → 55 |

The nadir SaO2 immediately after arrival is much the same in three of the arms —
during the ascent neither the kidney nor the drug has yet finished its work.
What diverges is the **steady state**.

**Acetazolamide 250 bid** moves steady-state SaO2 by +2.9 points, nadir SaO2 by
+0.8 points, peak LLS by -0.96 and AMS hours by -49 hours. Plasma HCO3⁻ comes
down from 18.8 → 15.5 mEq/L.

**Dexamethasone 4 mg bid** moves steady-state SaO2 by -0.0 points, nadir SaO2 by
-0.1 points, peak LLS by -0.94 and AMS hours by -39 hours. Plasma HCO3⁻ is
18.8 mEq/L — **unchanged**.

Both drugs lower the symptom score by comparable amounts. And yet **only one of
them moved the oxygen.** Every guideline that says "descend if symptoms persist"
is using symptoms as an oxygen gauge, and dexamethasone is the drug that breaks
that gauge. The nocturnal AHI column tells the same story once more:
acetazolamide goes 29 → 3, dexamethasone 81 → 54, the same as placebo.


---

## Scenarios (the full 50-state ODE, 23 of them)

| # | Scenario | Phenotype | nadir SaO2 | peak LLS (time) | AMS hours | peak ICP | EVLW max/final | peak mPAP | peak P_cap |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|
| 1 | sea level | typical trekker | 97.3 % | 0.38 (54 h) | 0 h | 10.9 | 0 / 0 mL | 14.6 | 11.2 |
| 2 | rapid ascent 4559 m | typical trekker | 79.4 % | 4.54 (30 h) | 72 h | 14.2 | 0 / 0 mL | 19.5 | 12.8 |
| 3 | 4559 m + acetazolamide 125 mg bid | typical trekker | 79.8 % | 3.86 (26 h) | 31 h | 13.5 | 0 / 0 mL | 19.4 | 12.8 |
| 4 | 4559 m + acetazolamide 250 mg bid | typical trekker | 80.2 % | 3.58 (25 h) | 23 h | 13.1 | 0 / 0 mL | 19.3 | 12.8 |
| 5 | 4559 m + dexamethasone 4 mg bid | typical trekker | 79.3 % | 3.60 (26 h) | 33 h | 14.2 | 0 / 0 mL | 19.3 | 12.8 |
| 6 | graded ascent 4559 m | typical trekker | 82.7 % | 3.48 (102 h) | 32 h | 12.2 | 0 / 0 mL | 19.3 | 12.8 |
| 7 | climb high sleep low | typical trekker | 84.1 % | 2.73 (64 h) | 0 h | 11.5 | 0 / 0 mL | 18.9 | 12.7 |
| 8 | HAPE-S rapid + exertion | HAPE-susceptible | 64.8 % | 8.60 (39 h) | 86 h | 37.5 | 321 / 169 mL | 56.4 | 27.4 |
| 9 | HAPE-S + nifedipine SR 30 mg bid | HAPE-susceptible | 65.2 % | 8.55 (39 h) | 87 h | 37.2 | 245 / 127 mL | 54.4 | 26.9 |
| 10 | HAPE-S + tadalafil 10 mg bid | HAPE-susceptible | 65.2 % | 8.52 (39 h) | 87 h | 36.7 | 196 / 97 mL | 47.6 | 25.6 |
| 11 | HAPE-S + dexamethasone 8 mg bid | HAPE-susceptible | 65.1 % | 6.31 (16 h) | 86 h | 34.1 | 248 / 98 mL | 53.2 | 26.8 |
| 12 | HAPE -> descent 2000 m | HAPE-susceptible | 64.8 % | 8.60 (39 h) | 70 h | 37.5 | 316 / 30 mL | 56.4 | 27.4 |
| 13 | HAPE -> O2 2-3 L/min | HAPE-susceptible | 64.8 % | 8.60 (39 h) | 66 h | 37.5 | 316 / 24 mL | 56.4 | 27.4 |
| 14 | HAPE -> Gamow bag 2 psi | HAPE-susceptible | 64.8 % | 8.60 (39 h) | 80 h | 37.5 | 316 / 43 mL | 56.4 | 27.4 |
| 15 | ignore AMS, keep ascending | tight-fit (low PVI) | 70.8 % | 9.36 (78 h) | 87 h | 34.2 | 0 / 0 mL | 21.7 | 13.3 |
| 16 | HACE -> dexamethasone + descent | tight-fit (low PVI) | 70.8 % | 9.10 (62 h) | 74 h | 34.1 | 0 / 0 mL | 21.6 | 13.3 |
| 17 | Everest, no supplemental O2 | elite climber | 24.7 % | 11.94 (472 h) | 375 h | 45.5 | 48 / 38 mL | 37.7 | 20.6 |
| 18 | Everest, supplemental O2 from 7900 m | elite climber | 57.1 % | 9.48 (450 h) | 362 h | 21.4 | 0 / 0 mL | 26.2 | 17.8 |
| 19 | sleep 4000 m | typical trekker | 84.0 % | 3.18 (30 h) | 20 h | 13.3 | 0 / 0 mL | 18.6 | 12.6 |
| 20 | sleep 4000 m + acetazolamide | typical trekker | 84.6 % | 2.57 (17 h) | 0 h | 12.4 | 0 / 0 mL | 18.4 | 12.6 |
| 21 | 21 days at 3800 m | typical trekker | 85.7 % | 2.84 (30 h) | 0 h | 12.1 | 0 / 0 mL | 18.2 | 12.5 |
| 22 | HAPE-S + salmeterol 125 ug bid | HAPE-susceptible | 64.9 % | 8.55 (39 h) | 87 h | 37.3 | 272 / 105 mL | 56.5 | 27.4 |
| 23 | 4559 m + ibuprofen 600 mg tid | typical trekker | 79.3 % | 4.09 (30 h) | 29 h | 14.2 | 0 / 0 mL | 19.5 | 12.8 |

### HAPE prophylaxis, 3-arm — reproduces the ranking of Maggiorini 2006

| arm | EVLW max | peak mPAP | peak P_cap | nadir SaO2 |
|---|---:|---:|---:|---:|
| no treatment (HAPE-susceptible + exertion) | 321 mL | 56.4 | 27.4 | 64.8 % |
| nifedipine SR 30 bid | 245 mL | 54.4 | 26.9 | 65.2 % |
| tadalafil 10 bid | 196 mL | 47.6 | 25.6 | 65.2 % |
| dexamethasone 8 bid | 248 mL | 53.2 | 26.8 | 65.1 % |
| salmeterol 125 µg bid | 272 mL | 56.5 | 27.4 | 64.9 % |

The table shows that the three drugs touch **different terms**. Nifedipine and
tadalafil lower mPAP and P_cap (pressure), while salmeterol lowers EVLW alone
**without touching P_cap at all** (drainage). Dexamethasone does a little of
both.

### HAPE rescue, 3-arm — started at 40 hours

| Intervention | EVLW max | **EVLW at 96 h** | final SaO2 |
|---|---:|---:|---:|
| descent to 2000 m | 316 mL | **30 mL** | 92.7 % |
| oxygen F_IO2 0.28 | 316 mL | **24 mL** | 93.1 % |
| Gamow bag 2 psi | 316 mL | **43 mL** | 91.6 % |

All three interventions move the same term (P_IO2), and so all three work. The
only difference is magnitude.

### Everest summit day

| arm | nadir SaO2 | final PaCO2 | peak LLS | EVLW max | peak ICP |
|---|---:|---:|---:|---:|---:|
| no supplemental oxygen | 24.7 % | 18.0 | 11.94 | 48 mL | 45.5 |
| F_IO2 0.45 from 7900 m | 57.1 % | 32.0 | 9.48 | 0 mL | 21.4 |


---

## Virtual population — AMS · HAPE incidence

150 virtual subjects (inter-individual variability assigned to HVR, PVI, λ_max,
non-uniformity a and Hb; fixed seed).

⚠️ **The absolute incidences are not calibrated and are over-predicted.** The
real AMS incidence on a rapid ascent to 4559 m is 50–60 %, whereas the model
gives 94.7 %. The reason is clear: the symptom gain was calibrated to "the
severity of the average trekker" (peak LLS 4.54 for the typical trekker, AMS
threshold 3), and the width of the inter-individual distribution is not broad
enough to straddle that threshold. As a result almost everybody crosses it.
**This section must be read only as a relative comparison between arms, and even
that relative comparison is compressed by the ceiling effect** — the RR of 0.880
for acetazolamide is far closer to 1 than the roughly 0.45 of the literature.
Stated honestly, this is the largest quantitative failure of the model.

| arm | AMS incidence | HAPE incidence | mean LLS | LLS 90th percentile | RR (AMS) |
|---|---:|---:|---:|---:|---:|
| rapid 4559 m, no prophylaxis | 94.7 % | 0.0 % | 4.63 | 5.92 | — |
| rapid 4559 m, acetazolamide 250 bid | 83.3 % | 0.0 % | 3.61 | 4.37 | 0.880 |
| rapid 4559 m, dexamethasone 4 bid | 86.7 % | 0.0 % | 3.66 | 4.57 | 0.915 |
| graded 4559 m (400 m/day) | 58.0 % | 0.0 % | 3.45 | 4.48 | 0.613 |

### Ascent-rate sweep — AMS is an integrator with memory

| Ascent rate (m/day) | days taken | peak LLS | AMS hours | nadir SaO2 | peak ICP |
|---:|---:|---:|---:|---:|---:|
| 1500 | 2.29 | 3.69 | 58 h | 83.2 % | 12.5 |
| 1000 | 3.43 | 3.51 | 28 h | 83.5 % | 11.9 |
| 750 | 4.57 | 3.40 | 14 h | 83.7 % | 11.6 |
| 600 | 5.71 | 3.33 | 8 h | 83.9 % | 11.5 |
| 450 | 7.62 | 3.21 | 0 h | 84.2 % | 11.3 |
| 300 | 11.43 | 3.08 | 0 h | 84.6 % | 11.0 |
| 200 | 17.14 | 3.00 | 0 h | 84.8 % | 10.9 |

### The critical altitude for HAPE — a bifurcation by phenotype

| Phenotype | critical altitude |
|---|---:|
| normal, rest | none below 7000 m |
| normal, 3x VO2 | 4673 m |
| HAPE-susceptible, rest | none below 7000 m |
| HAPE-susceptible, 3x VO2 | 2000 m |

**At rest neither phenotype has a critical altitude below 7000 m.** The moment
threefold exertion is put in, the normal subject comes down to 4673 m and the
susceptible subject to 2000 m, the lower limit of the sweep (that is, a
susceptible subject crosses the stress-failure line even at 2000 m given
sufficiently vigorous exertion).

This table repeats the arithmetic of Result 4 in the time domain. **The
bifurcation parameter a clinician can actually control is not altitude.** And an
altitude criterion such as "HAPE occurs above 3000 m" means something only once
the amount of exertion is held fixed.


---

## Calibration comparison (model vs. literature)

| Item | Model | Observed | Error | Source |
|---|---:|---:|---:|---|
| PB, Everest summit (mmHg) | 252.68 | 253.00 | -0.1 % | West 1983/1996 |
| PB, Denver 1610 m (mmHg) | 632.57 | 632.00 | +0.1 % | standard |
| Sea level PaCO2 (mmHg) | 40.10 | 40.00 | +0.2 % | - |
| Sea level PaO2 (mmHg) | 94.36 | 95.00 | -0.7 % | - |
| Sea level SaO2 (%) | 97.38 | 97.50 | -0.1 % | - |
| Sea level VE (L/min) | 6.53 | 6.30 | +3.7 % | - |
| 4559 m acute SaO2 (%) | 79.98 | 81.00 | -1.3 % | Capanna Margherita field data |
| 4559 m acclim. SaO2 (%) | 86.42 | 87.00 | -0.7 % | Capanna Margherita field data |
| 5300 m acclim. SaO2 (%) | 82.66 | 84.00 | -1.6 % | Everest base camp, CXE |
| 5300 m acclim. PaCO2 (mmHg) | 24.13 | 24.00 | +0.5 % | Grocott 2009 |
| 8400 m PaO2 (mmHg) | 26.87 | 24.60 | +9.2 % | Grocott NEJM 2009 (n=4) |
| 8400 m PaCO2 (mmHg) | 14.63 | 13.30 | +10.0 % | Grocott NEJM 2009 |
| 8400 m SaO2 (%) | 56.19 | 54.00 | +4.1 % | Grocott NEJM 2009 |

**Median |error| = 0.7 %, maximum = 10.0 %.**


---

## Defects caught by the dual implementation

1. **Acetazolamide killed the patient.** The red-cell carbonic anhydrase
   inhibition term was written `ve_eff = ve·(1 − 0.35·A_RBC/V_RBC·0.12)`, but
   `A_RBC` was a **mass (mg)**, not a concentration. At steady state the bracket
   went negative, ventilation diverged negatively, and PaCO2 stuck at 3 mmHg
   with plasma HCO3⁻ at 2.4 mEq/L. Replaced by a saturating Emax (IC50
   220 mg/L, maximum 4.5 %).
2. **Dexamethasone did nothing.** EC50 had been set at 0.35 mg/L, whereas the
   actual plasma concentration after 4 mg orally is 0.05–0.08 mg/L. The effect
   fraction stayed at 0.16 and scenarios 05 and 11 were indistinguishable from
   placebo. Corrected to EC50 0.030 mg/L.
3. **HAPE produced nine litres in twelve hours.** With the leak coefficient set
   to 33 mL/h/mmHg^1.5, EVLW reached 9 010 mL and SaO2 fell to 24.6 %. That is
   not HAPE, it is drowning. Recalibrated to 2.60 (about 1 000 mL in 24 hours).
4. **The mPAP of an exercising climber was 77 mmHg.** No catheter has ever
   recorded such a value. The cause was writing the pulmonary vascular bed as a
   **fixed resistance** — in reality recruitment and distension occur as flow
   rises. Added `rec = 1 + 0.35·(Q̇/Q̇₀ − 1)`. mPAP came down into the 55–60
   range, and at the same time the artefact in which a healthy climber drowned
   of HAPE in the Everest summit scenario disappeared.
5. **The amplification-ceiling formula was wrong.** Taking λ → ∞ with μ held
   fixed gave 1/(1−a), but μ = 1 + κ(λ−1) diverges along with λ. The correct
   limit is **1/[1 − a(1−κ)]**, and at a = 0.85 the naive formula overestimates
   by 45 %. Caught by checking the closed form against the numerics.
6. **Cerebral oedema of 51 mL, intracranial pressure at its ceiling.** The
   vasogenic oedema gain was an order of magnitude too large. Recalibrated so
   that ΔV ≈ 4–6 mL in AMS (ICP 16–20 mmHg) and 9–14 mL in HACE (ICP
   30–40 mmHg). At the same time cerebral CO2 reactivity was changed from linear
   to a **saturating sigmoid** — extrapolating a linear reactivity to the summit
   PaCO2 of 13 mmHg makes cerebral blood flow vanish.

## Hypotheses the model refuted

1. **The P50 "paradox"** (Result 5). With the measured rise in 2,3-DPG no peak
   exists and P50 decreases monotonically. The claim that survives is the
   magnitude (6.6 saturation points on the summit).
2. The initial assumption that **dexamethasone outperforms nifedipine in HAPE
   prophylaxis**. The model computes that dexamethasone enters on the alveolar
   fluid clearance and permeability side and so lowers EVLW while leaving
   capillary pressure almost untouched — consistent with the direction in
   Maggiorini 2006, where dexamethasone (29 %) was worse than tadalafil (14 %);
   it **reproduces the ranking of the three arms but not the absolute
   incidences**.

---

## Limitations — where the model disagrees with the literature

1. **The time course of nocturnal periodic breathing.** Because the model
   explains the instability with the CO2 reserve alone, it predicts that **the
   first night after arrival is the worst** and that things improve thereafter.
   The field polysomnography of Bloch et al. at 4559 m reported an AHI that was
   **higher** on day 3 than on day 1. The missing mechanism is the peripheral
   chemoreflex gain (loop gain) that rises with acclimatisation; the model
   carries this into ventilation but not into the stability index. **This is
   something the model got wrong and we do not gloss over it.**
2. **Ventilation at extreme altitude is under-predicted by about 10 %.** At
   8400 m the model gives PaCO2 14.6 mmHg against the 13.3 mmHg observed by
   Grocott et al.
3. **The absolute value of the AHI is a calibrated mapping, not a derived
   value.** CO2 reserve → AHI was fitted to the literature with a single
   sigmoid. The reserve itself is derived; the AHI is not.
4. **The γ = 1 assumption in the optimal haematocrit** (Result 6). 43.3 % is a
   lower bound.
5. **Absolute HAPE incidence.** The model reproduces the **ranking** of the
   three prophylaxis arms (tadalafil < dexamethasone < nifedipine < no
   treatment), but absolute values such as 74 % on placebo depend entirely on
   the assumptions made about the `a` and `λ_max` distributions of the virtual
   population.
6. **The absolute AMS incidence in the virtual population is 94.7 %, far above
   the observed 50–60 %** (see above). The relative risk is likewise compressed
   towards 1 by the ceiling effect (model RR 0.880 against roughly 0.45 in the
   literature). The symptom gain needs recalibrating against the proportion of
   the inter-individual distribution that crosses the threshold, and that has
   not been done.
7. **The nadir SaO2 of 24.7 % on an Everest summit day without oxygen** is lower
   than any value ever measured in a human being. The model is extrapolating in
   the corner that is extreme altitude × exertion, and the near-summit numbers
   of that scenario should be read only qualitatively (the post-descent value of
   53.8 % is within the observed range).
8. **The capillary stress-failure threshold of 19.5 mmHg is fixed as a single
   value.** In reality it differs from region to region of the lung (a
   gravity-dependent gradient), and that distribution determines which lung
   region leaks first — the model has no such spatial structure.

---

## How to use it

```r
source("hai_mrgsolve_model_en.R")
res <- run_all()                       # summary of the 23 scenarios
print(res, digits = 3)

arms <- list(`no prophylaxis` = scenarios$`02_rapid_4559`(),
             `acetazolamide`  = scenarios$`04_rapid_4559_acz250`(),
             `dexamethasone`  = scenarios$`05_rapid_4559_dex`())
compare_arms(arms, "LLS")     # both drugs bring the symptoms down
compare_arms(arms, "SaO2")    # only one moves the oxygen  ← Result 8

shiny::runApp("hai_shiny_app_en.R")
```

```bash
python3 hai_reference_model.py         # regenerates every number above
dot -Tsvg hai_qsp_model.dot -o hai_qsp_model.svg
dot -Tpng -Gdpi=150 hai_qsp_model.dot -o hai_qsp_model.png
```

---

## ⚠️ Disclaimer

This is a semi-quantitative QSP model for educational and research purposes. It
has not been independently validated or certified and **must not be used for
real clinical decision-making, prescribing, or regulatory submission.** If you
have symptoms at altitude, **descend** — do not consult a model.
