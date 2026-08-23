# Iron Deficiency Anaemia — QSP Model

> **One-line summary** — oral iron absorption is not a *rate* but **the product of three
> factors**, and the third factor (enterocyte export capacity) is set by **the hepcidin the
> immediately preceding dose produced.** In other words, the drug closes the very door it must
> pass through. And beneath that door, what limits haemoglobin is **an entirely different
> bottleneck** (bone-marrow iron supply, dominated by macrophage recycling), so absorbed amount
> and Hb response diverge from each other at high doses.

| Deliverable | File |
|---|---|
| 🗺️ Mechanistic map (150 nodes · 13 clusters · 225 edges) | [`ida_qsp_model.dot`](../../../iron-deficiency-anemia/ida_qsp_model.dot) · [SVG](../../../iron-deficiency-anemia/ida_qsp_model.svg) · [PNG](../../../iron-deficiency-anemia/ida_qsp_model.png) |
| ⚙️ mrgsolve ODE model (37 compartments · 122 parameters · 15 scenarios) | [`ida_mrgsolve_model.R`](ida_mrgsolve_model.R) |
| 📊 Shiny dashboard (8 tabs) | [`ida_shiny_app.R`](../../../iron-deficiency-anemia/ida_shiny_app.R) |
| 📚 References (69 PubMed links, all verified) | [`ida_references.md`](ida_references.md) |
| 🔁 Independent-verification Python reference implementation (no dependencies) | [`ida_reference_model.py`](../../../iron-deficiency-anemia/ida_reference_model.py) |

<p align="center">
  <a href="../../../iron-deficiency-anemia/ida_qsp_model.svg"><img src="../../../iron-deficiency-anemia/ida_qsp_model.png" width="760" alt="IDA QSP map"></a>
</p>

---

## 1. Why this disease is modelled this way

Iron deficiency anaemia looks like a disease where "iron is lacking, so give iron." But
clinical reality is not that simple. Giving the same total amount split daily versus every
other day gives different results; doubling the dose of intravenous iron does not make Hb rise
twice as fast; and when inflammation is present, Hb does not rise even though iron is
absorbed. All three of these "odd" facts come from **the fact that iron absorption is a
regulated gate**, and **the fact that the bottleneck limiting the rise in Hb sits somewhere
other than that gate.** This model writes each of those two constraints as an equation.

### Constraint 1 — the gate

```
absorption flux = A_LUM  ×  DMT1 dose  ×  FPN_ENT
                    │          │            └── enterocyte export capacity:
                    │          │                set by the hepcidin from the preceding dose
                    │          └── saturable (Km 11 mg)
                    └── luminal available iron
```

### Constraint 2 — the ceiling

The bone marrow uses **13.3 mg** of iron per day, of which about **85%** is recycled from
macrophages that have consumed senescent red cells. The share coming from the gut is, at best,
only around 10 mg/day. So even if iron supply is increased, it is the bone marrow that sets the
ceiling on the rate at which Hb can rise.

---

## 2. Baseline states — equilibria the model found for itself

The two reference states are **the model's equilibrium solutions, not input values.** Starting
from the normal state, raising only the blood-loss parameter (`VBLEED`) to 6.375 mL/day
(≈ 178 mL/cycle) and integrating to steady state gives the IDA state below.

| Measure | Normal (VBLEED 0.6) | IDA (VBLEED 6.375) | Reference normal range (female) |
|---|---|---|---|
| Haemoglobin | 13.86 g/dL | **9.02 g/dL** | 12.0–15.5 |
| MCH | 30.3 pg | **21.2 pg** | 27–33 |
| RBC | 4.57 ×10¹²/L | 4.26 ×10¹²/L | 4.0–5.2 |
| Serum iron | 110 µg/dL | **28 µg/dL** | 60–170 |
| TIBC | 316 µg/dL | **427 µg/dL** | 250–400 |
| TSAT | 34.8 % | **6.5 %** | 20–45 |
| Hepcidin | 6.43 ng/mL | **0.32 ng/mL** | 1–20 |
| Ferritin | 66.9 ng/mL | **5.8 ng/mL** | 15–150 |
| Storage iron | 463 mg | **15 mg** | 300–500 |
| Tissue (non-erythroid) iron | 386 mg | 316 mg | ~380 |
| Reticulocytes | 38 ×10⁹/L | 42 ×10⁹/L | 25–75 |
| CHr (reticulocyte Hb) | 30.3 pg | **20.5 pg** | > 28 |
| sTfR | 1.36 mg/L | **3.63 mg/L** | 0.8–1.8 |
| EPO | 8 mIU/mL | 94 mIU/mL | 4–25 |
| Iron absorption | 1.24 mg/day | 2.53 mg/day | 1–2 |
| Bone-marrow iron consumption | 16.3 mg/day | 13.3 mg/day | — |
| Total body iron | 2769 mg | 1585 mg | — |

Two things stand out. **First, reticulocytes are inappropriately low for the degree of
anaemia** (42 ×10⁹/L), because the bone marrow cannot expand even though EPO is raised to
94 mIU/mL — the model expresses this with a gate, `F_EXP = TSAT/(TSAT + 12)` (at TSAT 6.5%,
only 35% of the EPO expansion is realised). **Second, absorption is raised to 2.53 mg/day,
twice normal**, the result of the deficiency suppressing hepcidin and opening the gate. Yet
the disease persists because losses are still larger. The IDA equilibrium is "the point where
Hb falls, bleeding iron loss falls with it, and it finally equals absorption."

---

## 3. A single dose — the door opening and closing

A single oral 60 mg dose of elemental iron (IDA baseline state):

| Time | Event | Value |
|---|---|---|
| 0 h | Baseline | serum iron 27.9 µg/dL · hepcidin 0.318 ng/mL · FPN_ENT 0.210 |
| 5.75 h | Peak serum iron | **224 µg/dL** (TSAT 52.5 %) |
| 10.25 h | Peak hepcidin | **1.24 ng/mL (3.9× baseline)** |
| 12.2 h | Peak enterocyte ferritin | 0.70 → **5.38 mg** |
| 22 h | Nadir of export capacity | FPN_ENT 0.170 (**80.9 % of baseline**) |
| 24 h | — | hepcidin still 1.36× · export capacity 81.0 % |
| 48 h | — | export capacity 85.4 % |

**Fractional absorption 22.0 % (13.2 mg / 60 mg)** — matches the observed 20–22% in
iron-deficient women.

The key is that there are **two clocks.** Hepcidin itself clears fast, with a half-life of
2.5 h, but the enterocyte's export capacity returns on a much slower clock. The effective
recovery time constant is `1/(KSYN_FPE + KDEG_FPE·HEP^0.6)`, about 51 hours (≈2 days) at the
IDA baseline hepcidin. So **absorptive capacity is still low even after circulating hepcidin
has returned to normal.** What Hahn called the "mucosal block" in 1943 is the difference
between these two clocks.

### The refractory window — measured with a probe dose

The same 60 mg was given at various intervals after a conditioning dose, and the probe's own
absorption was measured as the difference from a control simulation without the probe.

| Interval | Probe absorption | Vs rested gut |
|---|---|---|
| +4 h | 9.94 mg | **75.4 %** |
| +8 h | 11.60 mg | 88.0 % |
| +12 h | 11.69 mg | 88.7 % |
| +24 h | 12.06 mg | **91.5 %** |
| +36 h | 12.16 mg | 92.3 % |
| +48 h | 12.42 mg | 94.2 % |
| +72 h | 12.63 mg | 95.8 % |

A second dose 4 hours later on the same day loses a quarter of its absorption (exactly the
phenomenon observed by Moretti 2015). It has still not fully recovered after 24 hours.

---

## 4. The key result — fractional absorption and total absorption rank in opposite orders

14 days of dosing, assuming full adherence, with elemental iron per dose fixed:

| Regimen | Doses | Total taken | Absorbed | Fractional absorption | mg/day | mg/dose |
|---|---|---|---|---|---|---|
| 60 mg q12h | 28 | 1680 mg | 234.0 mg | 13.9 % | **16.72** | 8.36 |
| 60 mg q24h | 14 | 840 mg | 147.7 mg | 17.6 % | 10.55 | 10.55 |
| 60 mg q48h | 7 | 420 mg | 84.4 mg | 20.1 % | 6.03 | 12.06 |
| 60 mg q72h | 5 | 300 mg | 62.8 mg | **20.9 %** | 4.48 | 12.55 |
| 120 mg q24h | 14 | 1680 mg | 188.7 mg | 11.2 % | 13.48 | 13.48 |
| 120 mg q48h | 7 | 840 mg | 113.3 mg | 13.5 % | 8.10 | 16.19 |
| 180 mg q48h | 7 | 1260 mg | 130.0 mg | 10.3 % | 9.28 | 18.57 |
| 30 mg q24h | 14 | 420 mg | 105.7 mg | **25.2 %** | 7.55 | 7.55 |

**The two metrics line up in opposite orders.** Lengthening the interval raises fractional
absorption (13.9 → 20.9 %) and lowers daily absorption (16.72 → 4.48 mg/day). In other words,
**alternate-day dosing is the optimum for efficiency and tolerance, not for total delivery.**
Per mg swallowed, alternate-day dosing yields 1.14× the iron of daily dosing, but the iron
actually delivered per day with daily dosing is 1.75× that of alternate-day dosing.

Comparing the same total amount split differently (1680 mg either way): 60 mg q12h absorbs
234.0 mg, while 120 mg q24h absorbs 188.7 mg. Here, **the gain from avoiding dose saturation
(Km 11 mg) outweighs the loss from the hepcidin refractory window.** That this balance flips
depending on dose and interval is the model's non-trivial prediction.

### Tolerability, a second currency

Turning on the terms by which GI symptoms erode adherence (`K_GI`, `EMAX_ADH`) penalises more
frequent dosing further. Actual absorption relative to full adherence:

| Regimen | Vs full adherence |
|---|---|
| 60 mg q12h | 86.2 % |
| 60 mg q24h | 87.9 % |
| 60 mg q48h | **93.2 %** |
| 60 mg q72h | 95.0 % |

The longer the interval, the smaller the loss. The reason clinical guidelines recommend
alternate-day dosing is that absorption efficiency and tolerability move in **the same
direction**, and this table also shows that delivered amount moves in the opposite direction.

---

## 5. 12-week treatment scenarios (baseline Hb 9.02 g/dL)

| # | Scenario | Hb week 2 | week 4 | week 8 | week 12 | ΔHb | TSAT | Ferritin | Storage iron | Time to +1 g/dL | +2 g/dL |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | No treatment | 8.99 | 8.95 | 8.88 | 8.81 | −0.21 | 6.2 | 5.6 | 14 | — | — |
| 2 | Oral 65 mg daily | 9.45 | 10.04 | 10.98 | 11.62 | +2.59 | 16.8 | 28.6 | 180 | 28 days | 58 days |
| 3 | Oral 65 mg twice daily | 9.75 | 10.62 | 11.70 | 12.26 | +3.23 | 32.3 | 45.7 | 309 | 18 days | 37 days |
| 4 | Oral 130 mg daily | 9.55 | 10.25 | 11.30 | 11.93 | +2.91 | 20.7 | 37.8 | 249 | 24 days | 48 days |
| 5 | Oral 130 mg on alternate days | 9.33 | 9.78 | 10.61 | 11.23 | +2.20 | 15.8 | 23.9 | 145 | 36 days | 74 days |
| 6 | Oral 65 mg on alternate days | 9.26 | 9.62 | 10.31 | 10.86 | +1.84 | 13.3 | 17.9 | 101 | 44 days | — |
| 7 | Oral 195 mg on alternate days | 9.36 | 9.87 | 10.76 | 11.40 | +2.38 | 17.3 | 27.5 | 172 | 33 days | 67 days |
| 8 | **IV FCM 1000 mg** | 10.28 | 11.44 | 12.05 | 12.22 | +3.20 | 24.4 | 46.8 | 314 | **12.5 days** | **21.5 days** |
| 9 | IV FCM 750 mg ×2 | 10.30 | 11.79 | 12.56 | 12.75 | +3.72 | 38.6 | 94.1 | 686 | 12.5 days | 19.5 days |
| 10 | IV derisomaltose 1000 mg | 10.28 | 11.44 | 12.05 | 12.22 | +3.20 | 24.4 | 46.8 | 314 | 12.5 days | 21.5 days |
| 11 | IV iron sucrose 200 mg ×5 | 9.97 | 11.24 | 12.29 | 12.38 | +3.36 | 25.7 | 50.7 | 343 | 15 days | 25.5 days |
| 12 | Oral daily + IL-6 20 | 9.02 | 9.11 | 9.32 | 9.48 | **+0.45** | 12.3 | 50.3 | 176 | — | — |
| 13 | IV FCM + IL-6 20 | 9.75 | 10.38 | 10.54 | 10.40 | +1.38 | 20.5 | 109.1 | 565 | 18.5 days | — |
| 14 | Oral daily + bleeding corrected | 9.63 | 10.38 | 11.59 | 12.38 | +3.36 | 21.1 | 34.0 | 219 | 21.5 days | 42 days |
| 15 | Bleeding correction alone | 9.17 | 9.31 | 9.57 | 9.83 | +0.81 | 7.8 | 6.4 | 20 | — | — |

A few things to note.

- **The benefit of IV iron is speed, not magnitude.** FCM 1000 mg produces a 12-week ΔHb of
  +3.20, essentially the same as oral 65 mg twice daily (+3.23). What differs is the time to
  get there: 21.5 vs 37 days to +2 g/dL. The reticulocyte peak also differs, 74 vs
  60 ×10⁹/L.
- **Correcting the bleeding alone is not enough (+0.81), and iron replacement alone is not
  enough either (+2.59).** Combining the two gives +3.36, approaching IV iron — the model
  quantifies the point that correcting the cause is a cheaper intervention than escalating the
  dose.
- **In inflammation, oral therapy essentially fails** (+0.45). Even IV iron manages only a
  partial response (+1.38). In scenario 12, ferritin rises to 50.3 while TSAT stays at 12.3% —
  the dissociation of markers in functional iron deficiency comes out of the mechanism, not
  from an assumption.

---

## 6. Escalating the IV iron dose buys storage iron, not speed

| IV dose | ΔHb week 1 | week 3 | week 12 | Max weekly slope | Storage iron wk 12 | Ferritin wk 12 |
|---|---|---|---|---|---|---|
| 200 mg | 0.29 | 0.87 | 0.88 | 0.42 g/dL/wk | 23 mg | 7 |
| 500 mg | 0.35 | 1.55 | 2.14 | 0.72 | 72 mg | 14 |
| 750 mg | 0.36 | 1.83 | 2.78 | 0.83 | 172 mg | 28 |
| **1000 mg** | 0.37 | 1.99 | 3.20 | 0.89 | 314 mg | 47 |
| 1500 mg | 0.38 | 2.16 | 3.68 | 0.93 | 662 mg | 91 |
| 2000 mg | 0.38 | 2.23 | 3.91 | 0.96 | 1021 mg | 135 |
| 3000 mg | 0.39 | 2.30 | 4.13 | 0.98 | 1720 mg | 219 |

Going from 500 mg to 3000 mg (6-fold), **the maximum weekly Hb slope rises only 36%**
(0.72 → 0.98), **the week-1 response rises only 11%** (0.35 → 0.39), yet **storage iron rises
24-fold** (72 → 1720 mg). What holds the clock is not the iron supply but the bone marrow.
This is why "filling the total iron deficit (Ganzoni)" and "correcting the anaemia quickly"
are different goals.

---

## 7. Inflammation — absorbed iron does not become haemoglobin

| IL-6 input | IL-6 steady state | Hepcidin | FPN_ENT | Absorbed mg/day | Lost mg/day | Net | ΔHb oral | ΔHb IV | Ferritin | TSAT |
|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 0 | 1.69 | 0.075 | 9.10 | 3.41 | +5.70 | +2.59 | +3.20 | 21.9 | 12.8 |
| 0.25 | 5 | 3.29 | 0.052 | 6.90 | 3.23 | +3.67 | +1.51 | +2.36 | 30.4 | 10.7 |
| 0.5 | 10 | 4.45 | 0.045 | 6.05 | 3.15 | +2.90 | +1.01 | +1.91 | 36.8 | 10.8 |
| 1.0 | 20 | 6.14 | 0.039 | 5.21 | 3.06 | +2.15 | +0.45 | +1.38 | 45.5 | 11.4 |
| 2.0 | 40 | 8.26 | 0.033 | 4.49 | 2.97 | +1.52 | **−0.07** | +0.83 | 55.0 | 12.4 |
| 4.0 | 80 | 10.41 | 0.030 | 3.96 | 2.90 | +1.07 | −0.50 | +0.37 | 63.4 | 13.3 |

At IL-6 of 40 pg/mL, the gut still **absorbs 4.5 mg/day** and the net iron balance is still
**positive at +1.52 mg/day**, yet Hb does not rise at all over 12 weeks (−0.07). This is not
because iron is lacking, but **because the iron is in the wrong place**: macrophage
ferroportin is shut, trapping it in storage (ferritin 55), while IL-6 simultaneously suppresses
EPO production, progenitors, and red-cell lifespan directly. So "giving more iron" is not the
answer in this situation — turning off the inflammation is.

## 8. IRIDA — the same microcytosis, the opposite prescription

A separate equilibrium state with `TMPRSS6` lowered to 0.30 (no bleeding):

| | Baseline | Oral 65 mg daily, 12 wk | IV FCM 1000 mg, 12 wk |
|---|---|---|---|
| Hb | 12.78 g/dL | +0.32 | +0.94 |
| Hepcidin | **11.51 ng/mL** (inappropriately high) | — | — |
| TSAT | 19.2 % | — | — |
| Ferritin | 56.9 ng/mL (normal) | — | — |
| Absorption | 0.93 mg/day | 291.6 mg (12 wk) | 44.1 mg |
| (for comparison) same regimen in acquired IDA | — | **+2.59** | **+3.20** |

Because the lesion is not "storage iron depletion" but "inappropriately high hepcidin",
ferritin is normal but oral iron does not work. Even having absorbed 292 mg over 12 weeks, Hb
rose only 0.32 g/dL. **The clinical recommendation to check hepcidin when oral iron fails in a
microcytic anaemia with normal ferritin is derived here.**

> Limitation: this model's IRIDA baseline Hb (12.78) is milder than the moderate anaemia of
> real IRIDA. The qualitative conclusion of oral refractoriness holds, but severity is
> underestimated.

---

## 9. Safety — equal efficacy, different costs

A comparison that changes only the carbohydrate shell of the IV iron. **Haematological
outcomes are identical** (scenario 8 vs 10: ΔHb +3.20 g/dL in both), but phosphate metabolism
differs entirely.

| Regimen | iFGF23 peak | Peak day | Phosphate nadir | Nadir day | Days <2.0 mg/dL | Days <2.5 | 1,25D nadir | PTH peak |
|---|---|---|---|---|---|---|---|---|
| FCM 1000 mg | **231 pg/mL** | 2.0 | **1.98 mg/dL** | 6.8 | **4.2** | 23.8 | 27.3 | 64 |
| FCM 750 mg ×2 | 234 | 8.8 | **1.92** | 11.2 | **9.2** | 31.8 | 27.1 | 64 |
| FCM 500 mg | 198 | 2.0 | 2.10 | 6.0 | 0.0 | 17.8 | 28.3 | 63 |
| Derisomaltose 1000 mg | 40 (no change) | — | 3.02 | — | **0.0** | 0.0 | 37.0 | 50 |
| Iron sucrose 200 mg ×5 | 69 | 29.5 | 2.75 | 31.2 | 0.0 | 0.0 | 34.5 | 53 |

Carboxymaltose suppresses FGF23 cleavage and raises intact FGF23 (→ phosphate loss via
NaPi-2a/2c inhibition, 1,25D falls via CYP27B1 inhibition, secondary PTH rise), whereas
derisomaltose does not. **Efficacy metrics alone cannot distinguish the two drugs** — only a
mechanistic model can show this asymmetry, and this is where QSP tells you more than a
clinical-trial summary table. Splitting FCM 750 mg into two doses raises the total to 1500 mg
and more than doubles the days of hypophosphataemia, from 4.2 to 9.2 days.

## 10. Tissue iron recovers later than haemoglobin

| Scenario | Hb wk 4 | Tissue iron wk 4 | IRLS wk 4 | Hb wk 12 | Tissue iron wk 12 | Tissue-iron recovery wk 12 | Hb recovery wk 12 |
|---|---|---|---|---|---|---|---|
| No treatment | 8.95 | 314 mg | 5.6 | 8.81 | 309 mg | −9 % | −4 % |
| Oral 65 mg daily | 10.01 | 320 mg | 5.1 | 11.56 | 329 mg | **18 %** | **53 %** |
| Oral 130 mg on alternate days | 9.77 | 320 mg | 5.1 | 11.19 | 327 mg | 16 % | 45 % |
| IV FCM 1000 mg | 11.44 | 395 mg | **0.0** | 12.22 | 390 mg | **105 %** | 66 % |

Oral treatment fills 53% of the Hb deficit by 12 weeks while filling only 18% of the
tissue-iron deficit. IV iron does the opposite, filling tissue iron first (105%), because the
bone marrow takes circulating iron preferentially. **Patients whose fatigue and restless legs
persist despite corrected anaemia**, and **patients whose fatigue improves with IV iron despite
having no anaemia** (Krayenbuehl 2011, FAIR-HF), both arise from the same structure.

---

## 11. Model structure

37 state variables:

| Compartment group | State variables |
|---|---|
| Gut (3) | `A_LUM` luminal available iron · `A_ENT` enterocyte available iron · `A_EFT` enterocyte ferritin iron (the mucosal block's memory) |
| Plasma and storage (6) | `A_TF` transferrin-bound iron · `A_NTBI` · `A_COL` IV colloid · `A_RES` macrophage · `A_LIV` hepatocyte · `A_TISS` non-erythroid tissue |
| Regulatory (5) | `TIBC` · `HEP` hepcidin · `FPN_ENT` enterocyte export capacity · `FPN_RES` macrophage/hepatocyte ferroportin · `ERFE` · `IL6` |
| Erythroid lineage (11) | `PROG` · `N1–N3` erythroblast cell counts · `H1–H3` those cells' Hb mass · `NRET`/`HRET` · `NRBC`/`HRBC` |
| Biomarkers (2) | `FERR` · `STFR` |
| Phosphate axis (5) | `CLV` cleavage-inhibition signal · `FGF23` · `PHOS` · `CTRIOL` · `PTH` |
| Tolerability and audit (4) | `GI` · `CUM_ABS` · `CUM_IV` · `CUM_LOSS` |

Three structural choices determine the character of this model.

1. **Cell count and Hb mass are tracked as two separate chains.** Only this way do MCH and CHr
   fall out of the calculation, and only this way can "a small number of normal cells" be
   distinguished from "a large number of hypochromic cells."
2. **Give ferroportin two clocks:** the enterocyte clock (`KSYN_FPE` 0.004/h, slow) creates the
   24–48 hour memory of the mucosal block, while the macrophage clock (`KSYN_FPR` 0.035/h,
   fast) creates the hours-scale response of inflammatory hypoferraemia. A single clock cannot
   reproduce both.
3. **Put an iron gate on EPO's marrow-expansion effect** (`F_EXP`). Without it, the model
   responds to anaemia with EPO by manufacturing cells without limit, and reticulocytes come
   out higher than reality.

### Mass conservation audit

`BODY_FE(t) = BODY_FE(0) + CUM_IV + CUM_ABS − CUM_LOSS` was verified across a 12-week
simulation. Error **~1×10⁻¹² %** (the limit of double precision). The enterocyte pools
(`A_ENT`, `A_EFT`) were deliberately excluded from body iron — iron that entered the enterocyte
and was shed with the cell never actually entered the body, so counting it would inflate both
sides together.

### Cross-validation with a dual implementation

Every result was computed with **two independent implementations** and compared: mrgsolve
(C++/LSODA) and a pure-Python fixed-step RK4 implementation in
[`ida_reference_model.py`](../../../iron-deficiency-anemia/ida_reference_model.py). The two
implementations share no code, so agreement between them is evidence that the results come
from **the equations themselves**, not from the behaviour of a particular solver. Example
(single 60 mg dose): peak serum iron 224.2 vs 224.2 µg/dL, peak hepcidin 1.243 vs 1.243 ng/mL,
fractional absorption 21.97 vs 21.97 %, +24 h probe 91.5 vs 91.5 %, 14-day q24h absorption
147.70 vs 147.70 mg.

---

## 12. Usage

```r
# Run the full ODE model + 15 scenarios + validation notes
Rscript -e 'options(ida.report=TRUE); source("ida_mrgsolve_model.R")'

# Interactive dashboard (8 tabs)
Rscript -e 'shiny::runApp("ida_shiny_app.R")'

# Reproduce without R (standard library only)
python3 ida_reference_model.py

# Re-render the map (needs newrank)
dot -Tsvg ida_qsp_model.dot -o ida_qsp_model.svg
dot -Tpng -Gdpi=150 ida_qsp_model.dot -o ida_qsp_model.png
```

When moving to a different patient, change `BV_L`, `PV_L`, `VBLEED`, `DIET`, `IL6_IN` and then
**you must re-equilibrate** (integrate for 2–3 years with no treatment). The default initial
values are the IDA equilibrium state in the table above, and that value is not an equilibrium
under other parameters. The Shiny app's "re-equilibrate baseline" checkbox performs this task.

## 13. Validation and limitations

The comparison against the literature and the explicit limitations are set out in the last two
sections of [`ida_references.md`](ida_references.md). In summary,

- What was reproduced: single-dose fractional absorption (22.0% vs the observed 20–22%), the
  24-hour-sustained hepcidin rise (1.36-fold), the reduced absorption of a same-day repeat dose
  (75.4%), the timing of the reticulocyte peak (12–14 days), the magnitude and timing of FCM
  hypophosphataemia and its absence with derisomaltose, and the normal/IDA reference values
  overall.
- **What was underestimated**: the alternate-day-to-daily fractional absorption ratio is 0.87
  in the model vs 0.75 observed — this model estimates the refractory-window penalty
  **conservatively**, because with the circulating hepcidin half-life fixed at 2.5 h, building
  the 24-hour memory from the enterocyte recovery clock alone forces a trade-off between depth
  and persistence.
- Otherwise: the serum-iron peak is somewhat late (5.8 h vs 2–5 h), the IRIDA phenotype is
  mild, the GI-adherence term is semi-quantitative, and the model is calibrated to a single
  60 kg female body size.

> ⚠️ A model for education and research purposes. It cannot be used for actual clinical
> practice, prescribing, or regulatory submission.
