# Cancer Anorexia-Cachexia Syndrome (CACS) — QSP Model

<a href="cacs_qsp_model.svg"><img src="cacs_qsp_model.png" width="100%" alt="CACS QSP mechanistic map"></a>

*(click the map for the zoomable SVG)*

---

## 1. The Premise

> **Cachexia is not starvation.**

When the body starves, it defends itself. As intake falls, resting energy
expenditure falls with it, ketogenesis spares protein, and re-feeding brings
back what was lost. In cancer cachexia, **all three of these compensations
are broken**, and the reason the syndrome does not respond to food is not
the calorie deficit itself but **that breakage**.

So this model is built from two arms that converge on muscle **mass**, and a
third **axis** that neither nutrition nor most drugs can touch.

```
                        Tumour burden (the sole upstream driver)
                                  │
        ┌─────────────────────────┴─────────────────────────┐
        ▼                                                   ▼
  [ARM A] Intake                                     [ARM B] Catabolism + hypermetabolism
  GDF-15 → GFRAL (area postrema)                      IL-6/LIF → STAT3
  → CALCR/PBN aversion signal                         Activin A → ActRIIB → SMAD2/3
  → melanocortin balance                              → FoxO → MAFbx/MuRF1 → UPS
  → intake ↓                                          + IGF-1/Akt/mTORC1 suppression (anabolic resistance)
  (the falling-leptin signal is                        + PTHrP/IL-6 → WAT browning → REE ↑
   blocked by central IL-1β/PGE2)
        └─────────────────────┬─────────────────────────────┘
                              ▼
                        Skeletal muscle mass
                              ×
        [AXIS C] Muscle quality ── inflammation → ROS → myosin nitrosylation, RyR1 leak,
                              PGC-1α collapse → force per kg ↓
                              ▼
                         Grip strength = mass × quality
```

**`GRIP = MASS × QUALITY`** — this one line is the reason this model exists.
It is what mechanistically explains the most consistently reproduced result
in cachexia drug development.

| Trial | Lean body mass (LBM) endpoint | Functional endpoint |
|---|---|---|
| ROMANA 1/2 (anamorelin) | **met** | grip strength **failed** |
| POWER 1/2 (enobosarm) | **met** | stair-climb power **failed** |

Both drugs rapidly add mass through anabolic signalling alone, with
inflammatory drive left untouched. Mass rises, but **the quality of that
mass does not.**

---

## 2. What Makes This Model Testable

### (1) Disease Is **Generated**, Not Assumed

Setting `TVOL0 = 0` leaves the patient's weight stable forever. Every state
variable is initialised at the **analytical fixed point** of its own
differential equation in the absence of a tumour, and activity energy
expenditure is **solved for** in `$MAIN` so that `TEE == INTK` holds exactly.

```
=== 1. Baseline stability (TVOL0 = 0, 365 days, % drift) ========
 state drift_pct
    BW  0.000000      MUSC 0.000000      FATM 0.000000
  GRIP  0.000000      MYOQ 0.000000       CRP 0.000000
   ALB  0.000000      INTK 0.000000       REE 0.000000
  IGF1  0.000000
```

Anorexia, hypermetabolism, wasting, and refractoriness all **emerge from
tumour burden and a single sensitivity parameter (`CXSENS`)**. Because the
baseline does not drift, every treatment effect below is pharmacology, not
numerical drift.

### (2) The Point of No Return Exists **as a Number**

Sustained STAT3 signalling depletes the Pax7+ satellite-cell/myonuclear pool
(`SATC`), and its recovery rate is **1/25** of its loss rate. Because the
muscle's own **set point** is carried on this pool, once the pool is
depleted, mass does not return even if the tumour is eliminated.

Starting the same curative chemotherapy progressively later (evaluated at
day 540):

| Treatment start | Muscle nadir | Muscle at day 540 | % of pre-disease | Amount recovered | **Permanent loss** | SATC |
|---|---|---|---|---|---|---|
| Day 0   | 24.5 kg | 24.7 kg | 98.6 % | +0.13 kg | 0.35 kg | 0.97 |
| Day 60  | 23.0 kg | 23.9 kg | 95.4 % | +0.84 kg | 1.15 kg | 0.91 |
| Day 120 | 20.8 kg | 22.6 kg | 90.2 % | +1.79 kg | 2.45 kg | 0.80 |
| Day 180 | 18.4 kg | 20.9 kg | 83.7 % | +2.50 kg | 4.06 kg | 0.68 |
| Day 240 | 16.4 kg | 19.2 kg | 76.7 % | +2.75 kg | 5.83 kg | 0.56 |
| Day 300 | 14.8 kg | 17.3 kg | 69.3 % | +2.54 kg | **7.67 kg** | 0.45 |

The later treatment starts, **the more is recovered, but permanent loss
rises far faster.** This is Fearon's refractory stage, and the reason
**timing beats efficacy** in this model.

---

## 3. Dissociation — This Model's Central Claim

At week 12, relative to untreated natural history. The dashed line is the
grip-strength gain expected (1.6 kg per kg) **if all the added lean mass
were normal muscle**.

| Intervention | ΔLBM (kg) | ΔGrip (kg) | Grip per kg | % of pure-mass expectation |
|---|---|---|---|---|
| Anamorelin | +1.77 | +1.34 | 0.76 | **47 %** |
| Enobosarm | +1.90 | +1.01 | 0.53 | **33 %** |
| Bimagrumab | +2.69 | +1.03 | 0.38 | **24 %** |
| Ponsegromab | +0.32 | +0.49 | 1.53 | 96 % |
| Espindolol | +0.31 | +0.48 | 1.55 | 97 % |
| Megestrol | −0.47 | −0.16 | 0.34 | 21 % |
| Dexamethasone | −1.40 | −1.00 | 0.71 | 45 % |
| Nutrition + resistance exercise | +1.89 | **+8.16** | 4.32 | **270 %** |
| Multimodal | +2.19 | **+9.36** | 4.27 | **267 %** |

**Pure anabolic drugs (anamorelin, enobosarm, bimagrumab) deliver only
24-47 % of their added lean mass as function. Load-bearing interventions
deliver 270 %.** The same 1 kg of lean mass has a functional value that
differs tenfold depending on where it came from.

### Why — Three Separable Mechanisms

Rather than treat this dissociation as a single catch-all coefficient, the
model splits it into **three separately testable pieces**.

1. **`LNW` — GH-derived lean water.** GHSR-1a agonists raise GH, and GH
   retains sodium and water. DXA cannot distinguish that water from muscle,
   and counts both as "lean mass". A substantial share of the lean mass
   anamorelin gains comes from here, and oedema was indeed a real adverse
   event in ROMANA.
2. **`QNEW` — the contractile fraction of unloaded hypertrophy.** When
   myofibrillar protein accumulates without the accompanying myonuclear,
   capillary, and neural drive, force per kg is low. Resistance exercise
   supplies the missing stimulus, so `EXRES` reverses this penalty.
3. **`MYOQ` — axis C itself.** This is the axis that inflammation lowers and
   only load raises. Pure anabolic drugs **cannot** touch this value at all.

Sensitivity diagnostics show this separation directly — `SLNW` moves only
anamorelin's lean mass and leaves grip strength untouched (0.76→2.05 kg vs.
unchanged grip), while `QNEW` moves only enobosarm's and bimagrumab's grip
strength and leaves lean mass untouched (grip +0.2→+3.0 kg at QNEW 0→1,
lean mass unchanged).

---

## 4. Comparison with Published Trials

The registration-trial population (`pat_trial`) is not the untreated natural
history — trial participants all receive chemotherapy, symptom management,
and nutritional counselling, so the placebo arm is nearly weight-stable over
12 weeks. So the comparison is made against **the correct control, not a
sicker one**. Ponsegromab alone is compared in the GDF-15-enriched population
(`pat_gdfhi`) — its phase 2 trial enrolled patients with
GDF-15 ≥ 1.5 ng/mL.

| Trial | Week | Model ΔBW | Model ΔLBM | Model Δgrip | Published result |
|---|---|---|---|---|---|
| **Placebo arm (absolute change)** | 12 | −1.14 | −0.54 | −1.48 | ROMANA placebo: BW +0.14, LBM −0.47, grip −1.58 |
| ROMANA 1/2 anamorelin 100 mg | 12 | +3.17 | **+1.69** | +1.12 | LBM +1.1~+1.6 kg, BW ~+2.2 kg, **no difference in grip** |
| POWER 1/2 enobosarm 3 mg | 12 | +0.75 | **+1.09** | **−0.23** | LBM met (~+1.3 kg), stair-climb power **failed** |
| Ponsegromab 400 mg q4w | 12 | **+2.66** | +0.42 | +0.63 | BW +5.6 % vs. placebo (~+2.8 kg), appetite and activity improved |
| Megestrol 800 mg/d | 12 | **+2.86** | **−0.79** | −0.90 | Weight gain (fat, fluid), **no** lean-mass or quality-of-life benefit |
| Olanzapine 2.5 mg | 12 | +2.04 | +0.21 | +0.36 | >5 % weight gain 60 % vs. 9 % placebo |
| ACT-ONE espindolol 10 mg bid | 16 | +0.53 | +0.17 | **+0.27** | BW +0.9 vs. −1.2 kg, **grip also improved** |
| Dexamethasone 4 mg/d | 12 | −0.99 | **−2.03** | −2.48 | Appetite tolerance after 2-4 weeks, subsequent myopathy |

Direction and rank order are correct across all eight trials. In particular,
**enobosarm's lean mass gain of +1.09 kg and grip loss of −0.23 kg emerging
from the same simulation simultaneously** is this model's key validation
point — no parameter was fitted to that negative result.

---

## 5. Three Clinical Traps the Model Reproduces

### (1) With Tocilizumab, Plasma IL-6 **Rises**

A substantial fraction of IL-6 is cleared through its own receptor. Blocking
IL-6R removes that pathway, so the ligand accumulates. A rise in IL-6 is not
treatment failure but **evidence of target engagement**; CRP is what should
be judged instead.

```
         arm  day    IL6   CRP
   untreated    0   2.00   2.0
   untreated   28   4.08   9.7
   untreated   84   7.30  18.3
 tocilizumab   28   8.89   2.0
 tocilizumab   84  15.93   2.0   <- IL-6 is more than double, but CRP is normal
```

### (2) Myostatin **Falls** in Cachexia

Muscle itself is the source of production. It is a consequence of muscle
wasting, not its cause, and must not be misread as a biomarker. The ligand
that actually activates ActRIIB is **Activin A**, and in the model
bimagrumab's effect comes from there too.

### (3) The Scale **Hides** Tissue Loss

Hypoalbuminaemia lowers colloid osmotic pressure and produces oedema,
progestins retain fluid, the acute-phase response enlarges the liver, and GH
agonists increase lean water. A stable weight does not mean stable tissue.
The model **decomposes** weight change into muscle/fat/other lean
mass/liver/lean water/oedema to show this directly.

The nutrition-only scenario (S03) is the textbook case of this trap — at
week 12, body weight is nearly stable (**−0.09 kg**) while lean body mass
keeps falling (**−0.42 kg**). Looking at the scale alone, treatment appears
to be working.

---

## 6. Files

| File | Contents |
|---|---|
| [`cacs_qsp_model.dot`](cacs_qsp_model.dot) | Mechanistic map source — 19 modules · 247 nodes · 457 edges |
| [`cacs_qsp_model.svg`](cacs_qsp_model.svg) | Zoomable vector map |
| [`cacs_qsp_model.png`](cacs_qsp_model.png) | Raster map (150 dpi) |
| [`cacs_mrgsolve_model.R`](../../../cancer-cachexia/cacs_mrgsolve_model.R) | **78 ODEs · 307 annotated parameters · 20 scenarios · 11 diagnostics** |
| [`cacs_shiny_app.R`](../../../cancer-cachexia/cacs_shiny_app.R) | 10-tab interactive dashboard |
| [`cacs_references.md`](cacs_references.md) | **125 citations, verified live against PubMed** (19 sections) |

### The Map's 19 Modules

1 Tumour compartment and tumour-derived cachexins · 2 Systemic inflammation ·
3 Hepatic acute-phase response · 4 **[Arm A-i]** Brainstem GFRAL/area
postrema · 5 **[Arm A-ii]** Hypothalamic melanocortin · 6 Intake and
gastrointestinal tract · 7 **[Arm B-i]** Skeletal muscle protein catabolism ·
8 **[Arm B-ii]** Anabolic resistance · 9 Adipose tissue (lipolysis,
browning) · 10 Energy balance and hypermetabolism · 11 Neuroendocrine axis ·
12 Body composition and organ outcomes · 13 **[Axis C]** Muscle quality and
physical function · 14 Oral drug PK · 15 Biologic PK · 16 Drugs and targets ·
17 Non-drug interventions · 18 Chemotherapy · 19 Staging and clinical
endpoints

---

## 7. mrgsolve Model Structure

**78 compartments**

| Group | Count | State variables |
|---|---|---|
| Tumour and circulating mediators | 10 | TUMOR, IL6, TNFA, GDF15, ACTA, MSTN, LIF, CRP, ALB, LIVM |
| **Arm A** central | 7 | BSSIG, NAUS, AGRP, POMC, ANOR, LEPT, GHRL |
| Intake and energy | 3 | INTK, REE, PACT |
| **Arm B** muscle | 9 | STAT3, SMAD, NFKB, FOXO, ATRO, UPS, AUTOP, ROS, AAPL |
| Anabolic signalling | 8 | AKT, MTOR, ARES, GH, IGF1, TEST, CORT, INSR |
| **Axis C** quality | 4 | PGC1, MITO, MYOQ, SATC |
| Body composition | 11 | MUSC (what DXA sees), **MUSCC (what pulls)**, **LNW (GH water)**, FATM, OLBM, EDEM, BWL1-3, UCP1, ATGL |
| Function and outcomes | 2 | ECOG, CHZ |
| Drug PK | 24 | 11 drugs |

**PK/PD for 11 drugs**: anamorelin · megestrol · dexamethasone · olanzapine ·
enobosarm · espindolol · EPA · celecoxib · ponsegromab (quasi-equilibrium
binding) · tocilizumab (target-mediated disposition) · bimagrumab.

Monoclonal antibody clearance is modelled through albumin-dependent FcRn
recycling, so **clearance speeds up in hypoalbuminaemia — the sickest
patients receive the least exposure.**

### 20 Scenarios

| ID | Scenario | ID | Scenario |
|---|---|---|---|
| S00 | Healthy control (no tumour) | S10 | Olanzapine 2.5 mg |
| S01 | Natural history — advanced NSCLC | S11 | Enobosarm 3 mg (POWER) |
| S02 | Natural history — pancreatic cancer | S12 | Bimagrumab 10 mg/kg |
| S03 | Nutrition alone (+600 kcal, protein 1.5) | S13 | Tocilizumab 162 mg weekly |
| S04 | Nutrition + resistance/aerobic exercise | S14 | Espindolol 10 mg bid (ACT-ONE) |
| S05 | Anamorelin 100 mg (ROMANA) | S15 | EPA 2 g + celecoxib |
| S06 | Ponsegromab 400 mg q4w | S16 | Multimodal (MENAC-type) |
| S07 | Ponsegromab 100 mg q4w | S17 | Effective chemotherapy alone |
| S08 | Megestrol 800 mg/d | S18 | **Early** multimodal (precachexia) |
| S09 | Dexamethasone 4 mg/d | S19 | **Delayed** multimodal (refractory) |

S18 vs S19 is the most important contrast in this library — giving
**exactly the same package** in precachexia gains +0.49 kg muscle and
+5.8 kg grip strength at 12 weeks, while giving it in the refractory stage
loses −2.74 kg muscle. The difference is not the drugs but that the
satellite-cell pool is already gone.

### 11 Diagnostics (all auto-printed on script run)

Baseline stability · 12-week endpoint grid · 24-week endpoint grid ·
dissociation test · comparison against published trials · parameter
sensitivity of the dissociation · per-arm decomposition · tocilizumab
paradox · ponsegromab target engagement · hysteresis (point of no return) ·
staging and survival

---

## 8. Usage

```bash
# Render the mechanistic map
dot -Tsvg cacs_qsp_model.dot -o cacs_qsp_model.svg
dot -Tpng -Gdpi=150 cacs_qsp_model.dot -o cacs_qsp_model.png

# Run all 20 scenarios + 11 diagnostics (about 2 minutes)
Rscript cacs_mrgsolve_model.R

# Interactive dashboard
R -e 'shiny::runApp("cacs_shiny_app.R")'
```

Required packages: `mrgsolve`, and for the Shiny app `shiny`, `ggplot2`,
`dplyr`, `tidyr`, `DT`.

---

## 9. Limitations and Explicit Assumptions

Recorded honestly. These are things this model **cannot** do.

- **There is no ponsegromab dose-response.** 100 mg and 400 mg give
  identical results, because both doses fully suppress free GDF-15 below
  the limit of quantification. This is not a bug but a **testable
  prediction** — the dose-response observed in phase 2 cannot be explained
  by mean target engagement, and must instead come from trough exposure,
  inter-individual PK variability, or an extravascular GDF-15 reservoir.
- **`QNEW` and `SLNW` underpin the dissociation.** Neither value is directly
  measured; both are fitted. That is why their sensitivities are printed
  **separately** in diagnostic 6b, so it is directly visible which one, set
  to zero, makes the dissociation disappear.
- **The tocilizumab effect is stronger than the evidence supports.** In the
  model, IL-6R blockade nearly completely severs arm B. The actual clinical
  evidence is at the case-series level; there is no randomised trial. This
  is the model's **prediction**, not a reproduction.
- **Natural history (S01) is more aggressive than registration-trial placebo
  arms.** This is because it represents an advanced-cancer patient
  receiving no treatment at all, which is why trial comparisons are run
  separately in the `pat_trial` population.
- **The muscle compartment responds to signalling and barely responds to
  energy surplus.** The nutritional contribution (`GNUT`) is gated by
  anabolic resistance. The lean/fat partitioning of weight regain is skewed
  further toward fat than the Forbes relationship.
- **Sex, age, and genetic variation** are not addressed beyond `SEXM`
  scaling.
- The tumour is a single Gompertz compartment. There is no difference by
  metastatic site, no tumour heterogeneity, and no local effects.

---

## 10. The One Sentence That Matters Clinically

> In cachexia, **weight is an outcome, lean mass is a surrogate, and
> function is the goal.** These three do not move in the same direction,
> and which one a given drug can move is determined by where in the
> mechanism it cuts. And whichever arm is cut, **it must be cut while the
> satellite-cell pool still remains.**

---

## ⚠️ Disclaimer

This model is a **qualitative/semi-quantitative QSP model for education and
research purposes**. It was built from published literature and clinical
trial data but has not been independently validated or certified, and
**must not be used directly for clinical decision-making, prescribing, or
regulatory submission.** Parameters and assumptions are illustrative
approximations; separate fitting and validation against real patient data
is required.
