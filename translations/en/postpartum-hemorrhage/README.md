# Primary Postpartum Haemorrhage (PPH) — QSP Model

**PPH is not a coagulation disease. It is a flux problem against a reservoir that
pregnancy has enlarged to look big, and the valve is mechanical.**

At term, uterine blood flow is 750 mL/min — 10-15% of a 7 L/min cardiac output —
pouring into the cut face of roughly 120 spiral arteries that trophoblastic
remodelling has stripped of their own smooth muscle, leaving them **unable to
contract.** The circulating-volume increase that pregnancy has banked is
2100 mL (70 → 100 mL/kg), so

> **2100 mL ÷ 750 mL/min = 2.8 minutes**

That is the entire time allowance pregnancy has purchased. The model's fully
atonic uterus (ATN = 1) reaches a peak of 516 mL/min, crosses the ATLS class III
threshold (2100 mL) at **6.0 minutes**, and reaches exsanguination at
**10 minutes**. Every treatment-related claim in this document is a claim about
how much of that clock it buys back.

| File | Contents |
|------|----------|
| [`pph_qsp_model.dot`](../../../postpartum-hemorrhage/pph_qsp_model.dot) · [SVG](../../../postpartum-hemorrhage/pph_qsp_model.svg) · [PNG](../../../postpartum-hemorrhage/pph_qsp_model.png) | Mechanistic map — 137 nodes / 158 edges / 12 clusters |
| [`pph_mrgsolve_model.R`](../../../postpartum-hemorrhage/pph_mrgsolve_model.R) | mrgsolve ODE model — 41 compartments (10 drug PK + 31 disease/physiology), 12 scenario blocks |
| [`pph_shiny_app.R`](../../../postpartum-hemorrhage/pph_shiny_app.R) | Shiny dashboard — 9 tabs |
| [`pph_references.md`](../../../postpartum-hemorrhage/pph_references.md) | 90 references (every PMID verified via the PubMed E-utilities) |

---

## What this model argues (nine claims, all derived from numbers)

### 1. The valve is Poiseuille, not linear

Haemostasis at the placental attachment site is neither suturing nor
coagulation — it is a mechanical process in which myometrial bundles (the
"**living ligature**") clamp the vessels shut. Flow through a compressed vessel
scales with the fourth power of its radius, so the model writes
`patency = (1 − TONE)^4`.

| Uterine tone | Patency | Leak |
|---|---|---|
| 0.91 (firm uterus) | 6×10⁻⁵ | 0.05 mL/min |
| 0.50 | 0.063 | 47 mL/min |
| 0.20 (boggy uterus) | 0.410 | 307 mL/min |

In other words, the bedside distinction of "firm vs boggy" spans only a few
percentage points of contractility, and the model is correspondingly
**bistable** across that same narrow band.

| ATN (atony severity) | Untreated outcome |
|---|---|
| 0.60 | 2218 mL, survives, MAP nadir 78 mmHg |
| **0.63** | **death at 70 minutes (3544 mL)** |

Three percentage points of contractility separate survival from death.

### 2. Coagulation stays locked out until the machine engages

Clot formation is gated by shear: `fSH = 1/(1 + (Q/80)²)`. And a clot forming
over a vessel that the myometrium is not clamping down can never occlude more
than 35% of patency (`CSEAL0`). Model values at 10 minutes:

| | Leak | Shear gate | Placental-bed clot |
|---|---|---|---|
| Normal | 0 mL/min | 1.00 | → 0.98, maturing → 1.00 |
| Moderate atony | 39 | 0.81 | 0.83 |
| **Severe, untreated** | **168** | **0.19** | **stalls at 0.22** |
| Severe, treated | 22 | 0.93 | 0.87 |

So "uterotonics and tamponade before blood products" is not a convention but
**arithmetic**. On top of that, durable haemostasis needs time for the clot to
**organise** under low shear (state variable `MAT`) — what every mechanical
intervention actually buys is that time.

### 3. Oxytocin erases its own target

The OTR internalises with agonist exposure (`KDES` 0.08/min, occupancy-driven)
and recovers with a half-life of 2.6 hours. Two consequences follow:

- After 8 hours of oxytocin-augmented labour alone, steady-state OTR reaches
  **0.39** (consistent with the reported 30-60% reduction in myometrial binding
  sites after augmented labour)
- **Treatment itself desensitises**: a standard infusion starting from OTR 1.00
  falls to **0.09** by 120 minutes

For a patient with severe atony after augmented labour (2605 mL on
oxytocin+massage), the two remaining options are not equivalent.

| Option | Blood loss | Reduction |
|---|---|---|
| Double the oxytocin | 2462 mL | −143 |
| Quadruple the oxytocin | 2328 mL | −277 |
| **Add carboprost (FP)** | **1680 mL** | **−925** |
| Add ergometrine (α₁/5-HT₂A) | 1781 mL | −824 |
| Add misoprostol (EP2/EP3) | 1953 mL | −652 |

**Switching receptors is worth 3-6 times more than raising the dose of a
receptor that has already vanished.**

### 4. One parameter builds the entire treatment ladder

A single `ATN` sets both the endogenous contractile drive (`DRVI`) and the
**ceiling** (`CAP`) that drugs can push against. Blood loss (mL); D@n = death at
minute n:

| ATN | Untreated | Massage | Oxytocin | Oxytocin+massage | 4-agent ladder |
|---|---|---|---|---|---|
| 0.60 | 2218 | 940 | 1078 | 743 | 641 |
| 0.63 | D@70 | 1170 | 1434 | 882 | 741 |
| 0.66 | D@44 | 1504 | 2081 | 1072 | 870 |
| 0.70 | D@30 | 2236 | **D@77** | 1446 | 1102 |
| 0.74 | D@22 | **D@58** | **D@51** | **2054** | 1427 |
| 0.78 | D@18 | D@32 | D@38 | **D@93** | **1877** |
| 0.84 | D@13 | D@19 | D@24 | D@37 | **D@57** |
| 0.90 | D@11 | D@13 | D@16 | D@20 | D@21 |

The ATN 0.74 row shows **superadditivity**: massage alone fails at 58 minutes
and oxytocin alone at 51 minutes, but **combining them survives at 2054 mL** —
because two sub-threshold drives sum past the threshold together. Above ATN
0.84 the entire drug ladder falls short of the threshold — this is the model's
internal definition of "uterotonic-refractory atony".

### 5. Fibrinogen crosses the critical line first — by ratio arithmetic

What determines which factor is depleted first in a diluting, bleeding patient
is simply the **ratio of threshold to baseline**.

| | Threshold | Term baseline | Ratio |
|---|---|---|---|
| **Fibrinogen** | 2.0 g/L | 4.5 g/L | **0.44 ← highest** |
| Coagulation factors (pooled) | 0.30 | 1.00 | 0.30 |
| Platelets | 50 ×10⁹/L | 250 ×10⁹/L | 0.20 |

In the crystalloid-first arm, the model confirms the order: fibrinogen crosses
2.0 g/L at a cumulative blood loss of **3476 mL**, Hb crosses 7 g/dL at
3555 mL, and platelets (nadir 71) and coagulation factors (nadir 0.33) **never
reach threshold at all.**

Corollary — also arithmetic: **FFP cannot raise a low fibrinogen.** FFP's own
fibrinogen concentration is 2 g/L, and raising 7 L of plasma by 1 g/L requires
7 g = **14 units = 3.5 L** of FFP. Because FFP's concentration is itself the
transfusion target, FFP always pulls fibrinogen **towards** 2 g/L, from either
side.

### 6. Resuscitation that saves the circulation can wreck the ligature

Calcium is both a coagulation cofactor and **the smooth-muscle contraction
ion**. So citrate-induced ionised hypocalcaemia is a double hit specific to
obstetric haemorrhage. Cold blood products chill enzymes and muscle together.
And because uterine circulation is pressure-passive, raising MAP with
crystalloid **raises the leak along with it.**

| Resuscitation strategy (ATN 0.88, 4-agent + TXA) | Blood loss | Side effects |
|---|---|---|
| Goal-directed 1:1 + fibrinogen + calcium + warming | **3489 mL** | iCa 1.00, 37.0 °C |
| Crystalloid-first, RBC only | 3953 mL | fibrinogen 1.37 g/L, 34.6 °C, pH 7.31, 6.2 units |
| No calcium replacement during massive transfusion | +922 mL | iCa 0.45 |
| No fluid warmer during massive transfusion | +2873 mL | 31.8 °C, 16 vs 11 units |

**A fluid warmer is a hemostatic drug.**

### 7. TXA's window is not 3 hours. It is the duration of the bleed

TXA does not build clot. It **protects** the clot the mechanical system has
managed to build. In a refractory patient bleeding for about 40 minutes:

| TXA timing | Blood loss | Benefit |
|---|---|---|
| Not given | 4099 mL | — |
| 5 minutes | 3400 mL | −699 |
| 15 minutes | 3645 mL | −454 |
| 30 minutes | 4030 mL | −69 |
| 45 / 90 / 180 minutes | 4099 mL | **0** |

The WOMAN trial's 3-hour boundary is a property **not of the molecule but of
that population's bleeding duration**. The model does not assume "10% loss of
benefit every 15 minutes" — it re-derives that from plasmin kinetics. (Peak
plasmin in the model: normal 1.0 · moderate 1.4 · severe treated 1.7 ·
**severe untreated 4.1** — hyperfibrinolysis is not a separate disease here but
**a consequence of flow.**)

### 8. Timing overwhelms the instrument

Catastrophic atony (ATN 0.94) + the full drug bundle + goal-directed
resuscitation:

| Additional intervention | Blood loss | Outcome | Oxygen debt |
|---|---|---|---|
| None | 10044 mL | **death at 80 minutes** | 1791 mL |
| Balloon tamponade at 20 minutes | **4528 mL** | survives | 152 |
| Balloon at 45 minutes | 7617 mL | survives | 1195 |
| Balloon at 75 minutes | 10156 mL | **death at 150 minutes** | 8260 |
| **Aortic compression 8-35 min → balloon at 35 min** | **2782 mL** | survives, **zero transfusion** | **0** |
| Uterine artery ligation at 50 minutes | 9624 mL | survives | 1527 |
| Hysterectomy at 60 minutes | 8451 mL | survives | 710 |
| Hysterectomy at 120 minutes | 10044 mL | dies | — |

**An early bridge intervention beats a late definitive operation.** Aortic
compression — needing no equipment and taking seconds — is the single best move
in this model.

### 9. Haemoglobin is the last thing to move

Whole blood leaves the body at a constant concentration, so the untreated
exsanguination arm **dies at Hb 11.2-11.4 g/dL.** A transfusion policy triggered
by haemoglobin is structurally late. That is why the model's massive
transfusion protocol is triggered by **ongoing loss** (>1500 mL with a leak
>25 mL/min) — which is how real protocols actually operate.

### Appendix: the other three T's are algebraically independent terms

The 4 T's are not a differential-diagnosis checklist but four terms invisible
to one another.

| Situation | Uterotonics alone | Definitive treatment |
|---|---|---|
| Cervical laceration, 260 mL/min | **death at 79 minutes**; even the full ladder still **dies at 80 minutes** | repair at 25 minutes → survives at 3794 mL / repair at 60 minutes → 8179 mL |
| Retained tissue, 35% | even the 4-agent ladder still **dies at 106 minutes** | manual removal at 30 minutes → survives at 4407 mL |
| Presenting fibrinogen 1.6 g/L | 2193 mL (nadir 0.89) | replacement to 2.5 g/L → 1677 mL |

**Failure to respond to uterotonics is not a dosing problem — it is diagnostic
information.**

---

## Model structure

41 ODE compartments, time unit = minutes, 70 kg term mother, starting
immediately after placental delivery.

| Block | State variables |
|---|---|
| Drug PK (10) | `OXY` `CBT` `MSD`/`MSC` `ERD`/`ERC` `PGD`/`PGC` `TX1`/`TX2` |
| Volume/flow (5) | `V` `CIV` `LOSS` `UO` `ODEBT` |
| Blood components (4) | `HBM` `FBM` `PLTM` `FCT` — all **masses**, so dilution follows automatically |
| Uterus (4) | `TONE` `OTR` `RET` `MECH` |
| Haemostasis (5) | `CLT` `CLTR` `MAT` `TPA` `PLS` |
| Internal environment (5) | `CIT` `CAI` `TMP` `LAC` `HCO` |
| Product counters (5) | `RBCU` `FFPU` `PLTU` `FGCU` `CAU` |
| Exposure integrals (3) | `TSEV` (MAP<50) `TAKI` (MAP<65) `XSEV` (Fib<2) |

Death is not adjudicated by an arbitrary blood-loss threshold but by whichever
comes first of a cumulative oxygen debt of **120 mL O₂/kg (= 8400 mL)** or a
45% loss of circulating volume.

### Calibrated parameters (6)

1. `UBF0` and the tone→patency exponent `NP` — normal third-stage blood loss
   of about 270 mL, haemostasis by about 5 minutes
2. The `ATN → (DRVI, CAP)` mapping — to reproduce the treatment ladder in §4
   above
3. `KFORM` `QS50` `CSEAL0` `KMAT` — so that a low-flow placental bed clots
   durably by 15-30 minutes and a high-flow placental bed never clots at all
4. `KDES`/`KREC` — OTR ≈ 0.39 after 8 hours of augmented labour
5. `KFDEG`/`KLYS`/`KTPA*` — fibrinogen <2 g/L and plasmin 3-4× in severe PPH
6. Citrate, calcium, and heat-capacity constants — ionised calcium
   0.5-0.9 mmol/L in obstetric massive transfusion, and about 2.8 °C/h of
   cooling per 100 mL/min of unwarmed 4 °C product

### Derived outputs (falsifiable)

The 2.8-minute reserve · reaching class III at 6 minutes · the ATN 0.60/0.63
knife-edge · the full treatment ladder · the futility of dose escalation versus
the benefit of switching receptors · fibrinogen crossing first and the OBS2
negative result · TXA's window equalling bleed duration · timing beating the
instrument · the fact that haemoglobin moves last.

### Reproduced trial results

- **WOMAN / Gayet-Ageron IPD** — benefit confined to early administration and
  vanishing with delay
- **OBS2 · FIB-PPH** — when fibrinogen is already up on the plateau
  (2.5-4 g/L), empirical 4 g dosing changes blood loss from 3489 to only
  3473 mL — a **16 mL** difference. The model reproduces this negative result
  and explains it: at that point `fFIB` is already 0.94 and clot formation is
  no longer fibrinogen-limited.
- **CHAMPION** — carbetocin and oxytocin's prophylactic effects sit close
  together (946 vs 1212 mL)
- **Hiippala** — the order of depletion on dilution (fibrinogen first)
- **Suarez et al.** — balloon tamponade's success rate and its dependence on
  timing

---

## How to run

```r
# Load the model and run scenarios
source("pph_mrgsolve_model.R")
print(S03_ladder())        # severity × treatment grid (Table 4 above)
print(S05_txa_timing())    # TXA timing series
print(S07_mechanical())    # mechanical/surgical escalation

# Interactive dashboard
shiny::runApp("pph_shiny_app.R")
```

### A note on where the numbers come from

Every figure cited in this document and in the `.R` comments was produced by an
**independent implementation of the same equations and parameters** (pure
Python RK4, dt = 0.01 min). The environment in which these files were written
had no R or mrgsolve, so `pph_mrgsolve_model.R` itself was not run here. If a
re-run gives different results, **trust the R output** and report the
discrepancy.

## The three assumptions that most deserve attack

1. **`NP = 4`** (the Poiseuille valve). Setting `NP = 2` erases the knife-edge
   and makes tone a forgiving continuous variable, and the model can no longer
   explain why "firm vs boggy" is a usable bedside distinction. `NP = 1` makes
   a normal third stage leak 3 L/h. No human measurement of the actual
   flow-tone relationship at the placental attachment site appears in the
   literature.
2. **`CSEAL0 = 0.35`** (the ceiling on clot over an unclamped vessel). Raising
   it to 1.0 lets severe atony achieve spontaneous haemostasis, and refractory
   PPH stops existing. Lowering it to 0 means coagulation replacement never
   helps at all. This model's entire claim that PPH is a mechanical disease
   hangs on this one parameter, so it is the first number that should be
   measured.
3. **`KREC` (t½ 2.6 h)** — if recovery is substantially faster, dose escalation
   should work and the advice to "switch receptors" is wrong. This is testable
   with myometrial sections from augmented versus unaugmented labour.

---

⚠️ **Disclaimer**: This is a qualitative/semi-quantitative QSP model for
educational and research purposes. It has not been independently validated or
certified and must not be used for actual clinical decision-making,
prescribing, or regulatory submission.
