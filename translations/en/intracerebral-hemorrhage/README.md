# Spontaneous Intracerebral Haemorrhage (ICH) — QSP Model

<p align="center">
  <a href="ich_qsp_model.svg"><img src="ich_qsp_model.png" width="720" alt="ICH QSP mechanistic map"></a>
</p>

A quantitative systems pharmacology (QSP) model of spontaneous
intracerebral haemorrhage, with hypertensive deep haemorrhage as its
archetype. It unites haematoma growth, haemostatic competence, the
Monro-Kellie volume balance, the three stages of perihaematomal oedema,
haem-iron chemical toxicity, neuroinflammation, and the PK/PD of
nicardipine · tranexamic acid · 4F-PCC · andexanet alfa · idarucizumab ·
deferoxamine · mannitol · and intraventricular alteplase into a single
system of differential equations.

| File | Contents |
|---|---|
| [`ich_qsp_model.dot`](ich_qsp_model.dot) | Mechanistic map source — 20 clusters, 219 nodes, 361 edges |
| [`ich_qsp_model.svg`](ich_qsp_model.svg) / [`ich_qsp_model.png`](ich_qsp_model.png) | Rendered map (Graphviz `dot`, PNG 150 dpi) |
| [`ich_mrgsolve_model.R`](../../../intracerebral-hemorrhage/ich_mrgsolve_model.R) | mrgsolve model — **60 ODE compartments**, 12 scenarios, including response surfaces and falsification switches |
| [`ich_shiny_app.R`](ich_shiny_app.R) | Shiny dashboard — 8 tabs |
| [`ich_references.md`](ich_references.md) | 112 PubMed-verified references (every PMID checked) |

---

## What This Model Claims (The three structural commitments)

The textbook ICH model ends at "the bigger the bleed, the worse the
outcome." That single picture cannot simultaneously satisfy what the
literature demands — facts such as **lowering blood pressure reduced
haematoma expansion but barely moved the 90-day outcome** (INTERACT2 ·
ATACH-2), **haematoma volume stayed the same yet the 90-day mRS moved**
(i-DEF), and **it was residual volume, not the surgery itself, that
mattered** (MISTIE III). So three structural choices were made.

### ① Haematoma Volume Is Not a State but an **Integral**

Nowhere in the model is there a term saying "drug X shrinks the
haematoma." There is a single bleeding flux, and it is **the product of
three factors**.

```
BLEEDR = KBLEED · NOPEN · (1 − CLOT) · max(0, MAP − PTISS)
              ↑        ↑          ↑
        open bleeding    haemostatic     driving
        points           incompetence    pressure
```

Antihypertensives act on the first factor, reversal agents and
tranexamic acid on the second, and mechanical avalanche and
procedure-related rebleeding on the third. Because these are **factors of
a single product, not additive terms**, it follows without any separate
assumption that a single agent's effect on volume is small (INTERACT2 ·
ATACH-2 · TICH-2 all reduced expansion but did not move the mRS), while
combination is synergistic (INTERACT3). Set any one factor to zero and the
entire flux vanishes — that is the falsifiable content of this structure.

Also, since `PTISS = ICP + KTAMP·VHEM`, **tamponade is endogenous**. The
larger the haematoma grows, the more it stops its own bleeding. This is
why untreated expansion is 20-38%, not 100%. Setting `KTAMP = 0` raises
expansion from 8.4 to 9.7 mL and abolishes the self-limiting behaviour.

### ② Outcome Is Determined by **Two Damage Channels Running on Different Clocks**

| Channel | Mechanism | Timescale | Reversed by evacuation? |
|---|---|---|---|
| **Mass effect** | `STRAIN` → tissue deformation → axonal displacement | Hours to days | Partially yes |
| **Haem-iron chemistry** | `HEME → HO-1 → Fe²⁺ → Fenton → LPO → ferroptosis` | Days to weeks | **No** |

Both channels erode `NEUR` (surviving neurons) and `WMI` (white matter
integrity), and appear in the mRS. This is why deferoxamine can move the
90-day mRS distribution **without touching haematoma volume at all** — in
the model, DFO has **not a single edge** into `VHEM`. And it is also why
MISTIE III showed benefit only in the ≤15 mL residual-volume subgroup:
evacuation acts only on Channel 1, and the mRS only moves once the
`STRTHR` threshold is crossed, while iron already liberated keeps flowing
down Channel 2.

### ③ MAP Appears **Twice, with Opposite Signs**

MAP is both the driving pressure of the bleed (lower is better) and,
simultaneously, the numerator of cerebral perfusion pressure
(`CPP = MAP − ICP`, where lower is worse — more so as mass effect and
oedema raise ICP, and as the lesion itself erodes autoregulatory capacity
`AUTOR`).

**No equation anywhere encodes a U-shaped curve.** And yet sweeping the
intensity of blood pressure lowering from 0 to 1 produces one.

| BP-lowering intensity | SBP nadir | Haematoma expansion (mL) | CPP<60 cumulative (h) | P(mRS 0-2) |
|---|---|---|---|---|
| 0.000 | 185 | 8.37 | 0.0 | 0.193 |
| 0.250 | 147 | 6.47 | 0.0 | 0.204 |
| 0.500 | 129 | 5.62 | 0.0 | 0.208 |
| **0.625** | **123** | **5.34** | **4.5** | **0.210 ← optimum** |
| 0.750 | 118 | 5.12 | 35.3 | 0.209 |
| 1.000 | 112 | **4.77 (minimum)** | 59.9 | 0.201 |

Haematoma expansion decreases **monotonically**, yet the outcome bends
near an SBP of 123 mmHg. This means the optimal SBP target is not a fixed
number but **state-dependent**. Scenario 4 was deliberately placed on that
harmful arm — of all the BP-lowering arms it has **the smallest bleed**
(4.68 mL) yet a 90-day outcome **worse** than Scenario 3's. Setting
`KISCH = 0` leaves the volume unchanged while removing only the harm,
confirming that the U is carried through perfusion.

---

## Model Structure (60 ODEs)

| Block | Compartments | Contents |
|---|---|---|
| Drug PK | 22 | Nicardipine (2 compartments + effect site), clevidipine, labetalol, TXA (2 compartments), 4F-PCC, vitamin K recovery, apixaban (2 compartments), andexanet + complex, dabigatran, idarucizumab + complex, deferoxamine, feroxamine, mannitol, intraventricular alteplase |
| Haemostasis | 6 | Prothrombin activity, fibrinogen, platelet function, plasmin, clot competence, open bleeding points |
| Haemodynamics | 3 | SBP, sympathetic tone (including the Cushing reflex), autoregulatory capacity |
| Intracranial | 5 | CSF volume, haematoma, intraventricular blood, early oedema, late oedema |
| Barrier · oedema mediators | 3 | Tissue thrombin, BBB permeability, MMP-9 |
| Haem-iron | 4 | Free haem, free Fe²⁺, ferritin sequestration capacity, lipid peroxidation |
| Inflammation | 5 | MG1/MG2 microglia, neutrophils, IL-6, IL-10 |
| Tissue · clinical | 6 | Surviving neurons, white matter integrity, plasticity reserve, NIHSS, temperature, glucose |
| Cumulative exposure | 6 | SBP>140 AUC, CPP<60 time, ICP>20 time, Fe²⁺ AUC, IL-6 AUC, total extravasated blood |

**ICP is not a state variable.** It is the residual of the Monro-Kellie
volume balance.

```
VADD = VHEM + VIVH + OEDE + OEDL + (VCSF − VCSF0) − VCOMPE      (0 if negative)
ICP  = ICP0 · exp(EELAST · VADD)
```

CSF absorption follows the Davson equation (absorption = (ICP − Pss) ×
conductance), and in the intact state is exactly balanced against
production of 21 mL/h. Because the displaceable CSF reserve is finite at
18 mL, the pressure-volume curve becomes **biphasic** rather than
exponential — flat while reserve remains, and steep once it is exhausted.
This is the actual physiology of compensatory exhaustion.

---

## 12 Scenarios and Verification Results

The system was **verified by actually integrating it in mrgsolve 2.0.1 / R
4.3.3.** The numbers below are run outputs, not claims.

| # | Scenario | Expansion (mL) | SBP nadir | Peak PHE (mL) | CPP<60 (h) | Fe AUC(14d) | NIHSS 90d | P(mRS 0-2) |
|---|---|---|---|---|---|---|---|---|
| 1 | Untreated natural history | 8.37 (27.9%) | 185 | 26.2 | 0.0 | 71.0 | 13.9 | 0.195 |
| 2 | Guideline BP lowering (<180) | 7.51 | 170 | 24.9 | 0.0 | 70.1 | 13.7 | 0.200 |
| 3 | Intensive BP lowering (<140, within 1 hour) | 6.24 | 145 | 23.1 | 0.0 | 68.8 | 13.5 | 0.207 |
| 4 | **Excessive BP lowering (harm arm)** | **4.68 (minimum)** | 109 | 20.9 | **64.3** | 67.1 | **13.6** | **0.204** |
| 5 | Tranexamic acid alone (TICH-2) | 7.18 | 185 | 24.6 | 0.0 | 69.9 | 13.7 | 0.201 |
| 6 | **Care bundle (INTERACT3)** | 5.44 | 145 | 20.6 | 0.0 | 68.0 | **11.4** | **0.295** |
| 7 | Warfarin ICH, no reversal | **14.59 (48.6%)** | 163 | 19.2 | 0.0 | 77.0 | 13.7 | 0.200 |
| 8 | Warfarin + 4F-PCC + vitamin K | 8.66 | 138 | 21.5 | 0.0 | 71.0 | 13.6 | 0.207 |
| 9 | Apixaban + andexanet | 6.84 | 137 | 22.8 | 0.0 | 69.2 | 13.5 | 0.207 |
| 10 | **Deferoxamine (i-DEF)** | **6.24 (identical to #3)** | 145 | 22.3 | 0.0 | **48.1 (−30%)** | 13.3 | 0.215 |
| 11 | **MISTIE III evacuation** | 6.24 | 145 | **12.4** | 0.0 | 50.6 | 11.9 | 0.273 |
| 12 | Large ICH + IVH + EVD + alteplase | 3.84 | 149 | 27.3 | 0.0 | 79.9 | 14.6 | 0.173 |

Two rows worth noting:

- **Scenario 10 vs 3** — the 24-hour haematoma volume is **identical to
  the decimal** (6.24 mL), yet the 14-day iron-exposure AUC alone falls
  from 68.8 to 48.1, and the outcome improves. Setting `KCHEL = 0` returns
  iron AUC to 68.8 and removes the effect entirely. This is how i-DEF's
  "mRS shift with no volume change" is structurally reproduced in the
  model.
- **Scenario 4 vs 3** — the excessive-lowering arm bleeds less, 4.68 vs
  6.24 mL, yet the 90-day NIHSS is worse, 13.6 vs 13.5. The only item
  making the difference is the 64.3 hours of CPP<60 exposure.

### Comparison Against Literature Anchors

| Anchor | Literature | This model |
|---|---|---|
| Untreated expansion | 20-38% | +8.4 mL (27.9%) |
| Peak perihaematomal oedema | 3-5 days | 26.2 mL, 4.3 days |
| Intensive vs guideline lowering | ~15-25% reduction in expansion | 6.2 vs 7.5 mL (−25%) |
| TXA (TICH-2) | ~10-15% reduction in expansion | 7.2 mL (−14%) |
| Warfarin INR 3.0, no reversal | ~2-fold | 1.74-fold |
| Warfarin + 4F-PCC | Near-normalisation | 8.66 (no anticoagulation: 8.37) |
| Andexanet | >90% reduction in anti-Xa, then rebound | 171 → 1.0 ng/mL, rebounding to 7.2 at 12 hours |

---

## Falsification Switches

Every structural claim has a parameter that kills it. Checked with the
`falsify()` function.

| Switch | Result | What it proves |
|---|---|---|
| `KBLEED = 0` | Expansion 8.37 → −0.77 mL | All volume change comes from this single flux |
| `KTAMP = 0` | Expansion 8.37 → 9.73 mL | Tamponade is the source of self-limitation |
| `KAVAL = 0` | Expansion 8.37 → 7.79 mL | Mechanical avalanche creates the early-time sensitivity |
| `KCHEL = 0` | Iron AUC 48.1 → 68.8 (volume unchanged) | DFO's effect passes **entirely** through the iron channel |
| `KFERRO2 = KOEDL = 0` | NIHSS 13.3 → 10.0 (volume unchanged) | The chemical channel really does carry the outcome |
| `KISCH = 0` | Scenario 4's harm disappears (volume unchanged) | The U-shaped curve is carried through perfusion |
| `KOFFA = 0` | Andexanet rebound disappears | The rebound is an emergent property of dissociation + redistribution, not a separate term |

---

## Usage

```r
# model · scenarios · dosing helpers
source("ich_mrgsolve_model.R")

sim <- run_all_scenarios()          # 12 scenarios (tens of seconds)
summarise_scenarios(sim)            # summarised using the metrics the trials actually reported

# the emergent U-shaped curve
bp_response_surface(end = 1440)

# falsification: does killing the iron channel remove the deferoxamine effect
falsify("10_deferoxamine_idef", list(KCHEL = 0))

# interactive dashboard (8 tabs)
shiny::runApp("ich_shiny_app.R")
```

Required packages: `mrgsolve`, `dplyr` (model), plus `shiny`, `ggplot2`,
`tidyr` (dashboard).

---

## Cases Where the Literature Forced a Model Fix During Development

A QSP model's value lies less in what it got right than in **what was
shown to be wrong and fixed.** Actually integrating the draft revealed the
following, which were then corrected.

1. **CSF diverged to negative values.** There was no pressure-driven
   absorption term to balance production (21 mL/h). Adding the Davson
   equation and a finite displaceable reserve (18 mL) made the
   pressure-volume curve biphasic, which is better physiology than
   originally intended.
2. **IL-6 ran away.** Writing the `IL6`↔`NEU` feedback linearly gives a
   loop gain exceeding 1, so the cytokine diverges to 18-fold by 24 hours.
   Both sides had to be changed to saturating functions before
   inflammation became a finite, self-resolving pulse.
3. **The blood pressure effect was calibrated larger than the literature
   supports.** The draft's comments recorded a target of "35-40%
   reduction in expansion," whereas INTERACT2's mean growth was 4.5 vs
   5.3 mL and ATACH-2's expansion was 18.9% vs 24.4% — actually 15-25%.
   Because this amounted to **calibrating to reproduce a result the
   trials never obtained**, the target was pulled back toward the
   literature.
4. **Bleeding ended before the drug could arrive.** At `KSEAL = 0.80`,
   bleeding terminated within 1.5 hours, leaving antihypertensives no time
   to act, and the BP-lowering arms were effectively indistinguishable.
   Extending the bleeding window to 6-12 hours (consistent with the
   literature) made time-to-arrival meaningful.
5. **Nicardipine was too potent.** At `EC5NIC = 0.035 mg/L`, 3 mg/h
   already saturated Emax, making the guideline, intensive, and excessive
   arms pharmacologically identical.
6. **Warfarin bleeding came out at 3.5-fold** (the literature says about
   2-fold). Writing `CLOT` as a single multiplicative term let INR erode
   the platelet arm along with it. Separating the platelet arm (`WPLT`)
   from the coagulation arm (`WFIB`) brought it to 1.74-fold, and as a
   side effect it became consistent with PATCH (platelet transfusion is
   harmful in antiplatelet-associated ICH), giving the model a reason why
   platelet transfusion and DDAVP are ineffective for VKA-associated
   bleeding.
7. **The U-shaped curve did not emerge.** At `CBFISC = 0.55`, the
   perfusion channel never activated, giving a monotonic "the lower, the
   better." Because perihaematomal tissue is already hypoperfused
   (Zazulia), the damage-onset threshold **must be higher** than the
   classical infarction threshold. Raising it to 0.70 produced the U.
8. **Mortality was scored on the recovered state.** The ICH score is
   designed as a baseline score, but it was being computed from the
   90-day state. Changed to the 24-hour score plus cumulative
   intracranial pressure exposure.

---

## Limitations

- **This is not for individual patient prediction.** It has not been
  fitted to patient-level data. "Verification" here means the equations
  integrate stably and reproduce trial-level direction and approximate
  magnitude — it does not mean fitness for clinical use.
- **Recurrence and long-term secondary prevention are not in the ODEs.**
  They lie outside the 90-day horizon and exist only on the map (cluster
  17).
- **The self-fulfilling prophecy of DNR/withdrawal of care** (Becker
  2001) confounds the interpretation of every ICH trial, but the model
  includes it only as a map node, not as an ODE implementation. The
  model's mortality reflects the biological trajectory alone.
- The distinction between CAA and hypertensive pathology exists only at
  the level of the `LOC`·`SVD` parameters; there is no Aβ dynamics.
- The perihaematomal neuronal survival rate (`NEUR`) falling to 0.12-0.15
  by 90 days is on the pessimistic side. The perihaematomal rim does
  indeed mostly die, but the NIHSS mapping weights (`WNEUR`·`WWM`·`KDEF`)
  would be worth refitting against patient data.

---

## ⚠️ Disclaimer

This is a QSP model for educational and research purposes. It was built
from public literature but has not been independently validated or
certified, and **must not be used for actual clinical decision-making,
prescribing, or regulatory submission.**
