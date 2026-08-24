# Controlled Ovarian Stimulation (COS) — QSP Model
## Controlled Ovarian Stimulation for IVF · Quantitative Systems Pharmacology

<p align="center">
  <a href="cos_qsp_model_en.svg">
    <img src="cos_qsp_model_en.png" width="900" alt="COS QSP mechanistic map">
  </a><br>
  <sub><a href="cos_qsp_model_en.svg">View full-resolution SVG</a> · 169 nodes · 215 edges · 17 clusters (16 mechanisms + 1 legend)</sub>
</p>

---

## The claim this model is built to make

> **Yield and OHSS are not two different variables but two readings of the
> same state variable — growing granulosa cell mass. So they cannot be
> separated by FSH dose. Separation comes only from the shape of the
> trigger — because oocyte maturation needs only the leading edge of the
> LHCGR signal, while OHSS is proportional to its area.**

As with the other models in this repository, what matters is not "what was
put in" but **"what came out despite not being put in."** The following are
**absent** from this model's list of parameters:

- The number of growing follicles, the number of oocytes retrieved
- Days of stimulation
- Any coefficient stating "hCG is riskier than an agonist"
- An OHSS severity scale, risk score, or threshold
- A classification stating "PCOS is a high-responder"
- A regression equation stating "AMH predicts oocyte number"

And yet every result below comes out anyway. What was put in is only **one
mass, one ligand pool, and a set of different kernels that read them**.

| What was put into the model (structure) | What the model computed (result) |
|---|---|
| A per-follicle FSH threshold distribution (log-normal T50 9 IU/L, σ 0.45) + E2/inhibin FSH feedback | **Single ovulation** in a natural cycle (1.2 follicles), 18 growing simultaneously in a stimulated cycle — the same equations |
| The threshold **falls** with size, and **switches over** to LH support past 11 mm | Selection of the dominant follicle and atresia of the rest (FSH 7.8 → 2.6 IU/L) |
| E2, P4, oocytes, and VEGF are all read from a single granulosa-cell mass | 3× the dose → oocytes **+15%**, ascites **+125%** (dose is not a way to separate them) |
| The P4 leak term is **linear** in mass (no saturation) | Trigger-day P4 becomes a follicle counter: 0.66 (AFC 12) / 1.31 (25) / 1.49 (32) ng/mL |
| The LHCGR signal is read through three kernels, K = 2 / 25 / 40 IU/L | Agonist trigger: the same maturation rate (0.77 vs 0.78), a **43-fold** reduction in ascites |
| Maturation = leading edge + an autonomous clock / rupture = clock × integrated signal | A 34–38-hour retrieval window is derived (61% loss at 44 hours) |
| Progesterone suppresses GnRH pulsatility | Multiple corpora lutea produce a feedback that **cuts off its own support** → luteal-phase deficiency |
| AMH raises the threshold, and LH tone arrests follicles before they acquire LHCGR | PCOS's 6–9 mm arrest (the polycystic appearance) and hyper-response to stimulation, together |
| Age → ploidy only | From age 28 to 44, oocyte number stays **unchanged** at 8.0, while CLBR falls from 76% to 21% |

Every figure was actually computed by
[`cos_reference_check_en.py`](cos_reference_check_en.py),
and the full output is in
[`cos_reference_output.txt`](cos_reference_output.txt).

---

## Deliverables

| File | Contents |
|------|------|
| [`cos_qsp_model_en.dot`](cos_qsp_model_en.dot) / [`.svg`](cos_qsp_model_en.svg) / [`.png`](cos_qsp_model_en.png) | Mechanistic map — 169 nodes · 215 edges · 17 clusters |
| [`cos_mrgsolve_model_en.R`](cos_mrgsolve_model_en.R) | **65-ODE** mrgsolve model · 22 scenarios · 4 sweeps |
| [`cos_shiny_app_en.R`](cos_shiny_app_en.R) | 11-tab Shiny dashboard |
| [`cos_references_en.md`](cos_references_en.md) | **87 references** (every PMID individually verified, 48 used directly for parameter calibration) |
| [`cos_reference_check_en.py`](cos_reference_check_en.py) | Independent Python/scipy re-implementation of the same equations (for validation) |
| [`cos_reference_output.txt`](cos_reference_output.txt) | Full output of the validation script (the source of every number in this README) |

---

## 1. Model structure (65 ODEs)

| Compartment group | Count | Contents |
|--------|------|------|
| Follicles | 30 | 10 threshold-quantile slots × (diameter D · granulosa mass G · viability V) |
| Drug PK | 12 | rFSH (depot/central), corifollitropin, GnRH antagonist, hCG (depot/central), GnRH agonist (depot/central), cabergoline, letrozole, vaginal P4 |
| Endocrine | 11 | Pituitary LH storage pool · LH · endogenous FSH · E2 · inhibin B · AMH · P4 · luteal mass · surge readiness · GnRH receptor availability |
| OHSS | 5 | VEGF · vascular permeability · third-space fluid · plasma volume · pregnancy hCG |
| Accumulators | 7 | Maturation-kernel exposure · autonomous maturation clock · rupture integral · oocyte count · MII count · E2 AUC · VEGF kernel area |

### Why follicles are discretised into "slots"

Instead of tracking the follicle cohort individually, **equal-probability
deciles of the FSH threshold distribution** were taken as slots, with each
slot representing `AFC/10` follicles. The consequences of this choice:

- Inter-patient variation enters through **AFC alone** (the number of slots
  and the equations stay fixed).
- AMH is not a separate input but a **computed output** of the
  small-follicle pool (AFC 12 → 2.38, AFC 32 → 6.34 ng/mL).
- "How many grow" is determined by the length of time `FSHTOT` stays above
  each slot's threshold — in other words, it is an **output**.

### One ligand, three kernels (the core of this model)

LH, hCG, and pregnancy hCG all use the same receptor. So the ligand is kept
as a single one (`LHEQ`), and the only thing added is the fact that **each
downstream process needs a different signal**.

| Process | Kernel | What it needs |
|------|------|-----------|
| Thecal androgen · luteal support | `LHEQ/(LHEQ+2)` | Baseline pulsatility is enough |
| Resumption of meiosis (maturation) | `Hill(LHEQ, 25, 6)` → **commits** at 0.06 d exposure → then a 0.85 d **autonomous clock** | Only the **leading edge** of the surge |
| Follicle-wall rupture (ovulation) | The above clock reaches 1.55 d **and** the signal integral ≥ 0.90 d | The leading edge **and** a sustained signal |
| VEGF induction | The **time integral** of `Hill(LHEQ, 40, 3)` | The **area**, over several days |

This single asymmetry derives the whole of the modern trigger literature
(Table B below).

---

## 2. Verification

Because R/mrgsolve could not be run in this environment,
**all 65 ODEs were independently re-implemented in Python/scipy**,
integrating 22 scenarios × 34 days and 6 parameter sweeps. This process
uncovered and fixed **8 structural defects**.

| # | Defect uncovered | Symptom | Fix |
|---|------------|------|------|
| i | The meiotic-commitment integral was driven by **total LHCGR occupancy** | Baseline LH of 5 IU/L committed oocytes on cycle day 1, so 4 mm follicles all ovulated by day 2 | Separating the kernel for each downstream process — this became the model's **central claim** |
| ii | The E2 positive-feedback switch's time constant was too fast (KR 3.0) | A partial surge at E2 148 pg/mL (stimulation day 2) → every follicle arrested at 7.1 mm | KR 1.20 · KROFF 0.45 · a dominant-size gate · the antagonist **blocking entry into surge mode** |
| iii | The pituitary LH storage pool had **infinite capacity** | The pool grew without bound under antagonist dosing → **the antagonist's steady-state LH suppression effect was zero** | A saturable pool (SLHMAX 12 000) → the agonist flare size falls out as a result |
| iv | The oocyte-counting window and the follicle-collapse window **overlapped** | Counting occurred in follicles that were already emptying → at AFC 12, oocyte count came out as **2.4** (correct value 8.0) | Separated into a counting window (1.2 h) followed by a collapse window (7.2 h) |
| v | The agonist bolus was dosed in the **morning** of trigger day | Given 9.6 h earlier than the actual trigger → the dual-trigger arm ovulated before retrieval (14.6 → 2.8 oocytes) | Separated into events integrated at the exact clock time |
| vi | Luteal P4 output was **linear** in mass | 2.3 ng/mL in a natural luteal phase (observed 10–20), 173 ng/mL in an IVF luteal phase with 15 corpora lutea (observed 25–60) | A saturating term `CL/(1+CL/1.46)` → reproduces both simultaneously |
| vii | Progesterone's GnRH feedback was **missing** | Under any trigger, LH recovered to 8–10 IU/L within 2 days → the corpus luteum was always rescued (agonist deficiency only −13%) | A `KP4LH` brake → luteal mass cuts off its own support, and dual-trigger rescue is reproduced |
| viii | The threshold distribution was centred at 6.5 IU/L | Endogenous FSH alone grew 6 slots → **the natural cycle became a triple ovulation**, with day-5 E2 of 300 pg/mL | Re-centred at 9.0 IU/L (σ 0.45) — the stimulated-arm calibration was left unchanged |

---

## 3. Scenario results (22 arms)

`dstim` days of stimulation · `E2trg` trigger-day E2 (pg/mL) · `P4trg`
trigger-day P4 (ng/mL) · `foll` follicles ≥11 mm · `OPU` oocytes retrieved ·
`MII` mature oocytes · `AUCveg` VEGF kernel area (days) · `Hct` peak
haematocrit (%) · `asc` peak ascites (L)

| Scenario | dstim | E2trg | P4trg | foll | OPU | MII | AUCveg | Hct | asc | OHSS |
|---|---|---|---|---|---|---|---|---|---|---|
| 0 Natural cycle (no drugs) | – | 430* | – | – | – | – | 0.49 | 38.1 | 0.00 | none |
| 1 Antagonist + rFSH 150 + hCG 10000 | 11.4 | 1615 | 0.66 | 9.3 | 8.0 | 6.3 | 10.84 | 39.3 | 0.27 | mild |
| 2 Individualised low dose (rFSH 100) | 11.4 | 1336 | 0.55 | 7.5 | 6.8 | 5.2 | 10.84 | 38.6 | 0.11 | mild |
| 3 High responder (AFC 25) + hCG | 12.4 | 3594 | 1.31 | 18.3 | 15.7 | 12.5 | 9.87 | 45.4 | 1.56 | **severe** |
| 4 High responder + **agonist trigger + freeze-all** | 12.4 | 3594 | 1.31 | 18.3 | 14.7 | 11.5 | **0.42** | 38.2 | **0.02** | none |
| 5 High responder + rFSH **75** + hCG (dose reduction) | 12.4 | 2416 | 0.89 | 12.5 | 11.3 | 8.7 | 9.88 | 41.7 | 0.78 | moderate |
| 6 High responder + hCG + cabergoline | 12.4 | 3594 | 1.31 | 18.3 | 15.7 | 12.5 | 9.87 | 43.9 | 1.29 | moderate |
| 7 High responder + coasting (stopped CD8) | 9.4 | 1503 | 0.69 | 14.7 | 12.4 | 8.1 | 12.89 | 39.0 | 0.20 | mild |
| 8 Low responder (AFC 5) + rFSH 300 | 12.4 | 889 | 0.35 | 4.3 | 3.7 | 3.1 | 9.86 | 38.1 | 0.01 | none |
| 9 Low responder + hp-hMG LH 75 | 12.4 | 955 | 0.33 | 4.2 | 3.6 | 2.9 | 9.87 | 38.1 | 0.01 | none |
| 10 Corifollitropin alfa 150 µg single dose | 11.4 | 2088 | 0.84 | 10.8 | 9.4 | 7.6 | 10.83 | 41.3 | 0.71 | moderate |
| 11 **No antagonist** (premature LH surge) | 16.4 | 139 | 18.75 | 0 | 0.2 | 0.0 | 6.20 | 38.1 | 0.00 | cancelled |
| 12 Dual trigger (agonist + hCG 1500) | 12.4 | 3594 | 1.31 | 18.3 | 14.7 | 11.5 | 5.85 | 38.3 | 0.05 | none |
| 13 Agonist trigger + **fresh transfer** | 12.4 | 3594 | 1.31 | 18.3 | 14.7 | 11.5 | 5.17 | 38.2 | 0.02 | none |
| 14 hCG, retrieval delayed to **40 hours** | 11.4 | 1615 | 0.66 | 9.3 | 7.6 | 5.3 | 10.66 | 39.9 | 0.40 | moderate |
| 15 Age 41, AFC 8, rFSH 300 | 11.4 | 1305 | 0.53 | 7.0 | 6.1 | 4.9 | 10.84 | 38.5 | 0.10 | none |
| 16 High responder + hCG + freeze-all | 12.4 | 3594 | 1.31 | 18.3 | 15.7 | 12.5 | 5.17 | 45.5 | 1.57 | **severe** |
| 17 High responder + PPOS + agonist trigger | 12.4 | 3866 | **5.22** | 18.2 | 14.7 | 11.4 | 0.36 | 38.1 | 0.01 | none |
| 18 AFC 32 + hCG 10000 + fresh transfer | 12.4 | 4026 | 1.49 | 20.7 | 18.5 | 14.4 | 9.86 | 46.2 | 1.74 | **severe** |
| 19 AFC 32 + **agonist trigger + freeze-all** | 12.4 | 4026 | 1.49 | 20.7 | 17.1 | 13.2 | **0.42** | 38.3 | **0.04** | none |
| 20 AFC 32 + hCG + cabergoline + freeze-all | 12.4 | 4026 | 1.49 | 20.7 | 18.5 | 14.4 | 5.16 | 44.5 | 1.44 | moderate |
| 21 AFC 32 + hCG + freeze-all | 12.4 | 4026 | 1.49 | 20.7 | 18.5 | 14.4 | 5.16 | 46.2 | 1.74 | **severe** |

\* Peak E2 of the natural cycle (no trigger). LH surge 55.7 IU/L (day
11.2), ovulation on day 12.5, mid-luteal P4 12.5 ng/mL, luteal-phase LH
2.7–3.7 IU/L.

### How to read this — three pairwise comparisons

- **3 vs 4** (same patient, same stimulation, only the trigger differs):
  oocytes 15.7 → 14.7 (**−6%**), ascites 1.56 → 0.02 L (**−99%**), VEGF area
  9.87 → 0.42 d. The maturation rate is 0.80 → 0.78.
- **3 vs 5** (same trigger, dose halved): oocytes 15.7 → 11.3 (**−28%**),
  ascites 1.56 → 0.78 L (−50%), **still moderate OHSS**. Dose only ever
  moves both values in the same direction.
- **18 vs 21** (fresh vs freeze-all, same hCG): early OHSS is **identical**
  (1.74 L). What freeze-all removes is only the **late OHSS** produced by
  pregnancy hCG.

---

## 4. Sweep A — does FSH dose separate yield from risk?

| Patient / dose | dstim | E2trg | Oocytes | MII | Euploid | Hct | Ascites (L) | OHSS |
|---|---|---|---|---|---|---|---|---|
| AFC 12 / 75 IU | 11.4 | 1121 | 5.8 | 4.4 | 1.22 | 38.3 | 0.05 | none |
| AFC 12 / 112 IU | 11.4 | 1420 | 7.1 | 5.5 | 1.54 | 38.8 | 0.15 | mild |
| AFC 12 / **150 IU** | 11.4 | 1615 | 8.0 | 6.3 | 1.75 | 39.3 | 0.28 | mild |
| AFC 12 / 225 IU | 11.4 | 1830 | 8.9 | 7.1 | 1.96 | 40.2 | 0.47 | moderate |
| AFC 12 / 300 IU | 11.4 | 1932 | 9.1 | 7.3 | 2.04 | 40.6 | 0.56 | moderate |
| AFC 12 / 450 IU | 11.4 | 2008 | 9.2 | 7.5 | 2.08 | 41.0 | 0.63 | moderate |
| AFC 25 / 75 IU | 12.4 | 2416 | 11.3 | 8.7 | 2.51 | 41.7 | 0.78 | moderate |
| AFC 25 / **150 IU** | 12.4 | 3594 | 15.7 | 12.5 | 3.61 | 45.4 | 1.56 | severe |
| AFC 25 / 450 IU | 12.4 | 4476 | 18.3 | 14.9 | 4.32 | 46.7 | 1.88 | severe |

**150 → 450 IU (3×): oocytes +15%, ascites +125%.** This is because
oocytes saturate at the cohort size (AFC), while VEGF area does not
saturate. The model has no "dose ceiling" or "saturation parameter" — once
a finite cohort has been fully recruited, there is simply nothing left to
give. This is the same structure as the flattening of the clinical
dose–response observed in practice (references 67 · 71).

## 5. Sweep B — one ligand, three kernels (AFC 32)

| Trigger | Peak ligand (IU/L) | VEGF area (d) | Maturation area (d) | Oocytes | MII | Maturation rate | Hct | Ascites (L) | Luteal day-7 P4 |
|---|---|---|---|---|---|---|---|---|---|
| hCG 10 000 IU | 199 | **5.16** | 6.21 | 18.5 | 14.4 | 0.78 | 46.2 | **1.74** | 37.4 |
| hCG 5 000 IU | 101 | 3.55 | 4.75 | 18.5 | 14.3 | 0.77 | 45.0 | 1.36 | 37.1 |
| hCG 2 500 IU | 51 | 1.76 | 3.12 | 18.5 | 13.7 | 0.74 | 41.0 | 0.57 | 36.5 |
| hCG 1 500 IU | 32 | 0.71 | 1.59 | 18.5 | 12.2 | **0.66** | 38.2 | 0.02 | 35.8 |
| Dual (agonist + hCG 1 500) | 159 | 1.13 | 1.97 | 17.2 | 13.2 | 0.77 | 38.5 | 0.09 | 34.8 |
| **Agonist 0.2 mg alone** | 152 | **0.41** | 0.52 | 17.1 | 13.2 | 0.77 | 38.3 | **0.04** | **22.2** |

This table is the whole model.

1. **Reducing the hCG dose and switching to an agonist are not the same
   thing.** hCG 1500 IU brings ascites down to 0.02 L, but its peak of
   only 32 IU/L barely clears the maturation kernel (K = 25, n = 6),
   dropping the maturation rate to 0.66. The agonist reaches a peak of
   152 IU/L (ample for maturation) with an area of 0.41 d (no OHSS). **A
   signal with a high leading edge and a narrow area** — that is what an
   agonist trigger actually is.
2. **The luteal phase pays a price.** Luteal day-7 P4 with the agonist
   alone is 22.2 vs 37.4 ng/mL with hCG (−41%). The dual trigger brings it
   back to 34.8. The model reproduces this direction, but at a smaller
   magnitude than clinical observations (references 29 · 30) — the
   remaining mechanism is presumed to be a qualitative deficit in
   luteinisation itself, and is recorded as a limitation.

## 6. Sweep C — the 34–38-hour retrieval window is the gap between two clocks

| Trigger→retrieval | Oocytes | MII | Maturation rate | Follicle loss |
|---|---|---|---|---|
| 28 h | 7.9 | 5.2 | 0.66 | 1% |
| 32 h | 8.0 | 5.9 | 0.75 | 0% |
| **34 h** | 8.0 | 6.1 | 0.77 | 0% |
| **36 h** | 8.0 | 6.3 | **0.79** | 0% |
| **38 h** | 8.0 | 6.3 | 0.79 | 0% |
| 40 h | 7.6 | 5.3 | 0.70 | 5% |
| 44 h | 3.1 | 0.8 | 0.24 | **61%** |

The leading edge is set by the autonomous maturation clock (`TMAT` 0.85
d), and the trailing edge by the rupture integral (`RUPX` 0.90 d × clock
1.55 d). There is no constant "36 hours" anywhere in the model.

## 7. Sweep D — ovarian reserve and cumulative live-birth rate

| AFC | AMH (computed) | Oocytes | MII | Blastocysts | Euploid | CLBR (%) | OHSS |
|---|---|---|---|---|---|---|---|
| 4 | 0.79 | 3.1 | 2.5 | 1.0 | 0.69 | 39.9 | none |
| 6 | 1.19 | 4.4 | 3.5 | 1.4 | 0.97 | 50.9 | none |
| 9 | 1.78 | 6.2 | 5.0 | 2.0 | 1.38 | 63.7 | none |
| 12 | 2.38 | 8.0 | 6.3 | 2.5 | 1.75 | 72.2 | mild |
| 16 | 3.17 | 10.1 | 7.8 | 3.1 | 2.17 | 79.6 | moderate |
| 20 | 3.96 | 11.9 | 9.1 | 3.6 | 2.52 | 84.3 | moderate |
| 25 | 4.95 | 13.8 | 10.4 | 4.1 | 2.89 | 88.0 | moderate |
| 32 | 6.34 | 16.1 | 11.9 | 4.7 | 3.30 | 91.1 | moderate |

The saturation of CLBR (+32 percentage points from AFC 4 to 12, +7
percentage points from AFC 20 to 32) is a consequence of the binomial
chain `1 − (1 − 0.52)^n`, and was not forced. AMH, too, is a value
computed from the small-follicle pool, not a regression equation.

## 8. Sweeps E · F — antagonist start day and age

| Antagonist start | Peak LH | Premature surge | dstim | Oocytes | MII |
|---|---|---|---|---|---|
| None | 74.0 | **YES** | 16.4 | 0.2 | 0.0 |
| CD4 | 6.0 | no | 11.4 | 7.9 | 6.2 |
| CD5 | 6.0 | no | 11.4 | 8.0 | 6.3 |
| CD6 | 6.0 | no | 11.4 | 8.0 | 6.3 |
| CD7 | 6.0 | no | 11.4 | 8.0 | 6.3 |
| CD8 | 60.0 | **YES** | 16.4 | 0.3 | 0.0 |

CD4–CD7 are equivalent, and CD8 is already too late — the model says the
surge window closes when the dominant follicle reaches 14 mm.

| Age | Oocytes | MII | Blastocysts | Euploid | CLBR (%) |
|---|---|---|---|---|---|
| 28 | 8.0 | 6.3 | 2.5 | 1.96 | 76.3 |
| 32 | 8.0 | 6.3 | 2.5 | 1.75 | 72.2 |
| 36 | 8.0 | 6.3 | 2.3 | 1.26 | 60.4 |
| 39 | 8.0 | 6.3 | 2.2 | 0.84 | 45.9 |
| 42 | 8.0 | 6.3 | 2.0 | 0.48 | 29.9 |
| 44 | 8.0 | 6.3 | 1.9 | 0.31 | 20.6 |

With AFC held fixed, **oocyte number is independent of age**, while CLBR
alone falls 3.7-fold. The clinical common sense that stimulation is a
matter of oocyte number and age is a matter of ploidy falls directly out
of which term in the chain depends on age.

---

## 9. How to reproduce

```bash
# Mechanistic map
dot -Tsvg cos_qsp_model_en.dot -o cos_qsp_model_en.svg
dot -Tpng -Gdpi=150 cos_qsp_model_en.dot -o cos_qsp_model_en.png

# Validation (Python/scipy — 22 scenarios + 6 sweeps, about 15 minutes)
python3 cos_reference_check_en.py           # full
python3 cos_reference_check_en.py --brief    # summary table only

# mrgsolve model and Shiny app (R required)
R -e "source('cos_mrgsolve_model_en.R')"
R -e "shiny::runApp('cos_shiny_app_en.R', port = 8080)"
```

---

## 10. Known limitations and falsifiable predictions

| Item | Model | Observed | Interpretation |
|---|---|---|---|
| Natural-cycle dominant follicle diameter | 16.9 mm | 20–22 mm | The E2 that fires the surge is reached before the diameter is. Short by 3–5 mm |
| Agonist-trigger luteal-phase deficiency | −41% (P4 22.2 vs 37.4) | Larger | The model cannot capture the qualitative deficit in luteinisation itself |
| Cabergoline | Ascites −17% | Meta-analysis RR 0.38 | The model is conservative |
| Luteal-phase letrozole | E2 −79%, ascites change <4% | Uncertain | **A falsifiable prediction**: the structure in which E2 is a marker and VEGF is the mediator |
| Spontaneous ovulation after an agonist trigger | Almost none even at 44 hours | Uncertain | **A falsifiable prediction**: the rupture integral is never reached |
| Coasting | OHSS reduction only comes with a −35% drop in oocytes | Uncertain | **A falsifiable prediction**: a direct consequence of the single-mass structure |
| Absolute inhibin B | Natural-cycle peak of 51 pg/mL | 80–150 | The small-follicle weighting term (`WSMIB`) is too low |
| LH supplementation (hp-hMG) | Per-follicle E2 +7%, oocytes unchanged | Meta-analyses also show no increase in oocytes | Consistent |

This model does not capture the **molecular mechanisms of endometrial gene
expression and the receptivity window**, **sperm/fertilisation-side
factors**, **embryo culture conditions**, **multiple pregnancy**, or
**long-term effects on the ovarian tissue itself**.

---

## ⚠️ Disclaimer

This is a **qualitative/semi-quantitative QSP model for educational and
research purposes.** It was built from published literature and clinical
trial data but has not been independently validated or certified, and
**must not be used directly for clinical decision-making, prescribing, or
regulatory submission.** Ovarian stimulation in particular is an area
directly tied to individual patient safety, and no number in this model
can serve as the basis for an actual dosing decision.
