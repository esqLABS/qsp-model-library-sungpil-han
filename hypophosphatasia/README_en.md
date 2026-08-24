# Hypophosphatasia (HPP) — QSP Model

**Hypophosphatasia · ALPL (1p36.12) loss of function → TNSALP deficiency →
accumulation of extracellular pyrophosphate (PPi) → competitive inhibition of
hydroxyapatite crystal growth**

This directory holds one complete set of quantitative systems pharmacology (QSP)
models for hypophosphatasia. It comprises the mechanistic map, an mrgsolve ODE
model, a Shiny dashboard, the references, and a dependency-free Python reference
implementation that **actually computes every number** cited in this README.

---

## 0. The question this model is built around

The textbook account is "ALPL loss of function → reduced alkaline phosphatase
(ALP) → soft bone". That account is qualitative, and it conceals the four
quantitative facts that actually determine the behaviour of the disease and the
dosing strategy. The model was built to expose those four, and every figure below
is taken verbatim from a run of `hpp_reference_model.py`
(`hpp_model_report.txt`).

### T1 — This is a threshold phenomenon, not a proportional relationship

Extracellular PPi is cleared **by TNSALP itself**, so `PPi ∝ 1/E`, and the
mineralisation rate is `∝ 1/(1 + PPi/Ki)`. Composing the two relations gives a
**hyperbola** in the local enzyme activity `E`, and the exact solution is the
positive root of the following quadratic (which the model solves in closed form
at every time point):

```
α·P² + [Vm·E + α·Km − S]·P − S·Km = 0
  where  α = kout·c₂/(c₁+c₂),  β = kout·Jsys/(c₁+c₂),  S = Jppi + β,
         c₁ = fvol·kout,  c₂ = kcatpp·E_plasma + CLppi,  Vm = kcat·Km
```

In the `PPi ≪ Km` limit the structure simplifies into a visible form:

```
Mrel(E) = (1 + PPi₀/Ki) / (1 + (Jppi + kout·PPi_plasma) / (Ki·(kcat·E + kout)))
```

The computed result — **a single curve generates the entire clinical spectrum:**

| Phenotype | Residual local activity | Extracellular (perivesicular) PPi | Mrel (1 = normal) | Osteoid volume fraction | Plasma PLP |
|--------|---------------|--------------------------|-----------------|------------------------|----------|
| Perinatal lethal (perinatal) | 0.010 | 47.3 µM | 0.120 | 10.6 % | 207 nM |
| Infantile | 0.030 | 40.1 µM | 0.137 | 9.4 % | 190 nM |
| Childhood | 0.120 | 21.3 µM | 0.229 | 6.0 % | 138 nM |
| Adult / odonto | 0.350 | 8.5 µM | 0.465 | 3.4 % | 85 nM |
| Heterozygous carrier | 0.500 | 6.0 µM | 0.604 | 2.8 % | 69 nM |
| Normal | 1.000 | 3.0 µM | 1.000 | 2.0 % | 45 nM |

- Raising the enzyme **11.7-fold** (3 % → 35 %) raises the mineralisation rate
  by only **3.39-fold**. Raising it a further 1.4-fold from there (35 % → 50 %)
  buys a mere 1.30-fold.
- The local activity required for half-maximal mineralisation is **E50 = 0.766**
  (exact solution) / 0.691 (linear approximation). The steep part of the curve is
  confined to **below 0.2** of normal activity. That is why a heterozygous
  carrier at 50 % activity is nearly normal, and an infant at 3 % is not.
- The osteoid volume fraction (normal < 2 %, osteomalacia > 10 %) reproduces the
  histomorphometric spectrum **even though it was never fitted**. `Jppi` was
  fixed by the single condition "PPi = 3 µM at normal activity", and was not
  fitted to clinical outcomes.

### T2 — The pharmacodynamic marker is being measured in the wrong compartment

Plasma PLP and plasma PPi are hydrolysed by the **plasma** enzyme. What
determines mineralisation is the enzyme **bound to hydroxyapatite** (which is
exactly what the deca-aspartate tail of asfotase alfa is for). The two pools do
not fill at the same rate.

- In the infantile form the **measurable** rise in plasma PPi is 2.57-fold, while
  the **causal** rise in extracellular PPi is 13.38-fold — a 5.2-fold
  underestimate.
- A 26-week dose-ranging exercise (infantile form, 5 kg, three times weekly):

| Weekly dose (mg/kg/wk) | Serum ALP (U/L) | HA occupancy | Local activity ELOC | Plasma PLP (nM) | Extracellular PPi (µM) | Mrel | RSS |
|---|---|---|---|---|---|---|---|
| 0 | 53 | 0.000 | 0.03 | 190 | 40.1 | 0.137 | 5.88 |
| 0.75 | 673 | 0.035 | 0.31 | 27 | 9.0 | 0.446 | 2.50 |
| 1.5 | 1252 | 0.062 | 0.53 | 20 | 5.4 | 0.657 | 1.44 |
| 3.0 | 2401 | 0.109 | 0.90 | 16 | 3.2 | 0.974 | 0.59 |
| **6.0 (label dose)** | **4768** | **0.193** | **1.57** | **14** | **1.81** | **1.205** | **0.55** |
| 12.0 | 9827 | 0.325 | 2.63 | 13 | 1.08 | 1.379 | 0.54 |
| 30.0 | 26800 | 0.562 | 4.52 | 12 | 0.63 | 1.513 | 0.53 |

The weekly dose at which each marker reaches a given percentage of **its own**
maximal change:

| Marker | 50 % | 80 % | 90 % |
|------|------|------|------|
| Plasma PLP | 0.22 | 0.35 | **0.64** |
| RSS (26 weeks) | 0.50 | 1.40 | 2.19 |
| Extracellular Mrel | 2.30 | 6.91 | **11.85** |

At the label dose of 6 mg/kg/week, PLP is already at **99 %** of its maximal
effect and RSS at 99 %, yet mineralising capacity (Mrel) is at **78 %**. The
ratio of doses required for a 90 % effect is **18.4-fold** (Mrel / PLP). In other
words, **a normalised PLP is not evidence that the skeletal dose is adequate.**
It must not be used as grounds for dose reduction.

### T3 — There are two clocks (reversible osteomalacia vs irreversible structural damage)

Same drug, same dose, **differing only in a 351-day start delay**. Both read out
at the 3-year time point:

| | Started at day 14 | Started at day 365 |
|---|---|---|
| Irreversible growth plate damage (GPD) | 0.060 | **0.295** (4.9-fold) |
| Change in height Z-score | −0.31 | **−1.58** (1.27 SD lost) |
| RSS | 0.24 | 1.18 |
| Osteoid volume fraction | 0.0182 | **0.0182** (identical) |
| Craniosynostosis index | 0.794 | 0.843 |
| Survival probability | 0.857 | 0.483 |
| Cumulative PPi exposure (µM·d) | 2 628 | 16 036 |

**The reversible markers (RSS, osteoid) converge, but height and
craniosynostosis do not.** The delay is paid for permanently. In the time course
of the label-dose treatment arm, extracellular PPi responds within days
(40.1 → 1.46 µM), the osteoid within a month (9.4 % → 2.3 %) and RSS over six
months (5.88 → 0.79), but GPD never decreases.

### T4 — Here an antiresorptive is mechanistically a poison

Bisphosphonates are **non-hydrolysable PPi analogues** with a skeletal half-life
of years. They enter through the very same `Ki` term that the disease is already
saturating. Adult HPP, 3 years:

| | Untreated | + bisphosphonate | + teriparatide |
|---|---|---|---|
| Bone-bound inhibitor (µM-eq) | 0 | 13.84 | 0 |
| Extracellular PPi (µM) | 8.48 | 8.48 | **4.80** |
| Local enzyme activity ELOC | 0.350 | 0.350 | **0.630** |
| Mrel | 0.465 | **0.319** (−31 %) | **0.712** (+53 %) |
| Osteoid volume fraction | 0.034 | 0.049 | 0.031 |
| Fracture burden | 0.094 | **0.775** (8.2-fold) | 0.052 |
| Pain (0–10) | 2.77 | 4.47 | 1.64 |
| Dental attachment (relative) | 0.511 | 0.361 | 0.748 |

The fracture burden becomes 8.2-fold higher **despite** resorption having been
suppressed. Conversely, teriparatide increases osteoblast numbers, raising the
**endogenous** enzyme supply 1.80-fold (a self-amplifying loop) and lowering
extracellular PPi from 8.48 to 4.80 µM — a route to more enzyme without
administering any.

### T5 — Split dosing is almost irrelevant; the weekly total is what matters

The HA-bound pool has a half-life of 3.5 days, which buffers the dosing interval.
The same 6 mg/kg/week (at the 1-year time point):

| | 2 mg/kg three times weekly | 1 mg/kg six times weekly | Ratio |
|---|---|---|---|
| HA occupancy | 0.181 | 0.180 | 0.995 |
| Extracellular PPi (µM) | 1.913 | 1.924 | 1.006 |
| Mrel | 1.184 | 1.182 | **0.998** |
| RSS | 0.254 | 0.255 | 1.007 |

By contrast, escalating from 2 to 3 mg/kg three times weekly raises Mrel from
1.184 to 1.297 (**+9.6 %**). What needs changing is not the schedule but the
total.

---

## 1. Deliverables

| File | Contents |
|------|------|
| [`hpp_qsp_model.dot`](hpp_qsp_model.dot) · [`.svg`](hpp_qsp_model.svg) · [`.png`](hpp_qsp_model.png) | Mechanistic map — **167 nodes, 224 edges, 16 clusters** |
| [`hpp_mrgsolve_model.R`](hpp_mrgsolve_model.R) | mrgsolve ODE model — **32 differential states + 4 quasi-steady-state (QSS) species**, 17 scenarios, calibration notes |
| [`hpp_reference_model.py`](hpp_reference_model.py) | Dependency-free Python reference implementation (**the numerical authority**). Its parameter block corresponds 1:1 with the R file |
| [`hpp_model_report.txt`](hpp_model_report.txt) | The full output of that script — the source of every number in this README |
| [`hpp_shiny_app.R`](hpp_shiny_app.R) | Shiny dashboard — **12 tabs** (profile · PK/target delivery · PPi economy · biochemical markers · skeleton · respiration/survival · mineral/kidney · CNS · vitamin B6 · scenario comparison · threshold explorer · dose-response · model notes) |
| [`hpp_references.md`](hpp_references.md) | **57 references** (sections 1–10 have PMIDs verified via the PubMed E-utilities, section 11 is separated out explicitly as unverified search links, and section 12 lists the parameters with no evidential basis) |

### Running it

```r
# mrgsolve model
source("hpp_mrgsolve_model.R")
HPP_report()            # summary table of the 17 scenarios
HPP_threshold_curve()   # T1: the hyperbola and its closed form
HPP_dose_ranging()      # T2: which marker saturates first

# Shiny
shiny::runApp("hpp_shiny_app.R")
```

```bash
# Python reference implementation (no external packages required)
python3 hpp_reference_model.py            # full report
python3 hpp_reference_model.py --brief    # omits the numerical verification section
```

---

## 2. Model structure

### The 32 differential states

| Group | States |
|------|------|
| **PK (4)** | `ASC` subcutaneous depot · `AC` plasma enzyme · `EB` **hydroxyapatite-bound enzyme (the effect compartment)** · `ADA` anti-drug antibody |
| **Cell and inhibitor pools (4)** | `OB` osteoblasts · `OC` osteoclasts · `BPB` bone-bound bisphosphonate · `OPN` phosphorylated osteopontin |
| **Vitamin B6 (1)** | `PLBR` central nervous system pyridoxal cofactor pool |
| **Mineral homeostasis (5)** | `CAS` serum Ca · `PIS` serum P · `PTHS` PTH · `VITD` 1,25(OH)₂D · `NEPH` nephrocalcinosis |
| **Bone matrix (2)** | `OST` unmineralised osteoid volume fraction · `BMIN` mineralised matrix |
| **Growth plate and height (4)** | `GP` growth plate integrity · `GPD` **irreversible** growth plate damage · `RSS` rickets severity score · `HTZ` change in height Z-score |
| **Thorax, respiration, survival (3)** | `RIB` rib cage mineralisation · `RESP` respiratory function index · `SURV` survival probability |
| **Symptoms and other organs (7)** | `MUS` muscle function (6MWT-like) · `PAIN` pain · `DENT` dental attachment · `CRAN` craniosynostosis · `SEIZ` seizure burden · `ECT` ectopic calcification · `FX` fracture burden |
| **Bookkeeping (2)** | `CUMHAZ` cumulative hazard · `AUCPPI` cumulative PPi exposure |

### The four species solved at quasi-steady state (QSS)

Extracellular PPi (τ ≈ 9 minutes), plasma PPi, plasma PLP and plasma PEA reach
equilibrium thousands of times faster than the disease time scale. These are
solved **in closed form** at every time point (the coupled PPi pair collapses
into the single quadratic above). The system thereby becomes non-stiff, so
fixed-step RK4 (Python) and LSODA (mrgsolve) give the same answer.

### The genotype enters as a single number

`FRACENZ` = the residual **local** TNSALP activity (as a fraction of normal). The
clinical phenotype is an **output, not an input**. Every scenario starts from a
burn-in at the untreated steady state for the genotype in question; at no point
is a "disease baseline" specified by hand.

---

## 3. The 17 scenarios

Each row reads out at the last time point of its own time course (S1 · S2 ·
S5–S8 · S11–S13 · S16 = 365 days, the rest = 1095 days). `ΔHTZ` is the change in
height Z-score relative to the starting point.

| Scenario | Extracellular PPi | Mrel | Osteoid | RSS | ΔHTZ | RESP | Survival | Muscle function | Seizures | Pain |
|---|---|---|---|---|---|---|---|---|---|---|
| S1 perinatal, untreated | 47.3 | 0.120 | 0.106 | 7.31 | −0.90 | 0.20 | **0.123** | 0.16 | 1.00 | 10.0 |
| S2 infantile, untreated | 40.1 | 0.137 | 0.094 | 6.88 | −0.86 | 0.39 | **0.512** | 0.22 | 0.99 | 10.0 |
| S3 childhood, untreated | 21.3 | 0.229 | 0.060 | 5.51 | −2.10 | 0.55 | 0.984 | 0.44 | 0.72 | 7.8 |
| S4 adult, untreated | 8.5 | 0.465 | 0.034 | 2.21 | 0.00 | 0.79 | 1.000 | 0.75 | 0.14 | 2.8 |
| S5 infantile + AA 2 mg/kg TIW | 1.91 | 1.184 | 0.018 | 0.25 | −0.18 | 0.95 | **0.887** | 0.80 | 1.00 | 0.0 |
| S6 infantile + AA 1 mg/kg six times weekly | 1.92 | 1.182 | 0.018 | 0.26 | −0.18 | 0.95 | 0.886 | 0.80 | 1.00 | 0.0 |
| S7 infantile + AA 3 mg/kg TIW | 1.40 | 1.297 | 0.017 | 0.25 | −0.18 | 0.95 | 0.888 | 0.80 | 1.00 | 0.0 |
| S8 infantile + AA 0.5 mg/kg TIW | 5.64 | 0.631 | 0.027 | 1.66 | −0.33 | 0.85 | 0.867 | 0.72 | 1.00 | 1.8 |
| S9 infantile + AA started at day 14 (3 years) | 1.95 | 1.176 | 0.018 | 0.24 | −0.31 | 0.95 | 0.857 | 0.80 | 1.00 | 0.0 |
| S10 infantile + AA started at day 365 (3 years) | 1.95 | 1.177 | 0.018 | 1.18 | **−1.58** | 0.95 | **0.483** | 0.80 | 1.00 | 0.0 |
| S11 infantile + AA, high-titre ADA | 2.96 | 1.006 | 0.020 | 0.26 | −0.18 | 0.95 | 0.887 | 0.80 | 1.00 | 0.0 |
| S12 perinatal + pyridoxine alone | 47.3 | 0.120 | 0.106 | 7.31 | −0.90 | 0.24 | 0.180 | 0.32 | **0.03** | 10.0 |
| S13 perinatal + AA + pyridoxine | 1.94 | 1.179 | 0.018 | 0.27 | −0.19 | 0.84 | **0.772** | 0.96 | **0.03** | 0.0 |
| S14 adult + bisphosphonate | 8.5 | **0.319** | 0.049 | 3.37 | 0.00 | 0.68 | 1.000 | 0.65 | 0.14 | 4.5 |
| S15 adult + teriparatide | 4.80 | **0.712** | 0.031 | 0.98 | 0.00 | 0.89 | 1.000 | 0.80 | 0.14 | 1.6 |
| S16 infantile + Ca/vitamin D restriction | 40.1 | 0.137 | **0.110** | 7.01 | −0.86 | 0.37 | 0.478 | 0.15 | 0.99 | 10.0 |
| S17 childhood + AA discontinued after 1 year | 21.3 | 0.229 | 0.060 | 5.16 | −1.36 | 0.56 | 0.992 | 0.45 | 0.72 | 7.8 |

A few things worth reading out of this:

- **Survival.** One-year survival is 12.3 % for untreated perinatal disease,
  51.2 % for untreated infantile disease, and 88.7 % at the label dose. Relative
  to the reported values (~95 % on treatment vs 27–42 % in historical controls),
  the treated arm has been left **deliberately conservative** — it was not
  adjusted to match.
- **The trade-off in S16.** Ca/vitamin D restriction lowers serum Ca from 2.60 to
  2.42 mmol/L and urinary Ca from 1.48 to 1.23 mmol/kg/day, and abolishes
  nephrocalcinosis. But it lowers the supersaturation along with them, so the
  osteoid volume fraction **worsens from 0.094 to 0.110.** Fixing the chemistry
  numbers and fixing the bone are different things.
- **Discontinuation in S17.** The plasma drug disappears within two weeks, but
  the skeletal effect takes two months to disappear (Mrel 1.198 at day 364 →
  0.390 at day 395 → 0.246 at day 425). What sets the rate of loss of effect is
  not the plasma PK but the HA-bound pool.
- **ADA in S11.** A high ADA titre cuts plasma exposure by 40 % (9.42 →
  5.67 mg/L) and lowers local activity from 1.48 to 0.96, yet **RSS is
  essentially unchanged at 0.254 vs 0.256.** At the label dose there is headroom
  in the radiographic endpoint, and it masks the immunogenicity.

---

## 4. The vitamin B6 axis — substrate excess and product deficiency at once

PLP does not cross the blood-brain barrier as such. TNSALP must first
dephosphorylate it to pyridoxal. Hence **plasma PLP is high while the cerebral
cofactor pool is low** — which is why the seizures respond to pyridoxine.

| | Plasma PLP (nM) | Cerebral cofactor (1 = normal) | Seizure burden | Muscle function | 1-year survival |
|---|---|---|---|---|---|
| S1 perinatal, untreated | 207 | **0.016** | 1.00 | 0.16 | 12.3 % |
| S12 perinatal + pyridoxine | 207 | **0.816** | **0.03** | 0.32 | 18.0 % |
| S13 perinatal + AA + pyridoxine | 13.5 | 0.805 | 0.03 | 0.96 | **77.2 %** |
| S2 infantile, untreated | 190 | 0.048 | 0.99 | 0.22 | 51.2 % |
| S5 infantile + AA | 13.5 | **0.014** | 1.00 | 0.80 | 88.7 % |

Pyridoxine restores the central cofactor without touching the skeleton — when the
thorax is the problem it **changes the seizures but not survival** (S12: 12.3 % →
18.0 %). Asfotase alfa, conversely, changes survival by fixing the skeleton. The
two are not substitutes for each other.

> ⚠️ **A hypothesis generated by the model, not something clinically
> established.** Asfotase alfa is bone-targeted and does not enter the central
> nervous system, so lowering plasma PLP amounts to lowering the substrate
> concentration gradient the brain was depending on (the cerebral cofactor of
> 0.014 in S5 is **lower** than the 0.048 of untreated S2). This is a mechanistic
> argument for continuing to monitor B6 status during ERT, and nothing more. This
> model contains no evidence that the effect has been observed in actual clinical
> practice.

---

## 5. Mineral homeostasis — paradoxical hypercalcaemia

Calcium that is not deposited in bone stays in the extracellular fluid. PTH is
suppressed, urinary calcium rises, and nephrocalcinosis accumulates.

| | Ca (mmol/L) | P (mmol/L) | PTH (pmol/L) | 1,25D (pmol/L) | Urinary Ca (mmol/kg/d) | Nephrocalcinosis |
|---|---|---|---|---|---|---|
| Normal reference | 2.40 | 1.80 | 3.50 | 100 | 1.200 | 0 |
| S2 infantile, untreated | **2.60** | **2.23** | **1.61** | 73 | **1.481** | **0.147** |
| S5 infantile + AA | 2.37 | 1.74 | 3.76 | 104 | 1.164 | 0.006 |
| S16 Ca/vitamin D restriction | 2.42 | 2.15 | 3.29 | 97 | 1.229 | 0.000 |

**Reversing the mineralisation defect is the way to correct the hypercalcaemia —
because bone is the calcium sink.** Dietary restriction works for the same
reason, but it cannot fix the bone (see the S16 trade-off in §3). Meanwhile,
extracellular PPi falls below normal during ERT (1.91 vs a normal 3.0 µM), so the
ectopic calcification term becomes active (index 0.0011 at 1 year) — which
connects to the point that ENPP1 deficiency (GACI) is the mirror-image disease on
the PPi **deficiency** side.

---

## 6. Numerical checks — discrepancies are reported, not absorbed

| Check | Result |
|------|------|
| Step halving (RSS at 182 d, h = 0.25 / 0.125 / 0.0625 d) | 0.774618 / 0.774612 / 0.774611 → \|Δ\| = 8.7 × 10⁻⁷. Discretisation error is not the issue here. |
| Exact quadratic QSS solution vs the 32-state integration (extracellular PPi) | relative deviation **0** (identical by construction) |
| The same comparison, Mrel | **22.6 %** — the whole of this gap is the osteopontin co-inhibition state (`OPN`), which the closed form, containing only PPi, does not have |
| **Linearised hyperbola** vs the integration | worst case **44.1 %** (E = 0.10), because in the severe phenotypes PPi reaches 43–95 % of Km. The hyperbola is therefore used **only to display the structure**, and every reported figure comes from the exact solution |
| PK mass balance (infantile, 2 mg/kg TIW, 182 days) | administered 720.0 mg · absorbable 329.8 mg · residual in depot 0.647 mg · plasma 2.694 mg · HA-bound 11.089 mg (capacity 57.6 mg, occupancy 0.193) |

---

## 7. What was calibrated to what

| Module | Basis | Status |
|------|------|------|
| asfotase alfa PK: F = 0.458, tmax 24–48 h, terminal t½ = 2.28 d | product label | label values (see references §11 — the label is not a PubMed record) |
| normal plasma PPi ≈ 3 µM, and a low-micromolar Ki of PPi for hydroxyapatite growth | PMID [13893487](https://pubmed.ncbi.nlm.nih.gov/13893487/), [6326671](https://pubmed.ncbi.nlm.nih.gov/6326671/) | literature-based |
| TNSALP–ENPP1–ANK antagonism sets the extracellular PPi balance | PMID [12082181](https://pubmed.ncbi.nlm.nih.gov/12082181/), [15039209](https://pubmed.ncbi.nlm.nih.gov/15039209/) | literature-based (structural) |
| `Jppi` | fixed by the condition "PPi = 3 µM at normal activity" | **not fitted to outcomes** |
| osteoid volume fraction 2 % → 9–11 % | — | **emergent result** (not fitted) |
| treating residual activity as a continuous variable | PMID [32160374](https://pubmed.ncbi.nlm.nih.gov/32160374/) | literature-based |
| the requirement for dephosphorylation of PLP at the BBB → pyridoxine-responsive seizures | PMID [7550313](https://pubmed.ncbi.nlm.nih.gov/7550313/), [17395561](https://pubmed.ncbi.nlm.nih.gov/17395561/) | literature-based |
| 1-year survival (12 % / 51 % / 89 %) | PMID [26529632](https://pubmed.ncbi.nlm.nih.gov/26529632/) | calibrated, with the treated arm kept conservative |
| 6MWT and pain scales | PMID [30576866](https://pubmed.ncbi.nlm.nih.gov/30576866/), [40138164](https://pubmed.ncbi.nlm.nih.gov/40138164/) | rough scaling |
| `KON` · `KOFF` · `KDEGB` · `BMAXKG` · `KBONE` | — | **assumed values**. HA occupancy must be read as a relative scale only, not as an absolute value |
| ADA module | — | **assumed values**. S11 is a sensitivity analysis, not a prediction |
| craniosynostosis (`KCS` · `TCS`) | — | **a phenomenological rate law**. It reproduces "common, and not reliably prevented by ERT" by construction; it does not explain it |

---

## 8. Limitations — do not cite these results without reading this list

1. **Body weight is fixed per scenario.** Since mg/kg dosing is handled without
   reflecting the change in body weight as the child grows, actual exposure in
   the multi-year paediatric scenarios declines more gently than in the model.
2. **The bone-binding PK is not identified** (see the table above). The absolute
   value of the occupancy, and the saturation statement that "three times the
   label dose raises occupancy by only 2.21-fold", depend directly on the choice
   of `BMAXKG`.
3. **The calcium homeostasis module operates near its excretory capacity**
   (normal urinary excretion is ~93 % of intake). That is reasonable as neonatal
   physiology, but serum Ca is **the most parameter-sensitive output** in this
   model. Read it as a direction, not as an absolute value.
4. **`HTZ` is a change relative to the starting point.** Because the burn-in
   resets height Z to 0, the baseline short stature an actual childhood-form
   patient already carries is not represented. The −2.10 over three years in S3
   must be read as "additional loss".
5. **Craniosynostosis, myopathy and pain** are rate laws and weighted sums, not
   mechanistic models. The mechanism of the myopathy in HPP in particular is not
   established, and the model assigns arbitrary weights to B6, bone pain and
   osteoid.
6. **The seizure burden is close to maximal even in the infantile form**
   (`FRACENZ` = 0.03). Clinically, seizures are largely confined to the perinatal
   and severe infantile forms, so this part discriminates severity more coarsely
   than reality does.
7. **The prediction that ERT may worsen the central B6 supply is a product of the
   model** and has no clinical evidence behind it (see the warning in §4).
8. **Immunogenicity, ectopic calcification and injection-site reactions** are
   represented at the level of expressing that a safety signal exists; they do
   not predict incidence rates.
9. This is a single-patient deterministic model. **There is no inter-individual
   variability (IIV), no residual error, and no uncertainty propagation.** If
   confidence intervals are needed, `HPP_run()` must be wrapped over a sample of
   parameters.

---

## 9. Falsifiable predictions

Set down in a form in which the model can be shown to be wrong:

1. **Even among patients whose PLP has been normalised on the standard dose**,
   the radiographic and growth responses will continue to diverge with the total
   weekly dose (§T2). A cohort that uses PLP normalisation as a dose-reduction
   criterion will lose out on the skeletal endpoints.
2. **Three times weekly versus six times weekly should make no difference at the
   same total dose** (predicted difference 0.2 %). If a significant difference is
   observed, the assumed half-life of the HA-bound pool (3.5 days) is wrong.
3. **A one-year delay in starting treatment permanently costs about 1.3 SD of
   final height while leaving no trace in RSS or osteoid** (§T3). The predicted
   pattern in a delayed cohort is that RSS catches up but height does not.
4. **Bisphosphonate exposure in adult HPP raises the fracture rate** (model:
   8.2-fold). This already agrees with clinical reports (PMID
   [41503045](https://pubmed.ncbi.nlm.nih.gov/41503045/)); the model quantifies
   the magnitude.
5. **The benefit of teriparatide should be proportional to the residual
   activity** — increasing osteoblast numbers has nothing to multiply if
   `FRACENZ` is close to 0. The prediction, in other words, is that it is
   effective in the adult form (0.35) and essentially ineffective in the
   infantile form (0.03).

---

## ⚠️ Disclaimer

This is a qualitative to semi-quantitative QSP model for the purposes of
education, research and hypothesis generation. It has not been independently
validated and **must not be used for clinical decision-making, prescribing, or
regulatory submission.** The parameters are approximations for illustrative
purposes, and transferring the numbers without also citing the assumptions and
limitations set out in §7 and §8 is a misuse of them.
