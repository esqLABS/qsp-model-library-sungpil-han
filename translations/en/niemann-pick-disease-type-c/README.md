# Niemann-Pick Disease Type C (NPC) — Quantitative Systems Pharmacology Model

<a href="../../../niemann-pick-disease-type-c/npc_qsp_model.svg"><img src="../../../niemann-pick-disease-type-c/npc_qsp_model.png" width="820" alt="NPC QSP mechanistic map"></a>

> **198 nodes · 287 edges · 22-cluster mechanistic map · 41-ODE mrgsolve model ·
> 18 scenarios · 9-tab Shiny app · 144 individually PubMed-verified references ·
> an independently written Python reference implementation (all 142 shared
> parameters agreeing)**

---

## 1. The question this model sets out to answer

In NPC there are four drugs acting at **four points along one causal chain**, and two of them
were **approved by the FDA in the same month (September 2024)**. Yet the two were tested with
entirely different designs and different endpoints.

| Drug | Mechanistic point | Trial design | Primary endpoint | Result |
|----|-----------|-----------|----------|------|
| Miglustat | UGCG inhibition (substrate influx↓) | randomised, 12 months | horizontal saccade velocity | improved |
| **Arimoclomol** | HSF1 co-induction → folding yield of mutant NPC1↑ | randomised parallel, **12 months** | **change in 5-domain NPCCSS** | −1.40 (p=0.046) |
| **Levacetylleucine** | cerebellar/vestibular metabolism and function (symptomatic) | randomised **crossover**, **12 weeks** | **change in SARA** | −1.28 (p<0.001) |
| Adrabetadex | direct extraction of cholesterol from the membrane (bypassing NPC1) | open-label → sham-controlled | NPC-NSS | positive open-label, **negative randomised** |

This model was built in order to **separate three things** that the literature habitually
lumps together.

### ① The compartment claim — the markers are not looking at where the disease is

The plasma markers (cholestane-3β,5α,6β-triol, lysoSM-509/PPCS, bile acid B) are oxidation
products of sterol stored in **liver and spleen macrophages**. That is, they read the
**visceral compartment**. The disease that kills people is in the **cerebellum**. The two
compartments differ in drug access (oral vs intrathecal), in turnover (days vs years), and in
reversibility.

### ② The reserve claim — function splits into two pools

```
FUNC = 1 − D_rev − D_irr
```
`D_rev` returns with a time constant of about two weeks; `D_irr` is **monotone
non-decreasing**. Two treatments with **exactly identical** effects on storage produce
completely different trajectories depending on which side of this split they act. This is the
only structure that separates a **symptomatic drug** from a **disease-modifying** one.

### ③ The design claim — each drug was tested with a design that can only see its own assumed mechanism

This is a necessary consequence of ②. A 12-week window can only see `D_rev`, and a 12-month
"change from baseline" **cannot distinguish a constant offset from a change in slope.**

---

## 2. The single structural choice that does the most work

`DAM` is **the integral of cerebellar stress, and it has no repair term.** Purkinje cell death
is gated on `DAM` exceeding `D_reserve`. From this one line **two published facts that were
never calibrated to fall out at once.**

* Because the delay to onset is **inversely proportional** to stress, the residual NPC1
  activity spreads the age of onset out over several decades.
* Because the death rate past the gate **saturates**, the slope after onset becomes almost
  independent of the age of onset.

That is exactly the title of Yanjanin 2010
([PMID 19415691](https://pubmed.ncbi.nlm.nih.gov/19415691/)) —
*"linear clinical progression, **independent of age of onset**"*.

> **The discarded alternative supports this choice.** The first version used
> `dDAM = k·stress − krep·DAM`, with a repair term. With a repair term the time at which the
> gate is passed depends on stress **only logarithmically**, so the range of ages of onset
> collapses to 2–3 years and **the adult-onset form cannot be produced at all.** In that
> version the mild/mild genotype in fact never reached the gate and was asymptomatic for life.

The age-of-onset forms the model reproduces (developmental vulnerability `v_dev` is a
patient-level covariate):

| Age form | Genotype | v_dev | Residual f_NPC1 | Cerebellar CHOL fold | Model onset | Model median survival |
|--------|--------|-------|-------------|----------------|-----------|-----------------|
| Perinatal | null/null | 2.5 | 0.000 | 27.9× | 4.4 y | **4.5 y** (liver failure) |
| Early infantile | I1061T/null | 2.5 | 0.040 | 30.2× | 5.0 y | 22.8 y |
| Late infantile | I1061T/null | 1.2 | 0.040 | 30.2× | 6.5 y | 24.4 y |
| Juvenile | I1061T/I1061T | 0.6 | 0.080 | 26.4× | 8.0 y | 26.0 y |
| Adolescent/adult | mild/mild | 0.2 | 0.416 | 4.6× | 29.8 y | 52.4 y |
| NPC2 disease | NPC2 | 1.0 | 1.000 | 31.2× | 8.4 y | 26.4 y |

The perinatal form **dies of liver failure (4.5 y) just after its neurological onset (4.4 y)**
— the model produces the "phenotype that dies of liver disease first" without any separate
rule. The `NPC2` row shows the same disease as the juvenile form arising even though the NPC1
protein is **completely normal**. That is because NPC1 and NPC2 act **in series** on one flux,
so their functional fractions multiply.

---

## 3. Deliverables

| File | Contents |
|------|------|
| [`npc_qsp_model.dot`](npc_qsp_model.dot) · [`.svg`](../../../niemann-pick-disease-type-c/npc_qsp_model.svg) · [`.png`](../../../niemann-pick-disease-type-c/npc_qsp_model.png) | 198 nodes · 287 edges · 21 mechanism clusters + legend |
| [`npc_mrgsolve_model.R`](../../../niemann-pick-disease-type-c/npc_mrgsolve_model.R) | 41-ODE mrgsolve model, 148 annotated parameters, 6 genotypes · 6 age forms · 18 scenarios · design-experiment functions |
| [`npc_reference_model.py`](../../../niemann-pick-disease-type-c/npc_reference_model.py) | Independently written Python twin. Prints every number in this README (`python3 npc_reference_model.py`) |
| [`npc_shiny_app.R`](../../../niemann-pick-disease-type-c/npc_shiny_app.R) | 9-tab dashboard (patient · PK · lysosome · reserve · endpoints · scenarios · biomarkers · design experiments · calibration comparison) |
| [`npc_references.md`](npc_references.md) | 13 sections, 144 references, every PMID individually looked up and confirmed |

The **142 parameters shared by the two implementations are synchronised by a script**
(`npc_sync_params.py`), and during development that script caught three values that had
drifted while being copied across by hand (`Emax_fold` 3.4453↔0.30, `n5_k` 58.4↔30.0,
`h_swal` 6.79e-4↔1.23e-4).

---

## 4. Summary of the model structure

### 41 state variables

| Group | State variables |
|------|----------|
| Drug PK (13) | `MIG_GUT/CEN/PER/BR` · `ARI_GUT/CEN/CSF` · `NAL_GUT/CEN/BR` · `CD_CSF/BR/SYS` |
| NPC1 protein (2) | `NPC1_ER` · `NPC1_L` |
| Visceral lipids (2) | `CHOL_V` · `GSL_V` |
| CNS lipids (4) | `CHOL_C` · `GSL_C` · `SPH` · `CA_LY` |
| Lysosomal / metabolic function (4) | `HYD` · `AUTOPH` · `MITO` · `ROS` |
| Cerebellum (7) | `PC` · `PC_S` · `PC_LOST` · `INFL` · `SYN` · `CBL` · **`DAM`** |
| Functional decomposition (2) | **`D_REV`** · **`D_IRR`** |
| Biomarkers (5) | `TRIOL` · `PPCS` · `TCG` · `NFL` · `CALB` |
| Cochlea · survival (2) | `OHC` · `CUMHAZ` |

### Separation of time constants — what the latency actually is

The CNS lipid compartment is **25 times slower** than the visceral compartment, at
`tau_cns = 25` (`Jin_c`, `Vmax_c` and `kbas_chol_c` are all divided by 25). The steady state
is **unchanged**; only the rate of approach to it slows. It is a direct reflection of the fact
that brain cholesterol turnover takes months to years while the hepatic pool takes days, and
this single factor generates the several-year latency between birth and neurological onset.

### Why `D_rev` was not written as a "fraction of stressed cells"

At first `D_rev ∝ PC_S` (the fraction of cells in the stressed state). With that, `PC_S`
saturates above 0.9 within a few weeks of birth and **NPCCSS5 reaches 21 points in the first
year of life**. No NPC patient progresses like that. So `D_rev` was rewritten as a
**saturating function of the current biochemical stress** — intensity, not cell count.

---

## 5. Calibration vs validation (7 against 8)

Of the 15 targets in `npc_references.md`, **only 7 were used for calibration**; the rest are
validation only.

| # | Target | Source | Published | Model | Role |
|---|------|------|--------|------|------|
| T1 | plasma C-triol, patients | PMID 33228797 | 88.31 ng/mL | **88.28** | calibration |
| T2 | plasma C-triol, controls | PMID 33228797 | 5.97 ng/mL | **5.95** | calibration |
| T3 | 5-domain NPCCSS progression | PMID 33228797 | 1.50 points/year | **1.500** | calibration |
| T4 | 17-domain NPCCSS progression | PMID 33228797 | 2.7–2.9 points/year | **2.925** | 🔍 validation |
| T8 | levacetylleucine 12-week SARA difference | PMID 38294974 | −1.28 points | **−1.27** | calibration |
| T9 | C-triol vs NPCCSS5 ρ | PMID 33228797 | 0.265 | **0.526** | 🔍 validation (partial) |
| T11 | baseline SARA (IB1001-301) | PMID 38294974 | 15.91 | **15.91** | calibration |
| T14 | independence of progression from age of onset | PMID 19415691 | — | **reproduced** | 🔍 structural validation |
| T5 | arimoclomol 12-month 5-domain difference | PMID 34418116 | −1.40 points | **−0.17** | ❌ **falls short** |
| T6 | arimoclomol 12-month 4-domain difference | PMID 40520915 | −1.70 points | **−0.15** | ❌ **falls short** |
| T12/13 | levacetylleucine 18/12-month change from baseline | PMID 40513057 | −1.64 / −1.88 | **+0.88 / +0.17** | ❌ **disagreement** |
| T15 | arimoclomol miglustat subgroup | PMID 34418116 | −2.06 points | **no interaction** | ❌ **not reproduced** |

**Reproducing T14 matters particularly.** It comes out of the damage-integral structure by
itself:

| Age form | Onset | Slope over the 5 years after onset |
|--------|------|--------------------|
| Late infantile | 6.5 y | 1.52 points/year |
| Juvenile | 8.0 y | 1.42 points/year |
| Adolescent/adult | 29.8 y | 0.49 points/year |

The onset times are 23 years apart and the slopes are within a factor of three.

### Conservation and monotonicity checks (the model checks itself)

```
Purkinje pool conservation   max|PC + PC_S + PC_LOST − 1| = 9.8e-15
D_IRR monotone               min step = 0.0
PC_LOST monotone             min step = 5.0e-19
DAM monotone                 min step = 3.1
healthy subject, 45 years    NPCCSS5 = 0.0000, PC_LOST = 0, DAM/reserve = 0.000
```

---

## 6. What the model actually said — four results

### Result 1. The reason C-triol cannot grade severity is saturation, not noise

Calibrating **on just the two points** of a patient value of 88.31 and a control value of 5.97
forces a saturation constant of `Ktri = 8.35`. But the patient's visceral pool is
`CHOL_V = 36.8` — that is, **triol production is already at 81.5% of its own ceiling**.

Build a virtual cohort from **the same population as the published study (2–18 years,
paediatric genotypes)** and triol is confined to a 2.3 ng/mL span from 88.3 to 90.6 ng/mL,
while NPCCSS5 moves from 0.6 all the way to 21.9. A deterministic model with **not one noise
term** in it yields ρ = 0.526. That falls short of the published 0.265 (the residual is
presumably assay error and biological variability), but it is plainly different from the
ρ ≈ 1 that an unsaturated marker would give.

Widen the cohort to include the adult forms and the triol range opens out to 27.3–90.6 and ρ
rises to 0.605. **The marker carries information where it is not saturated, and the patients
the trials enrol are not there.**

### Result 2. Each drug was tested with a design that can only see its own mechanism

Match a **purely symptomatic** drug and a **purely disease-modifying** drug so that they have
identical 12-month NPCCSS5 efficacy (−0.177 points), then run each through the two published
designs.

| Design | Purely symptomatic | Purely disease-modifying | Sensitivity ratio |
|------|-----------|---------------|-----------|
| **A: 12-week crossover, SARA** (IB1001-301) | **−1.271** | −0.221 | **5.8-fold** |
| **B: 12-month parallel, NPCCSS5** (arimoclomol) | −0.177 | −0.177 | **1.0-fold (indistinguishable)** |

Design A is 5.8 times more sensitive to the symptomatic component, and design B **cannot in
principle distinguish the two mechanisms** (a constant offset and a change of slope give the
same "change from baseline"). **Neither published design identifies this split.**

The asymmetry comes from the scales themselves. SARA measures what one can do under cerebellar
control at this moment, so it weights `D_rev` heavily (`sara_a = 1.8`), whereas the NPCCSS
domains (gait · speech · swallowing · fine motor · cognition) are milestone-like and so weight
`D_irr` heavily (`n5_a = 0.25`).

**The only design that separates the split is a randomised withdrawal.** In the model,
stopping after 12 months of treatment lets SARA **rebound by +1.63 points** within 90 days,
and the residual gap against the untreated arm after 12 months is only **−0.06 points** — that
is, virtually all of this drug's benefit is the reversible component.

### Result 3. Arimoclomol's effect size is set by reserve, not by dose

Sweep the folding-rescue strength (`Emax_fold`) from the physiological range up to the
theoretical ceiling (patient entering at age 13, 12-month NPCCSS5):

| `Emax_fold` | Change in f_NPC1 | Reduction in progression |
|-------------|-------------|-----------|
| 0.15 | 0.080 → 0.155 (1.9-fold) | 6.9% |
| **0.30** (adopted) | 0.080 → 0.227 (2.8-fold) | **11.2%** |
| 1.20 | 0.080 → 0.589 (7.3-fold) | 29.3% |
| 3.40 | 0.080 → 0.921 (11.5-fold) | 48.2% |
| 60.0 | 0.080 → **1.000** (complete normalisation) | **48.4%** |

**Restoring folding entirely to normal erases only 48% of the progression.** The published
value is 65% (−1.40 points, 95% CI −2.76 ~ −0.03).

The reason lies in reserve being an **integral**. `DAM` does not decrease, so once the gate
opens it does not close. **Repairing the primary defect after the reserve has been spent
cannot return the rate of progression to normal.** Now leave the drug alone and change only
**the age at entry**:

| Age at entry | DAM/reserve | Reduction in progression |
|-----------|------------|-----------|
| 5 y | 0.63 (before the gate) | **46.9%** |
| 8 y | 1.14 (just past the gate) | 14.1% |
| 13 y | 2.04 | 11.2% |
| 16 y | 2.59 | 12.2% |

**A falsifiable prediction:** arimoclomol's effect should be strongly modulated by the state
of reserve at entry — large before the gate and almost absent after it. The published trial
enrolled **ages 2–18**, so the point estimate is probably being driven by the youngest
patients, and a confidence interval running from −2.76 to −0.03 is exactly what that looks
like. The model value of −0.17 lies **inside** that interval, but is far smaller than the point
estimate.

### Result 4. Cyclodextrin cannot separate efficacy from ototoxicity

Cyclodextrin's efficacy **is the extraction of cholesterol from membranes**, and the outer hair
cell is a membrane that structurally requires cholesterol
([PMID 26903308](https://pubmed.ncbi.nlm.nih.gov/26903308/)).
Hearing falls even in normal animals
([PMID 20357695](https://pubmed.ncbi.nlm.nih.gov/20357695/)).
The 17-domain NPCCSS **contains a hearing domain.**

| 12-month measure | Untreated | IT-CD 900 mg q2wk |
|-------------|--------|-------------------|
| Cerebellar lysosomal cholesterol (× normal) | 27.30 | **1.79** (normalised) |
| CSF calbindin | 0.25 | 0.20 |
| Plasma C-triol | 88.28 | **88.24** (essentially unchanged — intrathecal does not reach the viscera) |
| Hearing threshold shift | 0 dB | **32 dB** |

| 5-domain NPCCSS gain (no hearing domain) | **−0.72 points** |
|---|---|
| 17-domain gain, scored **excluding hearing** | **−1.40 points** |
| 17-domain gain, scored **including hearing** | **+1.55 points** |

**The hearing domain does not merely erase the benefit, it reverses its sign.** This is why
Ory 2017 had to score *"NSS minus hearing"*
([PMID 28803710](https://pubmed.ncbi.nlm.nih.gov/28803710/)), and why the randomised
sham-controlled phase 2b/3 that put hearing into the composite and removed the open-label
expectation effect **found no difference on any endpoint.**

Meanwhile plasma triol barely moves, from 88.28 to 88.24. **Because intrathecal administration
does not reach the visceral compartment, the plasma marker reports nothing whatsoever about
what this drug did in the brain** — the compartment claim of Result 1 confirmed again, this
time from the treatment side.

---

## 7. Discrepancies reported without adjustment (misses)

Following this repository's principle, **what does not fit is written down as it is rather
than made to fit.**

### ① Arimoclomol's effect size falls short (T5, T6)

Result 3 above. −0.17 points at physiological strength, and even at the theoretical ceiling
−0.71 points. The published value is −1.40. The model value is inside the published confidence
interval but is one eighth of the point estimate. The model's reading is that "in a patient
whose reserve is spent, this mechanism cannot produce an effect of that size", and that is
itself a falsifiable prediction.

### ② The levacetylleucine long-term extension does not add up arithmetically (T12, T13)

The 12-week randomised difference (−1.28) is reproduced **exactly**. But the published
single-arm extension gives a change in SARA from baseline of −1.88 at 12 months and −1.64 at
18 months.

| Window | Model (vs baseline) | Model untreated | Model difference | Published |
|----|--------------------|-------------|-----------|--------|
| 12 weeks | −0.93 | +0.35 | **−1.28** | −1.28 (vs placebo) ✓ |
| 12 months | **+0.17** | +1.50 | −1.34 | **−1.88** (vs baseline) ✗ |
| 18 months | **+0.88** | +2.24 | −1.36 | **−1.64** (vs baseline) ✗ |

That is, with a symptomatic effect of that size the disease progresses underneath, so **by 12
months baseline should already have been exceeded**. The untreated SARA progression implied by
the published values is **0.09–0.22 points/year**, whereas the value the model *derives* from
the NPCCSS target is **1.50 points/year** — a 7–16-fold discrepancy.

There are three possible readings and the model does not choose between them.
(a) SARA progression in NPC is far slower than NPCCSS progression — that is, the two scales
are not measuring a disease progressing at the same rate. (b) Levacetylleucine has a
disease-modifying component larger than the model can produce (blocking CNS lipid influx 100%
was still not enough — the required `nal_Emax_dm` was 2.47, that is, a physically impossible
value). (c) Change from baseline in a single-arm extension with no concurrent control
overestimates the effect.
**There is no published untreated SARA slope for NPC**, so these three cannot be adjudicated.
This is the most specific data gap this model points to.

### ③ No miglustat × arimoclomol synergy appears (T15)

**A prediction made in the expectation that it would appear turned out to be wrong.** HSP70
raises the folding *yield* and miglustat lowers the *load* on the same flux, so they ought to
act multiplicatively. Yet in the model:

```
arimoclomol alone      0.165 points saved
miglustat alone        0.013 points saved
additive prediction    0.179 points
model (combination)    0.178 points  → no synergy (−0.001)
```

The model states the reason outright: miglustat's **CNS** UGCG inhibition is only 12.5%
(brain:plasma partition 0.45; visceral 23.6%), and that is worth 0.013 points/year. **There is
nothing to multiply with.** The published difference between the miglustat subgroup (−2.06) and
the overall result (−1.40) is not reproduced by the model, and that non-reproduction is written
down.

Incidentally this becomes a testable statement about miglustat: the model predicts that
miglustat's effect should be **about twice as large on the visceral and bulbar side as on the
cerebellar side**, which agrees in direction with Patterson 2007 reporting improved swallowing
and stabilised hearing while failing to halt the progression of ataxia itself.

### ④ The late infantile and juvenile forms cannot be separated by genotype alone

Because the storage phenotype saturates, I1061T/null and I1061T/I1061T differ two-fold in
residual activity yet have almost the same steady-state burden (30.2× vs 26.4×). The two age
forms are specified **separately**, through the developmental vulnerability `v_dev`, and are
not derived from genotype. It is a limitation consistent with the real observation that
siblings with identical genotypes can have onsets ten years apart, but **it is a limitation.**

### ⑤ Five real defects caught during development

Recorded here. These are the ones that could have survived as plausible curves.

1. **The final state was not recorded.** `simulate()` returned the last state on the recording
   grid, so a burn-in run with `record_every = 400 days` had a trial entry age of **12.05
   years** rather than 13.0. Every number in the treatment scenarios was shifted by up to a
   year.
2. **The folding yield exceeded 1.** Making `Emax_fold` large gave `theta > 1`, so f_NPC1 went
   above normal and CHOL_C went negative, whereupon a fractional power of a negative number
   produced a **complex number** and the computation died. Solved by imposing the physical
   upper bound that a yield is a fraction.
3. **`D_rev` was tied to a cell fraction that saturates.** NPCCSS5 of 21 points within weeks
   of birth. Rewritten as a saturating function of the current biochemical stress.
4. **A damage variable with a repair term.** The range of ages of onset collapsed and the adult
   form disappeared. Switching to a pure integral made the delay to onset inversely
   proportional to stress.
5. **Developmental vulnerability could not move the age of onset.** `vuln` multiplied only the
   death rate, so it could not change the time at which the gate was passed. Multiplying it
   into the damage integral as well opened the ages out properly, to early infantile 5.0 y /
   late infantile 6.5 y / juvenile 8.0 y.
6. **Two cyclodextrin constants were unphysiological.** `cd_kext = 0.013` emptied cerebellar
   cholesterol down to **0.06× normal** (arithmetically harmless, since stress bottoms out at
   zero, but biologically absurd), and `cd_koto = 1.65e-4` wiped out the outer hair cells
   within a few months. They were re-anchored to normalisation (1.79-fold) and to 32 dB at 12
   months respectively. Efficacy barely changed — because the stress term saturates within the
   normal range.

---

## 8. The value of starting early — a consequence of reserve being an integral

The same triple combination is started at four time points and followed to age 20.

| Age at start | DAM/reserve at start | Duration of dosing | NPCCSS5 @ age 20 | D_rev@20 | D_irr@20 | Points saved |
|-----------|--------------------|-----------|--------------|----------|----------|-----------|
| 2 y | 0.18 | 18 years | **11.46** | 0.049 | 0.370 | **8.10** |
| 5 y | 0.63 | 15 years | 15.62 | 0.048 | 0.509 | 3.94 |
| 8 y | 1.14 | 12 years | 17.65 | 0.048 | 0.576 | 1.91 |
| 12 y | 1.86 | 8 years | 18.55 | 0.050 | 0.606 | 1.01 |
| Untreated | — | 0 | 19.56 | 0.121 | 0.622 | 0 |

`D_rev` recovers to the same 0.048–0.050 whatever the starting point — **the reversible
component can be bought at any time.** The whole difference comes from `D_irr` (0.370 vs
0.606). The price of starting six years late (age 2 → age 8) is **6.2 points** at age 20. This
is the quantitative argument for newborn screening for NPC
([PMID 27147587](https://pubmed.ncbi.nlm.nih.gov/27147587/)), and it is a necessary consequence
of the model's structure — an integral cannot be undone.

---

## 9. Drug exposure verification

| Drug | Regimen | Model Cavg | Model Cmax | Literature |
|----|------|-----------|-----------|------|
| Miglustat | 200 mg tid | 11.54 μM | 12.61 μM | Css 11–14 μM |
| Arimoclomol | 124 mg tid | 358 ng/mL | 435 ng/mL | linear PK, t½ ~4 h |
| Levacetylleucine | 4 g/day in divided doses | 1.20 mg/L | 2.83 mg/L | label Cmax 8.3 μg/mL |
| Adrabetadex | 900 mg IT q2wk | 55.5 mg/L (CSF) | 5721 mg/L | CSF t½ ~4 h |

Miglustat's CNS UGCG inhibition of **12.5%** versus **23.6%** in the viscera is a
**prediction** arising from the brain:plasma partition of 0.45, and it explains the absence of
synergy in result ③.

---

## 10. How to run

```r
# mrgsolve
library(mrgsolve)
mod <- mread("npc_mrgsolve_model.R")
npc_calibration_table(mod)        # model vs published values
npc_scenario_summary(mod)         # summary of the 18 scenarios
npc_design_experiment(mod)        # the design experiment

# Shiny
shiny::runApp("npc_shiny_app.R")
```

```bash
# the Python reference implementation — prints every number in this README
python3 npc_reference_model.py
python3 npc_reference_model.py --quick        # coarse dt, about 4× faster
python3 npc_reference_model.py --emit-param   # dump the mrgsolve $PARAM values
```

---

## 11. Disclaimer

This is a semi-quantitative QSP model for education and research. It was built from the public
literature and clinical trial data but has not been independently validated or qualified, and
**must not be used for real clinical decisions, prescribing, or regulatory submission.** In
particular, the discrepancies of §7 above are the parts the model cannot yet explain, and its
numbers for those parts must not be relied on as evidence.

For the complete reference list and the basis of each parameter, see
[`npc_references.md`](npc_references.md).
