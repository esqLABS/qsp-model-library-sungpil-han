# Arrhythmogenic Right Ventricular Cardiomyopathy (ARVC) — QSP Model
### Arrhythmogenic Right Ventricular Cardiomyopathy · Arrhythmogenic Cardiomyopathy

**One-line summary — this model claims that the clock in ARVC is cumulative
mechanical work, not age.** A desmosomal variant does not destroy
myocardium. It only **lowers the load threshold** at which the
intercalated disc fails. So genotype does not set the amount of damage —
it sets only *the price paid per heartbeat*. Age is not a variable
anywhere in this model; it appears only as the integral of load.

| File | Contents |
|------|------|
| [`arvc_qsp_model.dot`](../../../arrhythmogenic-right-ventricular-cardiomyopathy/arvc_qsp_model.dot) · [SVG](../../../arrhythmogenic-right-ventricular-cardiomyopathy/arvc_qsp_model.svg) · [PNG](../../../arrhythmogenic-right-ventricular-cardiomyopathy/arvc_qsp_model.png) | mechanistic map — 199 nodes, 21 clusters, 321 edges |
| [`arvc_mrgsolve_model.R`](../../../arrhythmogenic-right-ventricular-cardiomyopathy/arvc_mrgsolve_model.R) | 54-compartment ODE model (mrgsolve), a 10-scenario runner |
| [`arvc_reference_impl.py`](../../../arrhythmogenic-right-ventricular-cardiomyopathy/arvc_reference_impl.py) | a dependency-free Python twin of the same equations (47 states) + a 26/26 self-check |
| [`arvc_shiny_app.R`](../../../arrhythmogenic-right-ventricular-cardiomyopathy/arvc_shiny_app.R) | 13-tab interactive dashboard |
| [`arvc_references.md`](../../../arrhythmogenic-right-ventricular-cardiomyopathy/arvc_references.md) | 84 entries, every PMID confirmed via the PubMed API |

```bash
python3 arvc_reference_impl.py            # all scenarios
python3 arvc_reference_impl.py --check    # PASS/FAIL self-check
python3 arvc_reference_impl.py --falsify  # the falsification switch (PHI_EX = 0)
```

---

## 1. Why a linear cascade will not do (Five awkward observations)

ARVC is usually drawn as a linear cascade: **desmosomal variant → junction
failure → cardiomyocyte death → fibrofatty replacement → arrhythmia.**
That picture is not wrong, but it cannot hold simultaneously with the
following five observations.

| | Observation | Problem for a linear cascade |
|---|------|------------------------|
| **A** | *PKP2* truncating variants are far more common in the population than ARVC prevalence. In unselected biobanks, most carriers of an ARVC gene loss-of-function variant have no phenotype <br>[PMID 31638835](https://pubmed.ncbi.nlm.nih.gov/31638835/), [25820315](https://pubmed.ncbi.nlm.nih.gov/25820315/) | A cascade that starts from genotype predicts **complete penetrance** |
| **B** | Exercise dose is graded and large. Desmosomal variant carriers in the top exercise tertile reach the phenotype at a much younger age, with a VT/death hazard ratio of about **3.16** <br>[PMID 23871885](https://pubmed.ncbi.nlm.nih.gov/23871885/), [25896080](https://pubmed.ncbi.nlm.nih.gov/25896080/) | There is no place for exercise in the cascade |
| **C** | The same phenotype occurs in athletes with **no desmosomal variant at all**, and they exercised **more**, not less, than genotype-positive patients <br>[PMID 25516436](https://pubmed.ncbi.nlm.nih.gov/25516436/) | There is an outcome with no starting point for the cascade |
| **D** | *PKP2* is expressed equally in both ventricles, yet **the right ventricle goes first.** Left-dominant disease is almost entirely a *DSP*/*FLNC* phenomenon <br>[PMID 37844667](https://pubmed.ncbi.nlm.nih.gov/37844667/), [37048743](https://pubmed.ncbi.nlm.nih.gov/37048743/) | A uniform defect cannot explain why only one ventricle is chosen |
| **E** | The treatment hierarchy runs backwards. Sotalol was **no different from nothing** in a North American registry; flecainide (a class IC agent in structural cardiomyopathy) reduces VT when added on; isolated endocardial ablation recurs while combined endo-epicardial ablation does not; and the **ICD, which touches no part of the cascade, is the only therapy with a clear mortality benefit** <br>[PMID 19660690](https://pubmed.ncbi.nlm.nih.gov/19660690/), [27939893](https://pubmed.ncbi.nlm.nih.gov/27939893/), [26546346](https://pubmed.ncbi.nlm.nih.gov/26546346/), [22205683](https://pubmed.ncbi.nlm.nih.gov/22205683/) | This ordering cannot arise from a single cascade producing a single arrhythmia |

---

## 2. Three structural commitments

### Commitment 1 — the clock is cumulative mechanical work (`.dot` clusters 3, 4)

Cardiomyocyte loss is written as a **fatigue failure** process.

```
loss rate = (beat count) × (wall stress / reserve)^NMECH ,   NMECH = 4
```

`NMECH = 4` was **fixed in advance, from the Basquin exponent range for
load-bearing biological tissue, and not fitted.** Genotype enters through
exactly one term: reserve capacity.

Wall stress depends on wall thickness, and wall thickness depends on
remaining myocardial mass. The process is therefore **autocatalytic**,
which is why decades of latency followed by rapid overt progression can
arise without any phase-transition parameter.

**This commitment is directly tested by an experiment the model was not
fitted to.** In plakoglobin-deficient mice, **load-reducing therapy
(furosemide + nitrate)** — a combination with no desmosomal, ion-channel,
or antifibrotic action whatsoever — prevented right-ventricular dilatation
and the arrhythmic phenotype and preserved conduction.
[Fabritz 2011 PMID 21292134](https://pubmed.ncbi.nlm.nih.gov/21292134/),
[Kirchhof 2006 PMID 17030684](https://pubmed.ncbi.nlm.nih.gov/17030684/)

**The falsification switch is a single parameter.** Setting `PHI_EX = 0`
makes damage independent of load, and the model becomes a calendar-clock
model. At that moment, A, B, C, and D **all invert simultaneously.**

### Commitment 2 — RV selectivity is Laplace, not biology (cluster 3)

**No ventricle-specific biology at all** is written for *PKP2* disease.
Two measured asymmetries are enough.

- the right-ventricular free wall is about 4 mm, the left about 9 mm →
  Laplace σ = P·r/2h
- at peak exercise, **right-ventricular** end-systolic wall stress rises
  by about **+125%**, while the **left ventricle** rises by only about
  **+14%**, because pulmonary vascular resistance does not fall as far as
  systemic resistance does.
  [La Gerche 2011 PMID 21085033](https://pubmed.ncbi.nlm.nih.gov/21085033/)

Raised to the fourth power, this becomes roughly a 15-fold difference in
per-beat fatigue load during exercise. That is why the right ventricle
goes first, and why **the RV/LV gap widens as training load rises** —
both of these are predictions, not inputs. Left-dominant disease requires
a separate genotype-specific term, `KAPPA_LV`, switched on only for
*DSP*/*FLNC*. It is the model's only ventricle-specific term.

### Commitment 3 — there are two arrhythmia generators, and no drug crosses the generator it does not occupy

| | **Generator I** (cluster 9, blue) | **Generator II** (cluster 10, brown) |
|---|---|---|
| Character | early-stage · catecholamine-dependent · **structure-independent** | late-stage · re-entrant · **scar-dependent** |
| Mechanism | PKP2 loss strips Nav1.5/Cx43 from the intercalated disc and destabilises RyR2 → diastolic Ca leak → triggered activity | patchy fibrofatty replacement. **Maximal when viable and replaced tissue interlock at roughly 50/50** — even a homogeneous scar has no conduction channels |
| Evidence | [19661460](https://pubmed.ncbi.nlm.nih.gov/19661460/) · [24352520](https://pubmed.ncbi.nlm.nih.gov/24352520/) · [31438494](https://pubmed.ncbi.nlm.nih.gov/31438494/) | [34883271](https://pubmed.ncbi.nlm.nih.gov/34883271/) · [22205683](https://pubmed.ncbi.nlm.nih.gov/22205683/) |
| Therapies that occupy it | beta-blockers, flecainide (RyR2), exercise restriction | epicardial ablation, amiodarone |
| Clinical meaning | why young carriers with a normal echocardiogram die during exercise | why sustained monomorphic VT appears in advanced disease |

Two consequences follow from this structure, and both are the model's most
falsifiable predictions.

- **Flecainide carries two opposite signs.** Blocking RyR2 removes
  Generator I, while blocking Nav1.5 slows conduction and **enlarges**
  Generator II. So the sign of the same drug depends on where the patient
  sits on the trajectory. The lesson of CAST falls out of the geometry
  rather than being a bolted-on warning.
- **The ICD attaches to the outcome, not to either generator.** That is
  why the ICD is the only therapy with a clear mortality benefit, and, at
  the same time, why nothing upstream can substitute for it.

---

## 3. The arithmetic of the clock — before solving any ODE

Using only `NMECH = 4` and two coefficients taken from exercise CMR, the
weekly-average fatigue load rate (a fully rested healthy ventricle = 1.00):

| Exercise | MET-h/week | `LOAD_RV` | `LOAD_LV` | RV/LV |
|------|---------:|----------:|----------:|------:|
| Sedentary | 6 | 1.287 | 1.015 | 1.27 |
| Guideline-level | 15 | 1.718 | 1.037 | 1.66 |
| Competitive athlete | 60 | 3.871 | 1.148 | 3.37 |
| Elite | 100 | 5.785 | 1.246 | 4.64 |

> **Predicted** ratio of right-ventricular fatigue rates (competitive :
> sedentary) = **3.01**
> **Observed** VT/death hazard ratio (highest : lowest exercise tertile)
> = **3.16**
> [James 2013 PMID 23871885](https://pubmed.ncbi.nlm.nih.gov/23871885/)

This value was not fitted. `NMECH` was fixed in advance at 4, and the two
stress coefficients were taken directly from exercise CMR. The RV/LV row
is likewise a pure output.

---

## 4. What was fitted and what was predicted

Of roughly 130 total parameters, **only 4 were fitted.**

| Parameter | Fitted to |
|----------|-----------|
| `K_INJ` | median age of definite Task Force diagnosis in *PKP2* carriers at ordinary leisure-activity levels ([PMID 25820315](https://pubmed.ncbi.nlm.nih.gov/25820315/)) |
| `H0_VA` | an annual sustained ventricular arrhythmia risk of about 10% in definite ARVC |
| `LAM2` | the efficacy gap between substrate-directed and trigger-directed therapy ([PMID 19660690](https://pubmed.ncbi.nlm.nih.gov/19660690/)) |
| `K_DIL` | the RVEDVi trajectory toward overt disease |

**Predicted (not fitted) — where the model puts its neck out:**

- an exercise hazard ratio of about 3-fold, and its dose-response shape
- incomplete penetrance, and how it depends on exercise stratum
- gene-elusive ARVC, appearing only in athletes with extreme load
- male predominance arising from a single multiplier (`SEX_K_FEMALE`)
- RV-before-LV with no ventricle-specific biology, and a widening gap with
  training
- load-reducing therapy preventing the phenotype (Fabritz 2011)
- **sotalol's null effect** — the arm that lengthens wavelength and the
  arm that widens repolarisation dispersion in a heterogeneous substrate
  (thereby promoting re-entry initiation) cancel each other out. The
  dispersion term's magnitude was fitted not to an ARVC registry but to
  SWORD ([PMID 8691967](https://pubmed.ncbi.nlm.nih.gov/8691967/)) — a
  randomised trial in which pure IKr block **raised** mortality in scarred
  ventricles — which makes the registry's null result a prediction
- flecainide helping early and harming late: one drug, two signs
- the gap between isolated endocardial and combined endo-epicardial
  ablation
- the ICD changing mortality alone, without touching the substrate at all
- **timing dominating dose** in AAV9-*PKP2* gene therapy

### A stated miss

This model predicts **a larger reduction in arrhythmia from beta-blocker
monotherapy than the North American ARVC registry observed** (the registry
found no significant beta-blocker benefit on the VT endpoint,
[PMID 19660690](https://pubmed.ncbi.nlm.nih.gov/19660690/)). This model
wires almost all of Generator I through beta-adrenergic drive. If the
registry is right, part of the RyR2 leak must be
**adrenaline-independent.** This is the cleanest point at which to
falsify Commitment 3.

**A second miss — but the same miss.** Ermakov 2017 reported reduced VT
when flecainide was **added** to a beta-blocker
([PMID 27939893](https://pubmed.ncbi.nlm.nih.gov/27939893/)). This model
instead makes that combination slightly **harmful** (+3 to +7%) at every
age, because beta-blockade has already driven Generator I to zero, leaving
flecainide with nothing but the cost of conduction delay. That is,
**this is the first miss seen from the other side**, and both point to
the same thing — part of the RyR2 leak must be adrenaline-independent.
Adding that term would stop beta-blockade from fully eliminating
Generator I (fixing miss 1), and simultaneously give flecainide-on-top-of-
beta-blockade something to act on again (fixing miss 2). **That two
independent discrepancies converge on a single missing term, and that the
fix is testable, is the most useful thing this model has produced about
its own structure.**

---

## 5. Penetrance is a property of a cohort

A single trajectory cannot answer "is this variant penetrant?" — penetrance
is a statement about a population. So the fatigue scale `K_INJ` carries an
inter-individual multiplier, `FRAILTY` (log-normal, SD 0.55, **fixed
quantiles, no random draw**). `SIGMA_FRAILTY` and the female multiplier of
0.72 were set in advance; what comes out is the prediction.

- a male *PKP2* carrier at guideline-level activity: **56%** definite
  diagnosis by age 40 (33% for females) → consistent with the observed
  "35–50% by middle age" range
- a sedentary male carrier: 33% by age 40, 67% by age 60 — a substantial
  fraction never reach the phenotype in a lifetime
- a competitive athlete: **100%** penetrance by age 40
- male > female ([PMID 28329361](https://pubmed.ncbi.nlm.nih.gov/28329361/))
- variant-negative, sedentary group: **0%** by age 60. Variant-negative
  elite athletes: **100%** by age 60 → gene-elusive ARVC occurs only under
  extreme load
  ([PMID 25516436](https://pubmed.ncbi.nlm.nih.gov/25516436/))

---

## 6. 54 ODE compartments (47 in the Python twin)

| Group | Compartments |
|------|------|
| Drug PK (13) | nadolol gut/central · flecainide gut/central · sotalol gut/central · amiodarone gut/central/peripheral · cumulative amiodarone toxicity · AAV vector · transduction fraction · transgene PKP2 |
| Desmosome · intercalated disc (6) | `PKP2_ID` · `DSP_ID` · `PG_NUC` (nuclear plakoglobin) · `WNT_ACT` · `HIPPO` · `NAV_ID` |
| Generator I (4) | `CX43_ID` · `RYR_LEAK` · `CA_SR` · `CA_DIA` |
| Right-ventricular tissue (7) | `D_RV` (cumulative mechanical load) · `MYO_RV` · `NEC_RV` · `INF_RV` · `FAP_RV` · `FIB_RV` · `FAT_RV` |
| Left-ventricular tissue (7) | identical equations, differing only in `KAPPA_LV` |
| Ventricular mechanics (4) | `RVEDV` · `RV_CONT` · `LVEDV` · `LV_CONT` |
| Electrophysiology (5) | `CV_RV` · `CV_LV` · `SCARHETRV` · `SCARHETLV` · `PVC24` |
| Neurohormonal · events (8) | `SNS` · `BETA1_D` · `NTBNP` · `ABL_HOM` · `H_VT` · `H_DEATH` · `H_HF` · `ICD_SHK` |

Total: 54. The Python twin has 47 states because it replaces 8 fast oral
drug absorption compartments with steady-state concentrations (that
approximation is verified numerically) and splits ICD shocks into
appropriate/inappropriate. Both differences are documented in each file's
header.

### Two constraints that hold the structure together

- **`MYO_FLOOR` (regional preservation).** This is a lumped free-wall
  model, but real disease is regional. The subtricuspid basal wall, apex,
  and outflow tract carry the highest local stress and fail first, while
  the peri-septal and outflow myocardium remains below threshold. So the
  lumped average cannot reach zero, and it does not in pathology studies
  either. Damage acts only on `(MYO − MYO_FLOOR)`, which bounds the
  autocatalytic runaway without needing a phase-transition parameter.
  ([PMID 9362410](https://pubmed.ncbi.nlm.nih.gov/9362410/))
- **`FF_THRESH` (circuit-size threshold).** Re-entry needs a path longer
  than the wavelength. A few percent of diffuse interstitial fibrosis
  cannot supply that. So ordinary age-related fibrosis is not
  arrhythmogenic, and Generator II has a genuine starting point rather
  than rising from zero together with the first cardiomyocyte loss.

---

## 7. Diagnosis was implemented, not summarised

The 2010 Task Force Criteria
([PMID 20172912](https://pubmed.ncbi.nlm.nih.gov/20172912/),
[20172911](https://pubmed.ncbi.nlm.nih.gov/20172911/)) are implemented as
written, down to the category I–VI counts and the definite/borderline/
possible rules. Every ECG and imaging item is **an output of a state
variable, never an input.**

| Item | How it arises in the model |
|------|----------------------|
| number of inverted T waves in V1–V3 | function of substrate extent (myocardial loss) |
| epsilon wave | `CV_RV < 22 cm/s` |
| terminal activation duration (TAD) | `32 ms × CV0/CV_RV`, minor at ≥55 ms |
| SAECG late potentials | TAD ≥ 50 ms |
| % residual myocardium | `MYO_RV / (MYO+FIB+FAT)`, major at <60% |
| PVC/24h | Generator I readout, minor at >500 |
| category VI (genetic) | **automatic major criterion** for a variant carrier — why genotyped relatives are diagnosed earlier than the proband |

The risk calculator
([Cadrin-Tourigny 2019 PMID 30915475](https://pubmed.ncbi.nlm.nih.gov/30915475/))
was **not reimplemented inside the model.** Its 7 inputs (age, sex,
syncope, NSVT, PVC/24h, number of leads with T-wave inversion, RVEF) are
instead **exported as outputs**, so that the externally validated score
can be applied to the model's own trajectories. Transcribing the published
linear predictor into the model would make the comparison circular.

---

## 8. The Shiny app (13 tabs)

Organised **around the three commitments**, not around organ systems —
hiding the commitments behind neatly organised organ-specific tabs would
erase the model's whole point.

1. Patient — genotype, sex, body size, and **exercise history as the
   disease's clock**
2. The clock — wall stress, fatigue load rate, cumulative load
3. Drug exposure — each drug's PK and **occupancy at its own target**
4. Structure — myocardial loss, fibro-fatty replacement, wall thickness,
   ventricular volumes/ejection fraction
5. The two generators — Generator I and II side by side on the same time
   axis
6. ECG/imaging — all outputs
7. Diagnosis — real-time Task Force Criteria scoring + risk calculator
   inputs
8. Outcomes — arrhythmia risk, event-free survival, survival, ICD shocks,
   cumulative amiodarone toxicity
9. Treatment bench — comparing standard therapy arms. **Which arm moves
   the substrate and which arm only moves the arrhythmia** is Commitment 3
10. RV vs LV — ventricular selectivity, what happens when `KAPPA_LV` is
    switched on
11. Penetrance — the cohort view
12. Falsification — `PHI_EX = 0` alongside `PHI_EX = 1`
13. Model notes — fitted / predicted / stated misses

---

## 9. Validation

The Python twin reproduces every number above.

```bash
python3 arvc_reference_impl.py --check
```

The checks include: the exercise hazard-ratio prediction, the training
dependence of the RV/LV gap, calibration of diagnostic age, incomplete
penetrance, the gene-elusive requirement, ventricular selectivity, the
magnitude of arrhythmia risk (overt/latent), the effect of load reduction,
the Generator-I specificity of beta-blockade, sotalol's null effect, the
sign-flip of flecainide, the superiority of epicardial ablation, the
outcome-only benefit of the ICD, the timing-dependence of gene therapy,
the actual reversal produced by the falsification switch, numerical
equivalence of the PK approximation, conservation of compositional mass,
and integrator convergence.

The R/mrgsolve file and the Python twin differ in **exactly two places**,
both documented in their headers: (i) fast oral drugs are replaced with a
steady-state average concentration (the approximation itself is verified
numerically), and (ii) the AAV dose is injected directly rather than
passing through the vector compartment's 9-day delay.

---

## ⚠️ Disclaimer

A semi-quantitative QSP model for educational and research purposes. It
was built from published literature but has not been independently
verified, and **must not be used for clinical decision-making.** The
AAV-*PKP2*, GSK-3β inhibition, IL-1 blockade, and load-reduction-therapy
scenarios are either in clinical trials or extrapolated from animal data,
and were simulated solely to explore the model's structure.
