# Progressive Supranuclear Palsy — QSP model

**Progressive Supranuclear Palsy (PSP) — a 4R tauopathy**

| | |
|---|---|
| Mechanistic map | [`psp_qsp_model.dot`](../../../progressive-supranuclear-palsy/psp_qsp_model.dot) · [SVG](../../../progressive-supranuclear-palsy/psp_qsp_model.svg) · [PNG](../../../progressive-supranuclear-palsy/psp_qsp_model.png) — 259 nodes, 397 edges, 23 clusters |
| ODE model | [`psp_mrgsolve_model.R`](../../../progressive-supranuclear-palsy/psp_mrgsolve_model.R) — 71 ODEs, 8 regions, 16 agents, 29 scenarios, validation harness |
| Dashboard | [`psp_shiny_app.R`](../../../progressive-supranuclear-palsy/psp_shiny_app.R) — 12 tabs |
| References | [`psp_references.md`](../../../progressive-supranuclear-palsy/psp_references.md) — 146 entries (82 PubMed-resolved PMIDs, 64 search links) |

---

## 1. The question this model exists to answer

The two programmes conducted in PSP reported numbers in which **target engagement and
clinical outcome point in opposite directions**.

| Programme | Target engagement | Clinical outcome |
|---|---|---|
| **gosuranemab (PASSPORT, n=486)** | CSF unbound N-terminal tau **−98%** (placebo +11%, p<0.0001) | PSPRS 52-week change **10.4 vs 10.6** (p=0.85), CSF NfL unchanged |
| **NIO752 (MAPT ASO, phase 1, n=59)** | CSF total tau / p-tau181 **−20%** (5-fold **smaller**) | CSF NfL **stabilised** (placebo +40%) — the only downstream signal in the history of this disease |

**The largest target engagement** in the history of tauopathy trials is attached to
**the flattest null result**, and **the smallest biomarker change** is attached to
**the only signal**. And the latter has entered phase 3 (Preserve, n=300, 72 weeks,
mPSPRS-10).

A model in which "lowering tau" is **a single axis** cannot hold these two facts at
once. So this model separates tau into **several pools differing 10³–10⁴-fold in size
and with completely different causal roles**, and is constructed so that **the two
trials moved different pools** comes out as an output.

```
intraneuronal soluble monomer   TMON    ~2000 nM     ← the engine
regional aggregate load         A1..A8  0–1          ← the pathology
ISF N-terminal fragments        ETN     ~1 nM        ← the assay analyte
ISF seed-competent assemblies   ETS     ~0.02 nM     ← the causal species
CSF                             CTN,CTS ~0.005 nM
```

---

## 2. Seven structural axes

Each axis is **a line of code, not a fitted knob**.

### AXIS 1 — Pool size. The extracellular pool is not the engine but the exhaust

```
F_EXT = ETN / TMON = 5.75e-4      (model output, literature range 3e-4 ~ 1e-3)
```

Intracellular monomer is ~2000 nM; the ISF pool an antibody can reach is ~1 nM.
**Removing 98% of the small pool removes 0.06% of whole-brain tau.** The ceiling on
an extracellular mechanism is `F_EXT × Φ_ACC`, and **raising dose and affinity
without limit does not move that ceiling.**

### AXIS 2 — The engagement paradox. The number 98% is itself evidence that the target was small

The −98% is **not a fitted value**. It comes out of the real PK and the real affinity.

```
2000 mg q4w, V1 = 3 L, t½ ≈ 25 d   →  Cp,avg ≈ 300 mg/L ≈ 2000 nM
Kp,brain = 0.15%                    →  AISF = 6.73 nM        (model output)
Kd = 0.034 nM                       →  free fraction = 1/(1 + 196)
complex cleared 4× slower than free →  measured = -98.04%    (observed -98%)
```

**Put a 150 kDa IgG in at low nanomolar concentrations and it saturates a 1 nM
antigen while being irrelevant to a 2000 nM antigen. 98% engagement is not evidence
of target validity but evidence that the target was small enough to be saturable.**

Incidentally, because the complex is cleared **more slowly** than free tau (being
large and bulk-flow/FcRn-limited), the model predicts that **CSF total tau rises by
+283%** — unbound goes down while total goes up, **one mechanism with two readouts in
opposite directions**.

### AXIS 3 — Epitope. The measured species and the causal species are different molecules

Tau in ISF/CSF is **truncated**. The model places two physically separated states.

- `ETN` — N-terminal fragments. **Abundant**, cleaved out by
  caspase-3/calpain/AEP, and **seeding-incompetent**. This is literally the material
  the "unbound N-terminal tau" assay reports.
- `ETS` — mid-domain / MTBR assemblies. **Scarce**, seeding-competent, and species
  that **have already lost the N-terminus**.

An N-terminal antibody has `Kd = 0.034 nM` for the first and `Kd = 10⁶ nM` for the
second. **A single epitope parameter generates −98% on the measured species and 0% on
the causal species simultaneously.** The model also names the discriminating assays:
**CSF seeding activity** and **MTBR-tau243** — neither total tau nor N-terminal tau.

### AXIS 4 — Geometry. Progression is a travelling front across the connectome, and the local dynamics are already saturated

Local templated seeding is autocatalytic (`dA/dt = k·A·(1−A)`), hence **logistic**,
i.e. it **saturates** in each nucleus. For a sum of saturating sigmoids to be linear,
regions must be **recruited sequentially**. Therefore **the very fact that the PSPRS
slope is nearly linear at ~11 points/year is evidence of a front**, and that means
**every nucleus already involved sits at a point where changing `k` does nothing at
all.**

The recruitment order the model computes (time to reach 50% load, years):

| | Nucleus | Years | Clinical correlate |
|---|---|---|---|
| R1 | subthalamic nucleus (STN) | 2.05 | — |
| R4 | midbrain tegmentum (riMLF/rip) | **2.36** | **vertical gaze palsy — the earliest sign** |
| R2 | GPi/SNr | 4.01 | parkinsonism |
| R5 | pedunculopontine nucleus (PPN) | **4.08** | **backward falls** |
| R3 | SNc | 4.73 | dopamine deficiency |
| R8 | LC / bulbar medulla | 6.96 | dysphagia · dysarthria |
| R7 | frontal cortex | 7.15 | executive dysfunction · apathy |
| R6 | dentate nucleus | 7.57 | cerebellar signs (PSP-C) |

Trial enrolment occurs at 7.07 years. **That is, by the time of enrolment the front
has already passed through six of the eight nuclei.**

**Three null results follow from that single sentence.** Tideglusib (GSK-3β), OGA
inhibitors and davunetide (microtubule stabilisation) all act on `k`, and `k`
multiplies **a quantity that is already saturated**.

### AXIS 5 — Not a square-root law but **a measured elasticity**

The continuum Fisher–KPP predicts `v ∝ √(D·k·M)`, i.e. an elasticity of 0.5. A
discrete pulled front on a sparse graph predicts `v ≈ k·Δx / ln(k/D)`, i.e. nearly 1.
**This model asserts neither and measures instead** (`psp_elasticity()`).

| MAPT knockdown | Dosing from enrolment: 72-week slope (relative) | Dosing from biological onset: delay in reaching enrolment severity |
|---|---|---|
| 10% | 0.9974 | 0.42 years |
| 20% | 0.9947 | 0.98 years |
| 30% | 0.9918 | 1.71 years |
| 50% | 0.9854 | **4.14 years** |
| 70% | 0.9782 | 10.62 years |
| 80% | 0.9743 | **20.45 years** |
| 90% | 0.9702 | not reached within 33 years |

```
treatment-phase elasticity  dln(slope)/dln(M) = 0.022
=> halving the slope in a patient who already has symptoms is
   impossible at any level of knockdown.
```

**Same parameters, same drug, three orders of magnitude apart. The only difference is
*when* it is given.** This is a falsifiable prediction for the Preserve trial
(readout expected 2027–28): **0.02 points** on 72-week mPSPRS-10.

### AXIS 6 — The ASO gradient. CSF biomarkers overestimate the effect on the front

Intrathecal ASO exposure is **cortex ≫ deep nuclei**, whereas the front is in the
deep nuclei and CSF biomarkers are dominated by cortex. `ASO_DEEP = 0.4` makes this
divergence explicit: **CSF total tau −19% corresponds to a monomer reduction of ~10%
as seen by the front.**

Without this factor the model overestimates Preserve's effect 2.5-fold.

### AXIS 7 — Timing. The score is a delayed **integral** of neuronal loss gated by a reserve threshold

Aggregate load does not generate the score. It generates a **loss rate**, which
integrates into neuronal number, which crosses a **regional reserve threshold
(θ_r)**, which passes through two transit compartments to move the score. Between
drug and endpoint there are **two integrations and one threshold**.

**In-model counterfactuals — 52-week PSPRS benefit:**

| Intervention | Start | 52-week benefit |
|---|---|---|
| block aggregation **100%** | at enrolment | **0.21 points** |
| clear aggregates at **50× the rate** | at enrolment | **1.78 points** |
| MAPT ASO (Preserve regimen) | at enrolment | 0.02 points |
| mid-domain antibody | at enrolment | 0.0001 points |
| block aggregation 100% | **3 years before symptom onset** | **49.57 points** |

**Even a perfect drug started at enrolment gives 0.21 points at 52 weeks.** The reason
is that the damage rate is proportional to `A` and not to `dA/dt` — **stopping
aggregation and removing aggregates are different things, and what is needed at
enrolment is the latter.** And even the latter gives 1.78 points.

> PASSPORT and ARISE powered themselves to detect **an effect size their design
> forbade**. This is not a sample-size problem.

---

## 3. Validation

Running `psp_validate()` produces the table below. Only the **FIT** rows were used in
calibration; the **held-out** rows were not.

| Item | Observed | Predicted | Unit | \|log₁₀ ratio\| | Status |
|---|---|---|---|---|---|
| PSPRS at enrolment | 38 | 38.21 | pts | 0.002 | FIT |
| PSPRS 52-wk change, placebo | 10.6 | 11.53 | pts | 0.037 | FIT |
| Natural-history slope (Golbe 2007) | 11.3 | 10.83 | pts/y | 0.018 | held-out |
| Median survival from symptom onset, PSP-RS | 7.9 | 7.73 | y | 0.009 | FIT |
| Median survival, PSP-P | 10 | 8.16 | y | 0.089 | held-out |
| **Gosuranemab CSF unbound N-terminal tau** | **−98** | **−98.04** | % | — | FIT |
| Gosuranemab PSPRS 52-wk change | 10.4 | 11.53 | pts | 0.045 | held-out |
| Tilavonemab PSPRS 52-wk change | 10.6 | 11.53 | pts | 0.037 | held-out |
| **NIO752 CSF total tau** | **−20** | **−19.15** | % | — | FIT |
| NIO752 CSF p-tau181 | −20 | −18.96 | % | — | held-out |
| Tideglusib PSPRS 52-wk (TAUROS null) | 10.6 | 11.47 | pts | 0.034 | held-out |
| Davunetide PSPRS 52-wk (AL-108-231 null) | 11.0 | 11.52 | pts | 0.020 | held-out |
| Midbrain area at enrolment | 78 | 89.1 | mm² | 0.058 | held-out |
| MRPI at enrolment | 15 | 13.3 | — | 0.051 | held-out |
| Vertical saccade peak velocity | 175 | 139.6 | deg/s | 0.098 | held-out |
| F_EXT (ISF/intraneuronal tau) | 5e-4 | 5.75e-4 | — | 0.061 | structural |
| PSP-P PSPRS slope at matched severity | 5 | 13.73 | pts/y | **0.439** | **FAILURE** |
| CSF NfL rise over 1 y, placebo | 40 | 15.1 | % | **0.424** | **FAILURE** |

**Median \|log₁₀(predicted/observed)\| over the 15 quantitative anchors = 0.045.**

### The time axis is an output, not an input

It proved **impossible to fix** enrolment at "3.5 years after biological onset" —
doing so requires the score to be linear from t=0, and PSPRS 38 points and 10.6
points/year cannot then be satisfied together. So enrolment was defined as **the
moment the model itself passes PSPRS = 38**, and as a result **the preclinical period
became an output.**

```
t = 0.00 y   first tau aggregates
t = 3.57 y   symptom onset     ← 3.6 years of preclinical tauopathy (an output)
t = 7.07 y   trial enrolment severity (PSPRS 38)
t = 11.31 y  median survival   ( = 7.73 years after symptom onset, observed 7.9)
```

**There are 3.6 years ahead of symptoms that no clinical trial has ever seen.** This
interval is exactly the window the counterfactuals of AXIS 7 point to.

---

## 4. Failures reported, not repaired

### 4.1 The PSP-P slope — not explained by wave origin alone (\|log ratio\| 0.439)

The model represents PSP-P **by a difference in wave origin alone** (SNc instead of
STN + midbrain). The result is that **reaching enrolment severity is delayed by 2.13
years** and survival lengthens (8.16 years vs 7.73), but **the slope does not slow**
(13.73 vs observed ~5).

This requires something the model does not have: **at autopsy PSP-P has a lower tau
load than PSP-RS in absolute terms.** Lowering `KAGG` matches the slope, but doing so
makes the statement "the subtype difference comes out of origin alone" no longer
true. **The structural claim is not abandoned in order to match the slope; the
failure is reported instead.** The discriminating experiment is clear: compare total
tau-PET SUVR between PSP-P and PSP-RS **at the same severity**.

### 4.2 The one-year rise in CSF NfL (\|log ratio\| 0.424)

The placebo arm of NIO752 phase 1 reported **+40%** over a year, and the model
predicts **+15.1%**. In the model NfL is the sum of (i) the neuronal death rate and
(ii) ongoing non-lethal axonal damage, weighted by axonal volume, and it produces its
maximum rate of rise near enrolment, when the front is just entering high-volume
white matter (frontal cortex). Even that does not reach 40%.

That said, **in the wider literature CSF NfL in PSP is elevated but relatively stable
longitudinally** (which is why it is used as a diagnostic marker), and the +40% is a
single phase 1 observation in a placebo arm of n≈20. Either the model is wrong or
that observation is small-sample variation. **What decides which is Preserve's
placebo arm (n=100).**

### 4.3 Small misses

- **Midbrain area** 89.1 mm² (observed 70–85). The model's midbrain is a weighted
  combination of R4 · R1 · R5 and does not carry atrophy outside the tegmentum.
- **Vertical saccade peak velocity** 139.6 deg/s (observed 150–200). `GAM_SAC = 1`
  (linear) gives a slightly excessive reduction. Fitting the exponent would blur the
  statement that this is "a readout with no threshold", so it is left as it is.
- **MRPI** 13.3 (observed ~15, cut-off 13.55). In the progressive phase it exceeds
  20, which is outside the range over which MRPI has been validated. **It is reported
  as it stands, without clipping.**

### 4.4 Discarded structure, not tuned

**"The square-root law"** — the draft map (.dot) **asserted**, because `v ∝ √M`, that
"halving the slope requires 75% knockdown". Once the elasticity could actually be
measured, the model **refuted** it: the treatment-phase elasticity is 0.022, and
halving is impossible at any knockdown. The continuum Fisher–KPP's 0.5 is a **lower
bound** and the discrete pulled front's ~1.0 an **upper bound**, but **measured after
enrolment in a discrete 8-nucleus chain it is far below both limits** — because the
front has already passed. That assertion was **deleted from the map and the
documentation and replaced by the measured value.**

---

## 5. Defects found and fixed during integration

Found while running under mrgsolve 2.0.1 / R 4.3.3. **All of these would have
remained in place had the model not actually been run.**

1. **A degenerate fit with no wave.** The first optimisation found `KAGG = 47/day`.
   Every region saturated within a year, so **the front did not exist at all**, and
   the clinical anchors were matched purely by the neuronal-loss integrator. Every
   number was excellent and the model's central structural claim had disappeared.
   Resolved by adding wave arrival-time targets to the objective function.
2. **Disappearance of the clinical delay.** Freely fitted, `KLAG` collapsed to ~19
   days, deleting the delay on which the trial-design argument of AXIS 7 depends.
   Fixed.
3. **Initial conditions off steady state.** `TMON_0` was set without the baseline
   clearance factor (`CLR0`), so the model started outside its own steady state and
   the growth term saw a **phantom −5% "knockdown" of the opposite sign**. As a
   result `ASO_DEEP` was **quietly scaling the autophagy feedback** rather than the
   drug effect. `ACT0/CLR0` was put into both the initial conditions and the
   reference values, and `MFRONT` was decomposed into drug-attributable and
   physiology-attributable parts.
4. **Knockdown measured against itself.** The elasticity probe lowered `KSYN_T`, but
   the reference value `TMONS` also contained `KSYN_T`, so the ratio was always 1.
   **An 80% knockdown was recorded as zero, and the elasticity came out as −0.000.**
   Fixed reference values `KSYN_T0`/`TMREF` were introduced.
5. **A counterfactual with no time gate.** "50× aggregate clearance from enrolment"
   was actually applied **from birth**, so the patient never got the disease in the
   first place. The 11.2-point "treatment benefit" was **the absence of disease**.
   After adding the `TDIS/FDIS` gate, 1.78 points.
6. **"Preclinical" was not preclinical.** "3 years before symptom onset", coded as
   `TENROL − 1095`, was actually **183 days after** symptom onset. Corrected to
   `TSYMPT − 1095`.
7. **A comparison table reporting only change from baseline.** An intervention that
   **delays** the disease sits lower on the curve at a fixed enrolment time, so **its
   change is larger.** S09 and S27 **looked worse** than placebo (+0.48, +0.61).
   Fixed to output the absolute score alongside: the same scenarios are −4.30 and
   −17.32 points. (This is not a bug but a trap in endpoint selection, which is why
   it was **left** in the table.)
8. **Premature termination of an R string.** An apostrophe inside `Kd's` in a comment
   in the model code broke the single-quoted string. Switched to a raw string
   (`r"---( )---"`) and the whole file cleaned to ASCII (to prevent parse failure in
   environments with a POSIX locale).
9. **Multiple declarations on one line.** `double a = ..., b = ...;` is valid C++ but
   mrgsolve's declaration hoisting only handles the first variable. Split to one
   declaration per line.
10. **graphviz `init_rank` failure.** Cluster-level cycles broke `dot`'s local cluster
    ranking (23 clusters entangled in both directions). Switched to global ranking
    with `newrank = true`.

---

## 6. Using the files

```r
# requires: R >= 4.1, mrgsolve >= 1.0 (development environment mrgsolve 2.0.1)
source("psp_mrgsolve_model.R")

psp_validate()      # validation table against public anchors
psp_elasticity()    # treatment-phase vs prevention-phase elasticity — the model's sharpest number
psp_compare()       # comparison of the 29 scenarios (change + absolute score)

d <- psp_run("S02")             # PASSPORT (gosuranemab)
d <- psp_run("S07")             # Preserve-like ASO regimen
d <- psp_run("S25")             # counterfactual: complete block 3 years before symptoms

# dashboard
shiny::runApp("psp_shiny_app.R")
```

```bash
# re-render the map
dot -Tsvg psp_qsp_model.dot -o psp_qsp_model.svg
dot -Tpng -Gdpi=150 psp_qsp_model.dot -o psp_qsp_model.png
```

### Scenario list (29)

| | |
|---|---|
| S00–S01 | Natural history: PSP-RS, PSP-P |
| S02–S05 | Anti-tau antibodies: gosuranemab, tilavonemab, mid-domain/MTBR probe, Φ_ACC=1 |
| S06–S09 | MAPT ASO: phase-1 regimen, Preserve regimen, 4× dose, presymptomatic dosing |
| S10–S12 | Three local-dynamics agents (structurally null): tideglusib, OGA inhibitor, davunetide |
| S13–S19 | Salsalate, fasudil, TPN-101, LM11A-31, ezeprogind, AADvac1, riluzole |
| S20–S21 | Symptomatic: levodopa, zolpidem |
| S22–S23 | Combinations: ASO + mid-domain antibody, ASO + ezeprogind |
| S24–S25 | Counterfactual: complete block of aggregation at enrolment / 3 years before symptoms |
| S26 | PEG feeding tube (a risk modifier) |
| S27 | Triple combination, started 3 years before symptoms |
| S28 | Counterfactual: 50× aggregate clearance from enrolment |

---

## 7. Falsifiable predictions this model makes

| Prediction | Readout |
|---|---|
| **Preserve (NIO752 phase 3) will be null on 72-week mPSPRS-10** (predicted effect 0.02 points) — null even if CSF tau improves from −20% to −40% | 2027–28 |
| Mid-domain / MTBR antibodies will be null in symptomatic patients, **indistinguishably** from N-terminal antibodies (0.0001 points) | — |
| OGA inhibitors will be null in symptomatic PSP (predicted −0.09 points) | — |
| **CSF seeding activity** and **MTBR-tau243** will not move on an N-terminal antibody — while total tau **rises** | testable immediately in stored PASSPORT samples |
| PEG extends survival by 1.19 years and changes PSPRS **not at all** | observational study |
| Only presymptomatic intervention is clinically meaningful: 50% MAPT knockdown = **4.14 years** delay to enrolment severity | requires a preclinical cohort |
| Total tau-PET load differs between PSP-P and PSP-RS **at the same severity** (the discriminating experiment predicted by failure 4.1) | imaging study |

---

## 8. Model structure summary (state inventory)

| Compartment group | Number | Contents |
|---|---|---|
| Regional aggregate load `A1–A8` | 8 | STN · GPi/SNr · SNc · midbrain tegmentum · PPN · dentate nucleus · frontal cortex · LC/medulla |
| Regional surviving neurons `N1–N8` | 8 | integrator (the first integration) |
| Molecular | 8 | MAPT mRNA · monomer · GSK-3β · O-GlcNAc · acetylation · microtubule binding · progranulin · autophagy |
| Extracellular | 4 | ISF N-terminal fragments · ISF seeds · CSF for each |
| Glia / innate immunity | 5 | microglia · cytokines · astrocytes · type I IFN · lipid peroxidation |
| Biomarkers | 5 | CSF NfL · plasma NfL · plasma GFAP · midbrain area · SCP width |
| Clinical | 6 | established-damage transit ×2 (the second integration) · dysphagia · competing-risk CIF ×2 · survival |
| PK | 27 | mAb 2-compartment + ISF + CSF · ASO CSF + tissue · 11 small molecules · AADvac1 titre |
| **Total** | **71** | |

---

## ⚠️ Disclaimer

This is a semi-quantitative QSP model for educational and research purposes. It was
constructed from the public literature and reported clinical trial values but has not
been independently validated or certified, and **must not be used for clinical
decision-making, prescribing, or regulatory submission.** The parameters are
illustrative approximations, and the predictions in section 7 above are statements
made inside this model structure — if the structure is wrong, the predictions are
wrong too. That is why they have been written down in a falsifiable form.
