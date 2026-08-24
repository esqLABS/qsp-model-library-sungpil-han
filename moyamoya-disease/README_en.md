# Moyamoya Disease (MMD) — QSP Model

> **Thesis.** In moyamoya disease the *measured* quantity is cerebral blood flow
> (CBF), but the *state variable* is how much dilatory reserve the cortical
> arteriole has left. The arteriole is a **variable resistor with a floor**, and
> almost everything clinically distinctive about this disease is derived, as
> arithmetic, from what happens when that floor is reached.

---

## 1. Model structure

A hemisphere is written as a **3-node resistance network**. Two native inflow
routes (the antegrade internal carotid · the moyamoya perforators), one dangerous
collateral (the periventricular/choroidal anastomosis), and the surgical graft
all feed the same cortical arteriolar bed **in parallel**.

```
      Pa --[ g_ICA·(1-STEN)^4 ]--+
      Pa --[ g_moya            ]--+--> P_A --[ g_art ]--> Pv     (territory)
      Pa --[ g_PVA             ]--+       ↑ FLOOR
                                  |    [ g_leak ]
      Pa --[ g_bypass          ]--+--> P_F --[ g_artF ]--> Pv    (peri-anastomotic)
                                  |
                           [ g_coll ]        ← the safe collateral donated by the PCA
                                  |
      Pa --[ g_PCA ]----------> P_B --[ g_artB ]--> Pv           (donor)
```

Node balance equations:

```
A: gSA(Pa−P_A) + gc(P_B−P_A) + gl(P_F−P_A) = gA(P_A−Pv)
F: (gSF+gb)(Pa−P_F) + gl(P_A−P_F)          = gF(P_F−Pv)
B: gp(Pa−P_B) + gc(P_A−P_B)                = gB(P_B−Pv)
```

Autoregulation chooses `gA, gF, gB` so as to meet the demand of each
compartment, but is **clamped to `[1/R_art_max_eff, 1/R_art_min]`.** The entire
model is this one clamp. The remaining 37 ODEs (the lesion · the two collaterals
· the graft · arteriolar remodelling · angiogenic signalling · injury ·
haemorrhage · the PK of seven agents) hang from this algebraic layer.

---

## 2. Deliverables

| File | Contents |
|------|------|
| [`mmd_qsp_model.dot`](mmd_qsp_model.dot) · [`.svg`](mmd_qsp_model.svg) · [`.png`](mmd_qsp_model.png) | Mechanistic map — **126 nodes / 11 clusters** |
| [`mmd_mrgsolve_model.R`](mmd_mrgsolve_model.R) | mrgsolve model — **37 ODE compartments · 24 scenarios**, including the C++ 3-node network solver |
| [`mmd_shiny_app.R`](mmd_shiny_app.R) | Shiny dashboard — **12 tabs** |
| [`mmd_references.md`](mmd_references.md) | **99 references** (every PMID verified through the NCBI E-utilities) |
| [`mmd_reference_model.py`](mmd_reference_model.py) | Independent Python/scipy reimplementation (the numerical ground truth) |
| [`mmd_reference_output.txt`](mmd_reference_output.txt) | The run output that generated every number below |
| [`mmd_cross_validation.txt`](mmd_cross_validation.txt) | mrgsolve ↔ Python cross-validation results (80 paired values) |
| [`mmd_cross_validation.py`](mmd_cross_validation.py) · [`.R`](mmd_cross_validation.R) | The two-stage script that reproduces the cross-validation above |
| [`mmd_r_validation.txt`](mmd_r_validation.txt) | R model self-verification · JAM reproduction · MAP optimisation · virtual population output |

### Verification

Every equation was implemented **twice**: once in mrgsolve/C++, once in
numpy/scipy LSODA. Over 80 paired values the two implementations agree to
**a median of 0.000 %, a maximum of 0.15 %, and 80/80 within 1 %**
(`mmd_cross_validation.txt`).

What was and was not actually run is set down exactly (R 4.3.3,
mrgsolve 2.0.1):

| Item | Status |
|---|---|
| `mmd_reference_model.py` | Run — generated all 15 sections of `mmd_reference_output.txt` |
| `mmd_mrgsolve_model.R` compilation | Run — the C++ model built successfully |
| R self-verification · JAM reproduction · MAP optimisation · virtual population | Run — `mmd_r_validation.txt` |
| Python ↔ mrgsolve cross-validation | Run — `mmd_cross_validation.txt` |
| `mmd_shiny_app.R` parsing | Run — parses as valid R code |
| Shiny **server-side computation for all 12 tabs** | Run — the simulations, sweeps, probes and population computations each tab uses were called directly and the values checked |
| Shiny **app launch (browser rendering)** | **Not done** — installation of the `shiny` package fails in this container (the dependency `fs` fails to configure). This is an environmental constraint and not a defect of the model. The UI layout has not been confirmed in a real browser. |

---

## 3. Results that follow as arithmetic

Every number below stands as it is in `mmd_reference_output.txt`.

### (1) The infarction threshold is not a parameter — it is derived

```
CBF_crit = CMRO2 / (CaO2 · OEF_max)
```

At Hb 15 g/dL, **19.71 mL/100g/min** — the textbook "penumbra 20" simply drops
out of the oxygen arithmetic. It is not a fitted value. And it is at the same
time the explanation of sickle-cell moyamoya.

| Hb (g/dL) | CaO₂ | CBF_crit | vs Hb 15 |
|---|---|---|---|
| 15 | 0.1970 | 19.71 | 1.00× |
| 11 | 0.1445 | 26.88 | 1.36× |
| 8 | 0.1051 | **36.96** | **1.88×** |
| 7 | 0.0919 | 42.23 | 2.14× |

Changing Hb alone in five patients whose vascular lesions are **completely
identical** (5 years):

| Management | Hb | CBF_crit | CBF_ws | OEF | 5-yr infarction | 5-yr ischaemic events |
|---|---|---|---|---|---|---|
| No transfusion | 8 | 37.0 | 19.6 | 0.740 | **28.6 %** | 0.681 |
| Up to Hb 9 | 9 | 32.8 | 19.6 | 0.659 | 26.4 % | 0.509 |
| Chronic transfusion 11 | 11 | 26.9 | 19.3 | 0.545 | 17.6 % | 0.316 |
| Aggressive transfusion 13 | 13 | 22.7 | 19.0 | 0.466 | 6.8 % | 0.229 |
| Normal-Hb control | 15 | 19.7 | 18.9 | 0.406 | **1.6 %** | 0.182 |

**Why a transfusion, which touches nothing vascular at all, is an effective
treatment** — it alone moves not the "flow" but the "threshold that flow is
compared against" (STOP/SIT).

### (2) A critical inflow conductance gS\* exists, and three things switch at once

`gS* = 2.1995 mL/min/mmHg`. Above it, CBF is 50 whatever the pressure does and
`dCBF/dMAP = 0`. Below it the arteriole is pinned to its floor, CBF becomes a
linear function of MAP, and **the measured acetazolamide response inverts to
negative.** Three things change at the same point.

| gS | Equivalent stenosis | P_A | CBF_A | Arteriolar position | Intrinsic reserve % | **Measured response %** | dCBF/dMAP |
|---|---|---|---|---|---|---|---|
| 6.25 | 0.000 | 70.0 | 50.0 | 0.400 | 84.7 | +45.0 | 0.000 |
| 3.00 | 0.168 | 48.7 | 50.0 | 0.620 | 24.4 | +22.4 | −0.000 |
| **2.20** | **0.230** | 34.0 | 50.0 | **0.999** | **0.0** | **−0.2** | 0.000 |
| 1.80 | 0.268 | 30.8 | 43.4 | 1.000 | 0.0 | −0.3 | **0.546** |
| 1.00 | 0.368 | 23.3 | 27.7 | 1.000 | 0.0 | −0.6 | 0.350 |

**gS\* is a property of the network, not a property of the stenosis.** Converted
into an isolated internal carotid stenosis it is a diameter of 0.230, but a
patient with good collaterals reaches the same gS only at a far more severe
stenosis. **The Suzuki grade cannot mark this point.**

### (3) The silent phase is a consequence of the fourth power — reserve disappears three years before flow does

Stenosis is linear in time, but conductance is its fourth power. The natural
history of the adult ischaemic type:

| Year | STEN | g_ICA | gS | CBF_A | CBF_ws | OEF | **Intrinsic reserve %** | Infarction % |
|---|---|---|---|---|---|---|---|---|
| 1 | 0.125 | 3.657 | 3.66 | 50.0 | 49.8 | 0.335 | **40.6** | 0.00 |
| 2 | 0.234 | 2.150 | 2.16 | 49.4 | 48.8 | 0.339 | **0.0** | 0.00 |
| 3 | 0.353 | 1.094 | 1.57 | 41.4 | 33.1 | 0.405 | 0.0 | 0.09 |
| 5 | 0.559 | 0.235 | 1.40 | 39.7 | 21.0 | 0.422 | 0.0 | 0.85 |
| 10 | 0.816 | 0.007 | 1.27 | 39.6 | 17.9 | 0.423 | 0.0 | 6.88 |

**Reserve becomes 0 in year 2, yet flow does not move until years 3–5. That gap
is this disease.** Measuring CBF in the asymptomatic stage cannot find this
patient.

**And the steady state is "flow nailed to the threshold":**

| Year | CBF_ws | CBF_crit | Infarction % (= demand removed) |
|---|---|---|---|
| 5 | 21.0 | 19.7 | 0.85 % |
| 7 | 18.5 | 19.7 | 2.27 % |
| 10 | 17.9 | 19.7 | 6.88 % |

Because infarction **removes demand**, the flow per surviving gram comes back to
near the threshold. That the SPECT of chronic moyamoya looks "only slightly
abnormal" is *because* tissue has been lost. **Reading CBF alone systematically
underestimates the severity of this disease.**

### (4) Acetazolamide and PaCO₂ are not the same test

Acetazolamide acts on the arteriole alone (tissue carbonic anhydrase). PaCO₂ acts
on the arteriole **and on the leptomeningeal collateral conduits**. In the normal
brain this distinction is invisible. In moyamoya **the conduit is the
circulation**, so the two tests give different answers.

| State | Intrinsic reserve % | CBF_A | ACZ 1 g → ΔA % | ΔB % | PaCO₂ 50 → ΔA % | PaCO₂ 25 → ΔA % |
|---|---|---|---|---|---|---|
| Normal (y1) | 40.6 | 50.0 | **+38.1** | +39.0 | +35.0 | −51.8 |
| Reserve exhausted (y4) | 0.0 | 40.0 | **−1.7** | +39.0 | **+4.4** | −17.5 |
| Decompensated (y6) | 0.0 | 39.5 | **−1.9** | +39.0 | +6.2 | −20.0 |
| End stage (y10) | 0.0 | 39.6 | **−2.0** | +39.0 | +6.8 | −20.5 |

The donor territory (B) answers every time. The diseased territory (A) stops
answering, and then answers **with the wrong sign** — since B is A's supply,
dilating B steals from A. **Steal is a consequence of the network, not an
additional assumption.** Hypercapnia, conversely, helps by opening the conduit
itself — which is why the anaesthetic target is normocapnia to mild hypercapnia,
and why an acetazolamide test and a breath-hold test cannot be swapped for one
another.

### (5) The crying child: the same change in PaCO₂, the opposite outcome

| Patient | CBF_ws @40 | CBF_ws @25 | Ratio | CBF_crit | Threshold crossed? |
|---|---|---|---|---|---|
| Normal | 50.0 | 26.3 | 0.526 | 19.7 | No |
| Reserve preserved | 49.8 | 23.7 | 0.476 | 19.7 | No |
| Adult, reserve exhausted | 24.8 | 21.2 | 0.852 | 19.7 | No |
| Adult, decompensated | 19.3 | 15.7 | 0.811 | 19.7 | **already below** |
| Child, decompensated | 18.7 | 14.9 | 0.798 | 19.7 | **already below** |

The normal brain surrenders **44 %** of its flow to a 15 mmHg fall in PaCO₂ and
is still entirely safe, because flow was never the constraint in the first place.
The decompensated hemisphere infarcts having lost **less in percentage terms**,
because it had no margin.

### (6) The blood-pressure paradox has a computable optimum, and it differs by phenotype

Ischaemic risk flows through the pressure-passive term, haemorrhagic risk through
the periventricular wall stress `σ ∝ P_perf · dilatation^1.5`. The two have
**opposite signs with respect to MAP**.

| MAP | Ischaemic type: ischaemia/yr | Haemorrhage/yr | Sum | Haemorrhagic type: ischaemia/yr | Haemorrhage/yr | Sum |
|---|---|---|---|---|---|---|
| 70 | 0.1350 | 0.0030 | 0.1381 | 0.1284 | 0.0478 | 0.1762 |
| 85 | 0.0951 | 0.0044 | 0.0995 | 0.0912 | 0.0685 | **0.1597** |
| **90** | 0.0855 | 0.0049 | 0.0903 | 0.0822 | 0.0762 | **0.1584 ← optimum** |
| 110 | 0.0592 | 0.0071 | 0.0663 | 0.0503 | 0.1137 | 0.1640 |
| 125 | 0.0478 | 0.0102 | **0.0580 ← optimum** | 0.0448 | 0.1591 | 0.2039 |

The optimal MAP for the haemorrhagic type is **90 mmHg**, for the ischaemic type
the upper end of the sweep (125 mmHg, the point at which CBF saturates at 50) —
that is, within the physiological range the ischaemic type becomes "the higher
the better", while the haemorrhagic type acquires a distinct interior optimum.
**"Control the blood pressure" is not a moyamoya instruction until it says which
hemisphere it is talking about.**

### (7) A bypass treats both phenotypes by one mechanism (the JAM reproduction)

| Group | CBF_A | CBF_ws | P_A | **Q_pva** | g_pva | σ_pva | 5-yr haemorrhage | 5-yr ischaemia |
|---|---|---|---|---|---|---|---|---|
| Haemorrhagic / conservative | 40.7 | 18.4 | 28.3 | **31.6** | 0.513 | 1.144 | **0.310** | 0.336 |
| Haemorrhagic / direct bypass | 41.2 | 23.2 | 29.0 | **23.1** | 0.379 | 1.088 | **0.165** | 0.224 |
| Ischaemic / conservative | 39.6 | 17.9 | 27.7 | 11.3 | 0.181 | 0.937 | 0.023 | 0.344 |
| Ischaemic / direct bypass | 40.2 | 22.7 | 28.8 | 9.8 | 0.160 | 0.924 | 0.016 | 0.223 |

The mechanism is visible in the `Q_pva` column: **as P_A rises, the gradient
across the periventricular anastomosis collapses, VEGF falls, and that bed is
pruned.** Surgery does not reinforce the vessel that was going to burst — it
**retires** it.

> **Reported as an honest miss.** In the JAM trial (Miyamoto 2014, Stroke, PMID
> [24668203](https://pubmed.ncbi.nlm.nih.gov/24668203/)) 5-year rebleeding was
> 31.6 % conservative vs 11.9 % bypass, HR 0.355. **Only the conservative arm's
> proportion** was used to set `HEM_HAZ0`, so the HR is a prediction, and the
> model gives 31.0 % → 16.5 %, **HR 0.485**, thereby **underpredicting** the
> benefit by about a third. The most plausible reason: the model routes the whole
> of the benefit through the *pruning* of the periventricular bed, whereas a real
> bypass will already be lowering the perforator pressure faster than that bed
> involutes.

### (8) The hyperperfusion syndrome is focal and *relative*, and it is priced on the preoperative ischaemia

**A lumped territory cannot hyperperfuse.** Even with a patent graft in place the
territorial arteriole is still *dilating*, so the mean territorial CBF cannot
exceed demand. The syndrome must therefore necessarily be **focal** (the cortex
the graft was sewn to, `P_F ≫ P_A`) and **relative** (against the flow that
cortex had adapted to). 50 mL/100g/min is normal, and it damages a barrier that
has lived at 18 for five years.

| Preoperative hemisphere | REMOD | Intrinsic reserve % | CBF_F before | CBF_F peak | Relative peak | Day of peak | Days > 1.35 | Oedema |
|---|---|---|---|---|---|---|---|---|
| Reserve preserved | 0.000 | 31.4 | 50.0 | 50.0 | **1.00** | — | **0.0** | 0.000 |
| Reserve just exhausted | 0.593 | 0.0 | 38.1 | 65.7 | 1.41 | 7.7 | 9.0 | 0.166 |
| Decompensated 2 years | 0.593 | 0.0 | 37.1 | 65.2 | 1.43 | 7.8 | 10.7 | 0.226 |
| Decompensated 5 years | 0.593 | 0.0 | 35.6 | 64.6 | **1.46** | 8.1 | **13.4** | 0.332 |

A hemisphere with reserve left **does not hyperperfuse at all.** The syndrome
belongs not to the operating surgeon but to the **history of the arteriole**.
Vasoparalysis (`REMOD`) sets the **height** of the surge, and the barrier's
re-adaptation time constant sets its **duration** (peak 7–8 days, lasting 9–13
days — consistent with the clinically reported peak at 2–7 days and resolution
over 2–3 weeks).

**And the incidence is set by something the surgeon cannot see.** Sweeping the
leptomeningeal coupling `g_leak` between the focal cortex and the rest of the
territory trades territorial gain against focal surge **in opposite
directions**:

| g_leak | Q_byp | P_F−P_A | Territorial CBF_ws gain | Relative surge | Days > 1.35 | Syndrome? |
|---|---|---|---|---|---|---|
| 0.25 | 65.0 | 37.3 | +2.0 | **1.70** | 19.9 | **Yes** |
| 0.55 | 64.8 | 31.6 | +3.8 | 1.55 | 18.5 | Yes |
| 0.80 | 64.7 | 28.0 | +5.0 | 1.46 | 14.2 | Yes |
| 1.20 | 64.6 | 23.6 | +6.3 | 1.36 | 2.8 | No |
| 1.80 | 64.5 | 19.1 | **+7.7** | 1.34 | **0.0** | No |

`g_leak` is a property of the patient's leptomeningeal mesh, not a property of
the operation. Which is why this is **the reason hyperperfusion arises in only a
minority of technically perfect bypasses**, and the reason its occurrence is not
a measure of surgical quality. And it predicts a testable inverse correlation:
**the hemisphere that hyperperfused is the hemisphere that gained least.**

### (9) An indirect bypass is a gamble on angiogenic capacity — one parameter

| Patient | ANGIO | Operation | g_byp d90 | g_byp d365 | CBF_ws d365 | Infarction % d1095 |
|---|---|---|---|---|---|---|
| Child | 1.00 | none | 0.000 | 0.000 | 18.5 | 4.91 |
| Child | 1.00 | indirect | **1.020** | 1.150 | 23.4 | 2.67 |
| Child | 1.00 | direct | 1.050 | 1.050 | 23.2 | 2.67 |
| Adult | 0.58 | none | 0.000 | 0.000 | 21.0 | 2.27 |
| Adult | 0.58 | indirect | **0.861** | 1.147 | 26.1 | 1.46 |
| Adult | 0.58 | direct | 1.050 | 1.050 | 25.9 | 1.47 |
| Adult | 0.58 | combined | 1.787 | 2.195 | 27.8 | 1.31 |

The indirect construct has to be **grown by the patient himself**, and that
patient is precisely the patient in whom RNF213 has restricted angiogenesis.

### (10) PCA involvement changes the phenotype — one flag, two opposite signs

PCA involvement **abolishes the safe donor (leptomeningeal) and forces the
dangerous route (periventricular).**

| PCA involvement | g_coll | g_pva | σ_pva | ANEU | 10-yr ischaemia | **10-yr haemorrhage** |
|---|---|---|---|---|---|---|
| 0.00 | 0.332 | 0.181 | 0.937 | 0.013 | 0.404 | **0.031** |
| 0.50 | 0.232 | 0.360 | 1.066 | 0.063 | 0.399 | 0.156 |
| 1.00 | 0.133 | 0.505 | 1.137 | 0.131 | 0.396 | **0.316** |

Ischaemic risk is left almost as it was while **haemorrhagic risk alone moves
tenfold**. Haemorrhagic risk is not a separate disease with a cause of its own.
**It is what this circulation does when the collateral available to it is the
weak one.**

### (11) RNF213 enters at two sites with opposite signs

| Genotype | RNF | STEN 10y | Collateral ceiling | g_moya | gS | CBF_ws | Infarction % |
|---|---|---|---|---|---|---|---|
| Wild type / quasi-MMD | 0.0 | 0.810 | **1.972** | 1.014 | 1.27 | 18.3 | **5.11** |
| R4810K heterozygous | 1.0 | 0.816 | 1.144 | 1.080 | 1.27 | 17.9 | **6.88** |
| R4810K homozygous | 1.6 | 0.944 | 0.647 | 0.647 | 0.75 | 15.5 | **29.66** |

A gene that creates the obstacle and then forbids the detour is worse than either
effect on its own — **which is why RNF213 status is prognostic information over
and above the angiogram.** (The difference between wild type and heterozygote is
shallow, infarction 5.11 % vs 6.88 %: a higher ceiling is self-limiting because
the better perfusion means a smaller VEGF drive; the large separation appears in
the homozygote.)

### (12) Drugs move the risk terms and do not move the geometry

| Group | STEN | gS | CBF_A | CBF_ws | Infarction % | 5-yr ischaemia | 5-yr haemorrhage |
|---|---|---|---|---|---|---|---|
| None | 0.559 | 1.40 | 39.7 | 21.0 | 0.85 | 0.092 | 0.008 |
| Aspirin 100 mg/day | 0.560 | 1.41 | 39.6 | 21.0 | 0.39 | 0.070 | 0.008 |
| Cilostazol 200 mg/day | 0.560 | 1.42 | 39.6 | 21.0 | 0.71 | 0.086 | 0.008 |
| Atorvastatin 20 mg/day | 0.438 | 1.46 | 40.3 | 26.5 | 0.41 | 0.039 | 0.005 |
| **Nifedipine GITS 60 mg/day** | 0.594 | 1.50 | **34.1** | **17.3** | **2.15** | **0.152** | 0.005 |
| Minocycline 200 mg/day | 0.559 | 1.40 | 39.7 | 21.0 | 0.85 | 0.092 | **0.006** |
| Antihypertensive alone | 0.586 | 1.51 | 36.3 | 18.5 | 1.61 | 0.141 | 0.006 |
| Direct bypass (reference) | 0.540 | **2.12** | 40.6 | 26.3 | 0.47 | 0.047 | 0.005 |

**The only thing that changes gS meaningfully is the bypass** (1.40 → 2.12).
Nifedipine is an instructive failure: a genuine cerebral vasodilator that cannot
dilate any further a vessel already pinned to its floor, which instead dilates
the **donor territory** and lowers the MAP the pressure-passive territory was
relying on — making infarction 2.5-fold.

> **The gap in the literature is reported as it stands.** In moyamoya disease
> there is **no** randomised evidence for antiplatelet agents · statins ·
> calcium-channel blockers · minocycline. The drug arms above are not the
> reproduction of a trial but a **mechanistic hypothesis**, and they are marked
> as such.

### (13) When you measure is everything — reserve saturates, margin does not

Stratifying the virtual population of 400 (5 years) by **reserve** fails:
**99.0 %** of the cohort is already at zero reserve, so the quartile cut-points
come out as 0.00 % / 0.00 %. Reserve is extremely sensitive *while it is being
exhausted* (see (3): zero in year 2, flow normal) and carries no information
*once it has been exhausted*.

What still separates patients in established disease is the **margin to the
threshold**, `margin = CBF_ws − CBF_crit` (P10 −7.1 / median −0.3 / P90
+12.8 mL/100g/min). Give the same direct bypass to the worst and best quartiles
at year 5 and read year 8:

| Group | n | margin | STEN | gS | Infarction % no surgery | Infarction % bypass | Gain | Ischaemia no surgery | Ischaemia bypass | **ARR** |
|---|---|---|---|---|---|---|---|---|---|---|
| Worst margin (Q1) | 36 | −7.21 | 0.870 | 0.82 | 26.44 | 18.62 | **7.82** | 0.546 | 0.409 | **0.137** |
| Best margin (Q4) | 36 | +12.51 | 0.543 | 1.39 | 1.95 | 0.89 | **1.07** | 0.147 | 0.093 | **0.055** |

The absolute benefit of the same operation differs **about sevenfold**. But read
the `STEN` column honestly and in this stratification the angiogram is **not**
blind (stenosis 0.870 vs 0.543) — in this population the lesion and the
physiology are correlated and the DSA does contain real prognostic information.
So the defensible claim is the weaker one: **stenosis is information, but it is
not sufficient.** Because the same stenosis maps to different margins according
to collateral conductance and haemoglobin (see (3) · (1) · (11)).

---

## 4. What verification of this model refuted

- **A lumped territory cannot hyperperfuse.** The first implementation was a
  2-node territorial model, and even with a patent graft inserted the relative
  CBF peak came out at **1.00** — that is, no effect whatsoever. The cause was
  not a wrong assumption but the **resolution**: even after a patent graft,
  gS (≈2.4) falls short of the normal 6.25, so the arteriole is still dilating.
  A focal peri-anastomotic compartment (a third node) therefore had to be added,
  and as a result "hyperperfusion can only be a focal phenomenon" became a
  *prediction* of the model.
- **The definition of hyperperfusion was wrong as well.** Measured against the
  textbook normal (50), still nothing at all happens. Clinical hyperperfusion is
  a value referred to **the flow that cortex had adapted to**, and so the
  adaptation set-point (`CBFAD`, τ = 22 d) entered as a state variable. It is
  this that sets the duration of the syndrome.
- **The `map_mult` scenario bug.** The early "postoperative blood-pressure
  policy" experiment applied the MAP multiplier from day 0, and so was quietly
  rewriting the ten years of preoperative history; lowering MAP further made the
  surge appear to *grow*. Fixed with a policy start-day gate (`MAP_MULT_T`).
- **"Stratify the patients by reserve" does not work at year 5.** An attempt to
  stratify the virtual population of 400 into quartiles by intrinsic reserve at
  the 5-year time point gave cut-points of **0.00 % / 0.00 %** — because
  **99.0 %** of the cohort is already at zero reserve. The two "quartiles" became
  effectively the same set of patients, the numbers for the two groups printed
  out entirely identical, and the conclusion that had originally been written up
  ("the benefit is not the same") was refuted by its own output. Reserve is a
  **saturating variable**: extremely sensitive while it is being exhausted
  ((3) above: zero in year 2 while flow is normal), and almost devoid of
  information once exhausted. So this section was rewritten — the variable that
  still separates patients in established disease is the **margin to the
  threshold (CBF_ws − CBF_crit)**, and because `CBF_crit` depends on haemoglobin,
  **the same CBF means a different thing in different patients**. That is, this
  finding is not "measuring reserve is meaningless" but **"when you measure is
  everything"**.
- **A citation error.** The JAM trial PMID was written from memory as 24788972,
  but verification showed that to be a different Stroke 2014 paper. The correct
  value is **24668203**. One candidate in the reference list was also found to
  carry a PMID belonging to an unrelated paper (a calcaneal tumour) and was
  discarded. All 99 PMIDs in `mmd_references.md` were confirmed by reading title,
  journal and year back from NCBI.

---

## 5. Calibration anchors

| Item | Value | Use |
|---|---|---|
| Normal cortical CBF / CMRO₂ / OEF | 50 / 3.30 / 0.335 | geometry and oxygen constants |
| Large-vessel share of cerebrovascular resistance | 25 % | `FRAC_PROX` |
| Rise in normal CBF with acetazolamide 1 g IV | 30–40 % | `AZ_EMAX` 0.45 (model +38.1 %) |
| CO₂ reactivity | ~3.5 %/mmHg | `K_CO2` 0.035 |
| CBF and OEF on the affected side in symptomatic adult MMD | 15–25 %↓ vs normal, OEF 0.42–0.50 | collateral ceiling (model 39.6 / 0.423) |
| **JAM conservative arm, 5-year rebleeding** | **31.6 %** | `HEM_HAZ0` 0.414 — **the HR of the bypass arm is a prediction** |
| Hyperperfusion peak / resolution | 2–7 days / 2–3 weeks | `TAU_ADAPT` 22 d (model 7–8 days / 9–13 days) |
| Direct bypass patency rate / graft flow | >95 % / 30–60 mL/min | `GBYP_DIR` 1.05 (model Q_byp 65 mL/min) |

The most weakly constrained parameter is **`g_leak`** (leptomeningeal coupling),
and for exactly that reason it is exposed as a slider in the Shiny app and
reported as a sweep in (8) above.

---

## 6. Running

```bash
# mechanistic map
dot -Tsvg mmd_qsp_model.dot -o mmd_qsp_model.svg
dot -Tpng -Gdpi=150 mmd_qsp_model.dot -o mmd_qsp_model.png

# numerical ground truth (numpy/scipy) — regenerates every table above
python3 mmd_reference_model.py

# mrgsolve model: verification · JAM reproduction · MAP optimisation · virtual population
Rscript mmd_mrgsolve_model.R

# Shiny dashboard (12 tabs)
Rscript -e 'shiny::runApp("mmd_shiny_app.R")'
```

---

## ⚠️ Disclaimer

This is a mechanistic model for teaching and research. It has not been validated
against patient-level data and **must not be used for clinical
decision-making.** In particular, the "blood-pressure optimum" of (6) above is a
property of the model's internal risk functions and not a recommendation for
care.
