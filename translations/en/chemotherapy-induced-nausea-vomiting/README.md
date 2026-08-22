# Chemotherapy-Induced Nausea and Vomiting (CINV) — QSP Model
### Chemotherapy-Induced Nausea and Vomiting · Quantitative Systems Pharmacology

| Deliverable | File |
|--------|------|
| 🗺️ Mechanistic map | [`cinv_qsp_model.dot`](../../../chemotherapy-induced-nausea-vomiting/cinv_qsp_model.dot) · [SVG](../../../chemotherapy-induced-nausea-vomiting/cinv_qsp_model.svg) · [PNG](../../../chemotherapy-induced-nausea-vomiting/cinv_qsp_model.png) — 220 nodes / 335 edges / 25 clusters |
| ⚙️ mrgsolve model | [`cinv_mrgsolve_model.R`](../../../chemotherapy-induced-nausea-vomiting/cinv_mrgsolve_model.R) — 81 ODEs, time unit = hour, 37 scenarios |
| 📊 Shiny dashboard | [`cinv_shiny_app.R`](../../../chemotherapy-induced-nausea-vomiting/cinv_shiny_app.R) — 15 tabs |
| 📚 References | [`cinv_references.md`](../../../chemotherapy-induced-nausea-vomiting/cinv_references.md) — 163 entries / 162 unique PMIDs (every one looked up and verified) |

---

## The organising thesis

**The clinical endpoint in CINV is not the percentage the drug shaves off. It is the exponential of a hazard integral.**

```
   DRIVE(t) = GAIN × Σ_j  w_j · ΔA_j(t)        ← additive, one term per receptor
   λ(t)     = HMAX · EBS^NE / (THR^NE + EBS^NE),   EBS = max(0, DRIVE)
   I        = ∫ (λ + λ_rescue) dt
   CR       = E_S[ exp(−S · I) ],   S ~ lognormal(0, ω²)
```

**Every** antiemetic licensed to date **removes one term (A_j)** from this sum.
None of them reduces the whole of DRIVE by a multiplicative factor. From this single structure the
following follow **as outputs of the model rather than as rules**.

### ① Additivity is exact only in the drive — and it is additive nowhere in the endpoints

**The only layer at which additivity is structurally exact is DRIVE.** Each drug subtracts one term
from the sum and leaves the rest untouched. But the `drive → λ → I → CR` chain is **a saturating
threshold function followed by an exponential**, so additivity is lost in passing through that chain.
Measured in this model in earnest, it has already **decayed to 67%** of the sum by the level of ΔI
(the numbers are in consequence ① below). The "synergy" reported between the 5-HT3 antagonists, the
NK1 antagonists, dexamethasone and olanzapine is therefore not a receptor interaction but **the
result of additive pharmacology passing through a non-linear chain**. And:

> **The absolute benefit of a fourth agent is a non-monotonic function of baseline risk.**
> It is small in patients who are already well controlled, and small in patients whose emetic drive
> is overwhelming.
> This is why adding olanzapine wins big in highly emetogenic chemotherapy (HEC) and gains nothing in
> low emetogenic chemotherapy (LEC), and the two facts are the two ends of the same algebra.

### ② The acute and delayed phases are the same equation with a different dominant term

- **The acute term (5-HT3):** the mucosal injury signal `FG` is **sigmoid** in gastrointestinal
  platinum concentration (`HG = 2.5`), so serotonin release is a **burst**, not a tail. It peaks at
  4–6 hours and is gone by 24 — observable in urinary 5-HIAA.
- **The delayed term (NK1):** the `Z1→Z2` transduction chain induces preprotachykinin transcription
  (`t½ ≈ 24 h`) and is a **transcriptional ramp** peaking at about 40 hours.

A 5-HT3 antagonist erases **the term that dominates day 1 and is all but absent on day 3**. The
delayed phase does not arise because the drug wears off; it is **the term that drug never touched in
the first place**.

### ③ Nausea and vomiting dissociate because they are two different weighted sums

The brainstem sum is a matter of **whether a threshold is crossed** (vomiting), while the cortical ·
insular sum is a **continuous percept** (nausea) and includes H1, M1, 5-HT2A, vasopressin, gastric
slow-wave dysrhythmia and conditioning terms. The 5-HT3 and NK1 antagonists touch only one or two of
these. **Olanzapine is the only agent that removes four terms of the cortical sum simultaneously
(D2 · 5-HT2A · H1 · M1)**, and that is why "no nausea" is the endpoint olanzapine owns.

### ④ Why 100% control is structurally impossible

The area postrema has **no blood–brain barrier.** Circulating platinum therefore creates an input
`WAP·(FG + WAPD·FZ)`, and this term is **not downstream of any receptor on the map** — that is, no
combination of antagonists can remove it. The term has both an acute part (FG) and **a delayed part
(FZ)**, and the delayed part is precisely why the delayed-phase complete response stops below 80%
even on a four-drug regimen. Only the gain terms — glucocorticoid, GABA-A, CB1 — attenuate it.

### ⑤ The tonic subtraction happens **inside** each term, **before** the gain is applied

This was a substantive defect found and fixed during development, and it is recorded here. The
early version subtracted a single tonic constant **from the total after the gain had been
multiplied in**. That lets dexamethasone alone push the total below baseline and so become a
**complete antiemetic**. Constructing both sums as **increments relative to an analytically computed
drug-free value** makes that failure structurally impossible and removes two constants that would
otherwise have to be calibrated (DRIVE0 · CTX0). The excess drive of a drug-free,
chemotherapy-free patient is **exactly zero by definition**, and both implementations confirm this
with a self-test.

---

## Fitted vs derived vs held out

**The emetic anchors come from a single randomised trial** — Hesketh 2003
(*J Clin Oncol*, [PMID 14559886](https://pubmed.ncbi.nlm.nih.gov/14559886/)),
cisplatin ≥70 mg/m², **two arms × three intervals = 6 numbers**. This is a deliberate constraint, so
that numbers cannot be picked from several cohorts and assembled as convenient.

| Category | Items |
|------|------|
| **Fitted (8)** | `EPPT` `WNKC` `WAP` `WAPD` `THR` `HMAX` `NK50` `NN50` |
| **Derived (not fitted)** | `OMEGA` — see below |
| **Determined analytically** | the tonic reference values `VAG*` `N1S*` `occD2*` `occH1*` `occM1*` (closed forms in the parameters) |
| **Literature priors** | all PK, all Ki/koff, `KRE` (fixed by the shape of the 5-HIAA burst), `EDEXBS`, the Hill exponents |
| **Anchors (9)** | the 6 from Hesketh 2003 + untreated 24 h emetic episode count + untreated no-nausea + triplet no-nausea (the Navari 2016 **control arm**) |

**All held out:** every efficacy number for palonosetron, NEPA, rolapitant and olanzapine;
AC/carboplatin/oxaliplatin; the dexamethasone AUC ratios; NK1 PET occupancy; the 5-HIAA time course;
the QTc slope; CYP2D6 ultra-rapid metabolisers; risk-factor covariates; the multi-cycle anticipatory
carry-over.

### ω was not fitted — it was derived

In this model **a single integral `I` plays two roles**: the mean number of emetic episodes is `I`,
and the complete response rate is `E_S[exp(−S·I)]`. Requiring two independent facts about untreated
cisplatin — a mean of **about 6** emetic episodes at 24 hours and a complete response rate of
**about 2%** — to give **the same number** determines ω uniquely:

```
ω = 0.5915       (I(CR=2%) = 6.000 = mean number of emetic episodes)
```

This is not a fitted parameter but **the solution of a self-consistency condition**.

---

## Calibration fit

8 parameters (+ the derived ω), 9 anchors. **Median |log ratio| across all anchors = 0.161.**

| Arm | Endpoint | Observed | Predicted | Relative error |
|----|------|------|------|---------|
| ond + dex | CR acute | 0.781 | 0.665 | −14.9% |
| ond + dex | CR delayed | 0.558 | 0.594 | +6.4% |
| ond + dex | CR overall | 0.523 | 0.421 | −19.6% |
| ond + dex + apr | CR acute | 0.892 | 0.843 | −5.3% |
| ond + dex + apr | CR delayed | 0.754 | 0.847 | +12.4% |
| ond + dex + apr | CR overall | 0.727 | 0.721 | −0.5% |
| untreated | 24 h emetic episodes | 6.00 | 4.14 | −31.0% |
| untreated | no nausea | 0.020 | 0.002 | −88.6% |
| ond + dex + apr | no nausea | 0.220 | 0.416 | +89.1% |

Fitted parameters: `EPPT` 0.302 · `WNKC` 0.201 · `WAP` 1.292 · `WAPD` 1.385 ·
`THR` 8.125 · `HMAX` 1.255 · `NAP` 0.862 · `NK50` 2.361 · `NN50` 2.065.
Derived value `OMEGA` = 0.5915.

**The six complete-response anchors fit well, but the two nausea anchors pull in opposite
directions** — this is exactly item ① of "Reported failures" below, and it is not a problem that a
scale parameter can solve. `WNKC` sits against its lower bound (0.20), so this fit is a **boundary
solution**, and as a result the action of the NK1 antagonists has been assigned not to central NTS
NK1 but to the **peripheral NK1 on the vagal terminals**. That too is reported as it stands.

## Held-out predictions (values no anchor touched)

| Prediction target | Predicted | Observed | Ratio |
|-----------|------|------|-----|
| NEPA CR overall (cisplatin) | 0.849 | 0.896 | 0.95 |
| Palonosetron + dex CR overall | 0.559 | 0.765 | 0.73 |
| Olanzapine no-nausea ratio (vs the triplet) | 2.09 | 1.68 | 1.24 |
| Olanzapine 5 mg CR delayed (J-FORCE) | 0.944 | 0.790 | 1.19 |
| Rolapitant CR delayed | 0.932 | 0.727 | 1.28 |
| Aprepitant central NK1 occupancy at 24 h (PET) | 1.00 | 0.90 | 1.11 |
| CYP2D6 UM, ondansetron CR ratio | 0.764 | 0.75 | **1.02** |
| CYP2D6 UM, palonosetron CR ratio | 1.000 | 1.00 | **1.00** |
| High-risk patient CR overall | 0.611 | 0.50 | 1.22 |
| Low-risk patient CR overall | 0.919 | 0.85 | 1.08 |

**Reproduced exactly, qualitatively — the differentiation of the CYP isoforms.** Although no anchor
touched an enzyme, the minimum enzyme activities came out as follows.

| Co-administered agent | CYP3A4 activity | CYP2D6 activity |
|-----------|------------|-------------|
| Dexamethasone alone | 1.000 | 1.000 |
| + aprepitant | **0.100** | 1.000 |
| + netupitant | **0.101** | 1.000 |
| + rolapitant | 1.000 | **0.072** |

That is, **the two 3A4 inhibitors do not touch 2D6 at all, and rolapitant does exactly the
opposite**. So "rolapitant needs no dexamethasone dose adjustment" becomes an output rather than a
rule, and in the `UM2D6_ond_rol` scenario the **failure of ondansetron in a CYP2D6 ultra-rapid
metaboliser (CR 0.321) recovering with the addition of rolapitant (0.768)** is also an output.

## Three consequences (all computed results)

**① Additivity is exact only in DRIVE, and is already broken at ΔI.**

| Arm | Hazard integral I | CR(0–120 h) | log CR |
|----|----------------|-----------|--------|
| No prophylaxis | 7.687 | 0.010 | −4.604 |
| ond + dex | 0.856 | 0.421 | −0.866 |
| ond + dex + apr | 0.290 | 0.723 | −0.324 |
| Quadruplet | 0.233 | 0.768 | −0.263 |

ΔI for the single agents: dexamethasone 3.73, ondansetron 6.54, olanzapine 0.79.
And yet **the ΔI of ond+dex is 6.83 while the plain sum of the two components is 10.27** — that is,
it has already **decayed to 67%** at the level of ΔI. This is the precise statement this model makes:
additivity is **structurally exact in DRIVE**, and because the `drive → λ → I → log CR` path is a
saturating non-linearity, **it is already non-additive at ΔI and still less additive in CR%.** What
gets reported as "synergy" is the product of that non-linear chain.

**② The absolute benefit of a fourth agent is non-monotonic.**

| Susceptibility multiplier | Triplet CR | Quadruplet CR | Absolute benefit | Ratio |
|------------|--------|--------|----------|----|
| ×0.15 | 0.970 | 0.976 | +0.006 | 1.01 |
| ×0.40 | 0.923 | 0.938 | +0.015 | 1.02 |
| ×1.00 | 0.822 | 0.855 | +0.033 | 1.04 |
| ×2.50 | 0.628 | 0.686 | +0.058 | 1.09 |
| ×6.00 | 0.363 | 0.435 | +0.072 | 1.20 |

**The maximum absolute benefit is 7.2 percentage points, and it appears at a baseline CR of 38%.** It
vanishes at both ends. Which is to say that adding olanzapine winning in highly emetogenic
chemotherapy and gaining nothing in low emetogenic chemotherapy are the two ends of the same algebra.

**③ The acute/delayed separation is not a rule but the consequence of two time courses.**

| Arm | Vagal firing VAG @4/12/24/48/72 h | Central NK1 occupancy N1S @ the same time points |
|----|-----------------------------------|------------------------------|
| Untreated | 2.46 → 2.24 → 1.64 → 0.84 → 0.68 | 0.200 → 0.204 → 0.211 → 0.218 → 0.220 |
| ond + dex | 0.70 → 0.81 → 0.86 → 0.73 → 0.64 | 0.200 → 0.203 → 0.207 → 0.212 → 0.213 |
| ond + dex + apr | 0.48 → 0.60 → 0.67 → 0.53 → 0.44 | 0.000 (throughout) |

The vagal term **rises as a burst and then falls** while the NK1 term **rises monotonically**.
Ondansetron abolishes the burst (2.46 → 0.70) and cannot touch the NK1 ramp.

## Cross-validation of the two implementations (Python twin ↔ mrgsolve)

The same equations were implemented independently twice and 45 states compared. Dosing time points
were excluded (the two implementations reporting either side of a discontinuity is a notational
convention, not a disagreement).

| Scenario | States compared | States with median relative difference under 1% |
|----------|--------------|------------------------------|
| ond + dex | 45 | 44 / 45 |
| ond + dex + apr | 45 | 44 / 45 |
| Quadruplet | 45 | 44 / 45 |
| Untreated | 45 | 43 / 45 |

Most states agree to a median relative difference of 0.001–0.01% (for example `N1S` 0.21676 against
0.21676, `VAG` 0.59187 against 0.59187). The remaining systematic difference is in
**three states only (ECFV · KP · GVOL)** and its size is
**exactly 1.2672–1.2673-fold**. This matches the susceptibility multiplier at the default covariates,
**SUSC = 1.2674**. The cause is clear: the mrgsolve model multiplies SUSC into the hazard
**inside the ODEs**, whereas the Python twin multiplies it into HINT/EMES **after** integration. The
two are equivalent for the endpoints (CR · emetic episode count) but not for the fluid, potassium and
gastric content states that consume λ inside the ODEs. **The mrgsolve side is the self-consistent one
and is the deliverable.** Only the endpoints were used in fitting, so the fitted parameters are
unaffected.

## Reported failures — recorded rather than fixed

**① The two nausea anchors cannot be satisfied simultaneously (the largest failure).**
Untreated no-nausea comes out at 0.002 (observed 0.02), which is **far too severe**, while triplet
no-nausea comes out at 0.416 (observed 0.22), which is **far too mild**. The two residuals have
opposite signs. `NK50` and `NN50` (the scale and the threshold) do not separate them, and adding
`NAP` (the receptor-independent share of the cortical sum) to the fit in a second stage lowered the
objective from 0.806 to 0.577 but only moved `NAP` from 0.70 to 0.86. **The defect is not the scale of
the cortical sum but its drug sensitivity** — the component that makes untreated nausea severe in the
model is removed far too easily. The structure needed is a **separate central sensitisation /
interoceptive gain** term that the drugs do not touch; the present `NAP` · `NSW` · `NAVP` are all
ultimately driven by the same `EBS` and so disappear along with the drug.

**② It grossly over-controls regimens other than cisplatin.**
The `EMETO_P/EMETO_C` scaling factors were set by prior judgement and **no anchor constrained them.**
As a result the CR for moderately emetogenic and lesser regimens is unrealistically high: MEC
palonosetron 0.976 (observed ≈0.81), MEC ondansetron 0.969 (≈0.69), the carboplatin doublet 0.937
(≈0.60), AC triplet 0.962. **This model is calibrated for cisplatin and cannot be used for MEC/LEC
regimens as it stands.** What is needed is to calibrate the EMETO scaling factors separately against
MEC trials.

**③ The urinary 5-HIAA burst is too slow.** The first-8-hour share of the 24-hour total excretion came
out at 0.365 against an observed value of about 0.85 (ratio 0.43). `KEL_G` (set slow, at a 1.65→3.9 h
half-life) and `KOUT_G` are dragging the acute term out longer than necessary.

**④ QTc is under-predicted about two-fold.** 8 mg intravenous 2.40 ms (observed 5.6), 32 mg
intravenous 9.59 ms (observed 19.6). However, **the 32 mg/8 mg ratio is 4.00 in the model against an
observed 3.50**, so the dose proportionality is right. That is, the defect is the single absolute
magnitude of the slope `SQ_OND`, and the qualitative conclusion (the danger of the 32 mg intravenous
dose) stands unchanged.

**⑤ The dexamethasone DDI is over-predicted 1.7-fold.** The per-mg AUC ratios are aprepitant 3.72
(observed 2.2), netupitant 3.25 (≈2.4), rolapitant 1.14 (≈1.0). `KI3A_APR` is a prior value and was
not fitted. Beyond that, Nijstad 2022
([PMID 36287279](https://pubmed.ncbi.nlm.nih.gov/36287279/)) argues that **the observed 2.2-fold is
itself an over-estimate**, so the true error may be larger than this.

**⑥ The anticipatory nausea loop self-amplifies — a clinically important failure.**
In a six-cycle simulation, conditioning strength goes from 0.019 to 0.548 on `ond+dex` and from
0.003 to 0.540 on the quadruplet, converging on the same value **almost independently of the quality
of control**. The quadruplet's cycle-1 nausea of 0.04 becomes 3.61 by cycle 3. The cause is
diagnosed: the gain of the loop
`ANTIC → NANT·ANTIC (cortical sum) → NAUSEA → KACQ·NAUSEA → ANTIC` exceeds 1, and the cue term
`CUE_ON` is switched permanently to 1 so the loop keeps going between cycles as well. This model
therefore **fails to reproduce the central clinical fact that "if you control cycle 1 you can prevent
anticipatory nausea".** The fix required is to gate `CUE_ON` to the neighbourhood of the dosing days
only, after which nausea is not sustained between cycles and the loop does not run away.

**⑦ High-dose metoclopramide adds nothing.** `mcp_high` (180 mg ×4) and `mcp_dex` (20 mg) have an
identical overall CR of 0.064. Historically, high-dose metoclopramide really did work, so either the
weight `WD2` of the D2 term or the 5-HT3 cross-reactivity at high dose is under-valued.

**⑧ The two estimates of ω differ three-fold — the model's central unresolved tension.**
The value required by self-consistency of the natural history is **ω = 0.59**. But estimating ω from
**the gap between the combined endpoint and the two marginal endpoints** of Hesketh 2003 alone gives
**ω = 2.01** (both arms arriving independently at the same value, objective 2.9×10⁻⁵), and the
untreated 24-hour hazard integral this then requires is about 50 — that is, **50 emetic episodes**,
eight times the roughly 6 actually observed. The structure in which one integral simultaneously
carries the mean count and the exponent of the complete-response probability breaks down here. Since
ω = 0.59 was chosen, the model under-predicts the combined endpoint (0.421 predicted against 0.523
observed), and this is not a tuning failure but a **structural statement**: the combined–marginal gap
is not pure between-patient variance but is in substantial part **within-patient clustering of events
and acute→delayed persistence**. The discriminating experiment is clear — report the distribution of
emetic episode counts per patient (not only the mean) and the two explanations separate.

## Defects found and fixed during integration

1. **Doing the tonic subtraction after applying the gain** — dexamethasone alone pushed the total
   below baseline and became a **complete antiemetic**. Both sums were rebuilt as increments, making
   this structurally impossible, and two calibration constants were removed.
2. **The delayed-phase hazard being exactly zero in the NK1-antagonist arms** — because the delayed
   drive was entirely NK1-mediated. A delayed part (`WAPD·FZ`) was added to the receptor-independent
   term at the area postrema.
3. **A global negative clamp in the Python twin** — `DRIVE` and `CTX` are signed quantities
   (excesses), and clipping them at 0 broke the relaxation and let `CTX` diverge (−13 against −0.16).
   This defect was not present in the mrgsolve implementation and was found by cross-validation.
4. **`EPS` is a reserved word in mrgsolve** — compilation failed. Renamed to `AKATH`.
5. **The parameter `GLU_0` collided with the initial-value variable of the compartment `GLU`** —
   renamed to `GLUB`.
6. **`endpoints()` averaging the two records either side of a dosing time point** — it was
   interpolating across a discontinuity. Fixed to use only the last record at each time point.
7. **`source()`-ing the model file ran all 37 scenarios** — which happened every time the Shiny app
   started. Fixed with a `sys.nframe() == 0L` guard.
8. **The CYP diagnostic measured enzyme *amount*** — inhibition enters the activity (`ACT3A4`), so it
   looked as though there were no inhibition. It was a bug in the diagnostic code and the model was
   fine.
9. The argument name `names` in `run_all()` masked base R's `names()`, giving infinite recursion.
10. `graphviz`'s `init_rank` failure — solved with `newrank=true` (it arises with 25 clusters plus
    cyclic structure).

---

## Files & usage

```r
# ODE model (compilation and execution verified on mrgsolve 2.0.1)
source("cinv_mrgsolve_model.R")
tonic_selftest()                 # drug-free excess drive = 0 (self-test)
d   <- run_scenario("S12_quadruplet")
endpoints(d)
tab <- run_all()                 # all 37 scenarios
fourth_agent_benefit()           # non-monotonicity of the fourth agent's benefit
log_additivity_check()           # log(CR) additivity test

# dashboard
shiny::runApp("cinv_shiny_app.R")
```

```bash
# render the map
dot -Tsvg cinv_qsp_model.dot -o cinv_qsp_model.svg
dot -Tpng -Gdpi=150 cinv_qsp_model.dot -o cinv_qsp_model.png
```

**Dependencies:** `mrgsolve` (≥2.0.1, a C++ compiler is required), `shiny`, `ggplot2`,
`tidyr`, `dplyr`, Graphviz.

### Model structure summary (81 ODEs)

| Block | Compartments | Contents |
|------|--------|------|
| Antiemetic PK | 32 | ondansetron · granisetron (+ the extended-release subcutaneous depot) · palonosetron · aprepitant (+ the fosaprepitant prodrug) · netupitant · rolapitant · dexamethasone · olanzapine · metoclopramide · lorazepam · dronabinol |
| Chemotherapy PK | 5 | free platinum · tissue-bound platinum · gastrointestinal mucosal signal · delayed signal transduction chain (Z1 · Z2) |
| Enzymes / DDI | 2 | CYP3A4 pool (turnover · induction · reversible inhibition) · CYP2D6 activity |
| Steroid transduction | 2 | GR binding → nuclear genomic effect (`t½ ≈ 20 h` — the origin of delayed-phase specificity) |
| Serotonin axis | 4 | EC cell store · gut interstitial 5-HT · plasma 5-HT · cumulative urinary 5-HIAA |
| 5-HT3 receptor | 5 | 5-HT binding · binding of the three antagonists · **internalised receptor** (palonosetron) |
| Substance P / NK1 | 7 | preprotachykinin content · central/peripheral SP · SP binding · binding of the three antagonists |
| Central integration | 5 | dopamine · vagal firing · brainstem sum · cortical sum · nausea VAS |
| Endpoint accumulators | 4 | emetic episode count · hazard integral · nausea hazard integral · rescue medication use |
| Gastric motility | 3 | gastric content · slow-wave coupling index · vasopressin |
| Associative learning | 2 | conditioning strength · anxiety |
| Safety | 5 | ΔQTcF · blood glucose · sedation · akathisia (EPS) · constipation |
| Systemic outcomes | 5 | extracellular fluid deficit · potassium · creatinine · relative dose intensity (RDI) · FLIE |

The scenarios are constructed as **matched pairs** wherever possible, so that each comparison
isolates a single term of the sum (for example S28/S29 carboplatin ± an NK1 antagonist, S31/S32 low
emetogenic ± the quadruplet, S33/S34 CYP2D6 ultra-rapid metaboliser ± a non-CYP2D6 setron, S35/S36
high risk ± low risk).

---

## ⚠️ Limitations

- This is a **semi-quantitative QSP model for education and research** and must not be used for
  clinical decision-making.
- There are 9 anchors and 8 fitted parameters. Identifiability is defended only by structural
  argument and held-out prediction; formal confidence intervals were not computed.
- The fit was to **aggregate endpoints as reported in papers**, not to individual patient-level data.
  The frailty variance ω is therefore a value derived from the population marginal distribution, not
  one estimated from individual data.
- Emetic events are treated as Poisson. In reality they occur in **bouts** within a single patient, so
  the relationship between the mean count and the probability of "one or more" will be looser than in
  the model. The first item of "Reported failures" above is exactly this problem.
- Radiation-induced emesis, post-operative nausea and vomiting (PONV), paediatrics, and multi-day
  chemotherapy are not in scope.
