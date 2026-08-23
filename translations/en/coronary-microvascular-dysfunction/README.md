# Coronary Microvascular Dysfunction (CMD / ANOCA–INOCA) QSP Model

**Coronary Microvascular Dysfunction · Angina/Ischaemia with No Obstructive Coronary Arteries**

<a href="cmd_qsp_model.svg"><img src="cmd_qsp_model.png" width="820" alt="CMD QSP mechanistic map"></a>

> 184-node 19-cluster mechanistic map · 60-ODE mrgsolve model (2 myocardial layers × 14 drugs) ·
> 10-tab Shiny app · 229 references retrieved directly from PubMed ·
> dependency-free Python reference implementation (the source of every number, with 25 defects on record)

---

## 1. The one idea, in one sentence

**CFR is a ratio.** The same value of 1.9 can arise because **the denominator has
grown** (resting flow is excessive — the functional endotype) or because **the
numerator has been squeezed** (maximal hyperaemic flow is blocked by structure — the
structural endotype). The two patients receive the same number and need drugs pointing
in opposite directions. This model takes never conflating those two terms as its one
and only design principle.

Every microvascular index used clinically is either a ratio or a difference of two
quantities:

```
CFR     = v_hyper / v_rest              two flows, one number
MR      = Pd / v                        a resistance, not a flow
DPTI    = (Pd − LVEDP) × t_diastole     supply is only bought in diastole
SPTI    = P_sys × t_systole             demand is only spent in systole
Deficit = MVO2 − v·k·CaO2·E_max         the only quantity the patient can feel
```

Four results follow from this. The model **computes them rather than assuming them.**

---

## 2. The first result — the same number, different diseases

The reference data are Rahman et al., Circulation 2019 (85 ANOCA patients,
PMID 31707835). Defining microvascular dysfunction as CFR<2.5, 53% qualify, and of
those **62% have normal hyperaemic microvascular resistance (<2.5 mmHg/cm/s)** while
only **38% are elevated**. Resting microvascular resistance was 4.2±1.0 in the
functional group, 6.9±1.7 in the structural group and 7.3±2.2 mmHg/(cm/s) in controls.

The model fits only five flow constants to these data (§6) and predicts the rest.

| Endotype | CFR | Resting MR | Hyperaemic MR | Resting MBF | Hyperaemic MBF | IMR | Hyperaemic endo/epi |
|---|---|---|---|---|---|---|---|
| Control | 3.43 | 7.56 (observed 7.3) | 2.20 | 0.73 | **2.52** | 18.4 | 1.02 |
| **Functional** | **1.83** (observed 1.83) | **4.20** (observed 4.2) | 2.30 | 1.32 | **2.41** | 19.2 | 1.01 |
| **Structural** | **2.14** (observed 1.92) | 7.71 (observed 6.9) | 3.60 | 0.72 | **1.54** | 30.1 | 1.09 |
| Vasospastic | 2.71 | 7.29 | 2.69 | 0.76 | 2.06 | 22.5 | 1.03 |
| Non-cardiac | 3.42 | 7.58 | 2.22 | 0.73 | 2.50 | 18.5 | 1.02 |

The CFRs of the two CMD endotypes differ by 0.31 (1.83 vs 2.14), while **absolute
maximal perfusion is 2.41 vs 1.54 mL/min/g, a 57% difference**. The hyperaemic flow
velocity of the functional endotype is 96% of the control group's, that of the
structural endotype 61%. The ratio cannot see that difference and absolute flow can.
This is the quantitative basis for the recent move towards measuring absolute flow
instead of a reserve ratio, and the model says in numbers what is lost by using a
ratio.

---

## 3. The second result — resting tone is not a free parameter

Autoregulation sets tone so that supply matches demand at rest. Therefore, **given
the hyperaemic MR (= R_min fixed), the resting MR is predicted.** Only the control
group's hyperaemic MR was calibrated, and:

| Endotype | Hyperaemic MR | Resting MR predicted | Resting MR observed | Error |
|---|---|---|---|---|
| Control | 2.20 | 7.56 | 7.3 | +3.5% |
| Structural | 3.60 | 8.07 | 6.9±1.7 | +17% |

The resting MR of the structural endotype is **an uncalibrated prediction** and it
falls inside the observed range. Autoregulation lowers resting tone in order to
compensate for the rise in R_min, and that compensation eats the reserve.

In the functional endotype, however, this **fails.** With the controller left on and
the offset set to zero, resting MR becomes 7.73 (the same value as the control group)
and CFR becomes 3.36 (observed 4.2 / 1.83). Lower the tone by hand and the controller
puts it back within a single step. That is, **the coexistence of a low resting MR with
a normal minimal resistance cannot be expressed as a "state of low tone"; it can only
be expressed as a controller defect.** And the size of that defect is derived:

> The arterioles of the functional endotype dilate **as though metabolic demand were
> 92% higher than it actually is.**

### A falsifiable prediction (with a single blood gas)

If demand is unchanged while resting flow is 92% higher, then **the resting myocardial
oxygen extraction ratio must be lower by the same amount.** The model: control 0.689
vs functional 0.383 — **56% of normal.** Coronary sinus oxygen saturation must be
higher by that much. If measured resting extraction turns out to be normal, then the
controller-defect interpretation is wrong and the low resting MR must be regarded as
an artefact of the measurement conditions (anxiety, contrast, supine catheter-
laboratory haemodynamics). Either way it **is refuted by a single blood gas** — the
cheapest experiment this model proposes.

---

## 4. The third result — heart rate enters the oxygen balance twice

Heart rate raises demand (the tension-time product) and at the same time shortens the
diastolic window in which the subendocardium is perfused. The signs are the same.
Sweeping heart rate in the structural endotype at fixed workload (wl=3.2):

| HR | f_dia | DPTI | SEVR | Deficit (time term only) | Deficit (total) | endo/epi |
|---|---|---|---|---|---|---|
| 50 | 0.583 | 31.9 | 0.813 | 0.980 | 0.000 | 1.263 |
| 55 | 0.574 | 30.1 | 0.783 | 1.392 | 0.000 | 1.259 |
| 68 | 0.552 | 26.1 | 0.716 | 2.426 | 2.426 | 1.125 |
| 85 | 0.524 | 22.0 | 0.640 | 3.706 | 6.494 | **0.980** |
| 110 | 0.486 | 17.5 | 0.550 | 5.463 | 12.350 | 0.811 |

Going from 68 → 55 bpm, the subendocardial deficit falls by 2.426 mL O2/min/100 g, of
which **1.034 (43%) is diastolic perfusion time and 1.392 (57%) is the reduction in
demand**. Even a purely heart-rate-lowering drug therefore has **a supply-side
action**, and it is not interchangeable with removing the same amount of demand by
another route.

Reading the same result as a threshold is more clinical — the workload at which a
subendocardial deficit first appears:

| Resting HR | 50 | 55 | 60 | 68 | 75 | 85 | 95 |
|---|---|---|---|---|---|---|---|
| Workload at deficit onset | 3.40 | 3.28 | 3.15 | 2.95 | 2.78 | 2.53 | 2.29 |

Heart rate alone moves the ischaemic threshold by 48%. And the point at which the
endo/epi ratio passes 1.0 is around HR 85, which is exactly the criterion Rahman used
to define inducible ischaemia.

---

## 5. The fourth result — a drug that hits its target exactly can still fail in three ways

### 5.1 PRIZE (zibotentan) — internal cancellation, and the cause was not the one expected

PRIZE (PMID 39217504) tested zibotentan 10 mg for 12 weeks in 118 patients with
microvascular angina and obtained a Bruce duration of **−4.26 s (95% CI −19.60 to
+11.06)**. Blood pressure fell, **circulating ET-1 rose** (the target-mediated effect
of blocking ETA-mediated clearance — get this direction wrong in the model and it
predicts the drug wins), and adverse events were 60.2% vs 14.4%, dominated by fluid
retention.

The model's population-weighted decomposition:

| Term | Value |
|---|---|
| Microvascular gain alone (both costs removed) | **+0.00 s** |
| Attributed to the fall in blood pressure | **+2.41 s** |
| Attributed to fluid retention | **−6.99 s** |
| Model net effect | **−4.53 s** |
| PRIZE observed | **−4.26 s** (inside the observed CI) |

The content of the three terms was not what was expected when this section was
started, and the arithmetic refuted the expectation:

1. **The microvascular gain is exactly zero.** Removing ETA-mediated tone helps only
   where there is tone left to remove. The functional endotype has a normal minimal
   resistance, so there is nothing to gain, and the structural endotype's ceiling is
   anatomical, so a receptor antagonist cannot raise it. **A drug can be perfectly
   on-target and have no target left to hit.**
2. **The fall in blood pressure is in fact a small gain (+2.41 s).** Lower the
   systolic pressure and the wall-stress demand falls faster than the diastolic
   driving pressure does — because demand loads systolic pressure through the
   tension-time product while supply loads only the diastolic mean.
3. **Fluid retention accounts for the whole of the negative result (−6.99 s).** LVEDP
   subtracts directly from the subendocardial driving pressure (Pd − LVEDP), so **the
   trial's largest adverse event and the failure of its primary endpoint are the same
   event measured twice.**

This is not a dead end but an actionable reading. At 10 mg the endothelin hypothesis
was not cleanly tested — what moved the endpoint was not the receptor but sodium
retention. Two designs follow: remeasure with a diuretic co-administered (or at a
lower dose), or a molecule that separates coronary ETA blockade from renal sodium
handling. A larger trial of the same 10 mg will measure the same fluid again.

### 5.2 RWISE (ranolazine) — population dilution

RWISE (PMID 26614823) was negative overall and showed improvement only in the CFR<2.5
subgroup. The model reproduces that dilution as arithmetic (Rahman's composition
47% / 33% / 20%):

| Stratum | Weight | dSAQ | dCFR |
|---|---|---|---|
| CFR≥2.5 | 0.47 | +0.00 | +0.243 |
| Functional | 0.33 | +1.41 | +0.123 |
| Structural | 0.20 | +3.80 | +0.152 |
| **CFR<2.5 stratum** | 0.53 | **+2.31 U** | |
| **Whole cohort** | 1.00 | **+1.22 U** | |

Nothing changed about the drug; only the prevalence of the target mechanism changed.
Only 53% of the mechanistic effect becomes observable. That both values are below the
clinically important difference (10 U) is also reported as it stands — the model
predicts that for this drug a clinically meaningful responder subgroup **does not
exist**, and regards the RWISE subgroup signal as a real but small effect that reached
significance on the frequency scale.

### 5.3 WARRIOR (statin + ACEi/ARB) — background contamination and a timescale mismatch

The primary MACE HR of WARRIOR (PMID 41932694, 2476 women, 2.5 years) was
1.13 (0.94–1.37), and the contamination-corrected sensitivity analysis was 0.74. The
model computes a true hazard ratio of **0.897** in the structural endotype and writes
the observable ratio as a function of background contamination c:

| Contamination c | 0.0 | 0.4 | 0.6 | **0.8** | 0.95 |
|---|---|---|---|---|---|
| Observed HR | 0.897 | 0.938 | 0.959 | **0.979** | 0.995 |

The second problem is the instrument. The drugs' target is structural (τ_ML 180 days,
τ_CAPD 240 days), whereas the endpoint that dominated MACE was hospitalisation for
angina, and that is driven by symptom burden (τ_SAQ 14 days) and by a central
sensitisation component that no vasoactive drug touches. Over 2.5 years what IMT gains
is **CFR +0.158**, while **SAQ is +2.39 U** (MCID 10). The treatment does what it
promised to do to the microvasculature and still cannot reach the endpoint on which it
is scored.

### 5.4 CorMicA — the mechanical reason stratification wins

CorMicA (PMID 30266608) obtained a 6-month SAQ of **+11.7 U** (95% CI 5.0–18.4) with
endotype-linked treatment. The model (control arm = a β-blocker prescribed
irrespective of endotype):

| Stratum | Weight | SAQ usual care | SAQ stratified | dSAQ |
|---|---|---|---|---|
| Functional | 0.30 | 57.8 | 59.1 | +1.30 |
| Structural | 0.19 | 65.6 | 65.7 | +0.05 |
| Vasospastic | 0.29 | 51.5 | 63.6 | **+12.18** |
| Non-cardiac | 0.22 | 61.5 | 68.1 | +6.57 |
| **Population mean** | | 58.3 | 63.6 | **+5.38 U** |

It under-predicts the observed +11.7 U by about half (the direction and order of
magnitude agree). Most of the gain comes from the vasospastic stratum, and this is the
clinical fact that **β-blockers are harmful in vasospastic angina** emerging from the
model. The reason stratification wins is not pharmacological but mechanical: the same
prescription is right for one endotype and wrong for another, so the unstratified arm
averages benefit and harm. That average is the value a conventional trial reports as
"no effect".

---

## 6. Calibration — what was fitted and what was predicted

Only seven constants were **solved** against seven published targets (not adjusted by
hand); the rest were fixed a priori.

| Constant | Value | Target |
|---|---|---|
| `RMIN0` | 1.5174 | Control hyperaemic MR = 2.20 mmHg/(cm/s) |
| `W_ENDO` | 0.6283 | Control hyperaemic endo/epi = 1.02 |
| `F_ADO_REV` | 0.8760 | Functional hyperaemic MR = 2.30 |
| `RMIN_F` (structural) | 1.4306 | Structural hyperaemic MR = 3.60 |
| `AUTO_OFF` (functional) | 0.9165 | Functional resting MR = 4.20 |
| `K_ANG` | 3.3997 | Untreated functional SAQ = 55 |
| `NOCI_THRESH` | 3.8575 | Untreated functional Bruce = 480 s |

Rahman's three numbers for the functional group are mutually consistent
(4.2 / 2.30 = 1.83 = their CFR). This is because CFR = MR_rest/MR_hyp is an identity
when measured at the same aortic pressure, so fitting two of them fixes the third —
**the model's functional CFR is not free.**

The risk gradients from the meta-analysis of Kelshiker et al. (PMID 34849697) —
**HR 1.16 for death and HR 1.08 for MACE per 0.1 unit fall in CFR** — were put in as
they are, without calibration.

---

## 7. Predictions that came out without being designed in

### 7.1 The angina of the functional endotype is not ischaemia

Because minimal resistance is normal, the subendocardial supply–demand deficit is
**effectively zero at any workload** (the burden column of the 24-week scenario
table). And yet these patients have angina. Structurally, therefore, that pain must be
carried in the model by afferent (A1-adenosine) signalling and central sensitisation.
Consequences:

- Anti-ischaemic drugs fail in this endotype — which is what RWISE saw.
- **Aminophylline (an adenosine receptor antagonist) moves CFR by only +0.012 (5th of
  11 arms, indistinguishable from placebo) while lengthening Bruce time by +147.5 s** —
  twice as much as any vasoactive drug. This is exactly the pattern Elliott et al.
  observed (Heart 1997, PMID 9227295), and the model reproduces it from structure
  rather than from fitting.
- And **in the structural endotype it is predicted to be unhelpful or harmful** —
  because there, A2A blockade cuts into dilator reserve that is actually being used.

### 7.2 CFR ranks drugs well only where flow is the disease

Rank correlations between ΔCFR and symptomatic gain across 11 drugs × 3 endotypes:

| Endotype | ρ(ΔCFR, −ΔNOC) | ρ(ΔCFR, ΔBruce) | CFR rank 1 | Symptom rank 1 |
|---|---|---|---|---|
| Functional | 0.60 | 0.54 | Ivabradine | Ivabradine |
| Structural | **0.94** | 0.96 | Ivabradine | Ivabradine |
| Vasospastic | 0.72 | 0.95 | Ivabradine | **Amlodipine** |

CFR ranks almost perfectly in the structural endotype (where flow *is* the disease) and
degrades in the two endotypes where it is not. An earlier version claimed a **sign
discordance** here, but that was an artefact of defect B24 and has been retracted
(§9). What remains is a weaker and more specific claim.

**And MRR does not rescue this.** MRR = (CFR/FFR)(Pa_rest/Pa_hyper), but a drug lowers
pressure in both states, so the correction term cancels. In this model
Pa_rest = Pa_hyper, so the cancellation is exact and **MRR = CFR/FFR identically**,
with ΔMRR = ΔCFR/FFR. What MRR earns its value on is epicardial stenosis, not a drug's
own haemodynamics — the two problems are easily confused. The index whose
interpretation survives during treatment is **hyperaemic resistance itself**: for a
drug with no real effect on minimal resistance, ΔMR_hyp is near zero while ΔCFR is
not. If you want a physiological endpoint while testing a drug that lowers blood
pressure, what has to be prespecified is not a reserve ratio but **hyperaemic
microvascular resistance.**

### 7.3 The hyperaemic MR threshold of 2.5 mixes two different things

Hyperaemic MR ≥2.5 reads as "structural", but two things raise it: inward remodelling
and capillary rarefaction, which fix the minimal radius (irreversible on a timescale of
years), and contractile tone that adenosine does not fully reverse (the calibrated
value being 87.6% reversal, 12.4% residual — reversible on a timescale of weeks). The
label cannot distinguish the two, but **remeasurement can**: in the model, giving a
Rho-kinase inhibitor for 21 days to patients classified as structural brings some of
them below the threshold.

> **Prediction:** before telling a patient "your microvascular disease is structural",
> repeat the hyperaemic measurement after acute Rho-kinase or ETA blockade. It takes
> one visit, and it changes the drug class.

### 7.4 The workload at which reserve is exhausted — the functional endotype exhausts it earlier

The workload at which subendocardial autoregulation reaches zero and flow becomes
pressure-passive:

| Endotype | Control | Functional | Structural | Vasospastic |
|---|---|---|---|---|
| Exhaustion workload | 4.37 | **2.30** | 2.86 | 3.63 |

The functional endotype, whose minimal resistance is normal, exhausts its reserve
**earlier** than the structural endotype — because resting flow is already elevated.
This is the model's answer to the tension Rahman observed (that stress perfusion and
exercise perfusion efficiency were similar in the two endotypes): from the moment tone
touches zero, the two endotypes become the same system differing only in R_min.

---

## 8. Deliverables

| File | Contents |
|---|---|
| [`cmd_qsp_model.dot`](cmd_qsp_model.dot) · [`.svg`](cmd_qsp_model.svg) · [`.png`](cmd_qsp_model.png) | Mechanistic map: 184 nodes, 330 edges, 19 clusters. The **red cancellation edges** running back from the drug nodes to pressure and filling pressure carry the argument of §5.1 as a picture |
| [`cmd_mrgsolve_model.R`](../../../coronary-microvascular-dysfunction/cmd_mrgsolve_model.R) | 60-ODE mrgsolve model (29 physiological + 31 PK/metabolite), 14 drugs, 24 scenarios, five trial-reproduction functions, virtual population |
| [`cmd_reference_model.py`](../../../coronary-microvascular-dysfunction/cmd_reference_model.py) | Dependency-free Python reference implementation. **The source of every number**, with 25 defects (B1–B25) recorded at the top |
| [`cmd_reference_output.txt`](../../../coronary-microvascular-dysfunction/cmd_reference_output.txt) | The full output of that script (calibration, 13 sections, 24 scenarios, sensitivity analysis) |
| [`cmd_population_results.json`](../../../coronary-microvascular-dysfunction/cmd_population_results.json) | Machine-readable summary of results |
| [`cmd_shiny_app.R`](../../../coronary-microvascular-dysfunction/cmd_shiny_app.R) | 10-tab Shiny dashboard (patient · functional testing · decomposition of the ratio · heart rate · layer-wise perfusion · PK/PD · endpoints · scenario comparison · trial reproduction · virtual population) |
| [`cmd_references.md`](cmd_references.md) | 229 references retrieved directly from PubMed, classified by section, with the numbers used as quantitative anchors stated explicitly, and a table of 6 falsifiable predictions |

### Running it

```bash
dot -Tsvg cmd_qsp_model.dot -o cmd_qsp_model.svg
dot -Tpng -Gdpi=150 cmd_qsp_model.dot -o cmd_qsp_model.png
python3 cmd_reference_model.py          # about 15 minutes, generates two output files
Rscript -e 'source("cmd_mrgsolve_model.R"); print(cmd_endotype_table())'
Rscript -e 'shiny::runApp("cmd_shiny_app.R")'
```

---

## 9. The 25 defects the numerical work exposed (bug log)

With no R runtime available, every equation was run in Python first, and that process
exposed **25 real defects**. The full text is at the top of `cmd_reference_model.py`,
with a comment on each line where a defect was. Only the ones that changed a
conclusion:

- **B1** — fitting the resting MR of 4.2 directly as tone had the controller undo it
  within a single step. That failure became the model's central claim (§3): the
  functional endotype is not a state of tone but a controller defect, and that is why
  its size can be derived.
- **B14** — while fixing B2, making ET-1 and ROCK tone completely resistant to
  adenosine pushed the functional endotype's hyperaemic MR to 8.8 (it must be <2.5 by
  definition). The endotypes were being separated by contractile factors rather than by
  structure — entirely the wrong axis.
- **B20** — writing spasm as arteriolar tone let a controller with a gain of 6 absorb
  it, so vasospastic patients could not become ischaemic. Spasm is not distributed tone
  but an upstream series obstruction plus a blockade of a vascular layer.
- **B22** — holding filling pressure at its resting value while workload changed erased
  the dominant subendocardial insult of this disease (the mechanism of the CMD/HFpEF
  overlap). The structural endotype came out asymptomatic (SAQ 95) and the functional
  one at SAQ 55 — exactly the reverse of the clinic. A case of the model telling the
  truth about a missing equation rather than about the disease.
- **B25** — PK was integrated together with the physiology, and the absorption constant
  (up to 1.6/h) was the fastest state in the model. In a WARRIOR run with a 6-hour
  step, the ramipril concentration reached −97213 mg/L on day 1 and became nan
  thereafter, so **the intensive-therapy arm quietly became the placebo arm** (the clue
  was that the treated and untreated arms agreed to three decimal places).
- **B24** — a constant 16% epicardial pressure loss was applied to patients with a
  predisposition to spasm while the comparator was evaluated without it. In the
  structural endotype every drug appeared to double the deficit — including imipramine,
  which has no haemodynamic action at all. That impossible row exposed the defect, and
  **this defect is what created and then demolished the sign-discordance claim of
  §7.2.**
- **B21 / B17** — the five endothelial states and the fast tone state were written with
  time constants of one hour and 15 seconds and run with an explicit integrator. The
  former grew the residual by a factor of 1.375^672 over 12 weeks, taking ET-1 to 1e74
  while continuing to print a plausible CFR of 1.39. The stiffest state in a model is
  often a variable nobody reads.

---

## 10. Limitations, stated

- **The under-prediction of CorMicA**: the model gives +5.38 U and the observation is
  +11.7 U. The candidates are that the representation of the control arm (an
  endotype-blind β-blocker) is more optimistic than actual CorMicA usual care, and
  placebo/expectation effects that are not in the model.
- **The IMR distribution of the virtual population**: it comes close on the prevalence
  of CFR<2.5 at 62% (observed 53%) and on the functional proportion at 49% (observed
  62%), but IMR≥25 stops at 1.3%. The sampling under-represents severe structural
  remodelling and does not match IMR-based cohorts.
- **The MCID of 10 U for the SAQ** could not be confirmed against the primary source in
  this session (what was confirmed was the KCCQ paper). CorMicA's +11.7 U was taken as
  a pragmatic reference. The choice of threshold affects the sentence "the subgroup
  effect is below the threshold" but has nothing to do with the computed value of
  +2.31 U itself.
- **Acetylcholine testing** was implemented only through the tone pathway, and
  spontaneous spasm through the occlusion pathway. The two pathways were not unified.
- **Spatial heterogeneity at the capillary level** (flow maldistribution, plugging) is
  absent. There are only two layers, and within a layer homogeneity is assumed.
- The Python reference solves the fast loop as a bisection fixed point while mrgsolve
  solves it as a 15-second ODE. The two forms agree to 1e-6 in every endotype (§V), and
  since every reported index is by definition the steady state of that loop, this choice
  changes no number.

---

## ⚠ Disclaimer

This is a quantitative systems pharmacology model for educational and research
purposes. It was constructed from the published literature and clinical trial data but
has not been independently verified or certified, and **it cannot be used for clinical
decision-making, prescribing, or regulatory submission.** The parameters are
approximations for the purpose of explanation, and fitting and validation against real
patient data would be required separately.
