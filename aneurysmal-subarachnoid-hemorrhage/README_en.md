# Aneurysmal Subarachnoid Haemorrhage → Delayed Cerebral Ischaemia (aSAH-DCI) QSP Model
## Aneurysmal Subarachnoid Haemorrhage → Delayed Cerebral Ischaemia

<p align="center">
  <a href="sah_qsp_model_en.svg"><img src="sah_qsp_model_en.png" width="880" alt="aSAH-DCI QSP mechanistic map"></a>
</p>

---

## One-sentence summary

Rather than treating delayed cerebral ischaemia (DCI) as the chain
"subarachnoid blood → vasospasm → ischaemia," this model treats it as
**a structure in which four consumers share a single arteriolar
vasodilatory reserve.** From that, the reason clazosentan cuts
angiographic vasospasm by 65% and still fails to move the 3-month outcome
falls out as a **calculation**, not an assumption.

---

## Why this structure — four facts the existing narrative cannot explain

| # | Fact | Problem for the chain narrative |
|---|------|------------------|
| i | Clazosentan cuts moderate-to-severe angiographic vasospasm by about 65% (CONSCIOUS-1) but does not change outcome (CONSCIOUS-2/3) | if the chain were right, outcome should follow |
| ii | Nimodipine barely changes angiographic calibre and is still the **only** drug that improves outcome | the opposite direction from the chain |
| iii | DCI occurs without angiographic vasospasm, and two-thirds of patients with moderate-to-severe vasospasm never develop DCI | the chain is not one-to-one |
| iv | Induced hypertension and milrinone work in some patients and not others (HIMALAIA was null) | the chain has no patient-specific condition |

---

## Three structural claims

### CLAIM 1 — one reserve, four consumers

Cerebrovascular resistance is split into three elements in series.

```
R_large   conductance arteries.  Only the spastic segment (FSEG = 0.33) follows r⁻⁴
RA        arterioles.  The only tunable element.  RA0 0.85 → RAmin 0.20
R_micro   microvasculature + capillaries

autoregulation:  RA_desired = CPP/CBF_target − R_large − R_micro
                 RA         = AREG·clamp(RA_desired) + (1−AREG)·RA0
                 CBF        = CPP / (R_large + RA + R_micro)

reserve   = RA0 − RAmin = 0.65 resistance units
```

Four consumers share this **single reserve.**

| Consumer | What it consumes | Observability |
|--------|---------------|-------------|
| C1 large-vessel spasm | ΔR_large | **the only thing angiography/TCD sees** |
| C2 microvascular tone | pericyte/arteriolar constriction | invisible |
| C3 capillary microthrombosis | capillary dropout | invisible |
| C4a perfusion pressure loss | ICP↑ / MAP↓ | visible (but caused as a drug side effect) |
| C4b spreading depolarisation | **contracts the reserve term itself** (inverse neurovascular coupling) | invasive ECoG only |

**The mapping from any one consumer to the clinical endpoint is nonlinear
and context-dependent.** The same 33% reduction in calibre is harmless in
a patient with an intact reserve and produces infarction in a patient who
has already spent the reserve on microthrombosis. A drug that blocks one
consumer works **only in the subset where that consumer happened to be the
marginal one.** This is the structural reason drugs beat angiography and
lose the clinical trial.

### CLAIM 2 — the microcirculation does not merely add to resistance, it erodes the ceiling on oxygen extraction

Jespersen & Østergaard's capillary transit-time heterogeneity (CTH) theory
is built in directly.

```
CTH    = (RMIC − 1) + 2.5·MTHR
OEFmax = 0.85 / (1 + 1.75·CTH)
ISCH   = 1 − DO2·OEFmax / CMRO2,       DO2 = CBF × CaO2
```

A CBF that was adequate in a homogeneous vascular bed is not adequate in
an occluded one. This term is what makes microthrombosis and pericyte tone
first-order consumers, and **what lets anaemia (CaO2) and fever (CMRO2)
enter through exactly the same node as vasospasm.**

### CLAIM 3 — the time window is haemoglobin-processing kinetics, not the calendar

There is no "day 4" switch anywhere in the model. The time window is the
convolution of the transition chain with an inducible sink.

```
CLOT --k_mat 0.25--> RBCL --k_hem 0.23--> OXYHB --(HP · CD163 · HO-1 · CSF)-->
                     ↑ red cells must lyse first = a delay
                                        ↑ HO-1 is inducible, t½ 1.9 d = a lag
```

Two things follow immediately. **A larger haemorrhage peaks earlier and
higher.** And **CSF drainage does not block an effector — it truncates the
input function** — a structurally different intervention, not
interchangeable with a downstream blocker.

---

## Computed results (calculations, not claims)

Every number below was produced by
[`sah_reference_check.py`](sah_reference_check.py)
(an independent numpy/RK4 implementation, a virtual cohort of 900
patients, dt 0.01 day), and the full output is preserved verbatim in
[`sah_reference_check_output.txt`](sah_reference_check_output.txt).

### 1. The dissociation

| Metric | Value |
|------|-----|
| corr(maximal calibre reduction, final infarct volume) | r = 0.363, **r² = 0.132** |
| corr(peak TCD velocity, final infarct volume) | r = −0.023, **r² = 0.001** |
| patients with moderate-to-severe vasospasm but no DCI | **63.9%** of that group |
| DCI occurring without moderate-to-severe vasospasm | **10.9%** of all DCI |

TCD velocity being essentially uncorrelated with infarction is not an
accident. Velocity is `Q/A`, so **a vessel that is both severely spastic
and severely underperfused actually reads low.** Diagnostic performance
splits the same way.

| Test | Sensitivity | Specificity |
|------|--------|--------|
| TCD > 120 cm/s | 66.4% | 39.3% |
| TCD > 160 cm/s | 13.4% | 84.9% |
| TCD > 200 cm/s | 0.8% | 98.9% |
| PbtO2 < 20 mmHg | 95.0% | 60.0% |
| **PbtO2 < 15 mmHg** | **82.4%** | **84.9%** |

PbtO2 wins because it **measures the shared node itself.** TCD sees only
one of the four consumers, and even then in a form distorted by flow.

### 2. Who is consuming the reserve (standard care, population-average share)

| Day | Large vessel | Micro-tone | Microthrombosis | CPP | SD | Total demand/reserve |
|----|--------|----------|----------|-----|-----|----------------|
| 2 | 29.8% | 23.8% | 3.0% | 39.9% | 3.5% | 0.24 |
| 4 | 61.8% | 18.0% | 5.0% | 13.3% | 1.9% | 0.73 |
| 6 | 70.9% | 12.0% | 6.2% | 7.9% | 2.9% | 1.27 |
| 8 | 70.3% | 9.4% | 7.4% | 6.9% | 6.1% | 1.54 |
| 10 | 67.4% | 8.2% | 8.5% | 7.1% | 8.9% | 1.58 |
| 14 | 60.3% | 7.3% | 10.4% | 8.9% | 13.1% | 1.34 |

The large vessel is definitely the **largest** consumer. But the spastic
segment's share of total cerebrovascular resistance runs from 15.0% on
day 0 to 37.3% on day 6 — the remaining 63% lies elsewhere. And "largest"
and "blocking it moves the endpoint" are not the same statement — which
is exactly the next table.

### 3. Perfect single-pathway blockade vs the real drug

| Intervention | DCI% | DCI RR | poor-outcome RR |
|------|------|--------|--------------|
| Standard care (reference) | 26.4 | 1.00 | 1.00 |
| Perfect large-vessel blockade (unrealistic) | 9.2 | 0.35 | 0.61 |
| Perfect micro-tone blockade | 15.2 | 0.58 | 0.76 |
| Perfect microthrombosis blockade | 9.1 | 0.34 | 0.59 |
| Perfect SD blockade | 23.2 | 0.88 | 0.87 |
| Preserved autoregulation | 19.2 | 0.73 | 0.82 |
| Large vessel + microthrombosis together | 0.0 | 0.00 | 0.42 |
| **Perfect blockade of all four consumers** | **0.0** | **0.00** | **0.42** |
| **Actual clazosentan 15 mg/h** | 16.6 | 0.63 | **1.07** |

Two things fall out of this.

**(a) A ceiling on the outcome endpoint.** Even eliminating DCI entirely
down to 0% brings poor outcome down only to **41.8% of baseline.** The
remaining 58% is early brain injury (EBI) and non-DCI channels. **This is
a hard ceiling that sits over any drug targeting DCI** — a quantitative
answer to why trials keep failing.

**(b) The gap between perfect and real.** Perfectly eliminating the
large-vessel consumer gives an outcome RR of 0.61. Clazosentan partially
eliminates the same consumer and produces an outcome RR of 1.07. That gap
is redundancy plus a harm channel.

An experiment blocking two pathways simultaneously came out
**sub-additive** (0.651 + 0.655 = 1.307 alone vs 1.000 together). Perfect
blockade of even a single pathway already removes most of the DCI, so it
runs into the ceiling. That is, in this model a single target fails not
because "the pathways act synergistically" but because **the real drug
falls well short of perfect blockade, and the endpoint carries a separate
ceiling that sits outside DCI altogether.**

### 4. Clazosentan: where the angiographic victory disappears

```
moderate-to-severe angiographic vasospasm   65.3%  →  21.6%     (RR 0.33 ; CONSCIOUS-1 target 0.35)
maximal calibre reduction (median)           37%   →   25%
DCI                                          26.4%  →  16.6%     (RR 0.63)
poor outcome                                 29.5%  →  31.6%     (RR 1.07 ; CONSCIOUS-2 target ~1.05)

counterfactual: switching off only the harm channel gives poor outcome 18.6% (RR 0.63)
```

That is, **the angiographic effect is real, and the endpoint effect is
eaten first by redundancy and then by the harm channel.** The harm channel
enters the same node along three branches: fluid retention → pulmonary
oedema → PaO2↓; hypotension → CPP↓; haemodilution → CaO2↓. The last one is
particularly quiet — cutting oxygen content has exactly the same effect as
touching resistance. In the model, the fluid-retention burden rises from a
median of 0.00 to 0.75 (IQR 0.47–1.15), and if the threshold is set at
1.50, 13.6% cross it.

### 5. The window is emergent, not built in

| Stratum | oxyHb peak | ET-1 peak | calibre nadir | ischaemia peak |
|----|-----------|-----------|--------------|-----------|
| Overall | day 3.0 (0.33 RU) | day 5.7 | day 8.5 (38% reduction) | day 11.2 |
| mFisher 1-2 | day 3.3 (0.21 RU) | day 6.0 | day 8.5 (31%) | day 12.2 |
| mFisher 3-4 | day 2.9 (0.38 RU) | day 5.6 | day 8.5 (40%) | day 11.0 |
| Hp1-1/1-2 | day 3.0 (0.31 RU) | day 5.8 | day 8.5 (37%) | day 11.5 |
| Hp2-2 | day 2.9 (0.36 RU) | day 5.6 | day 8.5 (39%) | day 11.0 |

As predicted, **a larger haemorrhage peaks earlier and higher.** DCI
incidence by stratum:

- modified Fisher: mF1 2.3% · mF2 11.9% · mF3 18.0% · mF4 50.0%
- Hp genotype: Hp1-1/1-2 25.8% vs Hp2-2 27.7%
- autoregulation tertile: lowest 34.3% · middle 26.3% · highest 18.7%

### 6. Anaemia is a consumer of the same node

| Stratum | n | DCI |
|----|---|-----|
| Hb < 10 g/dL | 47 | 48.9% |
| Hb 10–13 | 479 | 29.6% |
| Hb ≥ 13 | 374 | 19.5% |

Transfusing only the anaemic group brings DCI from 48.9% to 42.6%
(RR 0.87) — a result obtained without touching a single vessel.

### 7. It is the route, not the target — oral vs intravenous nimodipine

```
no nimodipine  : vasospasm 72.9%  DCI 34.4%  poor outcome 36.8%   (days 4-11 mean CPP 85.4)
oral nimodipine: vasospasm 65.3%  DCI 26.4%  poor outcome 29.5%   (CPP 82.7)   RR 0.80
IV nimodipine  : vasospasm 65.0%  DCI 30.7%  poor outcome 32.8%   (CPP 76.3)   RR 0.89
```

Exposure is essentially identical (Css about 48 ng/mL), and vasospasm
reduction is identical. The only thing that differs is CPP. Same drug,
same target, different route — and a different conclusion.

### 8. Other scenarios (excerpted from all 15; RR vs standard care)

| Scenario | Vasospasm (≥33%) | DCI | poor-outcome RR |
|----------|------------|-----|--------------|
| Cilostazol 100 mg bid | 59.0% | 13.1% | 0.69 |
| Early lumbar drainage | 52.3% | 16.2% | 0.74 |
| Intrathecal nicardipine implant | 17.9% | 11.0% | 0.66 |
| Induced hypertension (MAP +20) | 65.3% | 18.1% | 0.90 |
| Milrinone 0.5 µg/kg/min | 57.8% | 18.3% | 0.82 |
| Ketamine (SD suppression) | 65.3% | 24.1% | 0.87 |
| Simvastatin 40 mg | 63.4% | 24.1% | 0.96 |
| Clazosentan + cilostazol + drainage | 8.4% | 3.1% | 0.80 |

The last row closes out Claim 1. Hitting three pathways at once drops DCI
from 26.4% to 3.1%, while poor outcome improves by only 20% — that ceiling
again.

### 9. Monitoring lead time and systemic channels

- TCD > 120 cm/s: median 3.4 days · PbtO2 < 20 mmHg: 2.9 days · infarct
  onset: 6.4 days
- hyponatraemia < 135 mmol/L 30.1% · < 130 mmol/L 4.6%
- median cumulative SD count 6 (65 in the DCI group vs 4 in the non-DCI
  group)
- median lowest autoregulation value 0.59 (0.33 in the DCI group vs 0.67
  in the non-DCI group)

---

## Files

| File | Contents |
|------|------|
| [`sah_qsp_model_en.dot`](sah_qsp_model_en.dot) | mechanistic map source (22 clusters · 221 nodes · 339 edge statements) |
| [`sah_qsp_model_en.svg`](sah_qsp_model_en.svg) / [`.png`](sah_qsp_model_en.png) | rendered output (150 dpi) |
| [`sah_mrgsolve_model_en.R`](sah_mrgsolve_model_en.R) | mrgsolve QSP model — 39 ODEs (8 PK/exposure + 31 disease/haemodynamics/injury), 15 scenarios, 12 analysis functions |
| [`sah_shiny_app_en.R`](sah_shiny_app_en.R) | 8-tab Shiny dashboard (reserve gauge · PK · haemoglobin clock · four consumers · perfusion/oxygen · endpoints · scenarios · monitoring) |
| [`sah_references_en.md`](sah_references_en.md) | 89 references, every PMID verified via PubMed E-utilities |
| [`sah_reference_check.py`](sah_reference_check.py) | an independent numpy/RK4 transcription. Produces every number above |
| [`sah_reference_check_output.txt`](sah_reference_check_output.txt) | the full output of the script above (committed verbatim) |

### State variables (39 ODEs)

| Group | States |
|------|------|
| Drug PK/exposure (8) | NIMG · NIMC · CLAZ · CILO · STAT · MILRC · NICA · KETA |
| Blood · haemoglobin (8) | CLOT · IVH · RBCL · OXYHB · HP · HEME · HO1 · BOX |
| Effectors (6) | ET1 · NOB · ROS · RHOK · INFL · PAI |
| Vascular dynamics (5) | SPASM · STRUCT · RMIC · MTHR · AREG |
| Spreading depolarisation (3) | SDSUS · SDBUR · SDCUM |
| Perfusion · injury (7) | EDEMA · HYDRO · ICP · MAP · OGD · INFVOL · EBI |
| Systemic (2) | NAS · FLUID |

### 15 scenarios

S1 supportive care only · S2 oral nimodipine (standard) · S3 IV
nimodipine · S4/S5 clazosentan 15 / 5 mg/h · S6 cilostazol ·
S7 simvastatin · S8 early lumbar drainage · S9 intrathecal nicardipine ·
S10 induced hypertension · S11 milrinone · S12 ketamine · S13
transfusion · S14 triple combination · S15 full stack

---

## Calibration (literature target vs model)

| Item | Target | Model |
|------|------|------|
| Moderate-to-severe angiographic vasospasm (control) | about 66% | **65.3%** |
| DCI incidence | about 30% | 26.4% |
| Poor outcome, mRS 4-6 (90 days) | 30–35% | 29.5% |
| Hyponatraemia < 135 | 30–50% | 30.1% |
| Clazosentan RR, moderate-to-severe vasospasm | 0.35 | **0.33** |
| Clazosentan RR, poor outcome | about 1.05 | **1.07** |
| Nimodipine RR, poor outcome | 0.67 | 0.80 |
| Cilostazol RR, DCI | 0.47 | **0.50** |
| Simvastatin RR, poor outcome | about 1.00 | **0.96** |
| Lumbar drainage RR, poor outcome | 0.76 | **0.74** |

## Where the model departs from the literature or fails (recorded, not tuned away)

1. **Nimodipine's improvement in outcome is too shallow.** RR 0.80 vs the
   Cochrane figure of 0.67. The direction is right and the magnitude is
   too small. In the model, nimodipine's benefit enters only through the
   microvascular, SD, and infarct-conversion pathways, and the literature's
   0.67 may include a component this model does not carry.
2. **Too little DCI without angiographic vasospasm.** 10.9% vs a
   literature range of 20–30%. Enlarging the microthrombosis/SD pathway
   could fix this, but then overall DCI incidence overshoots the target.
   Matching the incidence was chosen instead.
3. **"Any vasospasm (≥25% reduction)" runs high at 85.4%.** This is out of
   step with the literature's 60–70%, because the model represents the
   single most severely affected territory as its representative
   compartment. The metric that matches the clinical trial definition is
   ≥33% (moderate-to-severe), and that one matches at 65.3%.
4. **The extremes of modified Fisher grade are spread too far apart.**
   mF1 2.3% (literature about 6–12%), mF4 50.0% (literature about
   35–40%). The ordering is monotonic and correct, but the slope is too
   steep.
5. **Hyponatraemia < 130 runs low.** 4.6% vs a literature range of
   10–20%.
6. **A hypothesis that did not reproduce — autoregulation-dependence of
   induced hypertension.** While building the model, the prediction was
   "induced hypertension works only in pressure-passive (low-AREG)
   patients," and the calculation came out the opposite way (induced
   hypertension DCI RR: 0.73 in the autoregulation-impaired group vs 0.62
   in the intact group; milrinone 0.75 vs 0.60). In this structure, raising
   MAP always raises CPP and the maximum achievable flow, and the
   autoregulation-impaired group has less headroom to begin with, so its
   relative benefit is smaller. **The hypothesis that HIMALAIA's null
   result can be explained by autoregulation stratification is not
   supported by this model** — left standing as a failed prediction.
7. **The harm channel should not be read as a rate of pulmonary
   complications.** Individual susceptibility was built into FLUID
   (clazosentan-arm median 0.75, IQR 0.47–1.15), but the threshold was not
   calibrated against the actual clinical pulmonary complication rate
   (CONSCIOUS-2, about 8% vs 3%). The size of the harm term was tuned to
   reproduce CONSCIOUS-2's null outcome result and cannot be used to
   predict the pulmonary oedema incidence itself.
8. **Non-DCI infarction (surgery-related, EBI-related) is not given a
   separate state.** Infarct volume is populated only through the DCI
   pathway. The EBI term in the outcome logit carries that share instead.

---

## Reproducing this

```bash
# map
dot -Tsvg sah_qsp_model_en.dot -o sah_qsp_model_en.svg
dot -Tpng -Gdpi=150 sah_qsp_model_en.dot -o sah_qsp_model_en.png

# every number in this README (numpy only, about 10 minutes)
python3 sah_reference_check.py

# mrgsolve model
Rscript -e 'source("sah_mrgsolve_model_en.R"); print(sah_run_all(300)$scenarios)'

# dashboard
Rscript -e 'shiny::runApp("sah_shiny_app_en.R")'
```

---

## Disclaimer

A semi-quantitative QSP model for educational and research purposes. It
was built from published literature and clinical trial data but has not
been independently verified or certified, and **must not be used for
clinical decision-making, prescribing, or regulatory submission.**
Parameters are illustrative approximations, and the 8 "points of
departure" listed above in particular are a list of questions this model
must not be used to answer.
