# Immune Checkpoint Inhibitor-Induced Colitis (ICI colitis / IMDC) — QSP Model

<p align="center">
  <a href="../../../immune-checkpoint-inhibitor-colitis/icic_qsp_model.svg">
    <img src="../../../immune-checkpoint-inhibitor-colitis/icic_qsp_model.png" width="820" alt="ICI colitis QSP map">
  </a>
</p>

> **In one sentence.** Checkpoint-inhibitor colitis is not "a state of excessively
> activated immunity" but **a single ratio whose denominator can be depleted**, and
> from the one fact that the two drug classes enter that fraction at different
> positions there follow the two things the clinical data demand — that CTLA-4
> inhibitors have a dose-response and PD-1 inhibitors do not.

---

## 1. Why this disease

There are two observations in checkpoint-inhibitor colitis that the standard
"inflammation model" cannot explain.

1. **ipilimumab has a dose-response and anti-PD-1 does not.** From ipilimumab
   3 → 10 mg/kg, G3-4 colitis rises from 3.2% → 6.7% (Ascierto 2017), whereas
   across nivolumab 0.1-10 mg/kg and pembrolizumab 2 vs 10 mg/kg there is **no
   dose difference** in immune-related adverse events.
2. **ipilimumab + nivolumab colitis is not the sum of the two single agents.** In
   CheckMate 067 the combination is 7.7%, ipilimumab alone 8.7%, nivolumab alone
   0.6%.

If the two drugs were merely stronger and weaker versions of each other on the same
axis, both observations would be strange. This model takes the two drugs to enter
at **structurally different positions**, and derives (1) from there. It could not
derive (2), and that failure is written up as it stands in §6 below.

---

## 2. The central identities

### 2.1 Colonic inflammatory drive as a ratio

```
Φ_col = (a_eff · Teff + a_trm · Trm) / (a_reg · Treg + κ)
```

* **The numerator — anti-PD-1.** Per-cell activity is multiplied by the occupancy
  `O = C_t/(C_t + K_D)`. With K_D ≈ 1 nM and colonic tissue concentrations of
  2-200 nM, the computed occupancies are as follows.

  | nivolumab dose | plasma C_avg | tissue C_t | tissue nM | occupancy O | numerator multiplier Ppd1 |
  |---|---|---|---|---|---|
  | 0.1 mg/kg | 2.19 µg/mL | 0.285 | 2.0 | 0.661 | 1.529 |
  | 0.3 mg/kg | 6.58 | 0.855 | 5.9 | 0.854 | 1.683 |
  | 1.0 mg/kg | 21.9 | 2.85 | 19.5 | 0.951 | 1.761 |
  | 3.0 mg/kg | 65.8 | 8.55 | 58.6 | 0.983 | 1.787 |
  | 10 mg/kg | 219 | 28.5 | 195 | 0.995 | 1.796 |

  **0.3 → 10 mg/kg is a 33-fold increase in dose, yet the numerator multiplier
  moves 1.683 → 1.796, that is by only 6.7%.** There is no dose-response to find.
  The drug is already on the plateau across the whole clinical range in which a
  patient will be placed.

* **The denominator — anti-CTLA-4.** CTLA-4 occupancy also saturates (0.603 →
  0.938 from 1 → 10 mg/kg). What does **not** saturate is FcγRIIIa-mediated ADCC.
  Fc-mediated killing demands not simple occupancy but **immune-complex density**,
  so its effective EC50 is 10-100 times the binding K_D, and free drug in tissue is
  only 13% of plasma. That is what places EC50_ADCC **in the middle** of the
  clinical exposure range.

  | ipilimumab | tissue C_t | CTLA-4 occupancy | ADCC drive | steady-state Treg | release factor G |
  |---|---|---|---|---|---|
  | 1 mg/kg | 1.18 µg/mL | 0.603 | 0.228 | 0.756 | 1.269 |
  | 3 mg/kg | 3.54 | 0.820 | 0.470 | 0.600 | 1.532 |
  | 6 mg/kg | 7.08 | 0.901 | 0.639 | 0.525 | 1.704 |
  | 10 mg/kg | 11.8 | 0.938 | 0.747 | 0.486 | 1.809 |

  Because G = Reg₀/Reg is a **hyperbola**, it goes on rising after every occupancy
  term has flattened. What sets the ceiling of the hyperbola is **κ (the non-Treg
  regulatory floor)**, with G_max = Reg₀/κ = 7.67. κ is the most influential
  parameter in the model and the hardest to measure directly.

**A ratio with a saturated numerator and a depletable denominator** — that is the
whole of why one class has a dose-toxicity curve and the other does not.

### 2.2 The colon exhausts its reserve before the first symptom

Stool water is **the difference of two large numbers**.

```
Stool = max(Stool_min,  L + Sec − A_max · S_eff)
S*    = L / A_max = 1750 / 4500 = 0.3889
```

Against an ileal effluent L ≈ 1750 mL/day, the maximal colonic absorptive capacity
is A_max ≈ 4500 mL/day — a reserve of **2.57-fold**. So **61.1% of absorptive
function can be lost together with an entirely normal stool diary.** Converting the
grade boundaries into loss of absorptive function:

| Grade | Increase in stools/day | Stool water | S_eff required | Fraction lost |
|---|---|---|---|---|
| G1 onset | 1.0 | 255 mL | 0.332 | 66.8% |
| G2 onset | 4.0 | 720 mL | 0.229 | 77.1% |
| G3 onset | 7.0 | 1185 mL | 0.126 | 87.4% |

Three conclusions follow, and all of them are clinical.
1. **Grade is a late instrument with a large loss.** Calprotectin and endoscopy read
   the lesion itself, so they are not censored, and they therefore move weeks
   earlier.
2. **Two patients at the same grade can have entirely different reserve left.**
   Grade is a poor measure of severity and a poor trigger for treatment.
3. **Symptom-triggered treatment is late structurally, not through clinical
   carelessness.**

### 2.3 Steroid refractoriness is a statement about the crypt, not about the drug

Once injury exceeds `inj_deep`, the stem/progenitor cells at the base of the crypt
die, and they regenerate with a **20-day time constant** rather than a 3-day one.
From then on, however completely the immune drive is switched off, the absorptive
surface is trapped at `ISC · Rep`. A patient whose crypts are gone cannot be
rescued by turning immunosuppression up — which is why deep ulceration at endoscopy
predicts steroid failure, and why a minority end up at colectomy.

---

## 3. Deliverables

| File | Contents |
|------|------|
| [`icic_qsp_model.dot`](../../../immune-checkpoint-inhibitor-colitis/icic_qsp_model.dot) · [`.svg`](../../../immune-checkpoint-inhibitor-colitis/icic_qsp_model.svg) · [`.png`](../../../immune-checkpoint-inhibitor-colitis/icic_qsp_model.png) | Mechanistic map — **146 nodes, 12 clusters, 238 edges** |
| [`icic_reference_model.py`](../../../immune-checkpoint-inhibitor-colitis/icic_reference_model.py) | **The reference implementation — the numerical ground truth.** The file that computes every number below |
| [`icic_mrgsolve_model.R`](../../../immune-checkpoint-inhibitor-colitis/icic_mrgsolve_model.R) | mrgsolve 45-state ODE model (a 1:1 translation of the file above) + a 13-scenario driver |
| [`icic_shiny_app.R`](../../../immune-checkpoint-inhibitor-colitis/icic_shiny_app.R) | 8-tab Shiny dashboard |
| [`icic_references.md`](icic_references.md) | 98 PubMed links, with the source marked for each parameter |
| [`icic_reference_output.txt`](../../../immune-checkpoint-inhibitor-colitis/icic_reference_output.txt) | The full analysis output (A-L) |
| [`icic_cross_validation.txt`](../../../immune-checkpoint-inhibitor-colitis/icic_cross_validation.txt) | Python ↔ mrgsolve comparison and structural verification |

### State variables (45 ODEs)

| Compartment | States |
|------|------|
| PK (14) | ipilimumab·anti-PD-1·infliximab·vedolizumab·tocilizumab, two compartments each; prednisolone gut/central/effect site; JAK inhibitor |
| Immune (8) | Tn · Nclone (priming integral) · Teff · **Trm** · **Treg** · Th17 · Mac · Neut |
| Cytokines (8) | IFN-γ · TNF-α · IL-6 · IL-23 · **IL-15** · CXCL10 · MAdCAM-1 · IL-10 |
| Epithelium (4) | **ISC** · Ent (absorptive mass) · TJ · Muc |
| Lumen (3) | diversity Dv · butyrate Bu · translocated LPS |
| Systemic · outcomes (8) | **albumin** · CRP · calprotectin · ulceration · tumour · intratumoural effectors · cumulative steroid · infection risk |

---

## 4. Computed results

Every number was computed by `icic_reference_model.py` and stands verbatim in
`icic_reference_output.txt`. **The deterministic scenarios are run not on the median
patient but on an explicitly stated susceptible phenotype (`P_case()`)** — only 7-12%
of patients develop colitis on ipilimumab 3 mg/kg, so a model in which the median
patient gets colitis is a wrong model.

### 4.1 Structural self-validation

The undosed mucosa is an exact steady state to **machine precision**
(scaled |dy/dt| ≈ 9×10⁻¹⁶, no drift over a 365-day drug-free run). The two
elimination rate constants are not free parameters: they are obtained by **solving**
for this condition.

What makes this possible is the **mucosal tolerance band**. Every amplification loop
(IFN-γ → CXCL10 → macrophage → TNF → MAdCAM → influx → more effectors) only
operates once its input leaves the tolerance band. Without that band the open-loop
gain of the resting mucosa exceeds 1 and **any fluctuation whatsoever ends in
colitis**, which is a wrong story both numerically and biologically.

### 4.2 The censoring ladder — susceptible patient, ipilimumab 3 mg/kg

Each column is the day on which that instrument first said something.

| Regimen | Lesion (S↓10%) | Calpro >150 | >200 | G1 | G2 | G3 | Lowest S_eff |
|---|---|---|---|---|---|---|---|
| ipi 3 mg/kg | day 2 | 44 | 45 | 45.5 | **53** | 71.5 | 0.045 |
| ipi 10 mg/kg | day 2 | 26 | 26.5 | 26 | **30** | 94 | 0.161 |
| ipi3 + nivo1 | day 2 | 23.5 | 24 | 23 | **26.5** | 79 | 0.134 |
| nivolumab 3 mg/kg | — | — | — | — | — | — | 0.789 |
| ipi 1 mg/kg | — | — | — | — | — | — | 0.789 |

The day G2 is reached on ipilimumab 3 mg/kg is **day 53** — consistent with the
reported median of 6-7 weeks, and a value that came out rather than one the model
was told. And as the dose rises, onset moves earlier (53 → 30 → 26.5 days).

Look at the detailed time course and the censoring is plain to see. Between days 28
and 42, S_eff falls 0.80 → 0.45 while **the stool diary is entirely normal**, and
calprotectin only begins to move on day 42, reaching 437 µg/g by day 49.

### 4.3 The depth of the step is the onset time

The model was never told the onset times. The agents differ only in their
**distance** from the terminal event, epithelial apoptosis.

| Rescue therapy | Site of action | Model onset | Reported median |
|---|---|---|---|
| infliximab | TNF-α neutralisation = **the executioner** | 0.5 days | 2-3 days |
| corticosteroid | transcriptional programme | 2.0 days | 4-6 days |
| vedolizumab | α4β7 **trafficking blockade** (cannot touch cells already resident) | **36.5 days** | 14-21 days |

**The order is reproduced but the absolute values are compressed.** infliximab and
the steroid are faster than reported and vedolizumab slower. This appears to be
because the fast component of the symptom (secretory · transporter) carries a larger
share in the model than it does in reality, and it is written up as a limitation in
§6.

### 4.4 Selectivity is the currency (the tumour cost of rescue therapy)

The clone that damages the colon is the clone that kills the tumour. Every rescue
therapy therefore withdraws from the tumour account **in proportion to its
non-gut-selectivity χ**.

| Regimen | d180 tumour burden | vs no rescue | Cumulative steroid | Infection risk |
|---|---|---|---|---|
| ICI alone (no rescue needed) | 33.7 | 100% | 0 mg | 0.000 |
| Steroid only | 72.8 | 216% | 1487 mg | 0.075 |
| Steroid + infliximab (χ=0.60) | 293.0 | 870% | 1487 mg | 0.297 |
| Steroid + **vedolizumab** (χ=0.10) | 87.0 | **258%** | 1487 mg | 0.132 |
| vedolizumab + steroid sparing (0.5 mg/kg, 2-week taper) | 57.3 | **170%** | **429 mg** | 0.103 |

α4β7:MAdCAM-1 is **a gut address the tumour does not use**. Buy the same control of
colitis at the gut address and the tumour cost is far cheaper.

### 4.5 The disease eats its own antidote

`CL_ifx = 0.40 · (ALB/4.0)^−0.9 · (1 + 0.006·CRP)` (Fasanmade 2009).
Severe colitis causes a protein-losing enteropathy and drives albumin down.

| Albumin | CRP | CL_ifx | t½ | d14 trough | Tissue C_t | TNF neutralisation |
|---|---|---|---|---|---|---|
| 4.2 | 5 | 0.394 L/day | 14.8 days | 19.8 µg/mL | 1.98 | 88.5% |
| 3.5 | 40 | 0.559 | 11.1 | 13.9 | 1.39 | 86.8% |
| 3.0 | 80 | 0.767 | 8.8 | 9.01 | 0.90 | 83.7% |
| 2.6 | 120 | 1.014 | 7.4 | 5.53 | 0.55 | 78.8% |
| 2.2 | 150 | 1.302 | 6.5 | 3.25 | 0.33 | **71.1%** |

A patient with an albumin of 2.2 g/dL clears infliximab **3.30 times faster** than
one at 4.2 g/dL. A fixed 5 mg/kg delivers **the lowest exposure to the sickest
patient**. This is not an argument for watchful waiting in severe steroid-refractory
colitis; it is the mechanistic argument for shortening the interval and escalating
the dose.

### 4.6 Hysteresis — why 30-40% flare during the taper

Trm self-renews on epithelial IL-15, and IL-15 is raised by injury. **The drug is
the trigger and Trm is the state.** Stop the drug and the state does not go away.

| Taper duration | Cumulative steroid | Peak grade | Trm at end of taper | d180 Trm | Flare after taper | Final grade |
|---|---|---|---|---|---|---|
| 2 weeks | 998 mg | G3 | 0.062 | 0.146 | **yes** | G2 |
| 4 weeks | 1487 mg | G3 | 0.050 | 0.136 | **yes** | G2 |
| 6 weeks | 1977 mg | G2 | 0.042 | 0.055 | no | G0 |
| 8 weeks | 2468 mg | G2 | 0.036 | 0.052 | no | G0 |
| 12 weeks | 3448 mg | G2 | 0.027 | 0.045 | no | G0 |

A long taper lowers the flare risk but buys cumulative steroid and infection risk
along with it. A gut-selective agent buys the same durability **without either**
(§4.4).

### 4.7 Starting on an uncensored instrument

| Starting criterion | Start day | Peak grade | Lowest S_eff | Peak ulceration | Peak calpro | d180 Trm |
|---|---|---|---|---|---|---|
| Symptom-based (G2) | 53 | G3 | 0.159 | 1.119 | 836 | 0.144 |
| Symptom-based (G1) | 50 | G3 | 0.171 | 1.099 | 824 | 0.142 |
| **Calprotectin >200** | **45** | **G0** | **0.484** | **0.221** | **227** | **0.076** |
| Calprotectin >150 | 44 | G2 | 0.188 | 1.073 | 809 | 0.140 |

Starting at calprotectin >200 — merely starting eight days earlier — turns the
outcome from G3 into G0, cuts ulceration to a fifth and halves the residual Trm.
**What is bought is not speed but durability.**

Note, however, that **the >150 criterion starts a day earlier than the >200
criterion and the outcome is worse** (G2 vs G0). This arises because the system is
bistable near the separatrix, and it shows that a single day can flip the category.
It is left in as it stands, as a warning that the results must not be read as a
smooth "the sooner the better".

### 4.8 Virtual population incidence (n = 60 per arm, identical parameter set)

| Regimen | Model G≥1 | G≥2 | G≥3 | Reported |
|---|---|---|---|---|
| ipilimumab 1 mg/kg | 0.0% | 0.0% | 0.0% | low |
| **ipilimumab 3 mg/kg** | 13.3% | 13.3% | **5.0%** | G3-4 colitis 3.2-8.7% |
| **ipilimumab 10 mg/kg** | 21.7% | 21.7% | **11.7%** | G3-4 colitis 6.7% |
| nivolumab 3 mg/kg | 0.0% | 0.0% | 0.0% | G3-4 colitis 0.6% |
| nivolumab 10 mg/kg | 1.7% | 1.7% | 1.7% | no dose effect |
| ipi 3 + nivo 1 | 48.3% | 48.3% | 35.0% | colitis 7.7% |

* **Success — the shape of the dose-response.** From 3 → 10 mg/kg, G≥3 rises
  5.0% → 11.7%, that is **2.3-fold**. The reported values go 3.2% → 6.7%, about
  2.1-fold. This is the model's headline quantitative result, and it is made **by
  the unsaturated ADCC term alone**, in a range where both occupancy terms are flat.
* **Success — the absolute level at 3 mg/kg.** G≥3 5.0% (reported 3.2-8.7%).
* **Success — the flatness of anti-PD-1.** 0.0% vs 1.7% at 3 vs 10 mg/kg.
* The **failures** are in §6.

### 4.9 Sensitivity — the model tests its own claim

The model claims that the PD-1 arm is **saturated**. That claim makes a strong
prediction about its own sensitivity analysis: shaking those parameters by ±25%
must leave the result **exactly** unchanged.

| Parameter | Lowest S_eff at −25% / +25% | Normalised sensitivity |
|---|---|---|
| `KD_PD1_nM` | 0.2699 / 0.2699 | **+0.000** |
| `Emax_pd1` | 0.2699 / 0.2699 | **+0.000** |
| `L_pres` (ileal effluent) | 0.1739 / 0.3284 | **+1.145** |
| `k_hom` | 0.2619 / 0.1173 | −1.071 |
| `BDC` (tissue penetration) | 0.2453 / 0.1039 | −1.048 |
| `EC50_ADCC` | 0.0968 / 0.2290 | **+0.980** |
| `A_max` (colonic absorptive capacity) | 0.3284 / 0.2170 | −0.825 |
| `b_but` (butyrate) | 0.2499 / 0.1466 | −0.766 |
| `κ` (non-Treg floor) | 0.1358 / 0.2358 | +0.741 |
| `Emax_ctla` | 0.2295 / 0.1344 | −0.705 |

**The influence of the two PD-1 binding parameters is exactly zero, while that of
`EC50_ADCC` is almost one.** The claim of "a saturated numerator, an unsaturated
denominator" has come out as a number the model was never told.

And look at the top of the list. **`L_pres` and `A_max` — two physiological
constants with nothing whatever to do with immunity — come ahead of every immune
parameter.** Because they set the censoring threshold S* = L/A_max. If that is
right, then whether a patient is **reported** as having colitis matters as much on
baseline bowel physiology as on immune phenotype.

---

## 5. Scenarios

Thirteen are implemented in the driver at the foot of `icic_mrgsolve_model.R`, and
interactively in the Shiny app.

1-3. ipilimumab 1 / 3 / 10 mg/kg q3w ×4 (dose-response)
4-5. anti-PD-1 3 / 10 mg/kg q2w (**demonstration of flatness**)
6. ipi 3 + nivo 1 q3w ×4 → nivolumab 480 mg q4w
7. prednisolone 1 mg/kg started at G2 + a 4-week taper
8. steroid-refractory → infliximab 5 mg/kg (weeks 0·2·6)
9. steroid-refractory → vedolizumab 300 mg (weeks 0·2·6)
10. calprotectin-based early start
11. FCGR3A V/V · V/F · F/F genotypes
12. antibiotic exposure before treatment (the diversity · butyrate axis)
13. comparison of the tumour cost of each rescue therapy

In addition, tocilizumab, tofacitinib, FMT, taper durations of 2-12 weeks and a
pre-existing microscopic colitis phenotype are supported as parameters.

---

## 6. Where the model fails

This section is mandatory. To believe the results above you must read the following
alongside them.

| # | Failure | What it means |
|---|------|----------------|
| 1 | **Combination therapy.** Observationally ipi3+nivo1 colitis ≈ ipilimumab alone (7.7% vs 8.7%), yet this model is plainly **supra-additive**. Because Ppd1 and G multiply, the non-additivity would have to emerge from threshold censoring alone, and this population variance is not enough for that | Either the per-cell contribution of PD-1 is smaller in already-primed colonic clones, or the two drugs do not act multiplicatively on the same pool |
| 2 | **The Monte Carlo has not converged.** At n=60 the standard error on a 13% estimate is about 4.3 percentage points, and an earlier n=120 run of the same model gave 9.2% and 39.2% at 3 and 10 mg/kg | The two runs differ by much more than the nominal standard error. The incidence figures should be read only at the level of an order of magnitude, and the structural results of §2 and §4.9 do not depend on them |
| 3 | **Colitis on anti-PD-1 alone is 0%.** The observation is 0.6-2% G3 colitis | In the model anti-PD-1 causes colitis only where a latent priming pool (Trm₀) exists. The direction is right (the association with pre-existing IBD) but the magnitude falls short |
| 4 | **The absolute onset times of the rescue therapies.** The order is reproduced, but infliximab and the steroid are too fast and vedolizumab too slow | The fast component of the symptom (secretory · transporter) carries too large a share |
| 5 | **The depth of Treg depletion is a latent variable.** There are reports that ipilimumab does not deplete Treg in human tumours (Sharma 2019) | The model back-calibrates the depth of colonic mucosal Treg depletion from the toxicity outcome, and this is directly testable through prediction 1 of §8 |
| 6 | **Bistability.** Outcomes jump categorically near the separatrix (§4.7) | Individual predictions must not be read as a smooth dose-response |

Item 1 is the sharpest of them. **Because it is testable.** What the model says the
non-additivity really is, is not a biological interaction but **the product of a
censored endpoint**. If that is right, then measuring the same patients on an
uncensored instrument (calprotectin, the Mayo endoscopic score) should make the
combination look additive again. If it is still non-additive on a continuous
measure, then the censoring explanation is wrong, and so is this part of the model.

The structural results of §2 and §4.9 (the flatness, the reserve arithmetic, the
order of the steps, the selectivity) do not depend on population variance or on the
incidence calibration.

---

## 7. Falsifiable predictions

The full table is at the end of `icic_references.md`. The three central ones:

1. **An Fc-silenced anti-CTLA-4 will abolish the dose-response for colitis** —
   because the dose-response comes from ADCC and not from occupancy. The antitumour
   effect should be largely retained.
2. **Patients with colitis on anti-PD-1 alone will have high colonic Trm/clonality
   from before treatment** — because the drug cannot licence new clones, there has
   to be something already present for it to release.
3. **The combination will look additive when measured on a continuous instrument** —
   if the non-additivity is due to threshold censoring by grade rather than to
   biology, then measured by calprotectin or the Mayo score it should come close to
   the sum. (This is also the experiment that decides failure 1 of §6.)

---

## 8. How to reproduce

```bash
# mechanistic map
dot -Tsvg icic_qsp_model.dot -o icic_qsp_model.svg
dot -Tpng -Gdpi=110 icic_qsp_model.dot -o icic_qsp_model.png

# the reference implementation (the source of every table)
python3 icic_reference_model.py > icic_reference_output.txt

# mrgsolve
R -e 'library(mrgsolve); mod <- mread("icic","icic_mrgsolve_model.R");
      mrgsim(mod, end=365, delta=5) |> tail()'   # confirm the undosed steady state

# Shiny dashboard
R -e 'shiny::runApp("icic_shiny_app.R")'
```

---

## ⚠️ Disclaimer

This is a **qualitative / semi-quantitative QSP model for educational and research
purposes**. It was built from the public literature and clinical-trial data but has
not been independently validated or certified, and **must not be used directly for
real clinical decision-making, prescribing, or regulatory submission.**
Quoting the numbers of §4 without reading the list of failures in §6 is a misuse of
this document.
