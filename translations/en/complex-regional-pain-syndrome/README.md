# Complex Regional Pain Syndrome (CRPS) — QSP Model

> A quantitative systems pharmacology model that treats the chronic pain syndrome
> arising in a single limb after a fracture or an operation **not as a list of
> mechanisms but as a topology**.
> Because it separates a rapidly extinguishing **peripheral inflammatory node** from
> two **positive feedback loops** with time constants of weeks to months (the
> behaviour-cortex loop and the glial latch), in this model a single parameter set has
> **three stable states**, and every treatment effect becomes a question of "which
> state does it go to / can it still be changed".

| Deliverable | File |
|---|---|
| 🗺️ Mechanistic map (DOT, 20 clusters · 178 nodes · 263 edges) | [`crps_qsp_model.dot`](crps_qsp_model.dot) |
| 🖼️ Map (SVG) | [`crps_qsp_model.svg`](crps_qsp_model.svg) |
| 🖼️ Map (PNG, 150 dpi) | [`crps_qsp_model.png`](crps_qsp_model.png) |
| ⚙️ mrgsolve ODE model (34 ODEs · 9 scenarios · 10 analysis functions) | [`crps_mrgsolve_model.R`](crps_mrgsolve_model.R) |
| 📊 Shiny dashboard (8 tabs) | [`crps_shiny_app.R`](../../../complex-regional-pain-syndrome/crps_shiny_app.R) |
| 🔬 Independent implementation for verification (numpy/scipy) | [`crps_reference_check.py`](../../../complex-regional-pain-syndrome/crps_reference_check.py) |
| 📚 References (70, with PubMed links) | [`crps_references.md`](crps_references.md) |

---

## 1. The disease in one paragraph

CRPS is a syndrome in which, after a usually trivial inciting event such as a distal
radius fracture, foot surgery or minor trauma, persistent pain out of proportion to the
degree of injury appears in a single limb together with autonomic features (temperature
and sweating asymmetry), oedema, motor impairment (reduced ROM, dystonia) and trophic
changes (the Budapest criteria). In the acute phase there is SP/CGRP-mediated
neurogenic inflammation with raised IL-6 and TNF-α, activation of keratinocytes and
mast cells, and vasodilatation, producing the **warm, oedematous** phenotype; in the
chronic phase it shifts to the **cold, ischaemic** form through up-regulation of
α1-adrenergic receptors, endothelial dysfunction (ET-1↑, NO↓) and central sympathetic
outflow. On top of that come functionally active autoantibodies against β2AR/M2/AT1R
(with causality suggested by passive transfer), NMDA-dependent spinal sensitisation
with glial activation, loss of descending inhibition (CPM/DNIC), shrinkage of the S1
somatosensory map with body-perception disturbance, and the disuse bone loss (patchy
osteopenia) that arises out of pain-fear-disuse. Treatment is built around early
multidisciplinary rehabilitation, with steroids, free radical scavengers,
bisphosphonates, gabapentinoids, ketamine and spinal cord stimulation (SCS) deployed
according to stage and phenotype — but the RCT results are notoriously inconsistent
with one another.
**The aim of this model is to explain that inconsistency in terms of differences in
system state rather than in the size of the drug effect.**

## 2. What the model asserts structurally

```
(A) fast peripheral node (feed-forward, time constant = days)
    INJ(t) = exp(-t/400h) → SP/CGRP → cytokines → NGF → ROS · oedema · α1↑
                                                   → peripheral sensitisation PSENS
    → once the input goes, it extinguishes itself. On its own it cannot become chronic.

(B1) behaviour-cortex loop (positive feedback, time constant = weeks to months)
    PAIN --(KFEAR, Hill n=3)--> DISUSE --(Hill n=3)--> CORTEX --(W_CTX_PAIN)--> PAIN
                                    └--> hypoxia / bone remodelling --> PSENS

(B2) glial latch (glial priming switch)
    SSENS --> GLIA --(Hill n=16, GLIA50=0.36)--> SSENS
    KOUT_GLIA = 8e-4 /h → a memory of about 5 weeks. Past the threshold, spinal
    sensitisation sustains itself with no afferent input at all.
```

Because of this topology, **a single parameter set has three stable states**.

| Stable state | NRS at 3 years | CSS | GLIA | Interpretation |
|---|---|---|---|---|
| attractor 0 | 0.0 | 0.0 | 0.00 | resolution (full recovery, BMD 1.00, ROM 1.00) |
| attractor 1 | 5.5 | 4.6 | 0.23 | a chronic state with only the behaviour-cortex loop latched |
| attractor 2 | 6.9 | 6.8 | 0.59 | a severe chronic state that has gone on to the glial latch |

## 3. The mechanistic map — 20 clusters (178 nodes · 263 edges)

1. Inciting event and risk factors (fracture, surgery, immobilisation, initial pain intensity, ACE inhibitors, HLA)
2. Tissue and nerve injury and subtype (CRPS I/II, small-fibre neuropathy, IENFD↓, TLR4)
3. Neurogenic inflammation (SP/NK1, CGRP/CLR-RAMP1, plasma extravasation, keratinocytes, mast cells, NGF-TrkA)
4. Innate immunity and cytokines (IL-1β/IL-6/TNF-α/IL-8, IL-10↓, monocytes, NF-κB, complement)
5. Autoimmunity (β2AR/M2/AT1R autoantibodies, passive transfer, IVIG, rituximab)
6. Oxidative stress (ROS, •OH, ONOO⁻, lipid peroxidation, antioxidant defence↓, mitochondria)
7. Sympathetic-afferent coupling (α1 up-regulation, sympathetic sprouting, DRG, sympathetic block)
8. Microvasculature and hypoxia (endothelial dysfunction, ET-1↑/NO↓, AV shunt, lactate, pH↓, warm↔cold)
9. Peripheral sensitisation (TRPV1/TRPA1, Nav1.7/1.8, ASIC3, P2X3, allodynia)
10. Spinal central sensitisation (glutamate, NMDA-NR2B, wind-up, WDR expansion, disinhibition)
11. Spinal neuroimmunity (microglial P2X4, BDNF-TrkB-KCC2, astrocytes, CCL2)
12. Descending pain modulation (PAG-RVM on/off cells, LC-NE, 5-HT, loss of CPM/DNIC)
13. Cortical reorganisation (S1 map shrinkage, two-point discrimination threshold↑, body-perception disturbance, neglect-like features, GMI/mirror therapy)
14. Bone remodelling (RANKL/OPG, CTX-I, patchy osteopenia, bone scan, bisphosphonates)
15. Motor and trophic changes (dystonia, tremor, ROM↓, contracture, BoNT-A, intrathecal baclofen)
16. The psycho-behavioural disuse loop (catastrophising, kinesiophobia TSK, avoidance, CBT and graded exposure)
17. Drug PK (ketamine/norketamine, prednisolone, neridronate, gabapentin, amitriptyline, NAC, IVIG, lidocaine, vitamin C)
18. Drug PD and targets (NMDA blockade, NF-κB inhibition, FPPS inhibition, α2δ, NE/5-HT reuptake, antioxidant, FcRn)
19. Interventions and neuromodulation (SCS, DRG stimulation, sympathetic block, intrathecal infusion, multidisciplinary rehabilitation)
20. Clinical endpoints and natural history (NRS, CSS 0-16, Budapest, ROM, temperature asymmetry, remission/chronicity/the therapeutic window)

![CRPS QSP map](crps_qsp_model.png)

## 4. The mrgsolve model — 34 ODEs

| Compartment group | Compartments |
|---|---|
| PK (14) | `KET_C1` `KET_C2` `NORKET` / `PRED_GUT` `PRED_C` / `NER_C` `NER_BONE` / `GBP_GUT` `GBP_C` / `AMT_GUT` `AMT_C` / `NAC_GUT` `NAC_C` / `IVIG_C` |
| Peripheral node (8) | `NP` `CYT` `NGF` `EDEMA` `ROS` `AAB` `ALPHA1` `PSENS` |
| Vascular and hypoxia (2) | `PERF` `HYPOX` |
| Central loops (5) | `SSENS` `GLIA` `DINH` `CORTEX` `DISUSE` |
| Endpoints and bone (5) | `PAIN` `ROM` `OC` `BMD` `CTXI` |

Derived outputs: `CSS` (0-16, the four Budapest domains), `TEMP_ASYM` (°C),
`ACTIVE_CRPS`, `REMISSION`, `PHENOTYPE` (±1 = warm/cold), `LATCHED` (whether the glial
latch has closed), `RING_GAIN` (the local gain of the pain→disuse→cortex→pain loop).

**Nine scenarios**: ① untreated natural history ② early (day 7) multidisciplinary
package ③ the same package late (day 240) ④ rehabilitation alone (day 30)
⑤ prednisolone alone (day 7) ⑥ neridronate 100 mg IV ×4 ⑦ a 100-hour ketamine
infusion ⑧ ketamine + rehabilitation ⑨ the cold, long-standing phenotype + SCS.

**Ten analysis functions**: `run_scenarios` `CRPS_trait_bifurcation` `CRPS_insult_scan`
`CRPS_window_scan` `CRPS_arm_decomposition` `CRPS_ketamine_washout`
`CRPS_dose_vs_timing` `CRPS_phenotype_ordering` `CRPS_bone_axis`
`CRPS_scs_habituation` (plus `CRPS_report`, which runs them all).

## 5. Computed results

> Every number below is **a value printed by the code in this repository**, not a
> literature value. It reproduces exactly with `crps_reference_check.py` (which needs
> only numpy/scipy).
> The literature was used only to anchor parameters, and the correspondence is set out
> in the last table of [`crps_references.md`](crps_references.md).

### 5.1 Natural history — a two-step escalation

| Day | NRS | CSS | Temp. asymmetry | PSENS | SSENS | GLIA | CORTEX | DISUSE | ROM | BMD |
|---|---|---|---|---|---|---|---|---|---|---|
| 3 | 2.55 | 4.25 | **+1.64** | 0.17 | 0.10 | 0.00 | 0.00 | 0.10 | 0.93 | 1.00 |
| 14 | 5.39 | 5.43 | +0.55 | 0.22 | 0.30 | 0.09 | 0.19 | 0.49 | 0.53 | 0.99 |
| 30 | 6.05 | 5.39 | −0.09 | 0.17 | 0.30 | 0.20 | 0.45 | 0.59 | 0.38 | 0.97 |
| 90 | 5.91 | 5.16 | −0.38 | 0.10 | 0.26 | 0.31 | 0.52 | 0.51 | 0.44 | 0.93 |
| 180 | **6.92** | 6.81 | **−1.09** | 0.12 | **0.71** | **0.57** | 0.60 | 0.55 | 0.37 | 0.90 |
| 365 | 6.92 | 6.81 | −1.09 | 0.12 | 0.71 | 0.59 | 0.60 | 0.55 | 0.37 | 0.87 |

The peripheral node extinguishes after day 30 (PSENS 0.17→0.12) and yet the pain rises
**a second time between days 100 and 150**. This was not put into the parameters; it is
an emergent result that occurs at the point where GLIA passes its threshold. The
transition from warm (+1.64 °C) to cold (−1.09 °C) comes out the same way.

### 5.2 What determines chronicity — trait versus the size of the injury

| KFEAR (fear-avoidance gain) | 0.20 | 0.50 | 0.60 | 0.80 | **0.90** | 1.00 | 1.40 |
|---|---|---|---|---|---|---|---|
| NRS at 3 years | 0.00 | **0.00** | **4.43** | 5.20 | **6.89** | 6.92 | 6.99 |

| INJ_AMP (strength of the inciting event) | 0.10 | 0.35 | 0.70 | 0.90 | **0.95** | 1.20 | 1.50 |
|---|---|---|---|---|---|---|---|
| NRS at 3 years | 5.50 | 5.50 | 5.50 | 5.50 | **6.92** | 6.92 | 6.92 |

**The size of the injury does not determine whether the condition becomes chronic** —
from 0.10 to 0.90 the NRS at 3 years is 5.50 throughout, and at 0.95 it steps up to
6.92 (the selection of a severity tier). In the psychological trait parameter KFEAR, by
contrast, **resolution versus chronicity** is decided between 0.50 and 0.60, and
**attractor 1 versus 2** between 0.80 and 0.90. The same fracture, the same
pharmacology, a different outcome.

### 5.3 The therapeutic window — the same prescription, only the start day different

A six-month package of a prednisolone taper + NAC + rehabilitation:

| Start day | 3 | 30 | 60 | 75 | **90** | **95** | 120 | 240 | 365 |
|---|---|---|---|---|---|---|---|---|---|
| NRS at 3 years | 0.00 | 0.00 | 0.00 | 0.00 | **0.00** | **6.92** | 6.92 | 6.92 | 6.92 |
| BMD at 3 years | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 0.85 | 0.85 | 0.85 | 0.85 |

**t\* = 90-95 days.** Not one character of the prescription, the dose or the duration
differs. All that differs is whether GLIA had already passed its threshold by the time
the package arrived.

### 5.4 Which arm creates the window (NRS at 3 years)

| Start day | Prednisolone alone | NAC alone | **Rehabilitation alone** | Pred+NAC | Everything |
|---|---|---|---|---|---|
| d7 | 5.50 | 5.50 | **0.00** | 5.50 | 0.00 |
| d30 | 5.50 | 5.50 | **0.00** | 5.50 | 0.00 |
| d60 | 6.92 | 6.92 | **0.00** | 6.92 | 0.00 |
| d120 | 6.92 | 6.92 | 6.92 | 6.92 | 6.92 |
| d240 | 6.92 | 6.92 | 6.92 | 6.92 | 6.92 |

The only arm that reaches resolution is **rehabilitation, the arm that targets the
loop**. The whole of the measurable benefit of the anti-inflammatory arms is
"**prevention of the late transition (attractor 1→2)**" (5.50 vs 6.92), and it exists
only if given within 30 days. The implication for trial design is direct: **enrol
together the patients who were never going to undergo the transition and those who
already have, and the result reads as null however powerful the drug.**

Lower the rehabilitation intensity (adherence) and the window shortens — at an
intensity of 0.25-0.75 it is 60 days, and only at an intensity of 1.0 does resolution
hold out to 90 days.

### 5.5 Ketamine — the dissociation of acute analgesia from disease modification

| Day | No treatment | Ketamine | Rehabilitation | Ketamine + rehabilitation |
|---|---|---|---|---|
| 29 (before the infusion) | 6.04 | 6.04 | 6.04 | 6.04 |
| 32 (during the infusion) | 6.04 | **3.68** | 4.82 | **2.90** |
| 34 (end of the infusion) | 6.03 | **3.02** | 4.05 | **1.89** |
| 37 (after washout) | 6.02 | 5.78 | 3.47 | 3.08 |
| 44 | 5.98 | 5.86 | 2.85 | 2.67 |
| 210 | 6.92 | 5.51 | 0.00 | 0.00 |
| 700 | 6.92 | **5.50** | **0.00** | **0.00** |

Ketamine Cmax 258 ng/mL, norketamine 439 ng/mL. During the infusion 6.04→3.02
(corresponding to Sigtermans 2009's 7.2→2.7), and it returns to where it was three days
after stopping. All that remains in the long run is **a drop of one tier
(6.92→5.50)**.

### 5.6 A result contrary to the hypothesis this model was written to express

The model was written to embody the hypothesis that "rehabilitation works inside the
analgesic window ketamine creates, and so the combination is **supra-additive**".
**It was not.** At day 700 the interaction term is **−1.42 NRS (sub-additive)** —
because rehabilitation alone has already reached resolution, so there is nothing the
drug can add at the endpoint. Early on, by contrast, there is a clear additional
benefit (day 34: 1.89 vs 4.05). That is, **the drug does not raise the ceiling; it buys
time.**

The dose axis points the same way (NRS at 2 years):

| Maximum ketamine rate | d14 | d30 | d60 | d120 | d240 | Cmax |
|---|---|---|---|---|---|---|
| 5.5 mg/h | 6.92 | 6.92 | 6.92 | 6.92 | 6.92 | 65 ng/mL |
| 11 mg/h | 6.92 | 5.50 | 5.50 | 6.92 | 6.92 | 129 |
| 22 mg/h | 5.50 | 5.50 | 5.50 | 6.92 | 6.92 | 258 |
| 44 mg/h | 5.50 | 5.50 | 5.50 | 6.92 | 6.92 | 517 |
| **88 mg/h** | **5.50** | 5.50 | 5.50 | 6.92 | 6.92 | **1034** |

Above the threshold (about 11-22 mg/h), **quadrupling the dose changes the NRS at 2
years by exactly 0.00**. At the same maximum dose, moving the start day from 60 to 120
changes it by 1.42. Being a stronger analgesic is not the same as being a stronger
disease modifier.

### 5.7 The ranking of the arms by phenotype and stage (mean improvement in NRS over the 180 days after treatment)

| Warm, early (d14) | Δ mean NRS | | Cold, long-standing (d240) | Δ mean NRS |
|---|---|---|---|---|
| Rehabilitation/GMI | **+4.68** | | SCS | **+2.91** |
| SCS | +4.39 | | Rehabilitation/GMI | +1.49 (benefit at 2 years 0.00) |
| Gabapentin | +0.60 | | Gabapentin | +0.50 |
| Amitriptyline | +0.41 | | Amitriptyline | +0.34 |
| Bisphosphonate | +0.37 | | Ketamine | +0.07 |
| Vasodilator | +0.10 | | Vasodilator, bisphosphonate | +0.03 |
| Antioxidant, steroid, IVIG, sympathetic block | ≤ +0.07 | | Steroid, antioxidant, IVIG, sympathetic block | ≈ 0.00 |

The ranking is a property not of the molecules but **of the state of the system**. In
the cold, long-standing phenotype every drug arm collapses to a symptomatic effect, the
only thing that changes the 2-year endpoint is a device (SCS), and even that
habituates.

### 5.8 The bone axis (neridronate 100 mg IV ×4, started on day 30)

| Day | NRS untreated | NRS neridronate | CTX-I untreated | CTX-I treated | BMD untreated | BMD treated |
|---|---|---|---|---|---|---|
| 40 | 6.00 | 5.82 | 0.434 | 0.239 | 0.959 | 0.968 |
| 120 | 6.73 | 5.18 | 0.409 | 0.159 | 0.922 | 0.978 |
| 365 | 6.92 | **5.09** | 0.412 | **0.161 (−61%)** | 0.867 | **0.994** |

In this model the analgesia of a bisphosphonate comes not from the bone endpoint but
from cutting **the intraosseous acidosis pathway** (`W_BONE_PS`). It therefore shows
the same time dependence as the other peripheral-node drugs, which fits the fact that
the actual RCT was positive in patients with a short history.

### 5.9 SCS habituation (implanted at 6 months)

| Month | 7 | 12 | 24 | 36 | 60 |
|---|---|---|---|---|---|
| Δ NRS (control − SCS) | **+3.27** | +2.44 | +1.20 | +0.60 | **+0.15** |

With the same device and the same settings, a large initial difference decays
exponentially. It reproduces exactly the shape of Kemler's positive result at 6 months
(NEJM 2000) and the loss of the between-group difference at 5 years (J Neurosurg 2008),
while the underlying attractor does not change at all in the interval.

## 6. How to run it

```r
# --- the mrgsolve model ---
install.packages(c("mrgsolve", "dplyr"))
library(mrgsolve)
mod <- mread("crps_mrgsolve_model.R")
e   <- mod@envir

e$run_scenarios(mod)             # summary table of the 9 scenarios
e$CRPS_trait_bifurcation(mod)    # the table in 5.2
e$CRPS_window_scan(mod)          # the table in 5.3
e$CRPS_arm_decomposition(mod)    # the table in 5.4
e$CRPS_ketamine_washout(mod)     # the tables in 5.5 / 5.6
e$CRPS_dose_vs_timing(mod)       # the dose-timing grid in 5.6
e$CRPS_phenotype_ordering(mod)   # the table in 5.7
e$CRPS_bone_axis(mod)            # the table in 5.8
e$CRPS_scs_habituation(mod)      # the table in 5.9
e$CRPS_report(mod)               # all of it

# --- the Shiny dashboard (8 tabs) ---
install.packages(c("shiny", "ggplot2", "tidyr"))
shiny::runApp("crps_shiny_app.R")
```

```bash
# --- reproduce the numbers alone, without R/mrgsolve (numpy + scipy) ---
python3 crps_reference_check.py            # everything
python3 crps_reference_check.py win ket    # a selection

# --- rendering the map ---
dot -Tsvg crps_qsp_model.dot -o crps_qsp_model.svg
dot -Tpng -Gdpi=150 crps_qsp_model.dot -o crps_qsp_model.png
```

The first thing to try in the Shiny app: leave the treatment arms switched on and move
**only the "treatment start day" slider**, watching whether the GLIA curve on tab 4
crosses the latch threshold line. The result of 5.3 reproduces itself on screen
immediately.

## 7. Limitations

- **This is a semi-quantitative model for teaching and research** and it has not been
  fitted to patient data. The parameters were anchored to the literature and then
  adjusted to satisfy the target behaviours of §5 (incidence, phenotype switching, and
  the direction and magnitude of the major RCTs). It cannot be used for clinical
  decision-making.
- `REHAB = 1` is **the ideal upper bound of a fully adhered six-month
  multidisciplinary programme**. The "complete resolution" of §5.4 is the value at that
  upper bound, and a gradient with reduced adherence is presented alongside it.
- The existence of the three attractors is robust, being a consequence of the topology
  (positive feedback plus threshold non-linearity), but **the exact positions of the
  thresholds (GLIA50 = 0.36, KFEAR ≈ 0.55, t\* ≈ 90 days) are calibrated values** and
  are not themselves measured clinical quantities. What is offered as a testable
  prediction is not the positions but the **structure**: (i) the effect of the
  anti-inflammatory arms should appear only as prevention of the late transition,
  (ii) the length of the window should move with rehabilitation intensity, and
  (iii) the dose-response of an NMDA blocker should be flat above the threshold.
- Inter-individual variability (IIV/ω), adverse drug reactions (steroid metabolic
  effects, the psychotomimetic and dissociative effects of ketamine, the bisphosphonate
  acute-phase reaction, SCS lead migration/infection), paediatric CRPS, and the
  explicit dynamics of contralateral spread and of dystonia are not included.
