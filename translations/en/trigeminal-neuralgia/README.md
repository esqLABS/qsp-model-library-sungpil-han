# Trigeminal Neuralgia (TN) — QSP Model

> Integrated Quantitative Systems Pharmacology model linking neurovascular
> compression (NVC) of the trigeminal root entry zone (REZ) to focal
> demyelination, voltage-gated sodium channel (Nav1.3/1.6/1.7)
> upregulation, ephaptic crosstalk, ectopic/afterdischarge generation, and
> central trigeminal-nucleus sensitization that together produce
> paroxysmal lancinating facial pain — coupled to anticonvulsant PK/PD
> (carbamazepine with autoinduction, oxcarbazepine prodrug/MHD, baclofen,
> gabapentin, pregabalin) and interventional therapy (microvascular
> decompression, percutaneous radiofrequency rhizotomy) with time-dependent
> recurrence.

| Deliverable | File |
|---|---|
| Mechanistic map (DOT) | [`tn_qsp_model.dot`](tn_qsp_model.dot) |
| Map (SVG)             | [`tn_qsp_model.svg`](tn_qsp_model.svg) |
| Map (PNG, 150 dpi)    | [`tn_qsp_model.png`](tn_qsp_model.png) |
| mrgsolve ODE model     | [`tn_mrgsolve_model.R`](../../../trigeminal-neuralgia/tn_mrgsolve_model.R) |
| Shiny dashboard        | [`tn_shiny_app.R`](tn_shiny_app.R) |
| References             | [`tn_references.md`](tn_references.md) |

---

## 1. Disease in One Paragraph

The core pathogenesis of trigeminal neuralgia is chronic compression of
the trigeminal root entry zone (REZ) by a pulsatile artery (usually the
superior cerebellar artery) or vein, which causes focal demyelination.
Abnormal expression and redistribution of Nav1.7/1.3/1.6 sodium channels
occurs on the exposed naked axon, generating ectopic spontaneous firing
and afterdischarge (Devor's "ignition hypothesis"). Ephaptic crosstalk
between adjacent axons at the demyelinated site cross-wires light-touch
(Aβ) signals into the pain (Aδ/C) pathway, so an innocuous stimulus such
as lightly touching the face, chewing, or brushing teeth triggers an
electric-shock-like pain paroxysm. Repetitive afferent input produces
NMDA-mediated wind-up and central sensitisation in the spinal trigeminal
nucleus caudalis (Vc), amplifying the intensity, persistence, and
allodynia of the pain. First-line pharmacotherapy is with the sodium
channel blockers carbamazepine and oxcarbazepine; when response is
insufficient, baclofen (GABA-B) or gabapentin/pregabalin (α2δ) combination
therapy is added. Patients refractory to medication are transitioned to
destructive or non-destructive procedures such as microvascular
decompression (MVD, which removes the causative lesion), Gamma Knife,
percutaneous radiofrequency rhizotomy, or balloon compression, and even
after a procedure recurrence can occur over time.

## 2. Mechanistic Map Clusters (19 clusters, 108 nodes)

1. Aetiology/risk factors (advanced age, hypertension, vascular dilation, genetic predisposition, MS/tumour/AVM comorbidity)
2. Neurovascular compression (SCA/AICA/venous compression, REZ, pulsatile compression)
3. Focal demyelination (Schwann cell injury, naked axon, MRI findings)
4. Nav channelopathy (Nav1.7/1.3/1.6/1.8, ectopic discharge, afterdischarge, lowered threshold)
5. Peripheral crosstalk/sensitisation (ephaptic crosstalk, Aβ→Aδ/C miswiring, local neuroinflammation)
6. Central pathway/transmission (trigeminal ganglion → Vc → trigeminothalamic tract → VPM → somatosensory cortex, descending modulation)
7. Central sensitisation (NMDA, wind-up, glial activation, allodynia)
8. Secondary TN aetiology (MS plaques, cerebellopontine angle tumours, bilateral/sensory-loss red flags)
9. Carbamazepine PK (autoinduction, active metabolite epoxide, HLA-B*1502)
10. Oxcarbazepine PK (prodrug → MHD, hyponatraemia)
11. Adjunct/second-line drug PK (baclofen, lamotrigine, gabapentin, pregabalin, BoNT-A)
12. Drug PD (Nav blockade, reduced synaptic glutamate, α2δ/GABA-B, SNARE cleavage)
13. Procedures/surgery (MVD, Gamma Knife, RF rhizotomy, balloon compression, glycerol)
14. Adverse effects/safety (drowsiness, hepatotoxicity, hyponatraemia, SJS/TEN, agranulocytosis, post-procedure anaesthesia dolorosa)
15. Clinical endpoints (attack frequency, NRS, BNI, complete remission rate, time to relapse, QoL)
16. Natural history (relapsing-remitting course, progressive worsening, atypical TN)
17. Trigeminal branches/clinical presentation (V1/V2/V3, ICHD-3, high-resolution MRI, epidemiology)
18. Mental health/functional comorbidity (depression, anxiety, sleep disturbance, functional impairment, caregiver burden)
19. Second-line pharmacotherapy (topiramate, phenytoin, levetiracetam, sumatriptan)

## 3. mrgsolve Model (17 ODE Compartments)

* **Drug PK (5 drugs, 12 compartments)** — Carbamazepine (gut/central/epoxide/enzyme-induction-state, 4 compartments,
  with a 1-compartment turnover model capturing autoinduction), oxcarbazepine (gut/MHD, 2 compartments),
  baclofen (gut/central, 2 compartments), gabapentin (saturable-absorption gut/central, 2 compartments),
  pregabalin (linear gut/central, 2 compartments).
* **Disease/PD (5 compartments)** — Nav channel upregulation index (NAV_UPREG, residual NVC drive),
  ectopic discharge (ECTOPIC), central sensitisation index (CENTSENS), attack frequency (PAROX, per day),
  pain intensity (PAIN, NRS 0-10).
* **Safety/procedure (5 compartments)** — plasma Na+ (hyponatraemia), drowsiness/ataxia score,
  MVD_STATE and RF_STATE (post-procedure remission-to-relapse state variables).
* Nav blockade combines an Emax model of CBZ (parent drug + epoxide combined exposure) and OXC (MHD)
  (`1-(1-block_cbz)(1-block_oxc)`), and central-sensitisation inhibition is modelled as the additive
  inhibition of GABA-B (baclofen) and α2δ (gabapentin + pregabalin) mechanisms.

### 7 Scenarios

| # | Scenario | Calibration basis |
|---|---|---|
| 1 | Untreated natural history | Maarbjerg 2014 Headache natural-history cohort |
| 2 | Carbamazepine monotherapy (200 mg BID → titrated to 400 mg TID) | Zakrzewska 1989 JNNP crossover study, Wiffen 2014 Cochrane |
| 3 | Oxcarbazepine monotherapy (300 mg BID → 1200 mg/day) | Zakrzewska 1989 JNNP, Besi 2015 |
| 4 | Carbamazepine + baclofen combination (refractory) | Fromm 1984 Ann Neurol RCT |
| 5 | MVD procedure (day 14, post-operative CBZ taper) | Barker 1996 NEJM, Sindou 2006 |
| 6 | CBZ intolerance → switch to gabapentin + pregabalin | Al-Quliti 2015 review |
| 7 | Percutaneous radiofrequency rhizotomy (RF, day 30), relapse follow-up | Tronnier 2001 Neurosurgery |

## 4. Shiny Dashboard (8 Tabs)

1. **Patient profile** — adjusts NVC severity, branch involvement, secondary (MS) status.
2. **PK** — blood concentrations for each drug (CBZ/epoxide/OXC-MHD/baclofen/gabapentin/pregabalin).
3. **Pathway PD** — Nav channel block fraction, central sensitisation index.
4. **Clinical endpoints** — attack frequency, pain NRS, approximate BNI grade.
5. **Scenario comparison** — overlaid comparison of multiple scenarios plus a summary table.
6. **Biomarkers** — Nav channel expression/ectopic-excitability index (linked conceptually to imaging findings).
7. **Safety** — plasma Na+ (hyponatraemia warning line), drowsiness/ataxia score.
8. **References** — the full reference list.

## 5. How to Run

```bash
# 1) Render the mechanistic map
dot -Tsvg tn_qsp_model.dot -o tn_qsp_model.svg
dot -Tpng -Gdpi=150 tn_qsp_model.dot -o tn_qsp_model.png
```

```r
# 2) R/mrgsolve simulation
install.packages(c("mrgsolve","dplyr","tidyr","ggplot2","shiny","DT"))
library(mrgsolve)
mod <- mread("tn_mrgsolve_model.R")
source("tn_mrgsolve_model.R")  # to load the run_scenarios() helper (optional)
results <- run_scenarios(mod)
plot(results$cbz_mono %>% mrgsolve::filter_sims(time <= 720), c("PAROX","PAIN","CBZ_conc"))

# 3) Run the Shiny dashboard
shiny::runApp("tn_shiny_app.R")
```

## 6. Key Clinical Calibration Basis

| Endpoint | Compared against | Basis |
|---|---|---|
| CBZ autoinduction (t1/2 36h → 12-17h) | Increased clearance after repeated dosing | Bertilsson & Tomson 1986 Clin Pharmacokinet |
| OXC hyponatraemia incidence ~2.7% | Relative risk versus CBZ | Dong 2005 Neurology |
| Baclofen combination effect | Attack reduction versus CBZ alone | Fromm 1984 Ann Neurol |
| MVD long-term pain-free survival | ~70% pain-free at 10 years, ~1-4% annual recurrence | Barker 1996 NEJM; Sindou 2006 |
| RF rhizotomy/Gamma Knife recurrence | Tendency to recur faster than MVD | Tronnier 2001; Kondziolka 1996 |
| HLA-B*1502 SJS/TEN risk | Testing recommended before starting CBZ in Asians | Chung 2004 Nature; Ferrell 2008 |

## 7. Model Verification Status

Because this container does not have an R/mrgsolve execution environment
installed (no `Rscript`), the mrgsolve model has been completed through the
stage of **literature-based parameter design and self-review of the code
(dimensional and boundary-value checks)**, but its numerical output has not
been verified by actually compiling and integrating it. The `.dot` file
was rendered with Graphviz `dot`, and the resulting SVG/PNG were actually
generated and checked. Where an mrgsolve/R environment is available, it is
recommended to run it per the "How to Run" section above and confirm that
trajectories such as `PAROX`, `PAIN`, and `NA_PLASMA` move in the expected
direction (reduced attack frequency and pain under treatment, reduced
hyponatraemia with OXC, etc.).

## 8. Caveats

* For research, education, and hypothesis generation only; must not be
  used for clinical decision-making.
* The EC50/Emax values for Nav channel blockade and central-sensitisation
  inhibition are surrogate parameters not directly estimated from the
  literature, intended to show relative direction/ranking rather than
  absolute magnitude.
* MVD/RF recurrence rates are an approximation that simplifies the
  literature's survival curves to a single exponential function.
* This is a typical-value model that does not include inter-individual
  variability (IIV).
