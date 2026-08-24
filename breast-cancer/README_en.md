# Breast Cancer — Quantitative Systems Pharmacology Model

---

## Disease Overview

Breast cancer is the most common malignancy in women worldwide, with an estimated 2.3 million new patients diagnosed annually as of 2020, accounting for 11.7% of all cancer diagnoses. It is broadly classified into four molecular subtypes: estrogen-receptor-positive (ER+)/HER2-negative, which accounts for approximately 65–70% of all cases; HER2-positive (HER2+), 15–20%; and triple-negative breast cancer (TNBC), approximately 15%. Each subtype exhibits distinct pathophysiological mechanisms, treatment responsiveness, and prognosis, making subtype-tailored treatment strategies essential. In Korea, as in Western countries, the incidence of breast cancer has been steadily rising, with prevalence particularly high among women in their 40s and 50s. The HR+ subtype is driven by estrogen receptor (ER) signalling, HER2+ by overexpression/amplification of the ERBB2 receptor tyrosine kinase, and TNBC lacks these established targets and is associated with the poorest prognosis.

The treatment landscape for breast cancer has been transformed over the past two decades and differs substantially by molecular subtype. In ER+ disease, standard treatment is endocrine therapy based on an aromatase inhibitor (AI) or tamoxifen, and combining this with a CDK4/6 inhibitor (palbociclib, ribociclib, abemaciclib) has become first-line standard of care in metastatic settings, significantly extending progression-free survival (mPFS 25–28 months vs. ~14–16 months for endocrine monotherapy). HER2-positive disease benefits from highly effective anti-HER2 therapies, including dual blockade with trastuzumab and pertuzumab and next-generation antibody-drug conjugates (ADCs) such as T-DM1 and T-DXd, which have substantially improved survival. Targeted agents for TNBC remain more limited, but treatment options are expanding with the emergence of immune checkpoint inhibitors (pembrolizumab), PARP inhibitors (olaparib), and ADCs. Despite these advances, resistance remains a major clinical challenge: ESR1 mutations develop under aromatase-inhibitor pressure, CDK4/6-inhibitor resistance emerges via RB loss and cyclin E amplification, and HER2+ resistance is mediated by PI3K/AKT pathway activation. The QSP model developed here provides a mechanistic framework to simulate these dynamics, quantify treatment effects, and explore combination strategies.

---

## Model Structure

| Component | Details |
|-----------|---------|
| Mechanistic Map Nodes | 110+ nodes, 9 subgraph clusters |
| ODE Compartments | 22 compartments |
| Drug Coverage | 14 drugs (6 drug classes) |
| Treatment Scenarios | 6 clinical regimens |
| Simulation Duration | Up to 1 year (8,760 hours) |
| Key Pathways Modeled | ER/PI3K/AKT/mTOR, HER2/MAPK, CDK4/6-RB, DNA damage/PARP, PD-1/PD-L1 |
| Biomarkers Simulated | Ki-67, CA15-3, RB phosphorylation, CD8+ TILs, PD-L1, ctDNA (ESR1) |
| Resistance Modules | ESR1 mutation, RB loss, PI3K activation, CDK4/6i bypass |

---

## Drug Coverage

| Drug | Class | Target | Key Trial |
|------|-------|--------|-----------|
| Palbociclib | CDK4/6 inhibitor | CDK4/CDK6 | PALOMA-2 |
| Ribociclib | CDK4/6 inhibitor | CDK4/CDK6 | MONALEESA-2 |
| Abemaciclib | CDK4/6 inhibitor | CDK4/CDK6 | MONARCH-3 |
| Letrozole | Aromatase inhibitor | CYP19A1 | PALOMA-2 |
| Anastrozole | Aromatase inhibitor | CYP19A1 | ATAC |
| Tamoxifen | SERM | ERα | EBCTCG meta-analysis |
| Trastuzumab | Anti-HER2 antibody | HER2/ERBB2 | CLEOPATRA |
| Pertuzumab | Anti-HER2 antibody | HER2 dimerization domain | CLEOPATRA |
| T-DM1 | ADC (HER2-targeted) | HER2 + tubulin | EMILIA |
| T-DXd | ADC (HER2-targeted) | HER2 + Topo I | DESTINY-Breast03 |
| Pembrolizumab | Anti-PD-1 mAb | PD-1 | KEYNOTE-522 |
| Olaparib | PARP inhibitor | PARP1/2 | OlympiAD |
| Alpelisib | PI3K inhibitor | PIK3Cα (PI3Kα) | SOLAR-1 |
| Everolimus | mTOR inhibitor | mTORC1 | BOLERO-2 |

---

## Key Clinical Trials Simulated

| Trial | Regimen | Subtype | Key Result | Model Calibration Target |
|-------|---------|---------|------------|--------------------------|
| PALOMA-2 | Palbociclib + Letrozole | ER+/HER2- | mPFS 27.6 vs 14.5 mo | CDK4/6 inhibition depth, Ki-67 ≥50% suppression |
| MONALEESA-2 | Ribociclib + Letrozole | ER+/HER2- | mPFS 25.3 vs 16.0 mo | CDK4/6 + ER signaling synergy |
| MONARCH-3 | Abemaciclib + AI | ER+/HER2- | mPFS 28.2 vs 14.8 mo | CDK4/6 sustained inhibition (continuous dosing) |
| CLEOPATRA | Trastuzumab + Pertuzumab + Docetaxel | HER2+ | mPFS 18.7 mo, HR=0.62 | HER2 dimerization blockade, MAPK suppression |
| KEYNOTE-522 | Pembrolizumab + Chemo (neoadjuvant) | TNBC | pCR 64.8% vs 51.2% | CD8 expansion, PD-1 checkpoint release |
| OlympiAD | Olaparib monotherapy | BRCAm, HER2- | mPFS 7.0 vs 4.2 mo | PARP inhibition, synthetic lethality in BRCA-null |

---

## How to Run

### Prerequisites

```r
install.packages(c("mrgsolve", "ggplot2", "dplyr", "tidyr", "shiny",
                   "shinydashboard", "plotly"))
```

### mrgsolve ODE Model

```r
source("bc_mrgsolve_model.R")
# Simulates all 6 treatment scenarios
# Outputs: tumor volume, Ki-67, PK profiles, biomarker dynamics
```

### Shiny Interactive App

```r
shiny::runApp("bc_shiny_app.R")
# Opens interactive dashboard with 6 tabs:
# Tab 1: Patient Profile | Tab 2: PK | Tab 3: PD | Tab 4: Clinical Endpoints
# Tab 5: Scenario Comparison | Tab 6: Biomarkers
```

### Generate Mechanistic Map

```bash
# SVG (vector, for web/report embedding)
dot -Tsvg bc_qsp_model.dot -o bc_qsp_model.svg

# PNG (150 dpi raster, for README thumbnail)
dot -Tpng -Gdpi=150 bc_qsp_model.dot -o bc_qsp_model.png
```

---

## File Descriptions

| File | Description |
|------|-------------|
| [`bc_qsp_model.dot`](bc_qsp_model.dot) | Graphviz DOT mechanistic map (110+ nodes, 9 clusters) |
| [`bc_qsp_model.svg`](bc_qsp_model.svg) | Rendered SVG of mechanistic map (vector format) |
| [`bc_qsp_model.png`](bc_qsp_model.png) | Rendered PNG of mechanistic map (150 dpi) |
| [`bc_mrgsolve_model.R`](bc_mrgsolve_model.R) | mrgsolve ODE model with 22 compartments and 6 treatment scenarios |
| [`bc_shiny_app.R`](bc_shiny_app.R) | Interactive Shiny dashboard with 6 tabs |
| [`bc_references_en.md`](bc_references_en.md) | 50 PubMed-linked references organized by 9 sections |
| `README.md` | This file — model overview and usage guide |

---

## Key Parameters

| Parameter | Symbol | Value | Unit | Source |
|-----------|--------|-------|------|--------|
| Tumor proliferation rate | kprol | 0.0008 | 1/hr | Calibrated to PALOMA-2 control arm |
| Tumor death rate | kdeath | 0.0002 | 1/hr | Literature (Simeoni model) |
| Tumor carrying capacity | Kmax | 1000 | relative units | Estimated |
| E2 EC50 (proliferation stimulation) | EC50_E2 | 50 | pmol/L | ER biology literature |
| Letrozole Emax (E2 suppression) | Emax_AI | 0.98 | fraction | ATAC/BIG 1-98 PD data |
| CDK4/6 inhibitor Emax | Emax_CDK | 0.85 | fraction | PALOMA-2 Ki-67 PD data |
| CDK4/6 inhibitor EC50 | EC50_CDK | 100 | ng/mL | Population PK/PD modeling |
| Palbociclib clearance | CL_palbo | 63 | L/hr | Population PK (Friberg) |
| Palbociclib volume of distribution | Vd_palbo | 2583 | L | Population PK (Friberg) |
| Trastuzumab clearance | CL_tras | 0.225 | L/day | Population PK (Bruno) |
| Trastuzumab Vd (central) | Vc_tras | 3.1 | L | Population PK (Bruno) |
| HER2 signaling Emax | Emax_HER2 | 0.70 | fraction | HER2 CLEOPATRA PD data |
| HER2 EC50 (trastuzumab) | EC50_HER2 | 50 | mg/L | PK/PD estimate |
| CD8+ T cell kill rate | kCD8_kill | 0.002 | 1/hr | Immunology model (Jiang) |
| PD-1 inhibition Emax | Emax_PD1 | 0.80 | fraction | KEYNOTE-522 pCR calibration |
| Olaparib clearance | CL_olap | 8.6 | L/hr | Population PK (Karlsson) |
| Olaparib Emax (PARP inhibition) | Emax_PARP | 0.95 | fraction | OlympiAD biomarker data |
| ESR1 mutation rate | mu_ESR1 | 1e-6 | 1/cell/hr | Estimated from ctDNA kinetics |
| RB loss rate (CDK4/6i resistance) | mu_RBlos | 5e-7 | 1/cell/hr | Turner 2019 (PALOMA-3 resistance) |
| PI3K activation rate (resistance) | kPI3K | 0.0003 | 1/hr | Juric 2015 (SOLAR-1) |

---

## Expected Simulation Outputs

1. **Tumor Volume Trajectories** — All 6 regimens plotted as % change from baseline over 52 weeks, with uncertainty bands calibrated to trial median PFS times.
2. **Ki-67 Proliferation Index** — Dynamic change from baseline (%) at C1D14 and C2D1, calibrated to PALOMA-2 and MONARCH-3 biopsy data.
3. **CDK4/6 Inhibition Depth** — RB phosphorylation inhibition (%) over the dosing cycle; sustained inhibition shown for abemaciclib (continuous) vs. palbociclib/ribociclib (intermittent).
4. **CA15-3 Tumor Marker** — Serum biomarker response kinetics correlating with tumor burden; waterfall plots at Week 12.
5. **PK Profiles** — Drug concentration-time curves showing Cmax, Cmin, and AUC for each agent; steady-state achieved within 5–8 days for oral CDK4/6i.
6. **Immune Cell Dynamics** — CD8+ effector T cell expansion, regulatory T cell (Treg) suppression, and PD-L1 upregulation trajectories for TNBC/pembrolizumab simulation.
7. **PFS Curves** — Kaplan-Meier-like time-to-event curves for all scenarios, reproduced from model output; hazard ratios vs. control computed from tumor growth inhibition model.

---

## Limitations & Assumptions

- **Spatial homogeneity**: Tumor spatial heterogeneity is not modeled; the well-mixed compartment assumption is used throughout (no intratumoral gradients in drug penetration or oxygen).
- **Binary resistance states**: Acquired resistance is modeled as deterministic transitions between sensitive and resistant subpopulations, not as continuous stochastic evolutionary processes.
- **Simplified immune microenvironment**: The immune module uses a 3-cell-type model (CD8+ effectors, Tregs, tumor cells) without dendritic cell, macrophage, or NK cell contributions.
- **DDI not modeled**: Drug-drug pharmacokinetic interactions (e.g., CYP3A4 induction/inhibition affecting CDK4/6i exposure) are not explicitly represented for combination regimens.
- **Cardiac toxicity**: LVEF decline (trastuzumab cardiotoxicity, ribociclib QTc prolongation) is modeled empirically as a dose-dependent function rather than mechanistically via cardiomyocyte biology.
- **Single tumor site**: Metastatic spread to multiple organ sites (bone, lung, liver, brain) is represented by a single virtual tumor compartment; site-specific PK differences are approximated by scaling parameters.
- **Population homogeneity**: The base model uses a single representative patient; population variability (BSV) must be introduced via IIV parameter distributions for population-level predictions.

---

## Development Information

| Item | Details |
|------|---------|
| Created | 2026-06-21 |
| Version | 1.0.0 |
| Author | Claude Code Routine (CCR) — Automated QSP Library |
| Software | R ≥ 4.0, mrgsolve ≥ 1.0, Graphviz ≥ 2.40, Shiny ≥ 1.7 |
| Framework | ODE-based, deterministic; Simeoni-type tumor growth backbone |
| Validation | Calibrated against PALOMA-2, MONALEESA-2, MONARCH-3, CLEOPATRA, KEYNOTE-522, OlympiAD primary endpoints |
| References | See [`bc_references_en.md`](bc_references_en.md) (50 PubMed-linked citations, 9 sections) |

---

*This model is part of the QSP Disease Model Library (CCR). See the main [README](../README_en.md) for the full library index.*
