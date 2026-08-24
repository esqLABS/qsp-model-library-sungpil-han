# Psoriasis — QSP Model

> **Category:** Chronic Autoimmune Skin Disease
> **Key Pathway:** IL-23/IL-17 axis · TNF-α · Keratinocyte hyperproliferation
> **Date Added:** 2026-06-20

---

## Overview

Psoriasis vulgaris is a chronic immune-mediated skin disease affecting about 2-3% of the global
population, in which hyperproliferation of keratinocytes and immune-cell infiltration form the
characteristic scaly erythematous plaques.

**Core Pathogenic Mechanism:**
- Dendritic cells (mDC) → secrete **IL-23** → drive **Th17** differentiation
- Th17 → secrete **IL-17A** (the principal effector cytokine)
- IL-17A → activates keratinocyte **NF-κB** → hyperproliferation (acanthosis)
- TNF-α → activates vascular endothelial cells → amplifies immune-cell influx into skin

---

## Mechanistic Map

[![Psoriasis QSP Mechanistic Map](pso_qsp_model.png)](pso_qsp_model.svg)

*Click to zoom into the SVG vector image*

### Key Cluster Composition
| Cluster | Components |
|----------|-----------|
| Genetic/environmental triggers | HLA-Cw6, IL23R, CARD14, streptococcus, physical trauma, drugs, obesity |
| Innate immune activation | pDC, mDC, Langerhans cells, macrophages, neutrophils, LL-37–DNA complexes, TLR7/9 |
| Adaptive immunity (T cells) | Th17 (★), Th1, Th22, Treg, ILC3, tissue-resident memory T cells |
| Key cytokines | IL-17A (★), IL-22, TNF-α (★), IFN-γ, IL-36, TSLP |
| Keratinocyte pathology | NF-κB/STAT3 activation, hyperproliferation, abnormal differentiation, impaired barrier function |
| Vascular/stromal remodelling | VEGF-A-induced angiogenesis, ICAM-1/VCAM-1, leukocyte extravasation |
| Clinical measures | PASI, BSA, PGA, DLQI, psoriatic arthritis |
| Biologic PK | Adalimumab · secukinumab · risankizumab · ustekinumab |
| Small-molecule drug PK | Apremilast · tofacitinib · methotrexate |

---

## mrgsolve ODE Model (Pharmacokinetic/Pharmacodynamic Model)

**File:** [`pso_mrgsolve_model.R`](pso_mrgsolve_model.R)

### Model Structure (25 ODE Compartments)

| Compartment Category | Number of Compartments | Key Compartments |
|-----------|---------|-----------|
| Immune/cellular | 8 | DC, IL-23, Th17, IL-17A, TNF-α, IFN-γ, KC index, PASI |
| Adalimumab PK | 4 | SC depot, central, peripheral, TNF-bound complex |
| Secukinumab PK | 4 | SC depot, central, peripheral, IL-17A-bound |
| Risankizumab PK | 3 | SC depot, central, IL-23-bound |
| Ustekinumab PK | 3 | SC depot, central, p40-bound |
| Apremilast PK | 3 | GI, central, peripheral |
| Tofacitinib PK | 2 | GI, central |
| Methotrexate PK | 3 | GI, central, polyglutamate |

### Simulation Scenarios (7)

| Scenario | Drug | Regimen | Mechanism |
|----------|------|------|------|
| 1 | Untreated | — | Natural history |
| 2 | Adalimumab | 40mg SC q2w | TNF-α neutralisation |
| 3 | Secukinumab | 300mg SC wk0-4, q4w | IL-17A neutralisation |
| 4 | Risankizumab | 150mg SC wk0,4, q12w | IL-23 p19 inhibition |
| 5 | Apremilast | 30mg PO BID | PDE4 inhibition → increased cAMP |
| 6 | Tofacitinib | 10mg PO BID | JAK1/3 inhibition → decreased STAT3/1 |
| 7 | Methotrexate | 20mg PO once weekly | DHFR inhibition → decreased folate |

### Key PK Parameters

| Drug | Bioavailability (F) | CL (L/h) | V₁ (L) | t½ (days) | Kd / IC50 |
|------|-------------|---------|-------|---------|-----------|
| Adalimumab | 0.64 | 0.247 | 7.0 | ~14d | Kd=0.1 nM (TNF-α) |
| Secukinumab | 0.73 | 0.191 | 7.1 | ~27d | Kd=0.08 nM (IL-17A) |
| Risankizumab | 0.89 | 0.078 | 11.2 | ~28d | Kd=0.06 nM (IL-23p19) |
| Ustekinumab | 0.57 | 0.252 | 15.1 | ~21d | Kd=0.9 nM (p40) |
| Apremilast | 0.73 | 9.5 | 86.6 | 9h | IC50=74 nM (PDE4) |
| Tofacitinib | 0.74 | 22.8 | 87.0 | 3h | IC50=1-5 nM (JAK1/3) |
| Methotrexate | 0.70 | 4.8 | 24.0 | 3-10h | IC50=1 nM (DHFR) |

### Clinical Calibration Data (PASI Response Rate at Week 16)

| Drug | PASI75 (Actual) | PASI90 (Actual) | Clinical Trial |
|------|-------------|-------------|---------|
| Adalimumab 40mg q2w | 71% | 45% | CHAMPION |
| Secukinumab 300mg q4w | 77-80% | 59-67% | FIXTURE |
| Risankizumab 150mg q12w | 88-91% | 72-75% | UltIMMa |
| Ustekinumab 45mg q12w | 67-71% | 42% | PHOENIX-1 |
| Apremilast 30mg BID | 33-40% | 18% | ESTEEM-1 |
| Tofacitinib 10mg BID | 39-46% | 22% | OPT Pivotal |
| Methotrexate 20mg qw | 26-36% | — | Heydendael 2003 |

---

## Shiny Interactive App (Interactive Dashboard)

**File:** [`pso_shiny_app.R`](pso_shiny_app.R)

### 6-Tab Structure

| Tab | Contents |
|----|------|
| **Patient profile** | Setting baseline PASI, IL-17A, TNF-α, Th17 cell counts; disease phenotype selection |
| **PK profile** | Blood concentration-time curves for biologics and small-molecule drugs |
| **Cytokine PD** | IL-17A, TNF-α, Th17, IL-23 dynamics |
| **PASI endpoint** | PASI score change, PASI75/90/100 response rates |
| **Scenario comparison** | Simultaneous comparison of 7 treatments |
| **Biomarker dashboard** | Response-rate value boxes, heatmap, biomarker time series |

---

## References

**File:** [`pso_references.md`](pso_references.md)

- 42 references in total (with PubMed links)
- Grouped by section: pathogenesis · IL-23/Th17 axis · IL-17A/keratinocytes · TNF-α · clinical
  trials · PK/PD models · biomarkers · TYK2 inhibitors (emerging)

---

## File List

| File | Description |
|------|------|
| [`pso_qsp_model.dot`](pso_qsp_model.dot) | Graphviz mechanistic map source (100+ nodes, 11 clusters) |
| [`pso_qsp_model.svg`](pso_qsp_model.svg) | SVG vector image (zoomable) |
| [`pso_qsp_model.png`](pso_qsp_model.png) | PNG raster image (150 dpi) |
| [`pso_mrgsolve_model.R`](pso_mrgsolve_model.R) | mrgsolve ODE model (25 compartments, 7 scenarios) |
| [`pso_shiny_app.R`](pso_shiny_app.R) | Shiny dashboard (6 tabs) |
| [`pso_references.md`](pso_references.md) | 42 references (PubMed links) |

---

## Key Clinical Context

- **Prevalence:** 2-3% worldwide, about 1.5% in Korea
- **Moderate-to-severe criteria:** PASI ≥ 10, or BSA ≥ 10%, or DLQI ≥ 10
- **Treatment goal:** PASI90 (90% skin improvement) or IGA 0/1 (complete/almost complete
  clearance)
- **Psoriatic arthritis comorbidity:** occurs in about 30% of patients; IL-17A and TNF-α
  inhibitors are effective when combined
- **Comorbidities:** increased risk of cardiovascular disease (CVD), metabolic syndrome,
  inflammatory bowel disease (IBD), and depression

---

## Emerging Targets

| Target | Drug | Development Stage |
|------|------|-----------|
| TYK2 (IL-23/IFN-γ pathway) | Deucravacitinib | FDA-approved (2022) |
| IL-17A + IL-17F | Bimekizumab | FDA-approved (2023) |
| IL-13 | Tralokinumab | Exploratory clinical trials |
| OX40L | Amlitelimab | Phase 3 ongoing |
| IL-4Rα | Dupilumab | Exploring a psoriasis indication |
