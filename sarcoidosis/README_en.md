# Sarcoidosis (SARC) — QSP Model

> **QSP Disease Model Library** · A Quantitative Systems Pharmacology (QSP) model automatically generated via Claude Code Routine.
> Parent library → [../README_en.md](../README_en.md) · Category: Respiratory

[![SARC QSP Model](sarc_qsp_model.png)](sarc_qsp_model.svg)

## Overview

Sarcoidosis is a systemic inflammatory disease of unknown cause in which
non-caseating granulomas form across multiple organs, including the
lungs, lymph nodes, skin, eyes, and heart. Its prevalence is 10-40 per
100,000 and it is more common in Black and Scandinavian populations. The
core pathogenic mechanism is a Th1/IFN-γ-skewed immune response and
activated macrophages that mutually stimulate each other to form and
sustain granulomas. About two-thirds of cases resolve spontaneously, but
the remainder progresses to pulmonary fibrosis, pulmonary hypertension,
and cardiac sarcoidosis, all of which are life-threatening.
Corticosteroids are first-line therapy, with methotrexate and
azathioprine, and anti-TNF agents for refractory disease.

## Key Pathophysiological Pathways

| Pathway | Key molecules/mechanism | Clinical result |
|------|----------------|-----------|
| Antigen presentation and Th1 skewing | Unknown antigen (mycobacterium-like?) → APC → CD4+ Th1 → excess IFN-γ/IL-2 | Granuloma initiation, pulmonary infiltration |
| Macrophage M1 activation | IFN-γ → M1 macrophage polarisation → TNF, IL-12, IL-18 → epithelioid giant cells | Formation of non-caseating granulomas |
| Granuloma formation and maintenance | TNF, IFN-γ → CXCL10, CCL2 → sustained recruitment of CD4+ T cells and monocytes | Granuloma expansion and persistence |
| TGF-β/fibrosis pathway | Chronic granuloma → TGF-β1 secretion → myofibroblast activation | Progressive pulmonary fibrosis |
| Abnormal calcium metabolism | Overactive 1α-hydroxylase in granuloma macrophages → elevated 1,25(OH)₂D₃ | Hypercalcaemia, nephrolithiasis |
| Cardiac/conduction sarcoidosis | Myocardial granuloma and fibrosis → conduction disturbance and ventricular arrhythmia | Risk of sudden cardiac death |
| NF-κB/mTOR pathway | TNF, IL-1 → NF-κB; nutrient signalling → mTOR → macrophage proliferation | Survival signal sustaining chronic granulomas |

## Key Drug Targets

- **Corticosteroids (prednisolone)**: broad anti-inflammatory/Th1 suppression; first-line for symptomatic pulmonary, neurological, and cardiac sarcoidosis
- **Methotrexate**: folate antagonism → suppresses T-cell and macrophage function; a steroid-sparing agent
- **Azathioprine/mycophenolate**: inhibits purine/pyrimidine synthesis; chronic maintenance therapy
- **Hydroxychloroquine**: inhibits TLR9, blocks lysosomal acidification; used in cutaneous and hypercalcaemic sarcoidosis
- **Anti-TNF agents (infliximab, adalimumab)**: neutralises TNF → blocks the granuloma-maintenance signal; for refractory pulmonary, neurological, and cardiac sarcoidosis
- **Tetracyclines/clofazimine**: limited evidence; adjunctive in some neurosarcoidosis

## Model Files

| File | Description |
|------|------|
| [sarc_qsp_model.dot](sarc_qsp_model.dot) | Graphviz mechanistic map source (about 100+ nodes / 12 clusters) |
| [sarc_qsp_model.svg](sarc_qsp_model.svg) | SVG vector image (zoomable) |
| [sarc_qsp_model.png](sarc_qsp_model.png) | PNG image (150 dpi) |
| [sarc_mrgsolve_model.R](sarc_mrgsolve_model.R) | mrgsolve ODE model (about 21 compartments / about 17 scenarios) |
| [sarc_shiny_app.R](sarc_shiny_app.R) | Shiny dashboard |
| [sarc_references.md](sarc_references.md) | References (about 50, PubMed links) |

## mrgsolve Model (ODE Model)

- **Compartment structure**: steroid/methotrexate/infliximab PK compartments plus Th1 activation and IFN-γ kinetics, a granuloma-size index, a TGF-β/fibrosis module, serum ACE and calcium prediction compartments, and a module simulating pulmonary function (FVC/DLCO)
- **Key treatment scenarios**: an untreated arm split between spontaneous remission and progression, prednisolone alone, methotrexate, azathioprine, prednisolone + methotrexate, infliximab, infliximab + methotrexate, hydroxychloroquine, and others
- **Calibration/basis**: based on the steroid RCT by Baughman et al. in *Am J Respir Crit Care Med*, the GRAPPA infliximab study, and clinical data from European sarcoidosis guidelines

## Shiny Dashboard

Structured into 6 tabs: (1) **Patient profile** — sets Scadding stage,
organ involvement, serum ACE, calcium; (2) **PK profile** — blood
concentrations of steroids/immunosuppressants/biologics; (3) **Key PD
measures** — IFN-γ/TNF inhibition rate, granuloma-size kinetics; (4)
**Clinical endpoints** — FVC, DLCO, 6MWD, serum ACE over time; (5)
**Scenario comparison** — remission rate and fibrosis progression compared
across treatment strategies; (6) **Biomarkers** — serum ACE, calcium,
sIL-2R, chest X-ray stage.

## Usage

```r
library(mrgsolve)
mod <- mread("sarc_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny dashboard:
shiny::runApp("sarc_shiny_app.R")
```
```bash
# Render the mechanistic map
dot -Tsvg sarc_qsp_model.dot -o sarc_qsp_model.svg
```

## References

For full citations, see [sarc_references.md](sarc_references.md) (about 50 entries).

---
*This model is a qualitative/semi-quantitative QSP model for educational and research purposes and must not be used directly for clinical decision-making.*
