# Chronic Obstructive Pulmonary Disease (COPD) — QSP Model

> **QSP Disease Model Library** · A quantitative systems pharmacology (QSP) model automatically generated via the Claude Code Routine.
> Parent library → [../README_en.md](../README_en.md) · Category: Respiratory

[![COPD QSP Model](copd_qsp.png)](copd_qsp.svg)

## Overview

COPD is a major chronic respiratory disease affecting approximately 390 million people worldwide, with smoking as its most common cause. Smoking-induced oxidative stress and overactivation of proteases (neutrophil elastase, MMPs) cause airway epithelial damage and alveolar destruction (emphysema), while persistent airway inflammation mediated by IL-8 and TNF-α drives progressive airway narrowing. Under the GOLD guidelines, LAMA/LABA bronchodilators are the foundation of treatment, with ICS added for the eosinophilic or frequent-exacerbator phenotype, and a PDE4 inhibitor (roflumilast) combined in severe emphysematous disease.

## Key Pathways

| Pathway | Key molecule/mechanism | Clinical outcome |
|------|----------------|-----------|
| Smoking/oxidative stress | ROS → NF-κB → IL-8/TNF-α secretion | Neutrophil/macrophage airway infiltration |
| Protease/antiprotease imbalance | NE/MMP-9 ↑, α1-AT/TIMP ↓ | Destruction of alveolar septa → emphysema |
| Airway mucus hypersecretion | Goblet cell hyperplasia, MUC5AC overexpression | Chronic cough and sputum |
| Airway remodelling | Fibroblast activation, airway wall thickening | Airflow limitation (reduced FEV₁) |
| Eosinophilic airway inflammation | IL-5/IL-13 → eosinophil accumulation | Predicts ICS response, exacerbation risk |
| Pulmonary vascular remodelling | Hypoxia → pulmonary vasoconstriction → pulmonary hypertension | Right heart failure (cor pulmonale) |
| Acute exacerbation (AECOPD) | Viral/bacterial infection plus underlying inflammation | Rapid decline in lung function, hospitalisation |

## Drug Targets

- **LAMA (long-acting muscarinic antagonist)**: blocks the M3 receptor → bronchodilation → improves FEV₁ (tiotropium, glycopyrronium)
- **LABA (long-acting beta₂ agonist)**: stimulates the β₂ receptor → relaxes smooth muscle (salmeterol, olodaterol)
- **ICS (inhaled corticosteroid)**: suppresses IL-8/eosinophilic inflammation, prevents exacerbation (fluticasone, budesonide)
- **PDE4 inhibitor (roflumilast)**: inhibits cAMP degradation → anti-inflammatory effect, improves FEV₁, reduces acute exacerbations
- **Antibiotics (azithromycin)**: long-term prophylactic use in the frequent-exacerbator phenotype, combining antibacterial and anti-inflammatory effects

## Model Files

| File | Description |
|------|------|
| [copd_qsp.dot](copd_qsp.dot) | Graphviz mechanistic map source (approximately 430 nodes / 12 clusters) |
| [copd_qsp.svg](copd_qsp.svg) | SVG vector image (scalable) |
| [copd_qsp.png](copd_qsp.png) | PNG image (150 dpi) |
| [copd_mrgsolve_model.R](copd_mrgsolve_model.R) | mrgsolve ODE model (approximately 26 compartments / 6 treatment scenarios) |
| [copd_shiny_app.R](copd_shiny_app.R) | Shiny dashboard |
| [copd_references.md](copd_references.md) | References (approximately 52 articles, PubMed links) |

## mrgsolve Model (ODE Model)

- **Compartment structure**: LAMA/LABA/ICS (2 compartments each: lung and central/peripheral) + PDE4i plasma PK compartment + PD compartments for IL-8, neutrophil elastase (NE), CRP, eosinophils, FEV₁, emphysema index (Emph), pulmonary vascular resistance (PVR), cumulative/annual exacerbation rate, and CAT score
- **Key treatment scenarios**: (1) no treatment, (2) LAMA monotherapy, (3) LAMA + LABA combination, (4) LAMA + LABA + ICS triple therapy, (5) LAMA + LABA + ICS + PDE4i quadruple therapy, (6) smoking cessation alone
- **Calibration/evidence**: FEV₁ change based on the UPLIFT (tiotropium) and TRILOGY/TRINITY (ICS-LABA-LAMA) clinical trial data

## Shiny Dashboard

Comprises 6 tabs: (1) Patient Profile — sets GOLD stage, smoking history, and eosinophil count; (2) PK tab — lung and plasma concentrations of LAMA/LABA/ICS/PDE4i; (3) Key PD Metrics — IL-8, neutrophil elastase, CRP, eosinophils; (4) Clinical Endpoints — FEV₁ trend, cumulative acute exacerbations, CAT score; (5) Scenario Comparison — 1-year outcomes across 6 treatment strategies; (6) Biomarkers — serum fibrinogen, blood eosinophils, emphysema progression index

## Usage

```r
library(mrgsolve)
mod <- mread("copd_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny dashboard:
shiny::runApp("copd_shiny_app.R")
```
```bash
# Render the mechanistic map
dot -Tsvg copd_qsp.dot -o copd_qsp.svg
```

## References

For detailed citations, see [copd_references.md](copd_references.md) (approximately 52 articles).

---
*This model is a qualitative/semi-quantitative QSP model intended for education and research, and must not be used directly for clinical decision-making.*
