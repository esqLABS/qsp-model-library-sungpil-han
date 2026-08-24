# Chronic Kidney Disease (CKD) — QSP Model

> **QSP Disease Model Library** · A quantitative systems pharmacology (QSP) model automatically generated via the Claude Code Routine.
> Parent library → [../README_en.md](../README_en.md) · Category: Renal/Urological

[![CKD QSP Model](ckd_qsp_model.png)](ckd_qsp_model.svg)

## Overview

Chronic kidney disease (CKD) occurs in approximately 10~13% of adults worldwide, with diabetic nephropathy and hypertension as the most common causes. A vicious cycle of glomerular hyperfiltration (early stage) → nephron loss → compensatory hypertrophy → proteinuria/fibrosis progressively reduces eGFR over years. Triple combination therapy with a RAAS inhibitor (ACEi/ARB), an SGLT2 inhibitor (dapagliflozin, DAPA-CKD), and a non-steroidal MRA (finerenone, FIDELIO-DKD) is now established as the most powerful renoprotective strategy.

## Key Pathways

| Pathway | Key molecule/mechanism | Clinical outcome |
|------|----------------|-----------|
| Glomerular hyperfiltration/proteinuria | Angiotensin II → efferent arteriolar constriction → raised intraglomerular pressure | Increased UACR, podocyte injury |
| RAAS overactivation | Ang II → aldosterone → TGF-β secretion | Renal fibrosis, raised blood pressure, inflammation |
| Oxidative stress/NF-κB | Increased ROS, TNF-α, IL-6 | Tubulointerstitial inflammation |
| Interstitial fibrosis | TGF-β → myofibroblasts → excess collagen deposition | Nephron replacement → reduced eGFR |
| CKD-MBD | FGF-23 ↑, klotho ↓, reduced vitamin D → PTH ↑ | Vascular calcification, secondary hyperparathyroidism |
| Renal anaemia | Reduced EPO production, hepcidin ↑ → reduced haemoglobin | Fatigue, increased cardiovascular load |
| Left ventricular hypertrophy (LVH) | Blood pressure, anaemia, and RAAS overactivation → myocardial hypertrophy | Increased cardiovascular mortality risk |

## Drug Targets

- **ACE inhibitors/ARBs**: ramipril, losartan — inhibit Ang II → reduce intraglomerular pressure, reduce proteinuria, suppress TGF-β
- **SGLT2 inhibitor (dapagliflozin)**: blocks proximal tubular sodium-glucose reabsorption → reduces intraglomerular pressure, renoprotective (DAPA-CKD)
- **Finerenone (non-steroidal MRA)**: blocks the mineralocorticoid receptor → suppresses fibrosis and inflammation (FIDELIO/FIGARO-DKD)
- **ESA (erythropoiesis-stimulating agent)/darbepoetin**: stimulates the EPO receptor → corrects anaemia
- **Phosphate binders/active vitamin D**: regulate phosphate/PTH → manage CKD-MBD

## Model Files

| File | Description |
|------|------|
| [ckd_qsp_model.dot](ckd_qsp_model.dot) | Graphviz mechanistic map source (approximately 484 nodes / 10 clusters) |
| [ckd_qsp_model.svg](ckd_qsp_model.svg) | SVG vector image (scalable) |
| [ckd_qsp_model.png](ckd_qsp_model.png) | PNG image (150 dpi) |
| [ckd_mrgsolve_model.R](ckd_mrgsolve_model.R) | mrgsolve ODE model (approximately 30 compartments / 5 treatment scenarios) |
| [ckd_shiny_app.R](ckd_shiny_app.R) | Shiny dashboard |
| [ckd_references_en.md](ckd_references_en.md) | References (approximately 37 articles, PubMed links) |

## mrgsolve Model (ODE Model)

- **Compartment structure**: PK compartments for ACEi/ARB, finerenone, SGLT2i, ESA, and phosphate binder + PD compartments for nephron number, eGFR, UACR, Ang II, aldosterone, blood pressure, macrophages, IL-6, TNF-α, TGF-β, collagen, phosphate, klotho, FGF-23, vitamin D, PTH, EPO, hepcidin, haemoglobin, left ventricular hypertrophy index, and vascular calcification
- **Key treatment scenarios**: (1) natural course (no treatment), (2) ACEi monotherapy (ramipril 10 mg), (3) ACEi + finerenone 20 mg, (4) ACEi + dapagliflozin 10 mg (DAPA-CKD regimen), (5) triple combination (ACEi + Dapa + finerenone)
- **Calibration/evidence**: rate of eGFR decline and UACR response based on the DAPA-CKD (NEJM 2020) and FIDELIO-DKD (NEJM 2020) clinical trial data

## Shiny Dashboard

Comprises 6 tabs: (1) Patient Profile — sets CKD stage, presence of diabetes, and baseline proteinuria; (2) PK tab — plasma concentration of ACEi/dapagliflozin/finerenone; (3) Key PD Metrics — eGFR trend, UACR, Ang II; (4) Clinical Endpoints — renal failure progression, haemoglobin, CKD-MBD markers; (5) Scenario Comparison — 5-year outcomes across 5 treatment strategies; (6) Biomarkers — FGF-23, klotho, PTH, vascular calcification score

## Usage

```r
library(mrgsolve)
mod <- mread("ckd_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny dashboard:
shiny::runApp("ckd_shiny_app.R")
```
```bash
# Render the mechanistic map
dot -Tsvg ckd_qsp_model.dot -o ckd_qsp_model.svg
```

## References

For detailed citations, see [ckd_references_en.md](ckd_references_en.md) (approximately 37 articles).

---
*This model is a qualitative/semi-quantitative QSP model intended for education and research, and must not be used directly for clinical decision-making.*
