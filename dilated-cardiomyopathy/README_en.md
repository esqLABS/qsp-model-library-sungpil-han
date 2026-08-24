# Dilated Cardiomyopathy (DCM) — QSP Model

> **QSP Disease Model Library** · A quantitative systems pharmacology (QSP) model automatically generated via the Claude Code Routine.
> Parent library → [../README_en.md](../README_en.md) · Category: Cardiovascular

[![DCM QSP Model](dcm_qsp_model.png)](dcm_qsp_model.svg)

## Overview

Dilated cardiomyopathy (DCM) is a myocardial disease characterised by ventricular dilation and reduced systolic function (EF < 40%), and is one of the most common causes of heart failure. Its prevalence is approximately 36 per 100,000 population, with diverse causes including genetic (e.g. titin variants), viral, toxic, and autoimmune aetiologies. The core pathophysiology is a vicious cycle in which neurohormonal (RAAS/sympathetic) activation following myocardial injury drives adverse ventricular remodelling; GDMT (ARNI, beta-blockers, MRA, SGLT2i) induces reverse remodelling by suppressing neurohormonal activation.

## Key Pathways

| Pathway | Key molecule/mechanism | Clinical outcome |
|------|----------------|-----------|
| RAAS activation | Angiotensin II → AT1R → myocardial hypertrophy/fibrosis | LV dilation, raised LVEDP |
| Sympathetic activation | Norepinephrine → β1AR → increased heart rate and contractility | Increased myocardial oxygen consumption, cell death |
| Aldosterone excess | Promotes sodium retention and myocardial fibrosis | Ventricular stiffness, arrhythmia risk |
| Natriuretic peptides | ANP/BNP → vasodilation, natriuresis | Elevated BNP (a marker of heart failure progression) |
| Myocardial fibrosis | TGF-β → collagen synthesis | Ventricular stiffness, reduced ejection function |
| Mitochondrial dysfunction | Abnormal energy metabolism, increased ROS | Cardiomyocyte apoptosis |
| Systemic inflammation | TNF-α, IL-6 → myocardial depressant effect | Accelerated heart failure progression |

## Drug Targets

- **ARNI** (sacubitril/valsartan): inhibits neprilysin (increases BNP) plus ARB action (blocks AT1R) → dual neurohormonal suppression
- **Beta-blockers** (carvedilol, metoprolol): block β1AR → reduce heart rate, protect the myocardium
- **MRA** (spironolactone, eplerenone): block the aldosterone receptor → suppress fibrosis
- **SGLT2 inhibitors** (dapagliflozin, empagliflozin): improve osmotic diuresis and energy metabolism
- **ACE inhibitors** (enalapril): suppress angiotensin II generation, accumulate bradykinin

## Model Files

| File | Description |
|------|------|
| [dcm_qsp_model.dot](dcm_qsp_model.dot) | Graphviz mechanistic map source (approximately 148 nodes / 11 clusters) |
| [dcm_qsp_model.svg](dcm_qsp_model.svg) | SVG vector image (scalable) |
| [dcm_qsp_model.png](dcm_qsp_model.png) | PNG image (150 dpi) |
| [dcm_mrgsolve_model.R](dcm_mrgsolve_model.R) | mrgsolve ODE model (approximately 24 compartments / 5 treatment scenarios) |
| [dcm_shiny_app.R](dcm_shiny_app.R) | Shiny dashboard |
| [dcm_references.md](dcm_references.md) | References (approximately 52 articles, PubMed links) |

## mrgsolve Model (ODE Model)

- **Compartment structure**: drug PK compartments (depot + central compartments each for enalapril, carvedilol, spironolactone, sacubitril, dapagliflozin) + PD compartments (angiotensin II, aldosterone, BNP, LVEF, LV volume, heart rate, fibrosis index)
- **Key treatment scenarios**: ① no treatment (natural course), ② ACE inhibitor monotherapy, ③ ACE inhibitor + beta-blocker, ④ ARNI + beta-blocker + MRA (triple therapy), ⑤ full four-drug GDMT (ARNI + beta-blocker + MRA + SGLT2i)
- **Calibration/evidence**: parameters referenced from the PARADIGM-HF (LCZ696), DAPA-HF, and EMPEROR-Reduced clinical trial data and the Konstam et al. reverse remodelling literature

## Shiny Dashboard

Comprises 6 tabs: ① Patient Profile (sets baseline LVEF, BNP, and NYHA class), ② PK tab (blood concentration trends for 5 drugs), ③ Cardiac PD tab (changes in LVEF, LV volume, BNP), ④ Clinical Endpoints (6-minute walk distance, NYHA improvement, estimated hospitalisation rate), ⑤ Scenario Comparison (simultaneous comparison of 5 treatment strategies), ⑥ Biomarkers (angiotensin II, aldosterone, and fibrosis index trends).

## Usage

```r
library(mrgsolve)
mod <- mread("dcm_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny dashboard:
shiny::runApp("dcm_shiny_app.R")
```
```bash
# Render the mechanistic map
dot -Tsvg dcm_qsp_model.dot -o dcm_qsp_model.svg
```

## References

For detailed citations, see [dcm_references.md](dcm_references.md) (approximately 52 articles).

---
*This model is a qualitative/semi-quantitative QSP model intended for education and research, and must not be used directly for clinical decision-making.*
