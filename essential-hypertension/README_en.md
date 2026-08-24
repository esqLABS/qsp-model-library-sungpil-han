# Essential Hypertension (EH) — QSP Model

> **QSP Disease Model Library** · A quantitative systems pharmacology (QSP) model automatically generated via the Claude Code Routine.
> Parent library → [../README_en.md](../README_en.md) · Category: Cardiovascular

[![EH QSP Model](eh_qsp_model.png)](eh_qsp_model.svg)

## Overview

Essential hypertension is a rise in blood pressure (systolic ≥ 130 mmHg or diastolic ≥ 80 mmHg) occurring without a secondary cause, affecting approximately 30~45% of adults worldwide. It is a major cause of myocardial infarction, stroke, and chronic renal failure, and the leading risk factor for death worldwide. The core pathophysiology is a combined mechanism in which RAAS overactivation, increased sympathetic tone, and abnormal renal sodium/fluid homeostasis raise both cardiac output (CO) and total peripheral resistance (TPR). ACEi/ARB, CCBs, diuretics, and beta-blockers are first-line therapies.

## Key Pathways

| Pathway | Key molecule/mechanism | Clinical outcome |
|------|----------------|-----------|
| RAAS axis | Renin → angiotensinogen → Ang I → Ang II → AT1R | Vasoconstriction, aldosterone secretion, Na retention |
| Sympathetic activation | Norepinephrine → α1R (vasoconstriction), β1R (heart rate/contractility) | Increased cardiac output and TPR |
| Renal sodium homeostasis | Ang II/aldosterone → collecting duct Na reabsorption | Fluid excess, raised blood pressure |
| Endothelial dysfunction | Ang II → ROS → reduced NO | Impaired vasodilatory capacity |
| Pressure-natriuresis relationship | Raised blood pressure → increased renal pressure → delayed natriuresis | Resetting of the blood pressure set-point |
| Vascular remodelling | Chronic hypertension → medial hypertrophy → increased TPR | Entrenchment of hypertension |
| End-organ damage | Ventricular hypertrophy, proteinuria, retinal lesions | Cardiovascular and renal complications |

## Drug Targets

- **ACE inhibitors** (ramipril, enalapril): inhibit ACE → reduce Ang II, accumulate bradykinin → vasodilation
- **ARBs** (losartan, valsartan): directly block AT1R → block the effects of Ang II
- **CCBs** (amlodipine): block L-type calcium channels → relax vascular smooth muscle, reduce TPR
- **Beta-blockers** (bisoprolol): block β1AR → reduce heart rate and contractility, suppress renin secretion
- **Thiazide diuretics** (HCTZ): inhibit the distal tubular Na-Cl cotransporter → increase natriuresis

## Model Files

| File | Description |
|------|------|
| [eh_qsp_model.dot](eh_qsp_model.dot) | Graphviz mechanistic map source (approximately 173 nodes / 10 clusters) |
| [eh_qsp_model.svg](eh_qsp_model.svg) | SVG vector image (scalable) |
| [eh_qsp_model.png](eh_qsp_model.png) | PNG image (150 dpi) |
| [eh_mrgsolve_model_en.R](eh_mrgsolve_model_en.R) | mrgsolve ODE model (approximately 22 compartments / 6 treatment scenarios) |
| [eh_shiny_app.R](eh_shiny_app.R) | Shiny dashboard |
| [eh_references_en.md](eh_references_en.md) | References (approximately 41 articles, PubMed links) |

## mrgsolve Model (ODE Model)

- **Compartment structure**: drug PK compartments (central/peripheral compartments each for ACEi, ARB, CCB, beta-blocker, HCTZ) + PD compartments (renin, angiotensin II, aldosterone, fluid volume, cardiac output, total peripheral resistance, systolic/diastolic blood pressure, natriuresis, heart rate)
- **Key treatment scenarios**: ① no treatment (untreated hypertension), ② ACE inhibitor monotherapy (ramipril 10 mg QD), ③ ARB monotherapy (losartan 100 mg QD), ④ CCB monotherapy (amlodipine 10 mg QD), ⑤ beta-blocker monotherapy (bisoprolol 10 mg QD), ⑥ triple combination (ACEi + CCB + thiazide) — the standard first-line combination therapy
- **Calibration/evidence**: parameters referenced from ALLHAT (amlodipine/lisinopril/chlortalidone), the HOT trial, and HOPE trial (ramipril) data

## Shiny Dashboard

Comprises 6 tabs: ① Patient Profile (sets baseline blood pressure, renal function, salt intake, and BMI), ② PK tab (blood concentration and active metabolites of 5 drugs), ③ RAAS/haemodynamic PD tab (Ang II, aldosterone, TPR, and CO trends), ④ Clinical Endpoints (systolic/diastolic blood pressure, heart rate, achievement of blood pressure target), ⑤ Scenario Comparison (simultaneous comparison of 6 treatment strategies), ⑥ Biomarkers (renin, BNP, renal function, and natriuresis trends).

## Usage

```r
library(mrgsolve)
mod <- mread("eh_mrgsolve_model_en.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny dashboard:
shiny::runApp("eh_shiny_app.R")
```
```bash
# Render the mechanistic map
dot -Tsvg eh_qsp_model.dot -o eh_qsp_model.svg
```

## References

For detailed citations, see [eh_references_en.md](eh_references_en.md) (approximately 41 articles).

---
*This model is a qualitative/semi-quantitative QSP model intended for education and research, and must not be used directly for clinical decision-making.*
