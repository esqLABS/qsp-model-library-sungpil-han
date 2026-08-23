# Hypertrophic Cardiomyopathy (HCM) — QSP Model

> **QSP Disease Model Library** · A quantitative systems pharmacology (QSP) model automatically generated via the Claude Code Routine.
> Parent library → [../README.md](../README.md) · Category: Cardiovascular

[![HCM QSP Model](../../../hypertrophic-cardiomyopathy/hcm_qsp_model.png)](../../../hypertrophic-cardiomyopathy/hcm_qsp_model.svg)

## Overview

Hypertrophic cardiomyopathy (HCM) is the most common inherited cardiomyopathy, with a prevalence of 1 in 500, and its core pathogenic mechanism is hypercontractility driven by mutations in sarcomere protein genes such as MYH7 and MYBPC3. The mutant sarcomere increases the transition from the slow-relaxed (SRX) state to the disordered-relaxed/fast-active (DRX) state during the contraction-relaxation cycle, causing increased ATP consumption, heat generation, myocardial hypertrophy, and fibrosis. Obstructive HCM, in which left ventricular outflow tract obstruction (LVOTO) is present, accounts for approximately 70% of cases and is a cause of exercise intolerance, syncope, atrial fibrillation, and sudden death. Mavacamten, the first cardiac myosin inhibitor, increases the SRX state to suppress hypercontractility.

## Key Pathways

| Pathway | Key molecule/mechanism | Clinical outcome |
|------|----------------|-----------|
| Sarcomere hypercontractility | Increased SRX→DRX transition → excess cross-bridge formation → hypercontractility | Left ventricular hypertrophy, impaired relaxation |
| Left ventricular outflow tract obstruction | Hypertrophied septum + systolic anterior motion (SAM) → dynamic LVOTO | Increased pressure gradient, dyspnoea |
| Calcium signalling abnormality | Intracellular Ca²⁺ overload → calcineurin-NFAT → activation of hypertrophic genes | Myocardial hypertrophic gene expression |
| Myocardial fibrosis | TGF-β → cardiac fibroblast activation → collagen deposition | Interstitial fibrosis, diastolic dysfunction |
| ERK/MAPK signalling | Mechanical stress → ERK1/2 → activation of cell growth pathways | Maintenance of hypertrophy |
| Atrial fibrillation/sudden death | Raised LVEDP and fibrosis → atrial enlargement, ventricular tachycardia | Stroke, risk of sudden death |

## Drug Targets

- **Mavacamten**: inhibits β-cardiac myosin ATPase → restores the SRX state → reduces hypercontractility and relieves LVOTO; supported by the EXPLORER-HCM and VALOR-HCM trials
- **Aficamten**: a next-generation myosin inhibitor; reduced LVOT pressure gradient similarly to EXPLORER-HCM in the phase 3 SEQUOIA-HCM trial
- **Beta-blockers (metoprolol, atenolol)**: reduce heart rate and contractility → relieve LVOTO and symptoms; first-line drug therapy
- **Verapamil/diltiazem**: negative lusitropic action → prolongs relaxation time; used in patients intolerant of beta-blockers
- **Disopyramide**: negative inotropic action plus antiarrhythmic effect; combined with a beta-blocker when LVOTO persists
- **Septal reduction therapy**: surgical septal myectomy or alcohol septal ablation (ASA); for drug-refractory LVOTO

## Model Files

| File | Description |
|------|------|
| [hcm_qsp_model.dot](../../../hypertrophic-cardiomyopathy/hcm_qsp_model.dot) | Graphviz mechanistic map source (approximately 110+ nodes / 12 clusters) |
| [hcm_qsp_model.svg](../../../hypertrophic-cardiomyopathy/hcm_qsp_model.svg) | SVG vector image (scalable) |
| [hcm_qsp_model.png](../../../hypertrophic-cardiomyopathy/hcm_qsp_model.png) | PNG image (150 dpi) |
| [hcm_mrgsolve_model.R](../../../hypertrophic-cardiomyopathy/hcm_mrgsolve_model.R) | mrgsolve ODE model (approximately 19 compartments / 5+ treatment scenarios) |
| [hcm_shiny_app.R](../../../hypertrophic-cardiomyopathy/hcm_shiny_app.R) | Shiny dashboard |
| [hcm_references.md](../../../hypertrophic-cardiomyopathy/hcm_references.md) | References (approximately 56 articles, PubMed links) |

## mrgsolve Model (ODE Model)

- **Compartment structure**: drug PK (mavacamten gut/central/peripheral, beta-blocker single compartment) + cellular signalling (intracellular Ca²⁺/SR Ca²⁺, calcineurin, NFAT nuclear translocation, ERK activation) + cardiac structure (IVS thickness, LV mass, TGF-β1, collagen, LVOT pressure gradient, LVEDP, HR, NT-proBNP, troponin I, atrial fibrillation risk) compartments
- **Key treatment scenarios**: ① untreated baseline, ② mavacamten monotherapy, ③ mavacamten + beta-blocker combination, ④ beta-blocker monotherapy, ⑤ aficamten-analogue scenario, ⑥+ additional dose-response analyses
- **Calibration/evidence**: qualitatively calibrated from EXPLORER-HCM (Olivotto 2020, Lancet) — reduction in LVOT pressure gradient and KCCQ improvement, VALOR-HCM (Desai 2023, JAMA Cardiol) — reduced eligibility for septal reduction, and SEQUOIA-HCM (Nagueh 2024, Lancet) data

## Shiny Dashboard

The dashboard comprises a patient profile tab (baseline LVOT pressure gradient, IVS thickness, NYHA class, choice of genetic variant), drug PK dynamics (mavacamten blood concentration), sarcomere/cellular signalling PD metrics (SRX fraction, Ca²⁺, calcineurin), cardiac structure clinical endpoints (LVOT pressure gradient, LV hypertrophy, fibrosis), NT-proBNP/troponin biomarkers, comparison of 5+ treatment scenarios, and atrial fibrillation/sudden death risk metrics.

## Usage

```r
library(mrgsolve)
mod <- mread("hcm_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny dashboard:
shiny::runApp("hcm_shiny_app.R")
```
```bash
# Render the mechanistic map
dot -Tsvg hcm_qsp_model.dot -o hcm_qsp_model.svg
```

## References

For detailed citations, see [hcm_references.md](../../../hypertrophic-cardiomyopathy/hcm_references.md) (approximately 56 articles).

---
*This model is a qualitative/semi-quantitative QSP model intended for education and research, and must not be used directly for clinical decision-making.*
