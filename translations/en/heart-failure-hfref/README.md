# Heart Failure with Reduced Ejection Fraction (HFrEF) — QSP Model

> **QSP Disease Model Library** · A quantitative systems pharmacology (QSP) model automatically generated via the Claude Code Routine.
> Parent library → [../README.md](../README.md) · Category: Cardiovascular

[![HFrEF QSP Model](../../../heart-failure-hfref/hfref_qsp_model.png)](../../../heart-failure-hfref/hfref_qsp_model.svg)

## Overview

Heart failure with reduced ejection fraction (HFrEF, EF < 40%) has as its core pathophysiology a vicious cycle in which neurohormonal (RAAS/sympathetic nervous system) overactivation accelerates ventricular remodelling. Approximately 26 million people worldwide are affected, and 5-year mortality reaches approximately 50%. AngII, aldosterone, and norepinephrine, activated after myocardial injury, initially maintain cardiac output, but with chronicity drive ventricular dilation, fibrosis, and hypertrophy. The four pillars of guideline-directed medical therapy (GDMT), namely ARNI (sacubitril-valsartan), beta-blockers, MRAs, and SGLT2 inhibitors, now markedly improve survival. Heart rate control (ivabradine) and ventricular assist devices (VAD) are additional treatment options.

## Key Pathways

| Pathway | Key molecule/mechanism | Clinical outcome |
|------|----------------|-----------|
| RAAS activation | Renal hypoperfusion → renin → AngI → ACE → AngII → AT1R → aldosterone | Sodium retention, vasoconstriction, myocardial hypertrophy |
| Sympathetic nervous system overactivation | Reduced baroreceptor response → increased NE secretion → β₁R → increased heart rate/contractility | Chronic ventricular overload, β1-receptor downregulation |
| Natriuretic peptide compensation | Raised LVEDP → BNP/ANP secretion → NEP degradation → diuresis, vasodilation | Elevated NT-proBNP |
| cGMP-PKG signalling | ANP/BNP → pGC, NO → sGC → cGMP → PKG → myocardial relaxation | cGMP amplified by ARNI |
| Ventricular remodelling | TGF-β → fibroblasts → ECM → ventricular fibrosis and dilation | Reduced LVEF, ventricular enlargement |
| Inflammatory activation | TNF-α/IL-6 → cardiomyocyte apoptosis, worsening ventricular function | Accelerated decline in cardiac function |

## Drug Targets

- **ARNI (sacubitril-valsartan)**: inhibits NEP (reduces BNP degradation → raises cGMP) plus AT1R blockade; reduced death and hospitalisation versus enalapril in PARADIGM-HF
- **Beta-blockers (carvedilol, metoprolol)**: block β₁/β₂R → reduce heart rate and energy expenditure → reverse remodelling; supported by COPERNICUS and MERIT-HF
- **MRAs (eplerenone, spironolactone)**: block the aldosterone receptor → antifibrotic and diuretic effect; supported by RALES and EMPHASIS-HF
- **SGLT2 inhibitors (dapagliflozin, empagliflozin)**: osmotic diuresis, optimisation of myocardial energetics, reduced ventricular load; supported by DAPA-HF and EMPEROR-Reduced
- **Ivabradine**: blocks HCN (If channel) → selectively reduces heart rate; supported by the SHIFT trial
- **Diuretics (furosemide)**: relieve congestive symptoms

## Model Files

| File | Description |
|------|------|
| [hfref_qsp_model.dot](../../../heart-failure-hfref/hfref_qsp_model.dot) | Graphviz mechanistic map source (approximately 100+ nodes / 9 clusters) |
| [hfref_qsp_model.svg](../../../heart-failure-hfref/hfref_qsp_model.svg) | SVG vector image (scalable) |
| [hfref_qsp_model.png](../../../heart-failure-hfref/hfref_qsp_model.png) | PNG image (150 dpi) |
| [hfref_mrgsolve_model.R](../../../heart-failure-hfref/hfref_mrgsolve_model.R) | mrgsolve ODE model (approximately 26 compartments / 5 treatment scenarios + dose-response analysis) |
| [hfref_shiny_app.R](../../../heart-failure-hfref/hfref_shiny_app.R) | Shiny dashboard |
| [hfref_references.md](../../../heart-failure-hfref/hfref_references.md) | References (approximately 62 articles, PubMed links) |

## mrgsolve Model (ODE Model)

- **Compartment structure**: includes RAAS (AngI, AngII, Ang1-7, aldosterone), sympathetic (NE), BNP/NT-proBNP/cGMP, haemodynamics (LVEDV, HR, SVR, LVEF), remodelling (TGF-β1, fibrosis, hypertrophy), inflammation (TNF-α, IL-6), and drug PK (LBQ657, valsartan, beta-blocker, MRA, SGLT2i, ivabradine) compartments (26 compartments in total)
- **Key treatment scenarios**: ① untreated baseline, ② ACEi + beta-blocker (former standard), ③ ARNI + BB + MRA (triple therapy), ④ ARNI + BB + MRA + SGLT2i (four-pillar GDMT), ⑤ maximal GDMT + ivabradine (five-drug regimen)
- **Calibration/evidence**: LVEF recovery, NT-proBNP reduction, and mortality qualitatively calibrated from PARADIGM-HF (McMurray 2014), DAPA-HF (McMurray 2019), EMPEROR-Reduced (Packer 2020), and SHIFT (Swedberg 2010) clinical trial data

## Shiny Dashboard

The dashboard comprises a patient profile tab (baseline LVEF, BNP, NYHA class, choice of aetiology), drug PK dynamics (ARNI/SGLT2i/BB concentrations), RAAS/sympathetic nervous system PD metrics, cardiac function clinical endpoints (LVEF, NT-proBNP, CO, HR), comparison of 5 treatment scenarios, dose-response analysis, and a biomarker tab (NT-proBNP, BNP, MAP).

## Usage

```r
library(mrgsolve)
mod <- mread("hfref_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny dashboard:
shiny::runApp("hfref_shiny_app.R")
```
```bash
# Render the mechanistic map
dot -Tsvg hfref_qsp_model.dot -o hfref_qsp_model.svg
```

## References

For detailed citations, see [hfref_references.md](../../../heart-failure-hfref/hfref_references.md) (approximately 62 articles).

---
*This model is a qualitative/semi-quantitative QSP model intended for education and research, and must not be used directly for clinical decision-making.*
