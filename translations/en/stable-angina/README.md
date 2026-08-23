# Stable Angina (SA) — QSP Model

> **QSP Disease Model Library** · A quantitative systems pharmacology (QSP) model automatically generated via the Claude Code Routine.
> Parent library → [../README.md](../README.md) · Category: Cardiovascular

[![SA QSP Model](../../../stable-angina/sa_qsp_model.png)](../../../stable-angina/sa_qsp_model.svg)

## Overview

Stable angina (chronic coronary syndrome, CCS) is predictable chest pain that occurs when atherosclerotic coronary stenosis causes myocardial oxygen demand to exceed supply, affecting approximately 110 million people worldwide. The core pathogenic mechanism is coronary stenosis driven by lipid deposition, inflammation, and plaque formation, with symptoms reproduced by an oxygen supply-demand imbalance during exercise or emotional stress. Anti-ischaemic therapy (beta-blockers, calcium channel blockers, nitrates, ivabradine, ranolazine), antiplatelet therapy (aspirin, clopidogrel), and statins/ACEIs/ARBs form the foundation of standard treatment.

## Key Pathways

| Pathway | Key molecule/mechanism | Clinical outcome |
|------|----------------|-----------|
| Atherosclerotic plaque formation | LDL oxidation → foam cells → lipid core/fibrous cap formation | Coronary stenosis, reduced coronary flow reserve |
| Myocardial oxygen supply-demand imbalance | Increased heart rate, contractility, and afterload → rising MVO₂ vs. supply limited by stenosis | Ischaemic chest pain, ST changes |
| Endothelial dysfunction | Decreased eNOS → NO deficiency → vasoconstriction, platelet activation | Increased coronary tone |
| Sympathetic activation | Catecholamines → β1 receptor → increased heart rate, contractility, and oxygen consumption | Lowered threshold for triggering anginal symptoms |
| Late sodium current (INaL) | Ischaemia → increased INaL → intracellular Ca²⁺ overload → diastolic tension | Increased myocardial oxygen consumption, arrhythmia |
| Platelet activation/aggregation | ADP, TXA2 → GP IIb/IIIa → risk of thrombus formation | Risk of conversion to acute coronary syndrome |
| Inflammation/lipid pathway | hsCRP, IL-6, ox-LDL → long-term plaque destabilisation mechanism | Increased risk of cardiovascular events |

## Drug Targets

- **Beta-blockers (metoprolol, bisoprolol, atenolol)**: β1 blockade → reduced heart rate and contractility → reduced MVO₂; improves symptoms and prognosis
- **Calcium channel blockers (amlodipine, diltiazem, verapamil)**: L-type Ca²⁺ channel blockade → vasodilation, reduced heart rate
- **Nitrates (isosorbide mononitrate/dinitrate, sublingual NTG)**: eNOS-independent NO donation → venous and coronary dilation
- **Ivabradine**: If channel blockade → pure heart rate reduction; an alternative to beta-blockers in patients maintaining sinus rhythm
- **Ranolazine**: late INaL inhibition → reduced Ca²⁺ overload; add-on anti-anginal therapy
- **Aspirin/clopidogrel**: COX-1/P2Y12 inhibition → suppressed platelet aggregation; secondary prevention of cardiovascular events
- **Statins (atorvastatin, rosuvastatin)**: HMG-CoA inhibition → reduced LDL and pleiotropic anti-inflammatory effects; plaque stabilisation

## Model Files

| File | Description |
|------|------|
| [sa_qsp_model.dot](../../../stable-angina/sa_qsp_model.dot) | Graphviz mechanistic map source (approximately 100+ nodes / 15 clusters) |
| [sa_qsp_model.svg](../../../stable-angina/sa_qsp_model.svg) | SVG vector image (scalable) |
| [sa_qsp_model.png](../../../stable-angina/sa_qsp_model.png) | PNG image (150 dpi) |
| [sa_mrgsolve_model.R](../../../stable-angina/sa_mrgsolve_model.R) | mrgsolve ODE model (approximately 22 compartments / approximately 11 scenarios) |
| [sa_shiny_app.R](../../../stable-angina/sa_shiny_app.R) | Shiny dashboard |
| [sa_references.md](../../../stable-angina/sa_references.md) | References (approximately 49 articles, PubMed links) |

## mrgsolve Model (ODE Model)

- **Compartment structure**: 1~2-compartment PK for beta-blockers/CCBs/nitrates/ivabradine/ranolazine + heart rate/blood pressure/MVO₂ dynamics module, an INaL inhibition effect compartment, anginal-episode frequency prediction, and a long-term LDL/plaque-area progression module
- **Key treatment scenarios**: no treatment, beta-blocker monotherapy, CCB monotherapy, nitrate, beta-blocker + CCB, ivabradine, ranolazine, aspirin + statin, combined optimal medical therapy (OMT), etc.
- **Calibration/evidence**: based on clinical data from COURAGE (OMT vs. PCI), BEAUTIFUL (ivabradine), MERLIN-TIMI 36 (ranolazine), and TNT (intensive atorvastatin therapy)

## Shiny Dashboard

Comprises 6 tabs: (1) **Patient Profile** — set coronary stenosis severity, baseline heart rate, blood pressure, and risk factors; (2) **PK Profile** — time course of anti-anginal drug blood concentrations; (3) **Key PD Measures** — dynamics of heart rate, blood pressure, and MVO₂ reduction; (4) **Clinical Endpoints** — changes in weekly anginal episode frequency and exercise tolerance duration; (5) **Scenario Comparison** — comparison of symptom control and LDL reduction across treatment strategies; (6) **Biomarkers** — LDL-C, hsCRP, blood glucose, NT-proBNP, ischaemic threshold.

## Usage

```r
library(mrgsolve)
mod <- mread("sa_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny dashboard:
shiny::runApp("sa_shiny_app.R")
```
```bash
# Render the mechanistic map
dot -Tsvg sa_qsp_model.dot -o sa_qsp_model.svg
```

## References

For detailed citations, see [sa_references.md](../../../stable-angina/sa_references.md) (approximately 49 articles).

---
*This model is a qualitative/semi-quantitative QSP model intended for education and research, and must not be used directly for clinical decision-making.*
