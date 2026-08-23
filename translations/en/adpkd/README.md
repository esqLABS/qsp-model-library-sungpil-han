# Autosomal Dominant Polycystic Kidney Disease (ADPKD) — QSP Model

> **QSP Disease Model Library** · A quantitative systems pharmacology (QSP) model automatically generated via the Claude Code Routine.
> Parent library → [../README.md](../README.md) · Category: Renal/Urological

[![ADPKD QSP Model](../../../adpkd/adpkd_qsp_model.png)](../../../adpkd/adpkd_qsp_model.svg)

## Overview

Autosomal dominant polycystic kidney disease (ADPKD) is the most common monogenic kidney disease, caused by mutations in PKD1 (polycystin-1) or PKD2 (polycystin-2), occurring in approximately 1 in 1,000 live births. Loss of polycystin function lowers intracellular Ca²⁺ and raises cAMP, which promotes proliferation of cystic epithelial cells and fluid secretion, so that cysts multiply and expand from tens to thousands in number. Total kidney volume (TKV) is the key surrogate marker of disease progression, increasing by approximately 5–6% per year, with TKV > 750 mL (kidney length > 16.5 cm) defining a rapid progressor. Untreated, the median time to progression to ESRD is 58 years for PKD1 mutations and 79 years for PKD2 mutations. Tolvaptan (a V2 receptor antagonist), which inhibits the cAMP pathway, is the first approved drug shown to significantly slow the rise in TKV and the decline in eGFR.

## Key Pathways

| Pathway | Key molecule/mechanism | Clinical outcome |
|------|----------------|-----------|
| cAMP excess | Ca²⁺↓ → AC activation → cAMP↑ | Cell proliferation, CFTR-mediated fluid secretion |
| mTOR pathway | Loss of PC1 → mTOR overactivation → cell growth | Cystic epithelial proliferation |
| RAAS activation | Cyst-induced ischaemia → increased renin | Hypertension, glomerular hypertension |
| AVP/V2 receptor | Collecting duct V2R → cAMP → AQP2 expression and cyst fluid secretion | Accelerated increase in TKV |
| Fibrosis | TGF-β/PDGF → periscytic fibrosis | Loss of functional nephrons |
| Intracystic pressure | Cyst expansion → compression of adjacent vessels/nephrons | Reduced eGFR |
| Angiotensin II | Promotes oxidative stress and apoptosis | Worsening renal function |

## Drug Targets

- **Tolvaptan (V2R antagonist)**: demonstrated an approximately 49% reduction in TKV growth rate and delayed eGFR decline in the TEMPO 3:4 and REPRISE trials; requires hepatotoxicity monitoring
- **Everolimus (mTOR inhibitor)**: has a TKV-suppressing effect but renoprotection has not been demonstrated; limited clinical use
- **Somatostatin analogue (octreotide LAR)**: suppresses growth of hepatic and renal cysts; more effective in polycystic liver disease
- **ACE inhibitors/ARBs**: blood pressure control and renoprotection, currently the standard of care
- **Targets under investigation**: MEK inhibitors, CFTR inhibitors, antifibrotic agents

## Model Files

| File | Description |
|------|------|
| [adpkd_qsp_model.dot](../../../adpkd/adpkd_qsp_model.dot) | Graphviz mechanistic map source (approximately 184 nodes / 10 clusters) |
| [adpkd_qsp_model.svg](../../../adpkd/adpkd_qsp_model.svg) | SVG vector image (scalable) |
| [adpkd_qsp_model.png](../../../adpkd/adpkd_qsp_model.png) | PNG image (150 dpi) |
| [adpkd_mrgsolve_model.R](../../../adpkd/adpkd_mrgsolve_model.R) | mrgsolve ODE model (approximately 19 compartments / multiple treatment scenarios) |
| [adpkd_shiny_app.R](../../../adpkd/adpkd_shiny_app.R) | Shiny dashboard |
| [adpkd_references.md](../../../adpkd/adpkd_references.md) | References (approximately 53 articles, PubMed links) |

## mrgsolve Model (ODE Model)

- **Compartment structure**: tolvaptan (oral, 3-compartment PK) + everolimus (2 compartments) + octreotide LAR depot + ACEi compartment + PD compartments for AVP level, collecting duct cAMP, mTOR activity, angiotensin II, and blood pressure + outcome compartments for TKV, eGFR, urine osmolality, and functional nephron fraction
- **Key treatment scenarios**: ① no treatment, ② tolvaptan monotherapy (low risk), ③ tolvaptan monotherapy (high risk/rapid progression), ④ everolimus + tolvaptan combination, ⑤ ACEi + tolvaptan, ⑥ octreotide combination
- **Calibration/evidence**: parameters referenced from the TEMPO 3:4 clinical trial (TKV, eGFR data), the REPRISE study, and Torres et al. NEJM

## Shiny Dashboard

The dashboard comprises a patient profile tab (PKD1/PKD2, baseline TKV, eGFR, age), a tolvaptan PK and V2R occupancy tab, a TKV growth curve tab, an eGFR decline prediction tab, a treatment scenario comparison tab, and a urine osmolality/biomarker tab.

## Usage

```r
library(mrgsolve)
mod <- mread("adpkd_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny dashboard:
shiny::runApp("adpkd_shiny_app.R")
```
```bash
# Render the mechanistic map
dot -Tsvg adpkd_qsp_model.dot -o adpkd_qsp_model.svg
```

## References

For detailed citations, see [adpkd_references.md](../../../adpkd/adpkd_references.md) (approximately 53 articles).

---
*This model is a qualitative/semi-quantitative QSP model intended for education and research, and must not be used directly for clinical decision-making.*
