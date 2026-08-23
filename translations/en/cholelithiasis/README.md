# Cholelithiasis (CHOL) — QSP Model

> **QSP Disease Model Library** · A quantitative systems pharmacology (QSP) model automatically generated via the Claude Code Routine.
> Parent library → [../README.md](../README.md) · Category: Gastroenterology/Hepatobiliary

[![CHOL QSP Model](../../../cholelithiasis/chol_qsp_model.png)](../../../cholelithiasis/chol_qsp_model.svg)

## Overview

Cholelithiasis is a common digestive disease occurring in approximately 10~20% of adults worldwide, with a particularly high prevalence in Western countries and older women. Most (~80%) are cholesterol gallstones, in which crystal nuclei form in cholesterol-supersaturated bile, and gallbladder motility disorder (stasis) is the key factor promoting stone growth. Ursodeoxycholic acid (UDCA) is the only oral drug therapy that lowers the biliary cholesterol saturation index and dissolves stones, with statins and ezetimibe used in combination to further suppress cholesterol synthesis and absorption.

## Key Pathways

| Pathway | Key molecule/mechanism | Clinical outcome |
|------|----------------|-----------|
| Biliary cholesterol supersaturation | ABCG5/G8 overexpression, increased HMGCR activity | Raised cholesterol saturation index (CSI) |
| Promotion of nucleation | Secretion of mucin glycoproteins and pronucleating factors | Cholesterol crystal formation |
| Gallbladder stasis/motility disorder | Reduced CCK secretion, reduced smooth muscle contraction | Prolonged crystal residence time → stone growth |
| Enterohepatic circulation disruption | Reduced bile acid reabsorption, increased cholesterol/phospholipid ratio | Contraction of the bile acid pool |
| Chronic inflammation/mucosal response | Raised IL-6 and CRP, gallbladder wall thickening | Induces cholecystitis, pain |
| Statin effect | HMGCR inhibition → reduced hepatic cholesterol synthesis | Lowered biliary cholesterol concentration |
| UDCA mechanism | Inhibits cholesterol crystallisation, increases bile hydrophilicity | Promotes dissolution of small stones |

## Drug Targets

- **Ursodeoxycholic acid (UDCA)**: normalises bile acid composition, reduces CSI, dissolves small radiolucent stones (e.g. Actigall)
- **Statins (HMG-CoA reductase inhibitors)**: rosuvastatin, simvastatin — suppress hepatic cholesterol synthesis
- **Ezetimibe**: blocks NPC1L1 → inhibits intestinal cholesterol absorption
- **Bile acid sequestrants (cholestyramine)**: interrupt enterohepatic circulation, induce bile acid loss
- **CCK receptor agonist (sincalide)**: promotes gallbladder contraction, improves bile emptying

## Model Files

| File | Description |
|------|------|
| [chol_qsp_model.dot](../../../cholelithiasis/chol_qsp_model.dot) | Graphviz mechanistic map source (approximately 431 nodes / 10 clusters) |
| [chol_qsp_model.svg](../../../cholelithiasis/chol_qsp_model.svg) | SVG vector image (scalable) |
| [chol_qsp_model.png](../../../cholelithiasis/chol_qsp_model.png) | PNG image (150 dpi) |
| [chol_mrgsolve_model.R](../../../cholelithiasis/chol_mrgsolve_model.R) | mrgsolve ODE model (approximately 25 compartments / 5 treatment scenarios) |
| [chol_shiny_app.R](../../../cholelithiasis/chol_shiny_app.R) | Shiny dashboard |
| [chol_references.md](../../../cholelithiasis/chol_references.md) | References (approximately 46 articles, PubMed links) |

## mrgsolve Model (ODE Model)

- **Compartment structure**: separate absorption/plasma/hepatic/biliary/gallbladder PK compartments for UDCA, statin, and ezetimibe, plus PD compartments for bile acid pool, hepatic cholesterol, biliary cholesterol/phospholipid, gallbladder volume, crystal mass, stone volume, IL-6, and CRP
- **Key treatment scenarios**: (1) natural course (no treatment), (2) UDCA monotherapy, (3) statin monotherapy, (4) UDCA + statin combination, (5) UDCA + ezetimibe combination
- **Calibration/evidence**: parameters calibrated from the UDCA dissolution effect reported by Portincasa et al. (Lancet 2006) and relevant RCT data on the CSI-lowering effect of statins

## Shiny Dashboard

Comprises 6 tabs: (1) Patient Profile — sets risk factors (obesity, age, female sex, diet); (2) PK tab — UDCA/statin/ezetimibe plasma concentration time series; (3) Key PD Metrics — CSI, bile acid pool, gallbladder motility; (4) Clinical Endpoints — change in stone volume, dissolution rate; (5) Scenario Comparison — 1-year outcome comparison across treatment strategies; (6) Biomarkers — IL-6, CRP, biliary phospholipid/cholesterol ratio

## Usage

```r
library(mrgsolve)
mod <- mread("chol_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny dashboard:
shiny::runApp("chol_shiny_app.R")
```
```bash
# Render the mechanistic map
dot -Tsvg chol_qsp_model.dot -o chol_qsp_model.svg
```

## References

For detailed citations, see [chol_references.md](../../../cholelithiasis/chol_references.md) (approximately 46 articles).

---
*This model is a qualitative/semi-quantitative QSP model intended for education and research, and must not be used directly for clinical decision-making.*
