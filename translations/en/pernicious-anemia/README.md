# Pernicious Anemia (PNA) — QSP Model

> **QSP Disease Model Library** · A quantitative systems pharmacology (QSP) model automatically generated via the Claude Code Routine.
> Parent library → [../README.md](../README.md) · Category: Haematology

[![PNA QSP Model](pna_qsp_model.png)](pna_qsp_model.svg)

## Overview
Pernicious anaemia is a megaloblastic anaemia arising when autoimmune gastritis causes intrinsic factor (IF) deficiency, blocking absorption of vitamin B12 (cobalamin). It affects approximately 0.1% of the general population, is more common in people over 60 and in women, and is associated with an HLA-DR3/DR4 predisposition. Anti-intrinsic factor antibodies (blocking antibodies, 70%) and anti-parietal cell antibodies (APC-Ab, 90%) are the core pathogenic mechanisms, leading from destruction of gastric body parietal cells → IF deficiency → failure of B12 absorption → haematopoietic impairment and neurological complications (subacute combined degeneration of the spinal cord). Intramuscular vitamin B12 injection or high-dose oral supplementation is the foundation of treatment.

## Key Pathways
| Pathway | Key molecule/mechanism | Clinical outcome |
|------|----------------|-----------|
| Autoimmune gastritis pathway | APC-Ab/IF-Ab, autoreactive CD4 Th1 cells | Parietal cell destruction, achlorhydria |
| Intrinsic factor deficiency pathway | Loss of IF secretion → cobalamin-IF complex cannot form | Blocked B12 absorption in the terminal ileum |
| B12 systemic distribution pathway | Transcobalamin II (TC-II)-mediated transport, enterohepatic circulation | Reduced serum B12/holotranscobalamin |
| Megaloblastic anaemia pathway | Impaired DNA synthesis (thymidylate pathway), abnormal cell division | Increased MCV, reduced Hb, reticulocyte response |
| Neurotoxicity pathway | Methylmalonic acid/homocysteine accumulation, myelin damage | Subacute combined degeneration, peripheral neuropathy |
| Oral passive absorption pathway | IF-independent passive diffusion (~1%/dose) | Basis for high-dose oral therapy efficacy |

## Drug Targets
- **Intramuscular vitamin B12 (IM)**: cyanocobalamin/hydroxocobalamin — rapid replacement, preferred when neurological recovery is a priority
- **High-dose oral vitamin B12**: 1,000 µg/day — exploits passive absorption (1%), good adherence
- **Sublingual/intranasal formulations**: alternative delivery routes — bypass IF deficiency
- **Immunosuppression (experimental)**: short-term corticosteroids — suppress antibodies, limited clinical use
- **Folate supplementation**: used to correct a concurrent deficiency (caution: given alone it can mask neurological complications)

## Model Files
| File | Description |
|------|------|
| [pna_qsp_model.dot](pna_qsp_model.dot) | Graphviz mechanistic map source (approximately 154 nodes / 11 clusters) |
| [pna_qsp_model.svg](pna_qsp_model.svg) | SVG vector image (scalable) |
| [pna_qsp_model.png](pna_qsp_model.png) | PNG image (150 dpi) |
| [pna_mrgsolve_model.R](pna_mrgsolve_model.R) | mrgsolve ODE model (approximately 18 compartments / 5 treatment scenarios) |
| [pna_shiny_app.R](pna_shiny_app.R) | Shiny dashboard |
| [pna_references.md](pna_references.md) | References (approximately 35 articles, PubMed links) |

## mrgsolve Model (ODE Model)
- **Compartment structure**: compartments for the IM depot, oral GI absorption, IF pool, portal B12, plasma, holotranscobalamin, liver, bone marrow, and neural tissue, plus PD compartments for autoantibodies, gastric parietal cells, Hb, MCV, reticulocytes, and a neurological complication index
- **Key treatment scenarios**: ① no treatment, ② standard IM regimen (1,000 µg/day for 1 week, then monthly), ③ IM maintenance therapy, ④ high-dose oral therapy, ⑤ intensive aggressive replacement therapy
- **Calibration/evidence**: serum B12/MCV normalisation time and the course of neurological complication recovery calibrated with reference to the literature (Stabler 2013, Carmel 2008)

## Shiny Dashboard
Comprises tabs for patient profile (antibody positivity, gastritis stage, presence of neurological symptoms), B12 pharmacokinetic PK profile, haematology PD metrics (Hb, MCV, reticulocytes), neurological complication progression, route-of-administration scenario comparison, and biomarkers (holotranscobalamin, methylmalonic acid, homocysteine).

## Usage
```r
library(mrgsolve)
mod <- mread("pna_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny dashboard:
shiny::runApp("pna_shiny_app.R")
```
```bash
# Render the mechanistic map
dot -Tsvg pna_qsp_model.dot -o pna_qsp_model.svg
```

## References
For detailed citations, see [pna_references.md](pna_references.md) (approximately 35 articles).

---
*This model is a qualitative/semi-quantitative QSP model intended for education and research, and must not be used directly for clinical decision-making.*
