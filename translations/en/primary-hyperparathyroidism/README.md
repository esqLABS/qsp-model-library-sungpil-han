# Primary Hyperparathyroidism (PHPT) — QSP Model

> **QSP Disease Model Library** · A quantitative systems pharmacology (QSP) model automatically generated via the Claude Code Routine.
> Parent library → [../README.md](../README.md) · Category: Endocrine/Metabolic

[![PHPT QSP Model](../../../primary-hyperparathyroidism/phpt_qsp_model.png)](../../../primary-hyperparathyroidism/phpt_qsp_model.svg)

## Overview
Primary hyperparathyroidism (PHPT) is an endocrine disease in which autonomous excess PTH secretion from a parathyroid adenoma (85%), hyperplasia, or, rarely, carcinoma causes chronic hypercalcaemia. Its prevalence is 1~2 per 1,000 population, and it is most common in women in their 50s~60s. Loss of the calcium-sensing receptor's (CaSR) PTH-calcium negative feedback allows sustained PTH secretion, and PTH drives bone loss (particularly cortical bone) and nephrocalcinosis/nephrolithiasis through renal calcium reabsorption, 1,25-vitamin D activation, and osteoclast stimulation (the RANK-RANKL pathway). The curative treatment is parathyroidectomy, while the calcimimetic cinacalcet and denosumab are used as drug therapy when surgery is not possible or in asymptomatic PHPT.

## Key Pathways
| Pathway | Key molecule/mechanism | Clinical outcome |
|------|----------------|-----------|
| PTH-CaSR feedback disruption | Reduced CaSR sensitivity, autonomous adenoma secretion | PTH excess, hypercalcaemia |
| PTH-bone pathway | Increased RANKL → osteoclast activation → cortical bone loss | Reduced BMD, fracture risk |
| PTH-renal pathway | Increased distal tubular calcium reabsorption, FGF23 regulation | Hypercalciuria, nephrolithiasis, reduced renal function |
| Vitamin D activation pathway | PTH → CYP27B1 activation → increased 1,25-(OH)₂D | Increased intestinal calcium absorption |
| Bone remodelling imbalance pathway | Uncoupled osteoblast (OB)/osteoclast (OC) activity | Osteitis fibrosa (in severe cases) |
| FGF23-Klotho axis | Indirect FGF23 regulation by PTH | Phosphate regulation abnormality |

## Drug Targets
- **Cinacalcet (calcimimetic)**: a CaSR positive allosteric modulator → suppresses PTH secretion, reduces serum calcium (FDA approved, for patients unfit for surgery)
- **Denosumab**: an anti-RANKL monoclonal antibody → suppresses osteoclast activity → preserves BMD
- **Bisphosphonates (alendronate, etc.)**: induce osteoclast apoptosis → protect cortical bone
- **Parathyroidectomy**: removal of the PTH-hypersecreting adenoma — the only curative treatment
- **Fluid administration and loop diuretics**: symptomatic management of acute hypercalcaemia (in the inpatient setting)

## Model Files
| File | Description |
|------|------|
| [phpt_qsp_model.dot](../../../primary-hyperparathyroidism/phpt_qsp_model.dot) | Graphviz mechanistic map source (approximately 128 nodes / 12 clusters) |
| [phpt_qsp_model.svg](../../../primary-hyperparathyroidism/phpt_qsp_model.svg) | SVG vector image (scalable) |
| [phpt_qsp_model.png](../../../primary-hyperparathyroidism/phpt_qsp_model.png) | PNG image (150 dpi) |
| [phpt_mrgsolve_model.R](../../../primary-hyperparathyroidism/phpt_mrgsolve_model.R) | mrgsolve ODE model (approximately 20 compartments / 8 treatment scenarios) |
| [phpt_shiny_app.R](../../../primary-hyperparathyroidism/phpt_shiny_app.R) | Shiny dashboard |
| [phpt_references.md](phpt_references.md) | References (approximately 62 articles, PubMed links) |

## mrgsolve Model (ODE Model)
- **Compartment structure**: PK compartments for cinacalcet (1 compartment), denosumab (2 compartments + RANKL binding), and alendronate (bone binding), plus PD compartments for PTH, serum Ca, serum PO4, vitamin D (25 and 1,25), osteoblasts, osteoclasts, RANKL, urinary Ca, GFR, lumbar spine BMD, and femoral neck BMD
- **Key treatment scenarios**: ① normal (healthy), ② untreated mild PHPT, ③ untreated severe PHPT, ④ cinacalcet 60 mg/day, ⑤ denosumab 60 mg q6mo, ⑥ parathyroidectomy (day 90), ⑦ cinacalcet + denosumab, ⑧ PHPT with concurrent CKD + cinacalcet 90 mg
- **Calibration/evidence**: serum Ca/PTH normalisation time calibrated with reference to the international PHPT guidelines (5th International Workshop, 2022) and cinacalcet clinical trial data (SHOPPE, PRIMARY, etc.)

## Shiny Dashboard
Comprises tabs for patient profile (presence of symptoms, osteoporosis grade, nephrolithiasis history), cinacalcet/denosumab PK profile, the PTH-calcium-vitamin D PD axis, BMD and fracture risk clinical endpoints, treatment scenario comparison (serum Ca/BMD change), and biomarkers (PTH, 24-hour urinary calcium, GFR).

## Usage
```r
library(mrgsolve)
mod <- mread("phpt_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny dashboard:
shiny::runApp("phpt_shiny_app.R")
```
```bash
# Render the mechanistic map
dot -Tsvg phpt_qsp_model.dot -o phpt_qsp_model.svg
```

## References
For detailed citations, see [phpt_references.md](phpt_references.md) (approximately 62 articles).

---
*This model is a qualitative/semi-quantitative QSP model intended for education and research, and must not be used directly for clinical decision-making.*
