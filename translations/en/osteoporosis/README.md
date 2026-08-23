# Osteoporosis (OP) — QSP Model

> **QSP Disease Model Library** · A quantitative systems pharmacology (QSP) model automatically generated via the Claude Code Routine.
> Parent library → [../README.md](../README.md) · Category: Endocrine/Metabolic

[![OP QSP Model](../../../osteoporosis/op_qsp_model.png)](../../../osteoporosis/op_qsp_model.svg)

## Overview

Osteoporosis is a chronic metabolic bone disease in which reduced bone mass (T-score ≤ −2.5) and abnormal bone microarchitecture raise fracture risk, affecting approximately 30% of women and 12% of men aged 50 and older. Oestrogen deficiency (menopause) lowers the OPG/RANKL ratio, sharply increasing osteoclast activity, and a persistent imbalance between bone formation and resorption reduces BMD. Glucocorticoid-induced osteoporosis (GIOP) and excess parathyroid hormone (PTH) secretion are also important secondary causes. Bisphosphonates (suppress bone resorption), denosumab (anti-RANKL), and teriparatide/romosozumab (promote bone formation) are the current principal therapies.

## Key Pathways

| Pathway | Key molecule/mechanism | Clinical outcome |
|------|----------------|-----------|
| RANKL-RANK-OPG axis | Oestrogen↓ → RANKL↑/OPG↓ → increased osteoclast differentiation/activity | Promotes bone resorption, reduces BMD |
| Oestrogen deficiency | Reduced oestrogen receptor signalling → increased osteoblast apoptosis | Reduced bone formation, rapid postmenopausal bone loss |
| PTH signalling | Intermittent PTH↑ → cAMP → enhanced Wnt signalling → osteoblast proliferation | Bone formation promoted by teriparatide |
| Wnt/sclerostin | Sclerostin (a Wnt antagonist) → blocks LRP5/6 → suppresses osteoblast function | Reduced bone formation; the target of romosozumab |
| Glucocorticoid excess | GC → osteoblast apoptosis, RANKL↑ → osteoblast/osteoclast imbalance | GIOP (acute bone loss) |
| Calcium-vitamin D deficiency | PTH↑ → renal Ca reabsorption, reduced intestinal Ca absorption → secondary hyperparathyroidism | Osteomalacia, increased fracture risk |
| Bone turnover markers | CTX (osteoclast), P1NP (osteoblast), bsALP | Monitoring treatment response |

## Drug Targets

- **Alendronate/zoledronate (bisphosphonates)**: inhibit the osteoclast mevalonate pathway → osteoclast apoptosis → 50–70% reduction in bone resorption
- **Denosumab**: an anti-RANKL monoclonal antibody → blocks osteoclast differentiation (subcutaneous injection every 6 months)
- **Teriparatide (PTH1–34)**: intermittent PTH action → Wnt signalling → osteoblast proliferation/activity; increases BMD by more than 10%
- **Romosozumab**: anti-sclerostin → disinhibits Wnt signalling → dual effect of increased bone formation/decreased bone resorption (ARCH, FRAME trials)

## Model Files

| File | Description |
|------|------|
| [op_qsp_model.dot](../../../osteoporosis/op_qsp_model.dot) | Graphviz mechanistic map source (approximately 189 nodes / 10 clusters) |
| [op_qsp_model.svg](../../../osteoporosis/op_qsp_model.svg) | SVG vector image (scalable) |
| [op_qsp_model.png](../../../osteoporosis/op_qsp_model.png) | PNG image (150 dpi) |
| [op_mrgsolve_model.R](../../../osteoporosis/op_mrgsolve_model.R) | mrgsolve ODE model (approximately 22 compartments / 6 treatment scenarios) |
| [op_shiny_app.R](op_shiny_app.R) | Shiny dashboard |
| [op_references.md](../../../osteoporosis/op_references.md) | References (approximately 54 articles, PubMed links) |

## mrgsolve Model (ODE Model)

- **Compartment structure**: alendronate bone-binding (1 compartment), zoledronate bone-binding (1 compartment), denosumab SC+central (2 compartments), teriparatide central (1 compartment), romosozumab SC+central+peripheral (3 compartments); oestrogen, PTH, calcium, RANKL, OPG, sclerostin, osteoblast precursors, osteoblasts, osteoclast precursors, osteoclasts, BMD, CTX, P1NP, bsALP, and 10-year fracture risk
- **Key treatment scenarios**: S1 untreated postmenopausal, S2 alendronate 70 mg/week, S3 zoledronate 5 mg/year, S4 denosumab 60 mg/6 months, S5 teriparatide 20 μg/day, S6 romosozumab→denosumab sequential therapy
- **Calibration/evidence**: parameters referenced from the BMD data of the FIT (alendronate), HORIZON (zoledronate), FREEDOM (denosumab), and ARCH (romosozumab) clinical trials

## Shiny Dashboard

Comprises 6 tabs: ① **Patient Profile** (sets menopausal status, baseline BMD, T-score, and GIOP status), ② **PK** (plasma/bone drug concentration-time curves), ③ **Key PD Metrics** (RANKL/OPG ratio, osteoblast/osteoclast activity trends), ④ **Clinical Endpoints** (BMD g/cm², T-score, change in 10-year fracture risk), ⑤ **Scenario Comparison** (direct comparison of 6 treatment strategies), ⑥ **Biomarkers** (CTX, P1NP, bsALP trends).

## Usage

```r
library(mrgsolve)
mod <- mread("op_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny dashboard:
shiny::runApp("op_shiny_app.R")
```
```bash
# Render the mechanistic map
dot -Tsvg op_qsp_model.dot -o op_qsp_model.svg
```

## References

For detailed citations, see [op_references.md](../../../osteoporosis/op_references.md) (approximately 54 articles).

---
*This model is a qualitative/semi-quantitative QSP model intended for education and research, and must not be used directly for clinical decision-making.*
