# Diabetic Retinopathy QSP Model

> **Directory:** `diabetic-retinopathy/` | **Abbreviation:** DR | **Date:** 2026-06-24
> **Category:** Chronic disease / Ophthalmology / Diabetic complication

---

[![DR QSP mechanistic map](dr_qsp_model.png)](dr_qsp_model.svg)

---

## Disease Overview

**Diabetic retinopathy (DR)** is the most common microvascular complication
of diabetes and a leading cause of blindness in adults worldwide.

| Item | Content |
|------|------|
| **Prevalence** | About 34.6% of all diabetic patients (~146 million people) |
| **DME** | The leading cause of blindness in working-age adults; ~7% of diabetic patients |
| **PDR risk** | ~50% of severe NPDR progresses to PDR within 1 year |
| **Global burden** | ~224 million people with DR projected by 2045 (Teo 2021) |
| **Clinical outcome measures** | ETDRS visual acuity (letters), OCT CRT (µm), DR severity stage |

---

## Key Pathophysiological Pathways

### Stage 1: Hyperglycaemia → four biochemical pathways
| Pathway | Key enzyme/molecule | Result |
|------|--------------|------|
| **Polyol pathway** | Aldose reductase (AR) | Sorbitol↑, NADPH depletion, pseudohypoxia |
| **Hexosamine pathway** | GFAT, O-GlcNAc | Sp1 glycosylation → TGF-β, PAI-1 ↑ |
| **PKC pathway** | PKCβ1/β2/δ | VEGF↑, NF-κB↑, eNOS↓, ET-1↑ |
| **AGE-RAGE pathway** | Methylglyoxal | Protein cross-linking, basement-membrane thickening, NF-κB↑ |

### Stage 2: Oxidative-nitrosative stress
- Mitochondrial ETC overload → O₂•⁻ generation → PARP activation → GAPDH inhibition → amplifying feedback into pathways 1–4
- BH4 depletion → eNOS uncoupling → reduced NO + increased ONOO⁻ → apoptosis

### Stage 3: VEGF/neovascularisation
- HIF-1α (hypoxia) + NF-κB (inflammation) → VEGF-A165 overexpression
- VEGFR2 signalling → PI3K/AKT (permeability↑) + ERK (proliferation↑) + PLCγ (PKC)
- Increased Ang2 / Tie-2 destabilisation → loss of pericyte support → vascular fragility

### Stage 4: Vascular structural lesions
| Lesion | Stage | Mechanism |
|------|------|------|
| Microaneurysm | Early NPDR | Pericyte loss → vessel dilation |
| Hard exudate | NPDR | Lipid leakage |
| Cotton-wool spot | Severe NPDR | Nerve-fibre-layer ischaemia |
| IRMA | Severe NPDR | Capillary occlusion → collateral circulation |
| Neovascularisation (NVE/NVD) | PDR | VEGF-driven angiogenesis |
| Diabetic macular oedema (DME) | Any stage | BRB breakdown → intramacular fluid |

---

## QSP Model Files

| File | Description |
|------|------|
| [`dr_qsp_model.dot`](dr_qsp_model.dot) | Graphviz mechanistic map (source) |
| [`dr_qsp_model.svg`](dr_qsp_model.svg) | SVG vector image |
| [`dr_qsp_model.png`](dr_qsp_model.png) | PNG image (150 dpi) |
| [`dr_mrgsolve_model_en.R`](dr_mrgsolve_model_en.R) | mrgsolve ODE model |
| [`dr_shiny_app_en.R`](dr_shiny_app_en.R) | Shiny dashboard |
| [`dr_references_en.md`](dr_references_en.md) | 57 references |

---

## Model Specification

### Mechanistic map (DOT)
- **Node count:** 210+ (9 clusters)
- **Clusters:**
  1. Systemic Risk Factors
  2. Hyperglycaemia-driven Biochemical Pathways
  3. Oxidative-nitrosative Stress
  4. VEGF/Angiogenesis Signalling
  5. Neuroinflammation Pathways
  6. Retinal Vascular Pathology
  7. Retinal Neurodegeneration
  8. Drug PK/PD
  9. Clinical Endpoints

### mrgsolve ODE model
- **Compartments:** 18 ODEs
  - Drug PK: `DRUG_VIT`, `DRUG_CENT`, `DRUG_PERIPH`, `CORT_VIT`
  - Glucose: `BG`, `HBA1C`
  - VEGF: `VEGF_FREE`, `VEGF_BOUND`, `VEGF_PLANT` (PlGF)
  - Oxidative stress: `ROS`, `AGE`
  - Inflammation: `CYT`, `ICAM`
  - Cells: `PERICYTE`, `EC_COUNT`
  - Structure: `PERM`, `NV`, `CRT`
  - Vision: `VA`
- **Treatment scenarios:** 6
  | Scenario | Description | Supporting trial |
  |---------|------|-------------|
  | S0 | Untreated (poor glycaemic control) | DCCT control arm |
  | S1 | Glycaemic control alone (HbA1c → 7%) | DCCT/EDIC |
  | S2 | Aflibercept 2 mg IVT q4w×5→q8w | PROTOCOL T, PANORAMA |
  | S3 | Ranibizumab 0.5 mg IVT q4w | RISE/RIDE |
  | S4 | Faricimab 6 mg IVT q4w×4→q16w | TENAYA/LUCERNE |
  | S5 | Aflibercept + glycaemic control combined | Based on CLARITY + DCCT |

### Shiny dashboard (8 tabs)
| Tab | Contents |
|----|------|
| (1) Patient Profile | Patient parameter settings, automatic DR-stage classification, risk estimation |
| (2) Drug PK | Intravitreal drug concentration, PK parameters, AUC/Cmax |
| (3) VEGF / Angiogenesis | Free VEGF, drug-VEGF complex, NV index, pericytes |
| (4) Oxidative/Inflammation | ROS, AGE, cytokines, ICAM-1 |
| (5) Retinal Structure | CRT (OCT), vascular permeability, endothelial cell count |
| (6) Visual Outcomes | BCVA (ETDRS), VA change, outcome summary |
| (7) Scenario Comparison | Simultaneous plot comparing all 6 scenarios + endpoint table |
| (8) Biomarkers & About | Biomarker profile, model overview, calibration data |

---

## How to Run

```bash
# 1) Re-render the mechanistic map (requires Graphviz)
dot -Tsvg dr_qsp_model.dot -o dr_qsp_model.svg
dot -Tpng -Gdpi=150 dr_qsp_model.dot -o dr_qsp_model.png
```

```r
# 2) Run the mrgsolve model (requires R packages)
install.packages(c("mrgsolve", "dplyr", "ggplot2", "tidyr"))
source("dr_mrgsolve_model_en.R")

# 3) Run the Shiny dashboard
install.packages(c("shiny", "shinydashboard", "plotly", "DT"))
shiny::runApp("dr_shiny_app_en.R")
```

---

## Summary of Key Clinical Trials

| Trial | Drug | N | Duration | Primary outcome |
|--------|------|---|------|---------|
| **PROTOCOL T** | AFL vs RBZ vs Bev | 660 | 1 year | VA: +13.3/+11.2/+9.7 letters |
| **RISE/RIDE** | Ranibizumab | 382+382 | 2 years | VA: +10.9 letters |
| **CLARITY** | Aflibercept vs PRP | 232 | 1 year | VA: +3.3 letters (AFL superior) |
| **PANORAMA** | Aflibercept | 402 | 2 years | 2-step improvement: 65% vs 15% |
| **TENAYA** | Faricimab | 331+327 | 1 year | VA: +5.8 letters; CRT −189µm |
| **LUCERNE** | Faricimab | 330+338 | 1 year | VA: +6.6 letters; CRT −194µm |
| **DCCT** | Intensive insulin therapy | 1,441 | 6.5 years | 76% reduction in new-onset DR |

---

## References
A total of **57** PubMed citations → [`dr_references_en.md`](dr_references_en.md)

- Epidemiology/classification (4)
- Hyperglycaemic pathways (6)
- VEGF/angiogenesis (6)
- Anti-VEGF clinical trials (8)
- Glycaemic control (4)
- Oxidative stress/AGE (4)
- Neurodegeneration (3)
- PK/QSP modelling (5)
- Inflammasome/inflammation (3)
- Ang/Tie2 signalling (4)
- OCT biomarkers (3)
- Steroid therapy (2)
- Neuroprotection/novel therapeutics (3)
- GLP-1 agonists (2)
