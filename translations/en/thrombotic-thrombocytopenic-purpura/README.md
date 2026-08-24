# Thrombotic Thrombocytopenic Purpura (TTP) — QSP Model

[![Disease](https://img.shields.io/badge/Disease-TTP-red)](.)
[![Category](https://img.shields.io/badge/Category-Autoimmune%2FHematology-orange)](.)
[![Model](https://img.shields.io/badge/Model-mrgsolve%20ODE-blue)](.)
[![Drugs](https://img.shields.io/badge/Drugs-Caplacizumab%20·%20Rituximab%20·%20TPE-green)](.)

---

## Disease Overview

**Thrombotic thrombocytopenic purpura (TTP)** is a rare thrombotic
microangiopathy (TMA) caused by autoantibodies against ADAMTS13 (a von
Willebrand factor-cleaving metalloproteinase).

- **Core pathophysiology**: autoantibody -> inhibition of ADAMTS13 activity
  -> accumulation of ultra-large VWF multimers (ULVWF) -> microthrombus
  formation -> thrombocytopenia + microangiopathic haemolytic anaemia
  (MAHA) + multi-organ ischaemia
- **Incidence**: roughly 3-7 per million per year; female:male ratio ~3:1;
  most common in the 30s-50s
- **Pre-treatment mortality**: ~90% -> **with current TPE-based
  treatment**: ~10-20%

### The Diagnostic Triad (Thrombotic Microangiopathy Triad)
1. **Thrombocytopenia** (PLT <30 x 10^9/L)
2. **MAHA** (Coombs-negative haemolytic anaemia: LDH up, schistocytes up,
   haptoglobin down)
3. **ADAMTS13 deficiency** (<10 U/dL = %)

---

## Mechanistic Map

[![TTP QSP Model](../../../thrombotic-thrombocytopenic-purpura/ttp_qsp_model.png)](../../../thrombotic-thrombocytopenic-purpura/ttp_qsp_model.svg)

*Click to open the vector SVG image*

### Cluster Structure (13 Subgraph Clusters, 130+ Nodes)

| Cluster | Key components |
|---------|-------------|
| ① ADAMTS13 biology | ADAMTS13 gene, hepatic synthesis, plasma activity, autoantibody, inhibitor titre |
| ② VWF biology | WPB secretion, ULVWF pool, A1/A2 domains, ADAMTS13 cleavage, normal multimers |
| ③ Platelet biology | Platelet count, GPIbα, GPIIb/IIIa, activation, microthrombus formation |
| ④ Endothelial cell activation | EC activation, WPB release, eNOS/PGI2, ET-1, tissue factor |
| ⑤ Immunology/inflammation | B cells, plasma cells, CD4 Tfh, germinal centre, IL-6, TNF-α, NF-κB |
| ⑥ Secondary coagulation | Tissue factor pathway, thrombin, fibrin, D-dimer, anticoagulant system |
| ⑦ Multi-organ injury | Kidney (Cr), brain (neuro score), heart (troponin), gastrointestinal |
| ⑧ MAHA | Mechanical haemolysis, schistocytes, LDH, haptoglobin, Hgb |
| ⑨ Therapeutic plasma exchange (TPE) | FFP/SD plasma, ADAMTS13 replacement, inhibitor removal, ULVWF removal |
| ⑩ Caplacizumab PK/PD | Two-compartment nanobody PK, VWF A1 blockade, platelet-recovery effect |
| ⑪ Rituximab PK/PD | TMDD two-compartment PK, CD20 binding, B-cell depletion, ADAMTS13 recovery |
| ⑫ Immunosuppressive therapy | Prednisolone, cyclosporine, MMF, bortezomib |
| ⑬ Clinical endpoints | PLASMIC score, response criteria, relapse risk, complete remission, mortality |

---

## mrgsolve ODE Model Specification

**File**: `ttp_mrgsolve_model.R`

### Compartment List (18 ODEs)

| # | Compartment | Biological meaning |
|---|------|-------------|
| 1-3 | `CAPLA_GUT`, `CAPLA_C`, `CAPLA_P` | Caplacizumab two-compartment PK (nanobody) |
| 4-5 | `RTX_C`, `RTX_P` | Rituximab two-compartment PK (TMDD) |
| 6 | `A13_ACT` | ADAMTS13 activity (U/dL = %) |
| 7 | `INH` | Inhibitor titre (Bethesda units, BU) |
| 8 | `ULVWF` | ULVWF pool (ng/mL) |
| 9 | `PLT` | Platelet count (x10^9/L) |
| 10 | `MT` | Microthrombus burden (AU) |
| 11 | `BC` | B cells (% of normal) |
| 12 | `PC` | Plasma cells (AU) |
| 13 | `AUTOAB` | Autoantibody (BU) |
| 14 | `LDH_AB` | LDH (IU/L) — haemolysis |
| 15 | `CREAT` | Creatinine (μmol/L) — renal function |
| 16 | `TROP` | Troponin I (ng/mL) — cardiac injury |
| 17 | `HGB` | Haemoglobin (g/dL) — MAHA anaemia |
| 18 | `PRED_C` | Prednisolone plasma concentration (ng/mL) |

### Key PK Parameters

| Drug | Model | CL (L/d) | V1 (L) | t½ |
|------|------|---------|--------|-----|
| Caplacizumab | Two-compartment | 0.50 | 3.0 | ~6-13h |
| Rituximab | Two-compartment + TMDD | 0.80 | 4.0 | ~14-20d |
| Prednisolone | One-compartment | 15.0 | 35.0 | ~3h |

### Treatment Scenarios (6)

| Scenario | Regimen | Clinical basis |
|---------|-------|----------|
| S0 | Untreated (natural history) | — (mortality ~90%) |
| S1 | **TPE alone** | Rock 1991 NEJM (pre-standard era) |
| S2 | **TPE + prednisolone** | Current standard of care |
| S3 | **TPE + caplacizumab + Pred** (HERCULES) | Scully 2019 NEJM |
| S4 | **TPE + rituximab + Pred** (TITAN) | Peyvandi 2016 NEJM |
| S5 | **Triple combination** (TPE + CAPLA + RTX + Pred) | Combined strategy |
| S6 | **Congenital TTP** (FFP prophylaxis q2wk) | Upshaw-Schulman syndrome |

---

## Shiny Dashboard Specification

**File**: `ttp_shiny_app.R`

| Tab | Function |
|----|------|
| ① Patient profile | Weight/height, biomarkers at presentation, treatment selection, PLASMIC score calculation |
| ② Pharmacokinetics (PK) | Caplacizumab/rituximab/prednisolone concentration-time curves |
| ③ ADAMTS13/VWF/PLT | Key PD variables: ADAMTS13 activity, ULVWF, platelets, microthrombi |
| ④ Clinical endpoints | LDH, Hgb, creatinine, troponin, composite TMA index |
| ⑤ Scenario comparison | Head-to-head comparison of 5 treatment strategies (with response table) |
| ⑥ Biomarker panel | Schistocyte %, autoantibody (BU), B-cell kinetics, summary table |

---

## References

**File**: `ttp_references.md` — **47** PubMed citations in total

### Key Clinical Trials

| Trial | Summary of results |
|---------|----------|
| **HERCULES** (Scully 2019 NEJM) | Caplacizumab + TPE: reduced time to platelet response (2.69 vs. 2.88 days), fewer exacerbations (3% vs. 28%) |
| **TITAN** (Peyvandi 2016 NEJM) | Rituximab x4: higher complete remission at 12 weeks (59% vs. 43%), fewer relapses |
| **Rock 1991 NEJM** | TPE vs. plasma infusion: established the superiority of TPE (mortality 22% vs. 37%) |

---

## Usage

```bash
# Graphviz rendering
dot -Tsvg ttp_qsp_model.dot -o ttp_qsp_model.svg
dot -Tpng -Gdpi=150 ttp_qsp_model.dot -o ttp_qsp_model.png
```

```r
# Run the mrgsolve model
install.packages(c("mrgsolve","dplyr","ggplot2","tidyr","patchwork"))
source("ttp_mrgsolve_model.R")

# Run the Shiny dashboard
install.packages(c("shiny","shinydashboard","plotly","DT"))
shiny::runApp("ttp_shiny_app.R")
```

---

## File List

| File | Description |
|------|------|
| `ttp_qsp_model.dot` | Graphviz mechanistic map source |
| `ttp_qsp_model.svg` | Vector mechanistic map |
| `ttp_qsp_model.png` | PNG mechanistic map (150 dpi) |
| `ttp_mrgsolve_model.R` | mrgsolve ODE model + 6 scenarios + visualisation |
| `ttp_shiny_app.R` | Interactive Shiny dashboard (6 tabs) |
| `ttp_references.md` | 47 PubMed references |
| `README.md` | This file |

---

*Generated by Claude Code QSP Routine | 2026-06-25*
