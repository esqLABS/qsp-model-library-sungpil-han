# Lambert-Eaton Myasthenic Syndrome (LEMS) QSP Model

> **Category**: Autoimmune Neuromuscular Disease  
> **Abbreviation**: LEMS  
> **Model version**: v1.0 | 2026-06-27

---

## Disease Overview

**Lambert-Eaton Myasthenic Syndrome** is a rare autoimmune disease caused by autoantibodies that target the P/Q-type voltage-gated calcium channel (VGCC) at the presynaptic neuromuscular junction (NMJ).

### Core pathophysiology

```
SCLC tumour expressing VGCC
        ↓  (molecular mimicry)
  Generation of anti-VGCC IgG autoantibodies
        ↓  (P/Q-type VGCC blockade and internalisation)
  Reduced presynaptic Ca²⁺ influx
        ↓
  Reduced ACh vesicle fusion and release
        ↓
  Reduced endplate potential (EPP) amplitude → Safety Factor < 1
        ↓
  Failed muscle action potential → weakness / reduced CMAP
```

### Key clinical features
| Feature | Description |
|------|------|
| **Proximal lower-limb weakness** | Predominates in the legs over the arms |
| **Areflexia/hyporeflexia** | Reduced or absent tendon reflexes |
| **Autonomic dysfunction** | Dry mouth, constipation, erectile dysfunction, orthostatic hypotension |
| **Facilitation** | Transient improvement in strength and CMAP after brief exercise |
| **Anti-VGCC antibody** | Positive autoantibody against P/Q-type VGCC (>85%) |

### Epidemiology
- Prevalence: approximately 3.4 per million
- Approximately 50-60%: associated with small cell lung cancer (SCLC) (paraneoplastic LEMS)
- The remainder: non-tumour autoimmune LEMS
- Associated with HLA-B8, DR3 (non-tumour LEMS)

---

## Mechanistic Map

[![LEMS QSP mechanistic map](lems_qsp_model.png)](lems_qsp_model.svg)

> Click to open the high-resolution vector-format (SVG) map.

### Subgraphs included (12 clusters, 100+ nodes)

| Cluster | Content |
|----------|------|
| SCLC tumour | Tumour growth, VGCC antigen expression, paraneoplastic mechanism |
| Oncology treatment | Chemotherapy, immuno-oncology agents, radiotherapy |
| Immune system activation | DC, CD4+ Th2, Tfh, germinal centre response, B cells, plasma cells |
| VGCC autoantibody | Antibody kinetics, FcRn recycling, complement activation |
| Presynaptic terminal | SNARE complex, ACh vesicle pool, Ca²⁺ dynamics, facilitation |
| NMJ/postsynapse | nAChR, EPP, safety factor, muscle AP, contraction |
| Autonomic nervous system | Sympathetic/parasympathetic VGCC blockade → autonomic dysfunction |
| Amifampridine PK | 2-compartment PK (absorption, central, peripheral) |
| Immunosuppressant PK | Prednisolone, azathioprine, MMF PK |
| Amifampridine PD | K⁺ channel blockade → AP prolongation → increased Ca²⁺ influx |
| Immunosuppression PD | GR activation, NF-κB inhibition, B-cell suppression |
| Clinical endpoints | CMAP, QMG score, MRC, facilitation ratio, safety factor |

---

## mrgsolve ODE Model (`lems_mrgsolve_model.R`)

### Model compartments (15 ODE compartments)

```r
CMT: A_gut, A_central, A_periph          # Amifampridine PK (2-compartment)
     P_gut, P_central, P_periph          # Prednisolone PK (2-compartment)
     Ab_VGCC                             # Anti-VGCC antibody kinetics
     VGCC_free                           # Functional VGCC fraction
     RRP                                 # ACh readily releasable vesicle pool
     EPP_amp                             # Endplate potential amplitude
     CMAP                                # CMAP amplitude
     QMG                                 # QMG score
     Bcell                               # B-cell dynamics
     Tumor                               # SCLC tumour size
     Facil                               # Facilitation state variable
```

### Summary of key equations

| Equation | Description |
|--------|------|
| `VGCC_free' = k_recov×(1-f) - k_block×Ab×f` | Antibody-dependent VGCC blockade |
| `Ab_VGCC' = kin×Bcell×(1-Imax_pred) - kout×Ab` | Antibody production/clearance |
| `ACh_release ∝ Ca_pre^Hill_Ca` | Ca²⁺-dependent vesicle release (Hill equation) |
| `CMAP = f(EPP/EPP_thresh)` | Safety-factor-based CMAP |
| `K_block = Emax×C^n/(EC50^n + C^n)` | Amifampridine PD (Emax model) |

### Treatment scenarios (6)

1. **Natural course (untreated)**: high anti-VGCC antibody levels → sustained CMAP decline
2. **Amifampridine monotherapy** (15 mg TID): rapid CMAP improvement (within 1-2 days)
3. **Prednisolone monotherapy** (40 mg/day): slow antibody reduction, long-term CMAP improvement
4. **Combination therapy** (amifampridine + prednisolone): early and long-term effects combined
5. **Plasma exchange + amifampridine**: rapid antibody removal, followed by amifampridine for maintenance
6. **Paraneoplastic LEMS — chemotherapy + amifampridine**: tumour shrinkage → reduced antigen → reduced antibody

---

## Shiny App (`lems_shiny_app.R`)

### Tab layout (7 tabs)

| Tab | Key features |
|----|-----------|
| **Patient profile** | Set initial anti-VGCC antibody level, paraneoplastic status, number of PE sessions |
| **Drug PK** | Amifampridine/prednisolone PK profiles, K⁺ channel blockade Emax curve |
| **VGCC & antibody** | Antibody kinetics, VGCC functional status, antibody-CMAP scatter plot |
| **NMJ / CMAP** | CMAP time series, safety factor, EPP amplitude |
| **Clinical endpoints** | QMG score, facilitation phenomenon, clinical response summary table |
| **Scenario comparison** | Simultaneous visual comparison of 6 treatment scenarios |
| **Biomarkers** | Antibody-CMAP relationship, dose-response, B-cell dynamics, biomarker summary table |

### How to run

```r
library(shiny)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(plotly)
library(DT)
library(shinydashboard)

shiny::runApp("lems_shiny_app.R")
```

---

## Parameter Calibration

| Parameter | Value | Reference |
|----------|-----|----------|
| Amifampridine Ka | 0.693 h⁻¹ (Tmax ~1h) | Haroldsen et al. (2015) |
| Amifampridine CL | 18.0 L/h | Oh et al. (2009) |
| Amifampridine EC50 (K⁺ blockade) | 120 ng/mL | Sanders et al. (2018) |
| IgG half-life | ~14 days (kout=0.05 h⁻¹) | Standard FcRn model value |
| VGCC blockade rate | 0.001 L/pmol/h | Nagel et al. (1988) |
| EPP safety-factor threshold | 1.0 (EPP/threshold) | Wood & Slater (2001) |
| CMAP (normal) | 5.0 mV | Standard electrodiagnostic range |
| QMG maximum | 39 points | LEMS clinical trial scale |

---

## File Index

| File | Description |
|------|------|
| `lems_qsp_model.dot` | Graphviz mechanistic map source (12 clusters, 100+ nodes) |
| `lems_qsp_model.svg` | Vector-format mechanistic map |
| `lems_qsp_model.png` | Raster-format mechanistic map (150 dpi) |
| `lems_mrgsolve_model.R` | mrgsolve ODE PK/PD model (15 compartments, 6 scenarios) |
| `lems_shiny_app.R` | Shiny interactive dashboard (7 tabs) |
| `lems_references_en.md` | References (56 PubMed links) |
| `README.md` | This file |

---

## Abbreviations

| Abbreviation | Description |
|------|------|
| LEMS | Lambert-Eaton Myasthenic Syndrome |
| VGCC | Voltage-Gated Calcium Channel |
| SCLC | Small Cell Lung Cancer |
| NMJ | Neuromuscular Junction |
| CMAP | Compound Muscle Action Potential |
| EPP | Endplate Potential |
| DAP | Diaminopyridine (3,4-DAP = amifampridine) |
| QMG | Quantitative Myasthenia Gravis score |
| RRP | Readily Releasable Pool (ACh vesicles) |
| PK/PD | Pharmacokinetics/Pharmacodynamics |
| PE | Plasma Exchange |
| IVIG | Intravenous Immunoglobulin |
| GR | Glucocorticoid Receptor |
| FcRn | Neonatal Fc Receptor |
| RNS | Repetitive Nerve Stimulation |

---

*Model generated: Claude Code Routine (CCR) | Date: 2026-06-27*
