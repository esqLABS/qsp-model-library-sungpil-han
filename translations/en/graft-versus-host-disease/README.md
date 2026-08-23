# Graft-versus-Host Disease (GvHD) QSP Model

[![Disease](https://img.shields.io/badge/Disease-GvHD-red)](../../../graft-versus-host-disease/) [![Category](https://img.shields.io/badge/Category-Transplant%20Immunology-blue)](../../../graft-versus-host-disease/) [![Drugs](https://img.shields.io/badge/Drugs-CsA%20%7C%20TAC%20%7C%20RUX%20%7C%20BELU%20%7C%20MMF-orange)](../../../graft-versus-host-disease/)

## Disease Overview

**Graft-versus-Host Disease (GvHD)** is a severe complication following allogeneic haematopoietic stem cell transplantation (Allo-HSCT) in which donor immune cells attack the recipient's tissues. Acute GvHD (aGvHD) typically occurs within 100 days of transplantation, while chronic GvHD (cGvHD) usually develops thereafter and can persist for years.

---

## Mechanistic Map

[![GvHD Mechanistic Map](../../../graft-versus-host-disease/gvhd_qsp_model.png)](../../../graft-versus-host-disease/gvhd_qsp_model.svg)

*Click the image to open the full SVG with interactive zoom*

### Key Clusters

| # | Cluster | Key content |
|---|---------|---------|
| 1 | **HSCT Context & Conditioning** | MAC/RIC conditioning, TBI, tissue injury, DAMP/PAMP release, barrier damage |
| 2 | **Antigen Presentation & Priming** | Host DC, MHC mismatch, direct/indirect alloreactivity, CD28-B7 co-stimulation |
| 3 | **Donor T Cell Differentiation** | Th1/Th17/Treg/CD8 differentiation, NFAT/NF-κB/JAK-STAT/ROCK2 signalling |
| 4 | **Cytokine Network** | TNF-α, IFN-γ, IL-6, IL-17A, IL-10, TGF-β, BAFF |
| 5 | **Target Organ: Skin** | Keratinocyte apoptosis, lichenoid/sclerotic lesions, mLSS score |
| 6 | **Target Organ: Gut** | Intestinal epithelial cell apoptosis, crypt damage, ST2/REG3α biomarkers |
| 7 | **Target Organ: Liver** | Bile duct injury, cholestasis, Glucksberg grading |
| 8 | **Target Organ: Lung** | Bronchiolitis obliterans syndrome (BOS), reduced FEV1, CLAD score |
| 9 | **B Cell Pathology** | Tfh-B germinal-centre response, autoantibodies, BTK pathway |
| 10 | **Fibrosis (TGF-β/ROCK2)** | SMAD2/3, EMT, myofibroblast activation, ROCK2/IRF4 |
| 11 | **Drug PK: CNI** | CsA 2-compartment, TAC 2-compartment, CYP3A4/5 metabolism |
| 12 | **Drug PK/PD: Ruxolitinib** | JAK1/2 inhibition, STAT3/5 blockade, Treg expansion |
| 13 | **Other Drugs** | Steroids (NF-κB), Belumosudil (ROCK2), MMF (IMPDH), Ibrutinib (BTK) |
| 14 | **Clinical Endpoints & Biomarkers** | Glucksberg/NIH scores, ORR, FFS, OS, NRM |

**Total nodes: ~130+ | Subgraph clusters: 14**

---

## mrgsolve ODE Model (Pharmacokinetic/Pharmacodynamic Model)

**File**: `gvhd_mrgsolve_model.R`

### Compartments (32 total)

**Pharmacokinetics (PK) — 16 compartments:**
| Drug | Compartments | Characteristics |
|------|------|------|
| Cyclosporine A (CsA) | 3 (Gut/Central/Peripheral) | F=30%, CYP3A4, C₀ target 100-300 ng/mL |
| Tacrolimus (TAC) | 3 (Gut/Central/Peripheral) | F=25%, CYP3A5, C₀ target 5-15 ng/mL |
| Prednisone (PRED) | 2 (Gut/Central) | F=99%, GRα binding |
| Ruxolitinib (RUX) | 3 (Gut/Central/Peripheral) | F=95%, T½~3h, V=72L |
| Belumosudil (BELU) | 2 (Gut/Central) | F=80%, T½~20h, selective ROCK2 inhibition |
| MMF/MPA | 2 (Gut/Central) | F=94%, IMPDH inhibition |

**Immunology/biology (PD) — 16 compartments:**
- T cells: Th1, Th17, Treg, CD8 effector
- B cell pool
- 6 cytokines: TNF-α, IFN-γ, IL-17A, IL-10, TGF-β, IL-6
- 5 organ-damage compartments: Skin, Gut, Liver, Lung, Fibrosis

### Drug Effect Models (PD Effects)

| Drug | Target | Model |
|------|------|------|
| CsA/TAC | Calcineurin → NFAT → IL-2 | Emax (Hill n=1.5) |
| Ruxolitinib | JAK1/2 → STAT3/5 | Emax, including Treg-expansion effect |
| Belumosudil | ROCK2 → IRF4/STAT3 → Th17↓/Treg↑ | Emax, fibrosis suppression |
| MMF/MPA | IMPDH → suppression of lymphocyte proliferation | Emax |
| Prednisone | GRα → NF-κB suppression → broad cytokine↓ | Emax |

### Treatment Scenarios (6 Treatment Scenarios)

1. **No Prophylaxis** — baseline state (historical control)
2. **CsA Monoprophylaxis** — CsA monotherapy prophylaxis
3. **CsA + MMF** — standard combination prophylaxis
4. **TAC + MMF** — currently the most common standard regimen (NMDP/EBMT recommended)
5. **CsA → Ruxolitinib** — steroid-refractory cGvHD (based on REACH3)
6. **CsA → Belumosudil** — second-line-or-later cGvHD (based on ROCKstar)

---

## Shiny Dashboard (Interactive Dashboard)

**File**: `gvhd_shiny_app.R`

### Tab Layout (8 Tabs)

| Tab | Content |
|----|------|
| 1. **Patient & HSCT Profile** | Patient profile, GvHD risk radar chart, drug target descriptions |
| 2. **Drug PK Dashboard** | CsA/TAC/RUX/BELU concentration-time curves, PK summary table, PD effects |
| 3. **Immune Cell Dynamics** | Th1/Th17/Treg/CD8 dynamics, Th17/Treg ratio, B cells |
| 4. **Cytokine Network** | Pro-inflammatory (TNF-α, IFN-γ, IL-17A, IL-6) vs anti-inflammatory (IL-10, TGF-β) |
| 5. **Organ Damage & Endpoints** | Skin/gut/liver/lung damage scores, aGvHD/cGvHD grading, FFS |
| 6. **Scenario Comparison** | Side-by-side comparison of 6 treatment scenarios |
| 7. **Biomarkers** | ST2, REG3α, sTNFR1 biomarker time series |
| 8. **Mechanistic Map** | Full mechanistic map, PNG/SVG |

---

## Clinical Context

### Epidemiology
- **aGvHD incidence** after allogeneic HSCT: 30-50% with sibling donors, 50-70% with unrelated donors
- **cGvHD**: occurs in 40-70% of aGvHD survivors
- GvHD is the leading cause of non-relapse mortality (NRM)

### Three-Phase Pathophysiology
1. **Phase 1 (Afferent)**: tissue injury from conditioning → DAMP/PAMP release → host DC activation
2. **Phase 2 (Efferent)**: donor T-cell recognition of alloantigen → Th1/Th17 polarisation, CD8 activation
3. **Phase 3 (Effector)**: TNF-α/IFN-γ/IL-17A-mediated target-organ damage

### FDA-Approved Drugs
| Drug | Indication | Supporting trial |
|------|--------|------------|
| **Ruxolitinib** (Jakafi) | Steroid-refractory aGvHD, cGvHD | REACH1, REACH2, REACH3 |
| **Belumosudil** (Rezurock) | Second-line-or-later cGvHD | ROCKstar |
| **Ibrutinib** (Imbruvica) | cGvHD after first-line therapy | Single-arm study |

---

## Key Biomarkers

| Biomarker | Normal value | GvHD prediction | Target organ |
|-----------|--------|-----------|---------|
| ST2 (sST2) | <33 ng/mL | >33 ng/mL → grade 3-4 GI GvHD | Gut |
| REG3α | <10 ng/mL | >23 ng/mL → poor prognosis | Intestinal epithelium |
| sTNFR1 | <2 ng/mL | Elevation → aGvHD severity | Systemic/skin |
| CXCL9 | <100 pg/mL | IFN-γ-induced; T-cell homing | Systemic |
| Elafin | <10 ng/mL | Marker of skin GvHD | Skin |

---

## How to Run

```bash
# 1. Graphviz rendering
dot -Tsvg gvhd_qsp_model.dot -o gvhd_qsp_model.svg
dot -Tpng -Gdpi=150 gvhd_qsp_model.dot -o gvhd_qsp_model.png
```

```r
# 2. Run the mrgsolve model (R)
install.packages(c("mrgsolve", "dplyr", "ggplot2", "tidyr"))
source("gvhd_mrgsolve_model.R")

# 3. Shiny dashboard
install.packages(c("shiny", "shinydashboard", "plotly", "DT"))
shiny::runApp("gvhd_shiny_app.R")
```

---

## File List

| File | Description |
|------|------|
| `gvhd_qsp_model.dot` | Graphviz mechanistic map source (130+ nodes, 14 clusters) |
| `gvhd_qsp_model.svg` | Vector-format map (interactive zoom in browser) |
| `gvhd_qsp_model.png` | Raster-format map (150 dpi) |
| `gvhd_mrgsolve_model.R` | mrgsolve ODE model + 6 treatment scenarios + visualisation |
| `gvhd_shiny_app.R` | Shiny interactive simulator (8 tabs) |
| `gvhd_references.md` | 60 references (classified by section) |
| `README.md` | This file |

---

## Key References

- Zeiser R et al. *NEJM* 2020;382:1800 — Ruxolitinib for SR aGvHD (REACH2)
- Zeiser R et al. *NEJM* 2021;385:228 — Ruxolitinib for SR cGvHD (REACH3)
- Cutler C et al. *Blood* 2021;138:2278 — Belumosudil ROCKstar trial
- Ferrara JL et al. *Lancet* 2009;373:1550 — GvHD pathophysiology review
- Vander Lugt MT et al. *NEJM* 2013;369:529 — ST2 biomarker
- Peled JU et al. *NEJM* 2020;382:822 — Microbiome & GvHD outcome

Full set of 60 references → [`gvhd_references.md`](../../../graft-versus-host-disease/gvhd_references.md)

---

*Generated by Claude Code Routine (CCR) — 2026-06-25*  
*Disease category: Transplant Immunology / HSCT Complication*
