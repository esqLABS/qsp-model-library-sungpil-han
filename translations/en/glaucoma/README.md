# Primary Open-Angle Glaucoma (POAG) — QSP Model

**Abbreviation:** POAG · **Category:** Ophthalmology / Chronic Disease
**Date Created:** 2026-06-27 · **Directory:** `glaucoma/`

---

## Overview

Primary open-angle glaucoma (POAG) is the leading cause of irreversible blindness worldwide.
This QSP model integrates aqueous humor dynamics (Goldmann equation), trabecular meshwork (TM)
biology (ECM accumulation, ROCK signaling), optic nerve head (ONH) biomechanics (lamina cribrosa
deformation, BDNF retrograde transport), and retinal ganglion cell (RGC) apoptosis to quantitatively
predict the long-term effects of five drug classes and surgical interventions on IOP and visual field
(VF-MD) preservation.

---

## Mechanistic Map

[![POAG QSP Model](poag_qsp_model.png)](../../../glaucoma/poag_qsp_model.svg)

> Click the image to open the interactive SVG map.
> Rendered with `sfdp` layout engine (Graphviz 2.42).

### 10 Subgraph Clusters:
| # | Cluster | Key Components |
|---|---------|--------------|
| 1 | Genetics · Risk Factors | MYOC, OPTN, CDKN2B-AS1, CAV1/2, LOXL1 (XFS) |
| 2 | Ciliary Body · Aqueous Production | CA-II/IV, Na⁺/K⁺-ATPase, β₂-AR, α₂-AR, AQP1 |
| 3 | Trabecular Meshwork · Conventional Outflow | RhoA→ROCK→actin, TGF-β2→Smad2/3→ECM, TM senescence |
| 4 | Uveoscleral Outflow | FP-R→PKC→MMP-1,3→ciliary ECM remodeling |
| 5 | Intraocular Pressure Dynamics | Goldmann equation, diurnal variation, OHT, NTG |
| 6 | Drug PK/PD | PGA, BB, CAI, A2A, ROCK-I, SLT, MIGS, Trabeculectomy |
| 7 | Optic Nerve Head Biology | LC biomechanics, TPG, axonal transport, BDNF, TrkB/p75NTR |
| 8 | RGC Apoptosis Pathway | NMDA-Ca²⁺-nNOS, ROS, Bax/Bcl-2→Casp9→Casp3, DLK→JNK |
| 9 | Clinical Endpoints · Monitoring | RNFL, GCL, VF-MD, VFI, CDR, GPA |
| 10 | Systemic · Vascular Factors | OPP, autoregulation, XFS, PDG, neuroprotection |

---

## mrgsolve ODE Model

**File:** [`poag_mrgsolve_model.R`](../../../glaucoma/poag_mrgsolve_model.R)

### Compartments (16 ODE Compartments):
| # | State Variable | Description |
|---|---------|------|
| 1–5 | C_PGA, C_BB, C_CAI, C_A2A, C_ROCK | Drug PK in aqueous (ng/mL) |
| 6 | F_aq | Aqueous production rate (μL/min) |
| 7 | F_uv | Uveoscleral outflow rate (μL/min) |
| 8 | C_tm | TM outflow facility (μL/min/mmHg) |
| 9 | IOP | Intraocular pressure (mmHg) — Goldmann dynamic |
| 10 | ECM_TM | ECM accumulation index in TM |
| 11 | BDNF | BDNF in ONH (pg/mL) |
| 12 | Casp3 | Caspase-3 apoptotic index (0–1) |
| 13 | RGC | Retinal ganglion cell count (millions) |
| 14 | RNFL | RNFL thickness (μm) |
| 15 | VF_MD | Visual field mean deviation (dB) |

### Treatment Scenarios (8):
1. **Untreated POAG** — natural history (IOP ≈ 24 mmHg)
2. **Latanoprost QD** — PGA, ↑F_uv +100%
3. **Timolol BID** — β₂-blocker, ↓F_prod 30%
4. **Dorzolamide TID** — CAI, ↓F_prod 25%
5. **Brimonidine BID** — A2A + neuroprotection (↑BDNF)
6. **Netarsudil QD** — ROCK-I, ↑C_tm 35% + ↓P_ep 3 mmHg
7. **Latanoprost + Timolol** — Fixed-dose combination
8. **Triple therapy** — PGA + BB + CAI

### Key Equations:
```
Goldmann:  IOP_eq = (F_aq - F_uv) / C_tm + P_ep
TM ECM:    dECM/dt = k_prod + k_IOP·max(IOP−18,0) − k_clear·ECM − ROCK_effect
BDNF:      dBDNF/dt = k_prod·(1+A2A_boost) − (k_deg + k_IOP·max(IOP−18,0))·BDNF
RGC loss:  dRGC/dt = −(k_base + k_Casp3·Casp3_eff) · RGC
VF-MD:     VF_MD_eq = −25·(1 − RGC/RGC₀)^2.5  dB
```

---

## Shiny App

**File:** [`poag_shiny_app.R`](../../../glaucoma/poag_shiny_app.R)

### 7 Tabs:
| Tab | Content |
|----|------|
| 1. Patient Profile | Demographics, diagnosis stage, baseline IOP/RNFL/VF-MD |
| 2. Drug PK | Aqueous concentration profiles for 5 drug classes |
| 3. IOP Dynamics | Real-time IOP, F_aq, F_uv, C_tm with treatment |
| 4. Clinical Endpoints | VF-MD, RNFL, RGC count over 10 years |
| 5. Scenario Comparison | 8-scenario comparison with IOP, VF, RNFL tables |
| 6. Biomarkers & Neuroprotection | BDNF, Caspase-3, ECM_TM dynamics |
| 7. Sensitivity Analysis | Tornado plot, IOP vs VF scatter, parameter sensitivity |

---

## References

**File:** [`poag_references.md`](../../../glaucoma/poag_references.md) — **66** references

### Key Clinical Trials:
- **AGIS** (2000) — IOP control prevents VF deterioration
- **OHTS** (2002) — Ocular Hypertension Treatment Study
- **EMGT** (2003) — Early Manifest Glaucoma Trial
- **CIGTS** (2001–2011) — Collaborative Initial Glaucoma Treatment Study
- **LiGHT** (2019) — SLT first-line vs. eye drops (Lancet)
- **UKGTS** (2015) — Latanoprost RCT (Lancet)
- **TVT** (2012) — Tube vs. Trabeculectomy

---

## Key Parameter Sources

| Parameter | Value | Source |
|---------|---|------|
| F_prod baseline | 2.5 μL/min | Brubaker 1991 IOVS |
| F_uv baseline | 0.4 μL/min | Brubaker 1991 IOVS |
| C_tm normal | 0.30 μL/min/mmHg | Goldmann 1951 |
| P_ep | 8.0 mmHg | Goldmann 1951 |
| PGA ↑F_uv | +80–120% | Stjernschantz 2001 IOVS |
| BB ↓F_prod | 25–35% | Liu 2003 IOVS |
| ROCK-I ↑C_tm | +30–40% | Tian 2005 Br J Ophthalmol |
| RGC normal loss | 0.5%/yr | Harwerth 2004 IOVS |

---

## Usage

```r
# mrgsolve simulation
source("poag_mrgsolve_model.R")

# Run the Shiny app
library(shiny)
runApp("poag_shiny_app.R")
```

### Required Packages:
```r
install.packages(c("mrgsolve", "dplyr", "ggplot2", "tidyr",
                   "shiny", "shinydashboard", "plotly", "DT"))
```

---

## Files

| File | Description |
|------|------|
| `poag_qsp_model.dot` | Graphviz DOT — 150+ nodes, 10 subgraphs |
| `poag_qsp_model.svg` | Vector image (interactive) |
| `poag_qsp_model.png` | Raster image (150 dpi) |
| `poag_mrgsolve_model.R` | ODE model + 8-scenario simulation |
| `poag_shiny_app.R` | 7-tab Shiny dashboard |
| `poag_references.md` | List of 66 references |
| `README.md` | This file |
