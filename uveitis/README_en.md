# Uveitis QSP Model

[![Disease](https://img.shields.io/badge/Disease-Uveitis-blue)]()
[![Category](https://img.shields.io/badge/Category-Autoimmune%20%7C%20Ocular-orange)]()
[![Compartments](https://img.shields.io/badge/ODE%20Compartments-20-green)]()
[![Scenarios](https://img.shields.io/badge/Treatment%20Scenarios-7-purple)]()
[![References](https://img.shields.io/badge/References-60-yellow)]()

## Overview

**Uveitis** is inflammation of the uveal tract (iris, ciliary body,
choroid), classified anatomically as anterior, intermediate, posterior, or
panuveitis. Non-infectious uveitis is driven by autoimmune mechanisms
involving Th1/Th17 cells, TNF-α, IL-6, and VEGF-mediated disruption of the
blood-ocular barrier (BOB), and if left untreated it progresses to
permanent vision loss through cystoid macular oedema (CME), secondary
glaucoma, cataract, and retinal detachment. The primary site-specific PK
challenge is the blood-aqueous barrier (BAB) and the blood-retinal barrier
(BRB).

---

## Mechanistic Map

[![Mechanistic Map](uvt_qsp_model.png)](uvt_qsp_model.svg)

*Click to open the high-resolution SVG file.*

### Key Clusters (12 Clusters)

| Cluster | Content |
|---------|------|
| Ocular Anatomy & Immune Privilege | Anterior/posterior ocular structures, mechanisms of loss of immune privilege |
| Triggering Mechanisms | HLA-B27, IRBP/S-Ag autoantigens, infectious triggers, VKH/Birdshot |
| Innate Immune Activation | DCs, macrophages, neutrophils, NLRP3 inflammasome, NF-κB |
| T Cell Adaptive Immunity | Th1/Th17, Treg, CD8+, TCR/CD28 co-stimulation, FoxP3 |
| B Cell & Humoral Immunity | Plasma cells, autoantibodies (anti-IRBP, anti-recoverin), immune complexes |
| Cytokine & Chemokine Network | TNF-α, IL-6/STAT3, IL-17A, VEGF-A, COX-2/PGE2 |
| Blood-Ocular Barrier Disruption | BAB/BRB breakdown, anterior chamber cells/flare, CME, neovascularisation |
| Clinical Manifestations | Visual acuity (logMAR), IOP, SUN grading, cataract, secondary glaucoma |
| Drug PK | Eye drops, periocular/intravitreal injection, systemic oral/IV/SC |
| Drug PD Mechanisms | GR occupancy, NF-κB inhibition, TNF neutralisation, VEGF blockade, calcineurin inhibition |
| Complications & Disease Course | Acute-to-chronic transition, relapse rate, structural damage |
| Biomarkers & Monitoring | SUN grade, OCT-CST, FFA, serum TNF/IL-6, drug concentration/ADA |

---

## Model Structure

### Compartments (20)

| Module | Compartments | Description |
|------|------|------|
| **Drug PK** | Cgut, Cp, Cperiph | Oral/systemic PK (1-compartment absorption + 2-compartment) |
| | C_ant, C_post | Anterior chamber and posterior segment drug concentration |
| | C_depot | Intravitreal sustained-release implant (Ozurdex) |
| **Immune** | T_eff, T_reg | Th1/Th17 effector T cells, Treg |
| | APC_act, Macro | Activated APC/DC, M1 macrophages |
| **Cytokines** | TNF, IL6, IL17, VEGF | The four core cytokines |
| **Barrier** | BAB_int, BRB_int | Blood-aqueous barrier, blood-retinal barrier integrity (0-1) |
| **Clinical PD** | Cells_AH, CME | Anterior chamber cells, cystoid macular oedema |
| | VA_def, IOP_e | Visual acuity deficit (logMAR), excess IOP |
| | GR_occ | Glucocorticoid receptor occupancy |

### Treatment Scenarios (7)

| ID | Treatment | Dose/regimen | Key mechanism | Supporting trial |
|----|--------|----------|----------|--------------|
| S1 | No treatment | — | Natural history | — |
| S2 | Prednisolone 1% eye drops | Every 6 hours | Local NF-κB inhibition | Cunningham 1990 |
| S3 | Periocular triamcinolone 40 mg | Single dose | Sustained corticosteroid | Sallam 2011 |
| S4 | Intravitreal dexamethasone implant (Ozurdex) | 0.7 mg sustained release | Restores posterior BRB, ↓CME | HURON trial (Lowder 2011) |
| S5 | Systemic prednisone 1 mg/kg/day (tapered) | Daily → maintenance | Systemic GR activation, T-cell suppression | MUST trial (Kempen 2011) |
| S6 | Adalimumab 40 mg SC q2w | Fortnightly subcutaneous injection | TNF-α neutralisation, ↓leukostasis | VISUAL I (Jaffe 2016 NEJM) |
| S7 | Combination (prednisone + adalimumab) | S5 + S6 | Dual mechanism, rapid barrier recovery | VISUAL I/II |

---

## Key Equations

### BAB integrity
```
dBAB_int/dt = k_BABrep x (1 + 2*GR_occ) x (1 - BAB_int)
            - k_BABdeg x (TNF + 0.5*IL6) x BAB_int
```

### BRB integrity
```
dBRB_int/dt = k_BRBrep x (1 + 1.5*GR_occ + 2*aVEGF_eff) x (1 - BRB_int)
            - k_BRBdeg x (VEGF + 0.3*TNF) x BRB_int
```

### Cystoid macular oedema (CME)
```
dCME/dt = k_CMEform x (1-BRB_int) x VEGF
        - k_CMEres x (1 + 2*aVEGF_eff + GR_occ) x CME
```

### GR occupancy (Hill function)
```
GR_occ = Emax_cs x [Drug]^H / (EC50_cs^H + [Drug]^H)
```

---

## Clinical Parameter Calibration

| Parameter | Model value | Clinical data | Source |
|---------|--------|-----------|------|
| Adalimumab onset of effect | ~2-4 weeks | Significant reduction in inflammation within 4 weeks | VISUAL I |
| Ozurdex CME resolution rate | 60-70% | 65% complete resolution at day 60 | HURON trial |
| IOP elevation in steroid responders | 30-40% | ~38% with a ≥5 mmHg rise | Jaffe 2006 |
| Treg/Teff ratio | 0.67 at flare | About 1/3 of healthy controls | Caspi 2010 |
| BAB breakdown → anterior chamber cells | ~48-72h | Flare onset within 1-3 days | Nussenblatt 1985 |

---

## Shiny App Structure (8 Tabs)

| Tab | Contents |
|----|------|
| 1. Patient Profile | Patient settings, disease severity, treatment selection, run |
| 2. Pharmacokinetics | Plasma/ocular drug concentrations, GR occupancy |
| 3. Immune & Cytokines | T-cell populations, APC/macrophages, TNF/IL-6/IL-17/VEGF |
| 4. Barrier Integrity | BAB/BRB integrity, anterior chamber cells (SUN grade) |
| 5. Clinical Endpoints | Visual acuity (logMAR), OCT-CST, IOP, clinical milestones |
| 6. Scenario Comparison | Comparison of the 7 treatments (visual acuity/CME/TNF/summary table) |
| 7. Virtual Patients | VP population simulation (n up to 500), response rate by subtype |
| 8. Biomarker Monitor | Tracking of cytokines/OCT/drug concentration, monitoring schedule |

---

## Files

| File | Description |
|------|------|
| [`uvt_qsp_model.dot`](uvt_qsp_model.dot) | Graphviz mechanistic map (130+ nodes, 12 clusters) |
| [`uvt_qsp_model.svg`](uvt_qsp_model.svg) | Vector-format map (high resolution) |
| [`uvt_qsp_model.png`](uvt_qsp_model.png) | Raster-format map (150 dpi) |
| [`uvt_mrgsolve_model.R`](uvt_mrgsolve_model.R) | mrgsolve ODE QSP model (20 compartments, 7 scenarios, VP n=200) |
| [`uvt_shiny_app.R`](uvt_shiny_app.R) | Shiny interactive dashboard (8 tabs) |
| [`uvt_references.md`](uvt_references.md) | 60 references (15 sections) |
| [`README.md`](README.md) | This file |

---

## Usage

```r
# 1. Run the mrgsolve model
install.packages(c("mrgsolve","dplyr","ggplot2","tidyr"))
source("uvt_mrgsolve_model.R")

# 2. Run the Shiny app
install.packages(c("shiny","shinydashboard","DT"))
shiny::runApp("uvt_shiny_app.R")
```

```bash
# Re-render the Graphviz map
dot -Tsvg uvt_qsp_model.dot -o uvt_qsp_model.svg
dot -Tpng -Gdpi=150 uvt_qsp_model.dot -o uvt_qsp_model.png
```

---

## Key Findings

1. **Adalimumab (S6)** is the most effective at restoring the BAB and
   suppressing anterior-chamber inflammation through TNF neutralisation,
   reproducing the results of the VISUAL I/II trials.
2. **The intravitreal implant (S4)** is excellent at resolving CME but
   carries a risk of long-term steroid-induced IOP elevation.
3. **Combination therapy (S7)** shows the fastest barrier recovery and is
   superior at preventing structural damage.
4. **Virtual-patient analysis** shows that the Th17-dominant subtype
   responds relatively slowly to anti-TNF monotherapy alone, suggesting
   that adding an IL-17 blocker should be considered.
5. **Steroid eye drops (S2)** are effective for anterior uveitis but reach
   insufficient intraocular concentration for posterior lesions.

---

*Generated: 2026-06-26 | QSP Disease Model Library | Claude Code Routine (CCR)*
