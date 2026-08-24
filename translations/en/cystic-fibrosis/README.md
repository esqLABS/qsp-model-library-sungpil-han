# Cystic Fibrosis (CF) — QSP Model

## Overview

Cystic fibrosis is an autosomal recessive genetic disease caused by mutations in the
**CFTR (Cystic Fibrosis Transmembrane conductance Regulator)** gene, a rare disease
affecting about 100,000 people worldwide. Dysfunction of the CFTR protein causes thick
mucus secretion from exocrine glands and affects multiple organs, including the lungs,
digestive tract, and reproductive organs. The recent emergence of CFTR modulators
(Trikafta = ELX/TEZ/IVA) has fundamentally changed the treatment paradigm.

---

## Pathophysiology

| Stage | Key mechanism | Clinical outcome |
|------|----------|---------|
| CFTR gene mutation | Six classes (I–VI) of dysfunction | Absent/dysfunctional CFTR protein |
| ΔF508 processing defect (Class II) | ER quality control → ERAD (~99% degraded) | <1% of CFTR reaches the membrane |
| ASL dehydration | ENaC hyperactivity + reduced CFTR secretion → osmotic imbalance | Collapse of the airway surface liquid (<4 μm) |
| Impaired mucus transport | Reduced mucociliary transport velocity, increased viscosity | Bacterial (Pseudomonas) colonization |
| Chronic infection | Pa biofilm formation, antibiotic resistance | Neutrophil infiltration, airway destruction |
| Systemic inflammation | Excess IL-8/TNF-α/NE secretion | Decreased FEV1, pulmonary fibrosis |

---

## Model Components

### 1. Mechanistic map (`cf_qsp_model.dot`)

| Cluster | Key content | Node count |
|---------|---------|--------|
| ① CFTR Gene & mRNA | Gene, mRNA, six mutation classes | 12 |
| ② CFTR Protein Processing | ER QC, chaperones, ERAD, corrector binding | 14 |
| ③ CFTR Channel (Membrane) | Golgi trafficking, membrane localization, channel gating | 15 |
| ④ ENaC & Ion Transport | ENaC, NKCC1, Na/K-ATPase, TMEM16A | 11 |
| ⑤ ASL & Mucociliary Clearance | PCL, mucus layer, ciliary beat, MCC | 13 |
| ⑥ Airway Inflammation | NF-κB, IL-8, TNF-α, neutrophils, NETs, ROS | 25 |
| ⑦ Infection & Biofilm | P. aeruginosa, S. aureus, biofilm, antibiotics | 15 |
| ⑧ Lung Function | FEV1, LCI, bronchiectasis, transplant criteria | 13 |
| ⑨ CFTR Modulator PK/PD | Ivacaftor/ELX/TEZ pharmacokinetic compartments | 16 |
| ⑩ Other Therapies | DNase, hypertonic saline, PERT, gene therapy | 10 |
| ⑪ Systemic Effects | CFRD, pancreatic function, liver, bone density | 11 |
| ⑫ Clinical Endpoints | Sweat chloride, NPD, ppFEV1, CFQ-R | 10 |
| **Total** | | **165 nodes, 247 edges** |

### 2. mrgsolve ODE model (`cf_mrgsolve_model.R`)

**Compartments (25 ODEs):**
- Drug PK: Ivacaftor (3 compartments) + Elexacaftor (2 compartments) + Tezacaftor (2 compartments)
- CFTR biology: Band B processing, membrane-localized CFTR
- Airway surface liquid (ASL)
- Infection: Pa planktonic bacteria + biofilm
- Inflammation: IL-8, neutrophils, cumulative damage score
- Lung function: ppFEV1, cumulative exacerbation count
- Systemic: BMI, pancreatic function

**Treatment scenarios (7):**

| Scenario | Drug | Mutation | Target clinical metric |
|---------|-----|---------|-------------|
| 1 | No treatment | ΔF508/ΔF508 | Baseline |
| 2 | Ivacaftor (Kalydeco) | G551D | ppFEV1 +10.6 pp |
| 3 | LUM/IVA (Orkambi) | ΔF508/ΔF508 | ppFEV1 +2.6 pp |
| 4 | TEZ/IVA (Symdeko) | ΔF508/ΔF508 | ppFEV1 +3.4 pp |
| 5 | ETI (Trikafta) | ΔF508/ΔF508 | ppFEV1 +**14.3 pp** |
| 6 | ETI + Tobramycin | ΔF508/ΔF508 | Combined infection control |
| 7 | Early ETI (age 6) | ΔF508/ΔF508 | Effect of early intervention |

### 3. Shiny app (`cf_shiny_app.R`)

6 tabs:
1. **Patient Profile** — mutation class selection, clinical trial benchmarks
2. **CFTR Modulator PK** — Ivacaftor/ELX/TEZ plasma concentration profiles
3. **CFTR Function** — correction rate, potentiation rate, sweat chloride
4. **Lung Function** — ppFEV1 trajectory, cumulative exacerbation count
5. **Scenario Comparison** — comparison of 5 treatments, 52-week endpoint table
6. **ASL & Infection** — airway surface liquid height, bacterial burden, inflammatory markers

---

## Summary of Drug PK/PD Parameters

| Drug | Route | F (%) | Vc (L) | t½ (h) | Key mechanism |
|-----|------|--------|--------|--------|------------|
| **Ivacaftor** (VX-770) | p.o. 150mg q12h | 67 | 97 | 12 | Stabilizes the CFTR open state (↑Po) |
| **Elexacaftor** (VX-445) | p.o. 200mg q24h | 80 | 193 | 27 | NBD1 binding → CFTR correction |
| **Tezacaftor** (VX-661) | p.o. 100mg q24h | 70 | 271 | 14 | MSD1 binding → trafficking correction |
| **Lumacaftor** (VX-809) | p.o. 200mg q12h | 65 | ~300 | 26 | MSD1 binding (limited effect alone) |
| **Tobramycin** (TOBI) | Inhaled 300mg bid | — | — | — | Pa bactericidal (aminoglycoside) |

---

## Clinical Trial Calibration Targets

| Clinical trial | Treatment | Key metric | Result |
|---------|-----|---------|-----|
| STRIVE (NEJM 2011) | Ivacaftor | ppFEV1 | +10.6 pp (G551D) |
| TRAFFIC/TRANSPORT (NEJM 2015) | LUM/IVA | ppFEV1 | +2.6–3.0 pp (ΔF508) |
| EVOLENT (NEJM 2017) | TEZ/IVA | ppFEV1 | +3.4 pp (ΔF508) |
| VX-445-102 (NEJM 2019) | **ETI (Trikafta)** | ppFEV1 | **+14.3 pp** |
| VX-445-103 (AURORA) | **ETI** | Sweat Cl | **-41.8 mmol/L** |

---

## Model File List

| File | Description |
|------|------|
| [cf_qsp_model.dot](../../../cystic-fibrosis/cf_qsp_model.dot) | Graphviz mechanistic map (165 nodes, 12 clusters) |
| [cf_qsp_model.svg](../../../cystic-fibrosis/cf_qsp_model.svg) | SVG vector image |
| [cf_qsp_model.png](../../../cystic-fibrosis/cf_qsp_model.png) | PNG 150 dpi |
| [cf_mrgsolve_model.R](../../../cystic-fibrosis/cf_mrgsolve_model.R) | mrgsolve ODE model (25 compartments, 7 scenarios) |
| [cf_shiny_app.R](../../../cystic-fibrosis/cf_shiny_app.R) | Shiny dashboard (6 tabs) |
| [cf_references.md](../../../cystic-fibrosis/cf_references.md) | References (53 articles, PubMed links) |

---

## Preview

[![CF QSP Model](../../../cystic-fibrosis/cf_qsp_model.png)](../../../cystic-fibrosis/cf_qsp_model.svg)
