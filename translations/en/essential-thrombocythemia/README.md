# Essential Thrombocythemia (ET) — QSP Model

> **Quantitative Systems Pharmacology model** for Essential Thrombocythemia
> Driver mutations: JAK2 V617F / CALR / MPL | Megakaryopoiesis | Thrombosis | Drug PK/PD

[![Model](../../../essential-thrombocythemia/et_qsp_model.png)](../../../essential-thrombocythemia/et_qsp_model.svg)

---

## Disease Overview

**Essential Thrombocythemia (ET)** is a Philadelphia chromosome-negative myeloproliferative neoplasm (MPN) characterized by sustained thrombocytosis (platelets ≥450 ×10⁹/L) due to clonal expansion of a mutated hematopoietic stem cell. ET is driven predominantly by somatic mutations in:

| Mutation | Frequency | Key Features |
|----------|-----------|--------------|
| **JAK2 V617F** | ~55–65% | Constitutive JAK2 activation; higher thrombosis risk |
| **CALR Type 1** | ~15–20% | Del52bp; activates MPL; more benign course |
| **CALR Type 2** | ~8–12%  | Ins5bp; activates MPL; similar to JAK2 |
| **MPL W515L/K** | ~3–5%   | Activates TPOR/JAK2; clinically similar to CALR |
| **Triple negative** | ~10%  | Unknown driver; lowest thrombosis risk |

Key clinical features:
- **Thrombocytosis** (primary criterion; median PLT ~700–900 ×10⁹/L)
- **Thrombotic events** (arterial: stroke/MI; venous: DVT/PE; microvascular: erythromelalgia)
- **Hemorrhage** at very high platelet counts (>1500 ×10⁹/L, acquired vWD)
- **Splenomegaly** (~20–40% of patients)
- **Transformation** to post-ET myelofibrosis (~10–15% / 10 years) or AML (~1–3%)

---

## Mechanistic Map

**File:** [`et_qsp_model.dot`](../../../essential-thrombocythemia/et_qsp_model.dot) · [`et_qsp_model.svg`](../../../essential-thrombocythemia/et_qsp_model.svg) · [`et_qsp_model.png`](../../../essential-thrombocythemia/et_qsp_model.png)

### Clusters (9 subgraphs, ~160 nodes)

| # | Cluster | Key Nodes |
|---|---------|-----------|
| ① | Driver Mutations | JAK2 V617F, CALR Type1/2, MPL W515L, TET2, ASXL1, TP53, IDH1/2 |
| ② | JAK-STAT Signaling | TPO·TPOR→JAK2→STAT5/STAT3→PI3K/AKT/mTOR/MAPK |
| ③ | Megakaryopoiesis | LT-HSC→MPP→CMP→MEP→MKP→BFU-MK→CFU-MK→MK→Proplatelet→PLT |
| ④ | Cytokine Milieu | IL-3, IL-6, IL-11, TGF-β1, IFN-α/γ, PDGF-AB, VEGF, SCF |
| ⑤ | Platelet Biology | GPIb-vWF, GPIIb/IIIa, COX-1/TXA2, P2Y12, ADP, aggregation |
| ⑥ | Thrombosis/Hemorrhage | Arterial/venous events, erythromelalgia, coagulation cascade, acquired vWD |
| ⑦ | Drug PK/PD | Hydroxyurea (RRM2), Anagrelide (PDE3A), Ruxolitinib (JAK1/2), Peg-IFN-α, Aspirin (COX-1) |
| ⑧ | Clinical Endpoints | CHR, PHR, CMR, IPSS-ET, ELN response, OS, thrombosis-free survival |
| ⑨ | BM Pathology | Reticulin/collagen fibrosis, post-ET MF, blast phase, AML transformation |

### Key Pathway Logic

```
JAK2 V617F (constitutive) ──► JAK2 activation ──► STAT5 phospho
                                                  ──► PI3K/AKT/mTOR
                                                  ──► BCL-XL / MCL1 (anti-apoptosis)

CALR Type1/2 ──► aberrant MPL activation ──► JAK2 recruitment ──► (same cascade)

LT-HSC ──► MKP ──► IMK ──► PMK ──► Mature MK (endomitosis 8×) ──► Proplatelet ──► PLT
                                                     ↑
                               JAK2/STAT5 drives proliferation & survival
                                                     ↓
                               TGF-β1 secretion ──► fibroblast activation ──► reticulin fibrosis

Platelet > 1500 ──► acquired vWD ──► ultra-large vWF multimers ──► hemorrhage
PLT activation ──► TXA2/ADP ──► aggregation ──► arterial/venous thrombosis
```

---

## mrgsolve ODE Model

**File:** [`et_mrgsolve_model.R`](../../../essential-thrombocythemia/et_mrgsolve_model.R)

### Compartments (17)

| # | State | Description |
|---|-------|-------------|
| 1 | HSC  | Hematopoietic stem cells (JAK2-mutant pool) |
| 2 | MKP  | Megakaryocyte progenitors |
| 3 | MK   | Mature megakaryocytes |
| 4 | PLT  | Circulating platelets (×10⁹/L) |
| 5 | TPO  | Serum thrombopoietin (pg/mL) |
| 6 | JAK2 | JAK2 V617F allele burden (fraction 0–1) |
| 7 | SPL  | Spleen size (cm below costal margin) |
| 8–9 | HU_C / HU_P | Hydroxyurea central/peripheral (µg/mL) |
| 10–11 | ANA_C / ANA_P | Anagrelide central/peripheral |
| 12–13 | RUX_C / RUX_P | Ruxolitinib central/peripheral |
| 14–15 | pIFN_C / pIFN_P | Peg-IFN-α2a central/peripheral |
| 16 | RISK_T | Cumulative thrombosis risk (AU) |
| 17 | RISK_MF | Cumulative MF transformation risk (AU) |

### Key ODEs

```r
dPLT/dt = k_PLT_prod × k_MK_mat × MK - k_PLT_destr × (1 + SPL/10) × PLT
dJAK2/dt = k_JAK2_exp × JAK2 × (1-JAK2) - E_pIFN × 0.015 × JAK2
dRISK_T/dt = λ_T × (PLT/PLT_norm)^α_T × (1 + δ_JAK2 × JAK2) × (1 - E_ASA × 0.4)
```

### Drug Effect Functions (Hill equation)

```r
E_drug = Emax × C^γ / (EC50^γ + C^γ)
```

| Drug | EC50 | Emax | Target |
|------|------|------|--------|
| Hydroxyurea | 3.5 µg/mL | 0.85 | MKP proliferation (RRM2) |
| Anagrelide | 25 ng/mL | 0.75 | MK maturation (PDE3A) |
| Ruxolitinib | 150 ng/mL | 0.80 | JAK2 → pSTAT5 pathway |
| Peg-IFN-α | 8 ng/mL | 0.70 | JAK2+ clone suppression |
| Aspirin | — | 0.80 | TXA2/COX-1 (binary) |

### Treatment Scenarios (7)

| Scenario | Drug | Expected CHR | JAK2 Reduction |
|----------|------|-------------|----------------|
| ① No treatment | — | 0% | 0% (slow natural expansion) |
| ② Aspirin only | ASA 81 mg/d | 0% | 0% |
| ③ Hydroxyurea 500 mg | HU + ASA | ~30% | <5% |
| ④ Hydroxyurea 1500 mg | HU + ASA | ~60% | <10% |
| ⑤ Anagrelide 2 mg | ANA + ASA | ~35–45% | <5% |
| ⑥ Ruxolitinib 20 mg | RUX + ASA | ~60% | ~15–20% |
| ⑦ Peg-IFN-α2a 90 µg/wk | pIFN + ASA | ~45–75% | ~30–50% (CMR possible) |

### Clinical Trial Calibration

| Trial | Intervention | CHR | Reference |
|-------|-------------|-----|-----------|
| PT-1 (Harrison 2005 NEJM) | HU 1500 mg → ~60% CHR | 59% at 1yr | NEJM 353:33 |
| PT-1 | ANA 2 mg → ~35% CHR | 36% at 1yr | NEJM 353:33 |
| ANAHYDRET (Gisslinger 2013) | HU vs ANA — non-inferior | ~65% vs ~62% | Blood 121:1720 |
| RESPONSE-2 (ruxolitinib) | RUX 10 mg BID | ~50–60% PLT norm | Leukemia 2020 |
| Kiladjian 2013 | Peg-IFN-α2a | 45% CHR; JAK2 MR 18% | Haematologica |

---

## Shiny Dashboard

**File:** [`et_shiny_app.R`](../../../essential-thrombocythemia/et_shiny_app.R)

### 8 Interactive Tabs

| Tab | Content |
|-----|---------|
| ① Patient Profile | Baseline PLT, JAK2 AB%, age, prior thrombosis, ASXL1; IPSS-ET score; WHO criteria |
| ② Drug PK | Plasma concentration profiles; dose sliders for all 4 cytoreductive agents + aspirin |
| ③ Platelet Dynamics | Time-course PLT with CHR/PHR thresholds; response summary table |
| ④ JAK2 Allele Burden | JAK2 AB over time; CMR/PMR thresholds; molecular response classification |
| ⑤ Thrombosis Risk | Annual hazard over time; cumulative risk area plot; risk factor table |
| ⑥ Scenario Compare | Side-by-side multi-scenario plots; checkboxes for any subset of 7 scenarios |
| ⑦ Biomarker Panel | 6-panel dashboard: PLT, JAK2 AB, Spleen, TPO, MK pool, Cum MF Risk |
| ⑧ BM & Progression | MF transformation risk; pathway summary; co-mutation impact |

### How to Run

```r
install.packages(c("shiny","shinydashboard","dplyr","ggplot2","plotly","DT","purrr"))
shiny::runApp("essential-thrombocythemia/et_shiny_app.R")
```

---

## References

**File:** [`et_references.md`](../../../essential-thrombocythemia/et_references.md)

45 PubMed references organized in 12 sections:
1. Disease Definition & Epidemiology
2. Molecular Pathogenesis (JAK2/CALR/MPL)
3. JAK-STAT Signaling
4. Megakaryopoiesis
5. Thrombotic & Hemorrhagic Complications
6. Hydroxyurea PK/PD & Clinical Data
7. Anagrelide PK/PD & Clinical Data
8. Ruxolitinib in ET
9. Interferon Therapy
10. Disease Progression & MF Transformation
11. QSP / Mathematical Modeling
12. Clinical Guidelines & Risk Stratification

---

## File Summary

| File | Description | Size |
|------|-------------|------|
| [`et_qsp_model.dot`](../../../essential-thrombocythemia/et_qsp_model.dot) | Graphviz mechanistic map source (~160 nodes, 9 clusters) | ~18 KB |
| [`et_qsp_model.svg`](../../../essential-thrombocythemia/et_qsp_model.svg) | Vector map (scalable, full detail) | ~120 KB |
| [`et_qsp_model.png`](../../../essential-thrombocythemia/et_qsp_model.png) | Raster map (150 dpi) | ~180 KB |
| [`et_mrgsolve_model.R`](../../../essential-thrombocythemia/et_mrgsolve_model.R) | mrgsolve ODE QSP model (17 compartments, 7 scenarios) | ~9 KB |
| [`et_shiny_app.R`](../../../essential-thrombocythemia/et_shiny_app.R) | Shiny interactive dashboard (8 tabs) | ~14 KB |
| [`et_references.md`](../../../essential-thrombocythemia/et_references.md) | 45 PubMed references (12 sections) | ~6 KB |

---

## Key Model Assumptions & Limitations

1. **TPO feedback**: modeled as platelet-mediated TPO absorption (Kaushansky model); assumes normal hepatic TPO production
2. **JAK2 clone dynamics**: logistic growth model; CALR/MPL variants mapped to equivalent JAK2 burden for modeling purposes
3. **Drug PK**: steady-state CSS approximation used in Shiny (rapid equilibrium); mrgsolve has 2-compartment ODE
4. **Thrombosis risk**: hazard model calibrated to IPSS-ET data; does not model individual clotting factor dynamics
5. **MF progression**: simplified risk accumulator; does not include BM biopsy-grade ODE (future extension)
6. **Aspirin effect**: modeled as binary 80% TXA2 suppression; does not distinguish COX-1 isoforms
