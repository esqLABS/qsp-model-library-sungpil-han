# Multiple Myeloma (MM) QSP Model

[![MM QSP Mechanistic Map](../../../multiple-myeloma/mm_qsp_model.png)](../../../multiple-myeloma/mm_qsp_model.svg)

## Overview

Multiple myeloma (MM) is a haematological malignancy driven by proliferation
of a malignant plasma-cell clone in the bone marrow, and is the second most
common blood cancer, with roughly 35,000 new diagnoses per year in the
United States. It is characterised by osteolytic bone lesions, excess
M-protein secretion, renal failure, anaemia, and hypercalcaemia (the CRAB
criteria).

This QSP model includes a 13-cluster mechanistic map and mrgsolve ODE
simulation integrating the **bone marrow microenvironment (BMME), tumour
cell dynamics, bone remodelling, and drug PK/PD**.

---

## Key Pathophysiology

| Pathway | Core mechanism | Clinical outcome |
|------|-------------|---------|
| Proteasome (UPS) overload | Massive immunoglobulin synthesis by plasma cells -> ER stress -> UPS hyperactivation | Drug target (BTZ/CFZ) |
| IL-6 / JAK-STAT3 | BMSC secretes IL-6 -> MM cell proliferation and anti-apoptosis | M-protein up, Hgb down |
| RANKL/OPG imbalance | MM cells induce RANKL up / OPG down -> OC hyperactivation | Osteolytic lesions, pathological fractures |
| DKK1 / Wnt inhibition | MM cells secrete DKK1/sclerostin -> OB function suppressed | Suppressed bone formation -> net bone loss |
| IKZF1/IKZF3 (Ikaros/Aiolos) | IMiD reroutes the cereblon (CRBN) E3 ligase -> IKZF1/3 degradation -> IRF4 down | LEN/POM efficacy mechanism |
| CD38 / immunity | Plasma-cell-specific CD38 overexpression -> DARA ADCC/CDC activity | Daratumumab target |
| BCL-2 / BH3 mimetics | t(11;14) patients have increased BCL-2 dependence -> VEN releases BIM | Venetoclax sensitivity |
| PD-1/PD-L1 axis | MM cell PD-L1 expression -> CTL exhaustion | Immune evasion |

---

## File List

| File | Description |
|------|------|
| [`mm_qsp_model.dot`](../../../multiple-myeloma/mm_qsp_model.dot) | Graphviz mechanistic map (13 clusters, 130+ nodes) |
| [`mm_qsp_model.svg`](../../../multiple-myeloma/mm_qsp_model.svg) | SVG vector image (click to zoom) |
| [`mm_qsp_model.png`](../../../multiple-myeloma/mm_qsp_model.png) | PNG raster image (150 dpi) |
| [`mm_mrgsolve_model.R`](../../../multiple-myeloma/mm_mrgsolve_model.R) | mrgsolve ODE model (30 compartments, 6 treatment scenarios) |
| [`mm_shiny_app.R`](../../../multiple-myeloma/mm_shiny_app.R) | Shiny dashboard (6 interactive tabs) |
| [`mm_references.md`](../../../multiple-myeloma/mm_references.md) | 50 references (with PubMed links) |

---

## Mechanistic Map Clusters

1. **Plasma Cell Differentiation** — HSC -> CLP -> B cell -> GC -> plasma cell
2. **Myeloma Cell Biology** — MGUS -> SMM -> MM cell; tumour genomic abnormalities
3. **Bone Marrow Microenvironment (BMME)** — BMSC, BMEC, osteoblast/osteoclast, MSC, adhesion molecules
4. **Cytokine & Growth Factor Signaling** — IL-6, APRIL/BAFF, RANKL/OPG, DKK1, MIP-1α
5. **Intracellular Signaling (MM Cell)** — JAK-STAT3, PI3K/AKT/mTOR, NF-κB, proteasome, BCL-2 family
6. **Immune Evasion** — CTL, NK cell, Treg, MDSC, PD-1/PD-L1, CD38, BCMA, SLAMF7
7. **Myeloma Bone Disease** — Osteolytic lesion, NTX/PINP biomarkers, hypercalcaemia, SRE
8. **M-Protein & Clinical Markers** — M-protein, sFLC, β₂M, LDH, albumin, BMPC%
9. **Drug PK** — Bortezomib (3-cmt), carfilzomib (2-cmt), lenalidomide, daratumumab (TMDD), dexamethasone, venetoclax, zoledronic acid
10. **Drug PD** — Proteasome inhibition, CRBN/IKZF1/3 degradation, ADCC/CDC, BCL-2 BH3 mimetics, CAR-T/BiTE
11. **Clinical Endpoints** — ISS/R-ISS staging, CR/VGPR/PR/SD/PD, MRD negativity, PFS, OS
12. **Renal Complications** — Cast nephropathy, AL amyloidosis, eGFR
13. **Hematologic Complications** — Anaemia, neutropenia, thrombocytopenia

---

## mrgsolve Model Compartments

| Category | Compartments | Count |
|---------|------|---|
| Tumour cells | MM_S (sensitive), MM_R (resistant) | 2 |
| Disease markers | MP (M-protein), FLC (sFLC), IL6, VEGF | 4 |
| Bone remodelling | OB, OC, BV, NTX, PINP | 5 |
| Systemic biomarkers | Hgb, B2M | 2 |
| Bortezomib PK | BTZ1, BTZ2, BTZ3 | 3 |
| Lenalidomide PK | LEN_gut, LEN1 | 2 |
| Daratumumab PK | DARA1, DARA2, DARA_CD38 | 3 |
| Dexamethasone PK | DEX_gut, DEX1 | 2 |
| Venetoclax PK | VEN_gut, VEN1 | 2 |
| Zoledronic Acid PK | ZOL1, ZOL2, ZOL_bone | 3 |
| **Total** | | **28** |

---

## Treatment Scenarios

| Scenario | Regimen | Clinical trial basis |
|---------|--------|------------|
| 1 | Untreated (natural history) | — |
| 2 | **VRd** (bortezomib + lenalidomide + dex) | SWOG S0777 (Durie 2017) |
| 3 | **DRd** (daratumumab + lenalidomide + dex) | MAIA (Facon 2019) |
| 4 | **KRd** (carfilzomib + lenalidomide + dex) | ASPIRE (Stewart 2015) |
| 5 | **VenDex** (venetoclax + dex, t(11;14)) | BELLINI (Kumar 2020) |
| 6 | **DVRd** (dara + VRd, high-risk group) | PERSEUS (Sonneveld 2024) |

---

## Drug PK Parameter Summary

| Drug | Route | CL (L/hr) | Vd (L) | Half-life | Notes |
|------|---------|-----------|--------|--------|--------|
| Bortezomib (BTZ) | IV/SC | 9.0 | 4.7 (central) | ~76h (deep) | 3-compartment, SC F=83% |
| Carfilzomib (CFZ) | IV | 245 | 8.0 | ~1h | Very rapid elimination, irreversible |
| Lenalidomide (LEN) | PO | 3.2 | 65 | ~3h | F=90%, renal excretion |
| Daratumumab (DARA) | IV/SC | 0.2 (linear) | 3.0 (70kg) | ~14d | TMDD model |
| Dexamethasone (DEX) | PO | 20 | 130 | ~4h | F=78% |
| Venetoclax (VEN) | PO | 14 | 256 | ~12h | F=55% (with high-fat meal) |
| Zoledronic acid (ZOL) | IV | 3.8 | 4.0 | ~105h (bone) | Slow release after bone binding |

---

## Response Criteria (IMWG 2016)

| Response | M-protein criterion | Other criteria |
|------|-------------|--------|
| **sCR** | Negative serum/urine M-protein | Normal FLC ratio + BM <5% PC + NGS MRD negative |
| **CR** | Negative serum/urine M-protein | BM <5% PC |
| **VGPR** | >=90% reduction in M-protein | Urine <100 mg/24hr |
| **PR** | >=50% reduction in M-protein | Urine >=90% reduction (or <200 mg/24hr) |
| **SD** | Does not meet PR criteria, has not reached PD criteria | — |
| **PD** | >=25% increase in M-protein (from nadir) | New lesion or 25% increase in an existing lesion |

---

## Shiny App Tabs

| Tab | Content |
|----|------|
| **Patient profile** | Input for weight, age, ISS stage, cytogenetic risk, baseline biomarkers |
| **PK profile** | Drug concentration-time curves by regimen, EC50 reference lines, visualisation of drug effect |
| **PD / tumour burden** | MM cell dynamics (sensitive/resistant clones), M-protein response |
| **Clinical endpoints** | Response classification, bone remodelling, haemoglobin trend, response table |
| **Scenario comparison** | Simultaneous comparison of multiple regimens, best-response waterfall plot |
| **Biomarkers** | Bone turnover markers, cytokines, β₂M vs. M-protein correlation |

---

## Parameter Calibration Based on Key Clinical Trial Data

| Parameter | Value | Basis |
|---------|------|------|
| MM cell doubling time | ~35 days (kgrow=0.02/day) | Turesson 2004 |
| M-protein half-life | ~9 days (kMP_elim=0.08) | Clinical observation |
| VRd ORR (>=PR) | ~82% | SWOG S0777 |
| DRd PFS improvement | Median 61.9 months vs. 34.4 months | MAIA trial |
| DVRd MRD-negative rate | ~60% (10^-5) | PERSEUS trial |
| ZOL suppression of bone resorption | ~60-70% | Clinical NTX reduction |

---

*Date created: 2026-06-21 | Auto-generated by Claude Code Routine (CCR)*
