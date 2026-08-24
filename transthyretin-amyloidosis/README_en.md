# Transthyretin Amyloidosis (ATTR) — QSP Model

> **Directory:** `transthyretin-amyloidosis/` | **Abbreviation:** ATTR | **Date:** 2026-06-24
> **Category:** Rare Disease / Protein Misfolding Disorder / Amyloidosis

[![ATTR QSP mechanistic map](attr_qsp_model_en.png)](attr_qsp_model_en.svg)

---

## Disease Overview

**Transthyretin amyloidosis (ATTR)** is a protein-misfolding disease in which **transthyretin
(TTR)**, a transport protein synthesised in the liver, dissociates from its normal tetramer
structure and assembles into misfolded monomers → toxic oligomers → mature amyloid fibrils.

| Characteristic | ATTRwt (Wild-type) | ATTRv (Hereditary) |
|------|----------------|----------------|
| Cause | Age-related TTR destabilisation | Missense mutation in the *TTR* gene |
| Age of onset | ≥60 years, male predominance | 30–70 years, varies by mutation |
| Major phenotype | Cardiomyopathy (HFpEF→HFrEF) | Polyneuropathy / mixed cardiac |
| Representative mutations | — | V30M (Portugal/Japan/Sweden), V122I (of African descent) |
| Prevalence | ~10–13% of elderly HF autopsies | ~50,000 people worldwide |

---

## Core Pathogenic Mechanism — Four Stages

| Stage | Mechanism | Key Mediators |
|------|------|------------|
| **1. TTR synthesis** | TTR monomer synthesis in hepatocytes → β-sheet folding → homotetramer assembly → plasma secretion | TTR mRNA, pre-TTR → mature TTR (14 kDa), tetramer (55 kDa) |
| **2. Tetramer dissociation** | Thermal/pH stress or an ATTRv mutation drives tetramer dissociation → generation of misfolded monomers → nucleation → oligomers | Dissociation rate (kdis), oligomers (toxic intermediates), protofibrils |
| **3. Amyloid deposition** | Mature fibrils selectively deposit in the heart (ATTRwt≫), peripheral nerves (ATTRv≫), kidney, spleen, and GI tract | FIB_HRT, FIB_NRV, FIB_SYS |
| **4. Organ damage** | Deposition → cytotoxicity (NLRP3·IL-1β·TNF-α) → cardiomyocyte apoptosis → ventricular hypertrophy/diastolic dysfunction / Schwann cell compression → axonal degeneration → polyneuropathy | LVEF↓, NT-proBNP↑, NIS↑, mBMI↓ |

---

## Treatment Mechanism Overview

| Drug | Mechanism | Target | Key Clinical Trial |
|------|------|------|-------------|
| **Tafamidis** (61mg PO QD) | Occupies the T4-binding site → kinetically stabilises the tetramer (Emax ~80%) | Decreased TTR tetramer dissociation | ATTR-ACT (Maurer 2018 NEJM) |
| **Patisiran** (0.3mg/kg IV Q3W) | LNP-siRNA → ApoE-LDLR-mediated hepatic uptake → RISC/Ago2 → TTR mRNA cleavage | TTR mRNA ↓80% | APOLLO (Adams 2018 NEJM) |
| **Vutrisiran** (25mg SC Q3M) | GalNAc-siRNA → ASGR1 → RISC → mRNA degradation | TTR mRNA ↓83% | HELIOS-A (Gillmore 2021 NEJM) |
| **Inotersen** (300mg SC QW) | 2'-MOE ASO → RNase H1 → TTR mRNA cleavage | TTR mRNA ↓72% | NEURO-TTR (Benson 2018 Lancet) |
| **Acoramidis** (800mg PO BID) | Highly selective T4-site binding → tetramer stabilisation | Decreased TTR tetramer dissociation | ATTRiBUTE-CM (Elliott 2023 NEJM) |
| **Diflunisal** (250mg PO BID) | NSAID + weak T4-site binding | TTR stabilisation (non-selective) | Berk 2013 JAMA |

---

## Model File Composition

| File | Specification | Description |
|------|------|------|
| [`attr_qsp_model_en.dot`](attr_qsp_model_en.dot) | **116 nodes, 10 clusters** | Graphviz mechanistic map (DOT source) |
| [`attr_qsp_model_en.svg`](attr_qsp_model_en.svg) | Vector graphic | Zoomable SVG |
| [`attr_qsp_model_en.png`](attr_qsp_model_en.png) | 150 dpi PNG | Raster image for preview |
| [`attr_mrgsolve_model.R`](attr_mrgsolve_model.R) | **25-compartment ODE, 7 treatment scenarios** | mrgsolve PK/PD model |
| [`attr_shiny_app.R`](attr_shiny_app.R) | **8-tab dashboard** | Shiny interactive app |
| [`attr_references_en.md`](attr_references_en.md) | **60 PubMed citations (10 sections)** | Supporting reference list |

---

## mrgsolve Model Detail — 25 ODE Compartments

### Drug PK Compartments (9)

| Compartment | Drug | Route | Key Parameters |
|------|------|------|-------------|
| A_TAF_GUT, A_TAF_C, A_TAF_P | Tafamidis | Oral, 2-compartment | ka=0.42/h, CL=0.96L/h, t½~55h |
| A_VUT_SC, A_VUT_C | Vutrisiran | Subcutaneous, 1-compartment | ka=0.08/h, F=82%, t½~5d |
| A_INO_SC, A_INO_C | Inotersen | Subcutaneous, 1-compartment | ka=0.10/h, F=70%, CL=0.04L/h, t½~30d |
| A_PAT_C, A_PAT_P | Patisiran | Intravenous, 2-compartment | CL=0.18L/h, Vc=3.3L, t½~3d |

### Disease PD Compartments (16)

| Compartment | Biological Meaning | Key Dynamics |
|------|-------------|-----------|
| TTR_MRNA | Hepatic TTR mRNA (normalised) | kin(1-ERNA) - kout·mRNA; E_RNA = E_VUT + E_PAT + E_INO |
| TTR_TET | Plasma TTR tetramer | ksyn·mRNA - kdis(1-Estab)·TET - kout·TET |
| TTR_MONO | Misfolded monomer | 2·kdis·TET - kagg·MONO - kdeg·MONO |
| TTR_OLIGO | Toxic oligomer | kagg·MONO - kfib·OLIGO - kdeg·OLIGO |
| FIB_HRT | Cardiac amyloid burden | kfib·OLIGO·f_heart - kdeg_FIB·FIB_HRT |
| FIB_NRV | Peripheral nerve amyloid | kfib·OLIGO·f_nerve - kdeg_FIB·FIB_NRV |
| FIB_SYS | Systemic amyloid | kfib·OLIGO·f_sys - kdeg_FIB·FIB_SYS |
| INFLAM | Cardiac inflammation index | kin·FIB_HRT/(FIB50+FIB_HRT) - kout·INFLAM |
| LVEF | Left ventricular ejection fraction (%) | krec·(EF_base-EF) - kdet·FIB·INFLAM·EF |
| NT_proBNP | NT-proBNP (pg/mL) | kin·(INFLAM + 1/EF) - kout·NT_proBNP |
| NIS | Neuropathy impairment score | kin_NIS·FIB_NRV - kout_NIS·NIS |
| mBMI | Modified body mass index | -kdet·NIS/(NIS+NIS50)·mBMI |
| eGFR | Glomerular filtration rate | -kdet_eGFR·FIB_SYS·eGFR |
| SYMP_CARD | Cumulative cardiac symptoms | Integrated NYHA-like index |
| SYMP_NEURO | Cumulative neurological symptoms | Integrated FAP-stage-like index |

---

## Seven Treatment Scenarios

| # | Scenario | Drug/Dose | Population | Clinical Evidence |
|---|----------|----------|------|----------|
| S1 | Natural history, ATTRwt | None | ATTRwt-CM | Maurer 2018 (placebo arm) |
| S2 | Natural history, ATTRv | None | ATTRv-NP | Adams 2018 (placebo arm) |
| S3 | Tafamidis | 61mg PO QD | ATTRwt-CM | ATTR-ACT: 30% reduction in CV death + hospitalisation |
| S4 | Patisiran | 0.3mg/kg IV Q3W | ATTRv-NP | APOLLO: 34-point difference in mNIS+7 |
| S5 | Vutrisiran | 25mg SC Q3M | ATTRv-NP | HELIOS-A: 17-point improvement in NIS |
| S6 | Inotersen | 300mg SC QW | ATTRv-NP | NEURO-TTR: 19-point difference in mNIS+7 |
| S7 | Tafamidis + vutrisiran | Combination | Hypothetical | Dual mechanism (exploratory) |

---

## Shiny App Tab Structure (8)

| Tab | Contents |
|----|------|
| 1. Patient profile | ATTR phenotype/treatment settings, quick biomarker summary |
| 2. Drug PK | Visualisation of blood concentrations and drug effects (Estab, ERNA) for 4 drugs |
| 3. TTR misfolding | mRNA reduction, tetramer, oligomer, and fibril accumulation dynamics |
| 4. Cardiac outcomes | LVEF, NT-proBNP, cardiac fibril burden, estimated NYHA |
| 5. Neurological outcomes | NIS, mBMI, nerve fibril burden, eGFR |
| 6. Scenario comparison | Direct comparison of 5 treatment options, 18-month summary table |
| 7. Biomarker dashboard | Integrated cardiac/neurological biomarker panel, summary by timepoint |
| 8. Model information | ODE structure, parameter rationale, key references |

---

## Usage

```r
# 1) Run the mrgsolve ODE model
install.packages(c("mrgsolve","dplyr","ggplot2","tidyr","patchwork"))
source("transthyretin-amyloidosis/attr_mrgsolve_model.R")

# 2) Launch the Shiny dashboard
install.packages(c("shiny","shinydashboard","plotly","DT"))
shiny::runApp("transthyretin-amyloidosis/attr_shiny_app.R")

# 3) Render the Graphviz map
dot -Tsvg attr_qsp_model_en.dot -o attr_qsp_model_en.svg
dot -Tpng -Gdpi=150 attr_qsp_model_en.dot -o attr_qsp_model_en.png
```

---

## Summary of Key Clinical Trial Results

| Clinical Trial | Drug | N | Primary Outcome | Effect Size |
|---------|------|---|---------|---------|
| **ATTR-ACT** (2018) | Tafamidis 61mg | 441 | CV death + HF hospitalisation | HR 0.70 (95%CI 0.51-0.96) |
| **APOLLO** (2018) | Patisiran 0.3mg/kg | 225 | mNIS+7 improvement | -34 points vs. placebo (p<0.001) |
| **HELIOS-A** (2021) | Vutrisiran 25mg SC | 164 | mNIS+7 improvement | -17 points vs. placebo, adaptive comparison |
| **NEURO-TTR** (2018) | Inotersen 300mg | 172 | mNIS+7 improvement | -19 points vs. placebo (p<0.001) |
| **ATTRiBUTE-CM** (2023) | Acoramidis 800mg | 632 | 6MWT+KCCQ | p=0.0066 (hierarchical composite) |

---

*QSP Disease Model Library — a new model is added daily by Claude Code Routine.*
