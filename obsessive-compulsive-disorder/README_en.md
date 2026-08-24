# OCD — Obsessive-Compulsive Disorder QSP Model

> A Quantitative Systems Pharmacology (QSP) model of **Obsessive-Compulsive Disorder (OCD)**: a simulation integrating CSTC circuit dynamics, SERT-occupancy-based SSRI/clomipramine PK/PD, ERP treatment effects, and Y-BOCS clinical endpoints

---

## Mechanistic Map

[![OCD QSP Mechanistic Map](ocd_qsp_model.png)](ocd_qsp_model.svg)

*Click to view an enlargeable SVG version.*

---

## Overview

**Obsessive-Compulsive Disorder (OCD)** is a chronic neuropsychiatric disorder occurring in approximately 2–3% of the global population, characterised by repeated, unwanted intrusive thoughts (obsessions) and repetitive behaviours performed to relieve them (compulsions).

### Core pathophysiology

| Component | Role |
|-----------|------|
| **OFC (orbitofrontal cortex) hyperactivity** | Intrusive thoughts, generation of error signals |
| **Caudate nucleus hypermetabolism** | Maintenance of compulsive urges |
| **Direct-pathway (Go) dominance** | Thalamic disinhibition → circuit overload |
| **Indirect-pathway (No-Go) weakening** | Failure of behavioural inhibition |
| **Reduced SERT function** | Decreased synaptic 5-HT |
| **OFC → caudate glutamate hyperactivity** | Increased cortico-striatal drive |

---

## File Structure

| File | Content |
|------|------|
| `ocd_qsp_model.dot` | Graphviz DOT mechanistic map (130+ nodes, 13 clusters) |
| `ocd_qsp_model.svg` | Vector graphic (scalable) |
| `ocd_qsp_model.png` | Raster image (150 dpi) |
| `ocd_mrgsolve_model.R` | mrgsolve ODE PK/PD model (6 treatment scenarios, 100 virtual patients) |
| `ocd_shiny_app.R` | Interactive Shiny dashboard (8 tabs) |
| `ocd_references_en.md` | 50 references (including PubMed links) |
| `README.md` | This file |

---

## Model Architecture

### Mechanistic map — 13 clusters

```
1. SSRI PK           ← sertraline 2-compartment + BBB transport + CYP2D6/3A4 metabolism
2. Clomipramine PK   ← first-pass effect + active metabolite (DCMI)
3. Augmentation PK   ← risperidone, aripiprazole, memantine, D-cycloserine
4. Serotonin system  ← TPH2/AADC synthesis → SERT → MAO-A → 5-HT1A/2A/2C receptors
5. Dopamine system   ← VTA/SNc → striatum/OFC → D1R (direct) / D2R (indirect)
6. Glutamate/GABA    ← OFC-striatal glutamate synapse → NMDA/AMPA/mGluR
7. CSTC circuit      ← OFC → caudate → GPi/GPe/STN → thalamus → OFC loop
8. HPA axis/stress   ← CRH → ACTH → cortisol → hippocampal atrophy → anxiety loop
9. Neuroinflammation ← microglia → TNF-α/IL-6 → IDO → kynurenine/QUIN
10. Neuroplasticity   ← BDNF → TrkB → ERK/AKT → dendritic remodelling → LTP/LTD
11. CBT/ERP           ← fear extinction → OFC normalisation → habit-circuit reformation
12. Clinical endpoints ← Y-BOCS total score/subscales → response/remission criteria
13. Genetics/biomarkers ← SLC6A4/HTR2A/COMT/SLC1A1 variants → imaging/PET biomarkers
```

### mrgsolve ODE model — 23 compartments

```
Pharmacokinetics (PK):
  AG_SSRI → A1_SSRI → A2_SSRI / A_CNS     (sertraline 2-compartment + CNS)
  AG_CMI  → A1_CMI  → A_DCMI               (clomipramine + active metabolite)
  AG_RISP → A1_RISP                         (risperidone)

Pharmacodynamics (PD):
  SERT_OCC   ← Emax model of SSRI/CMI CNS concentration
  HT5_SYN    ← synaptic 5-HT (release-reuptake-degradation dynamics)
  DES_5HT1   ← 5-HT1A autoreceptor desensitisation (delayed-response mechanism)
  OFC_ACT    ← OFC activity (CSTC loop + 5-HT modulation)
  CAUD_ACT   ← caudate nucleus activity
  THAL_ACT   ← thalamic activity
  DIR_PATH   ← direct-pathway activation
  IND_PATH   ← indirect-pathway activation
  D2R_OCC    ← antipsychotic D2R occupancy
  BDNF_LV    ← BDNF level (neuroplasticity marker)
  ERP_EFF    ← cumulative ERP effect (0–1)
  YBOCS      ← Y-BOCS score (0–40, clinical endpoint)
  ANXIETY    ← anxiety state (normalised)
```

---

## Treatment Scenarios (6)

| Scenario | Treatment | Characteristics |
|----------|--------|------|
| 1 | Untreated (baseline) | Y-BOCS = 28 sustained |
| 2 | **Sertraline 200 mg/day** | Target SERT ≥80%, 6–12-week delayed effect |
| 3 | **Clomipramine 250 mg/day** | Strongest SERT inhibition, increased↑ side effects |
| 4 | Sertraline + **Risperidone 1.5 mg** | D2 augmentation after 12 weeks of insufficient SSRI response |
| 5 | Sertraline + **ERP** (combination) | Highest efficacy, dual OFC-normalisation mechanism |
| 6 | **ERP alone** | Psychotherapy alone, reduced↓ remission rate |

---

## Key Parameters

| Parameter | Value | Basis |
|----------|----|----|
| SERT EC50 (sertraline) | 1.2 ng/mL | Zitterl et al. 2008 |
| SERT EC50 (clomipramine) | 0.4 ng/mL | Estimated value (CMI more potent) |
| SERT occupancy threshold required | ≥80% | OCD vs depression (60%) |
| 5-HT1A desensitisation t½ | ~4 weeks | Explains SSRI delayed response |
| Y-BOCS response criterion | ≥35% reduction | Goodman et al. 1989 |
| Remission criterion | Y-BOCS ≤12 | APA guidelines |
| SSRI response rate | 40–60% | Soomro et al. 2008 |
| Maximum ERP-induced OFC suppression | 45% | Based on Foa et al. 2005 |

---

## How to Run

### mrgsolve model

```r
library(mrgsolve)
library(dplyr)
library(ggplot2)

# Run the model (6 treatment scenarios)
source("ocd_mrgsolve_model.R")
```

### Shiny dashboard

```r
library(shiny)
library(shinydashboard)
library(mrgsolve)
library(plotly)
library(DT)

shiny::runApp("ocd_shiny_app.R")
```

### Graphviz rendering

```bash
dot -Tsvg ocd_qsp_model.dot -o ocd_qsp_model.svg
dot -Tpng -Gdpi=150 ocd_qsp_model.dot -o ocd_qsp_model.png
```

---

## Shiny Dashboard Tabs (8)

| Tab | Content |
|----|------|
| 🧑 **Patient Profile** | Patient characteristics, pharmacogenomics (CYP2D6/SLC6A4/COMT), severity |
| 💊 **PK — Drug Levels** | SSRI/clomipramine blood and brain concentrations, SERT occupancy curve |
| 🧠 **PD — Neurotransmitters** | Synaptic 5-HT dynamics, 5-HT1A desensitisation, SERT-5HT correlation |
| 🔄 **CSTC Circuit** | OFC/caudate/thalamic activity, direct/indirect pathway balance, BDNF |
| 📊 **Clinical Endpoints** | Y-BOCS trajectory, response/remission valueBox, anxiety state |
| ⚖️ **Scenario Comparison** | Side-by-side comparison of 6 scenarios, results table |
| 🔬 **Biomarkers** | Imaging biomarkers, pharmacogenomic effects, DBS simulation |
| ℹ️ **About** | Model description, limitations, key references |

---

## Clinical Significance

- **Need for high-dose SSRI**: unlike depression (SERT 60%), OCD requires ≥80% SERT occupancy → titrate up to the maximum tolerated dose
- **Delayed response mechanism**: 5-HT1A autoreceptor desensitisation takes 2–4 weeks, delaying clinical effect → caution against premature discontinuation
- **ERP + SSRI combination**: superior to either treatment alone, via distinct mechanisms (drug: ↑5-HT / ERP: OFC downregulation)
- **Augmentation therapy**: when SSRI monotherapy is insufficient, risperidone (D2 blockade) restores the indirect pathway

---

## Key References

- Soomro GM et al. (2008) Cochrane Review: SSRIs for OCD
- Foa EB et al. (2005) JAMA: ERP vs Clomipramine RCT
- Zitterl W et al. (2008) Neuropsychopharmacology: SERT occupancy
- Goodman WK et al. (1989) Arch Gen Psychiatry: Y-BOCS development
- Saxena S & Rauch SL (2004) Annu Rev Neurosci: CSTC circuit
- Bloch MH et al. (2006) Mol Psychiatry: augmentation meta-analysis

→ Full set of 50 references: [`ocd_references_en.md`](ocd_references_en.md)

---

*Generated by Claude Code Routine (CCR) — 2026-06-28*
