# Sepsis & Septic Shock — QSP Model

> **Directory:** `sepsis/` | **Abbreviation:** SEP | **Date:** 2026-06-24
> **Category:** Acute Critical Care / Infectious Disease
> **Global Burden:** ~49 million cases/year; ~11 million deaths (22% of all global deaths)

[![SEP QSP Mechanistic Map](../../../sepsis/sep_qsp_model.png)](../../../sepsis/sep_qsp_model.svg)

---

## Disease Overview

**Sepsis** is a disease in which a dysregulated host response to infection leads to life-threatening organ dysfunction (Sepsis-3, Singer et al. 2016). Its pathophysiology goes beyond a simple systemic inflammatory response (SIRS) and is characterised by the simultaneous coexistence of immune hyperactivation (a pro-inflammatory storm) and immunosuppression (CARS).

**Septic shock** is a subset of sepsis defined by persistent hypotension (MAP <65 mmHg) together with a serum lactate >2 mmol/L, a severe state with a 28-day mortality exceeding 40%.

---

## Pathophysiology Summary

| Stage | Key mechanism | Key nodes |
|------|----------|----------|
| **1. Pathogen recognition** | LPS/PGN/DAMP → TLR4/TLR2/NLRP3 → NF-κB activation | PAMPs, DAMPs, PRRs, MyD88, NFkB |
| **2. Cytokine storm** | Excess secretion of TNFα, IL-1β, IL-6, IL-8 → secondary cell activation | TNF, IL6, IL1B, IL10, IFNγ |
| **3. Innate immune hyperactivation** | Neutrophil tissue infiltration, NETs, ROS, MMP-9 | Neut_T, NET, ROS, Elastase |
| **4. Complement activation** | C3→C5a (anaphylatoxin)→leukocyte/vascular response | C5a, C5aR, MAC |
| **5. Coagulation/DIC** | ↑ TF expression→↑ thrombin→↑ fibrin + ↑ PAI-1→DIC | Thrombin, Fibrin, PAI1 |
| **6. Endothelial dysfunction** | VE-cadherin disassembly · ↑ vascular permeability · ↑ NO → vasodilatory shock | ENDOT, VascPerm, MAP_vasc |
| **7. Multiple organ dysfunction syndrome (MODS)** | Lung (ARDS), kidney (AKI), liver, brain, circulatory failure | SOFA = 0–24 |
| **8. Late immunosuppression** | T-cell apoptosis · ↑ PD-1 · ↑ MDSC → CARS | Immunosuppression, CARS |

---

## Treatment Scenarios

| Scenario | Treatment | Mechanism | Key clinical trial |
|---------|------|------|-------------|
| **S1** | No treatment (natural course) | — | Observational cohort |
| **S2** | Antibiotic alone (meropenem 1g q8h) | Achieving fT>MIC → bacterial killing | Craig 1998 CID |
| **S3** | Antibiotic + norepinephrine | α1-adrenoceptor→restoration of MAP | SOAP II NEJM 2010 |
| **S4** | Bundle (antibiotic + NE + fluids 30 mL/kg) | RAAS correction · preload restoration | EGDT NEJM 2001 |
| **S5** | Bundle + hydrocortisone 200 mg/day | GR→cytokine suppression · ↑ vasopressor sensitivity | ADRENAL/APROCCHSS NEJM 2018 |
| **S6** | Bundle + HC + tocilizumab 8 mg/kg | IL-6R blockade→STAT3 inhibition→↓ CRP · ↓ PCT | REMAP-CAP NEJM 2021 |
| **S7** | Immunocompromised patient (antibiotic+NE, high bacterial burden) | Abnormal immune response, high morbidity/mortality | Clinical cohort |

---

## Model Specifications

| Component | File | Specification |
|---------|------|-----|
| 🗺️ Mechanistic map | [`sep_qsp_model.dot`](../../../sepsis/sep_qsp_model.dot) / [`.svg`](../../../sepsis/sep_qsp_model.svg) / [`.png`](../../../sepsis/sep_qsp_model.png) | **130+ nodes, 11 clusters** |
| ⚙️ mrgsolve ODE model | [`sep_mrgsolve_model.R`](../../../sepsis/sep_mrgsolve_model.R) | **24-compartment ODE**, **7 treatment scenarios** |
| 📊 Shiny dashboard | [`sep_shiny_app.R`](../../../sepsis/sep_shiny_app.R) | **8 tabs** (patient profile · antibiotic PK · cytokine/immune · haemodynamics/SOFA · organ function · scenario comparison · biomarkers · About) |
| 📚 References | [`sep_references.md`](../../../sepsis/sep_references.md) | **55 PubMed citations** (14 sections) |

---

## Mechanistic Map Clusters

| # | Cluster | Node count | Key content |
|---|---------|---------|---------|
| 1 | Pathogen Recognition | 20 | LPS, PGN, TLR4/2/5/9, NLRP3, NF-κB |
| 2 | Innate Immunity | 15 | Neutrophils, NETs, macrophages, ROS |
| 3 | Cytokine Storm | 18 | TNFα, IL-1β, IL-6, IL-8, IL-10, HMGB1 |
| 4 | Complement | 12 | C3, C5a, MAC, C5aR |
| 5 | Coagulation | 15 | TF, thrombin, fibrin, PAI-1, APC |
| 6 | Endothelial | 16 | VE-cadherin, ICAM-1, NO, oedema |
| 7 | Organ Failure/MODS | 22 | ARDS, AKI, SOFA 6-domain |
| 8 | Drug PK/PD | 15 | Abx, NE, hydrocortisone, tocilizumab |
| 9 | Biomarkers | 13 | PCT, CRP, lactate, SOFA |
| 10 | Metabolic | 8 | Lactate, ROS, mitochondria |
| 11 | Adaptive/CARS | 11 | Treg, MDSC, PD-1, apoptosis |

---

## mrgsolve Model Compartments (24 ODE Compartments)

| Group | Compartments | Description |
|------|------|------|
| **Pathogen** | BACT | Bacterial burden (CFU/mL, logistic growth + antibiotic killing) |
| **Antibiotic PK** | ABX1, ABX2 | 2-compartment meropenem PK (central · tissue) |
| **Cytokines** | TNF, IL6, IL10, IL1B | The four major cytokines (Emax-inhibition feedback) |
| **Immune cells** | NEUT_B, NEUT_T, MACS | Blood/tissue neutrophils, activated macrophages |
| **Complement** | C5A | C5a effector |
| **Coagulation** | THROMBIN, FIBRIN, PAI1 | Thrombin · fibrin · PAI-1 |
| **Endothelium** | ENDOT | Endothelial injury index (0–1) |
| **Organs** | PF_RATIO, CREATININE, BILIRUBIN, LACTATE, MAP_val, PLT_COUNT | Based on the 6 SOFA domains |
| **Drug PK** | NE_C, HC_C, TOCI_C | Norepinephrine · hydrocortisone · tocilizumab |

---

## SOFA Calculation (Based on Sepsis-3)

| Domain | Variable | 0 points | 1 point | 2 points | 3 points | 4 points |
|--------|------|-----|-----|-----|-----|-----|
| Lung | PaO₂/FiO₂ | ≥400 | <400 | <300 | <200 | <100 |
| Kidney | Creatinine (mg/dL) | <1.2 | 1.2–2.0 | 2.0–3.5 | 3.5–5.0 | >5.0 |
| Liver | Bilirubin (mg/dL) | <1.2 | 1.2–2.0 | 2.0–6.0 | 6.0–12.0 | >12.0 |
| Cardiovascular | MAP + vasopressor | ≥70 | <70 | MAP<65 | NE low dose | NE high dose |
| Coagulation | Platelets (×10⁹/L) | ≥150 | <150 | <100 | <50 | <20 |
| CNS | Consciousness/encephalopathy proxy | Normal | Mild | Moderate | Severe | Coma |
| **Total** | **SOFA 0–24** | — | — | **≥2 = Sepsis** | — | **≥11 ~40% mortality** |

---

## How to Run

```r
# 1) Run the mrgsolve model
install.packages(c("mrgsolve", "dplyr", "ggplot2", "tidyr"))
source("sep_mrgsolve_model.R")   # automatically simulates the 7 scenarios
plot(p_sofa)                     # SOFA score comparison
plot(p_mort)                     # 28-day mortality
print(summary_72h)               # 72h summary table

# 2) Run the Shiny dashboard
install.packages("shiny", "shinydashboard")
shiny::runApp("sep_shiny_app.R")

# 3) Render the mechanistic map (Graphviz)
# dot -Tsvg sep_qsp_model.dot -o sep_qsp_model.svg
# dot -Tpng -Gdpi=150 sep_qsp_model.dot -o sep_qsp_model.png
```

---

## Key Clinical Parameter Calibration Notes

| Parameter | Value | Basis |
|---------|-----|------|
| Bacterial growth rate | kgrow = 1.2/h | E. coli doubling ~35 min in vivo |
| Antibiotic MIC | 0.5 mg/L (meropenem) | EUCAST breakpoint |
| TNFα t1/2 | ~70 min | Beutler 1985 (kdeg=0.6/h) |
| IL-6 t1/2 | ~4.6 h | Taniguchi 1999 (kdeg=0.15/h) |
| Meropenem CL | 10 L/h | Roberts 2009 critical care PK |
| Hydrocortisone Emax | 0.65 (65% cytokine suppression) | Annane 2002 in vitro |
| Tocilizumab EC50 | 1.5 mcg/mL | IL-6R Kd ~1-2 μg/mL |
| SOFA mortality | logit = -6.5 + 0.45×SOFA | Seymour 2017 JAMA external validation |
| Antibiotic delay effect | ~7%/h survival decrease | Kumar 2006 Crit Care Med |

---

*Generated by Claude Code Routine (CCR) — 2026-06-24*
*For educational and research purposes. Not for clinical use.*
