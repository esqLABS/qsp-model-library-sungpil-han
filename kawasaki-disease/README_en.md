# Kawasaki Disease (KD) — QSP Model

[![Model](kd_qsp_model.png)](kd_qsp_model.svg)

## Overview

**Kawasaki disease (KD)** is an acute systemic vasculitis of unknown cause that occurs
mainly in children under 5 years of age, and is **the most common cause of acquired
heart disease in children** in developed countries. If not treated adequately, it leads
to the formation of coronary artery aneurysms (CAA), causing a lifelong risk of
ischemic heart disease.

**Core pathophysiology:**  
An unidentified infectious trigger → activation of innate immunity (macrophages,
neutrophils) → NLRP3 inflammasome → cytokine storm (IL-1β, IL-6, TNF-α) → vascular
endothelial damage → destruction of the coronary artery media → aneurysm formation

---

## Disease Features

| Item | Content |
|------|------|
| **Classification** | Pediatric Systemic Vasculitis |
| **Prevalence** | Japan ~300/100,000, Korea ~200/100,000, US ~20/100,000 (under age 5) |
| **Peak age** | 6 months–5 years (median ~2 years) |
| **Diagnostic criteria** | AHA 2017: fever for ≥5 days + at least 4 of 5 clinical criteria |
| **CAA incidence** | ~25% without treatment, ~3–5% after IVIG |
| **IVIG resistance** | About 10–20% |
| **Recurrence rate** | About 3% |

### Key Pathways

```
Unknown Trigger
  ↓
TLR/NLR Activation → Macrophage Activation
  ↓
NLRP3 Inflammasome → Caspase-1 → IL-1β (mature)
  ↓
Cytokine Storm: IL-1β + IL-6 + TNF-α
  ↓
Endothelial Activation (VCAM-1, ICAM-1, Tissue Factor ↑)
  ↓
Coronary Artery Inflammation → Medial Destruction → Aneurysm
  ↓
Thrombocytosis (peak Week 2–3) → Thrombosis Risk
  ↓
CAA: small (z 2.5–5) → medium (z 5–10) → giant (z ≥ 10)
```

---

## Treatment Algorithm

```
Diagnosis confirmed
├── Standard-risk (Kobayashi score < 4)
│   └── IVIG 2 g/kg × 1 dose + Aspirin 80–100 mg/kg/day
│       ├── Responder (80–90%) → Switch to low-dose aspirin (3–5 mg/kg/day)
│       └── Non-responder → IVIG-resistant protocol
│
├── High-risk (Kobayashi score ≥ 4 or Egami score ≥ 3)
│   └── IVIG + Aspirin + Prednisolone 2 mg/kg/day × 5–15 days
│
└── IVIG-Resistant (~10–20%)
    ├── Option A: 2nd IVIG 2 g/kg
    ├── Option B: Infliximab 5 mg/kg IV (TNF-α blockade)
    ├── Option C: IV Methylprednisolone pulse
    └── Option D: Anakinra 4 mg/kg/day SC (IL-1 blockade)
```

---

## Model Structure

### Mechanistic Map

| Component | Content |
|-----------|------|
| **Total nodes** | 134 |
| **Number of clusters** | 14 |
| **Cluster list** | Infectious trigger, innate immunity, NLRP3 inflammasome, cytokine network, adaptive immunity, endothelial cell activation, coronary artery pathology, platelet biology, fever/acute-phase response, IVIG PK, aspirin PK, corticosteroid PK, biologic agent PK, clinical endpoints |

### mrgsolve ODE Model

| Compartment | Description |
|------|------|
| A_IVIG_c/p | IVIG central/peripheral compartments (2-compartment + FcRn recirculation) |
| A_ASA_gut/c, A_SA_c | Aspirin absorption/plasma/salicylic acid metabolite |
| A_MP_c/p | Methylprednisolone 2-compartment |
| A_IFX_c/p | Infliximab 2-compartment |
| A_ANK_gut/c | Anakinra SC absorption + plasma |
| IL1b, IL6, TNFa | Cytokine dynamics |
| Mac_act, EC_act | Macrophage/endothelial cell activation |
| Fever | Body temperature (fever dynamics) |
| CRP | C-reactive protein |
| PLT_c | Platelet count (thrombocytosis) |
| CAL_Z | Coronary artery Z-score |

**A total of 21 compartments, 5 treatment scenarios**

### 5 Treatment Scenarios

| Scenario | Description |
|----------|------|
| **S1** | IVIG 2 g/kg + high-dose aspirin → low-dose aspirin (standard treatment) |
| **S2** | IVIG + aspirin + methylprednisolone (Kobayashi high-risk group) |
| **S3** | IVIG resistance → second IVIG dose |
| **S4** | IVIG resistance → infliximab rescue therapy |
| **S5** | IVIG resistance → anakinra rescue therapy (IL-1 blockade) |

### Shiny Dashboard (6 tabs)

| Tab | Content |
|----|------|
| 1. Patient Profile | Diagnostic criteria, treatment algorithm, value boxes |
| 2. Drug Kinetics (PK) | IVIG / aspirin / steroid / biologic agent PK curves |
| 3. Cytokines/Inflammation | IL-1β, IL-6, TNF-α, macrophage/endothelial cell activation |
| 4. Clinical Endpoints | Fever, CRP, coronary artery Z-score, platelets |
| 5. Scenario Comparison | Comparison of 5 treatment scenarios |
| 6. Biomarkers/Risk | Kobayashi score calculator, CAA risk probability |

---

## Files

| File | Description |
|------|------|
| `kd_qsp_model.dot` | Graphviz mechanistic map source (134 nodes, 14 clusters) |
| `kd_qsp_model.svg` | Vector image |
| `kd_qsp_model.png` | Raster image (150 dpi) |
| `kd_mrgsolve_model.R` | mrgsolve ODE QSP model (5 scenarios) |
| `kd_shiny_app.R` | Shiny interactive dashboard (6 tabs) |
| `kd_references.md` | References (60 PubMed links, 14 sections) |

---

## Usage

```bash
# Render the mechanistic map
dot -Tsvg kd_qsp_model.dot -o kd_qsp_model.svg
dot -Tpng -Gdpi=150 kd_qsp_model.dot -o kd_qsp_model.png
```

```r
# Run the mrgsolve model
install.packages(c("mrgsolve", "dplyr", "ggplot2", "tidyr"))
source("kd_mrgsolve_model.R")
results <- run_all_scenarios(wt = 15)

# Run the Shiny app
install.packages(c("shiny", "shinydashboard", "plotly", "DT"))
shiny::runApp("kd_shiny_app.R")
```

---

## Key Parameters

| Parameter | Value | Source |
|----------|----|------|
| IVIG dose | 2 g/kg IV | AHA 2017 |
| IVIG t½ | 21–28 days | Tremoulet 2015 |
| Infliximab EC50 | 2.5 μg/mL | Tremoulet 2020 |
| Anakinra EC50 | 1.0 μg/mL | Ouldali 2019 |
| Baseline IL-6 | 2.5 ng/mL | Matsubara 2013 |
| Coronary artery Z-score threshold | ≥ 2.5 (CAL), ≥ 10 (giant aneurysm) | AHA 2017 |

---

## References Summary

60 references across 14 sections: epidemiology (5) · pathophysiology (5) · cytokines
(6) · coronary artery pathology (5) · IVIG pharmacology (6) · risk scores (4) ·
corticosteroids (3) · infliximab (3) · anakinra/biologic agents (3) ·
platelets/aspirin (3) · echocardiography/Z-score (3) · QSP modeling (4) ·
COVID-19/MIS-C (3) · long-term prognosis (5)

See [`kd_references.md`](kd_references.md) for details.
