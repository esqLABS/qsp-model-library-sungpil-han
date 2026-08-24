# Acute Kidney Injury (AKI) QSP Model

[![AKI Map](aki_qsp_model.png)](aki_qsp_model.svg)

**Category:** Nephrology / Critical Care  
**Directory:** [`acute-kidney-injury/`](.)

---

## Pathophysiology

Acute kidney injury (AKI) is a clinical syndrome defined as **a rise in serum creatinine of ≥0.3 mg/dL within 48 hours, or a rise to ≥1.5 times baseline within 7 days, or urine output <0.5 mL/kg/h for 6 hours or more** (KDIGO 2012). It occurs in 10–15% of hospitalised patients worldwide and in more than 50% of ICU patients, and increases in-hospital mortality by up to 5–10-fold.

### Three Subtypes

| Subtype | Key mechanism | Representative causes |
|------|----------|----------|
| **Ischaemia-reperfusion injury (IRI)** | ATP depletion → mitochondrial dysfunction → apoptosis/necrosis/ferroptosis | Major surgery, cardiogenic shock, hypovolaemia |
| **Nephrotoxicity (NTX)** | Direct tubular cell toxicity, oxidative stress, GSH depletion | Cisplatin, aminoglycosides, contrast media |
| **Sepsis-associated (SA-AKI)** | LPS/DAMPs → TLR4 → NF-κB → microvascular dysfunction + cytokine storm | Sepsis, septic shock |

### Key Pathophysiological Pathways

| Pathway | Mechanism |
|------|---------|
| **Microvascular dysfunction** | ET-1↑, Ang II↑, NO↓, TXA2↑ → afferent arteriolar constriction → GFR↓ → medullary hypoxia |
| **Tubular cell injury** | ATP depletion → mPTP opening → Cyt-c release → Casp-9/3 → apoptosis; GPx4↓ → lipid peroxidation → ferroptosis |
| **Oxidative stress** | ROS↑ (mitochondrial/NOX) → GSH depletion → insufficient Nrf2/HO-1 counter-response → cell death |
| **Inflammatory storm** | DAMPs → TLR4 → NF-κB → IL-6, TNF-α, IL-1β↑ → neutrophil infiltration → NET formation |
| **Tubular obstruction and backleak** | Brush-border loss → intraluminal casts → obstruction → filtrate backleak |
| **Failed tubular repair** | G2/M arrest → TGF-β1↑ → myofibroblast activation → fibrosis → AKI-to-CKD |

---

## KDIGO AKI Staging

| Stage | Serum creatinine | Urine output |
|-------|--------------|-------|
| **1** | ×1.5–1.9 baseline or +0.3 mg/dL | <0.5 mL/kg/h for 6–12h |
| **2** | ×2.0–2.9 | <0.5 mL/kg/h for ≥12h |
| **3** | ×3.0 or more, or dialysis | <0.3 mL/kg/h for ≥24h or anuria for ≥12h |

---

## Deliverables

| File | Description |
|------|------|
| `aki_qsp_model.dot` | Graphviz mechanistic map (12 clusters, 100+ nodes) |
| `aki_qsp_model.svg` | Vector graphic (interactive view) |
| `aki_qsp_model.png` | Raster image (150 dpi) |
| `aki_mrgsolve_model_en.R` | mrgsolve ODE model (20 compartments, 7 scenarios) |
| `aki_shiny_app_en.R` | Shiny interactive dashboard (6 tabs) |
| `aki_references.md` | 55 references (classified by section) |

---

## Model Specifications

### Mechanistic map
- **12 subgraph clusters**: ① triggering factors/risk factors, ② renal microvascular dysfunction, ③ tubular cell injury (apoptosis/necrosis/ferroptosis), ④ oxidative stress (ROS/GSH/Nrf2), ⑤ innate immunity/inflammation (NF-κB/NLRP3/IL-6/TNF-α), ⑥ glomerular filtration and urinary dynamics, ⑦ AKI biomarkers, ⑧ furosemide PK (OAT1/3/NKCC2), ⑨ vasopressors/NAC/CRRT PK/PD, ⑩ tubular repair/regeneration, ⑪ maladaptive repair → AKI-to-CKD, ⑫ clinical outcomes
- Includes **100+ nodes**

### mrgsolve ODE Model (20 Compartments)
| Module | Compartments |
|------|------|
| **Drug PK** | Furosemide central/peripheral/gut · Norepinephrine central · NAC central |
| **AKI pathophysiology** | ATP · ROS · GSH · tubular cell viability (TCV) · GFR · IL-6 · TNF-α |
| **Biomarkers** | NGAL · KIM-1 · Serum creatinine · Cystatin C |
| **Repair/fibrosis** | Repair Capacity · TGF-β1 · Myofibroblast · Fibrosis Index |

### Treatment Scenarios (7)
1. **IRI — untreated**: natural-history control group
2. **IRI + furosemide 40mg IV q12h**: diuretic strategy (NKCC2 inhibition)
3. **IRI + norepinephrine**: blood pressure support → restoration of RBF autoregulation
4. **IRI + prophylactic NAC**: GSH replenishment → ROS clearance → reduced↓ oxidative damage
5. **SA-AKI — untreated**: natural history of sepsis
6. **SA-AKI + NE + furosemide + CRRT**: combined intervention strategy
7. **Nephrotoxic AKI (cisplatin)**: long-term monitoring of AKI-to-CKD transition

### Clinical Trial Calibration Parameters
| Reference trial | Calibration data |
|-------------|-----------|
| Mehta et al. Lancet 2015 | KDIGO staging criteria (sCr×1.5/2.0/3.0) |
| Gaudry et al. NEJM 2016 | Early CRRT vs delayed — 90d mortality ~48% |
| Zarbock et al. JAMA 2016 | Early RRT initiation benefit in AKI-3 |
| Mishra et al. JASN 2003 | NGAL rise within 2h of IRI |
| Han et al. Kidney Int 2002 | KIM-1 shedding kinetics |
| Meersch et al. ICM 2017 | TIMP-2·IGFBP7 prediction cutoff |
| Felker et al. NEJM 2011 | Furosemide dose-response in decongestion |

### Shiny Dashboard (6 Tabs)
| Tab | Content |
|----|------|
| ① Patient Profile | AKI subtype selection, risk factors, KDIGO staging reference table, overview graph |
| ② Drug PK | Furosemide blood concentration, NKCC2 inhibition rate, NE/NAC concentration, diuretic response |
| ③ Renal Biomarkers | NGAL, KIM-1, Cystatin C vs Cr, tubular cell viability, ROS/GSH |
| ④ Clinical Endpoints | eGFR trajectory, AKI stage time series, urine output, IL-6/TNF-α |
| ⑤ Scenario Comparison | GFR/sCr/NGAL comparison across all 7 scenarios + clinical summary table |
| ⑥ AKI-to-CKD | TGF-β1, myofibroblasts, fibrosis index, long-term (30-day) GFR, CKD risk |

---

## Usage

```r
# 1) Run the mrgsolve model
install.packages(c("mrgsolve", "dplyr", "ggplot2", "tidyr"))
source("aki_mrgsolve_model_en.R")

# 2) Run the Shiny app
install.packages(c("shiny", "shinydashboard", "mrgsolve", "DT"))
shiny::runApp("aki_shiny_app_en.R")
```

```bash
# 3) Re-render the mechanistic map
sfdp -Tsvg -Goverlap=prism aki_qsp_model.dot -o aki_qsp_model.svg
sfdp -Tpng -Goverlap=prism -Gdpi=150 aki_qsp_model.dot -o aki_qsp_model.png
```

---

## Key Drug PK/PD Summary

| Drug | Model characteristics |
|------|---------|
| **Furosemide** | 2-compartment PK · OAT1/3 active secretion → determines intratubular concentration · NKCC2 IC50=0.5 mg/L · reduced↓ OAT expression in AKI → diuretic resistance |
| **Norepinephrine** | 1-compartment PK (t½ ~2.5 min) · α1-receptor → MAP↑ → afferent arteriolar autoregulation → RBF recovery |
| **N-Acetylcysteine** | 1-compartment PK · cysteine precursor → GSH replenishment · thiol-NO exchange → eNOS↑ |
| **CRRT** | Creatinine clearance of 3 L/h (~50 mL/min) · IL-6 adsorption (0.5 L/h) · fluid balance regulation |

---

## References Summary

55 PubMed references, classified into 16 sections:
- Epidemiology/burden (4), KDIGO criteria (2), IRI pathophysiology (5), nephrotoxicity (3), sepsis-AKI (3),
  oxidative stress/mitochondria (4), inflammation (3), biomarkers (5), furosemide (3),
  vasopressors (2), NAC (2), CRRT (3), AKI-to-CKD (5), QSP modelling (4),
  clinical trials/guidelines (4), biomarker-guided prevention (3)

Full list: [`aki_references.md`](aki_references.md)

---

## Disclaimer

This model is a qualitative/semi-quantitative QSP model intended for educational and research purposes, and must not be used directly for clinical decision-making.
