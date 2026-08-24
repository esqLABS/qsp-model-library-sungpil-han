# CKD-Mineral Bone Disorder (CKD-MBD) QSP Model

## Overview

Chronic Kidney Disease–Mineral Bone Disorder (CKD-MBD) is a systemic syndrome characterized by dysregulation of the FGF23–Klotho–PTH–vitamin D axis, renal osteodystrophy, and vascular calcification as GFR declines. This QSP model captures these mechanisms in 17 ODE compartments with 7 treatment scenarios.

---

## Files

| File | Description |
|------|------|
| `ckdmbd_qsp_model.dot` | Graphviz mechanistic map source (115+ nodes, 10 clusters) |
| `ckdmbd_qsp_model.svg` | Mechanistic map (vector) |
| `ckdmbd_qsp_model.png` | Mechanistic map (150 dpi) |
| `ckdmbd_mrgsolve_model.R` | mrgsolve ODE model (17 compartments, 7 scenarios) |
| `ckdmbd_shiny_app.R` | Shiny interactive dashboard (7 tabs) |
| `ckdmbd_references.md` | 42 references (9 sections) |

---

## Mechanistic Map

[![CKD-MBD QSP Mechanistic Map](../../../ckd-mineral-bone-disorder/ckdmbd_qsp_model.png)](../../../ckd-mineral-bone-disorder/ckdmbd_qsp_model.svg)

---

## Model Architecture

### 10 Subgraph Clusters

| Cluster | Key nodes |
|---------|-----------|
| ① Kidney Function | GFR, CKD_Stage, NaPi-IIa/IIc, CYP27B1, CYP24A1, Klotho |
| ② Phosphate & FGF23 | Pi_plasma, FGF23_plasma, FGF23_bone, PHEX, DMP1, GALNT3 |
| ③ Vitamin D Axis | VitD3_skin, 25-OH-D, 1,25-OH₂D, VDR, RXR, VDRE, CaBP9k |
| ④ Calcium Homeostasis | Ca_plasma, Ca_ionized, CaSR_kidney, TRPV5, TRPV6, calcitonin |
| ⑤ Parathyroid Gland | PTH_plasma, CaSR_ptg, PTG_mass, SecHPT, TerHPT, nodular |
| ⑥ Bone Remodeling | Osteoblast, osteoclast, RANKL, OPG, BMD, sclerostin, DKK1 |
| ⑦ Vascular Calcification | VSMC, MGP, ucMGP, fetuin-A, hydroxyapatite, carotid IMT |
| ⑧ Drug PK | Sevelamer, cinacalcet, paricalcitol, etelcalcetide, denosumab |
| ⑨ Drug PD | CaSR_act, VDR_act, Pi_bind_GI, RANKL_inh, Emax models |
| ⑩ Clinical Outcomes | iPTH_lab, Pi_lab, Ca_lab, FGF23_lab, DXA_BMD, CV_event |

---

## ODE Compartments (17 ODEs)

| # | Compartment | Description |
|---|------|------|
| 1 | Pi | Serum phosphate (mg/dL) |
| 2 | FGF23 | Plasma FGF23 (pg/mL) |
| 3 | Klotho | Soluble Klotho (relative units) |
| 4 | PTH | Intact PTH (pg/mL) |
| 5 | VitD25 | 25-OH-vitamin D (nmol/L) |
| 6 | VitD_act | 1,25-OH₂D calcitriol (pg/mL) |
| 7 | Ca | Serum calcium (mg/dL) |
| 8 | OB | Osteoblast activity |
| 9 | OC | Osteoclast activity |
| 10 | BMD | Bone mineral density (relative to normal) |
| 11 | VC | Vascular calcification score |
| 12 | CIN_GUT | Cinacalcet gut (mg) |
| 13 | CIN_PLASMA | Cinacalcet plasma (ng/mL) |
| 14 | PAR_PLASMA | Paricalcitol plasma (ng/mL) |
| 15 | ETEL_PLASMA | Etelcalcetide plasma (ng/mL) |
| 16 | DEN_DEPOT | Denosumab SC depot (mg) |
| 17 | DEN_PLASMA | Denosumab plasma (mg/L) |

---

## Treatment Scenarios

| Scenario | Treatment | Represents |
|---------|------|------|
| S1 | No treatment (natural course of CKD G5) | Baseline |
| S2 | Sevelamer 2400 mg/day | Phosphate binder |
| S3 | Cinacalcet 60 mg/day | Calcimimetic |
| S4 | Paricalcitol 4 mcg 3x/week IV | Vitamin D receptor agonist |
| S5 | Sevelamer + cinacalcet | Combination |
| S6 | Etelcalcetide 5 mg 3x/week IV | Next-generation calcimimetic |
| S7 | Sevelamer + etelcalcetide + denosumab | Triple therapy |

---

## KDIGO 2017 Treatment Targets

| Index | Target range |
|------|----------|
| Intact PTH (iPTH) | 150–600 pg/mL (CKD G5D) |
| Serum phosphate | < 5.5 mg/dL |
| Serum calcium | 8.4–10.2 mg/dL |
| Ca×Pi product | < 55 mg²/dL² |
| 25-OH-vitamin D | > 75 nmol/L |

---

## Shiny App Tabs

1. **Patient Profile** — CKD stage, baseline lab values, simulation duration
2. **Drug PK** — PK profiles for phosphate binders · calcimimetics · vitamin D · denosumab
3. **PTH & Minerals** — iPTH, Pi, Ca, Ca×Pi, FGF23, vitamin D, KDIGO target attainment
4. **Bone Disease** — bone mineral density, osteoblasts/osteoclasts, P1NP/CTX, fracture risk
5. **Cardiovascular** — vascular calcification, GFR trend, Klotho, CV risk assessment
6. **Scenario Comparison** — comparison of 5 standard treatment strategies (iPTH, Pi, BMD, VC)
7. **Biomarkers** — normalised dashboard, spider chart, FGF23-Klotho axis

---

## Key References

- KDIGO CKD-MBD Guidelines 2017 · Gutierrez 2008 (FGF23-mortality, NEJM) · Hu 2011 (Klotho-VC, JASN)
- Block 2004 (mineral mortality, JASN) · Tentori 2008 (DOPPS) · Chertow 2002 (cinacalcet ACHIEVE)
- Luo 1997 (MGP knockout calcification, Nature) · Ketteler 2003 (fetuin-A, Lancet)
- Peterson & Riggs 2010 (calcium-bone QSP model, Bone)

Full references: [`ckdmbd_references.md`](ckdmbd_references.md)
