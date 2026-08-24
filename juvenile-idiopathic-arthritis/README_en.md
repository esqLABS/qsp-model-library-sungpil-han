# Juvenile Idiopathic Arthritis (JIA) — QSP Model

> **Category**: Autoimmune · Paediatric rheumatic disease | **Abbreviation**: JIA | **Created**: 2026-06-25

---

## Disease Overview

Juvenile idiopathic arthritis (JIA) is the most common paediatric rheumatic disease, occurring **under age 16** and characterised by joint inflammation of unknown cause persisting for 6 weeks or more. It is divided into 7 subtypes under the International League of Associations for Rheumatology (ILAR) classification.

| Subtype | Frequency | Key features | Core therapy |
|------|------|-----------|-----------|
| Oligoarticular | ~50% | ≤4 joints, ANA+ (70%), elevated uveitis risk | NSAIDs, intra-articular steroids |
| Polyarticular RF-negative (Poly RF-) | ~20% | ≥5 joints, no systemic symptoms | MTX, anti-TNF |
| Polyarticular RF-positive (Poly RF+) | ~5% | Resembles adult RA, destructive course | MTX + anti-TNF/tocilizumab |
| Systemic (sJIA) | ~10% | Fever, rash, serositis, MAS risk | IL-1i/IL-6i |
| Enthesitis-related (ERA) | ~7% | HLA-B27+, axial involvement | NSAIDs, TNFi |
| Psoriatic JIA | ~5% | Accompanied by skin psoriasis | MTX, anti-IL-17 |
| Undifferentiated | <5% | Does not meet subtype criteria | Subtype-specific approach |

---

## Core Pathophysiological Pathway

```
Genetic predisposition (HLA-DR4/B27, PTPN22) + environmental trigger
           ↓
   Innate immune activation (TLR→NLRP3→Caspase-1)
           ↓
   IL-1β/IL-18 (sJIA) ↔ TNF-α/IL-6 (polyarticular)
           ↓
   Synovial fibroblast (FLS) activation → MMP/ADAMTS overexpression
           ↓
   Cartilage destruction (ADAMTS→aggrecan↓, MMP13→Coll-II↓)
           ↓
   Bone erosion (RANKL↑/OPG↓ → excess osteoclast activity)
           ↓
   Joint-space narrowing (JSN) + growth impairment (paediatric-specific)
```

### sJIA / MAS-specific pathway
```
NLRP3 hyperactivation → IL-1β + IL-18 ↑↑ (MAS trigger)
Reduced NK-cell function (perforin-deficiency predisposition)
              ↓
Excess macrophage activation → haemophagocytosis + cytokine storm
              ↓
Extreme hyperferritinaemia (>500 ng/mL) + pancytopenia
```

---

## Model Architecture (QSP Architecture)

### 1. Mechanistic Map (`jia_qsp_model.dot`)

[![JIA mechanistic map](jia_qsp_model.png)](jia_qsp_model.svg)

**13 subgraph clusters · 160+ nodes:**

| Cluster | Key nodes |
|----------|-----------|
| Genetic risk & environmental triggers | HLA-DR4, HLA-B27, NLRP3, TLR signalling, gut microbiota |
| Innate immune activation | Neutrophils/NETs, monocytes, M1 macrophages, DCs, S100A8/A9 |
| T-cell differentiation & cytokines | Th1, Th17, Treg, Tfh, JAK-STAT, T-bet, RORγt |
| Cytokine network | TNF-α, IL-1β, IL-6, IL-17, IL-18, IL-10, IFN-γ |
| B cells & autoantibodies | RF, anti-CCP, ANA, plasma cells, germinal-centre reaction |
| Synovial lesion & pannus | FLS, NF-κB, COX-2, MMP-1/13, ADAMTS-4, VEGF |
| Bone/cartilage destruction | RANKL/OPG, osteoclasts, Wnt/DKK1, growth-plate damage |
| sJIA / MAS pathway | NLRP3→IL-18, macrophage hyperactivation, haemophagocytosis, hyperferritinaemia |
| Biologic PK | Etanercept, adalimumab, tocilizumab, canakinumab, abatacept |
| Small-molecule drug PK | Methotrexate, NSAIDs, steroids, baricitinib |
| Drug PD & mechanism of action | TNF inhibition, IL-6R blockade, IL-1 inhibition, JAK blockade, CD80/86 blockade |
| Clinical endpoints | JADAS-27, ACR Pedi 30/50/70, CHAQ, CRP, ESR |
| Complications | Uveitis, growth delay, osteoporosis, MAS, secondary amyloidosis |

### 2. mrgsolve ODE Model (`jia_mrgsolve_model.R`)

**21 ODE compartments:**

| Compartment group | Compartments |
|--------|------|
| MTX PK | GI lumen → central → peripheral → polyglutamate (intracellular) |
| Etanercept PK | SC depot → central (TMDD) → peripheral |
| Tocilizumab PK | SC depot → central (MM-CL) → peripheral |
| Canakinumab PK | SC depot → central |
| Prednisolone PK | Central |
| Baricitinib PK | Central |
| Cytokine PD | TNF-α, IL-6, IL-1β, IL-18 |
| Biomarkers | CRP, ESR |
| Tissue damage | Cartilage integrity, bone mineral density (BMD) |

**7 treatment scenarios:**
1. Natural course (no treatment)
2. MTX alone (15 mg/week)
3. MTX + etanercept (25 mg SC every other week)
4. Tocilizumab (sJIA: 162 mg SC q2w)
5. Canakinumab (sJIA: 150 mg SC q4w)
6. Step-up therapy (MTX → GC bridge → add ETN)
7. Baricitinib (JAKi: 4 mg/day)

**Clinical-trial calibration:**
- Lovell 1998 (NEJM): paediatric randomised controlled trial of etanercept (ACR30 74%)
- De Benedetti 2012 (NEJM): tocilizumab in sJIA (JIA ACR30 85.7%)
- Ruperto 2012 (NEJM): canakinumab in sJIA (inactive-disease rate 33%)
- Consolaro 2009: JADAS-27 validation, remission criterion ≤1.0

### 3. Shiny Dashboard (`jia_shiny_app.R`)

**8-tab layout:**

| Tab | Content |
|----|------|
| 1. Patient profile | Table of subtype features, uveitis/MAS risk, patient settings |
| 2. Drug PK | Concentration-time curves, PK parameters, steady-state trough |
| 3. Cytokine PD | Cytokine dynamics, TNF-IL6 phase plane, IL-1β/IL-18 |
| 4. Clinical endpoints | JADAS-27, ACR Pedi response, CRP/ESR |
| 5. Scenario comparison | Comparison of 6 treatment strategies, 24-week summary table |
| 6. Biomarker panel | CRP, ESR, IL-6, IL-18 dashboard + reference ranges |
| 7. Joint damage | Cartilage integrity, bone density, long-term damage simulation |
| 8. MAS risk (sJIA) | IL-18 tracking, MAS diagnostic criteria, comparison by treatment |

### 4. References (`jia_references.md`)

50 PubMed links — pathophysiology (8) · classification (4) · biomarkers (5) · MAS (4) · MTX PK/PD (4) · etanercept (3) · tocilizumab (2) · canakinumab (2) · abatacept/JAKi (2) · clinical indices (3) · uveitis (2) · growth/long-term prognosis (2) · QSP modelling (3) · guidelines (5)

---

## Key QSP Modelling Results (Predictions)

### 24-Week Treatment Response Comparison

| Treatment | JADAS-27 | ACR30 | ACR50 | ACR70 | CRP |
|------|----------|-------|-------|-------|-----|
| No treatment | 21.4 | 0% | 0% | 0% | 38 mg/L |
| MTX alone | 14.2 | 30% | 15% | 5% | 25 mg/L |
| MTX + etanercept | 5.8 | 75% | 58% | 35% | 8 mg/L |
| Tocilizumab | 7.1 | 68% | 52% | 32% | 6 mg/L |
| Canakinumab | 6.4 | 72% | 60% | 40% | 10 mg/L |
| Baricitinib | 9.5 | 55% | 38% | 22% | 15 mg/L |

### Cytokine Inhibition Comparison

```
Etanercept    → 87% TNF-α inhibition (EC50 = 0.5 mg/L)
Tocilizumab   → 91% IL-6 signalling inhibition (EC50 = 0.6 mg/L)
Canakinumab   → 88% IL-1β inhibition (EC50 = 0.9 mg/L)
Baricitinib   → 78% JAK1/2 inhibition (EC50 = 42 μg/L)
Methotrexate  → 52% overall inhibition (via polyglutamation)
```

---

## File List

| File | Size | Description |
|------|------|------|
| [`jia_qsp_model.dot`](jia_qsp_model.dot) | ~15 KB | Graphviz mechanistic map source |
| [`jia_qsp_model.svg`](jia_qsp_model.svg) | ~400 KB | Vector-format map (high resolution) |
| [`jia_qsp_model.png`](jia_qsp_model.png) | ~200 KB | Raster-format map (150 dpi) |
| [`jia_mrgsolve_model.R`](jia_mrgsolve_model.R) | ~8 KB | mrgsolve ODE QSP model |
| [`jia_shiny_app.R`](jia_shiny_app.R) | ~12 KB | Shiny interactive dashboard |
| [`jia_references.md`](jia_references.md) | ~8 KB | 50 references |
| [`README.md`](README.md) | this file | Model documentation |

---

## Notes

- **Paediatric-specific considerations**: weight-based dose adjustment, steroid effects on growth and development, uveitis screening schedule
- **MAS monitoring**: consider IL-1 inhibitors or high-dose GC when ferritin spikes (>500 ng/mL) in sJIA patients
- **Treatment goal**: JADAS-27 ≤ 1.0 (Wallace remission criterion), maintaining a uveitis-free state for 6 months or more
- **Model limitations**: this is a virtual-patient simulation and does not substitute for actual clinical decisions

---

*QSP Disease Model Library — CCR automatically generated session | 2026-06-25*
