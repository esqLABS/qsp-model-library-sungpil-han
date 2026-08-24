# Recurrent Pericarditis (RP) — QSP Model

[![RP Mechanistic Map](rp_qsp_model.png)](rp_qsp_model.svg)

**Classification:** Cardiovascular · Autoinflammatory
**Directory:** `recurrent-pericarditis/`

---

## Overview

Recurrent pericarditis (RP) is an inflammatory pericardial disease that relapses after a symptom-free interval (≥4-6 weeks) following an initial episode of acute pericarditis; without treatment, it recurs in **30-50%** of patients. The central pathomechanism is the **NLRP3 inflammasome → IL-1β/IL-18 axis**, making it an excellent target for IL-1 blockers, and colchicine, NSAIDs, corticosteroids, anakinra (IL-1Ra), and rilonacept (IL-1 trap) are the principal treatments. This QSP model integrates all of the major pharmacological pathways to simulate disease dynamics and treatment response.

---

## File List

| File | Description |
|------|------|
| `rp_qsp_model.dot` | Graphviz mechanistic map (14 clusters, 162+ nodes) |
| `rp_qsp_model.svg` | Vector graphic (for interactive browser viewing) |
| `rp_qsp_model.png` | PNG thumbnail (150 dpi) |
| `rp_mrgsolve_model.R` | mrgsolve ODE PK/PD model (24 compartments, 7 scenarios) |
| `rp_shiny_app.R` | Shiny dashboard (7 tabs) |
| `rp_references.md` | 57 references |

---

## Pathophysiology

| Key pathological pathway | Mechanistic detail |
|--------------|-------------|
| **NLRP3 inflammasome** | Signal 1 (TLR→NF-κB→NLRP3/pro-IL-1β transcription) + Signal 2 (ATP, crystals, K⁺ efflux, ROS) → ASC speck formation → Casp-1 → IL-1β/IL-18 cleavage and secretion |
| **IL-1β amplification loop** | IL-1β → IL-1R1 → IRAK4 → NF-κB → TNF/IL-6/CXCL8 → further NLRP3 activation (positive feedback) |
| **Neutrophil infiltration** | CXCL8/LTB4 → pericardial neutrophil influx → NET formation → DAMP release → inflammasome reactivation |
| **Fibrotic pathway** | M2 macrophages → TGF-β1 → pericardial myofibroblasts → collagen I/III → pericardial thickening → constrictive pericarditis |
| **Adaptive immunity** | Memory T cells → reactivation (triggers relapse); B cells → anti-heart antibodies → immune complexes → complement activation |
| **Eicosanoids** | AA → COX-2 → PGE2 (fever, pain, increased vascular permeability) |

---

## Model Specifications

### Mechanistic Map
- **14 subgraph clusters**: etiology · DAMP-PAMP recognition · NLRP3 inflammasome · NF-κB · innate immune cells · adaptive immunity · eicosanoids · pericardial pathology · clinical biomarkers · colchicine PK/PD · NSAID PK/PD · corticosteroid PK/PD · biologic PK/PD · clinical outcomes
- **162+ nodes**: covers all major molecular, cellular, and clinical components

### mrgsolve ODE Model (24 Compartments)

| Compartment group | # compartments | Content |
|-----------|---------|------|
| Colchicine PK | 3 | Gut · central · peripheral |
| Ibuprofen PK | 2 | Gut · central |
| Prednisone PK | 2 | Gut · central |
| Anakinra PK | 2 | SC depot · central |
| Rilonacept PK | 2 | SC depot · central |
| Immune/inflammation | 7 | NLRP3_ACT · IL1B · IL18 · TNF · IL6 · NEUTRO · M1_MACRO |
| Pericardial pathology | 4 | INFLAM · EFFUSION · FIBRIN · FIBROSIS |
| Clinical biomarkers | 2 | CRP · PAIN |

### Treatment Scenarios (7)

| # | Scenario | Evidence base |
|---|---------|------------|
| 1 | No treatment (natural course) | — |
| 2 | Ibuprofen 600 mg TID × 4 weeks | ESC 2015 |
| 3 | Colchicine 0.5 mg BID × 3 months | COPE, ICAP |
| 4 | Colchicine + ibuprofen (combination) | COPE, ICAP |
| 5 | Prednisone 0.5 mg/kg/d → gradual taper | ESC 2015 |
| 6 | Anakinra 100 mg/d SC × 6 months | AIRTRIP 2016 |
| 7 | Rilonacept 320→160 mg SC qw | RHAPSODY 2021 |

### Clinical Trial Calibration Data

| Trial | Drug | Control-arm recurrence rate | Treatment-arm recurrence rate | RRR |
|---------|------|------------|------------|-----|
| COPE (2005) | Colchicine + ASA | 45% | 24% | 47% |
| ICAP (2013) | Colchicine 0.5 mg BID | 32.3% | 16.7% | 48% |
| CORP (2011) | Colchicine (2nd episode) | 45.5% | 19.2% | 58% |
| AIRTRIP (2016) | Anakinra 100 mg/d | 90.9% | 18.2% | 80% |
| RHAPSODY (2021) | Rilonacept 320→160 mg qw | 74.4% | 8.8% | HR 0.04 |

---

## Shiny App Tab Structure (7 Tabs)

| Tab | Content |
|----|------|
| ① Patient & scenario | Patient profile, treatment selection, ESC diagnostic criteria, overview plot |
| ② Drug PK | Concentration-time curves, PD inhibition, PK parameter table |
| ③ Inflammasome/cytokines | NLRP3 · IL-1β · IL-18 · TNF · IL-6 · immune-cell dynamics |
| ④ Pericardial pathology | Inflammation · effusion · fibrin · fibrosis · risk stratification table |
| ⑤ Clinical endpoints | Pain VAS · CRP · summary table · treatment goals |
| ⑥ Scenario comparison | Inflammation/CRP/IL-1β/effusion comparison across all 7 scenarios + clinical-trial benchmarks |
| ⑦ Biomarkers | Time-point snapshots, trajectories, Emax curves, risk radar chart |

---

## Drug Mechanism of Action Summary

```
Colchicine:   binds β-tubulin → microtubule disassembly
              → inhibits NLRP3 ASC assembly (IC50 ≈ 0.5 ng/mL)
              → inhibits neutrophil migration, inhibits L-selectin shedding

Ibuprofen:    COX-1/2 inhibition → reduced PGH2 → PGE2↓
              → reduced fever, pain, and vascular permeability

Prednisone:   GRα → transcriptional repression of NF-κB → IL-1β/TNF/IL-6↓
              → but risk of rebound relapse with rapid tapering

Anakinra:     competitive IL-1R1 blockade (IC50 ≈ 0.1 nM)
              → blocks downstream IL-1β signalling

Rilonacept:   dual IL-1α/β trap (KD < 1 pM)
              → neutralises circulating IL-1β → blocks the NLRP3 positive loop
```

---

## References

57 total — see [rp_references.md](../../../recurrent-pericarditis/rp_references.md)

Key references:
- Imazio M, et al. NEJM 2016 (AIRTRIP): https://pubmed.ncbi.nlm.nih.gov/27668557/
- Klein AL, et al. NEJM 2021 (RHAPSODY): https://pubmed.ncbi.nlm.nih.gov/33405895/
- Imazio M, et al. NEJM 2013 (ICAP): https://pubmed.ncbi.nlm.nih.gov/24131175/
- Adler Y, et al. Eur Heart J 2015 (ESC Guidelines): https://pubmed.ncbi.nlm.nih.gov/26320112/

---

## Disclaimer

This model is a QSP model intended for educational and research purposes, and must not be used directly for clinical decision-making.
