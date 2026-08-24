# ALS (Amyotrophic Lateral Sclerosis) — QSP Model

[![ALS Mechanistic Map](../../../amyotrophic-lateral-sclerosis/als_qsp_model.png)](../../../amyotrophic-lateral-sclerosis/als_qsp_model.svg)

---

## Overview

**Amyotrophic lateral sclerosis (ALS)** is a fatal progressive neurodegenerative disease affecting both upper and lower motor neurons. Median survival from symptom onset is 2–5 years, with respiratory failure as the leading cause of death (~65%). Global prevalence is ~4–6 per 100,000. ~90% are sporadic; ~10% familial (SOD1 2%, C9orf72 5–10%, FUS/TDP-43 others).

---

## Model Components

### Mechanistic Map

**196+ nodes · 10 subgraph clusters**

| Cluster | Components | Content |
|---------|------------|------|
| 1. Genetic & Molecular | 18 | SOD1, TDP-43, FUS, C9orf72, UBQLN2, VCP, NEK1, etc. |
| 2. Protein Quality Control | 18 | 26S proteasome, autophagy, ER stress (PERK/IRE1/ATF6) |
| 3. Oxidative Stress & Mitochondria | 20 | ROS, GSH, Nrf2, mitochondrial function, apoptosis |
| 4. Glutamate Excitotoxicity | 18 | EAAT2, AMPA/NMDA receptors, Ca²⁺, nNOS, ONOO⁻ |
| 5. Motor Neuron Survival | 23 | BDNF/TrkB, IGF-1, GDNF, PI3K/Akt/mTOR, ERK1/2 |
| 6. Neuroinflammation | 25 | Microglial M1/M2, TNF-α, IL-1β, NLRP3, Treg |
| 7. Axonal Transport | 18 | NF-L/M/H, kinesin, dynein, mitochondrial transport |
| 8. NMJ & Motor Unit | 18 | ACh, nAChR, agrin/MuSK, EMG, MUNIX |
| 9. Drug PK/PD | 21 | Riluzole, edaravone, tofersen, AMX0035, masitinib |
| 10. Clinical Endpoints | 17 | ALSFRS-R, FVC%, NfL, pNF-H, survival |

### mrgsolve ODE Model (ALS ODE Model)

**26 ODEs · 60+ parameters · 7 treatment scenarios**

```r
# List of key compartments
# PK : DEPOT_RIL, C1_RIL, C2_RIL, IV_EDA, DEPOT_TOF, C1_TOF, C2_TOF,
#      DEPOT_PB, C_PB  (9 compartments)
# Disease: MN_upper, MN_lower, SOD1_wt, SOD1_mis, TDP43_nuc, TDP43_cyto,
#          Glut_syn, Ca_i, ROS, GSH, Mito, Mic_act, TNFa, BDNF  (14 compartments)
# Clinical: NfL_CSF, ALSFRS, FVC  (3 compartments)
```

#### Key Pathogenesis Equation

```
# Motor neuron death rate (multifactorial model)
MN_death_rate = k_MN_death × (1 + 3·SOD1_burden + 2·TDP43_cyto_frac
                              + 1.5·ROS_norm + 1.2·Glu_excess + 0.8·TNFa
                              - BDNF_prot_wt × (BDNF/BDNF₀))
```

### Shiny Dashboard Tabs (6 Tabs)

| Tab | Content |
|---|------|
| Patient Profile | Patient parameters, disease overview, PK summary table, clinical milestones |
| Drug PK | Riluzole/edaravone/tofersen concentration-time curves and drug effects |
| Biomarkers | NfL, TNF-α, ROS/GSH, mitochondrial function, BDNF |
| Clinical Endpoints | ALSFRS-R, FVC%, motor-neuron survival %, glutamate toxicity |
| Scenario Comparison | Simultaneous comparison of the 7 treatment arms (ALSFRS-R · MN survival · NfL) |
| Mechanistic Pathways | TDP-43 nuclear/cytoplasmic dynamics, SOD1 protein misfolding |

---

## Drug PK/PD Parameters

| Drug | Regimen | F(%) | t½ | CL | Vd | Mechanism | IC50/EC50 |
|-----|-----|------|----|----|----|---------|----|
| **Riluzole** | 50 mg PO BID | 60 | 12 h | 28 L/h | 245 L | ↓ Glu release via Na⁺ channel | IC50=0.5 μg/mL |
| **Edaravone** | 60 mg IV/day (cycle) | ~100 | 4.5 h | 18 L/h | 120 L | Free-radical scavenger | IC50=1.2 μg/mL |
| **Tofersen** | 100 mg SC q4w | — | 7 days | 0.5 L/h | 20 L | SOD1 mRNA knockdown (ASO) | EC50=0.1 μg/mL CSF |
| **AMX0035** | PB 3g + TUDCA 1g PO BID | 85 | 3 h | 12 L/h | 50 L | ↓ ER stress / mitochondrial protection | EC50=50 μmol/L |
| **Masitinib** | 4.5 mg/kg/day PO | 58 | 40 h | 25 L/h | 190 L | c-Kit/PDGFR-β → ↓ microglia | Phase 3 ongoing |

---

## 7 Treatment Scenarios

| Scenario | Treatment | Expected ALSFRS-R slowdown |
|---------|------|-----------------|
| 1. Untreated | None | Baseline |
| 2. Riluzole | 50 mg BID | ~12% improvement |
| 3. Edaravone | 60 mg/day (cycle) | ~10% improvement |
| 4. Riluzole + Edaravone | Combination | ~20% improvement |
| 5. Tofersen | SC q4w (SOD1-ALS) | ~30-40% improvement (SOD1 subtype) |
| 6. AMX0035 | PB+TUDCA BID | ~25% improvement |
| 7. All Drugs | All combined | ~45-50% improvement |

---

## Pathophysiology Summary

```
Genetic predisposition (SOD1/TDP-43/FUS/C9orf72)
    ↓
Collapse of protein homeostasis → formation of toxic aggregates
    ↓
① Glutamate toxicity: ↓EAAT2 → ↑synaptic Glu → AMPA/NMDA hyperactivation → ↑Ca²⁺
② Oxidative stress: ↓SOD1 function → ↑ROS → ↓GSH → mitochondrial damage
③ Neuroinflammation: microglial M1 skewing → ↑TNF-α/IL-1β → neurotoxicity
④ Trophic factor deficiency: ↓BDNF / ↓GDNF → reduced survival signalling
⑤ Axonal transport impairment: NF accumulation → blocked mitochondrial transport
    ↓
Motor neuron loss → denervation atrophy → ALSFRS-R decline → respiratory failure → death
```

---

## File List

| File | Description |
|------|------|
| [als_qsp_model.dot](../../../amyotrophic-lateral-sclerosis/als_qsp_model.dot) | Graphviz mechanistic map source (196+ nodes, 10 clusters) |
| [als_qsp_model.svg](../../../amyotrophic-lateral-sclerosis/als_qsp_model.svg) | SVG vector image (zoomable) |
| [als_qsp_model.png](../../../amyotrophic-lateral-sclerosis/als_qsp_model.png) | PNG raster image (150 dpi) |
| [als_mrgsolve_model.R](../../../amyotrophic-lateral-sclerosis/als_mrgsolve_model.R) | mrgsolve ODE model (26 compartments, 60+ parameters, 7 scenarios) |
| [als_shiny_app.R](../../../amyotrophic-lateral-sclerosis/als_shiny_app.R) | Shiny dashboard (6 tabs, plotly, DT) |
| [als_references.md](../../../amyotrophic-lateral-sclerosis/als_references.md) | 40 references (with PubMed links) |

---

## Key Clinical Trials

| Trial | Drug | Result |
|-----|------|------|
| Bensimon 1994 (N Engl J Med) | Riluzole | Survival +3 months |
| Lacomblez 1996 (Lancet) | Riluzole | Confirmation + dose optimisation |
| Abe 2017 (Lancet Neurol) | Edaravone | ALSFRS-R decline improved by 33% |
| ATLAS/Miller 2022 (NEJM) | Tofersen | SOD1 protein ↓50%, NfL ↓60% |
| CENTAUR/Paganoni 2020 (NEJM) | AMX0035 | ALSFRS-R +2.3 points/yr (45 weeks) |
| AB Science Phase 3 | Masitinib | Phase 3 ALSFRS-R primary endpoint ongoing |

---

*Model written: 2026-06-21 | ALS QSP Library v1.0*
