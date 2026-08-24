# Vascular Dementia (VaD) — Quantitative Systems Pharmacology Model

> **Category**: Neurological / Cerebrovascular disease
> **Abbreviation**: VaD
> **Build date**: 2026-06-27
> **Model version**: 1.0.0

---

## Disease Overview

**Vascular dementia (VaD)** is a cognitive impairment arising from brain injury caused by cerebrovascular disease, and is the **second most common cause of dementia after Alzheimer's disease** (roughly 20–30% of all dementia). It encompasses several pathological subtypes, including white-matter hyperintensities (WMH) and lacunar infarcts from small vessel disease (SVD), large-vessel atherosclerosis, and cortical infarcts from cardioembolism.

### Key Pathophysiological Pathways

```
Hypertension/diabetes/dyslipidaemia → small vessel disease (SVD)
     ↓
White-matter hyperintensities (WMH) + lacunar infarct + cortical microinfarcts
     ↓
Blood-brain barrier (BBB) disruption → DAMP influx
     ↓
Microglial M1 activation → NF-κB → IL-1β/TNF-α/IL-6
     ↓
Oxidative stress (ROS) → neuronal injury
     ↓
Cholinergic hypofunction + synaptic loss
     ↓
Cognitive decline (MMSE↓, executive function↓, processing speed↓)
```

---

## Model Structure

### 1. Mechanistic Map

| Item | Detail |
|------|------|
| Node count | **114** |
| Subgraph clusters | **9** |
| Cluster list | Drug PK/PD · vascular risk factors · cerebrovascular pathology · blood-brain barrier · neuroinflammation · oxidative stress · neurotransmitters · brain structure · clinical outcomes |

[![Mechanistic map preview](vad_qsp_model.png)](vad_qsp_model.svg)

### 2. mrgsolve ODE Model

| Item | Detail |
|------|------|
| ODE compartments | **18** |
| Drug PK compartments | 7 (AHT depot/central, APT, statin, AChEI brain, memantine brain, cilostazol) |
| Physiology/PD compartments | 11 (BP, LDL, CBF, WMH, infarct, microglia, cytokine, ROS, ACh, synapse, MMSE) |
| Treatment scenarios | **6** |
| Clinical trial calibration | PROGRESS, SCOPE, SPRINT-MIND, SIGNAL2, Black et al. 2003 (donepezil) |

**6 Treatment Scenarios:**

| # | Scenario | Drug |
|---|---------|------|
| 1 | No treatment | — |
| 2 | Antihypertensive monotherapy | Antihypertensive (ARB/ACEi) |
| 3 | Combined vascular risk-factor therapy | Antihypertensive + antiplatelet + statin |
| 4 | Symptomatic treatment | AChEI (donepezil) + memantine |
| 5 | Comprehensive treatment | Vascular + symptomatic + cilostazol |
| 6 | Optimal+ (high-dose statin) | Scenario 5 + statin 80 mg/day |

### 3. Shiny Dashboard

| Tab | Content |
|----|------|
| 1. Patient profile | Baseline settings for MMSE, WMH, BP, CBF, LDL, ACh, synaptic density, and risk-factor summary |
| 2. Drug PK/PD | Plasma concentration-time profile (48 h), steady-state PD effect summary |
| 3. Vascular/cerebral perfusion | SBP, CBF, WMH progression, cumulative microinfarct simulation |
| 4. Neurobiological mechanisms | Microglia/cytokines, ROS, ACh tone, synaptic density |
| 5. Clinical endpoints | MMSE trajectory, cognitive domains, VaD stage classification |
| 6. Scenario comparison | 2-year comparison of MMSE/WMH/CBF across the 6 treatment strategies |
| 7. Biomarkers | MRI/CSF/blood biomarker panel (normal range vs VaD range) |

### 4. References

| Item | Detail |
|------|------|
| Total references | **70** |
| Categories | Diagnosis/epidemiology · small vessel disease · blood-brain barrier · neuroinflammation · oxidative stress · cholinergic · clinical trials (BP/statin/antiplatelet/cilostazol/donepezil/memantine) · QSP methodology · biomarkers · prevention/risk factors |

---

## Key Parameter Calibration

| Parameter | Value | Basis |
|---------|-----|------|
| Emax_AHT_BP | 25 mmHg | PROGRESS trial (perindopril-based Rx) |
| WMH progression rate | ~0.65 mL/yr (untreated) | PROGRESS MRI substudy, Dufouil 2005 |
| AChEI effect | +1.0 MMSE @ 24 weeks | Black et al. 2003 (donepezil in VaD) |
| Cilostazol WMH reduction | −11.5% | SIGNAL2 trial |
| E_ST_LDL (statin) | 50% LDL reduction | Multiple statin RCTs |
| Emax_AChEI | 75% AChE inhibition | Pharmacodynamic modelling |
| CBF baseline (VaD) | ~45-55 mL/100g/min | ASL-MRI, xenon-CT data |

---

## Files

| File | Description |
|------|------|
| `vad_qsp_model.dot` | Graphviz mechanistic map source |
| `vad_qsp_model.svg` | SVG vector image |
| `vad_qsp_model.png` | PNG raster image (150 dpi) |
| `vad_mrgsolve_model.R` | mrgsolve ODE model (18 compartments, 6 scenarios) |
| `vad_shiny_app.R` | Shiny interactive dashboard (7 tabs) |
| `vad_references.md` | References (70) |
| `README.md` | This file |

---

## How to Run

### Render the Mechanistic Map (requires Graphviz)
```bash
dot -Tsvg vad_qsp_model.dot -o vad_qsp_model.svg
dot -Tpng -Gdpi=150 vad_qsp_model.dot -o vad_qsp_model.png
```

### Run the mrgsolve Model (requires R)
```r
install.packages(c("mrgsolve", "dplyr", "ggplot2", "tidyr"))
source("vad_mrgsolve_model.R")
```

### Run the Shiny App
```r
install.packages(c("shiny", "shinydashboard", "ggplot2", "dplyr", "DT", "tidyr"))
shiny::runApp("vad_shiny_app.R")
```

---

## Drug Mechanism Summary

| Drug | Mechanism of action | Key effect |
|------|---------|----------|
| **Antihypertensive** (ARB/ACEi) | RAAS inhibition → SVD prevention, BBB protection | Suppresses WMH progression, improves cerebral perfusion |
| **Statin** | HMG-CoA inhibition + pleiotropic effects (eNOS↑, NF-κB↓, BDNF↑) | LDL↓, neuroinflammation↓, improved vascular function |
| **Antiplatelet** | TXA2 synthesis inhibition (COX-1) | Platelet aggregation↓, small-vessel thrombosis↓ |
| **AChEI (donepezil)** | AChE inhibition → synaptic ACh↑ | Restores cholinergic tone, stabilises MMSE |
| **Memantine** | Low-affinity NMDA-R blockade | Excitotoxicity↓, suppresses synaptic loss |
| **Cilostazol** | PDE3 inhibition → cAMP↑ → vasodilation | CBF↑, WMH progression↓, platelet aggregation↓ |

---

## Limitations

- This model is a semi-quantitative QSP model intended for **educational and research purposes**.
- It cannot be applied directly in clinical practice without individual calibration against patient data.
- Mixed-type VaD (coexisting with Alzheimer's disease) requires a separate model extension.
- Parameters are representative (mean) values drawn from published literature and may not fully capture individual variability.

---

## Related Models in This Repository

- [Alzheimer's Disease](../alzheimers-disease/) — model centred on the Aβ/tau cascade
- [Essential Hypertension](../essential-hypertension/) — RAAS/SNS-centred BP regulation model
- [Ischemic Stroke](../ischemic-stroke/) — acute ischaemia and reperfusion injury model
- [Diabetic Nephropathy](../diabetic-nephropathy/) — metabolic-vascular mechanism model
- [Multiple Sclerosis](../multiple-sclerosis/) — CNS neuroinflammation model
