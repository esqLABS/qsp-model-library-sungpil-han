# Migraine QSP Model

## Overview

Migraine is one of the most common neurological disorders, affecting approximately 15% of the global population (roughly 1 billion people).
It ranks second in the Global Burden of Disease (GBD) analysis as a cause of years lived with disability (YLD), with the highest prevalence seen particularly in women aged 20–50.

| Item | Value |
|------|------|
| Global prevalence | ~15% (women 20%, men 10%) |
| Chronic migraine (≥15 MHD) | ~8% of prevalence |
| Mean attack duration | 4–72 hours |
| Annual direct cost (US) | >$3.6 billion |
| DALY (global) | 45.1 million DALY (GBD 2016) |

---

## Key Pathophysiological Pathways

### 1. Trigeminovascular Pathway

| Component | Role |
|----------|------|
| Trigeminal ganglion (TG) | CGRP/SP/PACAP synthesis and release |
| Trigeminal nucleus caudalis (TNC) | Second-order neuron, initiates central sensitisation |
| Trigeminocervical complex (TCC) | C1–C2 convergence, neck/occipital pain |
| Thalamus | Third-order neuron, pain perception |
| Cerebral cortex | Processing of photophobia/phonophobia |

### 2. CGRP Signaling

During a migraine attack, excessive CGRP is released from trigeminal nerve terminals, causing:
- Activation of **CLR/RAMP1 receptors** → ↑cAMP → opening of KATP channels → **dural vasodilation**
- Degranulation of dural mast cells → release of histamine/bradykinin → peripheral trigeminal sensitisation
- The therapeutic target of CGRP antibodies (erenumab, fremanezumab, galcanezumab) and gepants (rimegepant, ubrogepant)

### 3. Cortical Spreading Depression (CSD)

CSD is the neurophysiological mechanism of aura (visual scintillating scotoma) and triggers headache attacks via peripheral trigeminal sensitisation:
- K⁺ efflux (30–80 mM) → NMDA receptor activation → Ca²⁺ influx → neuronal depolarisation
- Propagation velocity: 3–5 mm/min (occipital → frontal)
- **CSD → Panx1 channel activation → dural trigeminal activation → CGRP release**

### 4. Central Sensitization

More than 10 minutes after attack onset, TNC second-order neurons become hyperactive:
- Cutaneous allodynia — hypersensitivity to touch and temperature
- PKCε → MAPK/ERK → BDNF upregulation → progression to chronic migraine
- **Triptans: most effective when taken before central sensitisation sets in (within 10 minutes)**

---

## Drug PK/PD Parameters

### Acute Treatments

| Drug | Mechanism of action | Key PK | 2h pain-free rate |
|------|---------|---------|--------------|
| Sumatriptan SC 6 mg | 5-HT1B/1D agonist | F=97%, t½=2h, CL=72L/h | ~35–40% |
| Sumatriptan oral 100 mg | 5-HT1B/1D agonist | F=14%, t½=2h | ~25–30% |
| Lasmiditan 200 mg | 5-HT1F agonist (no vasoconstriction) | F=38%, t½=5h | ~32–39% |
| Rimegepant 75 mg | CGRP-R antagonist (gepant) | F=64%, t½=11h | ~21% |
| Ubrogepant 100 mg | CGRP-R antagonist (gepant) | F=44%, t½=7h | ~19% |

### Preventive Treatments

| Drug | Mechanism of action | Administration | MMD reduction | 50% responder rate |
|------|---------|------|----------|-----------|
| Erenumab 140 mg | Anti-CGRP-receptor mAb | SC once monthly | 3.7 days | 47% |
| Fremanezumab 225 mg | Anti-CGRP-ligand mAb | SC once monthly | 3.7 days | 43% |
| Galcanezumab 120 mg | Anti-CGRP-ligand mAb | SC once monthly | 4.7 days | 52% |
| Rimegepant 75 mg QOD | CGRP-R antagonist | Oral, alternate days | 1.75 days | 28% |
| Topiramate 100 mg | AMPA/CA inhibition | Oral, divided | 2.1 days | 37% |
| Propranolol 160 mg | β1/β2-blocker | Oral, divided | 1.8 days | 35% |
| Amitriptyline 50 mg | SERT/NET inhibition (TCA) | Oral, nightly | 1.6 days | 30% |

---

## Model Structure

### mrgsolve ODE compartments (18 state variables)

```
Drug PK compartments (10):
  DEPOT_SUM, CENT_SUM, PERI_SUM     — sumatriptan SC 2-compartment
  DEPOT_ERE, CENT_ERE, PERI_ERE     — erenumab SC 2-compartment
  CENT_RIM                           — rimegepant 1-compartment
  DEPOT_TOP, CENT_TOP                — topiramate 2-compartment
  CGRPR_FREE                         — free CGRP receptor

Disease PD compartments (8):
  CGRP_TG     — trigeminal ganglion/plasma CGRP (pmol/L)
  CSD_ACT     — CSD activity (0–1)
  TG_ACT      — trigeminal activation (0–1)
  CS_STATE    — central sensitisation state (0–1)
  PGE2_COMP   — tissue PGE2 (pg/mL)
  NO_COMP     — nitric oxide level (pmol/L)
  SEROTONIN   — platelet 5-HT (ng/mL)
  PAIN_SCORE  — VAS pain score (0–10)
```

### Treatment Scenarios (7)

1. **Untreated acute attack** — 24h pain tracking after CSD induction
2. **Sumatriptan SC 6 mg** — 5-HT1B/1D agonism → CGRP release inhibition + vasoconstriction
3. **Lasmiditan 200 mg** — selective 5-HT1F agonism (excellent CNS penetration)
4. **Rimegepant 75 mg** — dual acute + preventive effect (CGRP-R antagonism)
5. **Erenumab 140 mg once monthly × 3 months** — CGRP-R occupancy and MMD reduction
6. **Topiramate 100 mg/day SS** — CSD suppression + AMPA blockade
7. **Chronic migraine vs erenumab, 1 year** — long-term progression modelling

---

## File List

| File | Description |
|------|------|
| [mgr_qsp_model.dot](../../../migraine/mgr_qsp_model.dot) | Graphviz mechanistic map source (100+ nodes, 12 subgraph clusters) |
| [mgr_qsp_model.svg](../../../migraine/mgr_qsp_model.svg) | SVG vector image (scalable) |
| [mgr_qsp_model.png](../../../migraine/mgr_qsp_model.png) | PNG raster image (150 dpi) |
| [mgr_mrgsolve_model.R](../../../migraine/mgr_mrgsolve_model.R) | mrgsolve ODE model (18 compartments, 7 scenarios) |
| [mgr_shiny_app.R](../../../migraine/mgr_shiny_app.R) | Shiny dashboard (6 tabs) |
| [mgr_references.md](../../../migraine/mgr_references.md) | 50 references (including PubMed links) |

---

## Mechanistic Map Preview

[![Migraine QSP Model](../../../migraine/mgr_qsp_model.png)](../../../migraine/mgr_qsp_model.svg)

*Click to open an enlargeable SVG image.*

---

## Key Clinical Trial Benchmarks

| Trial | Drug | Primary endpoint | Result |
|----------|------|-------------|------|
| STRIVE | Erenumab 140mg | MMD reduction | -3.7 days (vs -1.8 placebo) |
| EVOLVE-1 | Galcanezumab 120mg | MMD reduction | -4.7 days (vs -2.8 placebo) |
| HALO-EM | Fremanezumab 225mg | MMD reduction | -3.7 days (vs -2.5 placebo) |
| SAMURAI | Lasmiditan 200mg | 2h pain-free | 32.2% (vs 15.3% placebo) |
| ARTISAN-EM | Rimegepant 75mg | 2h pain-free | 21.2% (vs 10.9% placebo) |
| ACHIEVE-I | Ubrogepant 100mg | 2h pain-free | 19.2% (vs 11.8% placebo) |

---

## Generated

2026-06-20 · QSP Disease Model Library (CCR)
