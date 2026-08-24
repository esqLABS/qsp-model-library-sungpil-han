# Hepatorenal Syndrome (HRS) — QSP Model

> **QSP Disease Model Library** · A Quantitative Systems Pharmacology (QSP) model automatically generated via Claude Code Routine.
> Parent library → [../README.md](../README.md) · Category: Gastroenterology & Hepatobiliary / Renal

[![HRS QSP Model](../../../hepatorenal-syndrome/hrs_qsp_model.png)](../../../hepatorenal-syndrome/hrs_qsp_model.svg)

## Overview
Hepatorenal syndrome (HRS) is a functional, potentially reversible form of
acute kidney injury that arises in advanced cirrhosis with portal
hypertension. Its core pathophysiology is that profound splanchnic
vasodilation causes a fall in effective arterial blood volume (EABV),
and the RAAS, sympathetic nervous system, and vasopressin (AVP) become
maximally activated in compensation, strongly constricting the renal
vasculature. The 2015 International Club of Ascites (ICA) criteria and
their 2019 revision reclassified the condition into HRS-AKI (acute) and
HRS-NAKI (non-AKI, subacute/chronic); it is a fatal complication with a
30-day mortality of 40–90%. Combination therapy with a vasoconstrictor —
including terlipressin, approved by the FDA in 2022 — plus albumin is the
standard of care, and the definitive treatment is liver transplantation
(± simultaneous liver-kidney transplant, SLK).

## Key Pathophysiological Pathways
| Pathway | Key molecules/mechanism | Clinical result |
|------|----------------|-----------|
| Splanchnic vasodilation | eNOS/iNOS, NO, CO, PGI₂, glucagon ↑ | SVR ↓, MAP ↓, splanchnic pooling |
| Reduced EABV | Lowered effective blood volume → baroreceptor unloading | Hyperdynamic circulation, RAP ↓ |
| Neurohormonal activation | RAAS, sympathetic, AVP ↑↑ | Renal vasoconstriction, sodium/water retention |
| Renal vasoconstriction | AngII, NE, ET-1, TXA₂, adenosine, endocannabinoids | RBF ↓, GFR ↓, sCr ↑ |
| Loss of renal prostaglandins | PGE₂/PGI₂ ↓ (vulnerable to NSAIDs) | Worsened afferent-arteriole constriction |
| Cirrhotic cardiomyopathy | β-adrenergic downregulation, diastolic dysfunction | Insufficient CO under stress, precipitates HRS |
| Systemic inflammation | Bacterial translocation, LPS, TNF-α, IL-6, CRP | SIRS/ACLF, multi-organ dysfunction |
| Precipitating factors | SBP, GI bleeding, LVP without albumin, NSAIDs | Acute EABV reduction, onset of HRS-AKI |

## Key Drug Targets
- **Terlipressin**: a lysyl-vasopressin prodrug, V1a agonist → splanchnic
  vasoconstriction, increases MAP, RBF, and GFR (CONFIRM/REVERSE/OT-0401
  RCTs). Given as a bolus (1 mg q6h) or continuous infusion (2 mg/24h).
- **Norepinephrine**: an α1-adrenergic agonist → systemic vasoconstriction,
  an alternative option in the ICU setting, for terlipressin
  non-responders or where it is unavailable.
- **Midodrine + octreotide**: midodrine (α1, PO) + octreotide (inhibits
  glucagon and VIP, SC) + albumin. Used in outpatient/non-ICU settings.
- **Albumin 25%**: osmotic expansion + anti-inflammatory action +
  endothelial stabilisation. 1 g/kg (day 1) → 20–40 g/day (ATTIRE, ANSWER,
  SBP-prevention data).
- **Antibiotics**: SBP treatment (ceftriaxone), SBP prophylaxis
  (norfloxacin), rifaximin (suppresses bacterial translocation).
- **TIPS · RRT · liver transplantation**: procedural interventions, under
  the bridge-to-transplant concept.

## Model Files
| File | Description |
|------|------|
| [hrs_qsp_model.dot](../../../hepatorenal-syndrome/hrs_qsp_model.dot) | Graphviz mechanistic map source (~148 nodes / 15 clusters) |
| [hrs_qsp_model.svg](../../../hepatorenal-syndrome/hrs_qsp_model.svg) | SVG vector image |
| [hrs_qsp_model.png](../../../hepatorenal-syndrome/hrs_qsp_model.png) | PNG image (150 dpi) |
| [hrs_mrgsolve_model.R](../../../hepatorenal-syndrome/hrs_mrgsolve_model.R) | mrgsolve ODE model (26 compartments / 10 scenarios) |
| [hrs_shiny_app.R](../../../hepatorenal-syndrome/hrs_shiny_app.R) | Shiny dashboard (8 tabs) |
| [hrs_references.md](../../../hepatorenal-syndrome/hrs_references.md) | References (82 entries, PubMed links) |

## mrgsolve Model (ODE Model)
- **Compartment structure**:
  - Drug PK (10): terlipressin 2-cpt + lysyl-vasopressin, norepinephrine 1-cpt, midodrine→desglymidodrine, octreotide SC/central, albumin
  - Neurohormonal (4): PRA, aldosterone, endogenous NE, AVP
  - Haemodynamics (4): MAP, SVR, RBF, cardiac output
  - Renal/clinical (6): GFR, sCr, urine Na, urine output, serum Na, systemic inflammation
  - Outcomes (3): 30-day risk integral, terlipressin ischaemic AUC, bilirubin
- **Key treatment scenarios (10)**:
  1. Natural history
  2. Terlipressin bolus (1 mg q6h) + albumin 40 g/day
  3. Terlipressin continuous infusion (2 mg/24h) + albumin
  4. Norepinephrine CIV + albumin
  5. Midodrine + octreotide + albumin
  6. Albumin alone (ATTIRE-like comparator arm)
  7. Terlipressin response after SBP precipitation
  8. LVP-without-albumin PICD
  9. NSAID-induced worsening
  10. TIPS haemodynamic surrogate scenario
- **Calibration/basis**: CONFIRM (Wong 2021, NEJM), REVERSE (Boyer 2016), OT-0401 (Sanyal 2008), ATTIRE (China 2021), ANSWER (Caraceni 2018), Cavallin 2015/2016 (CI vs bolus).

## Shiny Dashboard (8 Tabs)
1. **Overview** — MAP/RBF/GFR/sCr trends + model schematic
2. **Drug PK** — terlipressin/lysyl-VP, NE, midodrine, octreotide, albumin
3. **Neurohormonal** — RAAS (renin, aldosterone), SNS, AVP
4. **Renal & urinary** — urine Na, urine output, serum Na
5. **Clinical endpoints** — HRS response, 30/90-day survival, MELD trajectory
6. **Scenario comparison** — table/graph comparing 7 regimens
7. **Safety** — terlipressin ischaemic AUC, risk of MAP overshoot
8. **References** — links to key references

## Usage
```r
library(mrgsolve); library(shiny)
source("hrs_mrgsolve_model.R")   # builds hrs_mod
shiny::runApp("hrs_shiny_app.R")
```

## References
[hrs_references.md](../../../hepatorenal-syndrome/hrs_references.md) — 82 references including ICA-2015/2019, CONFIRM, REVERSE, OT-0401, ATTIRE, ANSWER, PREDICT, CANONIC.

## License · Disclaimer
- Governed by the parent repository's license ([../../../LICENSE](../../../LICENSE)).
- This model is for education and research and does not replace clinical
  decision-making. Actual patient care should be guided by current
  guidelines and expert judgement.
