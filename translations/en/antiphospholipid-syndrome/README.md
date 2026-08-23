# Antiphospholipid Syndrome (APS) — QSP Model

> **QSP Disease Model Library** · A quantitative systems pharmacology (QSP) model automatically generated via the Claude Code Routine.
> Parent library → [../README.md](../README.md) · Category: Autoimmune/Rheumatic

[![APS QSP Model](../../../antiphospholipid-syndrome/aps_qsp_model.png)](../../../antiphospholipid-syndrome/aps_qsp_model.svg)

## Overview

Antiphospholipid syndrome (APS) is an autoimmune thrombotic disorder in which antiphospholipid antibodies (aPL: lupus anticoagulant, anticardiolipin antibodies, and anti-β₂-glycoprotein I antibodies) cause recurrent arterial and venous thrombosis and pregnancy complications. It occurs in approximately 30% of patients with systemic lupus erythematosus (SLE) (secondary APS) and can also occur independently (primary APS). aPL binds β₂-GPI to activate endothelial cells, platelets, and monocytes, and promotes thrombus formation by driving complement activation (C5a) and tissue factor (TF) expression. Patients who are triple positive (positive for all three aPL antibodies) carry the highest risk of thrombotic recurrence. Anticoagulation (warfarin/LMWH) and hydroxychloroquine are the core therapies, with rituximab and belimumab attempted in refractory disease.

## Key Pathways

| Pathway | Key molecule/mechanism | Clinical outcome |
|------|----------------|-----------|
| aPL-β₂GPI complex | aPL binds β₂-GPI → activates endothelial cells and platelets | Establishes a prothrombotic state |
| Endothelial activation | NF-κB → TF, ICAM-1, VCAM-1 expression | Arterial and venous thrombosis |
| Complement activation | C5a → neutrophil activation, NET formation | Thrombosis, pregnancy complications |
| Platelet activation | aPL → GpIbα receptor → platelet aggregation | Arterial thrombosis, thrombocytopenia |
| Excess coagulation | Increased thrombin generation, inhibition of protein C/S | Amplified thrombotic risk |
| mTOR/renal | Endothelial mTOR activation → renal microthrombosis | APS nephropathy (CAPS risk) |
| Placental complement | Complement deposition on placental trophoblast cells → cell damage | Pregnancy loss, pre-eclampsia |

## Drug Targets

- **Warfarin (vitamin K antagonist)**: standard anticoagulation for venous thrombotic APS; target INR 2.0–3.0 (3.0–4.0 in high-risk patients)
- **Low molecular weight heparin (LMWH)**: the core therapy for APS in pregnancy, including a placental protective effect
- **Hydroxychloroquine (HCQ)**: inhibits TLR signalling → reduces aPL production and platelet aggregation; prevents thrombosis in SLE-APS
- **Rivaroxaban (direct factor Xa inhibitor)**: inferior to warfarin in triple-positive patients; an alternative only in uncomplicated venous thrombosis
- **Rituximab (anti-CD20 antibody)**: reduces aPL through B-cell depletion in refractory APS and CAPS (catastrophic APS)
- **Low-dose aspirin**: combined for prevention in high-risk arterial thrombosis and pregnancy-related APS

## Model Files

| File | Description |
|------|------|
| [aps_qsp_model.dot](../../../antiphospholipid-syndrome/aps_qsp_model.dot) | Graphviz mechanistic map source (approximately 186 nodes / 13 clusters) |
| [aps_qsp_model.svg](../../../antiphospholipid-syndrome/aps_qsp_model.svg) | SVG vector image (scalable) |
| [aps_qsp_model.png](../../../antiphospholipid-syndrome/aps_qsp_model.png) | PNG image (150 dpi) |
| [aps_mrgsolve_model.R](../../../antiphospholipid-syndrome/aps_mrgsolve_model.R) | mrgsolve ODE model (approximately 22 compartments / multiple treatment scenarios) |
| [aps_shiny_app.R](../../../antiphospholipid-syndrome/aps_shiny_app.R) | Shiny dashboard |
| [aps_references.md](../../../antiphospholipid-syndrome/aps_references.md) | References (approximately 58 articles, PubMed links) |

## mrgsolve Model (ODE Model)

- **Compartment structure**: warfarin (2 oral compartments + effect compartment) + LMWH + HCQ (2 compartments) + rivaroxaban (2 compartments) + aspirin PK compartment + rituximab (2 compartments) + aPL IgG, B cells, complement C5a, endothelial TF, platelet activation, thrombin generation, deep vein thrombosis risk, pregnancy survival rate, mTOR renal, and INR PD compartments
- **Key treatment scenarios**: ① no treatment, ② warfarin monotherapy, ③ LMWH monotherapy (pregnancy), ④ warfarin + HCQ, ⑤ rivaroxaban (low risk), ⑥ aspirin + HCQ, ⑦ rituximab + anticoagulation (refractory APS)
- **Calibration/evidence**: parameters referenced from the RAPS trial (rivaroxaban vs. warfarin), Crowther et al. NEJM warfarin INR targets, and the Beppu et al. CAPS registry

## Shiny Dashboard

The dashboard comprises a patient profile tab (aPL profile, thrombosis history, pregnancy status), an anticoagulant PK and INR/anti-Xa level tab, a thrombus formation risk dynamics tab, a clinical endpoints tab (thrombotic recurrence, pregnancy outcomes), a treatment scenario comparison tab, and a biomarker tab (aPL titre, complement, TF).

## Usage

```r
library(mrgsolve)
mod <- mread("aps_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny dashboard:
shiny::runApp("aps_shiny_app.R")
```
```bash
# Render the mechanistic map
dot -Tsvg aps_qsp_model.dot -o aps_qsp_model.svg
```

## References

For detailed citations, see [aps_references.md](../../../antiphospholipid-syndrome/aps_references.md) (approximately 58 articles).

---
*This model is a qualitative/semi-quantitative QSP model intended for education and research, and must not be used directly for clinical decision-making.*
