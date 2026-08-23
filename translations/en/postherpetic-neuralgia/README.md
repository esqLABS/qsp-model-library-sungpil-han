# Postherpetic Neuralgia (PHN) — QSP Model

This is a Quantitative Systems Pharmacology (QSP) model package for postherpetic
neuralgia (PHN). It integrates in one place the pathophysiology running from
reactivation of the varicella-zoster virus (VZV) → ganglion injury →
peripheral/central sensitisation → diverse pain phenotypes, together with the
pharmacology of antivirals/preventive vaccine (Shingrix), gabapentinoids,
TCAs/SNRIs, 5% lidocaine patch, 8% capsaicin patch, NMDA blockers, opioids, and
future targets such as NaV1.7/NGF/P2X3.

## Files

| File | Description |
|---|---|
| `phn_qsp_model.dot` | Graphviz mechanistic map — 13 clusters, 170+ nodes (VZV latency-reactivation, nerve injury, peripheral/central sensitisation, neuroinflammation, clinical phenotypes, vaccine/antiviral/multi-drug pharmacology, PK/safety/patient covariates) |
| `phn_qsp_model.svg` / `.png` | DOT rendering output (svg + 150 dpi png) |
| `phn_mrgsolve_model.R` | mrgsolve QSP model — 24 ODE compartments (PK: GBP · PGB · AMI/NOR · DLX · LIDO · CAP · VAL · TRA · RZV; PD: VZV_LOAD, GANG_INJ, IENF, NAV_ACT, CSEN, MICROG, KCC2, NMDA_TONE, NGF, CMI, PAIN, ALLO, SLEEP, MOOD, AE_SED). 10 pre-defined scenarios |
| `phn_shiny_app.R` | Shiny dashboard — 8 tabs (Patient · PK · PD physiology · Clinical endpoints · Scenarios · Vaccine · Safety · Biomarkers/QST) |
| `phn_references.md` | 80 PubMed references (natural history · mechanism · vaccine · antiviral · drug-specific RCTs · QSP/PK models · QST phenotypes) |

## How to Use

```r
# 1) Render the DOT (SVG/PNG already included):
# dot -Tsvg phn_qsp_model.dot -o phn_qsp_model.svg
# dot -Tpng -Gdpi=150 phn_qsp_model.dot -o phn_qsp_model.png

# 2) Load the mrgsolve model & simulate scenarios
source("phn_mrgsolve_model.R")
out <- mrgsim_e(mod, e_combo, end = 180*24, delta = 24) |> as.data.frame()

# 3) Shiny dashboard
shiny::runApp("phn_shiny_app.R")
```

## Scenario Library (10)

1. **Placebo** — natural history
2. **Valaciclovir 1 g q8h × 7 d** — antiviral for acute zoster
3. **RZV (Shingrix)** — prevention via 2-dose vaccination at 0, 2-6 months
4. **Gabapentin titration → 3600 mg/d**
5. **Pregabalin 75 → 150 bid (300 mg/d)**
6. **Amitriptyline 25 → 75 mg HS**
7. **Duloxetine 30 → 60 mg qd**
8. **Lidocaine 5% patch** applied once daily for 12 h
9. **Capsaicin 8% patch** single 60-min application, reapplied at 90-day intervals
10. **Combo** (AV + RZV + PGB + lidocaine + AMI) — standardised integrated management

## Calibration Anchors

- RZV efficacy: HZ 97% (ages 50-69), 91% (≥70), PHN 88-91% (ZOE-50/70)
- Valaciclovir 1 g tid → ~30% shortening of acute pain resolution time (Beutner 1995, Tyring 1995)
- Pregabalin 300-600 mg → ~35-50% responder rate for 50% pain reduction (Dworkin 2003)
- Gabapentin 1800-3600 mg → ~30% reduction in NRS (Rice 2001)
- Amitriptyline 25-100 mg HS → NNT ~2.7 (Watson 1982)
- Capsaicin 8%: ~30% reduction in NRS over 12 weeks after a single application (Backonja 2008 / STRIDE)
- Lidocaine 5% patch: 30-50% reduction in allodynia, very low systemic Css (Galer 2002)

## Caveats

- The pain ODE is phenomenological — for in-silico scenario comparison and trial design purposes
- PK is simplified to 1/2-compartment models; popPK variability is not included
- VZV reactivation drive is abstracted as a single first-order variable
- The effect of capsaicin on IENF is assumed to be the average pharmacological defunctionalization effect
