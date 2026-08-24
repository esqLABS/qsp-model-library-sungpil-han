# Chronic Hepatitis C (CHC) — QSP Model

> **QSP Disease Model Library** · A Quantitative Systems Pharmacology (QSP) model automatically generated via Claude Code Routine.
> Parent library → [../README_en.md](../README_en.md) · Category: Infectious/Hepatobiliary

[![HCV QSP Model](HCV_qsp_model.png)](HCV_qsp_model.svg)

## Overview

Hepatitis C virus (HCV) is a blood-borne flavivirus with approximately 58 million people chronically infected worldwide; if untreated, it progresses over 20–30 years from hepatic fibrosis → cirrhosis → hepatocellular carcinoma (HCC). With the advent of direct-acting antivirals (DAAs), 8–12 weeks of treatment can now achieve a sustained virologic response (SVR12) in more than 95% of patients. This QSP model integrates the Perelson/Neumann target-cell-limited viral dynamics framework with modern DAA PK/PD, host immune modules, and hepatic fibrosis dynamics to simulate virologic response and long-term hepatic prognosis across treatment strategies.

## Key Pathways

| Pathway | Key molecule/mechanism | Clinical outcome |
|------|----------------|-----------|
| Viral replication | NS5B RNA-dependent RNA polymerase → negative-strand replication | High viraemia (log10 6–7 IU/mL) |
| NS5A replication complex | NS5A protein → membranous web formation | RNA replication hub |
| NS3/4A protease | Polyprotein processing → maturation of structural/non-structural proteins | Immune evasion (MAVS/TRIF cleavage) |
| Target-cell infection | HCV E2–CD81·SR-BI–CLDN1–OCLN receptor complex | Spread of hepatocyte infection |
| Innate immune evasion | NS3/4A → MAVS/TRIF cleavage → IFN-β suppression | Establishment of persistent infection |
| T-cell exhaustion | Chronic antigen exposure → PD-1/Tim-3 upregulation | CD8+ CTL dysfunction |
| Hepatic fibrosis | Chronic inflammation → hepatic stellate cell (HSC) activation → TGF-β/collagen accumulation | F0→F4 progression, cirrhosis |
| HCC development | Oxidative stress, TP53/β-catenin mutations | 1–5% per year (in cirrhosis) |

## Key Drug Targets

- **Sofosbuvir (SOF, NS5B inhibitor)**: nucleoside analogue → activated intrahepatically to the triphosphate (SOF-TP) → terminates RNA replication
- **Ledipasvir (LED) / Velpatasvir (VEL) / Pibrentasvir (PIB) (NS5A inhibitors)**: inhibit replication-complex assembly and viral secretion; picomolar-level potency
- **Glecaprevir (GLE) (NS3/4A protease inhibitor)**: blocks polyprotein processing → inhibits viral maturation
- **Ribavirin (RBV)**: immunomodulation and induction of error-prone replication, IFN-boosting effect; used in combination for GT3/cirrhosis
- **Peginterferon-α (PEG-IFN)**: JAK-STAT → ISG induction → pan-antiviral effect; now used only in special situations of DAA intolerance

## Model Files

| File | Description |
|------|------|
| [HCV_qsp_model.dot](HCV_qsp_model.dot) | Graphviz mechanistic map source (130+ nodes / 12 clusters) |
| [HCV_qsp_model.svg](HCV_qsp_model.svg) | SVG vector image (scalable) |
| [HCV_qsp_model.png](HCV_qsp_model.png) | PNG image (150 dpi) |
| [HCV_mrgsolve_model.R](HCV_mrgsolve_model.R) | mrgsolve ODE model (20 compartments / 7 treatment scenarios) |
| [HCV_shiny_app.R](HCV_shiny_app.R) | Shiny dashboard (6 tabs) |
| [HCV_references.md](HCV_references.md) | References (65 articles, 17 sections) |

## Mechanistic Map Clusters

| Cluster | Content |
|----------|------|
| Virion Entry & Uncoating | CD81, SR-BI, CLDN1, OCLN receptors, IRES translation, endocytosis |
| HCV Replication Complex | NS5B RNA polymerase, replication complex, NS5A domains |
| Assembly & Secretion | Core/E1/E2 assembly, VLDL pathway, NS2/NS3 processing |
| Innate Immune Evasion | MDA5/RIG-I-MAVS pathway, NS3/4A cleavage, STAT1 suppression |
| Adaptive Immunity | CD4 Th1, CD8 CTL, B cells, PD-1/Tim-3 exhaustion |
| Hepatic Pathology | HSC activation, TGF-β, Metavir F-score, HCC risk |
| DAA PK Compartments | SOF intrahepatic triphosphate, NS5A inhibitor plasma, GLE/PIB PK |
| Drug PD Effects | εp (replication inhibition), εi (infectivity inhibition), combined efficacy |
| Perelson Viral Kinetics | T (target cells), I (infected cells), V (viral RNA) ODEs |
| Clinical Outcomes | SVR12, relapse, resistance-associated substitutions (RAS), extrahepatic complications |
| Viral Dynamics Parameters | β, δ, p, c parameter configuration |
| DAA Regimens | SOF/LED, SOF/VEL, GLE/PIB, PEG-IFN/RBV, etc. |

## mrgsolve ODE Model

### Compartment structure (20 Compartments)

| Category | Compartment | Description |
|------|------|------|
| DAA PK | `SOF_Tp` | Sofosbuvir intrahepatic active triphosphate |
| DAA PK | `LED_p` | Ledipasvir plasma concentration |
| DAA PK | `VEL_p` | Velpatasvir plasma concentration |
| DAA PK | `NS5A_i` | Combined NS5A inhibitor effect |
| DAA PK | `GLE_p` | Glecaprevir plasma concentration |
| DAA PK | `PIB_p` | Pibrentasvir plasma concentration |
| DAA PK | `RBV_p` | Ribavirin plasma concentration |
| DAA PK | `RBV_RBC` | Ribavirin accumulation in red blood cells |
| IFN | `PEGIFN_p` | Peginterferon-α plasma concentration |
| Viral dynamics | `T_cell` | Number of target hepatocytes |
| Viral dynamics | `I_cell` | Number of infected hepatocytes |
| Viral dynamics | `V_rna` | Plasma HCV RNA (IU/mL) |
| Viral dynamics | `V_def` | Defective viral particles |
| Immune | `CTL` | CD8+ CTL (with exhaustion modelling) |
| Immune | `NK_cell` | NK cell activity |
| Immune | `Treg_HCV` | HCV-specific regulatory T cells |
| Hepatic pathology | `ALT` | ALT (surrogate for hepatocyte damage) |
| Hepatic pathology | `Fibro_met` | Metavir fibrosis score (F0–F4) |
| Hepatic pathology | `HSC_act` | Hepatic stellate cell activation index |
| Hepatic pathology | `HCC_idx` | Cumulative HCC risk index |

### Key treatment scenarios

| Scenario | Regimen | Duration | Supporting trials |
|----------|------|------|---------------|
| 1 | SOF/LED combination | 12 weeks | ION-1/2/3 |
| 2 | SOF/VEL combination | 12 weeks | ASTRAL-1/2/3 |
| 3 | GLE/PIB combination | 8 weeks | ENDURANCE, EXPEDITION |
| 4 | PEG-IFN + RBV | 48 weeks | ADVANCE, ILLUMINATE |
| 5 | SOF/VEL (cirrhosis) | 24 weeks | ASTRAL-4 |
| 6 | SOF/LED + RBV (GT3/RAS) | 24 weeks | LONESTAR |
| 7 | Untreated (natural history) | — | Natural history cohort |

### Key parameters (Perelson Viral Kinetics)

| Parameter | Symbol | Value | Source |
|----------|------|-----|------|
| Viral clearance rate | c | 22 /day | Neumann et al., Science 1998 |
| Infected-hepatocyte death rate | δ | 0.08–0.15 /day | Perelson et al., J Theor Biol 2013 |
| Viral production rate | p | 100 virions/cell/day | Dahari et al., Hepatology 2007 |
| Infection transmission coefficient | β | 1.5×10⁻⁷ /virion/day | Rong et al., Biophys J 2010 |
| SOF production-inhibition rate | εp_SOF | 0.999 | INSPIRE model |
| NS5A inhibition rate | εp_NS5A | 0.9999 | Based on picomolar IC50 |

## Shiny Dashboard

Comprises 6 tabs:

1. **Patient profile**: set genotype (GT1–6), fibrosis grade (F0–F4), baseline viral load, IL-28B genotype
2. **Pharmacokinetics (PK)**: time series of DAA plasma concentrations and intrahepatic SOF-TP concentration; εp/εi dynamics
3. **PD key indicators**: plasma HCV RNA log10 IU/mL, ALT, number of infected hepatocytes
4. **Clinical endpoints**: virologic response timeline (RVR/EVR/SVR12), fibrosis progression, HCC risk
5. **Scenario comparison**: comparison of 1-year outcomes across 7 treatment strategies, DataTable summary
6. **Immune landscape**: CTL dynamics (including exhaustion), NK/Treg ratio, CTL:Treg ratio

## Usage

```r
library(mrgsolve)
mod <- mread("HCV_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny dashboard:
shiny::runApp("HCV_shiny_app.R")
```
```bash
# Render the mechanistic map
dot -Tsvg HCV_qsp_model.dot -o HCV_qsp_model.svg
dot -Tpng -Gdpi=150 HCV_qsp_model.dot -o HCV_qsp_model.png
```

## References

See [HCV_references.md](HCV_references.md) for full citations (65 articles, 17 sections):

- Epidemiology · virology · viral dynamics modelling
- Innate/adaptive immunity · CTL exhaustion · NK cells
- Each class of DAA (NS5B / NS5A / NS3 protease / RBV / PEG-IFN)
- Resistance-associated substitutions (RAS) · SVR outcomes
- Fibrosis dynamics · extrahepatic complications · WHO elimination goals · IL-28B genotype

---
*This model is a qualitative/semi-quantitative QSP model intended for educational and research purposes, and cannot be used directly for real clinical decision-making.*
