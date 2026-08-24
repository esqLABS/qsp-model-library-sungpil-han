# Chronic Hepatitis B (CHB) — QSP Model

> **QSP Disease Model Library** · A quantitative systems pharmacology (QSP) model automatically generated via the Claude Code Routine.
> Parent library → [../README_en.md](../README_en.md) · Category: Gastroenterology/Hepatobiliary

[![CHB QSP Model](chb_qsp_model.png)](chb_qsp_model.svg)

## Overview

Chronic hepatitis B (CHB) is a public health problem affecting approximately 296 million people worldwide and a major cause of cirrhosis and hepatocellular carcinoma (HCC). HBV forms cccDNA (covalently closed circular DNA) in the nucleus of infected hepatocytes, which persists even after antiviral therapy, and its natural history is complex, comprising immune-tolerant, immune-active, and inactive phases. Nucleos(t)ide analogues (entecavir, tenofovir) suppress viral replication and slow the progression of liver fibrosis, while pegylated interferon-α2a combines immune modulation with a direct antiviral effect.

## Key Pathways

| Pathway | Key molecule/mechanism | Clinical outcome |
|------|----------------|-----------|
| cccDNA persistence/transcription | Nuclear cccDNA → pgRNA → HBV DNA replication cycle | Viral persistence, immune evasion |
| Immune tolerance/T-cell exhaustion | CTL exhaustion, increased immune checkpoint (PD-1/Tim-3) expression | Chronicity, delayed transition to the immune-active phase |
| Innate immune suppression | HBV blockade of IFN signalling (cGAS-STING) | Failure of early viral clearance |
| Hepatocyte injury | CTL-mediated cytolysis, raised ALT | Risk of acute exacerbation of hepatitis |
| Hepatic stellate cell (HSC) activation | TGF-β, TNF-α → HSC → collagen synthesis | Progression of liver fibrosis (F0→F4) |
| HCC risk | High viraemia, cirrhosis, persistent HBsAg | Annual HCC risk of 0.5~3% |
| HBsAg production/secretion | cccDNA transcription plus integrated HBV DNA | Sustained serum HBsAg positivity |

## Drug Targets

- **Entecavir (ETV)**: dNTP competition → inhibits HBV reverse transcriptase, high genetic barrier to resistance
- **Tenofovir disoproxil (TDF)/TAF**: a nucleotide analogue → inhibits HBV DNA polymerase
- **Pegylated interferon-α2a (Peg-IFN)**: activates innate and adaptive immunity plus direct antiviral action, with potential for HBsAg loss
- **siRNA/ASO (under investigation)**: degrades HBsAg mRNA → reduces serum HBsAg concentration → induces immune recovery
- **Capsid assembly modulators (CAM)**: a novel mechanism that blocks cccDNA replenishment

## Model Files

| File | Description |
|------|------|
| [chb_qsp_model.dot](chb_qsp_model.dot) | Graphviz mechanistic map source (approximately 438 nodes / 10 clusters) |
| [chb_qsp_model.svg](chb_qsp_model.svg) | SVG vector image (scalable) |
| [chb_qsp_model.png](chb_qsp_model.png) | PNG image (150 dpi) |
| [chb_mrgsolve_model.R](chb_mrgsolve_model.R) | mrgsolve ODE model (approximately 22 compartments / 6 treatment scenarios) |
| [chb_shiny_app.R](chb_shiny_app.R) | Shiny dashboard |
| [chb_references.md](chb_references.md) | References (approximately 54 articles, PubMed links) |

## mrgsolve Model (ODE Model)

- **Compartment structure**: ETV/TDF PK (gut, plasma, intracellular triphosphate) + Peg-IFN SC/plasma + PD compartments for target cells (T), infected cells (I), HBV DNA (V), cccDNA, HBsAg, CTL response, T-cell exhaustion, innate IFN, ALT, HSC activation, liver fibrosis, and cumulative HCC risk
- **Key treatment scenarios**: (1) no treatment, (2) entecavir 0.5 mg QD, (3) TDF 300 mg QD, (4) Peg-IFN-α2a for 48 weeks, (5) ETV + Peg-IFN combination, (6) ETV + siRNA combination
- **Calibration/evidence**: viral suppression kinetics for entecavir and tenofovir based on the Perelson AS et al. viral dynamics model (Hepatology 2012); HCC risk calibrated from PAGE-B score data

## Shiny Dashboard

Comprises 6 tabs: (1) Patient Profile — sets immune phase, baseline viral load, and liver fibrosis stage; (2) PK tab — ETV/TDF/Peg-IFN concentration time series; (3) Key PD Metrics — HBV DNA (log10), cccDNA, HBsAg, CTL; (4) Clinical Endpoints — ALT normalisation, liver fibrosis, cumulative HCC risk; (5) Scenario Comparison — 5-year outcomes across 6 treatment strategies; (6) Biomarkers — quantitative HBsAg, HBeAg, platelets (a surrogate marker for cirrhosis)

## Usage

```r
library(mrgsolve)
mod <- mread("chb_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny dashboard:
shiny::runApp("chb_shiny_app.R")
```
```bash
# Render the mechanistic map
dot -Tsvg chb_qsp_model.dot -o chb_qsp_model.svg
```

## References

For detailed citations, see [chb_references.md](chb_references.md) (approximately 54 articles).

---
*This model is a qualitative/semi-quantitative QSP model intended for education and research, and must not be used directly for clinical decision-making.*
