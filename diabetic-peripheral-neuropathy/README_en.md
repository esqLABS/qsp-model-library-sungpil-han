# Diabetic Peripheral Neuropathy (DPN) — QSP Model

> A Quantitative Systems Pharmacology model integrating the causal chain of diabetic peripheral neuropathy: quadruple metabolic injury (polyol · AGE · PKC · hexosamine) → microvascular/endoneurial hypoperfusion → nerve fibre damage → pain · sensory loss →
> leading to foot ulceration.

| Item | Value |
|------|-----|
| Directory | `diabetic-peripheral-neuropathy/` |
| Abbreviation | `dpn` |
| Category | Chronic disease · Neurological/Endocrine |
| Core ODE compartments | 26 (PK 10 + disease/pain 16) |
| Drug scenarios | 8 (untreated · pregabalin · duloxetine · combination · α-lipoic acid · epalrestat · capsaicin 8% · intensive glycaemic control + combination) |
| Shiny tabs | 8 |
| Mechanistic map nodes | 130+ (10 clusters) |
| References | 60 (18 sections) |

## 1. Deliverables
| File | Description |
|------|------|
| [`dpn_qsp_model.dot`](dpn_qsp_model.dot) / [`.svg`](dpn_qsp_model.svg) / [`.png`](dpn_qsp_model.png) | Graphviz mechanistic map |
| [`dpn_mrgsolve_model.R`](dpn_mrgsolve_model.R) | mrgsolve QSP ODE model |
| [`dpn_shiny_app.R`](dpn_shiny_app.R) | Shiny dashboard (8 tabs) |
| [`dpn_references.md`](dpn_references.md) | 60 PubMed citations |

## 2. Pathophysiology Summary
1. **Hyperglycaemia → intracellular glucose accumulation**: influx via GLUT1 in insulin-independent neurons/Schwann cells/endothelial cells → simultaneous activation of four injury pathways.
2. **Polyol pathway**: aldose reductase reduces glucose→sorbitol→fructose; NADPH consumption causes GSH depletion and osmotic stress.
3. **AGE/RAGE**: non-enzymatic protein glycation → RAGE → NF-κB → TNF-α · IL-6 · VCAM-1.
4. **PKC-β**: DAG↑ → eNOS↓ · ET-1↑ · VEGF abnormality.
5. **Hexosamine**: UDP-GlcNAc generation via GFAT → O-GlcNAc modification → regulation of GAPDH · PAI-1.
6. **Mitochondrial superoxide · PARP**: NAD⁺ depletion further amplifies the four pathways (Brownlee unified hypothesis).
7. **Microvascular injury**: vasa nervorum stenosis · basement membrane thickening · endoneurial hypoperfusion → endoneurial hypoxia.
8. **Nerve fibre damage**: distal die-back, reduced IENFD, decreased NCV, NGF↓, Schwann cell damage.
9. **Pain mechanism**: Nav1.7/1.8 hyperexcitability, TRPV1/TRPA1 sensitisation, spinal dorsal horn NMDA · microglial activation → central sensitisation.
10. **Clinical outcomes**: worsening pain NRS · MNSI · TCNS · BPI, impaired balance/sleep, foot ulceration → amputation → mortality.

## 3. Drug PK/PD
| Drug | Target | Modelling Depth |
|------|------|------------|
| Pregabalin | Cav α2δ-1 → ↓ glutamate release | First-order absorption · renal CL adjusted for eGFR, EC50_pain=4 mg/L, Emax=0.45 |
| Duloxetine | SNRI → descending 5-HT/NE | CYP2D6 phenotype adjustment, EC50=0.08 mg/L, Emax=0.35 |
| α-Lipoic acid | ROS scavenging · GSH regeneration | Oral F=0.30, IV F=1.0, EC50=0.5 mg/L |
| Epalrestat | Aldose reductase | IC50=5 mg/L, polyol inhibition |
| Capsaicin 8% | TRPV1 desensitisation | Effect compartment, τ=90 days |
| Lidocaine 5% | Nav1.7/1.8 blockade | Local effect compartment |
| Ruboxistaurin*, Aminoguanidine* | Failed drugs (reference line) | Negative control |

(*) Ph3 negative; included for negative-control purposes.

## 4. Calibration Anchors
- **DCCT/EDIC**: intensive glycaemic control → 60% reduction in 5-year neuropathy incidence.
- **SENZA-PDN HF10 SCS**: ≥50% pain reduction at 3 months in 76% vs 5% with CMM.
- **ALADIN/NATHAN-1**: ALA 600 mg IV for 3 weeks → TSS Δ −2.7; the 4-year NATHAN-1 trial showed NIS-LL improvement in some subgroups.
- **Pregabalin 300-600 mg/d**: NRS Δ −1.3 vs placebo.
- **Duloxetine 60-120 mg/d**: NRS Δ −1.4.
- **Capsaicin 8% (STEP)**: NRS Δ −1.0, sustained for 12 weeks.
- **Epalrestat (ADCT)**: NCV stabilisation at 3 years.

## 5. Usage
```bash
# Render the mechanistic map
dot -Tsvg dpn_qsp_model.dot -o dpn_qsp_model.svg
dot -Tpng -Gdpi=150 dpn_qsp_model.dot -o dpn_qsp_model.png
```
```r
# mrgsolve model (R)
install.packages(c("mrgsolve","dplyr","ggplot2"))
library(mrgsolve)
mod <- mread("dpn_mrgsolve_model.R")
ev  <- ev(amt=150, ii=12, addl=730, cmt="GUT_PG")   # Pregabalin 300 mg BID
out <- mrgsim(mod, events=ev, end=365, delta=1)
plot(out, NRS+IENFD+NCV~time)

# Shiny dashboard
install.packages(c("shiny","shinydashboard","DT","tidyr","scales"))
shiny::runApp("dpn_shiny_app.R")
```

## 6. Limitations
- Based on an average patient phenotype — population simulation (virtual patient cohort) requires future extension.
- Autonomic neuropathy (CAN), painless large-fibre variants, and exercise/therapeutic neurostimulation are represented only qualitatively.
- Capsaicin and lidocaine patches are abstracted as simple effect compartments — pharmacokinetic local exposure is not represented.
- All parameters are starting points based on literature and clinical trial mean values, and require re-estimation when fitted to real data.

## 7. Change History
- **2026-06-30** v1.0 — Initial model (Claude Code Routine).
