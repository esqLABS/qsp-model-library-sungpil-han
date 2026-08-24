# Narcolepsy Type 1 — QSP Model

> **Disease category**: Neurological / Sleep Disorder  
> **Directory**: `narcolepsy/`  
> **Abbreviation**: `narc`

---

## Pathophysiology

Narcolepsy type 1 (NT1) is a chronic neurological disease that arises from the
**selective, autoimmune-mediated destruction of orexin (hypocretin) neurons in the
lateral hypothalamus**. Whereas healthy individuals have about 70,000 orexin neurons,
NT1 patients lose **85–95% of them**, and CSF hypocretin-1 (OXA) concentration falls
below 110 pg/mL.

### Key pathological pathways

| Pathway | Detailed mechanism |
|------|----------|
| **Autoimmune pathogenesis** | HLA-DQB1\*06:02 susceptibility gene + molecular mimicry with H1N1 influenza/Pandemrix vaccine → CD4+/CD8+ T cells and NK cells destroy orexin neurons |
| **Orexin receptor system** | OXA (OX1R≫OX2R), OXB (OX2R≥OX1R) → dual Gq/Gi signaling → activates wake-promoting aminergic systems |
| **Collapse of the wake-promoting system** | Instability in LC (NE) · TMN (histamine) · VTA (DA) · DRN (5-HT) · basal forebrain (ACh) activity |
| **Flip-flop switch instability** | Orexin deficiency → disruption of VLPO–wake-system mutual inhibition balance → frequent switching between wake/sleep states |
| **Cataplexy** | Amygdala emotional stimulation → involuntary activation of the REM atonia circuit (SubC glutamatergic system) |
| **Metabolic complications** | Orexin deficiency → leptin resistance and energy-balance abnormalities → increased risk of obesity, T2DM |

---

## Model Files

| File | Description |
|------|------|
| [`narc_qsp_model.dot`](../../../narcolepsy/narc_qsp_model.dot) | Graphviz mechanistic map (source) |
| [`narc_qsp_model.svg`](../../../narcolepsy/narc_qsp_model.svg) | Vector graphic (zoomable) |
| [`narc_qsp_model.png`](../../../narcolepsy/narc_qsp_model.png) | Raster image (150 dpi) |
| [`narc_mrgsolve_model.R`](../../../narcolepsy/narc_mrgsolve_model.R) | mrgsolve ODE QSP model |
| [`narc_shiny_app.R`](narc_shiny_app.R) | Shiny interactive dashboard |
| [`narc_references.md`](narc_references.md) | References (35+ PubMed citations) |

---

## Mechanistic Map

[![Narcolepsy QSP Map](../../../narcolepsy/narc_qsp_model.png)](../../../narcolepsy/narc_qsp_model.svg)

### Cluster structure (12 subgraphs)

1. **Autoimmune pathogenesis** — HLA-DQB1\*06:02, T cells, molecular mimicry, autoantibodies
2. **Hypothalamic orexin neurons** — HCRT gene, OXA/OXB, neuronal loss
3. **Orexin receptor system** — OX1R, OX2R, Gq/Gi signaling
4. **Wake-promoting monoamine systems** — LC, TMN, VTA, DRN, basal forebrain
5. **Sleep-promoting system** — VLPO, adenosine, circadian rhythm (SCN)
6. **REM sleep regulation and cataplexy** — PPT/LDT, SubC, amygdala
7. **The four cardinal clinical symptoms** — EDS, cataplexy, hypnagogic hallucinations, sleep paralysis
8. **Drug PK compartments** — sodium oxybate, modafinil, pitolisant, solriamfetol
9. **Pharmacological targets and effects** — GABA-B, DAT/NET, H3R, SERT, OX2R
10. **Biomarkers and endpoints** — CSF OXA, MSLT, PSG, ESS
11. **Flip-flop switch** — wake/sleep bistable state, NT1 instability
12. **Metabolic complications** — obesity, T2DM, depression/anxiety, ADHD-like symptoms

---

## mrgsolve ODE Model (Pharmacological Model)

### Compartments (23)

| Compartment | Variable name | Description |
|------|--------|------|
| **Sodium oxybate PK** | GUT_OXY, CENT_OXY, PERI_OXY | 3-compartment model, t½ ≈ 30–60 min |
| **Modafinil PK** | GUT_MOD, CENT_MOD | 2-compartment, t½ ≈ 15 hours |
| **Pitolisant PK** | GUT_PIT, CENT_PIT | 2-compartment, t½ ≈ 10–12 hours |
| **Solriamfetol PK** | GUT_SOL, CENT_SOL | 2-compartment, t½ ≈ 7.1 hours |
| **Venlafaxine PK** | GUT_VEN, CENT_VEN | 2-compartment, anti-cataplectic |
| **Wake-promoting system** | WAKE_LC, WAKE_TMN, WAKE_VTA, WAKE_DRN | Monoaminergic firing rate (0–1) |
| **Sleep system** | SLEEP_P, VLPO_ACT, ADENOSINE | Homeostatic sleep pressure · VLPO · adenosine |
| **State variables** | WAKE_STATE, REM_STATE, NREM_STATE | Sleep/wake/REM flip-flop |
| **Clinical accumulation** | EDS_ACC, CATAPLEXY_ACC | Cumulative sleepiness · cataplexy |

### Treatment Scenarios (7)

| Scenario | Drug | Dose |
|----------|------|------|
| 1 | Untreated NT1 baseline | — |
| 2 | Sodium oxybate | 4.5 g + 2.75 g split doses (at bedtime) |
| 3 | Modafinil | 200 mg once daily |
| 4 | Pitolisant | 18 mg once daily |
| 5 | Solriamfetol | 150 mg once daily |
| 6 | Sodium oxybate + pitolisant | Combination |
| 7 | Venlafaxine | 75 mg (focused on cataplexy) |

### Key clinical trial calibration

| Reference | Treatment | Result |
|------|------|------|
| Black 2010 (Sleep Med) | Sodium oxybate | Cataplexy -69 to 75% |
| HARMONY I (Szakacs 2017, Lancet Neurol) | Pitolisant | ESS -5.8 vs. placebo -3.4 |
| TONES 3 (Schweitzer 2019, Sleep) | Solriamfetol | ESS -7.7 points |
| US MMSG 2000 | Modafinil | ESS -4.3 points |
| Mignot 2002 (Lancet) | CSF orexin-1 | Diagnostic criterion <110 pg/mL |

---

## Shiny Dashboard (Interactive Dashboard)

```r
shiny::runApp("narcolepsy/narc_shiny_app.R")
```

### Tab structure (8)

| Tab | Content |
|----|------|
| 1. Patient Profile | Orexin neuron survival rate, ESS, cataplexy frequency, HLA status |
| 2. Drug PK | Drug concentration-time curves, Cmax/Tmax/AUC table |
| 3. Sleep-Wake Regulation | Monoamine system activity, flip-flop switch, 24-hour circadian cycle |
| 4. Clinical Endpoints | Tracking of ESS · cataplexy · MSLT · sleep architecture |
| 5. Treatment Scenario Comparison | Clinical-trial benchmark waterfall chart · forest plot |
| 6. Biomarkers | CSF hypocretin-1, MSLT, HLA diagnostic visualization |
| 7. Autoimmune Mechanism | HLA association, neuronal destruction timeline, T cells |
| 8. References | Summary table of key clinical trials |

---

## Usage

```bash
# Render the mechanistic map (requires Graphviz)
dot -Tsvg narcolepsy/narc_qsp_model.dot -o narc_qsp_model.svg
dot -Tpng -Gdpi=150 narcolepsy/narc_qsp_model.dot -o narc_qsp_model.png
```

```r
# Run the mrgsolve model (R ≥ 4.0)
install.packages(c("mrgsolve", "dplyr", "ggplot2"))
source("narcolepsy/narc_mrgsolve_model.R")

# Run the Shiny dashboard
install.packages(c("shiny", "shinydashboard", "plotly"))
shiny::runApp("narcolepsy/narc_shiny_app.R")
```

---

## Key Parameters

| Parameter | Value | Source |
|---------|-----|------|
| Normal CSF orexin-1 | 200–300 pg/mL | Mignot 2002 |
| NT1 diagnostic criterion | <110 pg/mL | ICSD-3 |
| Orexin neuron loss | 85–95% | Thannickal 2000 |
| HLA-DQB1\*06:02 prevalence | ~95% (NT1) vs. 25% (general population) | Mignot 1997 |
| Prevalence | 25–50/100,000 | Longstreth 2007 |
| Average diagnostic delay | 10–15 years | Thorpy 2015 |
| Emax (pitolisant → ESS) | -5.8 points | HARMONY I |
| Emax (solriamfetol → ESS) | -7.7 points | TONES 3 |

---

## Disclaimer

This model is a qualitative/semi-quantitative QSP model built for **educational and
research purposes**. It is constructed from the public literature and must not be used
directly for clinical decision-making.
