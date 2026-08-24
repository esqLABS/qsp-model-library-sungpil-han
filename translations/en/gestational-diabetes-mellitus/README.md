# Gestational Diabetes Mellitus (GDM) — QSP Model
## Gestational Diabetes Mellitus · Quantitative Systems Pharmacology Model

<a href="../../../gestational-diabetes-mellitus/gdm_qsp_model.svg"><img src="../../../gestational-diabetes-mellitus/gdm_qsp_model.png" width="720" alt="GDM QSP mechanistic map"></a>

*(Click the image to open a zoomable SVG — 208 nodes, 14 clusters, 316 interactions)*

---

## One-Line Summary

**GDM is not "high blood glucose during pregnancy" — it is a race between
two clocks that start running together at the moment pregnancy begins.**
The placenta grows at nearly the same rate in almost every pregnancy,
pouring out counter-regulatory hormones that drop insulin sensitivity by
50-60%. Against this, the beta cell compensates through expansion and a
left-shift, and it is only this compensatory capacity that differs from
person to person. So in this model, **whether GDM occurs is an output, not
an input.**

---

## Why This Structure (The Modelling Thesis)

Existing models of glycaemia in pregnancy generally treat "hyperglycaemia"
as an input and predict fetal outcomes from it. That falls into circular
reasoning whenever treatments are compared — of course outcomes improve
once you have specified that glucose is lowered.

This model reverses the direction of causation. There are only two
patient-level parameters:

| Parameter | Meaning | Range |
|---|---|---|
| `BCAP` | Beta-cell adaptive capacity | 0.25 – 1.20 |
| `BMI` | Pre-pregnancy body mass index → baseline insulin sensitivity | 18 – 42 |

Glucose, insulin, HbA1c, birth weight, LGA probability, and postpartum
diabetes risk are **all outputs.**

### Two Results That Emerge Here

**(1) The Powe physiological subtypes appear without being specified.**
The GW-40 analytical fixed point of the maternal glucose/insulin subsystem
(the derivation is in CALIBRATION §D of `gdm_mrgsolve_model.R`):

| BCAP | BMI | Fasting glucose | Fasting insulin | Phenotype |
|---|---|---|---|---|
| 1.00 | 22 | ~77 mg/dL | ~14 µU/mL | Normal (NGT) — fasting glucose **falls** |
| 0.85 | 24 | ~87 mg/dL | ~13 µU/mL | Secretory-deficient type — fasting near normal, abnormal only after a load |
| 1.00 | 34 | ~100 mg/dL | ~18 µU/mL | Insulin-resistant type — hyperglycaemia + hyperinsulinaemia |
| 0.55 | 31 | ~110 mg/dL | ~13 µU/mL | Mixed/severe type |

Three different fasting/insulin signatures emerge under the same
diagnostic label, and **drug response should differ too** — metformin acts
on the hepatic side (EGP), insulin on the whole system. That the
secretory-deficient type is missed by fasting glucose alone is not a
prediction the model was built to produce; it follows directly from the
structure.

**(2) The fall in fasting glucose seen in normal pregnancy emerges without
being forced.** Fasting glucose at term falling from 85 to 77 mg/dL,
despite worsening insulin resistance, happens because uteroplacental
glucose uptake (the "siphon") outpaces the gentle rise in maximal EGP.
Without hepatic glucose auto-regulation (`KGEGP`), the model predicts
obvious fasting hyperglycaemia in every obese pregnancy — which is simply
the wrong prediction.

---

## File Structure

| File | Contents |
|---|---|
| [`gdm_qsp_model.dot`](../../../gestational-diabetes-mellitus/gdm_qsp_model.dot) | Mechanistic map source (Graphviz) — 208 nodes · 14 clusters · 316 edges |
| [`gdm_qsp_model.svg`](../../../gestational-diabetes-mellitus/gdm_qsp_model.svg) | Vector map (zoomable, searchable) |
| [`gdm_qsp_model.png`](../../../gestational-diabetes-mellitus/gdm_qsp_model.png) | Raster map (150 dpi) |
| [`gdm_mrgsolve_model.R`](gdm_mrgsolve_model.R) | **31-ODE** mrgsolve model + 10 scenarios + calibration rationale |
| [`gdm_shiny_app.R`](gdm_shiny_app.R) | **10-tab** Shiny dashboard |
| [`gdm_references.md`](gdm_references.md) | **87** references (every PMID machine-verified) |

### The map's 14 clusters

1. Pre-pregnancy susceptibility (the determinants of `BCAP`)
2. Placental endocrine unit — the counter-regulatory generator
3. Maternal insulin-resistance machinery (post-receptor defects)
4. Beta-cell adaptation — compensation vs failure
5. Maternal glucose kinetics · screening/diagnosis
6. Placental transport — nutrients **and drugs**
7. Fetal compartment — the Pedersen cascade
8. Medical nutrition therapy · exercise
9. Insulin pharmacology
10. Metformin pharmacology (placental transfer)
11. Sulfonylurea pharmacology
12. Maternal clinical endpoints
13. Neonatal · lifetime child endpoints
14. QSP model state variables · observations · calibration targets

---

## mrgsolve Model (31 ODEs)

### State variables

| Block | Compartments |
|---|---|
| Placenta (7) | `PLAC` `HPL` `PROG` `E2` `TNFA` `ADIPO` `CLP` |
| Maternal metabolism (8) | `BCM` `GLU` `INS` `FFA` `MBG` `GUT1` `GUT2` `FATM` |
| Fetus (6) | `GLUF` `INSF` `BCF` `IGF1F` `LEANF` `FATF` |
| Metformin (4) | `MGUT` `MCEN` `MPER` `MFET` |
| Insulin (2) | `IBAS` `IBOL` |
| Glyburide (3) | `GGUT` `GCEN` `GFET` |
| Safety (1) | `HYPOAUC` (cumulative maternal hypoglycaemia exposure) |

The time axis is **gestational day** (GW 8 → 12 weeks postpartum). At
delivery (`TDEL`) the placenta regresses, the entire fetal block is
frozen, and the counter-regulatory drive disappears.

### Three key equations

**Counter-regulatory index → insulin sensitivity**
```
CID   = SHPL·hPL + SPROG·(PROG+CLP)/100 + STNF·TNFα + SCORT·PLAC
SI    = SIREF · 1/(1+CID) · exp(−KBMI·(BMI−22)) · adiponectin term · (1+exercise+metformin)
```
CID ≈ 0.22 at GW 8, ≈ 1.48 at GW 40 → SI falls by about 49% (consistent
with Catalano's observations).

**Beta-cell compensation — two axes, mass and sensitivity**
```
dBCM/dt = KPROL·BCAP·[hPL/(hPL+KMH)]·BCM·(1−BCM/BCMMAX) − loss term (glucotoxicity, lipotoxicity, SU load)
KG50E   = KG50·(1 − KSENS·BCAP·[hPL/(hPL+KMH)])      ← left-shift of the secretion curve
```
Because human beta-cell mass increases only about 1.4-fold (Butler 2010),
a substantial part of the compensation is separated out as **function**
(the left-shift).

**Fetus — fat is the elastic compartment, lean mass is inelastic**
```
dLEANF/dt = KLN·IGF1F·LEANF·(1−LEANF/LNMAX)                        AIGF = 0.20
dFATF /dt = KFT·PLAC³·(1 + AFAT·(INSF/INSFREF−1) + BFAT·(FFA/FFAREF−1))·(LEANF/1000)
                                                                    AFAT = 0.55
```
For each 1 mg/dL rise in mean maternal glucose, birth weight rises by
about +16 g, and about two-thirds of that increment is **fat.** This is
why overgrowth in GDM is asymmetric, and why the risk of shoulder dystocia
differs between GDM and non-GDM neonates of the same birth weight
(Catalano 2003).

### Calibration targets

| Target | Observed | What corresponds to it in the model |
|---|---|---|
| HAPO continuous gradient | Category 1 LGA 5.3% → category 7 26.3% | `P(LGA) = Φ((BWZ−1.282)/1.0)`; a z of about 0.9 is needed per 25 mg/dL of fasting glucose |
| MFMU treatment effect | LGA 14.5→7.1, shoulder dystocia 4.0→1.5, pre-eclampsia 13.6→8.6, caesarean 33.8→26.9 (%) | `PSD0/KSD`, `PCS0/KCS`, `PPE0/KPEG/KPEB` were solved analytically from these pairs |
| MiG | 46% of the metformin arm needed supplemental insulin | Scenario S8 is built so this is **predicted** (not specified) |
| Bellamy | Postpartum T2DM RR 7.43 | `KDI = 3.8`; RR ≈ 8 at DI/DIREF = 0.45 |
| Placental drug transfer | Metformin U:M ≈ 1.0 / glyburide cord:M ≈ 0.7 / insulin 0 | `KPLM/(KPLM+CLMF) = 0.94` · `KPLG/(KPLG+CLGF) = 0.70` · insulin has **no transfer term at all** |

The full derivation (unit conversions, hand-solved fixed points, the
source of each coefficient) is in the `CALIBRATION` section, parts A–K,
of `gdm_mrgsolve_model.R`.

### 10 Scenarios

| # | Scenario | What it is meant to show |
|---|---|---|
| S1 | Normal-pregnancy baseline | The normal trajectory in which fasting glucose **falls** |
| S2 | Untreated GDM · secretory-deficient type | Fasting normal, abnormal after a load — the failure of fasting-only screening |
| S3 | Untreated GDM · insulin-resistant type | Hyperglycaemia + hyperinsulinaemia + elevated FFA + a contribution to pre-eclampsia |
| S4 | MNT + exercise | Carbohydrate amount/distribution + the insulin-independent (AMPK) pathway |
| S5 | MNT + metformin | Effect and **fetal exposure** shown on the same screen |
| S6 | MNT + insulin (auto-titrated) | Dose as a titration outcome, not a specification |
| S7 | MNT + glyburide | Both maternal and fetal consequences of glucose-**independent** secretion |
| S8 | Metformin → insulin add-on | Reproducing MiG's 46% as a prediction |
| S9 | HAPO gradient sweep | Threshold-free continuity — the core validation of this model |
| S10 | Postpartum follow-up · 5-year T2DM | What remains once the placenta is gone |

### Running it

```r
# Compile the model + 8 treatment scenarios + endpoint table
source("gdm_mrgsolve_model.R")
sims <- run_all()
endpoint_table(sims)

# HAPO validation sweep (BCAP x BMI grid)
sweep <- hapo_sweep()
plot_hapo(sweep)

# Insulin requirement is determined by titration (not specified)
titrate_insulin(list(BCAP = 0.55, BMI = 31, EXEFF = 0.20), fpg_target = 95)
```

---

## Shiny Dashboard (10 Tabs)

```r
shiny::runApp("gdm_shiny_app.R")
```

| Tab | Contents |
|---|---|
| 1 | Patient profile · **Powe subtype coordinates** (current position on the BCAP x sensitivity plane) |
| 2 | Placental endocrine axis (a clock that runs nearly identically in every pregnancy) |
| 3 | Maternal glucose · insulin + a 48-hour CGM-like profile + a TIR table |
| 4 | Beta-cell compensation · disposition index · left-shift of the secretion curve |
| 5 | Drug PK and **fetal exposure** (maternal vs fetal concentrations, checked against reported literature ratios) |
| 6 | Fetal growth · fat vs lean mass · the Pedersen cascade |
| 7 | Clinical endpoints **shown alongside the MFMU observed values** |
| 8 | Scenario comparison (same patient physiology, different strategies) |
| 9 | Biomarkers · postpartum trajectory · 5-year T2DM |
| 10 | HAPO validation (the 7 observed HAPO categories overlaid on the model sweep) |

Showing model output and actual trial results **overlaid on the same
axes** in tabs 7 and 10 is deliberate. A dashboard that shows only its own
output invites mistaking calibration for validation.

---

## What This Model Cannot Do (Limitations — Read This Part)

1. **It has not been fitted.** There is no individual patient data, and
   the parameter-uncertainty block (`$OMEGA`) is deliberately absent. Every
   scenario is a deterministic single individual.
2. **The hormones have been lumped into one.** hPL, progesterone,
   cortisol, and TNF-α cannot be individually identified from clinical
   data, so they are combined into a single `CID`. Claims about the
   contribution of any individual hormone therefore cannot be tested with
   this model.
3. **Hyperglycaemia-induced GLUT1 downregulation is not implemented**
   (Hahn 1998). Fetal exposure may be over-predicted under extreme
   hyperglycaemia.
4. **Fetal metformin is tracked, but its growth effect is fixed at zero.**
   In MiG-TOFU, metformin-exposed infants had more fat, but no intrauterine
   mechanism contained in this model can predict that result. Turning on
   such a term would mean manipulating a mechanism that does not exist, so
   it is left visible in the code multiplied by zero.
5. **Endpoint probabilities are population regressions laid on top of the
   mechanistic drivers** — they are not a mechanistic model of labour
   dynamics, placentation, or neonatal adaptation itself.
6. **Pre-eclampsia is handled only by regression.** sFlt-1/PlGF are drawn
   on the map but have no differential equation.
7. Stillbirth, congenital malformation (more a problem of pre-existing
   diabetes than of GDM), and neonatal respiratory outcomes are on the map
   but not in the equations.
8. Because postprandial absorption is approximated with a two-compartment
   transit model and incretin effect is a constant amplification, the
   model cannot reproduce the difference between oral and intravenous
   glucose loading.

---

## Reproducing the Map

```bash
dot -Tsvg gdm_qsp_model.dot -o gdm_qsp_model.svg
dot -Tpng -Gdpi=150 gdm_qsp_model.dot -o gdm_qsp_model.png
```

---

## Disclaimer

A semi-quantitative QSP model for educational and research purposes. It
was built on published literature but has not been independently
validated or certified, and has not been fitted to actual patient data.
**It must not be used for clinical decision-making, prescribing, or
regulatory submission.**
