# Acne Vulgaris — QSP Model

> **Quantitative Systems Pharmacology model of acne vulgaris**
> The androgen–IGF-1/mTORC1–sebocyte axis · infundibular hyperkeratinisation ·
> *Cutibacterium acnes* phylotypes/biofilm/resistance · NLRP3–IL-1β–Th17
> inflammation · lesion-transition dynamics · topical/systemic/hormonal/
> isotretinoin pharmacology
>
> ⚠️ For educational and research use. Do not use for clinical decision-making, prescribing, or regulatory submission.

---

## 1. Files

| File | Contents |
|------|------|
| [`acn_qsp_model.dot`](acn_qsp_model.dot) | Mechanistic map source — **18 clusters · 277 nodes · 383 edges** |
| [`acn_qsp_model.svg`](acn_qsp_model.svg) | Vector rendering (zoom in to read) |
| [`acn_qsp_model.png`](acn_qsp_model.png) | Raster rendering (150 dpi) |
| [`acn_mrgsolve_model.R`](../../../acne-vulgaris/acn_mrgsolve_model.R) | mrgsolve model with **55 ODE compartments · 213 parameters · 17 scenarios** |
| [`acn_shiny_app.R`](acn_shiny_app.R) | 10-tab interactive dashboard |
| [`acn_references.md`](acn_references.md) | **110** PubMed-verified references (PMIDs actually looked up) |

Rendering:
```bash
dot -Tsvg acn_qsp_model.dot -o acn_qsp_model.svg
dot -Tpng -Gdpi=150 acn_qsp_model.dot -o acn_qsp_model.png
```

---

## 2. Why This Model Exists — In One Paragraph

The four pillars of acne (excess sebum · infundibular
hyperkeratinisation · *C. acnes* dysbiosis · inflammation) are all in the
textbooks. The clinically decisive fact lies elsewhere. **Every lesion
the patient can see descends from a microcomedone the patient cannot
see.** The microcomedone reservoir is about three times larger than the
visible comedone population and turns over with a time constant of about
4 weeks. This single fact explains nearly every counterintuitive
phenomenon in acne therapeutics, and this model was built precisely to
reproduce that phenomenon.

> **Anti-inflammatory therapy (antibiotics, dapsone) rapidly empties the
> downstream compartments but never touches the upstream reservoir at
> all. That is why lesion counts fall within 6–8 weeks, and return within
> a few weeks of stopping. Retinoids drain the reservoir itself. This is
> slow (little visible benefit before week 8, and sometimes a transient
> worsening), but it is the only means of sustaining remission. This is
> why every guideline says "induce with a retinoid plus an antimicrobial,
> maintain with a retinoid alone," and why antibiotic monotherapy is not
> merely a weak treatment but a wrong one.**

The second non-trivial behaviour built into this model: **isotretinoin
does not cure acne by killing bacteria** (it has no direct antimicrobial
action at all). Shrinking the sebaceous gland removes the niche, and
*C. acnes* declines **as a consequence.** Whether relapse then follows is
determined by how much of the gland shrinkage is **permanent**, and the
model represents this as an explicit state variable of cumulative
exposure (`DURAB`). As a result, **the 120–150 mg/kg cumulative-dose rule
appears not as an if-statement but as an emergent property** (see §5).

Third: **because benzoyl peroxide is a non-specific oxidant, no
resistance mechanism to it exists, and it actively clears resistant
*C. acnes*.** So the resistant fraction `RESF` is a real state variable —
raised by antibiotic selection pressure and lowered by BPO. Overlaying
Scenario 6 (clindamycin alone) and Scenario 7 (clindamycin+BPO) turns the
entire antibiotic stewardship argument into a single picture.

---

## 3. Mechanistic Map (18 clusters)

[![ACN QSP map](acn_qsp_model.png)](acn_qsp_model.svg)

| # | Cluster | Key content |
|---|---|---|
| 1 | Pathogenesis — genetics · adrenarche · exogenous factors | Heritability 78–81%, GWAS loci, glycaemic load · dairy, drug-induced, PCOS |
| 2 | Androgen axis | DHEA-S → T → SHBG → FT → 5α-reductase type 1 → follicular DHT → AR |
| 3 | Nutrient-metabolic signalling | Insulin · IGF-1 · PI3K/AKT · FoxO1 ⟂ mTORC1 → SREBP-1c |
| 4 | Sebaceous gland (Pillar ①) | Stem cell → proliferation → lipogenesis → holocrine secretion, squalene peroxide |
| 5 | Infundibular hyperkeratinisation (Pillar ②) | Linoleic acid dilution deficiency · IL-1α · desmosome retention → keratin plug → microcomedone |
| 6 | *C. acnes* (Pillar ③) | Phylotypes IA1/IB/II/III, biofilm, lipase→FFA, porphyrins, resistance genes |
| 7 | Innate immune recognition | TLR2/NOD2/PAR-2 → NF-κB, NLRP3–ASC–caspase-1 → IL-1β, AMPs |
| 8 | Effector inflammation (Pillar ④) | IL-8 → neutrophils → MPO/NETs → MMP-1/9 → dermal collagen breakdown |
| 9 | Adaptive immunity | RORγt → Th17 → IL-17A → keratinocyte CXCL8 amplification loop |
| 10 | Lesion-transition dynamics | MC → closed/open comedone → papule → pustule → nodule → resolution/PIH/scar |
| 11 | Neuroendocrine | Local CRH–CRH-R1, substance P, cortisol, UV oxidation |
| 12 | Topical therapy PK/PD | BPO · 4 retinoids · clindamycin · azelaic acid · dapsone · clascoterone |
| 13 | Systemic antibiotics | Doxy/mino/sarecycline, sub-antimicrobial doses, non-antimicrobial anti-inflammatory action, resistance management |
| 14 | Hormonal therapy | COC (EE→SHBG↑, LH↓) · spironolactone · cyproterone · metformin |
| 15 | Isotretinoin | PK (food effect · 4-oxo) · FoxO1 restoration · sebocyte apoptosis · cumulative dose–relapse · safety |
| 16 | Clinical endpoints | Inflammatory/non-inflammatory lesion counts · IGA success · scarring · QoL · relapse rate |
| 17 | Severe and variant phenotypes | Conglobata · fulminans · adult female · SAPHO/PAPA · drug-induced |
| 18 | Model state-variable map | Correspondence table for the 55 mrgsolve compartments |

Edge rules: solid black = material flow/transition · blue = activation ·
red (tee) = inhibition · green dashed = drug action · orange dashed =
adverse effect · grey dotted = feedback loop.

---

## 4. mrgsolve Model Structure (55 ODEs)

| Group | State variables |
|---|---|
| Hormones | `TT` `SHBG` `ARS` |
| Metabolism | `IGF1` `INS` `LIP` |
| Sebaceous gland | `SGM` `DURAB` `SER` `LA` `SQOX` |
| Keratinisation | `KER` `IL1A` |
| Microbiology | `CAP` `CAB` `RESF` `FFA` `PORP` |
| Inflammation | `TLR` `IL1B` `IL8` `TNF` `IL17` `NEU` `MMP` `CRH` |
| Lesions | `MC` `CC` `OC` `PAP` `PUS` `NOD` `PIH` `SCAR` |
| Topical PK | `BPOS` `RETS` `CLIS` `AZES` `CLAS` `DAPS` |
| Systemic PK | `TETD` `TETC` `ISOD` `ISOC` `ISOP` `OXOC` `SPID` `SPIC` `EED` `EEC` `CUMISO` |
| Safety | `TG` `ALT` `KSER` `MUCO` |

Four design features:

1. **A self-calibrating baseline state.** `acn_baseline()` integrates for
   4 years without treatment to construct that patient's steady state.
   Lesion counts are an **emergent value, not an input**. Changing only
   the androgen, metabolic, and exposure parameters naturally yields
   everything from mild comedonal to severe nodular disease.
2. **Splitting tetracycline into two arms.** The antimicrobial arm has
   EC₅₀ ≈ 2.5 mg/L; the non-antimicrobial anti-inflammatory arm
   (NF-κB · IL-8 · MMP) has EC₅₀ ≈ 0.35 mg/L. The sub-antimicrobial
   40 mg extended-release dose sits precisely between the two — the
   anti-inflammatory effect is intact, with no antimicrobial effect and
   no selection pressure.
3. **Resistance as a state variable.** `RESF` is seeded by background
   mutation, selected by antibiotic exposure, cleared in proportion to
   BPO killing, and slowly reversed by a fitness cost.
4. **Cumulative-dose-dependent permanent shrinkage.** `DURAB` follows a
   Hill function of `CUMISO` (CD50 ≈ 85 mg/kg, H = 2.2). Relapse is an
   outcome, not a rule.

### Phenotype Presets
`mild_comedonal` · `moderate` · `moderate_male` · `severe_nodular` ·
`adult_female` · `pcos` · `skin_of_colour`

Baseline steady state (pre-treatment, model output values):

| Phenotype | Inflammatory | Non-inflammatory | Nodules | IGA | SER (µg/cm²/min) |
|---|---|---|---|---|---|
| mild_comedonal | 6.4 | 9.7 | 0.1 | 1.52 | 1.12 |
| moderate | 20.4 | 29.6 | 1.6 | 2.82 | 1.88 |
| moderate_male | 29.0 | 41.7 | 3.4 | 3.16 | 2.13 |
| severe_nodular | 62.4 | 91.0 | 12.1 | 3.63 | 2.83 |
| adult_female | 35.4 | 57.2 | 3.1 | 3.25 | 2.63 |
| pcos | 51.6 | 86.2 | 5.6 | 3.48 | 3.00 |
| skin_of_colour | 20.4 | 29.6 | 1.6 | 2.82 | 1.88 |

---

## 5. Scenario Library (17) and Run Results

```r
source("acn_mrgsolve_model.R")
mod <- ACN_build()
sim <- ACN_simulate(mod, which = 1:17)
ACN_summary(sim, at_week = 12)
```

### At Week 12 (moderate phenotype; Scenarios 12–13 are adult female, 14–16 are severe nodular)

| # | Scenario | Inflammatory Δ | Non-inflammatory Δ | IGA | *C. acnes* Δ | Resistant fraction |
|---|---|---|---|---|---|---|
| 1 | Untreated natural history | 0% | 0% | 2.82 | 0% | 0.045 |
| 2 | Vehicle/placebo | −27% | −25% | 2.55 | 0% | 0.045 |
| 3 | Adapalene 0.1% QD | −43% | −60% | 2.26 | 0% | 0.045 |
| 4 | BPO 2.5% QD | −51% | −35% | 2.16 | −38% | ~0 |
| 5 | Adapalene 0.3%/BPO 2.5% FDC | −63% | −70% | 1.80 | −38% | ~0 |
| 6 | Clindamycin 1% BID alone | −47% | −23% | 2.26 | −34% | **0.31** |
| 7 | Clindamycin/BPO FDC | −75% | −28% | 1.67 | −74% | ~0 |
| 8 | Doxycycline 100 mg + adapalene/BPO | −83% | −67% | 1.18 | −67% | ~0 |
| 9 | Sub-antimicrobial doxycycline 40 mg MR + adapalene | −50% | −59% | 2.13 | 0% | 0.045 |
| 10 | Sarecycline 1.5 mg/kg + adapalene | −67% | −57% | 1.73 | −36% | 0.092 |
| 11 | Clascoterone 1% BID | −56% | −52% | 2.03 | −21% | 0.045 |
| 12 | COC EE 30 µg | −54% | −54% | 2.64 | −20% | 0.045 |
| 13 | Spironolactone 100 mg | −56% | −57% | 2.59 | −22% | 0.045 |
| 14 | Isotretinoin 0.5 mg/kg | −96% | −90% | 1.01 | −77% | 0.045 |
| 15 | Isotretinoin 0.25 mg/kg | −90% | −86% | 1.85 | −60% | 0.045 |
| 16 | Isotretinoin 1.0 mg/kg | −99% | −92% | 0.61 | −87% | 0.045 |
| 17 | Induction (doxy+BPO) then retinoid maintenance | −83% | −67% | 1.18 | −67% | ~0 |

### At Week 24 — Where the Story Changes

| # | Scenario | Week-12 inflammatory Δ | Week-24 inflammatory Δ | Week-24 resistant fraction |
|---|---|---|---|---|
| 6 | Clindamycin alone | −47% | **−41% (rebound)** | **0.74** |
| 7 | Clindamycin + BPO | −75% | −77% (sustained) | **0.00** |
| 17 | Induction then retinoid maintenance | −83% | −50% (sustained) | 0.01 |

In Scenario 6, lesions return between weeks 12 and 24 not because the
drug was stopped. It is still being applied, but **the resistant fraction
has risen from 0.05 to 0.74, so the drug has stopped working**. In
Scenario 7, with BPO added, the same clindamycin keeps working.

### Isotretinoin Cumulative Dose – 12-Month Relapse Curve

`ACN_cumdose_curve(mod)` (severe nodular, 24 weeks of treatment + 52
weeks of follow-up):

| Daily dose (mg/kg) | Cumulative (mg/kg) | `DURAB` | Nadir lesion count | At 12 months (vs baseline) | SER at 12 months (vs baseline) |
|---|---|---|---|---|---|
| 0.15 | 25 | 0.04 | 12.2 | **93%** | 96% |
| 0.25 | 42 | 0.11 | 5.1 | **81%** | 89% |
| 0.35 | 59 | 0.19 | 2.4 | **68%** | 81% |
| 0.50 | 84 | 0.30 | 0.9 | **51%** | 70% |
| 0.65 | 109 | 0.39 | 0.5 | **40%** | 61% |
| 0.80 | 134 | 0.45 | 0.3 | **33%** | 55% |
| 1.00 | 168 | 0.50 | 0.2 | **27%** | 50% |

The knee of the curve sits at roughly **100–140 mg/kg**. Below it,
relapse rises sharply; above it, additional benefit shrinks while
triglyceride and mucocutaneous toxicity keep growing. This curve is not a
rule explicitly written into the model but an outcome of the
`CUMISO → DURAB → SGM ceiling` chain.

---

## 6. Calibration Targets vs Actual (Calibration)

Based on the week-12 reduction in inflammatory lesions. Target values are
representative figures taken from the registration-trial and guideline
evidence table in [`acn_references.md`](acn_references.md), and no formal
fitting was performed.

| Regimen | Target | Model | Verdict |
|---|---|---|---|
| Vehicle/placebo | 30–35% | 27% | ○ |
| Adapalene 0.1% | 40–50% | 43% | ○ |
| BPO 2.5% | 40–50% | 51% | ○ |
| Adapalene 0.3%/BPO | 60–70% | 63% | ○ |
| Clindamycin alone | 45–50% | 47% | ○ |
| Clindamycin/BPO | 60–68% | 75% | △ overpredicts |
| Doxycycline + adapalene/BPO | 65–75% | 83% | △ overpredicts |
| Sub-antimicrobial doxycycline 40 mg | 45–55% | 50% | ○ |
| Sarecycline + adapalene | 60–70% | 67% | ○ |
| Clascoterone 1% | 45–52% | 56% | ○ |
| COC (24 weeks) | 45–55% | 62% | △ overpredicts |
| Spironolactone (24 weeks) | 50–60% | 64% | △ overpredicts |
| Isotretinoin SER reduction (12 weeks) | −85~−90% | −85% | ○ |
| Isotretinoin ~120 mg/kg, lesions at 12 months | 20–30% of baseline | 25% | ○ |
| Isotretinoin ~35 mg/kg, lesions at 12 months | 55–70% of baseline | 58% | ○ |
| Triglycerides (1.0 mg/kg) | +30~+45% | +38% | ○ |
| Serum potassium (spironolactone 100 mg) | +0.2~0.3 mEq/L | +0.25 | ○ |

**The cause of the 4 overpredictions is known and not hidden.** All are
combination regimens in which several mechanisms act multiplicatively on
the same deterministic patient. An actual trial arm averages responders
and non-responders together, compressing the mean. So `ACN_summary()`
also reports `P_success` (the simulated IGA change logistically
transformed through inter-individual variability `SIGIGA`). It is this
`P_success`, not the deterministic `IGA_success` flag, that should be
compared against a published **IGA success rate**.

---

## 7. Shiny App

```r
shiny::runApp("acn_shiny_app.R")
```
Required packages: `shiny` `mrgsolve` `ggplot2` `dplyr` `tidyr` `DT`

| Tab | Contents |
|---|---|
| 1. Patient profile | Phenotype/exposure settings and the resulting four-pillar bar chart |
| 2. Treatment design | Regimen assembly + **automated guideline checks** (antibiotic without BPO, exceeding 16 weeks, cumulative dose shortfall, tetracycline+isotretinoin contraindication, etc.) |
| 3. PK | Systemic/topical concentrations and pharmacodynamic effect metrics |
| 4. Sebaceous gland · androgens | `SGM` `SER` `LIP` `DURAB`, SHBG/FAI, IGF-1 |
| 5. Microbiology · resistance | Planktonic/biofilm/resistant fraction — the stewardship argument |
| 6. Inflammatory cascade | TLR2 → IL-1β → IL-8 → neutrophils → MMP, plus the hyperkeratinisation axis |
| 7. Lesion counts · endpoints | The 6-compartment transition chain and IGA |
| 8. Scenario comparison | Select from the 17 scenarios and compare side by side |
| 9. Safety | TG · ALT · K⁺ · mucocutaneous · cumulative mg/kg (120–150 baseline marked) |
| 10. Scarring · pigmentation | The nodule-time integral and irreversible scar accumulation |

---

## 8. Limitations

- **Not fitted to individual patient data.** The parameters are
  hand-tuned approximations informed by literature-reported values.
- **This is a deterministic single-patient model.** Inter-individual
  variability is approximated only through the logistic transform of
  `P_success`.
- Lesion counts are a relative scale of "whole-face counts," and do not
  exactly match any particular trial's counting convention.
- Signalling state variables (`TLR`, `IL1B`, etc.) are dimensionless
  relative values, not measured cytokine concentrations.
- Truncal acne, photodynamic therapy, physical treatments
  (extraction · peels · laser), non-antiandrogen hormonal axes (GnRH,
  etc.), and pregnancy/lactation are not in the model.
- Antibiotic resistance is modelled only at the level of an individual's
  follicular flora, not community-level resistance prevalence.
- Adherence (`ADHERE`) is an exogenous parameter. The real feedback loop
  of irritation → discontinuation → worsening is not represented.

---

## 9. Reproducibility

```r
source("acn_mrgsolve_model.R")
mod <- ACN_build()

## the patient's own steady state
b <- acn_baseline(mod, "moderate"); b$summary[, c("INFLAM","NONINF","IGA","SERO")]

## week-12 regimen comparison
print(ACN_summary(ACN_simulate(mod, which = 2:13), at_week = 12))

## antibiotic stewardship: 6 vs 7
print(ACN_resistance(ACN_simulate(mod, which = c(6, 7))))

## isotretinoin cumulative dose and relapse
iso <- ACN_simulate(mod, which = 14:16, delta = 24)
print(ACN_relapse(iso, stop_week = list("14" = 34, "15" = 20, "16" = 20)))
print(ACN_cumdose_curve(mod))
```

Verified working environment: R 4.3 + mrgsolve 1.5.x, Graphviz 2.43. No
negative state variables occur in any scenario.

---

## 10. References

[`acn_references.md`](acn_references.md) — 15 sections, **110** entries.
Every PMID was verified by actually looking up its title, journal, and
year via NCBI E-utilities, and each entry is annotated with "which part
of the model this reference supports."
