# Chronic Hypoparathyroidism (HypoPT) — QSP Model
### Post-surgical HypoPT as the index phenotype · ADH1 (activating CASR) as the contrast case

<p align="center">
  <a href="hypopt_qsp_model_en.svg">
    <img src="hypopt_qsp_model_en.png" width="900" alt="HypoPT QSP mechanistic map">
  </a><br>
  <sub>Click to open as a zoomable SVG · 202 nodes · 276 edges · 16 clusters</sub>
</p>

---

## 1. What this model explains

Extracellular ionised calcium is held within about a 5% band by a single fast
negative-feedback loop. The parathyroid chief cell reads ionised calcium
through **CaSR** and secretes PTH along an inverse sigmoid curve. PTH does
exactly three things.

1. Switches on distal-tubule calcium reabsorption (TRPV5 · calbindin-D28k · NCX1)
2. Switches on renal CYP27B1 to raise 1,25(OH)₂D and increase intestinal absorption
3. Raises bone remodelling so the skeletal reservoir can be drawn on

A fourth action rides along with these — internalising NPT2a to **waste
phosphate**. When PTH disappears, all four disappear **simultaneously**.

### The decisive point of this disease is not "calcium is low"

The only regulatory lever left is the **gut**, and that lever moves only in
response to exogenous calcitriol and oral calcium — **neither of which is
sensed by feedback, and neither of which restores the renal axis**. In other
words, conventional treatment raises serum calcium by **pushing more calcium
into a kidney that has lost the ability to reclaim it**. The result is the
core counter-intuitive behaviour this model is built to reproduce.

> **The filtered load rises, but the fraction reabsorbed does not. So urinary
> calcium rises steeply in step with serum calcium. The patient buys
> normocalcaemia at the cost of a 2–17-fold risk of hypercalciuria ·
> nephrocalcinosis · CKD. At the same serum calcium, the only intervention
> that lowers urinary calcium is PTH replacement — because it restores the
> renal axis rather than overwhelming it.**

Because of this mechanism, the primary endpoint in trials in this field is not
serum calcium alone but a **triple composite response** (normalised corrected
calcium **AND** oral calcium ≤600 mg/day **AND** discontinuation of active
vitamin D). A serum-calcium target alone can easily look like "success" while
the kidney is being damaged. `RESP3` in the mrgsolve model is that composite
indicator.

---

## 2. Deliverables

| File | Contents |
|------|------|
| [`hypopt_qsp_model_en.dot`](hypopt_qsp_model_en.dot) | Mechanistic map source (202 nodes · 276 edges · 16 clusters) |
| [`hypopt_qsp_model_en.svg`](hypopt_qsp_model_en.svg) | Vector rendering (zoomable · searchable) |
| [`hypopt_qsp_model_en.png`](hypopt_qsp_model_en.png) | Raster rendering (150 dpi) |
| [`hypopt_mrgsolve_model_en.R`](hypopt_mrgsolve_model_en.R) | 36-compartment ODE model + 15 scenarios + dose-titration helper |
| [`hypopt_shiny_app_en.R`](hypopt_shiny_app_en.R) | 9-tab interactive dashboard |
| [`hypopt_references_en.md`](hypopt_references_en.md) | 96 PubMed links + parameter-literature cross-reference table |

---

## 3. Mechanistic map — 16 clusters

| # | Cluster | Key content |
|---|----------|-----------|
| 1 | Aetiology | Post-surgical (75–80%) · APS-1/AIRE · anti-NALP5 · CaSR-activating autoantibodies · **ADH1 (CASR)** · ADH2 (GNA11) · GCM2 · PTH gene · HDR (GATA3) · DiGeorge (TBX1) · Kenny–Caffey (FAM111A) · mitochondrial disease · infiltration · radiation · **hypomagnesaemia (reversible)** |
| 2 | Parathyroid chief cell | CaSR–Gq/11–PLCβ–IP₃ · Gi–cAMP↓ · **set point and sigmoid curve** · PTH transcription→prepro→pro→PTH(1-84) · secretory granules · PTH(7-84) antagonism · VDR/FGF23 feedback · Mg-dependent Gα |
| 3 | PTH1R signalling | R0 (sustained) vs RG (transient) states · Gs–AC–cAMP–PKA · Gq–PKC · β-arrestin internalisation · **sustained endosomal cAMP signalling** · desensitisation/recycling |
| 4 | Kidney | Filtered calcium load · proximal 65% · TAL claudin-16/19 and **TAL CaSR** · DCT TRPV5/calbindin/NCX1 (PTH-dependent) · NPT2a → **TmP/GFR↑** · CYP27B1↓ · CYP24A1↑ · hypercalciuria → nephrocalcinosis → CKD vicious cycle |
| 5 | Bone | Osteocyte–sclerostin–Wnt · RANKL/OPG↓ · **low-turnover bone (BMU activation frequency↓)** · exchangeable calcium reservoir · P1NP/CTX/BSAP↓ · BMD Z +1~+2 · microstructural heterogeneity · site-specific heterogeneous fracture risk |
| 6 | Gut | Calcium carbonate (gastric pH-dependent) vs calcium citrate · TRPV6/calbindin-D9k/PMCA1b · paracellular diffusion · **dose-dependent saturation of absorption fraction** · NPT2b · phosphate binders · intestinal oxalate |
| 7 | Vitamin D | 7-DHC→cholecalciferol→CYP2R1→25(OH)D→**CYP27B1**→1,25D · DBP binding · CYP24A1 autocatabolism · calcitriol/alfacalcidol/DHT |
| 8 | Phosphate · FGF23 | Hyperphosphataemia · FGF23–αKlotho · Ca×P product · pyrophosphate/fetuin-A saturation · **ectopic calcification (basal ganglia · vessels · lens)** · VSMC osteogenic transition |
| 9 | Magnesium | Intake/TRPM6 · PPI · diuretic losses · **simultaneous secretory blockade (Gα) and end-organ resistance (AC)** · reversible with repletion |
| 10 | Calcium species | Total calcium = albumin-bound 40% + complexed 10% + **ionised 50%** · pH effect · albumin correction formula · citrate |
| 11 | Neuromuscular | Na⁺ channel threshold instability → hyperexcitability↑ → paraesthesia · tetany · Chvostek/Trousseau · laryngospasm · seizure · QTc prolongation · reversible cardiomyopathy |
| 12 | Chronic complications | Basal ganglia calcification (Fahr-type) · brain fog · anxiety/depression (2-fold) · cataract · enamel hypoplasia · **HPES quality of life** |
| 13 | Conventional treatment | Calcium carbonate/citrate · calcitriol · alfacalcidol · cholecalciferol · **thiazide (the only conventional agent acting on urinary Ca)** · phosphate binders · Mg · IV calcium · **limitation: urinary calcium not controlled** |
| 14 | PTH replacement therapy | rhPTH(1-84) **biphasic SC absorption** · teriparatide BID/pump · **palopegteriparatide (TransCon PTH) sustained-release prodrug** · eneboparatide (RG-selective) · LA-PTH · tapering of conventional treatment |
| 15 | CaSR-targeted therapy | **encaleret (CaSR NAM)** → rightward shift of set point → recovery of residual PTH + urinary Ca↓ · autotransplantation · iPSC cell therapy |
| 16 | Endpoints | Corrected calcium 8.0–9.0 · 24h urinary Ca <300 mg · TmP/GFR · eGFR slope · nephrocalcinosis · **triple composite response** · HPES · QTc · BMD Z |

---

## 4. mrgsolve model — 36 compartments

```
Calcium       GUTCA (intestinal lumen) · CAE (extracellular) · CABONE (exchangeable bone) · CAMIN (structural bone)
Phosphate     PIE · PIBONE
Magnesium     MGE
PTH           PTHE (endogenous) · PTHD1/PTHD2 (SC biphasic depot) · PTHX (exogenous central)
TransCon      TCDEP (SC prodrug) · TCPRO (circulating prodrug)
Vitamin D     CTRDEP · ALFDEP · ALFC · D125 · D3DEP · D25
Other drugs   THZDEP/THZC (thiazide) · ENCDEP/ENCC (encaleret)
Signalling    PSIG (sustained PTH1R signalling — endosomal cAMP)
Cells/markers OC · OB · PKS (anabolic peak memory) · FGF23
End-organ     NC (nephrocalcinosis) · GFRC · BMDZ · BGCALC
Symptoms      CAEFF (effect compartment) · QOL
Cumulative    UCACUM · UPICUM
```

### Five behaviours the model was designed to reproduce

1. **Normalising serum calcium does not normalise urinary calcium** — scenarios 3 vs 8
2. **Only PTH replacement lowers urinary calcium at the same serum calcium** — scenarios 8·9·15
3. **Pulsatile and flat exposure act differently on bone** — osteoclasts follow
   time-averaged PTH while osteoblasts follow peak memory (`PKS`), so
   teriparatide is anabolic while TransCon PTH is closer to normalising
   turnover — scenarios 7 vs 8
4. **Calcium · calcitriol are especially dangerous in ADH1** — because the
   renal tubular CaSR is simultaneously overactive, urinary calcium is higher
   at the same serum calcium. Encaleret, by resetting the set point, is the
   only intervention that produces serum calcium↑ and urinary calcium↓
   **simultaneously** — scenarios 10 vs 11
5. **Hypomagnesaemia blocks both secretion and action at once** — it is
   refractory to calcium · calcitriol and recovers sharply with magnesium
   repletion — scenario 12

### 15 scenarios

| # | Scenario | Key point |
|---|----------|------|
| 1 | Normal control | Baseline steady state (Ca 9.5 · PTH 35 · 1,25D 40 · urinary Ca 150 mg/day) |
| 2 | Untreated post-surgical HypoPT | Ca 6.3–6.8 · phosphate 4.7 · 1,25D ~17 · bone turnover halved · BMD Z +1.5 |
| 3 | Conventional treatment (Ca 1000 mg + calcitriol 0.5 µg) | Reaches serum target, urinary Ca 229 mg/day |
| 4 | Conventional treatment + thiazide (calcium-sparing) | Same serum calcium at half the calcium dose · urinary Ca 149 mg/day |
| 5 | Overtreatment (Ca 3000 mg + calcitriol 1.0 µg) | Ca 11.3 · Ca×P >50 · urinary Ca 496 mg/day |
| 6 | rhPTH(1-84) 50→100 µg QD + taper | Biphasic absorption → PTH peak/trough ratio 29-fold |
| 7 | Teriparatide 30 µg BID | Short, high peak (ratio 79-fold) → anabolic bone response |
| 8 | Palopegteriparatide 18 µg QD | Flat exposure (ratio 1.0) · full discontinuation of conventional treatment · phosphate normalised |
| 9 | Palopegteriparatide 30 µg QD | Uptitrated · achieves triple composite response |
| 10 | Untreated ADH1 | PTH "inappropriately normal" at low calcium + high urinary Ca |
| 11 | ADH1 + encaleret 60 mg BID | Set point restored → calcium↑ urinary Ca↓ simultaneously |
| 12 | Hypomagnesaemic functional HypoPT | Mg 0.78 → Ca 6.98 · sharp recovery to Ca 9.6 after oral Mg on day 30 |
| 13 | Acute hypocalcaemic crisis | Ca falls 9.5→7.3 over 72 hours, rescued with IV calcium after 27.5 hours of tetany |
| 14 | 10-year conventional treatment | Nephrocalcinosis accumulates to 0.32 → eGFR falls 1.31 mL/min/year |
| 15 | 10-year palopegteriparatide | Same patient, nephrocalcinosis 0 · eGFR 0.80 mL/min/year |

### Usage

```r
source("hypopt_mrgsolve_model_en.R")
mod <- HYPOPT_build()

# Short-term scenarios
sim <- HYPOPT_simulate_scenarios(mod, which = 1:13)
HYPOPT_summary(sim)                    # steady-state summary table
HYPOPT_swing(sim, from_day = 150)      # calcium swing within a dosing interval

# Long-term renal outcomes
long <- HYPOPT_simulate_scenarios(mod, which = 14:15)
HYPOPT_summary(long)

# Conventional-treatment dose grid — the renal cost of hitting the calcium target
HYPOPT_titrate_conventional(mod)
```

The disease state is not "written in" as an initial value — it is where the
system settles after the parameters are changed. `hypopt_baseline()` integrates
forward 400 days with no treatment to reach a state where the derivatives
vanish, and passes that state vector as the `init` for the treatment
simulations.

---

## 5. Shiny app — 9 tabs

| Tab | Contents |
|----|------|
| 1. Patient profile | Aetiology presets · chief-cell function · CaSR set point · before/after treatment lab comparison · **Ca–PTH set-point curve** |
| 2. PK | Exposure by PTH molecule (full course + last-7-day zoom) · calcitriol · drug concentrations |
| 3. Calcium · phosphate homeostasis | Corrected calcium · ionised calcium · phosphate · magnesium · 1,25D · FGF23 (with normal-range bands) |
| 4. Renal safety | **24h urinary calcium** · FE_Ca · TmP/GFR · Ca×P · nephrocalcinosis · eGFR |
| 5. Bone | Bone turnover index · osteoblasts/osteoclasts · P1NP·CTX·BSAP · BMD Z-score |
| 6. Clinical endpoints | Symptom score · QTc · QoL · basal ganglia calcification · **triple composite response verdict** |
| 7. Scenario comparison | Compare selected scenarios among the 15 + summary table + calcium-swing table |
| 8. Biomarker summary | Steady-state table of 29 metrics + **conventional-treatment dose grid** |
| 9. References · Help | Rendered references |

```r
shiny::runApp("hypopt_shiny_app_en.R")
```

---

## 6. Verified model output

The table below is obtained by actually running the code in this repository
(`HYPOPT_simulate_scenarios()` → `HYPOPT_summary()`, mrgsolve 2.0.1, R 4.3.3).

| # | Ca (mg/dL) | Phosphate | PTH | 24h urinary Ca | FE_Ca (%) | Bone turnover | BMD Z | Symptoms | Composite response |
|---|-----------|----|-----|------------|-----------|--------|-------|------|---------|
| 01 Normal | 9.48 | 3.50 | 35.2 | 147 | 1.50 | 1.00 | 0.01 | 7 | — |
| 02 Untreated | 6.03 | 4.45 | 3.3 | 51 | 0.82 | 0.44 | 1.28 | 95 | 0% |
| 03 Conventional treatment | 8.66 | 4.48 | 2.2 | **229** | 2.55 | 0.41 | 1.31 | 22 | 0% |
| 04 +thiazide | 8.64 | 4.48 | 2.2 | **149** | 1.67 | 0.41 | 1.31 | 23 | 0% |
| 05 Overtreatment | 11.25 | 4.53 | 1.0 | 496 | 4.27 | 0.38 | 1.35 | 1 | 0% |
| 06 rhPTH(1-84) | 7.82 | 3.61 | 25.6 | 58 | 0.72 | 1.07 | 0.96 | 51 | 32% |
| 07 Teriparatide | 8.16 | 3.09 | 51.4 | 43 | 0.51 | 1.52 | 0.71 | 38 | 100% |
| 08 TransCon 18 µg | 8.53 | 3.83 | 19.5 | **125** | 1.42 | 0.79 | 0.93 | 26 | 96% |
| 09 TransCon 30 µg | 9.26 | 3.58 | 30.5 | 143 | 1.49 | 0.95 | 0.77 | 10 | 100% |
| 10 Untreated ADH1 | 7.81 | 3.86 | 18.7 | 135 | **1.68** | 0.78 | 0.51 | 52 | 2% |
| 11 ADH1+encaleret | **8.83** | 3.77 | 21.9 | **131** | **1.44** | 0.83 | 0.45 | 18 | 100% |
| 12 Low Mg → repletion | 9.62 | 3.45 | 36.1 | 151 | 1.52 | 1.03 | 0.73 | 6 | 73% |
| 14 10-year conventional | 8.67 | 4.52 | 2.3 | 207 | 2.55 | 0.41 | 1.75 | 19 | 0% |
| 15 10-year TransCon | 8.98 | 3.72 | 24.4 | 132 | 1.50 | 0.87 | 0.38 | 15 | 100% |

Three ways to read this:

- **3 vs 8** — at essentially the same serum calcium (8.66 vs 8.53), urinary
  calcium goes from 229 → 125 mg/day. It is restoration of the renal axis, not
  conventional treatment, that protects the kidney.
- **10 vs 11** — encaleret raises serum calcium 7.81 → 8.83 while
  **simultaneously** lowering FE_Ca 1.68 → 1.44%. No other intervention moves
  both directions together.
- **14 vs 15** — 10 years in the same patient. eGFR decline 1.31 vs 0.80
  mL/min/year, nephrocalcinosis 0.32 vs 0.

The dose-titration grid (`HYPOPT_titrate_conventional()`) shows this trade-off
directly. Every combination of conventional treatment that lands within the
8.0–9.0 calcium target carries 24-hour urinary calcium of 191–249 mg — there
is no way to avoid this trade-off with conventional treatment.

---

## 7. Calibration anchors

| State | Model target value |
|------|-------------|
| Healthy adult | Ca 9.5 mg/dL · iCa 1.19 mmol/L · PTH 35 pg/mL · 1,25D 40 pg/mL · phosphate 3.5 mg/dL · TmP/GFR 3.1 · urinary Ca 150 mg/day · intestinal absorption fraction ~30% |
| Untreated HypoPT | Ca 6.3–6.8 · PTH <5 · 1,25D ~17 · phosphate 4.6–4.8 · bone turnover markers 40–55% of normal · BMD Z +1.3~+1.6 |
| Conventional treatment | Ca ~8.6 (oral Ca 1000 mg + calcitriol 0.5 µg/day) · 24h urinary Ca ~240 mg |
| PTH replacement | 24h urinary Ca ~150 mg at the same serum calcium · phosphate normalised · oral calcium discontinued |

All anchors are reproduced by the run results in Section 6 above.

For the literature basis of each anchor and the parameter cross-reference
table, see [`hypopt_references_en.md`](hypopt_references_en.md).

---

## ⚠️ Disclaimer

This model is a **qualitative · semi-quantitative QSP model for educational and
research purposes**. It was built from published literature and clinical
trial data but has not been independently validated · certified, and
**must not be used directly for actual clinical decision-making, prescribing,
or regulatory submission.**
