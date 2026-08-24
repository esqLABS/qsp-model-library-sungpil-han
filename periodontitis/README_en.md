# Periodontitis — Quantitative Systems Pharmacology Model

> **The bacteria start it, the host does the damage, and only the bone does not come back.**

Severe periodontitis afflicts about 11% of adults worldwide, making it the
sixth most common disease in humanity and the leading cause of adult tooth
loss. This model describes periodontitis not as "bacteria dissolving bone" but
as a **bistable host–microbe feedback loop**, and re-derives the entire logic
of treatment from the fact that the loop's only irreversible output is
alveolar bone.

| Deliverable | File |
|--------|------|
| 🗺️ Mechanistic map (225 nodes · 18 clusters · 369 edges) | [`pdt_qsp_model.dot`](pdt_qsp_model.dot) · [SVG](pdt_qsp_model.svg) · [PNG](pdt_qsp_model.png) |
| ⚙️ mrgsolve ODE model (58 compartments · 205 parameters) | [`pdt_mrgsolve_model.R`](pdt_mrgsolve_model.R) |
| 📊 Shiny dashboard (11 tabs) | [`pdt_shiny_app.R`](pdt_shiny_app.R) |
| 📚 References (93 · 92 PubMed links) | [`pdt_references_en.md`](pdt_references_en.md) |

<p align="center">
  <a href="pdt_qsp_model.svg"><img src="pdt_qsp_model.png" width="820" alt="Periodontitis QSP mechanistic map"></a>
</p>

---

## 1. The organising claims

The common framework in periodontitis modelling — "bacterial load →
inflammation → bone loss" — fails to explain three facts that govern
treatment. So this model starts from a different framework.

### ① The engine is not toxicity but **nutrition** → hence bistability

Inflammation raises gingival crevicular fluid (GCF) flow from roughly 0.06 to
1.2 µL/min per site. GCF is a serum exudate, and red-complex bacteria are
**asaccharolytic, proteolytic anaerobes** — they feed on exactly that albumin
and haem. In other words:

```
Dysbiosis → inflammation → GCF exudation → nutrition → dysbiosis
```

In deep periodontal pockets this loop's loop gain exceeds 1, and the system
has **two attractors**. Whether the model actually has two attractors was
confirmed from different initial conditions.

| | Healthy attractor | Disease attractor |
|---|---:|---:|
| Dysbiotic bacterial load | 0.002 | 2.10 |
| *P. gingivalis* (% of biofilm) | 0.05 | 3.7 |
| Inflammation intensity (healthy = 1) | 0.95 | 4.68 |
| Neutrophil killing capacity (relative to healthy) | 1.00 | 0.19 |
| GCF flow (µL/min) | 0.07 | 1.25 |
| Probing depth (mm) | 1.4 | 7.3 |
| BOP (%) | 14 | 98 |
| PISA (mm²) | 190 | 1239 |

`PDT_separatrix()` uses bisection to numerically locate the **separatrix**. At
a 7.4 mm site, avoiding relapse requires bringing the residual dysbiotic
bacterial load down to **7.6% or less** (about a 1.1 log reduction), yet
closed instrumentation only achieves 89.6% biofilm disruption — the model
judges that site *"will relapse without adjunctive therapy or an open
approach."* This is the quantitative content behind clinical guidelines that
recommend adjunctive therapy or surgery for residual pockets ≥6 mm.

### ② Neutrophils are recruited, then **disarmed** → killing capacity is a state variable

*P. gingivalis* gingipains bypass the C3 convertase and **generate C5a
directly**, and C5aR1 × TLR2 crosstalk degrades MyD88 via Smurf1 while
**selectively preserving** the Mal–PI3K–Akt branch, blocking only
RhoA-dependent phagocytosis. The result is *a neutrophil that cannot kill but
still produces inflammation*.

So in this model, bacterial clearance is always proportional to
`PMN × PKILL`, and **never proportional to `PMN` alone.** At a diseased site,
neutrophils increase 3.6-fold, but killing capacity falls to 0.19-fold, so
**effective clearance actually decreases**.

`PDT_subversion_counterfactual()` — same bacteria, same neutrophils, killing
capacity restored alone:

| Condition | Dysbiotic bacterial load | Inflammation intensity | 2-year bone loss (mm) |
|---|---:|---:|---:|
| Subversion as-is (E_SUBVERT = 4.0) | 2.10 | 4.68 | 0.64 |
| Subversion halved (2.0) | 1.72 | 4.15 | 0.48 |
| **Subversion removed (0)** | **0.001** | **0.94** | **0.00** |
| C3 inhibitor alone (continuous dosing) | 1.90 | 3.51 | 0.34 |
| Antibiotic alone (no host change) | 2.10 | 4.68 | 0.63 |

This one switch matters more than any antibiotic in the file. And the same
structure also predicts the C3 inhibitor's **limit** — because the gingipain
pathway is complement-independent, even fully blocking C3 only brings C5a
down to 47%. In other words, **C3 inhibition is structurally a ~50%
intervention**, which is the mechanistic explanation for the partial response
observed in the AMY-101 phase 2a trial.

### ③ Bone is a ratchet → the value of treatment is not the final inflammation but the **integral of inflammation avoided**

Clinical attachment level (CAL) partially recovers via a long junctional
epithelium, but bone sets the floor of that recovery: `CAL ≥ BLOSS +
LJE_GAP`. Attachment cannot sit coronal to the bone. So two treatments with
identical inflammation values at 12 months will still produce different
10-year outcomes if one reached that value six months earlier.

`PDT_timing()` — identical treatment, only the start time differs (10-year
horizon, starting from 2.5 mm bone loss):

| Treatment delay | 10-year bone loss (mm) | 10-year CAL (mm) | Tooth survival (%) |
|---:|---:|---:|---:|
| 0 years | 2.70 | 4.28 | 89.2 |
| 1 year | 2.95 | 4.42 | 88.8 |
| 3 years | 3.48 | 4.71 | 88.2 |
| 6 years | 4.31 | 5.17 | 87.4 |

---

## 2. What clinical papers report lumped together, this model decomposes

> **A "2 mm reduction in probing depth" is the sum of four things with
> completely different meanings.**

```
Measured PPD = CAL − gingival margin position + probe over-penetration (inflammation-dependent)
```

`PDT_decompose()` (6.5 mm bone-loss site, 6 months after SRP):

| Component | mm | Share |
|---|---:|---:|
| True attachment gain (CAL reduction) | 0.40 | 16% |
| Gingival recession (oedema resolution + structural recession) | 1.65 | 68% |
| Reduced probe over-penetration (measurement artefact) | 0.39 | 16% |
| **Total PPD reduction** | **2.44** | 100% |

**Two-thirds of the PPD reduction is recession, and one-sixth is a
measurement artefact.** This is why PPD is a lenient metric and CAL an honest
one. Because CAL is also measured by probing, the model separately outputs
`CAL_m = CAL + over-penetration` so that it inherits the same over-penetration
error.

---

## 3. Depth-dependence is not biology but **anatomy**

The reason SRP works less well in deep pockets is not that deep pockets are
biologically different, but that **the clinician cannot reach the calculus**.
The model feeds in Waerhaug's (1978) direct measurements on extracted teeth
as-is, via a logistic function.

| Probing depth | Model calculus-removal rate | Waerhaug's measured value |
|---:|---:|---:|
| 3.0 mm | 82% | 83% |
| 4.5 mm | 67% | 61% |
| 6.5 mm | 37% | 32% |

There is one more important distinction here. Instrumentation does **two
different things**, and lumping them together is the most common error in
periodontal modelling.

- **Biofilm disruption** — nearly complete wherever the instrument reaches,
  weakly depth-dependent → *a short-term effect*
- **Calculus removal** — strongly constrained by access (table above) →
  residual calculus **reseeds** the site over the following weeks

This is why SRP always produces good short-term results and then shows
**depth-dependent relapse**.

Surgery changes **none** of the host parameters in this model. It changes
exactly one thing — access (`PPD50` 5.8 → 9.5 mm) — and pays a fixed cost in
recession (0.90 mm) and attachment trauma (0.45 mm). The **critical probing
depth emerges from that single trade-off.**

| | Model | Lindhe 1982 |
|---|---:|---:|
| SRP critical probing depth (shallower than this and attachment is **lost**) | 2.37 mm | 2.9 mm |
| Surgery critical probing depth | **4.24 mm** | **4.2 mm** |
| Depth at which surgery overtakes SRP | 7.69 mm | ~5.4 mm |

**None of these values ever entered an objective function.**

---

## 4. Sites are bimodal, clinical averages are not — and that is a prediction

A single site is bistable, so its outcome after treatment is either **heals or
does not heal**, with nothing in between. But clinical trials report
**averages over depth bands**, not individual sites, and those bands straddle
the separatrix. In other words, a clinical average that looks like a smooth
dose–response is actually **a knife-edge averaged over a distribution**. There
is no need to smooth the underlying biology, and this model does not.

`PDT_site_population()` weights active (disease-attractor) sites and quiescent
(healthy-attractor) sites by depth-specific activity prevalence to produce the
same band averages as clinical trials.

> **A testable prediction:** within the same depth band, the **distribution of
> responses should be bimodal, not unimodal around the mean** (responders /
> non-responders). This is what site-level longitudinal studies actually
> observe, and it is also why the "burst" description of periodontitis
> progression was coined.

---

## 5. Calibration

Every target below comes from the literature, and **none of them was used as
an objective function.** Parameters were set from mechanism and independent
measurements, and clinical outcomes were checked against them afterwards. To
reproduce: `Rscript pdt_mrgsolve_model.R`

| Metric | Model | Literature target |
|---|---:|---|
| SRP, 1-3 mm band PPD change (mm) | +0.04 | −0.03 (Cobb) |
| SRP, 1-3 mm band CAL change (mm) | **−0.22** | −0.34 (Cobb) |
| SRP, 4-6 mm band PPD reduction (mm) | +1.53 | +1.29 (Cobb) |
| SRP, 4-6 mm band CAL gain (mm) | +0.36 | +0.55 (Cobb) |
| SRP, ≥7 mm band PPD reduction (mm) | +1.40 | +2.16 (Cobb) |
| SRP, ≥7 mm band CAL gain (mm) | +0.29 | +1.19 (Cobb) |
| Surgery, 1-3 mm band CAL change (mm) | −0.45 | −0.60 (Lindhe) |
| Surgery, ≥7 mm band PPD reduction (mm) | **+2.92** | +2.90 (Lindhe) |
| Surgery − SRP, 4-6 mm CAL (mm) | **−0.22** | negative: SRP superior |
| Surgery − SRP, ≥7 mm CAL (mm) | **+0.31** | positive: surgery superior |
| Amox+metro added CAL gain, ≥7 mm (mm) | **+0.53** | +0.40~+0.50 (Feres) |
| SDD doxycycline added CAL gain (mm) | +0.16 | +0.30~+0.40 (Caton) |
| Local minocycline added PPD reduction (mm) | **+0.35** | +0.25~+0.30 (Williams) |
| Chlorhexidine rinse added PPD reduction (mm) | **+0.00** | ~0.00 (cannot reach subgingivally) |
| SDD, GCF collagenase reduction (%) | **63** | 60~70 (Golub) |
| Untreated alveolar bone loss (mm/year) | **0.33** | 0.15~0.40 |
| Calculus removal rate at 3.0 / 4.5 / 6.5 mm (%) | **82 / 67 / 37** | 83 / 61 / 32 (Waerhaug) |
| Maximum C5a inhibition from C3 blockade (%) | 47 | <100 (gingipain floor) |
| HbA1c reduction after periodontal therapy (%) | 0.32 | 0.43 (Cochrane 2022) |
| CRP reduction after periodontal therapy (mg/L) | **0.70** | 0.4~0.9 (D'Aiuto) |
| 6-month FMD improvement (%) | **+1.85** | +2.00 (Tonetti) |
| FMD at 24 hours after intensive treatment (%) | **−1.57** | negative (Tonetti) |
| PISA, generalised stage III (mm²) | **1239** | 1000~2000 (Nesse) |
| Smoker BOP vs non-smoker (%) | **64 vs 98** | lower in smokers (masking) |
| Plaque front–bone crest distance (mm) | **1.12** | 1.0~2.0 (Waerhaug) |
| Untreated diseased-site neutrophil killing capacity (% of healthy) | 19 | 20~50 (Maekawa) |

### What did not fit (reported honestly)

This is stated without concealment. These items say more about the model than
the ones that fit well.

1. **Underpredicts the SRP response in the ≥7 mm band** (PPD +1.40 vs +2.16,
   CAL +0.29 vs +1.19). This band mixes in very deep sites where SRP never
   crosses the separatrix, and the model's knife-edge is sharper than
   clinical reality. This is the price of representing a site mosaic
   (surfaces the instrument reaches and surfaces it does not, coexisting
   within the same pocket) as a single compartment, and it is the next point
   for improvement.
2. **Underpredicts AMY-101's anti-inflammatory effect** (28-day inflammation
   −24% vs an observed −46%). C5a inhibition itself matches well, at −45%. In
   other words, the model *lowers C5a well but underestimates C5a's share of
   the contribution to inflammation.* This implies C5a has roles beyond
   macrophage activation (endothelial activation, mast cells, tissue factor,
   and so on), and it is a testable structural gap.
3. **The SDD doxycycline adjunctive effect is somewhat low**
   (+0.16 vs +0.30~0.40), and **the reduction in smokers' response is too
   weak** (−0.22 vs −0.40~−0.50).
4. **Surgery's ≥7 mm CAL gain is low** (+0.60 vs +1.60). Because the PPD
   reduction matches almost exactly (2.92 vs 2.90), the model reproduces
   surgery's recession well but underestimates the attachment-gain component.
5. **The surgery–SRP crossover point is too deep** (7.69 vs ~5.4 mm). The
   direction and ordering are correct.

---

## 6. Model structure

### 58 ODE compartments

| Group | Count | Contents |
|---|---:|---|
| PK · delivery | 20 | Doxycycline (gut–central–GCF effect compartment), amoxicillin, metronidazole, minocycline microspheres (local depot + pocket), chlorhexidine (supragingival), AMY-101 (local), resolvin E1 (local), anti-TNF (2 compartments), denosumab (subcutaneous + central), instrumentation/surgery/regeneration kernels, acute bacteraemia |
| Microbiology | 6 | Commensal bacteria, dysbiotic bacteria, *P. gingivalis*, gingipain activity, tissue PAMP load, subgingival calculus depot |
| Innate immunity | 5 | Local C5a, crevicular neutrophils, **neutrophil killing capacity**, M1, M2 |
| Adaptive immunity | 3 | Th17, Treg, B/plasma cells |
| Mediators | 8 | IL-1β, TNF-α, IL-17A, PGE2, MMP-8, TIMP-1, RANKL, OPG |
| Structure | 9 | Osteoclasts, osteoblasts, **alveolar bone loss**, **attachment level**, gingival margin, instrumentation scarring, surgical recession, junctional epithelium barrier, GCF flow |
| Systemic | 5 | Systemic IL-6, CRP, insulin resistance, HbA1c, flow-mediated dilation (FMD) |
| Risk functions | 2 | Cumulative tooth-loss risk, cumulative MRONJ risk |

### 13 modelled treatments

Mechanical: oral hygiene instruction · SRP · open-flap surgery · regeneration
(EMD/GTR) · 3-monthly maintenance
Antimicrobial: amoxicillin+metronidazole (Van Winkelhoff) · local minocycline
microspheres · chlorhexidine rinse
Host modulation: sub-antimicrobial-dose doxycycline (Periostat) · **AMY-101
(C3 inhibitor)** · resolvin E1 · anti-TNF · NSAID
Systemic: denosumab (+MRONJ risk) · glycaemic control · smoking cessation

### Built-in analysis functions

| Function | What it does |
|---|---|
| `PDT_calibration_report()` | Table comparing 31 literature targets |
| `PDT_site_population()` / `PDT_bands()` | Site population → trial-style depth-band averages |
| `PDT_critical_depth()` | Critical-probing-depth emergence curve |
| `PDT_separatrix()` | Bisection to locate the separatrix, and whether instrumentation can cross it |
| `PDT_decompose()` | Four-component decomposition of PPD reduction |
| `PDT_timing()` | The cost of treatment delay (ratchet) |
| `PDT_subversion_counterfactual()` | Killing-capacity subversion on/off |
| `PDT_simulate_scenarios()` | 13 predefined scenarios |

---

## 7. How to run

```r
install.packages(c("mrgsolve", "shiny"))

# Model + full analysis (takes a few minutes)
Rscript pdt_mrgsolve_model.R

# Interactive dashboard
R -e 'shiny::runApp("pdt_shiny_app.R", launch.browser = TRUE)'
```

```r
# Example of individual use
source("pdt_mrgsolve_model.R")
mod <- pdt_model()

# Stabilise a 6.5 mm bone-loss site at its own attractor, then SRP + low-dose doxycycline
pt  <- param(PDT_patient(mod, bloss0 = 6.5), OHI = 0.6)
out <- mrgsim_df(pt, events = c(ev_srp(), ev_sdd()), end = 365, delta = 1)

PDT_decompose(mod, bloss0 = 6.5)     # What the PPD reduction is made of
PDT_separatrix(mod, bloss0 = 6.5)    # Can instrumentation cross the separatrix?
```

Re-rendering the map:

```bash
dot -Tsvg pdt_qsp_model.dot -o pdt_qsp_model.svg
dot -Tpng -Gdpi=150 pdt_qsp_model.dot -o pdt_qsp_model.png
```

---

## 8. Mechanistic map clusters (18)

① Host risk modifiers · ② Subgingival biofilm ecology (PSD) · ③ Keystone
virulence · ④ **Inflammophilic nutrient feedback engine** · ⑤ Complement and
C5aR1×TLR2 subversion · ⑥ Junctional epithelium barrier and innate
recognition · ⑦ Neutrophil biology · ⑧ Monocytes/macrophages/dendritic cells ·
⑨ Adaptive immunity · ⑩ Cytokine · lipid mediators · ⑪ Resolution programme ·
⑫ Connective-tissue attachment destruction · ⑬ **Alveolar bone remodelling
ratchet** · ⑭ Mechanical therapy · ⑮ Antimicrobial pharmacology · ⑯
Host-modulation pharmacology · ⑰ Clinical indices and biomarkers · ⑱ Systemic
sequelae

---

## ⚠️ Disclaimer

This model is a **qualitative/semi-quantitative QSP model for educational,
research, and hypothesis-generating purposes.** It was built from published
literature and clinical trial data but has not been independently validated
or certified, and **must not be used directly for clinical decision-making,
prescribing, or regulatory submission.** The parameters and assumptions are
illustrative approximations, and separate fitting and validation against real
patient data are required. Be sure to also read "What did not fit" in
Section 5.
