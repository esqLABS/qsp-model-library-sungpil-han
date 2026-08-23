# Hereditary haemorrhagic telangiectasia (HHT) — QSP model

**Hereditary Hemorrhagic Telangiectasia · Osler–Weber–Rendu disease**

| Deliverable | File |
|--------|------|
| Mechanistic map (203 nodes · 22 clusters) | [`hht_qsp_model.dot`](hht_qsp_model.dot) · [`svg`](hht_qsp_model.svg) · [`png`](hht_qsp_model.png) |
| mrgsolve ODE model (57 states · 18 scenarios) | [`hht_mrgsolve_model.R`](hht_mrgsolve_model.R) |
| Shiny dashboard (10 tabs) | [`hht_shiny_app.R`](../../../hereditary-hemorrhagic-telangiectasia/hht_shiny_app.R) |
| References (96 papers · each PMID individually confirmed) | [`hht_references.md`](hht_references.md) |

---

## 1. What this model does differently

HHT is usually drawn as "loss of the angiogenesis brake". That picture cannot explain
**why the same germline variant produces a 1 mm punctate lesion in the nose and a 3 cm
arteriovenous malformation in the liver**, nor **why antiangiogenic drugs reverse the
former and barely scratch the latter**. This model is built on four different claims.

### Claim 1 — ALK1/ENG/SMAD4 are shear stress **set-point regulators**, and losing them **inverts the sign** of the feedback

In normal endothelium, when wall shear stress rises the vessel remodels **inwards**, the
lumen narrows, and shear stress returns to the set point — **negative** feedback, and
what transmits that feedback is the BMP9/BMP10–ALK1–endoglin → SMAD1/5/8 axis
(PMID [27646277](https://pubmed.ncbi.nlm.nih.gov/27646277/),
[32078368](https://pubmed.ncbi.nlm.nih.gov/32078368/),
[37490341](https://pubmed.ncbi.nlm.nih.gov/37490341/),
[38727966](https://pubmed.ncbi.nlm.nih.gov/38727966/)).

In a lesion that has taken a second hit, this pathway is no longer present. Then **the
same rise in shear stress now causes outward remodelling.** The model writes this in one
line:

```
REMOD = KSH · g(WSS) · (1 − 2·pSMAD_rel)
```

| pSMAD_rel | Sign of the term | Consequence |
|-----------|-----------|------|
| 1 (normal) | negative | inward remodelling · self-limiting · no lesion |
| ~0.03 (second hit) | positive | outward remodelling · self-amplifying · lesion enlargement |

That is, **HHT is not "a disease of too much angiogenesis" but a disease in which the
control loop runs with the opposite sign.**

### Claim 2 — Shear stress is **not monotonically increasing** in lumen size. This is what separates telangiectasia from AVM

Because `Q = ΔP/(R_supply + R_lesion)` and `R_lesion ∝ S⁻⁴`,

```
WSS ∝ Q/S³ ∝ S / (1 + (S/S_supply)⁴)
```

Shear stress rises while **the lesion** dominates the resistance and **falls** once
**the supply vessel** begins to dominate it. Nasal mucosal lesions are fed by arterioles
— shear stress saturates already at a small lumen, so the positive feedback branch is
trapped below the inward-remodelling sink, the lesion has a **stable lumen**, and drugs
can reverse it. Hepatic and pulmonary shunts are fed by conduit arteries — shear stress
does not saturate, so there is a **floor** that no dose can go below.

**The same equation, two basins.** In this model telangiectasia and AVM are not two
diseases.

### Claim 3 — Anaemia is not only a consequence but an **input** to lesion growth

Hb ↓ → cardiac output ↑ to defend oxygen delivery → perfusion pressure and shear stress
↑ → outward remodelling ↑ → bleeding ↑ → Hb ↓. This loop closes from cluster 17 to
cluster 5 on the map.

**Correcting anaemia is therefore not simply supportive care but a measure that removes
a mechanical drive growing the lesions.**

### Claim 4 — The clinically decisive variable is not a score but the **flux balance**

The iron content of 1 mL of blood is `0.0347 × Hb[g/dL]` mg. Steady state requires

```
absorption − obligate losses = 0.0347 × Hb × bleeding rate
```

and inverting this gives **the maximum bleeding rate each iron replacement strategy can
sustain**. The benefit is **not continuous but threshold-shaped** — a drug that halves
the bleeding changes the game if it crosses the ceiling and changes almost nothing if it
does not.

---

## 2. Principal results

### Result 1 — the ceiling for oral iron is 17.8 mL/day, and 200 mg is no better than 65 mg

| Strategy | Iron absorbed (mg/day) | Maximum sustainable bleeding |
|------|----------------|--------------------|
| Diet only | 5.2 | **8.4 mL/day** (252 mL/month) |
| Oral 65 mg/day | 10.0 | **17.8 mL/day** (533 mL/month) |
| Oral 200 mg/day | 10.0 | **17.8 mL/day** (533 mL/month) |
| + IV 1 g / 8 weeks | 27.9 | 53.0 mL/day |
| + IV 1 g / 4 weeks | 45.7 | 88.3 mL/day |
| + IV 1 g / 2 weeks | 81.4 | 158.8 mL/day |

Because hepcidin caps daily absorption at about 10 mg, **tripling the oral dose does not
move the ceiling by even 1 mL/day.** This is a conclusion in the same direction as the
alternate-day dosing literature, and in the model it is derived from a single absorption
ceiling parameter rather than assumed.

### Result 2 — refractoriness to drugs arises from **supply-vessel geometry**, not lesion size

With a perfect antiangiogenic drug removing the whole of the drive, is the lesion
sustained by shear stress alone (ratio > 1 = cannot be regressed at any dose)?

| CO_rel | Nasal | Gastrointestinal | Liver | Lung | Brain |
|--------|------|------|-----|-----|-----|
| 1.00 | suppressed (0.50) | suppressed (0.61) | **escape (4.78)** | **escape (4.58)** | **escape (4.55)** |
| 1.45 | suppressed (0.64) | suppressed (0.77) | **escape (4.89)** | **escape (4.75)** | **escape (4.72)** |
| 2.00 | suppressed (0.78) | suppressed (0.91) | **escape (4.96)** | **escape (4.86)** | **escape (4.84)** |
| 2.60 | suppressed (0.90) | **escape (1.03)** | **escape (5.01)** | **escape (4.92)** | **escape (4.91)** |

Hepatic, pulmonary and cerebral AVMs are **already in the escape range with 5-fold
margin even at normal cardiac output**. Nasal lesions escape only once cardiac output is
2.6-fold (a cardiac index of about 7.4). This is where the reason that pulmonary AVMs
need embolisation rather than drugs comes from.

### Result 3 — the gain of the anaemia-shear stress loop is about one third

Fixing Hb and re-equilibrating the structure in the severe phenotype:

| Fixed Hb | Cardiac index | CO_rel | Mean lumen S_n | Bleeding per event (mL) | Bleeding rate (mL/day) |
|-----------|--------|--------|---------------|------------------|----------------|
| 14.6 | 2.93 | 1.027 | 2.726 | 16.53 | **129.5** |
| 11.0 | 3.55 | 1.247 | 2.797 | 17.86 | 141.8 |
| 9.0 | 4.08 | 1.431 | 2.849 | 18.88 | 151.3 |
| 7.5 | 4.62 | 1.622 | 2.898 | 19.86 | 160.3 |
| 6.0 | 5.39 | 1.892 | 2.959 | 21.14 | **172.0** |

Dropping Hb from 14.6 to 6.0 **increases the volume of blood lost itself by 33%** — with
the number of lesions untouched, purely through the mechanical feedback. The loop gain is
less than 1 so it does not run away, but **anaemia worsens itself by a third.**

### Result 4 — a substantial part of the haemodynamic benefit of bevacizumab may be **correction of anaemia** rather than anti-VEGF

Decomposing the 3-month fall in cardiac index (0.456 L/min/m²) in the reproduction of
Dupuis-Girod 2012 (PMID 22396517):

| Component | Contribution |
|------|------|
| Shunt regression | +0.241 (52.9%) |
| Correction of anaemia | +0.237 (51.9%) |

And a **counterfactual simulation**: giving IV iron 1 g/month alone with no anti-VEGF
drops the cardiac index from 4.87 to 4.56 (bevacizumab, 4.87 to 4.41). That is, **iron
alone reproduces about two thirds of the haemodynamic benefit.**

> **A testable prediction**: in HHT patients with high-output cardiac failure, correcting
> Hb alone to the same level with no anti-VEGF should reproduce a substantial part of the
> fall in cardiac index. This has never been tested, which also means that the
> bevacizumab trials did not separate the two mechanisms.

### Result 5 — intranasal bevacizumab is **not underdosed; the route is wrong**

| Route | 3-month MED (min/month) | Change | Peak submucosal concentration | Free VEGF |
|------|-------------------|------|------------------|-----------|
| Placebo | 243.6 | +0.2% | 0 nM | 62.5 pg/mL |
| Intranasal 25 mg ×3 | 238.2 | −2.1% | 2.3 nM | 62.5 |
| Intranasal 75 mg ×3 | 236.6 | −2.7% | 7.0 nM | 62.5 |
| Intranasal 7,500 mg ×3 | 230.8 | −5.1% | 698 nM | 62.5 |
| **Intranasal 75,000 mg ×3** | 228.3 | **−6.1%** | **6,984 nM** | 62.5 |
| **Systemic 5 mg/kg ×3** | **71.0** | **−70.8%** | 62.6 nM | **13.4** |

The model reproduces the flat non-response plateau between 25/50/75 mg in the published
trial (PMID 27599328) and **goes further here**: raise the dose 1,000-fold so that the
peak submucosal concentration is 100 times that of the systemic route and the effect
still stays at −6%. The reason lies not in the peak concentration but in **the time above
threshold** — put a 149 kDa IgG on a surface dominated by mucociliary clearance (t½ ~15
min) and the exposure is instantaneous, while the lesion biology integrates over weeks.
**The intranasal route cannot be rescued.**

### Result 6 — that TXA reduces only duration and not frequency is a **structural prediction**, not a fit

| Treatment | Monthly duration | Change | Episodes per month | Change |
|------|---------------|------|---------|------|
| Placebo | 243.5 | +0.1% | 15.9 | +0.2% |
| Oral TXA 3 g/day | 208.8 | **−14.1%** | 15.9 | **+0.2%** |
| Oral TXA 6 g/day | 204.1 | −16.1% | 15.9 | +0.2% |
| Topical TXA 40 mg/day | 241.7 | −0.6% | — | — |

Published values: the ATERO trial −17.3% (duration), episodes 22.1 vs 23.3 (unchanged)
(PMID [25040799](https://pubmed.ncbi.nlm.nih.gov/25040799/));
topical TXA was ineffective (PMID [27599329](https://pubmed.ncbi.nlm.nih.gov/27599329/)).

In the model TXA touches **only the clot** and not the lesion. So duration and volume
move while frequency moves by exactly 0%. **The theoretical ceiling for a perfect
antifibrinolytic is also derived, at −20.4%** — TXA cannot in principle treat HHT.

### Result 7 — the dose-response of pazopanib **inverts**

| Dose | Change in ESS (12 months) | Change in Hb | Bleeding rate | Transfusions U/3 months | ALT | Mural cell coverage |
|------|-------------------|---------|--------|--------------|-----|-------------|
| 0 mg | +0.02 | −0.35 | 153.5 | 15.5 | 24 | 0.522 |
| 50 mg | −2.40 | +1.72 | 12.3 | 1.5 | 35 | 0.589 |
| **150 mg** | **−2.56** | **+2.55** | **10.1** | **2.0** | 48 | 0.550 |
| 400 mg | −2.31 | +2.03 | 12.5 | 3.0 | 62 | 0.500 |
| 800 mg | −2.06 | +0.72 | 15.3 | 3.5 | 72 | **0.471** |

The model has pazopanib's two targets working **in opposite directions**: VEGFR2 blockade
removes the angiogenic drive and is beneficial, while PDGFRβ blockade strips mural cells
and is harmful. The result is the inversion that **150 mg is optimal and 800 mg
inferior** — a mechanistic interpretation of the practice of using 50–150 mg rather than
oncology doses in HHT, and it was not put in as a parameter but derived from the
difference between the two IC50s.

---

## 3. Calibration and validation

**What was used for calibration is 5 trials, 10 endpoints in total.**

| Trial | Endpoint | Model | Published | Verdict |
|------|------|------|------|------|
| PATH-HHT (PMID 39292928) | ESS difference vs placebo (24 weeks) | **−0.83** | −0.94 (CI −1.57~−0.31) | ✅ within CI |
| PATH-HHT | HHT-QOL difference vs placebo | **−1.52** | −1.4 | ✅ |
| PATH-HHT | Baseline ESS | 5.34 | 5.0 ± 1.5 | ✅ |
| Dupuis-Girod 2012 (PMID 22396517) | Cardiac index (baseline → 3 months) | 4.87 → **4.41** | 5.05 → 4.20 | ✅ direction · magnitude |
| Dupuis-Girod 2012 | Monthly epistaxis (3 months) | 254 → **102** | 221 → 134 | ✅ |
| Dupuis-Girod 2016 (PMID 27599328) | Intranasal 25/50/75 mg | **flat (−2 to −3%)** | flat (worse than placebo) | ✅ |
| ATERO (PMID 25040799) | Monthly duration | **−14.1%** | −17.3% | ✅ |
| ATERO | Episodes per month | **+0.2%** | unchanged | ✅ |
| Whitehead 2016 (PMID 27599329) | Topical TXA 40 mg/day | **−0.6%** | ineffective | ✅ |
| Parambil 2022 (PMID 34292451) | Transfusions (U/3 months) | 15.5 → **2.0** | 16 → 0 | ⚠️ partial |
| Parambil 2022 | Change in ESS (12 months) | **−2.56** | −4.77 | ❌ **underpredicted** |
| Parambil 2022 | Change in Hb (12 months) | **+2.55** | +4.80 | ❌ **underpredicted** |

### Four disagreements reported as they stand, without adjustment

1. **It underpredicts the Parambil pazopanib effect size by about half** (ESS −2.56 vs
   −4.77, Hb +2.55 vs +4.80). The transfusion burden, however, is largely reproduced at
   15.5 → 2.0 U/3 months. **Worth noting**: the model's −2.56 lies between the randomised
   trial (PATH-HHT's −0.94) and the uncontrolled retrospective study (Parambil's −4.77).
   Parambil is a single-centre retrospective cohort and its effect size is one no
   randomised trial has ever shown. That the model matches the randomised effect size and
   not the uncontrolled one is, if part of the latter is selection and regression to the
   mean, rather the expected direction. That does not make it not an underprediction, and
   it is reported as a disagreement as it stands.

2. **At the 6-month time point the model relapses on bevacizumab** (MED 102 → 143
   min/month). The published trial continued to improve at 6 months (134 → 43 min/month).
   The model's relapse is because the drug disappears with a t½ of 20 days after the last
   injection, and it is rather consistent with the fact that relapse after a bevacizumab
   course is common in practice and maintenance regimens had to be designed
   (PMID 25751241). But it disagrees with the 6-month point of that particular trial.

3. **The ESS floor in the mild phenotype is about 2.8.** In reality an HHT carrier with
   almost no nasal lesions can have an ESS below 1, whereas the model's duration and
   intensity domains are tied to the mean lumen and do not come down that far.

4. **The bleeding rate in the severe phenotype comes out as 153 mL/day.** Back-calculating
   from the published transfusion volume (5.3 units/month) and IV iron (1.5 g/month)
   requires 250–300 mL/day. The model reaches the same transfusion burden at a lower
   bleeding rate, and this difference narrows if some of the transfused iron is taken to
   accumulate rather than be lost — the model cannot decide which it is.

### An honest declaration about the placebo term

In the PATH-HHT reproduction the placebo group's −0.90 was put in explicitly as an
**empirical term (`PBO_MAX`), not a mechanism**. The model does not explain the placebo
response. The drug's **mechanistic** effect with the placebo term switched off is −0.85,
and the −0.83 in the table above is the difference between the two groups. That **all four
groups improved significantly in ESS** in Whitehead 2016 shows that the phenomenon this
term is reconstructing is real.

---

## 4. Model structure

### Mechanistic map — 203 nodes · 22 mechanistic clusters

1. Genotype · heterozygosity · somatic second hit · 2. The BMP9/BMP10–ALK1–endoglin
receptor complex · 3. SMAD1/5/8 transcriptional output · 4. Shear stress mechanotransduction
(the normal set point) · **5. ★ Feedback sign inversion** · 6. VEGF–VEGFR2–PI3K/AKT ·
7. Mural cell coverage · 8. Arteriovenous identity · Notch · 9. Lesion population dynamics ·
10. Nasal mucosa · epistaxis · 11. Gastrointestinal angiodysplasia · 12. Haemostasis ·
fibrinolysis · 13. Hepatic vascular malformations · high output ·
14. Pulmonary and cerebral AVMs · **15. ★ Iron flux balance** · 16. Erythropoiesis ·
anaemia · **17. ★ Cardiac output (where the loop closes)** · 18. Bevacizumab PK/PD ·
19. Pazopanib · 20. Immunomodulators · 21. Tranexamic acid · topical agents · procedures ·
22. Clinical endpoints · adverse events

### mrgsolve model — 57 ODE states

- **Pharmacokinetics 10**: bevacizumab (central · peripheral · nasal surface · nasal submucosa), pazopanib, pomalidomide, tranexamic acid
- **Signalling 18**: pSMAD1/5/8 · VEGF · AKT (5 vascular beds each), ID1, PDGF-B ×2
- **Lesions 9**: number ×2, lumen ×5 (nasal · gastrointestinal · liver · lung · brain), mural cell coverage ×2
- **Haemostasis · iron · erythrocytes 6**: plasmin, iron stores, Hb, hepcidin, EPO, reticulocytes
- **Cardiovascular 3**: cardiac output, left ventricular remodelling, right atrial pressure
- **Cumulative · scores 6**: cumulative bleeding · transfusions · IV iron, ESS, QOL, placebo term
- **Safety 5**: systolic blood pressure, proteinuria, ALT, neutrophils, VTE risk

**The five vascular beds** (nasal · gastrointestinal · liver · lung · brain) share **the
same equations** and differ only in `WSS0` (perfusion pressure) and `SSAT` (the
supply-vessel limit) — the phenotypic divergence comes out of two parameters.

### Phenotypes

| Phenotype | SEV | ANGF | GIF | HEPF | Corresponds to |
|--------|-----|------|-----|------|------|
| mild | 0.22 | 0.55 | 0.10 | 0 | — |
| moderate | 1.05 | 1.02 | 0.55 | 0.20 | The PATH-HHT enrolled group |
| severe | 1.70 | 1.25 | 45.0 | 0.30 | Parambil's transfusion-dependent group |
| hepatic | 1.00 | 1.00 | 0.25 | 1.00 | Dupuis-Girod's high-output group |
| hht1_pav | 1.00 | 1.00 | 0.30 | 0.10 (PULF 0.55) | HHT1 pulmonary AVM |
| jphht | 1.30 | 1.10 | 3.00 | 0.35 | SMAD4 overlap syndrome |

---

## 5. Defects found during development

Every equation was first written and run as an independent Python implementation
**before** being ported to mrgsolve. **Six real defects** emerged in the process, each
left as a comment at the place it was fixed.

1. **A sign error grew a non-existent AVM exponentially.** If a vascular bed has no
   second-hit lesion (`drive = 0`) the ceiling becomes 0, and the exception handling of the
   logistic cap fell to −1, so an already negative bracket was multiplied by −1. The
   pulmonary shunt index rose to 10¹³ and the cardiac index to 10¹³⁴.

2. **Shear stress was written as monotonically increasing in lumen size.** This is
   physically wrong — because `Q = ΔP/(R_supply + R_lesion)`, once the lesion becomes wider
   than the supply vessel the shear stress **falls** as `S⁻³`. Under the monotonic form no
   telangiectasia could have a stable lumen; they all ran to the anatomical ceiling, and
   there was no way to distinguish telangiectasia from AVM structurally. Fixing it to
   `S/(1+(S/S_supply)⁴)` produced the two basins.

3. **Dietary iron absorption was fixed at a constant 1.4 mg/day.** In reality it is
   demand-dependent, and fractional absorption rises from 7% to 35% as hepcidin falls.
   With it fixed, **every phenotype converged to an Hb of 11.5** whether it lost 1 mL a day
   or 150 mL.

4. **Erythropoiesis had a 15% floor.** In the form
   `min(iron_lim, 0.15 + 0.85·store_lim)`, production ran at 97% of normal even with iron
   stores at zero, and Hb never fell below 11 no matter how fast the bleeding. Replacing
   this with an explicit **iron flux balance** (recycling + absorption − obligate losses +
   mobilisation) brought the severe phenotype down to an Hb of 7.

5. **A drug concentration conversion was wrong by 1000-fold.** `amount[mg]/volume[L]` is
   mg/L, which is the same as µg/mL, and it was then multiplied by 1000 again. Exposure was
   computed 1000-fold too high for all four drugs, so **every IC50/EC50 sat in complete
   saturation and the dose-response was structurally flat.**

6. **The VEGF molar conversion was wrong by 10⁶-fold.** The mg/L factor (1e−3) was used
   for pg/mL → nM. Tissue VEGF was overestimated 10⁶-fold, so **bevacizumab was in vast
   molar deficit against the ligand and effectively did nothing.** Because this error was
   partly cancelling defect 5 (the concentration being 1000-fold higher produced partial
   binding), it only came to light the moment defect 5 was fixed — the classic case of two
   errors concealing one another.

Besides these there is **one structure that was rejected**. Writing the lesion lumen sink
as `S·(maturation + decay)` made lesions disappear entirely in patients with low drive, so
the phenotype range became a step function. It was changed so that the sink acts only on
**the dilation above the lumen it was born with** (`(maturation + decay)·(S − S_birth)`) —
because a telangiectasia is an already dilated venule and cannot remodel below that.

---

## 6. How to run it

```r
# 1) render the mechanistic map
#    dot -Tsvg hht_qsp_model.dot -o hht_qsp_model.svg
#    dot -Tpng -Gdpi=150 hht_qsp_model.dot -o hht_qsp_model.png

# 2) the mrgsolve model
source("hht_mrgsolve_model.R")
iron_ceiling()              # Result 1 — the iron absorption ceiling table
run_path_hht()              # reproduction of PATH-HHT
run_nasal_dose_ranging()    # Result 5 — why the intranasal route cannot be rescued
run_loop_gain()             # Result 3 — the anaemia-shear stress loop
run_scenarios("moderate")   # the 18 treatment scenarios

# 3) Shiny dashboard (10 tabs)
shiny::runApp("hht_shiny_app.R")
```

---

## 7. Limitations

- **The ESS is an approximate reconstruction of the domain structure of the Hoag
  instrument.** The six domains (frequency · duration · intensity · healthcare use ·
  transfusion · anaemia) and the relative magnitudes of the weights follow the original
  instrument, but the exact coefficients must be taken from the source, and this must be
  borne in mind when comparing the model's ESS with published ESS values in absolute
  terms. Comparisons of change are relatively safe.
- **The placebo response is not explained mechanistically** (see §3).
- **Haemorrhagic events from cerebral and spinal AVMs are represented only as a hazard**
  and individual events are not simulated.
- **Pregnancy and paediatric physiology are not included.**
- The model is for **educational and research purposes** and must not be used directly for
  clinical decision-making.

---

## 8. References

96 papers, every PMID individually confirmed → [`hht_references.md`](hht_references.md)

The five used directly for quantitative calibration:

- Al-Samkari H, et al. *N Engl J Med* 2024;391:1015-27 — PATH-HHT pomalidomide [PMID 39292928](https://pubmed.ncbi.nlm.nih.gov/39292928/)
- Dupuis-Girod S, et al. *JAMA* 2012;307:948-55 — intravenous bevacizumab [PMID 22396517](https://pubmed.ncbi.nlm.nih.gov/22396517/)
- Dupuis-Girod S, et al. *JAMA* 2016;316:934-42 — intranasal bevacizumab [PMID 27599328](https://pubmed.ncbi.nlm.nih.gov/27599328/)
- Gaillard S, et al. *J Thromb Haemost* 2014;12:1494-502 — ATERO tranexamic acid [PMID 25040799](https://pubmed.ncbi.nlm.nih.gov/25040799/)
- Parambil JG, et al. *Angiogenesis* 2022;25:87-97 — pazopanib [PMID 34292451](https://pubmed.ncbi.nlm.nih.gov/34292451/)

The four underpinning the structural claims:

- Baeyens N, et al. *J Cell Biol* 2016 [PMID 27646277](https://pubmed.ncbi.nlm.nih.gov/27646277/)
- Peacock HM, et al. *Arterioscler Thromb Vasc Biol* 2020 [PMID 32078368](https://pubmed.ncbi.nlm.nih.gov/32078368/)
- Banerjee K, et al. *J Clin Invest* 2023 [PMID 37490341](https://pubmed.ncbi.nlm.nih.gov/37490341/)
- Anzell AR, et al. *Angiogenesis* 2024 [PMID 38727966](https://pubmed.ncbi.nlm.nih.gov/38727966/)
