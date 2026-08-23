# Progressive Paediatric Myopia — QSP Model

[![Disease](https://img.shields.io/badge/Disease-Progressive%20Myopia-blue)]()
[![Category](https://img.shields.io/badge/Category-Ophthalmic%20%7C%20Paediatric-orange)]()
[![Compartments](https://img.shields.io/badge/ODE%20Compartments-52-green)]()
[![Scenarios](https://img.shields.io/badge/Treatment%20Scenarios-28-purple)]()
[![References](https://img.shields.io/badge/References-219-yellow)]()
[![Map](https://img.shields.io/badge/Map-169%20nodes%20%C2%B7%2019%20clusters-lightgrey)]()

## Overview

**Progressive (school) myopia** is projected to affect roughly half the world's
population by 2050, and its clinical weight sits not on the spectacle prescription
but on **axial length**. Below an axial length of 24 mm the lifetime risk of
irreversible visual impairment is 1.2%, whereas above 30 mm it is 90%
(Tideman 2016). This model is built around one structural commitment that follows
from that fact.

---

## The one structural commitment of this model

**Refraction (SER) is not a state variable.** It is an **output**, computed as the
difference between three optical elements each of which is regulated independently of
the others:

```
SER = optics( axial length , corneal power , lens power )
              element 1      element 2       element 3
              grows by       static after    falls 0.2–0.4 D
              scleral        age 3           per year through
              creep                          childhood
```

The optical expression is an exact paraxial two-thin-lens vergence trace, not a
fitted "D/mm" constant. As a result the following items become **outputs rather than
assumptions**.

- **dSER/dAL** comes out at the spectacle plane as −2.87 D/mm at an axial length of
  23 mm and −2.58 D/mm at 30 mm. At the corneal plane the same eye goes
  −2.91 → −1.73, so between axial lengths of 24 and 29 mm the sensitivity **falls by
  31%**, but the spectacle vertex-distance conversion almost cancels this and **only
  7% is left.** So the rule of thumb "2.7 D/mm" fits well only in spectacle
  refraction and does not fit for anyone thinking in contact lens or ocular
  refraction.
- **Stability in the emmetropic eye is cancellation, not standstill.** The +0.10 mm
  of annual axial growth in a normal child is a −0.28 D myopic shift, and this is
  cancelled almost exactly by an annual −0.30 D fall in lens power. Take the lens
  term out and every emmetrope becomes myopic.
- **Every myopia-control treatment acts on element 1 alone.** Element 3 goes on
  moving regardless of treatment, and therefore **biases the refractive endpoint.**

---

## Three results the simulations produced

### 1. LAMP's two endpoints are consistent with each other; the ATOM1↔ATOM2 comparison is not

Applying the model's optical expression (no fitted constants) directly to LAMP's four
arms, computing the refractive change implied by each arm's axial length change, and
reading the difference from the reported refractive change as the "implied lens
compensation":

| LAMP arm | Axial length (mm/yr) | Reported SER (D/yr) | SER implied by axial length | Implied lens compensation |
|---|---|---|---|---|
| Placebo | 0.41 | −0.81 | −1.15 | **+0.34** |
| 0.01%   | 0.36 | −0.59 | −1.01 | **+0.42** |
| 0.025%  | 0.29 | −0.46 | −0.82 | **+0.36** |
| 0.05%   | 0.20 | −0.27 | −0.56 | **+0.29** |

The implied lens compensation is **+0.353 ± 0.053 D/yr (CV 15%)** and is
**independent of the treatment arm.** Since the lens is the element the drug does not
touch, this is precisely the internal consistency check LAMP's two endpoints have to
pass — and it passes. The model's own physiological lens term is −0.270 D/yr in the
same year, inside the literature range (−0.20 ~ −0.45 D/yr).

Apply the same check **between trials** and it fails. The ATOM1 placebo arm (0.38 mm,
−1.20 D over 2 years) and the ATOM2 0.01% arm (0.41 mm, −0.49 D over 2 years) have
**effectively the same axial length while the reported refractions differ by
0.71 D.** Reconciling that would require the lens in the 0.01% arm to have lost
0.40 D more power per year than in the placebo arm, and since 0.01% takes away only
0.8 D of accommodation, no lens change of that size can arise. **Conclusion: the
ATOM1↔ATOM2 refractive comparison is not internally consistent, and it is the axial
length endpoint, which does agree between them, that must be believed.**

### 2. Atropine strikes three sites with different apparent affinities

Three constants were fitted to three arms (0.01% / 0.05% / 1%) and **0.025% was held
out**. The held-out arm is a genuine prediction:

| Concentration | Model axial length (mm/yr) | Trial | Model % reduction | Trial % reduction | Role |
|---|---|---|---|---|---|
| Placebo | 0.410 | 0.41 | 0% | 0% | Reference |
| 0.01% | 0.360 | 0.36 | 12.2% | 12.2% | Fitted |
| **0.025%** | **0.286** | **0.29** | **30.3%** | **29.3%** | **Held out → prediction** |
| 0.05% | 0.200 | 0.20 | 51.2% | 51.2% | Fitted |
| 0.1% | 0.113 | – | 72.5% | – | Prediction |
| 1% | −0.010 | −0.01 | 102.4% | 102.4% | Fitted |

The error in the held-out 0.025% arm is **1 percentage point**. The Hill slope of the
axial response comes out at **1.247** (a simple hyperbola with Hill=1 could not fit
the three arms simultaneously — defect D4), and the dose that halves axial elongation
is **0.048% w/v**, effectively the same as the 0.05% LAMP recommended.

The anterior sites (pupil and accommodation) were fitted separately to the LAMP data,
with small residuals (pupil ≤0.09 mm, accommodation ≤0.40 D). Taking the ratio of the
two curves:

| Concentration | Axial efficacy | Mydriasis fraction | Efficacy/mydriasis ratio |
|---|---|---|---|
| 0.01% | 0.122 | 0.141 | **0.86** |
| 0.05% | 0.512 | 0.418 | 1.23 |
| 0.1% | 0.725 | 0.566 | **1.28 (maximum)** |
| 1% | 1.024 | 0.899 | 1.14 |

**0.01% is the point with the worst efficacy/side-effect ratio.** That is, the case
for low-concentration atropine rests on absolute tolerability rather than on the
ratio. Because iris pigment modulates the anterior sites only, the model predicts
that a child with a light iris has **greater side effects at the same prescribed
dose** (pupil 7.26 vs 6.19 mm) and, through the fall in adherence, **a smaller axial
benefit as well** (0.236 vs 0.210 mm/yr) — a testable prediction which, if true,
supports pigment-stratified dosing.

### 3. Two actuators sharing one measurement

The biometer measures from the cornea to the RPE. The choroid lies inside that path,
changes by tens of µm within **days**, and is fully reversible. The sclera changes on
a millimetre scale over **months** and is irreversible.

```
AL_measured = AL_sclera − (choroidal thickness − baseline)/1000
```

The reversible (choroidal) component of the measured one-year effect:

| Treatment | Measured reduction (mm) | Scleral component | Choroidal component | Reversible fraction |
|---|---|---|---|---|
| Atropine 0.01% | 0.0498 | 0.0467 | 0.0031 | 6.3% |
| Atropine 0.05% | 0.2100 | 0.1940 | 0.0160 | 7.6% |
| Atropine 1% | 0.4200 | 0.3625 | 0.0575 | 13.7% |
| DIMS spectacles | 0.2361 | 0.2163 | 0.0199 | 8.4% |
| Orthokeratology | 0.2538 | 0.2314 | 0.0225 | 8.8% |
| **Red light (RLRL)** | 0.2800 | 0.2359 | 0.0442 | **15.8%** |

This decomposition also explains ATOM1's famous result — that the measured axial
length over two years on 1% atropine was **−0.02 mm (a shortening)** is the sum of
choroidal thickening (+58 µm) and a genuine arrest of growth. Since suppressing the
drive alone cannot go below the emmetropic growth rate (about 0.047 mm/yr at this
age), high-concentration atropine must act **directly on the sclera** (the KDIR term,
fitted to the 1% arm only).

---

## Two places where the simulation corrected the author

Recorded honestly.

1. It was expected that **the dioptre-per-mm conversion would fall substantially in a
   long eye**, and that this would be why the refractive endpoint most
   underestimates the effect in high-risk children. At the corneal plane it does
   indeed fall by 31%, but **the spectacle-plane conversion almost cancels it** and
   only 7% is left. The real endpoint asymmetry lies not in the optics but in **the
   convexity of risk** — at an axial length of 23 mm, 1 mm is 1.1 percentage points
   of lifetime visual impairment risk, whereas at 28 mm it is 22.8 percentage points,
   a **20-fold** difference. So "% reduction" is the wrong currency.
2. It was expected that **the rebound on stopping a choroid-led treatment (red light)
   would be visibly faster than that of a receptor-mediated treatment (atropine)**.
   It is not — the first 60 days account for 21% of the washout progression for red
   light and 26% for atropine 0.5%. The choroidal step is far too small relative to
   the persistent scleral term. Only the **magnitude** differs; the **shape** is
   indistinguishable.

---

## Where the model disagrees with the literature

Stated rather than erased.

- **The outdoor-activity arm is the weakest.** KDA was set a priori and He 2015 was
  used as the check, and it fails — the model predicts a **26% reduction** in
  three-year axial elongation for an extra 40 minutes of outdoor activity a day,
  whereas that trial measured about **11%** on refraction. This arm has not been
  adjusted to fit; it is flagged as a weak part instead.
- **Orthokeratology's efficacy is overpredicted.** Computing all four optical
  treatments from the peripheral myopic defocus derived from the optical design plus
  **a single** defocus response curve (shared with every other part of the model),
  DIMS (model 60% vs trial 62%) and MiSight (55% vs 52%) agree to within 3 percentage
  points. But HAL (68% vs 56%) and especially **orthokeratology (66% vs 43%)** are
  overpredicted. Mean |error| 10 percentage points. Given that the control-arm
  progression rates of the four trials themselves differ by a factor of 1.4
  (0.55–0.78 mm/2yr), the between-device difference the model predicts (≤14
  percentage points) is smaller than the between-trial noise — the model's position
  is that **ranking devices without a head-to-head RCT has no basis**.
- **The washout progression of the ATOM2 0.01% arm is overpredicted.** Fitting KRUP
  to the 0.5% arm alone and predicting the rest, the 0.1% arm comes out excellently
  at −0.711 D against an observed −0.68 D, but the 0.01% arm is −0.592 D against an
  observed −0.28 D. Either that cohort progressed unusually slowly, or 0.01% leaves a
  residual effect after discontinuation that the model does not carry.
- **The model's refractive progression is about 13% steeper than LAMP's** (−0.915 vs
  −0.81 D/yr). This is the same story as the lens compensation residual above (model
  0.270 vs LAMP-implied 0.353 D/yr).

---

## Mechanistic Map

[![Mechanistic Map](../../../myopia-progression/myp_qsp_model.png)](../../../myopia-progression/myp_qsp_model.svg)

*Click to go to the high-resolution SVG file. **169 nodes · 19 clusters · 246 edges.***

| Cluster | Contents |
|---|---|
| 1. Genetic and demographic susceptibility | GWAS 450+ loci, polygenic risk score, parental myopia, ethnicity, age at onset, plasticity envelope Φ(age) |
| 2. Visual environment | outdoor activity and illuminance, near work, educational intensity, urbanisation, season, sleep and circadian rhythm, chromatic environment |
| 3. Retinal image formation and the defocus signal | corneal and lens power, axial length, central and peripheral defocus, prolate eye shape, accommodative lag, pupil, higher-order aberrations, LCA, signed defocus response |
| 4. Retinal detection and the STOP/GO circuitry | cones, ipRGC/melanopsin, ON/OFF bipolar cells, amacrine cells, dopamine, glucagon/ZENK, VIP, melatonin, retinal NO, GABA, muscarinic and adenosine receptors |
| 5. RPE relay | RALDH2 → all-trans retinoic acid, RPE dopamine and muscarinic receptors, barrier and ion transport, BMP2/4 |
| 6. Choroid (the fast actuator) | thickness, blood flow, NO/eNOS, osmotic water movement, non-vascular smooth muscle, hypoxic HIF-1α, growth factors, circadian rhythm, **measurement artefact** |
| 7. Scleral matrix | fibroblast→myofibroblast, TGF-β/SMAD, MMP-2/MMP-14, TIMP-2, MMP-2:TIMP-2 balance, collagen I, proteoglycans, lysyl oxidase, integrins, hypoxia |
| 8. Scleral biomechanics (the final common pathway) | creep rate, IOP loading, scleral thickness, Laplace wall stress, **the mechanical vicious cycle**, posterior staphyloma, optic nerve deformation |
| 9. Biometry → refraction | vitreous chamber, ACD, lens thickness, corneal curvature, AL/CR ratio, **the optical identity**, dSER/dAL, **lens compensation** |
| 10. The lens (the forgotten regulator) | equatorial growth, thinning and flattening, GRIN redistribution, cessation of compensation at onset, **emmetropisation = cancellation** |
| 11. Accommodation and the ciliary muscle | ciliary muscle (low pigment → low apparent Kd), amplitude of accommodation, AC/A, zonular tension |
| 12. Atropine PK | 30 µL instillation, tear film (t½ 2 min), cornea, aqueous humour, iris (melanin binding), vitreous, **the transscleral route**, choroid, sclera, nasolacrimal duct, plasma |
| 13. Atropine PD | 2 anterior sites + 1 posterior low-affinity site, non-muscarinic contribution, **receptor upregulation**, rebound, **therapeutic index** |
| 14. Optical treatments | single vision, DIMS, HAL, MiSight, orthokeratology, progressive addition, **undercorrection (the historical error)**, **the ortho-K confounder**, red light |
| 15. Other pharmacology | 7-methylxanthine, pirenzepine, levodopa, scleral cross-linking (genipin, riboflavin), NVK002, timolol, anti-hypoxia agents |
| 16. Tolerability → adherence | photophobia, near blur, compensating spectacles, **adherence → effective dose feedback** |
| 17. Clinical endpoints | refraction, **axial length (preferred)**, incident myopia, incidence of high myopia, **CARE**, the instability of % reduction, responder analysis, statistical power |
| 18. Pathologic myopia | myopic macular degeneration, lacquer cracks, myopic CNV, retinal detachment, tractional maculopathy, glaucoma, cataract, **lifetime visual impairment**, **the convexity of risk** |
| 19. Biomarkers | optical biometry, EDI-OCT choroid, axial length percentile, peripheral refraction map, pupillometry, amplitude of accommodation, cycloplegic refraction, circadian time point, electronic eye-drop monitor |

---

## mrgsolve Model

[`myp_mrgsolve_model.R`](../../../myopia-progression/myp_mrgsolve_model.R) — **52 ODE compartments**, time unit **days**.

| Block | Compartments |
|---|---|
| Atropine PK (11) | tears · cornea · aqueous · iris · vitreous · choroid · sclera · plasma · peripheral + iris/ciliary receptor binding |
| Signalling cascade (7) | dopamine, retinoic acid, NO, TGF-β, MMP-2, TIMP-2, scleral hypoxia |
| Sclera (7) | aggrecan, collagen I, cross-linking, myofibroblasts, **creep rate**, scleral thickness, staphyloma |
| Choroid (4) | lacquer cracks, **choroidal thickness**, perfusion, hypoxia |
| Biometry (6) | **scleral axial length**, ACD, lens thickness, **lens power**, corneal curvature, peripheral refraction |
| Receptor and behaviour (4) | muscarinic receptor upregulation, adherence, amplitude of accommodation, near load |
| Devices (5) | ortho-K corneal reshaping, photobiomodulation (fast/slow), 7-MX gut and plasma |
| Mechanics and cumulative (8) | IOP, cumulative elongation, integrated myopia exposure, 5 risk hazards |

**The span of time constants**: tear film t½ 2 minutes → choroid 5 days → scleral
creep 20 days → receptor upregulation 90 days → growth trajectory 10 years.

**28 scenarios** (at the foot of the file): 6 natural history, 11 atropine
(concentration · discontinuation · tapering · adherence · iris shade · delayed
start), 7 optical and device, 4 environment and pitfalls.

**Thirteen constants were fitted**, and the file header lists which data each was
pinned to. Every other value comes either from physiology or from a structural
choice.

---

## Shiny App

[`myp_shiny_app.R`](myp_shiny_app.R) — **10 tabs**

1. **Patient profile** — baseline optics (axial length back-solved from the refraction by bisection), AL/CR ratio, risk stratum
2. **Growth trajectory** — axial length · refraction · choroid · creep, year-by-year progression table
3. **Optical decomposition** — the contributions of the three optical elements separated in dioptres, dSER/dAL at the spectacle plane against the corneal plane
4. **Choroid against sclera** — the two actuators, the reversible component fraction
5. **Atropine PK** — ocular and systemic concentrations, plasma C_max
6. **Dose-response** — efficacy, mydriasis, loss of accommodation, and their ratios
7. **Scenario comparison** — up to six arms compared simultaneously
8. **Rebound** — abrupt discontinuation against tapering, receptor upregulation
9. **Scleral biology** — MMP-2:TIMP-2, creep, scleral thickness, staphyloma
10. **Clinical endpoints** — CARE, lifetime visual impairment risk, the price of 1 mm

---

## Verification

There was no R runtime in this environment. So **all 52 equations were first
implemented and fitted in dependency-free Python RK4** and then carried across into
the mrgsolve file.

| File | Role |
|---|---|
| [`myp_reference_model.py`](../../../myopia-progression/myp_reference_model.py) | Independent re-implementation of the slow (disease) part, RK4 |
| [`myp_atropine_pk.py`](../../../myopia-progression/myp_atropine_pk.py) | The fast atropine PK subsystem (linear and autonomous, so integrated separately) |
| [`myp_calibration.py`](../../../myopia-progression/myp_calibration.py) | Fitting of the 13 constants and a verification battery of 13 sections |
| [`myp_calibration_output.txt`](../../../myopia-progression/myp_calibration_output.txt) | The full output of the script above |
| [`myp_fetch_references.py`](../../../myopia-progression/myp_fetch_references.py) | The PubMed E-utilities lookup script |

Integration check: three-year axial elongation differs by **3.9 × 10⁻¹⁴ mm** between
dt=0.100 day and dt=0.025 day (biometer repeatability is 2 × 10⁻² mm).

**Four real defects exposed** (each left as a comment where it was fixed):

| # | Defect | Symptom | Fix |
|---|---|---|---|
| D1 | Wrote the effector cascade as a **linear** coupling | TGF-β → TIMP-2 went **negative**, so MMP2/TIMP2 was 49, creep 125, and **−246 D** at three years | Every coupling rewritten as a **power law** (the sign cannot flip) |
| D2 | Used the analytical baseline as the choroidal reference for the measured axial length | t=0 failed to reproduce the specified baseline refraction (−1.78 D vs −1.50 D) | Reference the choroid **after** settling |
| D3 | RK4 intermediate stages do not guarantee positivity | A negative base with a fractional exponent → Python returns a **complex number**, and an unrelated comparison raises | All powers routed through `pw()`, which clamps the base |
| D4 | Wrote the atropine dose-response as a simple hyperbola (Hill=1) | The three arms could not be fitted simultaneously — with two constants pinned at their bounds, the 0.05% arm was still 9 percentage points short | Made the Hill exponent explicit (fitted value **1.247**) and fixed EMAXATR at 1 |

**A caution about the PK block.** Compartments 1–11 are real topical and systemic
atropine PK, and are integrated for the exposure read-outs (plasma C_max: about
8 pg/mL at 0.01% and about 765 pg/mL at 1% — consistent with the measured range of
300–900 pg/mL after instilling 1%). But **the PD reads a dose-level, cycle-steady-
state summary value.** The side-effect read-out measured in the clinic is taken about
14 hours after the night-time instillation, so it is **a time average of a
non-linear occupancy**, and `mean(Occ(C)) ≠ Occ(mean C)`. The posterior (efficacy)
site sits far below its EC50, so occupancy there is almost linear in dose and the two
forms agree. Do not rewire the PD onto the PK compartments without re-deriving the
summary values.

---

## Clinical implications the model computes

### When you start matters more than what you use

Atropine 0.05% given **always through to age 18**, varying only the age at which it
starts:

| Age at start | Duration | Untreated axial length | Treated axial length | CARE | Untreated lifetime risk | Risk avoided |
|---|---|---|---|---|---|---|
| Age 7 | 11 years | 27.14 | 25.39 | 1.75 mm | 27.3% | **20.7 %p** |
| Age 9 | 9 years | 25.41 | 24.54 | 0.88 mm | 6.8% | 3.7 %p |
| Age 12 | 6 years | 24.24 | 23.93 | 0.31 mm | 2.3% | 0.58 %p |
| Age 14 | 4 years | 23.89 | 23.74 | 0.15 mm | 1.7% | 0.21 %p |

The difference between starting at 7 and starting at 12 is **36-fold**. Because both
arms share the plasticity envelope Φ(age), the integral that is available to be
protected is itself concentrated at the early end of the age range.

### The plateau in CARE comes from the age envelope, not from the drug

Annual increments of CARE on atropine 0.05%: 0.200 → 0.150 → 0.120 → 0.095 → 0.074 →
0.057 → 0.044 → 0.034 mm. By **the fourth year it has fallen to 47% of the first
year's value**. And that happens even though the drug effect (EATR) is constant —
because both arms pass through the same Φ(age). The plateau Brennan 2021 observed is
derived here.

### NNT to prevent one case of lifetime visual impairment (age 8→18, atropine 0.05%)

| Baseline refraction | Untreated axial length | Treated axial length | Untreated risk | Treated risk | NNT |
|---|---|---|---|---|---|
| −0.50 D | 25.85 | 24.66 | 10.0% | 3.4% | 15 |
| −2.00 D | 26.65 | 25.34 | 19.2% | 6.4% | 8 |
| −4.00 D | 27.60 | 26.18 | 36.8% | 13.1% | 4 |
| −6.00 D | 28.49 | 26.99 | 57.6% | 24.6% | **3** |

The same % reduction is not the same benefit, because risk is convex in axial length.

### Tapering reduces rebound

After two years on 0.5%: abrupt discontinuation gives a washout of −0.870 D, a 90-day
taper −0.702 D, a 180-day taper −0.518 D, and a one-year taper −0.118 D. Because
receptor upregulation decays with τ=90 days, tapering has its largest effect when its
duration is of the same order as that time constant.

---

## References

[`myp_references.md`](myp_references.md) — **219 papers**.
Every entry was looked up directly in PubMed through the NCBI E-utilities (`esearch`
→ `esummary`), and **only records that were returned** were recorded. Every title,
journal, year, and PMID is a value PubMed returned; no citation was written from
memory. The 33 search terms whose lookups failed were not filled in by guesswork but
left as they stand in the final section of the file.

---

## ⚠️ Disclaimer

This is a semi-quantitative QSP model for educational and research purposes. It was
assembled on the basis of the published literature and clinical trial data but has
not been independently verified or certified, and **must not be used directly for
actual clinical decision-making, prescribing, or regulatory submission.** In
particular, do not cite the results without having read the four disagreements set
out above (the size of the outdoor-activity effect, orthokeratology efficacy, the
ATOM2 0.01% washout, and the rate of refractive progression).
