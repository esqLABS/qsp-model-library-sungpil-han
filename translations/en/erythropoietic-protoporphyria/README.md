# Erythropoietic Protoporphyria (EPP) / X-Linked Protoporphyria (XLP) — QSP Model

<a href="../../../erythropoietic-protoporphyria/epp_qsp_model.svg"><img src="../../../erythropoietic-protoporphyria/epp_qsp_model.png" width="760" alt="EPP QSP mechanistic map"></a>

*(click the image for the zoomable SVG)*

---

## The one sentence

**Protoporphyrin IX is at once the *product* of this pathway and the *substrate of the
broken enzyme*.** Everything about this disease is therefore a **ratio** (supply ÷
processing), and the **sign** of any intervention is a property not of the drug but
of **which arm is rate-limiting**. From that single structural fact, four different
results are *computed* — rather than written into the code.

---

## Files

| File | Contents |
|------|------|
| [`epp_qsp_model.dot`](../../../erythropoietic-protoporphyria/epp_qsp_model.dot) | Source of the mechanistic map — 14 clusters, 180 nodes, 238 edges |
| [`epp_qsp_model.svg`](../../../erythropoietic-protoporphyria/epp_qsp_model.svg) · [`.png`](../../../erythropoietic-protoporphyria/epp_qsp_model.png) | The renders (`dot -Tsvg` / `dot -Tpng -Gdpi=150`) |
| [`epp_mrgsolve_model.R`](../../../erythropoietic-protoporphyria/epp_mrgsolve_model.R) | **37-ODE** mrgsolve model + 9 treatment scenarios |
| [`epp_shiny_app.R`](../../../erythropoietic-protoporphyria/epp_shiny_app.R) | **9-tab** interactive dashboard |
| [`epp_references.md`](../../../erythropoietic-protoporphyria/epp_references.md) | **78 references** (13 sections, title-based PubMed links) |
| [`epp_reference_check.py`](../../../erythropoietic-protoporphyria/epp_reference_check.py) | Dependency-free pure-Python verification implementation (a 1:1 mirror of the `$ODE` block) |
| [`epp_reference_output.txt`](../../../erythropoietic-protoporphyria/epp_reference_output.txt) | The complete output of that script — **the source of every number in this README** |

### A note on the numbers

R is not installed in this repository, so the mrgsolve model could not be run
directly. A **pure-Python implementation transcribing the `$ODE` / `$TABLE` blocks
line by line** (`epp_reference_check.py`, with no numpy or scipy) is therefore
included alongside it, and it was actually run. Every figure below is taken from the
output of that run (`epp_reference_output.txt`) and none of it was written in by
hand. To reproduce:

```bash
python3 epp_reference_check.py > epp_reference_output.txt   # about 3 minutes
```

To run the R side:

```r
source("epp_mrgsolve_model.R")   # run_all(mod) is executed automatically
shiny::runApp("epp_shiny_app.R")
```

---

## Four results the model computes rather than assumes

### 1. The ~35% ferrochelatase threshold is a **derived value**

EPP is in practice neither simply recessive nor simply dominant. The commonest
genotype is a severe FECH null allele on one side meeting, in *trans*, the
**low-expression allele IVS3-48T>C** that about 10% of Europeans carry. That is why
the same null allele is silent in the parent and manifests in the child — because
what is required is not an allele but a **ratio**.

The number "35%" appears nowhere in the model. What the model contains is only the
structure that *the normal chelation capacity is 2.9 times the normal ALAS2 flux*,
and the threshold comes out of that as **1/2.9 = 34.5%**.

| Residual FECH activity | Erythrocyte PPIX (µmol/L) | vs normal | Tolerance time (min) | Symptom-free sun (h/day) |
|---:|---:|---:|---:|---:|
| 100% | 0.13 | ×1 | 1683 | 8.00 |
| 50% | 0.53 | ×4 | 416 | 7.41 |
| 40% | 1.24 | ×9.5 | 177 | 5.36 |
| **35%** | **2.66** | **×20** | **83** | **3.23** |
| 30% | 6.12 | ×47 | 36 | 1.61 |
| 25% | 10.92 | ×84 | 20 | 0.94 |
| **15%** (typical EPP) | **21.42** | **×164** | **10.2** | **0.50** |
| 10% | 26.77 | ×205 | 8.2 | 0.40 |

The knee is near 35%, and below it the curve explodes.

### 2. The zinc protoporphyrin signal — the same overload, the opposite metal

Ferrochelatase inserts not only iron but **zinc as well**. That one line is what
separates EPP from XLP.

| | Metal-free PPIX | Zn-PP | Total EP | metal-free % |
|---|---:|---:|---:|---:|
| Normal | 0.13 | 0.85 | 0.98 | 13% |
| **EPP** (FECH 15%) | **21.42** (×164) | **1.42** (×1.7) | 22.83 | **94%** |
| **XLP** (ALAS2 GoF) | **17.49** (×134) | **8.13** (×9.6) | 25.62 | **68%** |

In both diseases free PPIX rises more than a hundredfold, but Zn-PP rises 1.7-fold in
EPP and 9.6-fold in XLP — a **5.7-fold disparity**. The reason is simple. EPP breaks
the very enzyme that would insert the zinc, so the substrate piles up and the metal
cannot get in. XLP leaves the enzyme intact and merely pours in substrate, so zinc
insertion rises along with the substrate. **There is not a single branch in the code
that distinguishes the two diseases** — all that differs between the two runs is the
two parameters `FRES` and `GOF`.

On the same plane, iron deficiency and lead poisoning raise *Zn-PP alone* and leave
metal-free PPIX normal (enzyme normal, substrate normal, only the metal wrong). Three
diseases lie in three directions on one plane — Shiny tab 9.

### 3. The two therapeutic axes **multiply rather than add**

The photodynamic dose is a **product**:
`D = irradiance × epidermal transmittance T(melanin) × [PPIX]skin × time`.
So the tolerance time is `t_tol = D_prodrome / (E × T × PPIX)`, and two drugs acting
on different factors combine **multiplicatively**.

| Regimen (mean over a 180-day season) | Erythrocyte PPIX | Melanin OD | T(405 nm) | Tolerance time | Fold | Symptom-free sun over the season |
|---|---:|---:|---:|---:|---:|---:|
| No treatment | 21.42 | 1.00 | 0.500 | 10.2 min | ×1.00 | 89 h |
| Afamelanotide 16 mg q60d | **21.42** | 1.92 | 0.270 | 19.7 min | ×1.93 | 166 h |
| Dersimelagon 300 mg | **21.42** | 1.54 | 0.350 | 15.0 min | ×1.47 | 129 h |
| Bitopertin 60 mg | 13.07 | **1.00** | **0.500** | 16.3 min | ×1.59 | 140 h |
| Beta-carotene 180 mg | 21.42 | 1.00 | 0.500 | 11.1 min | ×1.09 | 97 h |
| **Combination (afa + bito)** | 13.07 | 1.92 | 0.270 | **31.5 min** | **×3.08** | **256 h** |

```
afamelanotide alone              ×1.925
bitopertin alone                 ×1.594
additive expectation (fa + fb − 1)   ×2.518
multiplicative (Bliss) expectation (fa × fb)  ×3.067
model observed                   ×3.078      ← +22.2% over the additive expectation
```

The cells set in bold are the point. **The shielding axis does not move PPIX at all,
and the source axis does not move melanin at all.** Two axes that are entirely
orthogonal mechanistically multiply in their effect. And the statement that
"erythrocyte protoporphyrin does not change even while afamelanotide is being given"
is a **falsifiable prediction** of this model — and it was indeed what happened in the
registration trial.

### 4. Protoporphyric liver disease is not accumulation but a **saddle-node bifurcation**

PPIX is hydrophobic and **is not excreted by the kidney at all**. Its only exit is
the bile — and PPIX crystals precipitated in the bile ducts **reduce bile flow
itself**. That is positive feedback on a single excretory system — there is no reason
for there to be only one steady state.

Continuing along residual FECH activity gives **two** fold points.

```
FECH > 16.1%          monostable, only the healthy branch exists
16.1% ~ 7.4%          ★ bistable — the healthy branch and the cholestatic branch coexist
FECH < 7.4%           the healthy branch has vanished entirely, liver disease is inevitable

example at FECH 10%:  healthy branch    liver PPIX 0.88, cholestasis 0.032
                      saddle (threshold) liver PPIX 1.46, cholestasis 0.360   ← 67% above
                      disease branch    liver PPIX 5.32, cholestasis 0.714
```

That is, the typical EPP patient is biochemically stable for decades but **carries a
threshold overhead**. A transient shock such as intercurrent infection, fasting,
alcohol or haemolysis pushes the system over the saddle, and it falls onto the
cholestatic branch and stays there. This is why protoporphyric liver failure appears
not gradually but **suddenly, and usually together with a precipitant**, and also why
the lifetime risk is a few per cent rather than universal.

Rescue from the established cholestatic branch (90 days):

| Intervention | Liver PPIX | Cholestasis | ALT | Bilirubin | Erythrocyte PPIX |
|---|---:|---:|---:|---:|---:|
| No intervention | 6.72 | 0.714 | 233 | 7.75 | 28.92 |
| Colestyramine 16 g/day | 4.62 | 0.714 | 148 | 4.80 | 28.92 |
| Plasma exchange | 3.92 | 0.712 | 126 | 4.06 | 28.92 |
| **Red cell transfusion (erythropoiesis −60%)** | **0.17** | **0.000** | **25** | **0.60** | 5.96 |
| Colestyramine + transfusion | 0.15 | 0.000 | 25 | 0.60 | 5.96 |
| **Bone marrow transplantation** | **0.01** | **0.000** | **25** | **0.60** | 3.92 |

**The absence of liver transplantation from this table is not an omission.** A liver
transplant changes not one of the parameters above — because the source is the bone
marrow. That is why the disease recurs in the graft, and why **only bone marrow
transplantation, which reverses `FRES`, is curative**. The model does not learn this
separately.

---

## A fifth result: the patient is a feedback controller

The prodrome is the **sensor**, moving into the shade is the **actuator**, and a
person's reaction time is the **loop delay**. The same sunlight, the same PPIX, with
only the delay changed over a single day:

| Reaction delay | No treatment: peak dose · peak NRS · days with a reaction | Afamelanotide: peak dose · peak NRS · days with a reaction |
|---:|---|---|
| 2 min | 1.27 · **1.21** · 0.00 | 1.22 · 1.07 · 0.00 |
| 10 min | 1.72 · 1.56 · 0.00 | 1.39 · 1.13 · 0.00 |
| 20 min | 2.37 · **3.54** · 0.00 | 1.72 · 1.49 · 0.00 |
| 40 min | 3.58 · **8.28** · 0.11 | 2.35 · **3.72** · 0.00 |
| 90 min | 6.08 · 9.75 · 0.34 | 3.67 · 8.69 · 0.16 |

There are two things to read here.

**(i) The step.** Between 20 and 40 minutes the peak pain jumps from 3.5 to 8.3.
Because once the critical dose is exceeded, the mast cell and complement arms go into
self-amplification. This is why patients say either "it was a warning" or "boiling oil
for three days" and there is very little in between. The pain reaches its peak
**18 hours after** the exposure.

**(ii) The second efficacy of a shielding drug.** At a 40-minute delay, afamelanotide
turns an NRS of 8.28 into 3.72 — that is, it converts **a full phototoxic reaction
into a prodrome**. And the behaviour has not changed at all. Slow the *rate* at which
dose accumulates and **the same person's same reaction time stops being
rate-limiting**. This benefit is not captured at all by a "time spent in sunlight"
endpoint — one reason why in this disease a patient-reported outcome can move more
than the primary measure.

And one clue the model produced by itself: the time constant of the nociception state
variable is about 50 hours, so **even obeying the prodrome every single time**, daily
outdoor activity accumulates pain and reaches an equilibrium in the NRS 6 range. The
chronic low-grade burning that patients attempting a normal outdoor life complain of
is not any single reaction but **a predicted consequence of a slow state variable**.

---

## Why window glass, an overcast day and sunscreen are all useless

The action spectrum of this disease is not ultraviolet but **visible violet at
400-410 nm** (the Soret band of PPIX, ε ≈ 1.7×10⁵ M⁻¹cm⁻¹).

| Condition | Transmittance | Tolerance time | Symptom-free sun (h/day) |
|---|---:|---:|---:|
| Midsummer direct sunlight | 1.00 | 10.2 min | 0.50 |
| Behind ordinary window glass | 0.90 | 11.4 min | 0.55 |
| Overcast sky | 0.50 | 20.5 min | 0.96 |
| **SPF 50 chemical sunscreen** | **1.00** | **10.2 min** | **0.50** |
| Iron-oxide tinted sunscreen | 0.35 | 29.2 min | 1.34 |
| Long sleeves + a wide-brimmed hat | 0.15 | 68.2 min | 2.78 |
| Complete opaque cover | 0.05 | 204.7 min | 5.77 |

That SPF 50 is entered with a transmittance of 1.00 is deliberate. UV filters **do not
attenuate at all** the Soret band that drives this disease. The only things that can
block it are broad, opaque, reflective barriers — and, from inside the body,
**melanin** alone. Which is why the one drug axis that works is an MC1R agonist.

---

## Afamelanotide: a five-day drug with a sixty-day effect

| Day | Implant remaining (µg) | Plasma (µg/L) | Melanin OD | T(405) | Tolerance time |
|---:|---:|---:|---:|---:|---:|
| 0 (dosing) | 16000 | 0.00 | 1.000 | 0.500 | 10.2 min |
| 1 | 7788 | 8.62 | 1.011 | 0.496 | 10.3 min |
| 5 | 437 | 0.48 | 1.224 | 0.428 | 12.0 min |
| 10 | 12 | 0.01 | 1.583 | 0.334 | 15.3 min |
| **20** | **0** | **0.00** | **1.881** | **0.272** | **18.8 min** |
| 30 | 0 | 0.00 | 1.873 | 0.273 | 18.7 min |
| 59 (immediately before the next dose) | 0 | 0.00 | 1.555 | 0.340 | 15.0 min |
| 150 (after three repeat doses) | 0 | 0.00 | 2.270 | 0.207 | 24.7 min |

The drug in plasma has **effectively gone by day 5**, yet the protective effect
reaches its **peak around day 20, two to three weeks later**, and is still 1.5 times
baseline at day 60. The PK is more than an order of magnitude shorter than the PD.
What sets the 60-day dosing interval is not the pharmacokinetics but **the turnover of
epidermal melanin (half-life about 35 days)**. On repeat dosing the trough climbs
steadily, accumulating to a maximum of 2.27 by the third dose. Because of this
sawtooth, assessment has to be made on the seasonal mean — a single snapshot at day
180 happens to fall in a trough.

---

## Model structure

### The mechanistic map — 14 clusters

1. Genetic lesions (FECH low-expression allele · ALAS2 gain of function · CLPX)
2. The eight steps of erythroid haem biosynthesis
3. Iron homeostasis — the point where there are three signs
4. PPIX distribution — zero renal excretion, one exit
5. Photophysics — the Soret band and singlet oxygen
6. The cutaneous injury cascade — complement · mast cells · TRPA1
7. Clinical phenotype and trial endpoints
8. Hepatobiliary complications — the saddle-node
9. MC1R → eumelanin
10. Drug PK — four phases
11. Drug PD — two orthogonal axes
12. The control loop — the patient as a feedback controller
13. Differential diagnosis and the biochemical fingerprint
14. QSP outputs

### The mrgsolve model — 37 ODEs

| Compartment group | Number | State variables |
|---|---:|---|
| Haem biosynthesis | 9 | `GLY` `ALA` `PBG` `UPG` `CPG` `PPG` `PPIXE` `HEME` `FE` |
| Erythrocyte pools | 2 | `PRBC` `ZRBC` |
| Distribution · hepatobiliary | 6 | `PPL` `PSK` `PLIV` `PGUT` `CHOL` `LINJ` |
| Photobiology · injury | 7 | `OX` `MAST` `C5A` `EDEMA` `NOCI` `AVOID` `DBOUT` |
| MC1R · melanin | 2 | `TYR` `MEL` |
| Drug PK | 9 | `ADEP` `AC` `DG` `DC` `BG` `BC` `CG` `CC` `HEMC` |
| Cumulative measures | 2 | `SUNCUM` `RXNCUM` |

**Nine scenarios**: ① untreated EPP ② afamelanotide implant ③ oral dersimelagon
④ oral bitopertin ⑤ combination (Bliss test) ⑥ XLP ⑦ the sign reversal of iron
supplementation ⑧ hepatic continuation and rescue ⑨ behavioural control-loop delay
sweep.

### The Shiny app — 9 tabs

① patient · genotype ② the haem pathway ③ drug PK/PD ④ photobiology ⑤ clinical
endpoints ⑥ scenario comparison ⑦ hepatic bistability ⑧ the control loop
⑨ differential diagnosis.

---

## Calibration

| Target reproduced | Literature target | Model output |
|---|---|---|
| Normal erythrocyte protoporphyrin | < 1.5 µmol/L, mostly zinc-bound | total 0.98, Zn 87% |
| Manifest EPP | 10-40 µmol/L, metal-free > 85% | 21.4, metal-free 94% |
| FECH threshold | about 35% | 34.5% (derived) |
| XLP Zn-PP fraction | distinctly higher than in EPP | EPP 6% vs XLP 32% |
| Sun tolerance | on the order of minutes | 10.2 min |
| Afamelanotide | about double the hours of sunlight | 166 h vs 89 h (×1.87) |
| Afamelanotide and PPIX | no change | no change |
| Bitopertin PPIX reduction | about 40% | −39% (60 mg) |
| Beta-carotene | weak evidence | ×1.09 |
| Iron supplementation | beneficial in XLP, uncertain in EPP | XLP −48%, EPP +4% |
| Liver disease | 2-5% lifetime, abrupt | bistable only at FECH 7.4-16.1% |

That beta-carotene comes out weak is not a bug but a **result**. What a quencher at a
blood level of 1-3 mg/L can achieve against a singlet-oxygen lifetime of 3
microseconds is about that much, and that is the physical ceiling — and the clinical
evidence is exactly that weak. A model that made it work well would be a wrong model.

---

## Limitations

- Erythrocyte, liver and skin are each treated as a single homogeneous compartment.
  There is no PPIX distribution by erythroblast maturation stage, no difference in
  irradiation dose between skin sites, and no intralobular gradient in the liver.
- The irradiance in the photobiology is in normalised units (midsummer noon = 1.0).
  It has to be recalibrated against real 405 nm spectral irradiance by latitude,
  season and altitude before the predictions in minutes are absolute.
- The nociception and oedema cascades are phenomenological, expressed with time
  constants and thresholds, and have not been validated at the cellular level.
- The positions of the fold points of the hepatic bistability (7.4-16.1%) are
  sensitive to the biliary secretion term `VBILE` and to the crystal-cholestasis Hill
  coefficient. The qualitative conclusion (that a bistable window exists, and that it
  lies in the most severe tail) is robust; the boundary values themselves are not.
- The clinical data on bitopertin are still limited, and the IC50 used here is a value
  back-calculated so as to reproduce the reported PPIX reduction of about 40%.
- Gene therapy and antisense splicing correction are on the map but were not put into
  the ODEs.

---

## ⚠️ Disclaimer

This is a qualitative / semi-quantitative QSP model for educational and research
purposes. It was built from the public literature but has not been independently
validated or certified, and **must not be used for clinical decision-making,
prescribing, or regulatory submission.**

---

*QSP Disease Model Library · [top-level repository README](../README.md)*
