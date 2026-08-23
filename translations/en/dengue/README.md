# Dengue Fever / Severe Dengue (Dengue Shock Syndrome) — QSP Model

<sub>Dengue · Severe Dengue · DSS — one antibody response, two channels of opposite sign</sub>

<a href="../../../dengue/denv_qsp_model.svg"><img src="../../../dengue/denv_qsp_model.png" width="640" alt="Dengue QSP mechanistic map"></a>

---

## The one idea this model starts from

Dengue is written here as **a disease in which one antibody response is read through
two channels of opposite sign**.

| | What happens | What is observed |
|---|---|---|
| **Channel 1 (protective)** | neutralising IgG clears virions → viraemia ↓ → interferon ↓ | **defervescence** |
| **Channel 2 (destructive)** | **the same IgG** binds NS1 and virions → immune complexes → complement, mast cell chymase, cross-reactive T cell TNF → shedding of the endothelial glycocalyx → the reflection coefficient σ falls | **plasma leakage** |

Because the two channels share the same driving variable, "the patient deteriorates
after the fever comes down" is not two events but **one event seen twice**, and the
24–48 hour critical phase is the width of the antibody rise itself.

Inside the model, body temperature is computed only from **interferon → hypothalamic
PGE₂**, and leakage only from **NS1 · TNF · VEGF · chymase → glycocalyx**. Beyond the
antibody rise the two chains **share no equation at all.** And yet the pulse pressure
nadir leads defervescence by 21.5 hours and the haematocrit peak leads it by 20.5
hours — within a day. This is the central claim of this model, and it is a result of
computation rather than a rule.

The second axis is that **leakage stops itself only by way of shock**.

```
J = Kf · [ (Pc − Pi) − σ·(Πp − Πi) ]
```

Once σ has fallen, **Pc is the only term** that can come down fast enough to stop the
filtration, and Pc comes down because the plasma volume that was feeding it has gone.
So fluid resuscitation restores perfusion **by switching the leakage back on**, and
the fluid dose-response is not monotone but **a U with an interior optimum**.

---

## Files

| File | Contents |
|------|------|
| [`denv_qsp_model.dot`](../../../dengue/denv_qsp_model.dot) · [`.svg`](../../../dengue/denv_qsp_model.svg) · [`.png`](../../../dengue/denv_qsp_model.png) | Mechanistic map — 129 nodes · 17 clusters · 225 edges |
| [`denv_mrgsolve_model.R`](../../../dengue/denv_mrgsolve_model.R) | mrgsolve QSP model — **45 ODEs**, 18 treatment scenarios, the WHO fluid ladder |
| [`denv_reference_model.py`](../../../dengue/denv_reference_model.py) | Independent Python/scipy re-implementation (for verification) |
| [`denv_reference_output.txt`](../../../dengue/denv_reference_output.txt) | The full output of the script above — the source of every number below |
| [`denv_scenario_results.json`](../../../dengue/denv_scenario_results.json) | Scenario and sweep results (machine-readable) |
| [`denv_shiny_app.R`](../../../dengue/denv_shiny_app.R) | Shiny dashboard — 10 tabs |
| [`denv_references.md`](../../../dengue/denv_references.md) | 121 PubMed links, 14 sections |

---

## How this was checked

The mrgsolve model's 45 ODEs were re-implemented in Python/scipy **from the written-out
equations rather than by translating the code**, and solved with LSODA. That
re-implementation exposed five defects in the original draft, all of which were
corrected.

1. **The Landis–Pappenheimer colloid osmotic pressure expression was applied to the
   albumin concentration.** That expression is for *total protein*. As a result the
   baseline net filtration pressure came out at 13 mmHg rather than 1.5 mmHg, a
   healthy person leaked 5 L an hour, and every downstream number became meaningless.
2. **The anorexia term was draining 15 L over the course of a run.** The fall in
   plasma volume was coming from dehydration rather than from capillary leakage, so
   the pathophysiology was bogus in its entirety. Baseline intake and output were
   balanced and `ANOREX` was lowered from 0.42 to 0.12.
3. **There was no compartment in which effusion could accumulate.** Because lymph
   protects the systemic interstitium (140 → 1100 mL/h), no pleural effusion or
   ascites formed. The reason effusion does in fact accumulate is that the absorptive
   limit of parietal pleural lymph is ~0.65 mL/kg/h, one hundredth of the systemic
   lymph (Miserocchi 1997). **A second Starling compartment** with its own
   pressure-volume relation and its own drainage ceiling was added.
4. **Kf was set by hand and the pre- to post-capillary resistance ratio was then
   changed.** The baseline no longer matched and an uninfected person leaked +52 mL
   an hour. `calibrate()` now **derives** Kf and Kfs from the condition that "an
   uninfected host must be at a standstill" (Kf = 91.77, Kfs = 13.89 mL/h/mmHg — they
   are not free parameters).
5. **Primary infection became chronic.** The host's own antibody sat at an
   intermediate titre and enhanced itself. Putting in the avidity advantage of
   homotypic IgG for quaternary-structure epitopes (`AVIDN = 3`, de Alwis 2012) made
   the convalescent relapse disappear.

**Baseline self-test** (uninfected host, 14 days): ΔVp +1.9 mL · ΔVi −33.8 mL ·
ΔTP +0.004 g/dL · ΔHct −0.016 % · ΔMAP +0.012 mmHg. A standstill.

---

## Key results (all solver output, not assertions)

### 1. The enhancement curve is bell-shaped because it is a **product**

Two occupancies are read off the same antibody pool. Neutralisation requires
multivalent binding and so is a steep Hill function; opsonisation needs only one IgG
per virion to engage FcγR and so saturates quickly.

```
E(A) = (1 − N) · [ 1 + (Φ−1)·A/(A+K) ]
       ↑ falling Hill     ↑ rising and saturating
```

A falling term × a rising, saturating term = a bell. It was not assumed; it came out.

| Titre | 1:0 | 1:20 | 1:40 | **1:57** | 1:80 | 1:320 | 1:696 | 1:1280 | 1:2560 |
|---|---|---|---|---|---|---|---|---|---|
| E | 1.00 | 11.81 | 12.69 | **12.83 (peak)** | 12.65 | 4.85 | **1.00** | 0.231 | 0.042 |

**Salje 2018 (Nature 557:719)** measured raised risk at 1:21–1:80 and protection above
1:1280 in a cohort. The model's peak is at 1:57 and its crossover to protection at
1:696.

### 2. The age of maximum severity in infant dengue is **a division problem**

Maternal IgG decays from a cord blood titre of 1:1280 with a half-life of 43 days. The
number of half-lives needed to fall to the titre of maximum enhancement, 1:57, is
log₂(1280/57) = **4.49**.

> 4.49 × 43 days = **193 days = 6.3 months**
> Observed peak incidence of infant DHF: **6–8 months** (Kliks 1988; Chau 2009)

The inputs are just two, the cord blood titre and the IgG half-life, and **no infant
data went into the model at all.**

### 3. A 20 % rise in haematocrit means **800 mL** has already gone

| | Value |
|---|---|
| Baseline Hct (red cells 2000 mL, plasma 2800 mL) | 41.7 % |
| WHO haemoconcentration criterion (+20 %) | 50.0 % |
| Plasma volume at that Hct | 2000 mL |
| **Plasma already lost** | **800 mL = 28.6 % of the compartment** |

This is why the criterion rings late, and why it is erased in a patient who has
already been given fluid.

### 4. σ is the disease — lesion by lesion

| Variant | Plasma lost % | Total protein nadir | Pulse pressure nadir | Effusion mL | Severity |
|---|---|---|---|---|---|
| Full model | 30.6 | 5.21 | 10.8 | 1825 | 0.372 |
| **σ clamped at 0.97** | **6.6** | 6.94 | 30.7 | **77** | **0.091** |
| Kf clamped at baseline | 27.6 | 5.15 | 14.8 | 1768 | 0.262 |
| Protein sieving removed (PS=0) | 31.3 | 5.86 | 10.3 | 1420 | 0.337 |
| Serosal drainage unlimited | 28.8 | 5.29 | 13.3 | 1128 | 0.261 |
| Baroreflex removed | 29.8 | 5.51 | 9.7 | 1117 | 0.661 |

Clamp σ and the disease disappears; clamp Kf and almost nothing changes.
**Permeability is a side branch and sieving is the disease.** Switch the baroreflex
off and the effusion *falls* while the severity doubles — because vasoconstriction
lowers Pc and so blocks the leakage at the cost of perfusion.

### 5. Primary against secondary — the same virus, the same host, one titre apart

| | Primary (naive) | Secondary (1:55) | Ratio |
|---|---|---|---|
| Peak viraemia (log₁₀) | 7.56 | 7.86 | 1.04 |
| Peak NS1 (ng/mL) | 432 | 643 | 1.49 |
| Peak neutralising IgG (titre) | 182 | 620 | 3.40 |
| Peak TNF-α (pg/mL) | 27.7 | 181 | **6.55** |
| Peak chymase (ng/mL) | 7.8 | 15.5 | 2.00 |
| Glycocalyx nadir | 0.70 | 0.24 | 0.35 |
| σ nadir | 0.81 | 0.62 | 0.76 |
| Hct rise % | 13.2 | 21.7 | 1.65 |
| Plasma lost % | 20.0 | 30.6 | 1.53 |
| Total protein nadir (g/dL) | 6.77 | 5.21 | 0.77 |
| Pulse pressure nadir (mmHg) | 24.6 | **10.8** | 0.44 |
| Effusion (mL) | 480 | 1825 | 3.80 |
| Platelet nadir (×10⁹/L) | 77 | 32 | 0.42 |
| Peak AST (U/L) | 135 | 556 | 4.12 |
| Peak lactate (mmol/L) | 1.03 | 4.35 | 4.23 |
| Number of WHO warning signs | 2 | 6 | 3.00 |

The difference in peak viraemia is a mere 0.3 log — because both exhaust the same
target cell pool. **The difference comes from the immunopathology, not from the viral
load.** (On day 0.5 of illness, though, the viraemia of secondary infection is about
70 times that of primary — in agreement with Duyen 2011.)

### 6. The fluid dose-response has **an interior optimum**

| × the WHO ladder | Infused mL | Hours in shock | Effusion mL | Lactate | Severity |
|---|---|---|---|---|---|
| 0 | 0 | 66.0 | 1825 | 4.35 | 0.372 |
| 0.25 | 1671 | 47.5 | 2031 | 2.50 | 0.289 |
| 0.50 | 3342 | 34.5 | 2188 | 1.63 | 0.254 |
| **0.75** | **5014** | **24.0** | **2293** | **1.14** | **0.249** |
| 1.00 | 6685 | 15.0 | 2359 | 1.09 | 0.249 |
| 1.50 | 10028 | 4.5 | 2456 | 1.09 | 0.256 |
| 2.00 | 13370 | 4.5 | 2561 | 1.09 | 0.269 |
| 3.00 | 20055 | 4.5 | 2737 | 1.09 | 0.297 |

The curve is a U because **crystalloid raises Pc and dilutes Πp at the same time, and
both terms increase J**. That the optimum (0.75×, about 5.0 L) falls near the WHO
ladder is an independent confirmation of a protocol that was tuned on real patients.

**Resuscitation efficiency collapses three-fold**: at σ = 0.97 it takes **5.5 mL** of
crystalloid to leave 1 mL of plasma behind, but at the nadir of leakage it takes
**16.2 mL**, and on a 3× ladder **29.7 mL**.

### 7. The window for an antiviral **closes before the patient reaches hospital**

| Start of dosing (after fever onset) | Integrated NS1 exposure (% of untreated) | Hours in shock | Benefit % |
|---|---|---|---|
| 0 h | 6.2 | 0.0 | **25.0** |
| 12 h | 30.5 | 0.0 | 20.6 |
| 24 h | 83.6 | 7.0 | 6.2 |
| 36 h | 98.5 | 14.0 | 1.4 |
| 48 h | 99.9 | 15.0 | 0.2 |
| ≥ 60 h | 100.0 | 15.0 | 0.0 |

By 24 hours, 84 % of the integrated NS1 exposure has already been spent.
**Balapiravir enrolled within 72 hours of fever and celgosivir within 48 hours, and
both were negative.** The model explains those failures by **the arithmetic of the
integral** rather than by drug efficacy.

**Corticosteroids have the same cliff.** Given at the onset of fever they lower the
severity from 0.249 to 0.195, but they delay viral clearance and so raise the
integrated NS1 exposure from 15.6k to 21.1k. The value of dosing at 72 hours is 5 % —
Tam 2012 enrolled within 72 hours and found nothing.

### 8. Prophylactic platelet transfusion raises **the platelet count and not haemostasis**

Haemostasis is written as **the product** of four requirements — platelet count ×
fibrinogen × **vessel wall** × perfusion. Because it is a product, **the worst term
settles the answer.**

| | Untreated | 3-pool transfusion |
|---|---|---|
| Platelet nadir (×10⁹/L) | 32 | 44 (+12) |
| Bleeding index | 0.490 | 0.444 |
| Severity | 0.249 | 0.237 |

At the nadir the worst term is **the vessel wall**, and platelets are not the vessel
wall. **AAPT (Lye 2017, Lancet 389:1611)** reached the same conclusion — the model was
not fitted to the trial; it agrees for structural reasons.

### 9. The benefit of colloid is proportional to σ

| Degree of leakage | σ nadir | Crystalloid severity | Colloid severity | Difference |
|---|---|---|---|---|
| Moderate | 0.656 | 0.195 | 0.179 | −0.016 |
| Severe | 0.616 | 0.249 | 0.222 | −0.027 |
| Critical | 0.594 | 0.292 | 0.262 | −0.030 |

The colloid osmotic pressure term that colloid adds enters the Starling expression
**already multiplied by σ**. **Wills 2005 (NEJM 353:877)** found no overall difference
in 512 children and a colloid benefit only in the stratum with the narrowest pulse
pressure. The model reproduces that pattern without being told to.

### 10. Pre-infection titre → outcome (an **output** of the model, not an input)

Every host was given **the same mosquito inoculum** and only the antibody titre was
varied (this sweep alone is run from the moment of inoculation — because a protected
host has to be free not to become infected in the first place).

| Titre | E | Peak V (log₁₀) | Hct+% | Pulse pressure nadir | Platelet nadir | Severity |
|---|---|---|---|---|---|---|
| naive | 1.00 | 7.58 | 14.3 | 23.6 | 79 | 0.002 |
| 1:20 | 11.81 | 7.83 | 22.2 | 9.9 | 32 | 0.399 |
| 1:55 | 12.83 | 7.84 | 22.2 | 10.0 | 32 | 0.394 |
| 1:160 | 10.31 | 7.60 | 21.9 | 10.6 | 34 | 0.376 |
| 1:320 | 4.85 | 6.94 | 19.9 | 13.8 | 40 | 0.287 |
| 1:640 | 1.21 | 5.44 | 13.3 | 23.3 | 72 | 0.072 |
| 1:1280 | 0.23 | 4.33 | 0.3 | 37.1 | 250 | 0.000 |
| ≥1:2560 | ≤0.04 | 1.70 | 0.3 | 37.1 | 250 | 0.000 |

Three things follow from this single curve: **why DHF is a disease of the second
infection**, **the age distribution of infant dengue**, and the fact that **a vaccine
is not a new mechanism but merely a device for setting the titre, so that
pre-vaccination serostatus testing is arithmetically compulsory** (Sridhar 2018;
WHO 2018 position paper).

---

## The 18 scenarios

| # | Scenario | Hct+% | PP nadir | Shock h | PLT nadir | Effusion mL | Lactate | AST | Severity |
|---|---|---|---|---|---|---|---|---|---|
| S01 | Primary infection (naive) | 13.2 | 24.6 | 0.0 | 77 | 480 | 1.03 | 135 | 0.001 |
| S02 | Secondary infection, untreated | 21.7 | 10.8 | 66.0 | 32 | 1825 | 4.35 | 556 | 0.372 |
| S03 | Secondary + WHO fluids | 17.0 | 17.3 | 15.0 | 32 | 2359 | 1.09 | 184 | 0.249 |
| S04 | Under-resuscitation (50 %) | 17.7 | 16.6 | 34.5 | 32 | 2188 | 1.63 | 258 | 0.254 |
| S05 | Over-resuscitation (200 %) | 17.0 | 17.3 | 4.5 | 32 | 2561 | 1.09 | 176 | 0.269 |
| S06 | Colloid rescue | 17.0 | 17.3 | 4.5 | 32 | 2306 | 1.09 | 176 | 0.222 |
| S07 | Albumin rescue | 17.0 | 17.3 | 18.5 | 32 | 2132 | 1.09 | 187 | 0.226 |
| S08 | Prophylactic platelet transfusion | 17.0 | 17.3 | 15.0 | 44 | 2359 | 1.09 | 184 | 0.237 |
| S09 | Steroid (at presentation) | 17.0 | 17.3 | 7.5 | 35 | 2346 | 1.09 | 176 | 0.223 |
| S09b | Steroid (at fever onset) | 13.1 | 22.2 | 0.0 | 35 | 2320 | 1.04 | 193 | 0.195 |
| S10 | Antiviral at 48 h | 17.0 | 17.3 | 15.0 | 32 | 2359 | 1.09 | 184 | 0.249 |
| S11 | Antiviral at presentation | 17.0 | 17.3 | 15.0 | 33 | 2359 | 1.09 | 184 | 0.249 |
| S12 | Vaccinated seronegative | 21.7 | 10.8 | 66.0 | 32 | 1825 | 4.37 | 558 | 0.373 |
| S13 | Vaccinated seropositive | 6.3 | 31.0 | 0.0 | 130 | 566 | 1.00 | 25 | 0.001 |
| S14 | Tertiary infection (high titre) | 6.4 | 31.0 | 0.0 | 129 | 581 | 1.00 | 25 | 0.001 |
| S15 | High-dose paracetamol | 17.0 | 17.3 | 15.0 | 32 | 2361 | 1.09 | 184 | 0.249 |
| S16 | Anti-NS1 monoclonal antibody | 15.1 | 19.6 | 1.0 | 38 | 2348 | 1.05 | 170 | 0.212 |
| S17 | Early oral rehydration | 13.1 | 22.3 | 0.0 | 32 | 2533 | 1.05 | 170 | 0.251 |
| S18 | Presentation delayed by 8 hours | 19.9 | 13.4 | 27.0 | 32 | 2366 | 2.24 | 269 | 0.286 |

**S12 = S02**: a seronegative recipient of CYD-TDV is indistinguishable from a natural
secondary infection. The vaccine does not create a new mechanism; it **puts the titre
into the enhancing range**. In S13 the same vaccine pushes the titre above 1:696 and
protects.

**The time of presentation (t = 48 h) is not a moment picked by hand.** It is the
instant at which, in the course of an untreated secondary infection, the model itself
first satisfies three WHO warning signs.

---

## Model structure

**45 ODEs**, 17 functional blocks:

| Block | State variables |
|---|---|
| Virology | target cells · eclipse · productively infected cells · viraemia · NS1 |
| Antibody | pre-existing heterotypic IgG · de novo neutralising IgG · plasmablasts · antibody-secreting cells |
| Cellular immunity | activated CD8 T cells |
| Mediators | type I IFN · TNF-α · IL-10 · VEGF-A · mast cell chymase |
| Endothelium | glycocalyx · hydraulic conductivity multiplier |
| Fluid (systemic) | plasma volume · interstitial volume · plasma protein · interstitial protein |
| Fluid (serosal) | pleural and ascitic volume |
| Haematology | red cell volume · platelets · megakaryocyte capacity · antiplatelet antibody · fibrinogen · white cells |
| Liver | surviving hepatocyte fraction · AST · ALT |
| Temperature | hypothalamic PGE₂ · core temperature |
| Circulation | SVR multiplier · heart rate · lactate |
| Drugs | colloid · antiviral depot/plasma · paracetamol · methylprednisolone |
| Integrals | cumulative fluid in/out · NS1 exposure · immune complex exposure |

**Parameters that are not free.** Put in Guyton's textbook values (Pc 17.3, Pi −3.0,
Πp 27.4, Πi 8.0 mmHg, baseline lymph flow 140 mL/h) and the two filtration
coefficients are **not a matter of choice**: `calibrate()` derives Kf = 91.77 and
Kfs = 13.89 mL/h/mmHg from the condition that "an uninfected host must be standing
still".

---

## Running

```r
library(mrgsolve); library(dplyr); library(ggplot2)
mod <- mread("denv_mrgsolve_model.R")

# secondary infection, untreated
out <- mod |> param(ABH0 = 55) |> mrgsim(end = 336, delta = 0.5)
plot(out, LOG10V + NS1 + SIGMA + Hct + PP + EFFUS + PLT + TEMP ~ time/24)

# dashboard
shiny::runApp("denv_shiny_app.R")
```

```bash
# re-run the independent verification (numpy + scipy required)
python3 denv_reference_model.py
```

---

## What this model does not do

* **It does not distinguish serotypes.** DENV-1~4 differ in replicative fitness and
  in cross-reactive structure, but the model has only one scalar titre and one
  `FCROSS`. It cannot say which serotype sequence is the worst.
* **It does not handle the incubation period.** t = 0 is the onset of fever. Every
  scenario starts from the same viraemia, so antibody status is the only difference.
  Only the titre sweep is run from the moment of inoculation.
* **It does not handle age.** Children leak more than adults at the same mediator load
  (Gamble 2000). This model is built on a 70 kg adult.
* **It does not handle the spatial heterogeneity of leakage.** The pleural and
  peritoneal cavities are lumped into a single compartment with a single compliance.
* **There is no coagulation cascade.** Fibrinogen is a single state variable, and
  there is no thrombin, no protein C axis, and no fibrinolysis.

### The most exposed prediction

That **clamping σ at 0.97 drops the severity from 0.372 to 0.091** is the strongest
single claim in this model. Since there is no way of measuring the reflection
coefficient directly in a human being, it depends entirely on using syndecan-1 and
hyaluronan shedding markers as surrogates for σ. Suwarto 2017 supports that surrogacy,
but the quantitative mapping is this model's own and has not been independently
validated.

---

## ⚠️ Disclaimer

This is a qualitative and semi-quantitative QSP model for educational and research
purposes. It was assembled from the published literature and clinical trial data but
has not been independently verified or certified, and **must not be used for actual
clinical decision-making, prescribing, or regulatory submission.**
