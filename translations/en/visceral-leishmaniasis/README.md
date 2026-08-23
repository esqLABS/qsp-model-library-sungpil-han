# Visceral Leishmaniasis (Kala-azar) — Quantitative Systems Pharmacology Model

*Leishmania donovani* / *L. infantum* · 73-ODE mrgsolve model · 4 organ
compartments × 4 antileishmanial drugs · 20 treatment scenarios

---

## The one structural claim this model makes

> **The drug and the parasite never meet in plasma.**

Amastigotes live inside macrophages of spleen, liver, bone marrow and dermis.
So every kill term in this model is driven by an **intramacrophage**
concentration — while amphotericin B nephrotoxicity is driven by the **free
plasma** concentration of the same dose. Liposomal encapsulation moves those
two integrals in *opposite* directions, because the liposome is cleared by the
very cell lineage that harbours the parasite.

Everything else in the model follows from taking that seriously:

| # | Claim | Result the model computed (computed, not assumed) |
|---|---|---|
| 1 | The same mg/kg splits into two integrals | Liposomal encapsulation gives spleen exposure ×4.0 and free plasma exposure ×0.30 → therapeutic index **13.3-fold** |
| 2 | Paediatric miltefosine failure is arithmetic | Dose ∝ WT¹·⁰, clearance ∝ WT⁰·⁷⁵ ⇒ AUC ∝ WT⁰·²⁵. Exposure in a 10 kg child is **0.669 times** the adult's |
| 3 | Cure is not zero parasites but crossing a separatrix | That separatrix **disappears** between CD4 350 and 300. In HIV-VL the deciding variable comes out as ART, not the total amount of drug |
| 4 | Asymptomatic infection and kala-azar are the outcome of a race | A **2-fold** difference in priming rate decides whether disease occurs (the parasite is identical) |
| 5 | Antimony resistance is an efflux titration and cannot be beaten with dose | From 20→40 mg/kg, QTc 455→511 ms, lipase 508→985 U/L, and the resistant strain still fails |
| 6 | The synergy of combination therapy lies in time, not in concentration | In a patient with 5-fold burden, L-AmB 5 alone fails while + miltefosine for 7 days cures |
| 7 | PKDL is what is left in the organ the drug reaches least | Time above the EC50 in skin: miltefosine 1205 hours vs L-AmB 133 hours (**9-fold**) |

---

## Files

| File | Contents |
|------|------|
| [`vl_qsp_model.dot`](../../../visceral-leishmaniasis/vl_qsp_model.dot) | Mechanistic map source — **164 nodes · 18 clusters · 220 edges** |
| [`vl_qsp_model.svg`](../../../visceral-leishmaniasis/vl_qsp_model.svg) | Zoomable vector map |
| [`vl_qsp_model.png`](../../../visceral-leishmaniasis/vl_qsp_model.png) | 150 dpi raster (22087 × 4393 px) |
| [`vl_mrgsolve_model.R`](vl_mrgsolve_model.R) | 73-ODE mrgsolve model + 20 scenarios + 10 analysis functions |
| [`vl_shiny_app.R`](../../../visceral-leishmaniasis/vl_shiny_app.R) | 10-tab interactive dashboard |
| [`vl_reference_model.py`](../../../visceral-leishmaniasis/vl_reference_model.py) | Independent Python implementation (for verification, see below) |
| [`vl_reference_output.txt`](../../../visceral-leishmaniasis/vl_reference_output.txt) | Run log of that implementation — the source of every number in this README |
| [`vl_population_results.json`](../../../visceral-leishmaniasis/vl_population_results.json) | Population simulation results for 16 arms (machine-readable) |
| [`vl_references.md`](vl_references.md) | **223 references** looked up directly in PubMed, in 17 sections |

---

## Why a Python implementation ships alongside (Provenance)

The environment in which this repository was built **had no R runtime.** Rather
than commit an unverified ODE model, it seemed more honest to implement every
equation independently in Python, actually integrate it, and commit the result.
`vl_mrgsolve_model.R` and `vl_reference_model.py` are identical term for term,
and every number in this README is output from the latter
(`vl_reference_output.txt`).

Doing it that way exposed **eight real defects**. Each is marked with a comment
in both files — anyone reading a QSP model has a right to know which lines were
difficult:

1. **An extinction floor that stalls the integrator.** Putting
   `if (P < 1e-7 && dP > 0) dP = 0` *inside* the derivative makes it a
   discontinuous switch. Invisible while the burden is high, it then flickers on
   and off between solver probes the moment any compartment empties, wrecking
   LSODA's error estimate — a 28-day miltefosine regimen would not finish in
   over 100 seconds. The floor belongs *between* integration intervals, not
   inside f(t,y).
2. **CD4 without a set-point.** With a loss term and no homeostatic term, an
   immunocompetent patient drifted down to a CD4 of 4/µL over 18 months of
   follow-up and the entire host immunity axis quietly switched off. Five
   regimens that actually cure more than 90% came out as "failures".
3. **TGF-β driven by a saturating signal.** Driven by the saturating antigen
   signal AG, it holds its maximum all the way down to 1e4 parasites and locks
   macrophage activation off permanently, so no slow drug could ever hand the
   baton to the host.
4. **A hundred-fold renal amplification.** The form `KIN·Kp·C·V − KOUT·A`
   multiplies the free concentration by Kp and then by up to KIN/KOUT = 12.5 on
   top. Deoxycholate gave a peak creatinine of 32 mg/dL and a **negative serum
   potassium**.
5. **Initialising resting macrophages in the activated state.** MPHA is the
   *activated* fraction, so a healthy host is near 0. Set high, it (a) inflates
   the effect of fast-acting drugs and (b) sterilises a fresh sand-fly
   inoculation the instant it arrives, so every natural-history analysis came
   out as "failure to infect".
6. **A 1e-4 scaling of inter-organ trafficking.** An effective rate of
   2e-8 /h is not trafficking, it is zero. A dermal inoculation stayed in the
   skin for 18 months and never spread to the viscera.
7. **Symmetric inter-organ trafficking.** Written symmetrically, the
   immune-privileged skin becomes an untreatable reservoir that keeps re-seeding
   the viscera and **all 20 regimens relapsed**. Dermal macrophages are largely
   tissue-resident and do not recirculate.
8. **IL-10 driven by total burden.** A growing dermal reservoir pushed past
   K50_IL10 and immunosuppressed the host *again* — the direction is backwards.
   PKDL patients are not immunosuppressed, and the IL-10 of active kala-azar
   comes from spleen and marrow.

---

## Compartment structure (73 ODEs)

| Group | Count | State variables |
|------|-----|----------|
| Amphotericin B | 8 | Intact liposome · released free drug · peripheral · spleen/liver/marrow/skin macrophage · renal cortex |
| Miltefosine | 7 | Gut depot · central · peripheral · 4 macrophage compartments |
| Paromomycin | 9 | IM depot · central · peripheral · 4 macrophage compartments · renal cortex · cochlea |
| Antimony | 8 | Depot · plasma Sb(V) · reserve · intracellular Sb(III) in 4 compartments · deep tissue reservoir |
| Parasite | 8 | 4 organs × (replicating + quiescent) |
| Host immunity | 7 | TMEM · IFN-γ · IL-10 · TNF-α · TGF-β · CD4 · activated macrophage |
| Clinical | 10 | Spleen · liver · Hb · platelets · white cells · albumin · polyclonal IgG · body temperature · weight loss · PKDL lesions |
| Toxicity | 9 | Tubular injury · creatinine · K · Mg · hearing · QTc · lipase · ALT · gastrointestinal |
| Cumulative | 7 | 5 exposure integrals · cumulative mortality risk · cumulative mg/kg |

### Three thresholds, two sources

In the model the burden is read through three different signals. This is the
model's second structural choice:

| Signal | Half-maximal burden | Source | Meaning |
|------|------------|------|------|
| `AG` antigen availability | **0.01** units (1e4 parasites) | **Total** (skin included) | T-cell priming — why PKDL patients are leishmanin-positive |
| `AGH` IL-10 signal | **60** units (6e7 parasites) | **Visceral only** | The immunosuppressive switch that turns on only at high burden |
| `AGS` systemic symptoms | **300** units (3e8 parasites) | **Visceral only** | Fever · cachexia · hypoalbuminaemia |

`KAG ≪ K50_IL10` is not a convenience but a claim: **antigen is still abundant
across the range of burdens in which IL-10 is already switched off.** That
window, roughly 3.8 log wide from 1e4 to 6e7 parasites, is the interval in which
the host can be primed, and pushing the patient into that window and holding him
there while memory cells accumulate is the whole of what drug therapy does.

---

## Results — all from `vl_reference_output.txt`

### 1. One dose, two integrals moving in opposite directions

Giving the same total of 10 mg/kg as the liposomal formulation and as
deoxycholate:

| Measure | L-AmB | d-AmB | Ratio |
|------|------|------|----|
| Total plasma Cmax (mg/L) | 79.8 | 2.32 | 34.4 |
| **Free** plasma AUC (mg·h/L) — drives nephrotoxicity | **137.9** | **459.7** | **0.30** |
| **Spleen macrophage** AUC (mg·h/L) — drives killing | **25 688** | **6 435** | **3.99** |
| Peak creatinine (mg/dL) | 1.01 | 1.30 | 0.78 |
| Therapeutic-index surrogate (AUC_spleen / AUC_free) | **186.3** | **14.0** | **13.3-fold** |

Sweeping the liposome leakage fraction `FREL` reveals that this is an effect of
the **delivery step** and not of the formulation. And at `FREL = 1.00` (the
extreme in which the liposome leaks entirely) the therapeutic index is
**14.0** — exactly the deoxycholate value. An internal consistency test that
emerged with no free parameter:

| FREL | AUC spleen | AUC free | Therapeutic index |
|------|---------|---------|---------|
| 0.05 | 32 564 | 23.0 | **1 417** |
| 0.30 (default) | 25 688 | 137.9 | 186 |
| 0.70 | 14 687 | 321.8 | 45.6 |
| 1.00 | 6 435 | 459.7 | **14.0** ← identical to d-AmB |

### 2. Paediatric miltefosine: it is all arithmetic

Dose is prescribed in proportion to body weight (exponent 1.0) but clearance
scales with exponent 0.75. Hence **AUC ∝ WT⁰·²⁵**:

| Weight (kg) | Linear dose (mg/d) | AUC (full course) | vs 50 kg | WT⁰·²⁵ prediction |
|-----------|-----------------|-----------|-----------|------------|
| 10 | 25.0 | 18 142 | **0.669** | **0.669** |
| 20 | 50.0 | 21 576 | 0.795 | 0.795 |
| 50 | 125.0 | 27 131 | 1.000 | 1.000 |

The observed ratios match the prediction to three decimal places — not a
coincidence but an identity. Allometric dosing (a fixed dose per weight band)
reverses this deficit, but **overcorrects** by about 10% at the smallest weights
(at 7 kg the requirement is 1.48-fold and the actual is 1.64-fold). That too is
reported as it stands.

**In a single typical patient, however, this 33% difference does not change the
outcome.** All six paediatric scenarios are cured. The difference appears in the
tail of the population simulation — see the table below. This is not a
limitation of the model but a result: an exposure deficit changes the outcome
only in patients whose margin is narrow.

### 3. The cure threshold is a computable separatrix, and at CD4 it **disappears**

The residual burden a patient with memory cells TMEM at the end of treatment can
finish off unaided, without drug:

| CD4 | TMEM | Critical residual burden | = amastigotes |
|-----|------|---------------|---------------|
| 700 | 0.55 | 120.2 units | 1.2 × 10⁸ |
| 700 | 0.30 | 104.1 units | 1.0 × 10⁸ |
| 700 | **0.15** | **none** | cannot finish off any burden |
| 500 | 0.55 | 103.8 units | 1.0 × 10⁸ |
| 350 | 0.55 | 88.7 units | 8.9 × 10⁷ |
| **300 or below** | 0.55 | **none** | cannot finish off any burden |

The point is not that the threshold declines gently but that it **vanishes
between CD4 350 and 300**. Above that line the drug need only put the patient
below the line and the host finishes the job. Below it the drug has to sterilise
on its own.

And one unexpected result follows from that structure. Sweeping the HIV-VL
scenarios, **the deciding variable is not the total amount of antileishmanial
drug but ART**:

| Scenario | CD4 | nadir | Day-540 burden | Outcome | Mortality |
|---------|-----|-------|-----------|------|-------|
| Immunocompetent, L-AmB 10 mg/kg | 700 | 0.00025 | 0.00026 | cure | 4.3% |
| HIV, **on ART**, L-AmB 10 mg/kg | 90 | 0.0056 | 0.0056 | cure | 10.9% |
| HIV, **on ART**, L-AmB 30 mg/kg | 90 | 0.0017 | 0.0056 | cure | 10.9% |
| HIV, **on ART**, L-AmB 30 + miltefosine 28 days | 90 | 8.8 × 10⁻⁵ | 0.0056 | cure | 10.9% |
| HIV, **no ART**, L-AmB 30 + miltefosine 28 days | 90 | 0.034 | 1.59 × 10⁴ | **relapse** | **97.5%** |
| HIV, on ART, L-AmB 10 mg/kg | 250 | 0.0056 | 0.0056 | cure | 10.9% |

Without ART, a total of 30 mg/kg plus 28 days of miltefosine still relapses;
with ART, 10 mg/kg alone cures. That does not mean escalating the total dose is
pointless — the nadir does genuinely fall, 0.0056 → 0.0017 → 8.8 × 10⁻⁵, so in a
patient with a narrow margin that difference changes the outcome. It is that in
this model CD4 recovery is **a different axis**, and without that axis you can
push the drug axis as hard as you like and never reach the finish line. This
runs in the same direction as the WHO guideline emphasis on starting ART in
HIV-VL.

### 4. The synergy of combination therapy lies in time

Holding the total L-AmB dose at 5 mg/kg:

| Regimen | Typical patient | Patient with 5-fold burden · CD4 300 · malnourished |
|------|----------|------------------------------|
| L-AmB 5 alone | cure | **failure** (nadir 68.5) |
| L-AmB 5 + miltefosine 7 days | cure | **cure** |
| Miltefosine first, then L-AmB | cure | cure |
| L-AmB 5 + miltefosine 14 days | cure | cure (highest TMEM@EOT at 0.234) |
| L-AmB 5 + paromomycin 10 days | cure | **relapse** |

In the typical patient every arm cures — meaning that in an immunocompetent
patient the efficacy margin of amphotericin is wide enough to mask the
contribution of the partner drug, and that is itself a result worth reporting.
The reason combination therapy exists is the patient in the second column. And
the place to look at the ordering is not the nadir but the `TMEM@EOT` column:
miltefosine's 7-day tail (terminal half-life 7 days) sustains drug pressure over
the several weeks in which memory cells accumulate. Paromomycin has slow
intracellular uptake (KIN 0.020 /h) and a plasma half-life of 2–3 hours, so it
cannot generate the same tail.

### 5. Antimony: an efflux titration cannot be beaten with dose

Writing resistance as **efflux** (MRPA / transport of the Sb-trypanothione
conjugate) **rather than as target affinity** makes it structurally clear why
dose escalation fails — the efflux term is exactly proportional to the
intracellular burden the escalation creates, whereas QTc and lipase are not:

| RES_SB | Dose (mg/kg/d) | Sb(III) spleen AUC | log drop | Peak QTc | Peak lipase | Outcome |
|--------|---------------|-----------------|---------|---------|-------------|------|
| 1 (susceptible) | 20 | 827 | complete | 455 | 508 | cure |
| 3 | 20 | 276 | −1.92 | 455 | 508 | failure |
| 3 | **40** | 552 | complete | **511** | **985** | cure — at a price |
| 9 (Bihar) | 20 | 92 | −0.50 | 455 | 508 | failure |
| 9 | **40** | 184 | −1.15 | **511** | **985** | **still fails** |

### 6. The skin is a pharmacologically different organ

Macrophage exposure by organ after a single dose of L-AmB 10 mg/kg, and the time
spent above the EC50:

| Organ | MPS share | AUC (mg·h/L) | AUC/EC50 (hours) |
|------|-----------|-------------|----------------|
| Liver | 0.60 | 26 504 | 4 417 |
| Spleen | 0.10 | 25 688 | 4 279 |
| Marrow | 0.12 | 5 990 | 998 |
| **Skin** | **0.03** | **801** | **133** |

Miltefosine in that same skin compartment gives AUC/EC50 = **1 205 hours** —
9-fold. That is the quantitative content of the sentence "miltefosine is the
PKDL drug", and it is a pharmacokinetic difference, not a microbiological one.

**The PKDL risk a VL regimen leaves behind** (residual dermal burden at 180 days
after apparent cure):

| VL regimen | Dermal drug AUC | Day-180 dermal burden |
|---------|-------------|----------------|
| L-AmB 10 mg/kg single | 801 | **90.8** |
| Paromomycin 21 days | 4 607 | **93.0** |
| L-AmB 5 + miltefosine 7 days | 27 126 | 14.8 |
| L-AmB 21 mg/kg divided | 1 682 | 6.2 |
| Miltefosine 28 days | 108 476 | **0.008** |

**Treatment of PKDL itself** (starting from a patient whose viscera are cured and
only the skin remains):

| Regimen | Day-180 dermal burden | Peak lesions | Day-360 lesions |
|------|----------------|----------|-----------|
| Miltefosine 12 weeks | 0.005 | 3.91 | 0.002 |
| L-AmB 5 mg/kg weekly × 4 (total 20) | 92.0 | 6.88 | 5.90 |
| L-AmB 5 mg/kg weekly × 8 (total 40) | 1.08 | 3.78 | 0.48 |
| No treatment | 93.2 | 6.32 | 4.78 |

### 7. Asymptomatic infection and kala-azar are the outcome of a race

Starting from a small visceral focus and sweeping only the T-cell priming rate
`KTP`:

| KTP (1/h) | Peak burden | Day-540 burden | Peak spleen (cm) | Peak IL-10 | Outcome |
|-----------|----------|-----------|---------------|-----------|------|
| 1.0 × 10⁻⁵ | 1.60 × 10⁴ | 1.60 × 10⁴ | 11.4 | 1.50 | **clinical kala-azar** |
| 1.0 × 10⁻⁴ | 1.51 × 10⁴ | 1.51 × 10⁴ | 11.3 | 1.50 | **clinical kala-azar** |
| 2.0 × 10⁻⁴ | 56.6 | 51.2 | 0.51 | 0.00 | asymptomatic |
| 5.0 × 10⁻³ | 2.53 | 2.53 | 0.50 | 0.00 | asymptomatic |

**A 2-fold difference in a single host parameter separates fulminant kala-azar
from no disease at all — the parasite does not change in the slightest.** This is
the model-level explanation of why the asymptomatic:clinical infection ratio is
so sensitive to region, nutritional status and age.

### 8. The model does not support one piece of received wisdom about divided dosing

Divided L-AmB dosing with the total fixed at 10 mg/kg:

| Schedule | AUC spleen | AUC free | Peak creatinine | Outcome |
|------|---------|---------|---------------|------|
| 10 mg/kg × 1 | 25 688.2 | 137.9 | 1.01 | cure |
| 2 mg/kg × 5 | 25 688.2 | 137.9 | 1.01 | cure |
| 0.5 mg/kg × 20 | 25 688.2 | 137.9 | 0.98 | cure |
| 2 mg/kg weekly × 5 | 25 688.2 | 137.9 | 0.96 | cure |

Because amphotericin distribution is linear, **at a fixed total dose both
integrals are exactly invariant to schedule** (25 688.2 and 137.9 in every row).
What division actually changes is the **peak** alone, and the lower peak
concentration in the renal cortex brings creatinine down slightly, 1.01 →
0.96–0.98. So the honest model-based argument for single-dose administration is
not "division raises the toxicity integral" but the operational one:
**equal efficacy · equal total exposure · one visit**.

### 9. Twenty scenarios (typical patient)

| Scenario | mg/kg | log drop | Day-540 burden | Spleen (cm) | Hb | Outcome |
|---------|-------|---------|-----------|----------|----|------|
| Untreated natural history | 0 | −0.00 | 1.07 × 10⁴ | 11.0 | 3.1 | path to death |
| SSG 20 mg/kg × 30 days (East Africa) | 600 | −9.34 | 0.0002 | 0.50 | 13.7 | cure |
| SSG 20 mg/kg × 30 days (Bihar, resistant) | 600 | −0.49 | 1.07 × 10⁴ | 11.0 | 3.1 | **failure** |
| AmB deoxycholate 15 mg/kg | 15 | −7.58 | 0.0003 | 0.50 | 13.7 | cure |
| **L-AmB 10 mg/kg single** | 10 | −7.58 | 0.0003 | 0.50 | 13.7 | cure |
| L-AmB 21 mg/kg divided | 21 | −8.31 | 0.0003 | 0.50 | 13.7 | cure |
| L-AmB 3 mg/kg × 5 | 15 | −7.89 | 0.0003 | 0.50 | 13.7 | cure |
| Miltefosine 28 days (adult) | 70 | sterilising | 2.5 × 10⁻⁵ | 0.50 | 13.7 | cure |
| Miltefosine 28 days (child, linear) | 70 | sterilising | 0.0002 | 0.50 | 13.7 | cure |
| Miltefosine 28 days (child, allometric) | 105 | sterilising | 8.2 × 10⁻⁵ | 0.50 | 13.7 | cure |
| Paromomycin 21 days (India) | 315 | −7.58 | 0.0003 | 0.50 | 13.7 | cure |
| Paromomycin 21 days (East Africa) | 315 | −2.21 | 1.07 × 10⁴ | 11.0 | 3.1 | **failure** |
| L-AmB 5 + miltefosine 7 days | 22.5 | −8.09 | 0.0003 | 0.50 | 13.7 | cure |
| L-AmB 5 + paromomycin 10 days | 155 | −7.58 | 0.0003 | 0.50 | 13.7 | cure |
| Miltefosine 10 days + paromomycin 10 days | 175 | −8.70 | 0.0003 | 0.50 | 13.7 | cure |
| SSG + paromomycin 17 days | 595 | −8.75 | 0.0003 | 0.50 | 13.7 | cure |
| **HIV-VL, L-AmB 10 mg/kg, before ART** | 10 | −4.23 | 1.59 × 10⁴ | 11.4 | 2.9 | **relapse** |
| HIV-VL, L-AmB 30 + miltefosine 28 days + ART | 100 | −8.04 | 0.006 | 0.50 | 13.7 | cure |
| PKDL, miltefosine 12 weeks | 210 | −5.09 | skin 0.06 | — | — | near-complete clearance |
| PKDL, L-AmB 5 mg/kg weekly × 4 | 20 | −1.78 | skin 17.7 | — | — | **failure** |

### 10. The toxicity ledger (typical patient, attributed by regimen)

| Regimen | Peak SCr | Lowest K | Lowest Mg | Hearing (dB) | Peak QTc | Lipase | ALT |
|------|---------|-------|--------|----------|---------|---------|-----|
| L-AmB 10 mg/kg single | 1.01 | 3.03 | 1.56 | — | 400 | 30 | 25 |
| AmB deoxycholate 15 mg/kg | **1.43** | **3.08** | **1.55** | — | 400 | 30 | 25 |
| SSG 20 mg/kg × 30 days | 0.85 | 4.10 | 2.00 | — | **455** | **507** | **149** |
| Paromomycin 21 days | 0.85 | 4.10 | 2.00 | **14.7** | 400 | 30 | 25 |
| Miltefosine 28 days | 0.85 | 4.10 | 2.00 | — | 400 | 30 | 55 |
| SSG + paromomycin 17 days | 0.85 | 4.10 | 2.00 | 11.9 | 451 | 504 | 147 |
| HIV: L-AmB 30 + miltefosine | **1.29** | **2.89** | **1.47** | — | 400 | 30 | 55 |

The point is that each column belongs to a different drug. One of the real
arguments for combination therapy is not efficacy but that it **spreads the
toxicity axes**.

### 11. Population simulation — the cure rates the model predicts

60 virtual patients per arm, 12 months of follow-up. What was varied is what
actually determines the outcome: burden at presentation (splenic aspirate grade
differs between patients by orders of magnitude), drug potency, clearance,
immune competence, nutritional status, the size of the quiescent reservoir, and
**adherence to the oral drug**.

| Regimen | Cure % | Failure % | Death % | Reference (literature range) |
|------|-------|-------|-------|-----------------|
| SSG 20 mg/kg × 30 days (East Africa) | **98.3** | 1.7 | 4.5 | ~93–95% |
| SSG 20 mg/kg × 30 days (Bihar, resistant) | **0.0** | 100.0 | 69.5 | ~35–65%; the model underestimates |
| AmB deoxycholate 15 mg/kg | **100.0** | 0.0 | 4.0 | ~97% |
| L-AmB 10 mg/kg single | **100.0** | 0.0 | 3.0 | ~95–96% |
| L-AmB 21 mg/kg divided | **100.0** | 0.0 | 3.1 | ~95–98% |
| Miltefosine 28 days, adult | **91.7** | 8.3 | 10.2 | ~85–94% |
| Miltefosine 28 days, child, **linear mg/kg** | **81.7** | 18.3 | 16.1 | ~75–80% |
| Miltefosine 28 days, child, **allometric** | **96.7** | 3.3 | 6.5 | adult level once exposure is corrected |
| Paromomycin 15 mg/kg × 21 days (India) | **86.7** | 13.3 | 12.7 | ~93% |
| Paromomycin 15 mg/kg × 21 days (East Africa) | **45.0** | 55.0 | 38.8 | ~63–85%; the model underestimates |
| L-AmB 5 + miltefosine 7 days | **100.0** | 0.0 | 3.0 | ~98% |
| L-AmB 5 + paromomycin 10 days | **100.0** | 0.0 | 3.0 | ~97% |
| Miltefosine 10 days + paromomycin 10 days | **96.7** | 3.3 | 5.9 | ~97% |
| SSG + paromomycin 17 days (East Africa) | **98.3** | 0.0¹ | 4.3 | ~91% |
| HIV-VL, L-AmB 10 mg/kg, **no ART** | **0.0** | 100.0¹ | 93.2 | prognosis extremely poor |
| HIV-VL, L-AmB 30 mg/kg + miltefosine 28 days + ART | **55.0** | 45.0¹ | 7.9 | ~55–70% at 12 mo |

¹ Relapse rate. A 0.0 in the other arms means the failures were classified not
as relapse but as **primary non-response** (never reached a 3 log drop).

Two things stand out.

**First, in paediatric miltefosine, fixing a single dosing exponent gives
81.7% → 96.7%.** Fifteen points come out of arithmetic rather than biology — the
conclusion of claim 2, and it agrees with the paediatric miltefosine literature
in both direction and magnitude.

**Second, the 55% / 45% of the HIV-VL arm happens to agree with the
literature.** This arm was not a calibration target — the cure rate simply falls
out of the separatrix that disappears with CD4 (claim 3), and it overlaps the
actual observed values at 12 months (~55–70%). The 0% in the arm without ART
comes out of the same structure.

**Three discrepancies to note honestly:**

- For the most powerful regimens (L-AmB, d-AmB, L-AmB-based combinations) the
  model **overestimates** the cure rate (100% vs an actual 95–98%). That is
  because of failure mechanisms absent from the model: reinfection, loss to
  follow-up, misdiagnosis, and treatment discontinuation for reasons other than
  adherence. The parameters could have been twisted further to bring 100% down
  to 95%, but rather than absorb a missing mechanism into the EC50 the choice was
  to write it down as a limitation.
- For resistant SSG in Bihar the model **underestimates**, giving a cure rate of
  0% (the reality is 35–65%). This follows from giving efflux resistance to every
  patient uniformly through the single parameter `RES_SB = 9`, whereas a real
  population contains a mixture of susceptible and resistant strains.
- East African paromomycin is likewise **underestimated** at 45% (the reality is
  63–85%). This is the same kind of simplification, imposing the regional
  difference in susceptibility uniformly through the single factor
  `EC50_P × 1.60`.

All three discrepancies have the same character: **wherever heterogeneity within
the population has been replaced by a single parameter, the model behaves too
deterministically.**

### 12. Why relapse cannot be diagnosed serologically — a ratio of two time constants

Time for each measure to recover halfway after a curative single dose:

| Measure | Half-recovery time (days) |
|------|----------------|
| Parasite burden (log) | **7.3** |
| Platelets | 9.6 |
| Spleen size | 10.8 |
| Albumin | 11.6 |
| Haemoglobin | 14.7 |
| **Polyclonal IgG** | **93.8** |

At 180 days the visceral burden is **7.6 log** below its value at presentation
and the patient has recovered completely in clinical terms, spleen 0.50 cm ·
Hb 13.7 g/dL, yet polyclonal IgG is still **3.89 g/dL** (3.60 at presentation,
1.10 in a healthy person). The burden takes 7.3 days, the antibody 93.8 — about
a **13-fold difference**.

That ratio is not a parameter but a model output (the IgG half-life alone is a
parameter). And it is the quantitative content of why the rK39 rapid test stays
positive for years after cure and therefore **cannot be used to diagnose
relapse**. A verdict of relapse is possible only on parasitological or clinical
endpoints.

### 13. Local sensitivity — what governs the outcome

`d log₁₀P / d log p` for the day-180 burden after a single dose of L-AmB
10 mg/kg:

| Parameter | Sensitivity | |
|---------|-------|---|
| `KIMM` maximal immune kill rate | **−2.24** | host |
| `KCD4` half-maximal CD4 for T-cell help | **+1.21** | host |
| `EC50_A` amphotericin EC50 | +0.25 | drug |
| `EMAX_A` maximal amphotericin kill | −0.13 | drug |
| `KOUT_A` tissue elimination | +0.10 | drug |
| `FREL` · `CL_LIP` · `FSP_A` · `V_LIP` | ≤ 0.02 | drug |

**The day-180 burden after a curative dose has been given is determined by host
parameters, not by drug parameters.** This confirms the thesis of the whole model
in reverse — and at the same time says why matching the drug parameters
precisely is not the first priority in this model.

---

## What comes from the literature and what is a model-level hypothesis

**Directly from the literature:** plasma pharmacokinetics of each drug
(clearance · volume of distribution · absorption · half-life), the allometric
exponents 0.75/1.0, the exposure shortfall of miltefosine in children, the
frequency of antimony resistance in Bihar, clinical cure rates and adverse-event
frequencies by regimen, the fact that IL-10 is a marker of active VL, the fact
that the leishmanin response converts after treatment, and the fact that PKDL
appears after treatment.

**Model-level hypotheses (not observations):** the EC50 values in the
intracellular macrophage compartments (free intracellular concentrations cannot
be measured, so these are calibration parameters fitted to clinical cure rates),
the liposomal MPS distribution shares, the size · EC50 multiple · conversion rate
of the quiescent population, all the rate constants of the TMEM · IL-10 · MACT
modules, the dermal immune-privilege coefficient `FIMM_SK`, and the absolute
position of the separatrix. **The reported ratios, orderings and windows are far
more robust than the absolute values.**

### Known limitations

- The population simulation **overestimates** the cure rate of the two most
  powerful regimens (L-AmB, miltefosine), because failure mechanisms absent from
  the model exist: reinfection, loss to follow-up, misdiagnosis, and treatment
  discontinuation other than through non-adherence. Adherence went in as an
  explicit parameter; the rest did not.
- The first few days after the sand-fly bite — the **pre-visceralisation stage**,
  in which neutrophils and monocytes carry the parasite out of the skin — **is
  not modelled.** The natural-history analysis starts from a small visceral focus
  that has already formed.
- Parasite species (*L. donovani* vs *L. infantum*) and regional strain
  differences enter only as EC50 multiples; they are not distinguished at the
  metabolic or genetic level.
- Granuloma formation is on the map but has no explicit compartment in the
  ODEs — the efficient control by the liver and the failure of the spleen are
  expressed only through `PMAX` and a difference in growth rate.
- Paediatric scenarios were evaluated at the single point of 10 kg body weight,
  and maturation effects in neonates and infants were not included.

---

## How to run

```r
# mrgsolve model
source("vl_mrgsolve_model.R")
sim <- vl_simulate(dosing = vl_regimens(50)$S05_lamb_single10, days = 360)
vl_outcome(sim, eot_day = 1)

# reproduce every number in the README
Sys.setenv(VL_RUN_ALL = "1"); source("vl_mrgsolve_model.R")

# interactive dashboard (10 tabs)
shiny::runApp("vl_shiny_app.R")
```

```bash
# independent Python verification implementation (runs without R)
python3 vl_reference_model.py          # full analysis + population simulation
python3 vl_reference_model.py --quick  # skip the population simulation

# re-render the mechanistic map
dot -Tsvg vl_qsp_model.dot -o vl_qsp_model.svg
dot -Tpng -Gdpi=150 vl_qsp_model.dot -o vl_qsp_model.png
```

---

## Disclaimer

This model is for **research and education** and must not be used for clinical
decision-making. A substantial number of the parameters are not literature
values but values calibrated to reproduce observed cure rates, and intracellular
concentrations are unmeasurable in principle. Actual practice follows the WHO
guidelines and each country's national guidelines.

This model is for research and educational use only and must not be used for
clinical decisions.
