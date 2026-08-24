# Major Thermal Burn Injury — QSP Model

> **One injury, two clocks, and a controller whose gain falls off while it
> is working.**
>
> Burn resuscitation is usually taught as "give by formula and watch urine
> output." This model instead starts from **a single product.**
>
> ```
>        intravascular gain  =  infused volume × f_ret ,        f_ret = f_ret(Π_p)
> ```
>
> and adds only the fact that Π_p (plasma colloid osmotic pressure) is a
> **convex function** of protein concentration — Landis–Pappenheimer,
>
> ```
>        Π = 2.1 C + 0.16 C² + 0.009 C³        [C = total protein, g/dL]
> ```
>
> Diluting total protein **by half**, from 7.0 to 3.5 g/dL, drops colloid
> osmotic pressure from 25.6 to 9.7 mmHg — **62%** gone. So **every litre
> of protein-free crystalloid costs more than the litre before it.** f_ret
> keeps falling for as long as resuscitation continues.

---

## Three sentences this model argues

### ① The product

Intravascular volume is not **what was infused** — it is **infused volume
× retained fraction.** And the retained fraction is a decreasing function
of cumulative crystalloid, because of the convexity above. This is not a
metaphor; it is a value the simulation prints out every hour.

```
 t[h]   R[mL/h]  multiple   UO[mL/h]  VP[%]    Cp    COP
    0       900   1.00       80     100   7.00   25.6
    4      1288   1.43        6      62   4.64   14.1
    8      2359   2.62       56      90   2.31    5.8
   12      1109   2.46       77      99   2.58    6.6
   24       179   0.40      262     125   4.37   13.0
```

**Colloid osmotic pressure collapses in exactly the window where the
infusion-rate multiple (R) is rising.** The controller is pushing harder
on top of the very gain it is destroying.

### ② The controller

Urine-output-titrated resuscitation is a **closed loop.**

```
        dR/dt = g · R · (UO_target − UO_measured) / UO_target
        loop gain = (dUO/dVP) · f_ret
```

Because f_ret falls under ①, **loop gain falls.** To hold the same urine
output, the controller must keep raising the infusion rate.

> **Fluid creep is not a dosing mistake — it is the fixed point of a
> controller whose gain is falling.**

The model produces this as a calculation, not a prediction: following the
formula literally, open-loop, gives exactly 4.00 mL/kg/%TBSA, but urine
output never reaches target (24-hour urine 0.53 L); titrating to urine
output instead gives **5.97 mL/kg/%TBSA**, or **1.49 times** the
prescribed volume. Literature observations are 5.2–6.7, with an in:out
ratio of 1.2–1.6.

Colloid works not by **adding volume** but by **restoring f_ret.** That is
why, in trials, albumin's signal shows up not as improved haemodynamics
but as **a reduced in:out ratio.**

And intra-abdominal hypertension closes **a second loop**: intra-
abdominal pressure cuts renal perfusion, lowering urine output at the same
plasma volume, and the controller reads that as hypovolaemia and
**infuses more.** This loop is positive.

### ③ The two clocks

The same %TBSA drives two processes whose time constants differ by three
orders of magnitude.

| Process | Time constant | Duration |
|---|---|---|
| Capillary leak | τ ≈ 9.5 h | closes within 24–36 hours |
| Hypermetabolism | τ ≈ 104 d | persists 12–24 months |

And **the slow clock's driver is not the %TBSA at admission but today's
open wound area, A_open(t).** This single substitution explains why early
excision and grafting lower REE — **because they remove the driver
itself.** Propranolol instead blocks **the messenger** (the β-receptor).

Because the two act on different factors of the same product, they are
**not additive**:

```
                                        peak REE    day-14 lean mass
  delayed excision (day 14), no drug   :  160.5 %      −11.3 %
  delayed excision (day 14) + propranolol : 157.3 %     −0.6 %   messenger alone −3.2 pts
  early excision (day 3), no drug      :  150.6 %       −5.2 %   driver alone   −9.9 pts
  early excision (day 3) + propranolol :  150.4 %       +3.4 %   both           −10.1 pts

  sum of each alone −13.2 pts  vs  combined −10.1 pts   →  SUB-ADDITIVE
```

**Propranolol earns −3.2 points when closure is slow, and only −0.2
points when closure is fast.** A testable claim.

---

## Deliverables

| File | Contents |
|---|---|
| [`mbi_qsp_model.dot`](mbi_qsp_model.dot) · [`.svg`](mbi_qsp_model.svg) · [`.png`](mbi_qsp_model.png) | mechanistic map — 18 clusters · 197 nodes · 287 edges. Treatment nodes are coloured by **which term of the model they touch** |
| [`mbi_mrgsolve_model.R`](mbi_mrgsolve_model.R) | mrgsolve ODE model — **56 compartments · 22 scenarios**, surgery and dosing handled through an event table |
| [`mbi_shiny_app_en.R`](mbi_shiny_app_en.R) | Shiny dashboard — **12 tabs**, real-time manipulation of patient, protocol, and drugs |
| [`mbi_references_en.md`](mbi_references_en.md) | **352 items** — only records returned by an exhaustive NCBI E-utilities lookup are included |
| [`mbi_reference_python.py`](mbi_reference_python.py) | a dependency-free Python RK4 re-implementation (reference implementation for validation) |
| [`mbi_reference_output.txt`](mbi_reference_output.txt) | the full validation output of the script above |

---

## Validation: 24 real defects found

This repository's environment has no R runtime. Because **an `$ODE` block
that has never been integrated is a hypothesis, not a model**, all 56
equations were reimplemented and run in Python RK4 using only the standard
library. **24 real defects** surfaced in the process, each fixed and
marked with a `DEFECT n` comment at the site of the fix. Representative
examples:

| # | Defect | Symptom | Fix |
|---|---|---|---|
| 1 | Kf·ΔP product had no ceiling | at t=0, **7.3 L/h** filtered into burned tissue, emptying plasma in 30 minutes | filtration cannot exceed **what plasma perfusion can deliver.** This ceiling is not a numerical patch — it is the documented reason resuscitation itself increases burn oedema: the stagnant zone is underperfused, and restoring perfusion restores the delivery term |
| 3, 9 | urine-output curve was too shallow | 2.4 mL/kg/h excreted at normal volaemia → the controller reads "too much urine" from t=0 and never titrates up | steepening the curve automatically reproduces a clinically important fact: **a urine output of 0.5 mL/kg/h can still mean plasma volume is 18% short.** It is a permissive target, not a normal one |
| 5 | σ_burn = 0.15 | convection alone stripped 86 g/h of protein, exhausting a 224 g plasma pool in 3 hours | the measured albumin reflection coefficient in burned tissue is 0.3–0.5, not 0.1, and total protein including globulin leaks even less |
| 6, 14 | lymphatic return too small, with no floor | at 1400 mL/h filtration, lymph of only 101 mL/h → the quasi-steady-state Ci→(1−σ)Cp is unreachable. Conversely, while the interstitium was empty, lymph kept flowing at its baseline rate, drying the interstitium to 2% of baseline by day 5 and driving plasma COP to **196 mmHg** | systemic lymph flow is about 120 mL/h at rest and rises 10–20-fold in burned tissue. A (V/V₀) factor was also added so lymph flow stops once the interstitium is empty |
| 7, 16 | albumin was **added on top of** crystalloid | the colloid arm ended up receiving **more** total volume — the opposite of every trial | colloid is a **substitution**, not an addition. Part of the same infusion is converted to 5% albumin, applied only in the first 24 hours |
| 11, 15 | the maintenance phase replaced **measured** urine output | positive feedback (give more → more comes out → replace more) → plasma volume ran away to **183%** of baseline. Conversely, replacing only the nominal urine output produced a negative balance of 8 L/day for 4 days | once resuscitation ends, input is titrated to **the patient**, not the formula. The proportional term on plasma volume is stable and matches actual clinical practice |
| 17 | no set point on lean body mass | recovery overshot baseline, 64 → **89 kg** | anabolism is **deficit-driven.** It stops once the deficit is repaid |
| 18 | no donor-site constraint | a 45% burn closed by day 12, and an **80% burn closed on the same day 12** | the rate-limiting resource in a large burn is unburned skin, which must re-heal before it can be re-harvested. So closure time is severely nonlinear in %TBSA |
| 19 | late-phase integration step of 0.25 h | intra-abdominal pressure of **2.4 × 10⁷ mmHg** at day 12. 0.10 h and 0.04 h agreed to the decimal | 0.08 h adopted as default, with a convergence check included in the output |
| 20 | the REE decay branch inherited 6% of the rise rate | the slow clock's time constant became **26 days** instead of 104 — erasing the best-documented feature of this syndrome (hypermetabolism outlasting wound closure) | the decay branch was separated out |
| 21 | no saturation on systemic bacterial burden | burden reached 240, pinning the sepsis switch at 1 → **antibiotics could not arithmetically change the outcome** (halving 240 to 120 left the saturated logistic unchanged) | added a saturation term + made the switch **continuous** rather than a step. Vancomycin then moved mortality from 78.7% to 47.6% |
| 22 | the hazard function had **no burn-size term at all** | every scenario clustered at 16% | it is the **open wound** that carries the risk — by its area, for as long as it stays open. The revised Baux score treats one year of age the same as one %TBSA point, which a linear age term cannot reproduce, so an exponential form was used |
| 23 | the glucose block confused concentration with rate | blood glucose of **6,500 mg/dL** (18,500 in the septic group) → that value multiplied the bacterial-invasion term 83-fold, **making antibiotics look useless** | rewritten in mass units (volume of distribution, hepatic glucose production, insulin-dependent/independent uptake). Peak glucose is now 310, mean 170 mg/dL |
| 24 | bacterial growth was normalised to %TBSA | infection dynamics were identical for a 20% and a 70% burn → burn size failed to generate its own principal late cause of death | **colonisation density is independent of wound size, but the number of invasion portals is not.** Size was attached to the **invasion** term, not growth |

The complete list and the location of each fix are in the `DEFECT n`
comments in the two source files.

---

## Calibration results

Reference patient: **80 kg · age 35 · 45% TBSA · no inhalation injury**

| Metric | Literature | Model |
|---|---|---|
| UO-titrated 24-hour volume vs Parkland | 5.2–6.7 mL/kg/%TBSA | **5.97** |
| in:out ratio (infused ÷ Parkland-prescribed) | 1.2–1.6 | **1.49** |
| 24-hour total | IAH threshold about 250 mL/kg | **269 mL/kg → intra-abdominal pressure 14.2 mmHg** |
| plasma-volume nadir | 60–80% of baseline | **62%** |
| 24-hour plasma colloid osmotic pressure | 10–16 mmHg (baseline about 26) | **13.0** |
| peak weight gain | +15 to +30% | **+20%** |
| peak REE | 120–180% of predicted | **157%** |
| duration of elevated REE | months to years | **135% still at day 60** |
| lean body mass (day 14, no drug) | about −9% | **−7.1%** |
| lean body mass (day 14, propranolol) | about +9% | **+3.4%** |
| fluid reduction with high-dose ascorbate | −45% | **−36%** |
| invasive wound infection threshold | above 10⁵ CFU/g | model threshold 5 log₁₀ |
| burn sepsis mortality (treated) | 30–60% | **47.6%** (untreated 78.7%) |

### An external validator — the revised Baux score

Mortality is computed **mechanistically** (open wound area × time open,
age, inhalation injury, sepsis, compartment syndrome, ARDS, lean-mass
loss). rBaux is used only as **a scoring rubric, never as an input.**

| Patient | rBaux | rBaux prediction | Model |
|---|---|---|---|
| age 25, 20% TBSA | 45 | 0.2% | 1.3% |
| age 50, 30% TBSA | 80 | 5.9% | **7.8%** |
| age 35, 45% TBSA (reference) | 80 | 5.9% | **9.1%** |
| age 60, 40% TBSA | 100 | 29.8% | 33.8% |
| age 35, 45% TBSA + inhalation | 97 | 24.2% | 20.8% |
| age 30, 80% TBSA | 110 | 52.4% | 57.1% |
| age 45, 70% TBSA | 115 | 63.9% | 62.8% |
| age 70, 50% TBSA + inhalation | 137 | 93.5% | 99.4% |

Mean absolute error across 8 patients: **3.2 percentage points.**

Look at the two rows in bold. For **two patients with the identical rBaux
score of 80** (age 35/45% and age 50/30%), the model gives 9.1% and 7.8%.
The model was **never told** these two are equivalent — the
interchangeability of one year of age for one %TBSA point is, here, an
**output, not an assumption.**

---

## 22 scenarios

| # | Scenario | What it addresses |
|---|---|---|
| 1 | No resuscitation | historical control |
| 2 | Parkland 4 mL/kg/%, **open-loop** | what happens if the formula is followed literally |
| 3 | Parkland, titrated to UO 0.5 | reference resuscitation |
| 4 | Modified Brooke 2 mL/kg/% | a lower starting rate |
| 5 | ISBI/ABA start-low | gentler titration |
| 6 | **Fluid creep** — chasing UO 1.0, no ceiling | the controller's fixed point |
| 7 | **Opioid creep** — a 35% drop in UO at the same plasma volume | contamination of the signal itself |
| 8 | 5% albumin, starting at 8 hours | timing-dependence of colloid |
| 9 | 5% albumin, starting at 0 hours | same fraction, different timing |
| 10 | high-dose ascorbate 66 mg/kg/h × 24 h | the oxidative component of the Kf lesion |
| 11 | early excision (day 3) | removing the driver |
| 12 | delayed excision (day 14) | keeping the driver |
| 13 | standard care (day-5 excision) | reference |
| 14 | propranolol 4 mg/kg/day (early closure) | blocking the messenger, fast closure |
| 15 | **propranolol (delayed closure)** | **the interaction test** |
| 16 | oxandrolone 10 mg BID | the anabolic term |
| 17 | propranolol + oxandrolone | combining two different terms |
| 18 | insulin, target 145 mg/dL | ABA-recommended range |
| 19 | intensive insulin, target 100 mg/dL | the cost of tight control |
| 20 | invasive wound sepsis | the main cause of late death |
| 21 | above + vancomycin 1 g q12h | does the antibiotic actually move the outcome |
| 22 | modern comprehensive protocol | colloid + day-3 excision + propranolol + oxandrolone + insulin |

### Glucose control — the cost of a strict target

```
  Scenario                mean glucose  max   min   time <70   mortality
  standard care                 170     310    72        3.9 h      9.1 %
  insulin target 145             164     304    72        4.0 h      7.4 %
  intensive insulin target 100   160     304    68        8.8 h      7.4 %
```

The model was never told that "tight control fails." This result arises
because **mean glucose falls by only a few mg/dL while hypoglycaemic
exposure roughly doubles**, and because burn-related hypoglycaemia is not
an insulin-dosing problem but an **interruption of enteral feeding around
surgery.**

---

## Where the model is wrong (not hidden)

1. **Albumin started at 8 hours reduces fluid by only 5%.**
   Meta-analytic estimates run 20–40%. The model's reason is mechanistic
   and testable — under Parkland's front-loaded, first-8-hour-heavy
   infusion pattern, most of the colloid-osmotic dilution has already
   happened by the 8-hour mark. Starting the same colloid fraction at hour
   0 gives **−38%.** If timing matters this little in trials, **this
   structure is wrong.**
2. **Ascorbate reproducing Tanaka's values well is a problem, not a
   victory.** No subsequent trial has reproduced that result, so the
   model is fitted to a result that may not be real. This scenario should
   be read not as a recommendation but as "what follows if the oxidative
   hypothesis of the Kf lesion is correct."
3. **5% albumin was entered into the total-protein colloid osmotic
   pressure equation.** Albumin is more osmotically active per gram than
   globulin, so the colloid arm is **structurally conservative.**
4. **Mortality was fitted to the rBaux logistic via the open-wound hazard
   function.** So the age–%TBSA interchangeability result above is an
   internal-consistency check, not an independent validation.
5. **Propranolol's lean-mass effect is smaller than in the literature.**
   The magnitude of change vs control (about +10.5 percentage points) has
   the right direction and order of magnitude, but the control arm loses
   less than Herndon's reported −9%, so both arms sit higher than
   reported.
6. **What is missing:** coagulopathy, rhabdomyolysis kinetics,
   drug-specific nephrotoxicity. Inhalation injury enters only as a fluid
   multiplier and a risk multiplier; gas exchange is not modelled.

---

## Running

```bash
# render the mechanistic map
dot -Tsvg mbi_qsp_model.dot -o mbi_qsp_model.svg
dot -Tpng -Gdpi=150 mbi_qsp_model.dot -o mbi_qsp_model.png

# reference implementation for validation (runs without R, standard library only)
python3 mbi_reference_python.py            # full validation suite (60 days)
python3 mbi_reference_python.py --quick    # abbreviated 14-day version
```

```r
# mrgsolve model
source("mbi_mrgsolve_model.R")
out <- run_scenario("parkland_titrated")
summarise_run(out)
all <- run_all()          # all 22 scenarios

# Shiny dashboard
shiny::runApp("mbi_shiny_app_en.R")
```

---

## ⚠️ Disclaimer

This model is a **qualitative/semi-quantitative QSP model for educational
and research purposes.** It was built from published literature and
clinical trial data but has not been independently verified or certified,
and **must not be used directly for clinical decision-making,
prescribing, or regulatory submission.** Burn resuscitation in particular
requires individualised bedside judgement for every patient, and this
model's fluid-infusion outputs must never serve as grounds for a
prescription under any circumstances.
