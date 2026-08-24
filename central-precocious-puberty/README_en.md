# Central Precocious Puberty (CPP) — QSP Model

> **Oestradiol enters the adult-height integral twice, with opposite signs.**
> A GnRH agonist removes both arms at once. Whether treatment is a benefit is
> therefore decided **not by any property of the drug but by the growth-plate
> reserve remaining at the moment it is started.**

| Deliverable | File |
|--------|------|
| 🗺️ Mechanistic map | [`cpp_qsp_model.dot`](cpp_qsp_model.dot) · [SVG](cpp_qsp_model.svg) · [PNG](cpp_qsp_model.png) — 196 nodes / 280 edges / 20 clusters |
| ⚙️ mrgsolve ODE model | [`cpp_mrgsolve_model_en.R`](cpp_mrgsolve_model_en.R) — 44 ODEs, 16 scenarios, 12 analysis functions |
| 📊 Shiny dashboard | [`cpp_shiny_app_en.R`](cpp_shiny_app_en.R) — 10 tabs |
| 📚 References | [`cpp_references_en.md`](cpp_references_en.md) — 97 items, every PMID verified |
| 🔬 Independent cross-check code | [`cpp_reference_check.py`](cpp_reference_check.py) → [`cpp_reference_output.txt`](cpp_reference_output.txt) |

---

## 1. What this model claims

Write a precocious puberty model as "puberty is too early → suppress it" and it
cannot answer the question that matters most in the clinic (**whom should we
treat**). This model instead puts adult height as a single integral and arranges for
oestradiol to act inside that integral through **two pathways of opposite sign** at
the same time.

```
adult height = HT(t₀) + ∫ GV(E2, IGF-1, GPRES) dt        (over the interval in which the growth plate is open)

  (+)  E2 → GH pulse amplitude ×2.5–2.8 → IGF-1 → chondrocyte proliferation
       ⟹ E2 raises the integrand (the growth velocity)

  (−)  E2 → growth-plate ERα → irreversible consumption of proliferative reserve → end of skeletal maturation
       ⟹ E2 pulls the upper limit of the integral forward
```

A GnRH agonist removes both of these arms **at once**. The sign of
`d(adult height)/d(treatment)` therefore depends not on the drug but on "is the
(−) arm still dominant right now".

**Nowhere in the model is a rule such as** *"treat before age 8"* **coded.** That
rule is an output (§3).

### Three structural claims

**① Bone age is "the rate of maturation relative to a normal child of the same age".**
Because the Greulich–Pyle atlas was calibrated in normal children, a normal child
**by definition** has BA = CA at every age and ΔBA/ΔCA ≡ 1.0. The model uses this
directly.

```
dBA/dCA = Rmat(E2, androgens, IGF-1) / Rmat(normal reference value at the same CA)
```

Two results that an additive formulation does not produce follow automatically from
this.

- An untreated CPP girl has **ΔBA/ΔCA ≈ 1.76 at chronological age 7** — because the
  numerator is pubertal and the denominator prepubertal.
- A **suppressed** girl has **ΔBA/ΔCA ≈ 0.56–0.67** at chronological ages 10–12, i.e.
  **lower** even than the normal prepubertal rate — because the reference child of
  that age is already pubertal, so the denominator has grown. This is why skeletal
  maturation on treatment looks like an "arrest", and there is no arrest term inside
  the model.

**② The pituitary decodes pulse frequency, and a depot agonist does not block the receptor but destroys that code.**

```
S = RS · ( Sendo · ffree  +  AINT · fa )
```

`fa` = agonist occupancy, `ffree` = 1 − fa − fx (the share available to endogenous
pulses), `AINT` = 1.6 (the agonist is a **superagonist** at the receptor), `RS` =
sensitised receptor fraction. Because `AINT > 1`, the first dose **raises** the
stimulus above baseline — this is the flare. Suppression comes afterwards, as `RS`
collapses (`KDES` 0.45/day, `KREC` 0.035/day → at fa = 0.8, `RS`ss ≈ 0.09). A GnRH
**antagonist** enters the same equation as `fx` in place of `fa` and has no intrinsic
activity term: with no flare and no desensitisation, occupancy must be maintained
continuously. Opposite pharmacology at the same receptor comes out of a single
equation.

**③ Growth velocity is an unreliable monitoring marker — quantitatively so.**
Effective suppression drops the growth velocity from 9.6 → 5.1 cm/yr, i.e. **below**
the normal prepubertal value (5.6 cm/yr). This is because the sex steroids and their
GH/IGF-1 amplification are lost at the same time. §8 converts this into
discriminative power: between a child who is adequately suppressed and one who is
not, growth velocity separates by **only 20%** (no discrimination), while basal LH
separates by 345% and ΔBA/ΔCA by 31%.

---

## 2. Natural history — the age of pubertal onset is an output, not an input

There is **one patient-level parameter, `MK0`** (the residual MKRN3 brake). MKRN3 is
an imprinted, paternally expressed ubiquitin ligase that **inhibits** the KNDy pulse
generator, and its loss of function is the commonest cause of monogenic CPP. Change
`MK0` alone and breast development, menarche, bone age and final height all move
together.

| | Normal girl (`MK0` 1.00) | CPP girl (`MK0` 0.40) |
|---|---|---|
| Breast development (Tanner B2) | 10.34 yr | **6.45 yr** |
| Menarche | 12.97 yr | **9.17 yr** |
| Completion of skeletal maturation | 14.75 yr | 11.42 yr |
| **Final adult height** | **163.9 cm** | **154.1 cm** |
| Peak growth velocity | 7.40 cm/yr @ 11.25 yr | 9.55 cm/yr @ 7.38 yr |
| Bone age at chronological age 8 | 7.98 yr (BA−CA −0.02) | 9.37 yr (BA−CA **+1.37**) |
| ΔBA/ΔCA (7–9 yr) | **0.97** | **1.76** |
| LH / E2 / uterine volume at chronological age 8 | 0.04 IU/L / 1.7 pg/mL / 1.37 mL | 1.80 IU/L / 38.4 pg/mL / 7.69 mL |
| IGF-1 at chronological age 8 | 161 ng/mL | 362 ng/mL |
| Cumulative E2 exposure | 517 pg/mL·yr | 711 pg/mL·yr |
| Adult BMD Z | −0.04 | +0.14 |

> **Untreated height loss = 9.8 cm** (154.1 vs 163.9)

The normal child's ΔBA/ΔCA = 0.97 is a **self-consistency check** (structurally it
must be 1.00, and it departs by 3% because the simulated E2 does not match the
prescribed reference curve exactly). Boys: normal completion of skeletal maturation
16.44 yr, final height 176.5 cm, peak growth velocity at 13.08 yr.

---

## 3. The central result — the sign flip (computed, not assumed)

Leuprolide 11.25 mg q12wk for 6 years, varying only the moment of starting:

| Start CA | Start BA | GPRES | Untreated (cm) | Treated (cm) | **Gain (cm)** | GV on treatment | ΔBA/ΔCA on treatment |
|---|---|---|---|---|---|---|---|
| 6.6 | 6.95 | 0.801 | 154.1 | 162.5 | **+8.46** | 5.36 | 1.00 |
| 7.0 | 7.61 | 0.732 | 154.1 | 161.3 | **+7.21** | 5.14 | 0.99 |
| 7.5 | 8.47 | 0.641 | 154.1 | 159.6 | **+5.53** | 4.81 | 0.96 |
| 8.0 | 9.37 | 0.546 | 154.1 | 158.0 | **+3.91** | 4.45 | 0.91 |
| 8.5 | 10.26 | 0.451 | 154.1 | 156.6 | **+2.50** | 4.03 | 0.85 |
| 9.0 | 11.13 | 0.358 | 154.1 | 155.5 | **+1.44** | 3.57 | 0.78 |
| 9.5 | 11.95 | 0.270 | 154.1 | 154.7 | **+0.68** | 3.04 | 0.72 |
| 10.0 | 12.70 | 0.190 | 154.1 | 154.3 | **+0.22** | 2.40 | 0.67 |
| 10.5 | 13.38 | 0.118 | 154.1 | 154.1 | **+0.02** | 1.38 | 0.56 |
| 11.0 | 14.00 | 0.051 | 154.1 | 154.0 | **−0.02** | 0.27 | 0.18 |

Interpolating:

- The point at which the gain falls **below 4 cm**: bone age **9.32 yr**
  (chronological age 7.97 yr)
- **Below 2 cm**: bone age **10.67 yr** (chronological age 8.74 yr)
- **Below 1 cm**: bone age **11.60 yr** (chronological age 9.29 yr)
- **Zero**: bone age **13.68 yr** (chronological age 10.75 yr)

The international consensus statement (PMID 19332438) states that there is "no
height gain when started after a bone age of 12 years". The model reproduces that
threshold by computing it independently as a bone age of **11.6 years**.

**Note, however, that the sign does not actually flip negative (−).** The gain
**approaches** zero asymptotically (−0.02 cm at a bone age of 14). That is, in this
structure late treatment is not *harmful* but *pointless*. The intuition that "late
treatment actually loses height because of the fall in growth velocity" is not
supported by the model — because once the growth plate is nearly closed there is no
growth left to suppress.

---

## 4. What is flare actually worth?

At the first dose the agonist raises the stimulus above baseline. That price was
converted into centimetres (leuprolide 3.75 mg q28d vs a hypothetical GnRH
antagonist 18 mg q28d, both started at chronological age 7.5).

| | Agonist | Antagonist | Untreated |
|---|---|---|---|
| Peak LH in the first 14 days (IU/L) | **3.71** | 1.46 | 1.49 |
| Peak E2 in the first 14 days (pg/mL) | **45.6** | 33.2 | 33.7 |
| LH at day 7 | 0.137 | 0.039 | 1.472 |
| LH at day 28 | 0.075 | 0.072 | 1.514 |
| Bone age at chronological age 9 | 9.94 | 9.92 | 11.13 |
| Final adult height (cm) | 160.1 | 160.3 | 154.1 |

> **The price of flare = 0.22 cm of adult height.**

So the flare is a clinically visible phenomenon (a 2.5-fold LH surge, withdrawal
bleeding) but has **essentially no effect on the height outcome.** It is something to
counsel about; it is not a reason to switch to an antagonist for the sake of height.

**And the model refutes one piece of received wisdom about flare.** In a virtual
cohort of 250, the incidence of withdrawal bleeding after the first depot comes out
at **7.6%** (19/250, reported 5–10%), but treating **the same cohort with a
flare-free antagonist still gives 6.4%** (16/250). In this structure bleeding after
the first depot is not the flare itself but **the abrupt withdrawal of oestrogen from
a primed endometrium**, and so the expectation that switching to an antagonist would
avoid withdrawal bleeding is not supported.

---

## 5. Formulation — the deciding variable is not potency but trough coverage

Five years from chronological age 7.4, the same patient:

| Formulation | Mean Cp (ng/mL) | Mean occupancy fa | Time with LH>0.5 (%) | Time with E2>10 (%) | Final height | Gain |
|---|---|---|---|---|---|---|
| Leuprolide 3.75 mg q28d | 0.960 | 0.703 | 0.0 | 0.0 | 160.4 | +6.39 |
| Leuprolide 7.5 mg q28d | 1.920 | 0.822 | 0.0 | 0.0 | 160.5 | +6.46 |
| Leuprolide 11.25 mg q12wk | 0.960 | 0.577 | 1.0 | 1.9 | 160.0 | +5.99 |
| Leuprolide 30 mg q12wk | 2.559 | 0.749 | 0.0 | 0.0 | 160.4 | +6.39 |
| Triptorelin 11.25 mg q12wk | 1.069 | 0.767 | 0.0 | 0.0 | 160.4 | +6.39 |
| Triptorelin 22.5 mg q24wk | 1.052 | 0.711 | 0.0 | 0.0 | 160.4 | +6.35 |
| Histrelin implant q12mo | 0.500 | 0.806 | 0.0 | 0.0 | 160.3 | +6.28 |
| Nafarelin nasal spray (adherence 1.00) | 0.163 | 0.226 | 0.0 | 0.0 | **158.8** | **+4.72** |
| Nafarelin nasal spray (adherence 0.80) | 0.130 | 0.186 | **14.2** | **50.3** | **158.1** | **+4.05** |
| 3.75 mg every 42 days (delayed) | 0.641 | 0.591 | 0.0 | 0.0 | 160.3 | +6.24 |
| 3.75 mg every 56 days (delayed) | 0.480 | 0.498 | 0.0 | 0.0 | 160.0 | +5.93 |
| 1.875 mg q28d (low dose) | 0.480 | 0.549 | 0.0 | 0.0 | 160.2 | +6.19 |
| GnRH antagonist 18 mg q28d | — | — | 0.0 | 0.0 | 160.7 | +6.62 |

(The Cp in the antagonist row is 0 because that column is the agonist total — the
antagonist is tracked as a separate state variable.)

**The computed conclusion:** depot formulations given correctly are effectively
indistinguishable from one another, and **the therapeutic margin is very wide** —
delaying the injection to twice the interval, or halving the dose, only drops the
gain from 6.4 → 5.9 cm. The only formulation that genuinely fails to suppress is the
**nasal spray** (systemic bioavailability 2.1%, half-life 4 hours).

The PK profile of the histrelin implant has **Cmax = Cmin = 0.500 ng/mL** — with
zero-order release there is no fluctuation at all. This is the reason the formulation
exists.

**Adherence sweep (nasal spray):**

| Adherence | Time with LH>0.5 (%) | Time with E2>10 (%) | Final height | Gain |
|---|---|---|---|---|
| 1.00 | 0.0 | 0.0 | 158.8 | +4.72 |
| 0.90 | 0.8 | 29.7 | 158.5 | +4.40 |
| 0.80 | 13.0 | 50.5 | 158.1 | +4.06 |
| 0.70 | 38.0 | 63.9 | 157.7 | +3.66 |
| 0.60 | 59.6 | 76.9 | 157.3 | +3.24 |
| 0.50 | 78.0 | 88.0 | 156.8 | +2.71 |

---

## 6. GnRH-independent disease — when the target and the mechanism are mismatched the drug effect is zero

In McCune-Albright syndrome (mosaic activating GNAS) the gonad makes E2 by itself,
**without LH**. In the model this term (`AUTSET`) does not exist downstream of the
GnRH receptor, so an agonist cannot structurally reach it.

| Strategy | Mean E2 (7–11 yr) | Bone age at chronological age 10 | Final height | Gain |
|---|---|---|---|---|
| No treatment | 50.6 | 13.88 | 151.9 | +0.00 |
| **Leuprolide 11.25 q12wk** | 42.1 | 13.82 | 152.0 | **+0.03** |
| High-potency AI (letrozole-type) | 2.6 | 10.67 | 160.6 | **+8.63** |
| AI + leuprolide | 2.2 | 10.66 | 160.8 | +8.88 |
| AI + tamoxifen | 2.6 | 10.62 | 161.4 | **+9.45** |

The identical leuprolide regimen gives **+9.78 cm in genuinely central disease**.
Same drug, same dose, same age — and the result is 0.03 cm versus 9.78 cm. This is
how a QSP model states quantitatively whether "the mechanism of action matches the
target".

---

## 7. Boys — aromatase inhibition separates the two signs

In the model a boy's E2 comes **only from the peripheral aromatisation of
testosterone** (the ovarian aromatase term is blocked by sex). An AI therefore
removes, in principle, only the (−) arm and leaves virilisation intact.

| Strategy | Mean E2 | Mean T | Bone age at chronological age 11 | Final height | Gain | BMD Z nadir |
|---|---|---|---|---|---|---|
| Untreated CPP boy | 29.7 | 384 | 14.03 | 167.4 | +0.00 | +0.00 |
| Leuprolide 11.25 q12wk | 3.9 | **35.9** | 12.05 | 170.2 | +2.88 | −0.12 |
| Anastrozole-type AI alone | 4.0 | **402** | 12.40 | 170.9 | +3.54 | −0.12 |
| GnRHa + AI | 1.1 | 36.1 | 11.86 | 171.5 | +4.17 | −0.24 |

The mechanistic separation is reproduced cleanly: the AI alone drops E2 from
29.7 → 4.0 while **holding** testosterone at 402 (virilisation preserved). The
reference final height for a normal boy is 176.5 cm.

**This result does, however, disagree with the literature** — see §12.

---

## 8. Monitoring — which marker actually discriminates

Comparing an adequately suppressed child (leuprolide 11.25 q12wk) with an
inadequately suppressed one (nasal spray, adherence 0.55) every 90 days from 4 months
of treatment onwards:

| Marker | Adequate | Inadequate | Difference | Discriminates? |
|---|---|---|---|---|
| Basal LH (IU/L) | 0.14 | 0.61 | **345%** | ✅ |
| E2 (pg/mL) | 3.63 | 15.20 | **319%** | ✅ |
| Uterine volume (mL) | 1.75 | 4.33 | **148%** | ✅ |
| Tanner breast stage | 1.34 | 2.36 | 77% | ✅ |
| IGF-1 (ng/mL) | 167 | 238 | 42% | ✅ |
| ΔBA/ΔCA | 0.80 | 1.05 | 31% | ✅ |
| **Growth velocity (cm/yr)** | **4.39** | **5.27** | **20%** | ❌ |

Growth velocity has **the lowest discriminative power of every marker**, and on top
of that **its sign is the opposite of the intuition** (the better-suppressed child
grows more slowly). The actual difference in final height between the two children is
3.1 cm.

**And the model refutes the stronger claim I set out to verify.** In a virtual cohort
of 185, the correlation between growth velocity on treatment and the realised height
gain is **r = +0.997**, i.e. a strong **positive** correlation. The reason is worth
learning: **between** individuals, the two values move together because a third
variable (the growth-plate reserve at the moment of starting) drives both. The
monitoring trap is therefore purely a **within**-individual problem — the step-like
drop at the instant treatment starts (9.6 → 5.1 cm/yr) — and not a problem of ranking
patients against one another. This result is reported as it stands rather than tuned
away.

---

## 9. The slowly progressive variant — the quantitative form of the overtreatment problem

Raise `MK0` and puberty starts later and progresses less urgently. These children also
present with "early breast development".

| `MK0` | Breast development | Untreated menarche | Untreated final height | Treated final height | **Gain** |
|---|---|---|---|---|---|
| 0.34 | 5.81 yr | 8.50 | 152.8 | 161.5 | **+8.72** |
| 0.40 | 6.45 yr | 9.17 | 154.1 | 161.5 | **+7.44** |
| 0.46 | 7.05 yr | 9.75 | 155.6 | 161.7 | +6.10 |
| 0.52 | 7.57 yr | 10.25 | 157.1 | 161.9 | +4.82 |
| 0.58 | 8.03 yr | 10.71 | 158.4 | 162.1 | +3.65 |
| 0.66 | 8.58 yr | 11.24 | 160.0 | 162.4 | **+2.35** |
| 0.76 | 9.18 yr | 11.83 | 161.6 | 162.8 | **+1.17** |

What is gained by treating a child whose breast development began at 8.6 years is
**2.35 cm**, and at 9.2 years **1.17 cm**. The price is 3–5 years of injections, a
fall in BMD Z, hot flushes, and a psychological burden reaching 5.8/10.

---

## 10. A virtual cohort of 400 — the actual yield of a "treat everyone who presents early" policy

Four hundred children with randomised `MK0`, growth parameters, height at age 5 and
**referral delay (0.3–3.8 years)** were screened, and the 198 who presented with
breast development before age 8 were given leuprolide 11.25 mg q12wk for 5 years.
Without the referral delay every child is diagnosed at a young bone age and the
policy comparison becomes meaningless — in real practice the bone age at the moment
the decision is made is spread from 6.9 to 13.5 years (mean 10.47).

| | Value |
|---|---|
| Untreated final height | mean 154.8 cm (SD 4.2, range 144.1–167.8) |
| Treated final height | mean 157.5 cm (SD 4.3) |
| **Mean gain if everyone is treated** | **+2.61 cm** (SD 2.45, range −0.11 ~ +9.29) |
| Proportion with a gain below 1 cm | **37.9%** |
| Proportion with a gain below 0 | 9.6% |

| Policy | Proportion treated | Cohort mean gain | Mean gain among the treated children |
|---|---|---|---|
| Treat only if bone age < 10.0 | 39.9% | +2.08 cm | **+5.21 cm** |
| Bone age < 10.5 | 48.0% | +2.27 cm | +4.73 cm |
| Bone age < 11.0 | 56.1% | +2.42 cm | +4.31 cm |
| Bone age < 11.5 | 66.7% | +2.53 cm | +3.79 cm |
| Bone age < 12.0 | 73.2% | +2.57 cm | +3.51 cm |
| **Treat everyone** | 100.0% | +2.61 cm | +2.61 cm |

How to read it: **treat everyone and the total cohort gain is maximised (+2.61 cm),
while the yield per child who receives the injections is minimised (+2.61 cm).**
Restrict to a bone age of 10 and the number treated falls to 40% while the cohort
gain retains only +2.08 cm (80% of the maximum), but the gain per treated child
**doubles (+5.21 cm)**. This, then, is not a question of efficacy but an allocation
question: **on whom do we impose 3–5 years of injections.**

The correlation between bone age at diagnosis and the gain is **r = −0.973**, making
it effectively the only variable that predicts the gain in this cohort (age at breast
development gives r = −0.428, and the predicted untreated final height r = −0.267).

---

## 11. Other computed results

**Bone density — the "fall" is mostly regression of an advanced baseline.**

| Chronological age | Normal | CPP untreated | CPP + GnRHa | + Ca/VitD |
|---|---|---|---|---|
| 7.4 | −0.00 | +0.28 | +0.28 | +0.28 |
| 9.0 | −0.03 | +0.53 | +0.15 | +0.38 |
| 10.0 | −0.04 | +0.59 | +0.05 | +0.37 |
| 11.4 | −0.00 | +0.53 | **−0.18** | +0.22 |
| 14.0 | +0.01 | +0.33 | −0.07 | +0.15 |
| 20.4 | −0.04 | +0.14 | **−0.07** | +0.05 |

Untreated CPP has a bone density that is **raised** for its age (+0.59 SD at 10
years), and once treatment starts that advanced baseline regresses towards the
reference line. The nadir on treatment is **−0.21 SD**, and adult peak bone mass
recovers to **−0.08 SD**. Give calcium and vitamin D alongside and the nadir becomes
**+0.00 SD**, i.e. the fall itself disappears — the quantitative form of the finding
reported by Antoniazzi et al. (PMID 10372699) that "calcium supplementation prevents
bone demineralisation".

**Axis recovery (after the last depot):**

| Days elapsed | Cp (ng/mL) | Occupancy fa | RS | LH | E2 |
|---|---|---|---|---|---|
| 0 | 0.776 | 0.689 | 0.091 | 0.07 | 2.8 |
| 30 | 0.198 | 0.361 | 0.143 | 0.19 | 4.4 |
| 60 | 0.051 | 0.126 | 0.254 | 0.77 | 17.0 |
| 90 | 0.013 | 0.036 | 0.398 | 1.76 | 35.8 |
| 180 | 0.000 | 0.001 | 0.591 | 2.94 | 50.5 |

Recovery of the sensitised receptor pool `RS` (t½ 20 days) is **slower than the PK**,
so suppression persists for several weeks beyond the last injection.

**rhGH co-treatment — partly self-cancelling.**

| Start | Strategy | Mean IGF-1 | Mean GV | ΔBA/ΔCA | Final height | Gain |
|---|---|---|---|---|---|---|
| 8.6 yr | GnRHa alone | 168 | 3.75 | 0.77 | 156.4 | +2.38 |
| 8.6 yr | GnRHa + rhGH | 260 | 4.57 | **0.84** | 158.5 | +4.41 |
| 8.6 yr | rhGH alone | 479 | 5.46 | **1.30** | 155.4 | +1.32 |
| 7.2 yr | GnRHa alone | 165 | 4.86 | 0.93 | 160.7 | +6.65 |
| 7.2 yr | GnRHa + rhGH | 251 | 5.92 | 1.01 | 163.2 | +9.12 |
| 7.2 yr | rhGH alone | 460 | 8.78 | **1.71** | 156.3 | +2.22 |

rhGH raises the growth velocity but **accelerates skeletal maturation along with it**
(ΔBA/ΔCA 0.77 → 0.84, and 1.30 when used alone). The net gain is +2.0~+2.5 cm when
laid on top of a GnRHa and stops at +1.3 cm on its own — that is, giving growth
hormone without suppression eats into itself.

**Sensitivity (±20%, on the adult-height gain, top 6):**

| Parameter | −20% | +20% | Sensitivity (cm/10%) |
|---|---|---|---|
| `BAFUS` (bone age at completion of skeletal maturation) | +1.11 | +8.39 | **+1.82** |
| `ME` (E2 drive on skeletal maturation) | +3.33 | +7.58 | +1.06 |
| `M0` (steroid-independent maturation) | +7.36 | +4.10 | −0.82 |
| `PIGF` (IGF-1 exponent) | +6.83 | +4.36 | −0.62 |
| `GVBASE` (growth-velocity scale) | +4.53 | +6.79 | +0.57 |
| `PGP` (reserve exponent) | +6.48 | +4.96 | −0.38 |

What governs the gain is not a drug parameter but the size of the **growth window**
(`BAFUS`) and **how fast E2 closes that window** (`ME`). The sensitivities of
pharmacological parameters such as `EC50L`, `AINT` and `KDES` are 5–25 times smaller
than those. This is what explains the result of §5 (that there is almost no
difference between formulations).

---

## 12. Where the model disagrees with the literature (not tuned away)

1. **The peak growth velocity of normal children is too low.** The model gives
   7.40 cm/yr (11.25 yr) for girls and 7.16 cm/yr (13.08 yr) for boys. Population
   data are about 8.3 and 9.5 cm/yr respectively. The cause is structural: the
   growth-plate reserve exponent (`PGP` = 0.45) that makes the final height come out
   right flattens the pubertal spurt. **Final height and bone age are right; only the
   peakedness of the spurt is too low.**

2. **With an early start, ΔBA/ΔCA is 0.95–1.00, higher than the literature
   (0.4–0.7).** With a late start it falls to 0.56–0.72 and enters the published
   range. The model's interpretation is that the commonly quoted low values were
   **measured in children whose bone age was already advanced**. There is no
   mechanism in this structure that *arrests* skeletal maturation below the
   prepubertal level.

3. **Time from stopping treatment to menarche is 5.6–8.0 months**, shorter than the
   reported figure (about 12–18 months). This is because in the model the primed state
   of the endometrium (`ENDO`) is not completely lost even during suppression.

4. **The Bayley–Pinneau predicted adult height is 5–10 cm higher than the model's
   actual final height** (BP-PAH 165–171 cm in the treated group vs 159.7 cm actual).
   This is not a bug but a prediction, and it points in the same direction as the
   clinical observation that the BP method overestimates in children with an advanced
   bone age.

5. **In boys, an aromatase inhibitor alone gives +3.54 cm.** However, the randomised
   trial by Varimo et al. (PMID 31024444) reported that **letrozole alone did not
   increase adult height**. The model does contain the fact that removing E2 also
   removes the (+) arm (the GH/IGF-1 amplification), but that cancellation is not as
   complete as it is in the actual trial. This failure in this direction is stated
   explicitly.

6. **The referral delay dominates the result.** In the §10 cohort, the proportion with
   a gain below 1 cm (37.9%) depends entirely on how the referral-delay distribution
   (0.3–3.8 years) is set. That distribution was **assumed** rather than computed by
   the model, and it must be replaced with real referral-delay data from a particular
   health system before the policy table in §10 can be applied to that system.

7. **`ENDOBLEED` = 0.503 is the only parameter in the entire model fitted to a
   literature figure** (withdrawal bleeding after the first depot 7.2%, reported
   5–10%). Every other parameter is either a literature value or hand-derived from
   published steady-state observations.

---

## 13. How to reproduce

```bash
# mechanistic map (196 nodes / 280 edges / 20 clusters)
dot -Tsvg cpp_qsp_model.dot -o cpp_qsp_model.svg
dot -Tpng -Gdpi=150 cpp_qsp_model.dot -o cpp_qsp_model.png

# independent cross-check: reproduces every number above without R / mrgsolve (numpy + RK4, dt 0.25 d)
python3 cpp_reference_check.py          # -> cpp_reference_output.txt

# mrgsolve model (44 ODEs, 16 scenarios, 12 analysis functions)
Rscript cpp_mrgsolve_model_en.R

# Shiny dashboard (10 tabs)
Rscript -e 'shiny::runApp("cpp_shiny_app_en.R")'
```

`cpp_reference_check.py` is an independent implementation, with fixed-step RK4, of
**the same 44 ODEs and the same parameter values** as `cpp_mrgsolve_model_en.R`. If the
two implementations do not agree it means one of them contains a transcription error,
and that disagreement is itself the bug report. `cpp_reference_output.txt` is the
Python-side output committed verbatim.

**A note on the numerics:** LH, FSH, E2 and testosterone are deliberately
coarse-grained as **pulse-averaged concentrations** (relaxation time 6–12 hours)
rather than at their real turnover (20–60 minutes). The timescale the model addresses
is days to years, and pulse **frequency** enters explicitly as a separate state
variable (`PULS`). Nothing clinically important is faster than 6 hours: the flare is
set by receptor desensitisation (t½ about 1.5 days) and between-dose escape by depot
dissolution (t½ about 15 days).

---

## 14. Disclaimer

This model is a **quantitative to semi-quantitative QSP model for educational and
research purposes** and has not been fitted to individual patient data. It must not
be used directly for actual clinical decision-making, prescribing or regulatory
submission.
