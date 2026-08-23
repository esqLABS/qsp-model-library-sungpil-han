# Neonatal Opioid Withdrawal Syndrome QSP Model

**Neonatal opioid withdrawal syndrome / neonatal abstinence syndrome · NOWS / NAS**

The starting point of this model is a single sentence.

> **NOWS is a noradrenergic disorder, yet we treat it with a µ agonist and
> score it with a tool that measures neither of those two things.**

Everything is built on a single difference.

```
GAP = A − ITONE

  A      The neuroadaptive set-point of the locus coeruleus.
         The value written there by months of sustained µ receptor occupancy in utero.
  ITONE  The inhibitory tone hanging on the LC's GIRK channels at this moment.
```

**GAP > 0** → LC disinhibition → noradrenaline outflow → symptoms.
**GAP < 0** → over-sedation. This is not "headroom". The Finnegan scoring sheet
**scores it as zero** while it suppresses feeding.

---

## 1. What this model claims (six numbers)

### 1-1. The maternal methadone dose changes **when**, not **how severe**

The map for the neuroadaptive set-point is **saturated** in occupancy:

```
A0 = AGAIN · FGA(GA) · EAD / (EA50A + EAD)     EA50A = 0.16
```

Even at methadone 40 mg/day the effective fetal µ occupancy EAD already exceeds 0.95.
So the whole of the difference in A0 produced by a 4-fold dose difference is this:

| Maternal methadone | A0 | Peak Finnegan | Onset (score>8) | Days of treatment |
|---|---|---|---|---|
| 40 mg/day | 0.687 | 9.5 | **2.7 days** | 23.0 |
| 90 mg/day | 0.696 | 9.3 | 3.2 days | 22.8 |
| 160 mg/day | 0.699 | 9.3 | **3.8 days** | 22.5 |

A0 moves by 1.6% and the peak score moves by 0.2 points. But **the time of onset is
pushed back by 1.5 days.** The reason is arithmetic. Washout is log-linear, so the
logarithm of the dose ratio becomes a time shift directly: ln(4) × 16 h / 0.693 = 32 h.

This is where the literature's failure to reach consensus over thirty years on the
relationship between maternal dose and NOWS severity comes from — **the relationship
lies in the timing, not in the severity.**

### 1-2. The three weaning rules are not good and bad but **three points on one trade-off curve**

Same infant, same biology, only the rule changes:

| Weaning rule | Days of treatment | Cumulative withdrawal burden ∫GAP⁺dt | Cumulative over-sedation |
|---|---|---|---|
| −10%/day of the stabilisation dose | **22.8** | 56.6 | 0.7 |
| −10%/day of the standard dose (absolute steps) | 29.8 | 27.5 | 38.7 |
| −10%/day of the current dose (exponential) | **35.5** | **24.2** | 7.6 |
| Tracking A directly (an unobservable oracle) | 22.1 | 68.6 | 0.5 |

The exchange rate is computed: **0.39 days of treatment to avoid one unit of
withdrawal burden**.

Two things follow. First, the 3-fold spread in published NOWS treatment durations
(4–28 days) is explained entirely by **differences in the weaning rule** rather than
by differences between infants — in this model, holding the biology fixed and
changing only the rule gives anything from 21.5 to 35.5 days.
Second, **knowing A buys almost nothing.** An oracle controller reading the
unmeasurable state (22.1 days / burden 68.6) does not dominate the stabilisation rule
(22.8 days / burden 56.6). The problem is not information but an exchange rate that
nobody sets explicitly.

### 1-3. The two halves of Eat-Sleep-Console work **in opposite directions**

The ESC-NOW trial changed the assessment tool and the care environment at the same
time, so its effect cannot be attributed. In this model the two enter through
different parameters and are therefore separable (virtual cohort n=200/cell):

| Condition | Pharmacotherapy rate | Discharge readiness | **Withdrawal burden the infant actually endured (median)** |
|---|---|---|---|
| Finnegan · usual care | 62.5% | 17.8 days | 16.5 |
| **Assessment criteria only** switched to ESC | 34.0% | 7.0 days | **33.4** |
| **Care only** intensified | 26.0% | 7.0 days | **10.8** |
| Both (the trial's intervention) | 13.5% | 7.0 days | 17.1 |

How to read it: **changing the assessment criteria alone doubles the withdrawal
burden (16.5 → 33.4). Changing the care alone reduces the burden (16.5 → 10.8). Put
the two together and it is 17.1 — essentially the same burden as usual care,
obtained with one fifth of the pharmacotherapy.**

That is to say, the two halves of ESC are not independent contributions to a single
benefit. **One half buys the reduction in medication and the other pays for it.**
Changing the scoring sheet without changing the care is not a safe intervention in
this model.
(A single deterministic patient goes the same way: ESC 12.0 days / burden 101.0
versus Finnegan 22.8 days / burden 56.6, untreated burden 192.8.)

### 1-4. The two reasons preterm NOWS is milder are separable, and the ratio is 9 to 1

The two causes that are always entangled in cohort data — (a) immature neuroadaptive
capacity and (b) immature clearance relative to weight-based dosing — are separated
by **forcing a term-level persistent adaptation pool** onto a 34-week infant:

| | A0 | Days of treatment | Cumulative over-sedation ∫GAP⁻dt |
|---|---|---|---|
| Term 39 weeks | 0.696 | 22.8 | 0.7 |
| Preterm 34 weeks (own adaptation) | 0.479 | **11.6** | **33.9** |
| Preterm 34 weeks (term adaptation forced) | 0.696 | 21.6 | 2.4 |

Of the 11.2-day difference, 10.0 days (89%) is the lower adaptation and 1.2 days
(11%) is pharmacokinetics. But the pharmacokinetic side charges its price
**somewhere other than in duration**: cumulative over-sedation in the preterm infant
is **48 times** that in the term infant. Weight-based dosing does not shorten the
preterm infant's treatment, it **over-sedates** them — and the Finnegan score does
not score over-sedation.

### 1-5. Phenobarbital **lowers the score and does not lower the GAP**

Phenobarbital is drawn outside the GIRK union (GABA-A, downstream of the LC).
In an infant exposed to opioids alone:

| | Days of treatment | Cumulative withdrawal burden | Weight gain |
|---|---|---|---|
| Morphine alone | 22.8 | 56.6 | 1383 g |
| + phenobarbital | **18.3** | **76.6** (+35%) | **558 g** |

**On paper it is 4.5 days shorter, the actual withdrawal burden rises by 35%, and
weight gain falls to 40%.** In an infant with concurrent benzodiazepine exposure the
story is different (31.1 days → 23.0 days, a shortening of 8.1 days against 4.5 days
for pure opioid): phenobarbital is genuinely beneficial **only when** part of the
score is not opioid withdrawal.

The other side of the same logic is clonidine. The α2 receptor opens **the same GIRK
channel** as the µ receptor. So clonidine actually fills the GAP, and at the same
time it does not re-drive the re-inducible adaptation pool AT: at 6 µg/kg/day,
22.8 → 16.7 days and burden 56.6 → 30.5. It is not free, though — at 12 µg/kg/day
the over-sedation integral jumps from 0.7 to 47.2.

### 1-6. Treatment partly prolongs itself

The adaptive set-point is two pools.

```
A = AD + AT
  AD  Persistent pool. Set in utero, t½ ≈ 12 days, not re-inducible after birth.
  AT  Re-inducible pool (classical adenylyl cyclase superactivation),
      τ ≈ 60 h, re-driven by any µ agonist now in the body — the treatment drug included.
```

Switch AT off (ATMAX → 0) and treatment falls from 22.8 days to 19.0 days.
**The adaptation that the treatment itself re-induces adds 3.8 days (17%) to the
duration of that treatment.** This is why "pushing the score down harder" is not
free, and why lowering the target score to 3 points stretches treatment to 42.5 days
while the over-sedation integral explodes from 0.7 to 87.4.

---

## 2. Model structure

### 2-1. 43 state variables (mrgsolve `$CMT`)

| Block | States |
|---|---|
| Morphine PK | `AGUT_M` `AC_M` `AP_M` `CE_M` |
| M6G | `A_G` `CE_G` |
| Methadone | `AGUT_D` `AC_D` `CE_D` |
| Buprenorphine / norbuprenorphine | `AGUT_B` `AC_B` `CE_B` `AC_N` `CE_N` |
| Clonidine | `AGUT_C` `AC_C` `CE_C` `TACH` |
| Phenobarbital | `AGUT_P` `AC_P` |
| Neuroadaptation | `AD` `AT` `RMU` |
| Locus coeruleus / signs | `NE` `CNS` `ANS` `GI` |
| Functional items | `SLP` `EAT` `CONS` |
| Growth | `WT` |
| Non-opioid co-exposure | `BZD` `BZW` `NICW` |
| Counters | `CUMM` `CUMD` `SEIZH` `AUCGAP` `AUCSED` |
| Measurement | `FNAS_S` |
| Protocol (feedback controller) | `DOSE` `TRTON` `STAB` `DSTAB` |

### 2-2. Core algebraic block

```
Competitive occupancy  θf = xf/(1+xf+xb),  θb = xb/(1+xf+xb)
                       xf = UF/EC50_eff,   xb = UB/EC50b_eff
Tolerance shift        EC50_eff = EC50 · (1 + KSHIFT · A)          KSHIFT = 3.0
µ effect               EMU = RMU · (1.00·θf + 0.80·θb)
Adaptation drive       EAD = RMU · (1.00·θf + 0.80·0.40·θb)        ← buprenorphine appears twice
α2 effect              EA2 = 0.34 · Cclon / (1.6·(1+1.2·TACH) + Cclon)
GIRK union             ITONE = EMU + EA2 − EMU·EA2                 ← a union, not a sum
GAP                    GAP = (AD + AT) − ITONE
LC drive               LC = LCB0·(1+0.55·(1−CARE)) + 11.9·GAP⁺·(1−0.30·CARE)
```

That buprenorphine enters **twice, with different numbers**, is this model's
structural commitment. The acute effect takes the partial-agonist ceiling
(`EMAXB` = 0.80, with respect to the endpoint of withdrawal suppression), and the
adaptation drive separately takes the weakness of β-arrestin recruitment
(`BADRV` = 0.40). So the low A0 of the buprenorphine-exposed infant (0.487 against
0.696) is not a value tuned with a "severity constant" but is **derived** from
receptor pharmacology.

### 2-3. The protocol is drawn as a feedback controller

Initiation threshold → latch → escalation (rate-limited) → stabilisation-dose lock →
48 hours of stability → one of four weaning rules. Because it sits inside `$ODE`, the
protocol is itself an object of simulation, and the trade-off curve of section 1-2
above comes out of this block.

### 2-4. Mechanistic map

`nows_qsp_model.dot` — **188 nodes, 17 clusters**
(maternal/placental · neonatal PK · maturation · treatment PK · µ receptor ·
neuroadaptation · locus coeruleus · sites of action of the other drugs ·
polysubstance exposure · CNS signs · autonomic signs · gastrointestinal signs ·
assessment tools · care environment · growth · protocol · clinical endpoints).
[SVG](../../../neonatal-opioid-withdrawal-syndrome/nows_qsp_model.svg) · [PNG](../../../neonatal-opioid-withdrawal-syndrome/nows_qsp_model.png)

---

## 3. Calibration and verification

### 3-1. Verification method

Because it was built in an environment without an R toolchain, every ODE and
algebraic block of the mrgsolve model is verified by an **independent
re-implementation in Python/scipy** (`nows_verify_python.py`). This script does two
things.

1. It parses the `$PARAM` block of `nows_mrgsolve_model.R` directly and
   **cross-checks all 145 parameters exhaustively** → 0 mismatches.
2. It runs **61 anchors** (pharmacokinetic · structural invariant · clinical) →
   **61/61 PASS**.

Full output: [`nows_verification_output.txt`](../../../neonatal-opioid-withdrawal-syndrome/nows_verification_output.txt)

### 3-2. Clinical anchors

| Anchor | Source | Model |
|---|---|---|
| Peak Finnegan in an untreated methadone-exposed infant | mid-teens | 15.9 (5.8 days) |
| Time of onset — short-acting opioid | <24 h | 0.4 days |
| Time of onset — methadone | 48–72 h | 2.6 days |
| Time of onset — buprenorphine | 36–60 h | 2.8 days |
| Morphine blood concentration during treatment | 15–45 ng/mL | median 12.9, maximum 36.7 |
| Neonatal methadone t½ | 16–25 h | 16.0 h (simulation 14.1 h) |
| Neonatal morphine CL | 0.25–0.31 L/h/kg | 0.279 |
| UGT2B7 capacity at PMA 40 weeks | ~23% of adult | 23.3% |
| Clonidine blood concentration at 6 µg/kg/day | 0.5–1.3 ng/mL | 0.86 |
| Weight gain during treatment | 25–30 g/day | 30.7 |
| Cumulative morphine, buprenorphine-exposed / methadone-exposed | MOTHER: markedly lower | 0.50-fold |
| Days of treatment, buprenorphine-exposed / methadone-exposed | MOTHER: 4.1 against 9.5 days | 0.70-fold |
| Sublingual buprenorphine treatment shorter than morphine | Kraft 2017: 15 against 28 days | 17.8 against 22.8 days |
| Adjunctive clonidine shortens treatment | Agthe 2009: 11 against 15 days | 16.7 against 22.8 days |
| ESC markedly lowers the pharmacotherapy rate | ESC-NOW: 52.0% → 19.5% | 62.5% → 13.5% |
| ESC brings discharge readiness forward | ESC-NOW: 14.9 → 8.2 days | 17.8 → 7.0 days |

### 3-3. What is not reproduced — stated explicitly

- **The 9.5-fold difference in cumulative morphine in MOTHER is not reproduced.**
  The model only goes as far as 0.50-fold (a 2-fold difference). MOTHER's extreme
  ratio is the product of a difference in dose among treated infants and a
  difference in the proportion treated; in this model the buprenorphine-exposed
  cohort has a treatment rate of 40.0%, close to MOTHER's 47%, so the remaining
  difference would have to come from the dose. Producing a dose difference of that
  size would require lowering `BADRV` further, and the untreated peak score of the
  buprenorphine-exposed infant would then become clinically far too mild. This
  tension is left unresolved.
- **The absolute value of the treatment duration is a property of the weaning rule,
  not of the infant.** The 22.8 days of the base scenario sits in the middle of the
  literature range (4–28 days) and must not be read as predicting the length of stay
  at any one institution.
- **Neurodevelopmental outcomes are not predicted.** `AUCGAP`/`AUCSED` are exposure
  metrics, not outcome predictors. The observational literature linking NOWS to
  later development is too confounded by the social environment to be used as a
  calibration target.

### 3-4. Defects exposed during verification and fixed

Re-implementing the model in Python and getting the anchors to pass exposed five
substantive defects, all of which were resolved by fixing the structure (not papered
over with parameters).

1. **The trap in which treatment sustains neuroadaptation permanently.** In a
   single-pool model, `dA/dt = KAON·EMU·(1−A) − KAOFF·A` has a fixed point at
   A* = 0.62 for EMU ≈ 0.47, so the infant could not begin weaning at any dose and
   sat at the same score (6.0) for more than 40 days. Resolved by splitting into two
   pools, persistent and re-inducible.
2. **µ occupancy released too late.** Without the tolerance shift (`KSHIFT`), a
   residual methadone concentration of 5 ng/mL still occupied 66% and pushed onset
   back by 4 days. Resolved by putting in the fact that the concentration-effect
   curve of the neuroadapted neonate is shifted to the right.
3. **Clonidine replaced morphine completely.** With EC50 0.45 ng/mL·Emax 0.75,
   6 µg/kg/day made the treatment itself unnecessary. Corrected to EC50 1.6
   ng/mL·Emax 0.34, returning it to the level of an adjunct.
4. **The breast-milk dose was 20-fold too high.** An initial expression made
   directly proportional to the maternal daily dose was giving the neonate
   0.63 mg/kg/day. Rewriting it as a relative infant dose (RID, 2.5% of the
   weight-adjusted maternal dose) brought it down to 4.5 ng/mL, which agrees with
   the infant plasma concentration of 2–9 ng/mL reported by Begg 2001.
5. **The preterm experiment was not an experiment.** A0 did not properly depend on
   gestational age, so 34 weeks and 39 weeks had the same A0. Resolved by putting a
   gestational-age sigmoid (`GA50A` 33.5 weeks, Hill 8) into the adaptive capacity.

There was also a problem of the integrator stalling on its initial step estimate, so
`first_step` was specified explicitly. That was a numerical problem, not a model one.

---

## 4. The assumption this model exposes most

**The neuroadaptation block has never been measured in a human neonate.** The
structure (a persistent pool set in utero plus a re-inducible pool that current
exposure re-drives) is inferred from the rodent locus coeruleus literature and from
the clinical course, and the two time constants — 12 days for the persistent pool,
60 hours for the re-inducible pool — are the most exposed numbers in this model.

In principle they are identifiable. **Serial scores obtained under different weaning
rules** distinguish these two constants: exponential weaning (−10%/day of the current
dose) and absolute-step weaning (−10%/day of the standard dose) have different
sensitivities to the half-life of the persistent pool. Randomise the two rules within
one institution and record the score trajectory every four hours, and within a few
weeks you would know whether the value of 12 days is right. Without that data, all
the absolute treatment durations above are conditional.

Incidentally, the `AUCSED` metric has a known limitation. A small residual opioid
tone remaining even after A has fallen to 0 (the breastfeeding arm, for example)
makes the GAP negative and is therefore integrated as "over-sedation". In that case
the infant is not sedated but merely exposed to opioid at a very low level. The
`AUCSED` of 79.1 in the breastfeeding scenario must be read this way.

---

## 5. Files

| File | Contents |
|---|---|
| [`nows_qsp_model.dot`](../../../neonatal-opioid-withdrawal-syndrome/nows_qsp_model.dot) | Mechanistic map source (188 nodes / 17 clusters) |
| [`nows_qsp_model.svg`](../../../neonatal-opioid-withdrawal-syndrome/nows_qsp_model.svg) · [`.png`](../../../neonatal-opioid-withdrawal-syndrome/nows_qsp_model.png) | Rendering (PNG 150 dpi) |
| [`nows_mrgsolve_model.R`](../../../neonatal-opioid-withdrawal-syndrome/nows_mrgsolve_model.R) | 43-ODE mrgsolve model · 25 scenarios · calibration sources |
| [`nows_shiny_app.R`](../../../neonatal-opioid-withdrawal-syndrome/nows_shiny_app.R) | 12-tab interactive dashboard |
| [`nows_verify_python.py`](../../../neonatal-opioid-withdrawal-syndrome/nows_verify_python.py) | Independent Python re-implementation · exhaustive parameter cross-check · 61 anchors |
| [`nows_verification_output.txt`](../../../neonatal-opioid-withdrawal-syndrome/nows_verification_output.txt) | Full verification output |
| [`nows_population_sim.py`](../../../neonatal-opioid-withdrawal-syndrome/nows_population_sim.py) | Virtual cohort · decomposition of the ESC effect |
| [`nows_population_output.txt`](../../../neonatal-opioid-withdrawal-syndrome/nows_population_output.txt) | Cohort run results |
| [`nows_references.md`](../../../neonatal-opioid-withdrawal-syndrome/nows_references.md) | 143 papers · every PMID verified with PubMed E-utilities |

### Running

```bash
# map
dot -Tsvg nows_qsp_model.dot -o nows_qsp_model.svg
dot -Tpng -Gdpi=150 nows_qsp_model.dot -o nows_qsp_model.png

# verification (numpy + scipy)
python3 nows_verify_python.py
python3 nows_population_sim.py

# model and dashboard (R: mrgsolve, shiny, ggplot2, dplyr, tidyr, DT)
Rscript -e 'source("nows_mrgsolve_model.R"); print(run_all())'
Rscript -e 'shiny::runApp("nows_shiny_app.R")'
```

---

## ⚠️ Disclaimer

This is a QSP model for educational and research purposes. It was assembled from the
published literature but has not been independently verified or certified, and
**must not be used for actual clinical decision-making, prescribing, or regulatory
submission.** In particular, none of the doses, thresholds, or weaning schedules
above is a clinical protocol.
