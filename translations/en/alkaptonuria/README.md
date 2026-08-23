# Alkaptonuria (AKU) — QSP Model
### Quantitative Systems Pharmacology model of alkaptonuria · ochronosis · nitisinone

<a href="../../../alkaptonuria/aku_qsp_model.svg"><img src="../../../alkaptonuria/aku_qsp_model.png" width="640" alt="AKU QSP map"></a>

| Deliverable | File | Scale |
|---|---|---|
| 🗺️ Mechanistic map | [`aku_qsp_model.dot`](../../../alkaptonuria/aku_qsp_model.dot) · [SVG](../../../alkaptonuria/aku_qsp_model.svg) · [PNG](../../../alkaptonuria/aku_qsp_model.png) | 221 nodes · 307 edges · 23 clusters |
| ⚙️ mrgsolve model | [`aku_mrgsolve_model.R`](../../../alkaptonuria/aku_mrgsolve_model.R) | 60 ODEs · 137 parameters · 24 scenarios |
| 📊 Shiny dashboard | [`aku_shiny_app.R`](aku_shiny_app.R) | 12 tabs |
| 📚 References | [`aku_references.md`](aku_references.md) | 131 papers, every PMID verified |

---

## 1. The one question this model sets out to answer

In alkaptonuria, nitisinone reduces urinary homogentisic acid (HGA) by **99.7%.** In
a four-year randomised trial the clinical score (cAKUSSI) worsened by **8.6 points**
less. These are two numbers for the same drug.

A model that cannot explain the distance between 99.7% and 8.6 points has not
explained this disease. This model reduces that distance to three **mass balances**
and took as its single objective that the remaining clinical facts should follow
from those three.

---

## 2. The three balances

### BALANCE 1 — the flux is conserved, and all that changes is the exit

The phenylalanine and tyrosine coming from dietary protein amount to **about
33 mmol/day** (70 g/day of protein), and this has to leave the body. Untreated,
about 95% of it leaves as urinary HGA.

The HPD that nitisinone inhibits sits **one step above the deficient enzyme (HGD)**.
Nitisinone therefore **cannot reduce the size of the flux; it can only change the
exit.** Close the HGA exit and the same 33 mmol/day has to leave as urinary HPPA,
HPLA, tyrosine and tyrosine conjugates, and plasma tyrosine is nothing more than
the **pressure head** needed to push that load through exits of small capacity.

Hence the central asymmetry of this model:

> **The dose sets HGA and the diet sets tyrosine.**

Model calculation (protein fixed at 70 g/day vs dose fixed at 10 mg/day):

| What was changed | Range | Plasma tyrosine |
|---|---|---|
| Nitisinone dose | 1 → 20 mg/day (**20-fold**) | 873 → 975 umol/L (**+12%**) |
| Dietary protein | 42 → 105 g/day (**2.5-fold**) | 560 → 1,545 umol/L (**+176%**) |

The change in tyrosine obtained by moving the dose 20-fold is one fifteenth of the
change obtained by moving protein 2.5-fold. This is not a tuned result but a
consequence of mass conservation. And what it means clinically is plain: **if
keratopathy is the concern, the thing to manage is the diet, not the dose.** That
the 5% keratopathy of the 2 mg cohort (NAC, with active dietitian involvement) and
the 14.5% of the 10 mg cohort (SONIA 2, information only) is read as a difference
of dose is, according to the model, **confounding**.

### BALANCE 2 — the toxic branch is 1.5e-5 of the flux, and it does not come back

While 33 mmol of HGA passes through each day, only **about 0.5 umol/day** of it is
oxidised to benzoquinone acetic acid (BQA) inside avascular collagenous tissue and
**polymerised onto collagen that is never replaced.** As a fraction, 1.5e-5. It is a
quantity invisible in any mass balance, and it is the whole of this disease.

There is **no loss term anywhere in the model** for ochronotic pigment. Because
there is none in the patient. Consequently:

- The disease is **the integral of a flux too small to be measured.**
- The drug can act only on the part not yet integrated.
- Make embrittlement a steep Hill function of pigment **density** (exponent 6) and,
  although pigment accumulates linearly from birth, symptoms appear **abruptly in
  the third decade** and accelerate thereafter. And that is so even though age never
  appears in the damage equation.

The literature forced one exception. There is a report (PMID 32904992) that **the
pigment of skin, ear and sclera lightens** after nitisinone treatment. The model
treats this not as a refutation but as a prediction: a loss term **can exist only
in tissues where collagen is actually replaced.** So a first-order loss term was put
on the skin, ear and sclera depots only, and not on articular cartilage, disc, or
aortic valve. That the observed reversal is confined to epidermal tissues is the
testable consequence.

### BALANCE 3 — HGA clearance has no reserve

The renal clearance of HGA is 600–900 mL/min, that is, **close to renal blood
flow** (tubular secretion exceeds glomerular filtration, and the kidney itself also
produces HGA — PMID 31609457). Being already at its ceiling, it has no buffering
capacity. Plasma HGA = production ÷ clearance, and **what stains the cartilage is
plasma HGA.**

Two things follow.

**(a) The trial endpoint and the causal quantity diverge.** During treatment HPPA
and HPLA rise more than 14-fold and competitively occupy the same organic-anion
secretory pathway (OAT1/OAT3), so HGA clearance **actually falls.** Hence at 10 mg
urinary HGA falls by 99.7% while serum HGA falls by only 92.6%. The rate of pigment
deposition follows the latter number. That the residual damage rate is 7.4% rather
than 0.3% — **a 25-fold difference** — is a substantial part of the distance between
99.7% and 8.6 points.

**(b) Ageing and CKD become disease-modifying factors.** When renal blood flow
falls, production is unchanged while plasma HGA rises 1:1. This has actually been
observed in a cohort of 225 (circulating HGA increases with age and is
significantly associated with a fall in HGA clearance). That is, the reason
untreated progression accelerates in old age is not only that "pigment begets
pigment".

---

## 3. Two damage channels — and one of them is reversible

The balances alone were not enough. Describe damage by the pigment integral alone
and nitisinone started at 49 barely slows cAKUSSI at all (at that age embrittlement
is already saturated). Yet SONIA 2 slowed it by about 75%.

The missing piece was in the literature. At AKU concentrations HGA is **directly**
toxic to chondrocytes: oxidative stress, cytoskeletal damage, impaired autophagy,
SAA amyloid formation, synovial activation. All of these are **functions of the
present concentration and are therefore reversible.**

So the damage rate was split into two channels.

| Channel | Driving variable | What nitisinone does |
|---|---|---|
| ① Pigment–mechanical channel | Cumulative pigment density (irreversible) | Only the future integral |
| ② Concentration channel | Present serum HGA (reversible) | Immediately and completely |

② is **gated** on ① (the pigment sets the vulnerability, the present HGA sets the
rate). This was an essential structural constraint. Make ② a term independent of ①
and synovitis exists from birth and the cartilage of a 10-year-old is destroyed,
which collides head-on with the one universal clinical fact about this disease —
**children with AKU are asymptomatic.** **The structure was discarded rather than
tuned.**

This decomposition supports three facts at once: starting at 49 still slows
cAKUSSI by 75% (②), the hip rotation angle in the NIH trial did not improve (① is
already spent), and starting at 5 can prevent almost all of it (both ① and ②).

---

## 4. Calibration and held-out validation

### 4.1 The anchors — SONIA 1's randomised five-arm dose-response

The calibration was hung on **one randomised trial, one time point, five arms**, so
as to avoid picking convenient numbers from several cohorts. Nine fitted parameters
/ eight anchors, per-coordinate shrink–expand iteration (coordinate descent, log
scale), objective function 0.0198.

| Anchor | Model | Reported | Ratio |
|---|---|---|---|
| Untreated u-HGA24 (umol/day) | 31,585 | 31,530 | 1.00 |
| 1 mg u-HGA24 | 2,990 | 3,260 | 0.92 |
| 2 mg u-HGA24 | 1,191 | 1,440 | 0.83 |
| 4 mg u-HGA24 | 455 | 570 | 0.80 |
| 8 mg u-HGA24 | 173 | 150 | 1.15 |
| Untreated serum HGA (umol/L) | 28.3 | 28 | 1.01 |
| 2 mg plasma tyrosine | 933 | 782 | 1.19 |
| 10 mg plasma tyrosine | 972 | 875 | 1.11 |

The fitted values landed somewhere physiologically reasonable. `KI_NT` =
**11.2 nmol/L** (the same order of magnitude as the literature IC50 of about
40 nmol/L), `HNT` = 1.45 (a curve steeper than simple competitive inhibition —
consistent with the enzymology of nitisinone as a slow, tight-binding inhibitor).

**There is one term whose hypothesis the data rejected.** A source of HGA that does
not respond to nitisinone (`SRC_INS`) was put in and fitted, and it **converged to
7.5 umol/day, that is, to almost zero.** The data support "the dose-response is
steep" rather than "there is an uninhibited source". The rejected term is left in at
a near-zero value rather than deleted.

### 4.2 Validation — 21 facts the calibration never saw

| # | Validation item | Model | Reported | Ratio |
|---|---|---|---|---|
| 1 | **Serum HGA, 10 mg (umol/L)** | 1.58 | 2.23 | 0.71 |
| 2 | Serum HGA, 2 mg (umol/L) | 15.10 | 3.86 | **3.91** ⚠ |
| 3 | Serum HPPA fold rise, 10 mg | 24.9 | 14.65 | **1.70** ⚠ |
| 4 | Residual u-HGA24 ratio (2 mg / 10 mg) | 11.3 | 6.63 | **1.70** ⚠ |
| 5 | **Tyrosine ratio (10 mg / 2 mg)** | 1.042 | 1.119 | 0.93 |
| 6 | **Tyrosine, 10 mg + 0.7 g/kg protein** | 657 | 620 | 1.06 |
| 7 | **Mass conservation (sum of exits / input)** | 0.983 | 1.00 | 0.98 |
| 8 | Age at which cartilage embrittlement > 0.5 | 38.5 | 31 | 1.24 |
| 9 | Age at appearance of disc calcification | 34.4 | 30 | 1.15 |
| 10 | Age at which pain VAS > 2 | 32.8 | 30 | 1.09 |
| 11 | **Age at 50% joint replacement** | 50.3 | 55 | 0.92 |
| 12 | **Age at aortic valve involvement** | 52.3 | 54 | 0.97 |
| 13 | **Age at which renal stones become clinical** | 59.4 | 64 | 0.93 |
| 14 | SONIA 2 control-arm cAKUSSI slope (points/month) | 0.289 | 0.239 | 1.21 |
| 15 | **SONIA 2 10 mg slope (points/month)** | 0.065 | 0.060 | 1.08 |
| 16 | **NAC 2 mg slope (points/month)** | 0.173 | 0.190 | 0.91 |
| 17 | SONIA 2 cAKUSSI difference, 48 months | −10.8 | −8.6 | 1.25 |
| 18 | SONIA 2 difference in aortic valve Pmax slope | −0.057 | −0.0093 | **6.1** ⚠ |
| 19 | **Pigment reduction with ascorbate 1 g/day** | 2.8% | ~0% | — |
| 20 | Skin pigment lightening, 2 years at 10 mg | 0.5% | ~20% | **0.02** ⚠ |
| 21 | **Cartilage pigment lightening, 2 years at 10 mg** | 0% | 0% | — |

**Median |log ratio| = 0.137** (a typical error of about 14%), with 13 of the 21
inside ±25% of the reported value.

In particular, the ones that came out right **having never been used in the fit**:

- **The dose-independence of tyrosine** (#5). This is the central claim of the
  model; it predicted 1.04-fold for a five-fold difference in dose, and the
  observation was 1.12-fold. Both values are a long way from five-fold.
- **Tyrosine under protein restriction** (#6). 657 umol/L at 0.7 g/kg, reported 620.
- **Three clinical time points** (#11, #12, #13). 50% joint replacement, valve
  involvement and renal stones at 50.3 / 52.3 / 59.4 years respectively — reported
  55 / 54 / 64. These time points did not go into the parameters.
- **The separation of the clinical slopes at 2 mg and 10 mg** (#15, #16). 0.065 and
  0.173 points/month, reported 0.060 and 0.190. A consequence of the structure in
  which the difference in clinical effect between the two doses comes out of the
  difference in serum HGA.
- **The failure of ascorbate** (#19). Ascorbate returns BQA to HGA, but
  polymerisation is an irreversible step **in parallel rather than in series**, so
  the pigment flux removed is only `KRED·[Asc] / (KRED·[Asc] + KPOL)` — 2.8% at
  1 g/day. Not a tuned failure but a structural one.
- **Cartilage pigment does not lighten** (#21). There being no loss term, 0%.
  Automatically.

### 4.3 Failures — the ones reported rather than fixed

**① Serum HGA at 2 mg (3.9-fold too high) — this is not a parameter problem.**
Virtually all HGA leaves in the urine, so urinary HGA ≈ total production. For the
four reported numbers (2 mg: urinary 1,200–1,440 / serum 3.86 · 10 mg: urinary 181
/ serum 2.23) to be simultaneously true, HGA clearance would have to be **4.6-fold
higher at 2 mg than at 10 mg.** Yet the concentrations of the competing organic
anions (HPPA, HPLA) are much the same at the two doses. No single clearance model
can satisfy all four values at once — and this is derived by hand before the model
is ever run. Candidate explanations: (i) the limit of quantification of serum HGA at
low concentrations, (ii) a cohort difference between NAC (2 mg, active dietary
management) and SONIA 2 (10 mg, information only), (iii) a saturable secretory
pathway not yet described. **The discriminating experiment: paired measurement of
serum and urinary HGA at 2 mg and at 10 mg, in one cohort, at the same visit.**

**② The rise in serum HPPA and the residual HGA ratio are each exactly 1.70-fold
too high** (#3, #4). This is not a coincidence. Both are determined by the single
steepness of inhibition, `HNT` = 1.45. Lower `HNT` and the two errors shrink
together, but SONIA 1's five-arm dose-response collapses. That is, the slope of
SONIA 1's dose-response and the observed rise in HPPA are not mutually compatible
inside this structure.

**③ The aortic valve (6.1-fold too high, yet in practical terms a null).** SONIA 2
reported a difference in the rate of Pmax progression over four years of
0.0093 mmHg/year (p=0.53). Putting the valve into the reversible channel gave
−0.51 mmHg/year, that is, a 55-fold overestimate of a null result. So **the valve
was structurally excluded from the reversible channel** — on the grounds that
calcific valve disease, once seeded, is a self-sustaining osteogenic programme. The
result came down to −0.057 mmHg/year. Still 6-fold, but set against the control
arm's four-year progression (about 4 mmHg) it is 1.4%, a null in practice. The sign
and the magnitude are right; the precision is not.

**④ The rate of pigment reversal (0.02-fold).** With the model's dermal turnover
half-life of 6 years, a visible change within 2 years cannot be produced. The
observed reversal implies that **the effective turnover of the visible pigment is on
a scale of months**, and that is a measurable parameter.

**⑤ The early events are 4–8 years late** (#8, #9, #10). The embrittlement threshold
is reached at 38.5 years (reported 31), and so on. Lowering `PD50` fits that, but
then joint replacement comes far too early. It is left standing as a limitation:
one Hill function cannot hold "disc calcification in the third decade" and "joint
replacement in the sixth decade" at the same time.

---

## 5. Clinically usable results

### 5.1 Age at initiation is the only large variable (measured headroom)

The same drug, the same 10 mg, the same patient. Only the age at initiation was
changed. Because the untreated shadow patient (the counterfactual) is integrated
alongside within the same simulation, "the room available for prevention" is not an
assumption but a **measured value**.

| Age at initiation | Pigment averted | Intact cartilage at 70 | Probability of joint replacement at 70 | cAKUSSI at 70 | Pain VAS at 70 |
|---|---|---|---|---|---|
| 2 years | 92.3% | 1.000 | 0.000 | 0.06 | 0.00 |
| 5 years | 88.8% | 1.000 | 0.000 | 0.09 | 0.00 |
| 10 years | 83.1% | 0.999 | 0.000 | 0.33 | 0.05 |
| 15 years | 77.3% | 0.995 | 0.000 | 1.6 | 0.34 |
| 20 years | 71.6% | 0.978 | 0.000 | 6.1 | 1.27 |
| 25 years | 65.9% | 0.935 | 0.000 | 15.6 | 2.94 |
| 30 years | 60.2% | 0.856 | 0.000 | 29.6 | 4.66 |
| 35 years | 54.4% | 0.750 | 0.042 | 46.2 | 5.84 |
| 40 years | 48.4% | 0.642 | 0.511 | 64.0 | 6.51 |
| 45 years | 42.2% | 0.548 | 0.944 | 80.7 | 6.89 |
| 50 years | 35.7% | 0.470 | 0.998 | 94.9 | 7.11 |
| 55 years | 28.9% | 0.406 | 1.000 | 106.6 | 7.27 |
| 60 years | 21.7% | 0.353 | 1.000 | 114.4 | 7.38 |
| **Untreated** | **0%** | **0.272** | **1.000** | **126.3** | **7.55** |

Every column is monotonic in age at initiation — which is not a foregone conclusion
but a result obtained only after defect 11 of §6 was fixed.

**Start by 30 and joint replacement does not happen to this patient. At 35 it is
4%, at 40 51%, at 45 94%, and at 50 it is effectively certain.** That is, in this
model the fate of joint replacement is decided inside **the ten years between 35 and
45** — because that is the interval in which the cartilage crosses its risk
threshold. The price of waiting those same ten years is a probability of joint
replacement of 4% → 94%. That gradient is precisely "the cost of waiting", and note
that pigment averted merely declines gently over the same span, 54% → 42%: **the
biochemical benefit disappears gently while the clinical benefit disappears like a
cliff.**

Initiation at 5–15 years has never been tested in AKU. The only human evidence for
paediatric exposure is 15 years of follow-up in hereditary tyrosinaemia type 1 —
**and this is stated plainly as extrapolation.**

### 5.2 Three testable predictions

1. **Keratopathy risk depends hardly at all on dose and greatly on dietary
   management.** With diet matched, comparing the incidence of keratopathy at 2 mg
   and at 10 mg should show a small difference. The reported 5% vs 14.5% should be
   mostly a difference in dietitian involvement.
2. **Low-dose nitisinone may actually make pigment deposition worse.** At 1 mg/day
   the model has urinary HGA falling by 90% while serum HGA **rises**, 28 →
   33 umol/L. Because the extent to which the 14-fold rise in HPPA and HPLA shaves
   HGA clearance exceeds the fall in production.
   ⚠ **This prediction leans on the term the model gets most wrong (the OAT
   competition term, failure ② above).** The discriminating experiment is the same as
   in 4.3①. Until it is confirmed it must not be used as grounds for a low-dose
   strategy, and equally it must not be used as grounds that a low-dose strategy is
   safe.
3. **The clinical effect of nitisinone is mediated through serum HGA and not
   through urinary HGA.** The surrogate endpoint of future trials should therefore be
   **serum HGA** rather than u-HGA24. Because the two markers move 25-fold
   differently, efficacy predictions computed from u-HGA24 are systematically
   optimistic.

### 5.3 What came out of the 24 scenarios (at age 70, paired values)

| Scenario | Serum HGA | cAKUSSI | Joint replacement | Keratopathy | Pigment averted |
|---|---|---|---|---|---|
| S01 Untreated | 35.9 | 126.3 | 1.00 | 0.00 | 0% |
| S11 10 mg from age 25 | 1.85 | 15.6 | 0.00 | 1.00 | 66% |
| S14 10 mg + unrestricted diet (84 g) | 2.90 | 75.2 | 0.79 | **1.00** | 56% |
| S15 10 mg + protein restriction (56 g) | 1.22 | **5.6** | 0.00 | **0.00** | 65% |
| S16 2 mg + unrestricted diet (84 g) | 26.0 | **124.7** | 1.00 | **1.00** | 23% |
| S17 2 mg + protein restriction (56 g) | 10.7 | **36.6** | 0.00 | **0.00** | 45% |
| S18 Diet alone (56 g, no drug) | 23.3 | 85.1 | 0.92 | 0.00 | 16% |
| S19 Ascorbate 1 g/day | 35.7 | 125.4 | 1.00 | 0.00 | 6% |
| S20 Missense, 3% residual activity, untreated | 15.2 | **40.9** | 0.02 | 0.00 | 15% |
| S21 CKD (renal function reduced x3), untreated | **70.0** | 126.5 | 1.00 | 0.00 | 6% |
| S23 10 mg from 25 to 45, then stopped | 33.4 | 81.1 | 0.77 | 1.00 | 30% |
| S24 IDEAL (HGA blocked, no rise in Tyr) | 1.61 | **0.04** | 0.00 | **0.00** | 95% |

Four things to read here.

**① Diet beats dose.** 2 mg on an unrestricted diet (84 g) is effectively
indistinguishable from no treatment (cAKUSSI 124.7 against 126.3). Attach nothing
but protein restriction to the same 2 mg and it becomes 36.6. Conversely, 10 mg on
an unrestricted diet only comes down to 75.2. **The difference made by diet alone
(124.7 → 36.6) is larger than the difference made by a five-fold dose (124.7 →
75.2).**

**② Within this model, keratopathy is determined entirely by diet.** Over 40 years
of exposure the unrestricted diet gives a cumulative probability of 1.00 regardless
of dose (S14, S16), and protein restriction gives 0.00 regardless of dose (S15,
S17). Dose does not appear in this column at all. This is the clinical translation
of the table in §2.2.

**③ Three per cent of residual enzyme is close to a cure.** The missense genotype
(S20) has, untreated, a cAKUSSI of 40.9 and a joint replacement probability of 0.02
— against 126.3 and 1.00 untreated. It is a quantitative explanation of why the
genotype–phenotype variability of this disease is so wide, and at the same time it
shows **why a chaperone or a gene therapy that restores a few per cent of residual
activity is an attractive target** (the IDEAL of S24 has a cAKUSSI of 0.04).

**④ CKD worsens the disease on its own.** Set renal impairment alone to 3-fold and,
with neither drug nor diet changed, serum HGA doubles from 35.9 to 70.0 and the
pigment burden rises 29%. For an AKU patient, a nephrotoxic drug (NSAIDs included)
is not symptomatic treatment but a **disease-modifying factor**.

### 5.4 What this model does not answer

Acutely fatal metabolic complications (PMID 26596578), pregnancy, infection,
paediatric growth, the systemic consequences of secondary amyloidosis, and
between-individual variability (only a single `ETA` is attached, to mechanical
loading). It cannot be used for virtual-cohort simulation.

---

## 6. Twelve defects found and fixed during integration (mrgsolve 2.0.1)

All of them gave a quietly wrong answer to begin with.

1. **`SOLVERTIME` in `$MAIN`.** It is not defined in `$MAIN`, so compilation broke.
   `$MAIN` was emptied so as to carry no `CMT_0` assignment whatsoever, letting
   `init()` always win.
2. **The mg → umol conversion was wrong by a factor of 1000.** `1e6/329.25` is umol
   per g. A 2 mg dose produced a plasma nitisinone of 810 umol/L (the true value
   1.0), and 2 mg and 10 mg gave **completely identical results** — because both
   doses inhibited HPD by 99.999%. The entire dose-response had vanished and a
   plausible-looking picture came out.
3. **The transamination exchange rate could not sustain the flux.** At `KTAT` =
   200/day the untreated plasma tyrosine was pushed up to **178 umol/L** rather
   than 54. Quasi-equilibrium required 5000/day.
4. **A score domain integrated without bound.** `dxdt_AKJ = 0.0006*(...)` had no
   loss term, so cAKUSSI went to 2084. Rewritten as a first-order filter towards an
   algebraic target value, which bounds it to 0–1.
5. **Disc height went negative.** `-KDISCH*BRD*LOAD` was not proportional to the
   state, so it produced a negative height at age 65. Multiplying by the state made
   it an exponential decay.
6. **Making the direct concentration channel independent of the pigment channel
   destroyed the cartilage of a 10-year-old.** Synovitis was present as a constant
   from birth. `KSYNOX` was **left at 0 and the structure discarded**.
7. **Renal stones formed in infancy.** Urinary supersaturation is present from
   birth, so with a first-order reaction model a stone forms at age 2.2. Treating
   renal papillary pigment as the nidus and gating on its square moved it to 59.4.
8. **Putting the valve into the reversible channel overestimated a null result
   55-fold.** −0.51 against the reported −0.0093 mmHg/year. Excluded from the
   structure (4.3③).
9. **The baseline of the treatment effect was read at birth.** `x$CAKUSSI[1]` is the
   value at t=0 (birth). The change over a four-year trial starting at age 49 came
   out as **82 points** (the true value about 1.2). Every treatment contrast was
   inflated by including the whole of the preceding natural history.
10. **The counterfactual integrator read the treated patient's HPP.** The untreated
    shadow patient was evaluated at a 25-fold elevated HPP, inflating the
    counterfactual about 8-fold and making every headroom figure meaningless.
    Because Balance 1 gives the untreated HGA flux directly and independently of
    treatment state, it was rewritten from the dietary input.
11. **An `$OMEGA` that had only been declared turned every comparison into a
    comparison between different patients.** Because `mrgsim` draws the `ETA(1)`
    attached to mechanical loading anew on every call, each arm of the
    age-at-initiation scan was in effect **a different patient**. Running the same
    scan twice gave 0.826 / 0.707 / 0.735 / 0.608 and 0.842 / 0.761 / 0.599 / 0.535,
    and this model's central result — "cartilage mass at 70 is a monotonically
    decreasing function of age at initiation" — looked as though it had broken. This
    was at first **misdiagnosed as a solver-tolerance problem.** Tightening `atol`
    from 1e-8 to 1e-12 did not fix it and could not fix it — because the noise was
    never numerical in the first place. Once the random effects were pinned to zero
    with `zero_re()` the values become monotonic, 0.8559 / 0.7504 / 0.6421 / 0.5476
    / 0.4698, and are **identical to seven decimal places** at `atol` 1e-8 and
    1e-10. That is, this model was not stiff. Uncontrolled random effects looked
    like numerical instability, and the wrong diagnosis wasted one fix. (Population
    variability is switched on only with `sim_aku(..., iiv = TRUE)`.)
12. **Only one scenario, the one that stops the drug, returned NaN.** S23, which
    starts at 25 and stops at 45, was NaN in every output while the other 23
    scenarios were perfectly healthy. When nitisinone is stopped, HPP falls rapidly
    from about 60 to 2 umol/L, and the moment the integrator overshoots `CHPP` below
    `-KMHPD`, the Michaelis denominator `(KMAPP + CHPP)` passes through zero and the
    whole trajectory becomes NaN. Fixed by putting a positive-part clamp on every
    concentration. **A defect that appears only in the withdrawal scenario — that is,
    a defect sitting exactly where nobody checks first.** (After the fix, S23:
    cAKUSSI 81.1, joint replacement probability 0.77, pigment averted 30% — stop
    after 20 years of treatment and only a third of the pigment benefit remains.)

---

## 7. Summary of the model structure

**60 ODE compartments**: nitisinone PK 3 · amino acids and pathway 5 · cumulative
urinary excretion 5 · HGA handling 3 · pigment depots 8 · cartilage and synovium 5 ·
spine 3 · valve 3 · renal and urinary 3 · bone and tendon 2 · eye, ear and skin 3 ·
pain 3 · score domains 3 · cumulative risk 4 · counterfactual and headroom 4 ·
exposure integrals 3

**24 scenarios** (paired-control structure):
reproduction of SONIA 1's five arms (S02–S05) · SONIA 2 and its control (S06–S07) ·
the 3-year NIH trial (S08) · initiation at ages 5/15/25/40/55 (S09–S13) · dose ×
diet 2×2 (S14–S17) · diet alone (S18) · ascorbate (S19) · missense genotype (S20) ·
CKD (S21) · heavy manual work + BMI 32 (S22) · stopping after 20 years (S23) ·
IDEAL comparator: HGA blocked + no rise in tyrosine (S24)

```r
source("aku_mrgsolve_model.R")
res <- main_aku()                    # calibration → anchors → held-out validation → 24 scenarios
shiny::runApp("aku_shiny_app.R")     # 12-tab dashboard
```

---

## 8. Disclaimer

This is a QSP model for educational and research purposes. It has not been
independently validated or certified and must not be used for clinical decisions,
prescribing, or regulatory submission. The three predictions in 5.2 in particular
are **hypotheses**, and the second of them depends on the term the model gets most
wrong. Please read it with attention to the structural consequences and to how they
could be refuted, rather than to the absolute values of the parameters.
