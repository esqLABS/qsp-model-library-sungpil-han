# Opioid-Induced Constipation (OIC) — Quantitative Systems Pharmacology Model

51-ODE mrgsolve model · 4 colonic compartments × 6 antagonist/laxative classes · 16
treatment scenarios · 160-node 18-cluster mechanistic map · 110 references queried
directly on PubMed

---

## The single structural claim this model makes

> **The therapeutic window of a PAMORA is not potency. It is the ratio of two
> occupancies computed at the same plasma concentration, and what sets that ratio is
> P-glycoprotein at the blood–brain barrier.**

```
        SI  =  OCC_antagonist(enteric plexus)  /  OCC_antagonist(central)
```

The enteric plexus is **outside** the blood–brain barrier and the central compartment is
**inside** it. In both places **the same receptor** is competing with **the same
agonist**. The SI therefore reduces (in the low-occupancy limit) to **the reciprocal of
the antagonist's Kp,uu**, corrected for the agonist's own brain distribution. Nowhere in
this file is there a term that *assumes* "peripheral selectivity" — it is a division.

The second claim is less pretty, and it was obtained by **the author's hypothesis being
refuted** in the course of building this model: what amplifies the drug effect is not
the transit–water feedback loop but the **anorectal brakes** (claim 3 below).

---

## Computed results (computed, not assumed)

| # | Claim | What the model computed |
|---|---|---|
| 1 | Selectivity is the reciprocal of Kp,uu | SI: naloxone **12.3** · naloxegol **92.8** · methylnaltrexone **1211** · naldemedine **1841**. A **149-fold** range from Kp,uu alone — the binding Ki spans 25-fold and in the **opposite** order |
| 2 | P-gp inhibition rotates the ratio, not the exposure | CYP3A4 inhibition (pure exposure) is recovered by halving the dose (SI 17.5→31.3, COWS 1.79→0.33). P-gp inhibition is not recovered (SI 9.8→13.3) and only efficacy is lost (SBM 4.30→**3.11**, below the response criterion) |
| 3 | The amplifier is **not** a feedback loop | Remove the transit–water loop entirely and the drug-effect fold change goes 2.89 → **2.91** (no change). The real amplification is the rectal urge threshold (fold **1.60**) and anal sphincter tone (**1.88**) — the **anorectal** terms |
| 4 | The price of the tolerance asymmetry is not worse constipation | Over a 24-week titration the dose goes 60→200 mg/day and central availability 0.413→0.373 while the gut holds at 0.915. And yet the SBM goes 1.49→**1.34** (not worse). The price is that **pain does not budge from 5.04** |
| 5 | Lactulose is a prodrug for osmoles, PEG 3350 is not | PEG 17 g = 5.07 mmol = **17.5 mL** held. Lactulose 20 g = 58.4 mmol → fermentation ×3.6 → **725 mL**. The price is a bloating score of 1.28→2.95 |
| 6 | Methadone and lubiprostone collide at the same protein | ClC-2 activity 0.305→**0.103**, SBM 1.84→1.59. Linaclotide goes via GC-C and is unaffected (1.81 either way) |
| 7 | Cmax vs Cavg is **decided by the endpoint** | Split the dose and the weekly frequency actually **improves** (4.30→5.11) while the probability of laxation within 4 hours is **halved** (0.241→0.081). The latter is precisely the endpoint methylnaltrexone was approved on |

The source of every number is `oic_reference_output.txt` (see "Why a Python
implementation ships with this" below).

---

## Files

| File | Contents |
|------|------|
| [`oic_qsp_model.dot`](../../../opioid-induced-constipation/oic_qsp_model.dot) | Mechanistic map source — **160 nodes · 236 edges · 18 clusters** |
| [`oic_qsp_model.svg`](../../../opioid-induced-constipation/oic_qsp_model.svg) | Zoomable vector map (10579 × 3565) |
| [`oic_qsp_model.png`](../../../opioid-induced-constipation/oic_qsp_model.png) | 150 dpi raster |
| [`oic_mrgsolve_model.R`](../../../opioid-induced-constipation/oic_mrgsolve_model.R) | 51-ODE mrgsolve model + 16 scenarios + 5 analysis functions |
| [`oic_shiny_app.R`](../../../opioid-induced-constipation/oic_shiny_app.R) | 10-tab interactive dashboard |
| [`oic_reference_model.py`](../../../opioid-induced-constipation/oic_reference_model.py) | Independent Python/scipy implementation — **the only implementation that was actually integrated** |
| [`oic_reference_output.txt`](../../../opioid-induced-constipation/oic_reference_output.txt) | Run log of that implementation — the source of **every** number in this README |
| [`oic_population_results.json`](../../../opioid-induced-constipation/oic_population_results.json) | Virtual population results (machine-readable) |
| [`oic_references.md`](../../../opioid-induced-constipation/oic_references.md) | **110 references** queried directly on PubMed, 14 sections |

---

## Why a Python implementation ships with this (Provenance)

The container in which this repository was built **had no R runtime.** Rather than commit
an unverified ODE model, the judgement was that it is more honest to implement every
equation independently in Python, actually integrate it, and commit the result.
`oic_mrgsolve_model.R` was **transcribed term by term** from `oic_reference_model.py`,
and not the other way round.

The consequence a reader needs to know: **the science was verified by integration, and
the R syntax was verified only by eye.** Whoever runs the R file first should be
prepared to fix typographical errors, and should not take a compilation failure as
evidence about the model.

### Seven defects exposed by actually integrating it

Each is marked with a comment in both files. Someone reading a QSP model has a right to
know which lines were hard.

| # | Defect | How it showed up |
|---|---|---|
| 1 | Water absorption was first-order in luminal water, so a compartment could dry out completely | The model diverged in the **healthy** arm — 1064 g of solids, SBM 0.02/week. The mucosa cannot extract water bound to the solid phase, so absorption has to stop at w_min = 0.615. This term is what makes the transit–water loop an *amplifier* rather than a *divergence* |
| 2 | With no rescue laxative the colon loaded without limit | The OIC arm converged on more than 2 kg of solids. And it was **required by the definition of the endpoint** — an SBM is "a bowel movement with no rescue medication in the preceding 24 hours", so an evacuation driven by rescue medication must empty the colon without being counted |
| 3 | The gut MOR was set up to act on **motility** only | The defining result for methylnaltrexone (laxation within 4 hours) could not be reproduced **at any potency** — even with complete blockade it ceilings at P=0.19, because it is impossible if transit time is the only route by which faecal water can change within 4 hours. MOR is also on the submucosal secretomotor neurons and inhibits secretion |
| 4 | The anal sphincter brake alone could not bring the chronic arm below 2.3 SBM/week | The mechanism actually missing was **opioid-induced rectal hyposensitivity** — the distension threshold at which the urge is felt itself rises. This is exactly the "no urge to go" that patients report, and it is also what a PAMORA reverses |
| 5 | The luminal opioid term overwhelmed the plasma 80-fold | 1500 nM delivered to the plexus against 19 nM free in plasma. **Opioid dose dependence was lost completely** (30 mg and 60 mg giving 2.60 vs 2.58 SBM/week) |
| 6 | The rescue-medication trigger rule was piecewise linear, which put a **discontinuity** into the steady state | Between AANO 3.0 and 3.2 the SBM jumped from 2.44 to 1.57. Smoothed, so that the steepness that remains is biology and not a rule |
| 7 | The lumen→plexus access coefficient was shared by every drug | **Oral** methylnaltrexone 450 mg became a miracle drug (SBM 11.0/week, complete normalisation). A permanently charged quaternary ammonium cannot cross the mucosa like that, and the fact that it cannot is precisely **why** this drug is peripherally restricted |

---

## What was fitted (Calibration provenance)

**Three parameters were fitted in the whole of this model.** Everything else is either a
literature value or a stated structural assumption.

| Parameter | Drug | Fitted to | Value |
|---|---|---|---|
| `KIGUT` | Naloxegol | KODIAC-04/05 25 mg arm, mean SBM/week 4.3 | 2.685 nM |
| `KIGUT` | Naldemedine | COMPOSE-1/2 0.2 mg arm, mean SBM/week 5.0 | 0.00821 nM |
| `KIGUT` | Methylnaltrexone | Thomas 2008, placebo-adjusted 4-hour laxation difference 33 pp | 0.795 nM |

### And what that fit revealed the model cannot explain

All three values had to be **3–41-fold more potent than the published MOR binding Ki**
of the drug concerned:

| Drug | Literature Ki (nM) | Fitted KIGUT (nM) | Ratio |
|---|---|---|---|
| Naloxegol | 7.400 | 2.685 | 0.363 |
| Naldemedine | 0.340 | 0.00821 | 0.024 |
| Methylnaltrexone | 28.000 | 0.795 | 0.028 |

**The model does not explain this gap, and it states that it does not.** Two candidates
were actually ruled out:

- If it were **a common error in the agonist-side anchor (`KIOP`)** the three ratios
  would have to move by **the same factor**. They differ 15-fold, so it is not that.
- **Plasma protein binding** cannot get the ordering right either. Naloxegol (fu 0.96)
  requires the smallest correction, and yet methylnaltrexone (fu 0.885) requires as
  large a correction as naldemedine (fu 0.065).

The candidates that remain (not ruled out): accumulation in the gut wall beyond the
pre-systemic term modelled here, active uptake into enteric neurons, error in the CL/F
values used to set the plasma exposure, and the fact that the three trials enrolled
different populations against endpoints over different windows.

> **The per-drug SBM predictions inherit this uncertainty in full. The selectivity
> results (claims 1 and 2) do not — because they are a ratio in which `KIGUT` cancels
> out.**

---

## The main results in detail

### 1 · The selectivity index — one division

The same OIC patient, each drug at its label dose:

| Drug | Kp,uu | 1/Kp,uu | Gut occupancy | Central occupancy | **SI** | COWS | Δpain | SBM/week |
|---|---|---|---|---|---|---|---|---|
| Naloxone PO | 1.000 | 1.0 | 0.977 | 0.0792 | **12.3** | 8.34 | −0.25 | 9.75 |
| Naloxegol | 0.020 | 50.0 | 0.569 | 0.0061 | **92.8** | 0.01 | −0.19 | 4.30 |
| Methylnaltrexone SC | 0.005 | 200.0 | 0.840 | 0.0007 | **1211** | 0.00 | −0.21 | 7.16 |
| Naldemedine | 0.012 | 83.3 | 0.715 | 0.0004 | **1841** | 0.00 | −0.18 | 5.00 |

Naloxone has the highest gut occupancy of all (0.977) and yet cannot be used — because
8 % central occupancy comes along with it and a COWS of 8.34 comes out. That is the
whole of "potency is not the therapeutic window".

One **as yet untested prediction** derives from this: since the SI also depends on the
brain concentration of the agonist (the more agonist there is in the brain, the harder
the antagonist's central competition), **the central safety margin of a PAMORA should be
narrower in a patient on morphine (Kp,uu ≈ 0.3)** than in one on oxycodone
(Kp,uu ≈ 3), whose brain uptake is active. It has never been confirmed clinically and is
offered only as a prediction.

### 2 · P-gp inhibitors — when the dose reduction aims at the wrong variable

Naloxegol 25 mg, three perturbations × the label's dose-reduction instruction:

| Situation | Free plasma (nM) | Gut | Central | SI | COWS | SBM/week |
|---|---|---|---|---|---|---|
| No interaction | 6.9 | 0.569 | 0.0061 | 92.8 | 0.01 | 4.30 |
| CYP3A4 inhibition only (AUC ×3.4) | 67.9 | 0.908 | 0.0518 | 17.5 | 1.79 | 7.10 |
| CYP3A4 inhibition + reduction to 12.5 mg | 34.0 | 0.831 | 0.0266 | **31.3** | 0.33 | 5.26 |
| P-gp inhibition only (Kp,uu ×10) | 6.9 | 0.569 | 0.0581 | 9.8 | 2.61 | 4.30 |
| P-gp inhibition + reduction to 12.5 mg | 3.5 | 0.397 | 0.0299 | **13.3** | 0.50 | **3.11** |
| Both (verapamil/diltiazem type) | 67.9 | 0.908 | 0.3531 | 2.6 | **30.38** | 7.10 |
| Both + reduction to 12.5 mg | 34.0 | 0.831 | 0.2144 | 3.9 | **23.50** | 5.26 |

**Look at the P-gp inhibition row. The plasma concentration has not changed (6.9 nM, the
same) and neither has the gut occupancy (0.569, the same). The only thing that moved is
the brain.** So the dose reduction recovers almost none of the safety margin (SI
9.8→13.3, about a seventh of baseline) while taking away exactly the efficacy the
patient needs (SBM 4.30→3.11, below the response criterion of 3).

What is common in the clinic is the last two rows — verapamil, diltiazem and quinidine
inhibit both.

### 3 · Where the amplifier is (refutation of the author's hypothesis)

The naloxegol 0 → 25 mg comparison, repeated with one brake switched off at a time:

| Brake switched off | Untreated | Naloxegol 25 | **Fold** |
|---|---|---|---|
| (none — the full model) | 1.49 | 4.30 | **2.89** |
| Transit–water feedback loop (GW=GWD=0) | 1.32 | 3.85 | **2.91** |
| Shift in the rectal urge threshold (AVSENS=0) | 4.19 | 6.69 | **1.60** |
| Anal sphincter tone (AANO=0) | 2.94 | 5.54 | **1.88** |
| cAMP→ACh→HAPC propulsion (EMAXMOR=0) | 2.23 | 4.83 | **2.16** |
| Segmental tone (ATONE=0) | 1.49 | 4.44 | 2.98 |
| Inhibitory NO/VIP tone (ANO=0) | 1.57 | 4.43 | 2.81 |
| Secretory inhibition (ESECMOR=0) | 1.63 | 4.51 | 2.77 |

**This analysis refuted the hypothesis the author set out to confirm.** Remove the
transit–water loop in its entirety and the fold change is essentially unchanged, 2.89 →
2.91. That loop sets the absolute *level* of faecal water (untreated 1.49 → 1.32) but it
does not **carry** the drug effect.

What carries it are the **anorectal** terms — the urge threshold and sphincter tone.
Colonic propulsion comes third, and the segmental tone, inhibitory tone and secretory
terms contribute essentially nothing to the fold change.

So the amplification is not a feedback gain but **the product of several independent
brakes that a single receptor releases at once**, and that product is dominated at the
**end** of the colon. Two testable consequences:

1. A PAMORA should work far better in a patient whose problem is **loss of the urge and
   sphincter tone** than in one with **genuine colonic inertia**. The prediction is that
   responders separate along an **anorectal phenotype** rather than a transit-time
   phenotype, and it can be checked with anorectal manometry.
2. A drug that removes the cause acts on every brake at once and therefore
   **multiplies**, while a drug that compensates one downstream branch (a secretagogue,
   an osmotic agent) **adds**. This is the structural reason a PAMORA beats the laxatives
   by more than their respective differences in potency.

### 4 · The tolerance asymmetry — the price is analgesia, not constipation

After 12 weeks of equilibration, titrated weekly over 24 weeks to hold pain at NRS 4.0:

| Week | Oxycodone mg/day | Central occupancy | **Central availability** | Pain | Gut occupancy | **Gut availability** | SBM/week |
|---|---|---|---|---|---|---|---|
| 0 | 60 | 0.722 | 0.413 | 5.04 | 0.625 | 0.926 | 1.49 |
| 4 | 133 | 0.852 | 0.383 | 5.00 | 0.787 | 0.919 | 1.36 |
| 8 | 200 | 0.897 | 0.373 | 5.03 | 0.848 | 0.915 | 1.34 |
| 24 | 200 | 0.897 | 0.373 | 5.04 | 0.848 | 0.915 | 1.34 |

The asymmetry is real and it is **7-fold** (central loss 0.627 against gut loss 0.085).
And yet **the bowel endpoint does not get worse** — gut occupancy climbs from 0.625 to
0.848 while the SBM is essentially unchanged, 1.49 → 1.34. That is because gut
transduction is already saturated at ordinary analgesic doses.

So this model **does not support** the story that "dose escalation progressively worsens
constipation". It says something more specific and more testable instead: **the
constipation a patient has at 60 mg/day is essentially the constipation they will have
at 200 mg/day.** This is the model's explanation for an epidemiological fact that a
simple dose–response account cannot handle — that OIC prevalence correlates only weakly
with opioid dose.

The asymmetry that matters clinically is on the **other side** of the ledger. The
titration never reaches its target: pain stays at 5.04 throughout a 3.3-fold escalation
and improvement stops beyond ~200 mg/day. **The escalation buys no analgesia, and the
bowel receives all of it.** This is the situation a PAMORA exists to break.

### 5 · The arithmetic of osmoles

| Agent | mmol | Water held | After fermentation |
|---|---|---|---|
| PEG 3350 17 g | 5.07 | 17.5 mL | 17.5 mL |
| Lactulose 20 g | 58.43 | 201.5 mL | **725.3 mL** |
| Lactulose 40 g | 116.86 | 403.0 mL | **1450.6 mL** |

PEG is a 3350 Da polymer given in **grams**, so 17 g is a mere 5.1 mmol. Lactulose is
342 Da and **on top of that** the colonic bacteria cut it into osmotically active
fragments. In the same simulation lactulose 20 g takes the SBM from 1.49 to 6.73, but
the bloating score from 1.28 to 2.95 — the fermentation that makes the osmoles makes the
gas as well.

### 6 · Methadone and lubiprostone

| Arm | SBM/week | Bristol | ClC-2 activity |
|---|---|---|---|
| Oxycodone, untreated | 1.49 | 2.01 | 0.000 |
| Oxycodone + lubiprostone | 1.84 | 2.60 | **0.305** |
| Oxycodone + linaclotide | 1.81 | 2.52 | 0.000 |
| Methadone, untreated | 1.49 | 2.01 | 0.000 |
| Methadone + lubiprostone | **1.59** | 2.18 | **0.103** |
| Methadone + linaclotide | **1.81** | 2.52 | 0.000 |
| Methadone + naldemedine | **5.00** | 2.85 | 0.000 |

Linaclotide reaches the same CFTR-mediated secretion by way of GC-C/cGMP and so bypasses
this block. A PAMORA is doubly unaffected — because it removes the cause instead of
compensating downstream.

### 7 · Acute vs chronic — the endpoint is what makes the mechanism visible

The same total daily dose, differing only in how it is split:

| Regimen | SBM/week (chronic) | P(laxation within 4 h) |
|---|---|---|
| (no antagonist) | 1.49 | 0.032 |
| Naloxegol 25 mg qd | 4.30 | **0.241** |
| Naloxegol 12.5 mg bid | 4.59 | 0.141 |
| Naloxegol 6.25 mg qid | **5.11** | **0.081** |
| Methylnaltrexone 12 mg SC qd | 7.16 | **0.362** |
| Methylnaltrexone 6 mg SC bid | 7.83 | 0.265 |
| Methylnaltrexone 3 mg SC qid | **7.91** | **0.161** |

Weekly frequency is a mean, so the Cmax–Cavg distinction is **invisible** in it —
splitting the dose actually raises the mean occupancy and improves it. The distinction
lives in the acute column, and that column is precisely the endpoint methylnaltrexone
was approved on (laxation within 4 hours without rescue medication) and the endpoint
naloxegol and naldemedine were **not** approved on.

---

## Known limitations, and the discrepancies left in rather than fixed

This section is where the author **thinks the model is wrong**. It has not been erased
by tuning.

1. **PEG 3350 is essentially inert in the model** (SBM 1.49 → 1.49, and 1.50 even at
   34 g/day). This is a consequence of the arithmetic: 17 g holds 17.5 mL
   colligatively, and hydrogen-bonded hydration (~2.5 waters per ethylene oxide unit →
   1.02 g water per g of polymer) captures only a comparable amount. **The author
   regards this prediction as wrong** — PEG 3350 plainly does work clinically. Its
   working must therefore lie not in the physics of the polymer itself but in the fact
   that it is given in 240 mL of vehicle and in the fact that it is **not fermented**
   (unlike lactulose it is not converted into absorbable SCFAs and gas). The model does
   not contain that route.
2. **Do not use this to predict responder rates.** In the default virtual population
   (analysis J) the untreated arm has **0 %** responders against 91.7 % for naloxegol
   and 95.0 % for naldemedine, which does not match the trials' placebo-adjusted
   differences (+12.7 and +16.0 pp) at all. Rather than leave that as a bare failure,
   **how much dispersion it takes to reproduce the placebo arm** was measured (analysis
   J2): inflate every parameter CV **3.2-fold** and the untreated arm becomes
   **28.3 %**, matching the trials' placebo 29 %. Reading the active arms again in that
   state:

   | Arm | Model | Trial | Model Δ | Trial Δ |
   |---|---|---|---|---|
   | Untreated / placebo | 28.3% | 29.4% | — | — |
   | Naloxegol 25 mg | 66.7% | 44.4% | **+38.3 pp** | +12.7 pp |
   | Naldemedine 0.2 mg | 78.3% | 47.6% | **+50.0 pp** | +16.0 pp |

   That is, even with the dispersion matched it **overestimates the size of the drug
   effect about 3-fold**. And the 3.2-fold inflation is not itself a mechanistic finding
   but **a measurement of a deficiency** — an opioid dose distribution with a log-scale
   CV of 1.44 does not exist. What that number says is that the mechanistic parameters
   varied here carry only **about a third** of the between-patient variance the trials
   actually contain, and that the rest lies in things this model has no state variable
   for (adherence, diet, activity, coexisting anorectal disease, co-prescribed
   constipating drugs, and the placebo response itself).
3. **There is a structural ceiling on methylnaltrexone's acute endpoint.** Even with the
   gut blocked **completely**, the probability of laxation within 4 hours stops at 0.481
   and does not reach the trials' 62 % on the first dose (as a placebo-adjusted
   difference it is close, 45 pp against 48 pp). The potency of this drug is therefore
   **unidentifiable** from the first-dose endpoint, and it was fitted to the ≥2/4-doses
   endpoint instead.
4. **The premise that ClC-2 is lubiprostone's principal route is under dispute.** §8 of
   the references (Oak 2022) reports that the contribution of ClC-2 to lubiprostone's
   secretory effect is incidental. This model puts ClC-2 as the principal route, so
   **this model's account of the methadone–lubiprostone interaction collapses along with
   that premise if the premise is wrong.**
5. **There is no placebo mechanism.** The untreated arm is mechanistically untreated and
   is not comparable with the trials' placebo arms (KODIAC 29 %, COMPOSE 34 %).
6. **The R file has never been run.** See the "Provenance" section above.
7. The effect sizes for prucalopride, linaclotide and lubiprostone were not fitted and
   are order-of-magnitude estimates. Their SBM predictions are predictions, not
   reproductions.

---

## The 16 scenarios (all of them)

| Scenario | SBM/week | CSBM | Bristol | Straining | PAC-SYM | Pain | COWS | Rescue/week |
|---|---|---|---|---|---|---|---|---|
| S01 Healthy (no opioid) | 11.16 | 6.02 | 3.99 | 0.94 | 0.60 | 7.26 | 0.00 | 0.04 |
| S02 OIC untreated (oxycodone 60/day) | 1.49 | 0.27 | 2.01 | 2.22 | 2.57 | 5.04 | 0.00 | 1.52 |
| S03 + PEG 3350 17 g/day | 1.49 | 0.27 | 2.01 | 2.22 | 2.57 | 5.04 | 0.00 | 1.52 |
| S04 + lactulose 20 g/day | 6.73 | 2.80 | 6.05 | 0.26 | 1.05 | 5.42 | 0.00 | 0.14 |
| S05 + lubiprostone 24 µg bid | 1.84 | 0.38 | 2.60 | 1.71 | 2.23 | 5.05 | 0.00 | 1.30 |
| S06 + linaclotide 145 µg qd | 1.81 | 0.37 | 2.52 | 1.77 | 2.27 | 5.05 | 0.00 | 1.32 |
| S07 + prucalopride 2 mg qd | 1.92 | 0.66 | 2.18 | 2.06 | 2.30 | 4.95 | 0.00 | 1.25 |
| S08 + naloxegol 12.5 mg qd | 3.11 | 1.05 | 2.37 | 1.89 | 2.16 | 4.90 | 0.00 | 0.94 |
| S09 + naloxegol 25 mg qd | 4.30 | 1.75 | 2.67 | 1.66 | 1.92 | 4.85 | 0.01 | 0.83 |
| S10 + naldemedine 0.2 mg qd | 5.00 | 2.15 | 2.85 | 1.53 | 1.68 | 4.86 | 0.00 | 0.70 |
| S11 + methylnaltrexone 12 mg SC qod | 4.22 | 1.81 | 2.63 | 1.68 | 2.18 | 4.87 | 0.00 | 1.13 |
| S12 + methylnaltrexone 450 mg PO qd | 6.23 | 2.90 | 3.13 | 1.36 | 1.46 | 4.84 | 0.00 | 0.56 |
| S13 + naloxone 20 mg PO tid | 9.75 | 5.09 | 3.81 | 1.02 | 0.69 | 4.79 | **8.34** | 0.10 |
| S14 + naloxegol 25 + PEG 17 g | 4.31 | 1.75 | 2.67 | 1.65 | 1.92 | 4.85 | 0.01 | 0.83 |
| S15 Naloxegol 25 + strong P-gp inhibition | 7.10 | 3.40 | 3.34 | 1.25 | 1.07 | 4.82 | **30.38** | 0.33 |
| S16 Methadone 60/day + lubiprostone | 1.59 | 0.30 | 2.18 | 2.06 | 2.48 | 5.04 | 0.00 | 1.46 |

The COWS column of S13 (naloxone) and S15 (P-gp inhibition) is the point of this model:
in both cases **the bowel effect is excellent** and the patient falls into withdrawal.

---

## Model structure (51 ODEs)

| Block | ODEs | Contents |
|---|---|---|
| Opioid PK | 5 | depot · central · peripheral · free brain concentration · lumen |
| Antagonist PK | 5 | the same structure + the brain separated by Kp,uu |
| Concomitant drugs | 6 | PEG · lactulose · lubiprostone · linaclotide · prucalopride ×2 |
| Receptor trafficking | 3 | gut availability · central availability · β-arrestin |
| Enteric signal transduction | 5 | cAMP · ACh · NO/VIP · cGMP · ClC-2 |
| Motility | 2 | HAPC · segmental tone |
| Colonic contents | 12 | 4 compartments × (solids · water · osmoles) |
| Symptoms | 5 | Bristol · straining · bloating · PAC-SYM · PAC-QOL |
| Central | 4 | pain · counter-adaptation · withdrawal · nausea |
| Counters · safety | 4 | SBM · CSBM · rescue medication · impaction risk |

---

## Usage

```r
# R
source("oic_mrgsolve_model.R")
res <- run_all_scenarios()          # the 16 scenarios
selectivity_table()                 # selectivity index
ddi_table()                         # P-gp vs CYP3A4
brake_decomposition()               # which brake carries the drug effect
shiny::runApp("oic_shiny_app.R")    # 10-tab dashboard
```

```bash
# Python (the implementation that was actually verified)
python3 oic_reference_model.py      # regenerates the whole A–J analysis
```

---

## ⚠️ Disclaimer

This is a **qualitative to semi-quantitative QSP model for educational and research
purposes**. It was assembled from the published literature and clinical-trial data but
has not been independently validated or certified, and **must not be used directly for
real clinical decision-making, prescribing, or regulatory submission.** In particular,
do not quote a number from this model without having read the seven items in the "Known
limitations" section above.
