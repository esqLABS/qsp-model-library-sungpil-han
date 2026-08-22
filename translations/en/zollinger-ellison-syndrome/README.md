# Zollinger–Ellison Syndrome (Gastrinoma) — QSP Model
### Zollinger–Ellison Syndrome / Gastrinoma · Quantitative Systems Pharmacology

| Deliverable | File |
|--------|------|
| 🗺️ Mechanistic map | [`zes_qsp_model.dot`](../../../zollinger-ellison-syndrome/zes_qsp_model.dot) · [SVG](../../../zollinger-ellison-syndrome/zes_qsp_model.svg) · [PNG](../../../zollinger-ellison-syndrome/zes_qsp_model.png) — 188 nodes / 275 edges / 22 clusters |
| ⚙️ mrgsolve model | [`zes_mrgsolve_model.R`](../../../zollinger-ellison-syndrome/zes_mrgsolve_model.R) — 59 ODEs, time unit = hour, 37 scenarios |
| 📊 Shiny dashboard | [`zes_shiny_app.R`](../../../zollinger-ellison-syndrome/zes_shiny_app.R) — 12 tabs |
| 📚 References | [`zes_references.md`](../../../zollinger-ellison-syndrome/zes_references.md) — 105 PubMed entries (every PMID verified) |

---

## The organising thesis

**Gastric acid output is not a state variable. It is a product.**

```
BAO  =  BAOCAP  ×   PCM     ×   PUMPA    ×   ACTP
                  ─────────    ─────────     ─────────
                  FACTOR 1     FACTOR 2      FACTOR 3
                  Parietal     Fraction of   Activation per pump
                  cell mass    pumps in      by secretagogue
                  (weeks)      the membrane  stimulus
                               (hours~days)  (minutes~hours)
```

Hypergastrinaemia raises **factors 1 and 3**. Yet every acid-suppressing agent ever
licensed acts on **factor 2 only**, and on nothing else. This asymmetry — not any
special "resistance" of the gastrinoma patient's parietal cell — is the single
structure this file sets out to express, and as a result all of the following become
**outputs of the model rather than coded-in rules**.

- **Why the maintenance dose in ZES is several times the reflux-disease dose.**
  Remove the same *fraction* of the pump pool, and if factors 1 and 3 are both
  raised then the *absolute* acid output that remains is correspondingly larger.
- **Why the dosing interval fails before the daily total does.** Pump synthesis
  (`KSYNP`) is upregulated by the trophic gastrin signal, so the pump pool of a
  hypergastrinaemic stomach **refills faster.** Nocturnal acid breakthrough is not a
  separate phenomenon but a consequence of pump-pool dynamics plus the histamine
  circadian rhythm.
- **Why the H2 blocker is outcompeted.** `DRIVE` is the **product** of the three
  stimulatory pathways, so blocking one leaves the other two still multiplying. A
  covalent pump blocker lies downstream of all three pathways and cannot be pushed
  aside by competition.
- **Why somatostatin analogues, debulking surgery and PRRT are additive to a PPI
  rather than redundant with it.** They lower `SPHEN` (gastrin secretion per unit
  tumour volume) and so touch factors 1 and 3. The PPI cannot touch either of them
  at all.
- **Rebound hypersecretion on withdrawal.** An enlarged `PUMPI` pool pouring onto an
  enlarged `PCM`.

### Two sign inversions — drawn exactly as sign inversions

1. **Secretin.** `SECN = 1/(1 + CSEC/KSECN)` multiplies the set-point of the
   **normal antral G cell** (inhibition). `SECT = 1 + ESECT·CSEC/(KSECT+CSEC)`
   multiplies the secretion rate of the **tumour** (stimulation). Same compartment,
   opposite sign, two targets. So the secretin stimulation test is an output of the
   model, and so is the way that test breaks down.
2. **Intragastric acid.** `BRK = 1/(1+(HC/HBRK50)^HBRKN)` regulates only the antral
   G cell and regulates **nothing at all** in the tumour. It is not that the
   secretion rate per cell is high — it is this **one missing edge** that makes the
   hypersecretion unstoppable. And because of the **same missing edge**, giving a PPI
   to somebody with no gastrinoma at all raises gastrin (BRK approaches 1) and
   manufactures the commonest false-positive diagnosis in this disease.

### The MEN1 arm is not a label but a feedback loop

`GLANDM` (parathyroid functional mass) → `PTH` → `CAION` → `CAF` (the CaSR-mediated
facilitation factor applied to gastrin secretion in both the tumour and the antrum),
and → the calcium arm of `DRIVE`. Nowhere in the file is there a sentence saying
"parathyroidectomy lowers gastrin". Scenario 17 returns `GLANDM` to 1.0, scenario 17b
runs the same patient without surgery alongside it, and the model states the
difference.

---

## The distinction between inputs and outputs

**Inputs** (what the user sets): tumour volume `TUM0`/`MET0`, grade `GRADEF`, SSTR2
density `SSTR2D`, MEN1 status (`GLANDM`, `MENFLG`), CYP2C19 phenotype `CYPF`, renal
function `RENF`, NSAID and *H. pylori* exposure, drug dose and regimen.

**Outputs** (what the model has to produce, and is never told): basal and maximal
acid output and the BAO/MAO ratio, fasting serum gastrin, intragastric pH and the
24-hour pH>4 holding time, the secretin increment and its false positives, the
maintenance dose required, nocturnal acid breakthrough, ulcer and oesophageal
injury, diarrhoea, rebound after withdrawal, the development of ECL hyperplasia and
gastric neuroendocrine tumours, tumour progression, renal absorbed dose, and the
three headroom fractions.

---

## The one deliberately phenomenological relationship — declared, not hidden

Intragastric pH is **not** computed as the `-log10[H+]` of the luminal H⁺ pool. Do
that and the achlorhydric limit cannot be reached. With first-order emptying alone
the residual H⁺ pool cannot fall by the four orders of magnitude that separate pH 1.7
from pH 6.3, so a completely blocked stomach comes out at pH 3.3 (which is exactly
what happened in the first version). Gastric contents are a buffer system, so pH
follows a **titration curve**, and the model fits a single-line Hill titration curve

```
pH = PHMIN + (PHMAX − PHMIN) / (1 + (HC/HC50)^PHN)
```

to the two directly measured endpoints — normal fasting pH about 1.7 and complete
achlorhydria about pH 6.3. This is the only baldly empirical relationship in the
file, and rather than bury it I flag it here.

---

## Validation

`zes_validate()` runs the A0–A13 suite. None of the values below is coded into the
equations. All of them are results actually run in mrgsolve 2.0.1.

### Fitted anchors (9) · coordinate-wise fixed-point iteration

Each anchor is paired with the one parameter that has dominant leverage over it. Of
26 iterations, the one with the minimum objective (iteration 23) was adopted, after
which two structural parameters (`KMPA` 0.10 → 0.035, `HC50` 5.0 → 7.0) were changed
and a further 18 iterations run. The second run produced no improvement (iteration 1
was the best), but the `KMPA` change raised the A7 ceiling from 49.7% to 51.6% and
moved A8 from 16.7% to 21.9%, closer to target, so that configuration was kept.
**Objective of the shipping configuration = 0.196** (normalised sum of squares).
The "Model" column below is entirely measured values from the shipping configuration.

| # | Anchor | Target | Model | Fitted parameter |
|---|------|------|------|---------------|
| A1 | Normal fasting basal acid output (mEq/h) | 3.0 | **2.99** | `BAOCAP` = 78.9 |
| A2 | Normal MAO/BAO ratio (pentagastrin) | 7.67 | **7.50** | `KDRV` = 45.4 |
| A3 | Untreated sporadic ZES fasting acid output (mEq/h) | 36 | **33.7** | `EMAXPCM` = 2.07 |
| A4 | Untreated sporadic ZES fasting gastrin (pg/mL) | 900 | **900** | `KSECT0` = 562 |
| A5 | Normal fasting gastrin (pg/mL) | 30 | **30.9** | `KSECA` = 901 |
| A6 | Normal parietal cell mass (no drift) | 1.00 | **0.98** | `GBNORM` = 27.1 |
| A7 | Normal, omeprazole 20 mg qd for 6 days: fall in acid output (%) | 66 | **51.6** ✗ | `KBIND` = 40 |
| A8 | Normal untreated 24-hour pH>4 holding time (%) | 20 | **21.9** | `BUFMEAL` = 107 |
| A9 | Normal, omeprazole 20 mg qd: pH>4 holding time (%) | 45 | **61.5** ✗ | `KEMPT` = 2.62 |

Six fit well and **three (A7, A9, and one of the held-out predictions below) do
not.** Rather than absorb the cause into a parameter I have left it standing and
written it up below.

`BUFMEAL = 107 mEq/meal` is larger than the measured titratable buffering capacity
of a meal (20~40 mEq). This is not a good sign but **a sign that the pH map is
absorbing something the buffer model does not capture**, and I record it as a
limitation in its own right.

### Methodological decisions — why Nelder-Mead was abandoned, why GSYN was not fitted

- A six-parameter Nelder-Mead did not converge in 40 minutes. Since each anchor has
  one dominant parameter, a coordinate-wise damped multiplicative update completes an
  iteration in six simulations. The coordinates are **not independent** (`BAOCAP`,
  `EMAXPCM` and `KSECA` all move acid output, acid output moves the D-cell brake, and
  that moves gastrin and comes back round), so every damping exponent is set well
  below 1 and **the iteration with the minimum objective is returned rather than the
  last one**.
- A tenth anchor (acid output ≈ 6 mEq/h in ZES on omeprazole 30 mg bd) was to be
  paired with `GSYN`, but **the calibration collapsed.** `GSYN` has more leverage over
  *untreated* ZES acid output (A3) than over *treated* ZES, so `GSYN` and `EMAXPCM`
  fought over A3 and the iteration oscillated and then diverged (untreated ZES acid
  output 144 mEq/h, `EMAXPCM` down to 0.07 before hitting the bound). `GSYN` is left
  at its literature value of 0.70, and the over-suppression is **reported below as a
  failure**.

### Held-out predictions — those that passed

These are values that were never fitted. All are actual `zes_validate()` output.

**(1) The decomposition into three factors is actually measured (A1).** A
counterfactual integrator runs alongside the real trajectory and measures each
factor's share. That `HF_PUMP = 0.000` in the untreated state means that with no
drug present, counterfactual C converges exactly onto the real trajectory (in an
earlier version it came out at −0.67).

| Scenario | BAO | HF_PCM | HF_DRIVE | HF_PUMP | PCM | ECL |
|----------|-----|--------|----------|---------|-----|-----|
| 02 untreated sporadic ZES | 38.8 | 0.421 | 0.645 | **0.000** | 1.73 | 1.95 |
| 06 omeprazole 30 mg bd | 1.4 | 0.422 | 0.645 | **0.878** | 1.73 | 1.97 |
| 20 omeprazole bd + octreotide | 1.2 | 0.415 | **0.549** | 0.915 | 1.72 | **1.06** |

The PPI cannot move `HF_DRIVE` at all (0.645 → 0.645). Add a somatostatin analogue
and it falls to 0.549 while the ECL mass normalises from 1.97 to 1.06. **This is the
central claim of the file, and it is measured rather than coded.**

**(2) At an identical 60 mg/day, the dosing interval beats the total (A2).** The
total mg is the same in all three arms.

| Arm | Fasting BAO | Nocturnal BAO | pH>4 (%) | Nocturnal pH<4 (h) |
|----|----------|----------|----------|---------------|
| 60 mg once daily | **6.8** | 5.6 | 50 | 8.0 |
| 30 mg twice daily | **2.0** | 1.1 | 76 | 2.2 |
| 20 mg three times daily | **0.9** | 0.8 | 98 | 0.5 |

**(3) Same dose, differing only in whether a tumour is present (A3).** The *fraction*
of the pump pool inhibited is almost the same while the *absolute* acid output
remaining is five-fold — because factors 1 and 3 are raised.

| Arm | Inhibited pump fraction | Fasting BAO | pH>4 (%) |
|----|------------------|----------|----------|
| Normal stomach, omeprazole 20 mg qd | 0.834 | **1.5** | 61 |
| ZES, omeprazole 20 mg qd | 0.867 | **7.5** | 41 |

**(4) The CYP2C19 gradient (A5)** — the same 30 mg bd. Fasting BAO: UM **3.3** → NM
**2.0** → PM **0.4**. Monotonic.

**(5) The H2 blocker is outcompeted, and the receptor is upregulated (A6).**
Famotidine 40 mg four times daily brings fasting BAO down only from 33.7 to 27.9
(17% inhibition), and tripling it to 120 mg four times daily stops at 21.4 (36%). At
the same time H2 receptor density rises from 1.35 to 1.58 (tolerance). It is because
only one of the three stimulatory pathways has been blocked.

**(6) The secretin sign inversion and its false positive (A7)** — the cleanest
held-out result in this model.

| Subject | Baseline gastrin | Peak | Nadir | Δ | Test result |
|------|---------------|------|------|---|-----------|
| Normal | 51 | 51 | 30 | 0 (**−21**) | negative |
| **Normal + PPI for 4 weeks** | **125** | 125 | 54 | 0 (**−71**) | **negative** |
| Sporadic ZES | 896 | 1712 | 893 | **+817** | positive |
| MEN1-ZES | 678 | 1289 | 674 | **+612** | positive |

The PPI raises gastrin 2.45-fold (51 → 125) and so makes the patient **look like a
gastrinoma**, but the secretin test **stays negative.** The discrimination actually
used in the clinic falls straight out of two edges of opposite sign.

**(7) Parathyroidectomy — a difference of differences against a paired control (A8).**
Because the tumour grows in both arms over six months (gastrin rises in both arms), a
control is indispensable.

| Arm | Ca²⁺ after | PTH after | Gastrin before → after | BAO before → after |
|----|---------|--------|------------------|-------------|
| 17 parathyroidectomy | 1.25 | 4.0 | 830 → 985 | 45.2 → 37.1 |
| 17b control (no surgery) | 1.41 | 4.5 | 830 → 1229 | 45.3 → 52.9 |
| **Difference of differences** | **−0.16** | −0.5 | **−244 pg/mL** | **−15.8 mEq/h** |

**(8) The antitumour arms — an identical metastatic baseline (A10).**

| Arm | TTP (months) | Literature | Best RECIST (%) | Tumour dose (Gy) | Renal dose (Gy) |
|----|-----------|------|-----------------|---------------|---------------|
| PPI alone (natural history) | **4.0** | 4.6~5.5 (placebo arms) | 0.0 | – | – |
| Octreotide LAR | **7.2** | — | 0.0 | – | – |
| Everolimus 10 mg | **8.9** | 11.0 (RADIANT-3) | 0.0 | – | – |
| Sunitinib 37.5 mg | **11.4** | 11.4 | 0.0 | – | – |
| CAPTEM ×12 | **21.9** | 22.7 (E2211) | **−32.7** | – | – |
| ¹⁷⁷Lu-DOTATATE ×4 | **12.1** | ~20.7 (OCLURANDOM) | −8.8 | 79.9 | **18.0** |
| PRRT, SSTR2-low | **6.5** | (loss of effect expected) | −0.6 | 21.4 | 19.3 |
| PRRT, no amino acid | 11.9 | — | −8.0 | 77.8 | **39.0** ✗ |

Sunitinib and CAPTEM hit the literature values almost exactly. **Take away amino acid
renal protection and the renal absorbed dose rises from 18.0 to 39.0 Gy, crossing the
23 Gy limit** — a result that was not coded in. In an SSTR2-low tumour the PRRT effect
essentially disappears (TTP 12.1 → 6.5).

**(9) Untreated natural history versus treatment (A11).** After five years:

| Arm | pH | BAO | Ulcer index | Oesophageal injury | 5-year bleeding/perforation risk | B12 | Mg | Gastric NET (cm³) |
|----|----|-----|----------|----------|--------------------|-----|----|--------------|
| 5 y omeprazole 30 bd | 5.53 | 0.9 | 0.00 | 0.00 | **0.004** | 0.42 | 0.75 | 0.07 |
| 5 y MEN1 omeprazole | 5.17 | 1.0 | 0.00 | 0.00 | 0.004 | 0.46 | 0.75 | **0.19** |
| 5 y untreated | 1.28 | 62.7 | **0.62** | **0.60** | **0.861** | 1.01 | 0.85 | 0.07 |

Two facts come out together: that pre-PPI-era ZES deaths were mostly ulcer
complications, and that the type 2 gastric neuroendocrine tumour burden is 2.7-fold
greater only in MEN1.

**(10) A somatostatin analogue combined with a PPI is slightly subadditive (A12).**
Baseline 33.7 → PPI alone 2.0, SSA alone 20.7, combination 1.8. If the two effects
were independent 1.2 would be expected, so observed − expected = **+0.52 mEq/h**.

### Held-out predictions — **those that failed** (with causes stated)

1. **Treated ZES is over-suppressed (the largest failure).** On omeprazole 30 mg bd
   the model gives a fasting acid output of **2.0 mEq/h**. The actual value from
   maintenance-therapy trials is **about 6 mEq/h**. The cause is the
   **acid-activation trap**. The prodrug needs the acid that active pumps make in
   order to convert to the sulfenamide, so in a hypersecretory stomach with many
   active pumps the drug is **more** efficient per mg. Pharmacologically right and
   clinically wrong. An attempt to paper over it by fitting `GSYN` collapsed the
   calibration (see above), so **it is left unfixed.** Two knock-on effects are
   recorded with it.
   - **The A4 titration curve does not reproduce "ZES needs several times the
     reflux-disease dose" in absolute terms.** At 20 mg per day the target
     (<10 mEq/h) is already met. But **the ordering is right**: at the same dose
     fasting BAO is normal 0.8 < sporadic ZES 2.8 < MEN1-ZES 3.2, and nocturnal
     pH<4 time falls monotonically with dose, 4.2 → 3.0 → 2.0 → 0.8 h.
   - **The long-term harms in A11 are overestimated.** Because the treated stomach
     becomes almost achlorhydric at pH 5.5, five-year cobalamin falls as far as 0.42.
     Real ZES patients are not completely achlorhydric even on treatment.
2. **The P-CAB (vonoprazan) comes out far too weak.** In a normal subject the pH>4
   holding time on 20 mg once daily is **34%**, whereas the literature value is above
   90% and superior to esomeprazole (Kagami 2016). In ZES it gives a fasting BAO of
   19.0, far behind omeprazole 30 bd (2.0). **The cause is clear:** the model has no
   separate compartment for the P-CAB's **intracanalicular ion-trapping reservoir**,
   so the effect tracks the plasma concentration directly, and with once-daily dosing
   at a half-life of 7.75 hours the blockade all but disappears at the trough (20~24
   hours post-dose). The right fix is to put a canalicular accumulation compartment
   ahead of `PUMPR`.
3. **Almost no rebound appears after withdrawal.** The rebound ratio (fasting BAO
   after withdrawal / before withdrawal) is **1.07** in normals and **0.66** in ZES
   (1.14 at four weeks). The literature value is 1.5~2-fold. The direction is right
   but the magnitude is short. The pump pool at the moment of withdrawal is genuinely
   enlarged (1.18-fold in normals, 1.33-fold in ZES), but the covalently bound `PUMPB`
   has to degrade with a 50-hour half-life before it vacates its place, so **an
   explosive release is structurally impossible.** Real rebound appears to depend more
   heavily on ECL hyperplasia persisting for several weeks.
4. **The TTP of everolimus and of PRRT is underpredicted** (8.9 vs 11.0 months, 12.1
   vs ~20.7 months). Sunitinib (11.4 vs 11.4) and CAPTEM (21.9 vs 22.7) are right, so
   the problem is not the growth model but the effect size of these two agents. For
   PRRT it is most likely that the effective α of low-dose-rate delivery was left as a
   single constant (`KRADK`).
5. **The residuals of the A7 and A9 anchors themselves** — see the ✗ entries in the
   "Fitted anchors" table above.

### A13 — the virtual population · where the variability disappears

The result of running three random effects (parietal cell trophic gain · tumour
secretion rate · ECL trophic gain) together from the run-in onwards. n = 40.

| | Median BAO | IQR | Range | Median FSG (IQR) | PCM IQR |
|--|-----------|-----|------|------------------|---------|
| **Untreated** | 33.4 | 28.0–39.1 | 20.6–53.7 | 752 (618–1119) | 1.59–2.01 |
| Omeprazole 30 mg bd | 1.9 | 1.8–2.2 | — | 903 | — |

In the untreated state, between-individual differences span 2.6-fold (20.6 → 53.7
mEq/h). That the IQR collapses to 1.8–2.2 after treatment and the target attainment
rate becomes 100% is **not because there is no variability but because of the ceiling
effect created by the over-suppression failure above.** With the pump pool 94%
inhibited, wherever you start you arrive at the same place. (In an earlier version the
random effects were injected after a deterministic run-in and only 10 days simulated,
so there was no variability **even in the untreated state** — the time constants of
the trophic gains are in weeks, so the run-in and the treatment period have to be
combined into a single simulation.)

### A0 — reference states (all outputs)

| State | Fasting BAO | MAO | BAO/MAO | pH (24h mean) | Gastrin (24h mean) | PCM | ECL | Ca²⁺ | Diarrhoea index |
|------|----------|-----|---------|---------------|---------------------|-----|-----|------|----------|
| Normal | 3.0 | 22.5 | 0.13 | 2.74 | 63 | 0.98 | 0.98 | 1.25 | 0.00 |
| Sporadic ZES | 33.7 | 55.9 | **0.60** | 1.33 | 902 | 1.73 | 1.95 | 1.25 | 0.44 |
| MEN1-ZES | 39.5 | 66.0 | **0.60** | 1.32 | 683 | 1.68 | 1.89 | **1.41** | 0.50 |
| Metastatic | 58.6 | 63.0 | **0.93** | 1.30 | 8005 | 1.86 | 2.14 | 1.25 | 0.62 |

The `BAO/MAO > 0.6` diagnostic criterion was not fitted. Both sporadic and MEN1 come
out at exactly 0.60 and the metastatic case at 0.93. The normal is 0.13.

---

## Thirty-seven scenarios — a paired control for every claim

| Group | Scenarios | What is held fixed and what is varied |
|------|----------|---------------------------------|
| Reference states | 01–04 | Normal / sporadic ZES / MEN1-ZES / metastatic. No drug |
| Interval vs dose | 05–07 | **The same 60 mg per day** split qd / bd / tid |
| Tumour or no tumour | 08–09 | **The same omeprazole 20 mg qd** in a normal stomach and in ZES |
| CYP2C19 | 10–11 | **The same 30 mg bd** in UM / NM / PM |
| Drug class | 12, 12b, 13–15 | PPI / P-CAB / H2RA at licensed doses, with a normal control |
| Attacking factors 1 and 3 | 16–19 | Octreotide · parathyroidectomy · duodenotomy · netazepide |
| Parathyroidectomy control | 17 vs 17b | **The same patient**, differing only in surgery (difference of differences) |
| Combination | 20 | PPI + somatostatin analogue |
| Withdrawal and rebound | 21–22 | Withdrawal after 8 weeks of dosing, in a normal subject and in ZES |
| Long horizon | 23–25 | 5 years treated / 5 years MEN1 treated / 5 years untreated |
| Antitumour | 26–33 | Six arms + SSTR2-low + no amino acid, all from **the same metastatic baseline** |
| Concomitant exposure | 34–35 | NSAID / impaired renal function |

---

## Defects found and fixed during integration (a record)

A mechanism looking plausible and equations actually computing are two different
things. Here is what was found and fixed in this version. Each item is written up
together with **the wrong output** it produced.

1. `posf()` is not exposed in mrgsolve 2.0.1 → replaced by a local `pz()` helper.
2. `STRICT` is a DOT reserved word → node renamed (the map would not render).
3. `SETINIT` is an mrgsolve reserved word → changed to `ICFLAG`.
4. **Assigning `CMT_0` in $MAIN silently invalidated `init()`.** Every scenario
   restarted the trophic feedback from a healthy stomach. In a patient with gastrin
   900 pg/mL, parietal cell mass stayed at 1.08 and ECL at 1.14, untreated ZES acid
   output came out as **15 mEq/h** rather than 40, and the maximal acid output of a
   gastrinoma patient printed as **18.94, completely identical** to that of a healthy
   subject. Resolved with an `ICFLAG` guard.
5. **A 40-week run-in doubled the tumour.** A phenotype specified at 900 pg/mL came
   out at **2355 pg/mL**. Tumour growth was frozen during the run-in (`GRADEF = 0`) so
   that `TUM0` means "the volume at presentation".
6. **`c(pheno, extra)` produced duplicate parameter names.** Because `param()` silently
   discards duplicates, scenario 32 ("PRRT, SSTR2-low") ran with `SSTR2D = 1.0`, that
   is, as **exactly the same simulation** as scenario 31. Replaced with a
   `modifyList`-based `.pp()`.
7. **pH could not reach the achlorhydric limit** (pH 3.3 under complete blockade) →
   titration curve.
8. **Gastric volume ignored the secreted volume.** Luminal [H⁺] went as high as
   **178 mmol/L**, beyond the physiological ceiling. The isotonic secretion that enters
   along with the acid was included in `VG`.
9. **A bimolecular neutralisation term turned the remaining buffer into an infinitely
   strong base.** Postprandial pH in untreated ZES came out at **6.3**. The
   neutralisation flux was saturated with respect to H⁺ and limited to the buffer
   capacity.
10. **The trophic drive `EXC` was rectified through `pz()`.** `gTD` oscillates around 1
    over the meal cycle, but only the positive half was being counted, and a healthy
    stomach with no tumour at all drifted to parietal cell mass **1.22** and ECL
    **1.27**. Changed to a signed `EXC` (floor −0.9).
11. **`GBNORM` was not the model's own healthy gastrin** → TD became 1.33 in the
    healthy state and the same drift arose by a different route. Resolved with a fixed
    point (anchor A6).
12. **Counterfactual C reported a headroom of −0.67 for an untreated patient.** It
    redistributed the pump pool to the trafficking equilibrium rather than to the
    current active:inactive ratio, so it did not converge onto the real trajectory when
    no drug was present.
13. **Oesophageal injury saturated in every arm.** Hung on "pH < 4", `ESOPH ≈ 0.75`
    even in arms whose acid secretion was completely controlled (a treated ZES stomach
    is still acidic). Changed to an excess-acid-load criterion.
14. **A healthy subject's diarrhoea index was 0.32** → changed from absolute acid
    output to the excess over mucosal handling capacity.
15. **A healthy subject developed ulcers in long-term simulation** → the mucosal
    tolerance threshold `ATHR` was below the normal postprandial acid output (9 → 14
    mEq/h).
16. **The dosing schedule started at midnight.** The pre-breakfast fasting window
    (03:00–07:00) then fell 3–7 hours after a once-daily dose, that is at **the moment
    of maximum effect**, so at identical total mg **once daily beat twice daily** —
    precisely the opposite of the model's central claim, and an artefact of clock
    alignment rather than of pharmacology. Changed to a 07:00 start.
17. **The nocturnal breakthrough index saturated at 8 hours in every arm**, so regimens
    could not be ranked → the nocturnal mean acid output is now reported alongside it.
18. **The secretin test reported 0 in a healthy subject.** When the response is
    inhibitory, using `max()` alone returns the baseline value, so an "opposite
    response" looks like "no response" → the nadir is now reported too.
19. **Parathyroidectomy was reported as raising gastrin.** It was being compared with
    the patient's own pre-operative state across four months of tumour growth. A paired
    unoperated control (17b) was added and the result reported as a **difference of
    differences**.
20. **The virtual population had no variability.** The three random effects apply to
    trophic gains whose time constants are in weeks, but they were injected after a
    deterministic run-in and only 10 days simulated, giving an interquartile range of
    3.6–4.2 mEq/h and a target attainment rate of **100%**. The run-in and the
    treatment period were combined into **a single simulation** per subject.
21. **`$OMEGA` was declared and not used** (dead ETA). The three random effects were
    applied in $MAIN to parietal cell trophic gain · tumour secretion rate · ECL
    trophic gain.
22. Calibration collapse (see "Methodological decisions" above) — the `GSYN` pairing was
    withdrawn.

---

## How to run

```bash
# 1. Render the map
dot -Tsvg zes_qsp_model.dot -o zes_qsp_model.svg
dot -Tpng -Gdpi=150 zes_qsp_model.dot -o zes_qsp_model.png

# 2. Build the model + run all 37 scenarios and the whole A0-A13 suite
Rscript -e 'options(zes.run.scenarios=TRUE); source("zes_mrgsolve_model.R")'

# 3. Re-run the calibration (coordinate-wise fixed-point iteration)
Rscript -e 'source("zes_mrgsolve_model.R"); print(zes_calibrate())'

# 4. Dashboard
Rscript -e 'shiny::runApp("zes_shiny_app.R")'
```

Packages required: `mrgsolve` (built and verified on 2.0.1), `shiny`, `ggplot2`,
`dplyr`, `tidyr`, and optionally `DT`. **All 17 outputs** of the Shiny app were
confirmed to render headlessly with `shiny::testServer`.

Run time for the whole suite (4 cores): run-in + A0 about 12 s, the 37 scenarios about
40 s, A1–A12 immediate, the A4 titration about 20 s, A13 (n=40) about 3 min.

---

## ⚠️ Disclaimer

This is a semi-quantitative QSP model for education and research. It was built from
the published literature and clinical trial data, but it has not been independently
verified or certified, and **it must not be used for real clinical decision-making,
prescribing, or regulatory submission.**
