# Toxic Alcohol Poisoning (Methanol · Ethylene Glycol) — QSP Model
### Quantitative Systems Pharmacology

| Deliverable | File |
|--------|------|
| 🗺️ Mechanistic map | [`tap_qsp_model.dot`](../../../toxic-alcohol-poisoning/tap_qsp_model.dot) · [SVG](../../../toxic-alcohol-poisoning/tap_qsp_model.svg) · [PNG](../../../toxic-alcohol-poisoning/tap_qsp_model.png) — 162 nodes / 246 edges / 20 clusters |
| ⚙️ mrgsolve model | [`tap_mrgsolve_model.R`](tap_mrgsolve_model.R) — 50 ODEs, time unit = hours, 25 scenarios |
| 🐍 Independent verification implementation | [`tap_python_reference.py`](../../../toxic-alcohol-poisoning/tap_python_reference.py) · [full output](../../../toxic-alcohol-poisoning/tap_reference_output.txt) — 15 verification sections |
| 📊 Shiny dashboard | [`tap_shiny_app.R`](tap_shiny_app.R) — 12 tabs |
| 📚 References | [`tap_references.md`](tap_references.md) — 88 PubMed entries (`mkrefs.py` queries them live; no PMID written from memory) |

---

## The organising thesis

**The parent alcohol is not the poison.** Methanol and ethylene glycol are, to a
good approximation, ethanol-like sedatives with an unusually large osmotic
footprint, and in themselves they have almost no organ toxicity. Every lesion in
this disease lies in **the acidic metabolites the patient's own liver
manufactures**.

```
methanol           --ADH-->  formaldehyde     --ALDH-->  FORMATE
ethylene glycol    --ADH-->  glycolaldehyde   --ALDH-->  glycolate  -->  oxalate
```

So this disease is **the integral of a flux, not a concentration**, that flux is
a single Michaelis–Menten expression, and **every therapy ever used owns exactly
one term of that expression, one each.**

```
                       Vmax_i · (S_i/Km_i)
   v_i  =  ─────────────────────────────────────────      ① generation
            1 + Σ_j (S_j/Km_j) + [FOM]/Ki                    fomepizole and ethanol
                                                             own the denominator

        −  CL_ren · (1 − f_reabs(pH_urine)) · ACID           ② renal excretion
                                                             bicarbonate owns f_reabs

        −  CL_hd · ACID                                      ③ dialysis
                                                             the machine owns this term alone

        −  Vmax_THF · THF · ACID/(Km + ACID)                 ④ folate-dependent oxidation
                                                             folinic acid owns Vmax
```

And **the injury is the CNS integral, not the plasma integral**, and inside its
threshold sits a term nobody prescribes but everybody measures.

```
   d(ACID_cns)/dt = PS · ( f_HA(pH_plasma)·C_plasma − f_HA(pH_brain)·C_cns )
        f_HA(pH) = 1/(1 + 10^(pH − pKa))     formate pKa 3.75 · glycolate 3.83
```

From this one structure the following follows **as arithmetic, not as
assertion**. Every number below was printed by
[`tap_python_reference.py`](../../../toxic-alcohol-poisoning/tap_python_reference.py),
and the full text is in
[`tap_reference_output.txt`](../../../toxic-alcohol-poisoning/tap_reference_output.txt).

---

## ① The two gaps are the same reaction read twice, and their sum is the dose

When 1 mol of formate is generated, the 1 mol of protons it carries titrates
1 mol of bicarbonate and that bicarbonate leaves the body as CO₂. Sodium does
not move. So **the osmolal gap falls by as much parent alcohol as has
disappeared and the anion gap rises by as much acid as has been generated.** The
**ratio** of the two reads elapsed time rather than severity, and the **sum** of
the two is an estimate of the amount ingested that is independent of when the
patient arrived.

Methanol 0.7 g/kg (1529 mmol), untreated:

| t (h) | Methanol mg/dL | Formate mM | OG mOsm | ΔAG mEq | **OG/ΔAG (the clock)** | Sum | pH |
|---|---|---|---|---|---|---|---|
| 0.5 | 88.7 | 1.18 | 28.1 | 0.8 | **37.5** | 28.9 | 7.360 |
| 2 | 101.4 | 4.43 | 32.8 | 3.3 | **9.85** | 36.2 | 7.334 |
| 8 | 61.9 | 13.91 | 22.6 | 12.4 | **1.83** | 34.9 | 7.259 |
| 16 | 22.8 | 23.03 | 15.1 | 18.8 | **0.81** | 33.9 | 6.985 |
| 24 | 4.5 | 25.62 | 12.7 | 19.6 | **0.65** | 32.3 | 6.694 |
| 48 | 0.0 | 18.61 | 12.4 | 9.8 | **1.26** | 22.3 | 7.300 |

**Three things come out of this, and one of them is the opposite of the
textbook.**

1. **The clock works — but only while it is going down.** From 37.5 at 1 hour to
   0.65 at 24 hours the ratio reads elapsed time and nothing else. But once the
   reaction is over, formate is cleared, the anion gap closes and the ratio
   **rises again** (0.65 → 1.26). A patient whose ratio is near 1 is ambiguous
   between "about halfway through" and "finished and recovering", and what
   distinguishes them at that point is **the sum**.
2. **The sum is a lower bound on the dose, not an estimate.** At 1 hour it reads
   35.0 against an actual ingestion of 36.4 mmol/L — 96% recovery. At 1.9 hours
   it is 4% high (the bicarbonate-space effect of §4), and by 48 hours it has
   decayed to 31% of its initial value. Read early it is quantitative; read late
   it underestimates.
3. **The osmolal gap does not return to zero — and this is a derived
   conclusion.** At 24 hours, when methanol has effectively all gone, the OG
   remains at 12.7 mOsm. The textbook says "formate replaces bicarbonate 1:1, so
   the osmolal gap normalises", but that is true **only if the two substances
   share a volume of distribution**, and they do not. Formate spreads through
   35 L while the bicarbonate it titrates is buffered across 45–73 L. So the
   rise in plasma formate is **larger** than the fall in plasma bicarbonate (the
   same fact as a delta-ratio > 1 in an organic acidosis), and that difference
   appears as a residual osmolal gap. **A methanol patient 24 hours in is
   predicted to have a moderately raised osmolal gap and a wide anion gap at the
   same time**, and a clinician who reads that remaining gap as evidence of
   unabsorbed alcohol will choose the wrong treatment.

### Why the sum is not exactly conserved — the bicarbonate space

```
apparent bicarbonate space = (0.40 + 2.6/[HCO3]) L/kg
```

| HCO3 mM | Space (70 kg) | Space / Vd(methanol) 42 L |
|---|---|---|
| 24 | 35.6 L | 0.847 |
| 16 | 39.4 L | 0.938 |
| 12 | 43.2 L | 1.028 |
| 8 | 50.8 L | 1.208 |
| 4 | 73.5 L | 1.750 |

The OG falls divided by 42 L and the ΔAG rises divided by the bicarbonate space.
At a normal HCO3 the space is 35.6 L, so **the ΔAG moves about 18% faster than
the OG and the sum drifts upward**; as the acidosis deepens and the space exceeds
42 L, the rise in ΔAG decelerates and the sum comes back. The deviation is
self-limiting and its size can be read off the table above.

---

## ② Why pH is a better prognostic factor than concentration — pH sits inside the transfer term

What crosses the membrane is **un-ionised formic acid only.**

| pH | f_HA plasma | Influx rate vs 7.40 | Equilibrium CNS/plasma ratio | vs 7.40 |
|---|---|---|---|---|
| 7.40 | 2.238e−04 | 1.00 | 0.398 | 1.00 |
| 7.20 | 3.547e−04 | 1.58 | 0.501 | 1.26 |
| 7.00 | 5.620e−04 | 2.51 | 0.631 | 1.58 |
| 6.90 | 7.074e−04 | 3.16 | 0.708 | 1.78 |
| 6.80 | 8.905e−04 | 3.98 | 0.795 | 2.00 |

**One equation, two effects, the same sign.** pH 7.40 → 6.80 makes the **influx
rate constant into the brain 3.98-fold** and the **equilibrium CNS/plasma ratio
2.00-fold**. *At the same blood formate*, the CNS burden differs two- to
three-fold. The pH range a clinician treats as "a number to be corrected" is in
fact the mechanism.

Bicarbonate reverses those two terms and at the same time moves **a third**.
Tubular reabsorption is non-ionic diffusion of the un-ionised acid, so the
reabsorbed fraction is a function of **urine pH**.

| HCO3 mM | Urine pH | Reabsorption | Renal formate clearance |
|---|---|---|---|
| 4 | 5.50 | 97.6% | 0.17 L/h |
| 12 | 6.13 | 93.4% | 0.47 L/h |
| 24 | 6.73 | 79.7% | 1.46 L/h |
| 28 | 7.66 | 32.8% | 4.84 L/h |
| 32 | 8.00 | 18.2% | 5.89 L/h |

**Note the direction of the trap: the sicker the patient (low HCO3, acid urine),
the more formate the kidney gives back.** This is a second self-amplifying loop,
and it is why "the acidosis is just a number" is wrong.

---

## ③ The same antidote, opposite management — one difference in the blocked-state half-life

Blocking ADH does not remove the poison; it **turns a self-eliminating poison
into a slowly leaking reservoir**. Whether that amounts to cure or to buying
time is decided by **a single clearance term**.

| Complete ADH blockade at 2 h, no dialysis | Peak | Blocked-state half-life | Time to <20 mg/dL |
|---|---|---|---|
| Methanol 0.7 g/kg | 105.1 mg/dL | **46.6 h** | 113 h |
| Ethylene glycol 90 mL | 200.6 mg/dL | **17.3 h** | 59 h |
| Ethylene glycol 90 mL, GFR 30% | 207.6 mg/dL | **54.7 h** | > 160 h |

The literature values are 43–52 hours for methanol and about 17 hours for
ethylene glycol. **No parameter was tuned separately to hit both values** — the
difference comes entirely out of `CLRE` alone (the term by which 20–30% of
ethylene glycol leaves unchanged in the urine). That is why ADH blockade can be
**definitive monotherapy** in ethylene glycol with intact renal function and is
**buying time** in methanol. And the third row is **the entire reason for
dialysing an ethylene glycol patient who is already in AKI**: once the poison has
taken the kidney, ethylene glycol starts behaving like methanol.

---

## ④ Competitive blockade is diluted by its own substrate — read the table, not the slogan

`v_i = Vmax_i·(S_i/Km_i) / (1 + Σ_j S_j/Km_j + F/Ki)`
Km (mM): ethanol 1.0 < ethylene glycol 6.0 < methanol 8.0 · fomepizole Ki = 0.15 µM

| Condition | Inhibition factor | Remaining flux | Flux mmol/h |
|---|---|---|---|
| Methanol 100, nothing else | 1.0 | 100% | 95.5 |
| Methanol 100 + ethanol 100 mg/dL | 5.4 | **18.4%** | 17.6 |
| Methanol 100 + ethanol 150 mg/dL | 7.6 | 13.1% | 12.5 |
| Methanol 100 + fomepizole 10 µg/mL | 166.7 | **0.60%** | 0.57 |
| Methanol 20 + ethanol 100 mg/dL | 13.2 | 7.6% | 4.0 |
| Methanol 20 + fomepizole 10 µg/mL | 457.1 | 0.22% | 0.12 |
| Methanol 400 + ethanol 100 mg/dL | 2.3 | **43.3%** | 48.9 |
| Methanol 400 + fomepizole 10 µg/mL | 49.9 | **2.0%** | 2.26 |

**"Fomepizole is 1000 times ethanol" is a false statement at the bedside.** The
inhibition factor is `(1 + ΣS/Km + I/Ki)/(1 + ΣS/Km)`, so the more substrate
there is the smaller it becomes. At a methanol of 100 mg/dL, ethanol at its
target concentration still lets **about a fifth** of the acid-generating flux
run, while fomepizole leaves **less than 1%** — the honest ratio is about
**31-fold**. And at a methanol of 400 mg/dL even fomepizole leaves 2.0%. This is
the arithmetic of "blockade is never absolute".

---

## ⑤ The only threshold in this model — calcium oxalate, hence rate-dependent toxicity

Formate injures in proportion to an integral, so **delivering the same dose more
slowly gives the same injury.** Oxalate precipitates only when it exceeds the
solubility product, so **delivering the same dose more slowly gives less
injury** — because the same oxalate simply leaves in the urine.

Ethylene glycol 90 mL, total dose identical, only the flux differing:

| | Peak oxalate µM | Renal CaOx mmol | GFR nadir % | Ionised Ca nadir mM |
|---|---|---|---|---|
| No treatment | 315.3 | 44.6 | **26** | 0.821 |
| Fomepizole 2 h | 56.9 | 5.1 | **95** | 1.127 |
| Fomepizole 6 h | 114.8 | 15.2 | **67** | 1.019 |
| Fomepizole 12 h | 201.8 | 28.7 | **43** | 0.914 |

**There is no corresponding table for methanol, and that asymmetry is the
point.** The sentence "fomepizole bought time" is a far stronger claim in
ethylene glycol than in methanol.

Because 1 mmol of deposited crystal takes 1 mmol of calcium with it, **ionised
calcium is the stoichiometric shadow of the crystal burden** and the QT reads
that shadow.

---

## Verification refuted two of the author's own claims

As is this repository's practice, a refuted claim is not quietly deleted but left
in place **as refuted**.

### Refutation A — "dialysis removes fomepizole, so the flux jumps mid-session"

Cluster 17 of the map draws this loop, and the model was built to display it.
**The model does not say that.** What matters is not how fast dialysis removes
the antidote but **the margin C/K between the therapeutic concentration and the
inhibition constant**.

| Antidote | Therapeutic concentration | Inhibition constant | **Margin C/K** | Dialysis k | Dialysis time needed to spend the whole margin |
|---|---|---|---|---|---|
| Fomepizole | 0.1218 mM | 1.5e−04 mM | **812** | 0.408/h | **16.4 h** |
| Ethanol | 21.71 mM | 1.0 mM | **22** | 0.531/h | **5.8 h** |

A conventional 4–6 hour session cannot spend fomepizole's margin, and spends most
of ethanol's. Checking that against the model:

| Comparison | redosing → no redosing | Change in oxidised dose | Change in CNS exposure |
|---|---|---|---|
| Fomepizole, 6 h IHD | identical (to 3 figures) | **+0.3%** | +0.1% |
| Ethanol, 6 h IHD | differs | **+2.4%** | +1.4% |
| Fomepizole, 30 h CRRT | identical | **+0.0%** | +0.0% |

Look inside the session and the mechanism is visible. The fomepizole arm with
redosing omitted still has `C/Ki = 173` after six hours and the flux barely moves,
0.132 → 0.174 mmol/h (against 67 with no blockade). The ethanol arm falls from
121 to 27 mg/dL, and the flux at the end of the session differs
**3.74 versus 1.36 mmol/h, 2.8-fold**.

So the corrected conclusion: **"increase the antidote during dialysis" is a
generous safety margin in fomepizole and a quantitative necessity in ethanol** —
the opposite direction from the way the two drugs are usually taught.

And one observation that arrived without being designed for: the reason the
total-dose effect is small is **that dialysis partly defends against its own
removal of the antidote**. The same session is pulling the substrate out as well,
so two terms of the same formula pull in opposite directions. To add it
honestly, the probability of death was **slightly lower, if anything**, in the
non-escalated ethanol arm (0.143 versus 0.155). A 2.5-fold escalation buys a
little blockade while buying more sedation than that, and sedation sits inside
the respiratory-failure hazard function. The model is not saying "it is better to
underdose ethanol"; it is saying that **the benefit and the harm of ethanol live
inside the same dose, and fomepizole's do not**.

### Refutation B — "pH is the best prognostic factor"

Across 16 runs crossing dose with time of arrival:

| Variable at admission | Spearman ρ with P(death) |
|---|---|
| Methanol concentration at admission | **−0.135** |
| Formate concentration at admission | **+0.935** |
| pH at admission | **−0.841** |

**What was confirmed, and confirmed strongly:** the alcohol concentration at
admission is useless and even **runs with the wrong sign.** That is not noise but
structure — a high methanol means the reaction has not happened yet, and that is
exactly the patient the antidote can save. A rule that grades severity from the
alcohol concentration alone sorts patients **by salvageability** and then calls
that severity.

**What was not confirmed:** "pH is best" is not true. **Blood formate is
better** (ρ +0.94 versus −0.84). Naturally so — formate is the proximate cause
and pH is one step downstream of it. The reason pH rules the bedside is not that
it is the better variable but **that formate cannot be measured within the four
hours in which the decision has to be made**. pH's residual advantage over
formate is **confined to a single transfer term**: two patients with the same
formate are not poisoned to the same degree, and the more acidotic one is worse.
Stated that carefully it is defensible; said as "pH is the best prognostic
factor" it is an artefact of assay availability.

---

## A control for falsifiability — an alcohol whose metabolite is not an acid

Isopropanol is not simulated as a separate compartment. It is simulated **by
deleting one term** (scenario M17, parameter `NOACID = 1`). ADH turns over as
before and the metabolite becomes acetone (a ketone), so neither a proton nor an
anion is generated.

| | Acid-forming case (M1) | Ketone-forming case (M17) |
|---|---|---|
| **Methanol oxidised (mmol)** | **1353** | **1353** ← the flux is identical |
| Peak osmolal gap | 33.7 mOsm | **32.8 mOsm** |
| Peak anion gap | 31.8 mEq | **12.0 mEq** |
| pH nadir | 6.680 | **7.399** |
| HCO3 nadir | 4.2 mM | **24.0 mM** |
| Peak plasma formate | 25.62 mM | **0** |
| Peak CNS formate | 23.00 mM | **0** |
| Basal ganglia injury index | 0.930 | **0** |
| Final logMAR | 1.90 | **0** |
| P(death) | 1.000 | **0.001** |

The same alcohol, the same ADH flux, the same osmolal gap, one term deleted — and
the disease disappears. That is isopropanol poisoning: a large osmolal gap, a
drunk patient, ketones, and **a normal anion gap**, managed with a bed and time.
**A methanol model that does not do this when the acid term is removed is
modelling a different molecule.**

The second control is in cluster 20 of the map: **propylene glycol** — the
osmolal gap and the anion gap both rise, but the acid is ordinary L-lactate and
the cause is the patient's own lorazepam infusion. The gaps are not specific.

---

## An expected non-additivity — thiamine · pyridoxine versus fomepizole (holds only in part)

Ethylene glycol 90 mL:

| | Peak oxalate µM | Renal CaOx mmol | GFR nadir % |
|---|---|---|---|
| Nothing | 315.3 | 44.6 | 26 |
| Thiamine + pyridoxine only | 140.9 | 28.7 | 29 |
| Fomepizole 2 h | 56.9 | 5.1 | 95 |
| Fomepizole 2 h + both cofactors | 48.3 | 2.8 | 95 |
| Fomepizole 12 h | 201.8 | 28.7 | 43 |
| Fomepizole 12 h + both cofactors | 94.8 | 13.8 | 50 |

| Cofactor effect on crystal burden | Absolute (mmol) | Relative (%) |
|---|---|---|
| Alone | −16.0 | 35.8% |
| Added to fomepizole 2 h | **−2.4** | 46.1% |
| Added to fomepizole 12 h | **−14.9** | 51.9% |

**It holds in absolute terms and fails in relative terms, and the difference
matters.** Thiamine and pyridoxine widen two escape routes out of the glyoxylate
node, so what they own is **a proportion multiplying the flux**. Proportions are
largely scale-invariant, which is why the percentage column looks flat, and why
**percentages are the wrong way to read this table**. The absolute effect tracks
the residual flux exactly: 16 mmol with no antidote, 15 mmol together with a
**late** antidote, 2.4 mmol together with an **early** one.

The corrected statement is narrower than what the map implies, and more useful.
The cofactors are not redundant with fomepizole in general, they are redundant
**only with fomepizole given immediately** — because only immediate blockade
leaves no glyoxylate to redistribute. For the patient who came late — the patient
for whom the antidote has already become largely pointless — two cheap and
harmless drugs retain essentially their full solo effect. **A trial that gives
everyone fomepizole within 2 hours and then randomises thiamine and pyridoxine is
a trial designed to find nothing; the same trial in late patients is not.**

Look too at what the cofactors do **not** fix: in the untreated arm the GFR moves
only from 26 to 29%. Most of that injury is glycolaldehyde acting directly on the
tubule, and that is **upstream** of the branch point on which these two act.

---

## The price of delay — the antidote window is a gradient, not a rule

Methanol 0.7 g/kg, fomepizole only (no dialysis, no bicarbonate):

| Start h | pH nadir | Peak formate mM | Peak CNS formate | logMAR | P(death) | % of dose already oxidised |
|---|---|---|---|---|---|---|
| 0 | 7.357 | 0.13 | 0.04 | 0.00 | 0.00 | 12 |
| 2 | 7.334 | 4.43 | 1.00 | 0.00 | 0.00 | 21 |
| 6 | 7.290 | 10.98 | 3.81 | 0.05 | 0.01 | 38 |
| 8 | 7.259 | 13.91 | 5.33 | 0.68 | 0.38 | 46 |
| 12 | 7.169 | 19.12 | 8.69 | 1.38 | 0.82 | 61 |
| 14 | 7.094 | 21.27 | 10.79 | 1.50 | 0.88 | 67 |
| 18 | 6.839 | 24.31 | 17.09 | 1.64 | 1.00 | 77 |
| 24 | 6.680 | 25.62 | 23.00 | 1.73 | 1.00 | 85 |
| 30 | 6.680 | 25.62 | 23.00 | 1.75 | 1.00 | 88 |

**There is no cliff.** The outcome columns worsen smoothly, so the "antidote
window" is a marketing term attached to a gradient. And **the last column
saturates before the outcome columns do**: between 18 and 30 hours the fraction
of the dose already oxidised moves only from 77 to 88%, yet the difference
between those two patients is the difference between a survivor and a fatality.
What is still changing after the reaction is over is **not how much acid was
made but how long it stays**, and fomepizole does nothing about that. This is
why the late-arriving patient must be dialysed, and **why "the concentration was
only 40 mg/dL so the antidote was withheld" is a false sentence** — the
concentration at that point is the amount of the poison that is no longer going
to hurt the patient.

---

## Twenty-five scenarios (summary)

| Scenario | pH nadir | HCO3 | AG | Formate | CNS | Glycolate | CaOx | iCa | GFR% | PUT | logMAR | P(death) | P(blindness) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| M1 no treatment | 6.680 | 4.2 | 32 | 25.62 | 23.00 | — | — | 1.17 | 100 | 0.93 | 1.90 | 1.00 | 1.00 |
| M2 fomepizole 2 h | 7.334 | 20.7 | 15 | 4.43 | 1.00 | — | — | 1.18 | 100 | 0.00 | 0.00 | 0.00 | 0.00 |
| M3 fomepizole 14 h | 7.094 | 6.2 | 30 | 21.27 | 10.79 | — | — | 1.16 | 100 | 0.71 | 1.50 | 0.88 | 1.00 |
| M4 fomepizole 8 h + IHD | 7.259 | 11.6 | 24 | 13.91 | 5.03 | — | — | 1.15 | 100 | 0.04 | 0.01 | 0.11 | 0.00 |
| M5 IHD alone 15 h | 7.044 | 5.7 | 30 | 22.21 | 11.19 | — | — | 1.15 | 100 | 0.47 | 1.01 | 0.66 | 1.00 |
| M6 ethanol co-ingestion | 7.097 | 6.2 | 30 | 21.92 | 12.14 | — | — | 1.16 | 100 | 0.83 | 1.74 | 0.90 | 1.00 |
| M7 ethanol antidote 4 h | 7.288 | 17.0 | 19 | 7.86 | 2.70 | — | — | 1.18 | 100 | 0.00 | 0.00 | 0.21 | 0.00 |
| M8 bicarbonate alone | 7.259 | 11.6 | 32 | 24.79 | 10.93 | — | — | 1.13 | 100 | 0.81 | 1.70 | 0.77 | 1.00 |
| M9 full treatment at 6 h | 7.290 | 14.3 | 22 | 10.98 | 3.57 | — | — | 1.15 | 100 | 0.00 | 0.00 | 0.00 | 0.00 |
| M10 fomepizole 14 h + folinic acid | 7.094 | 6.2 | 30 | 21.27 | 10.77 | — | — | 1.16 | 100 | 0.67 | 1.42 | 0.87 | 1.00 |
| M11 no redosing during dialysis | 7.259 | 11.6 | 24 | 13.91 | 5.03 | — | — | 1.15 | 100 | 0.04 | 0.01 | 0.11 | 0.00 |
| M12 massive 1.5 g/kg + 2 × IHD | 7.273 | 16.2 | 20 | 9.31 | 2.70 | — | — | 1.15 | 100 | 0.00 | 0.00 | 0.02 | 0.00 |
| M13 CRRT 30 h | 7.259 | 11.6 | 24 | 13.91 | 5.06 | — | — | 1.15 | 100 | 0.06 | 0.10 | 0.22 | 0.00 |
| M14 ethanol + IHD, escalated | 7.259 | 11.6 | 24 | 13.91 | 5.20 | — | — | 1.15 | 100 | 0.04 | 0.02 | 0.15 | 0.00 |
| M15 ethanol + IHD, not escalated | 7.259 | 11.6 | 24 | 13.91 | 5.20 | — | — | 1.15 | 100 | 0.04 | 0.02 | 0.14 | 0.00 |
| M16 CRRT, no redosing | 7.259 | 11.6 | 24 | 13.91 | 5.06 | — | — | 1.15 | 100 | 0.06 | 0.10 | 0.22 | 0.00 |
| M17 control: metabolite is a ketone | 7.399 | 24.0 | 12 | 0.00 | 0.00 | — | — | 1.18 | 100 | 0.00 | 0.00 | 0.00 | 0.00 |
| E1 EG no treatment | 7.196 | 6.9 | 29 | 0.05 | 0.02 | 21.81 | 35.4 | 0.82 | **26** | 0.00 | 0.00 | 0.00 | 0.00 |
| E2 EG fomepizole 2 h | 7.343 | 21.9 | 14 | 0.04 | 0.01 | 2.91 | 5.1 | 1.13 | **95** | 0.00 | 0.00 | 0.00 | 0.00 |
| E3 EG fomepizole 12 h + IHD | 7.298 | 13.9 | 22 | 0.05 | 0.02 | 12.10 | 7.5 | 1.01 | **46** | 0.00 | 0.00 | 0.00 | 0.00 |
| E4 EG cofactors only | 7.213 | 7.5 | 28 | 0.05 | 0.02 | 21.54 | 23.0 | 0.94 | **29** | 0.00 | 0.00 | 0.00 | 0.00 |
| E5 EG fomepizole + cofactors | 7.343 | 21.9 | 14 | 0.04 | 0.01 | 2.91 | 2.8 | 1.15 | **95** | 0.00 | 0.00 | 0.00 | 0.00 |
| E6 EG in CKD (GFR 30%) | 7.343 | 21.8 | 14 | 0.04 | 0.01 | 3.02 | 8.3 | 1.09 | **95** | 0.00 | 0.00 | 0.00 | 0.00 |
| E7 EG ethanol antidote + IHD | 7.339 | 20.9 | 15 | 0.05 | 0.01 | 4.03 | 3.4 | 1.10 | **87** | 0.00 | 0.00 | 0.01 | 0.00 |
| E8 EG + empirical calcium loading | 7.298 | 13.9 | 22 | 0.05 | 0.02 | 12.10 | 7.6 | 1.10 | **46** | 0.00 | 0.00 | 0.00 | 0.00 |

Comparing E8 with E3, **empirical calcium loading does not meaningfully increase
the crystal burden** (7.6 versus 7.5 mmol). This differs from the author's prior
expectation and the reason is structural: the deposition rate is limited by
**oxalate delivery** and not by calcium. Drawing calcium supplementation as
"conditional" on the map is still correct, but the reason is not that it promotes
crystallisation — it is that **the indication is narrow**.

---

## Verification

**All 50 ODEs were independently re-implemented in Python/scipy**
([`tap_python_reference.py`](../../../toxic-alcohol-poisoning/tap_python_reference.py)).
All the two implementations share is the equations — different language,
different integrator, written separately. This work **found eight defects, three
of which changed a conclusion.**

| # | Defect | Symptom | Fix |
|---|---|---|---|
| 1 | Renal bicarbonate regeneration written as 0.10/h | a lethal dose of methanol gave pH 7.33 | it is an NH₄⁺ excretory capacity term, so 0.010/h |
| 2 | Tubular CaOx deposition written as a rate proportional to the excess over supersaturation | 50 mmol/h of crystal deposited, ionised calcium 0 in every EG patient | deposition cannot exceed oxalate **delivery** → a saturating fraction of the filtered load |
| 3 | Optic nerve injury driven linearly by vitreous formate | every patient with a long exposure to a low concentration went blind (the entire ethanol-antidote arm) | the retina fails for the same reason as the basal ganglia → drive it by a **retinal ATP deficit**, a threshold function |
| 4 | No flux threshold on proximal tubular injury | even the 2-hour fomepizole patient crashed to GFR 50% | the residual non-ADH flux never stops, so a Hill threshold is needed |
| 5 | Dialysis bicarbonate transfer written with the small-molecule clearance | 300+ mmol/h delivered, HCO3 47, PaCO2 79 | effective systemic clearance 3.0 L/h |
| 6 | Winter's formula extrapolated above a normal HCO3 | PaCO2 79 mmHg in an alkalotic patient | metabolic alkalosis is ~0.7 mmHg per mM |
| 7 | Ignoring that exogenous NaHCO₃ raises sodium | a calculated anion gap of −11 during a bicarbonate infusion | count the sodium and **bicarbonate therapy barely changes the anion gap** (clinically important behaviour) |
| 8 | Continuing to integrate below pH 6.6 | pH 5.55, CNS formate 120 mM (a cadaver's numbers) | truncate at the first non-survivable point |

### Mass conservation

| Scenario | Input (mmol) | Accounted (mmol) | Error |
|---|---|---|---|
| M1 no treatment | 1529.2 | 1529.2 | −0.001% |
| M4 fomepizole + IHD | 1529.2 | 1529.2 | −0.001% |
| M9 full treatment | 1529.2 | 1529.2 | −0.001% |
| M12 massive | 3276.9 | 3276.9 | −0.001% |
| E1 EG no treatment | 1614.2 | 1614.1 | −0.003% |
| E3 EG fomepizole + IHD | 1614.2 | 1614.1 | −0.003% |
| E6 EG in CKD | 1614.2 | 1614.2 | −0.001% |

`ALL CHECKS PASS` — no negative states, no leaks.

### R ↔ Python cross-verification

A 16-row reference trajectory (19 variables) for scenario M4 is embedded in
[`tap_mrgsolve_model.R`](tap_mrgsolve_model.R)
as `TAP_REFERENCE_M4`, and `tap_verify()` compares against it. The tolerance was
deliberately set tight at a relative 2e−3 (three significant figures) — because
the two implementations share only the equations, this is a real test and not a
tautology.

```r
source("tap_mrgsolve_model.R")
tap_verify()          # 16 x 19 = 304 comparisons
```

> **Note.** R is not installed in this container (package repository error), so
> the mrgsolve model was **not run**. The R file was written to correspond
> equation by equation with the Python reference implementation, and the
> acceptance test above is included so as to check that correspondence, but
> **execution verification has not yet been carried out.** Running
> `tap_verify()` first in an environment that has R is recommended. Every figure
> in this README comes from the Python implementation, which has been verified by
> execution.

---

## Calibration anchors

| Item | Literature | Model |
|---|---|---|
| Methanol elimination (at saturating concentration) | 8.5 mg/dL/h (4.4–25) | 7.3 mg/dL/h at 100 mg/dL |
| Methanol terminal half-life (unblocked) | ~2.5–3 h | 2.4 h |
| Methanol half-life (ADH blocked) | 43–52 h | **46.6 h** |
| Ethylene glycol half-life (ADH blocked, normal renal function) | ~17 h | **17.3 h** |
| Ethylene glycol renal excretion unchanged | 20–30% | 26% (FE) |
| Fomepizole therapeutic concentration | >8.2 µg/mL | 21 µg/mL after a 15 mg/kg loading dose |
| Fomepizole Ki | 0.1–0.5 µM | 0.15 µM |
| Ethanol target concentration | 100–150 mg/dL | inhibition factor 5.4 at 100 mg/dL |
| Dialysis clearance (small molecules) | 200–300 mL/min | 240 mL/min |
| CaOx monohydrate Ksp (37 °C) | 2.32e−9 M² | same |
| Apparent bicarbonate space | 0.40 + 2.6/[HCO3] L/kg | same |
| Glycolate in severe EG | 10–25 mM | 21.8 mM untreated |
| Formate in severe methanol | 10–20 mM (fatal range) | 25.6 mM untreated (0.7 g/kg) |

---

## Mechanistic map (20 clusters)

1 exposure · products · routes · absorption · 2 parent alcohol handling (where the
two diseases diverge) · 3 the ADH→ALDH oxidative cascade · 4 competitive ADH
blockade · 5 formate handling (folate) · 6 glyoxylate branching ratios ·
7 calcium oxalate supersaturation (the only threshold) · 8 mitochondrial toxicity
and lactate · 9 acid-base · 10 **the two-gap clock** ·
11 **pH-dependent ion trapping** · 12 CNS (basal ganglia · optic nerve) ·
13 kidney · 14 calcium and the heart · 15 fomepizole pharmacology · 16 ethanol
antidote therapy · 17 extracorporeal removal · 18 adjunctive therapy ·
19 clinical endpoints · 20 **differential diagnosis and the two controls**

Along the bottom of the map the six sentences of the thesis are drawn as bold
edges cutting across the clusters. The "the antidote is dialysed too" loop in
cluster 17 has been left exactly as it was, with **refutation A** above attached
to it.

---

## Reproduce

```bash
# map
dot -Tsvg tap_qsp_model.dot -o tap_qsp_model.svg
dot -Tpng -Gdpi=150 tap_qsp_model.dot -o tap_qsp_model.png

# verification (fully reproducible without R)
pip install numpy scipy
python3 tap_python_reference.py                 # all 15 sections
python3 tap_python_reference.py --section 9     # refutation A only

# re-query the references
python3 mkrefs.py --refresh

# R model and dashboard
Rscript -e 'source("tap_mrgsolve_model.R"); tap_verify(); print(tap_all())'
Rscript -e 'shiny::runApp("tap_shiny_app.R")'
```

---

## ⚠️ Disclaimer

This is a quantitative systems pharmacology model for educational and research
purposes. It was assembled from the open literature but has not been
independently verified or certified, and **must not be used for actual clinical
decision-making, prescribing, or regulatory submission.** In particular, several
of the conclusions in this document (refutations A and B, the scope of the
cofactor non-additivity, the residual osmolal gap, the ineffectiveness of
empirical calcium loading) are **predictions of the model and not facts
confirmed by clinical data.** They have been stated falsifiably in order to be
tested, not in order to be adopted.
