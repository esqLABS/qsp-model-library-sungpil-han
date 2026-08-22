# Primary Hyperoxaluria — QSP Model

**Primary Hyperoxaluria (PH1 / PH2 / PH3) — Quantitative Systems Pharmacology Model**
*Enteric (secondary) hyperoxaluria included as a different source term in the same equations*

| Deliverable | File | Scale |
|---|---|---|
| 🗺️ Mechanistic map | [`ph1_qsp_model.dot`](../../../primary-hyperoxaluria/ph1_qsp_model.dot) · [SVG](../../../primary-hyperoxaluria/ph1_qsp_model.svg) · [PNG](../../../primary-hyperoxaluria/ph1_qsp_model.png) | 238 nodes · 21 clusters · 292 edges |
| ⚙️ mrgsolve ODE model | [`ph1_mrgsolve_model.R`](../../../primary-hyperoxaluria/ph1_mrgsolve_model.R) | **73 ODEs** · 138 parameters · 34 scenarios · 16 diagnostics |
| 📊 Shiny dashboard | [`ph1_shiny_app.R`](../../../primary-hyperoxaluria/ph1_shiny_app.R) | 11 tabs |
| 📚 References | [`ph1_references.md`](../../../primary-hyperoxaluria/ph1_references.md) | **294 entries** (every PMID actually looked up and verified through the NCBI E-utilities) |

Verification environment: **mrgsolve 2.0.1 / R 4.3.3**. All 34 scenarios and all 16 diagnostics were actually run and passed.

---

## 1. The organising idea

> **Oxalate is a terminal metabolite. Humans have no oxalate-degrading enzyme.
> Primary hyperoxaluria is therefore not "a disease with a severity score" but
> **an unclosed mass balance**, and every clinical finding is arithmetic about
> where the oxalate that could not be excreted has gone.**

```
d(total body oxalate)/dt = generation − renal excretion − intestinal excretion − dialytic removal
```

No parameter named `severity`, `stage`, or `progression_rate` **exists** in this
model. Urinary oxalate, plasma oxalate, the eGFR slope, the age at which ESKD is
reached, the stone rate and systemic oxalosis are all **outputs** of the balance
above. There is exactly one input separating the lethal infantile form from the
adult stone-former found by accident — `FAGT`, the residual
alanine:glyoxylate aminotransferase activity.

Rearranged, this single equation produces the one ratio that governs the whole
disease.

```
        generation (set by genotype)
        ────────────────────────────
        renal clearance (set by GFR)
```

And the two terms of this ratio move in **opposite directions** over a lifetime.
Because the oxalate the kidney excretes is what destroys that same kidney. The
denominator is eaten by the numerator.

---

## 2. Eight things this model is designed to **generate** (rather than assume)

### (1) AGT non-linearity — recessive inheritance and chaperone rescue on one curve

AGT is a high-Vmax, low-Km enzyme, so the mapping from residual activity to
urinary oxalate is sharply non-linear. **A single Michaelis-Menten term** yields
together two facts that a linear model cannot satisfy at once (diagnostic `D03`).

| Residual AGT activity | Urinary oxalate (mmol/1.73m²/day) | vs normal | vs ULN (0.46) |
|---|---|---|---|
| 100 % | 0.270 | 1.00× | 0.59× |
| **50 %** (carrier) | **0.289** | **1.07×** | 0.63× |
| 10 % | 0.478 | 1.77× | 1.04× |
| **3 %** (after B6 rescue) | **0.971** | 3.60× | 2.11× |
| **0.5 %** (G170R) | **1.464** | 5.43× | 3.18× |
| 0 % | 1.586 | 5.88× | 3.45× |

50 % of the enzyme is indistinguishable from 100 % (→ why PH1 is recessive and
carriers are entirely healthy), and merely raising 0.5 % to 3 % — a liver that is
still severely deficient — makes **33.7 %** of the oxalate disappear (→ why
pharmacological chaperone rescue is worth pursuing).

### (2) A threshold disease out of a linear sink

While GFR is preserved the kidney removes almost the whole of what is generated,
so plasma oxalate is near normal and the disease is a stone disease. Renal
clearance is proportional to GFR, so once GFR falls below
`generation/Pox_crit` the remainder has nowhere to go but tissue. The model
**computes** the eGFR at which plasma oxalate exceeds the plasma CaOx solubility
limit (≈30 µmol/L) (diagnostic `D04`).

> **Computed result: eGFR = 27.8 mL/min/1.73 m²**

Systemic oxalosis is not a separate mechanism that switches on late but **the
root of this equation.** The literature's observation that "oxalosis appears
below an eGFR of 30–45" is reproduced with no parameter that mentions oxalosis at
all.

### (3) A runaway loop whose gain grows

Crystals destroy nephrons → GFR falls → clearance falls → plasma oxalate rises →
deposition increases → nephrons destroyed again. Nowhere in the model is the
eGFR trajectory prescribed. Diagnostic `D14` shows that the gain of this loop
increases monotonically as GFR falls.

| Age | eGFR | Pox | dGFR/dt | Loop gain |
|---|---|---|---|---|
| 5 y | 107.5 | 8.2 | −5.50 | 0.0512 |
| 15 y | 53.4 | 16.5 | −4.46 | 0.0836 |
| 20 y | 34.1 | 24.7 | −3.63 | 0.1065 |
| 24 y | 18.9 | 40.0 | −3.87 | 0.2049 |
| 25 y | 15.1 | 47.3 | −3.76 | 0.2490 |

No new mechanism switches on. It is only the same loop gaining more gain as the
sink gets smaller, and this is why end-stage deterioration is abrupt rather than
linear. **The argument for early diagnosis is not sentiment but loop-gain
arithmetic** (diagnostic `D16`).

### (4) The biomarker that inverts — why one drug needed two endpoints

`Uox = Cl_ox × Pox`. Below the deposition threshold Uox *is* the amount
generated and a clean efficacy biomarker. Above the threshold, the falling GFR
**drags Uox down** while the disease accelerates. Diagnostic `D05` holds
generation exactly fixed and changes only the number of nephrons.

| eGFR | Uox | vs eGFR 125 |
|---|---|---|
| 125.0 | 1.465 | 100.0 % |
| 45.3 | 1.450 | 99.0 % |
| 22.9 | 1.297 | 88.5 % |
| **13.9** | — | **75 %** (computed crossing point) |
| **5.8** | — | **50 %** (computed crossing point) |

Between an eGFR of 125 and 40, Uox barely moves. That is, **a fall in Uox in
advanced PH is not therapeutic success.** Which is why ILLUMINATE-A (eGFR ≥ 30)
could use 24-hour Uox and ILLUMINATE-C (advanced CKD/dialysis) had to switch to
plasma oxalate. **The change of endpoint is arithmetic, not regulatory taste.**

### (5) 🏆 The headline — why LDHA siRNA is weaker than HAO1 siRNA and fails in PH2

LDHA is the last step in PH1, PH2 and PH3 alike. Read naively, LDHA siRNA ought
to treat all three. It does not — nedosiran works in PH1 and failed in PH2. The
model generates this from two structural facts, **with no PH2-specific drug
parameter**.

**(a) Glycolate oxidase does not merely *make* glyoxylate, it also *oxidises*
it to oxalate.** Silencing HAO1 cuts off the substrate supply and one of the two
oxalate-generating reactions at the same time. Silencing LDHA cuts off only the
other one and leaves the GO route as an untouchable floor.

**(b) Silencing an enzyme reduces *flux* only when there is a parallel branch to
take the substrate.** In PH1 that branch is GRHPR, so LDHA knockdown
redistributes glyoxylate into harmless glycolate. **In PH2, GRHPR is precisely
the gene that is missing.** LDHA becomes the **only exit** from the cytosolic
glyoxylate pool, and at steady state the flux through a single exit equals its
input regardless of how much enzyme there is. Knockdown merely raises glyoxylate
and hardly moves oxalate at all.

Diagnostic `D07` quantifies the share of each explanation.

| Condition | Nedosiran Uox change |
|---|---|
| PH1 | **−30.3 %** |
| PH2, extrahepatic source **switched off** | **−9.4 %** |
| PH2 (as modelled, extrahepatic source present) | **−3.3 %** |

- Loss explained by flux partitioning (loss of the parallel GRHPR branch): **20.9 points**
- Loss explained by hepatocyte-restricted delivery: **6.1 points**
- → **Flux partitioning explains 77 % of the PH1→PH2 collapse.**

That is, an LDHA-targeted siRNA is bound to work structurally less well in PH2
**before** hepatocyte-restricted delivery is even considered.

> **Falsifiable prediction.** In PH2, nedosiran should **raise plasma and urinary
> glyoxylate substantially** while oxalate does not change (model: 0.04 → 0.22
> µmol/L). The standard explanation (hepatocyte-restricted delivery) predicts no
> rise in glyoxylate. **A single measurement separates the two explanations.**

The same logic in reverse: a small-molecule LDH inhibitor (stiripentol)
distributes into both liver **and kidney**, so it ought to be superior to
nedosiran in PH2 — model scenario S18 (−3.3 %) versus S19 (**−52 %**).

### (6) A therapy that changes not one flux (Layer B)

Hyperhydration and potassium citrate halve the supersaturation without touching
any term of the mass balance. Because the stone endpoint is **a ratio, not a
flux** (diagnostic `D09`).

| Measure | Untreated | 3 L fluid + citrate | Change |
|---|---|---|---|
| Total generation (µmol/day) | 1473.2 | 1473.2 | **−0.0 %** |
| 24-hour Uox (mmol/day) | 1.463 | 1.464 | **+0.1 %** |
| Urinary oxalate concentration (mmol/L) | 0.975 | 0.488 | −49.9 % |
| AP(CaOx) index | 3.926 | 1.587 | −59.6 % |
| Symptomatic stones (events/year) | 1.807 | 0.364 | −79.9 % |

Generation and 24-hour excretion are effectively unchanged while supersaturation
and the stone rate collapse. This is the clearest possible separation between
**protecting the kidney** and **closing the balance**, and it is why supportive
therapy alone leaves systemic oxalosis behind.

### (7) Dialysis as arithmetic (diagnostic `D08`)

The amount removed is computed by **integrating** the real intradialytic
clearance acting on a plasma concentration that genuinely falls and then rebounds
— it is not assumed. Units: mmol/**week**.

| Prescription | Generation | Dialysis | Urine | Gut | **Shortfall** | Pox |
|---|---|---|---|---|---|---|
| Conventional HD 3×4 h | 10.3 | 3.4 | 4.5 | 0.5 | **+1.9** | 75.3 |
| Daily extended HD 7×6 h | 10.3 | 6.7 | 2.9 | 0.3 | **+0.4** | 45.4 |
| Peritoneal dialysis | 10.3 | 3.8 | 4.3 | 0.5 | **+1.7** | 60.3 |
| HD 3×4 h **+ lumasiran** | 3.6 | 1.2 | 2.3 | 0.2 | **−0.1** | **27.2** |

Conventional dialysis removes only about a third of the weekly generation, so
several mmol accumulate in tissue every week and Pox stays far above the
solubility limit. Intensified dialysis roughly doubles the removal and still
fails to close the balance. **Only a Layer A drug that cuts generation** brings
Pox below the solubility limit and brings the weekly shortfall to within 0.1 mmol
of zero. This is why dialysis is a bridge and not a destination in PH1, and it is
how an inadequate dialysis prescription is turned into an adequate one **without
changing the prescription at all**.

### (8) Not relapse but a reservoir — oxalate release after transplantation (diagnostic `D10`)

Bone is modelled as two pools: a **surface pool** that exchanges over months
(the one released after transplantation) and a **deep crystal-bound pool** that
empties over years. This is the standard structure of calcium bone kinetics.

| Months after transplantation | Uox (CLKT) | Uox (high-dose steroid) | Bone burden (mmol) |
|---|---|---|---|
| 0 | 1.179 | 1.179 | 188.7 |
| 1 | 0.508 | 0.644 | 183.0 |
| 3 | 0.443 | 0.495 | 176.7 |
| 6 | 0.398 | 0.432 | 171.2 |
| 12 | 0.360 | 0.394 | 164.2 |
| 48 | 0.302 | 0.317 | 133.2 |

**Even with a normal liver installed and generation near zero**, urinary oxalate
stays above the normal range for months. The disease has not returned; **the bone
is emptying**. The model predicts that high-dose steroid therapy prolongs this
release, which supplies a testable reason to prefer steroid minimisation after
transplantation in PH1.

---

## 3. Architecture

### 3.1 Mechanistic map — 21 clusters, 238 nodes

| # | Cluster | # | Cluster |
|---|---|---|---|
| 1 | Genetics · genotype-phenotype | 12 | Positive feedback loop of nephron loss |
| 2 | Substrate supply into the glyoxylate pool | 13 | Gut-kidney oxalate axis (+ enteric hyperoxaluria) |
| 3 | Hepatocyte peroxisome — the AGT reaction | 14 | siRNA PK/PD — GalNAc/ASGPR/RISC |
| 4 | Hepatocyte cytosol — LDHA, the terminal step | 15 | Pyridoxine and chaperone rescue |
| 5 | Hepatocyte mitochondrion — HOGA1 | 16 | Small molecules · gene therapy · experimental approaches |
| 6 | Systemic oxalate distribution and the solubility limit | 17 | Supportive therapy (Layer B) |
| 7 | Bone — the reservoir that buffers and betrays | 18 | Dialysis kinetics |
| 8 | Soft-tissue oxalosis | 19 | Transplantation |
| 9 | Renal oxalate handling | 20 | Clinical endpoints |
| 10 | Crystallisation physical chemistry | 21 | **The measurement layer** (biomarker inversion, assay pitfalls) |
| 11 | Crystal-induced tubulointerstitial injury | | |

### 3.2 Three therapeutic layers — none substitutes for another

| Layer | What it changes | Agents |
|---|---|---|
| **A · Generation** | the mass balance itself | lumasiran, nedosiran, stiripentol, pyridoxine (genotype-gated), liver transplantation, gene therapy |
| **B · Crystallisation** | prevents stones **while changing the balance not at all** | hyperhydration, potassium citrate, magnesium, neutral phosphate |
| **C · Removal** | the excretion routes | dialysis, intestinal oxalate degradation (Oxalobacter), oral oxalate decarboxylase |

B protects the kidney but can never prevent systemic oxalosis. C is
quantitatively insufficient on its own (§2-(7)). **Only A closes the balance.**
That ordering is the therapeutic argument of this map.

### 3.3 mrgsolve model — 73 ODEs

| Block | Compartments |
|---|---|
| Hepatic intermediary metabolism (8) | Glycolate, peroxisomal · cytosolic · mitochondrial glyoxylate, HOG, DHG, hydroxypyruvate, L-glycerate |
| Oxalate distribution (9) | Plasma, ECF, bone surface, bone deep, soft tissue, renal parenchyma, plasma glycolate, plasma glyoxylate, intestinal lumen |
| Enzymes · mRNA (9) | HAO1 mRNA, GO protein, LDHA mRNA, hepatic LDHA, **renal LDHA (separate)**, apo/holo-AGT, GRHPR, HOGA1 |
| Drug PK (13) | Lumasiran SC/plasma/liver/RISC, nedosiran SC/plasma/liver/RISC, stiripentol 3 compartments, pyridoxine · PLP |
| Renal structure · injury (10) | Nephrons, crystal burden, NLRP3, IL-1β, TGF-β, fibrosis, tubular injury, acidosis, urine volume, transplanted kidney |
| Stone · organ endpoints (9) | Stone mass, cumulative risk, retina, myocardium, nerve, skin, bone disease, marrow, B6 neurotoxicity |
| Audit integrals (15) | Cumulative generation · urinary excretion · intestinal excretion · dialytic removal · deposition, Pox AUC, AP AUC, time above threshold, etc. |

**Algebraic baseline solve.** Every rate constant is **derived** in `$MAIN` from
the annotated flux targets, so that at `FAGT=1` the whole system is an exact
fixed point. Diagnostic `D01`:

```
Uox   1 y = 0.269697  ->  40 y = 0.269697   drift = +0.0000 %
Pox   1 y = 1.302884  ->  40 y = 1.302884   drift = +0.0000 %
eGFR  1 y = 125.0000  ->  40 y = 125.0000   drift = +0.0000 %
tissue burden after 40 years = 0.0000 mmol
```

The mass balance closure error (diagnostic `D02`) is **below |0.001 %|** in all
three states. The disease is generated from the genotype by a silent
natural-history pre-run, and is never imposed as an initial condition.

---

## 4. Clinical anchor comparison (diagnostic `D06`)

| Anchor | Observed | Model |
|---|---|---|
| Normal Uox (mmol/1.73m²/day) | 0.15–0.46 | **0.27** |
| Uox in untreated PH1 | 1.0–2.5 | **1.46** |
| PH1 Pox, eGFR preserved (µmol/L) | 5–15 | **8.2** |
| PH1 Pox, ESKD (µmol/L) | 60–120 | **75** |
| **ILLUMINATE-A lumasiran Uox change** | **−65.4 %** | **−65.4 %** |
| PHYOX2 nedosiran (PH1) Uox change | ~−39 % | −30.3 % ⚠️ |
| PHYOX2 nedosiran (PH2) Uox change | no response | **−3.3 %** |
| Pyridoxine, G170R homozygous | −30 ~ −50 % | **−34.6 %** |
| Pyridoxine, null allele | 0 % | **+0.0 %** (structural null) |
| Time to ESKD in untreated PH1 | ~25 years | **25.1 years** |
| Rise in urinary glycolate on lumasiran | 2–5-fold | **3.7-fold** |
| Effect of hyperhydration on AP(CaOx) | ~−50 % | **−51 %** |
| Oxalobacter, enteric hyperoxaluria | effective | **−54 %** |
| Oxalobacter, PH1 (ePHex primary endpoint) | failed | **−6 %** |

The last two rows are **the same equations** — only the source term differs. One
intervention gives −54 % in one disease and −6 % in the other. **One equation,
two diseases.**

> ⚠️ The −65.4 % agreement with ILLUMINATE-A is striking, but the substrate
> source split (GO route 75 % / hydroxyproline 15 % / DAO 10 %) was chosen from
> the tracer literature and the result compared afterwards, so it should be read
> as **"a landing" rather than an independent prediction**. The sensitivity to
> this split is in diagnostic `D12`.

---

## 5. Negative and self-refuting results — reported rather than removed

The value of a QSP model lies less in what it got right than in **how it handles
what it got wrong**. Seven failures of this model are left in place, together
with their arithmetic.

**① It cannot generate PH3 (diagnostic `D15`).** Lose the aldolase that
*generates* glyoxylate and oxalate cannot rise, and the model says so. Implement
only the literature's main hypothesis (accumulated HOG inhibiting GRHPR) and PH3
urinary oxalate comes out normal (0.265) or below. That holds even with HOG
accumulated 100-fold. **The model fails to reproduce PH3, and it fails to
reproduce that hypothesis as it stands.** Only by adding a further route in which
the accumulated HOG escapes the mitochondrion and is cleaved **in the cytosol**,
bypassing peroxisomal AGT completely, does the PH3 range (0.429) appear. That is
a specific and testable claim about "where the aldolase activity is", and it is
reported as an unresolved gap.

**② It underpredicts the PHYOX2 nedosiran response (−30.3 % versus ~−39 %).**
Rather than tune it, the discrepancy was **converted into a measurable claim**
(diagnostic `D12`). The single parameter that determines the floor is
`PHI_GOOX` (the fraction of glyoxylate oxidised to oxalate by glycolate oxidase
rather than by LDH).

| `PHI_GOOX` | PH1 baseline Uox | Nedosiran % |
|---|---|---|
| 0.02 | 1.414 | −48.0 % |
| 0.08 | 1.438 | −38.6 % |
| 0.15 (model default) | 1.463 | −30.3 % |
| 0.20 | 1.481 | −25.6 % |

The value consistent with the reported ~−39 % is **0.077**. That is, the model
back-calculates a trial result into a falsifiable claim about enzymology —
testable in a hepatocyte assay.

**③ It cannot reproduce the time course of true infantile oxalosis.** At
`FAGT=0` ESKD comes at 22.6 years, which is not the infantile phenotype that
reaches ESKD within the first months of life. Even adding a growth (collagen
turnover) multiplier and an immature kidney gives years, not months. The
infantile form needs something this model does not have — probably a combination
of far lower renal reserve and a far higher generation rate per unit body surface
area.

**④ It overpredicts end-stage plasma oxalate.** Left alone with neither dialysis
nor transplantation, the model's Pox rises higher than any observed value.
Because the model has neither death nor mandatory renal replacement therapy in
it. This should be read as the model stating **"this state is not
survivable"**, which is why every long-term scenario pairs ESKD with dialysis or
transplantation.

**⑤ Urinary glycolate in PH1 is lower than in the literature.** Model 0.39–0.40
mmol/day against a frequently reported 0.5–1.5. The direction and the 3.7-fold
rise on lumasiran are right, but the absolute value is low. `GLYCO_SUP` could
have been raised to match, but that would disturb other anchors, so it is left
as it is.

**⑥ Layer B drives the tissue burden to zero — which is too good.** Diagnostic
`D09` says the tissue burden at 25 years is 247.5 mmol untreated versus 0.0 mmol
on Layer B. That is because Layer B preserves GFR well enough that Pox never
crosses the deposition threshold, whereas in reality a substantial fraction of
PH1 patients on optimal supportive therapy do eventually progress. It suggests
that the crystal→nephron-loss cascade depends excessively on supersaturation and
that a supersaturation-independent component of progression is missing from the
model.

**⑦ The anti-IL-1 scenario (S33) has no clinical basis whatsoever.** The model
says it preserves eGFR at 45 years from 34 to 77 while changing urinary oxalate
not at all (1.459). This is no more than a falsifiable prediction of the map;
there is no approved therapy in this class. That the effect comes out larger than
Layer B is most likely an artefact of the cascade gain.

---

## 6. Thirty-four scenarios

| ID | Scenario | ID | Scenario |
|---|---|---|---|
| S01 | Healthy control | S18 | PH2 + nedosiran (structural non-response) |
| S02 | Heterozygous carrier (AGT 50 %) | S19 | PH2 + stiripentol (reaches the kidney) |
| S03 | PH1 typical, untreated | S20 | PH3 — GRHPR inhibition hypothesis only |
| S04 | PH1 complete AGT loss, untreated | S21 | PH3 — cytosolic HOG aldolase hypothesis |
| S05 | PH1 mild (AGT 10 %) | S22 | Enteric hyperoxaluria + Oxalobacter |
| S06 | PH1 + hyperhydration alone | S23 | Enteric hyperoxaluria, untreated |
| S07 | PH1 + potassium citrate alone | S24 | PH1 ESKD, conventional HD 3×4 h |
| S08 | PH1 + full Layer B | S25 | PH1 ESKD, daily extended HD |
| S09 | PH1 G170R + pyridoxine 10 mg/kg | S26 | PH1 ESKD, peritoneal dialysis |
| S10 | PH1 null + pyridoxine (futile) | S27 | PH1 ESKD, HD + lumasiran (A+C) |
| S11 | PH1 + lumasiran (from age 5) | S28 | Combined liver-kidney transplantation (CLKT) |
| S12 | PH1 + lumasiran + Layer B | S29 | Kidney transplantation alone (liver retained) |
| S13 | PH1 + nedosiran | S30 | Kidney transplantation alone + lumasiran |
| S14 | PH1 + lumasiran + nedosiran | S31 | CLKT, high-dose steroid |
| S15 | PH1 + stiripentol | S32 | PH1 + lumasiran (started at age 20) |
| S16 | PH1 + Oxalobacter colonisation | S33 | PH1 + anti-IL-1/NLRP3 blockade |
| S17 | PH2 (GRHPR null), untreated | S34 | PH1 + high-dose ascorbate |

Selected key results (at age 20; for S24–S31, at the end of follow-up):

| ID | Uox | Pox | eGFR | AP | Stones/year | ESKD age |
|---|---|---|---|---|---|---|
| S01 healthy | 0.270 | 1.30 | 125.0 | 0.72 | 0.00 | — |
| S03 PH1 untreated | 1.390 | 24.65 | 34.1 | 4.25 | 2.04 | **25.1** |
| S08 Layer B | 1.463 | 8.40 | 105.2 | 1.59 | 0.37 | — |
| S11 lumasiran | 0.505 | 3.49 | 87.3 | 1.37 | 0.23 | — |
| S12 lumasiran + Layer B | 0.507 | 2.55 | **120.1** | **0.55** | **0.00** | — |
| S13 nedosiran | 1.014 | 11.18 | 54.8 | 2.90 | 1.18 | 38.6 |
| S14 lumasiran + nedosiran | **0.296** | 1.84 | 97.0 | 0.80 | 0.00 | — |
| S18 PH2 + nedosiran | 0.881 | 6.99 | 76.1 | 2.42 | 0.89 | — |
| S19 PH2 + stiripentol | **0.440** | 2.47 | 107.5 | 1.18 | 0.11 | — |
| S29 kidney transplant alone | 1.483 | 15.32 | 58.4 | 4.21 | 1.96 | — |
| S30 kidney transplant alone + lumasiran | 0.533 | 4.43 | **72.5** | 1.45 | 0.36 | — |

**Combination therapy normalises.** S14 (lumasiran + nedosiran) gives a Uox of
0.296 — within the normal range. Because the GO route and the LDHA route are the
**two** reactions generating oxalate and each drug blocks only one of them,
blocking both together brings the pharmacological floor below normal. A testable
prediction.

**A paradigm shift the model predicts.** S29 versus S30 — close the balance with
drugs and kidney transplantation alone becomes reasonable, sparing the liver
(eGFR at 45 years 58.4 → 72.5, Pox 15.3 → 4.4).

**The price of early diagnosis (diagnostic `D16`).**

| Age at which lumasiran is started | eGFR at 45 years | ESKD age |
|---|---|---|
| 2 y | **84.8** | — |
| 5 y | 69.3 | — |
| 10 y | 47.1 | — |
| 15 y | 30.6 | — |
| 20 y | 17.9 | — |
| 25 y | 0.2 | 25.0 |
| Untreated | 0.2 | 25.1 |

The same drug buys an entirely different kidney at 2 years and at 20. There is no
separate early-treatment parameter; it is loop-gain arithmetic. **Nephrons
already lost do not come back when the balance is closed later.**

---

## 7. Sixteen diagnostics

| ID | Question |
|---|---|
| D01 | Is the healthy baseline an exact fixed point? (drift 0.0000 %) |
| D02 | Does the mass balance close? (error < 0.001 %) |
| D03 | AGT non-linearity — recessive inheritance and chaperone rescue |
| D04 | At what eGFR does Pox exceed 30 µmol/L? (**27.8**) |
| D05 | The eGFR at which Uox falls 25 %/50 % with generation held fixed (**13.9 / 5.8**) |
| D06 | Comparison table of 14 clinical anchors |
| D07 | **Decomposition of the PH2 collapse — flux partitioning 77 % versus delivery 23 %** |
| D08 | Dialysis arithmetic — the computed weekly shortfall |
| D09 | Proof that Layer B changes no flux whatever |
| D10 | Release of the bone reservoir after transplantation |
| D11 | siRNA PK/PD disconnect (plasma t½ 5.5 h, effect 662 days → ~2900-fold) |
| D12 | Estimating `PHI_GOOX` by back-calculation from PHYOX2 (0.077) |
| D13 | The pyridoxine genotype gate and the conversion ceiling |
| D14 | How the loop gain changes with GFR |
| D15 | **PH3 — the model and the hypothesis fail together** |
| D16 | Early versus delayed treatment |

---

## 8. How to run

```bash
# 1) render the mechanistic map
dot -Tsvg ph1_qsp_model.dot -o ph1_qsp_model.svg
dot -Tpng -Gdpi=150 ph1_qsp_model.dot -o ph1_qsp_model.png

# 2) model + 34 scenarios + 16 diagnostics (full output)
Rscript ph1_mrgsolve_model.R

# 3) Shiny dashboard (11 tabs)
Rscript -e 'shiny::runApp("ph1_shiny_app.R", port=7788, launch.browser=FALSE)'
```

Required packages: `mrgsolve` (≥2.0.1), `dplyr`, `tidyr`, `ggplot2`, `shiny`,
`DT`. Run under a UTF-8 locale (`LANG=C.UTF-8`).

### The eleven tabs of the Shiny dashboard

1. **Mass balance** — displays the balance sheet itself like a budget, closure error included
2. **Genotype** — the AGT non-linearity curve in real time
3. **Urinary oxalate** — the efficacy biomarker while eGFR is intact + the glycolate "receipt"
4. **Plasma oxalate** — the biomarker once the threshold has been crossed + plasma glyoxylate
5. **Biomarker inversion** — the eGFR at which Uox loses its meaning (generation fixed)
6. **Kidney** — the eGFR trajectory and the crystal-injury cascade
7. **Stones** — the AP(CaOx) index and the event rate (Layer B)
8. **Systemic oxalosis** — bone (surface/deep), retina, heart, nerve, skin, marrow
9. **Drug PK/PD** — siRNA plasma versus RISC versus enzyme versus effect, renal LDHA included
10. **Dialysis · transplantation** — the weekly removal arithmetic and post-transplant release
11. **Scenario comparison** — side-by-side comparison of any two arms

> This app has **no severity slider.** Only genotype and treatment are set;
> everything else is an output of the balance.

---

## 9. An honest division of the basis of the parameters

The values grounded directly in the literature and the estimated values are
separated into a table in the section "An honest division of the basis of the
parameters" of
[`ph1_references.md`](../../../primary-hyperoxaluria/ph1_references.md). In
summary:

- **Literature-based**: the normal/PH1 Uox ranges, the normal/ESKD Pox ranges,
  the plasma CaOx solubility limit, the renal clearance ratio for oxalate, the
  ILLUMINATE-A/PHYOX2 results, the dosing regimens and plasma t½ of lumasiran and
  nedosiran, the G170R pyridoxine response, the hyperhydration targets, the
  Tiselius AP(CaOx) formula, the proportion reaching ESKD by age 25, the
  hydroxyproline contribution, and ASGPR hepatocyte selectivity.
- **Estimates (labelled as such)**: `PHI_GOOX`, the cytosolic glyoxylate
  branching ratios, `GXP_REF`/`KM_AGT` (AGT saturation — the largest
  sensitivity), the tissue deposition clearance and threshold, `BMAX_BONE`, the
  two-compartment bone kinetic constants, every individual rate constant of the
  crystal→fibrosis cascade (only the composite outcome was calibrated),
  `OX_XHEP2`, and the `IL1_BLOCK` effect size.

The headline result of `D07` holds even with `OX_XHEP2` set to zero, so it does
not depend on that estimate.

---

## 10. Disclaimer

This model is a **quantitative QSP model for educational and research purposes**.
It was assembled from the open literature and clinical trial data but has not
been independently verified or certified, and **must not be used directly for
actual clinical decision-making, prescribing, or regulatory submission.** In
particular, do not quote its figures without reading the seven failures in §5 and
the list of estimated parameters in §9.
