# Short Bowel Syndrome · Chronic Intestinal Failure (SBS-IF) QSP Model
## Short Bowel Syndrome with Chronic Intestinal Failure — Quantitative Systems Pharmacology Model

<p align="center">
  <a href="sbs_qsp_model_en.svg">
    <img src="sbs_qsp_model_en.png" width="880" alt="SBS-IF QSP mechanistic map">
  </a>
  <br>
  <em>Click to open the zoomable SVG · 23 clusters · 276 nodes · 381 edges</em>
</p>

---

## The organising idea of this model

Most other models in this library are built around a single **disease-severity
state** that a drug pushes down while a clinical score follows it. SBS-IF is
structurally different, and that difference is the whole of this model.

> ### Parenteral nutrition (PN) volume is not an input but an **output**.
>
> It is the **arithmetic residual** left over by two conservation equations.
>
> ```
>   water   dTBW/dt = oral intake + PN − output − insensible loss − urine
>   sodium  dNa /dt = oral Na    + PN Na − output Na − urine Na
>   energy  dE  /dt = absorbed kcal + colonic SCFA salvage + PN kcal − total energy expenditure
> ```
>
> Every arrow on the map justifies its place by **which term of these three
> equations it changes**. Villus height is not the endpoint — **litres** are the
> endpoint.

A patient whose 70 cm of jejunum ends in a stoma is not ill in the sense in which
a lupus patient is "ill". Their mucosa is histologically near-normal. What has
collapsed is the **budget**.

---

## Six structural commitments and what they generate

### 1. Intraluminal fluid is not one stream but **two**

Fluid that arrives **with** a meal (stream M) and fluid drunk **between** meals
(stream D) differ completely, in sodium concentration and in glucose content
alike. The driving force `[Na⁺]lum − CEQ` is computed for the two streams
**separately**, and only then averaged by volume fraction. This ordering
**generates** the most counter-intuitive fact about this disease.

| Composition of the 2.5 L drunk between meals | drink [Na⁺] | stream D driving force | direction | output |
|---|---|---|---|---|
| plain water only (ORS 0) | 5 mmol/L | **−95.0** | **SECRETES** | 4.15 L/day |
| ORS 0.2 | 22 | −66.3 | secretes | 3.65 L/day |
| ORS 0.4 (default) | 39 | −43.9 | secretes | 3.30 L/day |
| ORS 0.7 | 65 | −14.3 | secretes | 2.89 L/day |
| ORS 0.9 | 82 | **+4.3** | **absorbs** | 2.67 L/day |

The driving force of the same stream M **does not change** through all five rows:
+44.4 throughout. The only place the sign turns over is stream D. A single-stream
model that averages first has to **declare** this rule, whereas this model
**produces** it as a sign change in the flux equation.

### 2. Glucose does not add an absorption term, it **shifts the zero point**

Net sodium movement in the jejunum is zero at an intraluminal
`[Na⁺] ≈ 90-100 mmol/L` and becomes negative below it (Fordtran's perfusion
studies). SGLT1 cotransport therefore enters as **a term that drags CEQ
downwards** (~60 mmol/L at saturation). Put in as an added absorption term, it
double-counts the sodium that NHE3 immediately recycles and produces absurd
fluxes. Put in as a shift of the zero point, it reproduces the ORS effect and its
**saturability** at the same time.

### 3. Each nutrient has **a different reference bowel length**

Carbohydrate and protein are largely absorbed in the first 100-150 cm, fat needs
200-300 cm, and bile acids and B12 **will take nothing but the terminal ileum.**
A single "absorptive surface area" scalar cannot reproduce the ranking observed in
the same patient (CHO 70-90%, protein 60-80%, fat 30-50%). Three reference lengths
can, and in addition they place each micronutrient deficiency **on a named segment
rather than on a generic severity**.

### 4. Adaptation is not a sum but a **product**

```
TROPHIC = intraluminal nutrient × (1 + GLP2R occupancy term) × (1 + IGF-1 term) × Zn sufficiency
```

If enteral intake is zero the product is zero, and a completely fasted patient
should neither adapt nor respond to a GLP-2 analogue. **Enteral nutrition is not
an adjunct but a multiplying factor.** (This design intention is, however, partly
refuted by the model itself — see D07 and D17 below.)

### 5. The prescription is **a closed-loop controller running the actual trial protocol**

In STEPS it was not the investigator but an algorithm that decided PN reduction:
**if 24-hour urine output rises more than 10% above baseline, reduce PN; if weight
or electrolytes fall, hold.** This algorithm is implemented as a rate-limited
bidirectional controller on the `PNVOL` state, and three things follow.

- The registered endpoint (**≥20% reduction in PN volume**) comes out not as the
  result of fitting a PN-weaning parameter to the trial, but as **the result of
  simulating the protocol**.
- The controller **stops itself** (reduce PN and urine output returns to baseline,
  so weaning halts) — the reason real PN weaning reaches a plateau.
- And **a placebo-arm response emerges.** The urine-output trigger weans anyone
  with absorptive reserve left, drug or no drug.

### 6. The ileal brake **positive feedback loop is closed**

L cells sit in the terminal ileum and the proximal colon — exactly the tissue that
gets resected. Lose them and the brake disappears, transit accelerates, contact
time shortens, absorption falls, and the unabsorbed load so created **passes L
cells that are already gone.** Because contact time is a **multiplying factor** on
every absorptive fraction, an opioid that does nothing but slow transit time has a
genuinely quantitative effect in this disease.

---

## The four deliverables

| Deliverable | File | Scale |
|---|---|---|
| 🗺️ Mechanistic map | [`sbs_qsp_model_en.dot`](sbs_qsp_model_en.dot) · [`.svg`](sbs_qsp_model_en.svg) · [`.png`](sbs_qsp_model_en.png) | **23 clusters · 276 nodes · 381 edges** |
| ⚙️ mrgsolve ODE model | [`sbs_mrgsolve_model_en.R`](sbs_mrgsolve_model_en.R) | **71 ODEs · 28 scenarios · 17 diagnostics** |
| 📊 Shiny dashboard | [`sbs_shiny_app_en.R`](sbs_shiny_app_en.R) | **10 tabs** |
| 📚 References | [`sbs_references_en.md`](sbs_references_en.md) | **236 papers · every PMID verified** |

### How to run

```bash
# all scenarios + diagnostics (verified on mrgsolve 2.0.1, R ≥ 4.3)
Rscript sbs_mrgsolve_model_en.R

# interactive dashboard
Rscript -e 'shiny::runApp("sbs_shiny_app_en.R", port = 8080)'

# re-render the map
dot -Tsvg sbs_qsp_model_en.dot -o sbs_qsp_model_en.svg
dot -Tpng -Gdpi=150 sbs_qsp_model_en.dot -o sbs_qsp_model_en.png
```

---

## Model structure

### The 71 ODE compartments

| Layer | Compartments | Contents |
|---|---|---|
| Drug PK/PD | 18 | teduglutide (SC depot + central) · apraglutide · glepaglutide · native GLP-2 · GLP-2 effect site · somatropin · IGF-1 · loperamide · octreotide · proton pump pool · colestyramine (luminal) · rifaximin (luminal) · glutamine (luminal) · anti-drug antibody |
| Mucosa · adaptation · motility | 12 | villus · crypt · slow structural remodelling · SGLT1 · NHE3 · brush-border enzymes · colonic fermentative capacity · L-cell mass · PYY/GLP-1 brake tone · gastric acid hypersecretion · contact time · citrulline |
| Luminal ecology · bile acids | 5 | bile-acid pool · SIBO · D-lactate · portal endotoxin · bacterial mucosal injury |
| Conservation · body composition · prescription | 11 | total body water · exchangeable sodium · lean mass · fat mass · hyperphagia multiplier · **PN volume** · **PN kcal** · smoothed urine output · cumulative energy balance · bicarbonate · thirst behaviour |
| Micronutrients · bone | 10 | serum Mg · Zn · B12 · vitamins D/A/E · Se · essential fatty acids · PTH · bone mineral density |
| Liver (IFALD) | 5 | hepatic phytosterol · choline deficiency · hepatic inflammation · bilirubin · fibrosis |
| Kidney · stones | 4 | urinary oxalate · stone burden · eGFR · acute kidney injury |
| Catheter · outcomes | 6 | cumulative catheter-days · cumulative CRBSI · usable veins remaining · colonic polyps · quality of life · cumulative mortality risk |

### The 23 clusters of the mechanistic map

Aetiology · anatomy (parameters) · absorption coefficients (Damköhler product) ·
secretory load · **jejunal Na⁺/water flux (where the sign turns over)** · **the
colon = a digestive organ** · **the L-cell axis (the broken loop)** · motility and
transit time · structural adaptation · functional adaptation · **bile acids and the
100 cm rule** · SIBO and D-lactic acidosis · **the energy conservation equation** ·
**the water and sodium conservation equations** · **PN prescription = a closed-loop
controller** · GLP-2 analogue pharmacology · adjunctive pharmacotherapy ·
**IFALD** · **central venous access = a consumable resource** ·
kidney · stones · bone · micronutrients · clinical endpoints · legend

The map marks **structural null / trap** nodes separately (double octagons). At
those points the effect must be zero or **of the opposite sign**, and the
diagnostics test whether the model actually behaves that way: the plain-water trap ·
ORS correction · zero salvage in an end jejunostomy · complete fasting · native
GLP-2 · the reversal of colestyramine · the price of octreotide · D-lactic acidosis
requires a colon · Mg-deficient refractory hypocalcaemia · fish-oil emulsion · the
plateau of the placebo response · the anatomical floor.

---

## Calibration anchors and what the model produces

The default patient has **80 cm of jejunum + half a colon, no ileocaecal valve, no
ileum** (Messing type 2); a 900-day burn-in brings it to a self-consistent steady
state, and every scenario starts from there.

### Baseline phenotype (D02)

| Metric | Model | Published SBS-IF |
|---|---|---|
| PN volume | **13.1 L/week** | STEPS baseline ~12.9 L/week |
| Stool/stoma output | 3.3 L/day | 2-4 L/day |
| Effluent [Na⁺] | 84 mmol/L | 90-100 (the model is about 10 low — D17f) |
| Plasma citrulline | 16.1 µmol/L | SBS-IF 10-20 (normal 30-40) |
| Carbohydrate absorption | 76.6 % | 70-90 % |
| Protein absorption | 67.1 % | 60-80 % |
| Fat absorption | 32.2 % | 30-50 % |
| Hyperphagia multiplier | 1.33 × | up to ~2 × |
| Colonic SCFA salvage | 272 kcal/day | up to ~1000 (in the worst malabsorption) |

### Reproducing STEPS (D03) — the endpoint the protocol generated

| | Model | Observed (STEPS) |
|---|---|---|
| Teduglutide, change in PN at 24 weeks | **−4.45 L/week** | −4.4 L/week |
| Placebo (no protocol), change in PN | +0.02 L/week | — (confirms baseline stability) |
| Placebo (with protocol, D10) | **−23.9 %** | 30 % responders |
| Teduglutide (with protocol, D10) | **−53.8 %** | 63 % responders |
| Change in citrulline | +2.3 µmol/L | significant rise |

### The difference anatomy makes (D06)

Jejunum is held at 80 cm and only the colon changes.

| Colon in continuity | PN volume | Output | SCFA salvage | Fat absorption |
|---|---|---|---|---|
| 0 (end jejunostomy) | **26.2 L/week** | 5.26 L/day | 0 kcal/day | 28.7 % |
| 0.5 | 13.1 L/week | 3.30 L/day | 272 kcal/day | 32.2 % |
| 1.0 | **3.1 L/week** | 1.82 L/day | 181 kcal/day | 39.5 % |

A colon on its own changes the PN requirement by more than eightfold. And **the
colon's energy salvage is largest where the malabsorption is worst** — because the
substrate is precisely the unabsorbed carbohydrate (556 kcal/day with 40 cm of
jejunum + whole colon + a high-carbohydrate diet).

### The anatomical floor (D13) — the line no drug crosses

End jejunostomy, two years of teduglutide:

| Residual small bowel | PN at 2 years | Enteral autonomy |
|---|---|---|
| 40 cm | 28.0 L/week | no |
| 90 cm | 17.6 L/week | no |
| 130 cm | 10.1 L/week | no |
| 180 cm | 2.9 L/week | no |
| 250 cm | 0.2 L/week | **yes** |

The drug effect is roughly **a constant absolute gain**, so it cannot convert a
very short bowel into autonomy. The length at which autonomy is reachable is set
not by dose but by **anatomy**.

### The Ala2→Gly substitution is the drug (D11)

| Regimen | Mean GLP2R occupancy | Change in PN at 24 weeks |
|---|---|---|
| Teduglutide 0.05 mg/kg once daily | 0.439 | **−4.45 L/week** |
| Native GLP-2, equimolar, once daily | 0.058 | −1.47 L/week |
| Native GLP-2, equimolar, twice daily | 0.102 | −1.54 L/week |
| Native GLP-2 at **8× the dose**, twice daily | 0.473 | −4.62 L/week |

Because DPP-4 cleaves at Ala2-Gly3 and t½ is about 7 minutes, the same number of
moles given once a day does not sustain occupancy. **A supraphysiological dose
twice a day** reaches a comparable effect, which matches the regimen the human
studies of native GLP-2 actually used.

### Long-acting analogues: not mean exposure but **the shape of the profile** (D12)

| Regimen | Mean occupancy | Trough/peak | Change in PN at 24 weeks |
|---|---|---|---|
| Teduglutide once daily | 0.439 | 1.00 | −4.45 L/week |
| Apraglutide once weekly | 0.904 | 0.86 | −6.05 L/week |
| Glepaglutide once weekly | 0.551 | **0.07** | −4.30 L/week |
| Glepaglutide twice weekly | 0.855 | 0.81 | −5.85 L/week |

That the trough/peak ratio of once-weekly glepaglutide is 0.07 — that is, an
extremely peaked profile — shows exactly the structural reason a **twice-weekly**
regimen was the one used in the clinical trials. (That apraglutide shows the larger
effect is a **prediction** of the model and not a validated result.)

### IFALD: the composition of the lipid emulsion is the dose (D15)

| At 2 years | Hepatic phytosterol | Bilirubin | Fibrosis stage |
|---|---|---|---|
| Soybean oil (phytosterol ~350 µg/mL) | 527 | **2.44 mg/dL** | 0.74 |
| SMOF blend | 150 | 1.25 | 0.36 |
| Fish-oil based (phytosterol 0) | 0 | 1.00 | 0.26 |
| Soybean oil for 2 years → switched to fish oil for 1 year | — | **2.44 → 1.00** | — |

### Access is a consumable resource (D16)

| Over 2 years | Nights per week | Catheter-days | CRBSI per year | Veins remaining |
|---|---|---|---|---|
| Standard care | 5.28 | 548 | 0.39 | 5.73 |
| Taurolidine lock + dedicated team | 5.28 | 548 | **0.07** | 5.95 |
| Teduglutide | **3.37** | **366** | 0.25 | 5.82 |

**Catheter-days fall only when nights fall, not litres.** In the model each
infusion night delivers at most 2.5 L, so once the prescribed volume drops the team
removes a whole night — and that is what changes infection and access risk. One of
the clinical arguments for the GLP-2 analogues lies here.

---

## Where the model refutes itself — negative and self-refuting results (D17)

Six results contradict either this model's design intention or the literature it
was built on. **They are reported rather than deleted.**

**(a) The colestyramine trap comes out in the opposite direction.** The textbook
rule is that a bile-acid sequestrant helps below 100 cm of ileal resection and is
**harmful** above it. D09 reproduces the pool gradient exactly (2.3 g with 90 cm of
ileum remaining, 0.8 g with none at all). Yet the loss of fat absorption is largest
where the pool is **intact** (−16.5 %) and almost absent where it is depleted
(−0.5 %) — for the unavoidable reason that a depleted pool has nothing left to
sequester. The clinical rule must therefore rest on something this model does not
have: **binding of fat-soluble vitamins and drugs** by the resin, and the fact that
in extensive resection the remaining benefit is zero, so **no cost at all is
tolerable**. The model reproduces "there is no useful benefit" but not "it is
actively harmful", and that gap is a real limitation.

**(b) The individual anchor is right but the population effect is over-predicted.**
D03 lands almost exactly on the STEPS mean (−4.45 against −4.4 L/week) — because
that is the point it was calibrated to — but the virtual population in D04 produces
**a larger drug−placebo difference than the trial**. Matching a mean is not
matching a distribution. The causes of non-response (stopping the drug,
intercurrent illness, cautious weaning not driven by urine output) are absent from
the model, so **the responder fraction must not be quoted.**

**(c) Nutrient gating weakens the drug effect but does not abolish it (D07).** The
multiplicative structure was there to make a completely fasted patient unresponsive
to a GLP-2 analogue. It does not do that — because only the trophic arm is gated
and the motility arm is not (−4.47 at 100 % enteral, −1.57 L/week at 2 %, fasted).
Whether the real drug retains its transit-time effect in a completely fasted bowel
is an **open experimental question** that this model now poses explicitly.

**(d) Octreotide looks better than it is.** Even after halving its antisecretory
effect and adding a penalty for pancreatic lipase inhibition, D08 shows a sustained
PN reduction larger than the clinical literature supports. Tachyphylaxis, gallstone
formation and poor durability are absent from the model, and they are most of the
reason octreotide is a last resort rather than a first-line agent.

**(e) The baseline drift is small but not zero** (D01, about 0.5 %/year in PN
volume). The burn-in converges the physiological states, but the time constant of
the slow structural remodelling state is about 250 days, so a residual trend
survives even a 900-day run-in. **Effects smaller than 0.5 % per year cannot be
resolved with this model.**

**(g) The SIBO pathway is alive but silent at the endpoint.** Losing the ileocaecal
valve raises the bacterial load from 0 to 0.62, and a 14-day course of rifaximin
cuts it to 0.35 with a measurable recovery of the villi. But the load regrows to
baseline by about day 60 and PN volume at 24 weeks does not change. It is a
faithful reproduction of a short antibiotic course, but at the same time it means
that **the model has no mechanism by which SIBO changes the endpoint the trials
measure**. Bacterial deconjugation of bile acids is on the map but not in the
equations, and wiring it in is the obvious next step.

**(f) Effluent sodium comes out low.** About 84 mmol/L in D02, against an observed
90-100. The two-stream construction dilutes the lumen more than a real intermittent
drinker does, and no calibration was done to force this into agreement.

---

## Honest limitations

- It is a **whole-body daily-average model**. There is no meal timing, no circadian
  structure and no spatial gradient along the lumen. The two-stream construction is
  a crude surrogate for meal timing and must not be read as real pharmacokinetic
  compartmentalisation inside the lumen.
- **The fluid model cannot be used for prescribing.** These are population
  parameters calibrated to reproduce plausible outputs and PN volumes for the
  stated anatomy, not a validated absorption model.
- Several coefficients (jejunal `FMAXJ`, colonic `KCOLMAX`, the adaptation time
  constants) are **identifiable only as a set.** Other combinations reproduce the
  same anchors.
- **Survival, transplantation and quality of life are crude risk accumulators**,
  put in so that the **structure** of the competing risks is visible. Their
  absolute values must not be quoted.
- Paediatric SBS (growth, intestinal lengthening, NEC aetiology) is reachable only
  through parameters and has not been separately calibrated.
- Autologous intestinal reconstruction (STEP, Bianchi) and intestinal
  transplantation are on the map but not implemented as ODEs.

---

## Provenance

[`sbs_references_en.md`](sbs_references_en.md) classifies
**236 papers** into 20 sections. Every PMID was looked up through the NCBI
E-utilities and its title, first author, year and journal checked against the actual
record. A PMID written from memory frequently points at a completely unrelated
paper (an accident that actually happened in this repository), so no citation in
this file was written by hand.

The front of the references file carries an **anchor map**, so that each
quantitative claim in the model can be checked against the parameter it sits in and
the section of literature it rests on.

---

## ⚠️ Disclaimer

This is a **qualitative / semi-quantitative QSP model for educational and research
purposes**. It was built from the public literature and clinical-trial data but has
not been independently validated or certified, and **must not be used directly for
real clinical decision-making, prescribing, or regulatory submission.** The
parameters and assumptions are illustrative approximations, and separate fitting
and validation against real patient data are required. The absorption coefficients
and the fluid-balance parameters in particular cannot be applied to an individual
patient.

Read the negative-results section (D17) before quoting any number from it.
