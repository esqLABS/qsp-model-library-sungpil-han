# Elevated Lipoprotein(a) QSP Model
## Elevated Lipoprotein(a) — ASCVD · Calcific Aortic Valve Disease

<p align="center">
  <a href="../../../elevated-lipoprotein-a/lpa_qsp_model.svg"><img src="../../../elevated-lipoprotein-a/lpa_qsp_model.png" width="900" alt="Lp(a) QSP mechanistic map"></a>
</p>

---

## The question this model tries to answer

Lp(a) is a strangely difficult risk factor to work with. Causality has been
established by Mendelian randomisation, the target is safe (people with
LPA loss-of-function are healthy), and drugs that lower it by more than
90% already exist. And yet the field is full of observations that look
contradictory.

- Statins raise LDL receptors, and **Lp(a) rises instead of falling.**
- PCSK9 inhibitors lower LDL-C by 60% but Lp(a) by only 25–30%. **Same
  receptor — why a different answer?**
- Niacin lowers Lp(a) by 20%, yet its clinical trials were flatly negative.
- Muvalaplin is reported as **−85.8% by one assay and −70% by another.**
- Lipid-lowering trials for aortic valve stenosis (SEAS, ASTRONOMER,
  SALTIRE) were **all negative**, while genetics shows a strong causal
  signal.
- The reduction Mendelian randomisation demands and the reduction clinical
  trials demand differ **3-fold.**
- The same patient's Lp(a) can fall below threshold in mg/dL and above
  threshold in nmol/L.

Textbooks list these as separate cautionary notes. This model was built to
show that they are **arithmetic consequences of a single structure.**

---

## Structure — four numbers

> Lp(a) is **one** particle, made at **one** rate-limiting step, causing
> disease through **three** arms of action, and read on **two**
> mutually incompatible scales.

The 50 differential equations are **identical across every scenario**. The
seven contradictions above are not reproduced by changing parameters —
they fall out of this structure by calculation.

---

## (A) One rate-limiting step — production

Lp(a)'s fractional catabolic rate (FCR) is about 0.25 pools/day and
**stays essentially unchanged across a 10-fold range of concentration**
(Rader 1993, 1994). That is, the difference between individuals lies
entirely on the production side. The model builds this in as a
**structure**, not an assumption: `KCAT0` is fixed, and the translation
rate `KTL` is back-calculated from the target baseline concentration.

There is exactly one point in the production chain where isoform size
intervenes: pre-secretory degradation in the endoplasmic reticulum.

```
SECEFF = KSZ³ / (KSZ³ + NKIV2³)
```

| KIV-2 repeat count | 8 | 12 | 22 | 30 | 35 |
|---|---|---|---|---|---|
| Secretion efficiency | 0.954 | 0.860 | 0.500 | 0.283 | 0.199 |

This single term accounts for about a 5-fold range; the rest is
allele-specific transcription.

---

## (B) Only **24%** of Lp(a) catabolism is LDL-receptor dependent

This is the whole story behind the clearance-targeting drugs, and the
model **derives it rather than asserting it.**

```
KCAT0 = KLDLR_LPA + KOTH_LPA = 0.060 + 0.190 = 0.250 /day   (24% LDLR-dependent)
LDL   = KLDLR_LDL + KOTH_LDL = 0.300 + 0.050 = 0.350 /day   (86% LDLR-dependent)
```

When evolocumab raises LDLR from 1.00 to 2.75, **the same receptor change**
gives two completely different answers for the two lipoproteins.

| | LDL clearance | Lp(a) clearance | Result |
|---|---|---|---|
| Baseline | 0.350 | 0.250 | — |
| Evolocumab | 0.875 (×2.50) | 0.355 (×1.42) | **LDL-P −59.0% · Lp(a) −30.1%** |

The observed 25–30% was **never entered into the model.** It falls out of
the 24% receptor biology.

### The statin paradox — two terms of opposite sign

The same arithmetic runs in the opposite direction for statins.

```
clearance term  :  LDLR 1.00 → 1.78  →  Lp(a) −15.7%
transcription term :  sterol response element  →  ×1.30
                                    ────────────
                        net effect   +8.1%     ← observed +8 to +20%
```

The transcription term `ESTA = 0.30` is the only value fitted to this
endpoint, but it is not arbitrary. Given the clearance arithmetic,
**no smaller value could produce a net increase at all.**

And this yields a testable prediction: in patients with no LDL receptor
(HoFH), there is no offsetting clearance gain, so statins should raise
Lp(a) by the full transcription term. **Model value: +25.9%.**

### Ezetimibe — an unfitted prediction that turned out right

A drug that raises LDLR alone, without touching LPA transcription, should
lower Lp(a) by about a quarter of what it does to LDL.

| | model (no post-hoc adjustment) | observed |
|---|---|---|
| LDL-P | −19.3% | — |
| Lp(a) | **−6.6%** | about −7% (Awad 2018 meta-analysis) |

Ezetimibe is not a negative control — it is a **quantitative test** of the
24% figure, and it passes.

This same data simultaneously forces the assembly step to be **saturated
with respect to LDL substrate.** If assembly were first-order in LDL,
ezetimibe should have lowered Lp(a) by nearly 20% more than it did. (This
was in fact a defect in the first-draft model — see §Validation below.)

---

## (C) A second target — assembly

Treating free apo(a) as a separate state variable is the model's second
structural decision. An assembly inhibitor (muvalaplin) does not touch
mRNA, so apo(a) keeps being produced, and the apo(a) that fails to become
a particle **accumulates in plasma.**

```
intact Lp(a)   250 → 39.0 nmol/L   (−84.4%)     ← KRAKEN −85.8%
free apo(a)    2.60 → 10.66 nmol/L (×4.09)
```

So an assay that also captures free apo(a) is **bound to under-report the
effect.** The direction is exactly right. The magnitude is not.

| Traditional apo(a) assay's molar response to free apo(a), `RFREE` | Reported reduction |
|---|---|
| 1.0 (default) | −80.4% |
| 3.0 | −72.5% |
| **3.8** | **−69.4%** ← the −70% reported by KRAKEN |
| 5.0 | −64.9% |

That is, the model says: **for the KRAKEN discrepancy to close, the
traditional assay would need a molar response to free apo(a) of about 3.8
times its response to apo(a) inside a particle.** This is a measurable
quantity, and, as far as could be found, it has not been reported. `RFREE`
is left at **1.0, not 3.8** — rather than tuning the table to look pretty,
the fact that the model under-explains KRAKEN is left visible.

---

## (D) Two scales — and the real problem was not the units

Every pathological cargo of Lp(a) is **exactly 1 per particle.**

| Cargo | Count | Arm responsible |
|---|---|---|
| apoB-100 | 1 | atherogenic |
| OxPL (covalently bound to KIV-10) | about 1 | inflammatory |
| strong lysine-binding site (KIV-10) | 1 | antifibrinolytic |
| cholesterol core | nearly constant | atherogenic |

The only thing that varies with isoform size is the **apo(a) kringle
protein mass, and that is pathologically inert.** So nmol/L is not a matter
of unit preference — it is the **mechanistically correct scale.**

But running the model showed that the usual suspect was not the culprit.

| Isoform 8 → 35 repeats | Range of change |
|---|---|
| **Chemical conversion factor** (nmol/L per mg/dL) | 2.788 → 2.522 · **about 10%** |
| **Antibody bias** (polyclonal mass assay) | 0.562 → 1.406 · **about 150%** |

Unit conversion is almost blameless. **The antibody is the problem.**

### The patients who get missed

Five patients, all with an identical TRUE mass of 50 mg/dL:

| KIV-2 | apo(a) MW | TRUE particle (nmol/L) | Antibody bias | **Reported mass** | Mass ≥50? | Molar ≥125? |
|---|---|---|---|---|---|---|
| 8 | 287 kDa | 139.4 | 0.562 | **28.1** | missed | positive |
| 12 | 343 kDa | 137.2 | 0.688 | **34.4** | missed | positive |
| 22 | 483 kDa | 132.2 | 1.000 | 50.0 | borderline | positive |
| 30 | 595 kDa | 128.4 | 1.250 | 62.5 | positive | positive |
| 35 | 665 kDa | 126.1 | 1.406 | **70.3** | positive | positive |

**The mass assay misses precisely the small-isoform patients — the
genetically highest-risk patients — and flags the large-isoform patients
most strongly.**

A further point: the paired guideline thresholds of 50 mg/dL and 125
nmol/L are not equivalent to each other. A TRUE mass of 50 mg/dL
corresponds to 126–139 nmol/L depending on isoform, so the molar threshold
is always the more inclusive one.

### What's hiding inside 'LDL-C'

For a patient with Lp(a) 250 nmol/L and a small isoform:

```
reported LDL-C   127.3 mg/dL
              = true LDL-C 100.0  +  Lp(a)-cholesterol 27.3
                                     (21% of the reported value is not LDL)

Dahlén correction applied to the (biased) mass assay value → 108.5 mg/dL
                                    +8.5 mg/dL under-corrected vs true LDL-C of 100.0
```

This error runs **in the same direction, in the same patient**, as the
threshold error.

---

## (E) Why Mendelian randomisation and 5-year trials differ 3-fold

Risk is carried by two components with different time constants. That is
the whole explanation.

| Component | Half-life | Contribution to excess risk | What sees it |
|---|---|---|---|
| Slow — atheroma burden | about 14 years | 70% | lifelong genetic exposure (MR) |
| Fast — OxPL/IL-6 vulnerability | about 2 months | 30% | a 5-year trial |

A 62-year simulation (lifelong 250 vs 25 nmol/L):

```
HR at age 65                              2.29   (PAV 44.0 vs 29.8%, VULN 1.448 vs 1.000)
80% reduction from age 60, 5-year integrated risk reduction     20.0%
same 5-year window vs lifelong low Lp(a)                        53.1%
                                    ─────────
                              ratio      2.65-fold
```

This 2.65-fold is **not** a correction factor entered by hand. It is what
falls out of the two time constants.

### So why did niacin fail

It is the **absolute reduction**, not the percentage, that moves risk.

| Baseline (nmol/L) | Niacin (−19%) | Evolocumab (−30%) | Pelacarsen (−79%) |
|---|---|---|---|
| 30 | −5.8 | −9.0 | −23.8 |
| 60 | −11.6 | −18.1 | −47.6 |
| 250 | −48.5 | −75.3 | **−198.3** |
| 400 | −77.6 | −120.4 | −317.2 |

At the median Lp(a) of AIM-HIGH/HPS2-THRIVE participants (about
30 mg/dL), a 20% reduction is barely 6 mg/dL. Niacin was **not an
underpowered Lp(a) drug — it was a sufficiently potent drug given to a
population no percentage reduction could rescue.**

The same table also explains why the enrolment criteria of
Lp(a)HORIZON (≥70 mg/dL) and OCEAN(a)-Outcomes (≥200 nmol/L) were set
where they were. **Baseline selection overwhelms percentage potency.**

---

## (F) The aortic valve — a timing problem, not a potency problem

Once valve calcification crosses a threshold (Hill n=4), it becomes
self-perpetuating independent of Lp(a).

```
dVCALC/dt = KCA · ( VIC_OST  +  KSELFR · VCALC⁴/(KSELFH⁴ + VCALC⁴) )
                   ↑ Lp(a)-dependent          ↑ Lp(a)-independent, threshold-type
```

A 62-year simulation (age 18 → 80):

| | AVC (AU) | AVA (cm²) | mean gradient | calcium prevented |
|---|---|---|---|---|
| Low Lp(a), 25 nmol/L | 43 | 3.21 | 3.9 | — |
| High Lp(a), untreated | 1558 | 0.82 | 58.9 | — |
| **Pelacarsen from age 30** | **206** | **2.45** | 6.7 | **87%** |
| **Pelacarsen from age 60** | **1327** | **0.93** | 46.3 | **15%** |

Starting at age 30 blocks 87% of the calcium. Starting at age 60 — even
with **20 years** of a drug that removes 80% of the particles — blocks
only 15%.

This is why SEAS, ASTRONOMER, and SALTIRE were all negative, and at the
same time a prediction that **no drug, however potent, can succeed against
the valve indication with a secondary-prevention design.**

---

## (G) The feed-forward loop is real but weak — and that is the conclusion

`Lp(a) → OxPL → monocyte NF-κB → IL-6 → IL-6 response element on the LPA
promoter → Lp(a)`

The model computes an **open-loop gain g = 0.019**, amplification 1.02 —
clinically negligible. Narratives that emphasise the loop are overstating
it.

What cannot be ignored is the **baseline IL-6 tone.** Even at an IL-6 of
only 2 pg/mL, the promoter multiplier is already 1.33, so **about 25% of
baseline LPA transcription is IL-6-driven.** What IL-6 blockade removes is
not the loop — it is this tone.

Changing only the IL-6 value in a single equation:

| | model | observed |
|---|---|---|
| Non-inflammatory patients (IL-6 2.7) | Lp(a) **−26.6%**, hsCRP **−92.5%** | RESCUE −16 to −25%, hsCRP −92% |
| Rheumatoid arthritis (IL-6 21) | Lp(a) **−33.5%** | MEASURE −37% |

---

## Validation — five defects found by an independent Python re-implementation

All 50 ODEs were independently re-implemented in Python with fixed-step
RK4 and cross-checked. The five items below are **real defects found and
fixed** in that process.

| # | Defect | Symptom | Fix |
|---|---|---|---|
| a | Assembly first-order in LDL | High-intensity statin drove Lp(a) to **−4.1%**, inverting the paradox | saturated with respect to LDL substrate (`KMLDL`) |
| b | PCSK9 turnover 5-fold too fast | Same steady state, but numerically stiff | `KSP=KDP=1.0`, `KONE` re-tuned |
| c | Anti-IL-6 modelled as partial competition | Rheumatoid arthritis showed a **smaller** Lp(a) reduction than non-inflammatory disease (backwards) | at clinical dose the antibody is present at roughly a 10⁴-fold molar excess over IL-6 → near-complete blockade |
| d | No IL-6-independent CRP production | predicted hsCRP −98% (RESCUE observed −92%) | added `FCRP0 = 0.06` |
| e | Valve self-perpetuation term was a ramp, not a threshold | **an 80-year-old with Lp(a) 25 nmol/L got AVC 4,425 AU · AVA 0.34 cm²** — severe aortic stenosis for everyone | replaced with a Hill n=4 threshold and recalibrated |

In addition, the risk-index term was bounded to a finite range, and the
mean-gradient output was protected at the limit where the valve is nearly
closed.

---

## Deliverables

| File | Contents |
|---|---|
| [`lpa_qsp_model.dot`](../../../elevated-lipoprotein-a/lpa_qsp_model.dot) · [`.svg`](../../../elevated-lipoprotein-a/lpa_qsp_model.svg) · [`.png`](../../../elevated-lipoprotein-a/lpa_qsp_model.png) | mechanistic map — **152 nodes · 17 clusters · 210 edges** |
| [`lpa_mrgsolve_model.R`](lpa_mrgsolve_model.R) | **50 ODEs** · 20 scenarios · 4 analysis runs · calibration notes |
| [`lpa_shiny_app.R`](lpa_shiny_app.R) | **10-tab** interactive dashboard |
| [`lpa_references.md`](lpa_references.md) | **112** PubMed citations · equation-to-literature mapping table |

### State variable composition (50 total)

| Module | Count | Contents |
|---|---|---|
| Pharmacokinetics | 17 | pelacarsen (hepatic compartment) · siRNA (RISC compartment) · muvalaplin · evolocumab · PCSK9 · statin · niacin · anti-IL-6 · obicetrapib |
| Production · assembly · catabolism | 7 | mRNA · ER apo(a) · free apo(a) · Lp(a) particle · LDL · VLDL · LDLR |
| Inflammation (ARM 2) | 4 | OxPL · trained monocytes · IL-6 · CRP |
| Vessel wall (ARM 1) | 8 | intimal retention ×2 · foam cells · necrotic core · fibrous cap · atheroma · CAC · vulnerability |
| Valve | 4 | autotaxin · LysoPA · osteogenic VIC · valve calcium |
| Fibrinolysis (ARM 3) | 3 | LBS occupancy · PAI-1 · lysis resistance |
| Metabolism · kidney | 2 | HDL-C · eGFR |
| Cumulative quantities | 5 | Lp(a) AUC · apoB AUC · MACE risk · AVS risk · urinary apo(a) fragment |

### Running it

```r
# model + 20 scenarios + 4 analysis runs
Rscript lpa_mrgsolve_model.R

# dashboard
shiny::runApp("lpa_shiny_app.R")

# regenerate the map
dot -Tsvg lpa_qsp_model.dot -o lpa_qsp_model.svg
dot -Tpng -Gdpi=150 lpa_qsp_model.dot -o lpa_qsp_model.png
```

---

## What this model cannot answer

Four items left honestly open. All four are exposed as parameters, so each
can be tested immediately once new data appears.

1. **The size of the KRAKEN assay discrepancy.** The direction is
   reproduced but the magnitude is not. Closing it would require
   `RFREE ≈ 3.8`, a measurable but unreported value.
2. **Where assembly happens.** The extracellular (White 1994) vs
   intracellular (Frischmann 2012) debate is unresolved. The model adopts
   the extracellular view; if intracellular assembly is correct, the
   predicted rise in free apo(a) under muvalaplin would change.
3. **The cholesterol fraction of Lp(a) mass.** Classic 0.30 vs more recent
   0.17–0.25. `FCHOL` leaves this uncertainty exposed as-is, and whether
   the Dahlén correction over- or under-corrects depends on it.
4. **Whether apoB spared by blocked assembly is recycled back into plasma
   LDL** (`FRECY`). The magnitude of apoB reduction in the pelacarsen
   trials suggests `FRECY` is close to zero.

And **what was deliberately left out**: no sex differences, no
triglyceride/LPL axis, a single expressed isoform assumed rather than two
alleles, no plaque geometry, no rs3798220–aspirin interaction. The
antifibrinolytic arm was given a **deliberately low weight**
(`WARM3 = 0.25`), however strong the in-vitro biology, because Mendelian
randomisation shows no association with venous thrombosis. The model
follows the evidence, not the biology.

---

## ⚠️ Disclaimer

A qualitative/semi-quantitative QSP model for educational and research
purposes. It was built from published literature and clinical trial data
but has not been independently verified or certified, and **must not be
used for clinical decision-making, prescribing, or regulatory
submission.**
