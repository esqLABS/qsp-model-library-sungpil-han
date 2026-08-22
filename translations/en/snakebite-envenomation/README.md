# Snakebite Envenoming (Snake Venom Poisoning) and Antivenom — Quantitative Systems Pharmacology Model

4 toxin classes × 3 antivenom fragment formats × 8 snake archetypes ·
50-ODE mrgsolve model · 27 treatment scenarios · 300-subject virtual population

---

## The single structural claim this model makes

> **Antivenom binds. It does not undo.**

Antivenom antibody can only remove toxin that is **still floating freely in plasma
or lymph**. It cannot pull a three-finger toxin off an acetylcholine receptor it has
already bound, it cannot splice back together fibrinogen that has already been
cleaved, and it cannot rebuild a motor nerve terminal that has already been
destroyed.

So this model is written as **two structurally separated clocks**.

| | Clock | Time unit | Does antivenom act here |
|---|---|---|---|
| **clock 1** | venom clock | **hours** | **Yes — only here** |
| **clock 2** | substrate clock | **days** | **No — it does not even appear in the equations** |

And **every clinical endpoint is read off the slow clock.** The time integral of the
fast clock drives the slow one.

Take that seriously and the following results follow *as arithmetic, not as claims*.

| # | Claim | What the model computed (computed, not assumed) |
|---|---|---|
| 1 | The timing of antivenom changes **the depth of the nadir** and not **the slope of recovery** | Fibrinogen nadir 0.06 → 2.40 g/L (**42-fold**). The rate of rise measured at the same fibrinogen value (1.0 g/L) is 0.46–0.96 g/L/day (**within 2-fold, all below the hepatic synthesis ceiling of 1.55 g/L/day**) |
| 2 | Recurrence is **the ratio of the fragment's half-life to the half-life of venom input**, not a vial count | At the same molar dose, swapping Fab → F(ab')₂ reduces late (>48 h) free-venom exposure **6873-fold** (1.7999 → 0.0003 mg·h/L) and detectable time from 221 h → 0 h |
| 3 | Presynaptic and postsynaptic paralysis **look the same and demand opposite answers** | Neostigmine: cobra (postsynaptic) mechanical ventilation 26 h → **4 h** (84% reduction) / krait (presynaptic) 58 h → **58 h** (0%) |
| 4 | The kidney is **an integral, not a number** | With antivenom at 4 hours, venom falls below the limit of detection by 9.5 h, yet creatinine goes on rising for **another 48 hours**, peaking at 58 h |
| 5 | There is **a compartment antivenom cannot reach** | Even with plasma free-venom AUC fully controlled at 13.10 mg·h/L, **18.1%** myonecrosis remains (k_b,loc/k_b0 = 0.037, a **27-fold penalty**) |
| 6 | The stoichiometric shortfall is **built into the label, not into the prescriber** | Ten vials of Indian ASV (60 mgNE) give a **0.98×** margin against a median Russell's viper bite and **0.70×** against a cobra. In the population simulation, **43% are underdosed at the label dose** |
| 7 | In cobra bite, **dose and timing are different questions** | 10 → 20 vials at 4 hours: mechanical ventilation 12 h → 11 h (no change), yet modelled mortality rises 0.5% → **1.5%** (only reaction risk has been added) |

---

## Files

| File | Contents |
|------|------|
| [`sbe_qsp_model.dot`](../../../snakebite-envenomation/sbe_qsp_model.dot) | Mechanistic map source — **145 nodes · 21 clusters · 245 edges** |
| [`sbe_qsp_model.svg`](../../../snakebite-envenomation/sbe_qsp_model.svg) | Zoomable vector map |
| [`sbe_qsp_model.png`](../../../snakebite-envenomation/sbe_qsp_model.png) | 150 dpi raster |
| [`sbe_mrgsolve_model.R`](../../../snakebite-envenomation/sbe_mrgsolve_model.R) | **50-ODE** mrgsolve model + 27 scenarios + 7 analysis functions + virtual population |
| [`sbe_shiny_app.R`](../../../snakebite-envenomation/sbe_shiny_app.R) | **15-tab** interactive dashboard |
| [`sbe_reference_model.py`](../../../snakebite-envenomation/sbe_reference_model.py) | Standalone Python implementation (for verification, see below) |
| [`sbe_reference_output.txt`](../../../snakebite-envenomation/sbe_reference_output.txt) | Run log of that implementation — **the source of every number in this README** |
| [`sbe_scenario_results.json`](../../../snakebite-envenomation/sbe_scenario_results.json) | Full results for all 27 scenarios (machine-readable) |
| [`sbe_population_results.json`](../../../snakebite-envenomation/sbe_population_results.json) | Results for the 6 population arms (machine-readable) |
| [`sbe_references.md`](../../../snakebite-envenomation/sbe_references.md) | **145 references** retrieved directly from PubMed, in 26 sections |

---

## Why a Python implementation ships alongside (Provenance)

The environment in which this repository was built **had no R runtime.** Rather than
commit an ODE model that had never been integrated, it seemed more honest to
implement every equation independently in Python, actually integrate it, and commit
those results. `sbe_mrgsolve_model.R` and `sbe_reference_model.py` are identical
term for term, and **all 166 parameters have been mechanically cross-checked** (every
value agrees). Every number in this README is the output of the latter
(`sbe_reference_output.txt`).

Doing it that way exposed **9 real defects**. Each is flagged in both files with a
`NOTE(defect n)` comment — whoever reads a QSP model has a right to know which lines
were hard.

| # | Defect | Why it looked plausible | How it surfaced |
|---|------|---------------------|-------------------|
| 1 | Setting the local necrosis rate constant `knec` to the same magnitude as the plasma rate constant | Both are "1/h per mg/L", so the units matched | The bite-site volume is 0.15 L, so local concentration is 100× plasma → every snake produced **more than 50% limb necrosis within an hour** |
| 2 | Making CK release proportional to the necrosis **flux** only | It is true that necrosis releases CK | CK peaked at **4–5 hours at 2×10⁵ U/L**. Reality is 12–48 hours. A necrotic fibre goes on releasing CK *for as long as it stays necrotic* → a **stock** term is needed |
| 3 | Sizing the fibrin term of renal loss `kn_fib` as though it were dimensionless | The coefficient was small, so it looked safe | The time integral of `FGloss` is **total g/L consumed** (≈5). Its contribution was 0.05 → effectively zero. **The mechanism of renal injury by intravascular coagulation was quietly missing from the model** |
| 4 | A refill term 10× too weak relative to capillary leak | Looking at the leak constant alone, it was defensible | Plasma volume 0.9 L, MAP 20 mmHg → **every untreated patient died of shock before there was even a chance to read a coagulopathy or an AKI** |
| 5 | **Writing neostigmine as an independent multiplier on the safety factor** (`SF = SF0·TERM·(1−BR)·neo`) | The expression was smooth and looked pharmacologically natural | It **"treated" presynaptic krait bite, taking mechanical ventilation from 28 h → 4 h** — destroying claim 3 exactly. The error was not numerical but **mechanistic**: raising acetylcholine acts by *displacing a competitive antagonist*, so it must go **inside** the occupancy term |
| 6 | Giving venom a conventional two-compartment **systemic** distribution | The computed terminal half-life of 49.5 h agreed with the literature | **Claim 2 did not reproduce.** The depot emptied within 24 hours and the deep compartment was nearly empty, so not even Fab could be exhausted while venom was still present. The defect was **structural**: antigenaemia lasting days is not a distribution phenomenon but an **absorption** one → the peripheral compartment was redesigned as a **slow depot** |
| 7 | Detecting recurrence by comparison with 5% of the acute peak | A peak-based criterion looked objective | It reported a Fab arm in which free venom had **risen 3.6-fold from its trough** as "no recurrence". The acute peak belongs to a different regime (pre-neutralisation) → the reference must be the **post-treatment trough** |
| 8 | No absolute floor on the recurrence fold-change | Fold-changes are scale-invariant, so it looked safe | A **2700-fold "recurrence"** from a trough of 10⁻⁹ mg/L — a value no instrument can see and no patient can feel → an **assay limit of detection (2 ng/mL) gate** was introduced |
| 9 | No baseline production term for myoglobin | With a release term and two elimination terms it looked complete | In a patient with no bite at all, **plasma myoglobin decayed to zero over 14 days** — meaning the model's baseline state was not a fixed point (the classic sign of a missing zero-order production term). For the same reason the pigment cast compartment had its initial value corrected to its own steady state |

**Defects 5 and 6 were destroying the central claims while looking plausible.**
Had the model not been integrated, both would have been committed.

---

## State variables (50 ODEs)

| # | State | Unit | Which clock |
|---|------|------|-----------|
| 1–4 | `D_i` fast bite-site depot (4 classes) | mg | venom |
| 5–8 | `P_i` **slow / sequestered depot** (4 classes) | mg | venom |
| 9–12 | `V_i` plasma free toxin | mg | venom |
| 13–16 | `C_i` antivenom–toxin complex | mg | venom |
| 17–19 | `A_c` `A_p` `A_t` antivenom (central · peripheral · bite-site tissue) | mgNE | venom |
| 20–24 | Neostigmine (2) · varespladib (2) · tranexamic acid (1) | mg | — |
| 25–29 | `FG` `FX` `PLT` `XDP` `SYNUP` coagulation | g/L, fraction, 10⁹/L, µg/mL, × | **substrate** |
| 30–31 | `BR` receptor occupancy · `TERM` terminal integrity | 0–1 | **substrate** |
| 32–35 | `NEC` `EDEMA` `CK` `MB` local · muscle | 0–1, L, U/L, µg/mL | **substrate** |
| 36–39 | `NEPH` `SCR` `CAST` `FIBR` kidney | 0–1, mg/dL, AU, 0–1 | **substrate** |
| 40–41 | `PV` `GLX` circulation | L, 0–1 | **substrate** |
| 42–43 | `IL6` `TNF` inflammation | pg/mL | **substrate** |
| 44–46 | `MCA` `IC` `SS` antivenom adverse reactions | AU | **substrate** |
| 47–50 | `HBLD` `HDTH` `VAUC` `AVCUM` cumulative quantities | — | — |

### Unit convention — why mgNE

Antivenom is handled in units of **labelled potency**.
`1 mgNE` = the capacity to neutralise 1 mg of venom.

Indian polyvalent antivenom (ASV) is, per its regulatory documentation, **capable of
neutralising 0.6 mg of *Daboia russelii* venom per mL**, i.e. **6 mgNE** per 10 mL
vial. Using this unit makes antivenom potency *a regulatory fact rather than a fitted
parameter*, and turns "the stoichiometric margin of the standard vial count" from an
assumption into **a computed result**.

---

## Computed results (computed, not asserted)

Every table below is transcribed verbatim from `sbe_reference_output.txt`.

### A. Half-lives are computed, not asserted

| Entity | t½ α | **t½ terminal** | Ratio to venom input (44.7 h) |
|---|---|---|---|
| ovine Fab (CroFab type) | 0.88 h | **21.8 h** | 0.49 → **loses the race** |
| equine F(ab')₂ (ANAVIP type) | 2.65 h | **126.8 h** | 2.84 → wins |
| equine whole IgG (Indian ASV type) | 3.25 h | **136.8 h** | 3.06 → wins |

On the venom side, after the correction of defect 6, the process is
**absorption-limited (flip-flop)**.

| Class | MW | f_slow | Fast absorption t½ | **Slow depot t½** | Systemic elimination t½ | **Apparent terminal t½** | ε |
|---|---|---|---|---|---|---|---|
| SVMP | 50 | 0.35 | 2.77 h | 57.8 h | 8.09 h | **44.7 h** | 1.20 |
| SVSP | 30 | 0.32 | 1.98 h | 46.2 h | 6.07 h | **36.5 h** | 1.00 |
| PLA2 | 14 | 0.15 | 1.39 h | 23.1 h | 3.08 h | **19.3 h** | 0.60 |
| 3FTx | 7 | 0.05 | 0.58 h | 6.9 h | 1.85 h | **6.0 h** | 0.35 |

> **Venom leaves the bite site more slowly than it leaves the body.**
> The apparent terminal half-life is 5.5× the elimination half-life, and this is
> the one antivenom has to outlast.

### B. "How many vials is enough" as arithmetic

Requirement = Σ (mass per class / ε)

| Species | Venom mg | mgNE required | IgG vials | Margin at 10 vials (60 mgNE) |
|---|---|---|---|---|
| *Daboia russelii* (Sri Lanka) | 63 | **61.0** | 10.2 | **0.98×** |
| *Echis ocellatus* | 22 | 18.0 | 3.0 | 3.33× |
| ***Naja naja*** | **40** | **85.7** | **14.3** | **0.70×** |
| *Bungarus caeruleus* | 10 | 17.1 | 2.8 | 3.51× |
| *Crotalus atrox* | 55 | 34.3 | 5.7 | 1.75× |
| *Bothrops asper* | 50 | 48.8 | 8.1 | 1.23× |
| *Bitis arietans* | 90 | 74.1 | 12.3 | 0.81× |

**Please read the Daboia row and the Naja row side by side.** The cobra has **less**
venom protein than the Russell's viper (40 mg vs 63 mg) and yet needs **more**
neutralising capacity (85.7 vs 61.0 mgNE). This is because 55% of the mass is the
class antivenom covers worst (3FTx, ε = 0.35). **The clinical practice of using
20–30 vials for a neurotoxic cobra bite is not caution, it is stoichiometry.**

### C. Claim 1 — depth vs slope

| Antivenom | FG nadir | Incoagulable h | **Rate of rise at FG 1.0 g/L** | FG day 7 | SCr peak |
|---|---|---|---|---|---|
| None | **0.06** | 39.2 | 0.46 g/L/day | 2.58 | 2.12 |
| 10 v @ 1 h | **2.40** | 0.0 | n/a (never drops below 1) | 2.93 | 0.97 |
| 10 v @ 4 h | **0.58** | 0.0 | 0.82 g/L/day | 3.24 | 1.31 |
| 10 v @ 12 h | **0.06** | 18.5 | 0.96 g/L/day | 3.27 | 1.60 |
| 20 v @ 4 h | **0.61** | 0.0 | 0.81 g/L/day | 3.24 | 1.30 |

The nadir spreads **42-fold**, while the rate of rise *measured at the same
fibrinogen value* is 0.46–0.96 g/L/day (within 2-fold), all of it below the hepatic
synthesis ceiling of 1.55 g/L/day. Even the remaining 2-fold difference is not due to
antivenom — it is because in the untreated and late-treated arms **residual venom goes
on consuming fibrinogen even while it is rising**, and that is why the fastest rise
belongs to the **12-hour arm** (the one with the deepest deficit → the largest
acute-phase SYNUP).

Doubling the vials at the same 4-hour delay: nadir 0.58 → 0.61 g/L. **Timing is not
dose.** Substituting for the liver (cryoprecipitate), by contrast, takes P(bleed) from
3.1% → **0.6%**.

Tranexamic acid: nadir 0.576 → 0.580 g/L, P(bleed) 3.1% → 3.1%. **Unchanged to three
significant figures** — because 95% of fibrinogen loss is direct enzymatic cleavage
that plasmin cannot touch.

### D. Claim 2 — recurrence is a margin, not a vial count

| Regimen | Margin | Rebound fold | Time of recurrence | **Late (>48 h) free-venom AUC** | Detectable h |
|---|---|---|---|---|---|
| Fab 6 v, 1× bite | 1.40× | 3.6× | **124 h** | **1.7999** | 221 |
| F(ab')₂ 6 v, 1× bite | 1.40× | <LOD | none | **0.0003** | 0 |
| Fab 12 v, 1× bite | 2.80× | <LOD | none | 0.0000 | 0 |
| **Fab 12 v, 2× bite** | **1.40×** | 5.8× | **135 h** | **2.7901** | 259 |
| F(ab')₂ 12 v, 2× bite | 1.40× | <LOD | none | 0.0001 | 0 |
| Fab 6 v + label maintenance dose | 2.80× | 1.1× | none | 0.0000 | 0 |

At **the same molar dose and the same margin**, replacing Fab with F(ab')₂ cuts late
free-venom exposure **6873-fold** and takes detectable time from 221 h → 0 h.

Doubling the Fab dose abolishes recurrence for a 1× bite (margin 2.80×). But
**apply that same doubled dose to a 2× bite and the margin is back to 1.40×, and
recurrence comes straight back** (135 h). At the same 1.40× margin, the long fragment
does not recur.

> **A fixed vial count is an incomplete definition.** Whether recurrence happens
> depends on how much venom that patient happened to receive, and that cannot be
> measured at the bedside.

### E. Claim 3 — two paralyses, the same findings, opposite pharmacology

| Case | AChR occupancy | Terminal integrity | Ptosis h | Ventilation h | P(death) |
|---|---|---|---|---|---|
| Cobra, no antivenom | 0.822 | 0.763 | 58 | 26 | 1.2% |
| Cobra, **neostigmine only** | 0.822 | 0.763 | 56 | **4** | 0.3% |
| Cobra, ASV 10 v @ 4 h | 0.779 | 0.892 | 37 | 12 | 0.5% |
| Cobra, ASV 10 v **+ neostigmine** | 0.779 | 0.892 | **19** | **2** | 0.1% |
| Cobra, ASV 20 v @ 4 h | 0.775 | 0.895 | 36 | 11 | **1.5%** |
| Krait, ASV 10 v @ 6 h | 0.164 | **0.146** | 143 | 58 | 2.3% |
| Krait, ASV 10 v **+ neostigmine** | 0.164 | **0.146** | **143** | **58** | 2.3% |
| Krait, ASV 10 v **@ 0.5 h** | 0.016 | **0.931** | **0** | **0** | 0.1% |
| Krait, ASV 10 v @ 6 h, **no mechanical ventilation** | 0.164 | 0.146 | 143 | 58 | **99.4%** |

The same drug, two cases:

- Cobra (postsynaptic), no antivenom: mechanical ventilation **26 h → 4 h (84% reduction)**
- Krait (presynaptic), with antivenom: mechanical ventilation **58 h → 58 h (0%)**, ptosis **143 h → 143 h (0%)**

The bedside findings are the same. The drug acts on **occupancy**, but in the krait
occupancy is 0.16 and the deficit is *the terminal having been destroyed to 0.15 of
normal*. There is nothing for acetylcholine to displace.

**Dose and timing are different questions.** 10 → 20 vials at 4 hours: mechanical
ventilation 12 h → 11 h. By 4 hours the receptors are already 78% occupied and
occupancy decays on its own 35-hour clock. Yet modelled mortality goes **0.5% →
1.5%** — because 20 vials of whole equine IgG cross the reaction threshold. This is
not model noise, it is **the fact that the ledger has two sides**.

The krait shows the same asymmetry from the other side. Same product, same vial
count, **5.5 hours apart**: dosing at 6.0 h → ventilation 58 h, terminal 0.146 /
dosing at 0.5 h → ventilation 0 h, terminal 0.931.

And what actually saves a krait patient is not antivenom:
**with mechanical ventilation 2.3% vs without 99.4%.**

### F. Claim 4 — the kidney is an integral

| Antivenom | Free-venom AUC | SCr peak | Time | eGFR nadir | Permanent scarring |
|---|---|---|---|---|---|
| None | 95.45 | 2.12 | 103 h | 38 | 12.1% |
| 10 v @ 1 h | 2.41 | 0.97 | 59 h | 92 | 1.5% |
| 10 v @ 4 h | 15.98 | 1.31 | 58 h | 62 | 5.1% |
| 10 v @ 12 h | 47.36 | 1.60 | 68 h | 49 | 7.5% |

Treat at 4 hours and free venom falls below the assay limit of detection (2 ng/mL) by
**9.5 h**, yet creatinine goes on rising for **another 48 hours**, peaking at 58 h.
There is nothing left in plasma to measure.

Untreated, the same assay says something different: because the slow depot keeps
releasing, free venom stays above the limit of detection **out to 332 hours**. So **a
positive antigen assay on day 5 does not mean treatment failure, and a negative assay
on day 1 does not mean the kidney is safe.**

### G. Claim 5 — where antivenom cannot go, and what can

| Regimen | Myonecrosis | CK peak | Compartment pressure | Free-venom AUC |
|---|---|---|---|---|
| F(ab')₂ 10 v @ 4 h | **18.1%** | 6509 U/L | 29 mmHg | 13.10 |
| F(ab')₂ 10 v @ 4 h **+ varespladib @ 0.5 h** | **13.8%** | 4170 U/L | 28 mmHg | 13.10 |

Intravenous antibody controls the plasma compartment **completely** (free-venom AUC
13.10) and still leaves 18.1% myonecrosis. `k_b,loc / k_b0 = 0.037` — a **27-fold
penalty** at the bite site, plus a further 4-fold penalty in the slow depot.

Start an oral PLA2 catalytic-site inhibitor with a volume of distribution of 12 L at
0.5 hours and necrosis becomes 13.8% (24% reduction), CK 6509 → 4170 U/L. **This is
the entire rationale for a prehospital small-molecule programme** — a small molecule
goes where antibody cannot, and it can be swallowed before transport.

### H. The cost side of the ledger — treating a dry bite as envenoming

- **Benefit**: free-venom AUC 0.000 mg·h/L (because there was no venom)
- **Cost**: anaphylactoid index 0.89, serum sickness 1.25 (peaking at day 6.1),
  modelled P(death) 0.06% — **all of it iatrogenic**

20–50% of venomous snakebites (10–30% for elapids) inject no venom.
The 20-minute whole blood clotting test needs nothing but a single glass tube, and
that is the gate.
**This row exists so that this model is not read as an argument for "treat every bite".**

---

## What this model fails to reproduce (and what would have to be true)

This is the most important section. **It is reported as a negative result rather than
parameterised away.**

The Fab arm reproduces **recurrent venom antigenaemia** cleanly (3.6-fold rebound from
the trough, peaking on day 5). But it does **not** reproduce the **severe recurrent
hypofibrinogenaemia** the literature reports. The reason is arithmetic.

| Quantity | Value |
|---|---|
| Hepatic fibrinogen synthesis ceiling | **1.55 g/L/day** |
| Maximum late (>72 h) fibrinogen consumption rate | **0.094 g/L/day** (6.0% of the ceiling) |
| Shortfall factor needed for a second descent | **about 17-fold** |

Two routes were swept separately.

**Route 1 — sequester more into the slow depot.** `f_slow` cannot exceed 1.0 (the
whole bite sequestered), so the headroom is bounded. **Even sequestering the entire
bite, the late flux only reaches 57% of the synthesis ceiling.** Depot release cannot
be the mechanism of clinically recurrent coagulopathy.

**Route 2 — the antivenom–toxin complex dissociates.** Raising the default k_off = 0
all the way to 0.15/h does not reproduce it, and **the reason is instructive**: the
complex is cleared with a half-life of 11.6 h, so only 1% remains at 72 hours.
Dissociation releases venom **early**, while the complex pool is still large,
prolonging the acute phase rather than creating a late event.

**Conclusion (stated as a negative result):** two independently testable
interpretations remain.

1. The recurrent hypofibrinogenaemia reported after ovine Fab is mostly a
   **laboratory event** (a dip to the 0.1–0.4 g/L level, which this model does in
   fact generate) and not a return to incoagulability; or
2. **there is a separate venom reservoir** that has no compartment in this model.

It would have been easier to raise `f_slow` to 1.0 and call it "calibration", but that
is burying a real discrepancy under a plausible number.

---

## Virtual population — why the median patient cannot reproduce trial incidences

Every endpoint in this model is a **threshold**, and **the mean of a threshold is not
the threshold of the mean.** Six arms of 300 subjects (the dominant drivers of
inter-individual variability are injected venom dose, CV 55%, and delay to
presentation):

| Arm | Incoagulable | Recurrence | Late detectable | AKI (≥1.5×) | AKI stage 3 | Ventilation | Bleeding | Death | **Label dose insufficient** |
|---|---|---|---|---|---|---|---|---|---|
| Daboia + ASV 10 v, delay ~4 h | 37.7% | 15.3% | 22.3% | 37.0% | 2.7% | 6.3% | 8.7% | 2.6% | **43.0%** |
| Daboia + ASV 10 v, delay ~1 h | 15.7% | 21.0% | 33.3% | 27.0% | 2.0% | 2.0% | 7.3% | 2.3% | 37.7% |
| Daboia + ASV 20 v, delay ~4 h | 27.0% | 1.3% | 1.7% | 23.7% | 0.0% | 2.3% | 4.0% | 2.0% | 5.7% |
| Daboia, no antivenom | 98.7% | 0.0% | 100.0% | 100.0% | 12.0% | 25.0% | 39.7% | 12.4% | 100.0% |
| *C. atrox* + Fab 6 v @ 2 h | 7.3% | **24.7%** | **49.7%** | 21.0% | 0.3% | 0.0% | 3.6% | 1.1% | 20.3% |
| *C. atrox* + F(ab')₂ 6 v @ 2 h | 10.3% | **12.0%** | **17.0%** | 12.7% | 0.7% | 0.0% | 4.3% | 1.3% | 19.3% |

The last column is this model's answer. **The label's 10 vials give a 0.98× margin
against the *median* bite, so close to half of the real distribution of bites is
structurally underdosed.** That 30–50% of Russell's viper bites require repeat dosing
is not antivenom failing — it is **because the label is set to the median**.

---

## Conservation and consistency checks

- Venom injected 48.51 mg → residual in the body after 14 days 0.0000 mg (0.000%).
  No mass is trapped in a non-physical sink.
- All bounded states (`BR` `TERM` `NEC` `NEPH` `GLX` `FIBR`) stay within [0, 1].
- **Baseline fixed-point check** (dry bite, 14 days): maximum drift below 0.6%.
  This check is what caught defect 9.
- **R↔Python parameter cross-check**: all 166 values agree.

---

## The clinical decision layer — the whole model as four questions

This is cluster 21 of the map, and it is the summary of every calculation above.

1. **Is this envenoming?** 20WBCT + check for ptosis + limb oedema. Treat a dry bite
   as envenoming and you take all of the antivenom risk with none of the benefit.
   **This question is free.**
2. **Has the substrate already gone?** Antivenom changes the depth of the nadir, not
   the slope of recovery. Late antivenom still prevents *the consumption of the next
   few hours*, but that is all.
3. **Will the cover outlast the venom?** Compare t½(fragment) with t½(venom input).
   If the fragment is the one that loses, **plan the maintenance schedule at the time
   of the first dose** — not after recurrence has been discovered.
4. **What is the act that actually saves the life?** Viperid coagulopathy →
   antivenom / elapid presynaptic paralysis → **mechanical ventilation** / capillary
   leak → fluid / venom AKI → dialysis. **Antivenom is necessary, but rarely
   sufficient.**

---

## How to run

```r
# mrgsolve model
source("sbe_mrgsolve_model.R")
res <- run_all()                       # 27 scenarios
print(analyse_depth_vs_slope(res))     # claim 1
print(analyse_recurrence(res))         # claim 2
print(analyse_two_paralyses(res))      # claim 3
print(analyse_kidney_integral(res))    # claim 4
print(analyse_local(res))              # claim 5
print(analyse_stoichiometry())         # how many vials is enough
print(analyse_halflives())             # the recurrence inequality
print(virtual_population(n = 300))     # population

# Shiny dashboard (15 tabs)
shiny::runApp("sbe_shiny_app.R")
```

```bash
# standalone Python verification implementation (full reproduction without R)
python3 sbe_reference_model.py     # → sbe_reference_output.txt

# render the mechanistic map
dot -Tsvg sbe_qsp_model.dot -o sbe_qsp_model.svg
dot -Tpng -Gdpi=150 sbe_qsp_model.dot -o sbe_qsp_model.png
```

---

## ⚠️ Disclaimer

This model is a **semi-quantitative QSP model for educational and research purposes**.
It was constructed from the published literature and clinical data, but it has not
been independently verified or certified, and **it must not be used for actual
clinical decision-making, prescribing, or regulatory submission.**

The following values in particular are **model parameters, not observations**: the ε
vector (SVMP 1.20 / SVSP 1.00 / PLA2 0.60 / 3FTx 0.35), the slow-depot
compartmentalisation (`f_slow`, `ka_s`), the magnitude of the local neutralisation
penalty at the bite site, the neuromuscular safety factor and its thresholds, every
rate constant in the renal loss terms, and the baseline hazard of every hazard
function. The full separation is in the "What is sourced, and what is not" section of
[`sbe_references.md`](../../../snakebite-envenomation/sbe_references.md).

**Snakebite is an emergency.** Follow local WHO or national treatment guidelines and
the antivenom manufacturer's label, not this model.
