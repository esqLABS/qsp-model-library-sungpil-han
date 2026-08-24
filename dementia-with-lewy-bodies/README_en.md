# Dementia with Lewy Bodies — QSP Model

<a href="dlb_qsp_model.svg"><img src="dlb_qsp_model.png" width="620" alt="DLB QSP mechanistic map"></a>

**66-ODE mrgsolve model · 21-cluster 190-node 311-edge mechanistic map · 23 treatment scenarios · 11-tab Shiny dashboard · 154 PubMed references**

| File | Contents |
|------|------|
| [`dlb_qsp_model.dot`](dlb_qsp_model.dot) · [`.svg`](dlb_qsp_model.svg) · [`.png`](dlb_qsp_model.png) | Mechanistic map (Graphviz, 21 clusters / 190 nodes / 311 edges) |
| [`dlb_mrgsolve_model.R`](dlb_mrgsolve_model.R) | 66-compartment ODE model + 23 scenarios + mechanism verification functions |
| [`dlb_shiny_app_en.R`](dlb_shiny_app_en.R) | 11-tab interactive dashboard |
| [`dlb_references_en.md`](dlb_references_en.md) | 154 PubMed links (classified by section, with the model location marked) |

---

## 1. What this model sets out to do

The pharmacotherapy of dementia with Lewy bodies (DLB) is a list of facts that look as though
they contradict one another.

- Cholinesterase inhibitors (ChEIs) **work better than they do in Alzheimer's disease.**
- Levodopa **works less well than it does in Parkinson's disease.**
- Antipsychotics **can be lethal** (severe neuroleptic sensitivity, 2–3-fold mortality).
- 5-HT2A inverse agonists (pimavanserin) **work on the visual hallucinations but do not abolish them.**
- ChEI trials report **larger effects on neuropsychiatric and cognitive-fluctuation measures** than on cognitive scores.

The textbooks list these as separate items. The claim of this model is that **the five are
consequences of one structural fact**.

> **The presynaptic lesion is severe in all three ascending transmitter systems, but the
> postsynaptic receptor landscape is remodelled in a different direction in each system.**

| System | Presynaptic | Postsynaptic | Consequence |
|--------|----------|----------|------|
| Cholinergic | nbM loss severe (greater than in AD) | **M1/M4 preserved** | ChEIs work better than in AD |
| Dopaminergic | SNc loss severe | **no D2 upregulation** | levodopa works less well than in PD, and there is no reserve against D2 blockade |
| Serotonergic | dorsal raphe loss mild | **5-HT2A upregulated** | hallucinations arise, and inverse agonists work |

The model uses the same three-factor transducer

```
DRIVE = presynaptic capacity × postsynaptic density × (1 − drug occupancy)
```

**three times**. The only thing that differs is the differential equation for the middle factor.

```c
dxdt_M1R  :  target = M1BASE − KTAUM1·PTAU              // responds to tau only → preserved in DLB
dxdt_D2R  :  target = 1 + KD2UP·denervation·UPCAP       // UPCAP is what limbic α-syn kills
dxdt_HT2A :  target = 1 + KH2UP·denervation
                        + KH2NEO·(neocortical fibrils)  // goes up
```

Inside `dlb_mrgsolve_model.R` there is **no** rule saying "ChEIs work better in DLB" or
"DLB patients are sensitive to neuroleptics." Those facts are outputs.

---

## 2. The four structural commitments

### (1) One transducer, three times, with different signs

Run scenarios 12 · 13 · 14 side by side. The risperidone 1 mg/day exposure is **exactly
identical**, the only thing that changes is `PHENO`, and `PHENO` sets nothing but the initial
pathology distribution (the starting values of brainstem/limbic/neocortical α-synuclein, tau
and amyloid).

| | D2 occupancy | Postsynaptic reserve | Drug-induced deficit | Sensitivity index | ΔUPDRS-III |
|---|---|---|---|---|---|
| DLB | 0.58 | 0.114 | 1.402 | **0.463** | **+21.4** |
| PD dementia | 0.58 | 0.093 | 1.422 | **0.504** | **+18.9** |
| AD | 0.58 | 0.759 | 0.804 | **0.006** | +7.1 |

The D2 occupancy is the same in all three groups. The only thing that differs is the
**reserve**, and the deficit is amplified the less reserve there is. That the model predicts
DLB and PD dementia to be similar and starkly different from AD is consistent with the
literature — neuroleptic sensitivity is not specific to DLB but **a feature of Lewy body
disease as a whole** (McKeith 1992; Aarsland 2005).

> **A side prediction.** In the unmedicated state `OCCD2 = 0`, so `DEFICIT = 0` and
> `NSENS = 0` — however far the disease has progressed. That is, in this model neuroleptic
> sensitivity is **a drug event, not a severity milestone**. This is not a post-hoc
> rationalisation but a result of the shape of the equations.

### (2) Cognitive fluctuation is variance, not severity

Attention (`ATTM`) is not a level that a drug lifts, but the state of a bistable system
sitting close to a saddle-node.

```c
dxdt_ATTM = KATT · ( −ATTM³ + ALPHAB·ATTM + DRIVEA − ATTOFF )
```

The "fluctuation" that gets scored in the clinic is a state **transition**. A transition
requires two states, so the right measure is not the curvature of the fixed point but **how
deep inside the bistable region the system sits** — for this cubic that band is
`|DRIVEA − ATTOFF| < 2(α/3)^1.5`, about 0.157 wide in units of drive.

```c
BISTAB  = max(1 − |DRIVEA − ATTOFF| / DCRITB, 0)
AROUSD  = arousal instability (LC + PPN loss)
FLUCTGT = SIGN0 · (FLBASE + BISTAB · (ANOISE0 + KNOISE·AROUSD))
```

Bistability supplies the **opportunity** for a transition, and the noise that actually causes
one comes from the brainstem arousal systems (locus coeruleus + PPN). α-Synuclein destroys
both of these and tau largely bypasses them — which is why fluctuation is a core feature of
DLB and incidental in AD. And that holds even though both systems pass through the same cubic.

**The testable prediction the model makes:** a ChEI improves the variance (CAF)
proportionally more than the mean (MMSE). Donepezil 10 mg, 20 weeks:

| | Absolute change | Relative change |
|---|---|---|
| MMSE | **+2.28** | +10.8% |
| CAF (cognitive fluctuation) | **−3.05** | −37.6% |
| | | **ratio 3.49-fold** |

This figure of 3.49 was not fitted to any data. It is a falsifiable claim.

Fluctuation is also **non-monotone** — greatest in mild-to-moderate disease (model: CAF 8.2 at
about 10 months after diagnosis) and diminishing again once the state has fallen completely
onto the lower branch. The first draft's linear model could not produce this (see §5).

### (3) Visual hallucinations are the product of three deficits, not their sum

```c
VHDRIVE = (1 − bottom-up evidence fidelity) × (1 − top-down attentional binding) × (VHFLOOR + 5HT2A signal)
```

All three factors have to go bad **at the same time**. Because it is a product, no single
mechanism explains the hallucinations and no single mechanism abolishes them. The
justification for this structure is that in the literature neither occipital hypometabolism
alone, nor cholinergic loss alone, nor 5-HT2A upregulation alone separates hallucinators from
non-hallucinators cleanly, whereas their **combination** does (the PAD model of Collerton
2005).

`VHFLOOR` exists by the same logic. Without it the model predicts that pimavanserin all but
abolishes the hallucinations, which is not what the trial found (SAPS-PD −3.06, about one
third). Model calculation: hallucination burden **−39.6%**.

### (4) The GBA1 loop is an amplifier, not a switch — retracting the first draft's claim

The first draft claimed that the GBA1 feedback loop (GCase↓ → GlcCer↑ → slower oligomer
degradation → oligomers↑ → GCase trafficking blocked → GCase↓) has a saddle-node bifurcation,
and that there is therefore a **deadline** for disease modification. Actually computing the
loop gain gave

```
d ln GlcCer/d ln GCase × d ln KOLDE/d ln GlcCer
                       × d ln OLIG/d ln KOLDE × d ln GCTGT/d ln OLIG  ≈  0.17–0.26
```

which is less than 1. That is, it is **monostable**. Raising the Hill coefficients anywhere
within the range the literature supports never took it past 1, and getting past 1 required
contrived parameters.

So the claim was changed, and a **falsification test** was left inside the model — scenario
20b: ambroxol for two years, then **stopped**. If it really were bistable, GCase should stay
on the rescued branch.

| Scenario | GCase 2 y | GCase 3 y | GCase 5 y | Neocortical fibrils 5 y | MMSE 3 y |
|---|---|---|---|---|---|
| 18 GBA1 carrier, untreated | 0.195 | 0.174 | 0.167 | 0.816 | 6.48 |
| 19 ambroxol from day 0 | 0.356 | 0.295 | 0.268 | 0.790 | **8.50** |
| 20 ambroxol from day 1095 | 0.195 | 0.174 | 0.266 | 0.796 | 6.48 |
| **20b ambroxol day 0–730, then stopped** | **0.356** | **0.181** | **0.168** | 0.814 | 8.07 |

Stop the drug and it returns **completely** to the untreated level (0.167). It is not bistable.

Early (19) and delayed (20) also have essentially the same GCase at five years (0.268 vs
0.266). What differs is **the integral of the time spent in the suppressed state**, and that
appears as the difference between MMSE 8.50 and 6.48 at three years. That is precisely what
is meant by calling this a question of lever length.

So the claim this model actually supports is **not "miss the deadline and it is over" but
"the earlier you start, the longer the lever"**. A weaker claim, but one that can be sustained.
(The dose-response to GBA1 genotype itself is robust — untreated neocortical fibril burden at
five years: GBAF 1.00 → 0.711, 0.75 → 0.788, 0.55 → 0.816.)

---

## 3. Natural-history calibration

The values below were **obtained by running the model** and are not assertions.
mrgsolve 2.0.1 / R 4.3.3, LSODA, rtol = atol = 1e-8.
Where the first draft disagreed with the literature, **the model** was fixed, not the target.

| Anchor | Literature | This model |
|--------|------|---------|
| MMSE at diagnosis | 18–24 (DLB trials) | 23.9 |
| Annual MMSE decline | −4 ~ −5 | **−3.84** |
| MDS-UPDRS III at diagnosis | 15–25 | 19.3 |
| Annual MDS-UPDRS III rise | +5 ~ +9 | **+7.38** |
| Median survival from diagnosis | 4–7 years | **4.59 years** |
| MIBG heart/mediastinum ratio (at diagnosis) | < 2.0 (abnormal) | **1.66** |
| DaTSCAN striatal binding ratio | about 40–60% of normal | **0.36** |
| EEG dominant frequency (at diagnosis → 5 years) | pre-alpha 5.6–7.9 Hz | **8.13 → 5.62 Hz** |
| Peak cognitive fluctuation | mild-to-moderate stage | **CAF 8.2 at 10 months**; CAF>5 sustained for 705 days |
| Time hallucinations cross the clinical threshold | commonly around diagnosis | 26 months (a weakness of the model, §5) |

### Profiles of the three phenotypes at diagnosis

The equations are identical; only the initial pathology distribution differs.

| | MMSE | UPDRS-III | CAF | Hallucinations (5 y) | MIBG | DaTSCAN | EEG (Hz) | D2 reserve |
|---|---|---|---|---|---|---|---|---|
| **DLB** | 23.9 | 19.3 | 5.23 | 12.55 | **1.66** | **0.36** | **8.13** | 0.225 |
| **PD dementia** | 23.9 | 26.4 | 7.62 | 5.95 | 1.58 | 0.19 | 7.74 | 0.128 |
| **AD** | 23.3 | 3.8 | 3.00 | 0.77 | **2.54** | **0.78** | **9.10** | 0.742 |

The three indicator biomarkers (MIBG · DaTSCAN · EEG) separate in exactly the diagnostically
correct directions. None of them is an assigned value; all are computed.

---

## 4. PK/PD calibration

Steady state (dosing days 350–400), 0.02-day grid.

| Drug | Quantity | Literature | This model |
|------|------|------|---------|
| Rivastigmine 6 mg BID oral | Cmax | 20–30 ng/mL | **19.4** |
| | Brain AChE inhibition | ~40–60% | 0.50–0.79 (oscillating) |
| | BuChE inhibition | co-inhibited | 0.62 |
| Rivastigmine 9.5 mg/24h patch | Css | 3–9 ng/mL | **2.1–4.2** |
| | AChE inhibition | ~60% | **0.52–0.61** |
| | GI adverse-effect index (Cmax-driven) | patch < capsule (IDEAL) | **23 vs 58** |
| Donepezil 10 mg qd | Css | 40–60 ng/mL | **41.0–49.7** |
| | Brain AChE inhibition | 60–70% | **0.647** |
| Pimavanserin 34 mg qd | Cavg | 40–60 ng/mL | **53.6** |
| | 5-HT2A occupancy | 85–95% | **0.896** |
| | QTc prolongation | +5 ~ +10 ms | **+8.0** |
| Levodopa/carbidopa 150 mg | Cmax | 1000–2000 ng/mL | **1618** |
| Risperidone 1 mg/day | Active-moiety Cavg | ~10 ng/mL | **10.0** |
| | Striatal D2 occupancy | 60–70% | **0.662** |
| Quetiapine 50 mg/day | D2 occupancy | very low, transient | **0.038** |
| Anti-α-syn antibody 4500 mg q4w | CSF:plasma ratio | 0.1–0.3% | **0.30%** |
| Ambroxol 1.26 g/day | Increase in GCase activity | +35% (Mullin 2020) | **+34%** |

### Rivastigmine: the plasma half-life separated from the pharmacodynamic half-life

The first draft wrote rivastigmine as acting through its plasma concentration. That makes the
duration of action 1.5 hours, which renders 12-hourly dosing meaningless. **The carbamylated
enzyme had to be made a separate state variable**, and its rate of recovery is set not by drug
elimination but by **enzyme resynthesis** (half-life about 9 hours).

```c
dxdt_CARBA = KCARB·CRIV·(1 − CARBA) − KDECARB·CARBA     // KDECARB = ln2 / (9 h)
```

A side consequence: because the GI adverse effects are driven by `Cmax` while efficacy is
driven by the exposure integrated by the carbamylated enzyme, the IDEAL study finding that
**the patch achieves the same inhibition at a far lower adverse-effect index** falls out
without being coded separately.

### Effects comparable with the clinical trials

| Scenario | ΔMMSE (12 weeks) | ΔMMSE (20 weeks) | ΔNPI | ΔCAF | ΔHallucinations (3 y) |
|---|---|---|---|---|---|
| 04 Rivastigmine 6 mg BID oral | +3.68 | +3.91 | −32.9% | −63.0% | −59.1% |
| 05 Rivastigmine 9.5 mg/24h patch | +3.42 | +3.56 | −30.7% | −61.4% | −53.4% |
| 06 Donepezil 10 mg qd | **+2.26** | +2.27 | −19.0% | −37.4% | −44.9% |
| 07 Donepezil + memantine | +2.31 | +2.36 | −19.0% | −37.4% | −45.6% |
| 09 Donepezil + oxybutynin | **+0.26** | +0.46 | −16.6% | −52.7% | +0.1% |
| 10 Pimavanserin 34 mg | 0.00 | 0.00 | 0.0% | 0.0% | **−39.6%** |
| 21 Anti-α-syn antibody | +0.02 | +0.03 | −0.1% | +0.4% | −7.8% |
| 22 Optimal combination | +3.46 | +3.61 | **−43.5%** | −63.0% | **−73.8%** |

- **Donepezil +2.26 (12 weeks)** is essentially identical to the +2.2 of Mori 2012.
- **Rivastigmine ΔNPI −32.9%** agrees with the roughly 30% of McKeith 2000.
- **Donepezil + an anticholinergic**: same drug, same dose, and the effect collapses from
  +2.26 to +0.26. This is the result of one term, `ACHSIG = ACHS × M1R × (1 − 0.18·ANTICH)`,
  and it is the quantitative expression of the commonest iatrogenic error in the clinic.
- **Phenotype dependence of the ChEI effect**: donepezil's 12-week ΔMMSE is **+2.26** in DLB
  and **+1.76** in AD (ratio 1.29). Once tau has taken the postsynaptic M1 away, the
  transducer has less left to multiply.

### Levodopa — the same exposure, a different receptor landscape

| | Untreated UPDRS-III | Levodopa UPDRS-III | Change |
|---|---|---|---|
| DLB | 27.4 | 22.3 | **−18.5%** |
| PD dementia | 28.0 | 21.9 | **−21.9%** |
| Responder rate (≥20% improvement, n=400) | | | **48.0% vs 48.2% — see below** |

In the model, three things act **at the same time** to blunt the response in DLB:
(a) limbic α-synuclein lowers postsynaptic striatal integrity (`UPCAP`),
(b) D2 upregulation does not happen,
(c) a substantial part of DLB parkinsonism is not dopaminergic in the first place —
levodopa cannot touch the `MOTNDA·(FIBN + 0.5·FIBL)` term.
At diagnosis this non-dopaminergic share is about 30% of UPDRS-III in DLB and about 7% in PD
dementia.

> **A negative result is reported as a negative result.** The model separates the two
> phenotypes on **mean change** (−11.9% vs −14.7%, n=400) but fails to separate them on
> **responder rate** (48.0% vs 48.2%). The frequently quoted contrast of "about one third in
> DLB versus about 90% in PD" is a comparison of DLB with **early Parkinson's disease**, and
> this model has no early-PD arm — levodopa responsiveness in PD dementia is itself already
> blunted. Reproducing that contrast would require adding a fourth phenotype, and that was not
> done. Do not cite 48%/48% as the basis for any claim.

---

## 5. What the first draft got wrong

This section is the most useful part of this README. Everything below was **found by running
the model and then fixed**.

1. **The aggregation loop diverged.** The first version had no capacity-limiting term, so the
   fibril burden ran away from 0.35 to 27 in six months and the system went `NaN`. Unbounded
   autocatalysis is not a disease model, it is a divergence. The `(1 − OLIG)` and `(1 − FIB)`
   saturation terms were added to make it sigmoidal.
2. **Fibril release was wrongly subtracted as a loss of mass.** Subtracting `− KREL·FIB` in
   the fibril equation made the whole pathology actually **regress** (0.35 → 0.066). Seed
   release is a trace secretory flux, not a mass sink.
3. **The attention model was linear.** That predicts fluctuation growing indefinitely in
   proportion to severity, whereas clinically fluctuation is most prominent in mild-to-moderate
   disease and diminishes in the severe stage. The cubic produces this non-monotonicity, and
   the linear model could not.
4. **Measuring fluctuation by curvature gave a peak that lasted only 3.5 months.**
   `1/√curvature` becomes singular only at the bifurcation point, so it is far too narrow on
   the time axis. Switching to the **depth** inside the bistable region gave a clinically
   correct window lasting about two years (705 days).
5. **The phenotype was a label.** At first `PHENO` changed only the connectivity weights,
   while the initial fibril burden was identical in all three arms. So DLB and PD dementia
   were **indistinguishable** in levodopa response and in neuroleptic sensitivity. Making the
   phenotype set the initial **regional distribution** solved both problems at once.
6. **Neuroleptic sensitivity was written as an absolute threshold on the residual signal.**
   That made an advanced DLB patient mount a sensitivity reaction **without taking any drug**.
   Separating reserve (postsynaptic capacity in the unmedicated state) from deficit (the block
   the drug creates, amplified by reserve) gave `OCCD2 = 0 ⇒ NSENS = 0`, which is correct.
7. **The motor score was `(1 − dopamine signal)²`.** With that, even a healthy striatum
   reaches UPDRS-III 24 once D2 is 58% blocked. In reality a healthy striatum has reserve and
   shows only mild extrapyramidal signs. Switching to a sigmoidal reserve curve plus a
   terminal-loss exponent (`TERMEXP = 1.6`, the same exponent as the DaTSCAN mapping)
   produced AD +7.1 versus DLB +21.4.
8. **The levodopa conversion gain was 10-fold too large.** In the first version 300 mg/day
   almost doubled striatal dopamine. It was necessary to build in the fact that AADC and
   VMAT2 sit **in the very terminals that have been lost** (`AADCF`, `VMATF`).
9. **The antibody's CSF concentration came out at 88 mg/L.** That was down to using CSF as a
   mass sink for plasma. Rewriting it as a **partition** gave a steady-state CSF/plasma of
   0.3%.
10. **`$OMEGA` was declared and never used.** Inter-individual variability was applied to no
    parameter, so the "virtual patient population" was all the same person and responder rates
    came out as either 0% or 100%. ETAs were put on the residual neuron counts and the
    pathology burden at diagnosis.
11. **Then the opposite problem appeared.** The random effects leaked into the deterministic
    scenario comparisons and made comparison between the `PHENO` control arms impossible.
    `run_scn()` was made to use `zero_re()`, and variability is switched on only in
    `responder_rates()`, **where the distribution is the object of interest**.
12. **The GBA1 switch was not a switch.** §2(4) above.
13. **The initial conditions were not at equilibrium.** The first six months of every scenario
    were the model settling rather than the disease progressing. The fast states (ACh, M1,
    attention state, fluctuation, cognition, motor, RBD, autonomic, somnolence) were placed at
    their own equilibria in `$MAIN` — the attention state is a cubic equation, so it is solved
    by Newton iteration.
14. **That Newton iteration initially converged onto the wrong branch.** Picking the starting
    value from the sign of the drive made a newly diagnosed patient **start on the lower
    branch**, giving an MMSE 4 points too low. A patient at diagnosis has come down from
    health and has not yet fallen — which branch they are on is a question of **history, not
    of the current drive**, and that is exactly what bistability means. Fixed so that it always
    starts on the upper branch.

---

## 6. Model structure

### 66 compartments

| Group | Compartments |
|------|------|
| PK (22) | rivastigmine oral/patch/central, donepezil absorption/central/peripheral/effect, memantine, pimavanserin + AC-279, levodopa absorption/central/brain, antipsychotic, zonisamide, ambroxol absorption/central/brain |
| Antibody (3) | central · peripheral · CSF |
| Enzyme (2) | carbamylated AChE · carbamylated BuChE |
| Proteostasis (10) | GCase, GlcCer, α-syn monomer, regional oligomers ×3, regional fibrils ×3, interstitial seed pool |
| Neuronal populations (9) | nbM, PPN/LDT, SNc, locus coeruleus, dorsal raphe, REM atonia circuit, orexin, cardiac sympathetic, cortical synaptic density |
| Transmitters · receptors (6) | ACh, striatal DA, cortical NE, M1/M4, D2, 5-HT2A |
| Inflammation · co-pathology (5) | microglial activation, microglial exhaustion, astrocytes, amyloid, tau |
| Clinical states (7) | cognition, attention state, fluctuation amplitude, hallucination burden, MDS-UPDRS III, RBD, autonomic, somnolence, cumulative hazard |

### 23 scenarios

The controls are laid out in pairs — so that **every claim reads as a difference**.

| # | Scenario | Paired control |
|---|---|---|
| 01–03 | Natural history: DLB · PD dementia · AD | — |
| 04–07 | Rivastigmine oral / patch / donepezil / + memantine | 01 |
| 08 | Donepezil (AD phenotype) | 03 — phenotype dependence of the ChEI response |
| 09 | Donepezil + anticholinergic | 06 |
| 10–11 | Pimavanserin / quetiapine | 01 |
| 12–14 | Risperidone 1 mg: DLB / PD dementia / AD | one another |
| 15–16 | Levodopa 150 mg TID: DLB / PD dementia | 01 / 02 |
| 17 | Zonisamide | 01 |
| 18–20, 20b | GBA1 carrier untreated / ambroxol early / delayed / **stopped** | one another |
| 21 | Anti-α-synuclein antibody | 01 |
| 22 | Optimal combination (patch + pimavanserin + melatonin + droxidopa) | 01 |

### Mechanism verification functions

```r
source("dlb_mrgsolve_model.R")

neuroleptic_demo()   # same exposure, three receptor landscapes (§2-1 table)
levodopa_demo()      # same exposure, two phenotypes
responder_rates()    # responder rate in a virtual patient population (uses inter-individual variability)
ambroxol_demo()      # the GBA1 loop: early vs delayed vs stopped
chei_ratio()         # fluctuation/mean improvement ratio — the model's falsifiable prediction
run_all()            # all 22 scenarios
summarise_scn(run_all())
```

---

## 7. How to run

```bash
# render the mechanistic map
dot -Tsvg dlb_qsp_model.dot -o dlb_qsp_model.svg
dot -Tpng -Gdpi=150 dlb_qsp_model.dot -o dlb_qsp_model.png
```

```r
# the ODE model
install.packages(c("mrgsolve", "dplyr", "tidyr", "ggplot2", "shiny"))
source("dlb_mrgsolve_model.R")
out <- mrgsim_df(zero_re(mod), end = 2920, delta = 1)

# the Shiny dashboard
shiny::runApp("dlb_shiny_app_en.R")
```

`mrgsolve` requires a C++ compiler (`r-base-dev` on Linux, Xcode CLT on macOS, Rtools on
Windows).

> **Use `zero_re()`.** `$OMEGA` is defined, so calling `mrgsim()` plainly draws different
> random effects every time and the deterministic comparisons will not reproduce. `run_scn()`
> already uses `zero_re()` internally.

---

## 8. Limitations — what the model does not support

1. **The GBA1 feedback is not bistable.** §2(4). The first draft's claim of a "deadline for
   disease modification" has been retracted.
2. **The functional form of `FLUCTGT` was not fitted.** That cognitive fluctuation is the
   variance of a bistable attention state is an **interpretation** consistent with the
   literature, and the specific functional form is a choice. The 3.49-fold prediction it
   produces is the falsifiable content of that choice.
3. **`KSUPP` (the slope by which limbic α-syn suppresses postsynaptic striatal integrity)** is
   not a measured value. It is a single parameter back-calculated to match the incidence of
   neuroleptic sensitivity.
4. **All the neuronal-loss rate constants are correlated with one another.** They were tuned
   to hit six natural-history anchors simultaneously, so no biological meaning should be
   attached to the individual values.
5. **Rivastigmine's ΔMMSE of +3.9 is quite likely an overprediction.** The primary endpoint of
   the DLB rivastigmine trial (McKeith 2000) was NPI-4, and the model matches that well
   (−32.9% vs about 30%). But no ChEI raises MMSE by 4 points. The benefit of BuChE
   co-inhibition (`BCHEFR`/`KBCHUP`) appears to be overestimated, and with no head-to-head
   data against donepezil there is no basis on which to recalibrate.
6. **The time at which hallucinations cross the clinical threshold (26 months) is late.** In
   real DLB, hallucinations are already common around the time of diagnosis. Raising the
   initial neocortical fibril burden (`FIBN_0`) would fit it, but that throws other anchors
   off, so it is reported as it stands rather than fitted.
7. **The levodopa responder-rate contrast is a negative result.** See §4. The model separates
   DLB from PD dementia on mean change but not on responder rate, and it has no early
   Parkinson's disease arm.
8. **There is no prospective validation.** None of this model's predictions has been confirmed
   against new data.

---

## 9. Disclaimer

This model is a **semi-quantitative QSP model for education and research**, based on the
public literature. It has not been independently validated or qualified, and **must not be
used directly for real clinical decisions, prescribing, or regulatory submission.** The
parameters and assumptions are illustrative approximations.

Licence: see the repository [LICENSE](../LICENSE).
