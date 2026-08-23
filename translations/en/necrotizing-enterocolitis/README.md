# Necrotising Enterocolitis of the Preterm Infant (NEC) — QSP Model

**Not written as inflammation that gets worse — written as a closed loop hanging on one variable, enterocyte and villous integrity.**

---

## 0. The one structural claim of this model

NEC is a **closed positive-feedback loop** on a single state variable — enterocyte
and villous integrity `E` — whose loop gain is set by the luminal pathobiont
load `B`:

```
   E ──▶ BI = E·(TJ)·(MUC)                       barrier integrity
     ──▶ Pb = Pmin + (Pmax−Pmin)·(1−BI)^hP       permeability
     ──▶ Jtr = Pb · B                            bacterial translocation flux
     ──▶ TLR4s = TLR4expr · Jtr/(Jtr + Ktlr)     innate immune signal
     ──▶ (a) cell loss  kA·INJ·E/(E+KEap)   rises
         (b) crypt proliferation  1/(1+(INJ/Ki)²)   falls
     ──▶ E falls  ... and round again
```

Both injury arms (a) and (b) act on `E`, and that `E` is in turn what sets `Pb`.
That is the loop. Reducing to one dimension (`TJ = MUC = 1`, `B` held) gives the
field

```
   g(E) = kE(1−E)·Ftroph/(1+(INJ/Ki)²) − kA·INJ·E/(E+KEap),   INJ = INJ(Jtr(E;B))
```

Because `g(0) > 0` and `g(1) < 0`, the number of roots in `(0,1)` is **always
odd**. At low `B` there is only one root (health); as `B` rises a saddle-node
bifurcation creates **a pair** (the separatrix `E*` plus a necrotic attractor);
higher still and the healthy root is annihilated.

Five results then follow as **arithmetic rather than assertion**.

> **Where every number comes from.** The build environment of this repository has
> no R runtime, so `nec_mrgsolve_model.R` cannot be executed. Committing an ODE
> model that has never been integrated is not honest, so every equation of the
> mrgsolve file was reimplemented term by term in
> `nec_reference_model.py` (pure Python, RK4) and **actually integrated**. Every
> table below is transcribed as it stands from that run log,
> `nec_reference_output.txt`.

---

## 1. There are two critical loads, not one

Counting the roots of `g(E)` while `B` is raised yields two saddle-node
bifurcation points.

* **`B_lo`** = the load at which bistability appears → NEC becomes **possible**
* **`B_hi`** = the load at which the healthy root is annihilated → NEC becomes **inevitable**

The colonisation niche has a ceiling of `Ktot = 26 ×10⁹ CFU/g`, so where
`B_lo > Ktot` NEC is **physically unreachable** in that infant.

| GA (wk) | Feed | Ftroph | Ktlr | TLR4expr | **B_lo** | **B_hi** | Window ratio | Reachable? |
|---|---|---|---|---|---|---|---|---|
| 25 | Formula | 0.948 | 60.1 | 0.821 | **11.3** | **29.8** | 2.6× | Yes |
| 25 | Donor milk | 0.991 | 69.7 | 0.821 | 13.5 | 36.3 | 2.7× | Yes |
| 25 | Own mother's milk | 1.121 | 73.3 | 0.821 | 15.7 | 43.8 | 2.8× | Yes |
| 28 | Formula | 0.948 | 60.2 | 0.739 | **14.9** | **38.1** | 2.6× | Yes |
| 28 | Donor milk | 0.991 | 69.7 | 0.739 | 17.9 | 46.4 | 2.6× | Yes |
| 28 | Own mother's milk | 1.121 | 73.3 | 0.739 | 20.8 | 56.0 | 2.7× | Yes |
| 32 | Formula | 0.948 | 60.3 | 0.638 | 20.9 | 50.2 | 2.4× | Yes |
| 32 | Donor milk | 0.991 | 69.8 | 0.638 | 25.2 | 61.2 | 2.4× | Yes |
| 32 | **Own mother's milk** | 1.121 | 73.3 | 0.638 | **29.5** | 74.1 | 2.5× | **No (> Ktot)** |

**The reading.** When the load lies between `B_lo` and `B_hi`, what decides the
outcome is **not the load but the initial condition**. The infant sits on the
healthy branch and needs **something to push `E` below `E*`** — a precipitating
insult — before it crosses over. That is why NEC is all-or-none rather than a
disease that "gradually gets worse", why it clusters, and why two infants on the
same unit receiving the same feeds do different things.

This window is why the model's **insult term** actually does any work:
apnoea-bradycardia, hypotension, PDA diastolic flow reversal, packed red-cell
transfusion, late-onset sepsis.

---

## 2. The separatrix `E*` and the logarithmic dependence of the transition time

The roots of `g(E) = 0` at 28 weeks · formula · `B = 18`:

| Root | Value | Character |
|---|---|---|
| Necrotic attractor `E_low` | **0.0530** | stable |
| **Separatrix `E*`** | **0.2889** | **unstable (repelling)** |
| Healthy attractor | **0.8603** | stable |

The local eigenvalue at `E*` is `λ = +0.354 /d`. In a linear field the transition
time is `t = (1/λ)·ln(span/δ)`, and so **depends on δ only logarithmically**.

| Displacement δ from `E*` | ln-law prediction (h) | Actual integration (h) |
|---|---|---|
| 0.0005 | 417.7 | **465.2** |
| 0.0020 | 323.6 | 371.3 |
| 0.0080 | 229.6 | 277.6 |
| 0.0300 | 139.9 | 187.5 |
| 0.1000 | 58.2 | **101.1** |

**A 200-fold range in δ is compressed into a 4.60-fold range in transition
time** — the fingerprint of a logarithm. (The ≈ 48 h gap between prediction and
integration is the time spent traversing the remainder of the trajectory after
leaving the neighbourhood of `E*`, and it is constant, independent of δ. That is
the evidence that the linearisation is right.)

**A caveat that has to be written down honestly.** The `λ` above is a
**quasi-static** eigenvalue with `B` and `NECa` held fixed. The full 39-state
system restores the two amplifiers the reduction threw away — malabsorption
(`E↓ → SUB↑ → B↑`) and the DAMP term (`wNEC·NECa`). So it traverses the same
logarithm **in hours rather than in days**:

| Case | Bell II | Bell III | Interval |
|---|---|---|---|
| S5 (25 wk · formula · antibiotics 7 days + indomethacin) | d10.50 | d11.33 | **19.9 h** |
| S11 (27 wk · severe · no rescue intervention) | d9.63 | d10.54 | **21.7 h** |

This is the arithmetical reason why "Bell II to III is always about a day", and
the reason **why a biomarker that rings at Bell II buys so little time**
(see §7).

---

## 3. The point of no return is a function of lesion size, not of time

Because necrotic tissue is itself an inflammatory stimulus (DAMPs), the lesion
**sustains the injury it created**:

```
   INJ = TLR4s + wNO·NOx + wISCH·ISCH + wNEC·NECa
```

For `INJ` to exceed `INJth` on the `NECa` term alone requires

```
   NECa_crit = INJth / wNEC = 0.35 / 1.80 = 0.194
```

That is: **once about 19 % of the bowel is transmurally necrotic, the lesion
grows by itself from then on, even if every other source of injury disappears.**
Nil-by-mouth and broad-spectrum antibiotics can arrest the lesion **only below**
this line.

So in this model **what separates medical from surgical NEC is not severity but
whether the lesion crosses this line during the diagnosis-to-treatment delay
`dx_lag`**. Re-epithelialisation occurs only below `NECrev = 0.10` (mucosal and
submucosal injury heals, transmural necrosis does not), and intramural gas
appears only once the mucosa has already been breached
(`Pb > PbPNE = 0.20`) — that is, **pneumatosis is modelled as a consequence, not
a cause**.

---

## 4. A drug that moves the threshold and a drug that moves the state multiply

The model has two kinds of intervention, and they are not substitutes for one
another.

| Kind | Term it enters | Examples |
|---|---|---|
| **State-mover** | `B` only | Probiotics |
| **Threshold-mover** | `Ktlr` only | Oral 2′-FL, anti-TLR4, recombinant PAF-AH |

In one and the same population (n = 260, GA 24–32 wk weighted towards the lower
gestations, 3 days of empirical antibiotics):

| Arm | NEC % | RR | Surgery % | Median day confirmed |
|---|---|---|---|---|
| Formula (reference) | 25.4 | 1.000 | 1.9 | 20.3 |
| + state-mover only (probiotics) | 2.3 | **0.091** | 0.4 | 19.9 |
| + threshold-mover only (`thr_boost` 1.6) | 6.9 | **0.273** | 1.2 | 21.1 |
| **+ both** | 1.2 | **0.045** | 0.0 | 16.4 |
| Own mother's milk | 5.4 | 0.212 | 1.2 | 21.4 |

```
   multiplicative prediction   RR_t × RR_p          = 0.273 × 0.091 = 0.025
   additive prediction         1−(1−RR_t)−(1−RR_p)  = −0.636        ← inadmissible from the outset
   observed combination                             = 0.045
   → multiplicative error −0.021,  additive error −0.682
```

**Because one changes the graduation of the axis and the other moves the position
along that axis, the two effects cannot add — they multiply.** The additive model
predicts a negative risk and does not hold up in the first place.

**And human milk is both.** It raises `Ktlr` through HMOs and sIgA (threshold),
raises `Ftroph` through EGF/HB-EGF, and lowers `B` by supplying HMOs as a
**private substrate that only the commensals can use** (state). Human milk is
therefore not something to be set against a single-target drug — it is
**already combination therapy**, and the fair comparator is the "both" row of
the table above.

---

## 5. Feed advancement — where our hypothesis failed

Enteral substrate enters in three places with different signs: `Ftroph`
(+, trophic stimulus), `dB/dt` (+, bacterial substrate), and the oxygen
supply-demand balance (−, `demand = 1 + aDEM·FEED/200`). We therefore expected
an **interior optimum**.

**The integration does not support that expectation.** NEC risk is monotonically
increasing in the rate of advancement.

**An exclusively formula-fed cohort** (a deliberate stress test, and not the
standard of care of any real unit; n = 170 per rate):

| Advancement rate (mL/kg/d²) | GA 25 wk NEC % | GA 25 wk surgery % | GA 30 wk NEC % | GA 30 wk surgery % | Full feeds reached (d) | Weight gain (g) |
|---|---|---|---|---|---|---|
| 0 (minimal enteral feeds maintained) | **0.0** | 0.0 | **0.0** | 0.0 | not reached | 210 / 217 |
| 10 | 37.6 | 5.3 | 10.0 | 0.0 | 16 | 213 / 244 |
| 20 | 47.6 | 6.5 | 17.6 | 0.6 | 9 | 208 / 235 |
| 30 | 55.9 | 8.2 | 22.4 | 2.9 | 7 | 198 / 220 |
| 40 | 55.9 | 14.7 | 18.8 | — | 6 | 184 / 252 |
| **Spread across rates** | **55.9 pp** | | **22.4 pp** | | | |

Two things have to be read together.

1. **The same equation is shallow at 30 weeks and steep at 25** (spread 22.4 pp
   vs 55.9 pp). The null result of the large feed-advancement trial (SIFT) and
   clinical wariness about the most premature infants are not in contradiction:
   they are **one equation evaluated at two gestational ages**.
2. **But on NEC alone, the slowest advancement is always the lowest.** The price
   of slow advancement appears not in NEC but in **days to full feeds · TPN
   exposure · cholestasis · growth** (the last two columns, and the cholestatic
   bilirubin of §8). The trade-off lies **between endpoints, not within a single
   endpoint**.

The place where the hypothesis was refuted is left standing rather than erased.

---

## 6. Empirical antibiotics change sign at a computable duration

**The mechanism is not resistance but luminal delivery.** What changes the
microbial community is not the plasma concentration but the **luminal
concentration**.

| Drug | Luminal delivered fraction `f_lum` | Principal target | Direction |
|---|---|---|---|
| Ampicillin | **0.45** (biliary and renal excretion) | commensal anaerobes `C` | **harmful** |
| Gentamicin (intravenous) | **0.05** (essentially no oral bioavailability) | pathobionts `B` | beneficial |
| Metronidazole | **0.80** (almost free distribution) | commensal anaerobes `C` | **harmful** |

That is, empirical ampicillin + gentamicin **scorches the anaerobes alone and
leaves the Enterobacteriaceae comparatively intact**. Lose `C` and the SCFAs
disappear with it (and hence the pH suppression of `B`), and the niche opens.
**The benefit is immediate and the harm is delayed**, so the two integrals must
cross somewhere.

We run it deliberately in milk-fed infants — since the harm consists of the loss
of `C`, it can only appear where `C` is present (n = 200 per duration, GA 26 wk):

| Empirical antibiotics (days) | NEC % | Surgery % | `B` (d21) | `C` (d21) | SCFA (d21) |
|---|---|---|---|---|---|
| 0 | 6.0 | 0.5 | 12.22 | **4.33** | 4.0 |
| 2 | 5.5 | 1.0 | 12.05 | 4.42 | 4.1 |
| 3 | 8.0 | 1.5 | 12.11 | 4.37 | 4.0 |
| 5 | 5.5 | 0.5 | 12.22 | 4.26 | 3.9 |
| **7** | **9.5** | 0.5 | 13.02 | 3.62 | 3.3 |
| **10** | **11.0** | 1.5 | 14.48 | 2.13 | 1.9 |
| **14** | 9.0 | 2.0 | **15.17** | **0.74** | **0.6** |

The ecological indices (`C` 4.33 → 0.74, `B` 12.22 → 15.17, SCFA 4.0 → 0.6) are
monotonic and unambiguous in duration. NEC incidence is noisy because of the
insult draw (the 8.0 % of the 3-day arm), but beyond 7–10 days it decisively
exceeds the reference arm.

**And this generates one falsifiable prediction of the model:** in
**exclusively formula-fed infants, who have no commensals to lose, antibiotics
may even appear protective** (which is exactly what scenarios S3 and S10 of the
model produce). It can be tested by stratifying an antibiotic-exposure cohort by
feed type.

---

## 7. Most clinical biomarkers move behind the switch

Reference case S5 (25 wk · formula · antibiotics 7 days + indomethacin),
Bell II = d10.50:

| Biomarker | Day of crossing | Lead time relative to Bell II |
|---|---|---|
| **`Jtr > 3.0` (a mechanistic quantity)** | **d9.00** | **+35.9 h (leads)** |
| Intramural gas (pneumatosis) | d10.50 | 0.0 h (simultaneous — by definition) |
| Platelets < 100 | d11.54 | **−25.0 h (lags)** |
| CRP > 12 | never crosses | — |
| Lactate > 4 | never crosses | — |

What is to be seen here is not that the biomarkers rise well but that **most of
them rise after the switch has already been thrown**. What moves ahead of it are
the state variables (`E`, `Jtr`). Combined with the logarithmic transition time
of §2 the conclusion follows: **a biomarker that rings at Bell II can in
principle buy almost no time, and that is the structural reason why prevention
comes before early diagnosis.**

---

## 8. Scenarios (11) — the high-risk index patient

> These scenarios are **not a population average but a high-risk index
> patient**: an infant with dysbiosis who receives **one deep hypotensive insult
> on day 14**, the point at which the luminal load has fully accumulated. The
> median patient on the same regimen does not cross over at all (the population
> results of §9), and comparing arms in a patient who never even approaches the
> separatrix reveals nothing.

| Scenario | Bell | Day confirmed | Surgery | `E` (end) | `B` (end) | SCFA | Cholestatic bili | Weight gain g |
|---|---|---|---|---|---|---|---|---|
| S1 · 28 wk · own mother's milk | 0 | — | No | 0.921 | 0.02 | 14.6 | 0.00 | **389** |
| S2 · 28 wk · formula | **2** | 17.1 | No | 0.836 | 15.40 | 0.3 | 0.00 | 297 |
| S3 · 28 wk · formula + antibiotics 7 days | 2 | 17.7 | No | 0.836 | 15.41 | 0.2 | 0.00 | 303 |
| S4 · 28 wk · formula + probiotics | 1 | — | No | 0.906 | 0.05 | 11.0 | 0.00 | 384 |
| **S5 · 25 wk · formula + antibiotics 7 days + indomethacin** | **3** | **10.5** | **Yes** | **0.000** | 14.45 | 3.2 | 0.00 | **−130** |
| S6 · 25 wk · human milk + probiotics | 0 | — | No | 0.920 | 0.02 | 14.5 | 0.00 | 386 |
| S7 · 32 wk · formula | 1 | — | No | 0.865 | 13.13 | 2.6 | 0.00 | 306 |
| S8 · 28 wk · donor milk (pasteurised) | 1 | — | No | 0.916 | 0.02 | 14.6 | 0.00 | 385 |
| S9 · 28 wk · complete nil-by-mouth + TPN | 0 | — | No | 0.999 | 0.03 | 0.0 | **3.18** | 243 |
| S10 · 28 wk · formula + antibiotics 10 days + metronidazole | 2 | 18.9 | No | 0.836 | 15.44 | 0.2 | 0.00 | 304 |
| **S11 · 27 wk · severe (no rescue intervention)** | **3** | **9.6** | **Yes** | **0.000** | 6.64 | 8.5 | 0.00 | **−156** |

**What S9 says.** Complete nil-by-mouth is in this model the **safest** thing
there is for the barrier (`E` = 0.999). With no substrate `B` does not grow and
the loop does not close. The price appears elsewhere — cholestatic bilirubin
3.18 mg/dL (TPN cholestasis) and the lowest weight-gain group. The hypothesis we
first put up, that "fasting is dangerous by way of atrophy", is **not
supported** by this model, and that fact is left rather than erased.

### Extra-intestinal outcomes

| Scenario | Cholestatic bilirubin | Neuroinflammation index (drives NDI) | Tubular injury | Weight gain g |
|---|---|---|---|---|
| S1 · human milk | 0.00 | 0.000 | 0.000 | 389 |
| S2 · formula | 0.00 | 0.102 | 0.005 | 297 |
| S5 · 25 wk severe | 0.00 | **13.748** | 0.002 | −130 |
| S9 · nil-by-mouth + TPN | **3.18** | 0.000 | 0.000 | 243 |
| S10 · antibiotics 10 days + metronidazole | 0.00 | 0.102 | **0.008** | 304 |

---

## 9. Calibration — and where it does not fit

A population with a realistic feeding mix (60 % human milk, 55 % early empirical
antibiotics), n = 170 per stratum:

| GA stratum | Model NEC % | Model surgery % | Median day confirmed | Literature target NEC % |
|---|---|---|---|---|
| 24–25 wk | **21.8** | 3.5 | 17.3 | 12–15 |
| 26–27 wk | **12.4** | 1.8 | 17.4 | 8–11 |
| 28–29 wk | **2.9** | 0.0 | 20.7 | 4–6 |
| 30–32 wk | **0.6** | 0.0 | 14.8 | 1–3 |

**What fits.** The direction and the order of magnitude of the gestational-age
gradient, the timing of onset (2–3 weeks after birth), the fact that fulminant
deterioration runs in hours, and the fact that pneumatosis appears only after
mucosal destruction.

**What does not fit — stated explicitly.**

1. **The gradient is too steep.** At the lowest gestational age it is higher than
   the literature (21.8 % vs 12–15 %) and at the highest it is lower (0.6 % vs
   1–3 %). In a threshold system a small parameter difference is amplified into a
   large difference in incidence, so when the gestational-age-dependent terms
   (`phiTLR`, PAF-AH, IL-10, initial `E`·`TJ`·`MUC`) all act together in the same
   direction, the gradient is exaggerated.
2. **The surgical proportion is low.** About 15–25 % of confirmed NEC reaches
   Bell III, which is below the one-third to one-half of the literature.
3. **The effect size of human milk is stronger than in the literature.** The
   model gives RR ≈ 0.21 while the pooled estimate from trial meta-analyses is
   roughly 0.4–0.6. This model attributes a substantial part of the milk effect
   to the **ecological route** (colonisation by commensals capable of
   metabolising HMOs), and as a consequence predicts an interaction in which the
   effect divides according to whether `binfantis` is carried — this is the
   falsifiable prediction of §11, but it is also what makes the pooled RR
   overestimated.

---

## 10. Model structure (39 ODEs)

| Module | State variables |
|---|---|
| Mucosal barrier | `E` `TJ` `MUC` `IgA` `DEF` |
| Luminal ecology | `B` `C` `P` `SUB` `SCFA` `GAS` `HMO` |
| Innate immune signalling | `TLR4s` `IL1B` `TNF` `IL8` `IL10` `PAF` `NEUT` `NOx` |
| Perfusion · lesion | `PERF` `NECa` `PNEU` |
| Systemic · remote organs | `LPSp` `PLT` `CRP` `LAC` `NIN` `BILI` `WT` `KIN` |
| Pharmacokinetics | `AMPc` `GENc` `GENp` `MTZc` `INDc` `IBUc` `DEXc` `PGE` |

The core nonlinear structures:

* **`Pb = Pmin + (Pmax−Pmin)(1−BI)³`** — permeability responds to barrier loss cubically (`Pmax/Pmin = 24`)
* **`block = 1/(1+(INJ/Ki)²)`** — TLR4 signalling does not merely increase injury, it switches crypt regeneration off as well
* **`kA·INJ·E/(E+KEap)`** — loss saturates in `E` (a denuded mucosa goes on losing)
* **`alphaX = 0.55`** — cross-competition between the two bacterial guilds is below 1, so coexistence is possible (this is not single-niche competition)
* **`fpH = 1/(1+(SCFA/KpH)²)`** — SCFAs / low pH suppress the Enterobacteriaceae
* **`ISCH = max(0, 1 − (PERF/demand)/SDRcrit)`** — there is an oxygen reserve, and feeding eats into it
* **`wNEC·NECa`** — the lesion sustains its own injury → `NECa_crit = INJth/wNEC`
* **`f_lum`** — luminal delivery: the concentration the microbes see is not the plasma concentration
* **Secondary feedback: `E↓ → absorption↓ → SUB↑ → B↑ → E↓`** — why feed intolerance precedes NEC

---

## 11. The falsifiable predictions this model generates

These come not out of the literature but **out of the structure of the model**.
They are written down so as to be falsifiable.

1. **The protective effect of human milk divides according to whether a strain
   capable of metabolising HMOs is carried.** A study that measures HMO
   concentrations alone will see the effect diluted.
2. **The harm of empirical antibiotics should be observable only in milk-fed
   infants.** In exclusively formula-fed infants, who have no commensals to lose,
   antibiotics may appear protective. Stratifying an existing antibiotic-exposure
   cohort by feed type tests this.
3. **The point of no return is a function of lesion size and not a function of
   time** (`NECa_crit ≈ 0.19`). An intervention that shortens the
   diagnosis-to-treatment delay should therefore lower the surgical proportion
   **irrespective of severity**, and on a unit where the delay is already short
   the additional benefit should fall away sharply.
4. **A biomarker that rings at Bell II can in principle buy almost no time.** A
   useful biomarker has to measure not the stage but the **translocation flux**
   (intestinal permeability · bacterial translocation).
5. **Feed-advancement trials should divide by gestational age** — null at 30
   weeks and above, positive at 25 weeks and below. A pooled analysis will cancel
   the two effects against one another.

---

## 12. Files

| File | Contents |
|---|---|
| [`nec_qsp_model.dot`](nec_qsp_model.dot) | Mechanistic map source (167 nodes · 15 clusters) |
| [`nec_qsp_model.svg`](nec_qsp_model.svg) | Vector render |
| [`nec_qsp_model.png`](nec_qsp_model.png) | 150 dpi raster |
| [`nec_mrgsolve_model.R`](nec_mrgsolve_model.R) | mrgsolve 39-ODE model · scenarios · bifurcation analysis · population simulation |
| [`nec_reference_model.py`](nec_reference_model.py) | Pure-Python RK4 reimplementation of the same equations (for verification; actually executed) |
| [`nec_reference_output.txt`](nec_reference_output.txt) | The run log of the file above — **the source of every number in this README** |
| [`nec_scenario_results.json`](nec_scenario_results.json) | Bifurcation points · separatrix · scenarios · feeding · antibiotic results (machine-readable) |
| [`nec_population_results.json`](../../../necrotizing-enterocolitis/nec_population_results.json) | Factorial design · calibration results |
| [`nec_shiny_app.R`](nec_shiny_app.R) | 9-tab interactive dashboard |
| [`nec_references.md`](nec_references.md) | 106 references — each section marked with **which term of the model it supports** |

Reproduction:

```bash
# render the map
dot -Tsvg nec_qsp_model.dot -o nec_qsp_model.svg
dot -Tpng -Gdpi=150 nec_qsp_model.dot -o nec_qsp_model.png

# regenerate every number (pure standard library, about 8 minutes)
python3 nec_reference_model.py > nec_reference_output.txt
```

```r
# in an environment that has an R runtime
source("nec_mrgsolve_model.R")
sim <- run_all_scenarios();  summarise_scenarios(sim)
critical_load_table();       point_of_no_return()
factorial_arms();            abx_duration_sweep()
shiny::runApp("nec_shiny_app.R")
```

> `nec_mrgsolve_model.R` and `nec_reference_model.py` are **two implementations
> of the same equations**. Fix one and the other must be fixed too, otherwise the
> two files diverge.

---

## ⚠️ Disclaimer

This is a **qualitative to semi-quantitative QSP model intended for teaching and
hypothesis generation**. The parameters were hand-calibrated against **aggregate
indices** from the published literature (incidence by gestational age, the
distribution of onset times, risk ratios for human milk and probiotics, surgical
proportions, neonatal pharmacokinetics); they are **not** fitted to individual
patient data. The covariance structure was set arbitrarily. `B` and `C` are each
lumped into a single population, "Enterobacteriaceae-like" and
"obligate-anaerobe-like" respectively, whereas real NEC dysbiosis is determined
at strain level. The calibration mismatches set out in §9 must be read alongside
this.

**It must not be used for actual clinical decision-making, prescribing, or
regulatory submission.**
