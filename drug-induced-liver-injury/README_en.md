# Drug-induced liver injury (DILI) — QSP model
## Drug-Induced Liver Injury · Quantitative Systems Pharmacology

<p align="center">
  <a href="dili_qsp_model_en.svg">
    <img src="dili_qsp_model_en.png" width="900" alt="DILI QSP mechanistic map">
  </a><br>
  <sub><a href="dili_qsp_model_en.svg">View the full-resolution SVG</a> · 222 nodes · 330 edges · 15 clusters (14 mechanisms + 1 legend)</sub>
</p>

---

## The claim this model is built to make

> **Liver injury is not a problem of dose but of rate, and because of the JNK–Sab
> positive feedback loop the system is bistable.**

As with the other models in this repository, what matters here is not "what was put in"
but **"what came out although it was not put in"**. The following do **not exist** in
this model's parameter list:

- a dose threshold parameter
- an antidote time-window constant such as an "8-hour rule"
- injury-type switches such as `pattern` / `cholestatic` / `hepatocellular`
- a severity scale

And yet all the results below emerge. Because what went into the model is only **a
competition of rates** and **one positive feedback loop**.

| What was put into the model (structure) | What the model computed (result) |
|---|---|
| Reactive metabolite formation flux vs GSH resynthesis flux (cysteine rate-limited) | A steep dose threshold between 280 and 290 mg/kg |
| p-JNK → Sab → ROS → p-JNK positive feedback (saturating gain) | Bistability, and the separatrix between "adaptation" and "necrosis" |
| Cofactor (PAPS · UDPGA) depletion | The **fraction** bioactivated rises by itself at high dose |
| A single drug parameter competitively inhibiting BSEP | A cholestatic R ratio of 0.65 (same code, same equations) |
| Treg / PD-1 tolerance terms | The same drug at the same dose causes hepatitis only in carriers |
| ALT is a release rate; bilirubin is a clearance capacity | Hy's Law is a joint test of **rate × reserve** |

Every number was actually computed by
[`dili_reference_check_en.py`](dili_reference_check_en.py), and
the full output is in
[`dili_reference_output_en.txt`](dili_reference_output_en.txt).
No number was written in by hand.

---

## Result 1 — the threshold is not a parameter but a computed outcome

A single acute ingestion in a normal host. There are only mechanistic parameters such as
`KMD_ADD` and `KJ_K`; no value called "threshold dose" exists anywhere in the code.

| Dose (mg/kg) | GSH nadir | JNK peak | Peak ALT | Hepatocytes lost | Total bilirubin | Hy's Law |
|---:|---:|---:|---:|---:|---:|:--:|
| 100 | 84.2% | 0.003 | 25 | 0.0% | 0.71 | no |
| 160 | 72.6% | 0.010 | 25 | 0.0% | 0.82 | no |
| 220 | 60.3% | 0.030 | 26 | 0.0% | 0.97 | no |
| 250 | 54.0% | 0.050 | 35 | 0.1% | 1.08 | no |
| **280** | 47.6% | 0.181 | **509** | **4.1%** | 1.31 | no |
| **310** | 41.3% | 0.384 | **4 150** | **37.5%** | 1.84 | no |
| 340 | 34.9% | 0.507 | 6 099 | 59.6% | 2.34 | no |
| 370 | 28.6% | 0.591 | 7 070 | 72.9% | 2.87 | **YES** |
| 400 | 22.2% | 0.646 | 7 596 | 81.3% | 3.41 | **YES** |

- The 5% mass loss point is crossed between **280 and 290 mg/kg**.
- There is a region (300 mg/kg) where a single 10 mg/kg increment destroys a further
  **12.1 percentage points** of the liver. This is the fingerprint of a bistable system.
- Clinical anchors: 150 mg/kg = the (conservative) treatment line, above 250 mg/kg = high
  risk (Rumack 2002). Because the model's "normal host" is an idealised individual with an
  intact GSH pool and no CYP2E1 induction, it is consistent for the threshold to sit above
  the clinical treatment line — the reason the treatment line is low is that **hosts
  vary**, and that point is quantified in result 5.

---

## Result 2 — the same total dose, a different rate (the model's central claim)

**Fixing the total at 350 mg/kg** (24.5 g in a 70 kg adult) and varying only the rate at
which it arrives. If liver injury were a problem of dose, every row below would be
identical.

| Ingestion time | Plasma AUC (µM·h) | GSH nadir | JNK peak | Peak ALT | Hepatocytes lost |
|---|---:|---:|---:|---:|---:|
| Single bolus | 21 118 | 32.8% | 0.539 | 6 495 | **64.8%** |
| 3 h | 20 746 | 33.8% | 0.525 | 6 327 | 62.5% |
| 6 h | 20 221 | 35.3% | 0.501 | 6 028 | 58.7% |
| 12 h | 19 189 | 38.9% | 0.433 | 5 032 | 47.0% |
| 18 h | 18 241 | 43.4% | 0.337 | 3 197 | 28.0% |
| 24 h | 17 392 | 48.2% | 0.177 | 477 | 3.9% |
| 36 h | 15 996 | 58.3% | 0.039 | 29 | **0.0%** |
| 48 h | 14 959 | 67.3% | 0.021 | 25 | **0.0%** |

**The same 24.5 g either destroys two thirds of the liver or is metabolised without a
trace.** The exposure (AUC) differs only 1.4-fold and the outcome is 0% against 65%.
The reason is simple — GSH resynthesis capacity is a finite **rate** (about 0.6 mM/h), and
reactive metabolite arriving more slowly than that is disposed of indefinitely, while
metabolite arriving fast depletes the pool and ignites the positive feedback loop.

---

## Result 3 — the antidote window is not a constant but the moment the trajectory crosses the separatrix

Nowhere in the model is there an "8 hours" or a "10 hours". NAC only increases the cysteine
supply (the `CYS` state variable) and scavenges ROS directly.

**APAP 350 mg/kg, normal host** (the Prescott 21-hour IV regimen)

| NAC start (h) | 2 | 4 | 6 | 8 | 10 | 12 | 14 | 16 | 20 | 24 | 32 | 48 | untreated |
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| Hepatocytes lost (%) | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.6 | 8.4 | 18.3 | 33.9 | 46.3 | 58.2 | 64.6 | 64.8 |
| Peak ALT | 26 | 25 | 25 | 26 | 30 | 88 | 981 | 2 055 | 3 628 | 4 875 | 6 035 | 6 495 | 6 495 |

Administration within 10 hours is **completely protective** (normal ALT, 0% lost); it
starts to waver at 12 hours and collapses abruptly between 14 and 24 hours. This
**reproduces without putting it in as a parameter** the observation of Smilkstein 1988
(NEJM) — hepatotoxicity near 0% when given within 10 hours, rising sharply at 16–24 hours.

And the window **moves**:

| Situation | When half the available benefit has disappeared |
|---|---|
| 350 mg/kg, normal host | between 16 and 20 h |
| 500 mg/kg, normal host | between 12 and 14 h |
| 150 mg/kg, chronic alcohol + fasted host | **between 8 and 10 h** |

---

## Result 4 — the 150 mg/kg treatment line is a line about the host, not about the dose

**Exactly the same 150 mg/kg ingestion**, different hosts:

| Host | GSH nadir | JNK peak | Peak ALT | Total bilirubin | INR | Hepatocytes lost | Hy's Law |
|---|---:|---:|---:|---:|---:|---:|:--:|
| Normal | 74.6% | 0.008 | 25 | 0.80 | 1.09 | **0.0%** | no |
| Chronic alcohol + fasted<br><sub>(CYP2E1 induced 2.2-fold, reduced cysteine supply, GSH pool 60%)</sub> | 11.9% | 0.702 | 7 911 | 3.86 | 6.03 | **86.7%** | **YES** |
| The same host + NAC at 8 h | 19.1% | 0.283 | 3 240 | 2.06 | 2.02 | 30.0% | no |

A dose that does nothing at all to a normal host causes fulminant hepatic failure in a
susceptible one. This is where the reason the clinical nomogram's treatment line is drawn
far below the actual toxic threshold is explained by computation.

---

## Result 5 — in idiosyncratic DILI the variable is not dose but tolerance

Drug C (strong bioactivation, slow adduct removal, no BSEP liability) given at 10 mg/kg/day
for 56 days. **The drug, the dose and the adduct burden are all completely identical**, and
the only difference is whether the immune tolerance axis (Treg / PD-1) is intact.

| Scenario | Peak ALT | Time of onset | Peak drug-specific T cells | Treg nadir | Hepatocytes lost | Total bilirubin |
|---|---:|---:|---:|---:|---:|---:|
| **Non-carrier** of the HLA risk allele | **25** | — | 0.000 | 1.000 | **0.0%** | 1.00 |
| HLA carrier | 725 | **9.1 weeks** | 0.200 | 1.000 | 7.4% | 1.30 |
| HLA carrier + checkpoint inhibitor | **2 796** | **4.2 weeks** | 0.494 | 0.012 | 40.2% | 2.37 |
| The above + steroid (from day 42) | 2 796 | 4.2 weeks | 0.478 | 0.000 | **36.8%** | 2.21 |

- **The non-carrier is completely asymptomatic.** Adducts are formed (a subclinical adduct
  burden) but there is no clinical injury — this is the model's representation of why
  idiosyncratic DILI is rare.
- Remove tolerance with a checkpoint inhibitor and the same latent burden converts into
  hepatitis that is **faster (4.2 weeks against 9.1) and more severe**. This reproduces the
  experiment of Metushi–Uetrecht (2015) that "unmasked" amodiaquine liver injury in
  PD-1⁻/⁻ + anti-CTLA-4 mice.
Sweeping the steroid start day produces **exactly the same cliff as NAC** (drug, dose and
tolerance state all identical, only the start day differing):

| Steroid start | none | day 14 | day 21 | day 28 | day 35 | day 42 | day 49 |
|---|--:|--:|--:|--:|--:|--:|--:|
| Hepatocytes lost (%) | 40.2 | **0.3** | **1.9** | 27.2 | 34.1 | 36.8 | 38.8 |
| Peak ALT | 2 796 | 57 | 214 | 2 570 | 2 796 | 2 796 | 2 796 |

**There is a cliff between day 21 and day 28.** Before it, protection is nearly complete;
after it, nearly meaningless. The mechanism is identical to the NAC window — once the injury
has mobilised innate immunity and ROS amplification, suppressing T cells does not stop a
loop that is already turning. Nowhere in the model is there a term saying "treat early".

---

## Result 6 — why a threshold exists: bistability of the JNK–Sab loop

Isolating the loop alone and finding the fixed points against the redox headroom
`g = (0.30 + 0.70 · GSH/GSH₀) · NRF2`:

| g | State |
|---:|---|
| 1.60 – 0.90 | **monostable OFF** (JNK → 0) — injury impossible |
| 0.85 | **BISTABLE** — separatrix JNK = 0.031, ON state 0.391 |
| 0.80 and below | **monostable ON** (JNK → 0.536 and above) — injury inevitable |

The critical headroom is around **g\* ≈ 0.87**. GSH consumption **lowers** g and Nrf2
induction **raises** it. Because the two forces fight each other during the run, the
separatrix does not sit in a fixed place, and that is why the antidote window in result 4 is
not a constant.

> **An honest caveat.** The bifurcation diagram above is for *the isolated loop*, an
> approximation that holds NRF2 fixed and ignores the additional ROS input arriving from
> mitochondria and Kupffer cells. In the full model the transition happens earlier because
> of that additional input. The bifurcation diagram is a tool for explaining the *cause* of
> the threshold, not a formula for predicting the threshold value.

---

## Result 7 — Hy's Law is a joint test of rate × reserve

Why does Hy's Law (ALT ≥ 3×ULN **and** total bilirubin ≥ 2×ULN, without biliary
obstruction) require both arms at once? In the model the two measures are **structurally
different quantities**.

- `dALT/dt = k_rel · (necrosis flux) − k_el · ALT` → ALT is a low-pass filter of **the rate
  of dying** (half-life 47 h).
- `dTBIL/dt = production − k_clr · (surviving hepatocyte mass) · TBIL` → bilirubin is
  determined by **the reserve remaining** (clearance being proportional to surviving mass).

Sweeping the dose to see when each arm is crossed:

| Arm | Dose at which it is first crossed | Hepatocyte mass lost at that point |
|---|---:|---:|
| ALT ≥ 3×ULN (120 U/L) | 280 mg/kg | **4.15%** |
| Total bilirubin ≥ 2×ULN (2.4 mg/dL) | 350 mg/kg | **64.75%** |

**The bilirubin arm requires 15.6 times as much liver to be lost as the ALT arm.** So the
conjunction of the two arms becomes a test requiring simultaneously that "the liver is dying
fast (rate)" **and** that "the reserve is exhausted (reserve)", and that is why neither arm
alone predicts death.

> Of the 51 doses simulated, 33 crossed the ALT arm and 26 crossed both (79%). This
> proportion **depends entirely on the dose distribution sampled** (100–600 mg/kg was swept
> uniformly, so severe doses are over-represented), so it must not be compared directly with
> the clinical "Hy's Law positive → about 10% mortality". What matters here is not the
> proportion but the structural fact that **the two arms measure different physical
> quantities**.

---

## Result 8 — the injury pattern (the R ratio) is not specified but emergent

The three runs differ only in `KI_BSEP` (the BSEP inhibition Ki), `VMAX_CYP` (the
bioactivation capacity) and `KP` (the hepatic partition coefficient). The equations are
completely identical.

| Drug | Peak ALT | Peak ALP | R ratio | Verdict | Serum bile acids | Hepatocytes lost |
|---|---:|---:|---:|---|---:|---:|
| A — APAP 350 mg/kg (bioactivation type) | 6 495 | 70 | **278.4** | hepatocellular | 8.3 µM | 64.8% |
| B — BSEP inhibitor for 28 days (Ki 0.5 µM) | 141 | 648 | **0.65** | cholestatic | 39.1 µM | 1.0% |
| A+B — dual liability | 2 352 | 278 | **25.4** | hepatocellular (tending to mixed) | 37.7 µM | 25.0% |

Because BSEP is an ABC transporter, its flux is scaled by `× ATP` in the model. So **the
fact that hepatocellular injury produces secondary cholestasis** also follows without a
separate term (bile acids 37.7 µM in the A+B row).

---

## Result 9 — the reason miR-122 rises before ALT is nothing but its half-life

APAP 350 mg/kg, untreated:

| Measure | Time at which the threshold is reached |
|---|---:|
| miR-122 ≥ 5× baseline | **12.5 h** |
| AST ≥ 120 U/L | 15.5 h |
| ALT ≥ 3×ULN | 16.0 h |
| Total bilirubin ≥ 2×ULN | 48.2 h |

miR-122 leads ALT by **3.5 hours**. The model has no separate "early biomarker" term — both
measures are released from **the same necrosis flux**, and the only difference is the
elimination half-life (miR-122 3 h against ALT 47 h). A short half-life just is a fast rise.

---

## Result 10 — a refutation test: does the JNK–Sab loop really hold the structure up?

Sab (SH3BP5)-deficient mice are protected from APAP hepatotoxicity (Win 2011). Setting
`KROS_SAB = 0` in the model cuts only the positive feedback loop and leaves everything else
in place — bioactivation, GSH depletion and adduct formation all happen identically.

| Dose (mg/kg) | Loop intact: lost (JNK peak) | **Loop cut**: lost (JNK peak) |
|---:|---:|---:|
| 250 | 0.1% (0.050) | 0.0% (0.035) |
| **310** | **37.5%** (0.384) | **1.7%** (0.115) |
| 350 | 64.8% (0.539) | 19.9% (0.248) |
| 400 | 81.3% (0.646) | 46.5% (0.347) |
| 500 | 92.1% (0.729) | 75.9% (0.512) |
| 700 | 96.7% (0.758) | 89.3% (0.609) |

At 310 mg/kg the injury **disappears from 37.5% to 1.7%** and the threshold is pushed far
to the right. The loop is not decoration but **the device that converts metabolic damage
into necrosis**, and this agrees with what the knockout experiment says. (At sufficiently
large doses injury remains even without the loop, because the adducts themselves destroy
mitochondria — the loop is the device that creates the threshold, not the only killing
route.)

---

## Honest limitations

1. **Non-monotonicity of the NAC time window.** At 500 mg/kg and in the susceptible host,
   starting NAC at 2 hours is *slightly worse* than starting at 6 hours. Two hypotheses were
   tested and one was rejected:
   - "the fixed-duration regimen ends before the metabolism is finished" → extending the
     maintenance infusion **at the same rate** reduces the loss greatly, from 49.3% to
     29.3%, but the ordering still does not invert. So coverage duration alone does not
     explain it.
   - "the 2-hour run simply takes more damage" → **rejected.** The 6-hour run actually has
     the lower GSH nadir (23.3% against 26.2%) and more adducts (23.5 against 21.0). And yet
     its injury is less.

   So this is a dynamic effect of **when the thiol supply arrives relative to the reactive
   metabolite pulse**, not of the total quantity of damage. It is reported as it stands
   rather than smoothed away, and nobody should believe it without experimental
   verification. The clinically usual 350 mg/kg curve is completely monotonic and "NAC as
   early as possible" remains the model's unvarying recommendation.

2. **The normal host's threshold (≈285 mg/kg) is higher than the clinical treatment line
   (150 mg/kg).** This is because the model's "normal host" is a single individual whose GSH
   pool, cysteine supply and CYP2E1 are all ideal. Real population variability is reproduced
   only by putting `FCYP` / `CYSBASE` / `FGSH` / `AGEF` as distributions, and the models in
   this repository are deterministic single-individual simulations.

3. **The bistability bifurcation diagram is an approximation for the isolated loop** (see the caveat in result 6).

4. **Drugs B and C are hypothetical drugs.** Each is a parameter set modelled on the
   *archetype* of troglitazone/bosentan and flucloxacillin/amoxicillin-clavulanate
   respectively, and neither is a validated model of a specific real drug. Molecular weight
   conversion uses APAP (151.16) for convenience.

5. **Calibration status.** The parameters were tuned by hand to the physiological ranges of
   the literature and the clinical anchors below; no formal fitting to patient data or
   uncertainty quantification was performed.

---

## Calibration anchors

| Anchor | Source | Model result |
|---|---|---|
| At therapeutic doses glucuronide:sulfate:CYP ≈ 55:30:6, total clearance ≈ 19–21 L/h | Prescott 1980 | structurally matched (Vmax/Km allocation) |
| Sulfation saturates first at high dose | Slattery 1987 | reproduced by the PAPS/UDPGA depletion ODEs |
| Hepatic GSH 5–10 mM, cysteine the rate-limiting substrate | Lu 1999/2013 | GSH₀ = 6.5 mM, NAC acts only through `CYS` |
| Apoptosis with ATP, necrosis without | Kon 2004 | the ATP gate on the MPT term |
| Injury abolished in Sab deficiency | Win 2011 | with `KROS_SAB = 0` the loop is cut and the threshold disappears |
| NAC within 10 hours → hepatotoxicity near 0%, sharp rise at 16–24 hours | Smilkstein 1988 | reproduced in result 3 |
| Drugs with a low BSEP IC50 = cholestatic risk | Morgan 2010 · Dawson 2012 | result 7 |
| ALT half-life 47 h, AST 17 h | Ozer 2008 | `KALT_EL`, `KAST_EL` |
| Hy's Law | Temple 2006 | result 6 |

For the full list of 108 references and which structure or parameter each one determined,
see [`dili_references_en.md`](dili_references_en.md).

---

## Files

| File | Contents |
|---|---|
| [`dili_qsp_model_en.dot`](dili_qsp_model_en.dot) | Mechanistic map source (222 nodes · 330 edges · 15 clusters (14 mechanisms + 1 legend)) |
| [`dili_qsp_model_en.svg`](dili_qsp_model_en.svg) / [`.png`](dili_qsp_model_en.png) | The rendered map (`dot -Tsvg` / `dot -Tpng -Gdpi=150`) |
| [`dili_mrgsolve_model.R`](dili_mrgsolve_model.R) | 33-ODE mrgsolve model + 13 scenarios + 5 analysis functions |
| [`dili_shiny_app_en.R`](dili_shiny_app_en.R) | 10-tab interactive dashboard |
| [`dili_references_en.md`](dili_references_en.md) | 108 annotated references (16 sections) |
| [`dili_reference_check_en.py`](dili_reference_check_en.py) | A Python/scipy re-implementation with equations and parameters identical to the R model |
| [`dili_reference_output_en.txt`](dili_reference_output_en.txt) | The full output of the above script = the source of every number in this README |

### The 15 clusters of the map

1. Pharmacokinetics · hepatic distribution · 2. Biotransformation (activation vs
detoxification) · 3. Glutathione · thiol defence ·
4. Covalent adducts · mitochondrial bioenergetics · 5. Oxidative stress · Nrf2 adaptation ·
6. **The JNK–Sab amplification loop and the cell death decision** · 7. Bile acid
homeostasis · cholestasis ·
8. Innate immunity (DAMPs · Kupffer cells) · 9. Adaptive immunity (HLA-restricted
idiosyncrasy) ·
10. Regeneration · reserve · 11. Host susceptibility factors · 12. Biomarkers ·
13. Clinical endpoints · causality · prognosis · 14. Treatment · clinical management · 15. Legend

---

## How to run

```bash
# 1) render the mechanistic map
dot -Tsvg dili_qsp_model_en.dot -o dili_qsp_model_en.svg
dot -Tpng -Gdpi=150 dili_qsp_model_en.dot -o dili_qsp_model_en.png
#    note: newrank=true is required. There are so many edges crossing clusters that
#    graphviz's older per-cluster ranker fails with "trouble in init_rank".

# 2) re-verify the model (possible without R, about 2 minutes)
pip install numpy scipy
python3 dili_reference_check_en.py        # regenerates dili_reference_output_en.txt

# 3) simulate with mrgsolve
R -e "source('dili_mrgsolve_model.R'); print(run_all())"
R -e "source('dili_mrgsolve_model.R'); print(analysis_nac_window())"

# 4) Shiny dashboard
R -e "shiny::runApp('dili_shiny_app_en.R', port=8080)"
```

`dili_reference_check_en.py` and `dili_mrgsolve_model.R` contain **the same equations and the
same parameters**. If you change one, you must change the other and re-run the Python
verification.

---

## ⚠️ Disclaimer

This model is a **qualitative / semi-quantitative QSP model for educational and research
purposes**. It was constructed from the public literature but has not been independently
validated or certified, and **must not be used for real clinical decision-making,
prescribing, poisoning management, or regulatory submission.** All real poisonings,
including paracetamol overdose, must be managed according to the nomograms, antidote
protocols and poisons-centre guidance. The parameters are illustrative approximations and
fitting and validation against real patient data are separately required.
