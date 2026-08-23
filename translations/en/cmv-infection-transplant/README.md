# Post-transplant cytomegalovirus (CMV) infection — QSP model
### Cytomegalovirus infection and disease in transplant recipients · Quantitative Systems Pharmacology model

<a href="cmv_qsp_model.svg"><img src="cmv_qsp_model.png" width="880" alt="CMV QSP mechanistic map"></a>

| File | Contents |
|---|---|
| [`cmv_qsp_model.dot`](cmv_qsp_model.dot) · [`.svg`](cmv_qsp_model.svg) · [`.png`](cmv_qsp_model.png) | Mechanistic map — 171 nodes · 19 clusters · 219 edges |
| [`cmv_mrgsolve_model.R`](cmv_mrgsolve_model.R) | 48-ODE mrgsolve model + clinical decision-rule R driver + 19 scenarios |
| [`cmv_python_reference.py`](../../../cmv-infection-transplant/cmv_python_reference.py) | Independent Python/scipy re-implementation of the same 48 ODEs (**the verification original**) |
| [`cmv_reference_output.txt`](../../../cmv-infection-transplant/cmv_reference_output.txt) | The output of **actually running** the above file — the source of every number in this README |
| [`cmv_shiny_app.R`](../../../cmv-infection-transplant/cmv_shiny_app.R) | 12-tab Shiny dashboard |
| [`cmv_references.md`](cmv_references.md) | 116 papers · 19 sections · live PubMed lookup |
| [`mkrefs.py`](mkrefs.py) | Reference generator (a device for not writing PMIDs from memory) |

---

## Organising thesis — there is exactly one number every anti-CMV drug has to clear, and it is not the EC50

If systemic virus is written as a target-cell-limited process, then at quasi-steady
state for the virion pool the **net exponent** of the infected-cell pool becomes this
single line.

```
r = KPROD·(1−e_pol)·(1−e_pack) − DELI − KE8·E8eff − KENK·NKA
```

Two numbers actually measured in the clinic fix `KPROD` completely.

| Measurement | Constant derived |
|---|---|
| Untreated DNAemia doubling time ≈ **1.2 days** | `r0 = ln2/1.2 = 0.5776 /d` |
| Decline half-life on treatment ≈ **2.4 days** | `DELI = ln2/2.4 = 0.2888 /d` |
| | `KPROD = r0 + DELI = 0.8664 /d` |

Every argument about this disease then reduces to **two thresholds and one clock**.

### Threshold 1 (drug) — `e* = r0/KPROD = 0.667`

The fraction of virion production a regimen must remove for the sign of the exponent to
flip. **Potency, resistance, renal dose adjustment and drug interactions are all the same
question about this one number.**
(All values computed by execution — `cmv_reference_output.txt` §B)

| Regimen / strain | e | e−e\* | Verdict |
|---|---|---|---|
| valGCV 900 mg od (prophylaxis) | 0.890 | +0.223 | control |
| valGCV 900 mg BID (treatment) | 0.970 | +0.303 | control |
| valGCV 450 mg q2d (CrCl 30) | 0.770 | +0.104 | control |
| valGCV 450 mg twice weekly (CrCl 20) | 0.657 | **−0.010** | **breakthrough** |
| letermovir 480 mg od | 0.969 | +0.303 | control |
| letermovir 240 mg od (+CsA) | 0.940 | +0.274 | control |
| maribavir 400 mg BID | 0.792 | +0.125 | control |
| foscarnet 90 mg/kg q12h | 0.719 | +0.053 | control |
| valGCV 900 od / **UL97 mutant** | 0.112 | −0.555 | **breakthrough** |
| valGCV 900 BID / **UL97 mutant** | 0.335 | −0.332 | **breakthrough** |
| letermovir / **UL56 C325Y** | 0.010 | −0.656 | **breakthrough** |

### Threshold 2 (immunity) — `E8* = r0/KE8 = 4.81 CMV-specific CD8/µL`

The CD8 count at which the exponent turns negative **with no drug at all**. This is the
real endpoint of the disease, and **no drug can move this value.** Drugs only hold the
line while E8 climbs to it.

### The clock — the PCR interval is a dose: `2^(Δt/1.2 d)`

Pre-emptive therapy does not act on e at all. It sets **the height at which the exponent
begins to be reversed**.

| Monitoring interval Δt | Rise factor between draws | Actual value of a nominal 1000 IU/mL threshold |
|---|---|---|
| 3.5 days | 7.6-fold | 7,551 IU/mL |
| 7 days | **57.0-fold** | **57,018 IU/mL** |
| 14 days | 3,251-fold | 3,250,997 IU/mL |

The run results follow this arithmetic exactly (§G2):
q3.5d / q7d / q14d peak DNAemia = **3.44 / 3.94 / 4.62** log10,
P(disease) = **0.12 / 0.17 / 0.34**.
That is, monitoring frequency and drug potency are **quantities convertible into the same
units**, and the model prices each against the other.

### The coupling term — why the two thresholds fight each other

`d(E8)/dt` is proportional to **antigen**, and antigen is **the infected-cell pool the drug
has just removed**. Therefore **clearing threshold 1 delays reaching threshold 2.**
Late-onset CMV is not a device bolted on separately; it is exactly what this coupling term
does when the drug is stopped. Run result: at the moment 200 days of valganciclovir
prophylaxis ends, `E8eff = 0.094 /µL` (**1/51** of E8\* = 4.81), and **19 days later**
DNAemia crosses 1000 IU/mL again.

---

## Four structural decisions (each with what it derives)

### [1] The plasma DNA being measured is of two kinds — and the drugs hit different ones

- `Vv` = **capsid-packaged** virion DNA. The only infectious species, and what `BETA`
  multiplies.
- `Vl` = **free DNA** released from dying infected cells. At untreated quasi-steady state
  `Vl:Vv = 2:1`, i.e. **about 67% of what the PCR reports is not virion.**
- Polymerase inhibitors (GCV-TP · foscarnet · CDV-PP) reduce **both**.
- Terminase / kinase inhibitors (letermovir · maribavir) reduce only `Vv`, and because
  unit-length genomes are not cleaved from the concatemer they **actually increase DNA per
  dying cell** (`AMPPK`).

Immediate change in measured DNAemia computed with the infected-cell pool held fixed (§E):

| Regimen | e_pol | e_pack | Relative Vmeas | log10 reduction |
|---|---|---|---|---|
| no drug | 0 | 0 | 1.000 | 0.00 |
| ganciclovir 900 BID | 0.970 | 0 | 0.030 | **1.52** |
| foscarnet 90 q12h | 0.719 | 0 | 0.281 | 0.55 |
| **letermovir 480 od** | 0 | 0.969 | 1.000 | **0.00** |
| maribavir 400 BID | 0 | 0.792 | 1.000 | 0.00 |

Shaking `AMPPK` across the whole 0→1 range keeps letermovir's immediate reduction within
**0.17 to −0.12 log10**. The conclusion is derived, not assumed: **even if letermovir works
perfectly, plasma CMV DNA barely falls.** So DNAemia on letermovir cannot be an early
efficacy marker, and letermovir is a prophylactic drug, not a treatment drug.

### [2] Ganciclovir's activation step is a viral gene product — so one drug switches another off

pUL97 phosphorylates GCV to GCV-MP. Maribavir inhibits pUL97. Therefore maribavir enters
**not the effect term but the activation term** of ganciclovir, and the two drugs are
**structurally** antagonistic. Run values at treatment doses:

| | e |
|---|---|
| ganciclovir 900 mg BID alone | **0.970** |
| maribavir 400 mg BID alone | 0.792 |
| both together | **0.913** |

The combination is **0.057 lower** in e than ganciclovir alone. Both clear e\*, so it is not
a catastrophe, but adding maribavir on top of a failing ganciclovir course is not justified
by this arithmetic.

Meanwhile UL97 resistance is modelled as a **loss of kinase efficiency** (`FK_A = 0.125`)
and not as a polymerase EC50 shift. So that mutation does not touch maribavir (maribavir
resistance sits at **different residues** — T409M · H411Y · C480F, etc.). In a UL97
activation mutant: GCV alone e = 0.335 (fails e\*) · **MBV alone e = 0.792 (passes e\*)**.

### [3] The neutropenia loop closes through the kidney

Because 90% of GCV clearance is renal, CrCl alone determines **antiviral effect and marrow
toxicity simultaneously**. And two **errors in opposite directions** exist (§D, §G3).

| CrCl | Label dose | GCV Cavg | e(WT) | ANC steady state | Verdict |
|---|---|---|---|---|---|
| 100 | 900 od | 4.88 µM | 0.880 | 1.99 | OK |
| 60 | 900 od | 7.63 | 0.947 | 1.49 | OK |
| 40 | 450 od | 5.31 | 0.897 | 1.84 | OK |
| 30 | 450 q2d | 3.30 | 0.770 | 2.27 | OK (margin +0.10) |
| 20 | 450 twice weekly | 2.49 | **0.657** | 2.49 | **breakthrough** |

**When the dose is not reduced** (the classic ganciclovir marrow disaster):

| CrCl | Dose | GCV Cavg | e(WT) | ANC steady state | Verdict |
|---|---|---|---|---|---|
| 30 | 900 od | 13.20 µM | 0.982 | 0.97 | **Grade 3** |
| 20 | 900 od | 17.44 | 0.989 | 0.72 | **Grade 3** |

In simulation, the arm without renal adjustment (S11) spends **53 days with ANC < 1.0**,
reaches an ANC nadir of 0.26, and the dose is repeatedly halved and restarted on a
two-week cycle — exactly the pattern seen in real practice. Conversely the correctly
reduced arm (S10) protects the marrow but its **margin above e\* thins to +0.10.**

### [4] Resistance is a **shape**, not a magnitude — and its risk is **duration**

The same "single mutation" operates in a completely different shape in the two drugs.

| | Fold shift | Change in e | Residual effect |
|---|---|---|---|
| letermovir · UL56 C325Y | 3000-fold | 0.969 → **0.010** | **1.1%** (complete loss) |
| ganciclovir · UL97 M460V | 8-fold | 0.890 → **0.112** | 12.6% (crippled, but not zero) |

Yet even the latter, where dose escalation might in principle recover ground: clearing e\*
in a UL97 mutant requires GCV-TP **264.7 µM-eq**, i.e. a plasma Cavg of **20.4 µM = 2.0
times the licensed treatment dose**. **Escalation is not the answer; the drug has to be
changed** — that is the conclusion of this table.

And **a deterministic mutation flux was not used** (see "defects revealed by execution",
number 3 below). Instead **the number of days spent inside the selection-pressure window**
is integrated: the state in which wild type is suppressed (e_W > e\*) while at the same time
the mutant's own R_eff exceeds 1 — the only state in which a minority variant becomes the
dominant species.

| Arm | UL97 selection window | UL56 selection window |
|---|---|---|
| pre-emptive q3.5d | 46 days | 0 |
| pre-emptive q7d | 61 days | 0 |
| pre-emptive q14d | 83 days | 0 |
| valGCV prophylaxis 200 days | **200 days** | 0 |
| letermovir prophylaxis 200 days | 0 | **202 days** |

**Prophylaxis keeps a patient under selection pressure 3–4 times longer than pre-emptive
therapy.** The monitoring interval alone moves it by a factor of 2.

---

## Scenarios (19 arms, 365 days post-transplant)

The full table is in [`cmv_reference_output.txt`](../../../cmv-infection-transplant/cmv_reference_output.txt) §F. Extract:

| Arm | Peak VL (log10) | First cs-CMVi | P(disease) | E8@200 | ANC nadir | Selection window | Cost ($k) |
|---|---|---|---|---|---|---|---|
| S1 D+/R− untreated (natural history) | 4.83 (d38) | — | **0.60** | 4.31 | 3.19 | 0 | 0 |
| S2 pre-emptive q7d | 3.94 | d21 | 0.17 | 1.32 | 0.65 | 61 | 49.8 |
| S3 pre-emptive q3.5d | 3.44 | d18 | **0.12** | 3.21 | 1.07 | 46 | 21.0 |
| S4 pre-emptive q14d | 4.62 | d28 | 0.34 | 2.09 | 0.74 | 83 | 39.7 |
| S5 valGCV 100 days, no monitoring | 4.73 (d138) | — | 0.42 | 3.72 | 1.24 | 100 | 4.2 |
| S6 valGCV 200 days, no monitoring | 4.67 (d238) | — | 0.36 | 0.09 | 1.24 | 200 | 8.4 |
| S7 valGCV 200 days + monitoring | 4.22 | d228 | 0.15 | 0.09 | 0.77 | 245 | 29.3 |
| S8 letermovir 200 days + monitoring | 4.06 | d228 | 0.12 | 0.09 | 0.75 | 202 | 65.9 |
| S10 CrCl 30, correct dose reduction | 4.37 | d228 | 0.22 | 0.09 | 1.11 | 214 | 12.6 |
| S11 CrCl 30, **reduction omitted** | 4.10 | d228 | 0.17 | 0.09 | **0.26** | 204 | 67.4 |
| S12 D+/R+ pre-emptive q7d | 2.97 | not reached | 0.06 | 4.83 | 3.19 | 0 | 8.2 |
| S13 D−/R+ basiliximab | 1.78 | — | 0.00 | 8.22 | 3.19 | 0 | 8.2 |
| S15 valGCV 200 days + mTORi switch | 3.99 | d228 | **0.10** | 0.09 | 1.05 | 205 | 15.0 |

Two things to note. **(a) S3 dominates S2 on every axis** — doubling the PCR frequency gives
a peak DNAemia 0.5 log10 lower, a selection window 15 days shorter, days with ANC < 1.0
falling from 25 to 0, and the total cost actually halves (because treatment days and G-CSF
fall). **(b) letermovir matches valganciclovir on disease rate** (0.12 vs 0.15) **and
definitely protects the marrow** (ANC nadir 0.75 vs 0.77 — both occurring during the
treatment-dose window) **but costs 2.2 times as much.**

---

## What this model produced on its own — and what it got wrong on its own

### (i) Extending the prophylaxis duration does not **prevent** the event, it **postpones** it

Sweeping valganciclovir prophylaxis duration from 0 to 365 days in an arm with no
monitoring after prophylaxis ends (§I):

| Prophylaxis days | ISI@stop | E8@stop | Day of relapse | P(disease) |
|---|---|---|---|---|
| 0 | 2.07 | 0.00 | 15 | 0.600 |
| 60 | 1.91 | 0.05 | 78 | 0.473 |
| 100 | 1.59 | 0.07 | 119 | 0.422 |
| 200 | 1.35 | 0.09 | 219 | 0.359 |
| 300 | 1.19 | 0.12 | 320 | 0.319 |

**The day of relapse tracks the day of discontinuation exactly (always 18–20 days later).**
There is only one reason P(disease) nonetheless falls — the *same* relapse meets a *less
immunosuppressed* host (ISI 1.91 → 1.19). E8 at discontinuation barely moves, from 0.05 to
0.12. Because the drug removed the antigen that would have made it.

**An honest miss.** The IMPACT trial reported CMV disease at 36.8% for 100 days vs 16.1% for
200 days. This model produces **42.2% vs 35.9%** — **the direction is right and the
magnitude is about 1/3 of the observed effect.** That is, something that makes long
prophylaxis work is missing from this model. The next section is the most plausible
candidate the model can test on itself.

### (ii) The training-dose U-curve — complete suppression may not be optimal

If E8 proliferation is antigen-dependent, then prophylaxis that suppresses **completely**
leaves nothing durable behind, while suppressing **not at all** produces disease. Sweeping
adherence (which maps monotonically onto e) in 200-day prophylaxis, the curve is not
monotonic (§J):

| Adherence | e(WT) | vs e\* | E8@200 | P(disease)@200 days |
|---|---|---|---|---|
| 1.00 | 0.890 | above | 0.09 | 0.359 |
| 0.50 | 0.668 | above | 0.10 | 0.347 |
| 0.45 | 0.620 | **below** | 0.17 | 0.242 |
| **0.40** | **0.563** | **below** | **0.38** | **0.163** ← minimum |
| 0.35 | 0.497 | below | 0.65 | 0.190 |
| 0.25 | 0.335 | below | 1.96 | 0.432 |
| 0.00 | 0 | — | 4.31 | 0.600 |

The minimum sits **below e\*** (e = 0.563). There is only one mechanism — partial suppression
holds the load below the disease threshold while leaving enough antigen to train E8.

**Do not read this as a prescribing recommendation.** The three reasons are set down as
they stand. (1) The selection-window column collapses to 2 days at the minimum, but that
is an artefact of the definition and no grounds for reassurance — selection requires wild
type to be suppressed, and below e\* it is not suppressed, so what partial suppression
actually buys is **continuous wild-type replication**, i.e. the substrate on which mutants
*arise*. (2) It is a single deterministic patient. A real cohort spreads out along this
curve, and some of it lands on the rising limb, worse than complete suppression. (3) The
position of the minimum is set by `KAG` (the antigen half-saturation constant for CD8
proliferation), which is the most poorly measured parameter in the model. Sensitivity:

| KAG | Position of minimum (adherence) | P(disease) at the minimum | P(disease) at adherence 1 |
|---|---|---|---|
| 0.02 | 0.40 | 0.054 | 0.307 |
| 0.05 | 0.40 | 0.163 | 0.359 |
| 0.12 | 0.40 | 0.326 | 0.420 |
| 0.30 | 0.40 | 0.524 | 0.526 |

The U shape survives across the whole range but **the depth of the minimum does not** — it
is a qualitative prediction, not a quantitative one.

### (iii) A refuted expectation — reported rather than deleted

S16–S19 were the arms built to show that "in ganciclovir-resistant breakthrough, switching
drugs beats continuing". **In P(CMV disease) it did not, and the arm that continued was the
lowest** (§G4):

| Arm | Day judged refractory | Peak VL | P(disease) | **Mutant peak** | eGFR nadir | Mg nadir |
|---|---|---|---|---|---|---|
| S16 → maribavir | d42 | 4.06 | 0.17 | **1785** | 82.9 | 2.00 |
| S17 → foscarnet | d42 | 4.01 | 0.17 | 1785 | **74.6** | **0.69** |
| S18 → GCV+MBV | d42 | 4.06 | 0.15 | 1785 | 82.7 | 2.00 |
| S19 continue GCV | d42 | 4.11 | 0.13 | **2487** | 83.3 | 2.00 |

The reason is visible inside the table. By the time the guideline definition of refractory
CMV is met (d42) the CD8 response is already rising, and in this single deterministic
patient **it is immunity, not the drug, that ends the event.** What actually separates the
arms is the very quantity the switch was aimed at — **peak virion DNA of the resistant
strain, 1785 on maribavir vs 2487 continuing GCV (1.4-fold)** — plus the price of the
foscarnet route (eGFR nadir 74.6 vs 82.9, magnesium 0.69 vs 2.00 mg/dL). So the claim this
model can defend is narrower than the one it set out to make: **switching reduces resistant
strain burden and the selection-pressure window; a disease-rate benefit is not demonstrated
in a host whose T cells are recovering anyway.**

### (iv) Letermovir's tacrolimus interaction eats part of its own benefit

Letermovir raises tacrolimus AUC 2.4-fold through CYP3A/OATP1B1 inhibition. Comparison of
the arm that anticipated this and halved the dose (S8) with the arm that left it alone (S9)
(§G5):

| Arm | TAC@d100 | ISI@d100 | P(disease) | P(rejection) |
|---|---|---|---|---|
| S8 tacrolimus reduced | 5.90 | 1.58 | 0.118 | 0.074 |
| S9 not reduced | **14.13** | **2.15** | 0.127 | 0.066 |

The effect is small but its sign is unambiguous, and **it shows exactly the trade the
clinician actually makes** — immunosuppression rises, CMV risk goes up by +0.9 percentage
points and rejection risk goes down by −0.8 percentage points.

---

## Eight defects revealed by execution (found by integration, not by inspection)

The 48 ODEs of this model were independently re-implemented in Python/scipy and actually
integrated. The defects that surfaced and were fixed in the process are left here rather
than deleted.

1. **The tacrolimus interaction was boolean.** Written as `if (letermovir present) clearance
   /= 2.4`, **solver rounding at the 1e-12 level** in the letermovir state variable latches
   the switch permanently. In arms that never received letermovir, tacrolimus had risen from
   9 to 32 ng/mL. Replaced with a concentration-dependent Emax term.
2. **Zero in the resistant-strain compartment was an unstable equilibrium.** While the
   strain's exponent is positive, solver rounding of 1e-16 grows at 0.4–0.5/d and becomes
   **macroscopic in about 70 days**. Every arm was "discovering" resistance out of
   floating-point noise. Fixed by putting an extinction floor at one systemic infected cell
   (= 1e-8).
3. **A deterministic mutation flux makes resistance inevitable.** Because MU·NCELL ≈ 3e3, any
   meaningful replication creates numerous mutant lineages and the probability saturates at
   1 in every arm (in flat contradiction with the roughly 5% observed under valGCV
   prophylaxis). Replaced with the integral of **days resident in the selection-pressure
   window** — because what actually distinguishes the arms is duration, not probability.
4. **Friberg γ = 0.17 is a value fitted to transient cytotoxic exposure.** Under chronic
   dosing the steady state becomes `ANC = CIRC0·(1−Edrug)^(1/γ)`, which dropped a patient on
   mycophenolate alone to an ANC of 0.83 and fired the neutropenia rule in every arm.
   γ · EMAX · E_MPA were refitted to chronic-exposure ANC (valGCV 900 od + MMF → 1.99,
   letermovir + MMF → 3.61).
5. **GFR0 was fixed at 55 while the initial value was CrCl/1.1.** Every run was quietly
   declining from 86 to 55 and reporting that as CMV nephropathy.
6. **The sanctuary compartment grew without bound in every arm** — including the untreated
   arms. No condition existed under which clearance exceeded the growth rate.
7. **Monitoring and disease costs were not being billed.** So the arm ordering twice as many
   PCRs looked cheaper.
8. **The refractory switch fired on a fixed date** — that is, before the resistant strain had
   done anything. So the four salvage arms were numerically identical. Replacing it with the
   guideline criterion (≥14 days of treatment, <0.5 log10 reduction) separated the arms, and
   that result is the refutation in (iii) above.

---

## Calibration ledger (what was fitted to what)

| Quantity | How it was set |
|---|---|
| `KPROD`, `DELI` | Doubling time 1.2 days · on-treatment decline half-life 2.4 days (2 measurements) |
| `LAMT = 2.00` | **One number** fitted to the observed peak DNAemia of untreated primary D+/R− infection (10^4.83) |
| `KDIS = 0.0300` | **One number** fitted to the 1-year P(CMV disease) = 0.60 of untreated D+/R− |
| `GAM` · `EMAXMYE` · `EMPA` | Refitted to chronic-exposure ANC (defect 4 above) |
| All EC50 values | Published in-vitro EC50, **corrected to free drug concentration** for highly protein-bound drugs |
| All PK parameters | Published population PK (F · V · CL · t½ · renal fraction) |
| `AMPPK = 0.50` | **An assumed value** — never directly measured in humans. A sweep across the whole 0–1 range is presented |
| `KAG = 0.05` | The **most fragile** parameter in the model. Governs the depth of the training-dose U-curve minimum |
| Costs | Illustrative list-price level; only the order of magnitude is meaningful |

## What was deliberately left out

- **UL54 polymerase mutations** appear only in the map and the cross-resistance table, not as
  a separate ODE strain.
- **There is no between-individual variability.** All the probabilities above are hazard
  function values for a single deterministic virtual patient and cannot be quoted as trial
  incidences.
- **Rejection and eGFR are reduced-form hazard surrogates**, not a mechanistic alloimmune
  model. They are included only because reducing the dose for neutropenia really is a driver
  of rejection, not because they are calibrated.
- **Adoptive T-cell therapy** has only an input term (`ACT`) and no dosing schedule. Vaccines
  appear only in the map.
- **Whole blood vs plasma conversion** and the 1–2 log10 inter-laboratory variation are not in
  the model.
- **HSCT-specific engraftment kinetics and GVHD** are absent. The reversal of risk in HSCT
  (recipient seropositivity becoming the high-risk group) is represented only through the
  initial memory T-cell pool.

---

## Reproduction

```bash
# mechanistic map
dot -Tsvg cmv_qsp_model.dot -o cmv_qsp_model.svg
dot -Tpng -Gdpi=150 cmv_qsp_model.dot -o cmv_qsp_model.png

# verification run (the source of every number in this README)
python3 cmv_python_reference.py > cmv_reference_output.txt

# re-query the references
python3 mkrefs.py --refresh

# mrgsolve / Shiny  (no R toolchain in this environment, so these were not run — see the caveat below)
Rscript -e 'library(mrgsolve); mod <- mread("cmv_mrgsolve_model.R")'
Rscript -e 'shiny::runApp("cmv_shiny_app.R")'
```

> **The caveat is stated explicitly.** The 48 ODEs were integrated and verified, but what did
> that was `cmv_python_reference.py`. There was no R/mrgsolve toolchain in this environment,
> so `cmv_mrgsolve_model.R` and `cmv_shiny_app.R` are equation-by-equation transcriptions of
> the Python reference implementation but **were not themselves run.**

---

## Disclaimer

This model is a mechanistic exploration tool for **educational and research purposes**. It
must not be used for clinical care, prescribing or diagnosis. In particular the training-dose
U-curve in (ii) above is **not a prescribing recommendation**; partially suppressive
antiviral therapy raises the risk of resistance and disease in real patients. All probability
values are hazard functions for a single virtual patient and are not clinical trial
incidences.
