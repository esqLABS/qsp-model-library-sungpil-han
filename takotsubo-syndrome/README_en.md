# Takotsubo Syndrome — QSP Model

**Takotsubo syndrome · stress-induced / catecholamine-mediated cardiomyopathy · TTS**

| Deliverable | File |
|---|---|
| Mechanistic map (218 nodes · 374 edges · 26 clusters) | [`tts_qsp_model.dot`](tts_qsp_model.dot) · [SVG](tts_qsp_model.svg) · [PNG](tts_qsp_model.png) |
| mrgsolve ODE model (64 compartments) | [`tts_mrgsolve_model.R`](tts_mrgsolve_model.R) |
| Shiny dashboard (14 tabs) | [`tts_shiny_app.R`](tts_shiny_app.R) |
| References (125, in 15 sections) | [`tts_references.md`](tts_references.md) |

---

## The one structural claim

> **"Apical ballooning" is not a location but a threshold.**
> The segment with the highest β2-adrenoceptor density is the first to cross the
> Gs→Gi switch threshold.

The model has three myocardial segments (apex · mid · base), and they differ in
**exactly two values**. Both are measurable quantities, and neither mentions
"contraction".

| | Apex | Mid | Base |
|---|---|---|---|
| `RHO` total β-AR density | 1.40 | 1.15 | 1.00 |
| `FB2` β2 fraction | 0.42 | 0.32 | 0.24 |
| (`INN` sympathetic innervation, the ¹²³I-MIBG reverse gradient) | 0.62 | 0.90 | 1.20 |
| (`TH` wall thickness) | 0.36 | 0.72 | 1.00 |

And **one** agonist-occupancy-dependent Gs→Gi switch function is applied
**identically** to all three segments.

```
SIG_x    = ((occ2_x − THR_SW)+)² / (… + KSW²) · (1 + WPKA·PKA_x) · (β2 density factor)
PHOS_x′  = KPH·SIG_x·(1 − PHOS_x) − KDEPH·PHOS_x
GI_x     = b2_x · occ2_x · PHOS_x                    ← where density enters
GS_x     = b1_x · occ1_x + b2_x · occ2_x · (1 − PHOS_x)
AC_x     = (AC_BAS + EMAX_AC·GS_x/(GS_x+KM_GS)) / (1 + GAM_GI·GI_x)
```

Nowhere in this file is it written that "the apex is akinetic". The word `APEX`
appears only in the four parameter values above and in output names. Apical
akinesis, basal hypercontractility, the LVOT gradient, and the ballooning index are
all **one equation read at three receptor densities**.

### A corollary worth writing down

**At rest the three segments are effectively indistinguishable.** Because the
apex's high receptor density is offset by the apex's low sympathetic innervation,
resting Gs coupling is almost the same across segments: 0.149 / 0.147 / 0.152
(apex / mid / base). The apex is **not an intrinsically weak segment.** It becomes
weak only when a β2 agonist arrives by way of the bloodstream.

That is why this disease requires **adrenaline** rather than sympathetic release
(noradrenaline). Noradrenaline is about 25-fold weaker at β2 and therefore cannot
recruit Gi in any segment — this is falsification test 1 (S06).

---

## The second, forced result: the benefit is a product

```
benefit = (added contractile drive) × (fraction delivered without raising receptor phosphorylation)
```

Every drug that raises cAMP (dobutamine · adrenaline · dopamine · **milrinone**)
raises PKA, PKA raises SIG, SIG raises PHOS, PHOS raises GI, and GI divides AC. The
second factor becomes **negative** once the apex has crossed the threshold — the
same feedback term read at a higher density.

Levosimendan/OR-1896 enters through `CASENS` (cTnC Ca²⁺ sensitisation) and does not
touch cAMP, so its second factor is structurally 1.

The two drugs are therefore compared **not at the same dose but at the same added
drive** (`matched_drive()`). Comparing them at the same dose is comparing potency,
and says nothing about the pathway.

**Milrinone is in the model for exactly one reason**: it raises cAMP *downstream* of
the receptor. If the harm were a side effect of receptor occupancy, milrinone should
be safe. The model says it is not (falsification test 3).

---

## The third result: the sign of β blockade is conditional, and the condition is a number inside the model

Lowering occupancy removes drive at the base (harm, if antegrade flow is the
problem) and removes phosphorylation at the apex (benefit, if the switch is the
problem). Which prevails is decided by the LVOT gradient term, and the existence of
that term is the patient parameter **`SEPT`** (basal septal morphology) multiplied
by basal hypercontractility.

The registry observations — that short-acting β blockade helps in obstructive TTS,
and that β blockade does not reduce recurrence while ACE inhibition does — are a
**consequence** of this pathway allocation, not a look-up table.

The reason basal hypercontractility persists for days is also inside the model. When
antegrade flow falls, the baroreceptor reflex sustains interstitial noradrenaline
(`REFL`). It is this loop that holds the base in a hypercontractile state after the
primary surge has gone, and short-acting β blockade helps the obstructive patient
because it breaks this loop.

---

## Fourth: four counterfactual integrators — no mechanism's share is assumed

The half-life of cAMP is about one minute, tens to thousands of times faster than
every other state variable. Its quasi-steady-state (QSS) value can therefore be
obtained in closed form, which makes exact counterfactual computation cheap.

| Integrator | Meaning |
|---|---|
| `HD_GI` | the EF that would be obtained had PHOS been 0 (= the pertussis toxin experiment) |
| `HD_OED` | the EF that would be obtained had oedema been 0 |
| `HD_ATP` | the EF that would be obtained had energy metabolism been normal |
| `HD_LVOT` | the forward stroke volume that would be obtained had the gradient been 0 |

The QSS assumption is not left as an article of faith. `CQ_AP` (QSS cAMP) is output
alongside the integrated `CAMP_AP`, so agreement between the two can be checked on
every run.

---

## Validation

`mrgsolve 2.0.1`, LSODA, atol/rtol 1e-8. Time unit = **hours (h)**. Adrenaline with
a 2-minute half-life, receptor dephosphorylation with a 4-day half-life, and 1-year
recurrence risk all have to live in the same system, so the unit is set by the
fastest state and the slow things are integrated out to 8760 h.

### The baseline state is derived

There are no hand-entered initial values. `$MAIN` computes the resting reference
values (per-segment Gs, cAMP, PKA, contractile drive, and the reference stroke
volume `SV0R`) from the parameters, and every normalisation divides by those values.
As a result, whatever value `RHO_AP` is changed to, at t=0 `CONT_x = 1`,
`CAL_x = 1`, `ATP_x = 1`, `GRK2 = 1`, and `RDN = 0` hold **exactly**. The
unstimulated control (S01) does not move to the decimal place over 90 days from
LVEF 62.0 %, CO 5.00 L/min, MAP 80 mmHg, NT-proBNP 90 pg/mL, QTc 418 ms.

### Defects exposed and fixed during integration (all of them structural errors)

1. **Using NT-proBNP as an absolute flux.** At a gain of `KBNP_I` = 900 pg/mL/h, a
   7 % fall in the transcriptional signal made the steady-state concentration
   **negative**. Rewritten as a fold-change on the elimination constant.
2. **Non-saturating natriuresis.** As a term proportional to BNP, a 20-fold rise in
   NT-proBNP drained 1.9 L per hour and dried the circulating blood volume down to
   0.8 L. Not a fine calibration error but a structural one.
3. **Perfusion linear in pressure.** Without coronary autoregulation, a small fall
   in MAP alone started an ATP → contractility → MAP downward spiral that killed
   every run within 48 hours. An autoregulation function in which flow is almost
   pressure-independent above 60 mmHg, plus coronary flow reserve (`CFR` = 2.5),
   were added.
4. **Unbounded oedema and stunning.** The relaxation term could not beat the
   generation term, so both variables saturated at 1.0 and contractility went
   negative. Changed to saturating forms, with a check that the feedback loop gain
   is less than 1 (oedema 0.19 + stunning 0.5 = 0.69).
5. **Unbounded systolic bulging.** `SH_AP` went to −1.5. A segmental shortening
   fraction physically cannot take that value. A saturating bulging term plus a
   floor of `SHMIN` = −0.25.
6. **A squared LVOT gradient.** Making the gradient proportional to the square of
   basal hypercontractility produced 472 mmHg. The dynamic gradient is limited by
   the pressure the basal segment can generate, so it was changed to a saturating
   form (`GMAXG` = 150 mmHg).
7. **Excessive α1 gain.** At `EA1_SVR` = 1.6, MAP was 146 mmHg. Lowered to 0.45
   instead.
8. **No deadband on the reflex feedback.** A numerical 0.1 % error in cardiac output
   switched on the reflex sympathetic drive and contaminated the baseline state
   (baseline troponin 0.18 ng/mL). Changed to a 3 % deadband plus a saturating form.
9. **No receptor-independent reserve in contractility.** With 92 % of resting cAMP
   receptor-driven, blocking β1 by 70 % alone took contractility to zero. `AC_BAS`
   was raised to 40 % of the resting flux, so that complete β blockade consumes
   about 30 % of contractility.
10. **A TdP hazard function of the wrong form.** Integrating an exponential hazard
    function with a constant term over 90 days gave **a cumulative probability of
    4.9 % in the unstimulated control**. At a normal QTc the hazard should not
    merely be small, it should be **zero**. Changed to an `exp(...) − 1` form so
    that the hazard vanishes below 450 ms, with the coefficient matched to the
    reported incidence of 2–5 % (the original form gave 17 % in the reference arm).
11. **Taking the carvedilol Ki from a binding assay (referenced to free
    concentration) and comparing it against total plasma concentration.**
    Carvedilol is 98 % protein bound, so this is a 50-fold error. At 12.5 mg twice
    daily it blocked β1 by 96 % and dropped the 1-year LVEF to 17 %. Resetting Ki to
    a total-concentration basis (0.0008 → 0.040 mg/L) gave 35 % blockade and a
    1-year LVEF of 52 %. Metoprolol was over-potent by 3-fold for the same reason
    (57 % blockade at 50 mg twice daily).
12. **`$TABLE` recomputing the algebraic layer.** This is this repository's 2026-07
    lesson. If `$ODE` and `$TABLE` compute the same flux with slightly different
    gating, then the trajectory reported is not the trajectory integrated. The whole
    algebraic layer was made into **a single macro** (`TTS_ALGEBRA`, 168 shared
    variables) that expands to the same text in both blocks. A macro cannot drift.

### Calibration

Seven parameters were fitted to eight public anchors by Nelder-Mead. Objective
function **6.5 × 10⁻⁵**, all eight anchors **within 1.6 %**.

| Anchor | Target | Predicted | Relative error |
|---|---|---|---|
| LVEF, presentation echocardiogram at 24 h (%) | 38 | 38.1 | +0.2% |
| LVEF at 30 days (%) | 58 | 58.9 | +1.6% |
| Peak ballooning index | 0.75 | 0.754 | +0.6% |
| Peak hs-cTnI (ng/mL) | 6.0 | 6.03 | +0.4% |
| Peak NT-proBNP (pg/mL) | 4200 | 4210 | +0.2% |
| NT-proBNP : troponin ratio | 700 | 699 | −0.2% |
| Peak QTc (ms) | 488 | 493 | +1.1% |
| Time of peak QTc (days) | 2.5 | 2.51 | +0.4% |

**What was fixed, and why.** The oedema kinetics (`KOED`, `KOEDR`), the stunning
kinetics (`KSTUN`, `KSTUNR`), and the oedema→QTc slope (`AQT_OED`) were fixed. Left
free, `KOED` and `AQT_OED` cancel each other indefinitely (what the QTc anchor sees
is their product). In the first fit that did leave them free, `AQT_OED` climbed to
25 856 and `KOED` fell to 0.0008 — the objective function improves and the model
becomes uninterpretable.

The reflex sympathetic gain `KREFL` was fixed too, and this one matters more. Left
free, the optimiser pushed it up to **70 000 nmol/L/h**. Because the reflex is
saturating (`CODEF/(CODEF+REFL50)`), a value that large turns a graded baroreflex
into an **almost binary switch** that trips on a 3 % fall in cardiac output.
Interstitial noradrenaline leaves the physiological range and β1 occupancy
saturates, with the result that **the size of the primary stimulus no longer
matters.** The visible symptom was that the premenopausal arm (whose surge is
2.9-fold smaller) stopped being milder — a reflex that does not care about the size
of the primary injury had washed out the oestradiol effect. `KREFL` was therefore
fixed at 3300, and at that value the peak interstitial noradrenaline stays at
**55.6 nmol/L** (about 7 times resting). A better objective function bought with a
non-physiological reflex is not a better model.

**A second calibration target: TdP.** Separately from the eight anchors above,
`H0_TDP` was matched to the reported TdP incidence of 2–5 %. That gives 2.4 %
cumulative over 90 days in the reference arm, 0.00 % in the unstimulated control,
and 28.3 % in the QT-prolonging co-medication arm.

**The LVEF anchor is the presentation echocardiogram (24 h), not the instantaneous
nadir.** With `min(LVEF)` as the anchor the objective function depends on the output
resolution: the reported nadir at delta = 1 h is 38 %, but at delta = 0.25 h it is
33 %. That is because the model dips transiently deeper at the catecholamine peak in
the first two hours, and no patient has an echocardiogram at that instant. Reporting
a resolution-dependent quantity as a fitted value is a quiet artefact.

### Reproduced without being anchored (held-out)

**1. The threshold dissociation comes from density, not occupancy.** Ramping
adrenaline slowly, at a plasma level of 1.35 nmol/L **72 %** of the apical β2 pool
has switched to Gi while the base is at **9 %** (at 1.88 nmol/L, 93 % versus 16 %).
Yet at the same moment the β2 **occupancy is higher at the base** (0.585 versus
0.486, because there is more noradrenaline at the base). That what creates the
separation is density and not occupancy is confirmed as an output of the model.

**2. Falsification 1 — a pure noradrenaline surge does not produce ballooning.**
Ballooning index 0.041 versus 0.754 (a 95 % reduction), LVEF nadir 57.9 % versus
37.1 %.

**3. Falsification 2 — pertussis toxin abolishes it.** Ballooning index 0.014, Gi
contribution to the EF deficit 0.003.

**4. Falsification 3 — raising cAMP downstream of the receptor produces the same
harm.** Dobutamine raises the ballooning index to 0.900 and milrinone to 0.910
(untreated 0.754), and under dobutamine the apical shortening fraction **falls to
−0.047, i.e. becomes dyskinesis.** The Gi contribution to the EF deficit goes
untreated 0.641 → dobutamine 0.735 → milrinone 0.792. The harm is not a side effect
of receptor occupancy; it comes from cAMP itself.

**5. Both inotropes worsen the integrated EF deficit.** The cumulative EF deficit
(EF·h) is 78.2 untreated, 83.0 with levosimendan, 88.9 with dobutamine. The
direction of the pathway argument (levosimendan < dobutamine) is confirmed, but
**both drugs are worse than no treatment.** Matched on added drive (dobutamine
113.7 mg/h versus levosimendan 0.260 mg/h), dobutamine gives the higher day-1 EF
(67.9 % versus 59.3 %) while levosimendan gives the higher day-7 EF (42.4 % versus
40.8 %). That is, dobutamine buys now and pays later.

**6. The oestradiol effect.** For the same trigger, 24-hour LVEF is 51.4 %
premenopausal (E2 = 100 pg/mL) versus 38.1 % postmenopausal, ballooning index 0.267
versus 0.754, troponin 4.65 versus 6.03 ng/mL.

**7. Catecholamine levels fall inside the reported range.** Peak plasma adrenaline
4.83 nmol/L (16 times resting; literature 2–34 times).

**8. The reverse gradient inverts the dispersion.** Inverting `RHO`/`FB2`/`INN`
makes the ballooning index's minimum **−0.514** (that is, the base contracts less
than the apex) and the LVEF nadir a milder 52.7 %.

**9. LVOTO prevalence falls inside the registry range — but it depends on which
gradient you count.** In a 150-patient virtual cohort, the proportion whose
**sustained** gradient (mean over 4–48 hours) exceeds 30 mmHg is 8.7 %, and the
proportion exceeding 25 mmHg is 22 %, bracketing the reported 10–25 %. But counting
the **peak** gradient gives 61.3 %, a gross over-prediction. Registry LVOTO is the
sustained gradient measured at the bedside and not the instantaneous spike at the
catecholamine peak, so which quantity is being compared changes the answer more
than 3-fold. The comparable quantity in this model is the sustained gradient.

**10. Inotropes enlarge the obstruction.** The 4–48 hour mean LVOT gradient rises
from 12.0 mmHg untreated to 33.3 with dobutamine, 37.1 with milrinone, and
39.2 mmHg with levosimendan. In the obstructive phenotype it goes from 41.5
untreated to 71.2 with dobutamine and 91.4 mmHg with levosimendan. Everything that
raises basal contractility pays the same price.

**11. Only ACE inhibition moves 1-year recurrence.** Untreated 2.40 %, ramipril
1.45 % (−40 %), metoprolol 2.37 % (−1.2 %), carvedilol 2.36 %. A consequence of the
pathway allocation, not a look-up table.

### What it fails to reproduce (written down rather than hidden)

**F1. It fails to reproduce the conditional benefit of β blockade — this is the
biggest failure.**
Sweeping the esmolol dose in the obstructive phenotype, the mean LVOT gradient does
fall clearly, 38.3 → 16.7 mmHg (at 105 mg/h), but the mean cardiac output during
infusion **actually worsens by 1.9 %**, 3.841 → 3.767 L/min, and the LVEF nadir
worsens from 30.0 to 24.7 %. At no dose is there a net gain in antegrade flow.
Cutting the gradient from 38 to 17 mmHg reduces the antegrade flow loss from 0.41 to
0.23 and so buys about 18 % of flow, but the same blockade pays more than that out
of basal contractility. The model **disagrees with** the clinical observation that
short-acting β blockade helps in obstructive TTS. The conditional sign structure is
inside the model, but the benefit/cost ratio is mis-set.

**F2. The isolated mid-ventricular ring cannot be produced by any monotone
apex > mid > base gradient.** The model produces only apical or apical-plus-mid
forms. This is not a problem to be removed by parameters but a limit of the
structure.

**F3. Reverse takotsubo is not derived.** It is reproduced only by putting in an
inverted gradient, so the model presents it not as a result but as a **prediction**:
reverse takotsubo must have an inverted receptor gradient. Section 14 of the
references supports the existence of the phenotype, but no paper has measured an
inverted receptor gradient.

**F4. It grossly over-predicts embolic risk.** In the reference arm without
anticoagulation the 90-day embolic probability comes out at 36.0 % (apixaban 5.7 %),
whereas the embolic rate in clinical TTS is of the order of 1–3 %. The TdP hazard was
matched to an incidence, but no anchor could be found for `H0_EMB` and it was
therefore left uncalibrated. So **the absolute values of the thrombus/embolism arm
must not be used, and only the between-arm comparison should be read** (all that can
be read is the direction — that apixaban reduces the peak thrombus burden
0.351 → 0.123 and the embolic probability 36.0 % → 5.7 %).

**F5. The oedema time constant was matched to the QTc time course (2–3 days), which
is faster than CMR T2 normalisation (weeks to months).** The two observations speak
to different time scales, and the model followed the electrophysiological one. It
therefore underestimates the persistence of tissue oedema.

**F6. There is no death in this model.** Below a cardiac output of about
1.5 L/min the haemodynamic algebra leaves the domain in which it means anything (MAP
collapses to the CVP floor and the integration fails). The SAH arm (`neuro`) was
placed **just inside** that boundary — the originally chosen `AMP_TOT` = 600 /
`TAU_SUR` = 48 h went past it. The boundary is written down here rather than hidden
behind solver tolerances.

**F7. The specific values of `RHO_AP` · `FB2_AP` · `INN_AP` · `TH_AP` carry no
PMIDs.** The existence and the direction of the gradients are supported by the
literature; their magnitudes are not. These four numbers are the model's principal
falsifiable inputs.

**F8. There are no atherosclerotic plaque, coronary thrombus, or infarct
compartments.** This model cannot distinguish TTS from coronary artery disease —
that distinction is a premise of the model, not an output.

### Numerical hygiene

- The unstimulated control (S01) does not move to the decimal place over 90 days
  from LVEF 62.0 %, CO 5.00 L/min, MAP 80 mmHg, PCWP 9.0 mmHg, NT-proBNP 90 pg/mL,
  cTnI 0.010 ng/mL, QTc 418 ms.
- The relative difference between the quasi-steady-state cAMP used for the
  counterfactuals and the integrated cAMP (`CQERR`) is output for every arm.
  Excluding the sharply changing interval at the surge peak it is under 1 %.
- The per-mechanism shares (`SHR_*`) are **ratios** and must be read together with
  the denominator (`EFdef_hours`), and because each counterfactual removes only one
  mechanism at a time **they do not sum to 1.**

### Principal arm results (delta = 0.5 h, 90 days)

| arm | LVEF 24 h | LVEF nadir | LVEF 7 d | LVEF 30 d | Ballooning peak | Gradient peak/48 h mean | CO nadir | cTnI | NT-proBNP | QTc (day) | TdP% | Gi share |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| S01 unstimulated | 62.0 | 62.0 | 62.0 | 62.0 | 0.000 | 0 / 0.0 | 5.00 | 0.01 | 90 | 418 (0.0) | 0.00 | – |
| S02 emotional, untreated | 38.1 | 37.4 | 45.1 | 58.9 | 0.754 | 39 / 12.0 | 4.24 | 6.02 | 4204 | 493 (2.5) | 2.36 | 0.56 |
| S03 physical (NE-dominant) | 29.7 | 29.5 | 41.8 | 57.8 | 0.798 | 37 / 2.0 | 4.59 | 14.69 | 6782 | 521 (2.3) | 9.43 | 0.57 |
| S04 neurological (SAH) | 19.3 | 19.0 | 35.7 | 56.4 | 0.884 | 38 / 0.9 | 4.21 | 16.77 | 11011 | 537 (2.5) | 19.6 | 0.73 |
| S05 phaeochromocytoma | 24.1 | 20.6 | 24.3 | 54.0 | 0.948 | 41 / 8.2 | 3.80 | 7.84 | 11722 | 523 (2.5) | 18.9 | 0.89 |
| **S06 pure NE (falsification 1)** | **57.9** | **57.9** | 60.1 | 61.9 | **0.041** | 35 / 0.0 | 4.85 | 6.93 | 282 | 427 (2.2) | 0.00 | 0.39 |
| **S07 pertussis (falsification 2)** | **60.4** | **60.4** | 61.5 | 62.0 | **0.014** | 34 / 0.1 | 4.85 | 5.71 | 160 | 423 (1.5) | 0.00 | **0.003** |
| S08 dobutamine | 67.5 | 33.2 | 40.8 | 57.8 | 0.900 | 40 / **33.3** | 4.24 | 7.77 | 6334 | 516 (3.1) | 8.12 | 0.62 |
| S09 levosimendan | 64.4 | 33.4 | 43.0 | 58.2 | 1.015 | 45 / **39.2** | 4.24 | 7.13 | 5731 | 524 (2.1) | 8.28 | 0.79 |
| **S10 milrinone (falsification 3)** | 65.4 | 33.5 | 40.4 | 57.7 | **0.910** | 42 / 37.1 | 4.24 | 8.17 | 5338 | 520 (2.5) | 9.08 | 0.67 |
| S11 adrenaline | 22.8 | 21.7 | 40.4 | 57.4 | 0.950 | 39 / 7.9 | 3.79 | 9.09 | 8280 | 515 (1.5) | 7.64 | 0.67 |
| S14 esmolol, no obstruction | 24.3 | 17.6 | 38.7 | 56.8 | 0.876 | 39 / 3.8 | 2.85 | 6.59 | 9446 | 516 (2.6) | 8.92 | 0.68 |
| S15 esmolol, with obstruction | 23.2 | 16.9 | 38.6 | 56.8 | 0.885 | 123 / 7.8 | 2.02 | 8.21 | 10041 | 516 (2.6) | 8.89 | 0.68 |
| S16 obstructive, untreated | 37.7 | 30.0 | 41.2 | 57.7 | 0.829 | 123 / 41.5 | 2.49 | 7.81 | 7738 | 510 (2.3) | 6.48 | 0.65 |
| S18 obstructive + IABP | 36.8 | 27.3 | 37.5 | 56.6 | 0.923 | 123 / **48.6** | 2.17 | 9.65 | 9092 | 519 (3.8) | 11.2 | 0.68 |
| S19 obstructive + Impella | 30.4 | 28.6 | 35.4 | 56.2 | 0.874 | 123 / **2.4** | 2.50 | 7.86 | 7938 | 530 (4.1) | 13.6 | 0.63 |
| S22 premenopausal (paired with S02) | **51.4** | 47.0 | 51.8 | 60.9 | **0.267** | 39 / 11.1 | 4.56 | 4.65 | 1445 | 456 (2.8) | 0.05 | 0.52 |
| **S25 reverse gradient** | 54.5 | 52.7 | 57.6 | 61.6 | **nadir −0.514** | 0 / 0.0 | 4.78 | 2.53 | 176 | 424 (5.2) | 0.00 | – |
| S26 ramipril | 38.1 | 34.5 | 44.0 | **60.8** | 0.754 | 39 / 11.8 | 4.24 | 6.02 | 3959 | 491 (2.2) | 2.02 | 0.66 |
| S27 metoprolol | 38.1 | 34.9 | 44.9 | 55.1 | 0.754 | 39 / 12.2 | 4.24 | 6.02 | 4498 | 497 (2.7) | 2.66 | 0.31 |

Units: LVEF %, gradient mmHg, CO L/min, cTnI ng/mL, NT-proBNP pg/mL, QTc ms.
**All 28 arms integrated successfully** (`delta` = 0.5 h, 90 days), and the maximum
relative error of the quasi-steady-state cAMP is at most 4.6 % even including the
surge peak.

Two rows worth noting: **S18 (IABP) raises the 48-hour mean gradient from 41.5 to
48.6 mmHg while S19 (Impella) lowers it from 41.5 to 2.4 mmHg.** Both act through
the same afterload term, and that the choice of mechanical support is decided by the
gradient term rather than by the severity of shock is an output of the model.

---

## Scenarios (28, designed in pairs)

| # | Scenario | What it tests |
|---|---|---|
| S01 | unstimulated control | whether the baseline state is an exact steady state |
| S02 | emotional trigger, untreated | the reference trajectory |
| S03 | physical trigger (NE-dominant) | the agonist mix selects the phenotype |
| S04 | neurological trigger (SAH) | maximal surge |
| S05 | adrenaline-secreting phaeochromocytoma | a natural experiment on agonist identity |
| **S06** | **pure NE surge** | **falsification 1: ballooning must not occur** |
| **S07** | **emotional trigger + pertussis toxin** | **falsification 2: Gi blockade must abolish it** |
| S08 | dobutamine 5 µg/kg/min | the cAMP pathway (harm axis) |
| S09 | levosimendan 0.1 µg/kg/min × 24 h | the non-cAMP pathway |
| **S10** | **milrinone 0.5 µg/kg/min** | **falsification 3: cAMP downstream of the receptor must be harmful too** |
| S11 | adrenaline 5 µg/min | the very agonist that makes the disease |
| S12 | noradrenaline 5 µg/min | SVR↑ via α1, but it adds β1 drive too |
| S13 | phenylephrine | afterload only, without touching β receptors |
| **S14** | **esmolol, no obstruction** | **sign pair (harm expected)** |
| **S15** | **esmolol, with obstruction** | **sign pair (benefit expected)** |
| S16 | obstructive TTS untreated | the control for S15 |
| S17 | phenylephrine in obstruction | reducing the gradient with afterload |
| S18 | IABP in obstructive TTS | afterload reduction → worse gradient (harm expected) |
| S19 | Impella in obstructive TTS | antegrade flow bypassing the apex |
| S20 | dobutamine in obstructive TTS | when two harm pathways overlap |
| S21 | levosimendan in obstructive TTS | sensitisation enlarges the gradient too |
| S22 | same stimulus, premenopausal | paired with S02 (the oestradiol benefit) |
| S23 | anxiolysis at arrival | blocking the afferent stimulus |
| S24 | QT-prolonging co-medication | repolarisation reserve |
| S25 | inverted receptor gradient | the reverse takotsubo prediction |
| S26 | ramipril 5 mg once daily | 1-year recurrence |
| S27 | metoprolol 50 mg twice daily | 1-year recurrence (paired) |
| S28 | apixaban 5 mg twice daily | apical thrombus |

---

## Usage

```r
# Required packages: mrgsolve (>= 1.0). shiny/ggplot2/DT are needed only for the dashboard.
setwd("takotsubo-syndrome")
source("tts_mrgsolve_model.R")

b  <- baseline("emo")                       # derive the baseline with a 14-day unstimulated burn-in
o  <- mrgsim_df(zero_re(b$mod), end = 2160, delta = 0.5)
summarise_arm(transform(o, scenario = "emo", label = "emotional trigger"))

summarise_all(run_all())      # all 28 scenarios
threshold_scan()              # where the threshold is (a model output)
falsification()               # the three falsification tests
bblocker_sign()               # the conditionality of the β-blockade sign
matched_drive()               # dobutamine vs levosimendan at matched drive
headroom()                    # per-mechanism shares of the EF deficit (measured)
recurrence_1y()               # ACEi vs β blockade, 1-year recurrence
thrombus_arms()               # with and without anticoagulation
population(200)               # virtual cohort
objective(log(fit_pars), TRUE) # per-anchor residuals
calibrate()                   # recalibration (takes several minutes)

# rendering the map
# dot -Tsvg tts_qsp_model.dot -o tts_qsp_model.svg
# dot -Tpng -Gdpi=150 tts_qsp_model.dot -o tts_qsp_model.png

# dashboard
# shiny::runApp("tts_shiny_app.R")
```

---

## Map structure (26 clusters)

1 triggers and susceptibility · 2 brain-heart axis · 3 sympathoadrenal output and
catecholamine kinetics · 4 receptor field (density gradient) ·
**5 threshold: β2 stimulus trafficking** · 6 cAMP/PKA ·
7 calcium and calcium overload · 8 myofilaments and segmental contractility ·
9 energy metabolism · mitochondria · 10 coronary microvasculature ·
11 membrane injury and troponin release · 12 myocardial oedema and transient
inflammation · **13 segmental mechanics → ballooning geometry (emergent)** ·
14 dynamic LVOT obstruction · 15 systemic haemodynamics · 16 neurohormones ·
17 electrophysiology (QT · TdP) · 18 apical thrombus and embolism · 19 recovery ·
20 drug A: cAMP-pathway inotropes (harm axis) · 21 drug B: non-cAMP contractility ·
22 drug C: β blockade and α1 agonism · 23 drug D: RAAS · diuresis · anticoagulation ·
24 mechanical circulatory support · 25 clinical endpoints · InterTAK · recurrence ·
**26 inference structure: counterfactual terms**

---

## ⚠️ Disclaimer

A semi-quantitative QSP model for education and research. It was assembled from the
public literature but has not been independently validated or certified, and **must
not be used for clinical decision-making, prescribing, or regulatory submission.**
The parameters are illustrative approximations, and fitting and validation against
real patient data would be needed separately. For the licence see the repository's
[`LICENSE`](../LICENSE).
