# Refactor notes — `refeeding-syndrome/rfs_mrgsolve_model.R`

Scope: the two rows `driver-patches/data/compound_perturbation_census.md`
carries for this file, both classified **"Delete PK compartment;
concentration is itself the state"**: **Magnesium (MG)** and
**Potassium**. Phosphate, thiamine, glucose/insulin, calcium/PTH/vitamin
D/FGF23, and every other state in this 48-compartment model are
untouched — they carry no census row and are out of scope.

## Is Magnesium/Potassium really dosed here? (the task's specific check)

Checked explicitly, per the task's instruction and the precedent set in
`wilsons-disease/wd_refactor_notes.md` (Zinc) and
`hypercalcemia-of-malignancy/mah_refactor_notes.md` (PTH, found to be
*not* dosed and corrected). Both electrolytes here are genuinely
dual-purpose: endogenous state **and** a real, exercised IV-infusion
dosing route, but with **no separate PK compartment of their own** — the
census classification's structural premise is correct, not a
misclassification to fix (contrast `toxic-alcohol-poisoning/
tap_refactor_notes.md`'s FOM row and `clostridioides-difficile-infection/
cdi_refactor_notes.md`'s BEZ row, both "Delete PK compartment" rows that
*were* wrong because a real separate PK backbone existed underneath).

- `MGE` (ECF magnesium, mmol) and `KE` (ECF potassium, mmol) are the
  model's own physiological ECF pools — fed by dietary absorption
  (`JMgabs`/`JKabs` from `MGGUT`/`KGUT`), cellular exchange
  (`JMgin`/`JMgout`, `JKin`/`JKout`), new-tissue uptake during refeeding
  (`JMgnew`/`JKnew`), renal excretion (`ExMg`/`ExK`), and GI losses
  (`gi_Mg`/`gi_K`) — **and** by a genuine continuous IV-infusion term,
  `ivMg`/`ivK`, computed from the `MG_DOSE`/`K_DOSE` `$PARAM`s (mmol/kg/d)
  and added directly into `dxdt_MGE`/`dxdt_KE`:
  ```c
  ivK = K_DOSE*BW/24.0;
  if(MG_BOLUS > 0.5) ivMg = (fmod(SOLVERTIME,24.0) < 2.0) ? MG_DOSE*BW/2.0 : 0.0;
  else               ivMg = MG_DOSE*BW/24.0;
  ...
  dxdt_KE  = JKabs  + ivK  + JKout  - JKin  - ExK  - gi_K  - JKnew;
  dxdt_MGE = JMgabs + ivMg + JMgout - JMgin - ExMg - gi_Mg - JMgnew;
  ```
  This is exactly the guide's "Edge case — continuous infusion input, not
  an event" (Part 2): no `ev()`/`data_set()` anywhere in the file for
  either electrolyte (confirmed by exact-name grep across the whole
  file — in fact **no compound in this file** uses event-based dosing;
  phosphate and thiamine are driven the identical way via `P_DOSE`/
  `TH_IV`/`TH_PO`, also out of scope here).
- `MG_DOSE`/`K_DOSE` are not dead parameters: they are varied across 6 of
  the file's own 16 scenarios (`01 NICE 2006 cautious ramp`,
  `02 ASPEN 2020 faster ramp`, `04 full feed + electrolytes only`,
  `06 full feed + both`, `14 Mg as a 2 h bolus`, `15 Mg as a 24 h
  infusion`, `16 potassium without magnesium` — the last built
  specifically to isolate one electrolyte's dosing from the other's).
- There is **no separate "drug" compartment or concentration** anywhere
  in the file distinct from `MGE`/`KE` themselves — unlike Zinc in
  Wilson's disease (its own `GUT_ZN`/`CENT_ZN` pair, structurally
  identical to the other three drugs in that file) or fomepizole/
  bezlotoxumab in the two corrected rows above (each had a real
  depot/central/peripheral PK backbone feeding a *separate* concentration
  variable). Here the exposed concentration (`MGser`/`Kser`, computed
  once per `$ODE` substep as `posv(MGE)/Vecf` / `posv(KE)/Vecf`) **is**
  the disease-state concentration used by every downstream physiological
  equation, dosed or not.

So both rows are **confirmed correct as classified** ("Delete PK
compartment; concentration is itself the state") — not corrected — and
both are genuinely, actively dosable. Recorded back into the census
unchanged in method, with the redirect-site/unit/summary columns filled
in (see bottom of this file).

## Archetype: bespoke ("None of these fit")

Neither compound fits Archetypes 1–4 as written — all four assume a
compound with its own PK compartment(s) separate from the disease state
it acts on (`$CMT CENT_<STEM>` etc., a volume, a clearance). Here `MGE`/
`KE` are the disease state, shared by every physiological flux in
sections 7–9 of `$ODE`, not a dedicated drug compartment. Per the guide's
"None of these fit" fallback: renamed to the convention, and every
already-distinct, already-named downstream physiological term each
electrolyte feeds is pulled out as its own `EFFECT_<STEM>_*` — mirroring
`diabetic-ketoacidosis/dka_refactor_notes.md`'s treatment of insulin
(four independently-named `EFFECT_INS_*` axes for one compound with
several genuine physiological effects, rather than one artificial
composite term). The self-referential homeostatic terms each electrolyte
uses to regulate **its own** compartment (`depMg`/`depK`, `JMgin`/`JKin`,
`JMgout`/`JKout`, `TmMg`, `romk`) are *not* renamed to `EFFECT_` — they
are the electrolyte's own physiology, not its effect on something else,
and were left exactly as the original computed them (only the bare
`Kser`/`MGser` reads inside them were renamed to `C_K`/`C_MG`).

### Naming applied

| Role | Original | Refactored |
|---|---|---|
| The exposed concentration (potassium) | `Kser` (`$ODE`-local, `posv(KE)/Vecf`) | `C_K` — same line, same formula |
| The exposed concentration (magnesium) | `MGser` (`$ODE`-local, `posv(MGE)/Vecf`) | `C_MG` — same line, same formula |
| Mg's effect on PTH secretion | `fMgPTH` (normalised Michaelis-Menten ratio) | `EFFECT_MG_PTH` — rewritten in canonical Hill form (see below) |
| Mg's effect on renal K wasting | `fMgK` (linear deficit multiplier) | `EFFECT_MG_KWASTE` — rename only, formula unchanged |
| Mg's effect on QTc | inline term inside the `QTc` sum | `EFFECT_MG_QTC` — extracted, formula unchanged |
| K's effect on QTc | inline term inside the `QTc` sum | `EFFECT_K_QTC` — extracted, formula unchanged |
| K's effect on aldosterone secretion | inline `SQ(Kser/KSER0)` inside `ald_t` | `EFFECT_K_ALD` — extracted, formula unchanged |
| Mg's own dosing input | `MG_DOSE`, `MG_BOLUS` | unchanged — a dosing-rate input, not one of the guide's PK roles (CL/V/KA/Q); same treatment as `INS_IV` in `dka_refactor_notes.md` |
| K's own dosing input | `K_DOSE` | unchanged, same reasoning |

`MGE`/`KE`/`MGGUT`/`KGUT`/`MGI`/`KI` (the `$CMT` compartments themselves)
are **not** renamed to `CENT_MG`/`CENT_K` etc. — they are shared disease
state read and written by many equations that have nothing to do with
dosing (dietary absorption, cellular exchange, growth uptake, GI loss,
renal excretion), and relabelling them as if they were a dedicated "drug
compartment" would misrepresent that. Only the *derived, `$ODE`-local*
concentration and effect terms — the genuine "pluggable interface"
points — are renamed, per the guide's "exposes exactly one concentration
variable that PD equations read" requirement.

## The Hill interface

### `EFFECT_MG_PTH`: already Hill-shaped once rewritten — a rename, not a refit

The original:
```c
double fMgPTH = (MGser/(MGPTH50 + MGser))/(MGSER0/(MGPTH50 + MGSER0));
if(fMgPTH > 1.05) fMgPTH = 1.05;
```
is a Michaelis-Menten ratio normalised to 1 at `MGSER0` (the reference
magnesium), not literally written as `Emax*C/(EC50+C)` — but it is
**algebraically identical** to that canonical form for a specific,
derivable `Emax`. Setting `EC50_MG_PTH = MGPTH50` (renamed, same value,
0.22 mmol/L — used nowhere else in the file, confirmed by exact-name
grep) and `GAMMA_MG_PTH = 1.0` (new, disclosed — no explicit Hill
coefficient in the original, same treatment the guide sanctions for
Zinc/TTM in `wilsons-disease/wd_refactor_notes.md`), the required `Emax`
is the value that makes `Emax*C/(EC50+C) = 1` at `C = MGSER0`:

```
Emax * MGSER0/(EC50_MG_PTH + MGSER0) = 1
=>  Emax = (EC50_MG_PTH + MGSER0) / MGSER0
```

Substituting back: `Emax*C/(EC50+C) = [(EC50+MGSER0)/MGSER0] * C/(EC50+C)
= [C/(EC50+C)] / [MGSER0/(EC50+MGSER0)]` — exactly the original's ratio.
`EMAX_MG_PTH` is therefore declared as a **derived constant** (computed
once in `$MAIN` from `EC50_MG_PTH`/`MGSER0`, matching this file's own
existing convention for every other derived constant — `TMPMAX`, `KBRES`,
`CA50`, etc. — rather than a hardcoded decimal literal, to avoid any
rounding risk):

```c
EMAX_MG_PTH = (EC50_MG_PTH + MGSER0)/MGSER0;   // = 1.258823529... in $MAIN
```
```c
double EFFECT_MG_PTH = EMAX_MG_PTH * pow(C_MG, GAMMA_MG_PTH)
                      / (pow(EC50_MG_PTH, GAMMA_MG_PTH) + pow(C_MG, GAMMA_MG_PTH));
if(EFFECT_MG_PTH > 1.05) EFFECT_MG_PTH = 1.05;
```
Verified numerically, not just algebraically: `MGSER`/`C_MG` and every
downstream output (`PSER`, `KSER`, `MGSER`, `CASER`, `CAXP`, `QTC`,
`MORT`, `TBP`, …) matched the untouched original to floating-point
precision across three scenarios (see Verification below) — the rewrite
introduced no numeric drift.

### `EFFECT_MG_KWASTE`, `EFFECT_MG_QTC`, `EFFECT_K_QTC`: genuinely linear, not Hill-shaped — preserved as linear, not force-fit

`fMgK = 1.0 + BETAMGK*posv(1.0 - MGser/MGSER0)` (renal-K-wasting
amplification) and the `AK*posv(KSER0-Kser)` / `AMG*posv(MGSER0-MGser)`
QTc-prolongation terms are all **first-order in deficit, unbounded, with
no saturation on the electrolyte side** — the same shape the guide
explicitly instructs not to force into `Emax*Cᵞ/(EC50ᵞ+Cᵞ)` (see
`wilsons-disease/wd_refactor_notes.md`'s treatment of D-Penicillamine/
Trientine). `EFFECT_MG_KWASTE`/`EFFECT_MG_QTC`/`EFFECT_K_QTC` are exposed
exactly as the original computed them: renamed (`EFFECT_MG_KWASTE`) or
extracted into a named local (`EFFECT_MG_QTC`, `EFFECT_K_QTC`) with the
formula byte-for-byte unchanged.

### `EFFECT_K_ALD`: genuinely quadratic, not Hill-shaped — preserved as quadratic

`SQ(Kser/KSER0)` inside `ald_t` (potassium's drive of aldosterone
secretion) is a plain squared ratio with no `EC50`/saturation term at
all — extracted into `EFFECT_K_ALD = SQ(C_K/KSER0)` unchanged, for the
same "don't invent a shape the original doesn't have" reason.

## `$PARAM` vs `$CAPTURE` for the interface variables

Tried first per the guide's literal instruction (`C_MG = 0` etc. directly
in `$PARAM`, assigned in `$ODE`). Confirmed via `POST /model_manifest`
that this does not compile under mrgsolve 2.0.1 — `$PARAM` values are
read-only inside `$ODE`, the same constraint documented in
`wilsons-disease/wd_refactor_notes.md` and every other model in this
refactor batch. Followed that established, actually-compiling precedent
instead: `C_MG`, `C_K`, `EFFECT_MG_PTH`, `EFFECT_MG_KWASTE`,
`EFFECT_MG_QTC`, `EFFECT_K_QTC`, `EFFECT_K_ALD` are plain `double` locals
computed in `$ODE`, listed in a new `$CAPTURE` block appended after the
original's own `$TABLE`. Confirmed via `POST /model_manifest` on the
refactored DSL: all seven appear in `outputPaths` (see Verification
below) — fully discoverable through the output side of the manifest, not
the overridable `parameters` side. `MG_DOSE`/`K_DOSE` (the actual
externally-drivable inputs) were already in `$PARAM` in the original and
needed no change to be manifest-discoverable.

## The dose-instant reporting artifact: does not apply here

Checked per the guide's dedicated section. That artifact requires a
discrete `ev()`/dataset bolus event whose instant produces a duplicate
report row computed from pre-dose state. **This file has no `ev()`/
`data_set()` dosing anywhere, for any compound** — every dose (`ivP`,
`ivK`, `ivMg`, `ivTh`, `poTh`) is a continuous, `$PARAM`-driven flux term
recomputed every `$ODE` substep from `SOLVERTIME`/`day`, including the
`MG_BOLUS` branch (`fmod(SOLVERTIME,24.0) < 2.0`, a continuous piecewise
function, not a discrete event). There is no dose instant to produce a
stale-state artifact at. Verified directly: comparing the `MG_BOLUS=1`
scenario (which exercises the 2-hour-on/22-hour-off piecewise branch)
point-by-point across the full grid showed no elevated deviation at any
`fmod(SOLVERTIME,24.0)` boundary crossing — see Verification below.

## When the original doesn't compile at all

Confirmed via `POST /model_manifest` on the untouched original's own DSL
(extracted verbatim from `rfs_code <- '...'`): **HTTP 500**, a build
failure with nothing to do with either electrolyte. mrgsolve 2.0.1's
`$ODE` preprocessor keeps `double` on only the first declarator of a
multi-declarator line and drops it from every subsequent
comma-separated name that also carries its own initializer — six
separate lines in this file's `$ODE` hit this:

```c
double fm = posv(FM), ffm = posv(FFM);
double Vg = FVG*BW, Vp = FVP*BW;
double g_cho = kcal_h*cho_f/4.0, g_fat = kcal_h*fat_f/9.0, g_pro = kcal_h*pro_f/4.0;
double qP = 0.40 + 0.60*quality, qK = 0.70 + 0.30*quality;
double qMg = 0.45 + 0.55*quality, qTh = quality;
double ivP = 0.0, ivK = 0.0, ivMg = 0.0, ivTh = 0.0, poTh = 0.0;
```

Same defect class as `diabetic-ketoacidosis/dka_mrgsolve_model.R`'s
issue #37, independently present here — logged as a new entry,
**UPSTREAM_ISSUES.md #124** (not fixed in the original, per the
never-edit-upstream rule). Per the guide's now-settled policy (rather
than #37's earlier in-memory-only workaround, superseded going forward):
the delivered `rfs_mrgsolve_model_refactored.R` carries the syntax-only
fix forward directly — each of the six lines split into one `double`
statement per declarator, same initializer expressions, same values,
same order, nothing else changed on those lines. Confirmed via
`POST /model_manifest` that the patched-and-renamed DSL compiles (HTTP
200), and the Verification below used an identically-patched
*(syntax-fix only, no renaming)* copy of the original as the comparison
baseline, exactly as the guide specifies.

Note the fix is narrower than issue #37's description might suggest:
several **pure declaration lists with no initializers** appear earlier in
the same `$ODE` block (`double kcal_d, cho_f, fat_f, pro_f, quality;`,
and the `$GLOBAL` derived-constants block, e.g. `double VECF0, VG0, VP0,
...;`) and compile cleanly — confirmed by their absence from the compiler
error list. Only lines mixing an initializer with a following
comma-separated declarator are affected.

## Verification

Extracted the bare DSL text from both `rfs_mrgsolve_model.R` (`rfs_code
<- '...'`, patched only with the six-line syntax fix above, no
electrolyte renaming) and `rfs_mrgsolve_model_refactored.R` (same
variable name, syntax fix + full renaming), and ran three of them through
qspserver `mrgsolve_api`, requests spaced ≥2.2 s apart:

1. **Scenario `01 NICE 2006 cautious ramp`'s own dosing/refeeding block**:
   `PHASE=2, KCAL_START=10, KCAL_GOAL=30, ADV_FRAC=0.20, ADV_START=2,
   P_DOSE=0.5, K_DOSE=2.5, MG_DOSE=0.3, TH_IV=200, TH_LEAD=0.5,
   TH_DAYS=10, RESCUE=1`.
2. **Scenario `16 potassium without magnesium`** (built by the file's own
   author specifically to isolate one electrolyte's dosing from the
   other's — a direct test of independent driveability): `PHASE=2,
   DIURETIC=0.6, KCAL_START=20, KCAL_GOAL=20, ADV_START=99, MG_DOSE=0.0,
   K_DOSE=2.0`.
3. **Scenario `14 Mg as a 2 h bolus`** (exercises the `MG_BOLUS`
   piecewise branch): `PHASE=2, DIURETIC=0.6, KCAL_START=20,
   KCAL_GOAL=20, ADV_START=99, MG_DOSE=0.4, MG_BOLUS=1, K_DOSE=2.0`.

**Adaptation disclosed:** the original's own R driver (`rfs_run()`)
simulates a starvation phase first, resets several cumulative-hazard
states, then re-initialises phase 2 from phase 1's ending state — a
two-call, state-carrying workflow the stateless `/run_simulation`
endpoint cannot reproduce in one request (no mechanism to inject a custom
full-state initial condition; the same constraint documented in
`diabetic-ketoacidosis/dka_refactor_notes.md` for its `run_protocol()`
feedback loop). All three runs above instead start from the model's own
built-in `$MAIN` "HEALTHY INITIAL CONDITIONS" with `PHASE=2` set
directly, then apply each scenario's dosing/refeeding parameters for the
full 14-day window (`end=336, delta=0.25`, matching the original's own
phase-2 `mrgsim()` call). This exercises `MG_DOSE`/`K_DOSE`'s dosing
mechanism and the complete downstream `$ODE` system exactly as the
model's own scenarios do; only the starting point (healthy vs.
pre-starved) differs, and it is applied identically on both sides of the
comparison, so it cannot itself be a source of any divergence between the
two models. (One consequence, visible in the sanity spot-check below: `C_K`
rises above the reference range under scenario 1, because K is dosed on
top of a healthy — not depleted — baseline; this is an expected artifact
of the adapted starting condition, not a defect in either model.)

Every requested `$CAPTURE`d output shared by both files (`PSER`, `KSER`,
`MGSER`, `CASER`, `CAXP`, `BWkg`, `BMI`, `PICF`, `KICF`, `MGICF`, `TKACT`,
`QTC`, `MORT`, `WERN`, `OEDEMA`, `TBP` — 16 columns) was compared
point-by-point across the full 1,345-point time grid (`t=0` to `t=336h`,
`delta=0.25h`) for all three scenarios.

**Result:**
- Scenario 1 (NICE): max abs diff = **1.0e-4**, at a single point on
  `TBP` (t=68.25h, 22577.2312 vs 22577.2313) — consistent with the API's
  own output rounding (4 decimal places on a ~22,577-mmol quantity, i.e.
  ~4×10⁻⁹ relative), not a computational divergence; every other point on
  every other column matched to `0.0`.
- Scenario 2 (K without Mg): max abs diff = **0.0** across all 16
  columns, all 1,345 points — bit-identical.
- Scenario 3 (Mg 2 h bolus, exercising the piecewise `MG_BOLUS` branch):
  max abs diff = **1.0e-4**, again a single `TBP` point (t=207.25h),
  same output-rounding pattern; every other point matched to `0.0`.

This is the expected outcome for a pure structural reorganisation (no
Hill-refitting on the linear/quadratic terms; the one genuinely-rewritten
term, `EFFECT_MG_PTH`, is algebraically identical to the original, as
derived above and confirmed numerically) — **no tolerance was loosened to
reach this result.**

`/model_manifest` on the refactored DSL additionally confirmed
`C_MG`/`C_K`/`EFFECT_MG_PTH`/`EFFECT_MG_KWASTE`/`EFFECT_MG_QTC`/
`EFFECT_K_QTC`/`EFFECT_K_ALD` all appear in `outputPaths`, and
`MG_DOSE`/`K_DOSE`/`EC50_MG_PTH`/`GAMMA_MG_PTH` all appear in
`parameters`. Spot-checked non-degenerate, correctly-shaped values from
scenario 1's own run (t=48h): `C_MG = 1.0445` (above `MGSER0 = 0.85`) ⇒
`EFFECT_MG_PTH = 1.0398` (> 1, as expected — more magnesium, more PTH
secretion per this term); `MGSER` and `C_MG` are exact pass-throughs of
each other at every timepoint (`max abs diff = 0.0`), likewise `KSER`/
`C_K` — confirming the rename introduced no numeric drift, independent of
the column-level comparison above.

No scratch/debug artifacts were left in the repository — the DSL
extraction, the syntax-fix patching, and the comparison scripts ran
entirely from a session-local scratchpad directory outside the repo tree.

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`, both
`refeeding-syndrome` rows: perturbation method confirmed unchanged
("Delete PK compartment; concentration is itself the state" — the
structural classification was correct, not corrected), target/pathway
and redirect-site/unit/summary columns filled in for Magnesium (MG) and
Potassium.
