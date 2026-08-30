# Refactor notes — `dengue/denv_mrgsolve_model.R`

Compound: the antiviral (AV), per the existing
`dengue | Drugs (AV) | Redirect concentration (clean single site)` row in
`driver-patches/data/compound_perturbation_census.md`. The model's own
comments and parameter descriptions identify it only generically as "the
antiviral" (`kaAV`/`kelAV`/`VdAV`/`EMAXAV`/`EC50AV` block under
`// ---------- drugs ----------------------------------------------------------`);
no specific drug name is given anywhere in the file, and the README was not
treated as authoritative for this (per the fork guide's "README is context,
not ground truth" rule). **Scope of this refactor**: only the AV compound's
PK compartments (`AVD`/`AVC`) and its effect term (`av_block`) were touched.
Paracetamol (APAP) and methylprednisolone (STER) — the model's other two
drugs — are byte-for-byte identical to the original in
`denv_mrgsolve_model_refactored.R`.

## Archetype determined

**Bespoke: depot + single compartment, where the compartment is integrated
directly as a concentration** — a hybrid of the guide's archetype 3 (depot +
central, linear) and the "compartment IS the concentration" quirk documented
for tocilizumab in `sepsis/sep_refactor_notes.md`.

The original:

```
$PARAM @annotated
kaAV     : 1.10    : antiviral absorption (1/h)
kelAV    : 0.077   : antiviral elimination (1/h)
VdAV     : 92.0    : antiviral volume of distribution (L)
EMAXAV   : 0.965   : maximal block of viral production (-)
EC50AV   : 0.42    : antiviral EC50 (mg/L)

$CMT @annotated
AVD  : antiviral gut depot (mg)
AVC  : antiviral plasma concentration (mg/L)

$ODE
double av_block = EMAXAV * AVC / (AVC + EC50AV);
...
dxdt_AVD  = -kaAV * AVD;
dxdt_AVC  = kaAV * AVD / VdAV - kelAV * AVC;
```

`AVC` is not an amount that PD divides by a volume — it is integrated
directly in mg/L, with `VdAV` appearing only inside the depot-to-
concentration conversion `kaAV*AVD/VdAV`, and `kelAV` acting as an
independent first-order rate constant on the concentration itself (not a
`CL/V` ratio). This does not match archetype 3's literal template
(`double C_TCZ = CENT_TCZ / V1_TCZ;`) — applying that division would double-
divide by volume and change the model's numeric output. It also isn't quite
archetype 1 (there is a real depot with its own absorption ODE, unlike
tocilizumab's direct-dose case in sepsis). So this is deliberately not
forced into either named archetype; it is renamed to convention and
isolated from PD, per the guide's "none of these fit" fallback.

## Renaming

| Original | Refactored |
|---|---|
| `AVD` | `GUT_AV` |
| `AVC` | `CENT_AV` |
| `kaAV` | `KA_AV` |
| `kelAV` | `KEL_AV` |
| `VdAV` | `V1_AV` |
| `EMAXAV` | `EMAX_AV` |
| `EC50AV` | `EC50_AV` |
| `av_block` | `EFFECT_AV` |

Added (new, not present in the original, purely to complete the named Hill
interface): `GAMMA_AV = 1` — the original's `av_block` had no explicit Hill
exponent (implicit gamma = 1).

All parameter *values* are copied verbatim from the original (`KA_AV=1.10`,
`KEL_AV=0.077`, `V1_AV=92.0`, `EMAX_AV=0.965`, `EC50_AV=0.42`) — nothing
invented, nothing dropped.

```
$PARAM @annotated
KA_AV    : 1.10    : antiviral gut absorption rate constant (1/h)
KEL_AV   : 0.077   : antiviral elimination rate constant (1/h)
V1_AV    : 92.0    : antiviral volume of distribution (L)
EMAX_AV  : 0.965   : maximal block of viral production and NS1 secretion (-)
EC50_AV  : 0.42    : antiviral EC50 (mg/L)
GAMMA_AV : 1.0      : antiviral Hill coefficient (-) [new; original had none]

$CMT @annotated
GUT_AV  : antiviral gut depot (mg) [was AVD]
CENT_AV : antiviral plasma concentration (mg/L) [was AVC; concentration by construction, not amount/volume]

$ODE
double C_AV = CENT_AV;   // exposed concentration = compartment, undivided
double EFFECT_AV = EMAX_AV * pow(C_AV, GAMMA_AV)
                  / (pow(EC50_AV, GAMMA_AV) + pow(C_AV, GAMMA_AV));
...
dxdt_GUT_AV  = -KA_AV * GUT_AV;
dxdt_CENT_AV = KA_AV * GUT_AV / V1_AV - KEL_AV * CENT_AV;
```

`EFFECT_AV` replaces `av_block` at both its original use sites: the viral
production term `p_eff = PPROD0 * (...) * (1.0 - EFFECT_AV)` and the NS1
secretion term `dxdt_NS1 = pNS1 * (1.0 - EFFECT_AV) * INF - ...`. No other
line in `$MAIN`/`$ODE`/`$TABLE` was touched.

`C_AV` and `EFFECT_AV` were additionally appended to `$CAPTURE` (the
original's `$CAPTURE ABH0 V0 RCRYS RCOLL` becomes
`$CAPTURE ABH0 V0 RCRYS RCOLL C_AV EFFECT_AV`) purely so they are
discoverable as outputs via the qspserver `/model_manifest` and
`/run_simulation` endpoints — see "qspserver compatibility" below.

Outside the DSL block, the R driver's dosing (`CMT_AVD <- 38`, used by
`build_events()` to dose the `s$av` scenarios) addresses the compartment by
**numeric position**, not by name, and compartment order is unchanged by
this refactor (`GUT_AV` is still compartment 38, same as `AVD` was) — so no
functional change was needed there. The one accompanying comment
(`## 24 = PLT, 38 = AVD (antiviral depot), ...`) was updated to say `GUT_AV
(antiviral depot, was AVD)` purely for readability; this is a comment only,
not a semantic change.

## Hill interface: rename, not a fit

The original's `av_block = EMAXAV * AVC / (AVC + EC50AV)` is already an
exact Hill ratio with an implicit `gamma = 1` (a plain Emax/EC50 form, no
ODE-solved receptor kinetics behind it) — so per the guide, this is a
rename, not a refit. `EMAX_AV` and `EC50_AV` carry the original's exact
values; `GAMMA_AV = 1` is added explicitly since the original had no Hill
exponent term. `EFFECT_AV` is computed from the exact same algebraic
expression as `av_block` (renamed variables only), so there is no
approximation to introduce error.

## qspserver compatibility

- `KA_AV`, `KEL_AV`, `V1_AV`, `EMAX_AV`, `EC50_AV`, `GAMMA_AV` are declared
  in `$PARAM`, so `/model_manifest` lists them with their original default
  values (confirmed: `KA_AV=1.1`, `KEL_AV=0.077`, `V1_AV=92`,
  `EMAX_AV=0.965`, `EC50_AV=0.42`, `GAMMA_AV=1` all appear in the manifest
  response).
- `C_AV` and `EFFECT_AV` are computed in `$ODE` (they depend on the
  compartment state, so — per the precedent set in
  `thyroid-eye-disease/ted_refactor_notes.md` — they cannot themselves be
  fixed `$PARAM` values). Both were added to `$CAPTURE` so they are
  discoverable as **outputs** via `/model_manifest`'s `outputPaths` and
  retrievable via `/run_simulation`'s `outputs` list — confirmed both ways
  (see Verification below). This is additive only; it does not appear in
  the original-vs-refactored diff below, which compares only the original's
  own 71 pre-existing `$CMT`/`$CAPTURE` names.
- The `.cpp` bare-DSL extraction step (`model_content` must be pure mrgsolve
  DSL, no R wrapper) **does** apply to this file — `denv_mrgsolve_model.R`
  uses the `code <- '...'; mcode_cache(...)` R-string-wrapper pattern. The
  quoted block was extracted mechanically (the text between the `code <- '`
  line and the closing `'` line, byte-for-byte) for both the original and
  the refactored file, used only as transient request payloads to the
  qspserver API, and deleted afterward — no `.cpp` file is left on disk in
  `dengue/`.

## A pre-existing defect found while verifying (not fixed, logged upstream)

Neither `denv_mrgsolve_model.R` nor the refactored file compiles as written
under the installed mrgsolve version (2.0.1). `$ODE` declares nine local
`double`s (`SV`, `CO`, `MAP`, `PP`, `Jv`, `Pser`, `Jser`, `Dser`, `Hct`) used
only for that section's own `dxdt_*` computations; `$TABLE` independently
recomputes the same physiology under distinctly `_o`-suffixed names
(`SV_o`, `CO_o`, ... `Hct_o`) and reports them via `capture SV = SV_o;` etc.
— i.e. the author already used different names specifically to avoid a
clash. Under mrgsolve 2.0.1 the clash happens anyway, because the compiler
auto-promotes every bare `$ODE` assignment to a reportable class member,
so `$TABLE`'s `capture SV = SV_o;` then tries to redeclare a member `$ODE`'s
own `double SV = ...;` already created:

```
130:11: error: redefinition of 'capture {anonymous}::Jv'
  130 |   capture Jv;
      |           ^~
66:10: note: 'double {anonymous}::Jv' previously declared here
   66 |   double Jv;
      |          ^~
```

(and identically for the other eight names.) Confirmed on the *original*,
untouched file, via the qspserver `mrgsolve_api` container — this is not
something the AV refactor introduced, and is completely unrelated to AV
(none of the nine colliding names are AV-related; they are all Starling-
exchange/haemodynamics locals). **For verification only**, both the
original and the refactored model text were given the identical,
minimal, non-numeric workaround of suffixing all nine `$ODE`-local names
`_ode` (and updating every downstream in-`$ODE` reference to match), leaving
`$TABLE`'s `_o`-suffixed recomputation and `capture` lines completely
untouched. This workaround was applied symmetrically to scratch copies of
both models and does not touch any numeric value or `$TABLE` logic, so it
cannot bias the comparison. Neither `denv_mrgsolve_model.R` nor
`denv_mrgsolve_model_refactored.R` on disk contains this workaround — it
existed only in temporary, unsaved request payloads sent to the qspserver
API. Logged as entry #35 in `translations/UPSTREAM_ISSUES.md` (note: a
concurrent agent's edit landed entry #34 in between while this one was being
appended, which is why this is #35 rather than #33 — see that file's own
entry ordering).

## Verification

**Method**: qspserver `mrgsolve_api` container at `http://localhost:8007`
(`POST /model_manifest`, `POST /run_simulation`), per the guide's mandatory
protocol — local R/mrgsolve was not used.

**Scenario tested**: the AV dosing block from the original file's own
`S10_antiviral_day2_illness`/`S11_antiviral_at_present` scenarios (both
resolve to the same dosing under the original's own definitions: `av = 48`
and `av = T_PRESENT = 48` are numerically identical), reproduced exactly —
`450` mg boluses into compartment 38 (`AVD`/`GUT_AV`), starting at t = 48 h,
repeated every 12 h for 14 doses total (`seq(0,13)*12`, i.e. `ii=12,
addl=13` in the API's NM-TRAN-style `DoseSpec`) — plus `ABH0 = 55`
(secondary-infection titre, also part of both scenarios), run 0–336 h at
`delta = 0.5` (674 time points), matching `run_scenario()`'s own defaults.

**One deliberate simplification, disclosed**: the qspserver API's `dosing`
field (mrgsolve engine) supports only NM-TRAN-style bolus/infusion
`DoseSpec` records; it has no mechanism for the WHO fluid-ladder's
time-varying `RCRYS` **parameter** step-schedule that `S10`/`S11` also set
via `crys = 1.00` (implemented in the original's own R driver as a
`data_set` of `evid=2` parameter-change rows, a pattern the generic API does
not expose). Since `RCRYS`/the fluid ladder is completely orthogonal
machinery — it enters only the Starling-exchange/haemodynamics equations in
section 7, nowhere near AV's own compartments or effect term — and AV's
effect enters far upstream (viral production and NS1 secretion, section 2),
this was run without the fluid ladder (`RCRYS` left at its `$PARAM` default
of 0) rather than force an approximation of it through the dosing API. This
still exercises AV's effect on the **entire** downstream model (virus, NS1,
IFN, TNF, immune complexes, glycocalyx, both Starling beds, haemodynamics,
liver, fever, haematology) exactly as the original's own AV scenarios do,
so it remains a complete test of the refactored block; the fluid ladder
itself is untouched code, run through neither original nor refactored, and
not part of what this refactor changed.

Every `$CAPTURE`d output shared between the two models — all 45 `$CMT`
compartments and 26 `$TABLE`-computed captures (`TGT, ECL, INF, V, NS1,
ABH, ABN, PBL, ASC, CTL, IFN, TNF, IL10, VEGF, CHYM, GLX, LPX, VP, VI, PRP,
PRI, VSER, RBCV, PLT, MKC, APLT, FIB, WBC, HEP, AST, ALT, PGE, TEMP, SVRR,
HR, LAC, COLL, AVD/GUT_AV, AVC/CENT_AV, APAP, STER, CIN, COUT, NAUC, ICAUC,
ILLDAY, EADE, NEUTF, SIGMA, Jv, Jlymph, Jser, Dser, NETSER, Pcap, Pint,
Pser, COPp, COPi, Hct, HCTRISE, SV, CO, MAP, SBP, DBP, PP, TPROT, ALBUM,
EFFUS, LOG10V, BLEEDIX, WARNSIGN, SHOCK, FLUIDBAL, DO2out` — 71 outputs in
total) was compared point-by-point across all 674 time points.

**Result: exact match, max absolute deviation = 0.0 and max relative
deviation = 0.0 for every output, every time point.** This is a rename with
no change to the underlying arithmetic, so it exceeds the guide's
"near-exact match" expectation for pure structural reorganization — a
bit-for-bit match, as expected. The antiviral effect itself is far from
trivial over this run: `EFFECT_AV` (checked in the refactored output) rises
from ~0 pre-dose to 0.80 by t=60h, 0.84 by t=72h, and ~0.90 by t=150h, and
viraemia (`V`) collapses from a peak of ~2.1e6 copies/mL at t=48h to
~99 copies/mL by t=96h under treatment — confirming the comparison actually
exercises the refactored PK/effect block rather than passing on a
degenerate no-op case.

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`,
`dengue | Drugs (AV)` row.
