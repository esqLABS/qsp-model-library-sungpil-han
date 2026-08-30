# Refactor notes — `distal-renal-tubular-acidosis/drta_mrgsolve_model.R`

Calibration run: hydrochlorothiazide (HCTZ), the model's only add-on drug
(the disease's own alkali therapy — bicarbonate/citrate/potassium salts —
is not a "compound" in the driver-patch sense; HCTZ is the one exogenous
agent with its own PK block per `driver-patches/data/compound_perturbation_census.md`,
classified "Redirect concentration (clean single site)"). Only HCTZ's PK
block and its two downstream usages (renal potassium secretion, urinary
calcium) were touched. Every other compartment, parameter, and equation —
including the entire alkali-therapy and vitamin-D PK — is byte-for-byte
identical to `drta_mrgsolve_model.R` (see the diff-scoped summary at the
end of this file).

## Archetype

**Archetype 3, without peripheral compartment and without a depot
bioavailability parameter** (depot + central, linear elimination), applied
as a genuine reparameterization, not just a rename:

- The original tracks HCTZ with two compartments, `AG_hctz` (gut depot,
  mg) and `C_hctz` — but `C_hctz` is itself an ODE state holding
  **concentration directly** (mg/L), not an amount:
  ```
  dxdt_AG_hctz = -ka_hctz * AG_hctz;
  dxdt_C_hctz  = ka_hctz * AG_hctz / Vd_hctz - ke_hctz * C_hctz;
  ```
  This is mathematically a standard one-compartment depot model (bolus into
  a depot, first-order absorption, first-order elimination) written in its
  concentration form instead of its amount form — `d(C)/dt = ka·GUT/V -
  ke·C` is `d(CENT/V)/dt` for `dCENT/dt = ka·GUT - CL·C` with `CL = ke·V`.
- Refactored to the guide's convention: `GUT_HCTZ` (was `AG_hctz`, unchanged
  units/behavior) and `CENT_HCTZ` (an **amount**, was `C_hctz`), with
  `C_HCTZ = CENT_HCTZ / V1_HCTZ` computed each step — matching every other
  archetype-3 compound in this fork's naming convention, where `C_<STEM>`
  is always a derived quantity, never itself a state.
- `V1_HCTZ = Vd_hctz` (renamed, same value) and `CL_HCTZ = ke_hctz *
  Vd_hctz = 0.088 * 3.60 = 0.3168` (a pure reparameterization of the
  original's micro-constant `ke_hctz`, following exactly the rule the guide
  gives for archetype 2's `k10 -> CL/V1` rewrite — not an invented number).
  `KA_HCTZ` is a straight rename of `ka_hctz` (1.10 /h, unchanged).
- No `F_HCTZ` was added: the original has no bioavailability parameter for
  HCTZ (unlike bicarbonate/citrate/potassium, which do have explicit `F_*`
  parameters), so absorption is implicitly complete. Per the guide's rule
  never to invent a parameter the original didn't have, `GUT_HCTZ` feeds
  `CENT_HCTZ` at the full `KA_HCTZ * GUT_HCTZ` rate, exactly matching the
  original's `ka_hctz * AG_hctz` term with no added `F`.

This is why the "no peripheral, no depot-F" archetype-3 variant fits exactly
and needed no bespoke structure.

## Hill interface: rename, not a fit

The original computes one ratio, `hz_ = C_hctz / (EC50_hctz + C_hctz)`, and
uses it **directly, unscaled, at two separate points** in the disease
equations:

1. Renal potassium secretion (`Ksec`): `+ kK_hctz * hz_ * BW/24 * 0.05`
   (thiazide kaliuresis).
2. Urinary calcium (`UCa_h`): `* (1.0 - kUCa_hctz * hz_)` (thiazide
   hypocalciuric effect).

`hz_` is already exactly `C / (EC50 + C)` — a Hill function with an
implicit `Emax = 1` and no explicit Hill exponent (`gamma = 1`). Per the
guide's rule ("if the original's effect term is already this shape... pull
out `EMAX_<STEM>`, `EC50_<STEM>`, `GAMMA_<STEM>` as explicit named
parameters with the original's values"), this is a rename:

```
EMAX_HCTZ  = 1.0   // the ratio's implicit ceiling
EC50_HCTZ  = 0.09  // was EC50_hctz, unchanged value
GAMMA_HCTZ = 1.0   // no exponent in the original -> 1
EFFECT_HCTZ = EMAX_HCTZ * pow(C_HCTZ, GAMMA_HCTZ)
             / (pow(EC50_HCTZ, GAMMA_HCTZ) + pow(C_HCTZ, GAMMA_HCTZ));
```

`EFFECT_HCTZ` then replaces `hz_` verbatim at both of its two original call
sites (`Ksec` and `UCa_h`); `kK_hctz` and `kUCa_hctz` (how much the disease
model weights the compound's own effect) are untouched, since they are
disease-side scaling coefficients, not part of HCTZ's own PK/effect
interface.

No fit was needed or performed — `EMAX_HCTZ=1`/`GAMMA_HCTZ=1` are the
values the original's arithmetic already had, not new numbers.

## A necessary technical fix specific to this file: no `double` on `$ODE`-scoped capture names

`C_HCTZ` and `EFFECT_HCTZ` are computed inside `$ODE` (they must be — the
disease equations that consume them, `Ksec` and `UCa_h`, are themselves
computed later in the same `$ODE` block, at every solver step, not once per
interval) and are also listed in `$CAPTURE` so they are discoverable per
the guide's qspserver-compatibility rule (see below). Under mrgsolve 2.0.1,
this collides: **any bare `double NAME = expr;` inside `$ODE` is
auto-promoted by mrgsolve to a class member, and `$CAPTURE` naming that
same identifier requires mrgsolve to *also* declare a member of that exact
name** — so writing `double C_HCTZ = ...;` inside `$ODE` while also listing
`C_HCTZ` in `$CAPTURE` is a duplicate declaration:

```
158:10: error: redefinition of 'double {anonymous}::C_HCTZ'
  158 |   double C_HCTZ;
      |          ^~~~~~
50:10: note: 'double {anonymous}::C_HCTZ' previously declared here
```

This is the same underlying mrgsolve-2.0.1 auto-promotion mechanism already
documented for a different file in `translations/UPSTREAM_ISSUES.md` #35
(`dengue`, `$ODE` locals colliding with `$TABLE` captures) — but here it is
**introduced by this refactor's own new `$CAPTURE` entries**, not a
pre-existing defect in the untouched original (the original never captured
`C_hctz`/`hz_`), so it is fixed directly in the delivered
`drta_mrgsolve_model_refactored.R` itself (not logged to
`UPSTREAM_ISSUES.md`, and not a scratch-only workaround): inside `$ODE`,
`C_HCTZ`/`EFFECT_HCTZ` are assigned with **no `double` prefix** (mrgsolve
already provides the declaration because they're in `$CAPTURE`), while the
separate `$TABLE`-scoped `double C_HCTZ = o_C_HCTZ;` /
`double EFFECT_HCTZ = o_EFFECT_HCTZ;` (which is how this file's own
existing convention re-exposes every other `$ODE`-computed quantity for
capture, via a `o_*`-prefixed global bridge) keeps its explicit `double`,
since `$TABLE` does not hit the same auto-promotion collision. Confirmed
empirically both ways via the qspserver `mrgsolve_api` `/model_manifest`
endpoint (fails with `double` in `$ODE`, builds cleanly without it).

## qspserver `/model_manifest` discoverability

Per the guide's qspserver-compatibility requirements, checked against this
file, following the precedent set for `thyroid-eye-disease` (see
`ted_refactor_notes.md`, "`/model_manifest` discoverability"):

- `KA_HCTZ`, `CL_HCTZ`, `V1_HCTZ`, `EC50_HCTZ`, `EMAX_HCTZ`, `GAMMA_HCTZ`
  are all literal, fixed-value `$PARAM` entries, so `/model_manifest`
  already lists them (confirmed: `KA_HCTZ=1.1, CL_HCTZ=0.3168,
  V1_HCTZ=3.6, EC50_HCTZ=0.09, EMAX_HCTZ=1, GAMMA_HCTZ=1` all appear in the
  manifest response).
- `C_HCTZ` and `EFFECT_HCTZ` themselves are state-dependent (computed from
  `CENT_HCTZ`), so they cannot be `$PARAM` defaults. Both are added to
  `$CAPTURE` (additive-only — the original's own 24-item `$CAPTURE` list is
  untouched, these two are appended) so they are discoverable as
  **outputs** via `/model_manifest`'s `outputPaths` and retrievable via
  `/run_simulation`'s `outputs` list — confirmed both ways.
- No separate `<abbr>_mrgsolve_model_refactored.cpp` extraction was
  produced. Like `thyroid-eye-disease/ted_mrgsolve_model.R` (see
  `ted_refactor_notes.md`) and unlike `rheumatoid-arthritis`/others,
  `drta_mrgsolve_model.R` is **not** the `code <- '...'; mcode(code)`
  R-string-wrapper pattern the `.cpp`-extraction requirement targets — its
  own usage docstring (`mod <- mread("drta_mrgsolve_model.R")`) confirms it
  is raw mrgsolve block-marker DSL text meant to be read directly. The
  whole DSL portion of the file (`$PROB` through `$CAPTURE`) was posted
  as-is to `/model_manifest` and `/run_simulation` and compiled/ran
  correctly (see Verification below); the trailing `$ENV` block (R-only
  scenario-driver code — `SUBJ`, `LESION`, `GENOTYPE`, `make_regimen()`,
  `sim_scenario()`, the 28 `SCENARIOS`, `run_all_scenarios()`,
  `run_B21CS()`) was excluded from the posted `model_content` since it is
  not part of the mrgsolve DSL itself (mirrors the same exclusion in the
  TED case).

## A pre-existing upstream defect found while verifying (not fixed here)

Logged as `translations/UPSTREAM_ISSUES.md` #38. Unrelated to HCTZ:
`$GLOBAL`'s first line, `#define _MRG_DRTA_`, is a valueless macro
(unreferenced anywhere else in the file) that crashes mrgsolve 2.0.1's
parameter-definition parser (`Error in FUN(X[[i]], ...) : subscript out of
bounds`, in `pp_defs -> s_pick -> nonull -> unlist -> sapply -> lapply`)
before compilation is even attempted — **the untouched original does not
build at all** under the currently installed mrgsolve, independent of any
HCTZ content. Bisected down to this single line (confirmed: removing it
alone fixes the crash; restoring it alone reproduces it, with the rest of
`$GLOBAL`/`$PARAM`/`$CMT` held constant). Worked around, for verification
only, by deleting that one line from scratch copies of both
`drta_mrgsolve_model.R` and `drta_mrgsolve_model_refactored.R` (plus
excluding the trailing `$ENV` block per above) — neither the tracked
original nor the delivered `_refactored.R` was edited; both still contain
the bare `#define` exactly as written.

**Also noted (an engine-configuration mismatch, not a file defect):** this
model's own R driver runs every one of its 28 scenarios with
`mrgsim(..., maxsteps = 5e6, hmax = 0.5)` — 250x mrgsolve's default
`maxsteps` of 20000, plus an explicit step-size cap that the qspserver
`mrgsolve_api` `/run_simulation` endpoint has no field to set. Confirmed on
the untouched (`#define`-patched) original alone: an **undosed** run
exceeds the default 20000-step budget somewhere between 24h and 48h of
simulated time (`[lsoda] 20000 steps taken before reaching tout`), and a
run carrying **any** dosing event (tested individually on the bicarbonate,
citrate, potassium, KCl, and HCTZ compartments — every one fails alone)
can exceed it within hours. Bisection (holding dose amounts/compartments
fixed, varying only the requested end time) found dosed runs reliably
succeed through this API up to a 12h window and reliably fail by 24h, with
an unstable, non-monotonic boundary in between (18h failed, 20h passed in
one comparison) — consistent with the model being right at this API's
fixed step-budget edge, exactly the situation its own author's
`maxsteps=5e6` anticipated. This is why the verification below uses a 12h
single-dose window rather than reproducing the full 150-365 day duration
any of the model's own named scenarios (`S01`-`S28`) actually run for.

## Verification

Per the guide's mandatory protocol, via the qspserver `mrgsolve_api`
service (`POST /model_manifest`, `POST /run_simulation`,
`http://localhost:8007`), on the `#define`-patched (verification-only)
scratch copies described above.

**`/model_manifest`**: both the patched original and the refactored DSL
compiled successfully and reliably (3/3 repeat calls each) after the two
workarounds above; the refactored manifest additionally lists
`KA_HCTZ`, `CL_HCTZ`, `V1_HCTZ`, `EC50_HCTZ`, `EMAX_HCTZ`, `GAMMA_HCTZ` in
`parameters` and `C_HCTZ`, `EFFECT_HCTZ` in `outputPaths`, as required.

**`/run_simulation`**: reproduced the dosing composition of the original
file's own `S24_thiazide_addon` scenario — the *only* one of the model's 28
named scenarios that doses HCTZ at all (`reg = list(kind = "ADV",
mEq_kg_day = 1.0, times = c(8, 20), hctz_mg = 25, KCl_mmol_day = 20)`,
`subj = "child"`, `lesion = "complete"`) — as a single simultaneous
administration at `t=0` (rather than S24's own twice-daily-for-240-days
repeat pattern, which the step-budget limitation above rules out through
this API): `AG_citPR`/`AG_bicPR` 1.75 mmol / 9.75 mEq (the ADV7103 fast/slow
split of `1.0 mEq/kg * 30 kg` at `f_cit_ADV=0.35`), `AG_K` 15 mEq,
`AG_KCl` 10 mmol, and `AG_hctz`/`GUT_HCTZ` 25 mg — into the same 1-based
compartment indices in both files (dosing is index-addressed, and HCTZ's
compartment position in `$CMT` is unchanged by the rename). Subject
parameters set via the request's `parameters` field: `GFR0=115` (child),
`LES=0.18`, `LES_grad=0.58` (complete lesion) — `BW`, `AGE`, `BSA`,
`NINTAKE`, `f_Kcat` already matched the model's own defaults for this
subject/regimen and needed no override. `time = {end: 12, delta: 0.12}`
(100 points), comfortably inside the feasible window found by bisection.

Every one of the original's 24 pre-existing `$CAPTURE` outputs was
compared, plus the raw states `HCO3_e`, `K_pl`, `Ca_pl`, `UCa_s`, `Cl_pl`,
`ALDO`, `PTH`, and the HCTZ PK states themselves (`AG_hctz` vs `GUT_HCTZ`,
`C_hctz` vs `C_HCTZ`) — 33 outputs total, across all 101 time points.

**Result: near-exact match.** `AG_hctz` vs `GUT_HCTZ`: max abs diff = 0.0
(identical, as expected — this compartment's equation is an unchanged
rename). Every other output's max absolute deviation was at or below the
qspserver API's own JSON response rounding (outputs are serialized to 4
decimal places), the largest being `Hdef_day`/`NAE_day`/`acid_gap` at
0.0181 (relative deviation ~0.3-0.9%, on quantities the HCTZ kaliuresis
term feeds into only indirectly via the shared potassium-secretion
equation) and `C_hctz` vs `C_HCTZ` at exactly 1.0e-4 — a single tick of
that rounding grid. `reserve`, `eGFR`, `responder`, `HCO3_target`,
`NEAP_day`, `UCit_day`, and `Ca_pl` — outputs with no path to HCTZ at
all — matched at exactly 0.0, as expected. This is consistent with the
guide's expectation for a pure structural reorganization ("same math,
reorganized... anything beyond floating-point-scale deviation means a
bug") — no Hill-fit was performed or needed, since `EFFECT_HCTZ` is an
exact algebraic rewrite of the original's `hz_` ratio, and the residuals
here are bounded entirely by the API's own display-precision rounding, not
by any structural difference between the two files.

`EFFECT_HCTZ` itself ranged 0 to 0.9841 over the dosed window (peak
`C_HCTZ` = 5.5749 mg/L against `EC50_HCTZ` = 0.09 mg/L gives
5.5749/(5.5749+0.09) = 0.9841 exactly, confirming the Hill arithmetic is
correct, not just self-consistent).

## Anything else flagged

- `AG_hctz`/`GUT_HCTZ`'s position in `$CMT` (compartment 8 of 56) is
  unchanged, so the `$ENV`-block R driver's `CMT["AG_hctz"]` lookup (used
  by `make_regimen()`'s `hctz_mg` handling) was updated to
  `CMT["GUT_HCTZ"]` to keep pointing at the same compartment under its new
  name — this is a rename of HCTZ's own dosing-scaffolding reference, not a
  behavioral change (same index, same dose amounts/timing).
- No other compound's compartments, parameters, or equations were
  touched. The entire alkali-therapy PK (bicarbonate/citrate/potassium
  absorption and release, `NaDC1`-mediated citrate reabsorption), the
  saturating V-ATPase actuator, the three-sink acid disposal, the
  warm-started urine-pH bisection, bone/growth/calcium/PTH dynamics, and
  vitamin-D PK are all byte-identical between the two files.

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`,
`distal-renal-tubular-acidosis | HCTZ` row.
