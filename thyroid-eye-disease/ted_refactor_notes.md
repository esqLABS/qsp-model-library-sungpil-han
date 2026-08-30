# Refactor notes — `thyroid-eye-disease/ted_mrgsolve_model_refactored.R`

Scope: tocilizumab (TCZ) PK block and its downstream effect equation only, per
`FORK_WORKFLOW_GUIDE.md` Part 2. Every other compound in this file
(teprotumumab, glucocorticoid, rituximab, selenium) is copied verbatim,
unchanged, including all 9 prebuilt R-side dosing scenarios and the
simulation-running helper code at the bottom.

## Archetype determined

**Archetype 2 — no depot, two compartments, linear elimination.**

The original models tocilizumab with two compartments (`CEN_TCZ`, `PER_TCZ`)
and no absorption depot:

```
dxdt_CEN_TCZ = -(CL_TCZ/V1_TCZ)*CEN_TCZ - (Q_TCZ/V1_TCZ)*CEN_TCZ + (Q_TCZ/V2_TCZ)*PER_TCZ;
dxdt_PER_TCZ =  (Q_TCZ/V1_TCZ)*CEN_TCZ - (Q_TCZ/V2_TCZ)*PER_TCZ;
```

This is already CL/Q/V-parameterized (not micro-constants), so the rewrite is
a pure rename/reorganization into the guide's convention, combining the two
`-(CL/V1)` and `-(Q/V1)` terms on `CENT_TCZ` into a single `-(CL+Q)/V1` term
(algebraically identical, matches the archetype's worked example):

```
dxdt_CENT_TCZ = -(CL_TCZ + Q_TCZ)/V1_TCZ * CENT_TCZ + (Q_TCZ/V2_TCZ)*PERI_TCZ;
dxdt_PERI_TCZ =  (Q_TCZ/V1_TCZ)*CENT_TCZ - (Q_TCZ/V2_TCZ)*PERI_TCZ;
```

No TMDD, no depot, no infusion-as-ODE-input edge case — a clean fit to
Archetype 2, no compartments added or removed.

Renames: `CEN_TCZ`→`CENT_TCZ`, `PER_TCZ`→`PERI_TCZ`, local `CTCZ`→`C_TCZ`.
`CL_TCZ`, `V1_TCZ`, `V2_TCZ`, `Q_TCZ`, `EC50_TCZ` already matched the
convention in the original and are unchanged in value.

## Hill interface

The original's only use of tocilizumab's concentration downstream is a plain
occupancy ratio, `RO_IL6R = CTCZ/(EC50_TCZ+CTCZ)`, which is then scaled by an
**inline literal `0.3`** at its single point of use:

```r
double cyt_signal = CYT * (1.0 - 0.3*RO_IL6R);
```

This is already the Hill shape (`C/(C+EC50)`, i.e. gamma=1 implicit, Emax=0.3
implicit) — per the guide, this is a rename, not a refit. Two things were
extracted as explicit named parameters, values taken directly from the
original (no invented or refit values):

- `EMAX_TCZ = 0.30` — pulled out of the inline `0.3` literal in `cyt_signal`
- `GAMMA_TCZ = 1` — the original had no explicit Hill exponent

`EC50_TCZ` was already a named parameter (`1.0`, unchanged).

`RO_IL6R` (the bare occupancy ratio) is kept as-is because it is separately
captured (`IL6R_occupancy`) and used nowhere else — only its *scaled* use in
`cyt_signal` is now expressed as the named `EFFECT_TCZ`:

```
double EFFECT_TCZ = EMAX_TCZ * pow(C_TCZ, GAMMA_TCZ) / (pow(EC50_TCZ, GAMMA_TCZ) + pow(C_TCZ, GAMMA_TCZ));
...
double cyt_signal = CYT * (1.0 - EFFECT_TCZ);
```

Since `GAMMA_TCZ=1` and `EC50_TCZ`/`EMAX_TCZ` match the original exactly,
`EFFECT_TCZ ≡ 0.3*RO_IL6R` algebraically — an identity, not an approximation.
No steady-state fitting was needed (this is not a TMDD/ODE-derived interface).

## Verification

**Method:** ran the original file's own scenario 6 — "Tocilizumab 8 mg/kg IV
q4w x4 (steroid-resistant TED)", `ev(amt=8*70, cmt="CEN_TCZ"/"CENT_TCZ",
ii=28*24, addl=3)` — through both the original and refactored models via
`mread()` + `mrgsim()` (mrgsolve 2.0.1, `end=17520` h, `delta=24` h, 732 time
points), and compared all 23 `$CAPTURE`d outputs pointwise.

**Result: pass.** All 23 captured outputs matched to floating-point/solver
tolerance:

| Output | max relative dev. | max absolute dev. |
|---|---|---|
| `conc_tocilizumab` | 1.033e-04 | 6.34e-11 |
| `IL6R_occupancy` | 1.033e-04 | 2.32e-12 |
| `Fibroblast_idx` | 3.0e-10 | 4.7e-11 |
| `Hyaluronan_idx` | 1.8e-08 | 2.5e-09 |
| `EOMVolume_mL` | 1.2e-09 | 5.4e-08 |
| `Proptosis_mm` | 3.6e-10 | 5.4e-09 |
| all other 17 outputs | ≤ 3.3e-11 (several exactly 0) | ≤ 1.4e-10 |

The nominally largest *relative* deviation (1.033e-04, on `conc_tocilizumab`
and `IL6R_occupancy`) is an artefact of dividing by a near-zero value at a
trough time point, not a structural mismatch — the matching *absolute*
deviation is 6.34e-11 against a peak concentration of 176.7 mg/L, i.e.
~15 significant figures of agreement. This is consistent with ODE-solver
step-size noise from reordering algebraically-identical floating-point
operations (`-(CL/V1)-(Q/V1)` combined into `-(CL+Q)/V1`), not a bug, and
matches the guide's expectation for Archetypes 1–3 ("near-exact match... same
math, reorganized"). No Hill-fit was performed (none was needed), so there is
no fit-quality tradeoff to report — `EFFECT_TCZ` is an exact algebraic
rewrite of the original's inline `0.3*RO_IL6R`.

**Cross-checked via a second engine:** the same scenario was re-run through
the local qspserver `mrgsolve_api` service (`POST /model_manifest` then
`POST /run_simulation` at `http://localhost:8007`, dosing `amt=560` into
compartment index 7 — `CEN_TCZ`/`CENT_TCZ` — `ii=672`, `addl=3`, `end=17520`,
`delta=24`), as an independent confirmation alongside the local
`mread()`/`mrgsim()` run above. Result agreed with the local run: identical
pattern of exact matches on every non-TCZ output, and the same
floating-point-scale deviation on the TCZ-derived outputs (max relative
1.426e-04 on `conc_tocilizumab`/`IL6R_occupancy`, but max **absolute**
deviation only 1.0e-11) — consistent with the local result, not a
new/different discrepancy. `CAS_score` showed a max relative deviation of
1.0e-11 (absolute 1.0e-20) through this engine, not visible in the local run
at printed precision — still floating-point noise, several orders of
magnitude below anything that would indicate a structural difference.

## Two pre-existing upstream defects found while verifying (not fixed here)

Neither is related to tocilizumab or to this refactor; both were present in
the untouched original and were only worked around **in scratch, in-memory
copies used for verification** so the models would build at all — neither the
tracked original nor the delivered `_refactored.R` was edited to work around
them, and both files as committed still contain the original code exactly as
written. Logged in `translations/UPSTREAM_ISSUES.md`.

1. **`[PARAM]` annotation format error.** Two annotated parameter lines,
   `HILL_HEAR : 1.8` and `HILL_CUSH : 2.0`, have no third (description)
   field. mrgsolve's annotated-block parser rejects this outright
   (`Error: improper annotation format`) — the model does not build at all,
   for either file, until a description is added.
2. **`TIME` is not usable inside `$ODE` under mrgsolve 2.0.1.** The original
   uses `TIME` twice inside `[ODE]` (`t_month = TIME/730.0 + DURATION_MO;`
   and inside `smoke_mult_now`'s `exp(-TIME/4380.0)`). In this mrgsolve
   version, `TIME` expands to `self.time`, but `MRGSOLVE_ODE_SIGNATURE` (the
   generated `$ODE` function's parameter list) does not include a `self`
   parameter — only `$MAIN`/`$TABLE`/`$EVENT`/`$PREAMBLE` do. The correct
   in-`$ODE` accessor in this version is `SOLVERTIME` (`_ODETIME_[0]`).
   Confirmed via the generated `.cpp`: `error: 'self' was not declared in
   this scope`, pointing exactly at the `TIME` macro expansion on the
   `t_month` line. Reproduces identically from the refactored file (copied
   verbatim) and is unrelated to any TCZ change.

Verification of the TCZ refactor itself was carried out after mechanically
patching **both** scratch copies identically with (a) an added description
string on the two malformed `HILL_*` lines and (b) `TIME`→`SOLVERTIME` at the
two ODE-block usages — neither change alters any parameter value or any
PK/PD equation, so the comparison above is unaffected by these workarounds.

## qspserver `/model_manifest` discoverability

`FORK_WORKFLOW_GUIDE.md` picked up a new "qspserver compatibility
requirements" section mid-session (added by a sibling calibration run on a
different file) asking that a compound's Hill parameters and named
covariates be discoverable via the qspserver API. Checked against this file:

- `EMAX_TCZ`, `EC50_TCZ`, and `GAMMA_TCZ` are declared in `[PARAM]`, so
  `/model_manifest` already lists them (confirmed: `EMAX_TCZ=0.3`,
  `EC50_TCZ=1`, `GAMMA_TCZ=1` all appear in the manifest response).
- `EFFECT_TCZ` itself is computed in `[ODE]`, not `[PARAM]` (it depends on
  the compartment state, so it cannot be a fixed parameter). Added a
  `capture TCZ_effect = EFFECT_TCZ;` line to `[TABLE]`/`[CAPTURE]` so the
  quantity is at least discoverable as an **output** via `/model_manifest`
  and retrievable via `/run_simulation`'s `outputs` list — confirmed both
  ways after the change. This is additive only (a new capture unique to the
  refactored file, since the original never named this quantity) and does
  not appear in the original-vs-refactored diff above, which compares only
  the original's own pre-existing 23 captures.
- The guide section's requirement to extract a bare-DSL `.cpp` sibling
  (because `model_content` "must be pure mrgsolve DSL, no R wrapper") does
  **not** apply to this file: `ted_mrgsolve_model.R` is already written
  directly in mrgsolve block-marker syntax (`[PROB]`/`[PARAM]`/.../`[CAPTURE]`
  with `##`-prefixed prose before/after, consumed via `mread()`), not the
  `code <- '...'; mcode(code)` R-string-wrapper pattern that requirement is
  aimed at. Confirmed empirically: the whole file content posted as-is to
  `/model_manifest` and `/run_simulation` compiled and ran correctly (see
  Verification above) — no separate `.cpp` extraction was necessary or
  produced.

## Anything else flagged

- No other compound's compartments, parameters, or effect equations were
  touched. `conc_teprotumumab`, `GC_signal`, `BCell_idx`, `TRAb_idx`, and all
  safety-tracker outputs matched **exactly** (0.000e+00 deviation), as
  expected since that code is byte-identical between the two files.
- The R-side scenario-helper comment block at the bottom of the file
  (`TED_simulate_scenarios()`) doses tocilizumab via `ev(cmt="CEN_TCZ", ...)`
  in scenario 6. Since this refactor renames that compartment to
  `CENT_TCZ`, a one-line note was added next to that comment block warning
  that `ev(cmt=...)` calls need updating to `CENT_TCZ` if dosing against the
  refactored file directly — the comment block itself is otherwise copied
  verbatim (per scope, it is R-side documentation, not tocilizumab PK code).
