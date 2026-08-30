# Cluster Headache (CH) — PK/PD refactor notes

Refactored `ch_mrgsolve_model.R` -> `ch_mrgsolve_model_refactored.R`, applying
the fork's pluggable-PK naming convention to all seven of the file's drugs:
**Galcanezumab (GALCA)**, **Lithium (LI)**, **Prednisolone (PRED)**,
**Sumatriptan (SUMA)**, **Topiramate (TOPI)**, **Verapamil PR (VERA)**,
**Zolmitriptan IN (ZOL)**. The original is untouched. All numeric parameter
values are copied verbatim; only PK structure/naming and the Hill-interface
expression are reorganized.

## Archetype per compound

All seven were already clean, single-site "Redirect concentration" cases per
the census -- no TMDD, no shared/combined Hill terms, no bespoke structure
needed for any of them.

| Compound | Original compartments | Archetype | Notes |
|---|---|---|---|
| SUMA (sumatriptan SC) | `SUMA_DEPOT`, `SUMA_CENT` | 3 minus peripheral (depot+central) | plain first-order absorption/elimination |
| ZOL (zolmitriptan IN) | `ZOL_DEPOT`, `ZOL_CENT`, `ZOL_PER` | 3 (depot+central+peripheral) | see "Zolmitriptan IN absorption" below -- checked, not biphasic |
| VERA (verapamil PR) | `VERA_DEPOT`, `VERA_CENT`, `VERA_PER` | 3 (depot+central+peripheral) | CL/Q/V convention already used in the original |
| LI (lithium) | `LI_DEPOT`, `LI_CENT` | 3 minus peripheral (depot+central) | renal covariate scaling preserved, see below |
| TOPI (topiramate) | `TOPI_DEPOT`, `TOPI_CENT` | 3 minus peripheral (depot+central) | pre-existing build defect in this block, see below |
| GALCA (galcanezumab) | `GALCA_SC`, `GALCA_CENT` | 3 minus peripheral (depot+central) | see "Galcanezumab: checked for TMDD" below -- linear, not TMDD |
| PRED (prednisolone) | `PRED_DEPOT`, `PRED_CENT` | 3 minus peripheral (depot+central) | dosed as prednisone prodrug, active moiety named PRED |

Renamed compartments: `GUT_<STEM>`/`CENT_<STEM>`/`PERI_<STEM>` per the guide.
The nine disease-side/device compartments (`HYPO_DRIVE`, `CGRP`, `PACAP`,
`PIAL`, `ATTACK_HZ`, `CUM_ATTACKS`, `BOUT_TIMER`, `O2_EFFECT`, `GON_EFF`) are
untouched -- they are not one of the seven compounds' own PK. Compartment
*order* (and therefore 1-based `$CMT` index) is identical between the
original and the refactored file, which is what let dosing events in
verification target the same numeric `cmt` index in both.

## Galcanezumab: checked for TMDD

The original's own header comment already says "galcanezumab PK (SC mAb,
FcRn recycled, linear)", and the `$CMT`/`$ODE` blocks confirm this: there is
no free-receptor or drug-receptor-complex compartment anywhere in the file
for galcanezumab, just `GALCA_SC` (depot) and `GALCA_CENT` (central) with
plain linear clearance (`dxdt_GALCA_CENT = KA_GALCA*GALCA_SC -
(CL_GALCA/V_GALCA)*GALCA_CENT`). Its effect on hypothalamic drive is a plain
`EMAX_GALCA * CP_GALCA / (CP_GALCA + EC50_GALCA)` ratio, not a fractional-
occupancy term over any receptor-complex state. So this is Archetype 3 minus
peripheral (depot+central), a rename, **not** TMDD and not a curve fit — the
Hill-interface rename applies directly to `C_GALCA`, exactly as it does for
every other non-TMDD compound below.

## Zolmitriptan IN absorption

Checked per the task's flag that intranasal absorption can be biphasic: the
original's `$ODE` models zolmitriptan absorption as a single first-order
depot, `dxdt_ZOL_DEPOT = -KA_ZOL * ZOL_DEPOT`, feeding a two-compartment
central+peripheral system — no dual-rate or delayed-release term anywhere.
So despite the IN route, this file's own PK is ordinary Archetype 3 (full
depot+central+peripheral), a rename with no special handling required.

## Lithium renal covariate scaling

The original scales lithium's clearance by creatinine clearance in `$MAIN`:
`double cl_li = CL_LI * (CrCL/100);`, then uses `cl_li` (not `CL_LI`
directly) inside `$ODE`'s lithium elimination term. This `$MAIN`-computed,
`$ODE`-consumed local variable pattern is preserved exactly, renamed
`CL_LI_eff` for clarity; `CL_LI` itself keeps the plain naming-convention
name and its original value.

## Hill interface: renames, not fits

Every one of the seven original effect terms was already written as a plain
`EMAX * C / (C + EC50)` ratio (implicit Hill coefficient of 1, never stated
explicitly in the original) — `E_vera`, `E_li`, `E_topi`, `E_galca`,
`E_pred`, `E_suma_acute`, `E_zol_acute`. Every one of `EFFECT_SUMA`,
`EFFECT_ZOL`, `EFFECT_VERA`, `EFFECT_LI`, `EFFECT_TOPI`, `EFFECT_GALCA`,
`EFFECT_PRED` is the identical expression with `GAMMA_<STEM> = 1` added
explicitly (matching the original's implicit exponent) and no other change.
Confirmed a pure rename by the exact-match verification below (not merely
"close"). No `nls()` fitting was needed or performed for any compound.

The disease-side combination logic is untouched: `preventive` (Bliss-style
combination of `EFFECT_VERA`/`EFFECT_LI`/`EFFECT_TOPI`/`EFFECT_GALCA`/
`EFFECT_PRED`/GON-block) feeds `HYPO_DRIVE`'s dampening exactly as the
original's `preventive` did; `acute_abort` (combination of `EFFECT_SUMA`/
`EFFECT_ZOL`/O2) feeds the attack hazard exactly as the original's
`acute_abort` did.

## `$PARAM` vs `$GLOBAL` for `C_<STEM>`/`EFFECT_<STEM>`

Per the guide's qspserver compatibility requirement #2, these should ideally
live in `$PARAM` (with a `= 0` default) for direct `/model_manifest`
discoverability. This was not done here, for the same reason already
documented in the AAA/AMD/breast-cancer refactors: mrgsolve 2.0.1 compiles
`$PARAM` members as **read-only references** inside `$ODE`, so a value that
must be recomputed every timestep from state (all fourteen of `C_SUMA`,
`C_ZOL`, `C_VERA`, `C_LI`, `C_TOPI`, `C_GALCA`, `C_PRED`, `EFFECT_SUMA`,
`EFFECT_ZOL`, `EFFECT_VERA`, `EFFECT_LI`, `EFFECT_TOPI`, `EFFECT_GALCA`,
`EFFECT_PRED`) cannot also be declared in `$PARAM`. Instead, all fourteen
are predeclared as `double`s in `$GLOBAL` and listed in `$CAPTURE`: visible
in every simulation's output columns and in `/model_manifest`'s
`outputPaths`, just not in its `parameters` list. The original's own
`$TABLE`-style `capture ConcSuma = C_SUMA;` etc. aliases are kept alongside,
unchanged in name and formula, for continuity with the original's own output
columns.

## Pre-existing upstream build defect (fixed syntax-only, inside topiramate's own block)

**The original does not compile under mrgsolve 2.0.1 as written.** Confirmed
via `POST /model_manifest` on the untouched original's own DSL (extracted
verbatim from `ch_code <- '...'`):

```
Error: improper annotation format
 input: EMAX_TOPI :  0.35
 context: parse annotated parameter block (PARAM)
```

The `$PARAM @annotated` line `EMAX_TOPI :  0.35` supplies only two of the
three required colon-delimited fields (`name : value`, missing
`: description`), so mrgsolve's annotation parser rejects the whole `$PARAM`
block at this line. This is the only such line in the file (checked by
scanning every other `$PARAM @annotated` line for the same one-colon-short
pattern — none found). Logged as
[`UPSTREAM_ISSUES.md` #95](../translations/UPSTREAM_ISSUES.md).

Because this defect sits inside topiramate's own PK/PD block — one of this
refactor's seven target compounds, not an incidental/unrelated line — the
fix is folded directly into topiramate's own naming-convention rename in
`ch_mrgsolve_model_refactored.R`, per the guide's point 4 ("if a defect is
genuinely inside the scope compound's own block, fix it as part of the
refactor itself"), rather than treated as a separate disclosed build-compat
patch: a purely descriptive third field was appended (`EMAX_TOPI :  0.35  :
Max fractional attack-rate reduction (preventive); ...`), changing nothing
about the parameter's name, value, or use in any equation. The checked-in
original `ch_mrgsolve_model.R` is left untouched and still carries the
defect forward, per the never-edit-upstream rule.

## Verification

**Method.** Both the original's own DSL (with only the one-line `EMAX_TOPI`
description-field fix applied to a throwaway verification copy, so it would
compile at all — never applied to the checked-in original) and
`ch_mrgsolve_model_refactored.R`'s embedded DSL were extracted as bare
mrgsolve DSL text (accounting for the R-string `\'` escape used once in the
`$GLOBAL` comment) and run through the qspserver `mrgsolve_api` container
(`POST /model_manifest`, `POST /run_simulation`) at `http://localhost:8007`,
requests spaced ~2s apart. `POST /model_manifest` confirmed the refactored
DSL compiles cleanly and exposes all 14 `C_<STEM>`/`EFFECT_<STEM>` terms plus
every original alias (`ConcSuma`, `HypoDrive`, `Preventive`, `AcuteAbort`,
etc.) in `outputPaths`.

All seven of the original's own dosing scenarios (`S0`-`S6`, defined in the
original's own `build_scenario()`/`simulate_all()`) were run, identical
dosing (same compartment index, amount, `ii`/`addl`, timing), over the same
12-week/`delta=1` window `simulate_all()` itself defaults to (2016 h, ~2017-
2129 output points depending on scenario). Every output shared between the
two files' `$CAPTURE` blocks was compared point-by-point: `HYPO_DRIVE`,
`CGRP`, `PACAP`, `PIAL`, `ATTACK_HZ`, `CUM_ATTACKS`, `BOUT_TIMER`,
`O2_EFFECT`, `GON_EFF` (compartments, identical names in both files) plus
`ConcSuma`, `ConcZol`, `ConcVera`, `ConcLi`, `ConcTopi`, `ConcGalca`,
`ConcPred`, `HypoDrive`, `CGRPtone`, `PACAPtone`, `Pialtone`, `HazardPerH`,
`AttacksWeek`, `Preventive`, `AcuteAbort` (the original's own `$TABLE`
aliases, kept unchanged in the refactored file).

**Result: exact match on all seven scenarios** (max abs diff 0.0, max rel
diff 0.0) — pure structural reorganization, as expected for seven Archetype-3
compounds with rename-only Hill terms, no numeric drift whatsoever:

| Scenario | Description | Points | Max abs diff |
|---|---|---|---|
| S0 | No treatment | 2017 | 0.0 |
| S1 | O2 + sumatriptan abortive | 2019 | 0.0 |
| S2 | Verapamil 240->480 mg/d preventive | 2045 | 0.0 |
| S3 | Verapamil + lithium (chronic CH) | 2129 | 0.0 |
| S4 | Galcanezumab 300 mg SC monthly | 2020 | 0.0 |
| S5 | Prednisone bridge + verapamil | 2079 | 0.0 |
| S6 | GON block + verapamil | 2046 | 0.0 |

**A solver-flakiness note, disclosed rather than hidden:** on the first pass
at the finest (`delta=1`) grid, the *original* build-fix copy transiently
failed on S1 (an `lsoda`/`intdy` internal-step warning cascading into
"simulation terminated") and on S3 (a worker crash, "malloc(): invalid size
(unsorted)", HTTP 500 "killed by signal 6"), while the refactored DSL
succeeded on the same requests. Both were re-run immediately afterward with
byte-identical requests and both succeeded on retry, matching the refactored
DSL exactly (max abs diff 0.0) as shown in the table above. This reproduces
the API's documented behavior under load ("has crashed under heavy
concurrent load before") rather than any real difference between the two
DSLs — the two files are otherwise structurally identical aside from the
renamed PK compartments/parameters and the Hill-term rename, and the
retried, successful runs prove the numerics match exactly. No scenario
window needed to be shortened; the full 12-week/`delta=1` grid matching
`simulate_all()`'s own defaults was used throughout.

`EMAX_TOPI`'s description-field fix, and the `$GLOBAL`-instead-of-`$PARAM`
placement for `C_<STEM>`/`EFFECT_<STEM>`, are the only differences between
what the original *would* compile to and what the refactored file compiles
to; both are confirmed non-numeric by the exact-match result above.

## Census

Recorded against the seven existing rows in
`driver-patches/data/compound_perturbation_census.md` (see that file) rather
than starting a second parallel tracker.
