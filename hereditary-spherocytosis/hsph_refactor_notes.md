# Refactor notes — `hsph_mrgsolve_model_refactored.R`

Scope: **mitapivat (MIT) only** — its own PK compartments, its own dosing
logic, and its own two disease-facing effect terms (ATP increase, 2,3-DPG
decrease). Every other equation in the file (the 9-cohort geometric
haemolysis model, erythropoiesis, spleen, haem catabolism, iron,
transfusion, parvovirus) is untouched apart from the mechanical,
non-numeric build-compatibility fixes described below, which were needed
to get the file to compile at all and are disclosed in full.

## The assigned defect: duplicate concentration sites

Per `driver-patches/data/compound_perturbation_census.md`, mitapivat here
was classified **"Normalize duplicate concentration sites, then
redirect."** Reading the original confirmed exactly that: mitapivat's own
plasma concentration (`CENT/Vc*1000`) was computed **independently, with
two different local names, in two different mrgsolve blocks**:

- `$ODE`: `double conc = MCEN/Vc_m*1000.0;` — feeds `drv = conc/(EC50_m+conc)`,
  which drives both `dxdt_ATP` and `dxdt_DPG`.
- `$TABLE`: `double MITA = MCEN/Vc_m*1000.0;` — reporting-only, captured as
  the file's only mitapivat output.

Same formula, same numeric value at every timepoint, but two textually
separate definitions with no shared name — exactly the kind of thing that
silently drifts apart the next time someone edits one copy and not the
other. Normalized to **one name, `C_MIT`**, predeclared once in `$GLOBAL`
and assigned from the identical formula at both call sites (see the
qspserver-compatibility section below for why it has to be reassigned at
both sites rather than computed once and shared).

## Archetype determination

Mitapivat's own block, read directly from the `dxdt_` lines:

- 3 compartments: `MGUT` (oral depot), `MCEN` (central), `MPER`
  (peripheral).
- Linear absorption (`ka_m`) into a linear two-compartment distribution
  model, already parameterized as `CL_m`/`Vc_m`/`Q_m`/`Vp_m` (i.e.
  CL/Q/V-style, not micro-constants) — no conversion needed, just renaming.
- No receptor/complex compartment anywhere — mitapivat's mechanism (PKR
  activation) is modeled purely via a plasma-concentration Hill-style
  ratio (`conc/(EC50+conc)`), not TMDD.

**Conclusion: Archetype 3 (depot + central + peripheral, linear
elimination)**, with no bioavailability parameter in the original (no
`F_MIT` — the guide's archetype 3 makes that optional, and no `F` term was
invented here).

One structural quirk, preserved unchanged: the original's rate constants
are stated in **per-hour** (`ka_m`, `CL_m`, `Q_m` — see their own `$PARAM`
descriptions) while the model integrates in **days**, so every ODE term
carries an explicit `*24.0` conversion factor. This is dimensionally
required, not a stylistic artifact, and is kept exactly as the original
wrote it, just with renamed symbols.

**Dosing edge case, preserved unchanged:** mitapivat is dosed as a
continuous, algebraically-pulsed zero-order input directly inside `$ODE`
(`if (DOSE_MIT > 0.0) { double ph = fmod(SOLVERTIME, 0.5); if (ph < 0.02)
dose_mit = DOSE_MIT/0.02; }` — a ~30-minute-wide pulse every 0.5 days,
i.e. BID), not an `ev()`/`data_set()` event. Per the guide's named edge
case, this is kept exactly as-is; converting it to event-based dosing
would change how the model's own scenarios must be run.

## Renaming applied (values unchanged from the original)

| Original | Refactored | Role |
|---|---|---|
| `MGUT` | `GUT_MIT` | absorption depot compartment |
| `MCEN` | `CENT_MIT` | central compartment |
| `MPER` | `PERI_MIT` | peripheral compartment |
| `ka_m` | `KA_MIT` | absorption rate (1/h) |
| `CL_m` | `CL_MIT` | clearance (L/h) |
| `Vc_m` | `V1_MIT` | central volume (L) |
| `Q_m` | `Q_MIT` | intercompartmental clearance (L/h) |
| `Vp_m` | `V2_MIT` | peripheral volume (L) |
| `EC50_m` | `EC50_MIT` | concentration for half-maximal PKR activation (ng/mL) |
| `Emax_atp` | `EMAX_MIT_ATP` | maximal fractional ATP increase (-) |
| `Emax_dpg` | `EMAX_MIT_DPG` | maximal fractional 2,3-DPG decrease (-) |
| `dose_m` | `DOSE_MIT` | dose per administration, BID (mg) |
| `conc` ($ODE local) / `MITA`'s inline formula ($TABLE local) | `C_MIT` (one $GLOBAL-cached site) | the normalized, single-defined exposed concentration (ng/mL) |
| `drv` | `OCC_MIT` | fractional PKR-activation occupancy, `C/(EC50+C)` |

`k_atp`/`k_dpg` (ATP/2,3-DPG turnover rates) and the `ATP`/`DPG`
compartments themselves are disease-side red-cell metabolic physiology,
not part of mitapivat's own PK/PD block — left exactly as named and
valued, per the guide's instruction not to touch other compounds/equations
outside scope. (There is only one compound in this file, so "other
compounds" here means "the disease model," which was otherwise untouched.)

**New parameter, not in the original:** `GAMMA_MIT = 1`. The original's
`conc/(EC50_m+conc)` has no Hill exponent — this makes that explicit as a
named parameter per the guide's Hill-interface convention. Since
`pow(x,1) = x`, this is a rename, not a refit.

## Hill interface

The original's PKR-activation term is already a plain ratio,
`drv = conc/(EC50_m+conc)`, feeding two separate downstream effects with
two separate `Emax` scalars — so per the guide this is **a rename, not a
refit**:

```c
OCC_MIT = pow(C_MIT, GAMMA_MIT) / (pow(EC50_MIT, GAMMA_MIT) + pow(C_MIT, GAMMA_MIT));
EFFECT_MIT_ATP = EMAX_MIT_ATP * OCC_MIT;
EFFECT_MIT_DPG = EMAX_MIT_DPG * OCC_MIT;
```

`dxdt_ATP` now reads `+ EFFECT_MIT_ATP` in place of `+ Emax_atp*drv`;
`dxdt_DPG` now reads `- EFFECT_MIT_DPG` in place of `- Emax_dpg*drv`. Both
effects share the same underlying occupancy (`OCC_MIT`) exactly as the
original's single `drv` did — this is the same pattern as one compound
driving two named `EFFECT_<STEM>_*` terms in the polymyalgia-rheumatica
and cervical-cancer refactors (shared occupancy computed once, consumed by
more than one named, separately-scoped effect). No refit was performed:
`GAMMA_MIT = 1` makes `pow(x,1) = x`, so the formula shape is identical to
the original's.

## Why `C_MIT`/`OCC_MIT`/`EFFECT_MIT_ATP`/`EFFECT_MIT_DPG` are `$GLOBAL`-predeclared, not `$PARAM`-declared

The guide's qspserver-compatibility section calls for every `C_<STEM>`/
`EFFECT_<STEM>` to live in `$PARAM` (with a `=0` default) so
`/model_manifest` can list it as a discoverable covariate. That is not
possible here: mrgsolve 2.0.1 compiles `$PARAM` members as **read-only
references** inside `$ODE`/`$TABLE`, so a quantity that must be
*recomputed* every timestep from live state (which is exactly what
normalizing the duplicate-site defect requires) cannot also be declared in
`$PARAM` — the same constraint already documented in the breast-cancer,
AMD, membranous-nephropathy, cervical-cancer, and stable-angina refactors
in this corpus. Instead, `C_MIT`, `OCC_MIT`, `EFFECT_MIT_ATP`, and
`EFFECT_MIT_DPG` are predeclared as plain `double`s in `$GLOBAL`, assigned
by plain statement (no re-declaring `double`) in `$ODE`, **and
re-assigned again, identically, in `$TABLE`** (matching the original's own
`$TABLE`-recomputes-`MITA`-independently behavior, and the same pattern
used in the cervical-cancer/breast-cancer/AMD refactors to avoid reading a
one-timestep-stale `$GLOBAL` double at a row where `$ODE` and `$TABLE`
evaluation could otherwise disagree). All four are listed bare in
`$CAPTURE` and confirmed present in `/model_manifest`'s `outputPaths` (see
Verification below) — discoverable as **outputs**, not as `$PARAM`
covariates, which is the same trade-off already accepted and disclosed in
every other refactor that hit this same mrgsolve constraint.

This is precisely what "normalize the duplicate site" means in practice:
before, `conc` and `MITA`'s formula were two independently-maintained
copies of the same expression with two different names; now there is one
name (`C_MIT`) and one formula, assigned at exactly the two points
mrgsolve's block model requires a value to be (re)computed at, rather than
two different names computing two textually-independent copies of the
same thing.

## Build-compatibility fixes (unrelated to mitapivat, logged upstream)

**The untouched original does not compile at all** under mrgsolve 2.0.1 —
confirmed via the qspserver `mrgsolve_api` container's `/model_manifest`,
independent of any refactor changes. Four separate, pre-existing defect
classes, none inside mitapivat's own block, fully written up as
`translations/UPSTREAM_ISSUES.md` **entry #107**:

1. A stray, duplicate `@annotated` marker inside `$PARAM` (`Error: Found
   duplicated block option names.`) — deleted.
2. `N_1`..`N_9` used in `$MAIN` with no declaration anywhere in the file
   (`'N_1' was not declared in this scope`) — declared, one `double`
   statement per name, before first use.
3. Fourteen multi-variable `double a = x, b = y, ...;` declarations across
   `$MAIN`/`$ODE`/`$TABLE` that lose the `double` keyword from every name
   after the first during mrgsolve's build (same mechanism as issues #37
   and #104) — each split into one `double NAME = value;` statement per
   variable, values unchanged.
4. mrgsolve 2.0.1 hoists every `$MAIN`/`$ODE`/`$TABLE` block-local
   `double`/`int` into one shared anonymous-namespace scope, so eleven
   scratch names reused across those three blocks collide
   (`redefinition of 'double {anonymous}::...'`): the for-loop counter `i`
   (five separate loops, one shared bare name), `ph` (transfusion-phase
   vs. mitapivat-dosing-phase), and `n`/`ex`/`vk`/`rho`/`ce`/`ig`/`fo`/`sl`/
   `Aent` (each independently recomputed once in `$ODE` and again in
   `$MAIN`/`$TABLE`). Renamed for uniqueness: `$MAIN`'s loop counter ->
   `im`; `$ODE`'s three loop counters -> `i` (first, unchanged), `i2`,
   `i3`; `$TABLE`'s loop counter -> `it`; `$ODE`'s transfusion-phase `ph`
   -> `ph_tx` (mitapivat's own kept as `ph`); `$TABLE`'s per-cohort-loop
   `n`/`ex`/`vk`/`rho`/`ce`/`ig`/`fo`/`sl` -> `n_t`/`ex_t`/`vk_t`/`rho_t`/
   `ce_t`/`ig_t`/`fo_t`/`sl_t`; `$ODE`'s `Aent` -> `Aent_ode` (`$MAIN`
   keeps `Aent`). Every rename is purely a symbol change with identical
   per-iteration arithmetic.

All four fixes are applied **directly in the delivered
`hsph_mrgsolve_model_refactored.R`**, disclosed here in full, confirmed to
change nothing numeric by the exact-match verification below. The
checked-in `hsph_mrgsolve_model.R` is completely untouched and still
carries every one of these defects exactly as written.

## qspserver compatibility

- **`model_content` / no `.cpp` extraction needed.** `hsph_mrgsolve_model.R`
  is **not** the `code <- '...'; mcode(code)` R-string-wrapper pattern the
  guide's `.cpp`-extraction step targets — its own usage docstring (`mod
  <- mread("hsph_mrgsolve_model", ".")`) confirms it is raw mrgsolve
  block-marker DSL text meant to be read directly from disk, with no R
  code wrapped around it (confirmed: zero occurrences of `code <- '`
  anywhere in the file). Same situation as
  `distal-renal-tubular-acidosis/drta_mrgsolve_model.R` and
  `thyroid-eye-disease/ted_mrgsolve_model.R` (see their own
  `*_refactor_notes.md`). The whole file, verbatim, **is** the DSL text
  that was posted to `/model_manifest`/`/run_simulation` as
  `model_content` for both the original (with the in-memory-only
  build-compat patch applied for verification) and the refactored
  deliverable — no separate `hsph_mrgsolve_model_refactored.cpp` was
  produced, and none is needed.
- **Covariates discoverable.** `/model_manifest` on the refactored DSL
  confirms `KA_MIT`, `CL_MIT`, `V1_MIT`, `Q_MIT`, `V2_MIT`, `EC50_MIT`,
  `GAMMA_MIT`, `EMAX_MIT_ATP`, `EMAX_MIT_DPG`, and `DOSE_MIT` all appear
  under `parameters`; `C_MIT`, `OCC_MIT`, `EFFECT_MIT_ATP`,
  `EFFECT_MIT_DPG`, `GUT_MIT`, `CENT_MIT`, `PERI_MIT` all appear under
  `outputPaths` (121 parameters / 99 outputs total on the untouched,
  build-patched original vs. the refactored file's own manifest — see raw
  responses referenced below).
- **`$SET` unaffected either way** — this file has no `$SET end/delta`
  block at all; the time grid is supplied entirely by the caller (`mrgsim
  end=/delta=` in local usage, or the request's own `time` field via the
  API), so there is nothing here that could be load-bearing on a
  request-supplied grid vs. a model-default one.
- **`$CAPTURE` complete.** `MITA` (legacy name, now literally `= C_MIT`),
  `C_MIT`, `OCC_MIT`, `EFFECT_MIT_ATP`, and `EFFECT_MIT_DPG` are all
  explicitly listed in `$CAPTURE`; `GUT_MIT`/`CENT_MIT`/`PERI_MIT`/`ATP`/
  `DPG` are ordinary `$CMT` compartments and therefore already present in
  every simulation's output regardless of `$CAPTURE`.
- **No R-only syntax inside the DSL.** Nothing in mitapivat's block (or
  anywhere else touched) depends on an R-side helper, constant, or
  string-interpolated value — everything needed is declared inside the
  DSL block itself (`$PARAM`, `$GLOBAL`, `$CMT`).

## Verification

**Method:** qspserver `mrgsolve_api` container (`http://localhost:8007`,
confirmed healthy via `GET /health`), `POST /model_manifest` and
`POST /run_simulation` — no local R/mrgsolve used, per this fork's
verification protocol. Requests spaced ~2 s apart per the API's stated
concurrency limit.

- `/model_manifest` on the refactored DSL: **HTTP 200**, compiles cleanly,
  every renamed/added parameter and output confirmed present (see above).
- `/model_manifest` on the untouched original: **HTTP 500** (`Error:
  Found duplicated block option names.`) — does not compile at all, per
  the build-compatibility section above. A scratch, in-memory-only copy of
  the original with the four syntax-only fixes applied (never written to
  `hsph_mrgsolve_model.R`) was then used for the actual comparison run —
  **HTTP 200**, 121 parameters / 99 output paths.

**Scenarios run — the file's own usage.** This file has no bundled R
scenario-runner script; its own `Usage` docstring shows `mrgsim(end=1400,
delta=1)` as the standard run length, and the accompanying
`hsph_shiny_app.R`'s own preset scenarios use `dose_m = 100` (mg BID,
e.g. `` `4 미타피밧 100 mg` = c(dose_m = 100) ``) as mitapivat's own dosing
level — so both were used exactly as found, not invented:

| Scenario | Parameters | Result |
|---|---|---|
| No drug | `dose_m`/`DOSE_MIT = 0` | exact match, every output, 0 deviation |
| Mitapivat 100 mg BID | `dose_m`/`DOSE_MIT = 100` | exact match, every output, 0 deviation |

Both runs used `time = {end: 1400, delta: 1}` (1401 points, t = 0…1400
days) — the file's own stated horizon; no shortening was needed, both
completed well inside the API's default solver-step budget.

Every one of the original's 29 `$CAPTURE`d/state outputs (`Hb`, `Hct`,
`RBC`, `MCV`, `MCH`, `MCHC`, `RETpct`, `RETabs`, `LIFE`, `AREA`, `VOLc`,
`DCRIT`, `SPHER`, `EMA`, `MCF`, `IGGC`, `RCORD`, `FGEOM`, `FOPS`, `FSEN`,
`FLIV`, `FLYS`, `TBIL`, `BRPROD`, `SPLEEN`, `STONEPCT`, `MITA`, `ATP`,
`DPG`) matched **exactly** (max abs diff = 0.0, max rel diff = 0.0)
between the build-patched original and the refactored file, at every one
of the 1401 timepoints, in both scenarios. This is the expected result for
a pure structural reorganization with no Hill-fitting (Archetype 3, plain
rename), per the guide's tolerance rule.

Sanity-checked the four new outputs directly, at the mitapivat-100-mg-BID
scenario's final timepoint (t = 1400 d): `MITA == C_MIT` exactly at every
timepoint (1177.0658 ng/mL at steady dosing); `OCC_MIT =
1177.0658/(1150+1177.0658) = 0.5058`; `EFFECT_MIT_ATP = 0.62*0.5058 =
0.3136`; `EFFECT_MIT_DPG = 0.45*0.5058 = 0.2276` — all algebraically
consistent with the renamed parameters' own values, confirming the
rewired Hill interface is wired correctly, not just "matches the
original's numbers by coincidence."

**Result: verification passed, exact match, as required for a pure
Archetype-3 rename with no Hill-fitting.**

## Other findings flagged for manual review

- The file's own header states "THIS FILE HAS NOT BEEN RUN" (no R
  toolchain was available when it was authored; `hsph_python_reference.py`
  is the actual model of record). That claim is confirmed accurate by this
  refactor: the untouched original does not compile under mrgsolve 2.0.1
  at all, for reasons having nothing to do with mitapivat (see
  `UPSTREAM_ISSUES.md` #107). Worth a maintainer pass to reconcile the R
  and Python implementations properly, independent of this refactor's
  narrow scope.

## Discoverability fix

A corpus-wide discoverability audit found `C_MIT` was not written as a
single contiguous `double C_<STEM> = <expr>;` statement anywhere in the
file. It was `$GLOBAL`-predeclared on a shared line (`double C_MIT,
OCC_MIT, EFFECT_MIT_ATP, EFFECT_MIT_DPG;`), bare-assigned once in `$ODE`
(`C_MIT = CENT_MIT/V1_MIT*1000.0;`), and bare-reassigned again in `$TABLE`
with the identical formula — this file is in fact the one the guide's own
worked discoverability example is based on (see the extensive in-file
comments at the `$GLOBAL` declare and at the `$TABLE` reassignment
explaining why the value is deliberately recomputed in both blocks: it is
the corpus's reference case for the "recompute in `$TABLE` to avoid a
one-step-stale read" fix, just not yet in the literal-text-discoverable
shape). `OCC_MIT`, `EFFECT_MIT_ATP`, and `EFFECT_MIT_DPG` on the same
`$GLOBAL` line are untouched — out of scope for this fix.

**Fix applied:**
1. Removed `C_MIT` from the shared `$GLOBAL` declare line (`OCC_MIT`,
   `EFFECT_MIT_ATP`, `EFFECT_MIT_DPG` left in place).
2. The file already had a bare `$TABLE`-side reassignment (`C_MIT =
   CENT_MIT/V1_MIT*1000.0;`, immediately before the `OCC_MIT`/
   `EFFECT_MIT_ATP`/`EFFECT_MIT_DPG`/`MITA` lines, ahead of `$CAPTURE`);
   per the guide's instruction for this case, that line was converted in
   place to `double C_MIT = CENT_MIT/V1_MIT*1000.0;` rather than adding a
   second, duplicate initializer. The `$ODE`-side bare assignment
   (`C_MIT = CENT_MIT/V1_MIT*1000.0;`) is untouched, since a `double`
   there would collide with the `$TABLE`-side `double` declaration for the
   same reason described in the guide.

**Verification:** `Rscript -e 'parse("hsph_mrgsolve_model_refactored.R")'`
succeeds with no error. Extracted DSL posted to qspserver's `mrgsolve_api`
`/model_manifest` compiled cleanly with `C_MIT` listed in `outputPaths`
and `EC50_MIT` in the parameter manifest (both pre-edit and post-edit DSL
compiled without error, and both already had `DOSE_MIT` in the parameter
manifest).

`/run_simulation` was run against both the pre-edit (`git show HEAD:...`)
and post-edit DSL over 30 days (`end=30, delta=0.25`) with `DOSE_MIT=50`
(mitapivat's own BID pulse-input dosing parameter, per the file's
"Edge case — continuous infusion input, not an event" archetype — this
compound has no `ev()`/`data_set()` dosing scenario in the file at all;
`DOSE_MIT` is a `$PARAM` that gates a zero-order pulse computed inside
`$ODE` from `SOLVERTIME`). **Result: exact match, max absolute diff = 0
for every point of `C_MIT`, `OCC_MIT`, `EFFECT_MIT_ATP`, `EFFECT_MIT_DPG`,
`CENT_MIT`, and `Hb`, including every synthetic duplicate row** — no
dose-instant reporting divergence here, because this compound's dosing is
a continuous solver-driven pulse into the depot compartment (`GUT_MIT`),
not a discrete NM-TRAN dose event applied to the compartment `C_MIT`
reads, so the pre-dose/post-dose row-duplication mechanism that produces
the artifact for directly-dosed compounds elsewhere in this batch
(`hypercalcemia-of-malignancy`) does not apply to this dosing style.
