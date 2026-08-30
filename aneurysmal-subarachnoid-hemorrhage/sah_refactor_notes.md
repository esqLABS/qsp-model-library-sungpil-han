# Refactor notes — `sah_mrgsolve_model.R` (calibration-follow-on run: nimodipine only)

**Scope of this pass.** Per the fork's PK/PD refactor spec
([`FORK_WORKFLOW_GUIDE.md`](../FORK_WORKFLOW_GUIDE.md), Part 2), only
nimodipine's ("NIM") PK block and its interface into the disease equations
were rewritten, per its existing row in
[`driver-patches/data/compound_perturbation_census.md`](../driver-patches/data/compound_perturbation_census.md)
("Redirect concentration (clean single site)"). Every other compound in the
file (clazosentan `CLAZ`, cilostazol `CILO`, simvastatin `STAT`, milrinone
`MILRC`, intrathecal nicardipine `NICA`, ketamine `KETA`) was copied
verbatim — confirmed by diff, see below.

## Archetype determination

**Archetype 3 — depot + central, linear elimination** (no peripheral
compartment), with one structural wrinkle already present in the original
and preserved unchanged: the "central" compartment stores the plasma
**concentration** directly (ng/mL), not an amount — the volume division is
folded into the depot→central absorption term and the IV-infusion term
instead of a separate `C = CENT/V1` step. This is the same pattern already
established for tocilizumab in `sepsis/sep_mrgsolve_model_refactored.R`
(`C_TCZ = CENT_TCZ`, undivided) — see that file's row in the census and its
`sep_refactor_notes.md` for the identical rationale.

- PK: `NIMG` (gut depot, mg) → absorbed at rate `KANIM`, bioavailability
  `FNIM`, into `NIMC` (plasma concentration, ng/mL), eliminated at rate
  constant `KELNIM` (t1/2 1.7 h) scaled by the shared subject covariate
  `CLF` (a clearance multiplier applied to several compounds' elimination in
  this file — not NIM-specific, left untouched). No peripheral compartment.
  No TMDD — the disease effect is a plain saturable ratio
  `NIMC / (NIMC + EC50NIM)`, reused (with four different `EMXNIM_*` ceilings)
  across four separate downstream sites. This matches Archetype 3 minus its
  peripheral compartment, i.e. the guide's stated "drop `GUT_TCZ`/`KA_TCZ`/
  `F_TCZ` for a 2-compartment-no-depot variant" collapsed the other way
  (depot + single "central" compartment, no peripheral).
- Dosing: oral nimodipine is a genuine mrgsolve dosing event (`ev(amt=...,
  cmt="NIMG", ii=1/6, addl=...)`, built by the file's own `sah_events()` — a
  real q4h regimen, not a smooth-window covariate). IV nimodipine, like the
  file's other IV/continuous agents (clazosentan, milrinone), is the
  guide's "continuous infusion input, not an event" edge case: `NIMIV` is a
  rate parameter added directly into `dxdt_CENT_NIM`, gated by `ON =
  (SOLVERTIME between TSTART and TSTOP)`. Both preserved as-is.

## Renaming applied (nimodipine only)

| Original | Refactored | Value |
|---|---|---|
| `NIMG` (compartment) | `GUT_NIM` | — |
| `NIMC` (compartment) | `CENT_NIM` | — |
| `KANIM` | `KA_NIM` | 36.0 (1/d) |
| `KELNIM` (rate constant, 1/d) + `VNIM` | `CL_NIM` (new, = `KELNIM x VNIM`) + `V1_NIM` | `CL_NIM = 979.0` (L/d), `V1_NIM = 100.0` (L) |
| `FNIM` | `F_NIM` | 0.13 |
| `EC50NIM` | `EC50_NIM` | 30.0 (ng/mL) |
| — (none; original had no Hill exponent) | `GAMMA_NIM` (new, `= 1.0`) | documented as "original had no explicit Hill term" |
| `EMXNIM_SP` | `EMAX_NIM_SP` | 0.16 |
| `EMXNIM_RM` | `EMAX_NIM_RM` | 0.45 |
| `EMXNIM_SD` | `EMAX_NIM_SD` | 0.55 |
| `EMXNIM_NP` | `EMAX_NIM_NP` | 0.50 |
| `ENIM` (local `$ODE` ratio) | `EFFECT_NIM` | — |
| `DMAP_NIM_PO`, `DMAP_NIM_IV`, `NIMPO`, `NIMIV` | unchanged | route-specific MAP-effect constants and dosing-driver switches, out of the naming convention's role table (no matching row), and touching `NIMPO`/`NIMIV` would have broken every scenario definition (`sah_scenario_defs`) and `sah_events()` that reference them by name |

**`CL_NIM` is a new derived parameter, not an invented one.** The original
had no explicit clearance/volume split for nimodipine — elimination was
coded as a bare rate constant `KELNIM` (t1/2-derived). Per the guide's
Archetype 2 note ("always rewrite to CL/Q/V ... the more interpretable, more
common convention across this corpus"), extended here to Archetype 3/1-style
single-compartment PK for naming consistency with the rest of the corpus:
`CL_NIM := KELNIM x VNIM = 9.79 x 100 = 979.0`, so `CL_NIM / V1_NIM` in the
refactored `$ODE` reduces to exactly `9.79` — bit-identical to the original
`KELNIM` term for every timestep. No parameter *value* was invented; `CL_NIM`
is an exact arithmetic recombination of the original's own `KELNIM` and
`VNIM`.

**One shared `EFFECT_NIM`, four site-specific ceilings — not four separate
Hill terms.** Nimodipine's original effect term (`ENIM = NIMC/(NIMC+EC50NIM)`)
is a single concentration-response fraction with no Hill exponent, reused
identically at four disease sites (`vasod_sp`, `vasod_rm`, `sdrate`, `npro`),
each multiplied by its own `EMXNIM_*` ceiling. Rather than manufacture four
distinct `EFFECT_NIM_SP/RM/SD/NP` variables (which would not match how the
original computes the ratio exactly once), `EFFECT_NIM` is kept as the single
shared fractional term — exactly mirroring how this file's own `ENIM` already
worked, and consistent with the guide's multi-drug rule ("keep each
compound's `EFFECT_<STEM>` separate; combine them only at the point the
disease equations actually use them") applied here to multiple *sites* of one
compound instead of multiple compounds:
```c
double C_NIM      = CENT_NIM;
double EFFECT_NIM = pow(C_NIM, GAMMA_NIM) / (pow(EC50_NIM, GAMMA_NIM) + pow(C_NIM, GAMMA_NIM));
```
with `GAMMA_NIM = 1.0` reproducing the original's plain ratio exactly. Each
of the four usage sites was updated in place (`EMXNIM_SP*ENIM` →
`EMAX_NIM_SP*EFFECT_NIM`, etc.), including the two `map_t` terms
(`DMAP_NIM_PO`/`DMAP_NIM_IV`) and `npro` (infarct-conversion block) — every
site that read `ENIM` in the original now reads `EFFECT_NIM`, and nothing
else on those lines changed.

## Diff scope confirmation

`diff -u sah_mrgsolve_model.R sah_mrgsolve_model_refactored.R` touches: the
four nimodipine `$PARAM` groups (PK block, PD block), the `$CMT` line pair,
the `ENIM`→`C_NIM`/`EFFECT_NIM` definition in `$ODE` section 0, the
`dxdt_NIMG`/`dxdt_NIMC` lines, the four downstream usage sites
(`vasod_sp`, `vasod_rm`, `sdrate`, the two `map_t` lines, `npro`), two new
`$TABLE` lines (recompute + capture, additive only — see below), and one
R-side line in `sah_events()` (`cmt = "NIMG"` → `cmt = "GUT_NIM"`, needed so
the refactored file's own scenario-running helpers still dose the right
compartment by name). Nothing else in the file differs — every other
compound's `$PARAM`/`$CMT`/`$ODE`/`$TABLE` lines, the haemodynamic/injury
equations, the virtual population, the 15 scenario definitions, the 10
knockout definitions, and all twelve `sah_*()` result functions are
byte-identical.

## Verification

**Method.** Both `sah_mrgsolve_model.R`'s and
`sah_mrgsolve_model_refactored.R`'s embedded mrgsolve `code <- '...'` DSL
blocks were mechanically extracted (regex on the `sah_code <- '...'`
assignment, verbatim quoted text, no character changes) and POSTed to the
local qspserver `mrgsolve_api` service at `http://localhost:8007`
(`POST /model_manifest` then `POST /run_simulation`), which compiles and
runs each DSL block directly with mrgsolve 2.0.1 server-side — no local
R/mrgsolve install used. No pre-existing upstream compile defect was found
in this file (unlike some of the other calibration files) — both DSL blocks
compiled and ran without any workaround.

Three of the file's own scenarios were run through both models, identical
dosing/parameters, full `$CAPTURE` comparison:

1. **S1 — "supportive only (no nimodipine)"** (all NIM parameters at
   default 0): 14 shared captures, 211 timepoints (`end=21, delta=0.1`).
2. **S2 — "SoC: oral nimodipine 60 mg q4h x21d"** (`NIMPO=60, TSTART=0.5,
   TSTOP=21`, dosing event `amt=60` into compartment 1 (`GUT_NIM`/`NIMG`),
   `ii=1/6, addl=123` — reproducing the file's own `sah_events()` exactly):
   29 shared captures, 212 timepoints.
3. **S3 — "IV nimodipine 2 mg/h x14d"** (`NIMIV=2, TSTART=0.5, TSTOP=14`,
   no dosing event — exercises the continuous-infusion-input path in
   `dxdt_CENT_NIM` and the `DMAP_NIM_IV` route): 29 shared captures, 211
   timepoints.

**Result: exact match in all three scenarios. Maximum relative deviation
observed = 0.0 (bit-identical) across every shared `$CAPTURE`d output over
the full time grid.** This is the expected outcome for Archetype 3 (pure
structural reorganization, no Hill-fitting) per the guide's tolerance table.

## qspserver `/model_manifest` discoverability

- `KA_NIM`, `CL_NIM`, `V1_NIM`, `F_NIM`, `EC50_NIM`, `GAMMA_NIM`,
  `EMAX_NIM_SP`, `EMAX_NIM_RM`, `EMAX_NIM_SD`, `EMAX_NIM_NP` are all declared
  in `$PARAM`, so `/model_manifest` lists them with their defaults —
  confirmed against the live manifest response (e.g. `CL_NIM=979`,
  `EC50_NIM=30`, `GAMMA_NIM=1`).
- `C_NIM` and `EFFECT_NIM` themselves are computed in `$ODE` from the
  compartment state, so (as in `ted_refactor_notes.md`'s identical finding
  for `EFFECT_TCZ`) they cannot also be fixed `$PARAM` entries.
  **First attempt failed to compile**: adding `capture C_NIM = ...;` /
  `capture EFFECT_NIM = ...;` to `$TABLE` (reusing the exact `$ODE` names)
  produced `error: redefinition of 'capture {anonymous}::C_NIM'` /
  `'EFFECT_NIM'` — mrgsolve hoists every `$ODE` `double` declaration *and*
  every `capture` into the same anonymous-namespace scope, so a `$TABLE`
  capture cannot reuse a name already declared as a `double` in `$ODE`.
  (This is also *why* this file's own author already used a systematic
  `o`-prefix in `$TABLE`, e.g. `oSPT` vs `$ODE`'s `SPT` — same collision,
  pre-existing pattern, now understood rather than just imitated.) Fixed by
  naming the captures `NIM_conc` / `NIM_effect` instead (recomputed from
  state via new `$TABLE` locals `oC_NIM`/`oEFFECT_NIM`), confirmed present
  in `/model_manifest`'s `outputPaths` and retrievable via `/run_simulation`.
  This mirrors `ted_refactor_notes.md`'s `capture TCZ_effect = EFFECT_TCZ;`
  workaround for the same underlying compile error.
- Both new capture lines are additive only (the original never named these
  quantities) and are excluded from the verification diff above, which
  compares only the original's own pre-existing 29 captures.
- `model_content` here is already pure DSL text extracted from the
  `code <- '...'` R-string wrapper (this file uses that pattern, not the
  `[PROB]`/`mread()` file-native pattern some other files in this corpus
  use) — no separate `.cpp` sibling was left behind; the extraction was
  in-memory only, used to build the verification requests above and then
  discarded, per the workflow guide.

## Anything else flagged

- No pre-existing upstream compile defect was found in this file requiring
  an `UPSTREAM_ISSUES.md` entry — both the original and refactored DSL
  compiled cleanly on the first attempt (after fixing the capture-naming
  collision described above, which is a property of the *refactor*, not a
  pre-existing defect in the original).
- `CLF` (a clearance multiplier covariate shared by nimodipine, clazosentan,
  cilostazol, and milrinone) was left completely untouched, exactly as
  `WT`/other subject-level covariates were left untouched in the
  `kidney-transplant-rejection` calibration run — it is not nimodipine-
  specific and renaming it was out of scope.
- `DMAP_NIM_PO`/`DMAP_NIM_IV` (nimodipine's induced-hypotension side effect
  on `map_t`) are outside the guide's naming-convention role table (no
  "MAP effect" role listed) and were left unchanged, consistent with how the
  `kidney-transplant-rejection` run left `TCZIN`/`TCZON`/`TCZSTART`/`TCZN`
  unchanged for the same reason.
