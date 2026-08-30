# Refactor notes — `ald_mrgsolve_model_refactored.R`

Scope: **G-CSF (GCSF), NAC, and Prednisolone (PRED) only** — their PK
blocks and their downstream disease effect equations. Pentoxifylline,
anakinra, and every disease equation not driven by one of the three
in-scope compounds are byte-for-byte copies of `ald_mrgsolve_model.R`,
apart from two disclosed, non-numeric mrgsolve-2.0.1 build-compatibility
fixes described below (needed just to get the *original* to compile at
all, unrelated to the three compounds' own math).

## Archetype determination

Read the actual `dxdt_` lines rather than assuming from the file's own
"Drug PK" section comments.

- **Prednisolone**: 3 compartments (`PRED_gut`, `PRED_C`, `PRED_P`),
  first-order absorption (`ka_pred`, `F_pred`) into a linear two-compartment
  distribution model (`CL_pred`, `Vc_pred`, `Q_pred`, `Vp_pred`), already
  written in direct CL/Q/V form (no micro-constant decomposition to
  convert). **Archetype 3** (depot + central + peripheral, linear).
- **NAC**: 1 compartment (`NAC_C`), linear elimination
  `-(CL_NAC/Vc_NAC)*NAC_C` plus one extra first-order loss term
  `-k_NAC_tissue*NAC_C` (an added "tissue distribution" rate with no
  corresponding tissue compartment — a lumped extra clearance, not a real
  second compartment). **Archetype 1** (single compartment, linear
  elimination) plus this one extra elimination term, kept verbatim.
- **G-CSF**: 1 compartment (`GCSF_C`), linear elimination
  `-(CL_GCSF/Vc_GCSF)*GCSF_C`. **Archetype 1**. `ka_GCSF` is declared in
  `$PARAM` but never referenced anywhere in `$ODE` — a genuine dead
  parameter in the original (dosing into `GCSF_C` happens directly via an
  infusion `rate`, not through an absorption compartment) — kept, renamed,
  and flagged, not removed, matching the corpus precedent of preserving
  dead parameters found in-scope (e.g. `pmr_refactor_notes.md`'s
  `KON_TCZ`/`KOFF_TCZ`).

None of the three needed TMDD (Archetype 4) or a Hill refit — see below.

## Naming applied

| Original | Refactored |
|---|---|
| `ka_pred`, `F_pred`, `CL_pred`, `Vc_pred`, `Q_pred`, `Vp_pred` | `KA_PRED`, `F_PRED`, `CL_PRED`, `V1_PRED`, `Q_PRED`, `V2_PRED` |
| `PRED_gut`, `PRED_C`, `PRED_P` | `GUT_PRED`, `CENT_PRED`, `PERI_PRED` |
| `Emax_pred`, `EC50_pred` | `EMAX_PRED`, `EC50_PRED`, + new `GAMMA_PRED = 1` |
| `CL_NAC`, `Vc_NAC`, `k_NAC_tissue` | `CL_NAC` (unchanged), `V1_NAC`, `K_TISSUE_NAC` |
| `NAC_C` | `CENT_NAC` |
| `Emax_NAC`, `EC50_NAC` | `EMAX_NAC`, `EC50_NAC` (unchanged), + new `GAMMA_NAC = 1` |
| `ka_GCSF`, `CL_GCSF`, `Vc_GCSF` | `KA_GCSF` (unused, kept), `CL_GCSF` (unchanged), `V1_GCSF` |
| `GCSF_C` | `CENT_GCSF` |
| `Emax_GCSF`, `EC50_GCSF` | `EMAX_GCSF`, `EC50_GCSF` (unchanged), + new `GAMMA_GCSF = 1` |

`K_TISSUE_NAC` is not one of the guide's named roles, but is kept as an
explicit, stem-suffixed parameter (it was already a named `$PARAM` in the
original, `k_NAC_tissue`) rather than folded silently into `CL_NAC` — a
rename, not a restructuring.

Pentoxifylline (`Emax_pento`, `EC50_pento`, `PTX_C`) and anakinra
(`Emax_anakin`, `EC50_anakin`, `ANK_C`) are completely untouched, including
their own compartment names, parameter names, and PK equations.

## `C_<STEM>` is an identity alias, not a `/V1` division — and why

The guide's Archetype 3 template computes the exposed concentration as
`C_TCZ = CENT_TCZ / V1_TCZ` (central compartment tracked as an *amount*,
concentration derived by dividing by volume). Prednisolone's original code
does **not** do this: `PRED_C` is dosed and read as already being the
plasma concentration directly (`ka_pred*F_pred*PRED_gut/Vc_pred*1000.0`
converts the incoming dose straight into concentration units, and the
downstream Hill term `Emax_pred*PRED_C/(EC50_pred+PRED_C)` reads `PRED_C`
with no further division). Working the algebra backward: if `PRED_C` were
truly `CENT/V1` and `PRED_P` were truly `CENT2/V2` for amounts `CENT`,
`CENT2`, the correct concentration-form ODE would need
`dPRED_C/dt`'s peripheral-return term to use `Q/V1` (not `Q/Vp` as
written) and `dPRED_P/dt`'s central-inflow term to use `Q/V2` (not `Q/Vc`
as written) — the original has these two coefficients swapped relative to
a unit-consistent amount/concentration pair. In other words, `PRED_C` and
`PRED_P` are structurally amount-like state variables (their own `dxdt_`
recursion is the textbook amount-based two-compartment form) that the
original simply *reads directly* as the PD-facing concentration, without
ever dividing by a volume. This is a pre-existing, unfixed unit-consistency
quirk in the original — flagged here, not corrected (correcting it would
change the compound's numeric behavior, out of scope for a rename-only
refactor). The refactor preserves it exactly: `dxdt_CENT_PRED`/
`dxdt_PERI_PRED` are renamed byte-for-byte copies of the original's
`dxdt_PRED_C`/`dxdt_PRED_P` (same swapped coefficients, same numbers), and
`C_PRED = CENT_PRED` is a plain identity, matching what the original's own
Hill term actually reads. The same identity choice is applied to `C_NAC`
and `C_GCSF` (both single-compartment, dosed and read directly as
concentrations, no volume division anywhere in the original either).

This exact situation — a "central compartment" already stored as a
concentration, not an amount — has precedent elsewhere in this fork, e.g.
`diabetic-ketoacidosis/dka_refactor_notes.md`'s insulin PK
(`double C_INS = CENT_INS;`) and `aneurysmal-subarachnoid-hemorrhage`'s
nimodipine.

## The Hill interface

All three compounds' original effect terms are already a plain Emax ratio
(`Emax*C/(EC50+C)`), so per the guide this is **a rename, not a refit** —
`GAMMA_*_ = 1` throughout, fixed because the original had no Hill
coefficient at all:

- **Prednisolone → NF-κB inhibition**: `eff_pred_NFkB` renamed
  `EFFECT_PRED`, read at every original use site (`KC_inhib`, `TNF_inhib`,
  `IL1B_inhib`) in place of `eff_pred_NFkB`, combined with pentoxifylline's
  and anakinra's own (untouched) effect terms exactly as the original did.
- **NAC → GSH restoration**: `eff_NAC_GSH` renamed `EFFECT_NAC`, read in
  `GSH_synth_eff = k_GSH_synth * (1.0 + EFFECT_NAC)`.
- **G-CSF → neutrophil stimulation**: the original's
  `eff_GCSF_neut = 1.0 + Emax_GCSF*GCSF_C/(EC50_GCSF+GCSF_C)` is "1 + a
  plain Hill ratio." Per the guide's instruction to keep each compound's
  `EFFECT_<STEM>` as just its own Hill term and combine only at the point
  of use, `EFFECT_GCSF` holds just the ratio (`EMAX_GCSF*C^γ/(EC50^γ+C^γ)`,
  range 0..`EMAX_GCSF`), and the disease equation re-applies the "+1"
  baseline fold at its one point of use: `dxdt_NEUT = neut_drive * (1.0 +
  EFFECT_GCSF) - k_neut_clear * NEUT` — an exact algebraic identity with
  the original (`eff_GCSF_neut - 1.0 ≡ EFFECT_GCSF` by construction), not a
  behavioral change. The same identity is used in `regen_rate`, which
  originally read `(eff_GCSF_neut - 1.0)/Emax_GCSF` and now reads
  `EFFECT_GCSF/EMAX_GCSF` — numerically identical.

No curve fit was performed anywhere; `GAMMA_PRED = GAMMA_NAC =
GAMMA_GCSF = 1` makes `pow(x,1) = x` in every Hill expression.

## `$MAIN`-vs-`$TABLE` capture timing (a finding, not a numeric change)

`EFFECT_PRED`/`EFFECT_NAC`/`EFFECT_GCSF` and their underlying `C_PRED`/
`C_NAC`/`C_GCSF` are computed in `$MAIN` (matching exactly where the
original computed `eff_pred_NFkB` etc.), and the disease equations in
`$ODE` (`KC_inhib`, `TNF_inhib`, `IL1B_inhib`, `GSH_synth_eff`,
`dxdt_NEUT`, `regen_rate`) read those `$MAIN`-timed values, exactly as the
original did — this part of the refactor changes nothing about *when*
these quantities are evaluated relative to the disease dynamics.

However, adding new `$TABLE` captures (`C_PRED_out`, `EFFECT_PRED_out`,
etc., needed for qspserver `/model_manifest` discoverability — see below)
by simply reading back the `$MAIN`-scope `C_PRED`/`EFFECT_PRED` variables
turned out to under-report: mrgsolve evaluates `$MAIN` once at the
**start** of each requested interval, before that interval's `$ODE`
integration runs, so a `$TABLE` capture that read the `$MAIN` member
directly reported the *previous* interval's value — confirmed empirically
against this exact model (with `delta=1h` dosing, `capture C_PRED_out =
C_PRED;` trailed `CENT_PRED` by exactly one 1h step; verified by comparing
the two directly in a scratch run before settling on the fix below).
**Fix**: `$TABLE` recomputes `C_PRED_tbl`/`C_NAC_tbl`/`C_GCSF_tbl` and
`EFFECT_PRED_tbl`/`EFFECT_NAC_tbl`/`EFFECT_GCSF_tbl` fresh, directly from
the `$CMT` state (`CENT_PRED` etc.) at `$TABLE`'s own evaluation time, and
captures those instead. This only affects the **new** reporting-only
captures added by this refactor for API discoverability — it does not
touch the `$MAIN`/`$ODE` computation the disease equations actually run
on, which is unchanged from the original. Distinct names (`_tbl` suffix)
are used to avoid mrgsolve's `capture NAME = NAME;` self-shadowing pitfall
(same precedent as `dka_refactor_notes.md`'s `_report` suffix and
`pmr_refactor_notes.md`'s discussion of the same class of naming pressure).

## Pre-existing upstream build defects (fixed syntax-only in the delivered file)

**The original does not compile under mrgsolve 2.0.1 at all**, for two
independent reasons unrelated to any of the three refactored compounds.
Logged as [`UPSTREAM_ISSUES.md` #50](../translations/UPSTREAM_ISSUES.md).
Per the guide's settled policy for this situation ("When the original
doesn't compile at all"), both fixes are applied **directly to the
delivered `ald_mrgsolve_model_refactored.R`**, not just to a scratch copy:

1. **A `$TABLE`-local `double prob_d90` collides with its own
   `capture prob_d90 = prob_d90;` line** (mrgsolve auto-promotes the bare
   `double` assignment to a class member, then the `capture` line tries to
   redeclare the same name):
   ```
   59:11: error: redefinition of 'capture {anonymous}::prob_d90'
   56:10: note: 'double {anonymous}::prob_d90' previously declared here
   ```
   Fixed by renaming the capture to `prob_d90_out`, matching every other
   line in the same `$TABLE` block (`MELD_out`, `DF_out`, `ALT_out`, ...).
   Purely a rename — the underlying `double prob_d90` computation is
   untouched, and the R-side plotting/summary code below the model block
   was updated to read the new column name (`prob_d90` → `prob_d90_out`,
   4 sites: the mortality-risk plot, the fibrosis/mortality facet plot, and
   the day-90 summary table — all in code paths that actually execute; the
   `run_sim()` function's own `outvars` list, dead code not called by the
   file's own scenario loop, was updated too for consistency).
2. **Fourteen `$ODE` lines across twelve compartments assign directly to
   their own compartment state to clamp it** (`if(AA < 0) AA = 0;`,
   `if(H < 0.01) H = 0.01;`, `if(PRED_C < 0) PRED_C = 0;`, etc. — the full
   list covers `AA`, `GSH`, `KC`, `IL1B`, `NEUT`, `H` (×2), `F` (×2),
   `PRED_C`, `NAC_C`, `GCSF_C`, `PTX_C`, `ANK_C`). mrgsolve 2.0.1 passes
   every `$CMT` state into `$ODE` as a `const double&`, so any of these is
   a hard compile error. Same defect class as `UPSTREAM_ISSUES.md` #36
   (`age-related-macular-degeneration`): confirmed none of these twelve
   names is read again later in the same `$ODE` evaluation after its own
   clamp line, so — following #36's reasoning — only `dxdt_*` feeds the
   solver, meaning the clamp assignment itself could never have altered
   the integrated trajectory even on a build where it compiled. Fixed by
   deleting all fourteen lines outright; no replacement logic was added.
   Two of the twelve affected names (`PRED_C`, `NAC_C`) and one compound's
   PK (`GCSF_C`'s own clamp) fall inside this refactor's own scope — their
   clamp lines are deleted as part of the renamed PK blocks
   (`CENT_PRED`/`CENT_NAC`/`CENT_GCSF` no longer have a clamp line at all,
   matching the twelve-compartment fix uniformly); the other nine
   (`AA`, `GSH`, `KC`, `IL1B`, `NEUT`, `H`×2, `F`×2, `PTX_C`, `ANK_C`) are
   disease-side or out-of-scope-compound state, fixed only because the
   whole file must compile as one unit — no other change was made to any
   of those nine compartments' equations.

Both fixes are syntax-only and non-numeric — proven by the verification
below, which compares the original (with only these same two patches
applied to an in-memory scratch copy, never to the checked-in
`ald_mrgsolve_model.R`) against the refactored file and finds exact
(0.0) agreement on every shared output. The checked-in original is
untouched and still carries both defects exactly as written.

## qspserver `/model_manifest` discoverability

Confirmed via `POST /model_manifest` on the refactored DSL
(`http://localhost:8007`): all 21 new/renamed parameters are listed with
their original values —
`KA_PRED=1.5, F_PRED=0.82, CL_PRED=8.5, V1_PRED=35, Q_PRED=15, V2_PRED=50,
EMAX_PRED=0.75, EC50_PRED=150, GAMMA_PRED=1, CL_NAC=12, V1_NAC=30,
K_TISSUE_NAC=0.4, EMAX_NAC=0.65, EC50_NAC=80, GAMMA_NAC=1, KA_GCSF=0.5,
CL_GCSF=0.8, V1_GCSF=4.5, EMAX_GCSF=1.8, EC50_GCSF=5, GAMMA_GCSF=1` — and
`outputPaths` lists `GUT_PRED`, `CENT_PRED`, `PERI_PRED`, `CENT_NAC`,
`CENT_GCSF` (compartments) plus `C_PRED_out`, `EFFECT_PRED_out`,
`C_NAC_out`, `EFFECT_NAC_out`, `C_GCSF_out`, `EFFECT_GCSF_out` (the six new
reporting captures). `C_<STEM>`/`EFFECT_<STEM>` themselves are
state-dependent (recomputed every step from the `$CMT` state), so — per the
precedent in `chronic-hypothyroidism/hypo_refactor_notes.md` and
`diabetic-ketoacidosis/dka_refactor_notes.md` — they are discoverable as
`$CAPTURE` outputs, not `$PARAM` values; attempting to also declare them in
`$PARAM` would collide with their own `$ODE`/`$MAIN`/`$TABLE`
recomputation (confirmed by that same precedent, not re-tested here since
the mechanism is identical).

## Verification

**Method.** Both `ald_mrgsolve_model.R`'s embedded DSL (with the two
build-compat patches above applied to an in-memory scratch copy only) and
`ald_mrgsolve_model_refactored.R`'s embedded DSL (patches applied directly,
as delivered) were extracted mechanically (find the `code <- '...'`/
`ald_code <- '...'` assignment, pull the quoted contents verbatim) and
POSTed to the qspserver `mrgsolve_api` container (`http://localhost:8007`,
confirmed healthy) — `/model_manifest` first (both compile), then
`/run_simulation`.

**Note on a mid-verification infrastructure outage**: partway through this
work the shared `mrgsolve_api` container started returning `500` errors
(`readRDS: ReadItem: unknown type 47, perhaps written by later version of
R`) for every request, including trivial throwaway models — a corrupted
entry in the container's persistent build-cache volume (`mrgsolve_soloc`,
mounted at `/soloc`, explicitly retained "across restarts" per
`docker-compose.yaml`), evidently written by a different, incompatible R/
mrgsolve version from another concurrent session sharing this container.
Restarting/force-recreating the container did not fix it, because the
corrupted cache entries persist in the volume across container restarts.
Once the container was confirmed healthy again, the specific content
hashes already submitted during the broken window remained permanently
poisoned (mrgsolve's cache tries to `readRDS` the existing, corrupt file
for that exact hash rather than rebuilding); appending a harmless
no-op comment (a "cache-bust nonce") to each DSL text before resubmitting
produced a fresh, never-cached hash and avoided the poisoned entries
without touching either file's actual content. This is disclosed here in
case it recurs for a later session sharing the same container; it is an
infrastructure artifact, not a finding about either model.

Ran **3 dosing scenarios**, reproducing the original file's own
`build_events()` definitions (identical amounts/timing; `build_events()`
is itself dead code never invoked by the original's own "Run all 7
scenarios" pipeline further down the file, but it is the only place in the
file with real dosed PK events for these three compounds — the executed
pipeline instead bypasses PK entirely via a fixed steady-state proxy
concentration and an Emax on/off switch, see "Other findings" below):

| Scenario | Description | Dosing (identical, both models) |
|---|---|---|
| S3 | Prednisolone 40 mg QD × 28 days | `amt=40` into the gut compartment (cmt 15) at `t = 0, 24, ..., 648` h |
| S4 | NAC IV loading + maintenance | `amt=350, rate=17.5` at t=0 (20 h infusion); `amt=117, rate=4.9` at t=24,48,72,96 (per-dose ~23.9 h infusions) into the central compartment (cmt 18) |
| S6 | G-CSF 350 μg SC × 5 days | see note below |

**S6 dosing note**: the original's own `build_events()` computes the
infusion rate as `rate = 350 / Vc_GCSF_dose` — `Vc_GCSF_dose` is not
defined anywhere in the file (a pre-existing dead-code bug; this branch is
never executed by the file's own pipeline, so the bug never surfaces at
runtime). Since this is an R-side dosing-helper bug, not a DSL defect,
substituted a well-defined bolus dose (`amt=350, rate=0`, once daily × 5
days) of the same documented total daily amount for verification purposes
only — both the original(-fixed) and refactored models received the
identical substituted dosing, so the head-to-head comparison remains
valid regardless of this substitution. Not fixed in either file; the
original's `build_events()` still contains the undefined-variable
reference exactly as written (only its `cmt=` string was updated in the
refactored sibling, from `"GCSF_C"` to `"CENT_GCSF"`, to track the renamed
compartment — the dead code otherwise untouched).

For each scenario, every disease-side `$CAPTURE`d output shared by both
models (`MELD_out`, `DF_out`, `prob_d90_out`, `ALT_out`, `BILI_out`,
`INR_out`, `H_out`, `GSH_out`, `ROS_out`, `KC_out`, `NEUT_out`, `F_out`,
`TNF_out`, `IL1B_out`) plus the dosed compound's own PK output were
compared across the full time grid:

| Scenario | Grid | Compound PK output | max abs diff | max rel diff |
|---|---|---|---|---|
| S3 (Prednisolone) | 749 pts, 0–720 h, Δ=1h | `PRED_C_out` vs `C_PRED_out` (peak 542.60) | 0.0 | 0.0 |
| S4 (NAC) | 342 pts, 0–168 h, Δ=0.5h | `NAC_C_out` vs `C_NAC_out` (peak 21.88) | 0.0 | 0.0 |
| S6 (G-CSF) | 198 pts, 0–192 h, Δ=1h | `GCSF_C_out` vs `C_GCSF_out` (peak 354.98) | 0.0 | 0.0 |

All fourteen shared disease-side outputs (`MELD_out` through `IL1B_out`)
also matched **exactly (0.0 deviation)** in all three scenarios — including
`KC_out`/`TNF_out`/`IL1B_out` (driven by `EFFECT_PRED`) in S3, `GSH_out`
(driven by `EFFECT_NAC`) in S4, and `NEUT_out` (driven by `EFFECT_GCSF`)
in S6 — confirming the renamed Hill interfaces reproduce the original's
disease-level impact exactly, not just the isolated PK compartments.
Sanity range check on the new `EFFECT_*_out` captures: `EFFECT_PRED_out`
0–0.588 (< `EMAX_PRED`=0.75), `EFFECT_NAC_out` 0–0.140 (< `EMAX_NAC`=0.65),
`EFFECT_GCSF_out` 0–1.775 (< `EMAX_GCSF`=1.8) — all within the expected
Hill-ratio bound, none degenerate.

**Result: verification passed, exact match as expected for a pure
rename/reorganization with no Hill-fitting** (Archetypes 1 and 3, gamma
fixed at 1 throughout — this is the "near-exact match expected" case per
the guide, and in practice came out bit-identical rather than merely
close).

## Other findings flagged for manual review (not fixed, out of scope)

1. **`build_events()` is dead code.** Neither `run_sim()` nor the file's
   own executed "Run all 7 scenarios" loop ever calls it. The executed
   pipeline instead sets a fixed steady-state proxy concentration as an
   *initial condition* (e.g. `CENT_PRED = 200` when prednisolone is "on")
   and never doses the PK compartments via events at all — so the file's
   own headline scenario results were never actually produced by
   integrating the PK ODEs, only by the disease equations reading a static
   assumed concentration. This applies to Prednisolone and NAC's own
   "on/off" switches (parameter override + fixed initial concentration).
2. **G-CSF's scenario "on/off" (S6) has no numeric effect in the executed
   pipeline**, a real modeling no-op independent of this refactor: `S6`
   sets `EMAX_GCSF=1.8` (vs `1.0` elsewhere) but never sets
   `CENT_GCSF`/`GCSF_C`'s initial value away from its `$INIT` default of
   `0`. Since `EFFECT_GCSF = EMAX_GCSF * C^γ/(EC50^γ+C^γ)` is `0` whenever
   `C=0` regardless of `EMAX_GCSF`, S6 produces identical
   `dxdt_NEUT`/`regen_rate` behavior to every other Emax-off scenario in
   the same pipeline. Confirmed by inspecting the actual `init()` call in
   the "Run all 7 scenarios" loop: `GCSF_C_ss <- 0` is computed but never
   passed to `init()`. Not fixed — out of scope for a rename-only refactor,
   and changing it would alter the original's numeric behavior.
3. **`ka_GCSF` is a dead parameter** in the original (declared, never used
   in any `$ODE`/`$MAIN` expression) — G-CSF is dosed directly into the
   central compartment via an infusion `rate`, with no absorption
   compartment. Kept, renamed `KA_GCSF`, and flagged rather than removed
   (pruning dead parameters was not part of the assigned scope).
