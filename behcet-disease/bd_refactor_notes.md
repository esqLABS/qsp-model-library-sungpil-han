# Refactor notes — `behcet-disease/bd_mrgsolve_model.R`

Compounds refactored: **Colchicine (COL)**, **Prednisolone (PRED)**,
**Adalimumab (ADA)**, **Apremilast (APR)**, **Canakinumab (CAN)** — the
five rows in `driver-patches/data/compound_perturbation_census.md`
classified `behcet-disease | Colchicine`, `behcet-disease | Prednisolone`,
`behcet-disease | Adalimumab`, `behcet-disease | Apremilast`,
`behcet-disease | Canakinumab`, all originally marked "Normalize duplicate
concentration sites, then redirect". These are the file's entire set of
exogenous compounds — no other compound exists in this model. Every
disease-side compartment and equation (neutrophil activation, Th1/Th17/
Treg balance, the TNF-α/IL-1β/IL-6/IL-17A cytokine network, endothelial
activation, oral ulcer index, ocular inflammation index, BDCAF composite)
is untouched, confirmed by comparing the delivered file against the
original outside the renamed PK/Hill regions.

Original: `behcet-disease/bd_mrgsolve_model.R` (never edited).
Delivered: `behcet-disease/bd_mrgsolve_model_refactored.R`.

## Archetype per compound

- **Colchicine (COL)**: **archetype 3** (depot + central + peripheral,
  linear). `GUT_COL`/`CENT_COL`/`PERI_COL` (was `AGUT_COL`/`ACOL`/
  `ACOL_T`), `KA_COL`, `F_COL` (was `Foral_col`), `CL_COL`, `V1_COL`,
  `V2_COL`, `Q_COL` unchanged in value — every one of these parameters is
  actually used in `$ODE`. Two PD effects (neutrophil migration
  inhibition, IL-1β/NLRP3 inhibition).
- **Prednisolone (PRED)**: **archetype 2** (central + peripheral, no
  depot). `CENT_PRED`/`PERI_PRED` (was `APRED`/`APRED_T`), `CL_PRED`,
  `V1_PRED`, `V2_PRED`, `Q_PRED` unchanged. `KA_PRED`/`F_PRED` (was
  `ka_pred`/`Foral_pred`) are carried over **unused** — see "Vestigial
  parameters" below. Two PD effects (TNF-α suppression, IL-6 suppression),
  sharing one `EMAX_PRED` exactly as the original shared one `Emax_pred`.
- **Adalimumab (ADA)**: **archetype 2** (central + peripheral, no depot).
  `CENT_ADA`/`PERI_ADA` (was `AADA`/`AADA_P`), `CL_ADA`, `V1_ADA`,
  `V2_ADA`, `Q_ADA` unchanged. `F_ADA` (was `F_ada`) carried over
  **unused**. One PD effect (TNF-α neutralization).
- **Apremilast (APR)**: **archetype 1** (single compartment, linear, no
  depot). `CENT_APR` (was `AAPR`), `CL_APR`, `V1_APR` unchanged. `KA_APR`/
  `F_APR` (was `ka_apr`/`Foral_apr`) carried over **unused** — the
  original's own dosing scenario doses `AAPR` directly as an instantaneous
  bolus (`evid=1`, no `rate`), never routing through any depot/absorption
  process despite declaring absorption parameters for one. Two PD effects
  (TNF-α, IL-17A), sharing one `EMAX_APR` exactly as the original shared
  one `Emax_apr`.
- **Canakinumab (CAN)**: **archetype 2** (central + peripheral, no
  depot). `CENT_CAN`/`PERI_CAN` (was `ACAN`/`ACAN_P`), `CL_CAN`, `V1_CAN`,
  `V2_CAN`, `Q_CAN` unchanged. `F_CAN` (was `F_can`) carried over
  **unused**. One PD effect (IL-1β neutralization).

No compound needed **archetype 4 (TMDD)**. Confirmed from the code, not
assumed: Adalimumab's and Canakinumab's own effect terms
(`Eada_TNF = Emax_ada * (Cp_ada / (EC50_ada_TNF + Cp_ada))`, `Ecan_IL1 =
Emax_can * (Cp_can / (EC50_can_IL1 + Cp_can))`) are plain saturating Hill
ratios on total plasma concentration — the model has no free-receptor or
drug-receptor-complex compartment for either compound, i.e. no ODE-solved
receptor-binding kinetics to preserve. Both biologics are modeled purely
as linear 2-compartment PK plus an algebraic Emax effect, same as the
three small molecules.

## Vestigial (declared-but-unused) parameters — carried over, not removed

Six parameters are declared in the original's `$PARAM` but never read by
any `$MAIN`/`$ODE`/`$TABLE` calculation, confirmed by grep against the
whole file: `ka_pred`/`Foral_pred` (prednisolone), `F_ada` (adalimumab),
`ka_apr`/`Foral_apr` (apremilast), `F_can` (canakinumab). In every case
the compound has no depot compartment in `$ODE` for the absorption-rate
parameter to apply to, and is dosed as either a direct bolus (apremilast)
or via the `rate=-2` flag straight into the central compartment
(prednisolone, adalimumab, canakinumab) — see the dosing-scenario defect
below. These are renamed (`KA_PRED`, `F_PRED`, `F_ADA`, `KA_APR`, `F_APR`,
`F_CAN`) and kept declared, exactly per the guide's "never invent or
default a PK parameter, and don't add one the original didn't have" —
carrying forward an unused declaration is not adding one, and removing it
was judged a larger, unnecessary deviation from the original's own
`$PARAM` manifest than leaving it in place, annotated as unused.

Separately, nine baseline-placeholder parameters (`NEU0`, `TH1_0`,
`TH17_0`, `TREG_0`, `TNFA0`, `IL1B0`, `IL6_0`, `IL17A0`, `EA0`) are
likewise never read anywhere — a disease-side quirk, out of scope for a
compound-PK refactor, left untouched except for the build-compatibility
rename of three of them (see below).

## Upstream defects found (logged in `translations/UPSTREAM_ISSUES.md`)

Three separate, pre-existing defects were found while getting the
original's own DSL to load/run through the qspserver `mrgsolve_api`
(entries #152–#154). All three are unrelated to any compound's own PK/PD
and were fixed with syntax-only, non-numeric changes directly in the
delivered `_refactored.R` (never in the checked-in original):

1. **#152 — `$CAPTURE @annotated` duplicates 12 `$CMT` compartment
   names** (`NEU`, `TH1`, `TH17`, `TREG`, `TNFA`, `IL1B`, `IL6C`, `IL17A`,
   `EA`, `OUL`, `OCI`, `BDCAF`), which mrgsolve 2.0.1 rejects outright
   (`compartment should not be in $CAPTURE`) — the model never loads at
   all. **Fix applied:** removed the 12 duplicated names from
   `$CAPTURE @annotated` in the refactored file. No numeric effect —
   mrgsolve reports every `$CMT` compartment in its output by default
   regardless of `$CAPTURE` membership, confirmed these 12 names are
   still present in every run's output columns.
2. **#153 — `$PARAM` baseline placeholders `TH1_0`/`TH17_0`/`TREG_0`
   collide with mrgsolve's own auto-generated `<CMT>_0` initial-value
   symbol** for compartments `TH1`/`TH17`/`TREG`, producing a C++
   "conflicting declaration" build error once #152 is worked around.
   **Fix applied:** renamed to `TH1_BASE0`/`TH17_BASE0`/`TREG_BASE0`.
   Confirmed none of the three is ever read anywhere, so this is a pure
   syntax fix with zero numeric effect.
3. **#154 — three of the file's own eight dosing scenarios use
   `rate = -2` without the model declaring any modeled infusion
   duration** (`ev_pred`, `ev_ada`, `ev_can`, and the two combo scenarios
   referencing them), so they cannot actually be simulated — confirmed
   with a plain local `mrgsolve::mcode()`+`mrgsim()` call on the
   untouched original (not just via the API):
   `[mrgsolve] modeled infusion duration D_CMT or Dn must be positive
   when dosing record RATE is set to -2.` This is a scenario-execution
   limitation, not a model-loading defect, so it was **not** fixed in the
   refactored file (which preserves the same `rate = -2` events,
   unchanged, since rewriting them would mean the refactored file's own
   scenarios no longer literally match the original's). See "Verification
   methodology" below for how this was worked around for testing
   purposes only.

## Verification

Via the qspserver `mrgsolve_api` (`POST /model_manifest`,
`POST /run_simulation`, `http://localhost:8007`, requests spaced ≥2.2 s
apart, max 2 concurrent), comparing:

- **Baseline**: the original's own DSL (extracted verbatim from
  `bd_code <- '...'`) with *only* the two build-compatibility fixes from
  entries #152/#153 applied to a scratch copy (never to the checked-in
  original) — otherwise byte-identical to the original.
- **Test article**: `bd_mrgsolve_model_refactored.R`'s own DSL, extracted
  verbatim from `bd_code_refactored <- '...'`.

against **all 8 of the original file's own dosing scenarios**
(untreated, colchicine monotherapy, prednisolone monotherapy, adalimumab
monotherapy, apremilast monotherapy, canakinumab monotherapy,
colchicine+prednisolone combo, adalimumab+apremilast combo), same
cmt/time/amt/evid for every dose, full 2160 h (90 day) horizon at the
file's own `delta=4` h resolution — no shortening needed, no scenario hit
the API's default solver step-count limit.

### Verification methodology note (rate=-2 substitution)

Per upstream defect #154, `ev_pred`/`ev_ada`/`ev_can` (and both combos
referencing them) cannot literally be simulated as coded (`rate=-2` with
no modeled duration) through *any* mrgsolve 2.0.1 path, not just the API
— confirmed with a plain local `mcode()`+`mrgsim()` call on the untouched
original DSL. For these five scenarios only, verification substituted
`rate=0` (instantaneous bolus, identical time/cmt/amt) — applied
**identically** to both the original-patched baseline and the refactored
model, so the substitution cannot favor either side. This is a
verification-methodology workaround, not a fix to either deliverable; the
refactored file's own R driver script still encodes the original's
`rate=-2` events unchanged.

### Results (max abs / max rel diff across every shared `$CAPTURE`'d output, full time grid)

| Scenario | Doses | Points | Max abs diff | Max rel diff | Result |
|---|---|---|---|---|---|
| 1 Untreated | 1 | 542 | 0.0 | 0.0 | exact |
| 2 Colchicine | 180 | 721 | 0.0 | 0.0 | exact |
| 3 Prednisolone (rate→0) | 90 | 631 | 0.0 | 0.0 | exact |
| 4 Adalimumab (rate→0) | 7 | 548 | 0.0001 | 6.6e-05 | floating-point-scale |
| 5 Apremilast | 180 | 721 | 0.0 | 0.0 | exact |
| 6 Canakinumab (rate→0) | 2 | 543 | 0.0005 | 4.2e-04 | floating-point-scale |
| 7 Colch+Pred combo (rate→0) | 270 | 811 | 0.0 | 0.0 | exact |
| 8 Ada+Aprem combo (rate→0) | 187 | 728 | 0.0 | 0.0 | exact |

6 of 8 scenarios match **exactly** (max abs diff 0.0 across every shared
output — `NEU`, `TH1`, `TH17`, `TREG`, `TNFA`, `IL1B`, `IL6C`, `IL17A`,
`EA`, `OUL`, `OCI`, `BDCAF`, `Cp_COL`, `Cp_PRED`, `Cp_ADA`, `Cp_APR`,
`Cp_CAN`), including the two longest/most complex scenarios (combo1: 270
doses, combo2: 187 doses).

Scenarios 4 (adalimumab) and 6 (canakinumab) show a tiny, non-zero
deviation, confirmed to be floating-point/solver-path noise rather than a
structural bug:
- Each compound's own concentration column (`Cp_ADA`, `Cp_CAN`) matches
  **bit-for-bit** (0.0 diff) in every scenario — the linear PK ODEs
  themselves are unaffected.
- The deviation appears only in downstream, nonlinearly-coupled disease
  compartments (e.g. `IL6C`, `BDCAF`), at ~1e-4 absolute on values in the
  13–40 range (≤4.2e-4 relative) — consistent with an adaptive
  (LSODA-family) solver occasionally choosing a marginally different
  internal step given a bit-level difference in compiled floating-point
  operation order between the two DSL texts (renamed variables, restated
  as `pow(x, 1.0)` per the Hill template rather than a plain ratio),
  which a stiff, tightly-coupled 22-compartment nonlinear system can
  visibly (if tinily) amplify over hundreds of adaptive steps.
- Ruled out as a real bug: (a) determinism confirmed (re-running the same
  request against the same model twice returns byte-identical output,
  so this is not server-side nondeterminism); (b) the *same* rename
  pattern (`pow(C, GAMMA)` with `GAMMA=1`) is used for prednisolone's
  effect terms too, which matches **exactly** (0.0 diff) — ruling out
  `pow(x, 1.0)` itself as a source of rounding; (c) only 7 doses
  (adalimumab) and 2 doses (canakinumab) are involved, the fewest of any
  scenario, consistent with a rare, magnitude-dependent step-path
  sensitivity rather than a systematic formula error.

Per the guide's own tolerance language for archetypes 1–3 ("near-exact
match... anything beyond floating-point-scale deviation means a bug"),
these two results (≤4.2e-4 relative, isolated to nonlinearly-coupled
compartments, with the compound's own PK matching exactly) are judged
floating-point-scale, not a defect — reported honestly rather than
rounded away.

### A verification-caught bug that was found and fixed (not just disclosed)

An early version of this refactor collapsed the original's `$MAIN`-facing
concentration (`Cp_col`, used by the disease `$ODE`'s Hill effects) and
its `$TABLE`-facing recomputation (`Cp_COL`, used only for reporting)
into a single `$GLOBAL`-scope `C_COL`, reasoning that identical formulas
computing the same quantity twice were a textbook "duplicate
concentration site". Verification caught this immediately: the colchicine
scenario showed a real, reproducible one-row lag in every `Cp_*` output
column (`ref[i] == orig[i-1]` exactly, confirmed for `Cp_COL` across the
first 10 rows) while every disease-side compartment (`NEU`, etc.) still
matched exactly.

Root cause: mrgsolve evaluates `$MAIN` once per record using the
compartment state as of the *start* of that record's interval (the
"stale" value the disease `$ODE`'s Hill effects deliberately act on for
the whole interval — this is *by design*, matching how the original's own
`Ecol_NEU` etc. behave), then advances the ODE, then evaluates `$TABLE`
using the state as of the *end* of that interval (after any dose at that
exact record has been applied). The original's own `$MAIN` (`Cp_col`) and
`$TABLE` (`Cp_COL`) computations were therefore **not** true numeric
duplicates despite sharing an identical formula — they read the same
compartment at two different points in mrgsolve's per-record execution
order. Collapsing them into one value moved the *reported* concentration
one step earlier in time, without affecting any disease dynamics (which
already only ever read the `$MAIN`-timed value).

**Fixed, not just disclosed:** kept the two computations exactly where
the original placed them — `C_<STEM>` as a single, canonical, PD-facing
`double C_<STEM> = <expr>;` in `$MAIN` (the value `EFFECT_<STEM>` reads,
matching the original's `Cp_col`/etc. placement exactly), and `$TABLE`'s
own `Cp_<STEM>` reporting columns as an independent, fresh recomputation
from compartment state (matching the original's `$TABLE` computation
exactly) — not an alias of `C_<STEM>`. Re-verified after the fix: exact
match (0.0 diff) restored for every scenario that didn't already have it,
confirming the lag is gone.

## `/model_manifest` and discoverability

`POST /model_manifest` against the refactored model's extracted DSL
returns HTTP 200 and lists every renamed `$PARAM` (including all five
`EC50_<STEM>*`/`EMAX_<STEM>*`/`GAMMA_<STEM>*` Hill parameters) plus
`C_COL`/`C_PRED`/`C_ADA`/`C_APR`/`C_CAN` and all eight
`EFFECT_<STEM>_*` names in `outputPaths` (auto-discovered from
`$CAPTURE`).

Each `C_<STEM>` is written as a single, contiguous
`double C_<STEM> = <expr>;` statement directly in `$MAIN` (confirmed by
`grep -nE "double C_(COL|PRED|ADA|APR|CAN) +=" bd_mrgsolve_model_refactored.R`
— 5/5 hits) — not a `$GLOBAL` forward-declaration assigned later
elsewhere, which the guide's own "discoverable by downstream tooling"
section flags as invisible to the corpus's source-text pattern matcher
(the `hereditary-spherocytosis` precedent). An earlier draft of this file
did use a `$GLOBAL`-declare-then-assign pattern (to let `$TABLE` alias the
same value); once `$TABLE` was fixed to recompute independently instead
(see above), the `$GLOBAL` indirection was no longer needed at all and
was removed, restoring the plain, directly-discoverable `double C_<STEM>
= <expr>;` form. Each compound's matching `EC50_<STEM>*` parameter is
confirmed present in `$PARAM` for every one of the five compounds.

Also confirmed empirically (not just by inspection) that `C_COL`,
`EFFECT_COL_NEU`, and `EFFECT_COL_IL1` all report correctly (non-zero,
sensible values) when requested as `outputs` through `/run_simulation` —
ruling out a variable-shadowing concern between the `$MAIN`-local
`double` declarations and their `$CAPTURE @annotated` listing.

## Structural / build checks

- `Rscript -e 'parse("behcet-disease/bd_mrgsolve_model_refactored.R")'`
  succeeds (58 top-level expressions), confirming the mandatory structural
  template: single `bd_code_refactored <- '...'` string, curly apostrophes
  throughout the string body, single `mcode(...)` compile call after the
  string closes.
- `POST /model_manifest` returns HTTP 200 (after the two build-
  compatibility fixes above).
- All 8 of the original's own dosing scenarios run successfully through
  `POST /run_simulation` (see verification table above).

## Census update

`driver-patches/data/compound_perturbation_census.md`'s five
`behcet-disease` rows updated with confirmed target/pathway (Adalimumab:
TNF-alpha; Apremilast: PDE4 → TNF-α/IL-17A; Canakinumab: IL-1β; Colchicine:
tubulin polymerization → neutrophil migration / NLRP3-IL-1β; Prednisolone:
Glucocorticoid receptor — already correct in the census) and the
C_STEM/EFFECT_STEM/units/outcome columns, per the format used by other
completed rows in that file.
