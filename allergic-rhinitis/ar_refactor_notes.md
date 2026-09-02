# Allergic Rhinitis (AR) — PK/PD refactor notes

Refactor of `ar_mrgsolve_model.R` -> `ar_mrgsolve_model_refactored.R`, per
`FORK_WORKFLOW_GUIDE.md` Part 2. Original untouched.

## Census correction

`driver-patches/data/compound_perturbation_census.md` carries four rows for
this directory, all labeled by the automated classifier as
`Symptom model (<TOKEN>)` — a classifier artifact pulled from a nearby
comment, not a real compound name (exactly the failure mode the task
warned about). Checked against the actual code (`$PARAM`/`$INIT`/`$ODE`
names, `build_doses()`, the scenario list) and the directory's own
`README.md` — all four are genuine, named, externally dosed drugs:

| Census label | True compound | Route | Original PK shape |
|---|---|---|---|
| `Symptom model (CETI)` | **Cetirizine** — oral H1-antihistamine | oral, QD | depot -> central, linear elim (`CETI_D`, `CETI_C`) |
| `Symptom model (FP_LOCAL)` | **Fluticasone propionate** — intranasal corticosteroid | intranasal, local | single local-mucosal compartment, no depot (`FP_LOC`) |
| `Symptom model (MLKT)` | **Montelukast** — oral CysLT1 receptor antagonist | oral, QD | depot -> central, linear elim (`MLKT_D`, `MLKT_C`) |
| `Symptom model (OMA)` | **Omalizumab** — SC anti-IgE monoclonal antibody | SC, q4w | depot -> central -> peripheral + IgE-capture (`OMA_D`, `OMA_C`, `OMA_P`, `OMA_IGE`) |

Census rows updated accordingly (see bottom of this file).

## Build-compatibility fix (logged as `translations/UPSTREAM_ISSUES.md` #138)

The original does **not** compile under mrgsolve 2.0.1, for two independent,
syntax-only reasons — full detail in UPSTREAM_ISSUES.md #138:

1. `$INIT @annotated` uses `NAME = value : description` (an `=` sign) for
   all 22 entries; mrgsolve's annotated-block parser requires the colon
   form `NAME : value : description`, the same convention this file's own
   `$PARAM @annotated` block already uses correctly. Fixed in the
   refactored file: all 22 `$INIT` lines switched to the colon form, same
   names/values/descriptions.
2. `$TABLE` ends in five bare, space-separated `capture X Y Z ...` lines
   (no parentheses, no semicolon) — not valid C++, which is what `$TABLE`
   compiles as. One of those lines also lists seven names that are
   already-declared `$CMT` compartments (`IL4 IL5 IL13 HISTAMINE CYS_LT
   MAST_ACT EOS_N`), which mrgsolve 2.0.1 separately rejects if re-listed
   in `$CAPTURE`. Fixed in the refactored file: replaced with a single
   proper `$CAPTURE @annotated` block listing the 19 non-compartment names
   (compartments are exposed as output columns automatically, no
   re-listing needed).

Both fixes are syntax-only/non-numeric — confirmed by the verification
below, which compares the refactored file against a **patched-original
scratch copy** (original renaming untouched, carrying only these same two
fixes, never committed) and gets an exact match.

## Compound-by-compound archetype

- **Cetirizine (CETI)** — Archetype 3 minus peripheral (depot -> central,
  linear elimination). Renamed `CETI_D`/`CETI_C` -> `GUT_CETI`/`CENT_CETI`,
  `VD_CETI` -> `V1_CETI`, `H1_IC50_CETI` -> `EC50_CETI`, `H1_HILL` ->
  `GAMMA_CETI`. `EMAX_CETI = 1.0` is math-implied (the original's `H1_RO`
  ratio saturates at 1) — a rename/promotion, not a refit; the original's
  effect term was already the plain Hill-ratio shape.
  `EFFECT_CETI = EMAX_CETI * C_CETI^GAMMA_CETI / (EC50_CETI^GAMMA_CETI +
  C_CETI^GAMMA_CETI)` is computed once in `$TABLE`, matching where the
  original computed the analogous `h1ro` — the original's `$ODE`-block
  `h1ro` local was computed but never consumed by any `dxdt_` line (its
  only real consumer, `H1_EFF_HIST`, lived in `$TABLE` via a *second*,
  independently recomputed `H1_RO`), so single-sourcing it in `$TABLE`
  changes nothing behaviorally and removes a duplicate calculation.

- **Fluticasone propionate (FP)** — bespoke single compartment, no depot
  (`FP_LOC` -> `CENT_FP`), preserving the original's constant, dose-
  independent zero-order input term `KA_FP * FP_DOSE_BIOCONV / VD_FP_LOCAL`
  verbatim (a pre-existing behavioral quirk of the original — that term is
  a fixed number, not actually driven by any depot/dose state; disclosed
  here, not "fixed," since correcting it would be a real behavioral change
  out of scope for a structural refactor). Renamed `CL_FP_LOCAL` ->
  `CL_FP`, `VD_FP_LOCAL` -> `V1_FP`, `FP_DOSE_BIOCONV` -> `BIOCONV_FP`,
  `GR_IC50_FP` -> `EC50_FP`, `GR_HILL` -> `GAMMA_FP`. `EMAX_FP = 1.0`
  math-implied, same rename-not-refit reasoning as Cetirizine.
  `EFFECT_FP` (via `C_FP`) **must** stay computed where it feeds
  `dxdt_IL4`/`dxdt_IL5`/`dxdt_IL13`/`dxdt_EOS_N` live, i.e. in the ODE
  system every substep — matching "keep a calculation in the block the
  original used it in." See the dose-instant fix below for why these two
  are `$GLOBAL` macros rather than `$ODE`-scope doubles.

- **Montelukast (MLKT)** — Archetype 3 minus peripheral (depot -> central,
  linear elimination), directly analogous to Cetirizine. Renamed `MLKT_D`/
  `MLKT_C` -> `GUT_MLKT`/`CENT_MLKT`, `VD_MLKT` -> `V1_MLKT`,
  `CYSLTR1_IC50` -> `EC50_MLKT`. `GAMMA_MLKT = 1.0` (no explicit Hill
  coefficient in the original — its `CYSLTR1_INH` ratio was a plain
  `C/(IC50+C)` with implicit exponent 1), `EMAX_MLKT = 1.0` math-implied.
  Same "single-sourced in `$TABLE`, not consumed by any `dxdt_`" situation
  as Cetirizine's `EFFECT_CETI`.

- **Omalizumab (OMA)** — Archetype 3 (depot -> central -> peripheral) for
  its own PK (`OMA_D`/`OMA_C`/`OMA_P` -> `GUT_OMA`/`CENT_OMA`/`PERI_OMA`,
  `VC_OMA`/`VP_OMA` -> `V1_OMA`/`V2_OMA`, `KON_IGE`/`KOFF_IGE` ->
  `KON_OMA`/`KOFF_OMA`), **plus a bespoke TMDD-style IgE-capture term**
  (`OMA_IGE` -> `COMPLEX_OMA`) that is deliberately *not* forced into the
  standard Archetype-4 TMDD template: omalizumab binds the disease model's
  own shared `IGE_FREE` state (itself a dynamic disease compartment with
  its own synthesis/mast-binding/degradation terms), not an independent,
  compound-exclusive receptor pool with a fixed `RTOT`. There is no
  meaningful `RTOT_OMA` to declare. `EFFECT_OMA` is defined as the
  instantaneous fractional bound-IgE occupancy `COMPLEX_OMA / (IGE_FREE +
  COMPLEX_OMA + 1e-9)` — the closest honest analogue of the
  `COMPLEX_<STEM>/RTOT_<STEM>` pattern Archetype 4 uses, substituting the
  live free+bound total where a fixed `RTOT` would normally go. This is a
  disclosed bespoke structure, not a Hill-ratio "rename" (there was no
  original Hill term to rename — the original's actual PD action is the
  `ige_capture` ODE term itself, already fully expressed; `EFFECT_OMA` is
  an added descriptive readout, not a new gate consumed elsewhere).
  `C_OMA`/`ige_capture` stay in `$ODE`'s IgE-dynamics block, exactly where
  the original computed the analogous `oma_cp`/`ige_capture` (needed live,
  every substep, by `dxdt_IGE_FREE`).

Parameter *values* are unchanged from the original throughout — every rename
above is name-only; no PK parameter was invented, defaulted, or refit.

## Dose-instant reporting artifact — found and fixed (Fluticasone only)

Verification (below) initially showed `FP_LOCAL_NM`/`GR_OCC_FP` diverging
from the patched-original at the `t=0` duplicate report row of *every*
Fluticasone dose (abs diff up to 450 — the full dose amount — first-pass
result before the fix, not the final result), while every other output
across all seven scenarios matched exactly.

Root cause: Fluticasone doses directly into `CENT_FP` with no depot
compartment, and `C_FP`/`EFFECT_FP` were declared as `$ODE`-scope `double`
locals. mrgsolve/qspserver emits a duplicate report row at the exact
instant of a dose, and the `$ODE`-scope double is (re)computed from state
that has not yet had that instant's dose applied on that specific row —
the same quirk independently found and fixed for
`clostridioides-difficile-infection`'s bezlotoxumab (also dosed directly,
no depot; see that file's own notes and `UPSTREAM_ISSUES.md`). Cetirizine's
`EFFECT_CETI`, Montelukast's `EFFECT_MLKT`, and Omalizumab's `C_OMA` do
**not** hit this: the first two are computed only in `$TABLE` (not
`$ODE`), and Omalizumab is dosed into `GUT_OMA` (a *different* compartment
from the one `C_OMA` reads, `CENT_OMA`), so there is no discontinuity in
the read state at the dose instant either way.

**Fix applied** (matching the `clostridioides-difficile-infection`
precedent, preferred by the guide over disclosing-and-moving-on):
`C_FP`/`EFFECT_FP` are now `$GLOBAL` `#define` macros, not `$ODE`-scope
doubles:

```c
$GLOBAL
#define C_FP      (CENT_FP / V1_FP)
#define EFFECT_FP (EMAX_FP * pow(C_FP, GAMMA_FP) / (pow(EC50_FP, GAMMA_FP) + pow(C_FP, GAMMA_FP)))
```

A macro re-expands textually at every point of use (including inside
`$TABLE`, at `$TABLE`'s own execution time), so it always reads the
current, already-dosed state — no possibility of a stale cross-call read.
Confirmed via re-run through the live qspserver mrgsolve API (see the
verification table below): this **eliminates** the artifact entirely
(exact match, not just a narrower window).

**Disclosed discoverability trade-off** (same as
`clostridioides-difficile-infection`'s `C_VAN`/`C_BEZ`): a `#define` macro
is structurally incompatible with *also* writing a literal
`double C_FP = <expr>;` statement afterward (the preprocessor would
text-substitute the macro's own name inside that declaration, producing
invalid C++). So a naive `double C_<STEM> = ` text grep will find three of
this file's four compounds (`C_CETI`, `C_MLKT`, `C_OMA` — all genuine
`double ... = ...;` statements) but not `C_FP`. `EC50_FP` is still a
completely normal `$PARAM` entry, and `C_FP`/`EFFECT_FP` are both still
real output columns, confirmed present in `/model_manifest`'s
`outputPaths` (see below) — just not findable by that specific literal
grep pattern. Flagged here for whoever owns the downstream discovery
tooling to decide whether `#define C_<STEM> (...)` should also count as a
discoverable form (same open question cdi's notes already raised).

## Verification

**Method.** Built a patched-original scratch copy (`ar_mrgsolve_model.R`'s
own DSL block, unrenamed, carrying only the two UPSTREAM_ISSUES.md #138
syntax fixes — never committed, exists only as a local verification
artifact) and ran all seven of the original file's own `scenarios` (from
its own `build_doses()`/scenario list) through both it and
`ar_mrgsolve_model_refactored.R` via the live qspserver mrgsolve API
(`POST /model_manifest` then `POST /run_simulation`, `http://localhost:8007`,
sequential requests ~2 s apart). Same dosing (translated to the API's
`dosing`/`DoseSpec` list; compartment numbers unchanged 1:1 between the two
files, since compartment declaration order was preserved), same
`time = {end: 2016, delta: 1}` (84 days, hourly), identical for both
models. Compared every output name common to both files' `$CAPTURE` blocks
(19 unchanged-name reporting variables plus the 13 disease-side state
names shared by both models — 32 columns total) across the full time
grid.

**Result: exact match on every scenario, after the dose-instant fix
above.**

| Scenario | Dosing | n points | max abs diff | max rel diff |
|---|---|---|---|---|
| 1. Natural History (allergen only) | allergen pulse only | 2018 | `0.0` | `0.0` |
| 2. Cetirizine 10 mg QD | CETI | 2019 | `0.0` | `0.0` |
| 3. Fluticasone FP 200 μg/d | FP | 2019 | `0.0` | `0.0` |
| 4. Montelukast 10 mg QD | MLKT | 2019 | `0.0` | `0.0` |
| 5. Cetirizine + Fluticasone | CETI + FP | 2020 | `0.0` | `0.0` |
| 6. Omalizumab 300 mg q4w | OMA | 2019 | `0.0` | `0.0` |
| 7. Triple Therapy | CETI + FP + MLKT | 2021 | `0.0` | `0.0` |

(Point counts differ slightly across scenarios because mrgsolve/qspserver
emits an extra duplicate report row per dose event, and scenarios dose a
different number of times; both models produce the same count for the
same scenario, and every one of those rows matches exactly.)

Before the dose-instant fix, scenarios 3/5/7 (the only ones dosing
Fluticasone) showed a max abs diff of `450` (= the FP dose amount) on
`FP_LOCAL_NM`/`GR_OCC_FP` only, confined to the `t=0` duplicate report row
of each dose — every other output, every other timestep, and every other
scenario already matched exactly at that point. This confirmed the
mismatch was the known reporting artifact (see above), not a genuine
numeric bug in the renaming/restructuring itself, before the `$GLOBAL`
macro fix was applied and reverified to the exact-match table above.

**Structural checks:**
- `Rscript -e 'parse("ar_mrgsolve_model_refactored.R")'` succeeds (28 top-level
  expressions), confirming the file is valid R.
- Zero straight apostrophes (`'`) anywhere inside the quoted DSL string body
  (checked programmatically) — only curly apostrophes (`’`) are used.
- `ar_mrgsolve_model_refactored.cpp` (bare DSL extraction for qspserver
  `model_content`) is a byte-for-byte extraction of the quoted string body
  from `ar_mrgsolve_model_refactored.R` (`AR_model_refactored_code <- '...'`),
  confirmed via `/model_manifest` (200 OK, 66 parameters, 49 output paths)
  and via the `/run_simulation` comparisons above.
- `/model_manifest` confirms every compound's `EC50_<STEM>` is a real
  `$PARAM` entry (`EC50_CETI`, `EC50_FP`, `EC50_MLKT` all present) and that
  `C_FP`/`EFFECT_FP` (the macro-based pair) are present in `outputPaths`
  despite not being `double`-declared.

## Census update

`driver-patches/data/compound_perturbation_census.md` rows updated in place
(true compound identity, target/pathway, redirect site, outcome) — see
that file's `allergic-rhinitis` rows.
