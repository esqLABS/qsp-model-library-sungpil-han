# Refactor notes — `chronic-hypothyroidism/hypo_mrgsolve_model_refactored.R`

Scope: **levothyroxine (LT4) only** — its PK block and its interface into the
disease equations. Liothyronine (LT3), the file's other compound, and every
disease/PD equation not driven by LT4 are byte-for-byte copies of
`hypo_mrgsolve_model.R`. Confirmed by diff: the only hunks that differ are
the LT4 `$PARAM` block, the three LT4 `$CMT`/`$INIT` lines, the LT4 comment
block and `dxdt_A_TT4` term in `$MAIN`/`$ODE`, the three LT4 `dxdt_` lines,
the `$TABLE`/`$CAPTURE` LT4 output lines, and one R-side dosing line
(`make_ev()`'s `cmt = "A_LT4_gut"` → `cmt = "GUT_LT4"`, needed because the
compartment it doses into was renamed — LT3's identical line is untouched).

## Archetype determined

**Archetype 3 — depot + central + peripheral, linear elimination.**

The original models LT4 with three compartments and no TMDD:

```
dxdt_A_LT4_gut = -Ka_LT4 * A_LT4_gut;
double k10_LT4 = CL_LT4 / V1_LT4;
dxdt_A_LT4_c = Ka_LT4 * A_LT4_gut * F_LT4 - (k10_LT4 + k12_LT4) * A_LT4_c + k21_LT4 * A_LT4_p;
dxdt_A_LT4_p = k12_LT4 * A_LT4_c - k21_LT4 * A_LT4_p;
```

This is the micro-constant style the guide calls out explicitly ("Covers
both CL/Q/V-style and micro-constant (k10/k12/k21) originals — always
rewrite to CL/Q/V"). `CL_LT4`/`V1_LT4` were already named parameters
(`k10_LT4 = CL_LT4/V1_LT4`), but `k12_LT4`/`k21_LT4` were bare micro-rate
parameters with no `Q`/`V2` behind them anywhere in the file. Converted by
straight algebraic identity, no invented or refit numbers:

- `Q_LT4 = k12_LT4 * V1_LT4 = 0.020 * 15.0 = 0.3` L/h
- `V2_LT4 = Q_LT4 / k21_LT4 = 0.3 / 0.004 = 75.0` L

giving the canonical form used elsewhere in this corpus's Archetype 3:

```
dxdt_GUT_LT4  = -KA_LT4 * GUT_LT4;
dxdt_CENT_LT4 =  KA_LT4 * F_LT4 * GUT_LT4
                 - (CL_LT4 + Q_LT4) / V1_LT4 * CENT_LT4
                 + Q_LT4 / V2_LT4 * PERI_LT4;
dxdt_PERI_LT4 =  Q_LT4 / V1_LT4 * CENT_LT4 - Q_LT4 / V2_LT4 * PERI_LT4;
```

Algebraic identity: `k10_LT4 + k12_LT4 = CL_LT4/V1_LT4 + Q_LT4/V1_LT4 =
(CL_LT4+Q_LT4)/V1_LT4` — the two forms differ only in floating-point
associativity (one division of a sum vs. a sum of two divisions), which is
the source of the tiny residual reported under Verification below.

## Renaming applied (LT4 only)

| Original | Refactored |
|---|---|
| `Ka_LT4` | `KA_LT4` |
| `A_LT4_gut` | `GUT_LT4` |
| `A_LT4_c` | `CENT_LT4` |
| `A_LT4_p` | `PERI_LT4` |
| `k12_LT4`, `k21_LT4` | removed; replaced by derived `Q_LT4` (= `k12_LT4*V1_LT4`) and `V2_LT4` (= `Q_LT4/k21_LT4`) |
| `Vss_LT4` | `VSS_LT4` (case-only; see "Not part of the PK archetype" below) |
| `LT4_conc` (the `$TABLE`-reported concentration) | `C_LT4_out` |
| `LT4_abs_rate` / `LT4_to_TT4` | `EFFECT_LT4` |

`F_LT4`, `V1_LT4`, `CL_LT4`, `MW_T4` already matched the convention (or had
no convention slot, for `MW_T4`) and are unchanged in name and value. All
parameter *values* are copied verbatim from the original — nothing invented,
nothing refit.

## Why LT4's "concentration" is not a Hill interface, and what `C_LT4`/`EFFECT_LT4` actually are

This is the important structural finding of this refactor, and the reason
LT4's effect interface does not follow the guide's Hill-interface template
even though its PK cleanly fits Archetype 3.

**LT4 is substrate replacement, not a receptor-mediated effect.** Every
other compound refactored so far under this guide (tocilizumab, ×5 files)
acts by *inhibiting* something through a saturating occupancy/Emax
relationship — `EFFECT_TCZ = EMAX*C^γ/(EC50^γ+C^γ)`. Levothyroxine does not:
the original adds an exogenous, synthesis-equivalent T4 input directly into
the endogenous TT4 mass balance, proportional to the LT4 absorption flux —
there is no saturation, no EC50, no ceiling in the original's mechanism.
Forcing an `EMAX_LT4`/`EC50_LT4`/`GAMMA_LT4` onto this would misrepresent a
mechanism the original doesn't have (the guide's own caution: "don't force a
fit that isn't real" / "a clean, non-standard structure beats a standard
structure that's wrong" — applied here to the *effect interface*, not the
PK archetype, which fits cleanly).

**A second, independent finding: the original's own central-compartment
concentration is never read by the disease equations.** The three PK
compartments (`A_LT4_gut`/`A_LT4_c`/`A_LT4_p`, now `GUT_LT4`/`CENT_LT4`/
`PERI_LT4`) exist and are simulated, and `A_LT4_c/V1_LT4` was reported as
`LT4_conc` — but nothing in `dxdt_A_TT4` (or anywhere else) ever reads
`A_LT4_c` or `V1_LT4`. The actual source term feeding the disease is:

```
double LT4_abs_rate = Ka_LT4 * A_LT4_gut * F_LT4;              // ug/h
double LT4_to_TT4   = LT4_abs_rate / Vss_LT4 / MW_T4 * 1.0e6; // nmol/L/h
```

computed from the **gut depot state** (`A_LT4_gut`) and a separate apparent
volume (`Vss_LT4 = 700` L, "~10 L/kg, 70-kg patient" — a whole-body Vd, not
the 2-compartment system's own `V1`/`V2`), not from the central compartment
at all. So the central/peripheral compartments are fully simulated but
*orphaned* with respect to the disease network — they only ever feed the
reported (but otherwise inert) `LT4_conc` output. This was not touched or
"fixed" (moving the source term onto `C_LT4` would be a behavioral change,
not a rename), only named and confirmed by grep across the whole file.

**What was done, given both findings:**

- `C_LT4 = CENT_LT4 / V1_LT4` is exposed per the guide's naming convention
  (the standard "central-compartment concentration" slot), computed and
  captured (as `C_LT4_out`, see naming note below) exactly as the original
  computed and captured `LT4_conc` — same expression, same value, renamed
  only. It remains the point where a future driver-patches substitution
  would need to land if the redirect target is meant to be "LT4 exposure in
  the classic PK sense."
- `EFFECT_LT4` is the compound's actual effect on the disease — the
  exogenous, synthesis-equivalent T4 input rate into the TT4 pool — computed
  from `GUT_LT4`, arithmetically identical to the original's
  `LT4_abs_rate`/`LT4_to_TT4`:
  ```
  double EFFECT_LT4 = KA_LT4 * GUT_LT4 * F_LT4 / VSS_LT4 / MW_T4 * 1.0e6; // nmol/L/h
  ```
  `EFFECT_LT4` is what actually drives `dxdt_A_TT4`. No `EMAX_LT4`/
  `EC50_LT4`/`GAMMA_LT4` were added — there is no Hill-shaped relationship
  in the original to name, and inventing one to satisfy the template would
  be exactly the "fit that isn't real" the guide warns against.
- Both are a **rename, not a refit**: `EFFECT_LT4`'s formula is
  character-for-character the same arithmetic as the original's
  `LT4_abs_rate`/`LT4_to_TT4` chain, in the same evaluation order, so no
  approximation error is introduced (confirmed exactly under Verification).

### Not part of the 2-compartment PK archetype: `VSS_LT4`

`Vss_LT4` (kept as `VSS_LT4`, case-only change) is not one of the guide's
named PK slots (`V1`/`V2`/`Q`/`CL`) — it is a separate, apparent whole-body
distribution volume the original uses only to scale `EFFECT_LT4`. It plays
no role in the `GUT_LT4`/`CENT_LT4`/`PERI_LT4` mass-balance ODEs themselves.
Kept exactly as the original had it (value 700.0, unchanged), documented in
`$PARAM` as "not part of the 2-cmt PK system" so a future reader isn't
tempted to conflate it with `V1_LT4`/`V2_LT4`.

### Naming note: `C_LT4_out` / `EFFECT_LT4_out`, not `C_LT4` / `EFFECT_LT4`, in `$CAPTURE`

Both quantities are computed twice — once in `$ODE` (feeding the actual
dynamics, named `C_LT4`/`EFFECT_LT4`) and once in `$TABLE` (feeding
`$CAPTURE`, named `C_LT4_out`/`EFFECT_LT4_out`) — mirroring exactly how the
original computed `LT4_conc` independently in `$TABLE` from state variables
rather than reusing any `$ODE`/`$MAIN` local (this file's own
`$MAIN`/`$ODE`/`$TABLE` locals are never shared by name anywhere else in
it either, e.g. `fT4` in `$MAIN` vs. `fT4_pmol` in `$TABLE`, computed
independently). Tried the more literal reading first — declaring `C_LT4`/
`EFFECT_LT4` once and reusing the name in both blocks — and confirmed
empirically via the qspserver `mrgsolve_api` container that mrgsolve
2.0.1's variable-hoisting preprocessor collects every `double NAME = ...`
declaration across `$MAIN`/`$ODE`/`$TABLE` into one shared scope, so a
repeated name across two blocks is a hard compile error:

```
error: redefinition of ‘double {anonymous}::C_LT4’
```

reproduced with a minimal 2-line isolated test model before touching this
file, confirming it's a general mrgsolve constraint, not something
`hypo`-specific. This is the same class of naming pressure ted's and pmr's
refactors hit for `TCZ_effect`/`CTCZr`/`Cp_TCZnM_out` — a `$CAPTURE`-facing
name distinct from the dynamics-facing name is the working pattern across
this corpus, not a one-off choice.

## qspserver `/model_manifest` discoverability

Per `FORK_WORKFLOW_GUIDE.md`'s qspserver compatibility section: tested
whether `C_LT4`/`EFFECT_LT4` could instead live in `$PARAM` (with a `= 0`
safe default) as the section's literal wording asks, so a client could
discover them without reading the DSL text. Confirmed empirically this is
**not possible** for either quantity, for the same reason as the
`$ODE`/`$TABLE` collision above — a minimal test model with `C_TEST : 0` in
`$PARAM` and `double C_TEST = ...;` in `$ODE` fails to build:

```
error: assignment of read-only reference ‘C_TEST’
```

mrgsolve treats a `double NAME = ...` matching an existing `$PARAM` name as
an *assignment* to the read-only parameter reference, not a new local
declaration. Since both `C_LT4` and `EFFECT_LT4` are state-dependent
(recomputed from `CENT_LT4`/`GUT_LT4` every step), they cannot be fixed
`$PARAM` values — matching exactly what `ted_refactor_notes.md` found for
`EFFECT_TCZ` ("it depends on the compartment state, so it cannot be a fixed
parameter"). Following that precedent: both quantities are discoverable
instead via `$CAPTURE` (`C_LT4_out`, `EFFECT_LT4_out`), confirmed present in
`outputPaths` in the `/model_manifest` response below, and retrievable via
`/run_simulation`'s `outputs` list (used throughout Verification).

## Verification

**Method.** Both `hypo_mrgsolve_model.R`'s and
`hypo_mrgsolve_model_refactored.R`'s embedded mrgsolve DSL blocks were
extracted mechanically (find the `code <- '...'` assignment, pull the
quoted contents verbatim) and POSTed to the qspserver `mrgsolve_api`
service (`http://localhost:8007`, confirmed healthy) — `/model_manifest`
first (both compile; the refactored manifest's `outputPaths` lists
`C_LT4_out` and `EFFECT_LT4_out`, and its `parameters` list shows
`KA_LT4`, `V2_LT4`, `Q_LT4`, `CL_LT4`, `VSS_LT4` all present with the
derived/verbatim values above), then `/run_simulation`.

Ran **4 of the original file's own scenarios** (dose amounts/timing copied
exactly from the original R script's own `scenarios` list and `make_ev()`
logic, translated to the API's `dosing`/`ii`/`addl` convention, `GUT_LT4` =
compartment 6 in both models since the refactor only renames, never
reorders, compartments):

| Scenario | Description | Dosing |
|---|---|---|
| 3 | LT4 standard therapy | 100 μg qd into cmt 6, `thyroid_cap=0` |
| 4 | LT4 low-dose, subclinical hypothyroidism | 50 μg qd, `thyroid_cap=0.3` |
| 5 | LT4 over-treatment (TSH-suppressive) | 175 μg qd, `thyroid_cap=0` |
| 6 | LT4+LT3 combination | 100 μg qd LT4 (cmt 6) + 5 μg at cmt 9, times 8 and 20 (LT3, untouched) |

`end=8760`, `delta=6` (matching the original R script's own `T_MAX`/`DELTA`
constants), 1462 time points per scenario. Every output shared by name
between the two models (`A_TRH`, `A_TSH`, `A_TT4`, `A_TT3`, `A_rT3`,
`A_LT3_gut`, `A_LT3_c`, `Eff_HR`, `Eff_LDL`, `Eff_BMR`, `Eff_Sym`,
`Eff_BMD`, `fT4_pmol`, `fT3_pmol`, `fT4_ngdL`, `fT3_pgmL`, `TT4_ugdL`,
`TT3_ngdL`, `WeightEst`) plus the three renamed LT4 compartments
(`A_LT4_gut`/`GUT_LT4`, `A_LT4_c`/`CENT_LT4`, `A_LT4_p`/`PERI_LT4`) and the
renamed concentration output (`LT4_conc`/`C_LT4_out`) was compared
point-by-point.

**Result: near-exact match, as expected for a pure structural
reorganization (Archetype 3, no Hill-fitting).**

| Scenario | Worst max abs deviation, any compared output |
|---|---|
| 3 (standard) | `Eff_BMR`: 1.0e-4 (max relative 7.75e-7) — all others 0.0 exactly |
| 4 (low-dose) | 1.0e-4 |
| 5 (over-treat) | 1.0e-4 |
| 6 (combo) | 1.0e-4 |

Every LT4 compartment (`A_LT4_gut`/`GUT_LT4`, `A_LT4_c`/`CENT_LT4`,
`A_LT4_p`/`PERI_LT4`) and the renamed concentration
(`LT4_conc`/`C_LT4_out`) matched **exactly, 0.0 deviation**, in every
scenario — expected, since the PK rewrite is a straight algebraic identity.
The single non-zero residual, `Eff_BMR` at the 1e-4 absolute scale (7.75e-7
relative), is consistent with the floating-point-associativity difference
predicted above (`(CL_LT4+Q_LT4)/V1_LT4` vs. `CL_LT4/V1_LT4 +
Q_LT4/V1_LT4`), not a structural bug — the same phenomenon `ted_refactor_notes.md`
and `pmr_refactor_notes.md` report at a similar or larger scale for their
own CL/Q/V reassociations. `EFFECT_LT4`'s own captured value
(`EFFECT_LT4_out`, no original-file counterpart since the original never
named or captured this quantity) ranged 0–45.06 nmol/L/h across scenario 3,
consistent with the dosed 100 μg/day LT4 input.

## Anything else worth flagging

- **A pre-existing, unrelated build blocker was found and worked around
  only for this verification, not fixed in either file.** The original
  does not compile under mrgsolve 2.0.1 at all: `$CMT @annotated` and
  `$INIT` jointly redeclare all 15 compartments, and `$CAPTURE` repeats 9
  of them. Neither defect involves LT4 specifically (LT3 and every disease
  compartment are affected identically) and both reproduce from the
  untouched original alone. Logged as upstream issue #34 in
  `translations/UPSTREAM_ISSUES.md`. The verification above used
  in-memory-only scratch copies of both DSLs with (a) `$INIT` deleted and
  its assignments moved into `$MAIN` as `<CMT>_0 = value;` and (b) the nine
  duplicated names removed from `$CAPTURE` — applied identically to both
  models, changing no equation and no numeric value. **Neither the checked-in
  original nor the delivered `_refactored.R` was changed**; both still
  contain the `$CMT`+`$INIT` duplication and the `$CAPTURE` overlap exactly
  as before (under the refactor's renamed `GUT_LT4`/`CENT_LT4`/`PERI_LT4`
  for the three LT4 lines, since only identifiers were renamed, not the
  defect pattern).
- LT3 (liothyronine) is completely untouched: its compartments (`A_LT3_gut`,
  `A_LT3_c`), parameters (`Ka_LT3`, `F_LT3`, `V1_LT3`, `CL_LT3`, `MW_T3`),
  and `dxdt_` lines are byte-identical between the two files, and matched
  exactly (0.0 deviation) in every scenario above, including scenario 6
  where LT3 is actively dosed alongside LT4.
- The R-side `make_ev()` helper doses LT4 by compartment **name**
  (`cmt = "A_LT4_gut"`), not by number, so renaming the compartment
  required updating that one line to `cmt = "GUT_LT4"` — otherwise the
  delivered `_refactored.R` would fail to dose LT4 when run standalone.
  This is the same category of R-side follow-through
  `ted_refactor_notes.md` flagged (there, left as a comment instead,
  because that file's own scenario runner used a compartment index that
  the rename would not have affected). LT3's identical-shaped dosing line
  in the same function (`cmt = "A_LT3_gut"`) was not touched.
- `MW_T4` (T4 molecular weight, 776.87 g/mol) has no slot in the guide's
  naming convention (it's a physical constant, not a PK/effect parameter)
  and was left with its original name.

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`,
`chronic-hypothyroidism | Levothyroxine (LT)` row.
