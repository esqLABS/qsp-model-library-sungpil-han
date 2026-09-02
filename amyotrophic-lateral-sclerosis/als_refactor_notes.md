# Refactor notes — `amyotrophic-lateral-sclerosis/als_mrgsolve_model.R`

## Scope

Four assigned compounds, all confirmed present and correctly identified in
the actual code (the census census-row identities were accurate — no
misclassified/renamed drug found here):

| Compound | Stem | Census classification | Confirmed? |
|---|---|---|---|
| Riluzole | `RIL` | Redirect concentration (clean single site) | Yes |
| Edaravone | `EDA` | Redirect concentration (clean single site) | Yes |
| Phenylbutyrate (PB) | `PB` | Normalize duplicate concentration sites, then redirect | Yes — see below |
| Tofersen | `TOF` | Normalize duplicate concentration sites, then redirect | Yes — see below |

No fifth compound was found in the code. The README/PK table also names
**Masitinib** ("Phase 3 ongoing") but it has no PK/PD block in the actual
`als_mrgsolve_model.R` — matches the guide's own warning that the README is
unreliable for structural detail; nothing to refactor there since it isn't
modeled at all.

### Why PB and TOF genuinely needed "normalize, then redirect" (not just "redirect")

- **PB**: the original's own PK compartment is literally named `C_PB` (an mg
  **amount**), which collides with the naming convention's own token for the
  **exposed concentration** (`C_<STEM>`, μmol/L). This is a real naming
  duplicate that had to be resolved by renaming the compartment (`C_PB` →
  `CENT_PB`) before `C_PB` could be given its intended meaning. The original's
  own `$TABLE` re-derivation also used a *different* name for the same
  quantity (`Cp_PB` in `$MAIN` vs `Cp_PB_uM` in `$TABLE`) — a second
  divergent-naming signal for the same compound.
- **Tofersen**: the original genuinely has *two* candidate "concentration"
  quantities for the same compound — `C1_TOF` (plasma) and `C2_TOF` (CSF,
  read by the model's own `EC50_TOF` effect term via `Ccsf_TOF`). Only one of
  these is PD-facing; picking the correct one and demoting the other to an
  informational, non-exposed diagnostic is exactly the same
  "duplicate-site-normalize" pattern already used for
  `abdominal-aortic-aneurysm`'s doxycycline (aortic-tissue-vs-plasma).

By contrast Riluzole and Edaravone each have exactly one PK-derived
concentration site with no competing name, matching "clean single site."

## Archetypes

### Riluzole (RIL) — Archetype 3 (depot + central + peripheral)

Straight rename, no structural change:

| Original | Refactored |
|---|---|
| `DEPOT_RIL`/`C1_RIL`/`C2_RIL` (cmt) | `GUT_RIL`/`CENT_RIL`/`PERI_RIL` |
| `ka_RIL`/`CL_RIL`/`V1_RIL`/`Q_RIL`/`V2_RIL`/`F_RIL` | `KA_RIL`/`CL_RIL`/`V1_RIL`/`Q_RIL`/`V2_RIL`/`F_RIL` |
| `IC50_RIL`/`Emax_RIL` | `EC50_RIL`/`EMAX_RIL` (+ new `GAMMA_RIL=1`, none in original) |

`C_RIL = CENT_RIL/V1_RIL`. `EFFECT_RIL` is a rename of the original's
`E_RIL = Emax_RIL*Cp_RIL/(IC50_RIL+Cp_RIL+1e-9)` — already the Hill-interface
shape with an implicit gamma of 1, so this is a rename, not a fit. Feeds
`dxdt_Glut_syn` exactly as `E_RIL` did.

### Edaravone (EDA) — Archetype 1 (single compartment, linear)

| Original | Refactored |
|---|---|
| `IV_EDA` (cmt) | `CENT_EDA` |
| `CL_EDA`/`V_EDA` | `CL_EDA`/`V1_EDA` |
| `IC50_EDA`/`Emax_EDA` | `EC50_EDA`/`EMAX_EDA` (+ `GAMMA_EDA=1`) |

Dosed via mrgsolve's own `rate=` infusion field (60 mg over 1 h) — the
guide's "continuous infusion input" edge case does **not** apply here; this
is a normal event-based zero-order infusion, not a manually-added ODE input
term, so no special handling was needed. `C_EDA = CENT_EDA/V1_EDA`.
`EFFECT_EDA` renames `E_EDA`, feeding the `ROS_scav` term exactly as before.

### Tofersen (TOF) — Archetype 3 (depot + central(plasma) + peripheral(CSF))

| Original | Refactored |
|---|---|
| `DEPOT_TOF`/`C1_TOF`/`C2_TOF` (cmt) | `GUT_TOF`/`CENT_TOF`/`PERI_TOF` |
| `ka_TOF`/`CL_TOF`/`V1_TOF`/`Q_TOF`/`V2_TOF` | `KA_TOF`/`CL_TOF`/`V1_TOF`/`Q_TOF`/`V2_TOF` |
| `EC50_TOF`/`Emax_TOF` | `EC50_TOF`/`EMAX_TOF` (+ `GAMMA_TOF=1`) |

`C_TOF = PERI_TOF/V2_TOF` — the **CSF** concentration, matching what the
original's own `EC50_TOF` effect term reads (`Ccsf_TOF`). Plasma
(`Cp_TOF_plasma = CENT_TOF/V1_TOF`) is kept as an informational,
**non**-exposed diagnostic (not read by any PD term), same pattern as
`abdominal-aortic-aneurysm`'s `Cp_doxy_plasma`.

**Is Tofersen's PK/PD structure "bespoke"?** The task flagged this as a real
possibility since tofersen is an antisense oligonucleotide (ASO) acting via
direct SOD1 mRNA knockdown, which often has no meaningful EC50 in a
mechanistic sense. Checked directly against the code: **no bespoke handling
was needed.** This particular model represents tofersen's effect as a plain
phenomenological `Emax*C/(EC50+C)` ratio on CSF concentration
(`E_TOF = Emax_TOF*Ccsf_TOF/(EC50_TOF+Ccsf_TOF+1e-9)`), not as an ODE-solved
mechanistic gene-knockdown/translation system. The author chose a standard
empirical PK/PD simplification for this compound despite its real mechanism
being direct antisense knockdown — that's a modeling choice already baked
into the original, not something this refactor invented. Because the
original's own effect term is already exactly the Hill-interface shape, this
is a rename (gamma=1, no original Hill coefficient), not a mechanistic fit,
and Tofersen fits Archetype 3 cleanly — nothing about its structure required
bespoke treatment here. (If a different ALS model in this corpus represents
tofersen via an explicit ODE-solved translation/knockdown system instead,
*that* file would need the bespoke/fitted treatment — this one does not.)

### Phenylbutyrate / AMX0035 (PB) — Archetype 3 minus peripheral (depot + central)

| Original | Refactored |
|---|---|
| `DEPOT_PB`/`C_PB` (cmt, mg amount) | `GUT_PB`/`CENT_PB` |
| `ka_PB`/`CL_PB`/`V_PB`/`F_PB` | `KA_PB`/`CL_PB`/`V1_PB`/`F_PB` |
| `EC50_PB`/`Emax_PB`/`Emax_mito_PB` | `EC50_PB`/`EMAX_PB`/`EMAX_PB_MITO` (+ `GAMMA_PB=1`, shared by both effects, matching the original's single shared `EC50_PB`) |

`C_PB = (CENT_PB/V1_PB)*1000.0` (μmol/L) — a rename of the original's
`Cp_PB`, now finally able to use the name `C_PB` because the compartment
that used to occupy that name was renamed to `CENT_PB` (see "duplicate
site" discussion above).

**Two effects, one shared concentration, matching the original exactly:**
`EFFECT_PB` (ER stress/CHOP reduction) and `EFFECT_PB_MITO` (mitochondrial
protection) both read `C_PB`, both use `GAMMA_PB=1`, and — matching the
original precisely — both reuse the **same single** `EC50_PB` (the original
never gave the ER-stress and mitochondrial effects separate potencies; only
`Emax` differs between them). Each has its own `EMAX_<STEM>`. This mirrors
the doxycycline MMP-9/MMP-2 precedent (`abdominal-aortic-aneurysm`) of one
concentration driving two named effects.

**Pre-existing dead-variable defect, disclosed, not "fixed":** the
original's `E_PB` (`Emax_PB*Cp_PB/(EC50_PB+Cp_PB+1e-9)`, the ER-stress/CHOP
effect) is computed in `$MAIN` but **never read by any `$ODE` line** —
grepped the original file directly: `E_PB` appears exactly once, at its own
definition. Only `E_mito_PB` is actually wired into `dxdt_Mito`. This is a
genuine incompleteness in the original PD model (an effect term computed
and then silently dropped), not introduced by this refactor. `EFFECT_PB` is
kept in the refactored file with the same "computed but not read by any
$ODE" status — same treatment this corpus already gives other
found-dead/unused parameters (e.g. `EMAX_BB_HR`/`EC50_BB_HR` in
`abdominal-aortic-aneurysm`) — disclosed here rather than silently wired up
(which would be a behavioral change, not a rename).

## Two pre-existing mrgsolve 2.0.1 build defects, fixed syntax-only in the sibling only

`POST /model_manifest` on the untouched original's own extracted DSL
(`als_code <- '...'`) fails before any compound-specific code is reached.
Logged as **`translations/UPSTREAM_ISSUES.md` entry — see the file's own
numbering; both are new entries added by this session**:

1. **`$CMT` (bare) + `$INIT` (`name = value`) jointly declare every
   compartment.** mrgsolve 2.0.1's `validObject()` rejects this outright:
   `invalid class "mrgmod" object: Duplicated model names: DEPOT_RIL C1_RIL
   ... FVC` (all 26 compartments). Same defect class as `sepsis`'s issue #30
   point 2 in `UPSTREAM_ISSUES.md`.
2. **Eleven `$PARAM` names collide with mrgsolve's own auto-generated
   `<compartment>_0` initial-value symbol**, once (1) is worked around by
   giving `$CMT` its own 3-field annotated value: `MN_upper_0`,
   `MN_lower_0`, `TDP43_nuc_0`, `Ca_i_0`, `ROS_0`, `GSH_0`, `Mito_0`,
   `BDNF_0`, `NfL_CSF_0`, `ALSFRS_0`, `FVC_0` are declared in `$PARAM` as
   plain baseline-reference constants (used inside `MN_death_rate`,
   `Mic_stim`, `EAAT2_eff`, `MN_frac`, etc.), and every one of them has the
   exact same name as the compartment it references (`MN_upper`, `ROS`,
   `BDNF`, ...). mrgsolve auto-generates a `const double& <CMT>_0` reference
   for every compartment that has an explicit non-default initial value;
   with the identically-named user parameter also present, gcc reports
   `conflicting declaration`/`redeclaration` errors for all eleven at the
   C++ compile stage. Same defect class as `sepsis`'s issue #30 point 3.

**Syntax-only fix applied in this sibling only** (never in the checked-in
original):

- The eleven colliding baseline-reference parameters are renamed
  `<name>_0` → `<name>_REF` (`MN_upper_REF`, `MN_lower_REF`,
  `TDP43_nuc_REF`, `Ca_i_REF`, `ROS_REF`, `GSH_REF`, `Mito_REF`,
  `BDNF_REF`, `NfL_CSF_REF`, `ALSFRS_REF`, `FVC_REF`), every usage site
  updated in lockstep. `Mic_0` did not collide (the compartment is
  `Mic_act`, not `Mic`) and is unchanged.
- `$CMT` is declared bare (2-field: name + description only, no inline
  value) and the separate `$INIT` block is removed. Nonzero initial
  conditions are instead set via the modern `<CMT>_0 = value;` idiom inside
  `$MAIN` (now free to use, since the colliding params were renamed away).
  **This was not an arbitrary choice** — see the next section for why.
- No numeric parameter value or model behavior changes; both fixes are
  purely structural/naming. Confirmed non-numeric by the verification
  results below (exact match, `0.0` on every compared output).

## A third defect found only once the first two were fixed: `$CAPTURE` cannot recompute a name already hoisted from `$MAIN`

Once defects 1–2 above were fixed (necessary just to get `/model_manifest`
to return successfully at all), a **third**, previously-masked compile
error surfaced: the original recomputes `Cp_RIL`/`Cp_EDA`/`Ccsf_TOF` in
`$TABLE` under the **same name** already declared as a `double` in `$MAIN`
(`capture Cp_RIL = C1_RIL/V1_RIL;` vs `double Cp_RIL = C1_RIL/V1_RIL;`).
mrgsolve 2.0.1 rejects this: `redefinition of 'capture {anonymous}::Cp_RIL'
... previously declared here: 'double {anonymous}::Cp_RIL'`. The original's
own `Cp_PB`/`Cp_PB_uM` pair and `E_RIL`/`E_RIL_cap`,
`E_EDA`/`E_EDA_cap`, `E_TOF`/`E_TOF_cap` pairs all *already* use a
different name for the `$TABLE` recompute than the `$MAIN` double, and
those four do **not** hit this — only `Cp_RIL`/`Cp_EDA`/`Ccsf_TOF` reuse the
exact same name in both blocks. Same defect class as
`UPSTREAM_ISSUES.md`'s "`$CAPTURE` duplicates compartment names" entries,
just against a `$MAIN`-declared double instead of a `$CMT` compartment.

This is disclosed here as **found only during verification, and is why the
in-memory-only original-side scratch copy used for verification (never
saved to the checked-in original) needed the same class of fix already
applied to the deliverable** — the *real* discovery this defect forced,
described next, is the one that actually mattered for correctness.

## The real finding: a genuine reporting-timing bug, not just a compile workaround, and why $TABLE reassignment (not $CAPTURE) is the fix

The first attempt to resolve the defect above moved
`Cp_RIL`/`Cp_EDA`/`Ccsf_TOF` (and, in the refactored sibling,
`C_RIL`/`C_EDA`/`C_TOF`/`C_PB`) out of `$TABLE` entirely and exposed the
`$MAIN` doubles directly via a bare `$CAPTURE` block. This **compiles** and
even looked correct in a no-dosing spot-check — but running the `amx0035`
scenario's own dosing through the qspserver `mrgsolve_api` exposed a real
bug: the reported `C_PB` value was **permanently offset by one reporting
row** against its own `CENT_PB/V1_PB` (e.g. at the first `t=12` row,
`CENT_PB=178.9278` but the naively-`$CAPTURE`d `C_PB` still read `0`,
matching what `CENT_PB` had been at the *previous* row, not the current
one). Isolated with a minimal single-compartment reproduction
(`$MAIN: double C_X = CENT_X/V1_X;` `$CAPTURE C_X`) dosed repeatedly: `C_X`
was offset by exactly one row from `CENT_X/V1_X` at **every** row, not just
at dose instants — confirming this is the guide's own documented `$MAIN`
semantics ("`$MAIN` evaluates once per interval using state from the
*start* of that interval") manifesting as a **permanent** one-step lag
whenever a `$MAIN`-computed double is read directly via `$CAPTURE`, not the
milder, self-healing "dose-instant reporting artifact" described elsewhere
in the guide for `$ODE` locals.

**This explains why the original itself recomputes these quantities in
`$TABLE`** (`Cp_PB_uM`, `E_RIL_cap`, `E_EDA_cap`, `E_TOF_cap`) rather than
just exposing its `$MAIN` locals — `$TABLE` runs after `$ODE` has
integrated up to the row actually being reported, so a value **recomputed**
there is contemporaneous. The fix used in the refactored deliverable: keep
`double C_<STEM> = ...`/`double EFFECT_<STEM> = ...` in `$MAIN` (needed
there so the disease ODEs read them with the same timing the original
used, and needed for the corpus-wide `double C_<STEM> = <expr>;`
discoverability grep), then in `$TABLE` **reassign** the same names with a
bare `NAME = expr;` statement — no `double`/`capture` keyword, so it is a
reassignment of the symbol mrgsolve already hoisted from `$MAIN`, not a
redeclaration, and compiles cleanly while reading fresh, correct,
per-row state. Verified directly: the minimal reproduction model with this
fix applied read `C_X` correctly at every row, and the full ALS model's
`amx0035` scenario changed from `max_abs_diff=3578.6` (the naive-`$CAPTURE`
version) to `0.0` (the `$TABLE`-reassignment version) against the
identically-fixed original.

**Disclosure:** the untouched original itself is unaffected by this bug for
`Cp_RIL`/`Cp_EDA`/`Ccsf_TOF` in normal use (it already recomputes them
under matching names in `$TABLE`) — the bug was only ever a risk in the
in-memory verification scratch copy and in an early draft of this
refactor's own deliverable, both fixed before delivery. Nothing about this
finding implicates the checked-in original.

## Verification

**Method.** Extracted the bare DSL text from
`als_mrgsolve_model_refactored.R` (`als_refactored_code <- '...'`, saved
verbatim as `als_mrgsolve_model_refactored.cpp`, confirmed byte-identical)
and, since the untouched original does not compile (three defects above),
an **in-memory-only** scratch copy of the original with the same class of
syntax-only fixes was used for the original side of the comparison — the
checked-in original was never modified. Both DSLs run through the local
qspserver `mrgsolve_api` (`http://localhost:8007`), `POST /model_manifest`
then `POST /run_simulation`, requests spaced ~2 s apart, one at a time (no
concurrency).

**Dosing.** Built from the original file's own `make_dosing()` R function,
translated into the API's `dosing` field (numeric 1-based `cmt`, from each
side's own `/model_manifest` `compartments` list — both sides share
identical compartment ordering since the refactor is a pure rename, so the
same indices apply to both) — riluzole 50 mg q12h (1096 doses over 548 d),
edaravone 60 mg IV over 1 h per the original's 14-on/14-off cycle schedule
(114 doses), tofersen 100 mg SC loading (d0/14/28) + q28d maintenance (21
doses), AMX0035/PB 3000 mg q12h (1096 doses). `combo_ril_eda` and
`all_drugs` required merging and **sorting by time** — the API's runner
requires a monotonically-sorted dosing table (`the data set is not sorted
by time` otherwise); this is a harness-construction detail, not a modeling
issue, and was applied identically to both sides.

**All 7 of the original's own scenarios run, exact match (`max_abs_diff =
0.0`) across every shared/renamed `$CAPTURE`d output, full time grid:**

| Scenario | Points compared | Window | Max abs diff |
|---|---|---|---|
| `no_treatment` | 549 | 548 d | 0.0 |
| `riluzole` | 1645 | 548 d | 0.0 |
| `edaravone` | 663 | 548 d | 0.0 |
| `combo_ril_eda` | 1759 | 548 d | 0.0 |
| `tofersen` | 17 | **shortened to 14 d** (see below) | 0.0 |
| `amx0035` | 1645 | 548 d | 0.0 |
| `all_drugs` | 2876 | 548 d | 0.0 |

38 output columns compared per scenario (26 renamed/unchanged compartments,
4 renamed concentration pairs, 3 renamed effect pairs, 5 unchanged
`$TABLE`-derived disease outputs). Every one of them is bit-identical
between the original (fixed only for the three build defects above) and
the refactored sibling — the expected result for archetypes 1 and 3 with a
pure rename Hill interface (no fitting was needed or performed for any of
the four compounds).

**`tofersen` scenario window — a genuine, pre-existing numerical-stability
defect, not a solver-step-count budget limit and not a refactor bug.** At
full 548-day duration with `f_SOD1_mut=1.0` (the SOD1-ALS subtype
parameterization the original's own `run_scenario()` selects specifically
for this scenario), lsoda fails identically on **both** the original and
the refactored sibling at `t≈384 h` (`repeated error test failures`, then
`negative istate: -4`). Traced the mechanism directly: `dxdt_Mito` has no
floor, and `ROS_prod = k_ROS_prod*(1+SOD1_burden)*Ca_i/Mito` divides by it;
with `f_SOD1_mut=1` the SOD1-misfolding/aggregation cascade drives `Ca_i`
and `ROS` up quickly enough that `GSH` and `Mito` are both driven toward
zero (`Mito=1.37e-6`, `GSH=2.29e-6` by `t=336 h` in a diagnostic run), at
which point the `Ca_i/Mito` division creates a positive-feedback numerical
singularity lsoda cannot integrate through. This reproduces **identically**
on both models (same failure time, same lsoda diagnostic), confirming it
is a pre-existing defect in the original's own ODE system (an unfloored
divisor combined with an aggressive familial-ALS parameterization), not
something this refactor introduced. Per the fork guide's precedent for
solver-budget limits, the verification window for this scenario alone was
shortened to 336 h (14 days, well before the blow-up) — same dosing, same
parameterization — rather than declaring a mismatch; this is disclosed
here rather than silently substituting a different scenario. Logged as a
new `translations/UPSTREAM_ISSUES.md` entry (numerical-stability class, not
the usual build-defect class, but the same "log what you find, don't fix
it upstream" rule).

**Hill-fit quality:** not applicable — all four compounds (including
Tofersen, see the archetype discussion above) are pure renames of an
already-Hill-shaped ratio in the original; no `nls()` fit was performed or
needed.

## Census

`driver-patches/data/compound_perturbation_census.md` rows for
`amyotrophic-lateral-sclerosis` updated: `Edaravone`/`Riluzole` confirmed
"clean single site"; `Phenylbutyrate (PB)`/`Tofersen` confirmed "duplicate
sites, normalize" for the reasons given above; target/pathway and the
`C_<STEM>`/`EFFECT_<STEM>` mapping filled in for all four; outcome column
records this refactor and its exact-match verification result.
