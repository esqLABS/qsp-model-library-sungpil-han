# Refactor notes — `primary-sclerosing-cholangitis/psc_mrgsolve_model.R`

**Scope of this pass.** Per the fork's PK/PD refactor spec
(`FORK_WORKFLOW_GUIDE.md`, Part 2), all **three** of this file's census-row
compounds were rewritten, per their existing rows in
`driver-patches/data/compound_perturbation_census.md` (all three classified
"Redirect concentration (clean single site)"): **Bezafibrate (BEZ)**,
**Obeticholic acid (OCA)**, and **Ursodeoxycholic acid (UDCA)**. NorUDCA
(`NorUDCA_bile` compartment, `Kbile_NorUDCA` parameter) has no census row
and is therefore not a refactor target -- it is byte-identical to the
original. Every disease-side equation (gut-liver axis, FXR/bile-acid pool,
immune dynamics, cholangiocyte health, fibrosis scaffolding, portal
pressure, CCA risk) is untouched except where it reads a renamed drug
concentration or effect term.

## Archetype per compound

**UDCA and OCA: the same bespoke shape -- Archetype 3 variant (depot +
central), plus a bespoke bile tissue-site compartment instead of a
standard `PERI_<STEM>`.** Both compounds already model a third compartment
in the original (`UDCA_bile`/`OCA_bile`) fed unidirectionally from the
central compartment (`Kbile_UDCA * UDCA_plasma`/`Kbile_OCA * OCA_plasma`,
with its own first-order elimination, no back-flow to central) -- this is
structurally *not* a two-way distribution peripheral compartment (Archetype
2/3's `PERI_<STEM>`, which exchanges both directions), so it was renamed
`BILE_UDCA`/`BILE_OCA` rather than `PERI_UDCA`/`PERI_OCA`, following the
same "bespoke role name, not a standard one forced to fit" precedent as
`LEU_COL` in `familial-mediterranean-fever/fmf_refactor_notes.md` and
`TISSUE_DOXY` in `abdominal-aortic-aneurysm/aaa_refactor_notes.md`. Every
PD equation for both compounds (FXR activation, hydrophobic-BA-index
clearance for UDCA, ALP reduction for both) reads only this bile
compartment's derived concentration, never the central/plasma
compartment -- exactly the guide's stated exception ("two [exposed
concentrations] only when a genuinely different tissue site matters"; here
only one is actually PD-facing, the plasma value is not exposed at all
since the original never used it for anything but feeding the bile
compartment).

**Bezafibrate (BEZ): bespoke -- single compartment, no depot, PLUS a
constant zero-order input term.** The original's own `dxdt_BEZ_plasma =
Ka_BEZ * 1.0 - (CL_BEZ / Vd_BEZ) * BEZ_plasma` does not fit Archetype 1
cleanly: `Ka_BEZ` is applied as a bare zero-order input magnitude
(`Ka_BEZ * 1.0`, mg/h) rather than a true first-order absorption rate from
a depot compartment (there is no `BEZ_gut`) -- so BEZ_plasma always
accumulates toward a background steady state (`Ka_BEZ * V1_BEZ / CL_BEZ =
1.2` mg) regardless of whether the scenario doses BEZ at all, on top of
whatever bolus dosing is given directly into the same compartment. This
looks like an artifact of copy-pasting the norUDCA pattern
(`dxdt_NorUDCA_bile = Kbile_NorUDCA * 1.0 - ...`, also untouched, also out
of scope) rather than deliberate zero-order-infusion PK design, but per the
never-invent/never-flatten rule this was kept exactly as coded, only
renamed (`KA_BEZ`, `CL_BEZ` already conformed, `Vd_BEZ` -> `V1_BEZ`,
`BEZ_plasma` -> `CENT_BEZ`) -- not restructured into a "cleaner" single
compartment, and not given an invented depot. This quirk is flagged here
for visibility, not logged as an `UPSTREAM_ISSUES.md` defect: the model
still compiles and runs exactly as intended for this term (it isn't a
build defect), just an unusual PK design choice.

## Renaming applied

| Original | Refactored | Value | Role |
|---|---|---|---|
| `Ka_UDCA` | `KA_UDCA` | 0.8 (1/h) | absorption rate |
| `F_UDCA` | unchanged | 0.50 | bioavailability (already conformed) |
| `CL_UDCA` | unchanged | 12.0 (L/h) | clearance (already conformed) |
| `Vd_UDCA` | `V1_UDCA` | 25.0 (L) | central volume |
| `Kbile_UDCA` | `K_BILE_UDCA` | 0.15 (1/h) | central-to-bile transfer (bespoke role) |
| -- (inline literal `2.0`) | `V_BILE_UDCA` (new) | 2.0 | bile mg-to-conc proxy divisor, promoted from the literal in `UDCA_bile_conc = UDCA_bile / 2.0` |
| `UDCA_gut` (cmt) | `GUT_UDCA` | -- | depot |
| `UDCA_plasma` (cmt) | `CENT_UDCA` | -- | central |
| `UDCA_bile` (cmt) | `BILE_UDCA` | -- | bespoke bile tissue site |
| `UDCA_bile_conc` (local) | `C_UDCA` | -- | **the exposed concentration** (bile pool) |
| `Ka_OCA` | `KA_OCA` | 0.6 (1/h) | absorption rate |
| `F_OCA` | unchanged | 0.60 | bioavailability (already conformed) |
| `CL_OCA` | unchanged | 8.0 (L/h) | clearance (already conformed) |
| `Vd_OCA` | `V1_OCA` | 20.0 (L) | central volume |
| `Kbile_OCA` | `K_BILE_OCA` | 0.12 (1/h) | central-to-bile transfer (bespoke role) |
| -- (inline literal `0.5`) | `V_BILE_OCA` (new) | 0.5 | bile mg-to-conc proxy divisor, promoted from the literal in `OCA_bile_conc = OCA_bile / 0.5` |
| `OCA_gut` (cmt) | `GUT_OCA` | -- | depot |
| `OCA_plasma` (cmt) | `CENT_OCA` | -- | central |
| `OCA_bile` (cmt) | `BILE_OCA` | -- | bespoke bile tissue site |
| `OCA_bile_conc` (local) | `C_OCA` | -- | **the exposed concentration** (bile pool) |
| `IC50_OCA_FXR` | `EC50_OCA_FXR` | 0.10 (μmol/L) | Hill EC50 |
| -- (implicit, ratio max=1) | `EMAX_OCA_FXR` (new) | 1.0 | Hill Emax, "original ratio `C/(EC50+C)` implies max=1" |
| -- (none) | `GAMMA_OCA_FXR` (new) | 1.0 | Hill exponent, "no explicit Hill term in original" |
| `FXR_OCA` (local) | `EFFECT_OCA_FXR` | -- | OCA effect on FXR activation |
| -- (inline literal `0.25`) | `EMAX_OCA_ALP` (new) | 0.25 | Hill Emax, promoted from `OCA_ALP_eff` |
| -- (inline literal `2.0`) | `EC50_OCA_ALP` (new) | 2.0 | Hill EC50, promoted from `OCA_ALP_eff` |
| -- (none) | `GAMMA_OCA_ALP` (new) | 1.0 | Hill exponent |
| `OCA_ALP_eff` (local) | `EFFECT_OCA_ALP` | -- | OCA effect on ALP reduction |
| `Ka_BEZ` | `KA_BEZ` | 1.0 (mg/h, zero-order magnitude) | see bespoke note above |
| `CL_BEZ` | unchanged | 15.0 (L/h) | clearance (already conformed) |
| `Vd_BEZ` | `V1_BEZ` | 18.0 (L) | central volume |
| `BEZ_plasma` (cmt) | `CENT_BEZ` | -- | central (no depot) |
| -- (identity) | `C_BEZ` (new local) | -- | **the exposed concentration** -- equals `CENT_BEZ` raw amount (mg); the original never divided BEZ by a volume before using it in an effect term, so this is preserved exactly, not newly volume-normalized |
| -- (inline literal `0.3`) | `EMAX_BEZ_FIB` (new) | 0.3 | Hill Emax, promoted from `BEZ_antifib` |
| -- (inline literal `5.0`) | `EC50_BEZ_FIB` (new) | 5.0 | Hill EC50, promoted from `BEZ_antifib` |
| -- (none) | `GAMMA_BEZ_FIB` (new) | 1.0 | Hill exponent |
| (derived from `BEZ_antifib`) | `EFFECT_BEZ_FIB` (new) | -- | Bezafibrate antifibrotic effect (`BEZ_antifib = 1.0 - EFFECT_BEZ_FIB`, same formula) |
| -- (inline literal `0.20`) | `EMAX_BEZ_ALP` (new) | 0.20 | Hill Emax, promoted from `BEZ_ALP_eff` |
| -- (inline literal `10.0`) | `EC50_BEZ_ALP` (new) | 10.0 | Hill EC50, promoted from `BEZ_ALP_eff` |
| -- (none) | `GAMMA_BEZ_ALP` (new) | 1.0 | Hill exponent |
| `BEZ_ALP_eff` (local) | `EFFECT_BEZ_ALP` | -- | Bezafibrate effect on ALP reduction |
| `UDCA_hydro_eff` | `EMAX_UDCA_HYDRO` | 0.50 | Hill Emax (was already a named fraction parameter) |
| -- (inline literal `2.0`) | `EC50_UDCA_HYDRO` (new) | 2.0 | Hill EC50, promoted from `hydro_decrease` |
| -- (none) | `GAMMA_UDCA_HYDRO` (new) | 1.0 | Hill exponent |
| (derived from `hydro_decrease`) | `EFFECT_UDCA_HYDRO` (new) | -- | UDCA effect on hydrophobic-BA-index clearance |
| -- (inline literal `0.05`) | `EMAX_UDCA_FXR` (new) | 0.05 | Hill Emax, promoted from `FXR_UDCA` |
| -- (inline literal `5.0`) | `EC50_UDCA_FXR` (new) | 5.0 | Hill EC50, promoted from `FXR_UDCA` |
| -- (none) | `GAMMA_UDCA_FXR` (new) | 1.0 | Hill exponent |
| `FXR_UDCA` (local) | `EFFECT_UDCA_FXR` | -- | UDCA effect on FXR activation |

All parameter *values* are copied verbatim from the original (a promoted
magic-number literal keeps the exact number it replaces). The bile
compartments' own elimination-rate literals (`0.05` for `BILE_UDCA`,
`0.04` for `BILE_OCA`) were deliberately **left as inline literals**, not
promoted to named parameters -- they are plain PK elimination rates, not
part of the Hill effect interface the naming convention's promotion
requirement targets, and the original never named them either; promoting
every incidental PK-rate literal would be scope creep beyond what the
guide asks for. `NorUDCA_bile`/`Kbile_NorUDCA` are untouched (no census
row, not in scope).

## Hill interface: rename for six of seven effect terms, one Emax promoted from an implicit ratio

All seven of this file's drug effect terms were already exactly
`Emax * C^1 / (EC50^1 + C^1)`-shaped ratios in the original (or, for
`FXR_OCA`, `C / (EC50 + C)` with an *implicit* Emax of 1 -- the raw ratio's
own asymptote). Per the guide's rename rule, `EMAX_<STEM>_<target>`/
`EC50_<STEM>_<target>` were pulled out with the original's own values (or,
for `EMAX_OCA_FXR`, the value the un-scaled ratio already implied),
`GAMMA_<STEM>_<target> = 1` was added throughout (no compound/target had an
explicit Hill exponent), and every ratio was rewritten with `pow(x, 1.0)`
in place of the bare value -- `pow(x, 1)` is a required exact-identity
special case for IEEE 754 conforming `pow` implementations (returns `x`
exactly, per C99 Annex F / IEEE 60559), confirmed bit-identical by the
verification results below. **No `nls()` fit was needed or performed for
any of the three compounds** -- every rewrite here is a pure rename/
promotion of an already Emax/EC50-shaped ratio, never a refit.

**Each compound has multiple, genuinely distinct downstream effect
targets, kept as separate named terms, never collapsed.** UDCA drives
three independent effects from the same bile concentration (FXR
activation, hydrophobic-BA-index reduction, ALP reduction); OCA drives two
(FXR activation, ALP reduction); BEZ drives two (antifibrotic
collagen-synthesis inhibition, ALP reduction). These are genuinely
different physiological actions of the same drug on different disease
variables (not "several drugs on a combined multi-drug expression," the
case the guide's "combine only where disease equations actually use them"
targets), so each got its own `EFFECT_<STEM>_<target>` name and its own
Emax/EC50/Gamma triple -- same pattern as Colchicine's two effect terms in
`familial-mediterranean-fever/fmf_refactor_notes.md` and Encaleret's two in
`hypoparathyroidism/hypopt_refactor_notes.md`. Where the disease equation
combines multiple *compounds'* effects on the same target (ALP:
`ALP_drug_red = (EFFECT_OCA_ALP + EFFECT_UDCA_ALP + EFFECT_BEZ_ALP) * ALP`),
that summation is unchanged from the original's own
`(OCA_ALP_eff + UDCA_ALP_eff + BEZ_ALP_eff) * ALP` -- each compound's term
stays independently named and independently driveable, only combined at
the exact point the original already combined them.

## A subtlety caught during verification: `$MAIN` vs `$ODE` evaluation timing

The first verification pass (see below) found a small but real, non-zero
mismatch (max abs diff 0.03 mmHg on `PortalPressure`, localized entirely to
`Col1a1`/`Fibroscan`/`PortalPressure`/`CCA_risk` -- i.e. exactly BEZ's
antifibrotic downstream chain, nothing else) after declaring `double C_BEZ
= CENT_BEZ;` in `$MAIN`, mirroring where `C_OCA`/`C_UDCA` are declared.
This was **not** a rename bug in the Hill formula itself (`pow(x,1)==x`
holds exactly, confirmed once the real cause was found) -- it was an
evaluation-timing bug introduced by *where* the assignment was placed.
mrgsolve's `$MAIN` block runs once per dosing/observation interval, not at
every internal ODE solver step, whereas `$ODE` runs at every step. The
original computed `OCA_bile_conc`/`UDCA_bile_conc`/`FXR_input`/etc. in
`$MAIN` and fed them into `$ODE` (so hoisting `C_OCA`/`C_UDCA` into `$MAIN`
is a faithful structural match, confirmed exact below) -- but the original
computed `BEZ_antifib`/`BEZ_ALP_eff` **directly inside `$ODE`**, reading
the always-current `BEZ_plasma` at every internal step. Declaring `C_BEZ`
in `$MAIN` made it lag behind the true evolving `CENT_BEZ` within an
interval, a small but genuine numerical regression, not a pure rename.
**Fix:** `double C_BEZ = CENT_BEZ;` was moved into `$ODE`, immediately
after `dxdt_CENT_BEZ`, matching exactly where and how often the original
read `BEZ_plasma`. This is disclosed here because it is exactly the kind
of subtlety the guide's "a gate failure is a question, not a verdict"
principle exists for -- the first failing comparison correctly identified
a real difference in the rewrite, not a tolerance to loosen.

## Verification

**Method.** Both files' embedded mrgsolve `psc_model_code <- '...'` DSL
blocks were mechanically extracted (regex on the assignment, verbatim
quoted text; no separate `.cpp` file left behind, extraction was in-memory/
scratch-only), the original's copy was additionally patched with the
build-compat workaround below (in-memory only, never written back to the
tracked `psc_mrgsolve_model.R`), and both were POSTed to the local
qspserver `mrgsolve_api` service at `http://localhost:8007`
(`/model_manifest` then `/run_simulation`), which compiles and runs each
DSL block directly with mrgsolve 2.0.1 server-side -- no local R/mrgsolve
install used. Requests were spaced out sequentially (never fired
concurrently), respecting the service's stated `max_concurrent_jobs: 2`
limit and recent-crash history.

**Pre-existing upstream build defects, unrelated to any of the three
refactored compounds, patched for verification only.** Neither
`psc_mrgsolve_model.R` nor a pure rename of it compiles under mrgsolve
2.0.1 (two layered defects, logged as **upstream issue #54** in
`translations/UPSTREAM_ISSUES.md`): (1) `$INIT @annotated` uses plain
`NAME = value` entries instead of the required annotated `NAME : value :
description` triples; (2) once that tag is dropped, `$CMT @annotated` and
`$INIT` jointly redeclare all 27 compartments (same defect class as
issues #27/#34/#36/#42/#48/#49/#51). Fixing (2) via the modern
`<CMT>_0 = value;` `$MAIN` idiom surfaces one further incidental collision:
`$PARAM LOXL2_0` (an unused baseline value, 0.20) collides with the
auto-generated `LOXL2_0` init symbol mrgsolve creates for compartment
`LOXL2` (whose own init value, 0.35, was never actually the same number --
confirming the parameter was dead code even in the original). All three
fixes were applied identically to **both** the delivered
`psc_mrgsolve_model_refactored.R` (per this fork's settled build-compat
policy, since a "verified" deliverable nobody can run afterwards defeats
the point) and, in-memory only, to the untouched original for this
verification run: `$INIT` deleted, its 27 assignments moved into `$MAIN`
via the `<CMT>_0` idiom (values copied verbatim, using each file's own
compartment names), and `LOXL2_0` renamed `LOXL2_0_BASE` in both. None of
these three fixes touches a single number, equation, or PK/PD behavior --
confirmed by the exact-match results below, which is the actual proof, not
just the description.

**Results: exact match on both of the file's own dosing scenarios run for
their full stated duration.**

1. **"UDCA + OCA Combination"** (`build_dosing("UDCA", 900, freq_hours=12)`
   + `build_dosing("OCA", 10, freq_hours=24)`, i.e. 450 mg q12h x730 doses
   into `GUT_UDCA`/cmt 1 plus 10 mg q24h x365 doses into `GUT_OCA`/cmt 4),
   run for the full stated 5-year window (`end=43800h`, `delta=168h`,
   matching the original's own `obs_times`, 263 timepoints, no shortening
   needed -- this model's step count stayed well within the API's default
   `maxsteps` budget): **max abs deviation 1.0e-17** across all 32 shared
   `$CMT`/`$CAPTURE` outputs (i.e. floating-point noise on values already
   at the ~1e-13 level; every output at physiologically meaningful
   magnitude matched to displayed precision, e.g. `Col1a1`/`PortalPressure`
   trajectories bit-identical after the `$MAIN`/`$ODE` fix above). This run
   also exercises BEZ's own background zero-order kinetics (undosed in
   this scenario, but always running per its `KA_BEZ * 1.0` term) -- also
   exact.
2. **"Bezafibrate 400 mg/day"** (`build_dosing("BEZ", 400, freq_hours=24)`,
   400 mg q24h x365 doses into `CENT_BEZ`/cmt 8), same 5-year window:
   **max abs deviation 0.0** (exact) across all 32 shared outputs,
   including `EFFECT_BEZ_ALP`/`EFFECT_BEZ_FIB` (captured, PD-facing) and
   `ALP`/`Col1a1`/`Fibroscan`/`PortalPressure`/`CCA_risk` (the downstream
   chain that had previously exposed the `$MAIN`/`$ODE` timing bug above).

Compartment indices are identical between the original and the refactored
model (`GUT_UDCA`=1 ... `CENT_BEZ`=8, `LPS`=9 ... `CCA_risk`=27), since
only compartment *names* changed, never their `$CMT` declaration order.
This is the expected, and achieved, result for three compounds that are
pure structural renames/promotions per the guide's tolerance table for
Archetypes 1-3 and their bespoke variants: "expect a near-exact match...
anything beyond floating-point-scale deviation means a bug" -- both
scenarios match at or below floating-point noise, not merely "near."
("Natural Progression (No Treatment)" and "OCA 10 mg/day" were not
separately POSTed: the former is a strict subset of what scenario 1 above
already exercises for BEZ's background kinetics and the untouched disease
network, and the latter's OCA dosing is a strict subset of scenario 1's
combination dosing.)

## qspserver `/model_manifest` discoverability

Confirmed via `POST /model_manifest` against the extracted, renamed DSL:
`KA_UDCA`, `F_UDCA`, `CL_UDCA`, `V1_UDCA`, `K_BILE_UDCA`, `V_BILE_UDCA`,
`KA_OCA`, `F_OCA`, `CL_OCA`, `V1_OCA`, `K_BILE_OCA`, `V_BILE_OCA`,
`EC50_OCA_FXR`, `EMAX_OCA_FXR`, `GAMMA_OCA_FXR`, `KA_BEZ`, `CL_BEZ`,
`V1_BEZ`, `EMAX_UDCA_HYDRO`, `EC50_UDCA_HYDRO`, `GAMMA_UDCA_HYDRO`,
`EMAX_UDCA_FXR`, `EC50_UDCA_FXR`, `GAMMA_UDCA_FXR`, `EMAX_BEZ_FIB`,
`EC50_BEZ_FIB`, `GAMMA_BEZ_FIB`, `EMAX_OCA_ALP`, `EC50_OCA_ALP`,
`GAMMA_OCA_ALP`, `EMAX_UDCA_ALP`, `EC50_UDCA_ALP`, `GAMMA_UDCA_ALP`,
`EMAX_BEZ_ALP`, `EC50_BEZ_ALP`, `GAMMA_BEZ_ALP` all appear in the
manifest's `parameters` with their original numeric defaults.
`C_OCA`/`C_UDCA`/`C_BEZ` and all seven `EFFECT_<STEM>_<target>` terms are
state-derived (computed in `$MAIN`/`$ODE` from compartment concentrations),
so -- per the same reasoning established for `EFFECT_COL_NEU`/`EFFECT_ANA`/
`EFFECT_CANA` in `familial-mediterranean-fever/fmf_refactor_notes.md` --
they cannot also be `$PARAM` entries; all ten appear in the manifest's
`outputPaths` via the file's own `$CAPTURE @annotated` list, confirmed
discoverable.

## Anything else flagged

- No compound other than Bezafibrate, Obeticholic acid, and Ursodeoxycholic
  acid was touched. NorUDCA's compartment (`NorUDCA_bile`) and parameter
  (`Kbile_NorUDCA`) are byte-identical to the original; it has no census
  row.
- The R-side `build_dosing()` helper's `cmt_map` was updated only for the
  three renamed compartments it references (`"UDCA_gut"` ->
  `"GUT_UDCA"`, `"OCA_gut"` -> `"GUT_OCA"`, `"BEZ"` ->
  `"CENT_BEZ"`; `"norUDCA"` -> `"NorUDCA_bile"` unchanged) -- same dosing
  amounts, same timing, same scenario list, same plots throughout. No
  other R-side code references a renamed identifier (checked by grepping
  the full original file for every old PK name outside the DSL block).

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`, the
`primary-sclerosing-cholangitis | BEZ`, `primary-sclerosing-cholangitis |
OCA`, and `primary-sclerosing-cholangitis | UDCA` rows.
