# Refactor notes — `vte_mrgsolve_model.R` (apixaban, dabigatran, enoxaparin, rivaroxaban, warfarin)

**Scope of this pass.** Per the fork's PK/PD refactor spec
([`FORK_WORKFLOW_GUIDE.md`](../FORK_WORKFLOW_GUIDE.md), Part 2), all five of
this file's compounds were renamed to the fork's pluggable-PK convention, per
their existing rows in
[`driver-patches/data/compound_perturbation_census.md`](../driver-patches/data/compound_perturbation_census.md)
(all five classified "Redirect concentration (clean single site)"): **APIX**
(apixaban, stem `APIX`), **DABI** (dabigatran, stem `DABI`), **ENOX**
(enoxaparin, stem `ENOX`), **RIV** (rivaroxaban, stem `RIV`), and **WARF**
(warfarin, stem `WARF`). The file models no other compounds, so nothing was
left untouched at the compound level — but the disease-side coagulation
cascade (FXa/FIIa activity, fibrin, clot burden, plasmin/fibrinolysis,
D-dimer, INR/aPTT reporting, the factor-synthesis pools) is unchanged, as is
every dosing scenario in the R driver code below the model string (only the
renamed compartment names and output columns needed to change there).

## Archetypes determined

Only rivaroxaban has a genuine absorption depot in this file; the other four
compounds are dosed directly into their own single compartment with two
first-order loss terms summed together (the original's own `KA_<STEM>` acts
as an *extra elimination term*, not an absorption rate from a separate
depot — there is no `*_GUT` compartment for APIX/DABI/WARF/ENOX anywhere in
the original's `$CMT` block).

### Rivaroxaban (`RIV`) — Archetype 3 (depot + central + peripheral, linear)

```
dxdt_GUT_RIV   = -KA_RIV * GUT_RIV;
dxdt_CENT_RIV  =  KA_RIV * GUT_RIV - (CL_RIV+Q_RIV)/V1_RIV*CENT_RIV + Q_RIV/V2_RIV*PERI_RIV;
dxdt_PERI_RIV  =  Q_RIV/V1_RIV*CENT_RIV - Q_RIV/V2_RIV*PERI_RIV;
```
Matches Archetype 3 exactly, already using CL/Q/V-style parameters. `EFFECT_RIV`
(renamed from `INH_FXa_RIV`) is already the guide's canonical Hill ratio —
rename only, `GAMMA_RIV = 1.3` (the original's own `HILL_RIV`, not 1).

### Apixaban / Dabigatran / Warfarin / Enoxaparin — Archetype 1 (no depot, single compartment)

All four share the identical structural quirk: dosed directly into their one
compartment (bioavailability applied by the R driver at the dosing event,
`amt = dose * F_<STEM>`, never inside the DSL), with two summed first-order
loss terms:

```
dxdt_CENT_APIX = -KA_APIX * CENT_APIX - (CL_APIX/V1_APIX) * CENT_APIX;
dxdt_CENT_DABI = -KA_DABI * CENT_DABI - (CL_DABI/V1_DABI)/RF_adj_DABI * CENT_DABI;
dxdt_CENT_WARF = -KA_WARF * CENT_WARF - (CL_WARF/V1_WARF) * CENT_WARF;
dxdt_CENT_ENOX = -KA_ENOX * CENT_ENOX - (CL_ENOX/V1_ENOX)/RF_adj_ENOX * CENT_ENOX;
```

`KA_<STEM>` here behaves as a second elimination rate constant, not an
absorption rate (mathematically equivalent to a single combined elimination
constant `KA_<STEM> + CL_<STEM>/V1_<STEM>`, since both terms act on the same
compartment) — this is exactly "no depot, single compartment, linear
elimination," Archetype 1, just with the elimination expressed as two summed
named rate constants rather than one. Preserved exactly as the original wrote
it; not restructured into a genuine depot model (that would be adding
pharmacology the original never had) and not collapsed into a single combined
rate constant either (that would lose the original's own parameter names).
`DABI` and `ENOX` additionally divide their clearance term by a renal-function
covariate (`RF_adj_DABI = pow(eGFR_pat/eGFR_ref, 0.85)`,
`RF_adj_ENOX = pow(eGFR_pat/eGFR_ref, 0.65)`) — untouched, not part of the
compound-identity naming convention.

**Renaming applied (compartments and PK parameters, all five compounds):**

| Original | Refactored |
|---|---|
| `RIV_GUT` / `RIV_CENT` / `RIV_PERIPH` | `GUT_RIV` / `CENT_RIV` / `PERI_RIV` |
| `APIX_CENT` | `CENT_APIX` |
| `DABI_CENT` | `CENT_DABI` |
| `WARF_CENT` | `CENT_WARF` |
| `ENOX_CENT` | `CENT_ENOX` |
| `Cp_RIV` / `Cp_APIX` / `Cp_DABI` / `Cp_WARF` / `Cp_ENOX` | `C_RIV` / `C_APIX` / `C_DABI` / `C_WARF` / `C_ENOX` |
| `INH_FXa_RIV` / `INH_FXa_APIX` / `INH_FXa_ENOX` / `INH_FIIa_DABI` | `EFFECT_RIV` / `EFFECT_APIX` / `EFFECT_ENOX` / `EFFECT_DABI` |
| `HILL_RIV` / `HILL_APX` / `HILL_DABI` / `HILL_W` | `GAMMA_RIV` / `GAMMA_APIX` / `GAMMA_DABI` / `GAMMA_WARF` |
| `IC50_WARF` | `EC50_WARF` |
| `KIN_VK` / `KOUT_VK` / `VK0_ox` / `VK0_red` | `KIN_VK_WARF` / `KOUT_VK_WARF` / `VK0_ox_WARF` / `VK0_red_WARF` |
| — (none; original had no Hill coefficient) | `GAMMA_ENOX` (new, `= 1`, documented below) |

`KA_<STEM>`, `F_<STEM>`, `CL_<STEM>`, `V1_<STEM>`, `V2_RIV`, `Q_RIV`,
`EMAX_RIV`/`EMAX_APIX`/`EMAX_DABI`/`EMAX_ENOX`, `EC50_RIV`/`EC50_APIX`/
`EC50_DABI`/`EC50_ENOX` were **already** named exactly per this guide's
convention in the original (a pleasant exception to the corpus's usual
naming chaos) — left unchanged.

`F_<STEM>` is declared in `$PARAM` for every compound but is **not** used
inside `$ODE` for any of them — the original applies bioavailability by
pre-multiplying the R driver's dose amount (`amt = 15 * 0.93`, etc.) rather
than via mrgsolve's `F_<cmt>` mechanism or an `F_<STEM>*GUT` term inside the
absorption flow. Preserved exactly as the original does it (same amounts used
in the verification dosing below); not converted to a DSL-side bioavailability
mechanism, since that would change how dosing must be constructed against
this model versus the original.

**GAMMA_ENOX (new parameter).** The original's enoxaparin effect term had no
explicit Hill exponent: `INH_FXa_ENOX = (EMAX_ENOX * Cp_ENOX) / (EC50_ENOX +
Cp_ENOX)`. `EFFECT_ENOX` is written in the guide's canonical
`EMAX*pow(X,GAMMA)/(pow(EC50,GAMMA)+pow(X,GAMMA))` form with `GAMMA_ENOX = 1`
— algebraically identical to the original's plain ratio for any concentration
(confirmed in verification below), a rename/canonicalization, not a refit.

## A genuine PD defect found and preserved-but-flagged (not a build defect)

The original's warfarin Hill term reused **rivaroxaban's own `EMAX_RIV`**
as warfarin's effect ceiling — no `EMAX_WARF` parameter existed anywhere in
the file:

```c
double WARF_INH = (EMAX_RIV * pow(Cp_WARF, HILL_W)) /
                  (pow(IC50_WARF, HILL_W) + pow(Cp_WARF, HILL_W));
```

This is a cross-compound parameter leak — exactly the kind of coupling this
refactor's "isolate each compound's own PK/PD block" goal is meant to
eliminate — and it also violates the guide's requirement that a compound's
own effect be expressed as "one named function of that concentration... not
buried inside a combined multi-drug expression." Rather than either (a)
literally reproducing the cross-reference in the refactored file (which would
defeat independent driveability: overriding `EMAX_RIV` for a rivaroxaban
sensitivity analysis would silently also change warfarin's effect) or (b)
inventing an arbitrary new value for `EMAX_WARF` (a refit, out of scope for a
rename-only pass), a new parameter `EMAX_WARF = 0.97` was added carrying the
**identical** value `EMAX_RIV` already defaults to. `EFFECT_WARF` is
therefore numerically identical to the original's `WARF_INH` for every
scenario in this file (none of which override `EMAX_RIV`'s default,
confirmed in the verification below) while being properly self-contained.
Logged as **upstream issue #84** in
[`translations/UPSTREAM_ISSUES.md`](../translations/UPSTREAM_ISSUES.md) —
not fixed in the checked-in original, per the shared never-edit-upstream
rule.

A second, purely cosmetic observation from the same block: `INH_FVII`,
`INH_FX`, `INH_FII` (computed from `VK_ratio` in `$MAIN`) are dead code in
the original — no `dxdt_F*_POOL` line reads them; those instead recompute
`STIM_VK = VK_RED/VK0_red` independently in `$ODE`. Preserved as-is (not
removed — deleting unused-but-harmless code is out of scope for a rename-only
refactor); noted here for the record.

## When the original doesn't compile at all

**Neither `vte_mrgsolve_model.R` nor a naive rename-only version of it
compiles under mrgsolve 2.0.1.** Three independent, pre-existing build
defects were found, all unrelated to any compound's own PK — logged as
**upstream issue #83** in `UPSTREAM_ISSUES.md` — and all three were
given the syntax-only, non-numeric fix the workflow guide's policy calls for,
applied directly in this `_refactored.R` sibling (never in the checked-in
original):

**1. `$PARAM`'s `VK0_red` computed its default from a sibling `$PARAM` name.**
```
VK0_red  = VK0_ox * KIN_VK / KOUT_VK  // steady-state ratio
```
mrgsolve evaluates every `$PARAM` block's defaults as a single R
`list(...)` call — list elements have no access to each other during that
evaluation — so this fails with `object 'VK0_ox' not found`. **Fix:**
replaced with the literal value that expression evaluates to,
`VK0_red_WARF = 0.7777777777777778` (`= 1.0 * 0.14 / 0.18`). No scenario in
this file (or in the verification below) ever overrides `VK0_ox`/`KIN_VK`/
`KOUT_VK`, so this is numerically identical, not just structurally
equivalent.

**2. `$CMT` and `$INIT` jointly redeclared all 18 compartments.** Every
single compartment named in `$CMT` was *also* assigned an initial value in a
separate `$INIT` block — valid in older mrgsolve, but this build rejects it
outright: `invalid class "mrgmod" object: Duplicated model names: RIV_GUT
RIV_CENT ... FII_POOL`. **Fix:** removed `$INIT` entirely; every initial
value now set in `$MAIN` via the `<cmt>_0` idiom, guarded by
`if (NEWIND <= 1) { ... }` (same pattern as
`von-willebrand-disease/vwd_mrgsolve_model_refactored.R`). Values are
byte-identical to the original's `$INIT` block, including the ones that
referenced other parameters (`CLOT_SIZE_0 = CLOT_init;`,
`PLASMIN_ACT_0 = PLASMIN_base;`, etc.) — `$MAIN`, unlike `$PARAM`/`$INIT`,
has full parameter scope, so these no longer need the `VK0_red`-style
literal workaround.

**3. Six `$MAIN` doubles were re-declared under the identical name by a
self-referential `$TABLE` `capture` line.** `Cp_RIV`, `Cp_APIX`, `Cp_DABI`,
`Cp_WARF`, `ANTI_XA`, and `INR` are each computed as a `double` in `$MAIN`,
then the original's `$TABLE` block re-declares `capture Cp_RIV = Cp_RIV;`
(and the same pattern for the other five) — mrgsolve auto-promotes every
`$MAIN` double into the same anonymous-namespace scope as `$TABLE` captures,
so a same-named capture collides:
```
67:11: error: redefinition of 'capture {anonymous}::Cp_RIV'
30:10: note: 'double {anonymous}::Cp_RIV' previously declared here
```
This is the same defect *class* as `UPSTREAM_ISSUES.md` #43 (DIC's
`FREETPA`), but at `$MAIN` scope rather than `$ODE`. **Fix:** the six
redundant `$TABLE` lines were dropped; the bare names were already
independently listed in `$CAPTURE`, and — since a `$MAIN` double is
auto-promoted regardless of whether a `$TABLE` line also references it —
they remain fully reportable without the redundant capture line (confirmed
present in `/model_manifest`'s `outputPaths` and retrievable via
`/run_simulation`, see below).

**None of these three fixes are numeric or behavioral** — see Verification
below, which proves it empirically (max relative deviation 0.0 across every
shared output, in every one of the file's own 7 dosing scenarios).

## qspserver `/model_manifest` discoverability

- `KA_RIV`, `F_RIV`, `CL_RIV`, `V1_RIV`, `Q_RIV`, `V2_RIV`, `EMAX_RIV`,
  `EC50_RIV`, `GAMMA_RIV`; the equivalent `_APIX`/`_DABI`/`_ENOX` sets;
  `EC50_WARF`, `GAMMA_WARF`, `EMAX_WARF` (new), `KIN_VK_WARF`,
  `KOUT_VK_WARF`, `VK0_ox_WARF`, `VK0_red_WARF`; and `GAMMA_ENOX` (new) are
  all declared in `$PARAM`, confirmed present with their defaults in the
  live `/model_manifest` response.
- `C_RIV`/`C_APIX`/`C_DABI`/`C_WARF`/`C_ENOX` and
  `EFFECT_RIV`/`EFFECT_APIX`/`EFFECT_DABI`/`EFFECT_WARF`/`EFFECT_ENOX` are
  computed in `$MAIN` from compartment state, so — as already found for
  `EFFECT_NIM` (`sah_refactor_notes.md`), `EFFECT_TCZ`
  (`ted_refactor_notes.md`), and `EFFECT_ATR`/`EFFECT_TXA`
  (`dic_refactor_notes.md`) — they cannot also be `$PARAM` entries (a
  `$PARAM` parameter is a static default; a `$MAIN`-computed quantity under
  the identical name is a redefinition, the same class of error as build
  defect #3 above). Exposed instead as bare `$CAPTURE` entries (no `$TABLE`
  line needed, for the same reason removing the six redundant lines above
  didn't lose discoverability) — confirmed present in `/model_manifest`'s
  `outputPaths` and retrievable via `/run_simulation`.
- `model_content` was extracted from the `vte_model_code <- '...'` R-string
  wrapper (straight text extraction via regex on the assignment, verbatim
  quoted text, no character changes) purely to build the verification
  requests below; no `.cpp` sibling was left behind, per the workflow guide.
- One apostrophe-hygiene note specific to this refactor: three of the new
  comments drafted for this pass originally used a plain apostrophe
  (`warfarin's`, `original's`) — inside the DSL's enclosing R single-quoted
  string (`vte_model_code <- '...'`), an apostrophe closes the string early,
  exactly the defect class in `UPSTREAM_ISSUES.md` #1
  (`acute-bacterial-meningitis`). Caught by `Rscript`'s own `parse()` on this
  file before verification (`67:11: unexpected symbol`, tracing to the
  apostrophe) and fixed by rewording (`warfarins`, `originals`) rather than
  switching the block to double quotes, to keep the diff against the
  original's own quoting style minimal. `parse()` on both files now succeeds
  with an identical expression count (63 each).

## Verification

**Method.** Both the original's and the refactored file's embedded mrgsolve
DSL blocks were mechanically extracted (regex on the `<name>_code <- '...'`
assignment, verbatim quoted text) and POSTed to the local qspserver
`mrgsolve_api` service at `http://localhost:8007` (`/model_manifest` then
`/run_simulation`), which compiles and runs each DSL block directly with
mrgsolve 2.0.1 server-side — no local R/mrgsolve install used for the model
comparison itself (a separate, syntax-only `Rscript parse()` check was run
directly on both `.R` files, per the note above, to catch the apostrophe
defect — not a model run). Requests were spaced ~2 s apart per the API's
`max_concurrent_jobs: 2` limit. The **original's own build defects (#83
above) were also applied to the comparison copy of the true original** —
identically, syntax-only, non-numeric — purely so it would compile
side-by-side with the refactored file; **the checked-in
`vte_mrgsolve_model.R` on disk contains none of these three fixes** and does
not currently build against mrgsolve 2.0.1 without them.

All **7 of this file's own dosing scenarios** (its own R driver code,
`make_dosing()` plus the two dabigatran renal variants) were run through
both models, full `end`/`delta` windows exactly as the original specifies —
no shortening was needed; every run completed well within the API's default
step budget:

| Scenario | Dosing | Window | Param overrides |
|---|---|---|---|
| S1 — DVT: rivaroxaban | 15 mg BID×21d (×0.93 F) → 20 mg QD | 90 d, δ=1 h | none (defaults) |
| S2 — PE: apixaban | 10 mg BID×7d (×0.50 F) → 5 mg BID | 90 d, δ=1 h | `PLASMIN_base=0.4` |
| S3 — warfarin+enoxaparin bridge | warfarin 5 mg QD + enoxaparin 1 mg/kg BID×10d (×0.92 F) | 90 d, δ=1 h | none (defaults) |
| S4 — enoxaparin prophylaxis | 40 mg QD×14d (×0.92 F) | 14 d, δ=0.5 h | `CLOT_init=0` |
| S5 — extended rivaroxaban | 10 mg QD×180d (×0.93 F) | 180 d, δ=2 h | `CLOT_init=0` |
| S6a — dabigatran, normal GFR | 110 mg BID×90d (×0.065 F) | 90 d, δ=1 h | `eGFR_pat=90` |
| S6b — dabigatran, CKD3 | same dosing as S6a | 90 d, δ=1 h | `eGFR_pat=30` |

16 shared `$CAPTURE` outputs were compared per scenario (`C_RIV`/`C_APIX`/
`C_DABI`/`C_WARF` vs. the original's `Cp_RIV`/`Cp_APIX`/`Cp_DABI`/`Cp_WARF`,
plus `ANTI_XA`, `INR`, `aPTT_out`, `INH_FXa`, `INH_FIIa`, `TG_ETP`, `DDIMER`,
`CLOT_PCT2`, `FIBRIN_mg`, `FVII_pct`, `FX_pct`, `FII_pct`), across the full
time grid, plus a spot-check of renamed compartment state
(`GUT_RIV`/`CENT_RIV`/`PERI_RIV`/`CENT_APIX`/`CENT_DABI`/`CENT_WARF`/
`CENT_ENOX` vs. the originals) for scenario S3.

**Results: exact match on every scenario.** Max relative deviation
**0.0** across all 16 shared outputs, in all 7 scenarios, including the
two scenarios (S4, S5) whose `CLOT_init=0` override drives `CLOT_PCT2`
(`= CLOT_SIZE/CLOT_init*100`) to `NaN`/`Inf` from time 0 — confirmed this
happens **identically** in the original (reproduced independently with the
build-fixed original alone before comparing), so it is a pre-existing
property of dividing by a zero `CLOT_init`, not something the refactor
introduced; both sides produce the identical `NaN`/`Inf` pattern at the
identical timepoints. Renamed-compartment state (S3: `CENT_WARF` vs.
`WARF_CENT`, `CENT_ENOX` vs. `ENOX_CENT`, `VK_RED`, `FVII_POOL`) matched to
machine precision as well. `EFFECT_WARF`'s tail value (0.434, non-zero) and
`EFFECT_ENOX`'s tail value (0, drug fully cleared by day 90) were spot-checked
to confirm the newly-exposed effect terms are genuinely active, not
degenerate no-ops.

Per the guide's tolerance table: Archetypes 1 and 3 are pure structural
reorganization with no Hill-fitting, and reproduced bit-for-bit — the
expected outcome, and the strongest possible confirmation that (a) the
renaming introduced no behavioral change and (b) the `EMAX_WARF` defect
recovery (same 0.97 value, correctly scoped) is truly numerically inert
under every one of this file's own scenarios.
