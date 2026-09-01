# Refactor notes — `x-linked-hypophosphatemia/xlh_mrgsolve_model.R`

All three compounds modeled in this file were in scope, per the existing
rows in `driver-patches/data/compound_perturbation_census.md` (all three
classified "Redirect concentration (clean single site)"): **Burosumab**,
**Oral calcitriol**, and **Oral phosphate (PHOSORAL)**. There is no
disease-model code outside these three compounds' own PK/PD blocks in this
file, so nothing else changed except the two build-compat fixes described
below.

## Archetype per compound

**Burosumab — Archetype 3** (depot + central + peripheral, linear), matched
exactly, no deviation:

```
dxdt_GUT_BURO  = -KA_BURO * GUT_BURO;
dxdt_CENT_BURO =  KA_BURO * GUT_BURO - CL_BURO/V1_BURO*CENT_BURO
                    - (Q_BURO/V1_BURO)*CENT_BURO + (Q_BURO/V2_BURO)*PERI_BURO;
dxdt_PERI_BURO =  (Q_BURO/V1_BURO)*CENT_BURO - (Q_BURO/V2_BURO)*PERI_BURO;
```

**Checked explicitly for TMDD/receptor-binding kinetics, per the task
brief, since burosumab is an anti-FGF23 IgG1 mAb** — a class that often
warrants archetype 4. It does not have any here: there is no free-FGF23 or
drug-FGF23-complex compartment anywhere in the file, and FGF23 neutralization
(`FGF23_NEUT`, renamed `EFFECT_BURO`) is computed as a plain algebraic Hill
function of concentration (`pow(C_BURO,...)/(pow(EC50_BURO,...)+pow(C_BURO,...))`),
not from a receptor-occupancy state solved by its own ODEs. The original
models burosumab's mechanism as direct concentration-driven neutralization,
not target-mediated disposition — confirmed by reading the actual `$ODE`
block, not assumed from the drug class or the file's header comment (which
says only "directly neutralizes circulating FGF23," consistent with what
the code does).

**Oral calcitriol and Oral phosphate — bespoke, same pattern for both.**
Neither fits archetypes 1–3 cleanly: each has an absorption depot (mass,
mg/µg) feeding directly into a *second* compartment that is not a
mass-over-volume central compartment but is itself the concentration/
exposure-signal state PD reads:

```
dxdt_GUT_CALC = -KA_CALC * GUT_CALC;
dxdt_C_CALC   =  KA_CALC * GUT_CALC / V_CALC - KE_CALC * C_CALC;
```

(Oral phosphate is structurally identical, `GUT_PHOSORAL` → `C_PHOSORAL`.)
This is a valid, common apparent-volume-normalized-absorption
simplification — the division by `V_CALC`/`V_PHOSORAL` happens once, inside
the absorption term, rather than as a separate `C = CENT/V` algebraic line
after a proper mass-balance central compartment. No compartment was added
or removed to force this into archetype 3; the structure is preserved
exactly, only renamed. Per the guide's "None of these fit" edge case, this
is flagged here as bespoke rather than mis-labeled as archetype 3.

## Renaming

| Original | Refactored |
|---|---|
| `BURO_DEPOT` | `GUT_BURO` |
| `BURO_CENT` | `CENT_BURO` |
| `BURO_PERIPH` | `PERI_BURO` |
| `BURO_CP` | `C_BURO` |
| `EC50_BURO_FGF23` | `EC50_BURO` |
| `HILL_BURO` | `GAMMA_BURO` |
| `FGF23_NEUT` | `EFFECT_BURO` |
| `PHOSORAL_GUT` | `GUT_PHOSORAL` |
| `PHOSORAL_SIG` | `C_PHOSORAL` |
| `EMAX_PHOS_ORALBOOST` | `EMAX_PHOSORAL` |
| `PHOS_ORALBOOST` | `EFFECT_PHOSORAL` |
| `CALC_GUT` | `GUT_CALC` |
| `CALC_CENT` | `C_CALC` |
| `EMAX_CALC_BOOST` | `EMAX_CALC` |
| `CALC_EXOG_BOOST` | `EFFECT_CALC` |

`KA_BURO`, `CL_BURO`, `V1_BURO`, `V2_BURO`, `Q_BURO`, `F_BURO`,
`KA_PHOSORAL`, `KE_PHOSORAL`, `V_PHOSORAL`, `F_PHOSORAL`, `KA_CALC`,
`KE_CALC`, `V_CALC`, `F_CALC`, and `EC50_PHOSORAL`/`EC50_CALC` already
matched the convention and are unchanged. All parameter *values* are copied
verbatim from the original — nothing invented or defaulted.

Two mechanical consequences of renaming the compartments, both required for
the model to still mean the same thing (not scope creep):
- The `$MAIN` bioavailability variables, which mrgsolve matches to a
  compartment by name (`F_<cmtname>`), were updated from
  `F_BURO_DEPOT`/`F_PHOSORAL_GUT`/`F_CALC_GUT` to
  `F_GUT_BURO`/`F_GUT_PHOSORAL`/`F_GUT_CALC` to track the renamed
  compartments — the underlying `F_BURO`/`F_PHOSORAL`/`F_CALC` user
  parameters they read from are untouched.
- The R-side `scenarios` list's `cmt =` dosing targets were updated to the
  new compartment names (`"GUT_BURO"`, `"GUT_PHOSORAL"`, `"GUT_CALC"`).
  Compartment order (and therefore 1-based compartment number, used by the
  qspserver API's `DoseSpec.cmt`) is unchanged.

Added (new, not present in the original, purely to complete the named Hill
interface — see below): `EMAX_BURO`, `GAMMA_PHOSORAL`, `GAMMA_CALC`.

## Hill interface: rename in all three cases, not a fit

All three compounds' disease effects were already written as plain
`Emax*X^gamma/(EC50^gamma+X^gamma)`-shaped ratios in the original — no
ODE-solved kinetics to approximate, so this is a pure rename per the guide's
"already this shape" rule, not the TMDD-style fitting procedure.

- **Burosumab**: `FGF23_NEUT = pow(BURO_CP, HILL_BURO) / (pow(EC50_BURO_FGF23,
  HILL_BURO) + pow(BURO_CP, HILL_BURO))`. This already has the exact target
  form with an implicit `Emax=1` (no separate multiplier in front). Renamed
  to `EFFECT_BURO = EMAX_BURO * pow(C_BURO, GAMMA_BURO) / (pow(EC50_BURO,
  GAMMA_BURO) + pow(C_BURO, GAMMA_BURO))` with `EMAX_BURO = 1` added
  explicitly (multiplying by 1 changes nothing) so the interface is
  complete without needing to consult the DSL text to know Emax is fixed
  at 1.
- **Oral phosphate**: `PHOS_ORALBOOST = EMAX_PHOS_ORALBOOST * PHOSORAL_SIG /
  (EC50_PHOSORAL + PHOSORAL_SIG)` — already Hill-shaped with `gamma=1`
  implicit (no `pow()` in the original since gamma=1 collapses it to a
  simple ratio). Renamed to `EFFECT_PHOSORAL = EMAX_PHOSORAL *
  pow(C_PHOSORAL, GAMMA_PHOSORAL) / (pow(EC50_PHOSORAL, GAMMA_PHOSORAL) +
  pow(C_PHOSORAL, GAMMA_PHOSORAL))` with `GAMMA_PHOSORAL = 1` added
  explicitly. `pow(x, 1.0) == x` exactly (IEEE-754 guarantee for `pow` with
  `y=1`), so this introduces no floating-point drift — confirmed
  empirically below (exact match, not merely close).
- **Oral calcitriol**: same pattern, `CALC_EXOG_BOOST = EMAX_CALC_BOOST *
  CALC_CENT / (EC50_CALC + CALC_CENT)` → `EFFECT_CALC = EMAX_CALC *
  pow(C_CALC, GAMMA_CALC) / (pow(EC50_CALC, GAMMA_CALC) + pow(C_CALC,
  GAMMA_CALC))`, `GAMMA_CALC = 1` added explicitly.

No curve-fitting (`nls()` or otherwise) was needed or performed for any of
the three compounds — unlike the tocilizumab/TMDD calibration case this
guide is built around, none of these three compounds has ODE-solved
receptor kinetics to approximate; each was already the target closed form.

## Two pre-existing build defects found while verifying (not fixed
upstream, logged)

Confirmed on the **untouched original**, independently, via `POST
/model_manifest` against the original file's own extracted DSL — neither
introduced by this refactor:

1. **`$CAPTURE` lists fourteen compartment names.** The original's
   `$CAPTURE` line repeats fourteen names that are also `$CMT` compartments
   (`NPT2 TMPGFR PHOS CALCITRIOL PTH BSAP RSS HEIGHTZ_XLH SIXMWT WOMAC UCACR
   NEPHROCALC PHOSORAL_SIG CALC_CENT`). mrgsolve 2.0.1 refuses to build any
   model whose `$CAPTURE` repeats a compartment name:
   ```
   Error in validObject(.Object) :
     invalid class "mrgmod" object: compartment should not be in $CAPTURE: NPT2,TMPGFR,PHOS,CALCITRIOL,PTH,BSAP,RSS,HEIGHTZ_XLH,SIXMWT,WOMAC,UCACR,NEPHROCALC,PHOSORAL_SIG,CALC_CENT
   ```
2. **`$MAIN`'s `if (NEWIND <= 1) { ... }` block assigns directly to
   compartment names** (`NPT2 = NPT2_BASE;`, `TMPGFR = TMPGFR0;`, ... 12
   assignments total). This second defect only surfaces once defect 1 is
   worked around — under `$PLUGIN autodec` in this mrgsolve build, a bare
   compartment name in `$MAIN` is a read-only reference, so the C++
   compiler rejects every one of these assignments:
   ```
   190:14: error: assignment of read-only reference 'NPT2'
   ```

Per the guide's settled policy for this situation ("When the original
doesn't compile at all"), both fixes are syntax-only and non-numeric, and
were applied **directly to the delivered `xlh_mrgsolve_model_refactored.R`**
(not just a scratch copy), while `xlh_mrgsolve_model.R` was left completely
untouched:

1. The fourteen compartment names were removed from `$CAPTURE`. Every
   compartment's state is included in mrgsolve output regardless of
   `$CAPTURE` membership — confirmed by diffing `/model_manifest`'s
   `outputPaths` before and after: all nineteen compartments (nee
   fourteen-plus-renamed) are still present. This changes nothing about
   what is reported, only what compiles.
2. The twelve direct compartment assignments in `$MAIN` were switched to
   mrgsolve's `<cmt>_0` initial-value idiom (`NPT2_0 = NPT2_BASE;`, etc.),
   numerically identical to the original's intent (set the compartment's
   initial value once per individual, `NEWIND <= 1`).

Logged as entry **#58** in `translations/UPSTREAM_ISSUES.md`.

## qspserver `/model_manifest` discoverability

Confirmed via `POST /model_manifest` against the extracted, renamed DSL —
all Hill-interface parameters appear in the manifest's `parameters` with
their stated defaults: `EC50_BURO`, `GAMMA_BURO`, `EMAX_BURO`,
`EMAX_PHOSORAL`, `EC50_PHOSORAL`, `GAMMA_PHOSORAL`, `EMAX_CALC`,
`EC50_CALC`, `GAMMA_CALC`, plus every PK parameter
(`KA_BURO`/`CL_BURO`/`V1_BURO`/`V2_BURO`/`Q_BURO`/`F_BURO`,
`KA_PHOSORAL`/`KE_PHOSORAL`/`V_PHOSORAL`/`F_PHOSORAL`,
`KA_CALC`/`KE_CALC`/`V_CALC`/`F_CALC`).

`C_BURO`, `EFFECT_BURO`, `EFFECT_PHOSORAL`, and `EFFECT_CALC` are
state-derived doubles (recomputed every `$ODE` step, not fixed covariates),
so — per the same reasoning established in `mucopolysaccharidosis-type-1/
mps1_refactor_notes.md` and `primary-sclerosing-cholangitis/
psc_refactor_notes.md` for this engine/plugin combination — they cannot
also be `$PARAM` entries; they appear in the manifest's `outputPaths` via
`$CAPTURE` instead, confirmed present. `C_PHOSORAL` and `C_CALC` are
themselves `$CMT` compartments (see "bespoke" discussion above), so they
are *not* in `$CAPTURE` at all (that would re-trigger defect 1) but appear
automatically in every run's output as compartment states — also confirmed
via `/model_manifest`'s `outputPaths`.

`model_content` is pure DSL text extracted from the `code <- '...'`
R-string wrapper — no `.cpp` sibling was left behind; the extraction was
in-memory only, used to build the verification requests below and then
discarded.

## Verification

Per the guide's mandatory protocol: ran the original file's own dosing
scenarios through both the (defect-patched, for buildability only)
original and the refactored model via the qspserver `mrgsolve_api` service
(`POST /model_manifest`, `POST /run_simulation`), 1-year horizon, daily
(`delta=24`) output grid — matching the original scenarios' own
`run_scenario()` defaults (`end = 8760, delta = 24`). No shortening of the
window was needed; none of these scenarios hit the API's default solver
step-count limit.

Scenarios run (dose amounts/timing copied exactly from the original file's
own `scenarios` list, translated to the API's `dosing` convention;
compartment numbers are 1-based, unchanged by the rename since compartment
declaration order is unchanged):

1. **`1_Untreated_NaturalHistory`** — no dosing. Sanity check that all
   nineteen states hold flat at baseline with all three compounds silent.
2. **`2_Conventional_Ped_PhosCalcitriol`** — oral phosphate 180 mg (10
   mg/kg × 18 kg) into `GUT_PHOSORAL` (cmt 4), q6h × 1461 doses, **and**
   oral calcitriol 0.54 µg (30 ng/kg × 18 kg) into `GUT_CALC` (cmt 6), q24h
   × 365 doses — exercises **both** oral phosphate and oral calcitriol
   simultaneously, the original's own combination scenario.
3. **`5_Burosumab_Adult_1p0mgkg_Q4W`** — burosumab 70 mg (1.0 mg/kg × 70
   kg) into `GUT_BURO` (cmt 1), q672h × 13 doses — exercises burosumab
   alone.

Every shared output was compared point-by-point across the full time grid
for all three scenarios: `C_BURO`/`BURO_CP`, `EFFECT_BURO`/`FGF23_NEUT`,
`NPT2`, `TMPGFR`, `PHOS`, `CALCITRIOL`, `PTH`, `BSAP`, `RSS`, `HEIGHTZ_XLH`,
`AGV_CALC_XLH`, `SIXMWT`, `WOMAC`, `UCACR`, `NEPHROCALC`,
`C_PHOSORAL`/`PHOSORAL_SIG`, `C_CALC`/`CALC_CENT`.

**Result: exact match, max abs diff = 0.0 for every output, every
scenario**, both against the intermediate tested DSL and, as a final check,
against the DSL text extracted directly from the delivered
`xlh_mrgsolve_model_refactored.R` (the two differ only in added comments,
confirmed by diff before re-verifying). This is the expected outcome per
the guide's tolerance table for a pure structural reorganization with
rename-only Hill terms (no fitting performed, none needed).

Sanity check (not required for verification, but confirms the new
`EFFECT_PHOSORAL`/`EFFECT_CALC` captures are wired correctly): in scenario
2, `EFFECT_PHOSORAL` rises from 0 to ~0.133 and `PHOS` rises from 2.2 to
~2.59 mg/dL as the oral-phosphate/calcitriol regimen equilibrates; in
scenario 3, `EFFECT_BURO` rises to ~0.95, `NPT2` rises from 0.40 toward
~0.91, and `RSS` falls from 5.5 toward ~2.0 over the year — all directionally
consistent with the model's own calibration notes.

## Anything else worth flagging

- The three compounds are fully independent in this file — no shared
  multi-drug expression like `DRUG_ANTI_INFLAM` in the RA example. Each
  compound's `EFFECT_<STEM>` feeds a distinct disease target (`EFFECT_BURO`
  → `NPT2_TARGET` and `CALCITRIOL_TARGET`; `EFFECT_PHOSORAL` →
  `PHOS_TARGET`; `EFFECT_CALC` → the `CALCITRIOL` ODE's exogenous term and
  `UCACR_TARGET`), so there was no risk of accidentally collapsing separate
  drugs into one shared Hill term.
- Compartment ordering is unchanged (`GUT_BURO` is still compartment 1,
  etc.) since renaming was done in place without reordering — any external
  code addressing compartments by 1-based number, not name, is unaffected.
- `$MAIN`'s use of `NEWIND <= 1` to set initial conditions is otherwise
  untouched behaviorally by the `<cmt>_0` idiom switch (see build-defect
  fix #2 above) — same condition, same values, same timing, only the
  left-hand-side syntax changed to compile under this build.

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`, the
`x-linked-hypophosphatemia` rows for Burosumab, Oral calcitriol, and Oral
phosphate (PHOSORAL).

## Discoverability fix

A follow-up pass found a genuine naming-convention bug in the "bespoke"
compounds above, distinct from (and more invasive than) the `$GLOBAL`
forward-declare discoverability fix applied elsewhere in this batch.

**The bug.** The "Renaming" table above and the file's own top-of-file
comment describe oral phosphate's and oral calcitriol's second (non-depot)
compartment as renamed to `C_PHOSORAL`/`C_CALC` — but `C_<STEM>` is the slot
this guide's naming convention reserves for the *exposed-concentration
double* (`double C_<STEM> = <expr>;`), not for the `$CMT`-declared ODE state
itself. The earlier rename put the compartment in that slot directly
(`$CMT ... C_PHOSORAL ... C_CALC ...`, with `dxdt_C_PHOSORAL =`/`dxdt_C_CALC
=` as the governing ODEs), and there was no separate `double C_<STEM> =
<expr>;` statement anywhere in the file. The compound was completely valid,
compiling mrgsolve — just invisible to the corpus-wide discoverability scan
(pattern-matching `double C_<STEM> = <expr>;`) that the driver-PK dashboard
and other downstream tooling rely on to find a compound's exposed
concentration.

**The fix.** Renamed the compartment itself, and every reference to it as a
state, from `C_PHOSORAL` → `CENT_PHOSORAL` and `C_CALC` → `CENT_CALC`:
the `$CMT @annotated` declaration line, the `dxdt_CENT_PHOSORAL =`/
`dxdt_CENT_CALC =` equations (including the state read on their own
right-hand side), and the two places the state was read inside a Hill
formula in `$ODE` (`EFFECT_PHOSORAL`'s and `EFFECT_CALC`'s `pow(...)` calls,
now reading `CENT_PHOSORAL`/`CENT_CALC`). None of `KA_CALC`, `KE_CALC`,
`V_CALC`, `EC50_CALC`, `GAMMA_CALC`, `EMAX_CALC` (and the `PHOSORAL`
equivalents) were touched — those were already correctly named parameters,
not the compartment.

A new `$TABLE` block (none existed before) was added between `$ODE` and
`$CAPTURE`, containing the actual exposed-concentration lines:

```
double C_PHOSORAL = CENT_PHOSORAL;
double C_CALC     = CENT_CALC;
```

**Why this is a plain identity, not a mass/volume division** (checked, not
assumed): both ODEs are `dxdt_CENT_<STEM> = KA_<STEM>*GUT_<STEM>/V_<STEM> -
KE_<STEM>*CENT_<STEM>`. `GUT_<STEM>` is a mass (mg for phosphate, µg for
calcitriol) and `KA_<STEM>` is `1/h`, so `KA*GUT` is a mass flux (mass/h);
dividing that flux by `V_<STEM>` (L) gives `mass/(L·h)` — concentration per
unit time. That means `dxdt_CENT_<STEM>` is already in concentration/time
units, i.e. `CENT_<STEM>` itself integrates directly in concentration units
(matching the `$CMT` annotations: "Oral phosphate transient exposure signal
(a.u.)" and "Oral calcitriol central concentration (ng/mL equiv)" — the
compartment was already documented as a concentration, just misnamed). No
further division by a volume is needed or correct here, unlike Archetype
1/2/3's `C_<STEM> = CENT_<STEM> / V_<STEM>` (those apply when the
compartment holds an *amount*, which this bespoke pattern deliberately does
not — see the "Oral calcitriol and Oral phosphate — bespoke" discussion
above). `C_PHOSORAL = CENT_PHOSORAL` and `C_CALC = CENT_CALC` are therefore
exact identities.

`$CAPTURE` was extended from
`C_BURO EFFECT_BURO AGV_CALC_XLH EFFECT_PHOSORAL EFFECT_CALC` to also list
`C_PHOSORAL C_CALC` — per the guide's qspserver-compatibility rule 4, every
discoverable concentration must be captured, and neither was previously
(only the `EFFECT_` terms were). The pre-existing build-compat comment
immediately above `$CAPTURE` (about the original's 14-name `$CAPTURE`/`$CMT`
collision) was updated to note that the compartments it refers to are now
named `CENT_PHOSORAL`/`CENT_CALC`, so the new `C_PHOSORAL`/`C_CALC` capture
entries are the `$TABLE` doubles, not a re-introduction of that collision.

All other comments referencing the old `C_PHOSORAL`/`C_CALC` identifiers as
a *state* (the top-of-file refactor summary, the `$CMT` annotation trailing
comments, the bottom calibration notes) were updated to `CENT_PHOSORAL`/
`CENT_CALC`, or to both names where the comment now needs to distinguish the
state from the newly exposed concentration.

**Verification.** `Rscript -e 'parse(...)'` succeeds with zero errors after
the fix. `POST /model_manifest` against the extracted, edited DSL confirms
`C_PHOSORAL` and `C_CALC` now appear in `outputPaths` (they did not before).
Ran all three of the same dosing scenarios used for the original refactor's
verification — `1_Untreated_NaturalHistory` (no dosing),
`2_Conventional_Ped_PhosCalcitriol` (oral phosphate + oral calcitriol
together), `5_Burosumab_Adult_1p0mgkg_Q4W` (burosumab alone) — through both
the pre-edit DSL (`git show HEAD:.../xlh_mrgsolve_model_refactored.R`,
before this fix) and the edited DSL, 1-year horizon, `delta=24` daily grid,
requests spaced ~2s apart. Compared all 24 output columns each run (17
`$CMT` states unaffected by the rename, `C_BURO`, `EFFECT_BURO`,
`AGV_CALC_XLH`, `EFFECT_PHOSORAL`, `EFFECT_CALC`, plus the renamed pair
`C_PHOSORAL`→`CENT_PHOSORAL` and `C_CALC`→`CENT_CALC`) — **every column
matched exactly, max abs diff = 0.0, in all three scenarios** (72 column
comparisons total, all exact). Additionally confirmed, within the edited
model's own output, that `C_PHOSORAL == CENT_PHOSORAL` and `C_CALC ==
CENT_CALC` exactly at every timepoint in all three scenarios — direct
confirmation the new exposed-concentration identity is wired correctly, not
just that nothing else broke. This is a rename plus a pure-identity
addition, not a refit, and the verification result confirms no numeric
value changed anywhere in the model.
