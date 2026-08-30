# Refactor notes — `diabetic-nephropathy/dn_mrgsolve_model.R`

Scope: all **four** compounds this file models — **ACE inhibitor (ACEI)**,
**ARB (ARB)**, **SGLT2 inhibitor (SGLT)**, and **Finerenone (FINE)** — per the
four existing rows in `driver-patches/data/compound_perturbation_census.md`,
all classified "Redirect concentration (clean single site)". There is no
fifth pharmacological input in this file to leave untouched — the model has
exactly these four drug PK blocks.

## Archetype per compound

All four compounds share the same structure and were handled identically:

**Archetype 3 variant — depot + central, no peripheral, no bioavailability
term.** The original models each compound as a 1-compartment central volume
fed by a first-order absorption depot (`GI_<drug>  -> CENT_<drug>`), with
plain linear elimination `CL/V1`. There is no peripheral compartment and no
explicit oral bioavailability (`F`) parameter in the original for any of the
four drugs, so none was invented — per the guide's "never invent or default a
PK parameter the original didn't have," `F_<STEM>` is simply absent (fully
equivalent to `F=1`), matching the original's own math exactly.

| Role | ACEI | ARB | SGLT | FINE |
|---|---|---|---|---|
| Depot | `GUT_ACEI` (was `GI_acei`) | `GUT_ARB` (was `GI_arb`) | `GUT_SGLT` (was `GI_sglt2`) | `GUT_FINE` (was `GI_fine`) |
| Central | `CENT_ACEI` (was `CENT_acei`) | `CENT_ARB` (was `CENT_arb`) | `CENT_SGLT` (was `CENT_sglt2`) | `CENT_FINE` (was `CENT_fine`) |
| `KA_/CL_/V1_` | was `ka_acei/CL_acei/Vd_acei` | was `ka_arb/CL_arb/Vd_arb` | was `ka_sglt2/CL_sglt2/Vd_sglt2` | was `ka_fine/CL_fine/Vd_fine` |
| Exposed concentration | `C_ACEI` (was `Cp_acei`, ng/mL) | `C_ARB` (was `Cp_arb`, ng/mL) | `C_SGLT` (was `Cp_sglt2`, ng/mL) | `C_FINE` (was `Cp_fine`, ng/mL) |

All parameter values are copied verbatim from the original (`ka_acei=1.2` ->
`KA_ACEI=1.2`, `Vd_acei=25` -> `V1_ACEI=25`, `CL_acei=5.5` -> `CL_ACEI=5.5`,
etc. for all four drugs) — nothing invented, nothing defaulted, nothing
dropped. The original was already in CL/V convention (no `k10`/`k12`/`k21`
micro-constant form to convert), so the PK ODEs are a pure identifier rename
with identical arithmetic. Compartment order and numbering are unchanged
(`GUT_ACEI` is still compartment 1, `CENT_ACEI` compartment 2, ...,
`CENT_FINE` compartment 8), so the original's own 1-based dosing compartment
numbers (`cmt=1,3,5,7` in `build_events()`) still point at the same four
depots.

`C_<STEM>` keeps the original's own "rough unit" convention — the central
amount divided by volume, multiplied by 1000, labeled ng/mL in the original's
comments (`CENT_acei / Vd_acei * 1000`) — unchanged, just renamed.

## The Hill interface

Every one of the four original effect terms was already exactly
`Emax*C/(EC50+C)` — a plain ratio, no Hill exponent — so per the guide's
rename rule each got `GAMMA_<STEM> = 1.0` as an explicit named parameter (no
original Hill coefficient existed), and the formula was rewritten
`EMAX_<STEM>*pow(C_<STEM>,GAMMA_<STEM>)/(pow(EC50_<STEM>,GAMMA_<STEM>)+pow(C_<STEM>,GAMMA_<STEM>))`.
With `GAMMA=1`, `pow(x,1)==x` exactly, so this is arithmetically identical to
the original — confirmed by the exact-match verification below. This is a
**rename, not a refit**, for all four compounds.

| Compound | Was | Now |
|---|---|---|
| ACEI | `E_acei = Emax_acei*Cp_acei/(EC50_acei+Cp_acei)` | `EFFECT_ACEI` (`EMAX_ACEI`, `EC50_ACEI`, `GAMMA_ACEI=1`) |
| ARB | `E_arb = Emax_arb*Cp_arb/(EC50_arb+Cp_arb)` | `EFFECT_ARB` (`EMAX_ARB`, `EC50_ARB`, `GAMMA_ARB=1`) |
| SGLT | `E_sglt2 = Emax_sglt2*Cp_sglt2/(EC50_sglt2+Cp_sglt2)` | `EFFECT_SGLT` (`EMAX_SGLT`, `EC50_SGLT`, `GAMMA_SGLT=1`) |
| FINE | `E_fine = Emax_fine*Cp_fine/(EC50_fine+Cp_fine)` | `EFFECT_FINE` (`EMAX_FINE`, `EC50_FINE`, `GAMMA_FINE=1`) |

**Multi-drug combination points untouched.** The original already combined
ACEI and ARB only at the point the disease equations use them
(`RAAS_block = 1-(1-E_acei)*(1-E_arb)`); this was left exactly as-is,
substituting only the two renamed inputs
(`RAAS_block = 1-(1-EFFECT_ACEI)*(1-EFFECT_ARB)`). SGLT and FINE were never
combined multiplicatively with anything — each enters the disease ODEs
(`BG`, `TGF_cmpt`, `ECM_cmpt`, `Tub_cmpt`, `Fib_cmpt`, `GFR` hemodynamic
term, and the `$TABLE` clinical outputs `HbA1c`/`SBP`/`DBP`) as its own
independent `EFFECT_SGLT`/`EFFECT_FINE` term, matching the guide's "never
collapse several drugs into one shared Hill term" rule (already true of the
original, preserved as-is).

## qspserver compatibility

`C_<STEM>` and `EFFECT_<STEM>` are kept as plain `double` locals computed in
`$ODE` (not declared in `$PARAM`), matching the precedent set by every prior
refactor in this repo that hit the same mrgsolve constraint
(`$PARAM`/`$PARAM @annotated` values compile to read-only `const double&`
references, so a `$PARAM`-declared `C_ACEI` cannot be assigned inside
`$ODE`). Discoverability is satisfied the same way: every `C_<STEM>` and
`EFFECT_<STEM>` this file creates (`C_ACEI`, `C_ARB`, `C_SGLT`, `C_FINE`,
`EFFECT_ACEI`, `EFFECT_ARB`, `EFFECT_SGLT`, `EFFECT_FINE`) is listed in
`$CAPTURE`, confirmed present in `/model_manifest`'s `outputPaths`. The
`EMAX_*`/`EC50_*`/`GAMMA_*` Hill parameters are ordinary `$PARAM` values and
are fully overridable, confirmed present in `/model_manifest`'s `parameters`.

One build-specific finding worth flagging for future refactors in this repo:
**`$ODE`-local `double`s are directly visible inside `$TABLE`** in this
mrgsolve 2.0.1 build (they compile into the same scope) — an initial draft
of this file recomputed `C_<STEM>`/`EFFECT_<STEM>` a second time inside
`$TABLE` (matching what looked like a defensive pattern in an unrelated
prior refactor's notes) and hit a hard `redefinition of 'double
{anonymous}::C_ACEI'` C++ compile error. Fixed by computing each
`C_<STEM>`/`EFFECT_<STEM>` exactly once, in `$ODE`, and letting `$TABLE`/
`$CAPTURE` read that same value directly — no recomputation needed.

## When the original doesn't compile at all

The untouched original does not build under mrgsolve 2.0.1, for three
independent, pre-existing reasons unrelated to any of the four target
compounds — logged as **issue #61** in `translations/UPSTREAM_ISSUES.md`:

1. `$INIT` packs two `name=value` assignments per line with no comma
   separator (`GI_acei=0 CENT_acei=0`), which fails to parse as an R
   `list(...)` call before any C++ compilation is even attempted.
2. Once (1) is worked around, `$CMT` (bare) and `$INIT` are revealed to
   jointly redeclare all 19 compartment names — mrgsolve 2.0.1 treats
   `$INIT` as its own compartment-declaring block, not a companion to
   `$CMT` (same defect class as issues #34/#36).
3. `$CAPTURE` repeats 8 names already declared as `$CMT` compartments
   (`BG`, `AGE_cmpt`, `AngII_cmpt`, `TGF_cmpt`, `ROS_cmpt`, `ECM_cmpt`,
   `Pod_cmpt`, `Fib_cmpt`) — mrgsolve 2.0.1 refuses to build any model
   whose `$CAPTURE` repeats a compartment name (same defect class as
   issues #31/#34/#36).

Per the guide's settled policy for this situation, all three are fixed
**directly in the delivered `dn_mrgsolve_model_refactored.R`** (syntax-only,
non-numeric, disclosed here and in the manifest of changes above):

- The `$INIT` block was deleted; its 19 assignments were moved into a new
  `$MAIN` block using the modern `<CMT>_0 = value;` idiom
  (`GUT_ACEI_0=0;`, `BG_0=8.5;`, ... `GFR_cmpt_0=55;`) — no compartment,
  order, or starting value changed. This one change fixes defects (1) and
  (2) together.
- The 8 duplicated names were removed from `$CAPTURE`. Compartment states
  are always present in mrgsolve's own simulation output regardless of
  `$CAPTURE`, so `BG`, `AGE_cmpt`, `AngII_cmpt`, `TGF_cmpt`, `ROS_cmpt`,
  `ECM_cmpt`, `Pod_cmpt`, and `Fib_cmpt` are still available in every
  `as.data.frame(out)` the R script downstream produces — nothing in the
  plotting/summary code needed to change.

`dn_mrgsolve_model.R` itself is completely untouched and still carries all
three defects exactly as written.

## Verification

Per the guide's mandatory protocol: ran the original file's own dosing
scenarios through both a build-compat-patched scratch copy of the original
(in-memory only, never written to `dn_mrgsolve_model.R`) and the delivered
`dn_mrgsolve_model_refactored.R`, via the qspserver `mrgsolve_api` service
(`POST /model_manifest`, `POST /run_simulation`), comparing every shared
`$CAPTURE`d output point-by-point.

**Scenario S6 (Triple therapy — ACEi 10 mg q12h + SGLT2i 25 mg q24h +
Finerenone 20 mg q24h), full 2-year horizon (`end=730`, `delta=7` days,
matching the original's own `build_events("S6", 730)` and the R script's
own weekly `STEP=7` output grid exactly):** ran to completion on both
models (109 time points). All 23 shared disease/clinical outputs
(`BG, AGE_cmpt, AngII_cmpt, TGF_cmpt, ROS_cmpt, ECM_cmpt, Pod_cmpt,
UACR_cmpt, Tub_cmpt, Fib_cmpt, GFR_cmpt, Inh_ACE, Inh_AT1R, Inh_SGLT2,
Inh_MR, eGFR, UACR_o, HbA1c, SBP, DBP, CKD_Stage_num, pct_UACR, pct_GFR`)
and all 4 PK/exposure crosswalk pairs (`Cp_acei_out`<->`C_ACEI`,
`Cp_arb_out`<->`C_ARB`, `Cp_sglt2_out`<->`C_SGLT`, `Cp_fine_out`<->`C_FINE`)
matched with **max abs diff = 0.0** across every one of the 109 time
points.

**Scenario S2 (ARB monotherapy, 100 mg q24h) — isolates ARB, the one
compound not exercised by S6:** the original's own 730-day window exceeds
the qspserver API's default `maxsteps=20000` on both models identically
(`[lsoda] 20000 steps taken before reaching tout`) — a solver-budget limit,
not a mismatch (per the guide's step-count caveat). Shortened to a **60-day
window** (`end=60`, `delta=6`, `addl` scaled down to keep the same 100 mg
q24h regimen), which both models solved without error. Result: **max abs
diff = 0.0** across all 23 shared outputs and all 4 crosswalk pairs, over
all 12 time points (e.g. `Cp_arb_out`/`C_ARB` = 12.3424 ng/mL identically
at t=30d on both models).

**Result: exact match (max abs diff = 0.0) for every output, both
scenarios, everywhere both models could be run to completion.** This is
the expected outcome for a pure PK-structure rename plus `gamma=1` Hill
renames — there is no approximation anywhere in this file to introduce
error.

## Anything else worth flagging

- The four drugs' `Cp_<drug>_out` `$TABLE` variables in the original are
  exact duplicates, in both formula and value, of the `Cp_<drug>` doubles
  already computed in `$ODE` (e.g. `Cp_acei_out = CENT_acei/Vd_acei*1e3`
  vs. `Cp_acei = CENT_acei/Vd_acei*1000` — same quantity, recomputed under
  a different name purely for `$CAPTURE`). The refactor collapses this:
  `C_<STEM>` is computed once, in `$ODE`, and captured directly — see the
  "$ODE-local doubles are visible in $TABLE" finding above.
  `Inh_ACE`/`Inh_AT1R`/`Inh_SGLT2`/`Inh_MR` (legacy `%`-scaled output
  names) are preserved unchanged, now reading `EFFECT_<STEM>*100` instead
  of `E_<drug>*100`.
- No dead/unused parameters or dosing-switch flags were found in this file
  (unlike the AAA/RA precedents) — every `$PARAM` value in the original is
  read by some `$ODE`/`$TABLE` expression.
- The `HbA1c_base`/`kHbA1c` and `AGE_thresh`/`BG_min` parameters are
  declared but only `HbA1c_base` (not `kHbA1c`) and only the *value* of
  `BG_min` (not `AGE_thresh`) are actually read by the `$ODE`/`$TABLE`
  body — this predates the refactor and is not specific to any of the
  four target compounds, so it was left as-is (not logged separately;
  it's a much smaller, more common pattern of unused-parameter slack than
  the dedicated dead-parameter defects logged elsewhere in this repo).

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`, the
`diabetic-nephropathy | ACEI`, `diabetic-nephropathy | ARB`,
`diabetic-nephropathy | FINE`, and `diabetic-nephropathy | SGLT` rows.
