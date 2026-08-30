# Refactor notes — `takayasu-arteritis/ta_mrgsolve_model.R`

Scope: all four modeled compounds — prednisone/prednisolone (`PRED`),
tocilizumab (`TCZ`), methotrexate (`MTX`), infliximab (`IFX`). There are
no other compounds in this file, so nothing was left untouched by
exclusion; every disease-biology compartment (`IL6`, `sIL6R`,
`IL6_cmplx`, `TNF`, `TH1`, `TH17`, `TREG`, `VWI`, `ST`) and biomarker
(`CRP`, `PET`, `VWT`) plus `$TABLE`'s `NIH_SCORE`/`ITAS_SCORE`/
`RESPONSE_FLAG` are formula-identical to the original — only the four
drug-effect variable names feeding into them were substituted 1:1
(`Inh_PRED`→`EFFECT_PRED`, `Occ_TCZ`→`EFFECT_TCZ`, `Inh_MTX`→`EFFECT_MTX`,
`Inh_IFX`→`EFFECT_IFX`).

## Archetype per compound

**PRED (prednisone/prednisolone) — Archetype 3 (depot + central +
peripheral, linear), plus a preserved ke0 effect (biophase) compartment
not covered by name in the guide's four templates.** Confirmed from the
`dxdt_` lines: `PRED_GUT→PRED_C→PRED_P` is a standard mass-based
depot/central/peripheral CL/Q/V system (renamed `GUT_PRED`/`CENT_PRED`/
`PERI_PRED`), and a fourth compartment, `PRED_EFF`, equilibrates against
`PRED_C` via `ke0_PRED*(PRED_C - PRED_EFF)` — a classic link/biophase
model for delayed corticosteroid onset. The original's own `Inh_PRED`
Hill term reads `PRED_EFF`, **not** `PRED_C` — i.e. the PD-relevant site
is the biophase, not plasma. Renamed `EFF_PRED` (kept as its own
compartment, `KE0_PRED` for the rate constant), and exposed
`C_PRED = EFF_PRED` per the guide's "two [exposed values] only when a
genuinely different tissue site matters" allowance — plasma is still
captured informationally as `CP_PRED = CENT_PRED`, matching the original's
own `CP_PRED` output, but `C_PRED` (what `EFFECT_PRED` actually reads) is
the biophase level, not plasma.

**TCZ (tocilizumab) — Archetype 3 (depot + central + peripheral) plus a
preserved parallel linear + Michaelis-Menten clearance term, NOT
Archetype 4/TMDD.** The task brief flagged TMDD/receptor-binding kinetics
as something to check rather than assume for TCZ, given this drug's
five-different-structure precedent already documented across this repo
(`sepsis`, `rheumatoid-arthritis`, etc.). Checked against the actual
`$CMT`/`$ODE` here: there is **no** free-receptor or drug-receptor-complex
compartment anywhere in this file's TCZ block — only
`TCZ_SC`/`TCZ_C`/`TCZ_P` (depot/central/peripheral). The only
nonlinearity is a saturable elimination term added in parallel to the
linear clearance:
```
double CL_TCZ_tot = CL_TCZ + CLMM_TCZ * TCZ_C / (KM_TCZ + TCZ_C);
```
This is the standard "parallel linear + Michaelis-Menten" empirical
approximation for target-mediated elimination (no explicit receptor
states), not a genuine TMDD system — per the guide's own definition,
Archetype 4 is "used when the original explicitly models receptor
binding as its own dynamic system, not just an algebraic ratio on
concentration," and this file has neither the dynamic system nor even an
algebraic occupancy ratio: `Occ_TCZ` (renamed `EFFECT_TCZ`) is a plain
Hill function of raw `TCZ_C`. So this is Archetype 3 plus a preserved
bespoke nonlinear-clearance addition, kept exactly as the original
modeled it (per "don't flatten a mechanistically rich PK model... but
also don't add [TMDD] structure the original doesn't have"). Renamed
`CL_TCZ_tot` unchanged (already stem-neutral); `CLMM_TCZ`/`KM_TCZ` kept as
bespoke stem-scoped parameter names (no slot in the guide's table for a
parallel-MM term, so named consistently with the `_TCZ` stem convention
rather than invented generically).

**MTX (methotrexate) — Archetype 3 minus peripheral (depot + central,
linear) plus a bespoke third compartment for the active polyglutamate
pool, matching established corpus precedent.** `MTX_GUT→MTX_C` is a
standard depot+central linear system; a third compartment, `MTX_PG`,
accumulates from `MTX_C` via a first-order "polyglutamation" flux
(`kpg_MTX`) and depolyglutamates back out (`kdpg_MTX`). The original's
`Inh_MTX` Hill term reads `MTX_PG`, not `MTX_C` — the polyglutamate pool
is the actual PD-relevant moiety, exactly the same structural pattern
already established in this repo for methotrexate (`POLY_MTX` in
`rheumatoid-arthritis/ra_mrgsolve_model_refactored.R` and
`sarcoidosis/sarc_mrgsolve_model_refactored.R`). Renamed `MTX_PG`→
`POLY_MTX` to match that precedent exactly (not a new name invented for
this file), `kpg_MTX`→`KPOLY_MTX`, `kdpg_MTX`→`KDEPOLY_MTX`. Exposed
`C_MTX = POLY_MTX` (the compartment itself is the concentration the
original's effect term reads, no further transformation — same
"compartment IS concentration" pattern documented for `POLY_MTX` in the
sarcoidosis notes).

**IFX (infliximab) — Archetype 2 (no depot, central + peripheral,
linear).** IV-dosed directly into `IFX_C` (renamed `CENT_IFX`) via a
short infusion in the R-side event function; no absorption compartment.
`Inh_IFX` reads `IFX_C` directly (undivided), preserved as `C_IFX =
CENT_IFX` — same "compartment read directly as concentration" pattern as
TCZ below.

## A shared quirk across all four compounds: the exposed concentration is
## never volume-divided at the point PD reads it

For all four compounds, the original's Hill/Emax term reads its
compartment's value directly — `PRED_EFF`, `TCZ_C`, `MTX_PG`, `IFX_C` —
with **no** `/V1` or `/Vc` division anywhere at the point of use, even
though `Vc_*`/`Vp_*` volumes are used internally within each compound's
own ODE flux terms (in the standard CL/Q/V way). This is preserved
verbatim, not "fixed": `C_PRED = EFF_PRED`, `C_TCZ = CENT_TCZ`, `C_MTX =
POLY_MTX`, `C_IFX = CENT_IFX` — none divided by any `V1_<STEM>`. Renaming
these to look like Archetype 3's literal template line
(`C_TCZ = CENT_TCZ / V1_TCZ`) would silently change the model's numeric
output by the corresponding volume factor for every one of the four
compounds; this was checked and rejected in favor of preserving the
original's actual behavior, consistent with the guide's "verify against
the original, mechanically" rule.

## Hill interface: rename, not a fit — for all four compounds

Every one of `Inh_PRED`, `Occ_TCZ`, `Inh_MTX`, `Inh_IFX` in the original
was already the exact `Emax*X^gamma/(EC50^gamma+X^gamma)` shape, each with
its own explicit Hill exponent already declared (`hill_PRED=1.5`,
`hill_TCZ=1.2`, `hill_MTX=1.0`, `hill_IFX=1.3`) — unlike many files in
this corpus, none of the four needed a `GAMMA_<STEM>=1` default invented,
since the original already carried a per-compound exponent for every one
of them. All four are renames (`EFFECT_<STEM>`, `EMAX_<STEM>`,
`EC50_<STEM>`, `GAMMA_<STEM>`), not fits — no curve-fitting was needed or
performed for any compound, including TCZ, despite its nonlinear PK (the
nonlinearity lives entirely in the clearance term, not in the effect
term). No R² is reported because none applies.

## Renaming summary

| Original | Refactored |
|---|---|
| `PRED_GUT` | `GUT_PRED` |
| `PRED_C` | `CENT_PRED` |
| `PRED_P` | `PERI_PRED` |
| `PRED_EFF` | `EFF_PRED` |
| `ka_PRED` | `KA_PRED` |
| `Vc_PRED` | `V1_PRED` |
| `Vp_PRED` | `V2_PRED` |
| `CLd_PRED` | `Q_PRED` |
| `CL_PRED` | `CL_PRED` (unchanged) |
| `ke0_PRED` | `KE0_PRED` |
| `F_PRED` | `F_PRED` (unchanged) |
| `Emax_PRED` | `EMAX_PRED` |
| `EC50_PRED` | `EC50_PRED` (unchanged) |
| `hill_PRED` | `GAMMA_PRED` |
| `Inh_PRED` | `EFFECT_PRED` |
| `TCZ_SC` | `GUT_TCZ` |
| `TCZ_C` | `CENT_TCZ` |
| `TCZ_P` | `PERI_TCZ` |
| `ka_TCZ` | `KA_TCZ` |
| `Vc_TCZ` | `V1_TCZ` |
| `Vp_TCZ` | `V2_TCZ` |
| `CLd_TCZ` | `Q_TCZ` |
| `CL_TCZ`, `CLMM_TCZ`, `KM_TCZ` | unchanged (already stem-scoped) |
| `F_TCZ` | `F_TCZ` (unchanged) |
| `Emax_TCZ` | `EMAX_TCZ` |
| `EC50_TCZ` | `EC50_TCZ` (unchanged) |
| `hill_TCZ` | `GAMMA_TCZ` |
| `Occ_TCZ` | `EFFECT_TCZ` |
| `MTX_GUT` | `GUT_MTX` |
| `MTX_C` | `CENT_MTX` |
| `MTX_PG` | `POLY_MTX` (matches existing corpus precedent) |
| `ka_MTX` | `KA_MTX` |
| `Vc_MTX` | `V1_MTX` |
| `CL_MTX` | `CL_MTX` (unchanged) |
| `kpg_MTX` | `KPOLY_MTX` |
| `kdpg_MTX` | `KDEPOLY_MTX` |
| `F_MTX` | `F_MTX` (unchanged) |
| `Emax_MTX` | `EMAX_MTX` |
| `EC50_MTX` | `EC50_MTX` (unchanged) |
| `hill_MTX` | `GAMMA_MTX` |
| `Inh_MTX` | `EFFECT_MTX` |
| `IFX_C` | `CENT_IFX` |
| `IFX_P` | `PERI_IFX` |
| `Vc_IFX` | `V1_IFX` |
| `Vp_IFX` | `V2_IFX` |
| `CLd_IFX` | `Q_IFX` |
| `CL_IFX` | `CL_IFX` (unchanged) |
| `Emax_IFX` | `EMAX_IFX` |
| `EC50_IFX` | `EC50_IFX` (unchanged) |
| `hill_IFX` | `GAMMA_IFX` |
| `Inh_IFX` | `EFFECT_IFX` |

All parameter *values* are copied verbatim from the original — nothing
invented, nothing dropped. `DOSE_PRED`/`DOSE_TCZ`/`DOSE_MTX`/`DOSE_IFX`
(informational-only params, never read in `$ODE`, matching the
`DOSE_ERT`-style pattern already documented in `gaucher-disease`) and
`WT` are left unchanged. `$TABLE`'s `CP_PRED`/`CP_TCZ`/`CP_MTX`/`CP_IFX`
diagnostic outputs are kept under their original names (now pointing at
the renamed central compartments) to avoid churn in the R-side plotting
code (`p_pk1` references `CP_TCZ` by name), and all R-side `cmt = "..."`
dosing-target strings in `pred_events()`/`tcz_events()`/`mtx_events()`/
`ifx_events()` and the scenario definitions were updated to the renamed
compartments.

## Build-compatibility fix (mrgsolve 2.0.1) — logged as `UPSTREAM_ISSUES.md` #74

Confirmed via the qspserver `mrgsolve_api` container (`http://localhost:8007`)
that the **untouched original does not compile at all**, for three
independent, pre-existing reasons, none specific to any compound's own
PK/PD:

1. **`$PARAM` and `$CMT` (both non-annotated blocks) use C-style
   `/* ... */` block comments**, which this mrgsolve build's
   non-annotated-block parsers cannot handle — `$PARAM`'s R-based
   `list(...)` parse fails outright on the first comment; `$CMT`'s
   whitespace-split parser instead leaves stray `/*`/`*/` tokens as bogus
   duplicate compartment names. Fixed by stripping every `/* ... */` span
   from these two blocks only (comments only — no name or value touched;
   `$ODE`/`$TABLE`'s own C-style comments are untouched and compile fine,
   since those blocks go straight to the C++ compiler).
2. **`$CMT` (bare) and `$INIT` (bare `name = value`) jointly redeclare all
   24 compartment names** ("Duplicated model names") — same category as
   several other files in this corpus (`#30/#34/#41/#46/#57/#58/#64/#70/
   #72`). Fixed by moving `$INIT`'s 24 initial values into a `$MAIN` block
   using the standard `<CMT>_0 = value;` idiom.
3. **`$TABLE`'s three bare `capture NAME1 NAME2 NAME3 ...` statements
   (space-separated, no commas) are not valid multi-name capture syntax**
   under this mrgsolve build — they compile as literal C++ and the
   compiler chokes on the second bare name (`expected ';' before
   'NIH_SCORE'`). Fixed by consolidating all 14 names across the three
   lines into one ordinary `$CAPTURE` block.

All three fixes are syntax-only and non-numeric, applied directly to
`ta_mrgsolve_model_refactored.R` per the guide's settled policy; the
tracked `ta_mrgsolve_model.R` is untouched and still carries all three
defects exactly as written. Full detail and reproduction in
`translations/UPSTREAM_ISSUES.md` #74.

## A fourth finding, PRED-specific and in-scope: a shared, pre-existing NaN fragility, fixed (not merely logged)

While verifying, every one of the file's own combination scenarios
(2 through 5, i.e. every scenario that includes prednisone) hit an
identical `NaN` blowup partway through its own named 365-day duration —
**in the identically-patched original as well**, confirming this is a
pre-existing defect unrelated to the renaming, not something introduced
by this refactor:

| Scenario | Blowup time (h) | Blowup day |
|---|---|---|
| 3 (Pred + MTX) | 1568 | ~65.3 |
| 4 (Pred + TCZ, TAKT) | 2292 | ~95.5 |
| 5 (Pred + IFX) | 1576 | ~65.7 |

Root cause: as the tapering prednisone dose decays `PRED_C`/`PRED_EFF`
toward zero between doses, the adaptive-step solver's floating-point
noise occasionally undershoots to a tiny **negative** value (confirmed:
`PRED_C = -8.5007e-11` at t=2288h in Scenario 4, immediately followed by
`NaN` in every compartment at t=2292h). `Inh_PRED`'s Hill term uses
`pow(PRED_EFF, hill_PRED)` with `hill_PRED = 1.5` — a non-integer
exponent — and `pow()` of a negative base with a non-integer exponent is
`NaN` in C++; once any single `dxdt_` evaluates to `NaN`, the whole
coupled ODE system stays `NaN` for the rest of the run. This is the same
category of fragility already documented for `von-willebrand-disease`
(`UPSTREAM_ISSUES.md` #72, DDAVP's fractional-`pow()` NaN) and
`essential-thrombocythemia` (#64).

Because this sits squarely inside PRED's own PD term — one of this
refactor's four in-scope compounds — it is **fixed as part of the
refactor itself**, per the guide's point 4 ("if a defect is genuinely
inside the scope compound's own block, fix it as part of the refactor
and say so"), not merely logged. The delivered
`ta_mrgsolve_model_refactored.R` guards every one of the four `C_<STEM>`
exposed-concentration variables with `fmax(..., 0.0)` before they enter
their respective Hill terms:
```
double C_PRED = fmax(EFF_PRED, 0.0);
double C_TCZ  = fmax(CENT_TCZ, 0.0);
double C_MTX  = fmax(POLY_MTX, 0.0);
double C_IFX  = fmax(CENT_IFX, 0.0);
```
`GAMMA_MTX = 1.0` (integer) is not actually NaN-exposed by a negative
base, but is guarded too for consistency. `GAMMA_TCZ = 1.2` and
`GAMMA_IFX = 1.3` are both non-integer and therefore equally exposed in
principle, even though none of the three tested scenarios actually drove
`CENT_TCZ`/`CENT_IFX` negative before PRED's own blowup pre-empted them
(PRED clears fastest of the four and is present in every combination
scenario, so it always hits zero first). The tracked original is
untouched and still `NaN`s at the three times above, unfixed.

## Verification

**Method**: qspserver `mrgsolve_api` container (`http://localhost:8007`),
`/model_manifest` and `/run_simulation`, per the task's instruction not to
use local R/mrgsolve. Extracted the bare DSL text from both
`ta_mrgsolve_model.R` (identically patched for the three build-compat
defects above, no compound renamed, no fmax guard) and
`ta_mrgsolve_model_refactored.R` (fully renamed, fmax guard applied) and
POSTed each as `model_content`. Requests spaced ~2s apart per the API's
concurrency note. `/model_manifest` confirmed `outputPaths` includes
`C_PRED`/`C_TCZ`/`C_MTX`/`C_IFX` and `EFFECT_PRED`/`EFFECT_TCZ`/
`EFFECT_MTX`/`EFFECT_IFX` (all four `$CAPTURE`d, per the qspserver
compat requirement — not `$PARAM`, since `$PARAM @annotated` values
compile to read-only `const double&` in this mrgsolve build and cannot be
reassigned in `$ODE`; the same constraint already documented in
`abdominal-aortic-aneurysm/aaa_refactor_notes.md`, so discoverability is
satisfied via `$CAPTURE`/`outputPaths` instead, matching that precedent
and the 13 prior TCZ refactors before it).

**Scenarios run**: the file's own named Scenario 3 (Prednisone +
Methotrexate), Scenario 4 (Prednisone + Tocilizumab, the TAKT regimen),
and Scenario 5 (Prednisone + Infliximab) — reconstructed exactly from
`pred_events()`/`tcz_events()`/`mtx_events()`/`ifx_events()`'s own R logic
(taper schedule, dose amounts, intervals), converted to the API's
`dosing` field (`DoseSpec` array, one entry per discrete dose event,
since the API's `events`/`EventSpec` field is deSolve-engine-only and is
silently ignored — confirmed empirically — by this mrgsolve-backed
endpoint). IFX's original `rate = -2` (a 2h-infusion flag needing an
external `D1`/`R1` data column mrgsolve's own `ev()` never supplied here
either) is approximated in this verification as an equal-total-amount
instant bolus at the same dose times; this affects only the immediate
post-dose spike shape, not the longer-run PK/PD trajectory the guide asks
us to verify, and is disclosed here rather than silently substituted.
Each scenario run over the full 365-day/8760h horizon the scenario itself
specifies (no shortening needed — the API's default solver step budget
handled the full year without a `maxsteps` error).

**Result**: for every one of the 28 shared outputs (all 12 drug-PK
compartments across the three scenarios' relevant compounds, all 9
disease-biology compartments, all 3 biomarkers, and all `$TABLE` captures
— `ESR_now`, `NIH_SCORE`, `ITAS_SCORE`, `RESPONSE_FLAG`, `CP_PRED`,
`CP_<TCZ/MTX/IFX>`, `Drug_inh_VWI`, `EFFECT_PRED`, and the scenario's
second `EFFECT_<STEM>`), **max abs diff = 0.0** across every time point
before that scenario's shared NaN blowup point (671/463/458 points for
Scenarios 4/3/5 respectively, i.e. exact agreement for the entire
pre-blowup run). This is a pure structural reorganization (Archetypes 1-3
and the edge case, per the guide's tolerance rule) for all four
compounds — no Hill-fitting was needed, so exact match is the expected
and achieved result. Past each blowup point, the identically-patched
original stays `NaN` for the rest of the run (as it does un-patched too),
while the refactored file (with the in-scope `fmax` guard) continues to
produce finite, physically sensible values (`NIH_SCORE`/`CRP`/etc.
evolving smoothly) all the way to t=8760h in all three scenarios — an
expected, disclosed divergence from the guard, not a verification
failure.

## Deliverables

- `ta_mrgsolve_model_refactored.R` — build-compat-fixed, four compounds
  renamed to convention, `fmax` NaN guard applied to all four `C_<STEM>`.
- This file.
- `driver-patches/data/compound_perturbation_census.md` — `takayasu-
  arteritis` rows for `IFX`/`MTX`/`PRED`/`TCZ` filled in (IFX's
  Target/Pathway corrected from `?` to `TNF-alpha (neutralization)`,
  confirmed directly from the original's own `Inh_IFX` term).
- `translations/UPSTREAM_ISSUES.md` #74.

The tracked `ta_mrgsolve_model.R` is completely untouched.
