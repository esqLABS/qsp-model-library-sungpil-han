# Refactor notes — `peripheral-arterial-disease/pad_mrgsolve_model.R`

Scope: **all six modeled compounds** — Clopidogrel (parent) + its own active
metabolite CLOPAM, Aspirin (ASP), Ticagrelor (TICA), Rivaroxaban (RIVA),
Cilostazol (CILO), Atorvastatin (ATST) — the model has no other drugs. All
nine disease-PD compartments (`Plt_agg`, `Thrombin`, `LDL_C`, `ABI`,
`WalkDist`, `Collat`, `EF_idx`, `PlaqueVol`, `hsCRP`) and their kinetics are
untouched except for one disclosed build-compatibility fix (see below),
which touches only the `$INIT`→`$MAIN` initial-condition idiom, not any
dynamics.

## Clopidogrel's active metabolite (census row "Clopidogrel (AM)") is a genuine parent-metabolite chain, not a second independently-dosed compound

The task brief flagged this as something to check rather than assume, per
the sarcoidosis PRED/PREDL and urolithiasis ALLO/OXP precedent. The code
confirms a genuine parent-metabolite chain, same pattern:

- `C_am` (the active thiol metabolite) has **no `GUT_`/depot compartment and
  no `ka_`/`F_` absorption parameters of its own** anywhere in the file, and
  no dose event targets it in any of the original's own seven scenarios.
- `C_am`'s only inflow term is a conversion computed directly from the
  parent's own plasma concentration:
  ```
  dxdt_C_clopi  = ka_clopi * A_clopi / Vc_clopi
                  - (CL_clopi / Vc_clopi) * C_clopi
                  - kact * C_clopi;
  dxdt_C_am     = kact * C_clopi * Vc_clopi / Vc_am
                  - (CL_am / Vc_am) * C_am;
  ```
  `kact` is literally named "CYP2C19 activation rate to thiol metabolite" in
  its own `$PARAM` comment, and is subtracted from `C_clopi`'s own mass
  balance (real mass diversion, not just an informational read).
- The P2Y12 disease-effect term (`Inh_AM`) reads only `C_am`, **never**
  `C_clopi` — parent clopidogrel itself has no direct pharmacodynamic effect
  in this model; only its active metabolite does.

**Handled as a genuine parent-metabolite chain:** clopidogrel keeps its own
PK block (`GUT_CLOP`/`CENT_CLOP`) which exists only to drive the conversion
into the metabolite (no `EFFECT_CLOP`, matching the original having no
direct clopidogrel-parent PD effect). The metabolite keeps its own
single-compartment PK (`CENT_CLOPAM`, fed by that conversion) and exposes
the `C_CLOPAM`/`EFFECT_CLOPAM` pair the disease equation (P2Y12 inhibition,
combined with ticagrelor's own effect) actually reads. The census row for
"Clopidogrel (AM)" (previously `Target/Pathway: ?`) is corrected below to
reflect the parent-metabolite relationship and the real target (P2Y12
receptor), not an independent compound with an unknown target.

## Every compartment already IS its own concentration — no `C_<STEM> = CENT_<STEM>/V1` division needed

Unlike a typical mass-in-a-volume compartment, every one of this file's six
PK blocks parameterizes its central compartment directly in concentration
units. Confirmed from the actual `dxdt_` lines — e.g. aspirin:
```
dxdt_C_asp = ka_asp * A_asp / Vc_asp - (CL_asp / Vc_asp) * C_asp;
```
This is dimensionally the derivative of `A_asp/Vc_asp` (chain rule, `Vc_asp`
constant), so `C_asp` integrates concentration directly — the same
"compartment IS concentration" pattern already documented for `HTBZ_brain`/
`POLY_MTX` in the huntingtons-disease and sarcoidosis refactor notes. All
six compounds (and the metabolite CLOPAM) follow this pattern identically.
Renamed `CENT_<STEM>` per the naming convention's "central compartment"
role, with a pass-through `double C_<STEM> = CENT_<STEM>;` exposing the
convention's own concentration name — same structure the huntingtons-disease
notes use for `C_HTBZ = CENT_HTBZ`.

## Archetype per compound

**Clopidogrel (parent) — Archetype 3 minus the peripheral compartment**
(depot + central, linear elimination) **plus a first-order conversion
pathway out of the central compartment** that forms CLOPAM — not flattened
to a plain archetype, since the conversion is a real part of the original's
own pharmacology (same "don't simplify away a mechanistically richer term"
principle the guide applies to TMDD):
```
dxdt_GUT_CLOP  = -KA_CLOP * GUT_CLOP;
dxdt_CENT_CLOP = KA_CLOP*GUT_CLOP/V1_CLOP - (CL_CLOP/V1_CLOP)*CENT_CLOP
                 - KCONV_CLOPAM*CENT_CLOP;
```
No direct disease effect (`EFFECT_CLOP` does not exist), matching the
original having no direct parent-clopidogrel PD effect.

**Clopidogrel active metabolite (CLOPAM) — Archetype 1** (no depot, single
compartment, linear elimination), fed by the conversion inflow above:
```
dxdt_CENT_CLOPAM = KCONV_CLOPAM*CENT_CLOP*V1_CLOP/V1_CLOPAM
                   - (CL_CLOPAM/V1_CLOPAM)*CENT_CLOPAM;
```

**Aspirin (ASP) — Archetype 3 minus the peripheral compartment** (depot +
central, linear elimination):
```
dxdt_GUT_ASP  = -KA_ASP * GUT_ASP;
dxdt_CENT_ASP = KA_ASP*GUT_ASP/V1_ASP - (CL_ASP/V1_ASP)*CENT_ASP;
```

**Ticagrelor (TICA) — Archetype 3 minus the peripheral compartment** (depot
+ central, linear elimination): same structure as aspirin.

**Rivaroxaban (RIVA) — Archetype 3 minus the peripheral compartment** (depot
+ central, linear elimination): same structure.

**Cilostazol (CILO) — Archetype 3 minus the peripheral compartment** (depot
+ central, linear elimination): same structure. `Emax_walk` (a downstream
disease-side scaling factor on the walking-distance equation, not a Hill
`Emax`) is kept as a separate, unrenamed parameter — it multiplies
`EFFECT_CILO` inside the disease equation rather than being part of
cilostazol's own Hill interface, so it does not fit the `EMAX_<STEM>` slot.

**Atorvastatin (ATST) — Archetype 3 minus the peripheral compartment**
(depot + central, linear elimination): same structure.

None of the six compounds have a peripheral compartment or TMDD in the
original — every PK block is exactly depot+central, so no archetype had to
be forced or flattened.

## Renaming

| Original | Refactored |
|---|---|
| `ka_clopi`, `CL_clopi`, `Vc_clopi` | `KA_CLOP`, `CL_CLOP`, `V1_CLOP` |
| `A_clopi`, `C_clopi` | `GUT_CLOP`, `CENT_CLOP` |
| `kact` | `KCONV_CLOPAM` (renamed for its actual role — a first-order conversion rate constant forming the metabolite — same naming logic as `k_Allo_Oxp`→`KCONV_OXP` in `urolithiasis/uri_refactor_notes.md`) |
| `CL_am`, `Vc_am`, `C_am` | `CL_CLOPAM`, `V1_CLOPAM`, `CENT_CLOPAM` |
| `EC50_P2Y12`, `Emax_P2Y12` | `EC50_CLOPAM`, `EMAX_CLOPAM` |
| — | `GAMMA_CLOPAM = 1.0` (new; implicit in the original's plain-ratio shape) |
| `ka_asp`, `CL_asp`, `Vc_asp` | `KA_ASP`, `CL_ASP`, `V1_ASP` |
| `A_asp`, `C_asp` | `GUT_ASP`, `CENT_ASP` |
| `EC50_COX1`, `Emax_COX1` | `EC50_ASP`, `EMAX_ASP` |
| — | `GAMMA_ASP = 1.0` (new) |
| `ka_tica`, `CL_tica`, `Vc_tica` | `KA_TICA`, `CL_TICA`, `V1_TICA` |
| `A_tica`, `C_tica` | `GUT_TICA`, `CENT_TICA` |
| `EC50_tica`, `Emax_tica` | `EC50_TICA`, `EMAX_TICA` |
| — | `GAMMA_TICA = 1.0` (new) |
| `ka_riva`, `CL_riva`, `Vc_riva` | `KA_RIVA`, `CL_RIVA`, `V1_RIVA` |
| `A_riva`, `C_riva` | `GUT_RIVA`, `CENT_RIVA` |
| `EC50_FXa`, `Emax_FXa` | `EC50_RIVA`, `EMAX_RIVA` |
| — | `GAMMA_RIVA = 1.0` (new) |
| `ka_cilo`, `CL_cilo`, `Vc_cilo` | `KA_CILO`, `CL_CILO`, `V1_CILO` |
| `A_cilo`, `C_cilo` | `GUT_CILO`, `CENT_CILO` |
| `EC50_PDE3`, `Emax_PDE3` | `EC50_CILO`, `EMAX_CILO` |
| — | `GAMMA_CILO = 1.0` (new) |
| `Emax_walk` | kept unchanged (disease-side scaling on `EFFECT_CILO`, not part of the Hill interface — see above) |
| `ka_atst`, `CL_atst`, `Vc_atst` | `KA_ATST`, `CL_ATST`, `V1_ATST` |
| `A_atst`, `C_atst` | `GUT_ATST`, `CENT_ATST` |
| `EC50_HMG`, `Emax_HMG` | `EC50_ATST`, `EMAX_ATST` |
| — | `GAMMA_ATST = 1.0` (new) |
| `Inh_AM` | `EFFECT_CLOPAM` |
| `Inh_tica` | `EFFECT_TICA` |
| `Inh_P2Y12` (`fmax(Inh_AM, Inh_tica)`) | `EFFECT_P2Y12` (`fmax(EFFECT_CLOPAM, EFFECT_TICA)`, kept as the same combination formula, computed only at the point of use) |
| `Inh_COX1` | `EFFECT_ASP` |
| `Inh_FXa` | `EFFECT_RIVA` |
| `Eff_PDE3` | `EFFECT_CILO` |
| `Eff_HMG` | `EFFECT_ATST` |

All parameter *values* are copied verbatim from the original — nothing
invented except the six `GAMMA_*` = 1.0 additions (each a rename target the
original's bare Emax-ratio shape mathematically already implied).

## Hill interface: rename, not a fit — for all six compounds

Every one of `Inh_AM`, `Inh_tica`, `Inh_COX1`, `Inh_FXa`, `Eff_PDE3`, and
`Eff_HMG` in the original was already the exact plain-ratio shape
`Emax*C/(EC50+C)` (implicit `gamma=1`):
```
double Inh_AM   = Emax_P2Y12 * C_am / (EC50_P2Y12 + C_am);
double Inh_tica = Emax_tica  * C_tica / (EC50_tica  + C_tica);
double Inh_COX1 = Emax_COX1 * C_asp / (EC50_COX1 + C_asp);
double Inh_FXa  = Emax_FXa  * C_riva / (EC50_FXa  + C_riva);
double Eff_PDE3 = Emax_PDE3 * C_cilo / (EC50_PDE3 + C_cilo);
double Eff_HMG  = Emax_HMG  * C_atst / (EC50_HMG  + C_atst);
```
Rewritten to the explicit `pow(C, GAMMA)/(pow(EC50,GAMMA)+pow(C,GAMMA))` form
with `GAMMA_*=1.0` — mathematically identical (`pow(x,1)=x`), a pure rename
with no curve fitting performed or needed. No R² is reported for any
compound because none applies.

The combined P2Y12 term (`fmax` of the clopidogrel-metabolite and ticagrelor
effects) is preserved exactly as the original's own combination rule — per
the guide's "combine each compound's `EFFECT_<STEM>` only at the point of
use" instruction, `EFFECT_CLOPAM` and `EFFECT_TICA` are kept as two
independently-driveable terms, combined into `EFFECT_P2Y12` only where the
platelet-aggregation and plaque-progression equations actually consume it
(same `fmax`, not rewritten to Bliss independence or anything else, since
the original's own combination rule is `fmax`, not a Bliss/multiplicative
form).

## `$CAPTURE` naming for the exposed concentrations and effect terms

`C_CLOP`, `C_CLOPAM`, `C_ASP`, `C_TICA`, `C_RIVA`, `C_CILO`, `C_ATST`, and
all six `EFFECT_<STEM>` plus the combined `EFFECT_P2Y12` are computed as
`$ODE`-local `double`s (needed there for the model's own dynamics), so they
are captured under an `_OUT` suffix (`C_CLOP_OUT`, `EFFECT_CLOPAM_OUT`, …)
rather than a bare self-referential `capture C_CLOP = C_CLOP;`, which would
be a `redefinition of capture` under mrgsolve 2.0.1 — same pattern as the
sarcoidosis/urolithiasis/huntingtons-disease precedents. Confirmed
discoverable via `/model_manifest`'s `outputPaths`. All of the original's
own pre-existing `$TABLE` captures (`Clopi_ng`, `AM_ng`, `Tica_ng`,
`Riva_ng`, `Cilo_ng`, `Atst_ng`, `Inh_P2Y12_pct`, `Inh_COX1_pct`,
`Inh_FXa_pct`, `Eff_PDE3_pct`, `LDL_reduction`, `ABI_out`, `Walk_m`,
`Plt_pct`, `Thrombin_out`, `Collat_idx`, `EF_pct`, `LDL_out`, `CRP_out`,
`Plaque_out`, `MACE_risk`, `MALE_risk`, `Rutherford`) are kept under their
original names — every one of these was already named independently of the
renamed compound internals, so the post-DSL R script's `outvars=` list and
all plotting code needed **no changes** to keep working.

## Build-compat fix (mrgsolve 2.0.1) — logged as `UPSTREAM_ISSUES.md` #96

Confirmed via the qspserver `mrgsolve_api` container (`http://localhost:8007`)
that the **untouched original does not compile at all**:

```
Error in validObject(.Object) :
  invalid class "mrgmod" object: Duplicated model names: A_clopi C_clopi
  C_am A_asp C_asp A_tica C_tica A_riva C_riva A_cilo C_cilo A_atst C_atst
  Plt_agg Thrombin LDL_C ABI WalkDist Collat EF_idx PlaqueVol hsCRP
```

Same defect class as issues #29/#34/#36/#59/#76/#81/#82/#83: `$CMT`+`$INIT`
jointly redeclare all 20 compartments. Not specific to any of the six
compounds — it's in the shared `$INIT` boilerplate.

**Fix applied directly to the delivered `pad_mrgsolve_model_refactored.R`**
(not just a scratch copy), per the guide's settled policy, syntax-only and
non-numeric: the `$INIT` block deleted, its 20 assignments moved into
`$MAIN` as `<CMT>_0 = value;` under the refactor's renamed compartments
(e.g. `GUT_CLOP_0 = 0;`, `CENT_CLOPAM_0 = 0;`), same values, no compartment
added or removed. The checked-in original (`pad_mrgsolve_model.R`) is
untouched and still carries the defect exactly as written; an identically
syntax-fixed scratch copy of the original (same fix, original
compartment/param names, never committed) was built in-memory purely to
construct the verification comparison target. Full disclosure and the exact
error trace are in `translations/UPSTREAM_ISSUES.md` #96.

## A second, independent defect found and logged, not fixed — `UPSTREAM_ISSUES.md` #97

While extending the verification window, the original's own LDL-C turnover
parameters were found to be grossly imbalanced: `k_ldl_prod * LDL_base`
(production, `0.312 mg/dL/h`) outweighs `k_ldl_cl * LDL_C` (clearance,
`0.0024 mg/dL/h` at the stated baseline) by roughly 130-fold, so `LDL_C`
climbs almost linearly rather than sitting near its own stated 130 mg/dL
baseline — reaching ~1005 mg/dL by day 121 and ~3627 mg/dL by day 528, at
which point the solver fails (`lsoda`, convergence error / `NaN`), well
inside the model's own stated 2-year default horizon and its own 180/365-day
summary-table checkpoints. `k_ldl_cl`'s own comment (`t1/2 ~55h`) implies a
value roughly 680x larger than what is actually declared — consistent with
a units/typo defect. This is a disease-PD-side parameter defect, unrelated
to any of the six compounds' own PK/PD (statin's `Eff_HMG`/`EFFECT_ATST`
can at most quadruple the already-too-small clearance term, nowhere near
enough to rescue the imbalance). **Not fixed** — `k_ldl_cl`'s value is
copied verbatim into the refactor per the never-invent-or-correct-a-
parameter-value rule; full detail in `translations/UPSTREAM_ISSUES.md` #97.

This defect is why the verification window below is shortened to 60 days
rather than treating a longer-window mismatch/solver failure as a refactor
bug — the original itself becomes numerically unstable well before either
the 2-year default or even the 365-day summary-table checkpoint, for every
scenario tested (with or without statin).

## qspserver compatibility checklist

- `model_content` is pure mrgsolve DSL text (no R wrapper) — confirmed by
  extracting the quoted `pad_code <- '...'` block byte-for-byte and building
  it standalone via `POST /model_manifest`; the extracted text is
  byte-identical (15,506 characters) to the block embedded in
  `pad_mrgsolve_model_refactored.R`.
- Every `KA_`/`CL_`/`V1_`/`KCONV_`/`EMAX_`/`EC50_`/`GAMMA_` parameter for all
  six compounds (plus CLOPAM) lives in `$PARAM` (confirmed present in the
  manifest's `parameters` list, 60 total).
- `C_CLOP_OUT`/`C_CLOPAM_OUT`/`C_ASP_OUT`/`C_TICA_OUT`/`C_RIVA_OUT`/
  `C_CILO_OUT`/`C_ATST_OUT`/`EFFECT_CLOPAM_OUT`/`EFFECT_TICA_OUT`/
  `EFFECT_ASP_OUT`/`EFFECT_RIVA_OUT`/`EFFECT_CILO_OUT`/`EFFECT_ATST_OUT`/
  `EFFECT_P2Y12_OUT` are all discoverable `outputPaths` (58 total outputs,
  see manifest above).
- `$SET end/delta` is not present in the original either — the R-side
  `run_scenario()` wrapper supplies `end`/`tgrid` explicitly; the
  verification below drives `time.end`/`time.delta` through the request, not
  through any model-internal default.
- No R-only syntax inside the DSL block — the extracted text compiled
  standalone with no surrounding R script.
- Solver step budget: not the limiting factor here — the limiting factor is
  the original's own numerical instability (issue #97 above), which forces
  a shorter window than the default 20000-maxsteps budget would otherwise
  allow.

## Verification

Per the guide's mandatory protocol: extracted the quoted DSL block from both
the build-compat-patched original (scratch copy only, never committed) and
the delivered `pad_mrgsolve_model_refactored.R`, confirmed both build via
`POST /model_manifest` on the qspserver `mrgsolve_api` container
(`http://localhost:8007`, confirmed healthy throughout, requests spaced ~2s
apart), then ran the original file's own seven named dosing scenarios
(`scenarios` list) through `POST /run_simulation`:

| Scenario | Dosing (from the original's own `scenarios` list) | Result |
|---|---|---|
| 1_Untreated | no dosing | exact match, max abs diff 0.0 |
| 2_Aspirin_100 | `GUT_ASP` (cmt 4): 100mg q24h x 59 addl | exact match, max abs diff 0.0 |
| 3_Clopidogrel_75 | `GUT_CLOP` (cmt 1): 75mg q24h x 59 addl | exact match, max abs diff 0.0 |
| 4_DAPT_Clopi_ASA | `GUT_CLOP` + `GUT_ASP` combined | exact match, max abs diff 0.0 |
| 5_COMPASS | `GUT_RIVA` (cmt 8): 2.5mg q12h + `GUT_ASP` 100mg q24h | exact match, max abs diff 0.0 |
| 6_Cilostazol_ASA | `GUT_CILO` (cmt 10): 100mg q12h + `GUT_ASP` 100mg q24h | exact match, max abs diff 0.0 |
| 7_Optimal | `GUT_CLOP`+`GUT_ASP`+`GUT_RIVA`+`GUT_ATST` (cmt 12) combined | exact match, max abs diff 0.0 |

**Verification window: 60 days at `delta=6h` (241 points), not the
original's own 730-day/2-year default** — per the guide's provision to
shorten the window rather than treat a solver failure as a mismatch, and
per issue #97 above: the original's own LDL-C/plaque runaway causes `lsoda`
convergence failures or `maxsteps` exhaustion starting around day 142-146
for several scenarios (confirmed by first attempting the full 365-day
window and observing `Aspirin_100`, `DAPT_Clopi_ASA`, `COMPASS`, and
`Cilostazol_ASA` all fail with `lsoda`/`istate -5` convergence errors around
t≈3400-3650h, and `Cilostazol_ASA` separately exhausting the API's default
20000-maxsteps budget) — 60 days is safely inside the stable regime for
every scenario tested (confirmed: `ABI_out` has already declined from its
0.70 baseline to ~0.43 by day 60, well along the same trajectory that
eventually diverges, but the ODE system itself is still numerically stable
there).

Every shared output was compared point-by-point across the full 241-point
time grid: all 20 raw compartments (mapped through the renames given above)
and all 23 pre-existing `$TABLE` captures (`Clopi_ng`, `AM_ng`, `Tica_ng`,
`Riva_ng`, `Cilo_ng`, `Atst_ng`, `Inh_P2Y12_pct`, `Inh_COX1_pct`,
`Inh_FXa_pct`, `Eff_PDE3_pct`, `LDL_reduction`, `ABI_out`, `Walk_m`,
`Plt_pct`, `Thrombin_out`, `Collat_idx`, `EF_pct`, `LDL_out`, `CRP_out`,
`Plaque_out`, `MACE_risk`, `MALE_risk`, `Rutherford`). **Result: exact
match, max abs diff = 0.0, for every output in all seven scenarios** — the
expected result for pure structural reorganization with no Hill-refitting
(every one of the six compounds' effect terms was already an exact
Hill/plain-ratio shape in the original, see Archetype section above).

Sanity-checked the new `_OUT` captures are non-trivial (not degenerately
zero): in a bespoke Clopidogrel+Rivaroxaban run (60 days), `C_CLOPAM_OUT`
peaks at 0.0004 mg/L and `EFFECT_CLOPAM_OUT`/`EFFECT_P2Y12_OUT` peak at
0.0552 (5.5% P2Y12 inhibition) — consistent with the original's own
"rapidly hydrolyzed" comment on clopidogrel's very fast parent clearance —
while `EFFECT_RIVA_OUT` reaches 0.44 (44% FXa inhibition), a materially
larger, non-degenerate effect.

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`, all seven
`peripheral-arterial-disease` rows (Aspirin, Atorvastatin (ATST), Cilostazol,
Clopidogrel, Clopidogrel (AM), Rivaroxaban, Ticagrelor) — Clopidogrel (AM)'s
`Target/Pathway` corrected from `?` to reflect the parent-metabolite finding
above (P2Y12 receptor, as clopidogrel's own active thiol metabolite); the
other six rows' `Target/Pathway` filled in from the actual code (P2Y12
receptor for Clopidogrel parent/Ticagrelor — parent has no direct effect,
only PK; COX-1 for Aspirin; Factor Xa for Rivaroxaban; PDE3 for Cilostazol;
HMG-CoA reductase for Atorvastatin).

## Anything else flagged

- The post-DSL R script required three small, disclosed edits to keep the
  file runnable standalone after the rename: `cmt_map`'s names (cosmetic
  only — confirmed by grep that `cmt_map` is never referenced anywhere else
  in the script), and the sensitivity-analysis `params_to_vary` vector
  (`"EC50_P2Y12"→"EC50_CLOPAM"`, `"EC50_FXa"→"EC50_RIVA"`,
  `"Emax_HMG"→"EMAX_ATST"`; `"k_plaque"`/`"ABI0"` are disease-side params,
  unchanged). Everything else — all seven `scenarios` event definitions
  (numeric `cmt=` targets are unchanged, since compartment declaration order
  is preserved), `run_scenario()`'s `outvars=` list, all eight figures, and
  the summary table — is unchanged, since every capture name the post-DSL
  script references was kept under its original name (see `$CAPTURE`
  naming section above).
- `mread_cache()`'s model identifier was changed from `"pad_qsp"` to
  `"pad_qsp_refactored"` purely to avoid an `mrgsolve` cache collision if
  both files are ever compiled in the same R session; non-behavioral.
