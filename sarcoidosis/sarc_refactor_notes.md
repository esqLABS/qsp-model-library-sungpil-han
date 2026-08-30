# Refactor notes — `sarcoidosis/sarc_mrgsolve_model.R`

Scope: **prednisone (PRED), prednisolone (PREDL), and methotrexate (MTX)
only** — their PK blocks and their downstream effect equations.
Azathioprine (AZA), hydroxychloroquine (HCQ), and infliximab (IFX) — all
modeled in the original as constant steady-state concentration
placeholders (`AZA_CSS`/`HCQ_CSS`/`IFX_CSS`), not real PK compartments —
are completely untouched: same params, same `Eaza`/`Ehcq`/`Eifx` formulas,
same combination structure. All 13 disease-PD compartments and their
kinetics are untouched except for a build-compatibility rename (see
below), which touches only baseline-parameter names and the
initial-condition idiom, not any dynamics.

## PREDL is prednisone's own active metabolite, not a second drug

The task brief flagged this as something to check rather than assume, and
the code confirms it cleanly. The original file's own header already says
so (`Drug PK: Prednisone/Prednisolone, Methotrexate ...` and, in the model
scope comment, `Prednisone (prodrug) --> Prednisolone (active)`), and the
`$ODE` structure backs it up:

- `PREDL_C`/`PREDL_P` have **no `GUT_`/depot compartment and no `ka_`/`F_`
  absorption parameters of their own** anywhere in the file.
- `PREDL_C`'s only inflow term is a conversion computed directly from
  `PRED_C`'s own plasma concentration:
  ```
  double pred_conc_ug_mL = PRED_C / Vc_pred_ode;
  dxdt_PREDL_C = k_predl_conv * CL_pred * pred_conc_ug_mL * Vc_predl_ode / 1000.0
                + k21_predl * PREDL_P - k_predl_el_ode * PREDL_C - k12_predl * PREDL_C;
  ```
  `k_predl_conv` is literally named "fraction prednisone converted to
  prednisolone" in its own `$PARAM` comment.
- The disease-effect parameters (`EC50_pred`/`Emax_pred`) are computed
  against `predl_c` (derived from `PREDL_C`), **never** against `PRED_C`
  — i.e. prednisone itself has no direct pharmacodynamic effect in this
  model; only its metabolite does.

**Handled as a genuine parent-metabolite chain, not two independently
dosed compounds:** `PRED` keeps its own PK block (`GUT_PRED`/`CENT_PRED`)
and its own `C_PRED`, which exists **only** to drive the conversion into
`PREDL` (prednisone has no `EFFECT_PRED`, matching the original having no
direct prednisone PD effect). `PREDL` keeps its own PK block
(`CENT_PREDL`/`PERI_PREDL`, fed by that conversion) and exposes the single
`C_PREDL`/`EFFECT_PREDL` pair the disease equations actually read. The
census row for PREDL (previously marked `Target/Pathway: ?`) is corrected
below to reflect that it is prednisone's metabolite sharing the
glucocorticoid-receptor mechanism, not an independent compound with an
unknown target.

## Archetype per compound

**PRED (prednisone) — Archetype 3 minus the peripheral compartment**
(depot + central, linear elimination). Confirmed from the actual `dxdt_`
lines:
```
dxdt_PRED_GUT = -ka_pred * PRED_GUT;
dxdt_PRED_C  = ka_pred * PRED_GUT - k_pred_el_ode * PRED_C;
```
No bioavailability factor in the original (`ka_pred*PRED_GUT` enters
`PRED_C` directly) — none invented here either.

**PREDL (prednisolone) — Archetype 2** (2-compartment, linear, CL/Q/V),
fed by the conversion inflow above. The original parameterizes the
central↔peripheral transfer as classic micro-constants (`k12_predl`,
`k21_predl`, applied directly to compartment masses, not divided by any
volume in the ODE) rather than CL/Q/V — per the guide, rewritten to
CL/Q/V. Because `V1_PREDL` (like `V1_PRED`, `V1_MTX`) is expressed
per-kg and multiplied by `WT` inside `$MAIN`/`$ODE`, a literal
`Q = k12*V1` computed once at a fixed body weight would only reproduce
the original at that one weight. Instead, `Q_PREDL`/`V2_PREDL` are
derived **per-kg**, matching `V1_PREDL`'s own convention, so the ratio is
weight-invariant:

```
Q_PREDL  = k12_predl * Vd_predl = 0.80 * 1.50 = 1.20   (L/h per kg)
V2_PREDL = Q_PREDL / k21_predl  = 1.20 / 0.40 = 3.00   (L/kg)
```

At runtime: `Q_PREDL_L = Q_PREDL*WT`, `V1_PREDL_L = V1_PREDL*WT`,
`V2_PREDL_L = V2_PREDL*WT`. Then for any `WT`:
`Q_PREDL_L/V1_PREDL_L = Q_PREDL/V1_PREDL = 1.20/1.50 = 0.80 = k12_predl`,
and `Q_PREDL_L/V2_PREDL_L = Q_PREDL/V2_PREDL = 1.20/3.00 = 0.40 = k21_predl`
— exactly the original's micro-constants, for every body weight, not just
the default `WT=70`. Combined with `CL_PREDL` (kept flat, exactly as the
original's `CL_predl` was — never itself WT-rescaled, only `V1_PREDL_L`
is), the refactored central-compartment elimination term
`(CL_PREDL + Q_PREDL_L)/V1_PREDL_L` algebraically equals the original's
`k_predl_el_ode + k12_predl`, and the peripheral term
`Q_PREDL_L/V1_PREDL_L` / `Q_PREDL_L/V2_PREDL_L` equals the original's
`k12_predl`/`k21_predl` term-for-term. This is a pure reparameterization
(derived from the original's own three values), not a new or fitted
value.

**MTX (methotrexate) — Archetype 3 minus the peripheral compartment**
(depot + central, linear) **plus a bespoke third compartment** (`POLY_MTX`,
name unchanged — already matched the convention) for the intracellular
active polyglutamate moiety that the original's own Hill effect term
actually reads:
```
dxdt_MTX_GUT = -ka_mtx * MTX_GUT;
dxdt_MTX_C   =  ka_mtx * MTX_GUT - k_mtx_el_ode * MTX_C - k_poly * MTX_C;
double mtx_nmol = MTX_C * 1e6 / 454.0;
dxdt_MTX_POLY = k_poly * mtx_nmol - k_poly_el * MTX_POLY;
```
Not a second independently-dosed compound — an accumulation/elimination
extension of MTX's own single stem. The identical `POLY_MTX` name and
structure is already used for methotrexate in
`rheumatoid-arthritis/ra_mrgsolve_model_refactored.R` (`GUT_MTX`/
`CENT_MTX`/`POLY_MTX`), so this refactor follows that established
corpus-wide precedent rather than inventing a new one.

## Renaming (PRED/PREDL/MTX only)

| Original | Refactored |
|---|---|
| `PRED_GUT` | `GUT_PRED` |
| `PRED_C` | `CENT_PRED` |
| `ka_pred` | `KA_PRED` |
| `Vd_pred` | `V1_PRED` |
| `CL_pred` | `CL_PRED` |
| `PREDL_C` | `CENT_PREDL` |
| `PREDL_P` | `PERI_PREDL` |
| `k_predl_conv` | `FCONV_PREDL` |
| `ka_predl` | `KA_PREDL` (unused in the original — see below) |
| `Vd_predl` | `V1_PREDL` |
| `CL_predl` | `CL_PREDL` |
| `k12_predl`, `k21_predl` | `Q_PREDL`, `V2_PREDL` (derived, see above) |
| `Emax_pred` | `EMAX_PREDL` |
| `EC50_pred` | `EC50_PREDL` |
| — | `GAMMA_PREDL = 2.000` (new; copied from the original's shared `HILL_n`) |
| `MTX_GUT` | `GUT_MTX` (name unchanged) |
| `MTX_C` | `CENT_MTX` |
| `ka_mtx` | `KA_MTX` |
| `Vd_mtx` | `V1_MTX` |
| `CL_mtx` | `CL_MTX` |
| `k_poly` | `KPOLY_MTX` |
| `k_poly_el` | `KDEPOLY_MTX` |
| `MTX_POLY` | `POLY_MTX` (name unchanged, already matched convention) |
| `Emax_mtx` | `EMAX_MTX` |
| `EC50_mtx_poly` | `EC50_MTX` |
| — | `GAMMA_MTX = 2.000` (new; copied from the original's shared `HILL_n`) |
| `Epred` | `EFFECT_PREDL` |
| `Emtx` | `EFFECT_MTX` |

`ka_predl`/`KA_PREDL` is declared in the original ("Prednisolone
absorption/intercompartment rate") but **never referenced by any `dxdt_`
line anywhere in the file** — preserved unused, renamed only, per the
"don't invent, don't drop the original's parameters" rule (same treatment
already given to `RTOT_TCZ` in `rheumatoid-arthritis/ra_refactor_notes.md`).

`CL_PRED`, `CL_PREDL`, `CL_MTX` carry one more preserved quirk: their own
`$PARAM` comments all read "L/h per 70 kg", but none of them are actually
rescaled by `WT/70` anywhere in the ODE — only each compound's own `V1_*`
is WT-scaled. This is a pre-existing property of the original's own math
(so e.g. elimination half-life is implicitly weight-dependent through
volume only, not through clearance), preserved exactly rather than
"fixed," since fixing it would be a numeric/behavioral change out of
scope for a rename.

All parameter *values* are copied verbatim from the original — nothing
invented except `GAMMA_PREDL`/`GAMMA_MTX` (both `2.000`, copied from the
original's own shared `HILL_n`, since the original had no per-compound
Hill exponent — `HILL_n` itself is untouched and still governs
`Eaza`/`Ehcq`/`Eifx`, out of this refactor's scope).

## Hill interface: rename, not a fit — for both compounds

Both `Epred` and `Emtx` in the original were already the exact
`Emax*C^n/(EC50^n+C^n)` shape (shared exponent `n = HILL_n = 2`):
```
double Epred = Emax_pred * pow(predl_c, n) / (pow(EC50_pred, n) + pow(predl_c, n) + 1e-9);
double Emtx  = Emax_mtx  * pow(mtxpg_c, n) / (pow(EC50_mtx_poly, n) + pow(mtxpg_c, n) + 1e-9);
```
`EFFECT_PREDL`/`EFFECT_MTX` use the identical formula with the renamed
per-compound `EMAX_`/`EC50_`/`GAMMA_` parameters — same values, same
`+1e-9` zero-division guard, same execution order. This is a rename, not
a refit; no curve fitting was needed or performed, and no R² is reported
because none applies.

Every direct downstream use of `Epred` was redirected to `EFFECT_PREDL`
(four sites: `dxdt_TREG`'s `(1+Epred*0.2)` term, `dxdt_IFNG`'s
`(1-Epred*0.6)` term, `dxdt_IL12`'s `(1-Epred*0.5)` term, `dxdt_CALIT`'s
`(1-Epred*0.4)` term) plus its use inside the combined `Elymph`/`Emac`/
`Etnf_inhib` expressions in `$ODE` Section 2. `Emtx` was used only inside
`Elymph`, redirected the same way. AZA's/HCQ's/IFX's own terms inside
these combined expressions (`Eaza`, `Ehcq`, `Eifx`) are byte-identical to
the original.

## `C_MTX` — the exposed concentration is the compartment itself

`POLY_MTX` (the intracellular active polyglutamate pool) **is** the
concentration the original's Hill term reads (`mtxpg_c = MTX_POLY`, no
further division) — the same "compartment IS concentration" pattern
already documented in `eosinophilic-esophagitis/eoe_refactor_notes.md` and
`dengue/denv_refactor_notes.md`. `C_MTX = POLY_MTX` needs no further
transformation. A separate, purely informational plasma-level
concentration (`Cp_MTX_plasma = CENT_MTX/V1_MTX_L`) is also captured but
is **not** the exposed effect-driving concentration — same "kept only as
an informational, non-exposed diagnostic" pattern already documented for
the tissue/plasma split in `abdominal-aortic-aneurysm/aaa_refactor_notes.md`.

Also preserved exactly (not "fixed"): the original computes the
mg/L→nmol/L conversion feeding `MTX_POLY` directly from `MTX_C` (a mass in
mg, per the `$CMT` comment: "MTX central plasma (mg/L equivalent)"),
i.e. `mtx_nmol = MTX_C * 1e6/454.0`, not from an actual mg/L
concentration. Renamed `MTXPG_input_nmol = CENT_MTX * 1e6/454.0` —
identical formula, identical (already slightly ambiguous) units handling.

## Build-compat fix (mrgsolve 2.0.1) — logged as `UPSTREAM_ISSUES.md` #60

Confirmed via the qspserver `mrgsolve_api` container (`http://localhost:8007`)
that the **untouched original does not compile at all**, entirely inside
the disease-PD initial-condition block — unrelated to PRED/PREDL/MTX:

1. All 13 disease-compartment initial conditions are set with the
   deprecated `_init_<CMT> = value;` idiom, not accepted by mrgsolve
   2.0.1 (`'_init_MAC_ACT' was not declared in this scope`, etc.). Same
   defect class as issues #29, #41, #42, #59.
2. `$PARAM` baseline names `TH1_0`, `TREG_0`, `IL12_0` collide with
   mrgsolve's own auto-reserved per-compartment `<CMT>_0` initial-value
   symbol for compartments `TH1`, `TREG`, `IL12`
   (`conflicting declaration 'double& TREG_0'`, etc.) — a distinct defect
   from (1), reproducing independently.
3. A pre-existing typo, independent of (1)/(2): the original declares
   `TREG_0` but reads it back as `TREG0` (no underscore) in three places
   (`_init_TREG = TREG0;`, and twice in the macrophage-suppression logic,
   `treg_ratio`/`treg_suppress`). `TREG0` is never declared anywhere, so
   as literally written the model could never have compiled under any
   mrgsolve version requiring every referenced symbol to exist.

**Fix applied directly to the delivered `sarc_mrgsolve_model_refactored.R`**
(not just a scratch copy), per the guide's settled policy, all
syntax-only and non-numeric: (a) `TH1_0`/`TREG_0`/`IL12_0` renamed to
`TH10`/`TREG0`/`IL120` — matching the no-underscore convention every
other baseline param in the same block already uses (`MAC_ACT0`, `TNF0`,
`IFNG0`, `GRAN0`, `FIBR0`, `ACE0`, `CALIT0`, `CA0`, `SIL2R0`, `FVC_P0`),
which also fixes the Defect-3 typo in the same stroke, since `TREG0` — the
name the code actually reads — becomes the declared name; (b) all 13
`_init_<CMT> = value;` lines converted to the modern
`<CMT>_0 = value;` idiom (no compartment, index, or starting value
changed). Full disclosure and the exact error traces are in
`translations/UPSTREAM_ISSUES.md` #60. The checked-in original
(`sarc_mrgsolve_model.R`) is untouched and still carries all three
defects exactly as written; the identical pair of fixes was applied to an
in-memory-only scratch copy of the original (never committed) purely to
build a working comparison target for verification.

## `$CAPTURE` naming for the state-dependent quantities

`C_PRED`, `C_PREDL`, `Cp_MTX_plasma`, `C_MTX`, `EFFECT_PREDL`,
`EFFECT_MTX` are computed every `$ODE` evaluation from compartment state,
so — per the precedent in `eosinophilic-esophagitis/eoe_refactor_notes.md`
and `thyroid-eye-disease/ted_refactor_notes.md` — they cannot themselves
be `$PARAM` (a `$PARAM` of the same name would be a redeclaration against
the `$ODE`-local `double`). They are captured in `$TABLE` instead, with an
`_OUT` suffix (`C_PRED_OUT`, `C_PREDL_OUT`, `Cp_MTX_plasma_OUT`,
`C_MTX_OUT`, `EFFECT_PREDL_OUT`, `EFFECT_MTX_OUT`) — a bare
`capture C_PRED = C_PRED;` self-reference is itself a
`redefinition of capture {anonymous}::C_PRED` under mrgsolve 2.0.1 (this
file's own first attempt hit exactly that error), the same collision
class as `UPSTREAM_ISSUES.md` #35/#43/#50/#59, here newly introduced by
this refactor's own capture choice rather than pre-existing in the
original (the original's captures used distinct names —
`PREDL_conc_ngmL`, `MTX_plasma_ugL`, `MTXPG_nmolL`, `E_pred_fx`,
`E_mtx_fx` — so this is *not* a new upstream defect to log, just an
authoring choice fixed here). Confirmed discoverable via
`/model_manifest`'s `outputPaths` (see Verification).

## qspserver compatibility checklist

- `model_content` is pure mrgsolve DSL text (no R wrapper) — confirmed by
  extracting the quoted `sarc_code <- '...'` block byte-for-byte and
  building it standalone via `POST /model_manifest`.
- Every `EMAX_`/`EC50_`/`GAMMA_` parameter for both compounds lives in
  `$PARAM` (confirmed present in the manifest's `parameters` list, 70
  total).
- `C_PRED_OUT`/`C_PREDL_OUT`/`Cp_MTX_plasma_OUT`/`C_MTX_OUT`/
  `EFFECT_PREDL_OUT`/`EFFECT_MTX_OUT` are discoverable `outputPaths` (see
  above).
- `$SET end/delta` is untouched (not present in the original either — the
  R-side `run_sim()` wrapper supplies `end`/`delta` explicitly); the
  verification below drives `time.end`/`time.delta` through the request,
  not through any model-internal default.
- No R-only syntax inside the DSL block — the extracted `.cpp` compiled
  standalone with no surrounding R script.
- Solver step budget: 104 weeks at `delta=6h` is 2913 points — well under
  the API's default 20000-step budget; no shortening was needed for
  either scenario.

## Verification

Per the guide's mandatory protocol: extracted the quoted DSL block from
both the build-compat-patched original (scratch copy only, never
committed) and the delivered `sarc_mrgsolve_model_refactored.R`, confirmed
both build via `POST /model_manifest` on the qspserver `mrgsolve_api`
container (`http://localhost:8007`, confirmed healthy throughout, requests
spaced ~2s apart), then ran the original file's own dosing scenarios that
exercise the three refactored compounds through `POST /run_simulation`:

| Scenario | Dosing (from the original's own `dose_pred()`/`dose_mtx()`/`taper_pred()` helpers) | Result |
|---|---|---|
| Scenario 2 — Prednisone 40mg taper | `GUT_PRED`/`PRED_GUT` (cmt 1): 40mg→20mg→10mg→5mg, q24h, 4 weeks each | exact match, max abs diff 3.73e-10 |
| Scenario 3 — Prednisone 20mg + MTX 10mg/wk | `GUT_PRED` (cmt 1): 20mg→10mg→5mg taper over 36 weeks; `GUT_MTX` (cmt 5): 10mg q168h × 104 weeks | exact match, max abs diff 1.00e-12 |

Both scenarios run the full 104-week (17472h) window at `delta=6h`
(2913 points) — the same window the original's own `run_sim()` default
uses — no shortening needed; neither approached mrgsolve's default
20000-step budget.

Every shared output was compared point-by-point across the full time
grid: all 20 raw compartments (mapped through the renames:
`PRED_GUT↔GUT_PRED`, `PRED_C↔CENT_PRED`, `PREDL_C↔CENT_PREDL`,
`PREDL_P↔PERI_PREDL`, `MTX_GUT↔GUT_MTX`, `MTX_C↔CENT_MTX`,
`MTX_POLY↔POLY_MTX` for the refactored compounds; the 13 disease
compartments unchanged), all 20 `$CAPTURE` outputs (mapped
`PREDL_conc_ngmL↔C_PREDL_OUT`, `MTX_plasma_ugL↔Cp_MTX_plasma_OUT`,
`MTXPG_nmolL↔C_MTX_OUT`, `E_pred_fx↔EFFECT_PREDL_OUT`,
`E_mtx_fx↔EFFECT_MTX_OUT`; the 15 disease-biomarker/flag captures
unchanged), and `C_PRED_OUT` (no original counterpart — prednisone plasma
concentration was never separately captured in the original; checked only
for internal self-consistency, `CENT_PRED/V1_PRED_L`).

**Result: exact match to floating-point noise for every output, both
scenarios** (max abs diff 3.73e-10 for Scenario 2, 1.00e-12 for
Scenario 3 — both scenarios' worst-case deviations occur on the PK
compartments themselves, not the disease outputs, and are consistent with
ordinary floating-point accumulation over ~2900 solver steps, not a
structural mismatch). This is the expected result for pure structural
reorganization with no Hill-refitting (both compounds' effect terms were
already exact Hill shapes, see above) — per the guide's own tolerance
rule, "anything beyond floating-point-scale deviation means a bug, not a
tolerance to loosen"; none was observed.

Sanity-checked the dynamics are non-trivial (not degenerately zero on
both sides): Scenario 2's `C_PREDL_OUT` reaches 7.35 ng/mL by t=24h then
continues rising with accumulation; Scenario 3's `C_MTX_OUT` (MTX
polyglutamates) reaches ~2790 nmol/L by week 104 — two to three orders of
magnitude above the original header's own stated "5-15 nmol/L" clinical
reference range. This is a pre-existing property of the original's own
arithmetic (confirmed identical, to floating-point noise, in both the
original and the refactored model — not something this refactor
introduced or could have introduced via a pure rename) and is flagged here
for whoever next works on this file's calibration, not logged as a build
defect since the model compiles and runs correctly with this value.

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`, the
three `sarcoidosis` rows (MTX, PRED, PREDL) — PREDL's `Target/Pathway`
corrected from `?` to `Glucocorticoid receptor (as prednisone's own active
metabolite)` per the parent-metabolite finding above.

## Anything else flagged

- AZA/HCQ/IFX (`Eaza`, `Ehcq`, `Eifx`, and every param/CSS-state they
  read) are byte-identical to the original in the refactored file —
  confirmed by the exact-match verification above, since Scenario 2 and
  Scenario 3 both hold `AZA_CSS=HCQ_CSS=IFX_CSS=0` and the disease
  outputs still matched exactly.
- The post-DSL R script (helper functions, scenario definitions, plots,
  sensitivity analysis, validation table) required five small, disclosed
  edits to keep the file runnable standalone after the rename: `dose_pred()`/
  `dose_mtx()`'s `cmt=` targets (`"PRED_GUT"→"GUT_PRED"`,
  `"MTX_GUT"→"GUT_MTX"`), the Scenario-1 null event's `cmt=`, Plot 5/6's
  y-aesthetic columns (`PREDL_conc_ngmL→C_PREDL_OUT`,
  `MTXPG_nmolL→C_MTX_OUT`), and the sensitivity-analysis param name
  (`EC50_pred→EC50_PREDL`, including its own R-side bookkeeping columns).
  Everything else in the post-DSL script — `run_sim()`, the six-scenario
  color palette, all plot styling, the summary tables, the clinical
  validation targets — is unchanged.
- `mcode()`'s model identifier was changed from `"sarcoidosis_qsp"` to
  `"sarcoidosis_qsp_refactored"` purely to avoid an `mrgsolve` `soloc`
  cache collision if both files are ever compiled in the same R session;
  non-behavioral.
