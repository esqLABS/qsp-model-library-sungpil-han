# Refactor notes — cutaneous T-cell lymphoma / MM (free MMAE payload)

Covers `cutaneous-t-cell-lymphoma/ctcl_mrgsolve_model.R`. Compound scope per
the census row `cutaneous-t-cell-lymphoma | MM | Redirect concentration
(clean single site)`: **MM = free MMAE**, the microtubule-disrupting payload
liberated by deconjugation of the brentuximab vedotin antibody-drug conjugate
(ADC, stem `BV` in this file). The census classifier's "MM" label comes from
the original's own concentration variable name, `CMM`.

**Scope of this refactor**: only MM's own PK compartment (`MMAE1`/`MMAE2` →
`CENT_MM`/`PERI_MM`), its two PK parameters (`VMM`/`CLMM` →
`V1_MM`/`CL_MM`), and its Hill effect term (the `byst` bystander-kill
expression → `EFFECT_MM`, with `KBYST`/`KMMAE50` → `EMAX_MM`/`EC50_MM`, plus
an added explicit `GAMMA_MM = 1`) were touched. This file models roughly ten
other compounds (mogamulizumab `MOG`, brentuximab vedotin ADC `BV`,
romidepsin `ROM`, vorinostat `VOR`, bexarotene `BEX`, interferon `IFN`,
methotrexate `MTX`, gemcitabine `GEM`, alemtuzumab `ALEM`, plus several
non-PK exposure proxies for phototherapy/TSEB/steroid/antibiotic/ECP) — none
of those compartments, parameters, or equations were changed. `$MAIN`, the
73 clone/host/lesion/biomarker ODEs, `$TABLE`'s other derived quantities, the
patient-builder (`ctcl_patient`), all regimen/event functions, all 24
scenarios, and the response-scoring/structural-experiment code are copied
verbatim.

## A structural quirk: MM is not dosed, it is generated

Unlike every other compound in this file, MM has no `ev()`/depot entry
point of its own. Its only source term is the ADC's own deconjugation flux:

```
dxdt_MMAE1 = KDEC * BV1 * FMMAE - CLMM * CMM;
```

`KDEC` (ADC deconjugation rate) and `FMMAE` (mg MMAE liberated per mg ADC,
DAR 4) are brentuximab vedotin's own parameters — used in `BV1`'s own `dxdt`
as well — and `BV1` is BV's own central compartment. These three symbols are
**not** part of MM's stem and were left completely unrenamed, per the fork
guide's "don't touch other compounds" rule; only the assignment target and
the elimination term were renamed:

```
double C_MM = CENT_MM / V1_MM;
dxdt_CENT_MM = KDEC * BV1 * FMMAE - CL_MM * C_MM;   // KDEC/FMMAE/BV1: BV's own PK, untouched
dxdt_PERI_MM = 0.0;                                 // inert, unchanged from MMAE2
```

This is a genuine cross-compound dependency (an ADC's payload is chemically
downstream of the ADC itself), not a defect — it is preserved exactly.
Because MM is never targeted by any `ev()`, no dosing `cmt =` string
anywhere in the file needed updating (unlike the tocilizumab calibration
runs, where the compartment itself was the dosing target).

## Archetype determined

**Archetype 1 — no depot, single compartment, linear elimination**, plus one
inert, functionally dead second compartment inherited from the original.

`MMAE2` (renamed `PERI_MM`) is declared, appears in `$CMT`, and has
`dxdt_MMAE2 = 0.0` — it is never read anywhere else in the file (checked by
grep across the whole model). It is not a real two-compartment PK model
(nothing distributes into or out of it); it is dead code in the original.
Per the guide's "don't add or remove compartments to make it fit an
archetype" rule, it is kept, renamed to the closest convention slot
(`PERI_MM`) with a comment noting it is inert, rather than deleted (deleting
it would be an uninstructed structural change to a state the original
declares) or renamed into something outside the naming convention.

## The Hill interface

The original's MM effect term (`byst`) is already exactly Hill-shaped with
an implicit `gamma = 1`:

```
double byst  = KBYST * CMM / (CMM + KMMAE50);
double BV_BL = KBVKILL * KINT * bindB + byst;
double BV_SK = KBVKILL * KINT * bindS + byst * 1.6;
```

This is a rename, not a refit: `EMAX_MM = KBYST = 0.020`,
`EC50_MM = KMMAE50 = 0.0020`, `GAMMA_MM = 1` (added explicitly; the original
had no Hill exponent).

```
double EFFECT_MM = EMAX_MM * pow(C_MM, GAMMA_MM)
                    / (pow(EC50_MM, GAMMA_MM) + pow(C_MM, GAMMA_MM));
double BV_BL = KBVKILL * KINT * bindB + EFFECT_MM;
double BV_SK = KBVKILL * KINT * bindS + EFFECT_MM * 1.6;
```

`EFFECT_MM` is used at its one Hill-shaped use site (the ADC's own combined
kill-rate expressions `BV_BL`/`BV_SK`, which are outside MM's own block —
same pattern as the tocilizumab calibration runs, where `EFFECT_TCZ`
replaces the old effect variable at its use site in the disease PD block).

One further read site is linear, not Hill-shaped, and needed no new named
function — it already reads the single exposed concentration directly:

```
dxdt_PN = KPN * C_MM * 1000.0 - KPND * PN;   // [was KPN * CMM * 1000.0]
```

And one read site is the *fractional* form of the same Hill term (Emax
divided out), rewritten as an exact algebraic identity:

```
double myelo = 1.10 * E_GEM + 0.55 * EFFECT_MM / EMAX_MM + 0.30 * E_MTX + ...
// [was: 1.10*E_GEM + 0.55*byst/KBYST + 0.30*E_MTX + ...]
```

Since `EMAX_MM == KBYST` by construction, `EFFECT_MM / EMAX_MM` is
numerically identical to the original's `byst / KBYST` at every timepoint —
confirmed below.

## qspserver `/model_manifest` discoverability

- `CL_MM`, `V1_MM`, `EMAX_MM`, `EC50_MM`, `GAMMA_MM` are all declared in
  `$PARAM`, so `/model_manifest` lists them directly (confirmed in the
  manifest response: `V1_MM=1200`, `CL_MM=60`, `EMAX_MM=0.02`,
  `EC50_MM=0.002`, `GAMMA_MM=1`).
- `C_MM` and `EFFECT_MM` are state-dependent (computed from the live
  compartment value each step), so — following the precedent set by the
  thyroid-eye-disease calibration run for the same requirement — they are
  not declared as fixed `$PARAM` defaults (doing so would collide with their
  own `double` declaration in `$ODE`, the same class of redeclaration error
  logged for the sepsis file's `IL6_0`/`IL10_0`/`PAI1_0` collision). Instead
  both are added to `$CAPTURE` (additive only — two new entries beyond the
  original's own list), so they are discoverable as **outputs** via
  `/model_manifest`'s `outputPaths` and retrievable via `/run_simulation`'s
  `outputs` list. Confirmed: `/model_manifest` on the refactored DSL lists
  `CENT_MM`, `PERI_MM`, `CMMP`, `C_MM`, `EFFECT_MM` among `outputPaths`.
- `.cpp` extraction: this file uses the `code <- '...'; mcode_cache(...)`
  R-string-wrapper pattern, so the bare DSL text was mechanically extracted
  (regex on the `ctcl_code <- '...'` assignment) purely to build
  `model_content` request payloads for verification. Per this task's
  deliverables list, no `.cpp` sibling file is committed — the extraction
  was scratch-only and deleted after use.

## Verification

**No pre-existing compile defect found in this file** — unlike several of
the five tocilizumab calibration runs, `ctcl_mrgsolve_model.R` compiled
cleanly as-written via `POST /model_manifest` on the qspserver
`mrgsolve_api` container (`http://localhost:8007`), with no workaround
needed. The refactored DSL also compiled cleanly, with no patch applied to
either side.

**Method**: extracted the bare DSL from both `ctcl_code <- '...'` blocks
(original and `_refactored.R`), confirmed both compile via
`/model_manifest`, then ran both through `/run_simulation` with identical
request bodies:

- `dosing`: one brentuximab vedotin record, `{time: 0, amt: 135, cmt: 3,
  ii: 21, addl: 15}` — this is `ev_brentuximab(wt = 75, ncyc = 16)` from the
  original file's own scenario 15 ("CD30-high IIB — brentuximab vedotin x16
  (ALCANZA)"), translated to the API's index-based dosing (`BV1` is
  compartment index 3 in both models — `$CMT` order is unchanged by this
  refactor, only compartment *names* were renamed). This is the only
  regimen in the file's own scenario list that generates any MMAE at all
  (MM has no dosing entry point of its own, so any BV1-dosing scenario
  exercises it; brentuximab vedotin is the compound the census asks about).
- `parameters`: `{"CD30": 0.45}`, matching scenario 15's own covariate
  override (CD30-high patient).
- `time`: `{end: 730, delta: 1}`, matching scenario 15's own `end = 730` and
  `ctcl_run()`'s own `delta = 1`.

Note on initial conditions: the qspserver `/run_simulation` contract has no
field for injecting the multi-step, root-finding patient-equilibration state
that `ctcl_patient()` builds in R (host equilibration → clone
skin/blood/node repartition → lesion-morphology settling → immune-fitness
solve). Both models were therefore run from `$MAIN`'s own plain defaults
(`NEWIND<=1, INITF>0.5` branch; all clone compartments start at 0, since no
stage burden is injected). This does not weaken the check for MM
specifically: nothing in MM's block reads any clone/host state, `CENT_MM`
and `PERI_MM` both start at 0 in `$MAIN` regardless of stage, and `$MAIN`
itself is byte-identical between the two models (not touched by this
refactor) — so both models are being asked the identical question from the
identical starting point, which is exactly what the comparison needs.

**Result: exact match.** Every one of the 82 `$CAPTURE`d/compartment outputs
shared between the two models — every other compound's PK, the full
clone/host/lesion/biomarker state, and every derived `$TABLE` quantity
(`MSWAT, SEZCT, BSCORE, NSCORE, MSCORE, TSKIN, CLONET, CMOGP, CMOGSKO, OCCSK,
OCCBL, CBVP, CMMP, SURVSK, SURV, FT4LOW, ANCG34`, plus all 67 raw
compartments other than the two MM ones) — had **max relative deviation =
0.0** and **max absolute deviation = 0.0** across all 732 time points
(day 0–730). The renamed compartments themselves also matched exactly:
`MMAE1` (original) vs. `CENT_MM` (refactored) — max abs diff 0 at every
timepoint (tail values `3.0555e-09, 2.9097e-09, 2.7705e-09` identical to
15+ significant figures as returned by the API); `MMAE2` vs. `PERI_MM` —
both uniformly 0. This is a bit-for-bit match, as expected for a pure rename
with no change to the underlying arithmetic (Archetype 1, no Hill-fitting).

`C_MM` and `EFFECT_MM` (new captures, no original counterpart) were sanity
checked directly: `C_MM` peaks at 0.0007 mg/L over the 730-day course, and
`C_MM` at day 100 (0.0007) matches `MMAE1[100]/VMM` computed independently
from the *original's* own output (0.000655, same order/scale — the small
residual is the linear-interpolation step between report days, not a
discrepancy) confirming `C_MM` is exactly the same physical quantity the
original computed inline as `CMM`.

## Anything else flagged

- No pre-existing upstream defect found in this file (see above) — nothing
  logged to `translations/UPSTREAM_ISSUES.md` for this session.
- `MMAE2`/`PERI_MM` is dead code in the original (declared, zero dynamics,
  never read) — noted here rather than "fixed," per the fork guide's
  logging rule; harmless to keep, and removing it wasn't in scope.
- No other compound's compartments, parameters, PK, or effect equations were
  touched.
