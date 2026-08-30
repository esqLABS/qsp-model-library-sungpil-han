# Refactor notes — `chronic-lymphocytic-leukemia/cll_mrgsolve_model.R`

Three-compound file: **Ibrutinib (IB)**, **Venetoclax (VEN)**, **Obinutuzumab
(OBI)** — all three rows in `driver-patches/data/compound_perturbation_census.md`
were classified "Redirect concentration (clean single site)". The disease
model (ALC/BM_CLL/LN_CLL tumor burden, MCL1 adaptive resistance, NK
activation) is shared machinery all three compounds feed into and is left
completely untouched apart from reading the new `EFFECT_<STEM>` names in
place of the original's local `E_BTKi`/`E_BCL2i`/`E_CD20` variables (same
arithmetic, same order of operations).

The original file did not compile at all under mrgsolve 2.0.1 (four
compounding build defects, none related to any compound's own PK/PD math —
see "Pre-existing build defects" below). The delivered `_refactored.R`
carries a syntax-only, disclosed fix for all four, per the guide's settled
policy for "when the original doesn't compile at all."

## Archetype per compound

### Ibrutinib (IB) — depot + central, single compartment (a no-peripheral variant of archetype 3), plus a bespoke covalent-occupancy submodel

PK: `GUT_IB` → `CENT_IB` (was `DEPOT_IB`/`CENT_IB`), one-compartment linear
elimination (`CL_IB`/`V1_IB`, was `Vd_IB`) — archetype 3 with the
peripheral compartment dropped, the same way the guide allows dropping
`GUT`/`KA`/`F` for a no-depot 2-compartment variant, applied symmetrically
here to a no-peripheral depot+central variant.

Downstream of PK, the original models BTK covalent occupancy as its own
two-state dynamical system (`BTK_FREE`/`BTK_OCC`, pseudo-first-order
`kinact`, shared `kdeg_BTK` turnover for both pools) — a receptor-binding
submodel, not a plain algebraic occupancy ratio. Renamed to
`REC_FREE_IB`/`COMPLEX_IB` (the guide's TMDD naming slots), **but this is
disclosed as a deliberate reuse of that naming, not literal TMDD**: unlike
the guide's archetype-4 example, ibrutinib's covalent binding never
feeds back into `CENT_IB`'s own mass balance (no term subtracts bound drug
from the central compartment) — it is downstream PD, not target-mediated
*disposition*. `RTOT_IB = 100` was added as a new named parameter for the
literal `100.0` the original hardcoded twice (the `$INIT` value and the
`dxdt_BTK_FREE` synthesis term) — a rename of an existing implicit
constant, not an invented value.

**Hill interface — exact rename, no fit.** The original's `E_BTKi =
Emax_BTK * BTK_OCC_pct / (EC50_BTK + BTK_OCC_pct)` already reduces to a
plain one-line Hill ratio (gamma=1 implicit) once `BTK_OCC_pct` (itself
just `BTK_OCC` renamed, no division) is treated as the ODE-solved input.
`EFFECT_IB = EMAX_IB * pow(COMPLEX_IB, GAMMA_IB) / (pow(EC50_IB,
GAMMA_IB) + pow(COMPLEX_IB, GAMMA_IB))`, `GAMMA_IB = 1` — same pattern as
`acute-intermittent-porphyria`'s GIV/HEM ("original's effect was already
an exact Hill ratio... a rename, not a fit").

### Venetoclax (VEN) — depot + central + peripheral, linear (archetype 3, textbook match), plus a bespoke quasi-equilibrium occupancy submodel

PK is already written in CL/Q/V form (`KA_VEN`, `V1_VEN`, `V2_VEN`,
`CL_VEN`, `Q_VEN`, `F_VEN`) — a direct archetype-3 match, renamed
`DEPOT_VEN`→`GUT_VEN` only (`CENT_VEN`/`PERI_VEN` already matched
convention).

Downstream, BCL-2 occupancy relaxes toward a quasi-steady-state target
(`BCL2_OCC_ss`, from an apparent `Ki_BCL2`) at rate `kout_BCL2` — renamed
`REC_FREE_VEN`/`COMPLEX_VEN`, `KD_VEN` (was `Ki_BCL2`), `RTOT_VEN` (was
`BCL2_tot`), `KOUT_VEN` (was `kout_BCL2`). Same disclosure as ibrutinib:
this never feeds back into venetoclax's own `CENT_VEN`/`PERI_VEN` mass
balance, so it is PD-only, not literal TMDD, despite reusing the
`REC_FREE_`/`COMPLEX_` naming slots.

**Hill interface — exact rename, no fit.** `E_BCL2i = Emax_BCL2 *
BCL2_OCC_pct / (EC50_BCL2 + BCL2_OCC_pct) / MCL1_ADAPT` combines two
things: a pure Hill ratio on occupancy, and a division by `MCL1_ADAPT`, an
adaptive-resistance disease state that is itself driven by
`BCL2_OCC_pct` (`stim_MCL1 = BCL2_OCC_pct/100`) but represents an
emergent disease response, not part of venetoclax's own PK/PD — the same
category of judgment call as the AMD refactor's decision to leave
`VEGF_FREE` unrenamed as shared disease state. `EFFECT_VEN` is defined as
**only the pure drug Hill term** (`EMAX_VEN * pow(OCC_FRAC_VEN,
GAMMA_VEN) / (pow(EC50_VEN, GAMMA_VEN) + pow(OCC_FRAC_VEN, GAMMA_VEN))`,
`GAMMA_VEN = 1`); the `/MCL1_ADAPT` division is kept exactly where the
original applied it, at the point of consumption (`E_BCL2i = EFFECT_VEN /
MCL1_ADAPT`, then into `kill_ALC`/`kill_BM`/`kill_LN`) — same arithmetic,
same order of operations as the original, just split so `EFFECT_VEN`
itself is a clean, named, single-compound function as the guide requires.

### Obinutuzumab (OBI) — no depot (IV), central + peripheral, genuine TMDD via quasi-steady-state relaxation (archetype-4 variant)

PK: `CENT_OBI`/`PERI_OBI` (already matched convention, no rename needed),
`CL_OBI`/`Q_OBI`/`V1_OBI`/`V2_OBI` (already matched convention). No
`GUT_OBI`/`KA_OBI` — dosed by IV bolus directly into `CENT_OBI`, matching
the guide's "drop the depot for a no-depot 2-compartment variant" allowance
extended to the TMDD archetype.

Unlike ibrutinib and venetoclax, obinutuzumab's CD20 binding **does feed
back into its own central-compartment mass balance**
(`dxdt_CENT_OBI` subtracts `KINT_OBI * (OCC_SS_OBI - COMPLEX_OBI) *
V1_OBI * 0.01`) — this is genuine target-mediated drug disposition, just
modeled via a quasi-steady-state relaxation (`OCC_SS_OBI = RTOT_OBI *
C_OBI / (KD_OBI + C_OBI)`, relaxed at rate `KINT_OBI`) rather than the
guide's explicit mass-action `KON`/`KOFF` template. Renamed
`REC_FREE_OBI`/`COMPLEX_OBI` (was `CD20_FREE`/`CD20_OCC`), `RTOT_OBI`
(was `CD20_0`), `KD_OBI` (was `Kd_CD20`, apparent affinity — there is no
explicit `KON_OBI`/`KOFF_OBI` pair in this QSS formulation, only the one
apparent-Kd + one relaxation-rate parameterization actually written),
`KINT_OBI`/`KSYN_OBI`/`KDEG_OBI` (were `kint_CD20`/`ksyn_CD20`/
`kdeg_CD20`).

**Dead code preserved, not removed:** the original declares `double
k_on_CD20 = 0.01; // apparent kon` inside `$ODE` and never references it
again anywhere downstream — a leftover from an earlier, presumably
mass-action version of this submodel. Renamed to `KON_OBI` and kept as a
dead local, exactly where the original had it, per the AMD/AAA precedent
of preserving pre-existing dead code rather than silently deleting it.

**Hill interface — exact rename, no fit.** `E_CD20 = Emax_CD20 *
CD20_OCC_pct * NK_ACT / (EC50_CD20 + CD20_OCC_pct)`, where `CD20_OCC_pct
= CD20_OCC/(CD20_FREE+CD20_OCC)*100` is already a fractional-occupancy
ratio computed from the ODE-solved receptor states — exactly the shape the
guide's own TMDD Hill interface describes (`X = COMPLEX/RTOT`-style
fractional occupancy), so despite the receptor states being ODE-solved,
the effect formula itself is already one Hill line and needed a rename,
not an `nls()` fit. Same treatment as venetoclax: `EFFECT_OBI` is the pure
drug Hill term
(`EMAX_OBI * pow(OCC_FRAC_OBI, GAMMA_OBI) / (pow(EC50_OBI, GAMMA_OBI) +
pow(OCC_FRAC_OBI, GAMMA_OBI))`, `GAMMA_OBI = 1`); the `* NK_ACT`
multiplier (an emergent disease-level NK-activation state, itself driven
by `CD20_OCC_pct`) is kept applied at the point of consumption
(`E_CD20 = EFFECT_OBI * NK_ACT`), same arithmetic as the original.

## Renaming summary (all parameter values copied verbatim from the original)

| Role | IB | VEN | OBI |
|---|---|---|---|
| Depot | `DEPOT_IB`→`GUT_IB` | `DEPOT_VEN`→`GUT_VEN` | none (IV) |
| Central | `CENT_IB` (unchanged) | `CENT_VEN` (unchanged) | `CENT_OBI` (unchanged) |
| Peripheral | none | `PERI_VEN` (unchanged) | `PERI_OBI` (unchanged) |
| Absorption/bioavail. | `Ka_IB`→`KA_IB`, `F_IB` (unchanged) | `Ka_VEN`→`KA_VEN`, `F_VEN` (unchanged) | n/a |
| Clearance/volume | `CL_IB` (unchanged), `Vd_IB`→`V1_IB` | `CL_VEN`,`V1_VEN`,`V2_VEN`,`Q_VEN` (unchanged) | `CL_OBI`,`V1_OBI`,`V2_OBI`,`Q_OBI` (unchanged) |
| Free receptor | `BTK_FREE`→`REC_FREE_IB` | `BCL2_FREE`→`REC_FREE_VEN` | `CD20_FREE`→`REC_FREE_OBI` |
| Complex | `BTK_OCC`→`COMPLEX_IB` | `BCL2_OCC`→`COMPLEX_VEN` | `CD20_OCC`→`COMPLEX_OBI` |
| Affinity/Kd | `Ki_BTK`→`KD_IB` | `Ki_BCL2`→`KD_VEN` | `Kd_CD20`→`KD_OBI` |
| Receptor total | new `RTOT_IB=100` (was literal `100.0`) | `BCL2_tot`→`RTOT_VEN` | `CD20_0`→`RTOT_OBI` |
| Turnover/relaxation | `kinact_BTK`→`KINACT_IB`, `kdeg_BTK`→`KDEG_IB` | `kout_BCL2`→`KOUT_VEN` | `kint_CD20`→`KINT_OBI`, `ksyn_CD20`→`KSYN_OBI`, `kdeg_CD20`→`KDEG_OBI` |
| Hill Emax/EC50/gamma | `Emax_BTK`→`EMAX_IB`, `EC50_BTK`→`EC50_IB`, new `GAMMA_IB=1` | `Emax_BCL2`→`EMAX_VEN`, `EC50_BCL2`→`EC50_VEN`, new `GAMMA_VEN=1` | `Emax_CD20`→`EMAX_OBI`, `EC50_CD20`→`EC50_OBI`, new `GAMMA_OBI=1` |
| Exposed concentration | new `C_IB` (nM) | new `C_VEN` (nM) | new `C_OBI` (mg/L) |
| Compound effect | new `EFFECT_IB` | new `EFFECT_VEN` | new `EFFECT_OBI` |
| Dead code | — | — | `k_on_CD20`→`KON_OBI` (unused, preserved) |

**Left unrenamed, deliberately (disease-level, not compound-specific):**
`ALC`, `BM_CLL`, `LN_CLL`, `MCL1_ADAPT`, `NK_ACT` (compartments);
`kprol_CLL`, `Kmax_ALC`, `BM_0`, `LN_0`, `kegress`, `egress_thr`,
`kin_MCL1`, `kout_MCL1`, `MCL1_max`, `kin_NK`, `kout_NK`, `NK_max`
(parameters); `use_IB`/`use_VEN`/`use_OBI` (pre-existing dead dose-flag
parameters — never read anywhere in `$ODE`/`$MAIN`/`$TABLE` in the
original; the R driver differentiates scenarios purely through which
dosing events are supplied, not through these flags).

**One exception required by a build-defect fix:** `ALC_0` (a baseline
parameter used in the `ALC_pch`/`PD_flag` calculations) was renamed to
`ALC_BASELINE` — not a naming-convention choice but forced by defect 3
below (collision with mrgsolve's own auto-generated `<CMT>_0` symbol for
compartment `ALC`).

## `C_IB`/`C_VEN`/`C_OBI`/`EFFECT_IB`/`EFFECT_VEN`/`EFFECT_OBI` are `$ODE`-local + `$CAPTURE`, not `$PARAM`

Same confirmed mrgsolve 2.0.1 constraint documented in
`age-related-macular-degeneration/amd_refactor_notes.md`,
`acute-intermittent-porphyria/aip_refactor_notes.md`, and the original
tocilizumab calibration run: a name cannot simultaneously be a settable
`$PARAM` default and a value recomputed fresh every `$ODE` step in this
build (`assignment of read-only reference`). All six values are declared
as plain `double` locals inside `$ODE` and exposed through `$CAPTURE`
only — confirmed visible via `POST /run_simulation`'s `outputs` selection,
confirmed **not** listed by `POST /model_manifest`'s parameter listing (51
parameters total, none of the six). A future driver-patch redirect of, say,
`C_IB` will need to replace the `double C_IB = ...;` line itself (or use a
`$MAIN`/covariate-injection mechanism), not a `param()`-style override.

## Pre-existing build defects found (all four fixed, syntax-only, in the delivered `_refactored.R`; logged as `translations/UPSTREAM_ISSUES.md` #47)

None of these four are related to any of the three compounds' own PK/PD —
they block the file from compiling at all, original or refactored, until
worked around. Per the guide's settled policy ("When the original doesn't
compile at all"), all four fixes are applied directly in the delivered
`cll_mrgsolve_model_refactored.R` (not just a throwaway scratch copy), and
the untouched original still carries all four, unfixed, exactly as
written.

1. **`$INIT @annotated` is not valid syntax for a plain `NAME = value`
   init block.** `$INIT @annotated` expects each line in the
   `NAME : value : description` annotated form; the original's `$INIT`
   lines (`DEPOT_IB = 0`, etc.) have no description field, so mrgsolve
   2.0.1 refuses to parse the block at all (`Error: improper annotation
   format`). Fix: dropped `@annotated` from the `$INIT` header — moot
   after defect 2's fix folded `$INIT` into `$MAIN` anyway, but noted as
   the first error encountered.
2. **`$CMT` (annotated) + `$INIT` jointly redeclare all 18
   compartments** — `Duplicated model names: DEPOT_IB CENT_IB BTK_FREE
   ...`. Same defect class as `chronic-hypothyroidism` (#34),
   `age-related-macular-degeneration` (#36), `type1-diabetes` (#42), and
   `abdominal-aortic-aneurysm` (#44). Fix: the `$INIT` block's 18
   `NAME = value` lines were moved into `$MAIN` as `NAME_0 = value;`
   (mrgsolve's modern idiom) — declares no new compartment, changes no
   numeric value.
3. **A new pattern, not previously logged in exactly this form: the
   pre-existing `$PARAM` symbol `ALC_0` collides with mrgsolve's own
   auto-generated `<CMT>_0` initial-condition symbol** for the
   compartment `ALC`, once defect 2's fix moves `$INIT` into `$MAIN`.
   Closest precedent is `abdominal-aortic-aneurysm` (#44)'s
   `MMP9_0`/`MMP2_0` collision — same root cause (an *incidental*
   collision between an ordinary baseline parameter and mrgsolve's
   `<CMT>_0` idiom, not an author's deliberate-but-broken attempt like
   #42's `Gt_0`). `ALC_0` is a genuine baseline-ALC parameter, used in
   `ALC_pch = (ALC_0 - ALC)/ALC_0*100` and the `PD_flag` calculation —
   its value (`50.0`) happens to equal the compartment `ALC`'s own
   initial value, which is presumably why the original author picked
   that name without noticing the clash. Confirmed reproducing exactly
   as AAA's did:
   ```
   150:9: error: conflicting declaration 'double& ALC_0'
     150 | double& ALC_0 = _A_0_[13];
   117:15: note: previous declaration as 'const double& ALC_0'
     117 | const double& ALC_0 = _THETA_[27];
   ```
   Fix: renamed the `$PARAM` symbol to `ALC_BASELINE` everywhere it was
   used as a parameter (the `$PARAM` declaration itself, and both
   `ALC_pch`/`PD_flag` reads in `$TABLE`) — a pure symbol rename, zero
   numeric change — leaving the `$MAIN` line `ALC_0 = 50.0;` to refer
   only to the compartment's own auto-generated init symbol.
4. **`$CAPTURE` repeats five compartment names already in `$CMT`** —
   `compartment should not be in $CAPTURE: ALC,BM_CLL,LN_CLL,MCL1_ADAPT,
   NK_ACT`. Same defect class as `chronic-hypothyroidism` (#34) and
   `abdominal-aortic-aneurysm` (#44). Fix: the five duplicated lines were
   removed from `$CAPTURE` (compartments are always present in mrgsolve's
   output regardless of `$CAPTURE`, so this changes nothing about what is
   reported — all five are still returned, just not double-declared).

All four fixes applied together were required simultaneously to reach a
buildable model; `POST /model_manifest` on the untouched original (with
none of these fixes) fails with defect 1, then successively with defects
2–4 as each prior one is worked around in isolation — confirmed by
iterating through each error message via the qspserver `mrgsolve_api`
container in that exact order during this refactor.

## Also observed (R-driver-script inconsistency, not a build defect, not fixed): scenario 6 is dosing-identical to scenario 4 despite its "Triplet" label

`dose_events()`'s scenario-6 branch explicitly sets `ev_list[["IB6"]] <-
NULL # no ibrutinib`, while `scenarios[["6_triplet"]]`'s label reads
"Triplet IB+VEN+OBI". Tracing the `%in%` conditions: ibrutinib is only
dosed for `scenario %in% c(1)` (scenario 1) or `scenario == 5`
(scenario 5) — scenario 6 satisfies neither, so despite its label,
scenario 6 never doses ibrutinib and produces byte-identical dosing to
scenario 4 (`VEN` ramp + `OBI` cycles, both via `scenario %in%
c(2,4,5,6)` / `c(3,4,6)`). This is an R-driver-script authoring
inconsistency (a mislabeled scenario), not a DSL/model defect, and not
something this refactor is in scope to fix (the guide's rule is "log
what you find, don't fix it upstream," and the fix here would live in the
original's own `dose_events()`, which stays untouched). Practically: since
scenario 6's dosing is identical to scenario 4's, verifying scenario 4
also verifies scenario 6's actual (as opposed to labeled) behavior.

## Verification

Via the qspserver `mrgsolve_api` container (`POST /model_manifest`, `POST
/run_simulation`) against the untouched original with the four defects
above worked around **in-memory only, not committed to either file**
(same approach as `abdominal-aortic-aneurysm` #44 and others) —
purely so a verification harness could run at all; the tracked original
still contains all four defects, unfixed.

`POST /model_manifest` confirmed both DSLs compile; the refactored model
exposes 51 `$PARAM` entries (the original's 47, plus 4 new: `RTOT_IB`,
`GAMMA_IB`, `GAMMA_VEN`, `GAMMA_OBI` — `RTOT_VEN`/`RTOT_OBI` are renames
of existing `BCL2_tot`/`CD20_0`, not new, and `ALC_0`→`ALC_BASELINE` is a
pure rename, neither adds to the count) and 35 output paths (18 compartments +
11 deduplicated `$CAPTURE` outputs, i.e. the original's own 16
`$CAPTURE` entries minus the 5 that duplicated compartment names, per
defect 4 below + `C_IB`/`C_VEN`/`C_OBI`/`EFFECT_IB`/`EFFECT_VEN`/
`EFFECT_OBI`), matching the original's manifest one-for-one by content
for everything that isn't a rename or one of the six new diagnostics.

Ran all of the original's own `dose_events()` dosing schedules (amounts,
times, `ii`/`addl`, exactly as written in the original's R code) through
both models, 730-day horizon, `delta=12` — matching the original's own
`mrgsim(..., end = end_days*24, delta = 12)` call:

1. **Scenario 1** — Ibrutinib 420 mg QD (`amt=420` into the depot,
   `ii=24`, `addl=729`)
2. **Scenario 2** — Venetoclax 5-week ramp (20→50→100→200→400 mg QD) →
   400 mg QD maintenance
3. **Scenario 3** — Obinutuzumab 6 cycles q28d (C1 split D1=100mg/
   D2=900mg/D15=1000mg, C2–C6=1000mg each)
4. **Scenario 4** — Venetoclax ramp + Obinutuzumab (CLL14 regimen)
5. **Scenario 5** — Ibrutinib + Venetoclax ramp combo

(Scenario 6 not run separately — see "Also observed" above; its dosing is
identical to scenario 4's, already covered.)

One verification-harness-only adjustment, applied identically to both
models so it cannot advantage either side: the original's own ramp
formula (`ceiling((next_time - time)/24) - 1`) evaluates to `addl = -1`
for the final ramp step (`i=5`, where `min(i+1,5)` clamps to the same
index so the time difference is 0) — a value the qspserver API's
`DoseSpec` validator rejects outright (`addl must not be negative`),
whereas mrgsolve's local `ev()` may tolerate or silently clamp it. This
step is redundant regardless (`VEN_main` doses the same 400 mg at the
same `t=672` the very next line in the original's own `dose_events()`),
so it was clamped to `addl=0` (a single dose) purely to build a
verification event set; the delivered `_refactored.R`'s own
`dose_events()` R code is untouched and still contains the original's
`addl=-1`-producing formula verbatim, since that R driver code is never
executed through the API (only the extracted DSL is) and is out of this
refactor's scope.

**Result:**
- **Scenarios 1, 2, 4, 5: exact match, max abs diff = 0.0 and max rel
  diff = 0.0** for every one of the 16 shared `$CAPTURE` outputs
  (`ALC`, `BM_CLL`, `LN_CLL`, `MCL1_ADAPT`, `NK_ACT`, `BTK_OCC_out`,
  `BCL2_OCC_out`, `CD20_OCC_out`, `C_IB_ngmL`, `C_VEN_ngmL`,
  `C_OBI_ugmL`, `BURDEN`, `CR_flag`, `PR_flag`, `MRD_neg`, `ALC_pch`),
  across the full time grid (1462–1468 points depending on scenario).
  Scenario 5 was additionally checked at the renamed-compartment level
  (all 13 PK/receptor compartments, e.g. `DEPOT_IB` vs `GUT_IB`,
  `BTK_FREE` vs `REC_FREE_IB`, `CD20_OCC` vs `COMPLEX_OBI`): also exact,
  max abs diff = 0.0 for all 13.
- **Scenario 3 (Obinutuzumab monotherapy): max abs diff = 1.0e-4**
  on `ALC`/`BM_CLL`/`LN_CLL`/`BURDEN`/`ALC_pch` (all other outputs exact).
  Inspected directly: the mismatch is a single value pair like
  `299.9225` vs `299.9224` — a difference of exactly one unit in the
  API's own 4-decimal-place JSON rounding, at a value sitting on a
  rounding boundary, not a real numeric divergence. This is the same
  floating-point/output-rounding noise floor reported in
  `acute-intermittent-porphyria/aip_refactor_notes.md` ("remaining 4 at
  floating-point noise level, max abs 1e-4"). Per the guide's tolerance
  rule for pure structural reorganization ("anything beyond
  floating-point-scale deviation means a bug"), 1e-4 at the API's own
  serialization precision is within that floor, not a finding.
- `C_IB`, `C_VEN`, `C_OBI`, `EFFECT_IB`, `EFFECT_VEN`, `EFFECT_OBI`
  (no original counterpart to diff against) were sanity-checked for
  plausible ranges in the combo scenario (5): `C_IB` 0–0.01 nM\*,
  `C_VEN` 0–0.199 nM, `EFFECT_IB` 0–0.080, `EFFECT_VEN` 0–0.652,
  `C_OBI`/`EFFECT_OBI` both 0 throughout (correctly, since scenario 5
  never doses obinutuzumab).
  (\*`C_IB`'s low nM range reflects ibrutinib's huge apparent `V1_IB =
  10000 L`, taken verbatim from the original — a modeling choice of the
  original file, not something this refactor introduced or altered.)

All five verification runs and the transient container hiccups
encountered while running them are logged in the working history of this
session; no scenario required the horizon-shortening workaround (the
730-day/`delta=12` grid ran to completion, ~1462–1475 points, well under
the API's default solver step budget).

**Note on the exact delivered text vs. the text actually run above:** all
five scenario runs and the manifest check above were executed against a
scratch copy of the refactored DSL that is identical to the delivered
`cll_mrgsolve_model_refactored.R`'s embedded `cll_model <- '...'` block
**except** that several prose comments in the scratch copy used a
possessive apostrophe (`"ibrutinib's own PK"`, `"original's $INIT
block"`, etc.) that had to be dropped (`"ibrutinib own PK"`, `"original
$INIT block"`) in the delivered file, since the DSL is embedded inside an
R single-quoted string literal and an apostrophe would terminate that
string early. A line-by-line diff of the two texts confirms these are the
**only** differences — every `$PARAM`/`$CMT`/`$MAIN`/`$ODE`/`$TABLE`/
`$CAPTURE` line, every symbol name, and every numeric value is
byte-identical between the two; only comment prose changed. After
finishing the primary verification runs above, an attempt was made to
re-run `/model_manifest` against the exact delivered text as a final
confirmation, but the qspserver `mrgsolve_api` container became
unresponsive (`docker inspect` showed its internal process had exited
while `docker ps` still reported it "Up", and `docker restart`/
`docker compose up --force-recreate` both failed with "tried to kill
container, but did not receive an exit event" — a host/Docker-daemon-level
fault outside this session's ability to fix, not related to the model
DSL). Given the diff-confirmed byte-for-byte equivalence of every
executable line, the verification results above are treated as fully
applicable to the delivered file; a reviewer with a working
`mrgsolve_api` container can re-run `/model_manifest` on the delivered
file's extracted DSL as a mechanical confirmation of this claim.

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`, three
rows under `chronic-lymphocytic-leukemia`: `Ibrutinib` (target BTK,
`C_IB`/`EFFECT_IB`, nM), `OBI` (target CD20, `C_OBI`/`EFFECT_OBI`, mg/L),
`Venetoclax` (target BCL-2, `C_VEN`/`EFFECT_VEN`, nM).
