# GERD PK/PD refactor notes

Source: `gerd/gerd_mrgsolve_model.R` (never edited).
Output: `gerd/gerd_mrgsolve_model_refactored.R` +
`gerd/gerd_mrgsolve_model_refactored.cpp` (bare DSL extraction, byte-identical
to the quoted block in the `.R` sibling -- confirmed by direct string
comparison).

## Confirmed compound identities (from the code, not the census label)

The census's generic "H2RA dose"/"P-CAB dose (PCAB)"/"PPI dose"/"Prokinetic
dose" labels are the classifier's dosing-variable-name artifacts, exactly the
failure mode the guide warns about. The original file's own `$PROB` header
comment and every `$PARAM` description name the real drug for all four
compounds, and the R-side scenario list dose them by name:

| Stem | Real drug | Evidence |
|---|---|---|
| `PPI`  | Omeprazole (default); Esomeprazole (one scenario overrides `F_PPI=0.73, CL_PPI=12`, same `GUT_PPI/CENT_PPI` compartments) | `$PROB` header, `$PARAM` comments, scenario names "Omeprazole 20 mg QD" / "Esomeprazole 40 mg QD" |
| `H2RA` | Famotidine | `$PARAM DOSE_H2RA : 40 : H2RA dose (mg, famotidine equivalent)`, scenario "Famotidine 40 mg BID" |
| `PCAB` | Vonoprazan | `$PARAM DOSE_PCAB : 20 : P-CAB dose (mg, vonoprazan)`, scenario "Vonoprazan 20 mg QD (P-CAB)" |
| `PROK` | Domperidone | `$PARAM DOSE_PROK : 10 : Prokinetic dose (mg, domperidone)`, scenario "...+ Domperidone 10 mg TID" |

No placeholder/generic compound here -- all four are real, named drugs.

## Archetype per compound

All four: **archetype 3 minus peripheral** (depot + central, linear
elimination, no peripheral compartment) -- matches the original exactly for
H2RA/P-CAB/Prokinetic. PPI additionally carries the original's own
`CYP2C19` genotype multiplier (scales `KA*F` up and `CL/V1` down by the same
factor) -- kept byte-for-byte, renamed variables only, no new pharmacology
added or removed.

- `PPI` (omeprazole/esomeprazole): `GUT_PPI`/`CENT_PPI` (was `PPI_GUT`/
  `PPI_CENT`), `V1_PPI` (was `V_PPI`). `EFFECT_PPI` is a rename of `INH_PPI`
  (already `Emax*C^gamma/(EC50^gamma+C^gamma)`; `EC50_PPI`=`IC50_PPI`,
  `GAMMA_PPI`=`HILL_PPI`=1.5, `EMAX_PPI`=1.0 new-explicit since the original
  had no separate scaling factor). Rename only, no refit.
- `H2RA` (famotidine): `GUT_H2RA`/`CENT_H2RA` (was `H2RA_GUT`/`H2RA_CENT`),
  `V1_H2RA` (was `V_H2RA`). `EFFECT_H2RA` renamed from `INH_H2RA`
  (`EMAX_H2RA*C/(EC50_H2RA+C)`, already Emax-scaled in the original;
  `EMAX_H2RA` keeps its original name+value, `EC50_H2RA`=`IC50_H2RA`,
  `GAMMA_H2RA`=1.0 new-explicit). Rename only.
- `PCAB` (vonoprazan): `GUT_PCAB`/`CENT_PCAB` (was `PCAB_GUT`/`PCAB_CENT`),
  `V1_PCAB` (was `V_PCAB`). `EFFECT_PCAB` renamed from `INH_PCAB` (plain
  `C/(IC50+C)` ratio; `EC50_PCAB`=`IC50_PCAB`, `EMAX_PCAB`=1.0 and
  `GAMMA_PCAB`=1.0 new-explicit). Rename only.
- `PROK` (domperidone): `GUT_PROK`/`CENT_PROK` (was `PROK_GUT`/`PROK_CENT`),
  `V1_PROK` (was `V_PROK`). **Two** named effects sharing one compound and
  one promoted `EC50_PROK`: `EFFECT_PROK_LES` (was the inline
  `K_PROK_LES*CP_PROK/(0.005+CP_PROK)` expression feeding `LES_P`) and
  `EFFECT_PROK_EMP` (was the inline `K_PROK_EMP*CP_PROK/(0.005+CP_PROK)`
  expression feeding `EMP_RATE`). `EMAX_PROK_LES`=`K_PROK_LES`(5.0),
  `EMAX_PROK_EMP`=`K_PROK_EMP`(0.30), `EC50_PROK`=0.005 (promoted from a
  literal hardcoded identically in *three* places in the original -- see
  below), `GAMMA_PROK`=1.0 new-explicit. Rename/promotion only, no refit
  (same shared-EC50 shape the original already had, just not named).

Each compound keeps its own separate `EFFECT_<STEM>`; they are combined only
where the disease equations actually read them (`ACID_INH` combines
`EFFECT_PPI`/`EFFECT_PCAB`/`EFFECT_H2RA`; `LES_P`/`EMP_RATE` each read only
their own prokinetic effect), per the guide's multi-drug rule.

## Real duplicate-concentration-site defect found -- corrects the census for 3 of 4 compounds

The census only flagged `PPI` as "Normalize duplicate concentration sites,
then redirect" (`H2RA`/`P-CAB`/`Prokinetic` were labelled "clean single
site"). The actual code shows **all four** compounds independently
compute their own concentration TWICE under the identical name:

```r
$MAIN
double CP_PPI  = PPI_CENT  / V_PPI;   # ... same for CP_H2RA, CP_PCAB, CP_PROK
...
$TABLE
capture CP_PPI   = PPI_CENT / V_PPI;  # ... same for CP_H2RA, CP_PCAB, CP_PROK
```

This is not just cosmetic duplication -- confirmed live against qspserver
(see "Upstream build defects" below), it is a genuine mrgsolve 2.0.1 build
failure: `$MAIN`'s `double CP_PPI` and `$TABLE`'s `capture CP_PPI = ...`
share the same compiled scope, so the second declaration is a hard
`redefinition of capture ::CP_PPI` compile error. Census census row for
PPI happened to be right about *something* being wrong here, but for the
wrong reason (it called it a duplicate-derivation semantic issue; it is
also a literal build defect), and it missed that H2RA/P-CAB/Prokinetic have
the exact same defect. **Census updated for all four rows** (see below).

`LES_pressure` in the original's own `$TABLE` is a **third**, independent
re-derivation of the prokinetic LES effect (`LES_P0 + K_PROK_LES*CP_PROK/
(0.005+CP_PROK)`, the 0.005 literal baked in yet again) -- a plain semantic
duplicate (different name from `PROK_EFF_LES`, so it did not itself trigger
a redefinition error), also consolidated.

### Fix in the refactored file

`C_PPI`/`C_H2RA`/`C_PCAB`/`C_PROK` and `EFFECT_PPI`/`EFFECT_PCAB`/
`EFFECT_H2RA`/`EFFECT_PROK_LES`/`EFFECT_PROK_EMP` are now computed **once**,
in `$MAIN` (the same block the original used them for its own PD -- per the
guide's "keep a calculation in the block the original used it in" rule, since
`$MAIN`-computed values feed directly into `$ODE`'s `dxdt_PUMP_ACT`/
`dxdt_PUMP_INH` and other terms in the original too). They are never
recomputed in `$TABLE`; `$TABLE`'s `LES_pressure` now reads
`LES_P0 + EFFECT_PROK_LES` instead of re-deriving it. All nine are exposed
via a bare `$CAPTURE` block (not `$PARAM` -- see below) so they remain
genuine downstream-discoverable outputs (confirmed present in
`/model_manifest`'s `outputPaths`, not silently dropped by removing the
`$TABLE` recompute).

## Upstream build defects found (pre-existing, unrelated to refactor scope)

Confirmed live against the qspserver mrgsolve API (`POST /model_manifest`)
on the **unmodified** original `gerd_mrgsolve_model.R`'s own extracted DSL
block:

1. **`$CMT @annotated` + a separate `$INIT` block jointly redeclaring the
   same 17 compartment names.** Compile error:
   `invalid class "mrgmod" object: Duplicated model names: PPI_GUT PPI_CENT
   H2RA_GUT H2RA_CENT PCAB_GUT PCAB_CENT PROK_GUT PROK_CENT PUMP_INACT
   PUMP_ACT PUMP_INH ACID_RATE GAS_pH AET MUC_DMG MUC_HEAL SYM_SCORE
   BE_RISK`. This is the exact "`$CMT`+`$INIT` jointly redeclaring
   compartments" defect class named in `FORK_WORKFLOW_GUIDE.md`. **Fixed
   in the refactored sibling only** by dropping `$INIT` entirely and
   setting the same non-zero initial values via mrgsolve's own
   `<cmt>_0 = value;` idiom inside `$MAIN` (this guide's own Archetype 4
   example uses the identical idiom) -- no numeric change; the PK
   compartments still start at 0 exactly as before, the 10 disease-state
   compartments start at the same values the original's `$INIT` gave them.
2. **`$MAIN`/`$TABLE` double-declaration of `CP_PPI`/`CP_H2RA`/`CP_PCAB`/
   `CP_PROK`**, described above -- confirmed as its own separate compile
   error once defect (1) above is fixed in isolation:
   `error: redefinition of 'capture {anonymous}::CP_PPI'` /
   `note: 'double {anonymous}::CP_PPI' previously declared here` (times 4,
   one per compound). Fixed in the refactored sibling as part of the
   normal refactor (single canonical `C_<STEM>` computed once, see above)
   -- this is simultaneously the required PK/PD refactor work AND the
   build-defect fix, so no separate throwaway patch was needed for this
   one.

Both defects logged as **`translations/UPSTREAM_ISSUES.md` #<TBD, see
current tail of that file at time of commit>** -- entry text: *"gerd:
`$CMT @annotated` + `$INIT` jointly redeclare all 17 compartment names
(`Duplicated model names` build failure); `$MAIN` `double CP_<STEM>` +
`$TABLE` `capture CP_<STEM> = ...` (PPI/H2RA/PCAB/PROK, 4x) is a second,
independent redefinition-of-capture build failure once (1) is patched in
isolation. Neither defect is specific to any one compound's own PK/PD;
both confirmed live via qspserver `/model_manifest` against the
unmodified original. Not fixed upstream; fixed syntax-only in
`gerd_mrgsolve_model_refactored.R`/`.cpp` only."*

Neither fix changes any parameter value, dosing, or model equation --
confirmed by the manifest's `parameters` list (identical names/values to
the original's own `$PARAM` block, aside from the renames documented above)
and by verification below.

## `C_<STEM>`/`EFFECT_<STEM>` are `$MAIN`-local + `$CAPTURE`, not `$PARAM`

Matches the corpus's established, disclosed exception (same precedent as
`acute-intermittent-porphyria`, `sepsis`, `age-related-macular-degeneration`):
declaring these as `$PARAM` would make them read-only where `$MAIN`
recomputes them every step, a confirmed mrgsolve 2.0.1 conflict. They are
`double`s declared once in `$MAIN` and listed in an annotated `$CAPTURE`
block, which the discoverability contract (`double C_<STEM> = <expr>;` +
a same-stem `EC50_<STEM>` parameter) does not require to be `$PARAM` --
confirmed present in both places by direct grep of the finished file:

```
double C_PPI  = CENT_PPI  / V1_PPI;      double C_H2RA = CENT_H2RA / V1_H2RA;
double C_PCAB = CENT_PCAB / V1_PCAB;     double C_PROK = CENT_PROK / V1_PROK;
EC50_PPI / EC50_H2RA / EC50_PCAB / EC50_PROK   (all in $PARAM)
```

## Verification

**STATUS: COMPLETE. Exact match on every shared output, every scenario.**

### Static checks (all passed)

- `Rscript -e 'parse("gerd_mrgsolve_model_refactored.R")'` succeeds (35
  top-level expressions parsed), no error.
- No straight apostrophes anywhere in the quoted DSL body (checked
  programmatically, 0 matches after the fix pass); exactly one `<stem>_code
  <- '` string open and one closing `'` immediately followed by
  `mod <- mcode(...)`.
- `gerd_mrgsolve_model_refactored.cpp` is byte-identical to the quoted
  string body in the `.R` sibling (direct Python string comparison,
  `True`).
- `POST /model_manifest` against the finished `gerd_mrgsolve_model_refactored.cpp`
  succeeds (HTTP 200) and returns exactly the expected `outputPaths`:
  all 18 renamed compartments, `C_PPI`/`C_H2RA`/`C_PCAB`/`C_PROK`,
  `EFFECT_PPI`/`EFFECT_PCAB`/`EFFECT_H2RA`/`EFFECT_PROK_LES`/
  `EFFECT_PROK_EMP`, and the 10 unchanged disease-side `$TABLE` outputs
  (`ACID_INH_pct`, `pH`, `AET_pct`, `DMG`, `HEAL`, `SYM`, `PUMP_active`,
  `PUMP_inhibited`, `LES_pressure`, `Barrett`) -- confirming both build
  defects are actually fixed and every promised covariate/effect is
  genuinely declared and capturable, not just present in the DSL text.
- A **reference-fixed copy of the original** (throwaway scratch `.cpp`,
  never committed -- same two syntax-only fixes as above: `$INIT`
  dropped in favour of `<cmt>_0` in `$MAIN`, and the `$TABLE` `CP_<STEM>`
  redeclarations replaced by a bare `$CAPTURE CP_PPI CP_H2RA CP_PCAB
  CP_PROK` block referencing the `$MAIN`-declared doubles directly) also
  compiles cleanly via `/model_manifest` (HTTP 200, correct
  `outputPaths` including `CP_PPI`/`CP_H2RA`/`CP_PCAB`/`CP_PROK`). This is
  the ground-truth baseline used for the side-by-side run below.
- Solver-step-budget check (per the guide's known API limitation):
  `end=1344` (the original's own full 56-day window) times out
  (`lsoda ... 20000 steps ... consider increasing maxsteps`) on the
  reference-fixed original at the default `maxsteps=20000`, confirmed
  independent of dosing (reproduces even for the zero-dose Control
  scenario, so it is inherent to the disease-state system's stiffness,
  not a refactor artifact). `end=1200` (50 days) was found to run
  reliably; `end=1250`+ intermittently hits a genuine `lsoda` corrector
  failure (`step size h_ = 2.4e-8 ... corrector convergence failed`)
  suggestive of real numerical stiffness near the 52-day mark (plausibly
  the `MUC_DMG`/`MUC_HEAL` hard floor/ceiling resets in `$ODE`) --
  **verification window shortened to `end=1200, delta=0.5` (50 days) for
  all six scenarios**, per the guide's explicit allowance to shorten
  rather than declare a mismatch on a solver-budget limit.

### Infrastructure blocker encountered and resolved (disclosed for the record)

`POST /model_manifest`/`/run_simulation` failed server-wide, twice, during
this session's testing, both times with a corrupted mrgsolve build-cache
file (`mrgmod_cache.RDS`) reported via `readRDS()` (`"embedded nul in
string"` the first time, `"ReadItem: unknown type 32, perhaps written by
later version of R"` the second time -- consistent with two different
concurrent sessions' writes colliding on the same shared, non-content-keyed
cache file for the qspserver mrgsolve API's `inline`-model project
directory). Both times reproduced with a trivial, completely unrelated
one-compartment model, confirming it was never specific to the GERD model
or this refactor. `/health` continued to report `"status":"healthy"`
throughout both incidents. The orchestrating session ran
`docker exec qspserver-mrgsolve_api-1 rm -f
/soloc/mrgsolve-so-2.0.1-x86_64-pc-linux-gnu/{inline,friberg_ro,pk1cmt,pkpd_indirect}/mrgmod_cache.RDS`
to clear the stale cache across all four affected project directories,
confirmed with a trivial test model that `/model_manifest` compiled again,
and the verification below was completed immediately afterward with no
further recurrence. (This session's own sandbox denied the same `docker
exec` command when attempted directly -- surfaced to the orchestrator
rather than worked around, and the orchestrator ran it instead.)

### Side-by-side numeric comparison -- RESULT: exact match, all scenarios

All six of the original's own named scenarios run through both the
reference-fixed original and the refactored file, identical dosing,
`end=1200, delta=0.5` (2402 points; 2403 for the combo scenario, which adds
one extra event-aligned row for the PROK q8h dosing -- time grids confirmed
identical between the two models in every scenario). Every one of the 14
outputs the two models share (`CP_PPI`/`C_PPI`, `CP_H2RA`/`C_H2RA`,
`CP_PCAB`/`C_PCAB`, `CP_PROK`/`C_PROK`, `ACID_INH_pct`, `pH`, `AET_pct`,
`DMG`, `HEAL`, `SYM`, `PUMP_active`, `PUMP_inhibited`, `LES_pressure`,
`Barrett`) compared at every time point (≈33,600 value pairs total):

| Scenario | Dosing | Points | Max abs diff | Max rel diff | NaNs (either side) |
|---|---|---|---|---|---|
| Control | none (amt=0 into GUT_PPI) | 2402 | 0.0 | 0.0 | 0 |
| Omeprazole 20 mg QD | GUT_PPI 20mg q24h x56 | 2402 | 0.0 | 0.0 | 0 |
| Esomeprazole 40 mg QD | GUT_PPI 40mg q24h x56, F_PPI=0.73/CL_PPI=12 | 2402 | 0.0 | 0.0 | 0 |
| Vonoprazan 20 mg QD | GUT_PCAB 20mg q24h x56 | 2402 | 0.0 | 0.0 | 0 |
| Famotidine 40 mg BID | GUT_H2RA 40mg q12h x112 | 2402 | 0.0 | 0.0 | 0 |
| Eso 40mg QD + Domperidone 10mg TID | GUT_PPI 40mg q24h + GUT_PROK 10mg q8h, F_PPI=0.73/CL_PPI=12 | 2403 | 0.0 | 0.0 | 0 |

**Every scenario, every shared output: max abs diff = 0.0, max rel diff =
0.0.** This is the expected result for a pure structural
reorganization/rename with no Hill-fitting (guide's Archetype 1-3
tolerance) -- confirms the refactor changed nothing numeric: the two
build-defect fixes are genuinely syntax-only, the four compounds'
concentration-site consolidation (single `$MAIN`-computed `C_<STEM>`
instead of the original's duplicate `$MAIN`+`$TABLE` declarations)
reproduces the original's own intended values exactly, and the
`LES_pressure` redirect to `EFFECT_PROK_LES` is a true algebraic identity
with the original's own three-times-duplicated formula, not an
approximation.

Verification script and raw data (session-local scratchpad, not
committed): `gerd_verify.py` (six-scenario harness), `gerd_verify_raw.json`
(full per-scenario time series from both models), `gerd_verify_summary.json`
(the per-output max-abs/max-rel table above), `gerd_ref_fixed2.cpp`
(reference-fixed original), `gerd_refactored.cpp` (extraction, identical
to the committed `.cpp`).

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`, the
`gerd | H2RA dose`, `gerd | P-CAB dose (PCAB)`, `gerd | PPI dose`,
`gerd | Prokinetic dose` rows -- corrected compound names (Famotidine,
Vonoprazan, Omeprazole/Esomeprazole, Domperidone), all four reclassified to
"Normalize duplicate concentration sites, then redirect" (census had only
PPI in that category; H2RA/P-CAB/Prokinetic share the identical
`$MAIN`+`$TABLE` duplicate-declaration defect, confirmed live, see above).
