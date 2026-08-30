# Refactor notes — `myotonic-dystrophy/dm1_mrgsolve_model.R`

Only **Mexiletine** was in scope, per its existing row in
`driver-patches/data/compound_perturbation_census.md` (classified
"Redirect concentration (clean single site)"). **ASO** (the file's second
compound) was left completely untouched — no renaming, no restructuring,
not even a variable rename — since it has no census row here and was not
part of this task.

## Archetype

**Mexiletine — Archetype 3** (depot + central + peripheral, linear),
matched exactly, no deviation:

```
dxdt_GUT_MEX  = -KA_MEX * GUT_MEX;
dxdt_CENT_MEX = KA_MEX * GUT_MEX
                - (CL_MEX / V1_MEX) * CENT_MEX
                - (Q_MEX / V1_MEX) * CENT_MEX
                + (Q_MEX / V2_MEX) * PERI_MEX;
dxdt_PERI_MEX = (Q_MEX / V1_MEX) * CENT_MEX - (Q_MEX / V2_MEX) * PERI_MEX;
```

## Renaming

| Original | Refactored |
|---|---|
| `MEX_GUT` | `GUT_MEX` |
| `MEX_CENT` | `CENT_MEX` |
| `MEX_PERI` | `PERI_MEX` |
| `Ka_mex` | `KA_MEX` |
| `F_mex` | `F_MEX` |
| `Vc_mex` | `V1_MEX` |
| `Vp_mex` | `V2_MEX` |
| `CL_mex` | `CL_MEX` |
| `Q_mex` | `Q_MEX` |
| `MW_mex` | `MW_MEX` |
| `fu_mex` | `FU_MEX` |
| `IC50_mex_nav` | `EC50_MEX` |
| `hill_mex_nav` | `GAMMA_MEX` |
| `Cp_mex` | `C_MEX_TOTAL` |
| `Cp_mex_free` | `C_MEX` |
| `nav_block` | `EFFECT_MEX` |

All parameter *values* are copied verbatim from the original — nothing
invented or defaulted. Added (new, not present in the original, purely to
complete the named Hill interface): `EMAX_MEX = 1.0`.

**A pre-existing dead parameter, preserved as-is, not fixed:** `F_mex`
(now `F_MEX`, "Bioavailability" = 0.85) is declared in `$PARAM` but is
**never applied anywhere in the original's own `$ODE`** — the absorption
term is `dxdt_MEX_CENT = Ka_mex * MEX_GUT - ...` with no `F_mex`
multiplier, and the R-side dosing event (`ev(cmt = "MEX_GUT", amt =
mex_dose, ...)`) doesn't apply it either. So the original's own simulated
behaviour is equivalent to `F=1` regardless of the declared 0.85 value.
Per the guide's "verify against the original, mechanically" rule, the
refactored ODE **reproduces this exactly** (no `F_MEX` multiplication
added) rather than "fixing" it to what the parameter's name/comment
implies — adding it now would silently change the compound's simulated
exposure by 15%, which is exactly the kind of unauthorized behavioural
change this refactor must not introduce. This is a modeling defect in the
original (not a build defect), logged separately below per the shared
"log what you find, don't fix it upstream" rule.

Two mechanical consequences of renaming the compartments, both required
for the model to still mean the same thing (not scope creep):
- The R-side `build_events()` helper's dosing target was updated from
  `cmt = "MEX_GUT"` to `cmt = "GUT_MEX"`.
- The R-side `plot_mex_pk` plot, which reads the captured concentration
  directly by name, was updated from `aes(Day, Cp_mex, ...)` to
  `aes(Day, C_MEX_TOTAL, ...)` — same data, renamed column.

Compartment declaration order (and therefore 1-based compartment number)
is unchanged: `GUT_MEX` is still compartment 1, `CENT_MEX` compartment 2,
`PERI_MEX` compartment 3, exactly as `MEX_GUT`/`MEX_CENT`/`MEX_PERI` were.

## The exposed concentration: `C_MEX`

The Nav1.4-blockade Hill term reads the **free** plasma concentration
(`Cp_mex_free` in the original), not the total. That is therefore the
single PD-facing signal, renamed `C_MEX`. The total (protein-bound+free)
concentration, `Cp_mex`, is kept as a diagnostic capture only — renamed
`C_MEX_TOTAL` — matching the original's own practice of capturing both,
but making explicit that only `C_MEX` is the point an external covariate
would later substitute into.

```
double C_MEX_TOTAL = (CENT_MEX / V1_MEX) * 1000.0 / MW_MEX;  // μg/mL → μM
double C_MEX = C_MEX_TOTAL * FU_MEX;                          // free, PD-facing
```

## Hill interface: rename, not a fit

The original's Nav1.4 blockade was already written as a plain
`pow(C,gamma)/(pow(EC50,gamma)+pow(C,gamma))` ratio with an implicit
`Emax=1` (no separate multiplier in front):

```
double nav_block = MEX_ON * pow(Cp_mex_free, hill_mex_nav) /
                   (pow(IC50_mex_nav, hill_mex_nav) + pow(Cp_mex_free, hill_mex_nav));
```

Renamed to `EFFECT_MEX = MEX_ON * EMAX_MEX * pow(C_MEX, GAMMA_MEX) /
(pow(EC50_MEX, GAMMA_MEX) + pow(C_MEX, GAMMA_MEX))` with `EMAX_MEX = 1.0`
added explicitly (multiplying by 1 changes nothing). No curve-fitting was
performed or needed — this is a pure rename of an already-Hill-shaped
term, per the guide's "already this shape" rule.

`MEX_ON` (the scenario on/off flag) is **not** part of the named Hill
interface — it is a pre-existing scenario switch (also present in the
original, multiplying the whole ratio), preserved unchanged and outside
the `EMAX`/`EC50`/`GAMMA` triple. It is functionally redundant with actual
dosing in every one of the original's own scenarios (if `MEX_ON=0` no
mexiletine is ever dosed, so `C_MEX` is already 0), but removing it would
be an unrequested behavioural change, so it stays exactly as the original
had it.

`EFFECT_MEX` feeds two disease targets, exactly as `nav_block` did:
`myotonia_target` (`1.0 - 0.80 * EFFECT_MEX`) and `QTc_mex_effect = 15.0 *
EFFECT_MEX`.

## Two pre-existing build defects found while verifying (not fixed
upstream, logged)

Confirmed on the **untouched original**, via `POST /model_manifest`
against the original file's own extracted DSL — neither introduced by
this refactor:

1. **`$CMT` and `$INIT` jointly redeclare all 19 compartments.** `$CMT`
   names all 19 compartments with inline comments, then a separate
   `$INIT` block assigns each of the same 19 names a starting value, in
   identical order. mrgsolve 2.0.1 refuses to build any model with
   duplicated compartment declarations:
   ```
   Error in validObject(.Object) :
     invalid class "mrgmod" object: Duplicated model names: MEX_GUT
     MEX_CENT MEX_PERI ASO_PLASMA ASO_MUSCLE ASO_NUCL CUG_FOCI MBNL1_FREE
     CUGBP1_ACT CLCN1_FETAL SERCA_FETAL INSR_FETAL MYOTONIA GRIP_STR
     MUSCLE_MASS PR_INT QTc_INT HOMA_IR FVC_PCT
   ```
2. **`$PARAM`'s `HOMA_IR_0` collides with mrgsolve's own auto-reserved
   `<CMT>_0` initial-value symbol** for the `HOMA_IR` compartment. This
   defect only surfaces once defect 1 is worked around:
   ```
   155:9: error: conflicting declaration 'double& HOMA_IR_0'
   131:15: note: previous declaration as 'const double& HOMA_IR_0'
   ```
   Checked: no other `$PARAM` name in this file collides with any of the
   other 18 compartments' `<CMT>_0` form (`Muscle_mass_0` and `FVC_0` are
   mixed/short-case names, not exact matches to `MUSCLE_MASS_0`/`FVC_PCT_0`).

Per the guide's settled policy for this situation ("When the original
doesn't compile at all"), both fixes are syntax-only and non-numeric, and
were applied **directly to the delivered
`dm1_mrgsolve_model_refactored.R`** (not just a scratch copy), while
`dm1_mrgsolve_model.R` was left completely untouched:

1. The `$CMT` block was removed; its per-compartment documentation was
   folded into a plain comment immediately above `$MAIN`. `$INIT` alone
   now declares and initializes every compartment — standard mrgsolve
   idiom, identical values, identical order, identical 1-based
   compartment numbers.
2. `HOMA_IR_0` was renamed `HOMA_IR_BASELINE`, and its one read site in
   `$MAIN` (`HOMA_IR_target = HOMA_IR_BASELINE * (1.0 + 1.5*(INSR_FETAL -
   INSR_fetal_norm))`) updated to match.

Logged as entry **#65** in `translations/UPSTREAM_ISSUES.md`.

## Verification

Per the guide's mandatory protocol: ran the original file's own dosing
regimen (mexiletine 300 mg TID / q8h, the only dose amount the original's
own `build_events()` ever actually issues — see below) through both the
(defect-patched, for buildability only) original and the refactored
model via the qspserver `mrgsolve_api` service (`POST /model_manifest`,
`POST /run_simulation`).

**Window shortened for verification-request convenience, not because of
a solver-step-count problem.** The original's own scenarios run a full
365-day (8760 h) horizon at `delta=6`. Verification here used a 30-day
(720 h) window at the same `delta=6` — chosen because most of this
model's disease-state rate constants are very slow (e.g.
`k_muscle_loss=0.0002/h`, `k_PR_prog=0.0001/h`) and a 30-day window
already fully exercises the fast mexiletine PK/PD pathway (`C_MEX` reaches
near-steady-state within ~3 days; `MYOTONIA`/`QTc_INT` respond within the
window) while keeping requests small. No step-count error or timeout was
encountered at either window length; this was a convenience choice, not
a workaround, and is disclosed here per the guide's requirement to be
honest about any shortening.

**Scenarios run** (compartment 1 = `GUT_MEX`/`CENT_MEX`, unchanged by the
rename since compartment order is unchanged):

1. **Natural history** — no dosing, `MEX_ON=0` (default). Sanity check
   that all 19 states hold flat at their `$INIT` baseline.
2. **Mexiletine 300 mg TID (q8h)** — `amt=300, cmt=1, ii=8, addl=89` (the
   original's own `build_events()` always doses 300 mg regardless of
   which scenario label is applied — see "other finding" below), with
   `MEX_ON=1` set explicitly as a parameter override (the original's own
   `idata` mechanism for turning the switch on per scenario).

Every shared output was compared point-by-point across the full time
grid for both scenarios: `MYOTONIA`, `GRIP_STR`, `MUSCLE_MASS`, `PR_INT`,
`QTc_INT`, `HOMA_IR`, `FVC_PCT`, `CLCN1_FETAL`, `SERCA_FETAL`,
`INSR_FETAL`, `CUG_FOCI`, `MBNL1_FREE`, `CUGBP1_ACT`, `myotonia_target`,
`HOMA_IR_target`, `FVC_target`, and the three renamed
Mexiletine-specific outputs (`Cp_mex`/`C_MEX_TOTAL`,
`Cp_mex_free`/`C_MEX`, `nav_block`/`EFFECT_MEX`).

**Result: exact match, max abs diff = 0.0 for every output, both
scenarios**, both against the intermediate tested DSL and, as a final
check, against the DSL text extracted directly from the delivered
`dm1_mrgsolve_model_refactored.R` (re-run after the file was finished, to
confirm the delivered artifact itself — not just a working copy —
verifies). This is the expected outcome per the guide's tolerance table
for a pure structural reorganization with rename-only Hill terms (no
fitting performed, none needed).

Sanity check (not required for verification, but confirms `EFFECT_MEX`
is wired correctly): with `MEX_ON=1` and the 300 mg TID regimen,
`EFFECT_MEX` rises from 0 to ~0.30 by day 3 of the 30-day window while
`MYOTONIA` falls from its 5.50 baseline to ~1.77 and `QTc_INT` falls
toward ~450.8 ms — directionally consistent with the model's own
Logigian/MELT-trial calibration notes (mexiletine reduces myotonia within
days, shortens QTc).

## Other finding (not a build defect, not fixed, disclosed only)

The original's R-side `build_events(mex_on, aso_on, gt_on, end_days)`
helper hardcodes `mex_dose <- mex_on * 300` regardless of which scenario
calls it — so **"Mexiletine 200 mg TID" (Scenario 2) and "Mexiletine 300
mg TID (MELT dose)" (Scenario 3) both actually dose 300 mg TID** in the
original's own `sim_all` run; the "200 mg" label is never realized as an
actual 200 mg dose anywhere in the simulated output. This is a
pre-existing authoring bug in the R wrapper, unrelated to mexiletine's
PK/PD math and unrelated to any of this refactor's renaming, reproduced
faithfully (not fixed) in the delivered `_refactored.R`. Not logged to
`UPSTREAM_ISSUES.md` as a numbered entry since it is an R-side scenario
labeling issue, not a DSL build defect or a discrepancy this refactor's
verification needed to resolve — flagged here for visibility only.

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`, the
`myotonic-dystrophy` row for Mexiletine.
