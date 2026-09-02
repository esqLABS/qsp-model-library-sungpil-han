# Essential hypertension (EH) — PK/PD refactor notes

Source: `eh_mrgsolve_model.R` -> `eh_mrgsolve_model_refactored.R`.
Guide: `FORK_WORKFLOW_GUIDE.md` Part 2. Census row:
`driver-patches/data/compound_perturbation_census.md`, `^| essential-hypertension |`.

## Compound identity (confirmed against the code and `README.md`, not just the census label)

The census lists these five rows under class-abbreviation labels (`ACEI`,
`ARB`, `BB`, `CCB`, `HCTZ`). All five are confirmed, real, specific
molecules from the model's own `$PARAM` comments and dosing scenarios
(`README.md`'s "주요 약물 표적" section lists the same class with several
example drugs each, e.g. "ACE억제제 (라미프릴, 에날라프릴)"; the code picks one
specific member of each class and gives it explicit PK):

| Census label | Confirmed compound | Evidence |
|---|---|---|
| ACEI | **Ramipril** (dosed prodrug) / **ramiprilat** (the active metabolite the ODE compartments actually track) | `$PARAM` comments name "ramipril"/"ramiprilat" on every ACEI line; dosing scenario 2 is literally `"S2: ACEI (Ramipril 10 mg)"`; calibration reference is Breslin et al. (ramipril PK) |
| ARB | **Losartan** (dosed prodrug) / **EXP3174** (the active metabolite the ODE compartments actually track) | `$PARAM` comments name "losartan"/"EXP3174"; scenario 3 is `"S3: ARB (Losartan 100 mg)"`; calibration reference is Gottwald et al. |
| BB | **Bisoprolol** | `$PARAM` comments name "bisoprolol" throughout; scenario 5 is `"S5: BB (Bisoprolol 10 mg)"`; calibration reference is Leopold et al. |
| CCB | **Amlodipine** | `$PARAM` comments name "amlodipine" throughout; scenario 4 is `"S4: CCB (Amlodipine 10 mg)"`; calibration reference is Faulkner et al. |
| HCTZ | **Hydrochlorothiazide** | already a specific real drug, not a class label; calibration reference is Beermann |

No compound in this file is a bespoke/generic class-level placeholder —
all five are named, specific real molecules with literature-cited PK
parameters.

## Archetype per compound

All five compounds are **archetype 2 (no depot, two compartments,
linear)**, except HCTZ, which is **archetype 1 (no depot, single
compartment, linear elimination)** — matching the original's own `$CMT`
block exactly (`HCTZ_C` has no peripheral compartment; the other four each
have a `_C`/`_P` pair). No compound needed a depot/absorption compartment
because **the original itself has none for any of the five** (see finding
below) — nothing was removed or added to fit the archetype; the rewrite
matches the original's actual structure.

Renames applied per the naming convention: `<STEM>_C -> CENT_<STEM>`,
`<STEM>_P -> PERI_<STEM>` (HCTZ: `HCTZ_C -> CENT_HCTZ`, no `PERI_HCTZ`).
`V1_<STEM>`, `V2_<STEM>`, `CL_<STEM>`, `Q_<STEM>` already matched the
convention in the original and are unchanged. All parameter *values* are
copied verbatim from the original.

## Found in passing: `KA_<STEM>`/`F_<STEM>`/`DOSE_<STEM>` are declared but never wired into the model

For every one of the 5 compounds, the original's `$PARAM` block declares an
absorption rate constant (`KA_ACEI`, `KA_ARB`, `KA_CCB`, `KA_BB`,
`KA_HCTZ`), a bioavailability fraction (`F_ACEI`, ... `F_HCTZ`), and a
`DOSE_<STEM>` on/off flag — but **none of these fifteen parameters is ever
referenced anywhere in `$MAIN`/`$ODE`/`$TABLE`**. There is no depot
compartment and no first-order absorption process anywhere in the file:
every drug is dosed as a direct bolus into its central compartment.
Bioavailability is instead applied *externally*, once, as a hardcoded
numeric literal that happens to equal `F_<STEM>`'s value, inside the R-side
dosing helper in section 3 of the file (e.g. `make_dose("ACEI_C", 10 *
0.28 * 1e3 * 1e-3)` — the `0.28` duplicates `F_ACEI`'s value as a magic
number rather than referencing the parameter). This is why the archetype
is 2 (central+peripheral, no depot) rather than 3 (depot+central(+
peripheral)) for every compound, despite `KA_<STEM>` existing in `$PARAM`.

Not a build defect (the file compiles and runs fine either way) and out of
scope to fix as part of a PK-reorganization refactor, so **`KA_<STEM>`,
`F_<STEM>`, and `DOSE_<STEM>` are preserved, unused, with their original
values, in the refactored file** — for parameter-value fidelity, not
because they do anything. Not logged as a separate `UPSTREAM_ISSUES.md`
entry (harmless dead parameters, not a defect that changes behavior or
blocks anything); disclosed here per the guide's "log what you find" rule.

## The Hill interface (rename, not a fit)

Every one of the original's 5 effect terms was already a plain Emax=1,
gamma=1 ratio (`C/(C+IC50)`), computed once in `$ODE` and used directly —
textbook case for a pure rename per the guide. Pulled out as explicit
`EMAX_<STEM>` (=1), `EC50_<STEM>`, `GAMMA_<STEM>` (=1) parameters, with the
Hill term itself written as
`EMAX*pow(C,GAMMA)/(pow(EC50,GAMMA)+pow(C,GAMMA))`:

| Stem | Original effect var | Original IC50 name/value | New EC50 name |
|---|---|---|---|
| ACEI | `ACE_inhib` | `IC50_ACEI = 0.005` (named param) | `EC50_ACEI` |
| ARB  | `AT1R_block` | `IC50_ARB = 0.02` (named param) | `EC50_ARB` |
| BB   | `BB_block` | `IC50_BB = 0.10` (named param) | `EC50_BB` |
| CCB  | `VGCC_block` | `0.003` (inline literal, **not** a named parameter in the original) | `EC50_CCB` (new named parameter, same value) |
| HCTZ | `NCC_inhib` | `0.02` (inline literal, **not** a named parameter in the original) | `EC50_HCTZ` (new named parameter, same value) |

Renamed to `EFFECT_<STEM>` throughout the RAAS/SNS/NO/TPR/cardiac/plasma
volume ODEs that read them (e.g. `ACE_activity = 1 - EFFECT_ACEI`,
`AngII_feedback = 1 + 0.8*EFFECT_ARB`, `TPR_CCB_effect = ALPHA_CCB *
EFFECT_CCB`), each compound's effect kept as its own separate named term —
none collapsed into a shared expression, per the guide.

## Normalizing duplicate concentration sites (the census's own classification for this file)

This is the actual defect the census's "Normalize duplicate concentration
sites, then redirect" label refers to. `$TABLE` in the original independently
recomputed `BB_C/V1_BB` **twice** (once inline inside `HR_out`'s expression,
again inline inside `BB_pct`'s expression) on top of the separate `$ODE`
computation (`BB_Cp`) — three total recomputations of the same physical
concentration for BB, and two each ($ODE + $TABLE) for the other four
compounds. Consolidated: one `double C_<STEM> = CENT_<STEM>/V1_<STEM>;` in
`$ODE` (the canonical site the PD math reads, satisfying the "exactly one
concentration variable" rule and the corpus's `double C_<STEM> = <expr>;`
discoverability contract), and — since mrgsolve does not share `$ODE`-local
variables with `$TABLE` (confirmed live below) — one more re-derivation in
`$TABLE` for the reporting-side quantities (`HR_out`, `*_pct`, `Cp_<STEM>`
captures), now computed once per compound there too and reused by every
`$TABLE` expression that needs it (`HR_out` and `BB_pct` both now read the
same `C_BB`, instead of each recomputing it independently).

## A live mrgsolve build gotcha hit and fixed during verification

First `/model_manifest` attempt failed to compile:
`redefinition of 'double {anonymous}::C_ACEI'` (and same for `C_ARB`,
`C_CCB`, `C_BB`, `C_HCTZ`, `EFFECT_ACEI`, `EFFECT_ARB`, `EFFECT_BB`,
`EFFECT_CCB`, `EFFECT_HCTZ`). mrgsolve forward-declares every
`double NAME = expr;` it finds anywhere in `$ODE`/`$TABLE` at file scope in
the generated C++, so declaring the same name with a second `double NAME =
...;` in `$TABLE` (after already declaring it in `$ODE`) is a genuine
redefinition, not two independent block-local scopes. Fixed by dropping
`double` from the `$TABLE`-side statements and using bare reassignment
(`C_ACEI = CENT_ACEI / V1_ACEI;`) instead, which reuses the forward
declaration mrgsolve already emitted from the `$ODE` initializer — this is
the same pattern `driver-patches/HANDOFF.md` documents for the
`$GLOBAL`-forward-declare case, except here mrgsolve does the forward
declaration automatically rather than it being written by hand. This is a
syntax-only fix (confirmed non-numeric: the resulting `/run_simulation`
comparison below is an exact match), not a change to any compound's PK/PD.
The single `double C_<STEM> = <expr>;` statement required for downstream
discoverability still exists exactly once, in `$ODE`; grepping the finished
file for `double C_<STEM> = ` and `EC50_<STEM>` both hit for all 5
compounds.

## Verification

Ran via the qspserver `mrgsolve_api` (`http://localhost:8007`), never local
R/mrgsolve, per the guide:

1. `POST /model_manifest` on both the extracted original DSL and the
   extracted refactored DSL — both compile cleanly (200). The original has
   **no pre-existing mrgsolve-2.0.1 build defect**; the only compile issue
   encountered was the `$ODE`/`$TABLE` redefinition above, introduced (and
   fixed) by this refactor itself, not inherited from the original.
2. `POST /run_simulation` for all 6 of the original file's own dosing
   scenarios (No Treatment; ACEI 10 mg QD; ARB 100 mg QD; CCB 10 mg QD; BB
   10 mg QD; Triple therapy ACEI 5 mg + CCB 5 mg + HCTZ 12.5 mg QD),
   identical dosing built from the same numbers the original's own
   `make_dose()` calls use, against both models, full 4032 h / 1 h grid
   (4033-4036 rows depending on scenario dose-time duplication), comparing
   all 23 outputs common to both `$CAPTURE` blocks (`MAP`, `SBP`, `DBP`,
   `PP`, `HR`, `CO`, `TPR`, `PV`, `AngII`, `Aldo`, `NO`, `LVM`, `eGFR`,
   `Cp_ACEI`, `Cp_ARB`, `Cp_CCB`, `Cp_BB`, `Cp_HCTZ`, `ACE_inhib_pct`,
   `AT1R_block_pct`, `BB_block_pct`, `VGCC_block_pct`, `NCC_inhib_pct`).
3. **Result: exact match. `maxabs = 0` and `maxrel = 0` for every one of
   the 23 shared outputs, at every timestep, in all 6 scenarios.** This is
   the expected outcome for a pure archetype-1/2 reorganization with no
   Hill-fitting (per the guide's tolerance rule) — same math, renamed and
   reorganized, requests spaced >=2s apart, one at a time (well under the
   API's 2-concurrent-job budget).
4. `Rscript -e 'parse(...)'` succeeds on the finished
   `eh_mrgsolve_model_refactored.R` (22 top-level expressions), after fixing
   one straight apostrophe in a `$ODE` comment ("original's" ->
   "original’s") that would otherwise have closed the R string early — the
   single-most-common defect class this guide's structural template
   section warns about. No other structural issues found.

## Found in passing: the model's own "No Treatment" baseline is not physiologically stable (logged as `UPSTREAM_ISSUES.md` #151)

`CO_L_0 = MAP0/80.0 = 1.25` L/min is inconsistent with the model's own
cardiac-output dynamics (`CO_target = HR0*SV0/1000 = 4.9` L/min at
baseline) — `CO_L` relaxes from 1.25 to ~4.84 L/min within ~10-15 hours in
every scenario including "No Treatment", dragging `MAP` from the documented
100 mmHg baseline up to a ~377-388 mmHg plateau (with `TPR_N` unchanged at
1.0) before any RAAS/TPR feedback has time to act. **Confirmed upstream,
not introduced by this refactor**: reproduced with zero difference by the
refactored file in the verification above (same values at every one of
4033+ timesteps, all 6 scenarios). Full account in `UPSTREAM_ISSUES.md`
#151. Not fixed here — parameter values and initial conditions are out of
scope for a PK-reorganization refactor — but flagged prominently since it
affects the interpretation of every scenario's absolute MAP/SBP/DBP level
(the relative direction and ranking of the treatment arms against "No
Treatment" is unaffected, since all 6 scenarios share the same `CO_L_0`
inconsistency equally).

## Documentation inconsistency noted, not fixed

The original's own header comment claims "22 ODE states"; the actual
`$CMT` block (both files) declares 18 (9 drug-PK + 2 RAAS + 5
haemodynamic + 2 remodeling — the block's own inline sub-counts also
don't sum consistently, e.g. "RAAS/PD compartments (6 compartments)" over
a block listing only `ANGII`/`ALDO`). Corrected in the refactored file's
own header comment to state the true count (18) and explain the
discrepancy; not changed in the untouched original.

## Census update

`driver-patches/data/compound_perturbation_census.md` rows for
`essential-hypertension` (ACEI, ARB, BB, CCB, HCTZ) updated in place with
confirmed compound identity, target/pathway, and outcome — see that file.
