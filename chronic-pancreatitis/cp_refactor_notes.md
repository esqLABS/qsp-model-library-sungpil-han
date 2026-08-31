# Refactor notes — `chronic-pancreatitis/cp_mrgsolve_model.R`

Compounds refactored: **OPIOID (Tramadol)** and **PERT (pancreatic enzyme
replacement therapy)** — the two rows in
`driver-patches/data/compound_perturbation_census.md` for
`chronic-pancreatitis`. Every disease-side compartment and equation
(inflammation `TNFa`/`IL6`/`TGFb`/`ROS`, PSC activation, fibrosis, exocrine
function, beta-cell mass, glucose, peripheral/central pain sensitisation)
is unchanged in structure and parameter values from `cp_mrgsolve_model.R`;
only the two compounds' own PK compartments, their exposed concentrations,
and their Hill effect terms were renamed/consolidated to this fork's
convention.

## Drug identity

- **OPIOID = Tramadol**, oral, per the original's own `$PARAM` comment
  ("Opioid PK (Tramadol oral model)") and its calibration note at the
  bottom of the file ("Opioid PK: 3-compartment model based on Tramadol
  data from Grond & Sablotzki (2004) Clin Pharmacokinet 43:879-923").
- **PERT = pancreatic enzyme replacement therapy** (generic lipase-based
  enteric-coated microspheres, e.g. pancrelipase), per the original's own
  `$PARAM` comment ("PERT PK (enteric-coated microspheres)") and dosing
  scenario (40,000 lipase units per main meal, 3x/day) — the model does
  not name a specific brand, so the compound identity carried forward is
  the generic "PERT" already used by the census row and the file's own
  variable names.

## Why OPIOID needed "normalize duplicate concentration sites" and PERT did not

The original computed the opioid plasma concentration **twice**, under two
independent names, that could drift apart under future edits:

- `Cop = CENT_OPIOID / V1_opioid` in `$ODE` (fed both Hill effect terms,
  `opi_pS`/`opi_cS`), and
- `Opioid_Cplasma = CENT_OPIOID / V1_opioid` in `$TABLE` (fed only the
  capture/plot output) — an independent re-derivation of the identical
  formula under a different name, with no code relationship to `Cop`.

Both were consolidated into one canonical `C_OPIOID` definition (assigned,
not merely renamed, at both the `$ODE` and `$TABLE` sites — see the
`$GLOBAL` note below for why mrgsolve needs the formula written twice).

PERT's concentration was already a single clean site (`Cpert = DUO_PERT`
in `$ODE`), used to compute a `pert_eff` local that turned out to be dead
code (see below); the *effect* itself was independently — but only
once — computed inline in `$TABLE`'s `EXO_FUN_PERT` (`PERT_ON * Emax_pert
* DUO_PERT / (EC50_pert + DUO_PERT) * (1.0 - EXO)`), which is why the
census classified PERT as "redirect (clean single site)" rather than
"normalize duplicate sites" — there was one real site to redirect, not two
to merge.

### `pert_eff` in the original was dead code

`double pert_eff = PERT_ON * Emax_pert * Cpert / (EC50_pert + Cpert);` was
declared in `$ODE` immediately before `dxdt_EXO`, but `dxdt_EXO` never
reads `pert_eff` (`dxdt_EXO = -k_exo_loss * FIB * EXO + 0.001 * (exo_0 -
EXO);`), and `$ODE`-local doubles are not visible in `$TABLE` (mrgsolve
block scoping), so `pert_eff` had zero effect on any output — the model's
actual PERT effect on exocrine function was always the independent inline
Hill term inside `$TABLE`'s `EXO_FUN_PERT`. Removed in the refactor (see
`EFFECT_PERT` below, which replaces both the dead `pert_eff` local and the
inline `$TABLE` term with one named, `$CAPTURE`d quantity). This changes
nothing about any output — confirmed by the exact-match verification
below.

## Archetype per compound

- **OPIOID (Tramadol)**: **archetype 3** (depot + central + peripheral,
  linear CL/Q/V). The original already used exactly this structure
  (`GUT_OPIOID`/`CENT_OPIOID`/`PERI_OPIOID`, `ka_opioid`/`CL_opioid`/
  `V1_opioid`/`Q_opioid`/`V2_opioid`/`F_opioid`) and even already used the
  guide's compartment-naming convention verbatim — only the `$PARAM` names
  were renamed to upper-case (`KA_OPIOID`, `CL_OPIOID`, `V1_OPIOID`,
  `Q_OPIOID`, `V2_OPIOID`, `F_OPIOID`), pure renames, same values, same
  `dxdt_` algebra (the guide's grouped `-(CL+Q)/V1*CENT` form is
  algebraically identical to the original's split `-(CL/V1)*CENT -
  (Q/V1)*CENT` and was used verbatim).
- **PERT**: **archetype 3 minus peripheral** (depot + central, linear),
  no bioavailability parameter in the original (none invented). The
  original's `GUT_PERT`/`DUO_PERT` dissolution-then-first-order-loss
  structure maps directly; `DUO_PERT` renamed `CENT_PERT` per the
  `CENT_<STEM>` convention. No volume in the original (`Cpert = DUO_PERT`
  directly, not divided by any V) — `C_PERT = CENT_PERT` preserves this
  exactly, same convention as `KE_LOP`/`C_LOP` in
  `bile-acid-diarrhea/bam_mrgsolve_model_refactored.R` and `C_LARA` in
  `celiac-disease/cd_mrgsolve_model_refactored.R`. The original's
  degradation-rate parameter `kd_pert` is likewise volume-less, so it was
  renamed `KE_PERT` (elimination rate, not a true CL/V clearance),
  following the same `KE_<STEM>` precedent as `bam`'s `KE_LOP`, rather
  than the guide's `CL_<STEM>` (which the guide's own template pairs with
  a volume).

## The Hill interface

Both compounds' original effect terms were already the guide's exact
Emax/EC50 ratio shape — every rename below is a rename, not a refit
(`GAMMA_*` is new only in the sense that the original had no explicit
exponent, i.e. implicitly 1):

- **OPIOID** drives *two* separate, independently named effects (matching
  the mitral-regurgitation `DOB` precedent of one compound with two named
  effect variables): `EFFECT_OPIOID_PS` (peripheral sensitisation, was
  `opi_pS = Imax_opioid_pS * Cop / (IC50_opioid_pS + Cop)`, renamed
  `EMAX_OPIOID_PS`/`EC50_OPIOID_PS`/`GAMMA_OPIOID_PS=1`) and
  `EFFECT_OPIOID_CS` (central sensitisation, was `opi_cS = Imax_opioid_cS
  * Cop / (IC50_opioid_cS + Cop)`, renamed `EMAX_OPIOID_CS`/
  `EC50_OPIOID_CS`/`GAMMA_OPIOID_CS=1`).
- **PERT** drives one effect, `EFFECT_PERT` (was the inline `Emax_pert *
  DUO_PERT / (EC50_pert + DUO_PERT)` term inside `$TABLE`'s
  `EXO_FUN_PERT`, renamed `EMAX_PERT`/`EC50_PERT`/`GAMMA_PERT=1`).
  `EXO_FUN_PERT = EXO + PERT_ON * EFFECT_PERT * (1.0 - EXO)` keeps the
  `PERT_ON` switch exactly where the original had it (gating the whole
  effect term, not the concentration).

## `$GLOBAL` predeclaration (mechanical necessity, not a design choice)

mrgsolve 2.0.1 hoists every `$ODE`- and `$TABLE`-local `double` declaration
into one shared anonymous-namespace scope. `C_OPIOID`/`C_PERT`/
`EFFECT_OPIOID_PS`/`EFFECT_OPIOID_CS`/`EFFECT_PERT` each need to exist in
*both* blocks (`$ODE`, because the disease-side pain/exocrine equations
read them; `$TABLE`, because `$ODE`-local doubles are not visible there,
and `EXO_FUN_PERT`/the `$CAPTURE`d covariates need them) — declaring
`double NAME` a second time in `$TABLE` for a name already declared in
`$ODE` is a compile error (`redefinition of double {anonymous}::NAME`),
confirmed directly against the qspserver API on a first attempt at this
refactor. Fix: predeclare the five names once in `$GLOBAL`
(`double C_OPIOID, C_PERT, EFFECT_OPIOID_PS, EFFECT_OPIOID_CS,
EFFECT_PERT;`) and *assign* (not redeclare) them at both sites, using the
identical formula in each — same pattern as `C_MIT` in
`hereditary-spherocytosis/hsph_mrgsolve_model_refactored.R`. This is a
structural/syntax necessity of the split, not a second "duplicate site" —
both assignments are the same canonical formula, kept in lockstep, as
disclosed in-line in the DSL.

## Pre-existing upstream defects found (not fixed upstream, logged)

The original **does not compile at all** under mrgsolve 2.0.1, and
separately has a runtime solver-stability defect unrelated to either
refactored compound. Both logged in `translations/UPSTREAM_ISSUES.md`
(re-checked the file's tail immediately before appending; two other
concurrent refactors had already claimed #108-#110 for unrelated files):

- **`#111`** — `$CAPTURE` re-lists nine `$INIT`-declared compartment names
  directly (`GLUC`, `TNFa`, `IL6`, `TGFb`, `ROS`, `PSC`, `FIB`, `EXO`,
  `BETA`), which mrgsolve 2.0.1 rejects outright (confirmed via
  `POST /model_manifest`: `compartment should not be in $CAPTURE:
  GLUC,TNFa,IL6,TGFb,ROS,PSC,FIB,EXO,BETA`). **Fix applied directly in
  `cp_mrgsolve_model_refactored.R`** (syntax-only, non-numeric, per the
  guide's settled policy): those nine names deleted from `$CAPTURE`. They
  remain present in every simulation output regardless (mrgsolve emits
  every compartment state automatically), so nothing is lost — confirmed
  by the exact-match verification below, where all nine appear
  byte-identical between original and refactored at every shared
  timepoint.
- **`#112`** — `PSC`'s hard ceiling clamp (`if (PSC > PSC_max) PSC_bounded
  = -kr_PSC * PSC;`) is a state-triggered discontinuity sitting directly
  inside `dxdt_PSC`. With the file's own `SEVERITY = 1.5` (used by every
  one of its five R-defined scenarios), `PSC` reaches the `PSC_max = 1.0`
  boundary at roughly hour 130-140 in every scenario tried (confirmed for
  "1. No Treatment", "3. Opioid + Gabapentin", and "4. Antifibrotic" —
  antifibrotic slows but does not prevent the rise), after which the
  derivative discontinuity forces the solver to keep halving its internal
  step, and the qspserver API's `run_simulation` fails with `[lsoda]
  20000 steps taken before reaching tout`. Binary-searched: 120h
  consistently succeeds, 140-150h consistently fails. **Not fixed** (pure
  disease-side defect, out of scope for an OPIOID/PERT PK refactor, and
  the never-edit-upstream rule applies regardless) — reproduced unchanged
  in `cp_mrgsolve_model_refactored.R`. This means the file's own R wrapper
  (`mrgsim(end = sim_end, delta = 24)`, `sim_end` = 2 years = 17,520h, no
  `maxsteps` override anywhere in the file) would itself be unable to
  complete any of its five advertised 2-year scenarios under a solver that
  does not silently raise `maxsteps` far past mrgsolve's low built-in
  default.

## Verification

Per the guide's known API limitation (`run_simulation`'s default
`maxsteps` well below what some scenarios need) *and* the runtime defect
above (which fails regardless of `maxsteps`, since the failure coincides
with a specific model state, not a fixed step budget), all five of the
original's own R-defined scenarios were run over a shortened 0-120h window
(hourly output, `delta = 1`) — the largest window confirmed clear of the
`PSC`-clamp instability in every scenario — rather than each scenario's
full 2-year horizon. Both DSLs (original with only the `#111` syntax fix
applied, and the delivered refactored file) were run through the
qspserver `mrgsolve_api` container (`POST /run_simulation`), ~2.5s apart
per call:

| Scenario (original's own params/dosing) | max abs diff | max rel diff |
|---|---|---|
| 1. No Treatment (`SEVERITY=1.5`, all switches off) | 0.0 | 0.0 |
| 2. PERT Only (`PERT_ON=1`, 40,000 units q8h into `GUT_PERT`) | 0.0 | 0.0 |
| 3. Opioid + Gabapentin (`E_gabapentin=0.35`, 50mg q8h into `GUT_OPIOID`) | 0.0 | 0.0 |
| 4. Antifibrotic (`PIRF_ON=1`, `LOSARTAN_ON=1`) | 0.0 | 0.0 |
| 5. Combination (all switches on, both compounds dosed q8h) | 0.0 | 0.0 |

Exact match (max abs/rel diff 0.0) on every shared output —
`PAIN_SCORE`, `FAT_MALAB`, `EXO_FUN_PERT`, `HBA1C`, `FIB_INDEX`,
`PSC_ACT`, `GLUC`, `TNFa`, `IL6`, `TGFb`, `ROS`, `PSC`, `FIB`, `EXO`,
`BETA`, `pSENS`, `cSENS`, all four PK compartments
(`GUT_OPIOID`/`CENT_OPIOID`/`PERI_OPIOID`/`GUT_PERT`), and the renamed
pair `Opioid_Cplasma`\<->`C_OPIOID` / `DUO_PERT`\<->`CENT_PERT` — at every
one of the 121-123 shared timepoints per scenario. This is exactly the
"near-exact match" the guide expects for a pure structural
reorganization (archetype 3 / archetype 3 minus peripheral, no
Hill-fitting on either compound).

`POST /model_manifest` against the extracted refactored DSL confirms
every covariate this guide's naming convention promises is discoverable
as an output: `outputPaths` includes `C_OPIOID`, `EFFECT_OPIOID_PS`,
`EFFECT_OPIOID_CS`, `C_PERT`, `EFFECT_PERT` (alongside every renamed PK
compartment and every disease-side state/capture), and `parameters`
includes every renamed `$PARAM` name (`KA_OPIOID`, `CL_OPIOID`,
`V1_OPIOID`, `Q_OPIOID`, `V2_OPIOID`, `F_OPIOID`, `KA_PERT`, `KE_PERT`,
`EMAX_PERT`, `EC50_PERT`, `GAMMA_PERT`, `EMAX_OPIOID_PS`,
`EC50_OPIOID_PS`, `GAMMA_OPIOID_PS`, `EMAX_OPIOID_CS`, `EC50_OPIOID_CS`,
`GAMMA_OPIOID_CS`) with its original value as default.

No R wrapper differences beyond the DSL and the necessary downstream
column rename: `p_opioid`'s ggplot call was updated from `Opioid_Cplasma`
to `C_OPIOID` (the new capture name) — cosmetic, not a numeric change.
