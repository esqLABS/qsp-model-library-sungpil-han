# Refactor notes — `nafld-masld/nafld_mrgsolve_model.R`

## Scope: two census rows checked, one corrected as a classifier artifact, one confirmed real and refactored; one more real compound found missing from the census and refactored too

This file carries two rows in `driver-patches/data/compound_perturbation_census.md`:
**"ER Stress (DEG)"** and **"GLP"**. Per the task's specific flag, both were
checked against the actual code rather than trusted at face value, and the
file was also checked for other compounds not carrying a census row — one
was found: the **FXR agonist (OCA-like)**, which is the file's primary,
most heavily-parameterized compound and the source of the file's own
"OCA 25mg/day" scenario.

| Row | Real drug? | Evidence |
|---|---|---|
| "ER Stress (DEG)" | **No — not a compound at all.** No compartment, parameter, or identifier named `DEG` exists anywhere in the file (confirmed by exact-name grep). "DEG" is a substring fragment shared by twelve unrelated first-order degradation-rate parameters scattered across the whole file (`kLTAG_deg`, `kROS_deg`, `kNRF2_deg`, `kER_deg`, `kTNF_deg`, `kIL6_deg`, `kIL1_deg`, `kTGFb_deg`, `kCol_deg`, `kFGF19_deg`, `kLPS_deg`, `kMCP1_deg`); the closest textual match to "ER Stress" is `kER_deg = 0.12 // ER stress resolution rate (h^-1)`, an endogenous clearance-rate constant, not a dosed compound. Same defect shape as the `EL` false positive in `drug-induced-liver-injury/dili_refactor_notes.md` (a shared-suffix fragment across several unrelated rate constants, mis-extracted as if it were a compound stem) and the precedent set for census mislabeling in `neonatal-hyperbilirubinemia/nhb_refactor_notes.md`. |
| "GLP" | **Yes.** A genuine GLP-1 receptor agonist (semaglutide-like), with its own two-compartment PK (`GUT_GLP1` depot, `CENT_GLP1` central; `ka_glp1`, `CL_glp1`, `V1_glp1`), its own three named Hill-shaped effects (`E_GLP1_IR`, `E_GLP1_DNL`, `E_GLP1_Kup`), and its own R-side dosing scenario ("Semaglutide 2.4mg/wk", weekly SC into `cmt=4`). The code's own stem is already `GLP1` (not `GLP`) — used consistently in `ka_glp1`, `CL_glp1`, `V1_glp1`, `GUT_GLP1`, `CENT_GLP1`, `Cp_GLP1`, `EC50_GLP1`, `Emax_GLP1_*`; the census's "Compound" column truncated it. Census's own perturbation-method classification, "Normalize duplicate concentration sites, then redirect," is also confirmed accurate: `Cp_GLP1` is computed once in `$MAIN` (`double Cp_GLP1 = CENT_GLP1 / V1_glp1;`) and a second time, independently, in `$TABLE` (`capture Cp_GLP1_ug = CENT_GLP1 / V1_glp1;`) — a genuine duplicate concentration site (see "Duplicate concentration sites" section below for how this was actually resolved). |
| FXR agonist (OCA-like) | **Yes — missing from the census entirely.** The file's *other* drug, and its most prominent: own three-compartment PK (`GUT`/`CENT`/`PERI`, depot+central+peripheral, `ka`/`CL`/`V1`/`V2`/`Q`/`F1`), own four named Hill-shaped effects (`E_FXR_DNL`, `E_FXR_TGF`, `E_FXR_LPS`, `E_FXR_FGF19`), and the file's own headline scenario ("OCA 25mg/day", `scen2_oca`, `cmt=1`). Confirmed via the file's header comment ("Drug PK Parameters (representative FXR agonist / OCA)"), the `EC50_FXR`/`Emax_FXR_*` parameter names used throughout `$MAIN`, and the REGENERATE-trial calibration note in the file's own top-of-file comment block. Same pattern as APAP in `drug-induced-liver-injury/dili_refactor_notes.md` (the file's own central compound, absent from its census rows). This compound also has the same duplicate-concentration-site pattern as GLP-1 RA (`Cp_FXR` in `$MAIN`, `Cp_FXR_ug` recomputed in `$TABLE`). |

**Conclusion:** "ER Stress (DEG)" is corrected in the census to "Not a
compound — classifier artifact, corrected" (left completely untouched in
the code, nothing to refactor). "GLP" is corrected to "GLP-1 RA
(semaglutide-like, GLP1)" and refactored. A new FXR-agonist/OCA row is
added to the census and refactored.

## What was *not* touched: Resmetirom (THRβ) and Pioglitazone (PPARγ) — dead `$PARAM` entries, no PK block

The original's `$PARAM` block also declares `EC50_THRb`, `Emax_THRb_DNL`,
`Emax_THRb_Box`, `EC50_Pio`, `Emax_Pio_IR`. An exact-name grep across the
whole file confirms none of the five appears anywhere in `$MAIN`, `$ODE`,
or `$TABLE` — no PK compartment, no effect expression, no dosing route.
The R-side "Resmetirom 80mg/day" scenario does not use them either: it
reuses the FXR/OCA depot compartment (`scen2_oca`) with *`Emax_FXR_DNL`/
`Emax_FXR_TGF`/`EC50_FXR` overridden*, never touching `Emax_THRb_*`.
These five parameters are carried over into the refactored file
byte-for-byte unchanged (still under their original names, not renamed to
`EMAX_THRB_*`/`EC50_THRB` convention) — there is no PK block to reorganize
and nothing to rename a compound stem *for*. This is disclosed here rather
than silently dropped or silently renamed.

## Archetype per compound

### FXR agonist (OCA-like), stem `FXR` — Archetype 3 (depot + central + peripheral, linear)

Clean match: `GUT`→`GUT_FXR` (depot), `CENT`→`CENT_FXR` (central),
`PERI`→`PERI_FXR` (peripheral), first-order absorption with
bioavailability, first-order + inter-compartmental clearance — exactly the
guide's Archetype 3 template shape, values unchanged.

| Original | Refactored | Value | Role |
|---|---|---|---|
| `GUT` (cmt) | `GUT_FXR` | -- | depot |
| `CENT` (cmt) | `CENT_FXR` | -- | central |
| `PERI` (cmt) | `PERI_FXR` | -- | peripheral |
| `ka` | `KA_FXR` | 0.8 (1/h) | absorption rate |
| `F1` | `F_FXR` | 0.65 | oral bioavailability |
| `CL` | `CL_FXR` | 4.5 (L/h) | clearance |
| `V1` | `V1_FXR` | 15.0 (L) | central volume |
| `V2` | `V2_FXR` | 30.0 (L) | peripheral volume |
| `Q` | `Q_FXR` | 2.0 (L/h) | inter-compartmental clearance |
| `Cp_FXR` (local) | `C_FXR` | -- | **the exposed concentration** |
| `EC50_FXR` | `EC50_FXR` (unchanged) | 0.5 (µg/mL) | Hill EC50, shared across all 4 FXR effects, same as the original |
| `Emax_FXR_DNL` | `EMAX_FXR_DNL` | 0.40 | Hill Emax, DNL pathway |
| `Emax_FXR_TGF` | `EMAX_FXR_TGF` | 0.35 | Hill Emax, TGF-β pathway |
| `Emax_FXR_LPS` | `EMAX_FXR_LPS` | 0.25 | Hill Emax, LPS pathway |
| `1.5` (hardcoded literal) | `EMAX_FXR_FGF19` (new name) | 1.5 | Hill Emax, FGF19 induction pathway — same value, now named |
| -- (none, implicit) | `GAMMA_FXR_DNL`/`_TGF`/`_LPS`/`_FGF19` (new, 4×) | 1.0 each | Hill exponent [original had no explicit Hill coefficient for any of the 4] |
| `E_FXR_DNL` | `EFFECT_FXR_DNL` | -- | fractional DNL reduction |
| `E_FXR_TGF` | `EFFECT_FXR_TGF` | -- | fractional TGF-β reduction |
| `E_FXR_LPS` | `EFFECT_FXR_LPS` | -- | fractional LPS reduction |
| `E_FXR_FGF19` | `EFFECT_FXR_FGF19` | -- | FGF19 induction factor |

### GLP-1 RA (semaglutide-like), stem `GLP1` — Archetype 3 minus peripheral (depot + central, linear)

Clean match to the guide's explicit "drop `GUT_TCZ`/`KA_TCZ`/`F_TCZ` for a
2-compartment-no-depot variant" allowance, in reverse: depot + central,
*no* peripheral compartment, no explicit bioavailability factor in the
original (subcutaneous dosing; no `F1_glp1` existed to rename, so none was
invented). `GUT_GLP1`/`CENT_GLP1` already matched the `GUT_<STEM>`/
`CENT_<STEM>` convention in the original — no compartment rename needed,
only the parameter names.

| Original | Refactored | Value | Role |
|---|---|---|---|
| `GUT_GLP1` (cmt, unchanged) | `GUT_GLP1` | -- | depot (SC) |
| `CENT_GLP1` (cmt, unchanged) | `CENT_GLP1` | -- | central |
| `ka_glp1` | `KA_GLP1` | 0.005 (1/h) | absorption rate |
| `CL_glp1` | `CL_GLP1` | 0.055 (L/h) | clearance |
| `V1_glp1` | `V1_GLP1` | 8.5 (L) | central volume |
| `Cp_GLP1` (local) | `C_GLP1` | -- | **the exposed concentration** |
| `EC50_GLP1` (unchanged) | `EC50_GLP1` | 0.3 (µg/mL) | Hill EC50, shared across all 3 GLP1 effects, same as the original |
| `Emax_GLP1_IR` | `EMAX_GLP1_IR` | 0.45 | Hill Emax, IR pathway |
| `Emax_GLP1_DNL` | `EMAX_GLP1_DNL` | 0.30 | Hill Emax, DNL pathway |
| `Emax_GLP1_Kup` | `EMAX_GLP1_Kup` | 0.35 | Hill Emax, Kupffer pathway |
| -- (none, implicit) | `GAMMA_GLP1_IR`/`_DNL`/`_Kup` (new, 3×) | 1.0 each | Hill exponent [original had no explicit Hill coefficient] |
| `E_GLP1_IR` | `EFFECT_GLP1_IR` | -- | fractional IR reduction |
| `E_GLP1_DNL` | `EFFECT_GLP1_DNL` | -- | fractional DNL reduction |
| `E_GLP1_Kup` | `EFFECT_GLP1_Kup` | -- | fractional Kupffer-activation reduction |

All parameter *values* are copied verbatim from the original.

## Hill interface: rename, not a fit — for both compounds

Every one of the seven effect terms (4 FXR + 3 GLP-1) was already exactly
the plain Hill/Emax ratio shape in the original (`Emax*C/(EC50+C)`, no
explicit exponent), so this is a rename to `EMAX_<STEM>_<pathway> *
pow(C_<STEM>, GAMMA_<STEM>_<pathway>) / (pow(EC50_<STEM>,
GAMMA_<STEM>_<pathway>) + pow(C_<STEM>, GAMMA_<STEM>_<pathway>))` with
`GAMMA = 1` throughout, not an `nls()` fit. The one deviation from a
literal 1:1 rename is `E_FXR_FGF19 = 1.5 * Cp_FXR / (EC50_FXR + Cp_FXR)`,
whose Emax (`1.5`, an induction factor greater than 1, reflecting FXR's
role inducing rather than suppressing FGF19) was a bare numeric literal in
the original with no name at all — promoted to `EMAX_FXR_FGF19 = 1.5`,
same value, same shape, now named per the guide's convention (same
pattern as the NF-κB/ROS magic-number promotions in
`abdominal-aortic-aneurysm/aaa_refactor_notes.md`'s statin row).

## Duplicate concentration sites — how they were actually resolved

The census correctly flagged both compounds' concentration as computed
twice: once in `$MAIN` (`Cp_FXR`, `Cp_GLP1`, feeding the disease-facing
effect terms) and once independently in `$TABLE` (`Cp_FXR_ug`,
`Cp_GLP1_ug`, feeding only the report). The literal reading of "normalize
duplicate sites, then redirect" would be to collapse these into one
site — but mrgsolve evaluates `$MAIN` once per reporting/dosing interval
using *start-of-interval* state, while `$TABLE` runs at the actual
requested output time; pointing `$TABLE`'s report values at `$MAIN`'s
`C_FXR`/`C_GLP1` would inject that timing lag into what was previously a
fresh, unlagged report value — exactly the same finding already made for
`idiopathic-pulmonary-fibrosis/ipf_refactor_notes.md`'s `Cp_pirf`/
`Cn_nint` (kept deliberately separate from `$MAIN`'s `C_PIR`/`C_NIN` for
the identical reason). Following that precedent: `$TABLE`'s
`Cp_FXR_ug`/`Cp_GLP1_ug` are kept as their own independent recomputation
(same formula, only the renamed compartments/params substituted), *not*
collapsed into a reference to `$MAIN`'s `C_FXR`/`C_GLP1`. What *was*
normalized: the `$TABLE` block's original `capture NAME = EXPR;` shorthand
syntax was converted to plain `double NAME = EXPR;` declarations plus one
consolidated `$CAPTURE` line at the end (needed anyway, since `$MAIN`-
declared `C_FXR`/`C_GLP1`/`EFFECT_FXR_*`/`EFFECT_GLP1_*` must also be
captured, and mixing the `capture NAME = EXPR;` inline shorthand with a
separate bare-name `$CAPTURE` list risks the double-registration collision
`neonatal-hyperbilirubinemia/nhb_refactor_notes.md` documented — the fix
used here instead is the working pattern already verified in
`idiopathic-pulmonary-fibrosis/ipf_mrgsolve_model_refactored.R`: `double`
declarations, wherever the original computed them, captured once via a
single bare-name `$CAPTURE` list at the very end). This is a syntax-only
reorganization of `$TABLE` — every captured value is the exact same
formula as the original computed, confirmed unchanged by the verification
below (`Cp_FXR_ug`/`Cp_GLP1_ug` match the original's own
`Cp_FXR_ug`/`Cp_GLP1_ug`, max abs diff 0.0, at every timepoint including
every dose-coincident reporting instant).

## When the original doesn't compile at all

Two pre-existing, unrelated build defects, both confirmed via
`POST /model_manifest` against the untouched original, logged as new
`translations/UPSTREAM_ISSUES.md` entries **#119** and **#120** (tail
re-checked immediately before appending — highest prior entry was #118).

1. **`$CMT` + `$INIT` jointly redeclare all 20 compartments**
   (`Duplicated model names: GUT CENT PERI GUT_GLP1 CENT_GLP1 ...`) — the
   same defect shape as `idiopathic-pulmonary-fibrosis` (#114) and
   `lymphangioleiomyomatosis` (#117). Fixed in the delivered
   `_refactored.R` by deleting `$INIT` and moving its 20 assignments into
   `$MAIN` via the `<CMT>_0 = value;` idiom (e.g. `GUT_FXR_0 = 0;`) —
   declares no new compartment, changes no numeric value, same order,
   same 1-based dosing indices.
2. **`dxdt_MCP1` uses an invalid `where` clause** — `mrgsolve` has no
   `where` keyword; the original wrote
   `dxdt_MCP1 = 0.15 * NFKB_drive - kMCP1_deg * MCP1 where NFKB_drive =
   KUP_ACT * (1 + TNF * 0.2);`, using `NFKB_drive` before it is declared.
   Fixed by declaring `double NFKB_drive = KUP_ACT * (1 + TNF * 0.2);` as
   its own line immediately before the `dxdt_MCP1` assignment and
   dropping the `where` clause — same formula, same value.

Both fixes are syntax-only, non-numeric, non-behavioral, applied directly
in `nafld_mrgsolve_model_refactored.R` (never in the checked-in original)
and disclosed here per the guide's settled policy. Confirmed both defects
are pre-existing and unrelated to either compound's own PK: reproduces
identically from the untouched original alone, via
`POST /model_manifest` against the qspserver `mrgsolve_api` container.

## Verification (qspserver `mrgsolve_api`, `http://localhost:8007`)

**Method.** Both files' embedded DSL blocks (`nafld_model_code`/
`nafld_refactored_code`) were mechanically extracted (regex on the
assignment, verbatim quoted text) and POSTed to the local qspserver
`mrgsolve_api` service (`/model_manifest` then `/run_simulation`), which
compiles and runs each DSL block directly with mrgsolve 2.0.1
server-side — no local R/mrgsolve install used. The **original** DSL, as
sent, was patched only with the two syntax-only fixes above (in-memory
only, matching the fixed content of the delivered `_refactored.R`, never
saved back over the checked-in original) so that a like-for-like build
could even be attempted — per the guide, "the delivered `_refactored.R`
should still actually compile," and the original's own untouched build
failure was independently confirmed first (see #119/#120 above) before
any patching. Requests were spaced several seconds apart, run
sequentially (never more than one in flight), respecting the service's
`max_concurrent_jobs: 2` limit. `POST /run_simulation`'s `dosing` field
addresses compartments by 1-based index; compartment order is identical
between the two files (only names changed) — confirmed from
`/model_manifest`'s `outputPaths`: `GUT`/`GUT_FXR`=1, `CENT`/`CENT_FXR`=2,
`PERI`/`PERI_FXR`=3, `GUT_GLP1`=4, `CENT_GLP1`=5 in both.

**Scenarios run — the file's own, not invented.** All three use the full
2-year duration the original's own `run_scenario()` wrapper uses
(`end=17520, delta=24`, 732–733 points) — no shortening was needed; none
of the three approached the API's default `maxsteps` budget:

1. **`scen2_oca`, "OCA 25mg/day"** (25 mg into `GUT`/`GUT_FXR`, cmt 1,
   q24h × 730 doses): **exact match, max abs diff 0.0** across all 36
   shared outputs (all 21 `$CMT` compartments plus all 15 pre-existing
   `$TABLE`/`$CAPTURE` names), at every one of 732 timepoints — including
   every reporting instant that coincides exactly with a daily dose (the
   guide's "dose-instant reporting artifact" scenario, since `ii=24`
   equals the reporting `delta=24`); no artifact of any size was
   observed.
2. **`scen3_sema`, "Semaglutide 2.4mg/wk"** (2.4 mg into `GUT_GLP1`, cmt
   4, q168h × 104 doses): **exact match, max abs diff 0.0** across all 36
   shared outputs, 732 timepoints.
3. **`scen4_combo`, "OCA + Semaglutide"** (both doses simultaneously):
   **exact match, max abs diff 0.0** across all 36 shared outputs, 733
   timepoints.

The "Resmetirom 80mg/day" and "No Treatment" scenarios were not run as
separate verification cases: per the original's own R code, "Resmetirom"
reuses `scen2_oca`'s exact dosing with parameter overrides (already
exercised structurally by scenario 1 above, with the same PK block), and
"No Treatment" is a zero-amount dummy event (trivial, exercises no drug
PK).

**Result: bit-exact across all three scenarios, consistent with the
guide's tolerance table for Archetype 3 ("pure structural
reorganization... expect a near-exact match").** No Hill fit was needed
(both compounds' effects were already the plain ratio shape); no
deviation of any kind was found, at any timepoint, in any of the 36
shared outputs, across either compound's own dosing regimen or the
combination.

## qspserver `/model_manifest` discoverability

Confirmed via `POST /model_manifest` against the extracted, renamed DSL:
`KA_FXR`, `F_FXR`, `CL_FXR`, `V1_FXR`, `V2_FXR`, `Q_FXR`, `EC50_FXR`,
`EMAX_FXR_DNL`, `GAMMA_FXR_DNL`, `EMAX_FXR_TGF`, `GAMMA_FXR_TGF`,
`EMAX_FXR_LPS`, `GAMMA_FXR_LPS`, `EMAX_FXR_FGF19`, `GAMMA_FXR_FGF19`,
`KA_GLP1`, `CL_GLP1`, `V1_GLP1`, `EC50_GLP1`, `EMAX_GLP1_IR`,
`GAMMA_GLP1_IR`, `EMAX_GLP1_DNL`, `GAMMA_GLP1_DNL`, `EMAX_GLP1_Kup`,
`GAMMA_GLP1_Kup` all appear in the manifest's `parameters` with their
original numeric defaults. `C_FXR`, `C_GLP1`, `EFFECT_FXR_DNL`,
`EFFECT_FXR_TGF`, `EFFECT_FXR_LPS`, `EFFECT_FXR_FGF19`, `EFFECT_GLP1_IR`,
`EFFECT_GLP1_DNL`, `EFFECT_GLP1_Kup` are state-derived (computed in
`$MAIN` from the compartment concentrations) — per the same reasoning
established for `C_COL`/`EFFECT_ANA`/etc. in
`familial-mediterranean-fever/fmf_refactor_notes.md` and `C_SNMP`/etc. in
`neonatal-hyperbilirubinemia/nhb_refactor_notes.md`, they cannot also be
`$PARAM` entries; all nine appear in the manifest's `outputPaths` via the
extended `$CAPTURE` list, confirmed discoverable. `GUT_FXR`, `CENT_FXR`,
`PERI_FXR`, `GUT_GLP1`, `CENT_GLP1` also appear in `outputPaths` as
ordinary compartments (indices 1–5), unchanged order from the original.

No `.cpp` extraction file was left behind — extraction was in-memory only,
used to build the verification requests above and then discarded, per the
workflow guide. All scratch files created during this pass were deleted
before delivery.

## Anything else flagged

- The R-side scenario list (`scen1_nodrug` … `scen5_resmet`) was updated
  only to reference the renamed depot compartments — switched from bare
  numeric `cmt=1`/`cmt=4` to name-based `cmt="GUT_FXR"`/`cmt="GUT_GLP1"`
  (same 1-based positions as the numeric original, but self-documenting
  after the rename; mirrors the style already used in
  `idiopathic-pulmonary-fibrosis/ipf_mrgsolve_model_refactored.R`) — same
  dosing amounts, same timing throughout. `params_override` lists in
  `run_scenario()` calls were updated to the new `EMAX_FXR_*`/`EC50_FXR`/
  `EMAX_GLP1_*` parameter names (values unchanged). No other R-side code
  (plotting, summary tables, responder analysis, sensitivity tornado)
  references any renamed compartment or parameter by name — all read from
  `$TABLE`/`$CAPTURE` output-column names, which are unchanged from the
  original (`Hepatic_TG`, `Kupffer_activation`, `Apoptosis`, `Collagen`,
  `TGFbeta1`, `HSC_act`, `ROS_level`, `ALT`, `AST`, etc.).
- Resmetirom and Pioglitazone `$PARAM` entries are unchanged/unrenamed —
  see the dedicated section above.

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`: the
`nafld-masld` "ER Stress (DEG)" row corrected to flag it as a classifier
artifact (not a compound); the "GLP" row corrected to "GLP-1 RA
(semaglutide-like, GLP1)" and refactored; a new row added for the FXR
agonist (OCA-like) and refactored.
