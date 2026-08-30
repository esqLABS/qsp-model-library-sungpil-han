# Refactor notes — `eoe_mrgsolve_model_refactored.R`

Scope: **mepolizumab (MEPO) and cendakimab (CENDA) only** — their PK blocks
and their downstream IL-5 / IL-13 effect equations. Budesonide, dupilumab,
and every disease equation not driven by mepolizumab/cendakimab
concentration are byte-for-byte identical to `eoe_mrgsolve_model.R` (the
post-DSL R script — `mcode()` call, `build_events()`, scenario list,
plotting, sensitivity analysis — is verified byte-identical to the original
in its entirety; see Verification below).

## Archetype determination

Read the actual `dxdt_` lines rather than assuming from the task brief's
"likely TMDD-flavored, similar to tocilizumab" hint — neither compound uses
receptor-binding kinetics at all.

**Mepolizumab:**
```
dxdt_MEPO_SC = -ka_mepo * MEPO_SC;
dxdt_MEPO_C  = ka_mepo * F_mepo * MEPO_SC / Vd_mepo
               - (CL_mepo / Vd_mepo) * MEPO_C;
```
**Cendakimab:**
```
dxdt_CENDA_GUT = -ka_cenda * CENDA_GUT;
dxdt_CENDA_C   = ka_cenda * F_cenda * CENDA_GUT / Vd_cenda
                 - (CL_cenda / Vd_cenda) * CENDA_C;
```

Both are a depot compartment feeding a single "central" compartment, no
peripheral compartment, no receptor-binding ODEs anywhere (confirmed by
grepping the whole file for `KON`/`KOFF`/`REC_FREE`/`COMPLEX` — none exist
for either compound). But neither central compartment is a plain
amount-in-a-volume: `MEPO_C`/`CENDA_C`'s own `$ODE` line divides the depot
input by `Vd_mepo`/`Vd_cenda` *on the way in*, and eliminates at
`(CL/Vd)*state` — a first-order rate applied directly to the *concentration*
state, not to an amount. So `MEPO_C`/`CENDA_C` already **are** the
concentrations (mg/L), exactly the same "compartment IS concentration"
quirk documented for tocilizumab in `sepsis/sep_refactor_notes.md` and for
the antiviral compound in `dengue/denv_refactor_notes.md`.

**Conclusion: Archetype 3 minus the peripheral compartment (depot + central,
linear, no peripheral), combined with Archetype 1's "compartment IS
concentration" behavior** for both compounds — not Archetype 4. This is not
a case of "none of these fit"; it is the corpus's third recurring variant of
Archetype 1/3 (bare concentration state, sometimes with a depot, sometimes
without), just not written down as its own numbered archetype in the guide.

## Renaming applied

Both compounds' stems (`MEPO`, `CENDA`) already matched the original's own
identifiers, so only the naming *pattern* changed, not the stems:

| Original | Refactored |
|---|---|
| `MEPO_SC` | `GUT_MEPO` |
| `MEPO_C` | `CENT_MEPO` |
| `ka_mepo` | `KA_MEPO` |
| `F_mepo` | `F_MEPO` |
| `CL_mepo` | `CL_MEPO` |
| `Vd_mepo` | `V1_MEPO` |
| `Emax_mepo_IL5` | `EMAX_MEPO` |
| `IC50_mepo_IL5` | `EC50_MEPO` |
| — | `GAMMA_MEPO = 1` (new; original had no explicit Hill exponent) |
| `CENDA_GUT` | `GUT_CENDA` |
| `CENDA_C` | `CENT_CENDA` |
| `ka_cenda` | `KA_CENDA` |
| `F_cenda` | `F_CENDA` |
| `CL_cenda` | `CL_CENDA` |
| `Vd_cenda` | `V1_CENDA` |
| `Emax_cenda_IL13` | `EMAX_CENDA` |
| `IC50_cenda_IL13` | `EC50_CENDA` |
| — | `GAMMA_CENDA = 1` (new) |

`GUT_MEPO`/`CENT_MEPO` and `GUT_CENDA`/`CENT_CENDA` keep the **same
positions** in `$CMT` (6th/7th and 8th/9th declared compartment,
respectively) as `MEPO_SC`/`MEPO_C` and `CENDA_GUT`/`CENDA_C` did, so the
numeric `cmt = 6` / `cmt = 8` dosing already used by this file's own
`build_events()` (mepolizumab and cendakimab scenarios) is unaffected —
confirmed empirically in Verification below, not just asserted.

`C_MEPO = CENT_MEPO` and `C_CENDA = CENT_CENDA`, both **undivided** (not
`C/V1` — `V1_MEPO`/`V1_CENDA` only ever appear inside each compound's own
depot-to-concentration conversion, exactly as `Vd_mepo`/`Vd_cenda` did in
the original), preserved unchanged for exact numerical equivalence.

## Hill interface

Both originals' effect terms are already a plain Emax ratio, so per the
guide this is **a rename, not a refit**:

- `Inh_mepo_IL5 = Emax_mepo_IL5 * MEPO_C / (MEPO_C + IC50_mepo_IL5)` →
  ```
  double C_MEPO = CENT_MEPO;
  double EFFECT_MEPO = EMAX_MEPO * pow(C_MEPO, GAMMA_MEPO) /
                       (pow(EC50_MEPO, GAMMA_MEPO) + pow(C_MEPO, GAMMA_MEPO));
  double Inh_mepo_IL5 = EFFECT_MEPO;
  ```
- `Inh_cenda_IL13 = Emax_cenda_IL13 * CENDA_C / (CENDA_C + IC50_cenda_IL13)` →
  ```
  double C_CENDA = CENT_CENDA;
  double EFFECT_CENDA = EMAX_CENDA * pow(C_CENDA, GAMMA_CENDA) /
                        (pow(EC50_CENDA, GAMMA_CENDA) + pow(C_CENDA, GAMMA_CENDA));
  double Inh_cenda_IL13 = EFFECT_CENDA;
  ```

Both `EMAX_*`/`EC50_*` carry the originals' exact literal values; `GAMMA_*
= 1` makes `pow(x,1) = x`, so both rewrites are algebraically identical to
the originals. `Inh_mepo_IL5` and `Inh_cenda_IL13` are kept as local aliases
so every downstream disease equation that reads them (`IL5_effective`,
`dxdt_IL5`, `IL13_signal`, `EPBAR_ss`, `dxdt_IL13`, `dxdt_FIBRO`) is
**untouched** — only the two effect-definition lines themselves changed.
No refit was performed for either compound; there is no fit-quality
tradeoff to report.

`C_MEPO`, `EFFECT_MEPO`, `C_CENDA`, `EFFECT_CENDA` are state-dependent
(computed from compartment state every `$MAIN` evaluation), so — per the
precedent set in `thyroid-eye-disease/ted_refactor_notes.md` and
`dengue/denv_refactor_notes.md` — they cannot themselves be `$PARAM` values.
`EMAX_MEPO`/`EC50_MEPO`/`GAMMA_MEPO`/`EMAX_CENDA`/`EC50_CENDA`/`GAMMA_CENDA`
**are** declared in `$PARAM` (confirmed discoverable via
`/model_manifest`); `C_MEPO_OUT`/`EFFECT_MEPO_OUT`/`C_CENDA_OUT`/
`EFFECT_CENDA_OUT` were added to `$TABLE`/`$CAPTURE` (additive only, a
capture unique to the refactored file) so the four state-dependent
quantities are at least discoverable as **outputs** via `/model_manifest`
and retrievable via `/run_simulation`'s `outputs` list — confirmed both ways
(see Verification).

## Upstream defect found and logged (not fixed here)

Confirmed via the qspserver `mrgsolve_api` container (`http://localhost:8007`,
healthy) that **the untouched original does not compile at all** under
mrgsolve 2.0.1 — three layered defects, none related to mepolizumab or
cendakimab:

1. `$INIT @annotated` uses `NAME = value : description` syntax throughout
   (all 18 lines, e.g. `BUD_ESO  = 0     : Initial budesonide esoph
   (mg/L)`). mrgsolve 2.0.1's annotated-block parser requires
   colon-delimited `NAME : value : description` (the format the same file's
   own `$PARAM @annotated` and `$CMT @annotated` blocks already use) —
   `Error: improper annotation format`. Reproduced in isolation with a
   minimal 1-compartment model using the identical `NAME = value :
   description` pattern; switching only the separator to `:` fixed the
   parse. Not specific to any one line — every `$INIT` line in the file
   uses this syntax.
2. Once (1) is worked around, `$CMT @annotated` and `$INIT` are revealed to
   jointly redeclare the same 18 compartment names — the same defect class
   as issue #34 (`chronic-hypothyroidism/hypo_mrgsolve_model.R`):
   `invalid class "mrgmod" object: Duplicated model names`.
3. Once (2) is worked around, `$CAPTURE` is revealed to list six
   compartment names directly (`IL13`, `IL5`, `EOTAX3`, `EOS_ESO`,
   `MAST_ESO`, `EPBAR`) — the same defect class as issue #34's second half
   and issue #30 (`sepsis`): `compartment should not be in $CAPTURE`.

**Verification workaround (in-memory only, not committed to either
file):** applied identically to scratch copies of `eoe_mrgsolve_model.R`
and `eoe_mrgsolve_model_refactored.R`: (a) the `$INIT @annotated` block
deleted and its 18 assignments moved into `$MAIN` using the `<CMT>_0 =
value;` idiom (same approach as issue #34), which resolves both (1) and (2)
in one step since it removes the malformed block entirely; (b) the six
compartment names removed from `$CAPTURE` (mrgsolve reports compartment
states regardless of `$CAPTURE`, so nothing is lost). No numeric value or
equation was touched by either change. Neither the tracked original nor the
delivered `_refactored.R` contains this workaround — both still contain the
`$INIT`/`$CMT`+`$INIT`/`$CAPTURE` defects exactly as written, so **neither
currently builds against mrgsolve 2.0.1** without the same workaround.
Logged as a new entry in `translations/UPSTREAM_ISSUES.md`.

Note on the refactored file specifically: renaming `MEPO_SC`/`MEPO_C`/
`CENDA_GUT`/`CENDA_C` required rewriting their `$CMT`/`$INIT` lines anyway.
The renamed `$INIT` lines keep the same malformed `NAME = value :
description` syntax as the rest of the (untouched) `$INIT` block — including
preserving the original's own missing space in `CENDA_GUT= 0` as
`GUT_CENDA= 0` — so the refactor does not incidentally fix a defect outside
its assigned scope, consistent with "log what you find, don't fix it
upstream" extending to the delivered sibling file as well as the original.

## Verification

Ran via the qspserver `mrgsolve_api` container (`POST /model_manifest` then
`POST /run_simulation`, `http://localhost:8007`, confirmed healthy) against
the **workaround-patched scratch copies** described above (identical
patch applied to both models; the tracked original and delivered
`_refactored.R` are unmodified).

`/model_manifest` confirmed both compile and that `KA_MEPO`, `F_MEPO`,
`CL_MEPO`, `V1_MEPO`, `EMAX_MEPO`, `EC50_MEPO`, `GAMMA_MEPO`, `KA_CENDA`,
`F_CENDA`, `CL_CENDA`, `V1_CENDA`, `EMAX_CENDA`, `EC50_CENDA`,
`GAMMA_CENDA` are all discoverable `$PARAM` entries in the refactored
manifest, and that `C_MEPO_OUT`/`EFFECT_MEPO_OUT`/`C_CENDA_OUT`/
`EFFECT_CENDA_OUT` are discoverable output paths.

Ran the original file's **own two named dosing scenarios** for these
compounds (from `build_events()`), reproduced as `dosing` records against
`/run_simulation` (identical dosing, both models):

| Scenario | Dosing | Result |
|---|---|---|
| `mepolizumab` | 300 mg SC q28d into `cmt=6` (`MEPO_SC`/`GUT_MEPO`), `end=364, delta=7` | exact match, 0 deviation |
| `cendakimab` | 160 mg PO QD into `cmt=8` (`CENDA_GUT`/`GUT_CENDA`), `end=364, delta=7` | exact match, 0 deviation |

Full `end=364` duration was used for both (52 weeks, matching `SIM_END` in
the original file) — no shortening was needed; neither scenario approached
mrgsolve's default 20000-step budget.

For each scenario, all **35 shared outputs** (every compartment and
`$CAPTURE` entry in the original, mapped through the compartment renames
`MEPO_SC→GUT_MEPO`, `MEPO_C→CENT_MEPO`, `CENDA_GUT→GUT_CENDA`,
`CENDA_C→CENT_CENDA`) were compared pointwise across the full time grid (67
points, t=0–365, for the mepolizumab scenario; 417 points, t=0–364, for the
cendakimab scenario). **Result: exact match, max absolute and max relative
deviation 0.000e+00 on every output, both scenarios** — including the
non-drug-specific outputs (`IL13`, `EOS_ESO`, `EREFS_SCORE`, etc.), the
other compound's untouched PK (`BUD_ESO`, `DUP_C`, `DUP_P`), and both
compounds' own renamed PK/effect outputs (`MEPO_TROUGH`, `CENDA_CONC`).
This is the expected result for a pure structural rename with no
Hill-fitting (Archetype 3-minus-peripheral/Archetype 1 hybrid, gamma=1
throughout).

Sanity-checked the dynamics are non-trivial (not just both sides
degenerately zero): mepolizumab dosing drives `IL5` from baseline 15 down
to ~4 pg/mL by day 21 and `MEPO_TROUGH` cycles through the expected
30–50 mg/L range between q28d doses; cendakimab dosing drives `IL13` from
baseline 80 down to ~52 pg/mL by day 21, consistent with each compound's
intended pharmacology.

## Census update

Recorded in `driver-patches/data/compound_perturbation_census.md`
(`eosinophilic-esophagitis` rows for Cendakimab and Mepolizumab).

## Anything else flagged

- No other compound's compartments, parameters, or effect equations were
  touched. Budesonide (`BUD_ESO`/`BUD_SYS`) and dupilumab (`DUP_SC`/
  `DUP_C`/`DUP_P`) matched **exactly** (0.000e+00 deviation) in both
  scenarios, as expected since that code is byte-identical between the two
  files.
- The post-DSL R script (compile call, `build_events()`, `SCENARIOS` list,
  plotting, sensitivity analysis) was confirmed **byte-identical** between
  the original and refactored files by direct string comparison — only the
  quoted DSL block itself and the added header comment differ.
- The extracted `.cpp` used for API verification (not a delivered
  artifact — deleted after use) was confirmed byte-identical to the quoted
  DSL string re-extracted from the finished `eoe_mrgsolve_model_refactored.R`,
  so no drift was introduced between the working copy and the delivered
  file.
