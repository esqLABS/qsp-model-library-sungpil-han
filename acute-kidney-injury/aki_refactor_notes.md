# AKI (`acute-kidney-injury/aki_mrgsolve_model.R`) — PK/PD refactor notes

Refactored: `aki_mrgsolve_model_refactored.R`. Original untouched.
Three compounds in scope, per the existing census rows
(`driver-patches/data/compound_perturbation_census.md`): **FUR**
(furosemide), **NAC** (N-acetylcysteine), **NE** (norepinephrine). All
three are real dosed compounds in this file (`FUR_GUT`/`FUR_CENT`/
`FUR_PERI`, `NE_CENT`, `NAC_CENT`) with their own `ev()` scenarios in the
original's R wrapper. No other compounds exist in this model.

## Census correction

The census's pre-existing classification for FUR in this file was
"Delete PK compartment; concentration is itself the state." That does
not match this file's actual code: FUR has a genuine depot + 2-compartment
linear PK system (`FUR_GUT` → `FUR_CENT` ⇄ `FUR_PERI`, `KA_FUR`/`F_fur`
absorption, `CL_fur`/`Q_fur`/`Vd_fur1`/`Vd_fur2` disposition), dosed via a
real `ev(cmt="FUR_GUT", ...)` oral/repeat-dose event in the original's own
Scenario 2 and Scenario 6. Per the guide's "match a compound to the
archetype closest to what the original already does" instruction, FUR is
classified here as **Archetype 3** (depot + central + peripheral, linear),
not the "delete PK compartment" pattern. The census row is updated
accordingly rather than forced to match the stale label.

## Archetype per compound

| Compound | Archetype | PK compartments (renamed) |
|---|---|---|
| FUR | 3 — depot + central + peripheral, linear | `GUT_FUR` → `CENT_FUR` ⇄ `PERI_FUR` |
| NE  | 1 — no depot, single compartment, linear | `CENT_NE` |
| NAC | 1 — no depot, single compartment, linear | `CENT_NAC` |

All three are pure structural reorganizations of the original's own
linear-PK math — no Hill-fitting was needed anywhere (see below), so the
guide's "near-exact match" tolerance applies to all three, and that is
what verification found (max abs diff 0.0 everywhere checked).

Compartment order in `$CMT` is preserved from the original (renamed in
place), so the original's own 1-based `ev(cmt=...)` dosing positions are
unchanged: `GUT_FUR`=1, `CENT_FUR`=2, `PERI_FUR`=3, `CENT_NE`=4,
`CENT_NAC`=5, followed by the file's 15 disease compartments in the
original's own order.

## Naming map

| Original | Refactored |
|---|---|
| `FUR_GUT` / `FUR_CENT` / `FUR_PERI` | `GUT_FUR` / `CENT_FUR` / `PERI_FUR` |
| `Vd_fur1` / `Vd_fur2` / `CL_fur` / `Q_fur` / `F_fur` | `V1_FUR` / `V2_FUR` / `CL_FUR` / `Q_FUR` / `F_FUR` |
| `NKCC2_IC50` | `EC50_FUR` |
| `E_NKCC2` | `EFFECT_FUR` |
| `NE_CENT` | `CENT_NE` |
| `CL_ne` / `Vd_ne` | `CL_NE` / `V1_NE` |
| `Emax_map` / `EC50_ne` | `EMAX_NE` / `EC50_NE` |
| `E_NE` | `EFFECT_NE` |
| `NAC_CENT` | `CENT_NAC` |
| `CL_nac` / `Vd_nac` | `CL_NAC` / `V1_NAC` |
| `Emax_gsh` / `EC50_nac` / `Emax_no` | `EMAX_NAC_GSH` / `EC50_NAC` / `EMAX_NAC_NO` |
| `E_NAC_gsh` / `E_NAC_no` | `EFFECT_NAC_GSH` / `EFFECT_NAC_NO` |

`OAT_eff`, `diuresis`, `RBF`, `E_RBF_gfr` are left unrenamed — they are
downstream disease-flux/physiology terms, not PK/Hill roles the naming
convention covers (same treatment as `FUR_DIU`/`FUR_NA` in
`hypercalcemia-of-malignancy/mah_refactor_notes.md`).

New, math-implied parameters (not in the original, added per the guide's
"`GAMMA_<STEM> = 1` if the original had no explicit Hill coefficient"
instruction): `GAMMA_FUR=1.0`, `GAMMA_NE=1.0`, `GAMMA_NAC=1.0`,
`EMAX_FUR=1.0` (the original's NKCC2-block ratio has no separate Emax —
its own maximum value is 1, made explicit here). `EC50_NAC` is
deliberately kept as one shared parameter for both `EFFECT_NAC_GSH` and
`EFFECT_NAC_NO`, exactly as the original's single `EC50_nac` was used for
both `E_NAC_gsh` and `E_NAC_no` — not split into two EC50s the original
never had.

## Finding, not fixed: all three compounds feed their Hill terms from raw
compartment content, not a true concentration

The original never divides any of `FUR_CENT`, `NE_CENT`, or `NAC_CENT` by
its own `Vd` before comparing it against `NKCC2_IC50` / `EC50_ne` /
`EC50_nac` in `$MAIN`, even though those parameters are documented in
concentration units (mg/L, mcg/mL, mg/L respectively). This is the same
pattern in all three compounds in this one file — consistent enough to
read as the author's actual (if dimensionally loose) modeling choice
rather than three independent typos, so it is preserved exactly rather
than "fixed":

```
double FUR_TUB = FUR_CENT * OAT_eff * TCV;             // raw amount, not FUR_CENT/Vd_fur1
double E_NKCC2 = FUR_TUB / (FUR_TUB + NKCC2_IC50);
double E_NE    = Emax_map * NE_CENT / (NE_CENT + EC50_ne);       // raw amount, not NE_CENT/Vd_ne
double E_NAC_gsh = Emax_gsh * NAC_CENT / (NAC_CENT + EC50_nac);  // raw amount, not NAC_CENT/Vd_nac
```

The refactored file exposes a genuine, standard concentration for each
compound per the naming convention — `C_FUR = CENT_FUR/V1_FUR`,
`C_NE = CENT_NE/V1_NE`, `C_NAC = CENT_NAC/V1_NAC` — for driveability and
manifest discoverability, **but** `EFFECT_FUR`, `EFFECT_NE`,
`EFFECT_NAC_GSH`, and `EFFECT_NAC_NO` all still take their Hill argument
from the raw compartment content (`CENT_FUR`/`CENT_NE`/`CENT_NAC`
directly, not `C_FUR`/`C_NE`/`C_NAC`), matching the original exactly.
This is a deliberate choice, not an oversight: substituting `C_<STEM>`
into those Hill terms would change the simulated trajectory (each
`Vd` is ≥8, so dividing by it moves every effect term far from where the
original always ran), which the guide's "never loosen a fit... to force a
pass" principle rules out. Disclosed here per the guide's "a poor fit is
a finding, not a bug to hide" instruction, generalized to this
unit-consistency quirk. A downstream driver wanting to substitute
`C_FUR`/`C_NE`/`C_NAC` directly into the disease effect should be aware
the original's own effect terms are calibrated to compartment *amount*,
not concentration.

## Build-compatibility fixes (all logged upstream, `translations/UPSTREAM_ISSUES.md` #125-#127)

Verifying against the original required patching a **scratch-only**
copy of the original's own extracted DSL (never the checked-in file) for
three independent, pre-existing defects, confirmed via `POST
/model_manifest`/`POST /run_simulation` against the untouched original
alone:

1. **`KA_FUR` used in `$ODE` but never declared in `$PARAM`** (#125). The
   raw `code <- '...'` string does not compile at all; the original's own
   R wrapper only ever runs a `gsub()`-patched copy (`code2`) that
   text-substitutes both `KA_FUR` occurrences for the literal `1.5`
   before `mcode()`. Because this is furosemide's own absorption-rate
   constant (in scope for this refactor, not an incidental defect), the
   refactored file resolves it as part of the refactor itself:
   `KA_FUR = 1.5` is a genuine, explicit `$PARAM` entry (matching the
   `KA_<STEM>` naming-convention role and the original's own hardcoded
   value) — no `gsub()` patch needed in the refactored file at all.
2. **`$CAPTURE` duplicates sixteen `$CMT` compartment names** (#126) —
   `FUR_CENT NE_CENT NAC_CENT ATP ROS GSH TCV IL6 TNFa NGAL KIM1 SCR CysC
   TGFb MYO FIBROSIS REPAIR_CAP` are all `$CMT` compartments, which
   mrgsolve 2.0.1 rejects outright (`compartment should not be in
   $CAPTURE`). The refactored file's `$CAPTURE` keeps only the file's
   genuinely non-compartment derived quantities; every removed name is
   still reported (unchanged) as a `$CMT` compartment, so nothing becomes
   unobservable via the API's `outputPaths`.
3. **The DSL alone never initializes any compartment away from 0**
   (#127, the most consequential of the three). `$MAIN` never assigns a
   single `<CMT>_0 = ...;`; every one of the original's own seven
   scenarios only starts from physiologically sensible values (GFR=100,
   TCV=1.0, etc.) because the R wrapper's `mod %>% init(...)` call runs
   *after* `mcode()`, entirely outside the quoted DSL string. A client
   driving only the bare `model_content` (as the qspserver API does)
   silently gets every compartment starting at 0 — confirmed: with only
   fixes #125-126 applied and no `init()` override, `GFR` stays at 0 for
   the entire 168 h horizon instead of declining from 100. The refactored
   file's `$MAIN` now assigns `<CMT>_0 = ...;` for all fifteen non-drug
   compartments, sourced from the original's own matching `$PARAM`
   baseline where one exists (`GFR_0=GFR0`, `ATP_0=ATP0`, `ROS_0=ROS0`,
   `GSH_0=GSH0`, `TCV_0=TCV0`, `IL6_0=IL6_base`, `TNFa_0=TNFa_base`,
   `SCR_0=SCR0`, `CysC_0=CysC0`, `REPAIR_CAP_0=EGF0`, `TGFb_0=TGFb0`) or
   hardcoded to the exact value the R wrapper's own `init()` call uses
   where no matching `$PARAM` exists (`NGAL_0=5`, `KIM1_0=0.5`,
   `MYO_0=0.0`, `FIBROSIS_0=0.0`). The five drug compartments correctly
   stay at mrgsolve's default 0, matching the R wrapper's own `init()`
   call there too. None of the three fixes touch any numeric value the
   original's own scenarios actually produce — proof is the verification
   result below.

All three fixes are syntax/wiring only — no compound's own PK/PD math
was changed by any of them; each is disclosed in the matching
`UPSTREAM_ISSUES.md` entry with the exact patch applied.

## Verification

Ran four of the original's own seven R-script scenarios (chosen to cover
each compound alone and two compounds combined) through both a
verification copy of the original (patched only for #125-127, no
renaming) and the refactored DSL, via the qspserver `mrgsolve_api`
container (`POST /model_manifest`, `POST /run_simulation`), full
168 h horizon at `delta=0.5` (the original's own grid) both times:

| Scenario | Dosing (original's own) | Compounds exercised |
|---|---|---|
| 2 | Furosemide 40 mg oral q12h ×14, `GUT_FUR`, start 2 h | FUR |
| 3 | NE continuous infusion, `CENT_NE`, rate 0.1, start 0 h | NE |
| 4 | NAC loading+maintenance infusion, `CENT_NAC`, rate 1000, start −2 h | NAC |
| 6 | NE infusion (start 0 h) + furosemide 40 mg q12h ×14 (`GUT_FUR`, start 6 h), sepsis-AKI + CRRT on | FUR + NE together |

All 28 shared outputs (every `$CMT` compartment plus `AKI_stage_out`,
`eGFR`, `UO_mLkgh`, `FE_Na`, `diuresis`, and the renamed
`E_NKCC2`→`EFFECT_FUR`, `E_NE`→`EFFECT_NE`, `E_NAC_gsh`→`EFFECT_NAC_GSH`)
matched **exactly** across all four scenarios: **max abs diff = 0.0** on
every output, every timepoint (339 points for scenarios 2/3/6, 338 for
scenario 4 which starts at t=−2). No tolerance was loosened to reach this
— it is a genuine bit-exact match, consistent with all three compounds
being pure structural reorganizations (Archetypes 1 and 3), no
Hill-fitting anywhere.

`POST /model_manifest` on the refactored DSL confirms every renamed PK
parameter (`KA_FUR, F_FUR, CL_FUR, V1_FUR, V2_FUR, Q_FUR, EC50_FUR,
EMAX_FUR, GAMMA_FUR, CL_NE, V1_NE, EMAX_NE, EC50_NE, GAMMA_NE, CL_NAC,
V1_NAC, EMAX_NAC_GSH, EC50_NAC, EMAX_NAC_NO, GAMMA_NAC`) appears in
`parameters` (77 total, fixed values unchanged from the original except
the three new `GAMMA_*=1.0`/`EMAX_FUR=1.0` and the newly-explicit
`KA_FUR=1.5`). `C_FUR, EFFECT_FUR, C_NE, EFFECT_NE, C_NAC,
EFFECT_NAC_GSH, EFFECT_NAC_NO` are state-dependent (derived inside
`$MAIN` from live compartment content each interval), so per this fork's
established precedent (`drug-induced-liver-injury/dili_refactor_notes.md`,
`acute-pancreatitis/ap_refactor_notes.md`) they are not `$PARAM`
defaults — confirmed discoverable instead in the manifest's
`outputPaths` (32 total: 5 renamed PK compartments + 15 disease
compartments + 5 original `$TABLE` outputs + 7 new
`C_<STEM>`/`EFFECT_<STEM>` quantities).

`$MAIN`-vs-`$ODE` placement: every derived quantity the original computed
in `$MAIN` (`inj_stim`, `FUR_TUB`/`E_NKCC2`/`diuresis`,
`E_NE`/`RBF`/`E_RBF_gfr`, `E_NAC_gsh`, `E_NAC_no`, `inj_composite`,
`nfkb`, etc.) stays in `$MAIN` in the refactored file, per the guide's
"keep a calculation in the block the original used it in" instruction —
none of it was moved into `$ODE`. `C_FUR`/`C_NE`/`C_NAC` are new
quantities (not present in the original at all), so placing them in
`$MAIN` alongside the effect terms they're grouped with is not a "moved"
calculation.

No separate `.cpp` extraction was left in the repository — the
deliverable list is just `aki_mrgsolve_model_refactored.R` and this notes
file plus the census update. The quoted DSL block was extracted to a
scratch file for every API call above (confirmed byte-identical to the
`code <- '...'` string in the delivered file, modulo R's own `\'`→`'`
string-escape unescaping, which is a no-op on the actual DSL text), then
the scratch copy was deleted along with the whole scratch directory.

## R syntax check

`Rscript -e 'cat(length(parse("aki_mrgsolve_model_refactored.R")))'`
succeeds (44 top-level expressions parsed) — an R-syntax check only (the
DSL itself is a quoted string, opaque to R's parser), not a local
mrgsolve build, per this task's "do not use local mrgsolve/R for
verification" instruction; all mrgsolve-level verification above was
done exclusively through the qspserver API.
