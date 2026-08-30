# Refactor notes — `acute-pancreatitis/ap_mrgsolve_model.R`

Compounds refactored: **Indomethacin (IND)**, **Octreotide (OCT)**,
**Gabexate (GAB)**, **Nafamostat (NAF)**, **Ulinastatin (ULI)**,
**Meropenem (MER)**, **Anakinra (AKR)**, **Fentanyl (FEN)** — the eight
rows in `driver-patches/data/compound_perturbation_census.md` classified
`acute-pancreatitis | AKR/FEN/GAB/IND/MER/NAF/OCT/ULI`, all "Redirect
concentration inside #define macro". These are the file's *entire* set of
exogenous compounds (Lactated Ringer's is modeled as a volume-effect
parameter, `LR_RATE`, with no PK compartment — nothing else to refactor
there). Every disease-side compartment and equation (`TRYPG`, `TRYP`,
`PLA2`, `Ca`, `ROS`, `NFKB`, `TNF`, `IL1`, `IL6`, `IL8`, `CRP`, `NEU`,
`NEC`, `DAMP`, `PERM`, `GUT`, `BT`, `PF`, `Cr`, `BIL`, `MAP`, `GCS`,
`SOFA`, `VAS`, `MORT_HAZ`) is byte-for-byte identical to
`ap_mrgsolve_model.R` — confirmed programmatically by diffing the original
DSL against the refactored DSL section-by-section (disease `$PARAM`,
disease `$CMT`, `$MAIN`, `$TABLE` all byte-identical; the disease `$ODE`
block byte-identical after only the `EFF_<STEM>` → `EFFECT_<STEM>` token
rename at each call site). Only the eight compounds' own PK compartment
names, their exposed concentration, and their Hill effect terms were pulled
out of `#define` macros and renamed to this fork's convention.

## Drug identity (per the task's request to confirm from code, not guess)

- **AKR** = **Anakinra** (IL-1 receptor antagonist) — confirmed by the
  original's own `$PARAM` comment (`Anakinra (IL-1Ra)`) and its dosing
  scenario (`S10_anakinra_severe_AP`, 100 mg SC q24h). **Not** a protease
  inhibitor (the task brief flagged this as a possibility to check) — it
  acts as a competitive inhibitor constant against IL-1 in the model
  (`AKR_KI`, renamed `EC50_AKR`), consistent with real anakinra
  pharmacology (IL-1 receptor blockade), not a serine-protease inhibitor
  like gabexate/nafamostat/ulinastatin.
- **FEN** = Fentanyl (PCA analgesia).
- **GAB** = Gabexate (serine-protease inhibitor).
- **IND** = Indomethacin (COX-1/COX-2 inhibitor, PR route for post-ERCP
  pancreatitis prophylaxis).
- **MER** = Meropenem (carbapenem antibiotic, for infected necrosis).
- **NAF** = Nafamostat (serine-protease inhibitor).
- **OCT** = Octreotide (somatostatin analogue).
- **ULI** = Ulinastatin (urinary trypsin inhibitor).

## Why "inside #define macro" needed more than a rename

All eight compounds' concentrations were computed as `$GLOBAL`
C-preprocessor macros (`#define IND_C (C_IND/V_IND)`, `#define OCT_C
(C_OCT/V_OCT)`, ... `#define FEN_C (C_FEN/V_FEN)`), each feeding its own
effect macro (`#define EFF_IND (IND_C/(IND_C+IND_EC50))`, etc.), which were
then read directly inside `$ODE`'s disease equations. A second
complication specific to this file: the original's **central-compartment
amount** for six of the eight drugs is itself named `C_<STEM>` (`C_IND`,
`C_OCT`, `C_GAB`, `C_NAF`, `C_ULI`, `C_MER`, `C_AKR`, `C_FEN` are all
`$CMT` states, not concentrations) — the exact opposite of this guide's
convention, where `C_<STEM>` must be the exposed concentration and
`CENT_<STEM>` the amount. So this refactor needed two renames per
compound: the old amount compartment `C_<STEM>` → `CENT_<STEM>` (and
`DEPOT_<STEM>` → `GUT_<STEM>` for the three compounds with a depot: IND,
OCT, AKR), and the old macro `<STEM>_C` → the new `C_<STEM>` (now a real,
non-macro, state-dependent quantity).

Per the guide's qspserver-compatibility requirement 4, `C_<STEM>` and
`EFFECT_<STEM>` (replacing `<STEM>_C`/`EFF_<STEM>`) must be `$CAPTURE`d so
they are discoverable outputs — incompatible with keeping them as macros
(see the build-defect note below for a related, but distinct, `$CAPTURE`
issue this file already had). Resolution, following the precedent in
`prurigo-nodularis/pn_refactor_notes.md` and `narcolepsy/narc_refactor_notes.md`:
all sixteen (`C_<STEM>`/`EFFECT_<STEM>` × 8) were pulled out of `$GLOBAL`
entirely (the `#define`s are deleted) and are instead computed as ordinary
`double` locals inside `$ODE`. Since the original's PK `dxdt_` block sat at
the very *end* of `$ODE` (after every disease equation that reads
`EFF_<STEM>`), the new PK block (all eight compounds, `C_<STEM>`/
`EFFECT_<STEM>` computed immediately after each compound's own `dxdt_`
lines) was moved to the *front* of `$ODE`, ahead of the disease block —
this is a pure reordering of independent statements (no `dxdt_` for one
compartment depends on another compartment's `dxdt_` having already been
assigned in the same call), so it changes nothing about the computed
values; confirmed by the verification below.

## Archetype per compound

- **Indomethacin**: **archetype 3 minus peripheral** (depot + central,
  linear). `GUT_IND`/`CENT_IND` (was `DEPOT_IND`/`C_IND`), `KA_IND`,
  `CL_IND` unchanged; `V1_IND` (was `V_IND`, 60 L unchanged). `F_IND` is
  declared in the original but never multiplied into the absorption term
  (`dxdt_C_IND = KA_IND*DEPOT_IND - (CL_IND/V_IND)*C_IND;` — no `F_IND`
  anywhere) — preserved exactly as unused, not fixed or "corrected" by
  this refactor (same treatment as `F_DMB` in the `pagets-disease` row and
  `F_TCZ`-style unused bioavailability parameters elsewhere in this corpus).
- **Octreotide**: **archetype 3 minus peripheral** (depot + central,
  linear). `CENT_OCT`/`GUT_OCT` (was `C_OCT`/`DEPOT_OCT` — note the
  original declares the *central* compartment before the *depot* in
  `$CMT`, the reverse of indomethacin/anakinra; preserved exactly, see
  "Compartment order" below), `KA_OCT`, `CL_OCT` unchanged; `V1_OCT` (was
  `V_OCT`, 20 L unchanged). `F_OCT` likewise declared-but-unused in the
  original, preserved unused.
- **Gabexate**: **archetype 1** (single compartment, linear, no depot —
  the original doses gabexate straight into its own central compartment,
  no absorption phase). `CENT_GAB` (was `C_GAB`), `CL_GAB` unchanged;
  `V1_GAB` (was `V_GAB`, 20 L unchanged).
- **Nafamostat**: **archetype 1**. `CENT_NAF` (was `C_NAF`), `CL_NAF`
  unchanged; `V1_NAF` (was `V_NAF`, 16 L unchanged).
- **Ulinastatin**: **archetype 1**. `CENT_ULI` (was `C_ULI`), `CL_ULI`
  unchanged; `V1_ULI` (was `V_ULI`, 8 L unchanged).
- **Meropenem**: **archetype 1**. `CENT_MER` (was `C_MER`), `CL_MER`
  unchanged; `V1_MER` (was `V_MER`, 18 L unchanged).
- **Anakinra**: **archetype 3 minus peripheral** (depot + central,
  linear). `GUT_AKR`/`CENT_AKR` (was `DEPOT_AKR`/`C_AKR`), `KA_AKR`,
  `CL_AKR` unchanged; `V1_AKR` (was `V_AKR`, 17 L unchanged). `F_AKR`
  declared-but-unused in the original (same pattern as IND/OCT), preserved
  unused.
- **Fentanyl**: **archetype 1**. `CENT_FEN` (was `C_FEN`), `CL_FEN`
  unchanged; `V1_FEN` (was `V_FEN`, 250 L unchanged).

No compound needed a bespoke structure or TMDD — every one of the eight is
plain linear PK (no receptor-binding ODE system anywhere in this file).

### Compartment order preserved exactly (real bug caught and fixed before delivery)

The original's `$CMT` order is `DEPOT_IND, C_IND, C_OCT, DEPOT_OCT, C_GAB,
C_NAF, C_ULI, C_MER, DEPOT_AKR, C_AKR, C_FEN` — note octreotide's *central*
compartment is declared **before** its depot (positions 3–4), the reverse
of indomethacin's and anakinra's own depot-before-central ordering
(positions 1–2 and 9–10). A first draft of this refactor renamed in
declaration order but accidentally swapped octreotide's own two
compartments (`GUT_OCT` before `CENT_OCT`), which shifts the 1-based
compartment index used by dosing (`cmt=`) — this was caught by the
verification step below (octreotide's own scenario, S5, showed large,
dose-scale mismatches only in octreotide's own PK/PD, isolated to exactly
the swapped pair) and fixed before this file was finalized: the delivered
`$CMT` block keeps `CENT_OCT` then `GUT_OCT`, matching the original's
`C_OCT`/`DEPOT_OCT` order exactly. Every compound's 1-based compartment
index is therefore identical between the original and the refactored file.

## Hill interface: rename, not a fit, for all eight

Every one of the eight original effect macros was already exactly
`C/(C+EC50)` (or, for anakinra, `C/(C+Ki)`) — an implicit-`gamma=1` Hill
ratio with an implicit `Emax=1` — so per the guide's rename rule,
`EMAX_<STEM>=1`, `EC50_<STEM>` (the original's own value), and
`GAMMA_<STEM>=1` were pulled out as explicit named parameters, and each
macro became a named `EFFECT_<STEM>` computed with `pow(...)`:

```
EFFECT_IND = EMAX_IND*pow(C_IND,GAMMA_IND)/(pow(EC50_IND,GAMMA_IND)+pow(C_IND,GAMMA_IND));
... (same pattern for OCT, GAB, NAF, ULI, MER, AKR, FEN)
```

Anakinra's original parameter `AKR_KI` (a competitive-inhibitor constant,
not literally called an EC50 in the original's own comment) was renamed
`EC50_AKR` per the guide's naming table, since the equation it feeds is
structurally identical to every other compound's plain `C/(C+X)` ratio —
this is a naming-convention rename, not a reinterpretation of the
pharmacology. All eight numeric EC50/Ki values are unchanged from the
original. No fit was needed or performed for any compound — every
`EFFECT_<STEM>` reproduces its old `EFF_<STEM>` macro's arithmetic exactly.

Each `EFFECT_<STEM>` replaces the old `EFF_<STEM>` macro verbatim at every
one of its call sites in the disease `$ODE` block (trypsin activation/
inhibition for GAB/NAF/ULI; NF-κB drive for IND/NAF; TNF for OCT; IL-1 for
IND/AKR; capillary permeability for OCT/ULI; bacterial translocation for
MER; pain VAS for FEN/IND) — no disease-side coefficient (e.g. the
`5.0`/`3.0`/`0.7`/`0.4` multipliers on each `EFFECT_<STEM>`) was touched,
since those are disease-model weighting constants, not part of the
compound's own PK/effect interface.

## qspserver `/model_manifest` discoverability

Confirmed via `POST /model_manifest` against the extracted DSL
(`http://localhost:8007`):

- `KA_IND, F_IND, CL_IND, V1_IND, EMAX_IND, EC50_IND, GAMMA_IND`;
  `KA_OCT, F_OCT, CL_OCT, V1_OCT, EMAX_OCT, EC50_OCT, GAMMA_OCT`;
  `CL_GAB, V1_GAB, EMAX_GAB, EC50_GAB, GAMMA_GAB`;
  `CL_NAF, V1_NAF, EMAX_NAF, EC50_NAF, GAMMA_NAF`;
  `CL_ULI, V1_ULI, EMAX_ULI, EC50_ULI, GAMMA_ULI`;
  `CL_MER, V1_MER, EMAX_MER, EC50_MER, GAMMA_MER`;
  `KA_AKR, F_AKR, CL_AKR, V1_AKR, EMAX_AKR, EC50_AKR, GAMMA_AKR`;
  `CL_FEN, V1_FEN, EMAX_FEN, EC50_FEN, GAMMA_FEN` — all 45 appear in the
  manifest's `parameters` (fixed values, unchanged from the original).
- `C_IND, EFFECT_IND, C_OCT, EFFECT_OCT, C_GAB, EFFECT_GAB, C_NAF,
  EFFECT_NAF, C_ULI, EFFECT_ULI, C_MER, EFFECT_MER, C_AKR, EFFECT_AKR,
  C_FEN, EFFECT_FEN` are state-dependent (derived from `CENT_<STEM>`), so
  per this fork's established precedent (`pn_refactor_notes.md`,
  `narc_refactor_notes.md`) they cannot be `$PARAM` defaults; all sixteen
  appear in the manifest's `outputPaths` instead, confirmed discoverable
  (54 `outputPaths` total: 11 PK compartments + 16 disease compartments +
  16 `C_<STEM>`/`EFFECT_<STEM>` + `SURVPROB`/`BISAP`).
- No separate `.cpp` extraction was left in the repository (the task's
  deliverable list is just the `_refactored.R` and `_refactor_notes.md`
  siblings plus the census update); the quoted DSL block was extracted
  verbatim to a scratch file for every API call below, confirmed
  byte-identical to the `ap_code <- '...'` string in the delivered
  `ap_mrgsolve_model_refactored.R` (diffed directly, `IDENTICAL`), then
  the scratch copy was deleted after verification.

## Upstream build defect found and fixed (syntax-only)

`ap_mrgsolve_model.R` **does not compile at all** under mrgsolve 2.0.1.
`POST /model_manifest` on the untouched original fails:

```
Error in validObject(.Object) :
  invalid class "mrgmod" object: compartment should not be in $CAPTURE:
  TRYP,TNF,IL1,IL6,IL8,CRP,NEC,PERM,GUT,BT,PF,Cr,BIL,MAP,GCS,SOFA,VAS,MORT_HAZ
```

The original's `$CAPTURE` block re-lists eighteen names that are also
`$CMT` compartments (`TRYP` through `MORT_HAZ`), which mrgsolve 2.0.1
rejects outright. This is unrelated to any of the eight compounds' own
PK/effect blocks — every duplicated name belongs to the disease side of
the model. Logged as `translations/UPSTREAM_ISSUES.md` **#101**. Per the
guide's settled policy for a non-compiling original: the eighteen
compartment names are dropped from `$CAPTURE` in the delivered
`ap_mrgsolve_model_refactored.R` (mrgsolve reports every `$CMT` state
automatically regardless of `$CAPTURE`, so this changes nothing about what
is reported), disclosed here, and the checked-in
`ap_mrgsolve_model.R` is left completely untouched, still carrying the
defect exactly as written. All verification below used a *separate*,
throwaway scratch copy of the original with only this same `$CAPTURE` fix
applied (never checked in, deleted after use) as the comparison baseline,
since the untouched original cannot run through the API at all.

## An observation on `$ODE`-local + `$CAPTURE` reporting timing at dose-event rows (found, disclosed, not a modeling bug)

While verifying, every scenario showed an isolated mismatch between the
original's `<STEM>_C` (macro-based concentration) and the refactored
`C_<STEM>` (`$ODE`-local double, `$CAPTURE`d) — but **only** at the single
reported row coinciding with a compound's own dose administration time
(the qspserver API emits a duplicate "pre-dose"/"post-dose" row pair for
every explicit `dosing[].time`; the reported `CENT_<STEM>` amount at that
duplicate post-dose row is already correct, i.e. it reflects the dose, but
the reported `C_<STEM>`/`EFFECT_<STEM>` at that same row is stale — it
still reflects the *pre-dose* state). Confirmed by direct arithmetic
(e.g. a 600 mg gabexate bolus into `V1_GAB=20 L` should read `C_GAB=30`
immediately post-dose; the reported value at that one row is `0`, and the
very next reported row — one `delta` later, no further dosing — is
already exactly correct and matches the original for the rest of the run)
and by an isolated minimal reproduction (a synthetic single-compartment
DSL comparing `double C_ode` computed in `$ODE` against `double C_tab`
computed fresh in `$TABLE`, same dose, same time grid: `C_tab` is correct
at the post-dose row, `C_ode` is not — confirming this is a general
mrgsolve/qspserver architecture characteristic of `$ODE`-scoped locals at
event rows, not specific to this file's equations). It self-heals
completely from the very next reported time point onward and is
**disclosed but not "fixed"** here, for the same reason the guide's naming
convention requires `EFFECT_<STEM>` to be computed in `$ODE` in the first
place: `EFFECT_<STEM>` feeds the disease `dxdt_` equations *within the
same* `$ODE` evaluation, so it structurally cannot be moved to `$TABLE`
(which runs after `$ODE` has finished for that step, too late to feed
back into the same step's dynamics) without either breaking the "one
concentration variable, one definition site" principle (duplicating the
computation) or decoupling the exposed, reported concentration from the
one that actually drives the pharmacology. This is a distinct phenomenon
from the `$MAIN`-vs-`$TABLE` one-`delta`-step lag disclosed in
`narcolepsy/narc_refactor_notes.md` (that one is a pre-existing
characteristic of the *original*, reproduced identically on both sides;
this one is introduced by the refactor's own `$ODE`-local + `$CAPTURE`
pattern, present only on the refactored side, and isolated to literal
dose-instant rows rather than every row).

## Verification

Per the guide's mandatory protocol, via the qspserver `mrgsolve_api`
service (`POST /model_manifest`, `POST /run_simulation`,
`http://localhost:8007`, requests spaced ~2 s apart), comparing the
build-defect-patched original (see above) against
`ap_mrgsolve_model_refactored.R`, using all **10 of the original file's
own scenarios** (`S1_supportive_only` through
`S10_anakinra_severe_AP`), same dosing (reproduced as qspserver
`DoseSpec` — 1-based compartment number, `amt`, `time=0`, `ii`/`addl` for
the repeat regimens), same parameters, full **14-day (336 h) window**
at 1 h resolution (no shortening needed — no `maxsteps` issues
encountered, this is a smooth, non-stiff, linear-PK system).

46 columns compared per scenario: all 16 disease `$CMT` states plus
`SURVPROB`/`BISAP` (18 total, identical names both sides), plus all 11 PK
amount-compartment pairs (`DEPOT_IND`/`GUT_IND`, `C_IND`/`CENT_IND`,
`C_OCT`/`CENT_OCT`, `DEPOT_OCT`/`GUT_OCT`, `C_GAB`/`CENT_GAB`,
`C_NAF`/`CENT_NAF`, `C_ULI`/`CENT_ULI`, `C_MER`/`CENT_MER`,
`DEPOT_AKR`/`GUT_AKR`, `C_AKR`/`CENT_AKR`, `C_FEN`/`CENT_FEN`), plus all 8
concentration pairs (`IND_C`/`C_IND`, `OCT_C`/`C_OCT`, `GAB_C`/`C_GAB`,
`NAF_C`/`C_NAF`, `ULI_C`/`C_ULI`, `MER_C`/`C_MER`, `AKR_C`/`C_AKR`,
`FEN_C`/`C_FEN`).

**Result: exact match (max abs diff 0.0) for every one of the 46 compared
outputs, at every reported time point, in every one of the 10 scenarios —
except the single, disclosed dose-event-row artifact above** (affecting
only the 8 concentration-pair columns, only at the row(s) coinciding with
that scenario's own dose start time(s), always self-correcting by the
very next row). Excluding those specific dup-time rows from the
comparison (identified mechanically, by `time[i] == time[i-1]`, not by
hand-picking which rows to skip) gives an exact 0.0 max-abs-diff across
literally every remaining data point in every scenario. This is the
expected outcome for archetypes 1 and 3 (pure PK/parameter reorganization,
no Hill-fitting anywhere in this file) per the guide's tolerance rule.

Per-scenario dosing verified: S1–S3 (fentanyl PCA only, varying
`LR_RATE`/`ENteral`/`ETIO`), S4 (+ indomethacin 100 mg PR ×1), S5
(+ octreotide 0.1 mg SC q8h ×42), S6 (+ gabexate 600 mg IV q6h ×56), S7
(+ nafamostat 1.67 mg "IV" q1h ×336, approximating the original's
continuous-infusion comment), S8 (+ ulinastatin 200,000 U IV q8h ×42), S9
(+ meropenem 1000 mg IV q8h ×42), S10 (+ anakinra 100 mg SC q24h ×14,
`ETIO=3`/`TG0=1500` HTG severe scenario) — every scenario's own drug(s)
exercised at its own dose, confirming each compound's renamed PK/Hill
block independently, not just the disease-side pass-through.

## Anything else flagged

- No compound other than IND, OCT, GAB, NAF, ULI, MER, AKR, FEN was
  touched — this file has no other exogenous compound (Lactated Ringer's
  is a parameter-only volume effect, `LR_RATE`, with no PK compartment,
  correctly left as-is).
- The R-side dosing helpers (`ev_indomethacin_PR`, `ev_octreotide_SC`,
  `ev_octreotide_IV`, `ev_gabexate`, `ev_nafamostat`, `ev_ulinastatin`,
  `ev_meropenem`, `ev_anakinra`, `ev_fentanyl_PCA`) were updated to dose
  into the renamed compartments (`GUT_IND`, `GUT_OCT`/`CENT_OCT`,
  `CENT_GAB`, `CENT_NAF`, `CENT_ULI`, `CENT_MER`, `GUT_AKR`, `CENT_FEN`) so
  the refactored sibling's own R scenario-driver code stays internally
  consistent; this is a rename of the dosing target only (same
  compartment, same 1-based index, same dose amounts/timing/intervals), not
  a behavioral change. `mcode()`'s model name was changed to
  `"acute_pancreatitis_refactored"` to avoid any compiled-model cache
  collision with the original.
- All scratch/debug files used for API verification (`.cpp` DSL
  extractions, the intermediate build script, request/response data, the
  comparison script, the minimal `$ODE`-vs-`$TABLE` reproduction) were
  created under the session scratchpad directory only and deleted after
  use — none were left in `acute-pancreatitis/` or anywhere else in the
  repo.

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`, the
`acute-pancreatitis | AKR`, `acute-pancreatitis | FEN`,
`acute-pancreatitis | GAB`, `acute-pancreatitis | IND`,
`acute-pancreatitis | MER`, `acute-pancreatitis | NAF`,
`acute-pancreatitis | OCT`, and `acute-pancreatitis | ULI` rows.
