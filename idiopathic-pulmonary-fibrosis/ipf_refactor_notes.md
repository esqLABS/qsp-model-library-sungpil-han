# Refactor notes — idiopathic-pulmonary-fibrosis / two compounds (PIR, NIN)

Covers `idiopathic-pulmonary-fibrosis/ipf_mrgsolve_model.R`. Both compounds
this file models were in scope — Pirfenidone and Nintedanib, per their
existing census rows (`N` and `P`) in
`driver-patches/data/compound_perturbation_census.md` — there are no other
compounds in this file to leave untouched. Every disease-side compartment
and equation (AEC-II, TGF-β1, M2 macrophage, ROS, fibroblast/myofibroblast,
collagen/MMP/TIMP, FVC, DLCO) is copied verbatim from the original into
`ipf_mrgsolve_model_refactored.R`; only the two drugs' own PK blocks and
their named Hill effect terms were renamed to the fork's convention.

## Compound identity and target correction (task-flagged)

The census listed both compounds only by single-letter code. Read from the
code directly:

- **`N` = Nintedanib** (`ka_N`/`F_N`/`CL_N`/`V1_N`/`V2_N`/`Q_N` match the
  file's own header comment "Nintedanib PK: Stopfer (2011) Clin
  Pharmacokinet; F=4.7%, t½=10h, ka=0.8/h" exactly — `F_N=0.047`,
  `ka_N=0.80`). Target in this file: fibroblast proliferation inhibition,
  modeled via an EC50/nM Hill term (`Emax_N_F`/`EC50_N_F=15.0 nM`); the
  file's own plot annotation labels this "IC50 ~20 nM (FGFR1)". Real-world
  nintedanib is a triple angiokinase inhibitor (VEGFR/FGFR/PDGFR); recorded
  in the census as **FGFR/PDGFR (fibroblast proliferation)** rather than a
  single receptor, since that is what the code and its own annotation
  actually represent.
- **`P` = Pirfenidone** (`ka_P=1.74`, `F_P=0.81` match the header's
  "Pirfenidone PK: Rubino (2009) Clin Pharmacokinet; F=81%, t½=2.4h,
  ka=1.74/h" exactly). The census's Target/Pathway column read **"PCSK9"**
  for this row — confirmed wrong: PCSK9 (a cholesterol/LDL-receptor
  regulator) appears nowhere in this file and has no established
  pharmacology in IPF or with pirfenidone. Corrected in the census to
  **TGF-β1 / anti-fibrotic-anti-inflammatory (no single defined
  receptor)** — pirfenidone's real, still-not-fully-elucidated mechanism,
  matching what the code models: three separate anti-fibrotic/
  anti-inflammatory effects (TGF-β production, M2/cytokine polarization,
  ROS/oxidative stress), not one receptor-ligand interaction. This looks
  like the same kind of cross-contaminated lookup-dictionary error the
  census's own limitations section warns about (Target/Pathway is `?` for
  ~70% of the corpus; where populated, apparently sometimes wrong).

## Archetype determined per compound

| Compound | Archetype | Structure |
|---|---|---|
| Pirfenidone (`PIR`) | 3 | depot + central + peripheral, linear CL/Q/V (already in that form in the original) |
| Nintedanib (`NIN`) | 3 | depot + central + peripheral, linear CL/Q/V (already in that form in the original) |

Both were already clean linear two-compartment-with-depot PK in the
original (`ka`/`F`/`CL`/`V1`/`V2`/`Q` — no micro-constant conversion
needed). Neither compound involves TMDD or any receptor-binding ODE
system.

## Naming map

| Original | Refactored | Value | Role |
|---|---|---|---|
| `DEPOT_P` / `CENT_P` / `PERI_P` | `GUT_PIR` / `CENT_PIR` / `PERI_PIR` | -- | depot/central/peripheral |
| `ka_P`, `F_P`, `CL_P`, `V1_P`, `V2_P`, `Q_P` | `KA_PIR`, `F_PIR`, `CL_PIR`, `V1_PIR`, `V2_PIR`, `Q_PIR` | 1.74, 0.81, 8.4, 20.4, 14.0, 3.6 | unchanged values, renamed |
| `DEPOT_N` / `CENT_N` / `PERI_N` | `GUT_NIN` / `CENT_NIN` / `PERI_NIN` | -- | depot/central/peripheral |
| `ka_N`, `F_N`, `CL_N`, `V1_N`, `V2_N`, `Q_N` | `KA_NIN`, `F_NIN`, `CL_NIN`, `V1_NIN`, `V2_NIN`, `Q_NIN` | 0.80, 0.047, 22.0, 730.0, 900.0, 15.0 | unchanged values, renamed |
| `Cp_P` | `C_PIR` | -- | µg/mL, the exposed PIR concentration |
| `Cn_nM` | `C_NIN` | -- | nM, the exposed NIN concentration |
| `Cn_N` | *(dropped)* | -- | dead code — computed in `$MAIN`, never read anywhere downstream in the original (only `Cn_nM` feeds the effect term); confirmed by grep, no other reference in the file |
| `EC50_P_TGFb` | `EC50_PIR_TGFB` | 30.0 (µg/mL) | Hill EC50, TGF-β term |
| `Emax_P_TGFb` | `EMAX_PIR_TGFB` | 0.65 | Hill Emax, TGF-β term (also reused for the M2 term, see below) |
| `EC50_P_M2` | `EC50_PIR_M2` | 25.0 (µg/mL) | Hill EC50, M2 term |
| `Emax_P_ROS` | `EMAX_PIR_ROS` | 0.45 | Hill Emax, ROS term |
| -- (none) | `GAMMA_PIR_TGFB`, `GAMMA_PIR_M2`, `GAMMA_PIR_ROS` (new) | 1.0 each | Hill coefficients [math-implied; the original had no explicit exponent for any of the three] |
| `EC50_N_F` | `EC50_NIN` | 15.0 (nM) | Hill EC50 |
| `Emax_N_F` | `EMAX_NIN` | 0.70 | Hill Emax |
| -- (none) | `GAMMA_NIN` (new) | 1.0 | Hill coefficient [math-implied] |
| `MW_P` | `MW_PIR` | 185.22 | unused in the original (declared, never referenced) — preserved unused, not repurposed |
| `MW_N` | `MW_NIN` | 539.63 | used in `C_NIN`'s nM conversion |
| `inh_P_TGFb` | `EFFECT_PIR_TGFB` | -- | drives `dxdt_TGFb` |
| `inh_P_M2` | `EFFECT_PIR_M2` | -- | drives `dxdt_M2` |
| `inh_P_ROS` | `EFFECT_PIR_ROS` | -- | drives `dxdt_ROS` |
| `inh_N_F` | `EFFECT_NIN` | -- | drives `dxdt_FIBRO` directly and `dxdt_MYOFIB` at half weight (`EFFECT_NIN * 0.5`, unchanged from the original's `inh_N_F * 0.5`) |

All disease-side parameters/compartments (AEC2, TGF-β1, M2, ROS,
fibroblast/myofibroblast, collagen/MMP/TIMP, FVC/DLCO kinetics) are
completely unchanged in name and value.

## Pirfenidone has three effect terms from one concentration, not one

`EFFECT_PIR_TGFB`, `EFFECT_PIR_M2`, and `EFFECT_PIR_ROS` all read the same
`C_PIR`, matching the "one drug, multiple independent downstream effects"
pattern already documented for dexamethasone
(`cytokine-release-syndrome/crs_refactor_notes.md`) and topiramate
(`essential-tremor/et_refactor_notes.md`). Handled as three named
`EFFECT_PIR_*` variables rather than one shared term, since collapsing
them would lose distinct Emax/EC50 ceilings the original computed
separately for each pathway.

**A parameter-sharing quirk in the original, preserved exactly, not
"fixed":** `inh_P_M2` reused `Emax_P_TGFb` (not its own Emax) with its own
`EC50_P_M2`, and `inh_P_ROS` reused `EC50_P_TGFb` (not its own EC50) with
its own `Emax_P_ROS` — i.e. the TGF-β and M2 terms share one Emax, and the
TGF-β and ROS terms share one EC50, and no separate `Emax_P_M2` or
`EC50_P_ROS` ever existed in the original. `EFFECT_PIR_M2` and
`EFFECT_PIR_ROS` reproduce this exactly (`EFFECT_PIR_M2` reads
`EMAX_PIR_TGFB`, `EFFECT_PIR_ROS` reads `EC50_PIR_TGFB`), rather than
inventing `EMAX_PIR_M2`/`EC50_PIR_ROS` the original never had.

## `$MAIN` vs `$ODE` for `C_<STEM>`/`EFFECT_<STEM>` — a real behavioral trap, not a style choice

The naming convention's own template computes `C_<STEM>` inside `$ODE`.
This file's original computes its PK-derived concentrations and effect
terms (`Cp_P`, `Cn_N`, `Cn_nM`, `inh_P_TGFb`, `inh_P_M2`, `inh_P_ROS`,
`inh_N_F`) in **`$MAIN`**, not `$ODE`. mrgsolve evaluates `$MAIN` once per
requested output/dosing interval using the *start-of-interval* state —
the same one-interval-late-lag behavior documented for
`rheumatoid-arthritis/ra_mrgsolve_model.R` (`UPSTREAM_ISSUES.md` #31 part
2). Confirmed empirically during development: an initial draft moved
`C_PIR`/`C_NIN`/`EFFECT_*` into `$ODE` (matching the template literally,
and matching how most other files in this corpus already had their
original concentration terms structured) — this compiled fine and looked
identical at a glance, but the verification comparison against the
original caught it immediately: `TGFb_norm` diverged by up to ~35% at the
96h/24h grid (e.g. 2.0068 vs 1.4616 at hour 24, Pirfenidone-alone
scenario) even though every parameter value and formula was identical.
Reverted to keep `C_PIR`, `C_NIN`, and all four `EFFECT_*` terms in
`$MAIN`, exactly where the original computed them, which restored an
exact match. Logged as `UPSTREAM_ISSUES.md` #114/#115's sibling
observation is issue #31's pattern recurring in a second file — noted here
rather than re-logged as a new numbered entry, since it is the same
underlying mrgsolve behavior, not a new defect.

**A related, deliberate non-normalization:** `$TABLE`'s own `Cp_pirf`/
`Cn_nint` independently recompute the identical `CENT_PIR/V1_PIR` and
`(CENT_NIN/V1_NIN)/MW_NIN*1e3` formulas a second time, rather than
referencing `$MAIN`'s `C_PIR`/`C_NIN`. This looks like exactly the
"duplicate concentration site" the census flagged for `P` ("Normalize
duplicate concentration sites, then redirect") — and, on inspection, `N`
has the identical structure (`Cn_nint` in `$TABLE` duplicates `Cn_nM` in
`$MAIN`), so the census's "clean single site" classification for `N`
appears to have missed this; the census's own limitations section
already flags this classifier as noisy near category boundaries, and this
looks like one more instance. However, given the `$MAIN`-lag finding
above, these are **not numerically redundant**: `$TABLE` runs at the
actual requested output time using the *current*, non-lagged state, while
`$MAIN`'s `C_PIR`/`C_NIN` are one interval stale. Redirecting `$TABLE`'s
`Cp_pirf`/`Cn_nint` to reference `$MAIN`'s `C_PIR`/`C_NIN` would inject
that lag into what was previously a fresh, unlagged report and silently
change those two outputs' numeric values — a behavior change disguised as
a "normalization". Kept `$TABLE`'s own independent recomputation exactly
as the original had it. The naming convention's "exactly one concentration
variable that PD equations read" is satisfied at the point that matters
(`$MAIN`'s `C_PIR`/`C_NIN`, the only site anything disease-side actually
reads) — `$TABLE`'s duplicate is reporting-only, feeds nothing, and
normalizing it away is not safe here.

## Build-compatibility fix (disclosed; syntax-only, non-numeric)

**The original `ipf_mrgsolve_model.R` does not compile under mrgsolve
2.0.1** — logged as `UPSTREAM_ISSUES.md` **#114**. `$CMT @annotated`
declares all 17 compartments and a separate `$INIT` block immediately
re-declares the same 17 with starting values — mrgsolve 2.0.1 rejects this
as `Duplicated model names`. Per the guide's settled policy, the fix —
deleting `$INIT` and moving its 17 assignments into `$MAIN` using the
`<CMT>_0 = value;` idiom — was applied **directly to the delivered
`ipf_mrgsolve_model_refactored.R`** (and, for comparison purposes only, to
an in-memory-only patched copy of the untouched original's own DSL sent to
the qspserver API; `ipf_mrgsolve_model.R` itself is untouched and still
carries the defect exactly as written).

## A second, more serious pre-existing defect: the disease system diverges within days, not weeks

**Independently of the `$CMT`/`$INIT` fix above, none of the file's own 5
scenarios can actually reach anywhere near its own claimed 52-week
horizon** — logged as `UPSTREAM_ISSUES.md` **#115**. `TGFb_prod` includes
uncapped positive-feedback terms from `M2` and `MYOFIB`, and `M2_prod`/
`M_diff` in turn scale with `TGFb`, closing two feedback loops with no
saturation anywhere. Confirmed with **no drug dosed at all** (Placebo arm
of the `$CMT`/`$INIT`-patched original): `TIMP` already exceeds 3x
baseline by hour ~29, and the mrgsolve/lsoda solver fails with the default
`maxsteps` around hour 112-114 as `TIMP`, `MMP`, and `COLLAGEN` diverge
toward astronomical values (confirmed via a fine probe to `end=113.8`:
`TIMP`~1.4e19). Binary-searched the solver-failure boundary for all 5 of
the original's own scenarios via the qspserver `mrgsolve_api` container:

| Scenario | Solver-failure boundary (h) |
|---|---|
| Placebo (Natural History) | 112-114 |
| Nintedanib 150 mg BID | 183-185 |
| Pirfenidone 267 mg TID (low dose) | 198-200 |
| Pirfenidone 801 mg TID | 202-204 |
| Combination (Pirf + Nint) | 296-298 |

This is entirely disease-side (`TGFb_prod`/`M2_prod`/`M_diff`, untouched
by this refactor) and reproduces identically whether or not either drug is
present — drug suppression of TGF-β/M2/fibroblast activity only delays
the divergence, it does not prevent it. Not something in scope for a
PK/PD refactor to fix; `ipf_mrgsolve_model_refactored.R` reproduces the
identical divergence, unchanged, at the identical time per scenario.

## Verification

**Method.** Extracted the bare DSL text from `ipf_mrgsolve_model.R`
(`code <- '...'`, with the `$INIT`→`$MAIN` build-compat fix above applied
in-memory only) and from the delivered `ipf_mrgsolve_model_refactored.R`
(`ipf_mrgsolve_model_refactored.cpp`, a byte-identical extraction of the
`code <- '...'` block, confirmed by direct comparison). Both run through
the qspserver `mrgsolve_api` container (`http://localhost:8007`),
`POST /model_manifest` (build-time check — confirmed all of `KA_/F_/CL_/
V1_/V2_/Q_` for both compounds, `EC50_/EMAX_/GAMMA_` for all four Hill
terms, `MW_PIR`/`MW_NIN` appear in `parameters`; `C_PIR`, `EFFECT_PIR_TGFB`,
`EFFECT_PIR_M2`, `EFFECT_PIR_ROS`, `C_NIN`, `EFFECT_NIN` all appear in
`outputPaths`) and `POST /run_simulation` (numeric comparison). Requests
were spaced ~2s apart per the API's stated concurrency limits; no failures
or crashes encountered once the shortened window (below) was adopted.

**Shortened verification window (disclosed):** per the guide's
solver-step-count policy, the original file's own 52-week (8736h)
scenarios do not complete under mrgsolve 2.0.1's default solver settings
for any of the 5 scenarios (see `UPSTREAM_ISSUES.md` #115 above) — this is
a pre-existing disease-side defect, not a mismatch between original and
refactored. Verification instead used a **0-96h window** (`end=96,
delta=24`, matching the original's own `delta=24` daily-observation
convention), safely inside every one of the 5 scenarios' own 112-298h
solver-success range (112h is the earliest boundary, Placebo).

**Dosing.** All 5 of the original file's own scenarios (from
`mk_pirf_ev()`/`mk_nint_ev()`/`mk_combo_ev()`/`scenarios`), unchanged
amounts/timing, truncated to the 96h window via `addl`:

| Scenario | Dosing |
|---|---|
| Placebo (Natural History) | no dosing |
| Pirfenidone 801 mg TID | 801,000 µg q8h, cmt 1 (`GUT_PIR`) |
| Nintedanib 150 mg BID | 150,000,000 ng q12h, cmt 4 (`GUT_NIN`) |
| Combination (Pirf + Nint) | both of the above together |
| Pirfenidone 267 mg TID (low dose) | 267,000 µg q8h, cmt 1 (`GUT_PIR`) |

**Result: exact match (max abs/rel deviation 0.0) for all 5 scenarios**,
across all 27 shared outputs (17 `$CMT` states common to both files —
disease-side compartments were never renamed — plus 10 `$TABLE`/`$CAPTURE`
outputs whose names are unchanged: `Cp_pirf`, `Cn_nint`, `FVC_pct`,
`DLCO_pct`, `FVC_decline_yr`, `Col_norm`, `AEC2_norm`, `TGFb_norm`,
`MMP_norm`, `TIMP_norm`, `MMP_TIMP_ratio`, `Myofib_norm`, `ROS_norm`,
`Periostin_proxy`, `MMP7_proxy`, `KL6_proxy`), 5-7 timepoints per scenario
(0-96h at 24h resolution, plus mrgsolve's own duplicate t=0 row):

| Scenario | max abs diff | max rel diff |
|---|---|---|
| Placebo (Natural History) | 0.0 | 0.0 |
| Pirfenidone 801 mg TID | 0.0 | 0.0 |
| Nintedanib 150 mg BID | 0.0 | 0.0 |
| Combination (Pirf + Nint) | 0.0 | 0.0 |
| Pirfenidone 267 mg TID (low dose) | 0.0 | 0.0 |

Consistent with the expectation for Archetype 3 with pure Hill-renames (no
fitting): both compounds' `inh_*` terms were already exact Hill ratios in
the original, so renaming, reorganizing into the `<ROLE>_<STEM>`
convention, and moving `$INIT` into `$MAIN` introduces no numeric
difference — this result is also the strongest available evidence that
keeping `C_PIR`/`C_NIN`/`EFFECT_*` in `$MAIN` (rather than `$ODE`) was the
correct call: had the lag not been reproduced faithfully, this match would
not have been exact (as the development-time `$ODE` draft demonstrated
directly, ~35% deviation on `TGFb_norm`).

No fitting was needed anywhere in this file (no TMDD, no ODE-solved
receptor kinetics) — the guide's "fit an approximation" path for
Archetype 4 does not apply here for either compound.

## Manifest check

`POST /model_manifest` on the refactored DSL confirms every
naming-convention parameter is discoverable in `parameters`: `KA_PIR,
F_PIR, CL_PIR, V1_PIR, V2_PIR, Q_PIR, KA_NIN, F_NIN, CL_NIN, V1_NIN,
V2_NIN, Q_NIN, EC50_PIR_TGFB, EMAX_PIR_TGFB, GAMMA_PIR_TGFB, EC50_PIR_M2,
GAMMA_PIR_M2, EMAX_PIR_ROS, GAMMA_PIR_ROS, EC50_NIN, EMAX_NIN, GAMMA_NIN,
MW_PIR, MW_NIN` (50 parameters total, including unchanged disease-side
ones). `C_PIR`, `EFFECT_PIR_TGFB`, `EFFECT_PIR_M2`, `EFFECT_PIR_ROS`,
`C_NIN`, `EFFECT_NIN` all appear in `outputPaths` (39 outputs total) —
computed in `$MAIN` as plain `double` locals (not `$PARAM`, per the
mrgsolve 2.0.1 read-only-reference constraint documented in
`pagets-disease/pbd_refactor_notes.md` and every refactor since) and
exposed via `$CAPTURE`.

## Census update

`driver-patches/data/compound_perturbation_census.md` rows for
`idiopathic-pulmonary-fibrosis` `N`/`P` updated with real compound names
(Nintedanib, Pirfenidone), the corrected Target/Pathway for `P` (was
`PCSK9`, now TGF-β1/anti-fibrotic-anti-inflammatory), compartment/output
names, and this outcome summary.
