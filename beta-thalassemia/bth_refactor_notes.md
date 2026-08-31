# Refactor notes — `bth_mrgsolve_model_refactored.R`

Scope: **all three of the file's PK compounds** — Luspatercept, Deferasirox
(DFX), Hydroxyurea (HU) — per the `beta-thalassemia | DFX`,
`beta-thalassemia | HU`, `beta-thalassemia | L` rows in
`driver-patches/data/compound_perturbation_census.md`, all three classified
"Normalize duplicate concentration sites, then redirect." The disease model
itself (erythropoiesis cascade, EPO feedback, ERFE/hepcidin axis, iron
metabolism, transfusion dosing) is untouched — no compartment added, removed,
or renamed outside the three compounds' own PK/PD blocks.

## What "L" turned out to be

The census row abbreviates the third compound as the bare letter `L`. The
file's own comments (`Luspatercept PK (2-cpt SC)`, `Ref: Platzbecker et al.
NEJM 2017 (MDS), Viprakasit NEJM 2020`) and its dosing/calibration references
(BELIEVE and BEYOND trials, both luspatercept trials in transfusion-dependent
and non-transfusion-dependent beta-thalassemia) confirm it is
**luspatercept**. Per the task instructions, a clearer stem than the bare
`L` was used throughout: **`LUSP`**.

## What the duplicate-definition problem actually was

Each compound's plasma concentration was computed **twice**, in two
independently-typed sites that had no relationship to each other in the
source text:

1. `$ODE`, cached once per compound: `double C_LUSPAT = LUSPAT_C1/Vc_L;`,
   `double C_DFX = DFX_CENT/V_DFX;`, `double C_HU = HU_CENT/V_HU;` — these
   fed the compounds' own disease-effect terms (`IE_effect_L`, `DFX_Fe_remov`,
   `HU_HbF`).
2. `$TABLE`, re-derived completely independently, under different names,
   purely for reporting: `double C_Luspa = LUSPAT_C1/Vc_L;`,
   `double C_Deferasirox = DFX_CENT/V_DFX;`, `double C_HU_calc = HU_CENT/V_HU;`.

Luspatercept's disease **effect** (not just its concentration) was
*additionally* re-derived a second, independent time in `$TABLE`:
`double IE_frac_eff = ie_frac * (1.0 - Emax_L*C_Luspa/(EC50_L+C_Luspa));`,
separate from the `$ODE`-computed `IE_effect_L`/`ie_eff` that the
erythropoiesis ODEs actually consume.

## Fix: `$GLOBAL` macros, not two independently-typed doubles

Rather than accepting "two names for one quantity, kept manually in sync"
(the pattern several earlier refactors in this fork settled for when the
duplicate spans the `$ODE`/`$TABLE` block boundary — see
`alkaptonuria/aku_refactor_notes.md` and
`benign-prostatic-hyperplasia/bph_refactor_notes.md`, both of which found
that mrgsolve 2.0.1 hard-errors on literally reusing the same variable name
in both blocks: `redefinition of 'double {anonymous}::C_NT'`), this file
follows the newer, stronger fix documented in FORK_WORKFLOW_GUIDE.md ("The
dose-instant reporting artifact — there is a real fix, not just a
disclosure"), which `clostridioides-difficile-infection`'s refactor
established: declare the concentration/effect quantities as **`$GLOBAL`
preprocessor macros** instead of `$ODE`-local (or `$TABLE`-local) `double`s.

```c
$GLOBAL
#define C_LUSP  (CENT_LUSP / V1_LUSP)
#define C_DFX   (CENT_DFX  / V1_DFX)
#define C_HU    (CENT_HU   / V1_HU)
#define EFFECT_LUSP (EMAX_LUSP * pow(C_LUSP, GAMMA_LUSP) / (pow(EC50_LUSP, GAMMA_LUSP) + pow(C_LUSP, GAMMA_LUSP)))
#define EFFECT_DFX  (KFE_DFX * C_DFX)
#define EFFECT_HU   (EMAX_HU * pow(C_HU, GAMMA_HU) / (pow(EC50_HU, GAMMA_HU) + pow(C_HU, GAMMA_HU)))
```

Because a `#define` is textual substitution, not a typed variable, the
*identical* definition is used at every read site in **both** `$ODE` and
`$TABLE` — this is genuine single-source-of-truth normalization (one
definition, many use sites), not two hand-written re-derivations that merely
happen to agree today. It also sidesteps the block-boundary redeclaration
error the `aku`/`bph` refactors hit, without needing their `CNTo`/`_tbl`-style
distinct-name workaround. As a side effect it also eliminates the
dose-instant stale-value reporting artifact described in
FORK_WORKFLOW_GUIDE.md by construction (a macro has no per-block-evaluation
cache to go stale) — confirmed empirically: the paired `/run_simulation`
comparisons below show `max abs diff 0.0` at every time point, including the
duplicate dose-instant row the API emits, for every scenario tested.

`$TABLE`'s old `C_Luspa`/`C_Deferasirox`/`C_HU_calc` local doubles are gone
entirely — superseded by directly `$CAPTURE`-ing the canonical
`C_LUSP`/`C_DFX`/`C_HU` macros (confirmed this works: mrgsolve's
`@`-list-style `$CAPTURE` accepts a bare macro name and expands it in place,
same as observed for `clostridioides-difficile-infection`'s `C_VAN`/
`EFFECT_VAN` macros). `IE_frac_eff` is **kept** as its own distinct
`$TABLE`-computed reporting quantity — it is not simply a duplicate of
`EFFECT_LUSP`, it *applies* `EFFECT_LUSP` to `ie_frac` — but its RHS now
reads the canonical `EFFECT_LUSP` macro instead of re-deriving the ratio
inline.

**Important, preserved-not-fixed caveat on `IE_frac_eff`:** in the
*original*, this `$TABLE` quantity was already a **luspatercept-only partial
recompute** — it never included hydroxyurea's additional reduction that
`$ODE`'s `ie_eff2` applies (`ie_eff2 = ie_eff * (1.0 - 0.3*HU_HbF)`).
`IE_frac_eff` in both the original and the refactored file therefore
diverges from the true simulated `ie_eff2` whenever HU is dosed alongside
luspatercept — this is a pre-existing quirk of the original's own reporting
code, reproduced exactly (not "fixed") by this refactor, since fixing it
would be a behavioral change beyond a rename/reorganization.

## Naming convention applied

### Luspatercept (archetype 3: depot + central + peripheral, linear)

| Original | Refactored | Original | Refactored |
|---|---|---|---|
| `LUSPAT_SC` (cmt) | `GUT_LUSP` | `Ka_L` | `KA_LUSP` |
| `LUSPAT_C1` (cmt) | `CENT_LUSP` | `F_L` | `F_LUSP` |
| `LUSPAT_C2` (cmt) | `PERI_LUSP` | `CL_L` | `CL_LUSP` |
| `Vc_L` | `V1_LUSP` | `Vp_L` | `V2_LUSP` |
| `Q_L` | `Q_LUSP` | `EC50_L` | `EC50_LUSP` |
| `Emax_L` | `EMAX_LUSP` | — | `GAMMA_LUSP` (new, =1.0) |
| `C_LUSPAT` (`$ODE` local) / `C_Luspa` (`$TABLE` local) | `C_LUSP` (`$GLOBAL` macro) | `IE_effect_L` | `EFFECT_LUSP` (`$GLOBAL` macro) |

### Deferasirox (archetype 3's no-peripheral variant: depot + central, linear)

| Original | Refactored | Original | Refactored |
|---|---|---|---|
| `DFX_GUT` (cmt) | `GUT_DFX` | `Ka_DFX` | `KA_DFX` |
| `DFX_CENT` (cmt) | `CENT_DFX` | `V_DFX` | `V1_DFX` |
| `k_DFX_Fe` | `KFE_DFX` | `C_DFX` (`$ODE` local) / `C_Deferasirox` (`$TABLE` local) | `C_DFX` (`$GLOBAL` macro) |
| `DFX_Fe_remov` | `EFFECT_DFX` (`$GLOBAL` macro) | | |

(`F_DFX`, `CL_DFX` kept their original names — already convention-compliant.)

### Hydroxyurea (archetype 3's no-peripheral variant: depot + central, linear)

| Original | Refactored | Original | Refactored |
|---|---|---|---|
| `HU_GUT` (cmt) | `GUT_HU` | `Ka_HU` | `KA_HU` |
| `HU_CENT` (cmt) | `CENT_HU` | `V_HU` | `V1_HU` |
| `Emax_HU` | `EMAX_HU` | — | `GAMMA_HU` (new, =1.0) |
| `C_HU` (`$ODE` local) / `C_HU_calc` (`$TABLE` local) | `C_HU` (`$GLOBAL` macro) | `HU_HbF` | `EFFECT_HU` (`$GLOBAL` macro) |

(`F_HU`, `CL_HU`, `EC50_HU` kept their original names — already
convention-compliant.)

All PK is a pure rename — no value changed, no compartment added, removed,
or reordered (compartment declaration order in `$INIT`, which fixes
`/run_simulation`'s 1-based `cmt` dosing index, is identical between the two
files; confirmed via the `outputPaths` comparison below).

## The Hill interface

**Luspatercept and hydroxyurea** were already bare `Emax*C/(EC50+C)` ratios
(`IE_effect_L`, `HU_HbF`) — promoting them to the guide's explicit Hill form
(`EMAX*C^GAMMA/(EC50^GAMMA+C^GAMMA)`) is a **rename, not a refit**:
`GAMMA_LUSP`/`GAMMA_HU` are new named parameters set to `1.0`, and
`pow(x,1)=x`, so the numeric value is unchanged (confirmed to floating-point
JSON-rounding precision, see the standalone macro-consistency check below).

**Deferasirox's** chelation effect (`k_DFX_Fe*C_DFX`, renamed
`EFFECT_DFX = KFE_DFX*C_DFX`) is **deliberately kept linear/unbounded**, not
forced into the saturating Hill template — the original's own mechanism has
no saturation term (it is a first-order removal-rate constant times
concentration, not a receptor-occupancy ratio), so imposing an `EMAX`/`EC50`
form would silently cap a mechanism the original leaves uncapped. Same
reasoning FORK_WORKFLOW_GUIDE.md documents for
`alkaptonuria/aku_refactor_notes.md`'s nitisinone term.

Each compound's `EFFECT_<STEM>` is combined into the disease equations
independently, exactly where the original did so (`ie_eff` reads
`EFFECT_LUSP` alone; `ie_eff2` then additionally applies `EFFECT_HU`;
`dxdt_FE_LIV` applies `EFFECT_DFX` alone) — no shared/collapsed Hill term
across compounds.

## Pre-existing build defects (unrelated to any compound's own PK) — logged as `translations/UPSTREAM_ISSUES.md` #131

The untouched original does not compile under this repo's mrgsolve 2.0.1
build, for three independent reasons, none of which touch any compound's own
PK/PD block:

1. `$PLUGIN autodiff nm-vars` — the `autodiff` plugin is not available
   (`Error: plugin autodiff could not be found`).
2. `$CAPTURE` lists ten names that are already `$INIT`-declared (i.e.
   implicit `$CMT`) compartments (`BFU_E CFU_E PRO_E BASO_E POLY_E ORTHO_E
   RETIC RBC_MAT FE_PL FE_CARD`) — mrgsolve 2.0.1 rejects this at
   model-object construction (`compartment should not be in $CAPTURE`).
3. `$ODE`'s cardiac-iron block calls `pmax(0.0, ...)` — `pmax` is R's
   vectorised max, not a C++ function available inside `$ODE`'s compiled
   scope (`'pmax' was not declared in this scope; did you mean 'fmax'?`).

Per FORK_WORKFLOW_GUIDE.md's "When the original doesn't compile at all"
policy: **not fixed in the checked-in original** (never-edit-upstream rule).
`bth_mrgsolve_model_refactored.R` carries the syntax-only fix forward
directly, disclosed inline at each site: (1) `$PLUGIN autodiff nm-vars` →
`$PLUGIN nm-vars`; (2) the ten compartment names dropped from `$CAPTURE`
(all ten remain reported, unchanged, as compartments — mrgsolve's
`outputPaths` concatenates compartments and captures, confirmed below); (3)
`pmax(...)` → `fmax(...)`, algebraically identical for two scalar arguments.
None of the three changes any numeric value or any compound's own PK/PD.
Full detail in `translations/UPSTREAM_ISSUES.md` #131.

## qspserver compatibility

Confirmed via `POST /model_manifest` against the qspserver `mrgsolve_api`
container (`http://localhost:8007`):

- The pristine original does **not** build (see the three defects above);
  the syntax-patched original (defects 1–3 fixed, no other change) builds
  cleanly, `outputPaths` count 34.
- The refactored DSL builds cleanly, `outputPaths` count 37 (three more than
  the patched original: `EFFECT_LUSP`, `EFFECT_DFX`, `EFFECT_HU`, the newly
  discoverable canonical effect terms — `IE_frac_eff` already existed in
  both).
- Every renamed compartment (`GUT_LUSP CENT_LUSP PERI_LUSP GUT_DFX CENT_DFX
  GUT_HU CENT_HU`) appears in `outputPaths` in the same 1-based position its
  original name held (`LUSPAT_SC LUSPAT_C1 LUSPAT_C2 DFX_GUT DFX_CENT HU_GUT
  HU_CENT`), so dosing into `cmt = 16/19/21` (the luspatercept/deferasirox/
  hydroxyurea depot compartments) is unaffected by the rename.
- All renamed/new `$PARAM` entries confirmed present:
  `KA_LUSP, F_LUSP, CL_LUSP, V1_LUSP, V2_LUSP, Q_LUSP, EC50_LUSP, EMAX_LUSP,
  GAMMA_LUSP, KA_DFX, F_DFX, CL_DFX, V1_DFX, KFE_DFX, KA_HU, F_HU, CL_HU,
  V1_HU, EC50_HU, EMAX_HU, GAMMA_HU` (61 total `$PARAM` entries).
- `C_LUSP, C_DFX, C_HU, EFFECT_LUSP, EFFECT_DFX, EFFECT_HU` are `$GLOBAL`
  macros, not `$PARAM` entries (a `#define`d name cannot also be a settable
  `$PARAM` default — same constraint documented in several other
  `*_refactor_notes.md` in this repo). They are discoverable via
  `outputPaths` instead (confirmed present).

## Verification

**Method.** Both DSL blocks were extracted verbatim from `code <- '...'`,
syntax-patched per the three build-compat fixes above (applied identically
to a scratch copy of the original, and as the delivered fix in the
refactored file), and POSTed to the qspserver `mrgsolve_api` container's
`/run_simulation`, spaced ~2.5 s apart. `outputs` requested the compound's
own-named columns from each model (e.g. `C_Luspa` from the original,
`C_LUSP` from the refactored file) plus every shared, unrenamed disease
output (`Hb_gdL, EPO_IUL, ERFE_rel, HEPC_nM, LIC_mgFe, FERR_ugL, CARD_T2star,
Retic_pct, IE_frac_eff`, and the eleven unrenamed compartments).

Dosing was built directly as `DoseSpec` records reproducing the magnitudes
the original file's own `ev()` objects define (`luspatercept_ev`,
`dfx_dose`, `hu_dose`, `tx_iron_ev`), using each object's own `amt`/`cmt`/
`ii`/`addl`. Five checks:

| # | Scenario | Dosing (from the original's own `ev()` objects) | Window | Params |
|---|---|---|---|---|
| A | Luspatercept only | `luspatercept_ev` (70000 µg → `GUT_LUSP`, q21d ×12) | 252 d (= scenario 4's own window) | scenario 4's own `ie_frac=0.55, EPO_base=30` |
| B | Deferasirox only | `dfx_dose` (2100 mg → `GUT_DFX`, qd) | 90 d (shortened from scenario 3's 365 d for solver budget) | scenario 3's own `ie_frac=0.80, EPO_base=50` |
| C | Hydroxyurea only | `hu_dose` (1400 mg → `GUT_HU`, qd) | 90 d | default `ie_frac=0.70` |
| D | Transfusion + Deferasirox, concurrent | `tx_iron_ev` (0.30 → `FE_LIV`, q21d) + `dfx_dose` | 90 d | scenario 3's own params |
| E | Luspatercept + Transfusion + Deferasirox, concurrent | all three, each at its own `ii`/`addl` | 252 d | scenario 5's own params |

**A note on scenario reconstruction (D/E vs the original's `ev_seq()`).**
The original file combines multiple `ev()` objects for scenarios 3 and 5 via
`ev_seq()`, whose exact time-offsetting semantics were not inspected here
(deliberately avoided invoking local R/mrgsolve for anything beyond
reading/writing files, per the task's constraint). Tests D and E instead
dose every compound concurrently starting at `time=0`, each at its own
`ii`/`addl` from the original's own event objects — a reasonable
reconstruction (concurrent transfusion+chelation, and luspatercept+
transfusion+chelation, are the clinically realistic reading of scenarios 3
and 5), but not a guaranteed byte-for-byte reproduction of `ev_seq()`'s
internal timing. This does not weaken the verification itself: both the
original and refactored DSL receive the *identical* reconstructed dosing in
each test, so any divergence between them still isolates a genuine
regression regardless of whether the reconstructed timing exactly matches
what `ev_seq()` would produce.

**`hu_dose` is defined in the file but never used by any of the six named
scenarios** (`sc1`...`sc6`) — grepped for `hu_dose` and confirmed its only
occurrence is its own definition. Test C runs it anyway, since it *is* one
of "the original file's own R code" objects, just not wired into a named
scenario — worth flagging as its own small finding, not logged as an
UPSTREAM_ISSUES.md defect (an unused-but-harmless R object is not a build or
correctness defect).

**Result: exact match, every time.** All five paired runs show **max abs
diff = 0.0, max rel diff = 0.0** across every shared/mapped output, at every
time point (including the API's duplicate dose-instant row) — pure
structural reorganization plus a rename-only Hill promotion (`GAMMA=1`), as
expected for archetype 3 and its no-peripheral variant.

**Additional macro-consistency check.** Recomputed `EFFECT_DFX`,
`EFFECT_LUSP`, `EFFECT_HU` from the refactored model's own returned
`C_DFX`/`C_LUSP`/`C_HU` and the known parameter values
(`KFE_DFX=0.12; EMAX_LUSP=0.65,EC50_LUSP=0.5; EMAX_HU=0.40,EC50_HU=10.0`) in
a 30-day, all-three-compounds-dosed run: residuals ≤ 4.9e-5 against the
API's own returned `EFFECT_*` values — consistent with the API's 4-decimal
JSON rounding on both the input (`C_*`) and output (`EFFECT_*`) sides, not a
calculation discrepancy (confirms the `$GLOBAL` macro formulas are wired up
correctly, independent of the exact-match original-vs-refactored comparison
above).

## Anything else worth flagging

- `WT = 70.0` (body weight) was left untouched — it is a shared, non-compound
  -specific parameter used only by the R-side dose-amount calculations for
  all three compounds, not part of any one compound's own PK block.
- The `IE_frac_eff` pre-existing partial-recompute caveat (see above) is
  worth a reviewer's attention if this model is ever used to validate a
  combined luspatercept+hydroxyurea regimen against `IE_frac_eff` specifically
  — the true simulated effective IE fraction is `ie_eff2`, not
  `IE_frac_eff`.
- `HGB_per_RBC = 28e-6` reports as `0` in `/model_manifest`'s parameter list
  for both the original and refactored model (a manifest-introspection
  rounding/formatting quirk on very small scientific-notation defaults, not
  a build or simulation defect — the parameter's actual runtime value is
  unaffected, confirmed by both models producing identical `Hb_gdL` output).
  Not logged as a new UPSTREAM_ISSUES.md entry: it does not affect this
  file's own use of the parameter (`HGB_per_RBC` is declared but never
  actually read anywhere in `$ODE`/`$TABLE` — `Hb_gdL` uses the separate
  simplified `RBC_MAT/RBC_ref*14.0` formula instead), and reproduces
  identically pre- and post-refactor.

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`, the
`beta-thalassemia | DFX`, `beta-thalassemia | HU`, `beta-thalassemia | L`
rows.
