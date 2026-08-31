# Refactor notes — `alpha1-antitrypsin-deficiency/aatd_mrgsolve_model.R`

**Scope of this refactor**: only the two compounds listed in
`driver-patches/data/compound_perturbation_census.md` for this file — AUG
(AAT augmentation therapy) and NEI (NE inhibitor, Alvelestat) — were
touched. Every other compound (Fazirsiran siRNA, gene therapy) and the
entire disease-mechanism network (hepatic Z-AAT polymer/fibrosis cascade,
pulmonary neutrophil/elastase/elastin/FEV1 cascade) is copied verbatim
from the original into `aatd_mrgsolve_model_refactored.R`.

## The original does not compile under mrgsolve 2.0.1 — logged, not fixed upstream

`POST /model_manifest` on the untouched original's own DSL (extracted
verbatim from `code <- '...'`) fails twice in sequence, both from the
same class of defect (a dead line using the removed `self` object):

1. `double dose_aug = self.mtime(0) ; // Handled via event table` —
   `self` undeclared in this scope. `dose_aug` is never read anywhere
   else in the file.
2. Once (1) is removed: `double FEV1_nat = FEV1_0 * (1.0 - FEV1_loss) -
   annual_decline_base * TIME / 365.0 ;` — mrgsolve's `TIME` macro
   expands to `self.time`, equally undeclared. `FEV1_nat` is likewise
   never read anywhere else; the model's real FEV1 trajectory comes
   entirely from `dxdt_FEV1_pct`, which never references `TIME`.

Both are logged as **`translations/UPSTREAM_ISSUES.md` #109**, not fixed
in the checked-in original, per the never-edit-upstream rule.

**Fix applied directly to the delivered `aatd_mrgsolve_model_refactored.R`**
(per the guide's settled policy for a non-compiling original — syntax-only,
non-numeric): both dead lines are deleted outright rather than patched,
since neither is read anywhere else in the file — deleting a provably-dead
assignment cannot change any numeric output. This is confirmed, not just
asserted: the verification below shows every `$CAPTURE`d output matching
exactly (max abs diff 0.0) between the syntax-fixed original and the
refactored model across three scenarios.

## A second, unrelated defect found in AUG's own block — also logged, not repaired

While isolating AUG's PK block, the augmentation-therapy pathway turned
out to be **completely disconnected from the disease state it is meant to
treat**. The original computes:

```
// Augmentation adds external M-AAT (already in AUG_C1 mg/dL)
double AAT_total_input = ZAAT_endogenous + Gene_MAAT + AUG_C1 ;

dxdt_AAT_C1 = ZAAT_endogenous + Gene_MAAT - k_AAT_elim * AAT_C1
              - k12_aug * AAT_C1 + k21_aug * AAT_C2 ;
```

`AAT_total_input` names the intent plainly ("Augmentation adds external
M-AAT") but is **never read again anywhere in the file** — the actual
`dxdt_AAT_C1` equation that drives every downstream output has no
`AUG_C1`/`AAT_total_input` term at all. `AUG_C1`/`AUG_C2` do have their
own working two-compartment linear PK and integrate correctly under
dosing (confirmed below, `C_AUG` rises and falls exactly as the dosing
schedule dictates) — the drug is delivered and cleared as intended, it
simply never touches the AAT/ELF/NE pathway. **Practical effect: dosing
the model's own "S2: AAT Augmentation" scenario changes zero captured
outputs relative to the untreated scenario**, confirmed empirically
below, not just by code inspection.

Logged as **`translations/UPSTREAM_ISSUES.md` #110**. Per the
never-edit-upstream rule and this fork's "log what you find, don't fix
it" principle, **the refactored model preserves this disconnection
bug-for-bug** — `C_AUG` is exposed as the redirect point per the naming
convention (so a future patch has somewhere to plug an external
covariate into), but no `EFFECT_AUG`/`EMAX_AUG`/`EC50_AUG` was invented
to paper over the gap, since the original computes no such term at all.
This is the same disclosed-disconnection pattern as
`vexas-syndrome/vexas_mrgsolve_model.R`'s azacitidine (UPSTREAM_ISSUES.md
#92) and `denosumab` in `hypercalcemia-of-malignancy` (structural
coupling preserved, not invented).

## Archetypes

### AUG (augmentation therapy) — Archetype 2 (two-compartment, no depot, linear), bespoke: micro-rate constants kept, no volume conversion

Census classification: "Delete PK compartment; concentration is itself
the state" — confirmed correct. `AUG_C1` (mg/dL) is dosed and read
directly as a concentration; `V1_aug` is declared in `$PARAM` but **never
referenced in any `$ODE` line** (checked by grep across the whole file).
Same "compartment already IS the concentration" pattern as `TOCI_C` in
`sepsis/sep_refactor_notes.md` and `CENT_HEM` in
`acute-intermittent-porphyria/aip_refactor_notes.md`.

| Original | Refactored | Role |
|---|---|---|
| `AUG_C1` (cmt) | `CENT_AUG` | central, itself a concentration (mg/dL), dosed directly |
| `AUG_C2` (cmt) | `PERI_AUG` | peripheral |
| `k10_aug` | `CL_AUG` | elimination rate constant, applied directly (no volume) |
| `k12_aug` | `K12_AUG` | central→peripheral micro-rate constant |
| `k21_aug` | `K21_AUG` | peripheral→central micro-rate constant |
| `V1_aug` | `V1_AUG` | **disclosed-dead**: unused by any `$ODE` line, original or refactored |
| `F_aug` | `F_AUG` | **disclosed-dead**: unused (no depot; IV dosed directly into `CENT_AUG`) |
| `ka_aug` | `KA_AUG` | **disclosed-dead**: unused (no depot) |

**Why `K12_AUG`/`K21_AUG` stay micro-rate constants rather than becoming
`Q_AUG`/`V1_AUG`/`V2_AUG`-style CL/Q/V** (the guide's default for
Archetype 2): that conversion presumes a real volume of distribution to
convert with, and this compound has none — `V1_aug` is dead, not merely
unused-but-present-in-the-math. Same disclosed deviation as `K12_ZOL`/
`K21_ZOL` in `hypercalcemia-of-malignancy/mah_refactor_notes.md`'s
zoledronate block, for the identical reason.

**A genuine, disclosed cross-coupling kept exactly as the original had
it**: `K12_AUG`/`K21_AUG` (AUG's own inter-compartmental transfer rates)
are *also* reused, by parameter value only, inside the disease-state
equations `dxdt_AAT_C1`/`dxdt_AAT_C2` (endogenous serum-AAT
central↔peripheral split) — a coincidental sharing of numeric parameter
values between two structurally unconnected compartment pairs in the
original, not a hidden functional coupling (`AUG_C1`/`AUG_C2` dosing does
not change `AAT_C1`/`AAT_C2`'s trajectory at all — see verification
below). Preserved via the renamed parameter at both use sites rather than
decoupled into two independent parameters, since decoupling would add a
parameter the original never had. `AAT_C1`/`AAT_C2` themselves are
disease state, not part of AUG's own block, and are otherwise untouched.

**No `EFFECT_AUG`**: the original computes no scalar effect term for
augmentation therapy at all (see the disconnection defect above) — same
"don't force a bad fit" / "don't invent an effect the original doesn't
have" precedent as denosumab in `hypercalcemia-of-malignancy` and
azacitidine in `vexas-syndrome`.

### NEI (NE inhibitor, Alvelestat) — Archetype 3 minus peripheral, bespoke identity concentration

Census classification: "Redirect concentration (clean single site)" —
confirmed correct. `NEi_C` has an absorption depot (`NEi_A`) but no
peripheral compartment and no volume parameter anywhere in the original
(`EC50_NEi` is itself denominated "µg/mL normalized", consistent with
`NEi_C` already being a concentration by construction). The original's
effect term is already a plain Emax ratio
(`Emax_NEi * NEi_C / (EC50_NEi + NEi_C)`) — a rename, not a refit.

| Original | Refactored | Role |
|---|---|---|
| `NEi_A` (cmt) | `GUT_NEI` | absorption depot |
| `NEi_C` (cmt) | `CENT_NEI` | central, itself a concentration, no volume in the original |
| `ka_NEi` | `KA_NEI` | absorption rate |
| `k_NEi_elim` | `CL_NEI` | elimination rate constant, applied directly (no volume) |
| `Emax_NEi` | `EMAX_NEI` | Hill Emax |
| `EC50_NEi` | `EC50_NEI` | Hill EC50 |
| — (none) | `GAMMA_NEI` (new) | Hill coefficient [math-implied =1; original had no explicit exponent] |

No `F_NEI` — the original has no bioavailability parameter for NEi at
all (implicitly F=1), so none is invented, same as `CTN`/`PRD`/`CIN`/
`FUR` in `hypercalcemia-of-malignancy/mah_refactor_notes.md`.

`EFFECT_NEI` replaces `NEi_effect` at its one use site (`NE_inhib_total`'s
`(1.0 - AAT_inhib_eff) * EFFECT_NEI` term).

## `$PARAM` vs `$CAPTURE` for `C_<STEM>`/`EFFECT_<STEM>`

Not attempted as `$PARAM = 0` declarations, per the mrgsolve 2.0.1
read-only-reference constraint already documented in
`pagets-disease/pbd_refactor_notes.md` and every refactor since (e.g.
`hypercalcemia-of-malignancy/mah_refactor_notes.md`). Followed the same
precedent directly: `C_AUG`, `C_NEI`, `EFFECT_NEI` are plain `double`
locals computed in `$ODE`, listed in `$CAPTURE`. Confirmed via `POST
/model_manifest` on the refactored DSL: all three appear in
`outputPaths`, and every renamed/new parameter (`CL_AUG`, `K12_AUG`,
`K21_AUG`, `V1_AUG`, `F_AUG`, `KA_AUG`, `KA_NEI`, `CL_NEI`, `EMAX_NEI`,
`EC50_NEI`, `GAMMA_NEI`) appears in `parameters` with its original
numeric default.

## Verification

**Method.** Extracted the bare DSL text from `aatd_mrgsolve_model.R`
(`code <- '...'`, with the two dead lines from UPSTREAM_ISSUES.md #109
removed — the minimal, symmetric, non-committed change needed to compile
it at all) and from `aatd_mrgsolve_model_refactored.R` (same variable,
refactored), and ran them through the local qspserver `mrgsolve_api`
(`http://localhost:8007`), `POST /model_manifest` then `POST
/run_simulation`, requests spaced ~2.2–2.5 s apart per the shared-service
note. The API's `dosing` field takes a literal numeric `rate` (0 = bolus)
per compartment number, not mrgsolve's own `rate = -2`/`D_<cmt>`-duration
convention that the original's own `aug_dose()` R helper uses (and for
which the original defines no `D_AUG_C1` parameter — a separate,
R-driver-only quirk, moot in practice since AUG has no effect on any
output regardless per UPSTREAM_ISSUES.md #110). Reproduced each
scenario's dosing intent as an equivalent bolus/repeated-bolus schedule
via the API's own `dosing` schema, applied identically to both models —
a verification-harness substitution, not a change to either model.

**Scenarios run** (compartment numbers 1-based, unchanged in position
between original and refactored since only names changed):

1. **AUG alone** — reproducing the file's own `S2` scenario (`aug_dose()`:
   55 mg/dL weekly IV): `dosing = [{amt:55, cmt:14 (AUG_C1/CENT_AUG),
   rate:0, ii:7, addl:25}]`, `time = {end:180, delta:1}` (180 days
   instead of the scenario's full 1825-day/5-year horizon — shortened per
   the guide's solver-step-budget note; 180 days is ample to demonstrate
   exact matching given AUG's own elimination/redistribution rates are
   all on a several-day-to-weeks time scale).
2. **NEI alone** — reproducing the file's own `S4` scenario (`nei_dose()`:
   0.30 units BID oral): `dosing = [{amt:0.30, cmt:16 (NEi_A/GUT_NEI),
   rate:0, ii:0.5, addl:59}]`, `time = {end:30, delta:0.5}` (30 days;
   NEI's own PK/PD equilibrates on an hours-to-day timescale, so 30 days
   is ample).
3. **AUG + NEI together** — reproducing the file's own `S6` scenario
   (both dosing schedules combined, to rule out any cross-block
   interference from the two independent renames): both `dosing` records
   above (AUG for 12 weeks, NEI BID) over `time = {end:84, delta:1}`.

**Result: exact match, max abs diff = 0.0, every one of the 14 shared
`$CAPTURE`d outputs (`AAT_mgdL, AAT_uMol, ELF_uMol, EmphIndex, FEV1,
SGRQ_score, Exacer_yr, GOLD_stage, Metavir_fib, ZPolymer, NEfree,
NE_inhib, HSC, LiverFib`), in all three scenarios** — compared
point-by-point across each scenario's full time grid (62–182 points
depending on scenario). This is the expected outcome for two pure
structural reorganizations (both AUG and NEI are renames, not fits —
`GAMMA_NEI = 1.0` reproduces the original's implicit-exponent-1 ratio
exactly).

**AUG's own PK integrates correctly** (confirmed independently of the
disease-state comparison above, since AUG has no disease-state effect to
compare against): `C_AUG` in the refactored model's own output rises and
falls in step with the weekly dosing schedule (e.g. samples `[0, 0,
50.59, 46.56, 42.89, ...]` mg/dL over the first several days post-dose,
each subsequent weekly dose producing the same rise-then-decay pattern) —
proof that `CENT_AUG`/`PERI_AUG`'s own two-compartment PK is intact and
correctly renamed, and that the zero effect on disease outputs is
specifically the pre-existing disconnection (UPSTREAM_ISSUES.md #110),
not a bug introduced by this refactor breaking AUG's PK itself.

**NEI's effect is correctly wired**: in the NEI-alone scenario, `NE_inhib`
(captured NE-inhibition percentage) drops from its untreated baseline
(0.9411) down through 0.891 → 0.845 → 0.804 → 0.765 as `EFFECT_NEI` rises
with accumulating drug exposure — matching exactly (0.0 max abs diff)
between the original (via its inline `NEi_effect` ratio) and the
refactored model (via the renamed, `pow()`/`GAMMA_NEI`-based
`EFFECT_NEI`), confirming the rename changed nothing about NEI's actual
pharmacology.

`/model_manifest` on the refactored DSL additionally confirmed all eleven
renamed/new parameters and seven new/renamed compartments/covariates
(`CENT_AUG, PERI_AUG, GUT_NEI, CENT_NEI, C_AUG, C_NEI, EFFECT_NEI`) appear
in `outputPaths`/`parameters` as required by the guide's qspserver
compatibility section.

No scratch/debug artifacts were left in the repository — DSL extraction,
the fixed-original copy, and the three-scenario comparison harness ran
entirely from a session-local scratchpad directory outside the repo
tree, deleted after use.

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`: the
`AUG` row confirmed "Delete PK compartment; concentration is itself the
state" (bespoke: micro-rate constants kept, disclosed pre-existing
disconnection from disease state — UPSTREAM_ISSUES.md #110); the `NEI`
row confirmed "Redirect concentration (clean single site)" (clean rename,
exact verification match).
