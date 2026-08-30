# Refactor notes — `acute-intermittent-porphyria/aip_mrgsolve_model.R`

Two compounds in this file, both refactored independently per the fork
guide's classification (`driver-patches/data/compound_perturbation_census.md`,
both rows tagged "Redirect concentration (clean single site)"): **givosiran
(GIV)** and **hemin (HEM)**. No other compound exists in this file.

## Archetype — GIV (givosiran)

**Archetype 3 (depot + central + peripheral, linear) plus one bespoke,
non-standard extra tissue compartment.** The original's `GIV_SC → GIV_C ⇄
GIV_P` block is a textbook depot/central/peripheral linear PK system,
renamed straight to convention: `GUT_GIV`/`CENT_GIV`/`PERI_GIV`. Beyond
that, the original adds a fourth compartment, `GIV_LIV` (liver tissue
concentration, ng/g), fed by a Michaelis–Menten-saturable hepatic uptake
term (`GALN_KM`, representing ASGPR/GalNAc receptor saturation) rather than
a full receptor-binding ODE pair — so this is *not* TMDD (archetype 4;
there is no separate free-receptor/complex state, just a saturable uptake
clearance into a fourth compartment) and doesn't fit archetypes 1–3 as a
whole file. Kept as its own tissue compartment rather than forced into
`PERI_GIV` (which is genuinely a different, linear, plasma-only
peripheral pool in the original) or dropped — renamed `LIV_GIV`, following
the same "extra role beyond the guide's table" precedent as `SYS_VIT` in
`age-related-macular-degeneration/amd_refactor_notes.md`.

`LIV_GIV` — not `CENT_GIV` — is the *only* compartment the disease
equations (ALAS1 mRNA synthesis) ever read, so it is the exposed
concentration: `C_GIV`. Like `TOCI_C` in `sepsis/sep_refactor_notes.md`,
`LIV_GIV` is already declared and dosed/used as a tissue *concentration*
(ng/g), not an amount divided by a volume, so `C_GIV = LIV_GIV` is a
straight alias, no division.

**A genuine dead-parameter observation, preserved as-is:** `V1_GIV` and
`V2_GIV` are declared in `$PARAM` but never referenced anywhere in the
original's `$ODE` — elimination and inter-compartmental transfer are
implemented as `CL_GIV * GIV_C` / `Q_GIV * GIV_C` directly (no `/V1_GIV`
anywhere), so both volumes are inert. Not fixed, not removed, per the
"don't invent, don't drop the original's parameters" rule — carried
through unchanged (still declared, still unused) in the refactored file.

## Archetype — HEM (hemin)

**Archetype 1 (no depot, single central compartment, linear elimination)
plus the same kind of bespoke extra liver tissue compartment as GIV.**
`HEM_C` (renamed `CENT_HEM`) is dosed directly and eliminated by
`-(CL_HEM + KUP_HEM*V_HEM)*HEM_C` — a single linear compartment, matching
archetype 1's shape (like `TOCI_C`/`CENT_TCZ` in the sepsis calibration
run, the compartment is dosed and used directly as a concentration, not
divided by a volume anywhere). `HEM_LIV` (renamed `LIV_HEM`) is a second,
non-standard tissue compartment fed by first-order hepatic uptake
(`KUP_HEM*HEM_C*V_HEM/(BW*0.025)`) and cleared by HO-1 catabolism
(`KHO1_HEM`) — again no receptor-binding ODE, so not TMDD.

`LIV_HEM` is the compartment the ALAS1-suppression term actually reads, so
`C_HEM = LIV_HEM` (straight alias, same "compartment is already a tissue
concentration" reasoning as GIV/TOCI_C).

## The Hill interface

### GIV: rename, not a fit

The original's siRNA knockdown term is already an exact Hill ratio:

```c
double siEFF = EMAX_siRNA * pow(GIV_LIV, HILL_siRNA) /
               (pow(EC50_siRNA, HILL_siRNA) + pow(GIV_LIV, HILL_siRNA));
```

Pulled out verbatim as `EMAX_GIV = 0.92`, `EC50_GIV = 450.0`,
`GAMMA_GIV = 1.50` (all three already had explicit values in the
original — no invented default needed here, unlike some other files'
implicit `gamma=1`). `EFFECT_GIV` is a straight rename of `siEFF`,
computed identically from `C_GIV` (= `LIV_GIV`, unchanged).

### HEM: rename via an algebraic identity, not a fit

The original's hemin suppression term is **not** written in Hill form —
it's an inverse-linear inhibitory expression:

```c
double HEM_EFF = 1.0 / (1.0 + (HEM_LIV / EC50_HEM) * KALAS1_HEM);
```

used directly as a multiplicative factor on ALAS1 synthesis (i.e., `1.0`
= no suppression, `→0` = full suppression as `HEM_LIV → ∞`). This is
exactly a Hill function once solved for its complement: let
`X = KALAS1_HEM * HEM_LIV`. Substituting into `HEM_EFF` and simplifying:

```
1 - HEM_EFF = (X/EC50_HEM) / (1 + X/EC50_HEM) = X / (EC50_HEM + X)
```

— an exact Hill function with `Emax = 1`, `gamma = 1`, potency
`EC50_HEM` acting on the *potency-scaled* variable `X = KALAS1_HEM *
C_HEM`, not on `C_HEM` directly. Verified by hand and by direct numeric
comparison in the run below (`HEM_EFF` in the original vs.
`1 - EFFECT_HEM` in the refactored file matched exactly at every
timepoint tested). `EMAX_HEM = 1.0` and `GAMMA_HEM = 1.0` are added as
explicit `$PARAM` values — both are forced by the original's own algebra
(the asymptote `HEM_EFF → 0` and the absence of any exponent on the
ratio), not invented free parameters. `EC50_HEM` and `KALAS1_HEM` both
keep their original names and values; `KALAS1_HEM` does not have a slot
in the guide's naming table (it's a potency-scaling factor, not one of
`EMAX`/`EC50`/`GAMMA`), so it was left as its original name rather than
folding it into a redefined `EC50_HEM` (which would have silently changed
the *numeric* value exposed under that name from the original's `4.50` to
a derived `4.50/3.50 = 1.2857` — kept both original parameters intact
instead, to avoid that).

In the refactored `$ODE`: `double X_HEM = KALAS1_HEM * C_HEM;` then
`EFFECT_HEM` is the standard Hill expression in `X_HEM`; the original's
`HEM_EFF` multiplicative factor is reconstructed as `1.0 - EFFECT_HEM` at
its one use site (`ALAS1_SYN`).

## `C_<STEM>`/`EFFECT_<STEM>` are `$ODE`-local + `$CAPTURE` only, not `$PARAM`

Attempted putting `C_GIV`, `EFFECT_GIV`, `C_HEM`, `EFFECT_HEM` in `$PARAM`
per the qspserver-compatibility section's literal request. Not possible
under the installed mrgsolve 2.0.1 build, for the same reason confirmed
and documented in `age-related-macular-degeneration/amd_refactor_notes.md`
(reproduced there with a minimal repro model): a `$PARAM` member is passed
into `$ODE` as a read-only reference, so a `double NAME = ...;` line
inside `$ODE` that reuses a `$PARAM` name is treated as an illegal
reassignment (`error: assignment of read-only reference`), not a fresh
per-step local. `$PARAM` and "recomputed every timestep" are mutually
exclusive for the same name in this build — the same conflict that made
`sepsis/sep_refactor_notes.md`'s `C_TCZ` and
`aneurysmal-subarachnoid-hemorrhage/sah_refactor_notes.md`'s `C_NIM`/
`EFFECT_NIM` land as `$CAPTURE`-only diagnostics rather than `$PARAM`
covariates too.
Given that conflict, all four are declared as plain `$ODE`-local `double`s
(`C_GIV`, `EFFECT_GIV`, `C_HEM`, `EFFECT_HEM`) and exposed via `$CAPTURE`
under a distinct `_out`-suffixed name (`C_GIV_out`, `EFFECT_GIV_out`,
`C_HEM_out`, `EFFECT_HEM_out`) so the capture statement itself doesn't
redeclare the same local it's reading — no actual name collision exists
in this file, the suffix is just the same defensive convention used
elsewhere in this fork's refactors. Fully visible via
`POST /run_simulation`'s `outputs` selection; **not** visible via `POST
/model_manifest`'s parameter listing — confirmed empirically (see
Verification below). A future driver-patch redirect of `C_GIV`/`C_HEM`
will need to replace the `double C_GIV = LIV_GIV;` / `double C_HEM =
LIV_HEM;` lines themselves, not override a `$PARAM` default.

## Two pre-existing upstream defects found while verifying (not fixed; logged as `translations/UPSTREAM_ISSUES.md` #39)

Neither is related to GIV or HEM specifically — both block the file from
compiling under mrgsolve 2.0.1 at all, original or refactored, until
worked around:

1. **`$PARAM @annotated` header, but the body uses plain (non-annotated)
   `NAME = value // comment` syntax**, not the colon-delimited annotated
   form mrgsolve's parser requires once `@annotated` is declared —
   `Error: improper annotation format` on the very first parameter line.
2. **`SOLVERTIME` (`_ODETIME_[0]`) used inside `$TABLE`**, where this
   mrgsolve build doesn't expose the macro — `'_ODETIME_' was not declared
   in this scope`, in the `HORM_TRIGGER` capture line, which redundantly
   recomputes the same value the `$ODE`-local `HORMTRIG` already holds.

**Verification-only workaround, applied identically to scratch copies of
both `aip_mrgsolve_model.R` and `aip_mrgsolve_model_refactored.R`, never
to either checked-in file:** dropped `@annotated` from the `$PARAM` header
(no parameter renamed or revalued), and replaced the `$TABLE` capture's
from-scratch `SOLVERTIME` recomputation with a direct reference to the
already-identical `HORMTRIG` local (`capture HORM_TRIGGER = HORMTRIG;`) —
a pure aliasing change, since the two expressions are the same formula.
Full detail in `translations/UPSTREAM_ISSUES.md` entry 39.

## Verification

**Method**: qspserver `mrgsolve_api` container, `POST /model_manifest`
and `POST /run_simulation`, `http://localhost:8007` — no local R/mrgsolve
used, per this fork's verification protocol.

**Scenario tested**: reproduces the original file's own Scenario 4
("Givosiran + Hemin (breakthrough)") — the only one of the original's six
named scenarios that doses both compounds together — via the API's
`dosing` field:

- Givosiran: 162.5 mg (65 kg × 2.5 mg/kg) into compartment 1
  (`GIV_SC`/`GUT_GIV`), `time=0`, `ii=28`, `addl=12` (13 doses total,
  matching the original's `make_givosiran_events(bw=65, dose_mgkg=2.5,
  n_months=13)`).
- Hemin: 195,000 µg (65 kg × 3 mg/kg × 1000) into compartment 12
  (`HEM_C`/`CENT_HEM`), `time=90`, `ii=1`, `addl=3` (4 daily doses,
  matching `make_hemin_events(bw=65, start_day=90, dose_mgkg=3, days=4)`
  called from Scenario 4). The original's `rate = -2` mrgsolve special
  code (its own comment: "zero-order 1h infusion") has no direct
  equivalent in the API's `DoseSpec.rate` field (a real infusion-rate
  number, not mrgsolve's `-1`/`-2` bioavailability-linked codes, and the
  original never defines the `R_HEM_C`/`D_HEM_C` parameter mrgsolve's
  `-2` convention would require to run for real). Translated faithfully
  to the *intended* behaviour instead: a literal zero-order infusion,
  `rate = amt / (1/24 day) = amt × 24`, delivering the same total dose
  over the same 1-hour duration the comment describes. Applied
  identically to both models, so it cannot bias the comparison.
- `parameters`: `CYCTRIG=1, PBGD_ACT=0.5` (matching the scenario's own
  `param()` override).
- `time`: `end=365, delta=0.5` — identical to the original's own
  `mrgsim(end=365, delta=0.5)` call for this scenario; no shortening was
  needed (733 output points ran well within the API's default step
  budget, no `maxsteps`/timeout issue encountered).

**Result: exact match to floating-point noise level.** All 35 shared
outputs (17 `$CAPTURE` diagnostics plus 18 raw compartment/state values,
comparing `GIV_SC`↔`GUT_GIV`, `GIV_C`↔`CENT_GIV`, `GIV_P`↔`PERI_GIV`,
`GIV_LIV`↔`LIV_GIV`, `HEM_C`↔`CENT_HEM`, `HEM_LIV`↔`LIV_HEM` under their
renamed identifiers) were compared point-by-point across all 733
timepoints:

- 29 of 35 outputs: **max abs diff = 0.0, max rel diff = 0.0** (bit-for-bit).
- `GIV_LIVER_ng`/`LIV_GIV`, `PBG_URINE_rate`, `ALA_URINE_rate`, `AUC_PBG`:
  max abs diff = 1e-4, max rel diff ≤ 1.3e-6 — floating-point/solver-step
  noise, not a structural difference (both far below any tolerance that
  would indicate a real bug).
- `HEMIN_PLASMA`/`CENT_HEM`: max abs diff = 3.6e-11 (i.e. zero to
  reporting precision); the reported max *relative* diff of ~1.97 is an
  artifact of dividing two numbers that are each individually ~1e-11 (long
  after the hemin dose has fully cleared, both trajectories are at the
  numerical noise floor around zero) — not a real discrepancy.

This is a pure structural reorganization for both compounds (GIV: exact
Hill rename; HEM: exact algebraic-identity rename, not a curve fit), so
the guide's "near-exact match expected" standard for archetypes 1/3 is
satisfied — in practice the match came out exact to solver-noise level,
same as every prior rename-only refactor in this fork.

**Diagnostic sanity check (no original counterpart to diff against):**
`C_GIV_out` ranges 0–1958 ng/g, `EFFECT_GIV_out` (= `siEFF_frac`, confirmed
identical at every timepoint) rises to ~83% knockdown at steady trough,
consistent with the 92% `EMAX_GIV` asymptote and the achieved liver
concentration. `C_HEM_out` peaks ~37,759 nmol/g during the 4-day hemin
course; `EFFECT_HEM_out` reaches exactly 1.0 at the dosing peak (full
ALAS1 suppression) and decays back toward 0 afterward — both ranges are
physiologically sensible given the model's own parameters.

`/model_manifest` was also called on both models as build-time insurance:
confirmed `EMAX_HEM`, `GAMMA_HEM`, `EC50_HEM`, `KALAS1_HEM` (and every
other renamed `$PARAM`) are listed with their correct values; confirmed
(as expected, see the `$PARAM`-vs-`$CAPTURE` section above) that
`C_GIV`/`EFFECT_GIV`/`C_HEM`/`EFFECT_HEM` do **not** appear in the
manifest's parameter list but do appear as `C_GIV_out`/`EFFECT_GIV_out`/
`C_HEM_out`/`EFFECT_HEM_out` in `outputPaths`.

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`, the
`acute-intermittent-porphyria | GIV` and `acute-intermittent-porphyria |
HEM` rows: target compartments `LIV_GIV`/`LIV_HEM` respectively; exposed
variables `C_GIV`/`EFFECT_GIV` (ng/g liver; fraction 0–1, EMAX 0.92) and
`C_HEM`/`EFFECT_HEM` (nmol/g liver; fraction 0–1, EMAX 1.0).
