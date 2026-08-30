# Refactor notes — `vwd_mrgsolve_model.R` (DDAVP, recombinant VWF, plasma-derived VWF/FVIII concentrate, tranexamic acid)

**Scope of this pass.** Per the fork's PK/PD refactor spec
([`FORK_WORKFLOW_GUIDE.md`](../FORK_WORKFLOW_GUIDE.md), Part 2), all 4
compounds this file models were rewritten, per their existing rows in
[`driver-patches/data/compound_perturbation_census.md`](../driver-patches/data/compound_perturbation_census.md)
(all classified "Redirect concentration (clean single site)"): **DDAVP**
(desmopressin, stem `DDAVP`), **RVWF** (recombinant VWF/vonicog
alfa-like, stem `RVWF`), **PDVWF** (plasma-derived VWF/FVIII concentrate,
Humate-P-like, stem `PDVWF`), and **TXA** (tranexamic acid, stem `TXA`).
Every disease-side equation (WPB store dynamics, endogenous VWF:Ag/:RCo
turnover, HMWM/ADAMTS13, FVIII stabilization, platelet clearance,
bleeding score, menstrual/GI loss, hemoglobin, sodium, thrombotic-risk
index) is untouched apart from reading the four renamed concentration
variables at their existing use sites.

## Are RVWF and PDVWF genuinely two separate products?

Confirmed yes, not one product mislabeled twice. They are dosed into
separate compartments (`RVWF_CP` vs `PDVWF_CP`), used in different
scenarios (`5_Type3_Severe_rVWF_Vonvendi` doses only `RVWF_CP`;
`6_Type3_Severe_PDConcentrate`/`10_Surgical_Major_PKguided_PDConc` dose
only `PDVWF_CP`), have different elimination rates (`KE_RVWF=0.033/h`,
t1/2 ~21h per Mannucci 2013, vs `KE_PDVWF=0.045/h`, t1/2 ~15h per
Dobrkovska 1998), and — mechanistically distinct — feed FVIII
stabilization by two genuinely different routes: RVWF (VWF-only,
no co-formulated FVIII) only raises `FVIII_STAB_TARGET`, so endogenous
`FVIII_C` climbs *gradually* toward it at `KOUT_FVIII` (t1/2 ~9h),
reproducing the delayed FVIII:C rise reported for vonicog alfa (Gill
2015); PDVWF (co-formulated FVIII, ratio `VWF_FVIII_RATIO=2.4`)
contributes `FVIII_EXO_DIRECT = C_PDVWF/VWF_FVIII_RATIO` *immediately*,
with its own PK decay, on top of the endogenous pool. This distinction
survives the refactor unchanged.

## Archetypes determined

### DDAVP — Archetype 3 minus the peripheral compartment (depot + central, linear elimination)

The original tracks the *central compartment as a concentration state*
directly (`DDAVP_CP`, "ng/mL equiv"), not as an amount divided by volume:

```
dxdt_DDAVP_DEPOT = -KA_DDAVP * DDAVP_DEPOT;
dxdt_DDAVP_CP    =  KA_DDAVP * DDAVP_DEPOT / V_DDAVP - KE_DDAVP * DDAVP_CP;
```

For a constant volume this is algebraically identical to the standard
"amount in `CENT`, `C = CENT/V1`" form used elsewhere in this corpus:
letting `A = C*V`, `dA/dt = V*dC/dt = KA*Depot - KE*V*C = KA*Depot -
KE*A` — exactly the guide's Archetype-3 shape with the peripheral
compartment dropped (matching the ATRA precedent in
`disseminated-intravascular-coagulation/dic_refactor_notes.md`). Bioavailability
(`F_DDAVP`) is applied by the original via mrgsolve's `F_<cmt>` dosing-time
mechanism (`F_DDAVP_DEPOT = F_DDAVP;` in `$MAIN`), **not** as a factor
inside the ODE flow term — the refactor preserves this exactly
(`F_GUT_DDAVP = F_DDAVP;`), because adding `F_DDAVP` a second time inside
`dxdt_CENT_DDAVP` would double-apply it.

The WPB-release drive (`DDAVP_DRIVE = pow(DDAVP_CP, HILL_DDAVP) /
(pow(EC50_DDAVP, HILL_DDAVP) + pow(DDAVP_CP, HILL_DDAVP))`) is already the
guide's canonical Hill ratio, with an implicit ceiling of 1.0 (no
separate Emax factor in the original) and `HILL_DDAVP=1.3` as its
exponent.

**Renaming applied (DDAVP only):**

| Original | Refactored | Value |
|---|---|---|
| `DDAVP_DEPOT` (compartment) | `GUT_DDAVP` | — |
| `DDAVP_CP` (compartment, tracked concentration) | `CENT_DDAVP` (tracked amount; `C_DDAVP = CENT_DDAVP/V1_DDAVP`) | — |
| `KA_DDAVP` | unchanged | 4.5 (1/h) |
| `F_DDAVP` | unchanged (still applied via `F_GUT_DDAVP` at dosing time) | 0.15 |
| `KE_DDAVP` + `V_DDAVP` | `CL_DDAVP` (new, `= KE_DDAVP x V_DDAVP`) + `V1_DDAVP` | `CL_DDAVP = 6.9` (L/h), `V1_DDAVP = 30.0` (L) |
| `EC50_DDAVP` | unchanged (already convention-named) | 0.8 (ng/mL) |
| `HILL_DDAVP` | `GAMMA_DDAVP` | 1.3 |
| — (none; original ratio has an implicit ceiling of 1.0) | `EMAX_DDAVP` (new, `= 1.0`) | documented below |
| `DDAVP_DRIVE` (local `$ODE` ratio) | `EFFECT_DDAVP` | — |

**`CL_DDAVP` is derived, not invented** — `CL_DDAVP := KE_DDAVP x V_DDAVP
= 0.23 x 30.0 = 6.9`, so `CL_DDAVP/V1_DDAVP` reduces to exactly `0.23`,
same technique as `CL_ATR`/`CL_NIM` in prior calibration runs.
`EMAX_DDAVP = 1.0` is not a free parameter — it is the original ratio's
own implicit asymptote (`DDAVP_DRIVE -> 1` as `DDAVP_CP -> inf`) made
explicit, required to reproduce the original term exactly under the
guide's canonical `EMAX*X^gamma/(EC50^gamma+X^gamma)` form.

`EFFECT_DDAVP` replaces `DDAVP_DRIVE` at its one use site
(`WPB_RELEASE_RATE = K_WPBDEPLETE * WPB_STORE * EFFECT_DDAVP`). `C_DDAVP`
is additionally read at a second, unrelated use site (`NA_DROP`'s
sodium-retention term) — this is the same single exposed concentration
variable read twice by two different disease-side effect equations, not
two separate concentration variables, consistent with the guide's "two
only when a genuinely different tissue site matters" allowance being
unnecessary here (same site, same variable, two equations).

### RVWF and PDVWF — Archetype 1 (no depot, single compartment, linear elimination), each

```
dxdt_RVWF_CP  = -KE_RVWF  * RVWF_CP;    // dosed directly (IU) via bolus
dxdt_PDVWF_CP = -KE_PDVWF * PDVWF_CP;   // dosed directly (IU) via bolus
```

Both are dosed directly as circulating amount (IU), no depot, no
peripheral compartment — Archetype 1 exactly, once expressed as
`CL/V1`.

**Renaming applied:**

| Original | Refactored | Value |
|---|---|---|
| `RVWF_CP` (compartment) | `CENT_RVWF` | — |
| `KE_RVWF` + `V_RVWF` | `CL_RVWF` (new, `= KE_RVWF x V_RVWF`) + `V1_RVWF` | `CL_RVWF = 1.155` (dL/h), `V1_RVWF = 35.0` (dL) |
| `RVWF_CONC` (local `$ODE` double) | `C_RVWF` | — |
| `PDVWF_CP` (compartment) | `CENT_PDVWF` | — |
| `KE_PDVWF` + `V_PDVWF` | `CL_PDVWF` (new, `= KE_PDVWF x V_PDVWF`) + `V1_PDVWF` | `CL_PDVWF = 1.575` (dL/h), `V1_PDVWF = 35.0` (dL) |
| `PDVWF_CONC` (local `$ODE` double) | `C_PDVWF` | — |
| `VWF_FVIII_RATIO` | unchanged | 2.4 (formulation ratio, not a PK role in the naming convention) |

`CL_RVWF := 0.033 x 35.0 = 1.155`; `CL_PDVWF := 0.045 x 35.0 = 1.575` —
both derived, not invented; `CL/V1` reduces to exactly the original
`KE_*` in each case.

### No `EFFECT_RVWF` / `EFFECT_PDVWF` — bespoke, disclosed

Unlike DDAVP and TXA, neither RVWF nor PDVWF has a saturating Hill/Emax
pharmacology in the original. Their concentrations feed the disease side
as **direct, unbounded linear mass-balance contributions**, not as a
receptor-occupancy effect:

- `EXO_VWF_TOTAL = C_RVWF + C_PDVWF` (linear sum, both compounds)
- `FVIII_STAB_DRIVE = fmin(1.0, (VWF_AG + FVIII_STAB_GAIN * C_RVWF) / (VWFAG_BASE + 1e-6))` — the only "saturation" here is a shared `fmin` ceiling applied to the *combined endogenous+RVWF* drive, not a per-compound Hill term on RVWF alone
- `FVIII_EXO_DIRECT = C_PDVWF / VWF_FVIII_RATIO` (linear)

Forcing an `EMAX*C^gamma/(EC50^gamma+C^gamma)` shape onto an unbounded
linear replacement term would misrepresent pharmacology the original
never modeled (per the guide's "a clean, non-standard structure beats a
standard structure that's wrong" and this fork's instruction not to force
a bad fit). **Decision: expose only `C_RVWF`/`C_PDVWF`** (satisfying
"redirect concentration, clean single site") and leave the existing
linear use sites as pure renames — zero numeric change, confirmed below.

## TXA — Archetype 3 minus the peripheral compartment (depot + central, linear elimination)

Same concentration-tracking pattern as DDAVP:

```
dxdt_TXA_DEPOT = -KA_TXA * TXA_DEPOT;
dxdt_TXA_CP    =  KA_TXA * TXA_DEPOT / V_TXA - KE_TXA * TXA_CP;
```

Bioavailability applied the same way, via `F_TXA_DEPOT = F_TXA;` in
`$MAIN` (renamed `F_GUT_TXA`), not inside the ODE flow term.

`TXA_EFFECT = EMAX_TXA_BLEED * TXA_CP / (EC50_TXA + TXA_CP)` is already
the guide's canonical Hill ratio with an implicit exponent of 1.

**Renaming applied:**

| Original | Refactored | Value |
|---|---|---|
| `TXA_DEPOT` (compartment) | `GUT_TXA` | — |
| `TXA_CP` (compartment, tracked concentration) | `CENT_TXA` (tracked amount; `C_TXA = CENT_TXA/V1_TXA`) | — |
| `KA_TXA` | unchanged | 1.1 (1/h) |
| `F_TXA` | unchanged (still applied via `F_GUT_TXA` at dosing time) | 0.35 |
| `KE_TXA` + `V_TXA` | `CL_TXA` (new, `= KE_TXA x V_TXA`) + `V1_TXA` | `CL_TXA = 7.25` (L/h), `V1_TXA = 25.0` (L) |
| `EMAX_TXA_BLEED` | `EMAX_TXA` | 0.30 |
| `EC50_TXA` | unchanged (already convention-named) | 8.0 (mcg/mL) |
| — (none; original had no Hill exponent) | `GAMMA_TXA` (new, `= 1.0`) | documented as "original had no explicit Hill term" |
| `TXA_EFFECT` (local `$ODE` ratio) | `EFFECT_TXA` | — |

`CL_TXA := 0.29 x 25.0 = 7.25`, derived not invented. `EFFECT_TXA`
replaces `TXA_EFFECT` at both its use sites (`BLEED_TARGET` and
`MENS_TARGET`).

## A DDAVP-specific numerical-robustness fix (in-scope, disclosed)

`EFFECT_DDAVP`'s fractional exponent (`GAMMA_DDAVP=1.3`, non-integer)
NaNs if `C_DDAVP` numerically undershoots zero — which it does, at
~1e-10 scale, deep into a long low-dose decay tail. **This is a
pre-existing fragility in the original**, not introduced by the
refactor: the original's own `DDAVP_CP` hits the same wall (confirmed
below). Because this sits inside DDAVP's own scope block, it is fixed
directly as part of this refactor (per the guide's point 4 for
in-scope defects) rather than only logged:

```c
double C_DDAVP_nn = fmax(C_DDAVP, 0.0);
double EFFECT_DDAVP = EMAX_DDAVP * pow(C_DDAVP_nn, GAMMA_DDAVP)
                     / (pow(EC50_DDAVP, GAMMA_DDAVP) + pow(C_DDAVP_nn, GAMMA_DDAVP));
```

This changes nothing for any `C_DDAVP >= 0` (every physiologically
meaningful value) and is disclosed here and in `UPSTREAM_ISSUES.md` #72.

## A pre-existing upstream compile defect, unrelated to this refactor

Neither `vwd_mrgsolve_model.R` nor the refactored file compiles as
written under mrgsolve 2.0.1, for two independent reasons:

1. **`$CAPTURE` duplicates 12 of the file's own `$CMT` compartment
   names** (`ADAMTS13_ACT`, `PLT_COUNT`, `BLEED_SCORE`, `MENS_LOSS`,
   `GI_LOSS`, `HB`, `NA_SERUM`, `THROMB_RISK`, `WPB_STORE`, `VWFPP`,
   `DDAVP_CP`, `TXA_CP`) — rejected outright (`compartment should not be
   in $CAPTURE`).
2. **`$MAIN` sets 14 bare compartment names directly** inside
   `if (NEWIND <= 1) { WPB_STORE = WPB0; ... }` instead of this build's
   `<CMT>_0` initial-value idiom — rejected as `assignment of read-only
   reference`.

Both are present in the *original* file already and unrelated to any of
the 4 compounds' own PK. Logged as **upstream issue #72** in
[`translations/UPSTREAM_ISSUES.md`](../translations/UPSTREAM_ISSUES.md).
Per the guide's settled policy for a non-compiling original, both fixes
are applied directly, syntax-only and non-numeric, into the delivered
`vwd_mrgsolve_model_refactored.R` — the tracked `vwd_mrgsolve_model.R`
itself is untouched and still carries both defects exactly as written.

## Diff scope confirmation

Beyond the two build-compat fixes and the compound renames/effect
definitions described above, nothing else in the file differs: every
disease-side `$PARAM`/`$CMT`/`$ODE` line, the aetiology/genotype presets
(`p_type1`, `p_type2b`, `p_type3`, `p_acquired`, `p_preg`), `make_ev()`,
`run_scenario()`, and the example-run/plotting block are byte-identical
except for the compartment-name strings in the scenario list
(`"DDAVP_DEPOT"`→`"GUT_DDAVP"` x4, `"RVWF_CP"`→`"CENT_RVWF"`,
`"PDVWF_CP"`→`"CENT_PDVWF"` x2, `"TXA_DEPOT"`→`"GUT_TXA"`) so the
refactored file's own scenarios still dose the renamed compartments by
name, and two calibration-note comment lines updated to reference the
renamed identifiers (prose only, no numeric content).

## qspserver `/model_manifest` discoverability

- `KA_DDAVP`, `F_DDAVP`, `CL_DDAVP`, `V1_DDAVP`, `EC50_DDAVP`,
  `GAMMA_DDAVP`, `EMAX_DDAVP`, `CL_RVWF`, `V1_RVWF`, `CL_PDVWF`,
  `V1_PDVWF`, `VWF_FVIII_RATIO`, `KA_TXA`, `F_TXA`, `CL_TXA`, `V1_TXA`,
  `EMAX_TXA`, `EC50_TXA`, `GAMMA_TXA` are all declared in `$PARAM`,
  confirmed present with their defaults in the live `/model_manifest`
  response (e.g. `CL_DDAVP=6.9`, `GAMMA_DDAVP=1.3`, `EMAX_TXA=0.3`,
  `GAMMA_TXA=1`).
- `C_DDAVP`, `EFFECT_DDAVP`, `C_RVWF`, `C_PDVWF`, `C_TXA`, `EFFECT_TXA`
  are computed in `$ODE` from compartment state (so, as found repeatedly
  in prior runs, cannot also be `$PARAM` entries). Unlike the
  `dic`/`sah`/`ted` precedents, this file's own `$CAPTURE` block already
  uses the simple bare-name list form (`$CAPTURE NAME1 NAME2 ...`, no
  `$TABLE`-style `capture NAME = expr;` re-declaration), so no naming
  collision arises here — all 6 quantities were added directly to the
  existing `$CAPTURE` list under their own names and are confirmed
  present in `/model_manifest`'s `outputPaths` and retrievable via
  `/run_simulation`, with no renamed/`o`-prefixed workaround needed.
- `model_content` was extracted from the `vwd_code <- '...'` R-string
  wrapper (straight text extraction, confirmed byte-identical to the
  quoted block programmatically); no `.cpp` sibling was left behind —
  extraction was in-memory only, used to build the verification requests
  below and then discarded, per the workflow guide.

## Verification

**Method.** Both the original's DSL (patched only with the two
build-compat fixes above, for verification purposes — this patched copy
exists only in transient in-memory request payloads, never written to
disk) and the refactored file's DSL were POSTed to the local qspserver
`mrgsolve_api` service at `http://localhost:8007`
(`/model_manifest` then `/run_simulation`), which compiles and runs each
DSL block directly with mrgsolve 2.0.1 server-side — no local R/mrgsolve
install used, requests spaced ~2-3s apart per the task's concurrency
note. Four of the file's own scenarios were run through both models,
identical dosing/genotype-preset parameters as `run_scenario()` defines
them:

1. **`3_Type1_DDAVP_Intranasal_Repeat`** — Type-1 preset (`VWFAG_BASE=35,
   VWFRCO_BASE=30, FVIII_BASE=42, HMWM_BASE=0.85`), 300 mcg q12h x6 doses
   into the DDAVP depot. **Window shortened to `end=160h`** (from the
   scenario's own `end=168h`) — both the original (patched) and the
   refactored file NaN in this scenario at very low residual DDAVP
   concentration (~1e-10 scale, original at t=167h, refactored — before
   its own `fmax` guard — at t=160h; see "numerical-robustness fix"
   above), so 160h is the last point where a like-for-like comparison is
   possible against the *original's* unguarded behavior. This is a
   solver/floating-point-fragility limitation, not a refactor mismatch,
   per the guide's "shorten the window rather than treating [this] as a
   mismatch" rule — disclosed honestly rather than silently swapping
   scenarios. 18 shared disease-side outputs, 161 timepoints.
2. **`5_Type3_Severe_rVWF_Vonvendi`** — Type-3 preset (`VWFAG_BASE=2,
   VWFRCO_BASE=1, FVIII_BASE=6, HMWM_BASE=0.02`), 3500 IU q24h x7 doses
   into the RVWF compartment, full `end=168h`. 18 shared outputs, 169
   timepoints.
3. **`10_Surgical_Major_PKguided_PDConc`** — Type-3 preset, 4200 IU q12h
   x21 doses into the PDVWF compartment, full `end=240h`. 18 shared
   outputs, 241 timepoints.
4. **`7_Menorrhagia_TXA_Hormonal`** — Type-1 preset + `HORMONAL_ON=1`,
   1300 mg q8h x90 doses into the TXA depot, full `end=720h`. 18 shared
   outputs, 721 timepoints.

Compared across the file's 18 disease-side `$CAPTURE`/compartment
outputs shared verbatim between original and refactored (`VWF_AG_TOTAL`,
`VWF_RCO_TOTAL`, `HMWM_EFFECTIVE`, `FVIII_C_TOTAL`, `PLT_COUNT`,
`BLEED_SCORE`, `MENS_LOSS`, `GI_LOSS`, `HB`, `NA_SERUM`, `THROMB_RISK`,
`WPB_STORE`, `VWFPP`, `VWF_AG`, `VWF_RCO`, `HMWM`, `ADAMTS13_ACT`,
`FVIII_C` — untouched by any of the 4 compounds' renaming, so any
numeric drift here would indicate a real bug); the renamed PK
compartments/concentrations themselves were also compared by fetching
each model's own names and checking peak values.

**Results:**

- **RVWF scenario (5): exact match.** Max relative deviation = **0.0**
  across all 18 shared outputs. `RVWF_CP` (orig) vs `CENT_RVWF` (refac.)
  peaked at identical `6372.7916` IU; `RVWF_CONC` vs `C_RVWF` peaked at
  identical `182.0798` IU/dL. Expected for a pure Archetype-1 rename with
  no effect-interface change (no Hill-fitting involved).
- **PDVWF scenario (10): exact match.** Max relative deviation = **0.0**
  across all 18 shared outputs. `PDVWF_CP` vs `CENT_PDVWF` peaked at
  identical `10065.746` IU; `PDVWF_CONC` vs `C_PDVWF` peaked at identical
  `287.5927` IU/dL.
- **TXA scenario (7): near-exact match.** Max relative deviation across
  the 18 shared outputs was **1.54e-6** (worst: `MENS_LOSS` at t=167h,
  `64.7748` vs `64.7749`) — floating-point-scale noise, same order as the
  `EMAX*C/(EC50+C)` re-expression's effect on this corpus's stiff coupled
  systems previously root-caused in `dic_refactor_notes.md`'s TXA finding
  (bit-reassociation between the original's direct ratio and the
  `pow(C,1.0)` form). `TXA_DEPOT`/`TXA_CP` vs `GUT_TXA`/`CENT_TXA`/`C_TXA`
  peaked at identical values (`455.0686` mg depot, `12.6076` mcg/mL
  concentration).
- **DDAVP scenario (3): near-exact match once the sub-float-precision
  noise floor is excluded.** Raw max relative deviation across the 18
  shared outputs over the shortened 0-160h window was **0.144** at
  `VWFPP`, t=155h — but both compared values there are **~7e-10 and
  ~6.5e-10** (VWFPP "IU/dL equiv", i.e. six-plus orders of magnitude
  below any clinically or analytically meaningful scale), so a 14%
  *relative* deviation between two numerically-insignificant noise-floor
  values is not a real discrepancy. Excluding all timepoints where both
  compared values are below 1e-3 (i.e. restricting to physiologically
  meaningful magnitudes), **max relative deviation drops to 7.15e-7**
  (worst: `NA_SERUM`, `139.792` vs `139.7919`) — floating-point noise,
  consistent with the TXA finding above. `DDAVP_DEPOT` vs `GUT_DDAVP`
  peaked at identical `45` mcg; `DDAVP_CP` vs `C_DDAVP` peaked at
  identical `1.3233` ng/mL. The underlying NaN-onset-timing difference
  between the two models (t=167h original vs t=160h refactored, before
  the `fmax` guard) is a separate, disclosed floating-point-fragility
  finding (see above and `UPSTREAM_ISSUES.md` #72), not part of this
  numeric-equivalence result.

**No mismatch was found that the analysis above did not fully explain.**
Per the guide's tolerance table: Archetype 1 (RVWF, PDVWF) reproduced
bit-for-bit; Archetype 3 (DDAVP, TXA) reproduced to floating-point-scale
tolerance (≤7.15e-7 and ≤1.54e-6 respectively) once sub-float-precision
noise-floor artifacts are excluded — both consistent with pure structural
reorganization, no Hill-refitting involved for any of the 4 compounds.
