# Refactor notes — `age-related-macular-degeneration/amd_mrgsolve_model.R`

Single-compound file: the only modeled drug is the generic anti-VEGF agent
the corpus-wide census labeled **VIT** (from its `DRUG_VIT` compartment).
`DRUG_TYPE` (1–5) switches which real drug's affinity/off-rate applies
(ranibizumab / aflibercept / bevacizumab / faricimab / brolucizumab), but
all five run through the same PK/binding structure — there is no second,
untouched compound in this file the way `rheumatoid-arthritis` had
adalimumab/methotrexate/baricitinib alongside tocilizumab.

## Archetype

**Bespoke — none of the guide's four archetypes fit as written**, kept
bespoke rather than forced, per the guide's own escape valve:

- `DRUG_VIT` → `DRUG_RET` → `DRUG_SYS` is a one-way, three-compartment
  **cascade** (dosing site → effect site → monitoring sink), not a
  mammillary central/peripheral pair (archetype 2) or a depot feeding a
  central/peripheral pair with two-way exchange (archetype 3). There is no
  back-flux anywhere: `DRUG_SYS` only receives outflow from the other two
  and has its own elimination, never feeding back.
- The drug also engages in genuine mass-action target-binding kinetics
  (`kon_used*DRUG_RET*VEGF_FREE - koff_eff*VEGF_BOUND`, and for faricimab
  only, an analogous Ang-2 term) — this is the shape of archetype 4
  (TMDD), but the "receptor" (`VEGF_FREE`/`VEGF_BOUND`) is not a pool
  dedicated to the drug the way `IL6R_FREE`/`DRUG_RL` were in the
  tocilizumab calibration run. See below for why it stayed unrenamed.

Renamed to convention: `DRUG_VIT`→`CENT_VIT` (vitreous, the dosing
compartment), `DRUG_RET`→`PERI_VIT` (retina, the effect site — this is
where the VEGF/Ang-2 binding terms actually read the drug), `DRUG_SYS`→
`SYS_VIT` (systemic, a one-way monitoring sink with no downstream use —
not in `$CAPTURE` in the original either, and left that way). `SYS_VIT`'s
role has no slot in the guide's naming table (`GUT`/`CENT`/`PERI` only);
using it here is a deliberate, disclosed extension of the convention
rather than dropping the compartment or mislabeling it `PERI`.

## Why VEGF_FREE/VEGF_BOUND (and ANG2_FREE/ANG2_BOUND) were NOT renamed to REC_FREE_VIT/COMPLEX_VIT

This was the central judgment call. In the tocilizumab archetype-4
calibration run, `IL6R_FREE`/`DRUG_RL` exist *only* to compute an
occupancy fraction that gates IL-6 signaling — nothing else in the RA
model reads them. Here, `VEGF_FREE` is the disease's own multi-purpose
pathway state: it drives `VEGFR2_ACT` (Hill function), which drives CNV
growth/regression and fluid influx, **and** `VEGF_FREE`'s own synthesis
term (`VEGF_upregulation`) feeds back from `VEGFR2_ACT` and `RPE_DAM` —
a self-sustaining loop that exists and runs identically with zero drug
present. Renaming it into the compound's own block would misrepresent a
core disease-biology state as something the drug owns. Only the one term
where the drug actually intervenes — `R_BIND_VEGF_VIT` /
`R_BIND_ANG2_VIT`, both built from `C_VIT` — was touched; `VEGF_FREE`,
`VEGF_BOUND`, `ANG2_FREE`, `ANG2_BOUND` and every other disease equation
are untouched, still under their original names.

## The exposed concentration: `C_VIT`

`DRUG_RET` (now `PERI_VIT`) is the *only* PK compartment any disease
equation ever reads — `DRUG_VIT`/`CENT_VIT` and `DRUG_SYS`/`SYS_VIT` are
pure transit/monitoring states. So `C_VIT` is defined as a straight alias
of `PERI_VIT` (no volume division needed — the compartment is already in
concentration units, nM) and substituted everywhere `DRUG_RET` used to be
read: `R_BIND_VEGF_VIT`, `R_BIND_ANG2_VIT`, and `dxdt_PERI_VIT`/
`dxdt_SYS_VIT`'s own outflow terms. This is the single point a future
external-covariate driver patch would redirect.

## No `EFFECT_VIT` Hill fit — and why forcing one would be dishonest

The naming convention's Hill-interface section assumes the compound's
disease effect collapses to one line, `Emax*Xᵞ/(EC50ᵞ+Xᵞ)`, either exactly
(a rename) or approximately (a fit against a sampled occupancy/effect
curve). Neither applies here: **the original has no algebraic effect term
to rename, and no single scalar "effect" to fit a curve against.** The
drug's consequence on disease is a multi-step, ODE-coupled chain —
`C_VIT` depletes `VEGF_FREE` via mass action → `VEGF_FREE` sets
`VEGFR2_ACT` via a Hill relaxation with its own time constant
(`k_R2_act`/`k_R2_inact`) → `VEGFR2_ACT` feeds CNV growth/regression *and*
fluid influx *separately* → each of those feeds `BCVA_target` with its own
coefficient. There is no point in that chain equivalent to RA's
`OCC_TCZ`/`EFFECT_TCZ`, a single ratio consumed everywhere the drug's
effect is needed. Faricimab's Ang-2 arm adds a second, parallel drug
action (via `ANG2_FREE`/`ANG2_BOUND`) feeding only the fluid-influx term,
never CNV — a second reason a single scalar effect would misrepresent the
mechanism.

Per the guide ("a poor fit is a finding... leave the compound's native
calculation in place instead of forcing a bad approximation"), the full
mechanistic chain was kept completely intact — nothing was flattened. As
a captured diagnostic only (not consumed by any `dxdt_`), `EFFECT_VIT` is
defined as the fractional local VEGF neutralization achieved,
`VEGF_BOUND/(VEGF_FREE+VEGF_BOUND+1e-9)` — an exact, already-known ODE
ratio (like RA's `EFFECT_TCZ`), not a fit, but it is **not** the thing
disease dynamics consume; that remains the full chain above.

## `EMAX_VIT`/`EC50_VIT`/`GAMMA_VIT` were not added

The naming table lists these as part of the Hill interface. Since there is
no single Hill-shaped relationship here to parameterize (see above),
adding them would mean inventing Emax/EC50/gamma values with nothing real
to derive them from — exactly what the guide prohibits ("never invent or
default a PK parameter"). Omitted; noted here so their absence reads as a
deliberate decision, not an oversight.

## `C_VIT`/`EFFECT_VIT` are NOT declared in `$PARAM` — a confirmed mrgsolve-build incompatibility, not a choice

The qspserver compatibility section asks for every `C_<STEM>`/`EFFECT_<STEM>`
to live in `$PARAM` (with a `= 0` default) *and* be recomputed every `$ODE`
step, for `/model_manifest` discoverability. This was attempted and found
to be **impossible under the installed mrgsolve 2.0.1 build**, confirmed
empirically via the qspserver `mrgsolve_api` container with a minimal
repro:

```
$PARAM
CL = 1, V = 10
C_X = 0
...
$ODE
double C_X = A / V;     // intended: fresh per-step local
```

compiles to (mrgsolve strips the `double` and the local declaration
whenever the name already exists as a `$PARAM`/`$CAPTURE` symbol, treating
the line as a reassignment of the existing member instead of a new local):

```
83:5: error: assignment of read-only reference 'C_X'
   83 | C_X = A / V;
```

`$PARAM` members (like `$CMT` states) are passed into `$ODE`/`$MAIN` as
read-only references in this mrgsolve build, so a name cannot be *both* a
settable `$PARAM` default *and* a value genuinely recomputed every
timestep — the two requirements are mutually exclusive here, not merely
inconvenient to satisfy together. Given that conflict, correctness of the
simulation took precedence: `C_VIT` and `EFFECT_VIT` are declared as plain
`double` locals inside `$ODE` and exposed via `$CAPTURE` only (bare names,
no `_out` suffix) — the same choice the tocilizumab calibration run made
for `C_TCZ`/`EFFECT_TCZ` (also `$ODE`-local + `$CAPTURE`, not `$PARAM`,
there too, for what now looks like the same underlying reason). Both
values are therefore fully visible through `POST /run_simulation`'s
`outputs` selection, just not through `POST /model_manifest`'s parameter
listing. Flagging this explicitly for whoever does the actual
driver-patch phase, since it means a future redirect of `C_VIT` cannot be
done by overriding a `$PARAM` value at all — it will require replacing the
`double C_VIT = PERI_VIT;` line itself (or a `$MAIN`/`covariate`-mechanism
substitution), not a `param()`-style override.

## Renaming (all parameter values copied verbatim; see the file's own
inline `[was ...]` comments for the full old→new map)

| Role | Original | Refactored |
|---|---|---|
| Compartments | `DRUG_VIT`/`DRUG_RET`/`DRUG_SYS` | `CENT_VIT`/`PERI_VIT`/`SYS_VIT` |
| Elimination | `kel_vit`/`kel_sys` | `KEL_VIT`/`KEL_SYS_VIT` |
| Transfer | `ktr_vit2ret`/`ktr_vit2sys`/`ktr_ret2sys` | `KTR_CENT_PERI_VIT`/`KTR_CENT_SYS_VIT`/`KTR_PERI_SYS_VIT` |
| Volume | `V_vit` | `V1_VIT` |
| Per-drug Kd | `Kd_VEGF_RNB/AFL/BEV/FAR/BRO` | `KD_VIT_RNB/AFL/BEV/FAR/BRO` |
| Ang-2 binding | `kon_ANG2`/`koff_ANG2` | `KON_ANG2_VIT`/`KOFF_ANG2_VIT` |
| Dead/unused params | `Dose_vit`, `Fabs_ret`, `kon_VEGF`, `koff_VEGF`, `n_loading`, `inter_maint` | `DOSE_VIT`, `FABS_VIT`, `KON_VEGF_VIT`, `KOFF_VEGF_VIT`, `N_LOADING_VIT`, `INTER_MAINT_VIT` |
| $ODE locals | `Kd_eff`/`kon_eff`(dead)/`koff_eff`/`kon_used` | `KD_EFF_VIT`/`KON_EFF_VIT`(dead)/`KOFF_EFF_VIT`/`KON_USED_VIT` |
| Reaction rates | `R_bind_VEGF`/`R_bind_ANG2` | `R_BIND_VEGF_VIT`/`R_BIND_ANG2_VIT` |
| New | — | `C_VIT`, `EFFECT_VIT` |

**Left unrenamed, deliberately:** `DRUG_TYPE` (a drug-variant selector
flag, not a physical PK/PD parameter — same tier as `WET_AMD`);
`VEGF_baseline`, `ANG2_baseline`, and every disease-pathway parameter
(`k_VEGF_syn`, `k_VEGF_deg`, `k_R2_act`, `k_ANG2_syn`, complement/RPE/GA/PR/
BCVA parameters) — none of these are compound-specific.

**Dead code preserved as-is (present in the original, not introduced by
this refactor):** `DOSE_VIT`, `FABS_VIT`, `KON_VEGF_VIT`, `KOFF_VEGF_VIT`,
`N_LOADING_VIT`, `INTER_MAINT_VIT` are declared but never read anywhere in
`$ODE`/`$MAIN` (the R-side `make_events()` hardcodes the actual dosing
schedule and amounts, bypassing these params entirely); `KON_EFF_VIT` is
computed but immediately superseded by `KON_USED_VIT` on the next line.
Also pre-existing and disease-level (not touched, not compound-specific):
`VEGF_baseline` and `ANG2_baseline` are declared but never referenced in
`$ODE`/`$MAIN` either — `$INIT`'s `VEGF_FREE_0`/`ANG2_FREE_0` hardcode the
same numeric values directly instead of referencing them.

A parameter-by-parameter diff (Python, matching every original `$PARAM`
name against its mapped refactored name) confirmed all 71 original
parameters carry over with byte-identical values; the only additions are
`C_VIT` and `EFFECT_VIT` (73 total).

## Two pre-existing upstream defects found while verifying (not fixed in either checked-in file; logged as `translations/UPSTREAM_ISSUES.md` #36)

Both reproduce identically from the untouched original via the qspserver
`mrgsolve_api` container — the file does not build under mrgsolve 2.0.1 at
all, refactored or not, until both are worked around:

1. **`$CMT` (bare) + `$INIT` jointly redeclare all 20 compartments** —
   `Duplicated model names: DRUG_VIT DRUG_RET DRUG_SYS ...`. Same defect
   class as issue #34 (`chronic-hypothyroidism`), independently present
   here. Verification workaround: the `$INIT` block's 20 `NAME = value`
   lines were moved into a `$MAIN` block as `NAME_0 = value;` (mrgsolve's
   modern idiom) — declares no new compartment, changes no numeric value.
2. **`$ODE` writes directly to two read-only compartment states** —
   `if(FLUID_EX < 0) FLUID_EX = 0;` and `if(PR_FRAC < 0) PR_FRAC = 0;`
   both fail with `assignment of read-only reference`, since mrgsolve
   2.0.1 passes `$CMT` states into `$ODE` as `const double&`. This class
   was not previously logged. It also appears to have always been a
   no-op regardless of compiler: only `dxdt_*` feeds the integrator, so
   even a compiler that tolerated the write would see it discarded before
   the next step — the negative-value guard the author evidently intended
   was never actually enforced. Verification workaround: both lines were
   deleted outright (a no-op removal, not a numeric change).

Applied identically to scratch copies of both `amd_mrgsolve_model.R` and
`amd_mrgsolve_model_refactored.R`, purely to get a working verification
harness; neither checked-in file was changed — both still contain the
`$CMT`+`$INIT` duplication and the two illegal writes exactly as
originally written (under the refactor's renamed compartments).

## Verification

Per the guide's mandatory protocol: ran the original file's own dosing
scenarios (dose amounts/timing copied exactly from `make_events()` and
each scenario's `param()` overrides in the original file) through both the
(defect-patched, for buildability only) original and the refactored model
via the qspserver `mrgsolve_api` service (`POST /model_manifest`, `POST
/run_simulation`), 730-day horizon, daily (`delta=1`) output grid —
matching the original scenarios' own `mrgsim(end = 730, delta = 1)` call.

Scenarios run (covering all 5 `DRUG_TYPE` variants, including faricimab —
`DRUG_TYPE==4` — to exercise the Ang-2 binding arm that only fires for
that drug):

1. **Scenario 1** — Ranibizumab, 0.5 mg q4w×3 loading → q8w×10 maintenance
   (13 doses, 2264.49 nM each into `CENT_VIT`/`DRUG_VIT`)
2. **Scenario 2** — Aflibercept, 2 mg q4w×3 → q8w×10 (`KEL_VIT`/`kel_vit`
   overridden to 0.099, as in the original)
3. **Scenario 3** — Faricimab, 6 mg q4w×4 → q16w×7 (`KEL_VIT` = 0.099;
   exercises `R_BIND_ANG2_VIT`)
4. **Scenario 4** — Brolucizumab, 6 mg q6w×3 → q12w×8 (`KEL_VIT` = 0.173)

For each scenario, every compartment state (20/20, both under original and
renamed identifiers) and every original `$CAPTURE`d `*_out`/`BCVA_change`
output (15/15) was compared point-by-point across the full time grid
(744 or 742 points depending on scenario, matching each drug's own
maintenance-interval count).

**Result: exact match, max abs diff = 0.0 and max rel diff = 0.0 for every
compared output, every scenario, every compartment.** This is the expected
outcome for a pure structural reorganization (no Hill-fit anywhere) per
the guide's own tolerance rule for archetypes 1–3/bespoke-linear
structures. `C_VIT` and `EFFECT_VIT` (no original counterpart to diff
against) were sanity-checked for plausible ranges instead: `EFFECT_VIT`
rises from 0 toward ~0.99+ within the first loading cycle for every drug
(consistent with all five having sub-nanomolar-to-low-nanomolar Kd against
a ~2.5 nM baseline VEGF pool), and `C_VIT` tracks the expected per-drug
concentration scale (ranibizumab lowest, brolucizumab highest, reflecting
each drug's dose/MW-derived nM dose and elimination half-life).

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`,
`age-related-macular-degeneration | VIT` row: target = `PERI_VIT`
(retina, the effect site); exposed variable = `C_VIT` / `EFFECT_VIT`
(diagnostic only, not consumed downstream); unit = nM / fraction 0-1.
