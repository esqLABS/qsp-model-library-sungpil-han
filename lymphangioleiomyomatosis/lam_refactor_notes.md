# Refactor notes — `lymphangioleiomyomatosis/lam_mrgsolve_model.R`

**Scope of this pass.** This file has two rows in
`driver-patches/data/compound_perturbation_census.md`, both classified
"Normalize duplicate concentration sites, then redirect":

- **Sirolimus** (stem `SIRO`)
- **Everolimus** (stem `EVER`)

Every disease-side equation (Rheb-GTP, mTORC1/S6K1/4E-BP1 signaling, LAM
cell dynamics, VEGF-D, MMP activity, estrogen, lung cyst volume, FEV1,
DLCO, renal AML volume) is untouched.

## What "duplicate concentration sites" turned out to be

For **each** compound, the original computed its own central-compartment
concentration under two different names, at two different sites, that
never talk to each other:

1. Inside `$ODE`:
   ```
   double Cp_siro  = SIRO_C  / V1_siro * 1000;
   double Cp_ever  = EVER_C  / V1_ever * 1000;
   ```
   feeding the mTOR-inhibition Hill terms (`I_siro`, `I_ever`).
2. Inside `$TABLE`, a second, independent re-read of the same
   `SIRO_C`/`EVER_C` amounts, exposed only for reporting:
   ```
   capture Cp_siro_ngml = SIRO_C / V1_siro * 1000;
   capture Cp_ever_ngml = EVER_C / V1_ever * 1000;
   ```

Both compounds' pair is normalized into one file-scope variable each
(`C_SIRO`, `C_EVER`), declared once in a new `$GLOBAL` block (the original
had no `$GLOBAL` block at all) and assigned once inside `$ODE` — the same
collision-avoidance pattern already established in
`celiac-disease/cd_refactor_notes.md` and
`neonatal-hyperbilirubinemia/nhb_refactor_notes.md` (mrgsolve 2.0.1 hoists
every block-local `double NAME` into one shared anonymous C++ namespace,
so declaring the same name as a local in two different DSL blocks
collides even though each block is textually separate). `Cp_siro_ngml`/
`Cp_ever_ngml` are removed; `C_SIRO`/`C_EVER` are captured directly
instead — a genuine single-site concentration definition per compound
where the original had two.

Each compound also has an unused, single-site peripheral-tissue
concentration (`Cpp_siro`/`Cpp_ever`, `PERI_*/V2_* * 1000`) that is
computed once and never read anywhere else in the file (not captured, not
fed into any Hill term) — this is not a duplicate site and is out of this
census row's scope; kept exactly as-is (informational, unreferenced
local), only its operand names updated to the renamed compartments/
parameters, same pattern as the doxycycline plasma diagnostic in
`abdominal-aortic-aneurysm/aaa_refactor_notes.md`.

## Archetype

**Archetype 3 (depot + central + peripheral, linear)** for both compounds
independently — the original's own 2-compartment-oral structure
(`*_GUT`/`*_C`/`*_P`, CL/Q/V parameterization) is already exactly this
shape; no compartments added or removed.

## Renaming applied

| Original | Refactored | Value | Role |
|---|---|---|---|
| `SIRO_GUT` (cmt) | `GUT_SIRO` | — | sirolimus depot (oral) |
| `SIRO_C` (cmt) | `CENT_SIRO` | — | sirolimus central (amount, mg) |
| `SIRO_P` (cmt) | `PERI_SIRO` | — | sirolimus peripheral (amount, mg) |
| `ka_siro` | `KA_SIRO` | 0.5 (h-1) | absorption rate |
| `F_siro` | `F_SIRO` | 0.15 | bioavailability |
| `CL_siro` | `CL_SIRO` | 10.0 (L/h) | clearance |
| `V1_siro` | `V1_SIRO` | 62.5 (L) | central volume |
| `Q_siro` | `Q_SIRO` | 8.0 (L/h) | inter-compartmental clearance |
| `V2_siro` | `V2_SIRO` | 1400 (L) | peripheral volume |
| `Imax_siro` | `EMAX_SIRO` | 0.90 | Hill Emax |
| `IC50_siro` | `EC50_SIRO` | 2.5 (ng/mL) | Hill EC50 |
| `hill_siro` | `GAMMA_SIRO` | 1.5 | Hill coefficient |
| `Cp_siro` / `Cp_siro_ngml` (two duplicate sites) | `C_SIRO` | — | **the exposed concentration** (single, normalized site, ng/mL) |
| `I_siro` | `EFFECT_SIRO` | — | rename, not a fit — this census row's own named Hill interface |
| `EVER_GUT` (cmt) | `GUT_EVER` | — | everolimus depot (oral) |
| `EVER_C` (cmt) | `CENT_EVER` | — | everolimus central (amount, mg) |
| `EVER_P` (cmt) | `PERI_EVER` | — | everolimus peripheral (amount, mg) |
| `ka_ever` | `KA_EVER` | 1.8 (h-1) | absorption rate |
| `F_ever` | `F_EVER` | 0.30 | bioavailability |
| `CL_ever` | `CL_EVER` | 20.0 (L/h) | clearance |
| `V1_ever` | `V1_EVER` | 90.0 (L) | central volume |
| `Q_ever` | `Q_EVER` | 15.0 (L/h) | inter-compartmental clearance |
| `V2_ever` | `V2_EVER` | 500 (L) | peripheral volume |
| `Imax_ever` | `EMAX_EVER` | 0.88 | Hill Emax |
| `IC50_ever` | `EC50_EVER` | 3.0 (ng/mL) | Hill EC50 |
| `hill_ever` | `GAMMA_EVER` | 1.4 | Hill coefficient |
| `Cp_ever` / `Cp_ever_ngml` (two duplicate sites) | `C_EVER` | — | **the exposed concentration** (single, normalized site, ng/mL) |
| `I_ever` | `EFFECT_EVER` | — | rename, not a fit — this census row's own named Hill interface |
| `I_drug` | `I_drug` (unchanged name) | — | **disease-level combination term**, not a per-compound name in the guide's convention; `1-(1-EFFECT_SIRO)*(1-EFFECT_EVER)`, same shape as the original's `1-(1-I_siro)*(1-I_ever)`, only its two inputs redirected |
| `Cpp_siro` (unchanged name) | `Cpp_siro` | — | peripheral-tissue conc.; single site, unused elsewhere in the original — **not** part of this census row, left as an informational local (operands redirected to `PERI_SIRO`/`V2_SIRO`) |
| `Cpp_ever` (unchanged name) | `Cpp_ever` | — | same status as `Cpp_siro` above |
| `Imax_mTORC2`, `IC50_mTORC2` (unchanged names) | unchanged | 0.15, 20.0 | declared but never referenced anywhere in `$ODE`/`$TABLE` in the original (dead parameters) — left unrenamed and unused; out of scope |

All parameter *values* are copied verbatim from the original.

## Hill interface: rename, not a fit

Both `I_siro = Imax_siro*C^h/(IC50^h+C^h)` and `I_ever = Imax_ever*C^h/
(IC50^h+C^h)` **are already** the canonical
`EMAX*X^gamma/(EC50^gamma+X^gamma)` shape with an explicit Hill
coefficient already present in the original (`hill_siro=1.5`,
`hill_ever=1.4` — renamed `GAMMA_SIRO`/`GAMMA_EVER`, not defaulted to 1).
`EFFECT_SIRO`/`EFFECT_EVER` are one-to-one renames. No `nls()` fit was
needed or performed.

## When the original doesn't compile at all

The untouched original fails to compile under mrgsolve 2.0.1, with two
independent, unrelated defects:

**1. `$CMT` and `$INIT` jointly redeclare all 18 compartments**
```
Error in validObject(.Object) :
  invalid class "mrgmod" object: Duplicated model names: SIRO_GUT SIRO_C SIRO_P
  EVER_GUT EVER_C EVER_P RHEB_GTP MTORC1 S6K1_P EBPP1 LAM_CELLS VEGFD MMP_ACT
  ESTROGEN CYST_VOL FEV1_PCT DLCO_PCT AML_VOL
```
The original's separate `$CMT` block names all 18 compartments (with
descriptive comments), and a following `$INIT` block re-declares the same
18 names with starting values — mrgsolve 2.0.1 treats `$INIT` as its own
compartment-declaring block, the same defect shape as
`idiopathic-pulmonary-fibrosis/ipf_mrgsolve_model.R`
(`translations/UPSTREAM_ISSUES.md` #114).

**2. `$PARAM` name `VEGFD_0` collides with mrgsolve's auto-generated
`<compartment>_0` symbol** for the identically-named compartment `VEGFD`
(only surfaces once defect 1 above is worked around):
```
200:15: error: redeclaration of 'const double& VEGFD_0'
343:15: error: conflicting declaration 'const double& VEGFD_0'
```
`VEGFD_0` (baseline VEGF-D, 400 pg/mL) is declared in `$PARAM` but never
referenced anywhere in `$ODE`/`$TABLE` — dead code, same as
`VEGFD_LAM` (1500 pg/mL) beside it — so this collision was invisible until
defect 1 was worked around and the compartment's own `VEGFD_0` init
symbol was generated. Same defect shape as `sepsis/sep_mrgsolve_model.R`
issue #30 (`$PARAM` name colliding with `<compartment>_0`).

Both are logged as new entries in `translations/UPSTREAM_ISSUES.md`
(**#117**, **#118**); the checked-in `lam_mrgsolve_model.R` is completely
untouched and still carries both defects exactly as written.

**Fix applied directly to the delivered `lam_mrgsolve_model_refactored.R`**
(per the guide's settled policy — syntax-only, non-numeric):
1. The `$CMT` block is dropped; `$INIT` alone declares all 18
   compartments, in the same order as the original's `$CMT` list, so
   1-based dosing indices are unchanged (`cmt=1` is still the sirolimus
   depot, `cmt=4` still the everolimus depot).
2. `VEGFD_0` is renamed `VEGFD_BASELINE`. Since the parameter is unused in
   the original, this changes no numeric value and no model behaviour.

## Verification

**Method.** Both files' embedded DSL blocks were mechanically extracted
(find the quoted-string assignment, pull the contents verbatim — 12,002
chars original, 14,175 chars refactored) and POSTed to the local
qspserver `mrgsolve_api` service at `http://localhost:8007`
(`/model_manifest` then `/run_simulation`), which compiles and runs each
DSL block directly with mrgsolve 2.0.1 server-side — no local R/mrgsolve
install used for either compile or run (a local `Rscript parse()` check
was used only to confirm the refactored file's *R wrapper syntax*, 67
top-level expressions in both the original and refactored `.R`, matching
exactly — not a model-compile check). For the "original" side of the
comparison, the same two disclosed build-compat fixes described above
were applied to a scratch copy only (never to the checked-in
`lam_mrgsolve_model.R`) so that a compiled baseline exists to verify
against at all. Requests were spaced ~2.5s apart and run sequentially
(never more than one in flight), respecting the service's
`max_concurrent_jobs: 2` limit and its history of crashing under
concurrent load. `/model_manifest` confirmed compartment order and
1-based dosing indices are identical between both files (`GUT_SIRO`/
`SIRO_GUT` both index 1, `GUT_EVER`/`EVER_GUT` both index 4).

**Scenarios run — all five of the file's own named scenarios, not
invented, full 2-year duration (`end=17520h, delta=168h`, 105-106 points;
no shortening needed, none approached the API's default `maxsteps`
budget):**

1. `1_Untreated` (no dosing)
2. `2_Sirolimus_2mgQD` (2 mg PO q24h, `addl=729`)
3. `3_Everolimus_10mgQD` (10 mg PO q24h, `addl=729`)
4. `4_Sirolimus_12mo_then_stop` (2 mg PO q24h, `addl=364` — 12 months on,
   then no further dosing through month 24)
5. `5_Everolimus_GnRH` (10 mg PO q24h, `addl=729`, plus the file's own
   `param(E2_basal=0.05, E2_LAMstim=1.0)` GnRH-ablation override)

All dosing amounts/intervals/durations reproduce the file's own
`dose_siro`/`dose_ever`/`dose_siro_12mo` event objects and `lam_GnRH`
parameter override exactly.

**Result: exact match, max abs diff 0.0**, across all 30 shared outputs
(18 compartments + 12 `$CAPTURE`d derived outputs, mapping the original's
`Cp_siro_ngml`/`Cp_ever_ngml` onto the refactored `C_SIRO`/`C_EVER`) at
every timepoint (105 points for the untreated scenario, 106 for each
dosed scenario — the extra row is mrgsolve's own duplicate report row at
the first dose instant, present identically in both files) in **all five
scenarios**. A supplementary dense-PK check (sirolimus q24h, `end=336h,
delta=1h`, matching the original's own `pk_dense` window) also matched
exactly (max abs diff 0.0) and showed no dose-instant reporting artifact
on `C_SIRO` (the `$GLOBAL`-declared, assign-in-`$ODE` pattern used here
follows the same precedent established in
`clostridioides-difficile-infection/cdi_refactor_notes.md` and
`celiac-disease/cd_refactor_notes.md` for avoiding this).

`EFFECT_SIRO`/`EFFECT_EVER` (new, not present in the original) were
sanity-checked directly: `EFFECT_SIRO` is nonzero only in scenarios 2 and
4 (peak 0.0864, sirolimus dosed) and exactly zero in scenarios 1/3/5;
`EFFECT_EVER` is nonzero only in scenarios 3 and 5 (peak 0.3719,
everolimus dosed) and exactly zero in scenarios 1/2/4 — confirming each
Hill interface fires only for its own compound, per the guide's
naming-convention item 4.

This is a pure structural reorganization (rename + duplicate-site
normalization per compound, no Hill-fitting), consistent with the guide's
tolerance table for Archetype 3: "expect a near-exact match... anything
beyond floating-point-scale deviation means a bug." Here the match is
bit-exact, not merely floating-point-scale.

## qspserver `/model_manifest` discoverability

Confirmed via `POST /model_manifest` against the extracted, renamed DSL:
`KA_SIRO`, `F_SIRO`, `CL_SIRO`, `V1_SIRO`, `Q_SIRO`, `V2_SIRO`,
`EMAX_SIRO`, `EC50_SIRO`, `GAMMA_SIRO`, `KA_EVER`, `F_EVER`, `CL_EVER`,
`V1_EVER`, `Q_EVER`, `V2_EVER`, `EMAX_EVER`, `EC50_EVER`, `GAMMA_EVER`
all appear in the manifest's `parameters` with their original numeric
defaults (60 parameters total). `C_SIRO`/`EFFECT_SIRO`/`C_EVER`/
`EFFECT_EVER` are state-derived (computed in `$ODE` from `CENT_SIRO`/
`CENT_EVER` each step), so — per the same reasoning established for
`C_LARA`/`EFFECT_LARA` in `celiac-disease/cd_refactor_notes.md` — they
cannot also be `$PARAM` entries (mrgsolve 2.0.1 treats `$PARAM` members as
read-only in `$ODE`); all four appear in the manifest's `outputPaths` via
the new `$CAPTURE` block, confirmed discoverable. `GUT_SIRO`/`CENT_SIRO`/
`PERI_SIRO`/`GUT_EVER`/`CENT_EVER`/`PERI_EVER` also appear in
`outputPaths` as ordinary compartments (indices 1-6), unchanged position
from the original's `$CMT` order.

No `.cpp` extraction file was left behind — extraction was scratch-only,
used to build the verification requests above and then discarded.

## Anything else flagged

- `Imax_mTORC2`/`IC50_mTORC2` are declared in the original's `$PARAM`
  block with sirolimus-specific comments ("max 15% mTORC2 inhibition by
  sirolimus") but are never referenced anywhere in `$ODE`/`$TABLE` — dead
  parameters describing a pathway (mTORC2/rapalogue resistance) the model
  never actually implements. Not fixed or removed (never-edit-upstream
  rule; and this refactor only touches PK/PD read/write sites, not
  unrelated dead code); left exactly as in the original, unrenamed. Worth
  a look if a future pass wants to actually wire up partial mTORC2
  inhibition.
- The R-side `dose_siro`/`dose_ever`/`dose_siro_12mo` event objects and
  `lam_GnRH` parameter override are unchanged in substance (same amounts,
  same `cmt=1`/`cmt=4` indices, same `ii`/`addl`/`time`); only their
  in-code comments were updated to name the renamed compartments
  (`GUT_SIRO`/`GUT_EVER`). Plot `p3` (drug concentration-time profile)
  now reads `C_SIRO + C_EVER` instead of the removed
  `Cp_siro_ngml + Cp_ever_ngml`. No other R-side code references any
  renamed parameter or compartment by name.
- `lam_shiny_app.R` and `lam_qsp_model.dot`/`.svg`/`.png` are untouched —
  out of scope for this refactor (only the mrgsolve model file and its
  siblings named in the task are deliverables).

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`: both
the `lymphangioleiomyomatosis` / Sirolimus and `lymphangioleiomyomatosis`
/ Everolimus rows' target/pathway (mTORC1, via Rheb-GTP/TSC2), exposed
variable names, units, and notes columns filled in per the above.
