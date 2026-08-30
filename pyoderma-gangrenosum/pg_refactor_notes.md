# Pyoderma Gangrenosum (PG) — PK/PD refactor notes

Refactor of `pg_mrgsolve_model.R` per `FORK_WORKFLOW_GUIDE.md` Part 2
("pluggable PK, named Hill interface"). Original untouched; new sibling
`pg_mrgsolve_model_refactored.R` added. Scope: all six compounds this file
models — Adalimumab (ADA), Anakinra (ANA), Cyclosporine (CSA), Infliximab
(IFX), Prednisone (PRED), Ustekinumab (UST) — there is no other compound
in this file to leave untouched.

## Archetype per compound

| Compound | Archetype | Notes |
|---|---|---|
| Adalimumab (ADA) | 3 (depot + central + peripheral, linear) | SC dosing, two-compartment disposition |
| Infliximab (IFX) | 2 (no depot [IV], 2-cmt, linear) | IV dosing directly into central |
| Anakinra (ANA)   | 3 minus peripheral (depot + central, linear) | SC dosing, one-compartment disposition |
| Cyclosporine (CSA) | 3 minus peripheral (depot + central, linear) | Oral dosing, one-compartment disposition |
| Ustekinumab (UST) | 3 minus peripheral (depot + central, linear) | See TMDD check below |
| Prednisone (PRED) | 3 minus peripheral, no bioavailability term | Oral dosing; original declares no F for PRED |

**Ustekinumab TMDD check** (called out explicitly in the task, since UST is
an anti-IL-12/23 antibody where receptor-binding kinetics would be
plausible): this file's own `$ODE` for UST is

```c
dxdt_UST_SC   = -KA_UST * UST_SC;
dxdt_UST_CENT =  KA_UST * F1_UST * UST_SC - (CL_UST/V_UST) * UST_CENT;
```

— a plain one-depot, linear-elimination compartment. There is no receptor
pool, on/off-rate parameter, or complex compartment for UST (or for any of
the other five compounds) anywhere in the file, and its effect on disease
(`EFF_UST_IL23 = Emax_drug * CONC_UST / (IC50_UST_IL23 + CONC_UST)`) is a
plain algebraic Emax ratio on concentration, not on a receptor-occupancy
state. Confirmed **not TMDD** in this file; kept as the plain-PK archetype
it actually is, per the guide's "match the archetype closest to what the
original already does" rule.

## Naming map

| Role | Original | Refactored |
|---|---|---|
| ADA depot/central/peripheral | `ADA_SC` / `ADA_CENT` / `ADA_PERI` | `GUT_ADA` / `CENT_ADA` / `PERI_ADA` |
| ADA bioavailability | `F1_ADA` | `F_ADA` |
| ADA conc / effect | `CONC_ADA` / `EFF_ADA_TNF` | `C_ADA` / `EFFECT_ADA` |
| IFX central/peripheral | `IFX_CENT` / `IFX_PERI` | `CENT_IFX` / `PERI_IFX` |
| IFX conc / effect | `CONC_IFX` / `EFF_IFX_TNF` | `C_IFX` / `EFFECT_IFX` |
| ANA depot/central | `ANA_SC` / `ANA_CENT` | `GUT_ANA` / `CENT_ANA` |
| ANA bioavailability / volume | `F1_ANA` / `V_ANA` | `F_ANA` / `V1_ANA` |
| ANA conc / effect | `CONC_ANA` / `EFF_ANA_IL1` | `C_ANA` / `EFFECT_ANA` |
| CSA depot/central | `CSA_GUT` / `CSA_CENT` | `GUT_CSA` / `CENT_CSA` |
| CSA bioavailability / volume | `F1_CSA` / `V_CSA` | `F_CSA` / `V1_CSA` |
| CSA conc / effect | `CONC_CSA` / `EFF_CSA_Th17` | `C_CSA` / `EFFECT_CSA` |
| UST depot/central | `UST_SC` / `UST_CENT` | `GUT_UST` / `CENT_UST` |
| UST bioavailability / volume | `F1_UST` / `V_UST` | `F_UST` / `V1_UST` |
| UST conc / effect | `CONC_UST` / `EFF_UST_IL23` | `C_UST` / `EFFECT_UST` |
| PRED depot/central | `PRED_GUT` / `PRED_CENT` | `GUT_PRED` / `CENT_PRED` |
| PRED volume | `V_PRED` | `V1_PRED` |
| PRED conc / effect | `CONC_PRED` / `EFF_PRED` | `C_PRED` / `EFFECT_PRED` |
| Hill EC50 (all six) | `IC50_ADA_TNF`, `IC50_IFX_TNF`, `IC50_ANA_IL1`, `IC50_UST_IL23`, `IC50_CSA_Th17`, `IC50_PRED` | `EC50_ADA`, `EC50_IFX`, `EC50_ANA`, `EC50_UST`, `EC50_CSA`, `EC50_PRED` |
| Hill Emax (all six) | shared `Emax_drug` | `EMAX_ADA`, `EMAX_IFX`, `EMAX_ANA`, `EMAX_UST`, `EMAX_CSA`, `EMAX_PRED` (see shared-EMAX finding below) |
| Hill gamma (all six) | not present | `GAMMA_ADA` = `GAMMA_IFX` = `GAMMA_ANA` = `GAMMA_UST` = `GAMMA_CSA` = `GAMMA_PRED` = 1 (added; original had no explicit Hill exponent, ratio was already gamma=1 shaped) |

`CL_ADA`, `V1_ADA`, `V2_ADA`, `Q_ADA`, `KA_ADA`, `CL_IFX`, `V1_IFX`,
`V2_IFX`, `Q_IFX`, `KA_ANA`, `CL_ANA`, `KA_CSA`, `CL_CSA`, `KA_UST`,
`CL_UST`, `KA_PRED`, `CL_PRED` were already spelled exactly per the fork
convention in the original and are untouched.

All six Hill effect terms were already exact ratios
(`Emax*C/(IC50+C)`, implicit gamma=1) in the original — every `EFFECT_<STEM>`
is a **rename, not a refit**.

## The six-way shared `Emax_drug` — a found, disclosed, and independently-fixed leak

The original reads a single `Emax_drug = 0.92` parameter as the ceiling for
all six compounds' effect terms (`EFF_ADA_TNF`, `EFF_IFX_TNF`,
`EFF_ANA_IL1`, `EFF_UST_IL23`, `EFF_CSA_Th17`, `EFF_PRED`) — six
pharmacologically unrelated mechanisms (two anti-TNF biologics, IL-1 axis
blockade, IL-23 blockade, calcineurin-mediated Th17 suppression, and
glucocorticoid-receptor-mediated NF-kB suppression) sharing one ceiling.
This is the same defect class as `venous-thromboembolism/vte_refactor_notes.md`'s
`EMAX_RIV`/`EMAX_WARF` finding (`translations/UPSTREAM_ISSUES.md` #84) —
logged here as a new entry (`translations/UPSTREAM_ISSUES.md`, "shared
`Emax_drug`" entry, appended after re-checking the file's tail immediately
before writing; the entry may carry a number that duplicates another
agent's concurrent entry for a different file, per the file's own
top-of-file note that entry numbers are "not meaningful beyond 'later
than'" under concurrent writes).

**Not fixed in the original**, per the never-edit-upstream rule. The
refactored sibling declares six independent `EMAX_<STEM>` parameters, each
carrying the identical `0.92` default — so every `EFFECT_<STEM>` is
numerically unchanged (no scenario in the file overrides `Emax_drug`'s
default; confirmed by the exact-match verification below), while each
compound becomes independently driveable per the fork convention (e.g. a
virtual-population override intended to lower only Adalimumab's ceiling
no longer silently changes the other five compounds' ceilings too).

## Build-compat fix required to compile this file at all

The untouched original does not compile under mrgsolve 2.0.1 — confirmed
via the qspserver `mrgsolve_api` container's `/model_manifest`:

```
Error: improper annotation format
 input: kout_Th17    : 0.05
 context: parse annotated parameter block (PARAM)
```

Twelve `$PARAM @annotated` lines are missing the required third
(description) field: `kout_Th17`, `kin_Treg`, `kout_Treg`, `kin_M1`,
`kout_M1`, `kout_ROS`, `kout_CRP`, `k_Calp_syn`, `kout_Calp`,
`IC50_IFX_TNF`, `IC50_UST_IL23`, `Emax_drug`. Logged as
`translations/UPSTREAM_ISSUES.md` "twelve `$PARAM @annotated` lines
missing a description field" entry.

Per the guide's settled policy ("When the original doesn't compile at
all"), a syntax-only, non-numeric fix was applied directly to the
delivered `pg_mrgsolve_model_refactored.R`:

- The three lines inside compounds already being renamed
  (`IC50_IFX_TNF`, `IC50_UST_IL23`, `Emax_drug`) simply gained a real
  description as part of their normal rename to
  `EC50_IFX`/`EC50_UST`/`EMAX_<STEM>` — no separate fix needed.
- The nine untouched disease-side lines (`kout_Th17`, `kin_Treg`,
  `kout_Treg`, `kin_M1`, `kout_M1`, `kout_ROS`, `kout_CRP`, `k_Calp_syn`,
  `kout_Calp`) each had a short, purely descriptive third field appended
  (e.g. `kout_Th17    : 0.05    : Th17 clearance (1/d)`), changing
  nothing about the parameter's name, value, or use in any equation.

None of this touches a single numeric value, parameter, compartment
initial condition, or equation — confirmed by the exact-match
verification below, which is the actual proof. The untouched original
`pg_mrgsolve_model.R` still carries the defect forward unfixed, per the
never-edit-upstream rule.

## Verification

**Method.** Both files' embedded mrgsolve DSL blocks (`pg_qsp_code <-
'...'`) were mechanically extracted (regex on the assignment, verbatim
quoted text) into scratch `.cpp` files; the original's copy was
additionally patched in-memory with the twelve-line build-compat fix
above (never written back to the tracked `pg_mrgsolve_model.R`). Both
were POSTed to the local qspserver `mrgsolve_api` container at
`http://localhost:8007` (`/model_manifest` then `/run_simulation`),
which compiles and runs each DSL block directly with mrgsolve 2.0.1
server-side — no local R/mrgsolve install used. Requests were spaced
~2s apart per the API's `max_concurrent_jobs: 2` limit; one transient
`worker_killed` / `corrupted size vs. prev_size` container crash occurred
mid-run (an infrastructure hiccup unrelated to model content — confirmed
by a clean immediate retry succeeding), recorded here for transparency.

**Results: exact match on all 8 of the file's own dosing scenarios**, run
for their own full stated duration (`end = 168 d`, `delta = 0.5 d`, 337
timepoints — no shortening needed, well within the API's default
`maxsteps` budget):

1. **Prednisone_SOC** (60→40→20→10 mg/d oral taper): exact match, max
   abs/rel deviation **0.0** across all 41 shared `$CAPTURE`/compartment
   outputs.
2. **Cyclosporine** (140 mg BID oral x6mo): exact match, max abs/rel
   deviation **0.0** across all 41 shared outputs.
3. **Infliximab** (350 mg IV at 0/2/6wk then q8w): exact match, max
   abs/rel deviation **0.0** across all 41 shared outputs.
4. **Adalimumab** (80 mg SC load then 40 mg q1w): exact match, max
   abs/rel deviation **0.0** across all 41 shared outputs.
5. **Anakinra** (100 mg/d SC): exact match, max abs/rel deviation
   **0.0** across all 41 shared outputs.
6. **Ustekinumab** (90 mg SC at 0/4/16wk): exact match, max abs/rel
   deviation **0.0** across all 41 shared outputs.
7. **Combo_PRED_CSA** (Prednisone taper + Cyclosporine maintenance,
   tests combined TNF/IL-1/IL-6 suppression logic across two
   simultaneously-dosed compounds): exact match, max abs/rel deviation
   **0.0** across all 41 shared outputs.
8. **Combo_IFX_low_CsA** (Infliximab + low-dose Cyclosporine): exact
   match, max abs/rel deviation **0.0** across all 41 shared outputs.

(Scenario 9, "NoTreatment," is trivially identical by construction — no
dosing at all — and was not separately POSTed.) The 41 shared outputs
compared per scenario: 22 identically-named disease-side compartments/
captures (`TNFa`, `IL1b`, `IL17A`, `IL6`, `IL23`, `IL8`, `Neutroph`,
`NET`, `Th17`, `Treg`, `M1`, `MMP9_act`, `ROS_les`, `Ulcer`, `Healed`,
`CRP`, `Calprot`, `Pain_VAS`, `PARACELSUS`, `DLQI`, `HiSCRpseudo`,
`CompleteHeal`), 6 renamed concentration captures (`CONC_<X>` vs
`C_<X>`), and 13 renamed PK compartments (`ADA_SC` vs `GUT_ADA`, etc.),
matched by identical `$CMT` declaration order (only compartment *names*
changed, never their order, so indices 1-13 line up 1:1). This is the
expected, and achieved, result for six compounds that are all pure
structural renames (Archetypes 2/3/3-minus-peripheral, no Hill-fitting)
per the guide's tolerance table: "expect a near-exact match... anything
beyond floating-point-scale deviation means a bug" — here the match is
exact, not merely near-exact.

## qspserver `/model_manifest` discoverability

Confirmed via `POST /model_manifest` against the extracted, renamed DSL
(86 parameters, 47 outputs total): every PK parameter (`KA_ADA`, `F_ADA`,
`CL_ADA`, `V1_ADA`, `Q_ADA`, `V2_ADA`, `CL_IFX`, `V1_IFX`, `Q_IFX`,
`V2_IFX`, `KA_ANA`, `F_ANA`, `CL_ANA`, `V1_ANA`, `KA_CSA`, `F_CSA`,
`CL_CSA`, `V1_CSA`, `KA_UST`, `F_UST`, `CL_UST`, `V1_UST`, `KA_PRED`,
`CL_PRED`, `V1_PRED`) and every Hill parameter (`EC50_<STEM>`,
`EMAX_<STEM>`, `GAMMA_<STEM>` for all six stems) appears in the
manifest's `parameters` with its original numeric default. `C_ADA`,
`C_IFX`, `C_ANA`, `C_CSA`, `C_UST`, `C_PRED`, `EFFECT_ADA`, `EFFECT_IFX`,
`EFFECT_ANA`, `EFFECT_UST`, `EFFECT_CSA`, `EFFECT_PRED` are state-derived
(computed in `$ODE` from compartment concentrations), so — per the same
reasoning established for `EFFECT_TOCI`/`EFFECT_ANA`/etc. in prior
refactors in this fork (e.g. `cytokine-release-syndrome/crs_refactor_notes.md`,
`familial-mediterranean-fever/fmf_refactor_notes.md`) — they cannot also
be `$PARAM` entries; all twelve appear in the manifest's `outputPaths`
via the renamed `$CAPTURE @annotated` list, confirmed discoverable.

No separate `.cpp` extraction file was left behind — extraction was
in-memory/scratch only, used to build the verification requests above
and then discarded, per the workflow guide.

## Anything else flagged

- No compound outside the six was touched — this file models exactly
  Adalimumab, Anakinra, Cyclosporine, Infliximab, Prednisone, and
  Ustekinumab, all six in scope.
- `SUP_TNF` is computed but never read anywhere else in the original's
  `$ODE` (`dxdt_TNFa` recomputes the same product inline instead of
  reading it) — dead code, found and preserved verbatim (renamed
  references only), not a defect worth a build-compat fix since it
  compiles and runs fine as dead code.
- The R-side `simulate_pg()` scenario builder, `make_vpop()`, and
  `run_vpop_scenario()` were updated only where they name a renamed
  compartment (`"ADA_SC"` -> `"GUT_ADA"`, `"IFX_CENT"` -> `"CENT_IFX"`,
  `"ANA_SC"` -> `"GUT_ANA"`, `"CSA_GUT"` -> `"GUT_CSA"`, `"UST_SC"` ->
  `"GUT_UST"`, `"PRED_GUT"` -> `"GUT_PRED"`) — same dosing amounts, same
  timing throughout. `make_vpop()`'s sampled parameter names (`CL_ADA`,
  `CL_IFX`, `CL_CSA`, `UlcerArea_0`, `DSEV`, `F_pathergy`,
  `F_comorb_IBD`) were already spelled identically in both files, so no
  change was needed there.

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`, the
six `pyoderma-gangrenosum | Adalimumab/Anakinra/Cyclosporine (CSA)/
Infliximab (IFX)/Prednisone/Ustekinumab` rows.
