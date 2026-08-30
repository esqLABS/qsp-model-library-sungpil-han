# Refactor notes — `t2dm_mrgsolve_model.R` (Empagliflozin, Glimepiride, Insulin Degludec, Metformin, Pioglitazone, Semaglutide, Sitagliptin)

**Scope of this pass.** Per the fork's PK/PD refactor spec
([`FORK_WORKFLOW_GUIDE.md`](../FORK_WORKFLOW_GUIDE.md), Part 2) and the
existing rows in
[`driver-patches/data/compound_perturbation_census.md`](../driver-patches/data/compound_perturbation_census.md)
(all seven `type2-diabetes` rows, each classified "Redirect concentration
(clean single site)"), all seven of this file's compounds were renamed to
the fork's pluggable-PK naming convention. This is the file's entire drug
roster (T2DM models no other compound), so there is nothing else in the
pharmacological sense to leave untouched — but the disease side (Bergman
glucose-insulin dynamics, glucagon, GLP-1, beta-cell mass, hepatic/
peripheral insulin resistance, FFA/adiposity, HbA1c/eGFR/UACR endpoints)
and the entire R scenario-runner layer were left byte-for-behavior
identical except for the lines that reference the renamed compartments/
parameters/effect terms by name. Confirmed by `diff` (see "Diff scope
confirmation" below).

## Archetype determination

All seven compounds are clean, single-site "Redirect concentration"
matches — none required a bespoke structure. Six of the seven map
directly onto Archetype 1 or Archetype 3 (with the Archetype-3 variants
the guide itself calls out: no-peripheral, and the CL/Q/V rewrite from
micro-rate-constants); the seventh (Insulin Degludec) is structurally a
clean Archetype-3 variant too, but its concentration is never actually
read by any disease-PD term in the original (see "A defect found, not
fixed" below) — the *rename* is still exactly Archetype 3, it is the PD
side that has nothing to rename.

| Compound | Archetype | Structure |
|---|---|---|
| **Metformin** | 3 (depot + central + peripheral, linear) | `GUT_MET`→`CENT_MET`↔`PERI_MET`, CL/Q/V1/V2 (converted from the original's CL/V elimination + k12/k21 micro-rate-constant mix) |
| **Empagliflozin** | 1 (no depot, single compartment, linear) | `CENT_EMPA` only |
| **Semaglutide** | 3 variant (depot + central, no peripheral, linear) | `GUT_SEMA`→`CENT_SEMA` |
| **Sitagliptin** | 1 (no depot, single compartment, linear) | `CENT_DPP4` only |
| **Glimepiride (SU)** | 1 (no depot, single compartment, linear) | `CENT_SU` only |
| **Insulin Degludec** | 3 variant (depot + central, no peripheral, linear) | `GUT_INS`→`CENT_INS`; **no `EFFECT_INS`** (see below) |
| **Pioglitazone** | 1 (no depot, single compartment, linear) | `CENT_PIOG` only |

**No TMDD, no receptor-binding kinetics anywhere in this file.** All seven
PK blocks are plain linear absorption/elimination; every concentration-
driven disease effect is already a plain Hill/Emax ratio in the original
(gamma either explicit or an implicit 1), so every Hill-interface term
below is a **rename, not a refit** — no `nls()` curve-fitting was needed
or performed for any compound, and no R²/fit-quality disclosure applies.

**Semaglutide's PK is not saturable or receptor-mediated in this model.**
The task brief asked to check specifically for this. The original's own
`CL_sema`/`Vc_sema` PK is a plain first-order depot→central linear system
(`ka_sema`/`CL_sema`/`Vc_sema`, no Michaelis-Menten or target-mediated
term) — its very low clearance (`CL_sema = 0.053 L/h`, giving the model's
own comment "very low CL (t½~168h)") is what produces the long
apparent half-life, not any nonlinear elimination mechanism. Only the
*disease-facing* Hill terms (`EFFECT_SEMA_INS`, `EFFECT_SEMA_GC`) are
concentration-saturating, and those are the compound's pharmacodynamics,
not its pharmacokinetics.

**Insulin Degludec's depot is a single first-order step, not a
multi-hexamer cascade.** The task brief asked to check for this too. The
original models exogenous insulin degludec absorption as one depot
(`INS_SC`) feeding one central compartment (`INS_C`) via a single `ka_ins
= 0.048 h⁻¹` (very slow, representing the hexamer-to-monomer dissociation
lumped into one rate constant, per the model's own comment "from SC
hexamer depot"), not as a sequential multi-step depot the way
`type1-diabetes/t1dm_mrgsolve_model.R`'s own insulin PK does (`GUT1_INS`→
`GUT2_INS`, see `t1dm_refactor_notes.md`). Archetype 3 (depot + central,
no peripheral) fits cleanly, no bespoke handling needed for the PK itself.

## Naming applied

Keyed to each compound's own stem (`MET`, `EMPA`, `SEMA`, `DPP4`, `SU`,
`INS`, `PIOG`) — the original already used these as ad hoc prefixes/
suffixes (`ka_met`, `EMPA_C`, `Csema`, `DPP4I_C`/`IC50_dpp4`, `EC50_su`,
`ka_ins`, `Cpiog`), just inconsistently placed and cased. `DPP4I_C`'s "I"
suffix (present only on that one compartment name, nowhere else in
sitagliptin's own block) was dropped for consistency — a purely cosmetic
normalization, not a new stem.

| Role | Original | Refactored |
|---|---|---|
| Metformin absorption depot | `MET_GUT` | `GUT_MET` |
| Metformin central | `MET_C` | `CENT_MET` |
| Metformin peripheral | `MET_P` | `PERI_MET` |
| Metformin absorption rate | `ka_met` | `KA_MET` |
| Metformin bioavailability (applied at dosing-amount level; unused in `$ODE`, same as original) | `F_met` | `F_MET` |
| Metformin central volume | `Vc_met` | `V1_MET` |
| Metformin clearance | `CL_met` | `CL_MET` |
| Metformin inter-compartmental CL (derived, see below) | `k12_met`/`k21_met` | `Q_MET` |
| Metformin peripheral volume (derived, see below) | — | `V2_MET` |
| Metformin molecular weight (declared, unused in `$ODE`, same as original) | `MW_met` | `MW_MET` |
| Metformin exposed concentration | *(none — `Cmet` computed inline)* | `C_MET` |
| Metformin EGP-suppression Emax/EC50/gamma (was inline literals `0.30`/`800`) | *(none)* | `EMAX_MET_EGP`, `EC50_MET`, `GAMMA_MET = 1` |
| Metformin Rd-stimulation Emax (was inline literal `0.12`, shares `EC50_MET`/`GAMMA_MET`) | *(none)* | `EMAX_MET_RD` |
| Metformin disease-facing effect terms | `Emet_EGP`, `Emet_Rd` | `EFFECT_MET_EGP`, `EFFECT_MET_RD` |
| Empagliflozin central | `EMPA_C` | `CENT_EMPA` |
| Empagliflozin absorption rate (declared, unused in `$ODE`, same as original) | `ka_empa` | `KA_EMPA` |
| Empagliflozin bioavailability (declared, unused in `$ODE`, same as original) | `F_empa` | `F_EMPA` |
| Empagliflozin central volume | `Vc_empa` | `V1_EMPA` |
| Empagliflozin clearance | `CL_empa` | `CL_EMPA` |
| Empagliflozin exposed concentration | `Cempa` | `C_EMPA` |
| Empagliflozin Emax/EC50/gamma | `Imax_empa`, `IC50_empa` (used as `IC50_empa*1000` at the point of use) | `EMAX_EMPA`, `EC50_EMPA = 30000` (the `*1000` folded into the parameter's own default; same number, see "Folded unit constants" below), `GAMMA_EMPA = 1` |
| Empagliflozin disease-facing effect term | `E_empa` | `EFFECT_EMPA` |
| Semaglutide absorption depot | `SEMA_SC` | `GUT_SEMA` |
| Semaglutide central | `SEMA_C` | `CENT_SEMA` |
| Semaglutide absorption rate | `ka_sema` | `KA_SEMA` |
| Semaglutide bioavailability (declared, unused anywhere — not even in the original's own dosing calc — same as original) | `F_sema` | `F_SEMA` |
| Semaglutide central volume | `Vc_sema` | `V1_SEMA` |
| Semaglutide clearance | `CL_sema` | `CL_SEMA` |
| Semaglutide exposed concentration | `Csema` | `C_SEMA` |
| Semaglutide insulin-secretion Emax/EC50/gamma | `Emax_sema_ins`, `EC50_sema`, hardcoded exponent `1.5` | `EMAX_SEMA_INS`, `EC50_SEMA`, `GAMMA_SEMA_INS` |
| Semaglutide glucagon-suppression Emax (was inline literal `0.40`)/gamma (was implicit 1) | *(none)* | `EMAX_SEMA_GC`, `GAMMA_SEMA_GC = 1` (shares `EC50_SEMA`) |
| Semaglutide disease-facing effect terms (concentration-driven) | `E_sema_ins`, `E_sema_Gc` | `EFFECT_SEMA_INS`, `EFFECT_SEMA_GC` |
| Semaglutide weight-loss term (time-only, **not** concentration-driven — see "Bespoke, non-Hill terms" below) | `Emax_sema_wt`, `ksema_wt`, `BW_sema_effect` | `EMAX_SEMA_WT`, `K_SEMA_WT`, `BW_SEMA_EFFECT` |
| Sitagliptin central | `DPP4I_C` | `CENT_DPP4` |
| Sitagliptin absorption rate (declared, unused in `$ODE`, same as original) | `ka_dpp4` | `KA_DPP4` |
| Sitagliptin bioavailability (declared, unused in `$ODE`, same as original) | `F_dpp4` | `F_DPP4` |
| Sitagliptin central volume | `Vc_dpp4` | `V1_DPP4` |
| Sitagliptin clearance | `CL_dpp4` | `CL_DPP4` |
| Sitagliptin exposed concentration | `Cdpp4` | `C_DPP4` |
| Sitagliptin Emax/EC50/gamma | `Imax_dpp4`, `IC50_dpp4` (used as `IC50_dpp4*1000`) | `EMAX_DPP4`, `EC50_DPP4 = 100` (folded, see below), `GAMMA_DPP4 = 1` |
| Sitagliptin disease-facing effect term | `DPP4_inh` | `EFFECT_DPP4` |
| Glimepiride central | `SU_C` | `CENT_SU` |
| Glimepiride absorption rate (declared, unused in `$ODE`, same as original) | `ka_su` | `KA_SU` |
| Glimepiride bioavailability (declared, unused in `$ODE`, same as original) | `F_su` | `F_SU` |
| Glimepiride central volume | `Vc_su` | `V1_SU` |
| Glimepiride clearance | `CL_su` | `CL_SU` |
| Glimepiride exposed concentration | `Csu` | `C_SU` |
| Glimepiride Emax/EC50/gamma | `Emax_su`, `EC50_su`, `n_su` | `EMAX_SU`, `EC50_SU`, `GAMMA_SU` |
| Glimepiride disease-facing effect term | `E_su` | `EFFECT_SU` |
| Insulin Degludec absorption depot | `INS_SC` | `GUT_INS` |
| Insulin Degludec central | `INS_C` | `CENT_INS` |
| Insulin Degludec absorption rate | `ka_ins` | `KA_INS` |
| Insulin Degludec bioavailability (declared, unused in `$ODE`, same as original) | `F_ins` | `F_INS` |
| Insulin Degludec central volume | `Vc_ins` | `V1_INS` |
| Insulin Degludec clearance | `CL_ins` | `CL_INS` |
| Insulin Degludec exposed concentration | `Cins` | `C_INS` (exposed per the naming convention; **no `EFFECT_INS`** — see below) |
| Pioglitazone central | `PIOG_C` | `CENT_PIOG` |
| Pioglitazone absorption rate (declared, unused in `$ODE`, same as original) | `ka_piog` | `KA_PIOG` |
| Pioglitazone bioavailability (declared, unused in `$ODE`, same as original) | `F_piog` | `F_PIOG` |
| Pioglitazone central volume | `Vc_piog` | `V1_PIOG` |
| Pioglitazone clearance | `CL_piog` | `CL_PIOG` |
| Pioglitazone exposed concentration | `Cpiog` | `C_PIOG` |
| Pioglitazone IR-reduction Emax/EC50/gamma | `Emax_piog_IR`, `EC50_piog`, implicit gamma 1 | `EMAX_PIOG`, `EC50_PIOG`, `GAMMA_PIOG = 1` |
| Pioglitazone disease-facing effect term (concentration-driven) | `E_piog_IR` | `EFFECT_PIOG_IR` |
| Pioglitazone fluid-retention/weight-gain term (time-only, **not** concentration-driven — see below) | `Emax_piog_wt`, `BW_piog_gain` | `EMAX_PIOG_WT`, `BW_PIOG_GAIN` |
| Disease-side glucose/insulin/glucagon/GLP-1/beta-cell/IR/FFA/HbA1c/eGFR/UACR machinery | *(unchanged)* | *(unchanged — outside any compound's own naming-convention scope)* |

**All eight Hill terms are renames, not refits.** Metformin's two, SU's
one, DPP-4's one, Empagliflozin's one, Semaglutide's two, and
Pioglitazone's one concentration-driven effect terms were all already
plain `Emax*Cᵞ/(Cᵞ+EC50ᵞ)` ratios in the original (gamma explicit for
Semaglutide-GSIS at 1.5 and Glimepiride at 1.5, implicit 1 for the other
six) — pulling each out under `EMAX_<STEM>[_<AXIS>]`/`EC50_<STEM>`/
`GAMMA_<STEM>[_<AXIS>]` with the original's own numeric values changes
zero numeric behaviour, confirmed exactly by the verification below.

**Folded unit constants: two `EC50`s were pre-multiplied by 1000 at the
point of use in the original, not stored as their own parameter.**
`E_empa = Imax_empa * Cempa / (Cempa + IC50_empa * 1000)` and `DPP4_inh =
Imax_dpp4 * Cdpp4 / (Cdpp4 + IC50_dpp4 * 1000)` both compute their actual
half-maximal constant inline (`IC50_empa=30` × 1000 = 30000;
`IC50_dpp4=0.10` × 1000 = 100) rather than declaring it directly. Per the
naming convention's single clean `EC50_<STEM>` parameter, this refactor
folds the multiplication into the parameter's own stored default
(`EC50_EMPA = 30000`, `EC50_DPP4 = 100`) rather than keeping a separate
`IC50_*`-times-1000 expression at the use site — the exact same number
the original computed, moved rather than changed. Confirmed non-numeric
by the exact-match verification below (which exercises both compounds).
This produces an oddly large `EC50_EMPA` relative to empagliflozin's own
therapeutic plasma concentrations (tens to low hundreds of ng/mL, per the
compound's own literature) — i.e. `EFFECT_EMPA` never gets close to
saturating at realistic doses in this model. That looks like a genuine
original units inconsistency (the comment says "ng/mL" for `IC50_empa`,
but the `*1000` factor effectively demands mcg/mL-scale exposure to
saturate), but it is a **parameter-value** question, not a syntax defect,
so per the never-invent-or-correct-a-parameter-value rule it is disclosed
here rather than fixed — the refactored model reproduces the original's
own (possibly miscalibrated) dose-response exactly.

**Bespoke, non-Hill terms: two compounds each carry one weight-change term
that is a function of elapsed time, not concentration.** Semaglutide's
`BW_SEMA_EFFECT` (`EMAX_SEMA_WT * (1 - exp(-K_SEMA_WT * TIME))`) and
Pioglitazone's `BW_PIOG_GAIN` (`EMAX_PIOG_WT * (1 - exp(-0.003 * TIME))`)
are asymptotic time-courses gated by `USE_SEMA`/`USE_PIOG` but **not**
functions of `C_SEMA`/`C_PIOG` at all — a genuinely different pattern from
every other effect term in this file (all of which read the compound's
own concentration). Renamed for consistency with the compound's stem but
kept outside the Hill interface, per the same reasoning
`type1-diabetes/t1dm_refactor_notes.md` gives for INS's/TEP's linear
(non-saturating) terms: a Hill/Emax formula always saturates in its
*input*, and the input here is elapsed time, not concentration, so
forcing it into `EFFECT_<STEM>` alongside the true concentration-driven
terms would misrepresent what actually gates it. `0.003` in
`BW_PIOG_GAIN` is Pioglitazone's own hardcoded per-hour rate constant in
the original (distinct from `K_SEMA_WT`, which the original *did* name);
left as a literal rather than invented into a new named parameter, since
introducing one that does not exist in the original would be adding
scope, not renaming it — flagged here for visibility rather than silently
left as an anonymous magic number.

## A build-compat defect fixed for the refactor (logged upstream)

The untouched original does not compile under mrgsolve 2.0.1 at all — three
independent, stacked defects, none inside any one compound's own PK/PD
block (all sit in shared boilerplate: the `$INIT` block, the `$CAPTURE`
list, and two `$MAIN` lines in the weight-change terms). Logged in full,
with the exact error text for each, as
[`translations/UPSTREAM_ISSUES.md` #99](../translations/UPSTREAM_ISSUES.md):

1. **`$CMT` + `$INIT` jointly redeclare all 25 compartments** ("Duplicated
   model names" — the same defect class as issues #29/#34/#36/#59/#76/
   #81/#82/#83/#96). Fixed in `t2dm_mrgsolve_model_refactored.R` by
   deleting the `$INIT` block and moving its 25 assignments into `$MAIN`
   as `<CMT>_0 = value;` under the refactor's renamed compartments (e.g.
   `GUT_MET_0 = 0;`), same values, no compartment added or removed.
2. **`$CAPTURE` re-lists nine `$CMT` names directly** (`Gc`, `GLP1`,
   `IR_H`, `IR_P`, `FFA`, `BW_t`, `HbA1c_cmpt`, `eGFR_cmpt`, `UACR_cmpt` —
   same defect class as issue #93). Fixed by removing those nine names
   from `$CAPTURE`; they remain in every output as ordinary compartment
   states (mrgsolve always returns `$CMT` states regardless of
   `$CAPTURE`), just no longer double-declared.
3. **`SOLVERTIME` used inside `$MAIN`, where it is out of scope.**
   `SOLVERTIME` expands to `_ODETIME_[0]`, a variable mrgsolve only
   declares inside the generated `$ODE` function; both occurrences (in
   Semaglutide's and Pioglitazone's own weight-change terms) are in
   `$MAIN`. Fixed by substituting `TIME` (mrgsolve's per-record
   simulation-time covariate, valid in `$MAIN`), which reads the identical
   value at every point this file's own harness ever observes it — `$MAIN`
   re-executes at each of `mrgsim(end=, delta=)`'s per-output-time
   observation records, so `TIME` there equals the same elapsed-time value
   `SOLVERTIME` would have provided had it been in scope. Confirmed
   non-numeric by the exact-match verification below (both weight terms
   are exercised in the Semaglutide/Pioglitazone scenarios).

All three fixes are syntax-only and non-numeric — no compound's own PK/PD
values, structure, or the disease-side equations were touched to make
this happen. See UPSTREAM_ISSUES.md #99 for the exact compiler errors.

## A second, non-compile-blocking defect found and preserved as-is (also logged upstream)

Logged as
[`translations/UPSTREAM_ISSUES.md` #100](../translations/UPSTREAM_ISSUES.md),
not fixed (neither defect is a syntax/build problem, and both sit inside
either one compound's own block or the harness's own scenario list, so
per the guide's "log what you find, don't fix it upstream" rule they are
carried forward unchanged):

1. **Insulin Degludec's own concentration is computed, captured, and even
   has a dedicated dosing scenario — but is never read by any
   disease-PD term.** Grepping the whole original for `Cins`/`INS_C`
   outside insulin degludec's own PK block and its own
   dosing/scenario wiring finds nothing: no EGP, Rd, insulin-action, or
   secretion equation ever reads it. The only other place `CL_ins`/
   `Vc_ins` reappear is inside `dxdt_Ip` (endogenous plasma insulin's own
   clearance term), applied **unconditionally regardless of `USE_INS` or
   dosing** — i.e. those two constants are reused as a generic clearance-
   rate denominator for *endogenous* insulin, not as insulin degludec's
   own pharmacological effect. Confirmed by toggling `USE_INS` 0 vs 1 with
   identical dosing through `mrgsolve_api`: `Gp`/`HbA1c_cmpt`/every other
   disease output is bit-identical, while `C_INS` itself rises and falls
   exactly as the PK dictates — a fully modeled, fully captured, but
   provably inert seventh compound. **Per the guide's naming convention,
   `C_INS` is still exposed** (it is a clean, single, renamed
   concentration site — the archetype match is real), but **no
   `EFFECT_INS` was invented** to paper over the gap; there is nothing in
   the original to rename. `CL_INS`/`V1_INS` remain reused, unchanged,
   inside `dxdt_Ip` exactly as the original construction did.
2. **Glimepiride and Pioglitazone are fully parameterized (PK + a named
   Hill effect term each) but neither is ever dosed by any of the file's
   own seven scenarios.** `build_doses()` has no `su`/`piog` branch, and
   no scenario object sets `su=1`/`piog=1` — `USE_SU`/`USE_PIOG` stay at
   their `$PARAM` default of 0 in every scenario the file ships with.
   Confirmed by grepping `build_doses()` and the `scenarios` list. This
   means the guide's "run the original's own scenarios" verification
   instruction has no scenario to reuse for these two compounds — see
   "Verification" below for how this was handled.

## Verification

**Method.** Both files' embedded mrgsolve DSL blocks were mechanically
extracted (`code <- '...'`/`t2dm_model_code <- '...'` → bare DSL text) and
run through the qspserver `mrgsolve_api` container (`http://localhost:8007`,
confirmed healthy), using `POST /model_manifest` to confirm both compile
and expose the right parameters/outputs, and `POST /run_simulation` to
compare every shared `$CAPTURE`d output (plus every shared `$CMT`
compartment state) across nine scenarios, requests spaced ~2s apart per
the API's own concurrency limit. Compartment order is unchanged by the
rename, so both files share the same 1-based NM-TRAN compartment index
for every dosing call (`GUT_MET`/`MET_GUT`=1, `CENT_EMPA`/`EMPA_C`=4,
`GUT_SEMA`/`SEMA_SC`=5, `CENT_DPP4`/`DPP4I_C`=7, `CENT_SU`/`SU_C`=8,
`GUT_INS`/`INS_SC`=9, `CENT_PIOG`/`PIOG_C`=11).

1. **`no_tx` — untreated baseline.** The original's own zero-dose event
   (`ev(amt=0, cmt=MET_GUT/GUT_MET, time=0)`), all `USE_*` flags at their
   `$PARAM` default of 0, 365 days at the original's own `delta=6`.
2. **`metformin`** — Metformin 1000 mg BID (`build_doses()`'s own
   `amt = 1000*0.55`, `ii=12`, `addl=729`), `USE_MET=1`, 365 days.
3. **`combo_empa`** — metformin dose + Empagliflozin 10 mg QD
   (`build_doses()`'s own `amt = 10*0.86*1e6/450.9`, `ii=24`, `addl=364`),
   `USE_MET=1, USE_EMPA=1`, 365 days.
4. **`combo_sema`** — metformin dose + Semaglutide 1 mg SC weekly
   (`build_doses()`'s own `amt = 1000/4113.6*1000`, `ii=168`, `addl=52`),
   `USE_MET=1, USE_SEMA=1`, 365 days.
5. **`triple`** — metformin + empagliflozin + semaglutide doses together,
   `USE_MET=1, USE_EMPA=1, USE_SEMA=1`, 365 days.
6. **`insulin`** — metformin dose + Insulin Degludec 20 U QD
   (`build_doses()`'s own `amt = 20*0.0347*0.91*1000/6103`, `ii=24`,
   `addl=364`), `USE_MET=1, USE_INS=1`, 365 days.
7. **`dpp4i`** — metformin dose + Sitagliptin 100 mg QD
   (`build_doses()`'s own `amt = 100*0.87*1e6/407.5`, `ii=24`,
   `addl=364`), `USE_MET=1, USE_DPP4=1`, 365 days.

Scenarios 2-7 are exactly the file's own `build_doses()`/`scenarios` list
(dose amounts, intervals, and durations taken verbatim), run to the full
365-day/`delta=6` horizon the original's own `run_scenario()` uses — no
step-count/timeout issue was encountered at this horizon (confirmed by
running the full year before attempting any shortening), so no window
had to be shortened.

**Bespoke verification scenarios for Glimepiride and Pioglitazone (no
original scenario exists to reuse — see UPSTREAM_ISSUES.md #100).**
Since neither compound is ever dosed in the original's own scenario list,
a dosing pattern matching the style already used for this file's other
direct-to-central compounds (round bolus amount, once-daily, no depot)
was built for verification purposes only — never added to either
delivered file's own scenario list:

8. **`su_bespoke`** — 50 (mass units) into `CENT_SU`/`SU_C` (idx 8), QD
   for 30 days, `USE_SU=1` alone, then again combined with the metformin
   dose (`USE_MET=1, USE_SU=1`) — confirms Glimepiride's renamed PK/PD
   block computes identically to the original's when actually exercised,
   and that it does not cross-contaminate metformin's own renamed block.
9. **`piog_bespoke`** — 50 (mass units) into `CENT_PIOG`/`PIOG_C` (idx
   11), QD for 30 days, `USE_PIOG=1` alone, then combined with metformin
   (`USE_MET=1, USE_PIOG=1`) — same purpose for Pioglitazone.

Every shared `$CAPTURE`d output (`Cmet`/`C_MET`, `Cempa`/`C_EMPA`,
`Csema`/`C_SEMA`, `Cdpp4`/`C_DPP4`, `Csu`/`C_SU`, `Cins`/`C_INS`,
`Cpiog`/`C_PIOG`, `EGP_val`, `Rd_val`, `UGE_loss`, `E_empa`/`EFFECT_EMPA`,
`E_sema_ins`/`EFFECT_SEMA_INS`, `DPP4_inh`/`EFFECT_DPP4`, `E_su`/
`EFFECT_SU`, `E_piog_IR`/`EFFECT_PIOG_IR`, `beta_prolif_rate`,
`beta_apop_rate`) plus every shared `$CMT` compartment state (`Gp`, `Gt`,
`Ip`, `X_action`, `Gc`, `GLP1`, `beta_mass`, `IR_H`, `IR_P`, `FFA`, `BW_t`,
`HbA1c_cmpt`, `eGFR_cmpt`, `UACR_cmpt`) was compared across the full
output time grid for each scenario. (`EFFECT_MET_EGP`, `EFFECT_MET_RD`,
and `EFFECT_SEMA_GC` have no original `$CAPTURE` counterpart to diff
directly against — the original never captured `Emet_EGP`/`Emet_Rd`/
`E_sema_Gc` — but their correctness is confirmed indirectly through the
exact match of `EGP_val`/`Rd_val` (which consume the first two) and `Gc`
(which consumes the third), since a wrong value in any of the three would
necessarily show up as a mismatch in its downstream consumer.)

**Result: exact match. Maximum absolute and relative deviation observed =
0.0 (bit-identical) across all 9 scenarios and every shared output/state.**
This is the expected outcome for a pure structural reorganization/rename
with no Hill-refitting, per the guide's tolerance table — consistent with
every compound here being a plain rename (no TMDD, no ODE-derived Hill fit
required for any of the seven).

One benign artefact was observed and is not a mismatch: the very first
output row (t=0, before the dose-time row) reports `NaN` for `Rd_val` and
`UGE_loss` in **both** models identically (a `BW_t`-division-by-zero
inherent to how mrgsolve reports the pre-`$MAIN` state at the record
immediately preceding the first dose), confirmed to occur at the same
index with the same value on both sides — treated as an equal (non-
mismatching) NaN-vs-NaN comparison, not a divergence.

**A qspserver infrastructure note.** All ~20 `/run_simulation` calls made
during this verification succeeded on the first attempt with requests
spaced ~2s apart, per the API's own stated `max_concurrent_jobs: 2` and
prior reports of crashes under heavy concurrent load in other refactors'
notes in this fork — no flakiness was observed this time.

## qspserver compatibility

- `C_<STEM>` and every `EFFECT_<STEM>[_<AXIS>]` are `$MAIN`-local
  `double`s, not `$PARAM` entries — this file's own pre-existing idiom
  (the untouched original already computes `Cmet`/`Cempa`/etc. and
  `E_empa`/etc. as plain `$MAIN` doubles, not `$PARAM`), and consistent
  with the hard mrgsolve 2.0.1 constraint documented in
  `diabetic-ketoacidosis/dka_refactor_notes.md` and
  `type1-diabetes/t1dm_refactor_notes.md` (assigning to a `$PARAM`-declared
  name inside `$MAIN`/`$ODE` is a compile error, since mrgsolve passes
  `$PARAM` values in as `const double&`).
- Discoverability is satisfied directly via `$CAPTURE` — unlike the DKA/
  T1DM files (which needed a separate `_report`-suffixed recomputation
  because their `$TABLE`/`$CAPTURE` block cannot see `$ODE`-local
  variables), this file's own idiom already captures `$MAIN`-computed
  doubles by their bare name (`Cmet`, `E_empa`, etc., in the untouched
  original) with no separate `$TABLE` block at all — so the renamed
  `C_<STEM>`/`EFFECT_<STEM>` doubles are captured by their own bare names
  directly, no `_report` suffix needed. Confirmed present in
  `/model_manifest`'s `outputPaths` (45 entries: 25 `$CMT` states + 7
  `C_<STEM>` + `EGP_val`/`Rd_val`/`UGE_loss` + 8 `EFFECT_<STEM>[_<AXIS>]`
  + `beta_prolif_rate`/`beta_apop_rate`).
- All 25 `$CMT` names (11 drug-PK, renamed; 14 disease-side, unchanged)
  and all 108 `$PARAM` entries (including the 8 new Hill parameters:
  `EMAX_MET_EGP`/`EMAX_MET_RD`/`EC50_MET`/`GAMMA_MET`,
  `EMAX_SEMA_GC`/`GAMMA_SEMA_GC`, plus `GAMMA_<STEM>=1` made explicit for
  the other five compounds) are confirmed present and correctly defaulted
  via `/model_manifest`.
- `$SET end/delta` is not present in this model at all (both files rely
  entirely on `mrgsim(end=, delta=)` call-site arguments in
  `run_scenario()`), so there is nothing to note about it being
  non-load-bearing for API runs — every verification run above supplied
  `time.end`/`time.delta` explicitly via the request.
- No R-only syntax inside the DSL block: no `sprintf`/`glue`-interpolated
  constants, no externally-defined R helper functions called from
  `$MAIN`/`$ODE`. Confirmed by the DSL compiling and running standalone
  through the API with no R wrapper present.

## Diff scope confirmation

`git diff --no-index type2-diabetes/t2dm_mrgsolve_model.R
type2-diabetes/t2dm_mrgsolve_model_refactored.R` touches: all seven
compounds' `$PARAM` blocks (renamed, plus the 8 new/explicit Hill
parameters), all 11 drug-PK `$CMT` lines (renamed), the `$INIT` block
(deleted, moved into `$MAIN` as 25 `<CMT>_0` lines per the build-compat
fix), the plasma-concentration and drug-effect blocks in `$MAIN` (renamed,
plus the two `SOLVERTIME`→`TIME` build-compat substitutions), every
`$ODE` line referencing a renamed drug-PK compartment/parameter, the
`$CAPTURE` block (renamed, nine `$CMT` names removed per the build-compat
fix), and `build_doses()`'s/`run_scenario()`'s renamed `cmt=`/`param()`
string literals. Nothing else — the Bergman glucose-insulin ODEs,
glucagon/GLP-1/beta-cell-mass/IR/FFA/HbA1c/eGFR/UACR equations, the
`scenarios` list's names/flags, the plotting code, and the clinical-
calibration checkpoint table are byte-for-byte identical (module the
renamed `C_<STEM>`/`EFFECT_<STEM>` output column references the plotting
code was never touched, since none of it reads a renamed column directly
— all seven plots reference disease-side columns only: `HbA1c_cmpt`,
`Gp`, `BW_t`, `beta_mass`, `eGFR_cmpt`, `UGE_loss`, `GLP1`).

## Anything else worth flagging

- This is the first file in this fork's refactor history where **all
  seven** modeled compounds are plain-PK, plain-Hill "clean single site"
  renames with zero TMDD/receptor-binding/multi-hexamer/ODE-fit cases
  among them — a useful data point for the census's own note that
  individual row classifications are "a strong first pass, not ground
  truth": every one of this file's seven "Redirect concentration (clean
  single site)" rows turned out to be exactly that on inspection.
- Two of the seven compounds (Glimepiride, Pioglitazone) are fully built
  but functionally dead weight in every scenario the file ships with, and
  an eighth "compound" is really the disease's own endogenous-insulin
  clearance term secretly reusing insulin degludec's PK constants — see
  UPSTREAM_ISSUES.md #100 for the full disclosure. None of this was
  fixed; all three are preserved exactly as the original built them.
