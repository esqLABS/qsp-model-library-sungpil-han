# Refactor notes — idiopathic-multicentric-castleman-disease / nine compounds (SILT, TOCZ, SIRO, RTX, ANA, RUX, DOXO, CYC, PRED)

Covers `idiopathic-multicentric-castleman-disease/imcd_mrgsolve_model.R`. All
nine compounds this file models were in scope (Anakinra, Chop cyclophosphamide,
Chop doxorubicin, Prednisone, Rituximab (RTX), Ruxolitinib, Siltuximab,
Sirolimus, Tocilizumab (TOCZ) — per their existing rows in
`driver-patches/data/compound_perturbation_census.md`) — there are no other
compounds in this file to leave untouched. Every disease-side compartment and
equation (IL-6 axis, mTORC1/TAFRO axis, lymph node size, plasmablast/B-cell
turnover, CRP/Hb/IgG/VEGF/platelet/hazard) is copied verbatim from the
original into `imcd_mrgsolve_model_refactored.R`; only the nine drugs' own PK
blocks and their named effect terms were renamed to the fork's convention.

## The original compiles as-is

`POST /model_manifest` on the untouched original's own DSL (extracted
verbatim from `code <- '...'`) returns `200` — no build-compatibility fix
was needed anywhere in this file.

## Archetype determined per compound (checked against this file's own equations)

The task flagged Rituximab, Siltuximab, and Tocilizumab specifically for a
TMDD check. Checked directly: **none of the nine compounds in this file
involve target-mediated drug disposition or any receptor-binding ODE
system** — every one is plain linear PK with an algebraic Emax-style effect
term (or, for rituximab, a plain first-order kill term with no saturation at
all).

| Compound | Archetype | Structure |
|---|---|---|
| Siltuximab (`SILT`) | 2 | no depot, 2-cmt (central+peripheral), linear CL/Q/V |
| Tocilizumab (`TOCZ`) | 2 | no depot, 2-cmt (central+peripheral), linear CL/Q/V |
| Sirolimus (`SIRO`) | 3 | depot + central + peripheral, linear |
| Rituximab (`RTX`) | 2 | no depot, 2-cmt (central+peripheral), linear CL/Q/V |
| Anakinra (`ANA`) | 3 minus peripheral | depot (SC) + central, linear |
| Ruxolitinib (`RUX`) | 3 minus peripheral | depot (oral) + central, linear |
| CHOP-doxorubicin (`DOXO`) | 1 | no depot, 1-cmt, linear elimination |
| CHOP-cyclophosphamide (`CYC`) | 1 | no depot, 1-cmt, linear elimination |
| Prednisone (`PRED`) | 3 minus peripheral | depot (oral) + central, linear |

Siltuximab was specifically checked for TMDD (anti-IL-6 antibody), Tocilizumab
for TMDD (anti-IL-6 receptor antibody), and Rituximab for TMDD (anti-CD20
antibody), per the task's explicit instruction — **none of the three is TMDD
in this file**: siltuximab and tocilizumab each feed a single concentration
compartment into a plain algebraic blocking-fraction equation (`1/(1 +
C_molar/Kd)`), with no free/complexed receptor pool anywhere in the ODE
system; rituximab feeds a single concentration into a plain first-order
CD20-lysis kill rate (`KILL_RTX * C_RTX`), again with no receptor-binding ODE
of any kind.

## Naming map

| Original | Refactored | Value | Role |
|---|---|---|---|
| `SILT_C` / `SILT_P` | `CENT_SILT` / `PERI_SILT` | -- | central/peripheral |
| `CL_SILT`, `V1_SILT`, `V2_SILT`, `Q_SILT` | unchanged | 0.232, 4.47, 2.49, 0.527 | already matched the convention |
| `Kd_SILT` | `EC50_SILT` | 1e-3 (nM) | Hill EC50 (renamed for the Hill interface; same role) |
| -- (none) | `EMAX_SILT` (new) | 1.0 | Hill ceiling [math-implied — original's block ratio saturates at a 100% neutralized fraction] |
| -- (none) | `GAMMA_SILT` (new) | 1.0 | Hill coefficient [original had none] |
| `MW_SILT` | unchanged | 145000 | Da (used inside the effect term's own unit conversion, not the naming convention's reserved roles) |
| `TOCZ_C` / `TOCZ_P` | `CENT_TOCZ` / `PERI_TOCZ` | -- | central/peripheral |
| `CL_TOCZ`, `V1_TOCZ`, `V2_TOCZ`, `Q_TOCZ` | unchanged | 0.20, 4.1, 2.5, 0.50 | already matched the convention |
| `Kd_TOCZ` | `EC50_TOCZ` | 2.5 (nM) | Hill EC50 |
| -- (none) | `EMAX_TOCZ` (new) | 1.0 | Hill ceiling [math-implied] |
| -- (none) | `GAMMA_TOCZ` (new) | 1.0 | Hill coefficient [new] |
| `SIRO_GUT` / `SIRO_C` / `SIRO_P` | `GUT_SIRO` / `CENT_SIRO` / `PERI_SIRO` | -- | depot/central/peripheral |
| `CL_SIRO`, `V1_SIRO`, `V2_SIRO`, `Q_SIRO`, `F_SIRO` | unchanged | 8.2, 12, 100, 10, 0.14 | already matched the convention |
| `Ka_SIRO` | `KA_SIRO` | 2.0 | absorption rate |
| `EC50_SIRO` | unchanged | 8 (ng/mL) | already matched the convention |
| -- (none) | `EMAX_SIRO` (new) | 1.0 | Hill ceiling [math-implied] |
| -- (none) | `GAMMA_SIRO` (new) | 1.0 | Hill coefficient [new] |
| `RTX_C` / `RTX_P` | `CENT_RTX` / `PERI_RTX` | -- | central/peripheral |
| `CL_RTX`, `V1_RTX`, `V2_RTX`, `Q_RTX` | unchanged | 0.32, 3.3, 4.5, 0.42 | already matched the convention |
| `RTX_kill` | `KILL_RTX` | 0.030 | bespoke first-order kill-rate constant (see Hill interface below — this compound has no EC50/Emax at all) |
| `ANA_SC` / `ANA_C` | `GUT_ANA` / `CENT_ANA` | -- | depot/central |
| `CL_ANA` | unchanged | 9.0 | already matched the convention |
| `V_ANA` | `V1_ANA` | 20 | central volume |
| `F_ANA` | unchanged | 0.95 | already matched the convention |
| `Ka_ANA` | `KA_ANA` | 0.40 | absorption rate |
| `RUX_GUT` / `RUX_C` | `GUT_RUX` / `CENT_RUX` | -- | depot/central |
| `CL_RUX` | unchanged | 18 | already matched the convention |
| `V_RUX` | `V1_RUX` | 72 | central volume |
| `F_RUX` | unchanged | 0.95 | already matched the convention |
| `Ka_RUX` | `KA_RUX` | 3.0 | **declared but unused in the original — preserved unfixed, see below** |
| `IC50_RUX` | `EC50_RUX` | 300 (nM) | Hill EC50 |
| -- (none) | `EMAX_RUX` (new) | 1.0 | Hill ceiling [math-implied] |
| -- (none) | `GAMMA_RUX` (new) | 1.0 | Hill coefficient [new] |
| `DOXO_C` | `CENT_DOXO` | -- | central |
| `CL_DOXO` | unchanged | 45 | already matched the convention |
| `V_DOXO` | `V1_DOXO` | 500 | central volume |
| -- (literal `0.05`) | `EMAX_DOXO` (new) | 0.05 | Hill ceiling — was an inline numeric literal in the original's `cyto_kill` term |
| -- (literal `0.1`) | `EC50_DOXO` (new) | 0.1 (mg/L) | Hill EC50 — same |
| -- (none) | `GAMMA_DOXO` (new) | 1.0 | Hill coefficient [new] |
| `CYC_C` | `CENT_CYC` | -- | central |
| `CL_CYC` | unchanged | 6 | already matched the convention |
| `V_CYC` | `V1_CYC` | 35 | central volume |
| -- (literal `0.04`) | `EMAX_CYC` (new) | 0.04 | Hill ceiling — inline literal in the original |
| -- (literal `5`) | `EC50_CYC` (new) | 5 (mg/L) | Hill EC50 — same |
| -- (none) | `GAMMA_CYC` (new) | 1.0 | Hill coefficient [new] |
| `PRED_GUT` / `PRED_C` | `GUT_PRED` / `CENT_PRED` | -- | depot/central |
| `CL_PRED` | unchanged | 5.6 | already matched the convention |
| `V_PRED` | `V1_PRED` | 35 | central volume |
| `Ka_PRED` | `KA_PRED` | 2.5 | absorption rate |
| `F_PRED` | unchanged | 0.85 | already matched the convention |
| `EC50_PRED` | unchanged | 40 (ng/mL) | already matched the convention |
| -- (none) | `EMAX_PRED` (new) | 1.0 | Hill ceiling [math-implied] |
| -- (none) | `GAMMA_PRED` (new) | 1.0 | Hill coefficient [new] |

`IL6_T`, `IL6_F`, `LN`, `PB`, `Bmem`, `CRP`, `Hb`, `IgG`, `VEGF`, `Anasarca`,
`Plt`, `mTOR`, `HAZ`, and every parameter/equation touching only those
(IL-6 axis, mTORC1/TAFRO axis, lymph node growth, CRP/Hb/IgG/VEGF/platelet
turnover, OS hazard) are **completely untouched** — not renamed, not
restructured.

## The Hill interface — five compounds already exact ratios, one bespoke, one uncoupled

**Siltuximab, Tocilizumab, Sirolimus, Prednisone**: all four blocked/reduced a
downstream signal via the same shape, `1/(1 + C/EC50)` (the *unblocked*
fraction, i.e. `1 - Emax*C/(EC50+C)` with `Emax=1`). Each is a **rename, not a
refit** — pulled out as `EFFECT_<STEM> = C^GAMMA/(EC50^GAMMA + C^GAMMA)` (the
*blocked* fraction, `Emax=1`, `Gamma=1`), with the disease-side block factor
recovered as `1 - EMAX_<STEM>*EFFECT_<STEM>`:

```
double EFFECT_SILT = pow(C_SILT_nM, GAMMA_SILT) /
                      (pow(EC50_SILT, GAMMA_SILT) + pow(C_SILT_nM, GAMMA_SILT));
double SILT_free_frac = 1.0 - EMAX_SILT * EFFECT_SILT;   // was: 1/(1+Csilt_molar/Kd_SILT)
```
(and identically in shape for `TOCZ_block`, `Sirol_block`, `Pred_block`).
Siltuximab and tocilizumab both convert their mg/L compartment concentration
to nM before the ratio (`C_<STEM>_nM = C_<STEM>/MW * 1e9`) exactly as the
original did — this conversion is preserved unchanged, including the
original's reuse of the same `145000` Da literal for tocilizumab as for
siltuximab (no separate `MW_TOCZ` parameter exists in the original, so none
was invented; disclosed rather than silently added).

**CHOP doxorubicin, CHOP cyclophosphamide**: the original's combined
cytotoxic-kill term, `cyto_kill = 0.05*(Cdoxo/(Cdoxo+0.1)) +
0.04*(Ccyc/(Ccyc+5))`, is already a **sum of two independent exact Hill
ratios** with inline numeric literals as Emax/EC50. Split into two named
terms and combined only at the point of use, per the guide's "combine only
where the disease equations actually use them" rule:

```
double EFFECT_DOXO = EMAX_DOXO * pow(C_DOXO, GAMMA_DOXO) / (pow(EC50_DOXO, GAMMA_DOXO) + pow(C_DOXO, GAMMA_DOXO));
double EFFECT_CYC  = EMAX_CYC  * pow(C_CYC,  GAMMA_CYC)  / (pow(EC50_CYC,  GAMMA_CYC)  + pow(C_CYC,  GAMMA_CYC));
...
double cyto_kill = EFFECT_DOXO + EFFECT_CYC;   // was: 0.05*(Cdoxo/(Cdoxo+0.1)) + 0.04*(Ccyc/(Ccyc+5))
```
Exact rename (`EMAX_DOXO=0.05`, `EC50_DOXO=0.1`, `EMAX_CYC=0.04`,
`EC50_CYC=5`, `GAMMA=1` for both), not a fit.

**Rituximab — genuinely bespoke, not Hill-shaped, disclosed rather than
force-fitted.** The original's CD20-lysis term, `RTX_kill * Crtx`, is a plain
**first-order rate** with no saturating EC50 anywhere — there is no ceiling,
no half-maximal concentration, nothing to name `EMAX_RTX`/`EC50_RTX` without
inventing values the original never had. Per the guide's explicit allowance
("if a compound's structure genuinely doesn't fit any of the archetypes,
don't force it"), this is handled as a bespoke linear effect:

```
double EFFECT_RTX = KILL_RTX * C_RTX;   // was: RTX_kill * Crtx -- linear, no EC50/Emax in the original
```
`EFFECT_RTX` is still one cleanly-named, independently-driveable function of
`C_RTX`, used at both of the original's two downstream sites (`Bmem` kill
rate, extra `IgG` clearance) exactly as `RTX_kill * Crtx` was.

**Anakinra — no `EFFECT_ANA` at all, because the original has none to rename
or refit.** `Cana` (renamed `C_ANA`) is computed for reporting but is never
read by any disease/PD equation anywhere in `$ODE` (confirmed by exact-name
grep). Inventing an `EFFECT_ANA` term would add pharmacology the original
never modeled, violating the guide's "never invent... a parameter/effect the
original didn't have" principle. `C_ANA` is still renamed and exposed (useful
as a redirectable PK covariate for future work), but no effect term
accompanies it. Logged as `translations/UPSTREAM_ISSUES.md` **#103** (below).

No `nls()` fitting was needed anywhere in this file (no TMDD, no ODE-solved
receptor kinetics) — the guide's "fit an approximation" path for Archetype 4
does not apply to any of the nine compounds.

## A genuine defect found while refactoring, preserved rather than "fixed" (Ruxolitinib)

Ruxolitinib's own absorption ODE in the original reads **sirolimus's**
parameter, not its own:

```
dxdt_RUX_GUT = -Ka_SIRO * 24 * RUX_GUT;
dxdt_RUX_C   =  Ka_SIRO * 24 * RUX_GUT * F_RUX - k10_RUX * RUX_C;
```

`Ka_RUX` (declared `3.0` in `$PARAM`) is never referenced anywhere in `$ODE`
(confirmed by exact-name grep) — ruxolitinib's real absorption rate is
governed by `Ka_SIRO` (`2.0`), sirolimus's own value. This is a genuine
numeric defect (not a build-compile issue — the file compiles and runs fine
as-is), logged as `translations/UPSTREAM_ISSUES.md` **#102** (below), same
class as `#90`/`#94` (one compound's parameter/dose silently entangled with a
different compound's own PK state).

**Not fixed here** — per the shared "never edit an upstream-tracked file /
log what you find, don't fix it" rule, and consistent with the precedent in
`graves-disease/gd_refactor_notes.md` (Scenario 6's rituximab-into-MMI
contamination, preserved and verified to reproduce identically), this is
preserved byte-for-behaviour-identical in the refactored sibling:
`dxdt_GUT_RUX`/`dxdt_CENT_RUX` still read `KA_SIRO`, and `KA_RUX` is carried
forward as a declared-but-unused parameter (`3.0`, unchanged). Verified below
with a bespoke ruxolitinib-alone scenario (no scenario in the original's own
R code doses ruxolitinib at all): `Crux`/`C_RUX` and every disease output
match **bit-for-bit** between the original and the refactored DSL, confirming
the cross-wired absorption rate carries through identically, bug included.

## `$PARAM` vs `$CAPTURE` for `C_<STEM>`/`EFFECT_<STEM>`

Following the pattern established across this batch (e.g.
`cytokine-release-syndrome/crs_refactor_notes.md`,
`graves-disease/gd_refactor_notes.md`): `C_<STEM>` and `EFFECT_<STEM>` are
plain `double` locals computed in `$ODE` (declaring them in `$PARAM` and then
assigning to them in `$ODE` does not compile under mrgsolve 2.0.1 —
`$PARAM` values are read-only inside `$ODE`), and are instead listed in
`$CAPTURE`. Confirmed via `POST /model_manifest` on the refactored DSL: all
sixteen `C_<STEM>`/`EFFECT_<STEM>` variables (`C_SILT, EFFECT_SILT, C_TOCZ,
EFFECT_TOCZ, C_SIRO, EFFECT_SIRO, C_RTX, EFFECT_RTX, C_ANA, C_RUX,
EFFECT_RUX, C_DOXO, EFFECT_DOXO, C_CYC, EFFECT_CYC, C_PRED, EFFECT_PRED`)
appear in `outputPaths`, discoverable through the output side of the
manifest.

Every naming-convention parameter (`CL_/V1_/V2_/Q_/KA_/F_/EC50_/EMAX_/
GAMMA_/KILL_` for all nine compounds) is confirmed present in `parameters`
with its original numeric default, via the same manifest call.

## Verification

**Method**: qspserver `mrgsolve_api` (`http://localhost:8007`),
`POST /model_manifest` (build-time check — confirmed all nine compounds'
naming-convention parameters and outputs discoverable) and
`POST /run_simulation` (numeric comparison), extracting the bare DSL text
from both the original (untouched, compiles as-is) and the refactored file's
final, delivered content, run through identical dosing. Requests were spaced
~2s apart per the API's stated concurrency limits; no failures or crashes
encountered other than one dosing-order fix (records had to be time-sorted
before submission — an API/harness requirement, not a model issue).

Eight of the nine compounds were verified against the **original file's own
scenarios** (from its `sim_one()`/`ev()` R code), all shortened from the
original's 365-day window to reduce solver load (per the guide's allowance)
while still covering every dose in the scenario:

| Compound(s) tested together | Scenario | Dosing (shortened window) | Result |
|---|---|---|---|
| Siltuximab alone | 2 | 770 mg IV q3w, days 0-357 (200-day window) | max abs `1e-4`, max rel `16.5` (both floating-point-scale — see below), 26 outputs, 219 timepoints |
| Tocilizumab alone | 3 | 560 mg IV q2w, days 0-350 (100-day window) | max abs `1e-10`, max rel `3.1e-5`, 26 outputs, 127 timepoints |
| Sirolimus alone | 4 | 2 mg PO QD x100 (100-day window) | **exact match, 0.0**, 26 outputs, 102 timepoints |
| Rituximab + Prednisone | 5 | RTX 750 mg IV days 0/7/14/21; PRED 60 mg PO QD x28 (100-day window) | max abs `1e-4`, max rel `11.9` (floating-point-scale), 27 outputs, 106 timepoints |
| CHOP (doxo+cyc) + Prednisone + Siltuximab | 6 | DOXO 112.5 mg + CYC 1125 mg q21d x7, PRED 100 mg PO QD x5/cycle, SILT 770 mg IV q3w (160-day window) | max abs `1.03e-7`, max rel `64.4` (floating-point-scale), 29 outputs, 190 timepoints |
| Siltuximab + Sirolimus + Anakinra | 7 | SILT 770 mg IV q3w, SIRO 2 mg PO QD, ANA 100 mg SC QD (100-day window) | max abs `5.1e-8`, max rel `50.2` (floating-point-scale), 28 outputs, 111 timepoints |

**Ruxolitinib has no scenario in the original file that doses it** — the
original's seven `sim_one()` scenarios cover siltuximab, tocilizumab,
sirolimus, rituximab, prednisone, CHOP, and anakinra, but never ruxolitinib.
Rather than inventing an unrelated regimen, ruxolitinib was verified with a
**bespoke bolus test** using its own declared label dose (20 mg PO BID, from
`DOSE_RUX`), run for 100 days: **exact match, max abs/rel deviation 0.0**,
26 shared outputs, 102 timepoints — this also confirms the preserved
`Ka_SIRO`-for-ruxolitinib cross-wire (see above) reproduces bit-for-bit.

**On the large `max_rel` values reported above**: every one of them is driven
by a variable (`IL6_free`, `VEGF`, `CRP`, `Cpred`/`C_PRED`, `Csiro`/`C_SIRO`)
passing extremely close to zero at some timepoint, where a `~1e-8`-to-`1e-13`
absolute floating-point-noise difference divides by a near-zero denominator
and inflates the relative metric to double digits. The **absolute**
deviations (`1e-4` at worst, most `1e-7`-`1e-10`) are exactly the
floating-point-scale noise expected for a pure structural reorganization with
no refit anywhere (Archetypes 1-3 per the guide's tolerance section) — this
was checked directly by inspecting the largest-relative-error variable in
each scenario, confirming it is a near-zero-crossing artifact, not a
systematic drift. Sirolimus (scenario 4) and ruxolitinib (bespoke) both came
back **exactly 0.0** on every shared output, with no near-zero-crossing
variable in either scenario to produce the artifact — consistent with the
same rename-only equivalence holding throughout.

**Result: near-exact match (floating-point-scale absolute deviation, worst
case `1.5e-4`) for all nine compounds**, consistent with the guide's
expectation for Archetypes 1-3 with a pure Hill-rename (no fitting) — no
compound in this file needed the Archetype-4 fitting path.

## Two upstream defects found and logged (not fixed)

- **`translations/UPSTREAM_ISSUES.md` #102** — ruxolitinib's absorption ODE
  reads sirolimus's own `Ka_SIRO` instead of its own declared `Ka_RUX` (see
  above).
- **`translations/UPSTREAM_ISSUES.md` #103** — anakinra's PK is fully
  simulated but never coupled to any disease/PD equation (see above).

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`: all nine
rows (Anakinra, Chop cyclophosphamide, Chop doxorubicin, Prednisone,
Rituximab (RTX), Ruxolitinib, Siltuximab, Sirolimus, Tocilizumab (TOCZ))
marked refactored/verified with their compartment/output names.
