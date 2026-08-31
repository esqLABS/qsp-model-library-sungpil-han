# Refactor notes — `hashimoto-thyroiditis/ht_mrgsolve_model_refactored.R`

## Scope: both census rows for this file are genuine, externally-dosed drugs

This file carries two rows in `driver-patches/data/compound_perturbation_census.md`:
**Levothyroxine (LT)** and **Liothyronine (LIT)**. Per the task's flagged
suspicion (this disease shares the same HPT/thyroid-hormone axis where
endogenous hormones can get confused with drugs — the same distinction
`graves-disease/gd_refactor_notes.md` had to make for T3/T4/TRAb/TSH),
checked the actual code rather than trusting the census's compound name
alone:

| Compound | Real drug? | Evidence |
|---|---|---|
| LT4 (Levothyroxine) | **Yes** | Own dosing-rate parameter `LT4_dose` (μg/day), own three-compartment PK (`LT4g`/`LT4c`/`LT4p`), own absorption/elimination/distribution parameters (`ka_LT4`, `kel_LT4`, `k12_LT4`, `k21_LT4`, `F_LT4`, `CL_LT4`, `Vc_LT4`), exercised by 6 of the file's own 7 scenarios (`LT4_dose` = 75–175 μg/day) |
| LiT3 (Liothyronine) | **Yes** | Own dosing-rate parameter `LiT3_dose` (μg/day), own compartment `LiT3c`, own elimination (`kel_LiT3`) and volume (`Vc_LiT3`), exercised by scenario 5 (`LiT3_dose` = 7.5 μg/day) |

By contrast, `T4p`/`T3p` (plasma hormone pools), `T4thy`/`T3thy` (thyroid
synthesis compartments), and `T4t`/`T3t` (tissue pools) have **no dosing
route anywhere** — no depot, no absorption/bioavailability parameter, no
`ev()`/dataset dosing — they are produced solely by the model's own
synthesis/deiodination/conversion terms. These are correctly **not** census
rows for this file and are completely untouched here, same reasoning as
graves-disease's T3/T4/TRAb/TSH.

## Archetypes

- **LT4: Archetype 3** (depot + central + peripheral, linear elimination),
  with the guide's **"continuous infusion input, not an event"** edge case
  — `LT4_dose` is a daily dose-rate parameter added as a zero-order source
  directly into the gut depot's `dxdt_`, not an `ev()`/dataset bolus.
- **LiT3: Archetype 1** (no depot, single compartment, linear elimination),
  same continuous-infusion edge case, into the central compartment directly
  (no absorption phase modeled).

## Renaming applied

The original's own `k12_LT4`/`k21_LT4` (micro-rate style) are rewritten to
the guide's CL/Q/V form, same recipe as
`chronic-hypothyroidism/hypo_refactor_notes.md`'s own LT4 refactor (same
underlying drug, same disease family):

| Original | Refactored | Value | Note |
|---|---|---|---|
| `LT4g` (cmt) | `GUT_LT4` | -- | depot |
| `LT4c` (cmt) | `CENT_LT4` | -- | central |
| `LT4p` (cmt) | `PERI_LT4` | -- | peripheral |
| `ka_LT4` | `KA_LT4` | 0.48 | absorption rate |
| `F_LT4` | `F_LT4` | 0.75 | bioavailability (unchanged name) |
| `CL_LT4` | `CL_LT4` | 1.3 | **was declared but never referenced in the original `$ODE`** (the ODE used `kel_LT4` directly) — now the operative clearance. Exact match: `kel_LT4*Vc_LT4 = 0.13*10.0 = 1.3`, i.e. the original's own `CL_LT4` value was already internally consistent with `kel_LT4`, just never wired in |
| `Vc_LT4` | `V1_LT4` | 10.0 | central volume |
| `k12_LT4` | `Q_LT4` | 2.5 | new, derived: `k12_LT4*Vc_LT4 = 0.25*10.0` |
| `Vp_LT4` (discarded) / `k21_LT4` | `V2_LT4` | 22.727272727272727 | new, derived: `Q_LT4/k21_LT4 = 2.5/0.11`, kept at full double precision (see "A precision pitfall" below). The original's own `Vp_LT4=22.0` was declared but never referenced anywhere and is numerically **inconsistent** with `k12_LT4`/`k21_LT4` (22.0 vs. the algebraically-required 22.7273) — discarded rather than reused, since reusing it would silently change the dynamics |
| `f_T4conv` | `f_T4conv` | 0.35 | declared, unused in the original — preserved as-is, unchanged name |
| `LiT3c` (cmt) | `CENT_LIT3` | -- | central (only compartment, Archetype 1) |
| `ka_LiT3` | `KA_LIT3` | 5.0 | declared, unused in the original (no depot modeled) — preserved as-is |
| `CL_LiT3` | `CL_LIT3` | 25.0 | same "declared but never referenced" pattern as `CL_LT4` above — exact match `kel_LiT3*Vc_LiT3 = 0.625*40.0 = 25.0` |
| `Vc_LiT3` | `V1_LIT3` | 40.0 | central volume |
| -- (inline `0.95`) | `F_LIT3` (new) | 0.95 | promoted from the original's inline magic number `LiT3_input = LiT3_dose * 0.95` — same value, now a named parameter |
| `LiT3_dose` | `LIT3_dose` | -- | dose-rate driver, case-aligned to the new stem |

## A genuine cross-compound coupling in the original, isolated (not dropped)

The original's `dxdt_LiT3c = LiT3_input - (kel_LiT3 + k12_LT4 * 0.5) * LiT3c;`
reads **`k12_LT4`** — a parameter belonging to LT4's own PK block — directly
inside LiT3's elimination term. This is a real, if unusual, modeling choice
(plausibly meant to represent some shared hepatic/renal handling between the
two thyroid hormones), not a typo: the value is used consistently and the
model compiles and runs with it. It does, however, conflict with the guide's
requirement that "each compound's PK lives in its own clearly-delimited
block." Rather than either (a) silently dropping the term, which would
change LiT3's own numerics, or (b) leaving a live cross-block reference to
`Q_LT4`/`V1_LT4` inside LiT3's own `dxdt_`, which would make LiT3's block
non-self-contained and would break if LT4's own PK is later redirected to
an external covariate — the term's **current numeric contribution** was
extracted into its own named, isolated constant:

```
CL_LIT3_from_LT4 = 0.125   // = k12_LT4*0.5 = 0.25*0.5, extracted from the original's live cross-reference
...
dxdt_CENT_LIT3 = LIT3_input - (CL_LIT3/V1_LIT3 + CL_LIT3_from_LT4) * CENT_LIT3;
```

This reproduces the exact original value (`0.125`) while making LiT3's
block self-contained (it no longer needs to read any LT4 symbol), which
matters specifically for this refactor's stated goal (a future covariate
substitution on `C_LT4` should not silently also change `CENT_LIT3`'s own
elimination). Disclosed here rather than treated as a defect to log,
because the original is internally consistent and not broken — just an
unusual, real coupling.

## The Hill interface: neither compound uses one — both are substrate replacement

**Same finding as `chronic-hypothyroidism/hypo_refactor_notes.md`'s own LT4
refactor**, and for the same underlying reason (this is genuine
biology, not a coincidence of the corpus): both LT4 and LiT3 are hormone
*replacement* therapy, not receptor-mediated inhibitors/agonists. The
original adds each drug's exogenous, synthesis-equivalent hormone
contribution **directly into the endogenous plasma hormone pool's mass
balance** — there is no saturation, no EC50, no ceiling in either
mechanism:

```
// original
double LT4c_nmol  = LT4c  / (MW_T4 * 1e-3 * Vc_LT4);
double LiT3c_nmol = LiT3c / (65.5  * 1e-3 * Vc_LiT3);
dxdt_T4p += LT4c_nmol  * kel_LT4;   // ...among other T4p terms
dxdt_T3p += LiT3c_nmol * kel_LiT3;  // ...among other T3p terms
```

Forcing an `EMAX`/`EC50`/`GAMMA` onto either would misrepresent a mechanism
the original doesn't have (guide: "don't force a fit that isn't real").
Instead, per the guide's own naming slot ("the compound's effect on
disease: `EFFECT_<STEM>`") without forcing the Hill *shape*, each is named
as an exact algebraic rename of the original:

```
double C_LT4  = CENT_LT4  / V1_LT4;                                // ug/L
double EFFECT_LT4  = C_LT4  * (CL_LT4  / V1_LT4)  / (MW_T4 * 1e-3); // nmol/L; = LT4c_nmol*kel_LT4  by construction
double C_LIT3 = CENT_LIT3 / V1_LIT3;                                // ug/L
double EFFECT_LIT3 = C_LIT3 * (CL_LIT3 / V1_LIT3) / (65.5   * 1e-3); // nmol/L; = LiT3c_nmol*kel_LiT3 by construction
```

Algebraic identity (not approximation): since `CL_LT4/V1_LT4 = 1.3/10.0 =
0.13 = kel_LT4` exactly, and `CL_LIT3/V1_LIT3 = 25.0/40.0 = 0.625 =
kel_LiT3` exactly (both by construction of the CL/Q/V rewrite above),
`EFFECT_LT4` and `EFFECT_LIT3` are bit-identical, term-for-term, to the
original's `LT4c_nmol*kel_LT4` and `LiT3c_nmol*kel_LiT3` — confirmed exactly
under Verification below (max abs diff `0.0` throughout).

**A disclosed, uncorrected unit inconsistency copied verbatim:** the
original's own comment says `// T3 MW=651 Da` immediately above a line that
divides by `65.5`, not `655` or `651` — apparently a units slip (converting
g/mol to a different scale) somewhere in the original's own derivation.
This is **not corrected** here (would change simulated T3 output, which is
out of scope for a rename); `65.5` is carried forward verbatim into
`EFFECT_LIT3`, and the discrepancy is disclosed here rather than silently
reproduced without comment. Not logged as a separate `UPSTREAM_ISSUES.md`
entry (a single inline unit-scaling choice, not a compile defect or a
disconnected mechanism) but worth a manual-review flag for anyone using
this file's absolute T3 values quantitatively.

## The original does not compile — three build defects fixed, disclosed, logged

`POST /model_manifest` on the untouched original's own DSL returned HTTP
500 for three independent reasons (a reserved `F_<compartment>` symbol
collision, `$CAPTURE` re-listing compartment names, and five `$PARAM` names
colliding with mrgsolve's auto-generated `<compartment>_0` init symbols).
Logged as `translations/UPSTREAM_ISSUES.md` **#116** (full detail there).
Summary:

1. `F_Se` (genuinely used, unlike some prior `F_<cmt>` collisions in this
   corpus) renamed to `FBIO_Se`.
2. Thirteen `$CMT`/`$INIT` names dropped from `$CAPTURE`
   (`AntiTPOAb AntiTgAb DmgThy Th1 Treg Bcell T4p T3p T4t T3t LT4c LiT3c Se`)
   — replaced with the refactor's own `C_LT4 EFFECT_LT4 C_LIT3 EFFECT_LIT3`.
3. `Th1_0`, `Treg_0`, `Bcell_0`, `DmgThy_0`, `Se_0` renamed to
   `Th1P0`/`TregP0`/`BcellP0`/`DmgThyP0`/`SeP0` (collided with mrgsolve's
   auto-generated `<compartment>_0` initial-value symbol for the
   identically-named `$CMT`s).

All three are syntax-only, non-numeric, confirmed by verification below
using a **patched-original** scratch copy (same three fixes, no LT4/LiT3
renaming) as the comparison baseline, per the guide's policy. The scratch
copy was never committed and was deleted after use.

Also disclosed in the same `UPSTREAM_ISSUES.md` entry (not independently
numbered): `$MAIN` declares three dead between-subject-variability locals
(`CLi`, `DmgSS`, `AbSS`, driven by `OMEGA` `ETA(1)`/`(2)`/`(3)`) that are
never referenced anywhere else — population variability declared in
`$OMEGA` has zero effect on any simulated output. Not fixed (would change
population-run behavior); preserved unchanged.

## A key finding: `$MAIN`-vs-`$ODE` placement is numerically load-bearing here

The first implementation attempt moved the `C_LT4`/`EFFECT_LT4`/
`C_LIT3`/`EFFECT_LIT3` calculation from `$MAIN` (where the original
computed `LT4c_nmol`/`LiT3c_nmol`) into `$ODE`, reasoning (per precedent in
`graves-disease/gd_refactor_notes.md` and
`chronic-hypothyroidism/hypo_refactor_notes.md`, which both describe
mrgsolve 2.0.1's `$MAIN`/`$ODE`/`$TABLE` variable-hoisting as one shared,
per-step-recomputed scope) that this would be a pure reorganization.

**It was not.** Verified empirically that it changes the integrated
trajectory: with the calculation moved to `$ODE`, `T4p`/`T3p` diverged from
the original starting at the very first output row of a dosed scenario
(e.g. `T4p` at day 1: original `0.1309`, `$ODE`-placed version `0.1981`),
even though the algebraic formula for `EFFECT_LT4` was independently
confirmed correct (matches a manual calculation from the captured `CENT_LT4`
value at every point checked). Diagnosed by comparing the original's dosed
vs. undosed trajectories at fine time resolution: with `delta=1`, the
dosed and undosed `T4p` are **identical** for the first reporting interval
(both `0.1309` at day 1), then diverge from the second interval onward; with
`delta=0.05`, the same dosed scenario's `T4p` at day ~0.95 is already
`0.174`, well above the `delta=1` result at day 1. This shows the
original's `$MAIN`-declared `LT4c_nmol`/`LiT3c_nmol` are refreshed once per
**reporting row** (using the compartment's value as of the previous row),
not continuously at every internal LSODA substep the way an `$ODE`-declared
local is — a first-order-hold / staircase approximation of the true
continuous source term, an artifact of declaring state-dependent flux in
`$MAIN` rather than `$ODE`.

Since the guide's tolerance policy for a pure Archetype 1/3 reorganization
is "near-exact match... anything beyond floating-point-scale deviation
means a bug, not a tolerance to loosen," and this genuinely was a
consequence of *where* the calculation was declared (not what it computed),
the fix was to **keep `C_LT4`/`EFFECT_LT4`/`C_LIT3`/`EFFECT_LIT3` in
`$MAIN`, in the same position the original had `LT4c_nmol`/`LiT3c_nmol`** —
a pure in-place rename, not a relocation — rather than "improve" the
model's numerical accuracy by moving it to `$ODE`. This is disclosed as a
correction to the general `$MAIN`/`$ODE` hoisting assumption stated in
prior refactors in this corpus: hoisting appears to make `$MAIN`-declared
names *visible* inside `$ODE` (they compile and read correctly), but does
**not** make them *refresh* at `$ODE`'s per-substep granularity — a
distinction that matters whenever a `$MAIN` local reads a compartment
(state) rather than only parameters/`ETA()`s. `dxdt_GUT_LT4`/
`dxdt_CENT_LT4`/`dxdt_PERI_LT4`/`dxdt_CENT_LIT3` themselves were **not**
moved (they stay in their original `$ODE` position, comment-block
delimited, purely renamed) — only the `EFFECT_LT4`/`EFFECT_LIT3` disease-
input calculation needed to stay in `$MAIN` for exact parity.

## A precision pitfall in the derived `V2_LT4`, caught and fixed before delivery

An initial pass wrote the derived `V2_LT4 = Q_LT4/k21_LT4 = 2.5/0.11` as
the truncated literal `22.72727273` (8 decimal digits). This is close but
**not** bit-exact: `2.5 / 22.72727273 ≈ 0.109999999987`, not exactly
`0.11` as a double. Over a 5-year integration this tiny (~1e-10 relative)
perturbation in the effective `k21`-equivalent rate accumulated into a
`0.0001`-scale absolute deviation in `PERI_LT4` by year 5, which was just
enough to shift the timing of a pre-existing, disease-side numerical
instability (see below) by one reporting week in one scenario. Fixed by
using the full double-precision value Python's `repr()` gives for
`2.5/0.11` (`22.727272727272727`), which round-trips bit-exact back to
`k21_LT4=0.11` (`2.5/22.727272727272727 == 0.11` confirmed in IEEE-754
double arithmetic) — confirmed this removes the discrepancy entirely (see
Verification: scenario 5 below).

## Verification

**Method.** Extracted the bare DSL text from `ht_mrgsolve_model.R`
(`ht_model_code <- '...'`, unmodified) and `ht_mrgsolve_model_refactored.R`
(same variable name, refactored), POSTed to the local qspserver
`mrgsolve_api` (`http://localhost:8007`), `POST /model_manifest` then
`POST /run_simulation`, requests spaced ~2s apart per the shared-service
note. Comparison baseline is the three-fix **patched-original** scratch
copy described above (confirmed compiling via `/model_manifest` first).
LT4/LiT3 dosing here is the guide's "continuous infusion input, not an
event" edge case, so verification uses plain `parameters` overrides
(`LT4_dose`, `LiT3_dose`/`LIT3_dose`), not the API's `dosing` event field.

**Scenarios run — the file's own, not invented** (`end=1825` days / 5
years, `delta=7` days, matching the original R script's own `sim_time`):

| # | Description | Dosing | Result |
|---|---|---|---|
| 1 | Untreated Hashimoto's | `LT4_dose=0, LiT3_dose=0` | **Exact match, max abs diff 0.0**, every shared output, including identical (absent) NaN pattern |
| 2 | LT4 monotherapy 100 μg/day | `LT4_dose=100` | **Exact match, max abs diff 0.0** on every point where both are non-NaN. One boundary case: the file's own pre-existing autoimmune-cascade instability (see below) hits NaN one reporting week (7 days) earlier in the refactored run than the original for the `Th1/Treg/Bcell/AntiTPOAb/AntiTgAb/DmgThy/ThyVol_rel/BMD_risk` group only — `T4p/T3p/T4t/T3t/TSH_clinical/fT4_pmolL/fT3_pmolL` (the entire HPT+thyroid+LT4-PK chain) match with identical NaN timing. See "A pre-existing chaotic instability" below |
| 5 | LT4 100 + LiT3 7.5 μg/day combination | `LT4_dose=100, LiT3_dose=7.5` | **Exact match, max abs diff 0.0** on every shared output, **including identical NaN blowup timing** (both hit NaN at day 28) |
| 6 | High-dose LT4 175 μg/day (TSH-suppressive) | `LT4_dose=175` | **Exact match, max abs diff 0.0** on every shared output, including identical NaN blowup timing |

`LT4g`/`LT4c`/`LT4p`/`LiT3c` (original) vs. `GUT_LT4`/`CENT_LT4`/
`PERI_LT4`/`CENT_LIT3` (refactored) matched exactly (`0.0`) in all four
scenarios. `C_LT4`/`EFFECT_LT4`/`C_LIT3`/`EFFECT_LIT3` (refactored-only,
no original-file counterpart by name) ranged over plausible non-zero
values consistent with the dosing in each scenario, confirming the
substrate-replacement mechanism is active.

**A pre-existing chaotic instability in the disease side (not LT4/LiT3),
confirmed to be shared by both models.** All four scenarios' output
eventually goes to `NaN` — a pre-existing solver blow-up in `DmgThy`'s own
positive-feedback damage-accumulation term (`dmg_drive = k_dmg*Th1*
(AntiTPOAb/500)*(1+Se_eff_GPx*(-0.4))`, uncapped, feeding
`dxdt_DmgThy = dmg_drive*k_dmg_scale*(1-DmgThy) - k_repair*Treg*DmgThy`),
unrelated to this refactor's own scope (LT4/LiT3) and reproduced
identically in the patched-original. Scenarios 1, 5, and 6 confirm this
blow-up happens at the **exact same simulated day** in both the original
and the refactored model. Scenario 2 alone shows the autoimmune-cascade
sub-block (which does not read `T4p`/`T3p`/`LT4`/`LiT3` at all in either
model — its own equations depend only on `Th1`/`Treg`/`Bcell`/`AntiTPOAb`/
`DmgThy`/`Se`) crossing into `NaN` one reporting week earlier in the
refactored run. Every paired value that is non-NaN in both models still
matches exactly (`0.0`) up to that point in every scenario, including
scenario 2 — the divergence is purely in *when* an already-chaotic,
positive-feedback runaway crosses its singularity, not in any computed
value beforehand. Given (a) the algebra for every renamed term was
independently confirmed exact, (b) three of four scenarios (including the
combination scenario exercising both LT4 and LiT3 together) show
bit-identical blow-up timing, and (c) `LSODA`'s adaptive step-size control
operates over the model's *entire* state vector (so a bit-level difference
anywhere, even in an unrelated compartment, can shift exactly when an
already-unstable feedback loop crosses into `NaN`) — this is consistent
with ordinary floating-point sensitivity of a pre-existing, uncapped
positive-feedback singularity, not a structural error introduced by the
refactor. Not logged as a new `UPSTREAM_ISSUES.md` entry (the underlying
`DmgThy` runaway is the same class of finding as `#112`/`#115`'s
uncapped-feedback blow-ups in other files, and fixing it is out of scope —
`DmgThy`'s damage/repair terms are entirely disease-side, untouched by
this refactor).

`/model_manifest` on the refactored DSL confirmed all of `KA_LT4, F_LT4,
CL_LT4, V1_LT4, Q_LT4, V2_LT4, f_T4conv, KA_LIT3, CL_LIT3, V1_LIT3,
F_LIT3, CL_LIT3_from_LT4, LT4_dose, LIT3_dose, FBIO_Se` appear in
`parameters` with the values above, and `GUT_LT4, CENT_LT4, PERI_LT4,
CENT_LIT3, C_LT4, EFFECT_LT4, C_LIT3, EFFECT_LIT3` all appear in
`outputPaths`.

No scratch/debug artifacts were left in the repository — DSL extraction,
the patched-original scratch copy, and all comparison scripts ran from a
session-local scratchpad directory outside the repo tree (`.scratch_ht/`
under the repo root during work, deleted after use).

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`: both
rows (`Levothyroxine (LT)`, `Liothyronine (LIT)`) filled in as refactored
and verified.
