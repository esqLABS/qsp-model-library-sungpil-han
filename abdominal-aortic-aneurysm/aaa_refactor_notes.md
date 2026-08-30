# Refactor notes — `abdominal-aortic-aneurysm/aaa_mrgsolve_model.R`

Scope: all **three** compounds this file models — **Doxycycline (DOXY)**,
**Propranolol/beta-blocker (BB)**, and **Statin (STAT)** — per the three
existing rows in `driver-patches/data/compound_perturbation_census.md`, all
classified "Redirect concentration (clean single site)". There is no fourth
compound to leave untouched: the file's only other pharmacological input,
ACE-inhibitor ("Perindopril"), has no PK compartments of its own (it is a
fixed `ACEi_dose` level with a simple additive Emax/EC50 term) and is
therefore not one of the three targets — it is byte-identical to the
original everywhere it appears.

## Headline finding: Doxycycline's effect reads AORTIC TISSUE, not plasma

This file was flagged going in as modeling "one compound acting on tissue
not plasma concentration," and that is exactly what the code does. Doxycycline
has **four** PK compartments, not three:

```
DGUT   : gut depot
DCENT  : central (plasma)
DPERIPH: peripheral
DAORTA : "aortic tissue (mg/L equivalent)" -- fed FROM plasma, not a straight
         partition of it: dxdt_DAORTA = kta_d*(DCENT/Vd_d) - kte_d*DAORTA
```

Both original effect terms that doxycycline drives —
`InhMMP9_doxy = Imax_d9 * Ct_doxy / (IC50_d9 + Ct_doxy)` and
`InhMMP2_doxy = Imax_d2 * Ct_doxy / (IC50_d2 + Ct_doxy)` — read `Ct_doxy`,
which is set to `DAORTA` (the tissue compartment), **never** `Cp_doxy`
(`DCENT/Vd_d`, plasma). Plasma doxycycline is computed and was captured in
the original purely for reporting/plotting; no PD equation reads it.

Per the task's explicit instruction, the exposed concentration is still
named `C_DOXY` (following the naming convention's `C_<STEM>` slot
regardless of which physical quantity it represents) — but it is defined as
`double C_DOXY = TISSUE_DOXY;` (the renamed `DAORTA`), not from
`CENT_DOXY/V1_DOXY`. This is documented at every place that matters so a
future maintainer cannot miss it:
- the `$PARAM` header comment block for the DOXY section,
- the `$CMT` description of `TISSUE_DOXY`,
- inline at the `C_DOXY` assignment itself (`/* the exposed, PD-facing
  concentration: AORTIC TISSUE level, NOT plasma */`),
- and again at both `EFFECT_DOXY_MMP9`/`EFFECT_DOXY_MMP2`.

The plasma level is retained under a deliberately non-convention name,
`Cp_doxy_plasma`, so it is visually obvious it is *not* the blessed
`C_<STEM>` — informational only (kept in `$CAPTURE` since the original
captured and plotted it), read by nothing.

## Archetype per compound

**Doxycycline (DOXY): Archetype 3 (depot + central + peripheral), plus a
genuinely distinct 4th tissue compartment** — the guide's own carve-out
("two [concentration variables] only when a genuinely different tissue site
matters") applies exactly here. Nothing was flattened or added beyond what
the original already modeled; the 4th compartment was renamed
(`DAORTA` → `TISSUE_DOXY`), not invented.

**Propranolol/beta-blocker (BB): Archetype 3 (depot + central + peripheral),
linear, no tissue distinction.** Single exposed concentration `C_BB =
CENT_BB/V1_BB` (plasma), matching what the original's `Cp_prop` already fed
into the SBP effect.

**Statin (STAT): Archetype 3 (depot + central + peripheral), linear, no
tissue distinction.** Single exposed concentration `C_STAT =
CENT_STAT/V1_STAT` (plasma). Statin is the one compound here with **three**
separate downstream effects rather than one (see below).

## Naming (all three compounds)

| Role | DOXY | STAT | BB |
|---|---|---|---|
| Depot | `GUT_DOXY` (was `DGUT`) | `GUT_STAT` (was `SGUT`) | `GUT_BB` (was `BGUT`) |
| Central | `CENT_DOXY` (was `DCENT`) | `CENT_STAT` (was `SCENT`) | `CENT_BB` (was `BCENT`) |
| Peripheral | `PERI_DOXY` (was `DPERIPH`) | `PERI_STAT` (was `SPERIPH`) | `PERI_BB` (was `BPERIPH`) |
| Tissue (DOXY only) | `TISSUE_DOXY` (was `DAORTA`) | — | — |
| `KA_/F_/CL_/V1_/Q_/V2_` | was `ka_d/F_d/CL_d/Vd_d/Q_d/Vp_d` | was `ka_s/F_s/CL_s/Vd_s/Q_s/Vp_s` | was `ka_bb/F_bb/CL_bb/Vd_bb/Q_bb/Vp_bb` |
| Tissue distribution/elimination (DOXY only) | `KTA_DOXY`/`KTE_DOXY` (was `kta_d`/`kte_d`) | — | — |
| Exposed concentration | `C_DOXY` (**tissue**, was `Ct_doxy`=`DAORTA`) | `C_STAT` (plasma, was `Cp_stat`) | `C_BB` (plasma, was `Cp_prop`) |

All parameter values are copied verbatim from the original (`Vd_d=110.0` →
`V1_DOXY=110.0`, etc.) — nothing invented, nothing defaulted, nothing
dropped. All were already in the CL/Q/V convention in the original (no
`k10`/`k12`/`k21` micro-constant form to convert), so the PK ODEs are a pure
identifier rename with identical arithmetic.

## The Hill interface: one compound, several targets

The guide's example shows one `EFFECT_<STEM>` per compound. Two of the three
compounds here genuinely drive more than one downstream target with
*independent* Emax/EC50 pairs, so each got its own suffixed effect name
rather than being forced into a single shared term (which the guide
explicitly warns against: "never collapse several drugs into one shared
Hill term"):

- **Doxycycline** → `EFFECT_DOXY_MMP9` (was `InhMMP9_doxy`, params
  `Imax_d9`/`IC50_d9` → `EMAX_DOXY_MMP9`/`EC50_DOXY_MMP9`) and
  `EFFECT_DOXY_MMP2` (was `InhMMP2_doxy`, params `Imax_d2`/`IC50_d2` →
  `EMAX_DOXY_MMP2`/`EC50_DOXY_MMP2`). Both read `C_DOXY` (tissue).
- **Statin** → `EFFECT_STAT_MMP9` (was `InhMMP9_stat`, params
  `Imax_s9`/`IC50_s9` → `EMAX_STAT_MMP9`/`EC50_STAT_MMP9`),
  `EFFECT_STAT_NFKB` (TNF suppression; was the **inline hardcoded**
  `0.30 * Cp_stat / (0.05 + Cp_stat)` with no named parameters at all —
  promoted to `EMAX_STAT_NFKB=0.30`/`EC50_STAT_NFKB=0.05`), and
  `EFFECT_STAT_ROS` (ROS suppression; was inline hardcoded
  `0.35 * Cp_stat / (0.08 + Cp_stat)` — promoted to
  `EMAX_STAT_ROS=0.35`/`EC50_STAT_ROS=0.08`). All three read `C_STAT`
  (plasma). Promoting the two hardcoded pairs into named `$PARAM` Hill
  parameters is a **rename/promotion, not a refit** — same values, same
  ratio shape, now discoverable as parameters instead of buried literals.
- **Propranolol** → `EFFECT_BB` (was `SBP_red_bb`, params
  `Imax_bb`/`IC50_bb` → `EMAX_BB`/`EC50_BB`). Reads `C_BB` (plasma).

Every one of these was already exactly `Emax*X/(EC50+X)` in the original —
a plain ratio, no Hill exponent — so per the guide's rename rule each got
`GAMMA_<STEM>_<TARGET> = 1.0` as an explicit named parameter (no original
Hill coefficient existed) and the formula was rewritten
`EMAX*pow(X,GAMMA)/(pow(EC50,GAMMA)+pow(X,GAMMA))`. With `GAMMA=1`,
`pow(x,1)==x` exactly, so this is arithmetically identical to the original,
confirmed by the exact-match verification below.

Where the original combined two drugs on the same target
(`InhMMP9_total = 1 - (1-InhMMP9_doxy)*(1-InhMMP9_stat)`, multiplicative
independence on MMP-9), that combination is untouched — only the two inputs
were renamed to `EFFECT_DOXY_MMP9`/`EFFECT_STAT_MMP9`. This already matched
the guide's "combine only at the point disease equations actually use them"
rule before the refactor.

## qspserver compatibility

Attempted the guide's literal instruction to declare `C_<STEM>`/`EFFECT_<STEM>`
in `$PARAM` with a `0.0` default for `/model_manifest` discoverability, and
hit a hard mrgsolve constraint: `$PARAM @annotated` values compile to
`const double&` references, which are **read-only** inside `$ODE`/`$MAIN`.
Attempting `C_DOXY = TISSUE_DOXY;` against a `$PARAM`-declared `C_DOXY`
fails to build (`error: assignment of read-only reference 'C_DOXY'`),
confirmed on the qspserver `mrgsolve_api` container. This is why every one
of the 13 prior refactors in this repo (TCZ, LT4, MX, ...) keeps
`C_<STEM>`/`EFFECT_<STEM>` as plain `double` locals computed in `$ODE`
instead — that pattern was followed here too. Discoverability is satisfied
the same way the `myopia-progression` refactor did it: every `C_<STEM>` and
`EFFECT_<STEM>` this file creates (`C_DOXY`, `C_STAT`, `C_BB`,
`EFFECT_DOXY_MMP9`, `EFFECT_DOXY_MMP2`, `EFFECT_STAT_MMP9`,
`EFFECT_STAT_NFKB`, `EFFECT_STAT_ROS`, `EFFECT_BB`) is listed in `$CAPTURE`,
which populates `/model_manifest`'s `outputPaths` — so a client can still
discover every one of these names from the manifest without reading the DSL
text, just not as an overridable `$PARAM`. The Hill parameters themselves
(`EMAX_*`, `EC50_*`, `GAMMA_*`) **are** ordinary `$PARAM` values and so are
fully overridable, matching the TCZ precedent.

Confirmed via `POST /model_manifest` on the qspserver `mrgsolve_api`
container (`http://localhost:8007`) that the refactored DSL compiles and
that `outputPaths` includes all 9 exposed names and `parameters` includes
every renamed `$PARAM`.

## A pre-existing defect found while verifying (not fixed, logged upstream)

Two independent, pre-existing build breaks, neither introduced by this
refactor, both reproduced from the untouched original alone via the
qspserver `mrgsolve_api` container:

1. **`$CAPTURE` repeats ten compartment names already in `$CMT`**
   (`MMP9, MMP2, ELAST, COLLAG, VSMC, ILT, DIAM, TNF, ROSO, MAC`). mrgsolve
   2.0.1 refuses to build any model whose `$CAPTURE` repeats a compartment
   name (`invalid class "mrgmod" object: compartment should not be in
   $CAPTURE: ...`) — the same defect class as several other files' entries
   in `translations/UPSTREAM_ISSUES.md` (e.g. issue affecting
   `rheumatoid-arthritis`).
2. **`$PARAM` names `MMP9_0`/`MMP2_0` collide with mrgsolve's own
   auto-generated `<CMT>_0` initial-condition symbol** for compartments
   `MMP9`/`MMP2`. The original's `$MAIN` line `MMP9_0 = MMP9_0;` (likewise
   `MMP2_0`) is a self-assignment that only "works" if the two `MMP9_0`
   symbols were distinct — they are not, so the C++ build fails with
   `conflicting declaration 'double& MMP9_0'` / `previous declaration as
   'const double& MMP9_0'`. This is a **new** defect pattern not previously
   logged (distinct from the `$CMT`+`$INIT` redeclaration pattern in issues
   #34/#36).

**Verification-only workaround** (in-memory scratch copies of both DSLs,
never applied to either checked-in file): (a) stripped the ten duplicated
compartment names from `$CAPTURE`; (b) renamed the `$PARAM` symbols
`MMP9_0`/`MMP2_0` to `MMP9_0_BASE`/`MMP2_0_BASE` everywhere **except** the
one `$MAIN` line that assigns the compartment's own auto-generated init
symbol (`MMP9_0 = MMP9_0_BASE;`) — a pure symbol rename that changes no
equation and no numeric value, applied identically to the original and the
refactored DSL. Neither `aaa_mrgsolve_model.R` nor
`aaa_mrgsolve_model_refactored.R` contains this workaround; both still
carry the original defect pattern exactly (the refactored file's
disease-PD block, including these exact two lines, is byte-for-byte
untouched — this was not one of the 3 target compounds). Logged as a new
numbered entry in `translations/UPSTREAM_ISSUES.md`.

## Verification

Per the guide's mandatory protocol: ran the original file's own dosing
scenarios through both the (workaround-patched, buildability-only) original
and the refactored model via the qspserver `mrgsolve_api` service
(`POST /model_manifest`, `POST /run_simulation`), comparing every shared
`$CAPTURE`d output point-by-point.

**Scenario 6 (Triple therapy — Doxycycline 100 mg q12h + Statin 40 mg q24h +
Propranolol 40 mg q8h), full 1-year horizon (`end=8760`, `delta=24`,
matching the original's own `scenario6` exactly):** ran to completion on
both models (369 time points). All 14 shared disease-PD outputs
(`MMP9, MMP2, ELAST, COLLAG, VSMC, ILT, DIAM, TNF, ROSO, MAC, SBP_mmHg,
Wall_Stress_idx, Rupture_P, MMP9_reduction`) and all 4 PK/exposure crosswalk
pairs (`Cp_doxy`↔`Cp_doxy_plasma`, `Ct_doxy`↔`C_DOXY`, `Cp_stat`↔`C_STAT`,
`Cp_prop`↔`C_BB`) matched with **max abs diff = 0.0** across the first 298
of 369 points (see the NaN note below for the remainder). This is the exact
match the guide anticipates for a pure structural reorganization plus a
Hill-parameter rename with `gamma=1` (no fitting was needed anywhere in
this file).

**Scenario 2 (Doxycycline monotherapy, 100 mg q12h) — isolates the
tissue-vs-plasma finding:** the original's own 8760h window for this
scenario (and even a 480h/20-day window) exceeds the qspserver API's
default `maxsteps=20000` on **both** models identically
(`[lsoda] 20000 steps taken before reaching tout`) — a solver-budget limit,
not a mismatch (see the guide's step-count caveat). Doxycycline alone
suppresses the disease system's runaway feedback loop (see below) far less
than the triple-therapy combination, so its own ODE gets numerically
stiffer sooner. Shortened to a **60h window** (`end=60`, `delta=12`,
`addl` scaled down to keep the same 100 mg q12h regimen), which both models
solved without error. Result: **max abs diff = 0.0** across all 14 shared
disease-PD outputs and both PK crosswalk pairs, over all 6 time points.
Concretely: at t=48h, `Cp_doxy`/`Cp_doxy_plasma` = 0.7491 mg/L (plasma)
while `Ct_doxy`/`C_DOXY` = 0.8807 mg/L-equivalent (aortic tissue) — visibly
different values, both reproduced exactly, confirming the tissue-vs-plasma
distinction survived the refactor unchanged and that `EFFECT_DOXY_MMP9`/
`EFFECT_DOXY_MMP2` are driven from the tissue value as intended.

**Result: exact match (max abs diff = 0.0) for every output, every
scenario, everywhere both models could be run to completion.** This is the
expected outcome for archetype 3 (pure PK reorganization) plus `gamma=1`
Hill renames/promotions — there is no approximation anywhere in this file
to introduce error.

## Anything else worth flagging

- **Both models diverge to numeric overflow (`NaN`) partway through the
  1-year horizon, identically.** In scenario 6, both original and
  refactored produce identical finite values through t=6960h (day 290) and
  then both go to `NaN` from t=6984h onward, at the exact same index. In
  the undosed natural-history case (no compounds at all) the same
  `MAC`→`MMP9`→`ELAST`/`DIAM`→`WALL_STR`→`MAC` positive-feedback loop blows
  up far faster (by ~t=480h). This looks like a genuine pre-existing
  modeling issue in the original — the PK rate constants are declared in
  h⁻¹ (`ka_d`, `CL_d`, etc.) while every disease-PD rate constant is
  declared in day⁻¹ (`k_MAC_in`, `k_MMP9_syn`, `k_diam_grow`, ...) but both
  sets of ODEs are integrated on the **same** time axis (hours — confirmed
  by the original R script's own `sim_time <- seq(0, 8760, by=24)` and
  `Day = time/24` conversions, and by the PK rate units only making sense
  per hour). None of the day⁻¹ disease rates are divided by 24 before use,
  so the disease system runs roughly 24× faster than its own labeled units
  imply, which is consistent with how quickly it blows up relative to the
  documented "~2 mm/yr" baseline expansion rate. This is **not** something
  introduced by the refactor (confirmed identical in the untouched,
  workaround-patched original) and **not** specific to the three refactored
  compounds (it reproduces with zero drugs dosed) — flagged here, and
  logged as an upstream defect, rather than fixed.
- The beta-blocker's heart-rate parameters (`Imax_bb_hr`/`IC50_bb_hr`,
  renamed `EMAX_BB_HR`/`EC50_BB_HR`) and the three dosing-switch flags
  (`DOXY_ON`/`STAT_ON`/`PROP_ON`) are declared in the original `$PARAM`
  block but **never referenced by any `$ODE`/`$TABLE` expression** —
  dead/unused, same treatment as the `RTOT_TCZ` precedent noted in the
  `rheumatoid-arthritis` refactor: preserved as-is (unused), not removed,
  not wired up, per "don't invent, don't drop the original's parameters."
  (Scenario on/off in this file's own R script is controlled entirely by
  which `events` a scenario attaches, not by these flags.)
- Compartment ordering is unchanged (`GUT_DOXY` is still compartment 1,
  etc.) — renaming was done in place without reordering, so 1-based
  compartment-number dosing (as used for `POST /run_simulation`'s `dosing`
  field) is unaffected between the original and refactored files.
- The ACE-inhibitor term (`ACEi_dose`, `Imax_acei`, `IC50_acei`,
  `SBP_red_acei`) is completely untouched — it has no PK compartments and
  is not one of the three target compounds.

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`, the
`abdominal-aortic-aneurysm | Doxycycline`, `abdominal-aortic-aneurysm |
Propranolol (BB)`, and `abdominal-aortic-aneurysm | Statin` rows.
