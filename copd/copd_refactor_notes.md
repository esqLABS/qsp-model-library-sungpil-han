# Refactor notes — `copd/copd_mrgsolve_model.R`

Four compounds refactored, per the existing rows in
`driver-patches/data/compound_perturbation_census.md`, all classified
"Redirect concentration (clean single site)": **ICS lung (Budesonide)**,
**LABA lung (Salmeterol)**, **LAMA lung (Tiotropium)**, **PDE4I
(Roflumilast)**. These are the only compounds modeled in this file —
nothing else was touched.

## Archetypes

- **LAMA, LABA, ICS** — all three inhaled compounds share the identical
  shape: **Archetype 3** (depot + central + peripheral, linear
  elimination) for the PK skeleton (`GUT_<STEM>`/`CENT_<STEM>`/
  `PERI_<STEM>`), **but with a deliberate deviation from the archetype's
  canonical `C_<STEM> = CENT/V1`**. In the original, the compound's *only*
  effect on disease reads a **lung (epithelial-lining-fluid) concentration**
  computed straight from the depot state (`LAMA_conc_lung =
  LAMA_LUNG/(F_lung_LAMA*200.0)`, etc.) — the systemic/plasma concentration
  (`Cp_LAMA = LAMA_C/Vc_LAMA`, etc.) is computed and was captured in the
  original's own output, but **never read by any effect term**. This is
  confirmed by direct inspection of every `E_LAMA`/`E_LABA`/`E_ICS`
  expression in the original, all of which read only the `*_conc_lung`
  variable. Per the guide's "two only when a genuinely different tissue
  site matters" allowance (normally invoked for a bespoke 4th tissue
  compartment), the lung fluid concentration and the systemic plasma
  concentration genuinely are two different tissue sites here — except no
  new compartment was needed, since the existing lung depot already *is*
  the tissue of effect. So: `C_<STEM> = GUT_<STEM>/(F_LUNG_<STEM>*200.0)`
  (the exposed, PD-driving concentration), and `C_<STEM>_PLASMA =
  CENT_<STEM>/V1_<STEM>` (the systemic concentration, kept only as an
  informational, non-driving diagnostic — same treatment as
  `Cp_doxy_plasma` in `abdominal-aortic-aneurysm/aaa_refactor_notes.md`).
- **PDE4I** — **Archetype 1** (no depot, single compartment, linear
  elimination). The original declares `ka_PDE4i`/`F_PDE4i` (oral
  absorption rate/bioavailability) in `$PARAM`, but **neither is ever
  referenced in `$ODE`** — dosing goes directly into `PDE4i_C` as a bolus
  (`dose_PDE4i <- ev(cmt = "PDE4i_C", amt = 400000, ...)`), bypassing any
  absorption compartment entirely. `dxdt_PDE4i_C = -CL_PDE4i * Cp_PDE4i`
  is a plain single-compartment linear-elimination equation. `KA_PDE4I`/
  `F_PDE4I` are carried over unchanged as disclosed dead parameters (same
  treatment as `KA_MEM` in `alzheimers-disease/ad_refactor_notes.md`).
  Because PDE4I is oral (not inhaled), its exposed concentration is the
  ordinary central/plasma concentration — no separate lung tissue site,
  unlike the three inhaled compounds above.

## Renaming (values unchanged from the original)

| Original | Refactored |
|---|---|
| `LAMA_LUNG` / `LAMA_C` / `LAMA_P` | `GUT_LAMA` / `CENT_LAMA` / `PERI_LAMA` |
| `LABA_LUNG` / `LABA_C` / `LABA_P` | `GUT_LABA` / `CENT_LABA` / `PERI_LABA` |
| `ICS_LUNG` / `ICS_C` / `ICS_P` | `GUT_ICS` / `CENT_ICS` / `PERI_ICS` |
| `PDE4i_C` | `CENT_PDE4I` |
| `ka_LAMA` / `Vc_LAMA` / `Vp_LAMA` / `F_sys_LAMA` | `KA_LAMA` / `V1_LAMA` / `V2_LAMA` / `F_LAMA` |
| `ka_LABA` / `Vc_LABA` / `Vp_LABA` / `F_sys_LABA` | `KA_LABA` / `V1_LABA` / `V2_LABA` / `F_LABA` |
| `ka_ICS` / `Vc_ICS` / `Vp_ICS` / `F_sys_ICS` | `KA_ICS` / `V1_ICS` / `V2_ICS` / `F_ICS` |
| `ka_PDE4i` / `Vc_PDE4i` / `F_PDE4i` | `KA_PDE4I` / `V1_PDE4I` / `F_PDE4I` (KA/F both dead, see above) |
| `Emax_LAMA` / `Emax_LABA` / `Emax_ICS` / `Emax_PDE4i` | `EMAX_LAMA` / `EMAX_LABA` / `EMAX_ICS` / `EMAX_PDE4I` |
| `Hill_ICS` | `GAMMA_ICS` (value unchanged, 1.5) |
| `LAMA_conc_lung` / `LABA_conc_lung` / `ICS_conc_lung` | `C_LAMA` / `C_LABA` / `C_ICS` |
| `Cp_LAMA` / `Cp_LABA` / `Cp_ICS` (informational only) | `C_LAMA_PLASMA` / `C_LABA_PLASMA` / `C_ICS_PLASMA` |
| `Cp_PDE4i` | `C_PDE4I` |
| `E_LAMA` / `E_LABA` / `E_ICS` / `E_PDE4i` | `EFFECT_LAMA` / `EFFECT_LABA` / `EFFECT_ICS` / `EFFECT_PDE4I` |
| `DOSE_PDE4i` | `DOSE_PDE4I` (case-only, for consistency with the renamed compartment/stem) |
| `F_lung_LAMA`/`F_lung_LABA`/`F_lung_ICS` | `F_LUNG_LAMA`/`F_LUNG_LABA`/`F_LUNG_ICS` (case-only; bespoke, no convention-table role) |

New (not present in the original, for the named Hill interface):
`GAMMA_LAMA = 1`, `GAMMA_LABA = 1`, `GAMMA_PDE4I = 1` — these make explicit
a shape the original already had implicitly as a plain ratio (no `pow()`,
no Hill exponent term at all). `GAMMA_ICS` is a rename of the original's
own explicit `Hill_ICS = 1.5`, not a new value.

`Q_LAMA`/`Q_LABA`/`Q_ICS`, `EC50_LAMA`/`EC50_LABA`/`EC50_ICS`/`EC50_PDE4I`,
`DOSE_LAMA`/`DOSE_LABA`/`DOSE_ICS`, and `Eos_thresh` already matched the
naming convention (or have no convention-table role) and are unchanged.

## Hill interface: renames, not fits

All four compounds' effect terms were already written as a plain Emax/Hill
ratio in the original — no ODE-solved receptor kinetics to approximate, so
every one of these is a rename, confirmed by the exact verification match
below (not merely "close"):

- **LAMA**: `E_LAMA = DOSE_LAMA * Emax_LAMA * LAMA_conc_lung / (EC50_LAMA +
  LAMA_conc_lung)` is `Emax*C^1/(EC50^1+C^1)` with implicit `gamma=1`.
  `EFFECT_LAMA` uses `pow(C_LAMA, GAMMA_LAMA)` with `GAMMA_LAMA=1` added
  explicit — `pow(x,1)=x`, so this changes nothing numerically (same
  precedent as `alzheimers-disease/ad_refactor_notes.md`, LEC/MEM).
- **LABA**: same shape, same treatment (`EFFECT_LABA`, `GAMMA_LABA=1`).
- **ICS**: `E_ICS_lung = DOSE_ICS * Emax_ICS * pow(ICS_conc_lung,Hill_ICS) /
  (pow(EC50_ICS,Hill_ICS) + pow(ICS_conc_lung,Hill_ICS))`, then
  `E_ICS = E_ICS_lung * Eos_factor` (an eosinophil-count-gated multiplier,
  `Eos_factor = (Eos>Eos_thresh) ? 1.0 : (Eos/Eos_thresh)`). `EFFECT_ICS`
  reproduces the identical two-step expression verbatim, with `Hill_ICS`
  renamed to `GAMMA_ICS` (value unchanged, 1.5) — a rename, not a refit.
  `Eos_factor` is disease/phenotype-side gating on a single compound's own
  effect (not a multi-drug combination), so it stays inside `EFFECT_ICS`
  per the guide's "expressed as one named function...not buried inside a
  combined multi-drug expression" — nothing here combines with another
  drug.
- **PDE4I**: `E_PDE4i = DOSE_PDE4i * Emax_PDE4i * Cp_PDE4i / (EC50_PDE4i +
  Cp_PDE4i)` — same `Emax*C^1/(EC50^1+C^1)` shape, same treatment
  (`EFFECT_PDE4I`, `GAMMA_PDE4I=1`).

No `nls()` fitting was needed or performed for any of the four compounds.

## A second, non-Hill ICS effect that stays outside `EFFECT_ICS`

The original has **two independent ICS-driven terms**, not one: besides
the Hill-shaped `E_ICS` above (used in the IL-8/NE/CRP/emphysema/
exacerbation-rate equations), `dxdt_Eos` uses a **separate, flag-only**
multiplier — `Eos_stim = (DOSE_ICS == 1.0) ? 0.3 : 1.0` — that depends only
on the on/off dosing switch, not on `C_ICS` at all (no Hill shape, no
concentration term anywhere in it). Since it is not a function of
concentration, it does not belong in the guide's concentration-driven Hill
interface and was **kept exactly as the original**, referencing the same
(renamed) `DOSE_ICS` flag. Flagged here, not changed.

## Verified dead parameters, kept and disclosed (not removed)

- **`KA_PDE4I`/`F_PDE4I`**: declared in `$PARAM`, never referenced in
  `$ODE`/`$TABLE` — dosing bypasses absorption entirely (see Archetypes,
  above). Confirmed by grep across the whole original file and by the
  exact-match verification below (removing neither param from the
  computation changes nothing, since neither was ever read).
- **`CRP0`/`Eos0`/`Emph0`** (renamed from `CRP_0`/`Eos_0`/`Emph_0`, see the
  upstream-defect section below): declared as "baseline" documentation
  params but never read anywhere in the DSL body. `Eos_0` is referenced
  from the original's own R-side `dose_response` scenario as a parameter
  override — but that override has **zero effect** on the simulated
  trajectory, because `$INIT` hardcodes `Eos = 200` as a literal,
  independent of the `Eos_0` param. This is a pre-existing authoring bug
  in the checked-in original, reproduced faithfully (not fixed) on both
  sides of the verification below.

## Pre-existing upstream build defects (fixed syntax-only in the delivered file)

**The original does not compile under mrgsolve 2.0.1 at all**, for three
reasons unrelated to any of the four refactored compounds' own PK,
confirmed via `POST /model_manifest` on the untouched original alone:

1. `$CMT` and `$INIT` jointly redeclare all 19 compartments (`$INIT`'s bare
   `NAME = value` assignments duplicate what `$CMT` already declared —
   mrgsolve 2.0.1 treats `$INIT` as an alternative compartment-declaring
   block, not a companion to `$CMT`). Same defect class as
   `UPSTREAM_ISSUES.md` #38/#42/#51/#61.
2. `$CAPTURE` lists 9 names that are already `$CMT` compartments (`IL8
   NE_sput CRP Eos FEV1 Emph PVR AE_cum AE_rate_ann`). Same defect class as
   #56/#57/#61.
3. Six `$PARAM` baseline names (`IL8_0`/`CRP_0`/`Eos_0`/`FEV1_0`/`Emph_0`/
   `PVR_0`) collide with mrgsolve's own auto-reserved per-compartment
   `<CMT>_0` initial-value symbol for the identically-named compartments
   `IL8`/`CRP`/`Eos`/`FEV1`/`Emph`/`PVR`. Same defect class as #60
   (sarcoidosis), here with six collisions instead of three.

```
Error in validObject(.Object) :
  invalid class "mrgmod" object: 1: Duplicated model names: LAMA_LUNG
  LAMA_C LAMA_P LABA_LUNG LABA_C LABA_P ICS_LUNG ICS_C ICS_P PDE4i_C IL8
  NE_sput CRP Eos FEV1 Emph PVR AE_cum AE_rate_ann
  invalid class "mrgmod" object: 2: compartment should not be in
  $CAPTURE: IL8,NE_sput,CRP,Eos,FEV1,Emph,PVR,AE_cum,AE_rate_ann
```

Logged as [`UPSTREAM_ISSUES.md` #63](../translations/UPSTREAM_ISSUES.md).
Per the guide's settled policy, the fix is applied **directly to the
delivered `copd_mrgsolve_model_refactored.R`**, not just to a scratch
copy, all three syntax-only and non-numeric:

- (a) the `$INIT` block was deleted; its 8 non-default initial values
  (`IL8`, `NE_sput`, `CRP`, `Eos`, `FEV1`, `Emph`, `PVR`, `AE_rate_ann`)
  were moved into a new `$MAIN` block using the modern
  `<CMT>_0 = value;` idiom — the all-zero PK compartments (`GUT_*`,
  `CENT_*`, `PERI_*`, `AE_cum`) need no explicit statement, matching their
  implicit-0 default, same as the original's own `$INIT` values for those.
- (b) the 9 duplicated names were removed from `$CAPTURE` (mrgsolve always
  reports every compartment's state regardless of `$CAPTURE`, confirmed
  via `/model_manifest`'s `outputPaths`, which still lists all 19
  compartments).
- (c) the six colliding baseline params were renamed `IL80`/`CRP0`/`Eos0`/
  `FEV10`/`Emph0`/`PVR0` (dropping the underscore, same convention as
  #60's fix) and every read site inside `$ODE` updated to match (3 sites:
  `dxdt_NE_sput`, `dxdt_CRP`, `inflam_penalty` for `IL80`; `FEV1_target`,
  `FEV1_frac` for `FEV10`; `dxdt_PVR` x2 for `PVR0`) — this frees the
  `<CMT>_0` symbols for their mrgsolve-reserved initial-value use. The
  R-side `dose_response` block's `Eos_0=350` override was updated to
  `Eos0=350` to match (still a no-op on the simulated trajectory, per the
  dead-parameter note above — reproduced faithfully, not fixed).

The checked-in original (`copd_mrgsolve_model.R`) was left untouched and
still carries all three defects.

## Verification

**Method.** Both the original's own model code (with the three
build-compat fixes above applied to an in-memory-only scratch copy, so it
would compile at all — never applied to the checked-in
`copd_mrgsolve_model.R`) and `copd_mrgsolve_model_refactored.R`'s embedded
DSL were extracted as bare mrgsolve DSL text and run through the qspserver
`mrgsolve_api` container (`POST /model_manifest`, `POST /run_simulation`)
at `http://localhost:8007`, spaced ~2s apart per request. Both sides were
confirmed to compile cleanly via `/model_manifest` before any simulation
was run.

All six of the original file's own dosing scenarios were run, exactly as
defined in its own `ev()`/`scenarios` list (dose amounts, `ii`, `addl`,
`time`, and `DOSE_*` flags all copied verbatim; compartment numbers
unchanged, since no compartment was added, dropped, or reordered), over
the model's own full 1-year horizon (`end = 8760`, `delta = 24`, 366-370
points depending on scenario — no shortening needed, no solver step-count
issues encountered):

1. **Scenario 1** — Placebo (Natural History), no dosing, all `DOSE_*=0`
2. **Scenario 2** — LAMA Monotherapy (Tiotropium 18µg qd — UPLIFT)
3. **Scenario 3** — LABA/LAMA Dual (Salmeterol + Tiotropium — FLAME/POET-COPD)
4. **Scenario 4** — ICS/LABA Dual (Budesonide + Salmeterol — TORCH)
5. **Scenario 5** — Triple Therapy (LAMA/LABA/ICS — IMPACT/ETHOS)
6. **Scenario 6** — Triple + Roflumilast (all four compounds active)

Plus one extra check: ICS monotherapy with the R script's own
`dose_response` block's `Eos_0=350`/`Eos0=350` override (confirmed a
no-op on both sides, per the dead-parameter note above, so this is
effectively a repeat of the default-`Eos` ICS-alone case — included for
completeness, not additional branch coverage).

Every shared output was compared point-by-point across the full time
grid for all runs: `Cp_LAMA/C_LAMA_PLASMA`, `Cp_LABA/C_LABA_PLASMA`,
`Cp_ICS/C_ICS_PLASMA`, `Cp_PDE4i/C_PDE4I`, `LAMA_conc_lung/C_LAMA`,
`LABA_conc_lung/C_LABA`, `ICS_conc_lung/C_ICS`, `IL8`, `NE_sput`, `CRP`,
`Eos`, `FEV1`, `Emph`, `PVR`, `AE_cum`, `AE_rate_ann`, `mPAP`,
`GOLD_stage`, `CAT_approx`, `SpO2`, `E_LAMA/EFFECT_LAMA`,
`E_LABA/EFFECT_LABA`, `E_ICS/EFFECT_ICS`, `E_PDE4i/EFFECT_PDE4I`, `E_BD`.

**Result: exact match, max absolute diff = 0.0 and max relative diff = 0.0
for every output, every scenario (all six named scenarios plus the extra
check), no solver timeouts or step-count issues at the model's own full
1-year horizon.** This is the expected outcome per the guide's tolerance
rule for pure structural reorganization with no Hill-fitting — every
renamed effect term computes the identical arithmetic on the identical
(renamed) inputs, so there was no source of numerical divergence to begin
with.

## `$PARAM` vs `$GLOBAL` for `C_<STEM>`/`EFFECT_<STEM>`

Per the guide's qspserver compatibility requirement #2, these should
ideally live in `$PARAM` (with a `= 0` default) for direct
`/model_manifest` discoverability. This was not done here, for the same
reason already documented in several sibling refactors (AMD, membranous
nephropathy, breast-cancer, sepsis, SAH): mrgsolve 2.0.1 compiles `$PARAM`
members as **read-only references** inside `$ODE`, so a value that must be
recomputed every timestep from state (`C_LAMA`, `C_LABA`, `C_ICS`,
`C_PDE4I`, `EFFECT_LAMA`, `EFFECT_LABA`, `EFFECT_ICS`, `EFFECT_PDE4I` —
all of which read a compartment) cannot also be declared in `$PARAM`.
Independently re-confirmed against the live `mrgsolve_api` container with
a minimal reproduction targeted at this exact file (`$PARAM C_X = 0.0;
$ODE C_X = CENT/V1;` fails to build with `error: assignment of read-only
reference 'C_X'`). Instead, all eleven (`C_LAMA`, `C_LAMA_PLASMA`,
`EFFECT_LAMA`, `C_LABA`, `C_LABA_PLASMA`, `EFFECT_LABA`, `C_ICS`,
`C_ICS_PLASMA`, `EFFECT_ICS`, `C_PDE4I`, `EFFECT_PDE4I`) are predeclared
as `double`s in `$GLOBAL` and listed in `$CAPTURE`: visible in every
simulation's output columns and in `/model_manifest`'s `outputPaths`, just
not in its `parameters` list.

## Anything else worth flagging

- The doc-comment header's compartment count was corrected from "26 ODE
  compartments" (already inaccurate in the original, which actually has 19)
  to "19 ODE states" — a cosmetic comment/print-statement correction, not a
  gated/verified output, same treatment as the compartment-count fix in
  `breast-cancer/bc_refactor_notes.md`.
- Compartment order and 1-based indices are fully preserved (no
  compartment added, dropped, or reordered), so no dosing-event
  compartment-number remapping was needed anywhere in the verification or
  the R-side scenario code.
- `k_FEV1_dec` (0.000055) and `k_emph_prog`/`k_PVR_prog` (0.000020/
  0.000015) are outside the scope of this refactor (disease-side natural
  history rates, not owned by any of the four compounds) and were left
  exactly as the original had them.

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`, the
`copd | ICS lung`, `copd | Laba lung`, `copd | Lama lung`, and
`copd | PDE4I` rows.

## Discoverability fix

A corpus-wide discoverability audit found `C_ICS`, `C_LABA`, `C_LAMA`, and
`C_PDE4I` were not written as a single contiguous `double C_<STEM> =
<expr>;` statement anywhere in the file. Each was individually
`$GLOBAL`-predeclared (`double C_LAMA;`, `double C_LABA;`, `double C_ICS;`,
`double C_PDE4I;`, each on its own line) and bare-assigned once in `$ODE`
— a legitimate, working pattern (not a bug), but not literal-text-
discoverable by tooling that regexes for `double C_<STEM> = ...;`. There
was no existing `$TABLE`-side reassignment for any of the four.
`C_LAMA_PLASMA`/`C_LABA_PLASMA`/`C_ICS_PLASMA` and the `EFFECT_<STEM>`
predeclares are untouched — out of scope for this fix.

**Fix applied**, identical for all four compounds:
1. Removed the four single-name `double C_<STEM>;` lines from `$GLOBAL`.
2. Added a new block of four `double C_<STEM> = <expr>;` lines to `$TABLE`,
   immediately before `$CAPTURE`, reusing the exact formulas already used
   in `$ODE` verbatim: `C_LAMA = GUT_LAMA/(F_LUNG_LAMA*200.0)`,
   `C_LABA = GUT_LABA/(F_LUNG_LABA*200.0)`, `C_ICS = GUT_ICS/(F_LUNG_ICS*
   200.0)`, `C_PDE4I = CENT_PDE4I/V1_PDE4I`.

**Verification:** `Rscript -e 'parse("copd_mrgsolve_model_refactored.R")'`
succeeds with no error. Extracted DSL posted to qspserver's `mrgsolve_api`
`/model_manifest` compiled cleanly with all four `C_<STEM>` names listed
in `outputPaths` and all four `EC50_<STEM>` in the parameter manifest.

**Result: exact match at every real timestep, one disclosed synthetic-row
divergence at t=0.** `/run_simulation` reproducing this file's own
scenario 6 ("Triple + Roflumilast", `c(dose_LAMA, dose_LABA, dose_ICS,
dose_PDE4i)`, `DOSE_LAMA=DOSE_LABA=DOSE_ICS=DOSE_PDE4I=1`) over 0–48 h
(delta 2), and separately scenario 2 ("LAMA Monotherapy" alone, and a q24h
x3 repeat-dose LAMA-only run to 50 h), run against both the pre-edit
(`git show HEAD:...`) and post-edit DSL, found:
- Every real timestep (`t>=1` onward, including every subsequent q24h
  repeat dose in the multi-cycle run) is **bit-identical** for `C_LAMA`,
  `C_LABA`, `C_ICS`, `C_PDE4I`, `EFFECT_LAMA`, `EFFECT_LABA`, `EFFECT_ICS`,
  `EFFECT_PDE4I`, `FEV1`, and `CAT_approx` (max abs diff = 0).
- At `t=0`, mrgsolve's own duplicate-row reporting (an implicit pre-dose
  baseline observation and the explicit `time=0` dose record sharing the
  same nominal time — this is the guide's documented "dose-instant
  reporting artifact") produces one synthetic extra row per compound's own
  dose event. On that one row, the **pre-edit** file shows a stale
  pre-dose value (`C_LAMA=0` where the true post-dose value is `90`, etc.
  — reproduced for `C_LABA`/`C_ICS`/`C_PDE4I` too, each on the duplicate
  row coincident with its own dose), because `C_<STEM>` in the original is
  only ever written inside `$ODE`, which for this row is evaluated before
  that row's own dose has been applied to the state — even though it is
  `$GLOBAL`-declared, matching the guide's documented anti-artifact
  convention. **This artifact already exists in the unmodified original**
  (confirmed directly against `git show HEAD:...`, not introduced by this
  fix) — the fix here happens to eliminate it as a side effect, because the
  new `$TABLE`-side recompute runs after `$ODE` has already applied the
  dose for that row, so it reports the correct, non-stale value one row
  earlier than the original does. This is a genuine, disclosed **numeric
  difference from the original at exactly that one synthetic duplicate row
  per dose**, not floating-point noise — flagged per the guide's "report
  the mismatch, don't adjust the comparison to make it pass" instruction.
  It is confined entirely to this synthetic zero-duration duplicate row
  (confirmed only 1 mismatched point out of 52 across a 3-cycle repeat-dose
  run — it does not recur at the `t=24`/`t=48` repeat doses, only at the
  very first `t=0` dose, where mrgsolve's baseline-observation-vs-dose-
  record duplication actually occurs) and self-heals by the very next row;
  it does not reflect any change to the model's integrated ODE trajectory.
  Contrast with `cluster-headache` (see that model's own notes, fixed in
  the same batch): its seven compounds' exposed concentrations are all
  central-compartment ratios reached only after `KA_<STEM>`-mediated
  absorption, so no compound is dosed directly into the same compartment
  its `C_<STEM>` reads, and no such divergence appears there at all.
  `C_PERT` in `chronic-pancreatitis` (same batch) is `$GLOBAL`-only in both
  old and new (no `$TABLE`-side reassignment was added there beyond
  upgrading the existing bare line to `double`), so it does not exhibit
  this either.
