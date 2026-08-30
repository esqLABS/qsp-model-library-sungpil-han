# Refactor notes — `gcd_mrgsolve_model.R`

**Scope of this pass.** Per the fork's PK/PD refactor spec
([`FORK_WORKFLOW_GUIDE.md`](../FORK_WORKFLOW_GUIDE.md), Part 2), only the
four compounds flagged in
[`driver-patches/data/compound_perturbation_census.md`](../driver-patches/data/compound_perturbation_census.md)
as "Redirect concentration (clean single site)" were rewritten: **ERT**
(imiglucerase/velaglucerase enzyme-replacement class), **ELIS** (eliglustat),
**MIGS** (miglustat), and **VENG** (venglustat). All four abbreviations were
checked directly against the code before starting: each has its own
dedicated `$CMT`/PK block, its own `USE_<STEM>` scenario flag, and its own
Hill-shaped effect term on GCS/GBA — none is a mislabeled compartment or a
process description masquerading as a compound name, so no census-row
identity correction was needed (only the blank Target/Pathway/notes columns
were filled in). No other compound or mechanism exists in this file (the
disease side — GBA/glucocerebroside/biomarker/organ/hematology/bone/
inflammation submodel — has no drug of its own), so nothing else was
touched.

## Archetype determination

### ERT — Archetype 2 (two-compartment, linear, bidirectional Q/V), bespoke exposed-concentration site

The original already parametrizes ERT's PK in the corpus's preferred
`CL`/`Q`/`V1`/`V2` macro-constant form (`CL_ERT`, `V1_ERT`, `Q_ERT`,
`V2_ERT`), converted internally to micro-constants
(`k10_ert=(CL_ERT*BW)/(V1_ERT*BW)`, `k12_ert=Q_ERT/(V1_ERT*BW)`,
`k21_ert=Q_ERT/(V2_ERT*BW)`) with genuine bidirectional exchange (`ERT_T`
flows back into `ERT_C` via `k21_ert`) — this is Archetype 2, not the
guide's "None of these fit" bespoke case, and not laronidase's one-way
irreversible-transfer shape from `mucopolysaccharidosis-type-1`. The
`BW`-normalized micro-constant algebra is preserved exactly as written
(renamed only) — not "fixed" to a cleaner form, since that would change
numeric behavior the original author chose (even though `k10_ert`'s `BW`
cancels algebraically while `k12_ert`/`k21_ert`'s does not, an asymmetry
that looks like it might be an authoring quirk but is not this refactor's
place to correct).

**The bespoke part:** the effect term reads the *peripheral/tissue* pool
(`ERT_T`, renamed `PERI_ERT`), not the central plasma compartment
(`ERT_C`/`CENT_ERT`) — `ERT_effect = USE_ERT*EMAX_ERT*VELA_MOD*ERT_T/
(EC50_ERT+ERT_T)`. This reflects the real pharmacology (M6P-receptor-
mediated cellular uptake into macrophages is what drives GBA restoration,
not circulating plasma levels) and is the same rationale documented for
laronidase's `C_LARO = TISSUE_LARO` in
`mucopolysaccharidosis-type-1/mps1_refactor_notes.md` — so `C_ERT` is
defined as `PERI_ERT` directly, the single point where an external
covariate would substitute in. `CENT_ERT` (plasma) is not separately
exposed as a second named concentration (unlike laronidase's ADA-exposure
signal) because nothing else in this file reads plasma ERT levels — the
original's `ERT_Cplasma` (`= CENT_ERT * BW`, informational "U total") is
purely a reporting variable, not a PD input, so it was left as an
unrenamed ancillary `$TABLE` variable (its one reference to the renamed
compartment was updated).

### ELIS, MIGS, VENG — Archetype 3 minus peripheral compartment (identical shape, three independent compounds)

All three are oral, one-compartment-plus-depot, first-order
absorption/elimination, no peripheral compartment: `dxdt_GUT = -KA*GUT;
dxdt_CENT = KA*GUT*F/V1 - (CL/V1)*CENT`. As in
`mucopolysaccharidosis-type-1/mps1_refactor_notes.md`'s genistein and
`aneurysmal-subarachnoid-hemorrhage/sah_refactor_notes.md`'s nimodipine,
each central compartment (`ELIS_C`, `MIGS_C`, `VENG_C`) already stores a
**concentration** directly (μg/mL, dosed via `F/V` and eliminated via
`CL/V` with no separate amount-to-concentration division downstream) —
preserved unchanged; `C_ELIS`/`C_MIGS`/`C_VENG` are straight aliases
(`= CENT_ELIS` etc.), not new computations.

Each compound's effect on glucosylceramide synthase (GCS) was already an
explicit `Emax * C/(IC50+C)` ratio (implicit Hill coefficient of 1) — the
guide's "already this shape: this is a rename, not a refit" case for all
three, no fitting performed or needed.

## Renaming applied

| Original | Refactored | Value | Role |
|---|---|---|---|
| `ERT_C` (cmt) | `CENT_ERT` | — | central (plasma), stored per-kg [U/kg] |
| `ERT_T` (cmt) | `PERI_ERT` | — | peripheral / M6P-receptor tissue-macrophage pool |
| `CL_ERT`, `V1_ERT`, `Q_ERT`, `V2_ERT` | unchanged | 1.4, 0.18, 0.35, 0.55 | already conformant to convention |
| `KM6P` | unchanged | 0.006 | **unused in original `$ODE`/`$TABLE`** — dead parameter, preserved as-is, flagged |
| `UPTK_ERT` | unchanged | 0.45 | **unused in original `$ODE`/`$TABLE`** — dead parameter, preserved as-is, flagged |
| `EMAX_ERT` | unchanged | 0.85 | already conformant Hill Emax |
| `EC50_ERT` | unchanged | 0.6 | already conformant Hill EC50 |
| — (none; no explicit Hill term) | `GAMMA_ERT` (new) | 1 | Hill exponent, "original had no explicit Hill term" |
| `VELA_MOD` | unchanged | 1.0 (1.05 for velaglucerase) | non-drug-specific class modifier, not a naming-table role — kept as-is |
| — (was the bare `ERT_T` state) | `C_ERT` (new alias) | — | **the exposed concentration** PD actually reads (= `PERI_ERT`) |
| `ERT_effect` | `EFFECT_ERT` | — | compound's effect on disease (same formula, `pow(C_ERT,GAMMA_ERT)/(...)` form) |
| `ELIS_GUT` (cmt) | `GUT_ELIS` | — | oral depot |
| `ELIS_C` (cmt) | `CENT_ELIS` | — | central (stores concentration directly, μg/mL) |
| `KA_ELIS`, `F_ELIS`, `CL_ELIS` | unchanged | 0.80, 0.20, 38.0 | already conformant |
| `V_ELIS` | `V1_ELIS` | 106 | central volume, renamed to convention |
| `IC50_ELIS` | `EC50_ELIS` | 0.010 | Hill EC50, renamed to convention |
| `IMAX_ELIS` | `EMAX_ELIS` | 0.95 | Hill Emax, renamed to convention |
| — (none) | `GAMMA_ELIS` (new) | 1 | Hill exponent, "original had no explicit Hill term" |
| — (was the bare `ELIS_C` state) | `C_ELIS` (new alias) | — | **the exposed concentration** |
| `ELIS_inh` | `EFFECT_ELIS` | — | compound's effect on disease |
| `MIGS_GUT` (cmt) | `GUT_MIGS` | — | oral depot |
| `MIGS_C` (cmt) | `CENT_MIGS` | — | central (concentration directly, μg/mL) |
| `KA_MIGS`, `F_MIGS`, `CL_MIGS` | unchanged | 0.60, 0.97, 4.5 | already conformant |
| `V_MIGS` | `V1_MIGS` | 28 | central volume, renamed |
| `IC50_MIGS` | `EC50_MIGS` | 50.0 | Hill EC50, renamed |
| `IMAX_MIGS` | `EMAX_MIGS` | 0.80 | Hill Emax, renamed |
| — (none) | `GAMMA_MIGS` (new) | 1 | Hill exponent, "original had no explicit Hill term" |
| — (was the bare `MIGS_C` state) | `C_MIGS` (new alias) | — | **the exposed concentration** |
| `MIGS_inh` | `EFFECT_MIGS` | — | compound's effect on disease |
| `VENG_GUT` (cmt) | `GUT_VENG` | — | oral depot |
| `VENG_C` (cmt) | `CENT_VENG` | — | central (concentration directly, μg/mL) |
| `KA_VENG`, `F_VENG`, `CL_VENG` | unchanged | 0.50, 0.70, 22.0 | already conformant |
| `V_VENG` | `V1_VENG` | 480 | central volume, renamed |
| `Kp_CNS` | unchanged | 0.30 | **unused in original `$ODE`/`$TABLE`** — dead parameter, preserved as-is, flagged |
| `IC50_VENG` | `EC50_VENG` | 0.015 | Hill EC50, renamed |
| `IMAX_VENG` | `EMAX_VENG` | 0.92 | Hill Emax, renamed |
| — (none) | `GAMMA_VENG` (new) | 1 | Hill exponent, "original had no explicit Hill term" |
| — (was the bare `VENG_C` state) | `C_VENG` (new alias) | — | **the exposed concentration** |
| `VENG_inh` | `EFFECT_VENG` | — | compound's effect on disease |
| `GL1_0` (`$PARAM`) | `GL1_SS0` | 8.5 | build-compat rename only, see "Build defect" below — unused dead parameter either way |
| `IL6_0` (`$PARAM`) | `IL6_SS0` | 12 | build-compat rename only, see "Build defect" below — unused dead parameter either way |

`SRT_inh = 1.0 - (1.0-EFFECT_ELIS)*(1.0-EFFECT_MIGS)*(1.0-EFFECT_VENG)` is
the one point the three SRT compounds' effects are combined (unchanged
formula, renamed inputs) — each compound's `EFFECT_<STEM>` stays
independently driveable up to that single combination point, per the
guide's "never collapse several drugs into one shared Hill term" rule.
All disease-side compartments/parameters (`GBA`, `GC_MAC/SP/LV/BM`, `GL1`,
`LYSOGL1`, `CHITR`, `FERRIT`, `SV`, `LV`, `HGB`, `PLT`, `BMD`, `OC`, `OB`,
`IL6`, `NFKB`, and every `k_*` rate constant) are byte-identical to the
original — none is a compound's own PK, confirmed by the diff scope below.

## Build defect found and fixed (disclosed per the guide's settled policy)

**The untouched original does not compile under mrgsolve 2.0.1**, for
three independent reasons — logged as
[`UPSTREAM_ISSUES.md` #69](../translations/UPSTREAM_ISSUES.md):

1. `$CMT` (bare names) and `$INIT` (`name = value`) jointly redeclare all
   26 compartments — `Duplicated model names`.
2. `$CAPTURE` lists 21 of those same 26 compartment names directly
   (everything except the 4 genuine `$TABLE` doubles) — `compartment
   should not be in $CAPTURE`.
3. Once (1) is worked around with the modern `<CMT>_0=` `$MAIN` idiom,
   the unused `$PARAM` names `GL1_0`/`IL6_0` collide with mrgsolve's own
   auto-generated `<CMT>_0` init symbols for the identically-named
   compartments `GL1`/`IL6` — the exact same collision, on the exact same
   `IL6_0` name, already logged for `sepsis` (#30).

**Per the guide's settled policy** ("When the original doesn't compile at
all"), all three fixes were applied directly to the delivered
`gcd_mrgsolve_model_refactored.R` (not just a scratch copy) — syntax-only,
non-numeric, disclosed here and in `UPSTREAM_ISSUES.md` #69:
- `$INIT` dropped, replaced by `$MAIN` using `<CMT>_0 = value;` (same 26
  values).
- `$CAPTURE` trimmed to the 4 genuine `$TABLE` doubles
  (`GC_TOTAL`, `GBA_PCT_NORMAL`, `ERT_Cplasma`, `SRT_GCS_inh`) plus the
  refactor's own new `C_<STEM>`/`EFFECT_<STEM>` doubles (in scope anyway)
  — every renamed compartment is still reported automatically via
  `/model_manifest`'s `outputPaths`, confirmed directly (see below).
- `GL1_0`→`GL1_SS0`, `IL6_0`→`IL6_SS0` (both dead parameters, zero
  numeric effect either way); the R-side virtual-population tibble
  column that also carried the name `GL1_0` was renamed to match so
  `param(mod, ...)` does not try to set a parameter that no longer
  exists.
- The untouched `gcd_mrgsolve_model.R` still carries all three defects
  exactly as written and does not build against mrgsolve 2.0.1.

## Diff scope confirmation

`diff gcd_mrgsolve_model.R gcd_mrgsolve_model_refactored.R` touches
exactly: the `GL1_0`/`IL6_0` rename (with an added explanatory comment),
the ERT/eliglustat/miglustat/venglustat `$PARAM` groups (renames +
`GAMMA_<STEM>` additions + inline archetype/build-fix comments), the 8
renamed `$CMT` lines, the `$INIT`→`$MAIN` block replacement, the four PK
`$ODE` blocks plus the new `C_<STEM>`/`EFFECT_<STEM>` block, the
`$TABLE`/`$CAPTURE` blocks, and five R-side strings (`"ERT_C"`→
`"CENT_ERT"` x1, `"ELIS_GUT"`→`"GUT_ELIS"` x3, `GL1_0`→`GL1_SS0` x2 in the
virtual-population tibble/`init()` call) so the refactored file's own
scenarios and VP loop still dose/read the right compartments by name.
Nothing else — the header comment/citation block, all disease-side
`$PARAM`/`$CMT`/`$ODE` content, all 6 scenario definitions' other
parameter overrides, `run_scenario()`, the plotting helpers, and the
Monte Carlo virtual-population code are otherwise byte-identical.

## Verification

**Method.** Both files' embedded `code <- '...'` DSL blocks were
mechanically extracted (regex on the `gcd_model <- '...'` assignment,
verbatim quoted text) and POSTed to the local qspserver `mrgsolve_api`
service at `http://localhost:8007` (`POST /model_manifest` then `POST
/run_simulation`), ~2s apart per the API's stated concurrency limit. The
build defect above was applied identically, in-memory only, to the
**original's** extracted text (never saved to `gcd_mrgsolve_model.R`) so
it could compile for comparison; the extracted refactored `.cpp` text was
confirmed byte-identical to the `_refactored.R`'s own inline copy before
use, then discarded (no `.cpp` sibling left behind).

Six dosing regimens were run through both models, identical dosing
(compartment index, amount, interval, `addl`), 1-year horizon (`end=8760
h, delta=24 h`, 366-368 timepoints depending on `addl` alignment — no
step-count issue encountered, well under the API's default solver
budget):

1. **S1 — Natural History** (no dosing, all `USE_*=0`): pure disease
   natural-history equations, all four drug compartments at zero.
2. **S2 — Imiglucerase 60 U/kg Q2W** (`USE_ERT=1`, dosing into
   `CENT_ERT`/`ERT_C`, `amt=333.33`, `ii=336h`, `addl=51`): ERT alone.
3. **S4 — Eliglustat 84 mg BID (CYP2D6 EM)** (`USE_ELIS=1, F_ELIS=0.20`,
   dosing into `GUT_ELIS`/`ELIS_GUT`, `amt=84`, `ii=12h`, `addl` for
   730 days): eliglustat alone, extensive-metabolizer parameters.
4. **S5 — Eliglustat 84 mg QD (CYP2D6 PM)** (`USE_ELIS=1, F_ELIS=0.35,
   CL_ELIS=8.0`, `ii=24h`): eliglustat alone, poor-metabolizer
   parameters — exercises the same PK block at different parameter
   values.
5. **S6 — Low-dose ERT (30 U/kg) + Eliglustat (84 mg BID)** (both
   dosed simultaneously): both compounds active together, exercising the
   `EFFECT_ERT`/`EFFECT_ELIS` combination through `GBA_eff`/`GC_clearance`.
6. **Bespoke miglustat/venglustat run** — **no scenario in the original
   R script actually doses miglustat or venglustat**: `USE_MIGS` and
   `USE_VENG` are declared as scenario flags and each compound has a full
   PK/PD block, but none of the file's 6 named `scenarios` list entries
   ever sets either flag to 1 or calls `make_oral_doses()` for
   `MIGS_GUT`/`VENG_GUT`. Rather than inventing a dosing scheme, this
   verification reused the original's own general-purpose
   `make_oral_doses(cmt, dose_mg, freq_h, n_days)` helper (already
   defined and used for eliglustat) together with the file's own stated
   defaults (`DOSE_MIGS=100` mg, comment "TID" → `freq_h=8`; `DOSE_VENG=15`
   mg, comment "QD" → `freq_h=24`; both ×730 days), applied identically to
   both models with `USE_MIGS=1`/`USE_VENG=1`. This is disclosed
   explicitly as a bespoke, non-named-scenario verification, not one of
   the file's own six labeled scenarios.

**One dosing simplification, disclosed:** the original's `make_ert_doses()`
helper passes `rate = -2` to mrgsolve's `ev()` (comment: "infuse over 1 h"),
a NONMEM-style special rate flag whose behavior through this API's simple
`DoseSpec` (a positional, non-`$ODE`-aware dosing record) was untested and
not attempted here. Verification instead used `rate = 0` (bolus) for both
the original and the refactored model identically — since both sides
receive exactly the same dosing input, this still isolates and confirms
the ODE reorganization is correct; it does not confirm the original
R script's own `rate=-2` infusion semantics reproduce through this
particular API (an orthogonal, dosing-mechanics question, not a
structural-refactor one).

**Result: exact match in all six runs. Maximum absolute deviation
observed = 0.0 (bit-identical)** across every shared disease-side output
(`GBA`, `GBA_PCT_NORMAL`, `GC_MAC`, `GC_SP`, `GC_LV`, `GC_BM`, `GC_TOTAL`,
`GL1`, `LYSOGL1`, `CHITR`, `FERRIT`, `SV`, `LV`, `HGB`, `PLT`, `BMD`, `OC`,
`OB`, `IL6`, `NFKB`, `ERT_Cplasma`, `SRT_GCS_inh`) and all 8 renamed/mapped
PK outputs (`ERT_C`/`CENT_ERT`, `ERT_T`/`PERI_ERT`, `ELIS_GUT`/`GUT_ELIS`,
`ELIS_C`/`CENT_ELIS`, `MIGS_GUT`/`GUT_MIGS`, `MIGS_C`/`CENT_MIGS`,
`VENG_GUT`/`GUT_VENG`, `VENG_C`/`CENT_VENG`) over the full shared time
grid in every scenario. This is the expected outcome for a pure structural
reorganization with rename-only Hill terms (no fitting performed, none
needed) per the guide's tolerance table. A separate spot-check of the
refactored model's own new `C_ERT`/`EFFECT_ERT` outputs under S2 dosing
confirmed non-trivial, sensibly-shaped values (e.g. `C_ERT≈0.42→0.06 U/kg`
and `EFFECT_ERT≈0.35→0.08` declining between two dosing intervals), as a
sanity check beyond the pure numeric-match requirement.

## qspserver `/model_manifest` discoverability

- `CL_ERT`, `V1_ERT`, `Q_ERT`, `V2_ERT`, `EMAX_ERT`, `EC50_ERT`,
  `GAMMA_ERT`, `KA_ELIS`, `F_ELIS`, `CL_ELIS`, `V1_ELIS`, `EMAX_ELIS`,
  `EC50_ELIS`, `GAMMA_ELIS`, `KA_MIGS`, `F_MIGS`, `CL_MIGS`, `V1_MIGS`,
  `EMAX_MIGS`, `EC50_MIGS`, `GAMMA_MIGS`, `KA_VENG`, `F_VENG`, `CL_VENG`,
  `V1_VENG`, `EMAX_VENG`, `EC50_VENG`, `GAMMA_VENG` are all declared in
  `$PARAM` and confirmed present with their defaults in the live
  `/model_manifest` response (85 total parameters).
- `C_ERT`, `EFFECT_ERT`, `C_ELIS`, `EFFECT_ELIS`, `C_MIGS`, `EFFECT_MIGS`,
  `C_VENG`, `EFFECT_VENG` are all state-derived (recomputed every `$ODE`
  step) and are exposed via `$CAPTURE` only, **not** `$PARAM`. Verified
  directly: a scratch variant adding `C_ERT : 0` to `$PARAM` alongside the
  existing `double C_ERT = PERI_ERT;` in `$ODE` failed to compile under
  mrgsolve 2.0.1 (`error: assignment of read-only reference 'C_ERT'`) —
  confirming the same `$PARAM`-vs-`$ODE`-local mutual exclusion already
  documented in `amd_refactor_notes.md`, `sah_refactor_notes.md`,
  `ted_refactor_notes.md`, and `mps1_refactor_notes.md` for this
  engine/plugin combination. This is a general mrgsolve/autodec behavior,
  not a defect specific to this file, so no additional `UPSTREAM_ISSUES.md`
  entry was added for it beyond noting it here. All eight remain fully
  discoverable via `/run_simulation`'s `outputs` selection through
  `$CAPTURE`, just not via `/model_manifest`'s parameter list.
- `model_content` is pure DSL text extracted from the `code <- '...'`
  R-string wrapper (this file uses the `gcd_model <- '...'` variant of
  that pattern) — no `.cpp` sibling was left behind; the extraction was
  in-memory only, used to build the verification requests above and then
  discarded.

## Anything else flagged

- `KM6P` (M6P receptor Km) and `UPTK_ERT` (tissue uptake rate) are
  declared in the original's `$PARAM` but never referenced anywhere in
  `$ODE`/`$TABLE` — dead/vestigial parameters, likely left over from an
  earlier authoring pass where tissue uptake may have been modeled with
  explicit Michaelis-Menten kinetics before being simplified to the
  current linear `k12_ert`/`k21_ert` exchange. Preserved unchanged, per
  the guide's "never invent or remove a parameter the original had" rule.
- `Kp_CNS` (venglustat brain/plasma ratio) is likewise declared but
  unused — the model has no separate CNS/brain compartment for
  venglustat despite the header comment describing it as "CNS
  penetration" and "high Vd." Preserved unchanged, flagged for whoever
  reviews this file's CNS pharmacology claims next.
- `DOSE_ERT`/`DOSE_ELIS`/`DOSE_MIGS`/`DOSE_VENG` are `$PARAM`-declared but
  are purely informational/scenario-labeling values — the actual dose
  amounts are hardcoded into each scenario's `make_ert_doses()`/
  `make_oral_doses()` call in the R-side scenario list, not read from
  these parameters inside the DSL. Left unchanged (not a naming-table
  role, not compound-specific PK).
- No compound in this file needed a bespoke, "none of the archetypes
  fit" structure — ERT is Archetype 2 with a bespoke *exposed-site*
  choice (documented above), and ELIS/MIGS/VENG are all a clean
  Archetype 3 minus peripheral.
