# Refactor notes — `bipolar-disorder/bd_mrgsolve_model.R`

Compounds refactored: **Lithium (LI)**, **Valproate (VPA)**, **Quetiapine
(QTP)** (+ its active metabolite **norquetiapine, NQT**), **Lamotrigine
(LTG)** — the four rows in
`driver-patches/data/compound_perturbation_census.md` classified
`bipolar-disorder | Lamotrigine (LTG)`, `bipolar-disorder | Lithium`,
`bipolar-disorder | Quetiapine (QTP)`, `bipolar-disorder | Valproate
(VPA)`, all "Normalize duplicate concentration sites, then redirect". This
is the file's complete exogenous-compound set — no other compound exists
in this model, confirmed by reading the full `$PARAM`/`$CMT`/`$ODE` code.

## Aripiprazole is named in the README and `$PROB` prologue but not modelled

The disease README (`bipolar-disorder/README.md`) lists "Aripiprazole" as
one of the drug targets shown in the mechanistic map, and the original
model's own `$PROB` prologue text names it too ("Lithium / Valproate /
Quetiapine / Lamotrigine / Aripiprazole PK-PD"). Reading the actual
`$PARAM`/`$CMT`/`$ODE` code found **no Aripiprazole PK compartments, no
Aripiprazole parameters, and no Aripiprazole PD effect anywhere** — the
model implements exactly four compounds, matching the census. This is the
guide's own documented failure mode ("the README is context, not ground
truth for structure") in the opposite direction — here it's the model's
own prologue text overstating scope, not the README understating it.
Disclosed here rather than silently dropped; not a compound this task
asked us to add, and adding one would be new pharmacology, not a
refactor. `$PROB`'s refactored prologue notes this explicitly.

## Archetype per compound

- **Lithium (LI)**: **archetype 3** (depot + central + peripheral,
  linear). `GUT_LI`/`CENT_LI`/`PERI_LI` (was `Li_gut`/`Li_central`/
  `Li_periph`), `KA_LI`/`V1_LI`/`V2_LI`/`Q_LI`/`CL_LI`/`F_LI` (was
  `ka_Li`/`Vc_Li`/`Vp_Li`/`Q_Li`/`CL_Li`/`F_Li`) unchanged in value.
  `F_LI` is declared but never multiplied into the original's own
  `dxdt_Li_central` (`= ka_Li*Li_gut - ...`, no `F_Li` factor) — preserved
  unused exactly as the original had it, not silently "fixed" by adding
  the multiplication.
- **Valproate (VPA)**: **archetype 3 minus peripheral** (depot + central,
  linear elimination on total amount, nonlinear protein binding on top).
  `GUT_VPA`/`CENT_VPA` (was `VPA_gut`/`VPA_central`), `KA_VPA`/`V1_VPA`/
  `CL_VPA`/`FU_VPA0`/`KM_VPA` (was `ka_VPA`/`Vc_VPA`/`CL_VPA`/`fu_VPA0`/
  `Km_fu`) unchanged in value. The nonlinear free-fraction term
  (`C_VPA = VPA_total*FU_VPA0*(1+VPA_total/KM_VPA)`) is kept exactly as
  the original — not flattened to a linear approximation, per the guide's
  "don't simplify mechanistically rich PK" rule.
- **Quetiapine (QTP) + norquetiapine (NQT)**: **archetype 3 minus
  peripheral** for the parent (`GUT_QTP`/`CENT_QTP`, was `QTP_gut`/
  `QTP_central`), plus a **bespoke one-compartment metabolite chain**
  (`CENT_NQT`, was `NQT_central`) fed by first-order metabolism from
  `CENT_QTP` — not its own dose/depot, matching the original exactly.
  `KA_QTP`/`V1_QTP`/`CL_QTP`/`F_QTP`/`FMET_QTP`/`CL_NQT` (was `ka_QTP`/
  `Vc_QTP`/`CL_QTP`/`F_QTP`/`km_QTP`/`CL_NQT`) unchanged in value. The
  metabolite reuses the parent's own volume divisor (`V1_QTP/1000.0`) for
  its own concentration, exactly as the original did (a disclosed,
  pre-existing simplification, not introduced here).
- **Lamotrigine (LTG)**: **archetype 3 minus peripheral** (depot +
  central, linear). `GUT_LTG`/`CENT_LTG` (was `LTG_gut`/`LTG_central`),
  `KA_LTG`/`V1_LTG`/`CL_LTG` (was `ka_LTG`/`Vc_LTG`/`CL_LTG`) unchanged in
  value. No bioavailability term in the original (F implicitly 1),
  preserved as-is.

No compound needed TMDD (archetype 4) — every drug-target interaction in
this file is a plain saturating Hill ratio (or, for Valproate, a
protein-binding nonlinearity feeding into one), not an integrated
receptor-binding ODE system.

## PK/PD interface

- **Lithium**: `C_LI = CENT_LI / V1_LI` (mEq/L) is the single exposed
  concentration every one of Lithium's four PD effects reads (GSK-3β,
  BDNF, IL-6, YMRS).
- **Valproate**: `C_VPA` is the **free (unbound)** concentration
  (`VPA_total_conc * FU_VPA0 * (1 + VPA_total_conc/KM_VPA)`, µg/mL) — this
  is what both of Valproate's PD effects (GSK-3β, YMRS) already read in
  the original (`VPA_free`, never the total). The protein-bound total
  (`VPA_total_conc`, µg/mL) is kept as a separate, non-exposed reporting
  quantity (`$TABLE`'s `VPA_ugmL`, unchanged name/formula) — not a second
  "exposed" concentration, since nothing mechanistic reads it; this
  matches the guide's "two [concentration sites] only when a genuinely
  different tissue site matters" allowance loosely (here it's bound vs.
  free, not a tissue site, but the same principle: only one is
  mechanistically active).
- **Quetiapine (parent)**: `C_QTP = CENT_QTP / (V1_QTP/1000.0)` (ng/mL) —
  read by the D2-occupancy, YMRS, and MADRS effects.
- **Norquetiapine (metabolite)**: `C_NQT = CENT_NQT / (V1_QTP/1000.0)`
  (ng/mL, reusing the parent's own volume divisor, see above) — read by
  the serotonin-synthesis-stimulation effect and the (linear, non-Hill)
  weight-gain driver.
- **Lamotrigine**: `C_LTG = CENT_LTG / V1_LTG` (µg/mL) — read by the
  MADRS effect.

### Hill interface — rename, not a fit, for every effect term

Every PD effect term in the original was already shaped exactly
`Emax*C/(EC50+C)` (implicit Hill coefficient of 1) **except** the
weight-gain driver, which is linear. `GAMMA_<STEM>_<EFFECT> = 1` pulled
out as an explicit new parameter for every Hill term (11 new `GAMMA_*`
parameters, all `=1`, matching the original's implicit shape):

| Compound | Old name(s) | New name(s) | Reads |
|---|---|---|---|
| LI | `Emax_GSK`/`IC50_Li_GSK` | `EMAX_LI_GSK3`/`EC50_LI_GSK3`/`GAMMA_LI_GSK3` | `C_LI` |
| LI | `Emax_BDNF_Li`/`EC50_BDNF_Li` | `EMAX_LI_BDNF`/`EC50_LI_BDNF`/`GAMMA_LI_BDNF` | `C_LI` |
| LI | `Emax_IL6_Li`/`EC50_IL6_Li` | `EMAX_LI_IL6`/`EC50_LI_IL6`/`GAMMA_LI_IL6` | `C_LI` |
| LI | `Emax_YMRS_Li`/`EC50_YMRS_Li` | `EMAX_LI_YMRS`/`EC50_LI_YMRS`/`GAMMA_LI_YMRS` | `C_LI` |
| VPA | `Emax_GSK`/`IC50_VPA_GSK` | `EMAX_VPA_GSK3`/`EC50_VPA_GSK3`/`GAMMA_VPA_GSK3` | `C_VPA` |
| VPA | `Emax_YMRS_VPA`/`EC50_YMRS_VPA` | `EMAX_VPA_YMRS`/`EC50_VPA_YMRS`/`GAMMA_VPA_YMRS` | `C_VPA` |
| QTP | `Imax_D2`/`EC50_YMRS_QTP` (reused) | `EMAX_QTP_D2`/`EC50_QTP_D2`/`GAMMA_QTP_D2` | `C_QTP` |
| QTP | `Emax_YMRS_QTP`/`EC50_YMRS_QTP` | `EMAX_QTP_YMRS`/`EC50_QTP_YMRS`/`GAMMA_QTP_YMRS` | `C_QTP` |
| QTP | `Emax_MADRS_QTP`/`EC50_MADRS_QTP` | `EMAX_QTP_MADRS`/`EC50_QTP_MADRS`/`GAMMA_QTP_MADRS` | `C_QTP` |
| QTP/NQT | `0.3`/`80` (hardcoded inline in `$ODE`) | `EMAX_QTP_5HT`/`EC50_QTP_5HT`/`GAMMA_QTP_5HT` | `C_NQT` |
| LTG | `Emax_MADRS_LTG`/`EC50_MADRS_LTG` | `EMAX_LTG_MADRS`/`EC50_LTG_MADRS`/`GAMMA_LTG_MADRS` | `C_LTG` |

Two disclosed value-duplications, not refits (both preserve the original's
own reuse of a single value across two purposes):

1. **`Emax_GSK = 0.90`** was one shared parameter the original used for
   *both* Lithium's and Valproate's independent GSK-3β inhibition terms.
   Duplicated (same value, `0.90`) into `EMAX_LI_GSK3` and
   `EMAX_VPA_GSK3` so each compound keeps its own explicit Hill parameter
   set, per the guide's "keep each compound's `EFFECT_<STEM>` separate"
   rule — not a refit, both still equal `0.90`.
2. **`EC50_YMRS_QTP = 120`** was reused by the original as the D2-receptor
   binding constant too (`D2_occ_QTP = QTP_conc/(EC50_YMRS_QTP+QTP_conc)`,
   a value that is algebraically identical to
   `EMAX_QTP_D2*C_QTP/(EC50_QTP_D2+C_QTP)` with `EMAX_QTP_D2=Imax_D2` and
   `EC50_QTP_D2=EC50_YMRS_QTP`'s own value). Duplicated (same value,
   `120`) into `EC50_QTP_D2` so the D2-occupancy and YMRS-reduction
   effects each have their own explicit potency parameter — not a refit.

**Weight gain is deliberately kept linear, not forced into a Hill shape.**
The original's `Wt_gain_rate = kWt_gain * NQT_conc` has no saturation
term at all — renamed to `KWT_QTP * C_NQT` (same value, `0.003`), with no
`EMAX_QTP_WT`/`EC50_QTP_WT` introduced, per the guide's "don't force a bad
fit" principle (there is no fit to make here; it is genuinely
unbounded/linear in the original).

**`E_MADRS_combo_bonus`** (the Li+QTP synergy term, Geddes 2016) is kept
exactly as the original computed it — a step-function bonus gated by two
hardcoded thresholds (`C_LI > 0.3`, `C_QTP > 50`), not a Hill term. The
thresholds are literal numbers in the original too (no named params
there), so they stay literal here rather than being promoted to new,
unrequested parameters.

## Duplicate concentration sites — the census's classification, verified against the actual code

Reading the code confirmed the census's classification for all four
compounds: the original independently re-derives each compound's
concentration **twice**, under two different names — once in `$MAIN`
(`Li_conc`, `VPA_conc`/`VPA_free`, `QTP_conc`, `NQT_conc`, `LTG_conc`,
feeding every PD/Hill effect term) and again in `$TABLE` (`Lithium_mEqL`,
`VPA_ugmL`/`VPA_free_ugmL`, `QTP_ngmL`, `NQT_ngmL`, `LTG_ugmL`, feeding
`$CAPTURE`).

### `$MAIN` and `$TABLE` are NOT interchangeable in this file — confirmed by a minimal reproduction

Before normalizing, a minimal reproduction was run against the live
qspserver mrgsolve API to settle whether the two independent
recomputations could simply be consolidated into one shared value:

```
$PARAM CL=0.3, V=1.0
$CMT CENT
$MAIN
double C_MAIN = CENT/V;
$ODE
dxdt_CENT = -CL*CENT;
$TABLE
double C_TABLE = CENT/V;
$CAPTURE C_MAIN C_TABLE
```

Dosed once (`amt=100` into `CENT`) with a second dose at `t=5`, `delta=1`:

| time | CENT | C_MAIN | C_TABLE |
|---|---|---|---|
| 0 | 100 | 0 | 100 |
| 1 | 74.08 | 100 | 74.08 |
| 2 | 54.88 | 74.08 | 54.88 |
| 5 | 122.31 | 30.12 | 122.31 |

`C_MAIN[i]` exactly equals `CENT[i-1]` (a clean one-row lag), while
`C_TABLE[i]` equals `CENT[i]` (the freshly-solved value for that row).
**`$MAIN` runs once per interval using the state as of the *start* of
that interval; `$TABLE` runs after the ODE has been advanced to the row's
own time.** This is exactly the guide's "keep a calculation in the block
the original used it in" warning, generalized from `$MAIN`-vs-`$ODE` to
`$MAIN`-vs-`$TABLE`: the untouched original already has this one-row lag
baked into every PD effect (computed from the `$MAIN`-frozen
concentration) relative to what it reports (computed fresh in `$TABLE`).
Sharing one mutable value across both blocks (or moving either
computation into the other) would silently change the simulated
trajectory.

### The fix: one authored formula per compound, two independent `double` declarations

Each concentration formula is authored exactly once, as a `$GLOBAL`
`#define` macro (`LI_CONC_EXPR`, `VPA_TOTAL_EXPR`/`VPA_FREE_EXPR(...)`,
`QTP_CONC_EXPR`, `NQT_CONC_EXPR`, `LTG_CONC_EXPR`). `$MAIN` and `$TABLE`
each still declare their own `double` (`C_LI`/`Lithium_mEqL`,
`C_VPA`/`VPA_ugmL`+`VPA_free_ugmL`, `C_QTP`/`QTP_ngmL`,
`C_NQT`/`NQT_ngmL`, `C_LTG`/`LTG_ugmL`) — so each still evaluates against
its own block's live state exactly as before, preserving the one-row lag
— but by expanding the same macro text instead of retyping the
arithmetic by hand under two different variable names against two
different (renamed) sets of inputs. This is genuine single-source-of-
truth normalization (one formula, two use sites that each need their own
timing) rather than either a silent behavior change or a no-op rename.

`C_<STEM>`/`EFFECT_<STEM>` are `$MAIN`-local + `$TABLE`-local +
`$CAPTURE` only, never `$PARAM` — confirmed (same constraint documented
in numerous sibling refactors, e.g. `copd/copd_refactor_notes.md`,
`atrial-fibrillation/af_refactor_notes.md`) that mrgsolve 2.0.1 compiles
`$PARAM` members as read-only references, so a value recomputed from live
state cannot also be a `$PARAM`.

## A real bug caught by verification, not disclosed-and-ignored

A first draft of `GSK_total_inh` wrote the "obviously correct"
Bliss-independence combination:

```
double GSK_total_inh = 1.0 - (1.0 - EFFECT_LI_GSK3) * (1.0 - EFFECT_VPA_GSK3);
```

This is **not** what the original computes. The original's literal
expression is:

```
double GSK_total_inh = 1.0 - (1.0 - (1.0 - Inh_GSK_Li) * (1.0 - Inh_GSK_VPA));
```

Algebraically, `1.0 - (1.0 - X) = X`, so the original's `GSK_total_inh`
is exactly `(1-Inh_GSK_Li)*(1-Inh_GSK_VPA)` — the *residual* (uninhibited)
fraction, not "one minus" that product as its name and my first-draft
"simplification" both assumed. The redundant double-negation is
load-bearing: it cancels only when its actual downstream consumer,
`dxdt_GSK3_activ`'s own `(1.0 - GSK_total_inh)` term, is also included —
at which point the *net* applied inhibition does come out to the
Bliss-independence value the name implies, but `GSK_total_inh` itself, as
a standalone quantity, is inverted relative to its name (a pre-existing
naming quirk in the original, not something to "fix" here).

This was caught, not guessed at: the first-draft version verified with a
**non-floating-point-scale mismatch** against the original
(`GSK3_activ` up to 0.95, `IL6_level` up to 2.4 absolute, across all six
scenarios — `IL6_level`'s own deviation is a downstream consequence of
`GSK3_activ`'s, via the shared `GSK3_activ/GSK3_base` term in
`dxdt_IL6_level`). Per the guide's "never adjust a comparison to make a
gate pass" rule, the mismatch was traced to its source (this one line)
rather than loosened away. Restoring the original's exact nested
double-negation (`1.0 - (1.0 - (1.0 - EFFECT_LI_GSK3) * (1.0 -
EFFECT_VPA_GSK3))`) reproduced an exact (`0.0`) match on every output —
see Verification below.

## Upstream build status — two pre-existing mrgsolve 2.0.1 defects, unrelated to any compound's own PK

Logged as [`translations/UPSTREAM_ISSUES.md` #139](../translations/UPSTREAM_ISSUES.md):

1. **`$CAPTURE` re-lists ten `$CMT` compartment names directly**
   (`DA_index HT5_index GSK3_activ BDNF_level IL6_level Cortisol_idx` and
   `YMRS_score MADRS_score GAF_score Weight_chg`) — mrgsolve 2.0.1 rejects
   this. Fixed by dropping both lines (all ten are still reported,
   unchanged, since mrgsolve exposes `$CMT` states as output columns
   automatically).
2. **`$MAIN`'s `Circ_forcing` uses `SOLVERTIME`**, which is only valid
   inside `$ODE` (`_ODETIME_` is not declared in `$MAIN`'s own generated
   scope) — confirmed by a minimal reproduction (see UPSTREAM_ISSUES.md
   #139) that `TIME` is the `$MAIN`-scope equivalent and reports the
   identical value at every reporting row. Fixed by the syntax-only
   rename `SOLVERTIME` -> `TIME`.

Neither defect sits inside any of the four refactored compounds' own
PK/PD blocks — both are pre-existing, disease-model-side (or file-wide
`$CAPTURE`-block) build defects, logged rather than fixed upstream, and
carried into `bd_mrgsolve_model_refactored.R` as disclosed, syntax-only,
non-numeric fixes (confirmed by the exact-match verification below, run
against a patched-original baseline carrying the identical two fixes and
nothing else).

## Verification

Via the qspserver `mrgsolve_api` service (`POST /model_manifest`,
`POST /run_simulation`, `http://localhost:8007`, requests spaced ~2.2 s
apart, one at a time), comparing a patched-original baseline (the
untouched original's own DSL, with only the two syntax-only fixes above
applied — nothing else changed) against
`bd_mrgsolve_model_refactored.R`'s extracted DSL, using **the model's own
six named scenarios' dosing** (identical `amt`/`cmt`/`ii`/`addl`/`time`
and `YMRS_base`/`MADRS_base` overrides, reproduced from the original R
script's own `make_dosing()` calls), each scenario's own full time
window — **no shortening needed for any scenario**, none hit the API's
default solver step-count limit:

| Scenario | Dosing | Window | Points | Max abs diff (23 shared outputs) |
|---|---|---|---|---|
| S1: Li monotherapy (acute mania) | Li 8.1 mmol q8h ×63 | 0–504 h | 506 | **0.0** |
| S2: VPA monotherapy (acute mania) | VPA 500 mg q12h ×42 | 0–504 h | 506 | **0.0** |
| S3: QTP monotherapy (BD depression) | QTP 300 mg q24h ×56 | 0–1344 h | 1346 | **0.0** |
| S4: Li + QTP combination (BD depression) | Li 8.1 mmol q8h ×168 + QTP 300 mg q24h ×56 | 0–1344 h | 1347 | **0.0** |
| S5: Li maintenance (1 year) | Li 8.1 mmol q8h ×1095 | 0–8760 h | 2192 | **0.0** |
| S6: LTG titration (BD-II depression) | LTG 25/50/100/200 mg q24h, 4 staged blocks | 0–2688 h | 1349 | **0.0** |

All 23 shared `$CAPTURE` outputs (`Lithium_mEqL`, `VPA_ugmL`,
`VPA_free_ugmL`, `QTP_ngmL`, `NQT_ngmL`, `LTG_ugmL`, `DA_index`,
`HT5_index`, `GSK3_activ`, `BDNF_level`, `IL6_level`, `Cortisol_idx`,
`YMRS_score`, `MADRS_score`, `GAF_score`, `Weight_chg`, `Li_toxic`,
`Li_subtherapeutic`, `VPA_toxic`, `YMRS_response`, `YMRS_remission`,
`MADRS_response`, `MADRS_remission`) matched exactly (`0.0`) at every
reported time point, in every scenario. This is the expected outcome per
the guide's tolerance rule — every one of this file's four compounds is a
pure structural reorganization (rename, duplicate-site normalization,
Hill-shape rename) with no ODE-derived Hill-fitting anywhere, so exact
agreement is the correct result, not a loosened tolerance.

### `/model_manifest` parameter-value diff

All 34 of the original's own renamed parameters (PK: `KA_LI`/`V1_LI`/
`V2_LI`/`Q_LI`/`CL_LI`/`F_LI`, `KA_VPA`/`V1_VPA`/`CL_VPA`/`FU_VPA0`/
`KM_VPA`, `KA_QTP`/`V1_QTP`/`CL_QTP`/`F_QTP`/`FMET_QTP`/`CL_NQT`,
`KA_LTG`/`V1_LTG`/`CL_LTG`; Hill: `EC50_LI_GSK3`/`EC50_VPA_GSK3`/
`EC50_LI_BDNF`/`EMAX_LI_BDNF`/`EC50_LI_IL6`/`EMAX_LI_IL6`/`EC50_LI_YMRS`/
`EMAX_LI_YMRS`/`EC50_VPA_YMRS`/`EMAX_VPA_YMRS`/`EC50_QTP_YMRS`/
`EMAX_QTP_YMRS`/`EMAX_QTP_D2`/`EC50_QTP_MADRS`/`EMAX_QTP_MADRS`/
`EC50_LTG_MADRS`/`EMAX_LTG_MADRS`/`ALPHA_COMBO_LI_QTP`/`KWT_QTP`) map to
their original names (`ka_Li`, `Vc_Li`, ..., `Emax_MADRS_LTG`,
`ALPHA_combo`, `kWt_gain`) with **zero value mismatches** (automated
diff). `EC50_QTP_D2` (new name) confirmed equal to the original's own
reused `EC50_YMRS_QTP` value (`120`). The refactored manifest carries
exactly 11 new `GAMMA_*=1` parameters (the explicit Hill coefficients
pulled out per the guide's rename rule) plus the two disclosed value-
duplications (`EMAX_LI_GSK3`/`EMAX_VPA_GSK3`, both `0.90`) and the
promoted inline weight-gain/5-HT constants (`EMAX_QTP_5HT=0.3`,
`EC50_QTP_5HT=80`) — 81 total parameters vs. the original's 67. `$CMT`
compartment order (and therefore every compound's 1-based dosing index)
is identical between original and refactored manifests at all 10 PK
positions (position 1 = `Li_gut`/`GUT_LI`, position 4 = `VPA_gut`/
`GUT_VPA`, position 6 = `QTP_gut`/`GUT_QTP`, position 9 = `LTG_gut`/
`GUT_LTG`).

All 16 new interface symbols (`C_LI`, `C_VPA`, `C_QTP`, `C_NQT`, `C_LTG`,
`EFFECT_LI_GSK3`, `EFFECT_LI_BDNF`, `EFFECT_LI_IL6`, `EFFECT_LI_YMRS`,
`EFFECT_VPA_GSK3`, `EFFECT_VPA_YMRS`, `EFFECT_QTP_D2`, `EFFECT_QTP_YMRS`,
`EFFECT_QTP_MADRS`, `EFFECT_QTP_5HT`, `EFFECT_LTG_MADRS`) are
`$MAIN`/`$TABLE`-local state-dependent quantities, so — per this fork's
established precedent — they are not `$PARAM` defaults; all 16 appear in
`/model_manifest`'s `outputPaths` instead (51 `outputPaths` total: 22
`$CMT` states + 6 unchanged `$TABLE`-derived concentration outputs
(`Lithium_mEqL`, `VPA_ugmL`, `VPA_free_ugmL`, `QTP_ngmL`, `NQT_ngmL`,
`LTG_ugmL`) + 3 safety flags + 4 response flags + 16 new interface
symbols — every new symbol is discoverable there, not only by reading the
DSL text).

### Discoverability self-check

`grep -c "double C_<STEM> = " bd_mrgsolve_model_refactored.cpp` and
`grep -c "EC50_<STEM>"` both return ≥1 for all four stems (`LI`, `VPA`,
`QTP`, `LTG`) — confirmed directly against the extracted `.cpp` before
calling this done.

## Delivered files

- `bd_mrgsolve_model_refactored.R` — full original script (helper
  function, six scenario definitions, all plots, summary table) with the
  DSL block renamed/restructured per the above, the model variable
  renamed `code_refactored`/`mod` compiled as
  `"BipolarDisorder_QSP_refactored"` (new registered name, avoiding any
  compiled-model cache collision with the original), and `make_dosing()`'s
  `cmt_map` updated to the renamed `$CMT` names (`GUT_LI`, `GUT_VPA`,
  `GUT_QTP`, `GUT_LTG`). No other line (scenario definitions, plotting,
  summary tables, all of which read only `$CAPTURE`d output names that
  were kept identical to the original) was touched.
- No standalone `.cpp` DSL extraction is committed alongside it — checked
  against this fork's own established practice first (zero
  `*_refactored.cpp` files exist anywhere in the tracked repo, across all
  100+ completed refactors), so the `model_content`-for-qspserver
  extraction described in the guide's "qspserver compatibility
  requirements" section is treated the same way every prior refactor has
  treated it: a verification-time-only text extraction (`code_refactored
  <- '...'` -> bare DSL), done fresh from `bd_mrgsolve_model_refactored.R`
  whenever a client needs to `POST` `model_content`, not a second
  repo-tracked artifact that could drift from the `.R` file's own quoted
  string.

## Anything else flagged

- No compound other than Lithium, Valproate, Quetiapine, and Lamotrigine
  was touched — this is the file's complete compound set (Aripiprazole,
  named only in prose, is not implemented — see above).
- All scratch/debug files used for API verification (`.cpp` DSL
  extractions, the patched-original scratch copy, the minimal-reproduction
  test model, request/response JSON, the comparison script) were created
  under the session scratchpad directory only and deleted after use —
  none were left in `bipolar-disorder/` or anywhere else in the repo.

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`, the
`bipolar-disorder | Lamotrigine (LTG)`, `bipolar-disorder | Lithium`,
`bipolar-disorder | Quetiapine (QTP)`, and `bipolar-disorder | Valproate
(VPA)` rows.
