# Refactor notes — `chronic-myeloid-leukemia/cml_mrgsolve_model.R`

Scope: all **four** TKI compounds this file models — **Asciminib (ASC)**,
**Dasatinib (DAS)**, **Imatinib (IMT)**, **Nilotinib (NIL)** — per the four
existing rows in `driver-patches/data/compound_perturbation_census.md`, all
classified "Redirect concentration (clean single site)". There is no other
pharmacological compound in this file (the immune/NK/CTL kill terms and the
Michor-hierarchy hematopoiesis/leukemia dynamics are not drug-driven) —
every equation outside the four TKIs' own PK/PD blocks is byte-identical to
the original, except the two mechanical build-compatibility fixes below,
which touch the file globally (they were pre-existing and unrelated to any
compound).

## Archetype per compound

All four share a genuine common structural quirk in the original, worth
stating once before the per-compound detail: **PD reads an intracellular
concentration, not plasma**, for three of the four compounds (IMT/DAS/NIL).
The original models a "biophase"/intracellular accumulation compartment fed
one-way from plasma (`Kin`/`Kout`- or `Kic`/`Kout`-style rate constants, no
back-flux to plasma), and every one of those three compounds' effect terms
(`E_imt`, `E_das`, `E_nil`) reads that intracellular state (`Cic_*`), never
the plasma one (`Cp_*`). This is exactly the guide's "two [concentration
variables] only when a genuinely different tissue site matters" carve-out —
so each of these three compounds gets **two** named concentration-like
states: `CENT_<STEM>` (plasma, informational) and the actual PD-facing
`C_<STEM>` (aliased to a new bespoke compartment, `IC_<STEM>`, since the
guide's naming table has no ready-made slot for "intracellular biophase").
Asciminib is the exception: its effect term reads `Cp_asc` (plasma) directly
— no intracellular compartment exists for it in the original.

**Imatinib (IMT) — Archetype 3 (depot + central, linear elimination, no
peripheral) plus the bespoke `IC_IMT` biophase compartment:**

```
dxdt_GUT_IMT  = F_IMT * dose_rate_imt - KA_IMT * GUT_IMT;
dxdt_CENT_IMT = KA_IMT * GUT_IMT / V1_IMT - (CL_IMT / V1_IMT) * CENT_IMT;
dxdt_IC_IMT   = KIN_IMT * CENT_IMT - KOUT_IMT * IC_IMT;
```
`C_IMT = IC_IMT` (straight alias — PD reads intracellular).

**Dasatinib (DAS) — no depot; continuous zero-order infusion approximation
straight into `CENT_DAS`, plus the same bespoke `IC_DAS` biophase
compartment.** This is the guide's own "edge case — continuous infusion
input, not an event": the original doses via
`dose_rate_das = dose_das/24/488*1000` added directly into the plasma
ODE, never through `ev()`/a depot; preserved exactly, unconverted to a
bolus-dosed depot (that would change how the model's own scenarios must be
run). `KA_DAS` exists as a declared `$PARAM` in the original but is **never
used anywhere in the original's `$ODE`** — a pre-existing dead parameter,
not something this refactor introduced; carried over verbatim (same value,
same name pattern, flagged inline in `$PARAM` as unused) rather than
silently dropped, per "never invent or default a PK parameter" (dropping
a parameter the original had would be the mirror-image mistake).

```
dxdt_CENT_DAS = F_DAS * dose_rate_das / V1_DAS - (CL_DAS / V1_DAS) * CENT_DAS;
dxdt_IC_DAS   = KIC_DAS * CENT_DAS - KOUT_DAS * IC_DAS;
```
`C_DAS = IC_DAS` (straight alias).

**Nilotinib (NIL) — Archetype 3 (depot + central, linear), plus the same
bespoke `IC_NIL` biophase compartment.** Structurally identical to
imatinib:

```
dxdt_GUT_NIL  = F_NIL * dose_rate_nil - KA_NIL * GUT_NIL;
dxdt_CENT_NIL = KA_NIL * GUT_NIL / V1_NIL - (CL_NIL / V1_NIL) * CENT_NIL;
dxdt_IC_NIL   = KIC_NIL * CENT_NIL - KOUT_NIL * IC_NIL;
```
`C_NIL = IC_NIL` (straight alias).

**Asciminib (ASC) — Archetype 3 (depot + central, linear), no peripheral,
no intracellular compartment.** The one compound of the four with a clean
textbook fit and no bespoke addition:

```
dxdt_GUT_ASC  = F_ASC * dose_rate_asc - KA_ASC * GUT_ASC;
dxdt_CENT_ASC = KA_ASC * GUT_ASC / V1_ASC - (CL_ASC / V1_ASC) * CENT_ASC;
```
`C_ASC = CENT_ASC` (plasma, PD-facing directly — no tissue distinction, per
the original).

**A note on `CENT_IMT`/`CENT_DAS`/`CENT_NIL`/`CENT_ASC` already being
concentration-valued, not amount-valued.** Unlike the `TOCI_C` quirk
documented in `sepsis/sep_refactor_notes.md` (a compartment dosed and read
directly as a concentration with *no* volume division anywhere), these four
compartments *do* divide by volume once, at the point mass enters from the
depot (`KA_* * GUT_* / V1_*`) and once again at elimination
(`CL_*/V1_* * CENT_*`) — this is algebraically identical to the standard
amount-based one-compartment model divided through by `V1_*`
(`d(amt)/dt = KA*GUT - CL/V1*amt`, then `C = amt/V1` substituted
throughout), so the compartment already *is* the concentration by
construction, consistently, not a shortcut that skips a volume term. No
renormalization was needed or performed; `CENT_<STEM>` carries the
original's exact numeric trajectory.

## Naming

| Role | IMT | DAS | NIL | ASC |
|---|---|---|---|---|
| Depot | `GUT_IMT` (was `Gut_imt`) | — (no depot; continuous rate) | `GUT_NIL` (was `Gut_nil`) | `GUT_ASC` (was `Gut_asc`) |
| Central (plasma) | `CENT_IMT` (was `Cp_imt`) | `CENT_DAS` (was `Cp_das`) | `CENT_NIL` (was `Cp_nil`) | `CENT_ASC` (was `Cp_asc`) |
| Intracellular (bespoke) | `IC_IMT` (was `Cic_imt`) | `IC_DAS` (was `Cic_das`) | `IC_NIL` (was `Cic_nil`) | — |
| Absorption/bioavail. | `KA_IMT`/`F_IMT` (was `ka_imt`/`F_imt`) | `KA_DAS` unused (was `ka_das`), `F_DAS` (was `F_das`) | `KA_NIL`/`F_NIL` (was `ka_nil`/`F_nil`) | `KA_ASC`/`F_ASC` (was `ka_asc`/`F_asc`) |
| Volume/clearance | `V1_IMT`/`CL_IMT` (was `Vd_imt`/`CL_imt`) | `V1_DAS`/`CL_DAS` (was `Vd_das`/`CL_das`) | `V1_NIL`/`CL_NIL` (was `Vd_nil`/`CL_nil`) | `V1_ASC`/`CL_ASC` (was `Vd_asc`/`CL_asc`) |
| Intracellular in/out rate | `KIN_IMT`/`KOUT_IMT` (was `Kin_imt`/`Kout_imt`) | `KIC_DAS`/`KOUT_DAS` (was `Kic_das`/`Kout_das`) | `KIC_NIL`/`KOUT_NIL` (was `Kic_nil`/`Kout_nil`) | — |
| Exposed concentration | `C_IMT` (= `IC_IMT`) | `C_DAS` (= `IC_DAS`) | `C_NIL` (= `IC_NIL`) | `C_ASC` (= `CENT_ASC`) |
| Hill params | `EMAX_IMT`/`EC50_IMT`/`GAMMA_IMT` (was `Emax_imt`/`IC50_imt`/`Hill_imt`) | `EMAX_DAS`/`EC50_DAS`/`GAMMA_DAS` (was `Emax_das`/`IC50_das`/`Hill_das`) | `EMAX_NIL`/`EC50_NIL`/`GAMMA_NIL` (was `Emax_nil`/`IC50_nil`/`Hill_nil`) | `EMAX_ASC`/`EC50_ASC`/`GAMMA_ASC` (was `Emax_asc`/`IC50_asc`/`Hill_asc`) |
| Effect | `EFFECT_IMT` (was `E_imt`) | `EFFECT_DAS` (was `E_das`) | `EFFECT_NIL` (was `E_nil`) | `EFFECT_ASC` (was `E_asc`) + `EFFECT_ASC_T315I` (was `E_asc_r`, params `EMAX_ASC_T315I`/`EC50_ASC_T315I`/`GAMMA_ASC_T315I`, was `Emax_asc_r`/`IC50_asc_r`/shared `Hill_asc`) |

All parameter *values* are copied verbatim from the original — nothing
invented, nothing defaulted, nothing dropped (including the dead `KA_DAS`).
`GAMMA_ASC_T315I` is a new independently-named parameter carrying the same
numeric value as `GAMMA_ASC` (`1.3`), since the original used one shared
`Hill_asc` for both of asciminib's effect terms; splitting it into two
named parameters (one per effect target, per the "never collapse several
drugs into one shared Hill term" spirit, applied here to one drug with two
targets) changes nothing numerically.

`dose_imt`/`dose_das`/`dose_nil`/`dose_asc` (the external mg/day switches
the R-side scenario code sets via `param()`) are **not** renamed — they are
covariate inputs, not part of the guide's PK-role naming table, and keeping
them identical means the original's own `scen1..scen7` blocks work
unmodified against the refactored model.

Two mechanical consequences of the compartment renames, both required for
the file to still mean the same thing (not scope creep):
- `init_tfr`'s named vector (Scenario 7, TFR attempt) was updated from
  `Gut_imt=0, Cp_imt=0, Cic_imt=0, Cp_das=0, Cic_das=0, Gut_nil=0, Cp_nil=0,
  Cic_nil=0, Gut_asc=0, Cp_asc=0` to the renamed compartments
  (`GUT_IMT=0, CENT_IMT=0, IC_IMT=0, CENT_DAS=0, IC_DAS=0, GUT_NIL=0,
  CENT_NIL=0, IC_NIL=0, GUT_ASC=0, CENT_ASC=0`), since `init()` matches by
  compartment name.
- Plot 3's `pk_data` column selections (`select(time_months, Cp_imt,
  Cic_imt)` etc.) were updated to the renamed columns
  (`CENT_IMT`/`IC_IMT`, `CENT_DAS`/`IC_DAS`, `CENT_NIL`/`IC_NIL`).

## The Hill interface: rename, not a fit

All five effect terms (`E_imt`, `E_das`, `E_nil`, `E_asc`, `E_asc_r`) were
already written as plain `Emax*pow(C,Hill)/(pow(IC50,Hill)+pow(C,Hill))` in
the original — the exact target shape, `Hill_*` already an explicit
parameter (not defaulted to 1). This is a pure rename in every case:
`EMAX_<STEM>`/`EC50_<STEM>`/`GAMMA_<STEM>` carry the original's
`Emax_*`/`IC50_*`/`Hill_*` values unchanged, no `nls()` fitting was needed
or performed. Where the original already combined each compound's own
effect only at the point of use (`E_total_sens = fmin(0.98, E_imt + E_das +
E_nil + E_asc)`, `E_total_res = fmin(0.98, E_asc_r)`), that combining logic
is untouched — it already matched the guide's "combine only at the point
disease equations actually use them" rule before this refactor; only the
five inputs were renamed to `EFFECT_IMT`/`EFFECT_DAS`/`EFFECT_NIL`/
`EFFECT_ASC`/`EFFECT_ASC_T315I`.

## qspserver compatibility

Following the precedent set by every prior refactor in this repo that hit
the same constraint (`sepsis`, `age-related-macular-degeneration`,
`acute-intermittent-porphyria`, `abdominal-aortic-aneurysm`, ...): declaring
`C_<STEM>`/`EFFECT_<STEM>` in `$PARAM` fails to build under mrgsolve 2.0.1
(`$PARAM` members are passed into `$ODE` as read-only references, and these
values are recomputed every timestep — the two requirements are mutually
exclusive in this build). All nine (`C_IMT`, `C_DAS`, `C_NIL`, `C_ASC`,
`EFFECT_IMT`, `EFFECT_DAS`, `EFFECT_NIL`, `EFFECT_ASC`,
`EFFECT_ASC_T315I`) are therefore plain `$ODE`-local `double`s, listed
(bare names, no collision so no `_out` suffix needed) in `$CAPTURE` — fully
discoverable via `POST /run_simulation`'s `outputs` selection and via
`POST /model_manifest`'s `outputPaths` (confirmed empirically, see
verification below), just not as an overridable `$PARAM`. The Hill
parameters themselves (`EMAX_*`, `EC50_*`, `GAMMA_*`) and every PK constant
**are** ordinary `$PARAM` values, fully overridable.

## Four pre-existing build defects found while verifying (not fixed
upstream, logged as `translations/UPSTREAM_ISSUES.md` #62)

Confirmed on the **untouched original** via `POST /model_manifest` — none
introduced by this refactor, all four already present in
`cml_mrgsolve_model.R`:

1. **`$INIT @annotated` uses `NAME = value : description`, not the required
   `NAME : value : description`.** Fails at the R-level annotation parser
   before any C++ compilation (`Error: improper annotation format`).
2. **Once (1) is worked around, `$CMT @annotated` and `$INIT @annotated`
   jointly redeclare all 22 compartments** (`Duplicated model names`).
3. **`$CAPTURE` repeats seven compartment names** already in `$CMT`
   (`Cp_imt`, `Cic_imt`, `Cp_das`, `Cic_das`, `Cp_nil`, `Cic_nil`,
   `Cp_asc`) — `compartment should not be in $CAPTURE`.
4. **`$TABLE` declares `double BCRPCT`/`double LOG_IS`/`double
   Resist_frac` as locals *and* has an explicit `capture NAME = NAME;` line
   for each, while the same bare names are also listed in `$CAPTURE`** —
   mrgsolve 2.0.1 auto-promotes the bare `$TABLE` local to a class member,
   which then collides with the `capture`-driven redeclaration
   (`redefinition of 'capture ...::BCRPCT'`). Every *other* `$TABLE`
   capture in this file already uses a distinct source name (e.g.
   `capture CHR = CHR_reached;`); only these three reuse the same name for
   both.

**Fix applied directly to the delivered `cml_mrgsolve_model_refactored.R`**
(per the guide's settled policy for this situation), all syntax-only and
non-numeric: (1) `$INIT`'s `=` changed to `:`; (2) the standalone `$CMT`
block was removed, its per-compartment descriptions folded into `$INIT`
(which already carries every value) — no compartment, order, or starting
value changed; (3) the seven duplicated names were dropped from
`$CAPTURE` (compartment states are always present in mrgsolve's output
regardless of `$CAPTURE` membership — confirmed by diffing
`/model_manifest`'s `outputPaths` before/after: every renamed compartment
is still listed); (4) the three colliding `$TABLE` locals were renamed
`BCRPCT_calc`/`LOG_IS_calc`/`Resist_frac_calc` (every internal reference to
them updated to match), while the `capture BCRPCT = BCRPCT_calc;`-style
lines and the external `$CAPTURE` names (`BCRPCT`, `LOG_IS`, `Resist_frac`)
are exactly as before — a pure local-variable rename, the reported values
are numerically identical to the original. `cml_mrgsolve_model.R` itself
was never touched and still carries all four defects exactly as written;
the identical four fixes were applied to an in-memory scratch copy of the
original purely so it could build for the `/run_simulation` comparison
below.

## A fifth finding: unbounded growth from a dead density-dependence term
(behavioral, not a build blocker — logged, not fixed, in either file)

The `$ODE` block computes `double suppress = fmax(0.0, 1.0 - Kcomp *
Ntotal);`, evidently meant to cap total marrow cellularity against
`K_total`, but it is only ever used as `suppress * 0.0` inside `dxdt_y0` —
multiplied by a literal zero, a complete no-op. Nothing in the model
actually enforces the stated carrying capacity. This is why, in three of
the four single-TKI monotherapy scenarios the original's own R script
defines (imatinib 400mg/day, dasatinib 100mg/day, nilotinib 300mg BID —
run in isolation, no combination), leukemic/normal cell counts grow without
bound and the shared simulation reaches `Inf` then `NaN` partway through
the original's own 10-year horizon. Confirmed identical in both models (see
verification below): this is a pre-existing dynamical defect, faithfully
reproduced by the refactor, not introduced by it.

## Verification

**Method**: qspserver `mrgsolve_api` container at `http://localhost:8007`
(`POST /model_manifest`, `POST /run_simulation`), requests spaced ~2s apart
per the API's stated concurrency limit. `model_content` was the bare DSL
extracted from each file's quoted `cml_code <- '...'` string (verbatim
text extraction, confirmed byte-identical to the in-file string).

**Scenarios tested**: exactly the four single-TKI dosing scenarios already
defined in the original file's own R code — Scenario 2 (imatinib
400mg/day), Scenario 3 (dasatinib 100mg/day), Scenario 4 (nilotinib 300mg
BID = `dose_nil=600`), Scenario 5 (asciminib 40mg BID = `dose_asc=80`) —
each run 0–87660h (the original's own `t_end_hr = 10yr`) at `delta=24`
(daily, matching `mrgsim(end = t_end_hr, delta = dt)` in the original),
3653 output rows per run, through both the (defect-patched, for
buildability only) original and the refactored model.

**Result: exact match up to the shared blow-up, for all four compounds.**

| Compound | Scenario | Shared outputs compared | Result |
|---|---|---|---|
| Imatinib | 400mg/day | `BCRPCT`, `LOG_IS`, `CHR`, `MMR`, `MR4`, `MR45`, `WBC`, `LSC_total`, `LSC_quies`, `Resist_frac`, `E_imt_out`/`EFFECT_IMT`, `Cp_imt`/`CENT_IMT`, `Cic_imt`/`IC_IMT` | Max abs diff = 0.0 across all 1711 rows before both models reach `Inf`/`NaN` at the identical row (t=41064h / ~4.68yr); `IC_IMT` vs `Cic_imt` also matched exactly (max abs diff = 0.0) across the full 1722 numeric rows available for that pair (the raw compartment, not the `>0`-guarded `C_IMT`, which reads `0` instead of `NaN` post-blowup by the guard's own design — noted, not a discrepancy) |
| Dasatinib | 100mg/day | same 13-output set (`E_das_out`/`EFFECT_DAS`, `Cp_das`/`CENT_DAS`, `Cic_das`/`IC_DAS`) | Max abs diff = 0.0 across all 1711 rows before the identical shared blow-up at the same row |
| Nilotinib | 300mg BID | same 13-output set (`E_nil_out`/`EFFECT_NIL`, `Cp_nil`/`CENT_NIL`, `Cic_nil`/`IC_NIL`) | Max abs diff = 0.0 across all 1711 rows before the identical shared blow-up at the same row |
| Asciminib | 40mg BID | same set (`E_asc_out`/`EFFECT_ASC`, `Cp_asc`/`CENT_ASC`) plus the new `EFFECT_ASC_T315I` (no original counterpart, since the original's own asciminib-only scenario never exercises `E_asc_r` as a captured output) | Max abs diff = 0.0 across **all 3653 rows, the full 10-year horizon** — this scenario alone does not trigger the shared blow-up within the tested window. `EFFECT_ASC_T315I` ranges 0 -> 0.819 (approaching `EMAX_ASC_T315I=0.82` as `C_ASC` saturates), behaving sanely |

This exceeds the guide's "near-exact match" expectation for a pure
structural reorganization with no Hill-fitting for all four compounds —
every comparable output is a bit-for-bit match, as expected for a rename
with no change to the underlying arithmetic. No scenario needed shortening
for solver step-count reasons; all four ran to completion (or to the
shared, pre-existing blow-up) at the original's own `delta=24`.

## Deliverables

- `cml_mrgsolve_model_refactored.R` (this refactor)
- This notes file
- `driver-patches/data/compound_perturbation_census.md` rows for
  Asciminib, Dasatinib, Imatinib (IMT), Nilotinib filled in
- `translations/UPSTREAM_ISSUES.md` #62 (four build defects + the dead
  `suppress` finding)

The original `cml_mrgsolve_model.R` is untouched.
