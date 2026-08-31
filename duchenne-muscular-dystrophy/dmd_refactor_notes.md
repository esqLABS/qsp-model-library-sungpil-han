# Refactor notes — `duchenne-muscular-dystrophy/dmd_mrgsolve_model.R`

## Scope: three census rows, real identities confirmed from the code itself

`driver-patches/data/compound_perturbation_census.md` carries three rows
for this file: **ASO**, **Corticosteroid (CS)**, **Gene therapy (DYS)**.
All three real identities are directly evidenced by the original file's
own scenario labels and doses, matching known real-world drug profiles
exactly — no external lookup needed:

| Row | Real drug | Evidence |
|---|---|---|
| ASO | **Eteplirsen** (EXONDYS 51, exon-51-skipping phosphorodiamidate morpholino oligomer) | Scenario 4 is literally named `"4. Eteplirsen 30 mg/kg/wk"`, IV, weekly — the real approved EXONDYS 51 dosing regimen. |
| Corticosteroid (CS) | **Deflazacort** (base PK, `F_CS = 0.72`) and **Prednisone** (same shared PK block, `F_CS` override `0.82` in its own scenario) | Scenario 2 `"2. Deflazacort 0.9 mg/kg/d"`, scenario 3 `"3. Prednisone 0.75 mg/kg/d"` — both real CINRG-era DMD corticosteroid regimens, sharing one PK backbone as the original wrote it (differentiated only by bioavailability). |
| Gene therapy (DYS) | **Delandistrogene moxeparvovec** (Elevidys / SRP-9001, Sarepta) | Scenario 5 `"5. Gene Therapy (Elevidys)"`, single IV dose `1.33e14 vg/kg` — matches the real, FDA-approved Elevidys dose exactly. |

## Stem-collision correction: "DYS" cannot be gene therapy's stem

The census pre-filled the gene-therapy row's compound name as `"Gene
therapy (DYS)"`. Checking the actual code: `$CMT DYS` already exists as a
**disease-side** compartment ("Dystrophin level, % of normal" — see the
`$CMT @annotated` block, "Disease pathophysiology compartments" section),
tracked by its own `dxdt_DYS` ODE and read throughout the disease model
(`DYS_TOTAL = DYS_BASE + DYS + ...`). Using `DYS` as the gene-therapy PK
stem would collide with this pre-existing, untouched compartment. Used
**`AAV`** instead — this matches the original's own PK compartment prefix
(`AAV_CIRC`, `AAV_MUS`) and avoids the collision entirely. This is a
naming-convention correction, not a compound-identity correction: the
compound itself (delandistrogene moxeparvovec/Elevidys) is exactly what
the census row was pointing at.

## Census classification checked against the actual code, corrected where wrong

Per the fork's standing practice (`achondroplasia/acdp_refactor_notes.md`,
`clostridioides-difficile-infection/cdi_refactor_notes.md`), the census's
automated classifier labels are a first pass, not ground truth:

- **ASO: "Redirect concentration (clean single site)" — accurate, with a
  caveat.** The original's own PD term (`DYS_from_ASO = EFF_ASO *
  IC_ASO`) already reads a single, clean concentration site — but it is
  **not** the plasma-level `C_ASO` the original also computes in `$MAIN`
  (`C_ASO = CENT_ASO / (V1_ASO/1000.0 * WT)`); that plasma value is
  computed and then never read anywhere else in the model, not even
  captured. The real PD-driving site is the bespoke third PK compartment,
  `IC_ASO` (intracellular active pool). "Redirect" is applied literally:
  the exposed `C_ASO` in the refactor is redirected to alias `IC_ASO` (the
  site PD actually reads), not the dead plasma calculation. See "ASO
  archetype" below.
- **Corticosteroid (CS): "Normalize duplicate concentration sites, then
  redirect" — half right, correction below.** The original does compute a
  concentration formula twice (`C_CS` in `$MAIN`, `C_CS_ngmL` in `$TABLE`,
  same formula `CENT_CS/V1_CS*1000.0`), plus a third, entirely unused
  `C_CS_PERIPH`. But collapsing the first two into one site is **not**
  safe: `$MAIN` and `$TABLE` are separate mrgsolve execution scopes, and
  (confirmed empirically, see "The duplicate is architecturally required"
  below) a `$MAIN`-computed local and a `$TABLE`-computed local reading the
  *same* compartment produce genuinely different captured values at every
  reported row. The correct "normalize" step here is: drop the
  genuinely-unused `C_CS_PERIPH` (pure dead code), keep `C_CS` in `$MAIN`
  (feeds PD, unchanged), keep `C_CS_ngmL` in `$TABLE` (reporting-only,
  unchanged) — not collapse them.
- **Gene therapy (AAV): "Delete PK compartment; concentration is itself
  the state" — half right, correction below**, same pattern as the
  `clostridioides-difficile-infection`/Bezlotoxumab precedent. The
  "concentration is itself the state" half is correct: `MUS_AAV` (was
  `AAV_MUS`) needs no volume division — `C_AAV = MUS_AAV` is a plain
  identity, exactly as the original's own `DYS_from_AAV = EFF_AAV *
  AAV_MUS * 1e-12 * 1e3` already used it directly. But "Delete PK
  compartment" is wrong: there are **two** genuine, dosed, first-order PK
  compartments (`AAV_CIRC`→`AAV_MUS`, renamed `CENT_AAV`→`MUS_AAV`), with
  real one-way transduction kinetics (`KAAV`, `DEGRAD`) — nothing was
  deleted or collapsed. See "AAV archetype" below.

## Archetypes

### Corticosteroid (CS): Archetype 3 (depot + central + peripheral, linear)

Already an exact structural match to the guide's Archetype 3 template —
only the compartment names needed renaming (`DEPOT_CS`→`GUT_CS`,
`PERIPH_CS`→`PERI_CS`; `CENT_CS` already matched). All six PK parameter
names (`KA_CS`, `F_CS`, `CL_CS`, `V1_CS`, `V2_CS`, `Q_CS`) already matched
the convention exactly — no renaming needed, values unchanged.

`C_CS = CENT_CS / V1_CS * 1000.0` (ng/mL), computed in `$MAIN`, matching
the original's own block placement exactly (see "Keep a calculation in
the block the original used it in" below for why this matters).

**Two named effects, not one** — the original's corticosteroid term
touches the disease model in two structurally distinct ways, matching the
precedent already used for two-mechanism single compounds (e.g.
`drug-induced-liver-injury`'s APAP → `EFFECT_APAP` +
`EFFECT_APAP_BSEP`):

1. **`EFFECT_CS`** — NF-κB inhibition. The original's `CS_NF_INH =
   C_CS / (EC50_CS_NF + C_CS)` is already an exact Emax-model ratio.
   Renamed to the Hill interface with `EMAX_CS_NF = 1.0` and
   `GAMMA_CS_NF = 1.0`, both **new, math-implied parameters, not fit** —
   the original ratio already saturates at 1 as `C_CS → ∞`, and had no
   explicit Hill exponent (`γ = 1` implied). Confirmed exact rename:
   `pow(x, 1.0) == x` bit-for-bit (verified via qspserver on both a static
   probe and full live integration).
2. **The membrane-stabilization term** — linear in `C_CS`, **not**
   Hill-shaped, and **kept fully inline in `$ODE`, not extracted into a
   named `EFFECT_CS_MEM` double.** See "A reassociation bug caught before
   delivery" below for why — this is the one place in this refactor where
   the guide's naming convention (one named function per compound effect)
   was deliberately *not* applied, and the reason is empirical, not
   stylistic.

Also disclosed: the original's own four-term membrane expression —
`- KD_CS_MEM*C_CS*MEMI*(-1.0) + KD_CS_MEM*C_CS*(1.0-MEMI)` — algebraically
reduces to `+ KD_CS_MEM*C_CS` independent of `MEMI` (the two `MEMI` terms
cancel). Not simplified in the refactor (would change nothing numerically
here, but *would* have changed the floating-point operation order — see
below), just noted as a curiosity in the original's own arithmetic.

**Dead code removed:** `C_CS_PERIPH = PERIPH_CS/V2_CS*1000.0`, computed in
the original's `$MAIN` and never referenced anywhere else in the model
(not in any `dxdt_`, not in `$TABLE`, not captured). Confirmed its removal
changes nothing (see Verification).

### ASO (Eteplirsen): bespoke 3-site (central/peripheral linear PK + intracellular active pool)

`CENT_ASO` already matched the convention. `MUS_ASO`→`PERI_ASO` (tissue
distribution site, matches `PERI_<STEM>` semantics). `IC_ASO` is kept
under its original name — a genuine, mechanistically real third site
(intracellular active pool, distinct pharmacokinetics: slow uptake
`KUP_ASO = 0.0015/h`, slow elimination `KOUT_ASO = 0.0008/h`) that doesn't
map to `GUT`/`CENT`/`PERI`; disclosed as bespoke. This is real,
well-supported PMO/eteplirsen pharmacology (renal-cleared from plasma,
slow muscle uptake, slower intracellular retention) — not a convenience
shortcut. The guide's naming convention explicitly allows "two [exposed
concentrations] only when a genuinely different tissue site matters" —
here that is exactly the case: `C_ASO` (informational plasma level, dead
in the original, dropped — see below) and the real PD-driving
intracellular site.

`C_ASO := IC_ASO` (plain identity, `$MAIN`, matching the original's own
block placement of `DYS_from_ASO = EFF_ASO * IC_ASO`).
`EFFECT_ASO = EFF_ASO * C_ASO` — an exact rename of the original's own
linear (non-Hill) term; **not** force-fit into an Emax/Hill shape the
original never had, per the guide's explicit instruction to leave a
native, non-saturating calculation in place rather than force a bad fit.

**Dead code removed:** the original's own `C_ASO = CENT_ASO /
(V1_ASO/1000.0 * WT)` (plasma concentration) is computed in `$MAIN` and
never read anywhere else in the model — not by any `dxdt_`, not by
`$TABLE`, not even captured. Confirmed its removal changes nothing (see
Verification). The name `C_ASO` is reused in the refactor for the
*actually-PD-driving* quantity (`IC_ASO`), per the redirect classification
above — this is a deliberate redirect, not an accidental collision.

### Gene therapy (AAV, delandistrogene moxeparvovec/Elevidys): bespoke 2-site, irreversible one-way transduction

`AAV_CIRC`→`CENT_AAV` (circulating vector, dosed directly — the file's
own scenario 5 is a single IV bolus into this compartment, matching real
Elevidys administration). `AAV_MUS`→`MUS_AAV` (muscle-transduced vector,
the actual PD-driving site) — kept as a bespoke name rather than
`PERI_AAV`: the transfer `CENT_AAV → MUS_AAV` is **irreversible** (vector
never returns from muscle to circulation), unlike the guide's Archetype 2
`Q_<STEM>` pattern, which is bidirectional by construction. Renaming it
`PERI_AAV` would misleadingly imply a reversible peripheral-distribution
compartment the original does not have.

Parameters renamed to the stem: `KAAV`→`KTRANSD_AAV` (transduction rate,
0.003/h), `DEGRAD`→`CL_AAV` (circulating vector degradation/clearance,
5e-5/h). `KDEC_DYS` (muscle-vector decay, 0.00008/h, "very slow") is
**split**: the disease-side `dxdt_DYS` still uses the original
`KDEC_DYS` name/value untouched (out of this compound's scope — see
below), while the AAV compartment's own decay uses a new,
identically-valued `KDEC_AAV` — a pure naming split (same number, two
names) so the AAV PK block does not share a parameter object with
untouched disease-side code, per the guide's "own clearly-delimited
block" requirement.

`C_AAV = MUS_AAV` — a plain identity (no volume division), matching the
original's own direct use of `AAV_MUS` inside `DYS_from_AAV`. This is the
one respect in which the census's "concentration is itself the state"
classification was correct.
`EFFECT_AAV = EFF_AAV * C_AAV * 1e-12 * 1e3` — exact rename of the
original's own `DYS_from_AAV`, same "simplified scaling" constant, linear
(non-Hill), not force-fit.

## Pre-existing defect, disease-side, left untouched: `KEXP_DYS` is dead code

`dxdt_DYS = KEXP_DYS * (DYS_TOTAL > 0 ? 0.0 : 0.0) - KDEC_DYS * DYS + 0.0;`
— both branches of the ternary are the literal `0.0`, so `KEXP_DYS`
("Micro-dystrophin expression rate") never contributes anything under any
value of `DYS_TOTAL`. The `DYS` compartment only ever decays from its
initial `0.1`; despite the adjacent comment claiming it is "replenished by
ASO effect ... and AAV," the actual drug-driven dystrophin restoration is
entirely algebraic in `$MAIN` (`DYS_TOTAL = DYS_BASE + DYS + DYS_from_ASO
+ DYS_from_AAV`), never fed back into `DYS`'s own ODE. `KEXP_DYS` is
declared inside the gene-therapy `$PARAM` section by original placement,
but its only (non-functional) use site is the disease-side `dxdt_DYS`
line — it is not functionally part of the AAV compound's own PK block at
all, so it was left completely untouched (name, value, and the dead
ternary) rather than folded into the AAV renaming. Logged as
`UPSTREAM_ISSUES.md` **#137**, not fixed, per "log what you find, don't
fix it" — fixing this would be a real behavioural change to the disease
model's own dynamics, out of scope for a structural PK/PD refactor.

## Pre-existing defect, disease-side, left untouched: `Dystrophin_pct` undercounts the gene-therapy scenario

`$TABLE`'s `Dystrophin_pct = DYS_BASE + DYS + EFF_ASO * IC_ASO` omits the
AAV term that `$MAIN`'s `DYS_TOTAL` (which drives every actual PD dynamic
— `DYS_EFF`, and through it `MEMI`, `SC`, `MF`, `SWD`, `FVC_pct`,
`LVEF_pct`) includes. Scenario 5 ("Gene Therapy (Elevidys)") therefore
correctly simulates AAV's disease-modifying effect internally, but the
file's own plotted output (`p2`, "Dystrophin Level Over Time," built from
`Dystrophin_pct`) silently undercounts it for that scenario specifically.
Confirmed this is not a `$MAIN`/`$TABLE` visibility artifact of this
refactor (see below) — the original *chose* to recompute
`Dystrophin_pct` from scratch in `$TABLE`, and that recomputation itself
omits the AAV term; a correct recomputation was always available (`+
EFF_AAV * AAV_MUS * 1e-12 * 1e3`) and simply wasn't written. Logged as
`UPSTREAM_ISSUES.md` **#137** (same entry as the `KEXP_DYS` finding — both
are dystrophin-accounting inconsistencies in the original's own code).
Preserved byte-for-byte identical in the refactor (confirmed by exact
match on `Dystrophin_pct` in the gene-therapy scenario, see Verification).

## Keep a calculation in the block the original used it in

Empirically confirmed via qspserver (minimal reproduction, then
reconfirmed against this file's own scenarios) exactly what the guide's
corresponding section warns about, with concrete numbers:

- A `double` declared in `$MAIN` and then listed in `$CAPTURE` **does**
  compile and **is** discoverable via `/model_manifest`'s `outputPaths` —
  no need to move it into `$ODE`/`$TABLE` for capturability.
- But its *reported/captured value at each output row is one full
  reporting interval behind the live state* — because `$MAIN` runs once
  per interval, using state from the *start* of that interval, and that
  value is then held fixed through all of `$ODE`'s internal solver
  substeps for the interval (this is the standard, deliberate mrgsolve
  covariate-update pattern, not a bug). A minimal test model confirmed
  this precisely: with a dose at `t=0` into a one-compartment model, a
  `$MAIN`-computed `C_TEST = CENT/V1` reports `[0, 0, 20, 18, 16.0, ...]`
  while a `$TABLE`-computed `OTHER = CENT/V1` (same formula, same state)
  reports `[0, 20, 18, 16.0, 14.4, ...]` — `C_TEST` is systematically one
  row behind `OTHER`, at every single reported row, not just at the dose
  instant.
- This directly explains — and justifies keeping — the original's own
  "duplicate" `C_CS` (`$MAIN`, PD-driving, one-interval-lagged) /
  `C_CS_ngmL` (`$TABLE`, reporting-only, live) pair: they are **not**
  numerically interchangeable, and the census's "normalize duplicate
  concentration sites" instruction, applied naively, would have either
  changed the corticosteroid PD trajectory (if `dxdt_MEMI`/`dxdt_NFkB`
  were switched to the live value) or changed the reported concentration
  column's own values (if `$TABLE` were switched to reference `$MAIN`'s
  lagged local) — both would be real behavioural changes disguised as
  cleanup. Kept as two independent, identically-formulaed computations,
  one in each original block, exactly as the original had it.
- `C_ASO`, `C_AAV`, `EFFECT_ASO`, `EFFECT_AAV`, `EFFECT_CS`, `DYS_TOTAL`,
  `DYS_EFF` all stayed in `$MAIN`, matching every one of their original
  computation sites (`DYS_from_ASO`, `DYS_from_AAV`, `DYS_TOTAL`,
  `DYS_EFF`, `CS_NF_INH` were all `$MAIN`-declared in the original).

## A reassociation bug caught before delivery

An early draft of this refactor extracted the corticosteroid
membrane-stabilization arithmetic into its own named `EFFECT_CS_MEM`
double:

```
double EFFECT_CS_MEM = - KD_CS_MEM * C_CS * MEMI * (-1.0)
                        + KD_CS_MEM * C_CS * (1.0 - MEMI);
dxdt_MEMI = KR_MEM * DYS_EFF * (1.0 - MEMI) * SC_EFF
            - KD_MEM * (1.0 - DYS_EFF) * CAI
            + EFFECT_CS_MEM;
```

This is algebraically identical to the original's own single four-term
chain, and a static, single-point probe (fixed inputs, no ODE
integration) confirmed `diff == 0.0` between the two forms. But run
through a real, live simulation (deflazacort scenario, 240h, qspserver
`/run_simulation`), this draft diverged from the patched-original baseline
by **up to 5.4×10¹¹** on `MEMI` and NF-κB pathway outputs by hour 240 — a
massive, clearly-a-bug-not-noise mismatch.

Root cause, isolated by a controlled A/B test (four variants: original
`EFFECT_CS` formula × original inline `dxdt_MEMI` vs. the Hill-rewritten
`EFFECT_CS` × the `EFFECT_CS_MEM`-extracted `dxdt_MEMI`, tested pairwise):
the `EFFECT_CS` Hill rewrite (`pow(C_CS, 1.0)` etc.) was **not** the
problem — reverting only that still showed the full divergence. Reverting
only the `EFFECT_CS_MEM` extraction (keeping the original inline
`- term*(-1.0) + term` chain, in the original's exact operation order)
eliminated the divergence completely (`max abs diff = 0.0`). The
extraction reorders the underlying IEEE754 floating-point operations
(`(A - B) - C + D` vs. `(A - B) + (C + D)`, not guaranteed bit-identical
under non-associative floating-point arithmetic) — normally an
imperceptible, sub-ULP effect, but this file's disease model has a
genuine, pre-existing runaway positive-feedback instability (see
`UPSTREAM_ISSUES.md` **#136**: `ROS → NFkB → M1 → ROS`, unbounded, no
ceiling/floor anywhere in the original), which amplifies that sub-ULP
seed by roughly eleven orders of magnitude within one 240-hour dosed
scenario.

**Fixed** by keeping the CS membrane term fully inline in `$ODE`, in the
original's exact four-term operation order, with only identifier renames
(`C_CS`, `MEMI`, `KD_CS_MEM` are all already the same names as the
original — the only names that changed anywhere nearby are the
compartment renames `DEPOT_CS`→`GUT_CS`/`PERIPH_CS`→`PERI_CS`, which do
not appear in this expression at all). No named `EFFECT_CS_MEM` variable
exists in the delivered file. This is disclosed here specifically because
it is the one place this refactor's own naming convention (one named
function per compound effect) was **not** applied to a real, distinct
compound-effect term — `EFFECT_CS` already satisfies the convention for
corticosteroid's primary (NF-κB) effect, and forcing a second explicit
name onto the membrane term turned out to be unsafe in this particular
file, not merely a stylistic choice foregone.

## Numerical fragility, pre-existing, disclosed, not fixed

Same class as `essential-thrombocythemia` (#64), `von-willebrand-disease`
(#72), `takayasu-arteritis` (#74), `achondroplasia` (#125) — logged as
`UPSTREAM_ISSUES.md` **#136**. The untreated natural-history scenario's
own `ROS → NFkB → M1 → ROS` positive-feedback loop has no ceiling or floor
anywhere in the original's equations. Run far enough (`t > 276.035h`),
`lsoda`'s adaptive step size collapses to its floor trying to resolve an
ever-steeper trajectory, and the solve fails outright (`negative istate:
-4`). Confirmed **identical** — same trigger time, `t = 276.035`, to the
reported precision — in three independently-tested DSL variants: the
patched-original baseline, the pre-fix refactor draft (before the
reassociation fix above), and the final delivered refactor. This is not
introduced or altered by the refactor in any way.

## Naming reference

| Original | Refactored | Value | Role |
|---|---|---|---|
| `DEPOT_CS` (cmt) | `GUT_CS` | — | CS absorption depot |
| `CENT_CS` (cmt) | `CENT_CS` | — | CS central (unchanged) |
| `PERIPH_CS` (cmt) | `PERI_CS` | — | CS peripheral |
| `KA_CS`,`F_CS`,`CL_CS`,`V1_CS`,`V2_CS`,`Q_CS` | unchanged | 1.2, 0.72, 25.0, 35.0, 80.0, 8.0 | already matched convention |
| — | `EMAX_CS_NF` (new) | 1.0 | Hill ceiling [math-implied] |
| `EC50_CS_NF` | unchanged | 15.0 | Hill EC50 |
| — | `GAMMA_CS_NF` (new) | 1.0 | Hill coefficient [explicit, original had none] |
| `CS_NF_INH` | `EFFECT_CS` | — | NF-κB inhibition (exact rename) |
| `CENT_ASO` (cmt) | `CENT_ASO` | — | ASO plasma (unchanged) |
| `MUS_ASO` (cmt) | `PERI_ASO` | — | ASO tissue |
| `IC_ASO` (cmt) | `IC_ASO` | — | ASO intracellular active (unchanged, bespoke) |
| `CL_ASO`,`V1_ASO`,`V2_ASO`,`Q_ASO`,`KUP_ASO`,`KOUT_ASO` | unchanged | 210, 290, 6500, 18, 0.0015, 0.0008 | already matched convention |
| — (dead, plasma-only, unused) | dropped | — | original's own unused `C_ASO` local |
| `DYS_from_ASO` | `EFFECT_ASO` | — | exact rename, linear |
| `AAV_CIRC` (cmt) | `CENT_AAV` | — | AAV circulating |
| `AAV_MUS` (cmt) | `MUS_AAV` | — | AAV muscle-transduced (bespoke, irreversible sink) |
| `KAAV` | `KTRANSD_AAV` | 0.003 | transduction rate |
| `DEGRAD` | `CL_AAV` | 5e-5 | vector degradation |
| `KDEC_DYS` (AAV_MUS use only) | `KDEC_AAV` (new) | 0.00008 | split, same value; disease-side `KDEC_DYS` untouched |
| `KEXP_DYS` | unchanged | 0.004 | dead/unused, disclosed, not renamed (see above) |
| `DYS_from_AAV` | `EFFECT_AAV` | — | exact rename, linear |
| — (dead) | dropped | — | original's own unused `C_CS_PERIPH` local |
| `M1_0` | `M1_BASELINE` | 15.0 | build-compat rename (used in `NF_INPUT`) |
| `M2_0`,`TGFb_0`,`FIB_0`,`SC_0`,`MF_0`,`SWD_0` | `*_BASELINE` | unchanged | build-compat rename (all dead/unused) |
| `FVC_0`,`LVEF_0` | unchanged | 95.0, 62.0 | no collision (compartments are `FVC_pct`/`LVEF_pct`), dead/unused, left as-is |

All disease-side compartments (`DYS MEMI CAI ROS NFkB M1 M2 TGFb FIB SC MF
SWD FVC_pct LVEF_pct`) and every disease-side parameter not listed above
are completely untouched — not renamed, not restructured.

## Build-compatibility fixes (the original does not compile at all)

Logged as `UPSTREAM_ISSUES.md` **#135**. `POST /model_manifest` on the
untouched original's own extracted `code_dmd <- '...'` string returns
HTTP 500 for three independent, pre-existing reasons:

1. `$CMT @annotated` and a separate `$INIT` block jointly redeclare all 22
   compartments (`Duplicated model names`). Fixed: `$INIT` dropped, the
   same 14 non-zero initial values set via `if (NEWIND <= 1) { <cmt>_0 =
   value; }` inside `$MAIN`.
2. `$CAPTURE` re-lists 14 of those same compartment names directly
   (`compartment should not be in $CAPTURE`). Fixed: dropped — they remain
   reported automatically.
3. Once (1) and (2) are worked around, seven baseline `$PARAM` names
   (`M1_0`, `M2_0`, `TGFb_0`, `FIB_0`, `SC_0`, `MF_0`, `SWD_0`) collide
   with mrgsolve's own auto-generated `<compartment>_0` initial-value
   symbol for the identically-named compartment (a real C++ redeclaration
   error, confirmed via a minimal reproduction before touching the full
   file). Fixed: renamed `<CMT>_BASELINE`, values unchanged — six of the
   seven are otherwise dead/unused, the seventh (`M1_0`) is used in
   `NF_INPUT` and its one use site was updated to `M1_BASELINE` too.

All three fixes are syntax-only and non-numeric. Applied directly to the
delivered `dmd_mrgsolve_model_refactored.R` (not just a scratch copy), per
the guide's settled policy. The same three fixes were applied to an
in-memory-only **patched-original** scratch copy (original DSL, otherwise
byte-identical) used solely as the verification baseline below — never
written into the repo tree, not part of the deliverables, deleted after
use.

## Verification

**Method.** Extracted the bare DSL text from both `dmd_mrgsolve_model.R`
(`code_dmd <- '...'`, unmodified) and `dmd_mrgsolve_model_refactored.R`
(`code_dmd_refactored <- '...'`, refactored) and ran them through the
local qspserver `mrgsolve_api` (`http://localhost:8007`), `POST
/model_manifest` then `POST /run_simulation`, requests spaced ~2.2s apart
per the shared-service note. Since the untouched original does not
compile at all (#135), the comparison baseline is the patched-original
scratch copy described above.

**Window.** The original's own scenarios run 6 years (`end = 52560,
delta = 12`). This file's disease model has a genuine, pre-existing
positive-feedback numerical instability (#136) that makes both the
original and the refactored model fail to solve past simulated hour
`276.035` on the untreated/corticosteroid arms. Per the guide's documented
solver-budget/stability allowance, verification uses a shortened window,
`end = 240, delta = 12` (matching the original's own `delta`), well inside
the stable region for every scenario tested — 10 or 12 doses per regimen,
enough to exercise absorption, distribution, and multi-dose accumulation
for every compound.

**Scenarios run — the file's own six, not invented, dosing amounts/timing
unchanged, only the `cmt` target renamed:**

1. **Natural History** (no dosing): exact match, max abs diff `0.0`, on
   every one of 26 shared outputs (8 renamed PK compartments/sites, 14
   disease compartments, `C_CS_ngmL`, `CK_serum`, `Dystrophin_pct`,
   `NSAA`), 21 rows.
2. **Deflazacort 0.9 mg/kg/d** (daily × 10): exact match, max abs diff
   `0.0`, all 26 outputs, 22 rows. (This scenario is what first exposed
   the reassociation bug above, before the fix — see that section for the
   pre-fix numbers.)
3. **Prednisone 0.75 mg/kg/d** (daily × 10, `F_CS = 0.82` override): exact
   match, max abs diff `0.0`, all 26 outputs, 22 rows — confirms the
   parameter-override path (shared PK block, different `F_CS`) is
   unaffected.
4. **Eteplirsen 30 mg/kg/wk** (weekly × 2 IV): exact match, max abs diff
   `0.0`, all 26 outputs, 22 rows.
5. **Gene Therapy (Elevidys)** (single IV bolus, `1.33e14 vg/kg`): exact
   match, max abs diff `0.0`, all 26 outputs including `Dystrophin_pct`
   (confirming the AAV-omission defect above reproduces identically, not
   just structurally), 22 rows.
6. **Deflazacort + Eteplirsen** (combination): exact match, max abs diff
   `0.0`, all 26 outputs, 23 rows — confirms the two independently-dosed
   compounds combine correctly with no cross-contamination.

This is the expected, and achieved, outcome for a pure structural
reorganization (Archetype 3 for CS, bespoke non-Hill archetypes for ASO
and AAV, no `nls()` fitting anywhere — every effect term in this file was
already either an exact Emax/Hill shape or a plain linear proportional
term): genuinely `0.0` at every point, over the full stable window, for
all six of the file's own scenarios.

`/model_manifest` on the refactored DSL additionally confirmed:
- **`parameters`** includes every renamed PK parameter (`KTRANSD_AAV`,
  `CL_AAV`, `KDEC_AAV`, `EMAX_CS_NF`, `GAMMA_CS_NF`, `M1_BASELINE`,
  `M2_BASELINE`, `TGFb_BASELINE`, `FIB_BASELINE`, `SC_BASELINE`,
  `MF_BASELINE`, `SWD_BASELINE`) with unchanged numeric defaults.
- **`outputPaths`** includes every `C_<STEM>`/`EFFECT_<STEM>` this
  refactor's naming convention creates (`C_CS`, `EFFECT_CS`, `C_ASO`,
  `EFFECT_ASO`, `C_AAV`, `EFFECT_AAV`), alongside every renamed
  compartment (`GUT_CS`, `PERI_CS`, `PERI_ASO`, `CENT_AAV`, `MUS_AAV`) and
  every original `$TABLE` output (`C_CS_ngmL`, `CK_serum`,
  `Dystrophin_pct`, `NSAA`).

No scratch/debug artifacts were left in the repository — DSL extraction,
the patched-original scratch copy, the A/B isolation tests, and the
comparison scripts ran entirely from a session-local scratchpad directory
outside the repo tree, deleted after use.

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`: all
three rows (ASO, Corticosteroid, Gene therapy) filled in with real drug
identities, target/pathway, renamed compartments, exposed
concentration/effect variable names, and this file's notes — including
the two classification corrections above (CS's normalize instruction was
half-safe, AAV's "delete PK compartment" was wrong about the compartment
count though right about the no-volume-division point).
