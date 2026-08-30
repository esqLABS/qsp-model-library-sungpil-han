# Refactor notes — `acn_mrgsolve_model.R` (ethinylestradiol, isotretinoin, 4-oxo-isotretinoin, spironolactone, tetracycline class)

**Scope of this pass.** Per the fork's PK/PD refactor spec
([`FORK_WORKFLOW_GUIDE.md`](../FORK_WORKFLOW_GUIDE.md), Part 2), 5
compounds were rewritten, per their existing rows in
[`driver-patches/data/compound_perturbation_census.md`](../driver-patches/data/compound_perturbation_census.md)
(all classified "Redirect concentration (clean single site)"): **EE**
(ethinylestradiol, stem `EE`), **ISO** (isotretinoin parent, stem `ISO`),
**OXO** (4-oxo-isotretinoin, the active metabolite of ISO, stem `OXO`),
**SPI** (spironolactone/canrenone, stem `SPI`), and **TET** (tetracycline
class — doxycycline default, minocycline/sarecycline via presets, stem
`TET`). Every other compound this 55-compartment file models — benzoyl
peroxide, topical retinoid, clindamycin, azelaic acid, dapsone,
clascoterone — is completely untouched: different compartments, different
`$PARAM` names, never referenced by any of the 5 renamed blocks except at
the pre-existing shared use sites noted below (which are read-only reads
of a renamed variable, not edits to the other compound's own logic).

## Is OXO really a 6th "compound", or is it ISO's metabolite?

Confirmed: OXO is 4-oxo-isotretinoin, ISO's own active metabolite (no
separate dose, no `GUT_OXO`; it is produced from `CENT_ISO` by first-order
metabolism, `KMET_ISO`). The census classifier found it because the
original file gives it its own compartment (`OXOC`) and its own PK
parameters (`CLOXO`, `VOXO`) — a real, independently-redirectable
concentration site — even though pharmacologically it is downstream of
ISO. Both are handled below, each on their own terms.

## Archetypes determined

### EE — archetype 3 minus the peripheral compartment (depot + central, linear elimination)

```
dxdt_GUT_EE  = -KA_EE * GUT_EE;
dxdt_CENT_EE =  KA_EE * GUT_EE * F_EE - (CL_EE / V1_EE) * CENT_EE;
```

Exactly the guide's Archetype-3 shape with the peripheral compartment
dropped. Two named Hill effects, matching two genuinely independent
pharmacological actions already present in the original as plain ratios:

| Original | Refactored | Value |
|---|---|---|
| `EED`/`EEC` (compartments) | `GUT_EE`/`CENT_EE` | — |
| `KAEE`/`FEE`/`VEE`/`CLEE` | `KA_EE`/`F_EE`/`V1_EE`/`CL_EE` | 1.20 (1/h) / 0.43 / 260 (L) / 12.0 (L/h) |
| `CEE` (local `$ODE` ratio) | `C_EE` | mg/L |
| `EMXEESH`/`EC5EESH` | `EMAX_EE_SHBG`/`EC50_EE_SHBG` | 1.70 / 2.5e-5 (mg/L) |
| — (none; implicit linear ratio) | `GAMMA_EE_SHBG` (new, `=1.0`) | rename, not a fit |
| `EMXLH`/`EC5LH` | `EMAX_EE_LH`/`EC50_EE_LH` | 0.42 / 1.8e-5 (mg/L) |
| — (none; implicit linear ratio) | `GAMMA_EE_LH` (new, `=1.0`) | rename, not a fit |
| `EESHBG` (local `$ODE` ratio) | `EFFECT_EE_SHBG` | hepatic SHBG-induction effect, used at `SHBGIN`'s one use site |
| `EELH` (local `$ODE` ratio) | `EFFECT_EE_LH` | LH-suppression effect, used at `TTIN`'s one use site |

`CEEO` (a pg/mL-scaled display duplicate the original computed
table-side, `POS(EEC)/VEE*1e6`) is kept unchanged in name and role,
formula updated to read `CENT_EE`/`V1_EE`. See the "$TABLE quirk" section
below for why this duplicate has to stay.

## OXO and ISO — a genuinely shared Hill, disclosed as bespoke

### ISO — archetype 3 (depot + central + peripheral), unchanged shape

```
dxdt_GUT_ISO  = -KA_ISO * GUT_ISO;
dxdt_CENT_ISO =  KA_ISO * GUT_ISO * F_ISO - (CL_ISO/V1_ISO)*CENT_ISO
                 - KMET_ISO*CENT_ISO - Q_ISO*(CENT_ISO/V1_ISO - PERI_ISO/V2_ISO);
dxdt_PERI_ISO =  Q_ISO * (CENT_ISO/V1_ISO - PERI_ISO/V2_ISO);
```

`F_ISO` (bioavailability) is food-dependent in the original
(`FISOF*(1+FOODEF*(LIDOSE>0.5?1:FOOD))`) — kept as a per-step `$ODE`
local (`F_ISO`), NOT a fixed `$PARAM`, exactly reproducing the original's
own mechanism; its 4 inputs are renamed `F0_ISO`/`FOODEF_ISO`/`FOOD_ISO`/
`LIDOSE_ISO`.

### OXO — bespoke: archetype-1-like single compartment, fed by metabolism not dosing

```
dxdt_CENT_OXO = KMET_ISO * CENT_ISO - (CL_OXO / V1_OXO) * CENT_OXO;
```

No `GUT_OXO` — this doesn't fit Archetype 1 as a *dosed* single
compartment; it's a metabolite pool with an inflow term instead of a
depot. Still a clean, independently-redirectable `C_OXO = CENT_OXO/V1_OXO`.

### Why one shared `EFFECT_ISO`, not two independent Hills

The original computes ISO's disease effect from a **potency-weighted sum
of parent + metabolite concentration**, then applies ONE Hill to that sum:

```
ISOEQ  = CISO + POTOXO * COXO;
ISOEFF = EMXISO * ISOEQ / (EC5ISO + ISOEQ);
```

This is real, single-drug pharmacology (an active metabolite pooled with
its parent by relative potency — the textbook way to handle an active
metabolite), **not** two independent drugs sharing a pathway. The guide's
"never collapse several drugs into one shared Hill term" rule targets the
latter (e.g. two different biologics hitting the same receptor); forcing
ISO's parent+metabolite sum apart into two separately-saturating Hills
would require inventing a decomposition the original's algebra doesn't
support (`Emax*(a+b)/(EC50+a+b)` is not separable into `f(a)+g(b)`) — an
actual refit changing the model's behavior, not a rename. Per the guide's
"a clean, non-standard structure beats a standard structure that's wrong"
and the RVWF/PDVWF and OLAP precedents (`von-willebrand-disease/vwd_refactor_notes.md`,
`breast-cancer/bc_refactor_notes.md`) of not inventing an effect term the
original never had, the decision here is: **keep the combined Hill
exactly as-is, renamed**, and disclose it plainly rather than force a fit.

```
double X_ISO      = C_ISO + POT_OXO * C_OXO;   // parent-equivalent composite
double EFFECT_ISO = EMAX_ISO * pow(X_ISO, GAMMA_ISO)
                     / (pow(EC50_ISO, GAMMA_ISO) + pow(X_ISO, GAMMA_ISO));
```

**Renaming applied (ISO + OXO):**

| Original | Refactored | Value |
|---|---|---|
| `ISOD`/`ISOC`/`ISOP` (compartments) | `GUT_ISO`/`CENT_ISO`/`PERI_ISO` | — |
| `OXOC` (compartment) | `CENT_OXO` | — |
| `KAISO`/`VISO`/`CLISO`/`QISO`/`VISOP` | `KA_ISO`/`V1_ISO`/`CL_ISO`/`Q_ISO`/`V2_ISO` | 0.45 (1/h) / 75.0 (L) / 2.00 (L/h) / 3.00 (L/h) / 90.0 (L) |
| `FISOF`/`FOODEF`/`FOOD`/`LIDOSE` | `F0_ISO`/`FOODEF_ISO`/`FOOD_ISO`/`LIDOSE_ISO` | 0.25 / 1.20 / 1 / 0 |
| `KMET` | `KMET_ISO` | 0.018 (1/h), shared inflow term for OXO |
| `VOXO`/`CLOXO` | `V1_OXO`/`CL_OXO` | 45.0 (L) / 1.10 (L/h) |
| `POTOXO` | `POT_OXO` | 0.35 |
| `CISO`/`COXO` (local `$ODE` ratios) | `C_ISO`/`C_OXO` | mg/L, each independently redirectable |
| `ISOEQ` (local `$ODE`) | `X_ISO` | composite parent-equivalent input to the Hill |
| `EMXISO`/`EC5ISO` | `EMAX_ISO`/`EC50_ISO` | 1.00 / 0.30 (mg/L) |
| — (none; implicit linear ratio) | `GAMMA_ISO` (new, `=1.0`) | rename, not a fit |
| `ISOEFF` (local `$ODE`) | `EFFECT_ISO` | used at `RETEFF`, `LIP`, `SGM` apoptosis, `TG`, `ALT`, `MUCO` |
| `CUMISO` (compartment) | `CUM_ISO` | cumulative mg/kg tracker, feeds the durable-downsizing Hill (`DURTGT`), untouched besides the rename |

No `EFFECT_OXO`: as shown above, OXO's only route to the disease is
through the shared `X_ISO`/`EFFECT_ISO` term — there is nothing
independent to name. `C_OXO` is still exposed as its own clean,
redirectable concentration per the census row.

## SPI — archetype 3 minus the peripheral compartment (depot + central, linear elimination)

```
dxdt_GUT_SPI  = -KA_SPI * GUT_SPI;
dxdt_CENT_SPI =  KA_SPI * GUT_SPI * F_SPI - (CL_SPI / V1_SPI) * CENT_SPI;
```

Three distinct pharmacological actions in the original, not one:

1. **17,20-lyase inhibition** (reduces gonadal androgen synthesis) — a
   plain Hill ratio, renamed:
   `EFFECT_SPI_LYASE = EMAX_SPI_LYASE*pow(C_SPI,g)/(pow(EC50_SPI_LYASE,g)+pow(C_SPI,g))`,
   used at `TTIN`'s one use site.
2. **AR competitive antagonism** — the original's own shape is
   `KDAR*(1 + SCLA/KICLAS + CSPI/KISPI)`: a **linear** term
   (`C_SPI/KI_SPI_AR`) that raises the apparent AR dissociation constant
   proportionally to spironolactone (canrenone) concentration — the
   textbook competitive-antagonist shape, genuinely **not** a saturating
   Hill (no `+C_SPI` in its own denominator; it only ever grows). Renamed
   `KISPI` -> `KI_SPI_AR` and left as the same inline linear expression
   (`C_SPI / KI_SPI_AR`) at its one use site inside `KDEFF` — **not**
   wrapped in a fabricated `EFFECT_SPI_AR` Hill name, because doing so
   would misrepresent an unbounded linear term as a bounded/saturating
   one (same reasoning as the RVWF/PDVWF precedent in
   `von-willebrand-disease/vwd_refactor_notes.md`).
3. **Serum-potassium safety effect** — a plain Hill ratio, renamed:
   `EFFECT_SPI_KAL = EMAX_SPI_KAL*pow(C_SPI,g)/(pow(EC50_SPI_KAL,g)+pow(C_SPI,g))`,
   used at `dxdt_KSER`'s one use site.

**Renaming applied:**

| Original | Refactored | Value |
|---|---|---|
| `SPID`/`SPIC` (compartments) | `GUT_SPI`/`CENT_SPI` | — |
| `KASPI`/`FSPI`/`VSPI`/`CLSPI` | `KA_SPI`/`F_SPI`/`V1_SPI`/`CL_SPI` | 1.00 (1/h) / 0.70 / 60.0 (L) / 2.60 (L/h) |
| `CSPI` (local `$ODE` ratio) | `C_SPI` | mg/L |
| `ELYASE`/`KILYA` | `EMAX_SPI_LYASE`/`EC50_SPI_LYASE` | 0.40 / 1.50 (mg/L) |
| — (none; implicit linear ratio) | `GAMMA_SPI_LYASE` (new, `=1.0`) | rename, not a fit |
| `LYASE` (local `$ODE`) | `EFFECT_SPI_LYASE` | — |
| `KISPI` | `KI_SPI_AR` | 0.95 (mg/L); kept as a linear competitive term, not a Hill |
| `EKAL`/`KIKAL` | `EMAX_SPI_KAL`/`EC50_SPI_KAL` | 0.13 / 1.20 (mg/L) |
| — (none; implicit linear ratio) | `GAMMA_SPI_KAL` (new, `=1.0`) | rename, not a fit |
| (inline, unnamed) | `EFFECT_SPI_KAL` | newly named at its one use site (`dxdt_KSER`) for consistency/discoverability; same math |

## TET — archetype 3 minus the peripheral compartment (depot + central, linear elimination)

```
dxdt_GUT_TET  = -KA_TET * GUT_TET;
dxdt_CENT_TET =  KA_TET * GUT_TET * F_TET - (CL_TET / V1_TET) * CENT_TET;
```

The model's own header comment calls out "TWO separate
concentration-effect relationships" for tetracycline — an antimicrobial
arm and a non-antimicrobial anti-inflammatory arm — so that is exactly
what got two names:

| Original | Refactored | Value |
|---|---|---|
| `TETD`/`TETC` (compartments) | `GUT_TET`/`CENT_TET` | — |
| `KATET`/`FTET`/`VTET`/`CLTET` | `KA_TET`/`F_TET`/`V1_TET`/`CL_TET` | 0.90 (1/h) / 0.90 / 52.0 (L) / 2.20 (L/h) |
| `CTET` (local `$ODE` ratio) | `C_TET` | mg/L |
| `EMXTET`/`EC5TET` | `EMAX_TET_AMR`/`EC50_TET_AMR` | 0.050 (1/h) / 2.50 (mg/L) |
| — (none; implicit linear ratio) | `GAMMA_TET_AMR` (new, `=1.0`) | rename, not a fit |
| `TETKIL` (local `$ODE`, SUBANTI-gated) | `EFFECT_TET_AMR` | bactericidal-rate contribution to `KILL`/`TETKL2` |
| `EMXAI`/`EC5AI` | `EMAX_TET_AI`/`EC50_TET_AI` | 0.22 / 0.35 (mg/L) |
| — (none; implicit linear ratio) | `GAMMA_TET_AI` (new, `=1.0`) | rename, not a fit |
| `AIEFF` (local `$ODE`) | `EFFECT_TET_AI` | anti-inflammatory effect |
| `AI` (capped copy, `min(AIEFF,0.70)`) | `AI_TET` | used at `dxdt_TLR` and `RESACC` |
| `EMMP` | `EMAX_TET_MMP` | 1.40; MMP-degradation-enhancement Emax |

Two further internal use-sites reuse `C_TET` directly with the renamed
EC50 constants, exactly as the original did, rather than being wrapped in
a third/fourth `EFFECT_TET_*` name: the MMP-degradation enhancement term
(`EMAX_TET_MMP * C_TET/(EC50_TET_AI + C_TET)`, deliberately reusing the
anti-inflammatory arm's own EC50, same as the original's `EMMP*CTET/(EC5AI+CTET)`)
and the antibiotic-selection-pressure term inside `PRESS`
(`C_TET/(EC50_TET_AMR + C_TET)`, reusing the antimicrobial arm's EC50,
SUBANTI/NARROW-gated exactly as before) — neither was an independently
named effect in the original, so neither invents one here.

## The `$TABLE` scoping constraint this file's own structure forced

This file (unlike several prior refactors in this batch, e.g.
`von-willebrand-disease/vwd_mrgsolve_model_refactored.R`, which have no
`$TABLE` block at all) has a genuine `$TABLE` block computing `INFLAM`,
`NONINF`, `IGA`, etc. **First attempt:** redeclare `C_TET`/`C_ISO`/
`C_OXO`/`C_SPI`/`C_EE`/`X_ISO`/`EFFECT_ISO`/`EFFECT_TET_AI`/
`EFFECT_TET_AMR`/`EFFECT_EE_SHBG`/`EFFECT_EE_LH`/`EFFECT_SPI_LYASE`/
`EFFECT_SPI_KAL` a second time inside `$TABLE`, mirroring the `$ODE`-side
definitions — **this failed to compile**: mrgsolve concatenates `$ODE`
and `$TABLE` into the *same* generated C++ function, so every one of
those names was a `redefinition of '{anonymous}::NAME'` error (confirmed
via the qspserver `mrgsolve_api` compiler output). This is exactly why
the *original* file's own author picked different, "O"-suffixed names for
its `$TABLE`-side recomputations (`CTETO`, `CISOO`, `COXOO`, `CSPIO`,
`CEEO`, `ISOEQT`, `ISOEFFT`, `AIEFFT`) instead of reusing `CTET`/`CISO`/
etc. directly — not a stylistic choice, a required workaround for this
build's `$ODE`+`$TABLE` scoping.

**Resolution:** left the original's own `CTETO`/`CISOO`/`COXOO`/`CSPIO`/
`CEEO`/`ISOEQT`/`ISOEFFT`/`AIEFFT` table-side recomputations exactly as
they were (only updating the identifiers they read, e.g. `TETC`->
`CENT_TET`), and captured the canonical `C_TET`/`C_ISO`/`C_OXO`/`C_SPI`/
`C_EE`/`EFFECT_ISO`/`EFFECT_TET_AI`/`EFFECT_TET_AMR`/`EFFECT_EE_SHBG`/
`EFFECT_EE_LH`/`EFFECT_SPI_LYASE`/`EFFECT_SPI_KAL` **directly from their
`$ODE`-scope declarations** (no `$TABLE` redeclaration needed or
possible — they are already in scope for `$CAPTURE` once `$TABLE` doesn't
redefine them). Confirmed via `/model_manifest`: all of the above appear
in `outputPaths` alongside the untouched `CTETO`/`CISOO`/`COXOO`/`CSPIO`/
`CEEO`/`ISOEFFT`/`AIEFFT`.

## A display-precision artifact, not a modeling bug (found during verification)

Comparing `C_EE` (mg/L) against the model's own `CEEO` (pg/mL,
`C_EE*1e6`) inside one verification run showed `C_EE` reporting as
exactly `0` at many timepoints where `CEEO` showed a sensible nonzero
curve (e.g. `CEEO=25.4518`, `C_EE` should be `~2.545e-5`). This is **not**
a computation error: `EFFECT_EE_SHBG`/`EFFECT_EE_LH` at the same
timepoints were correctly nonzero and Hill-shaped (e.g. `0.857`, `0.246`),
which is only possible if the underlying `C_EE` used internally was
correctly `~1e-5`-scale, not zero. The qspserver API's JSON output
appears to round/truncate very small captured values (ethinylestradiol
is dosed in micrograms, so its own mg/L concentration is natively
`~1e-5`) below its display resolution — precisely why the *original*
author scaled `CEEO` to pg/mL for reporting in the first place. This is a
reporting-precision artifact of very small native units, not a mismatch
between the original and refactored models (the original-vs-refactored
comparison below never relies on raw `C_EE`; it compares `CEEO`, which is
unaffected).

## Verification

**Method.** Both DSLs (the checked-in original's `code`, extracted
verbatim from its `code <- '...'` R-string, and the refactored file's
equivalent) were POSTed to the local qspserver `mrgsolve_api` service at
`http://localhost:8007` (`/model_manifest` then `/run_simulation`), which
compiles and runs each DSL block directly with mrgsolve 2.0.1
server-side — no local R/mrgsolve install used. The original compiled
cleanly on its own (no pre-existing build defect for this file — unlike
several other files in this batch, no `UPSTREAM_ISSUES.md` entry is
needed here). Requests were spaced ~2.2s apart per the task's concurrency
note.

Five of the file's own scenarios (`ACN_scenarios()`) were reproduced,
using the model's own dosing amounts/intervals/durations and phenotype
parameter presets, dosing directly by compartment index (`GUT_TET`=41,
`RETS`=36, `BPOS`=35, `GUT_ISO`=43, `GUT_SPI`=47, `GUT_EE`=49 — confirmed
identical indices in both models via `/model_manifest`'s `outputPaths`
order, since compartment declaration order is unchanged):

1. **Scenario 8** (`Doxycycline 100 mg + adapalene 0.3%/BPO 2.5%, 12 wk`)
   — full 3-drug dosing (`GUT_TET` 100 mg q24h x84d, `RETS` 0.3 q24h x84d,
   `BPOS` 2.5 q24h x84d), moderate phenotype (defaults), `end=2016h`,
   `delta=24h`. Exercises TET's antimicrobial arm (`SUBANTI=0`) alongside
   the untouched BPO/retinoid blocks.
2. **Scenario 9** (`Sub-antimicrobial doxycycline 40 mg MR + adapalene,
   12 wk`) — `GUT_TET` 40 mg + `RETS` 0.1, `SUBANTI=1`, same window.
   Exercises the `SUBANTI`-gated branch (`EFFECT_TET_AMR` forced to 0).
3. **Scenario 13** (`Spironolactone 100 mg QD, adult female, 24 wk`) —
   `GUT_SPI` 100 mg q24h x168d, adult-female phenotype parameters,
   `end=4032h`, `delta=24h`.
4. **Scenario 12** (`COC EE 30 ug x 6 cycles, adult female`) — 6 separate
   `GUT_EE` dosing blocks (0.03 mg q24h, 21 doses/cycle, cycles 672h
   apart, reproducing the 21-active/7-placebo pattern), adult-female
   phenotype, `end=4032h`, `delta=24h`.
5. **Scenario 14, shortened** (`Isotretinoin 0.5 mg/kg x 34 wk`) — `end`
   shortened from the scenario's own 602 days (including a 12-month
   post-course follow-up) to **238 days (5712h)**, i.e. the dosing period
   only, per the guide's "shorten the window rather than treating [a
   step-count/runtime constraint] as a mismatch" allowance — the 12-month
   follow-up mainly re-tests slow `DURAB` reversal over a very long
   horizon, not additional PK/PD structure, so the shortened window still
   fully exercises `GUT_ISO`/`CENT_ISO`/`PERI_ISO`/`CENT_OXO`/`CUM_ISO`,
   `EFFECT_ISO`, and the `TG`/`ALT`/`MUCO` safety trackers it drives.
   `GUT_ISO` 36 mg (0.5 mg/kg x72 kg) q24h x238d, severe-nodular phenotype,
   `delta=24h`.

Each scenario's own baseline burn-in (`acn_baseline()`, a 4-year
steady-state equilibration the original's `ACN_simulate()` runs before
applying dosing) was **not** reproduced — both models start from the raw
`$INIT` defaults instead. This is a disclosed simplification, not a
scope reduction: the burn-in touches no renamed identifier, both models
receive *identical* treatment (same skip), and every `$PARAM` default is
byte-identical between the two files, so comparing from the same
un-equilibrated starting point remains a fair, like-for-like test of the
renamed PK/PD blocks — it does not test the burn-in itself, which was
never in scope.

**Results.** Compared across every shared `$CAPTURE`/compartment output
(58 disease-side names untouched by any of the 5 compounds' renaming —
`TT` through `FAIO`'s full 55-compartment + `$TABLE` set minus the 11
renamed compartments — plus the 11 renamed compartments themselves,
matched by position: `TETD`/`TETC`->`GUT_TET`/`CENT_TET`,
`ISOD`/`ISOC`/`ISOP`->`GUT_ISO`/`CENT_ISO`/`PERI_ISO`,
`OXOC`->`CENT_OXO`, `SPID`/`SPIC`->`GUT_SPI`/`CENT_SPI`,
`EED`/`EEC`->`GUT_EE`/`CENT_EE`, `CUMISO`->`CUM_ISO`):

- **All 5 scenarios: exact match.** Maximum relative deviation across
  every shared output, every timepoint, was **0.0** in all five runs —
  including the untouched compartments (proving nothing outside the 5
  compounds' own blocks drifted), the renamed PK compartments themselves,
  and the pre-existing `CTETO`/`CISOO`/`COXOO`/`CSPIO`/`CEEO`/`ISOEFFT`/
  `AIEFFT` table-side recomputations (proving the renamed identifiers
  they now read produce identical values to the originals they replaced).
- An additional internal self-consistency check (refactored model only)
  confirmed `C_TET`==`CTETO`, `C_ISO`==`CISOO`, `C_OXO`==`COXOO`,
  `C_SPI`==`CSPIO` exactly (0.0 relative deviation) across all 5
  scenarios, and `EFFECT_ISO`==`ISOEFFT`, `EFFECT_TET_AI`==`AIEFFT`
  exactly as well — confirming the new canonical `$ODE`-scope-captured
  names and the pre-existing `$TABLE`-side duplicates agree perfectly
  (`C_EE` vs `CEEO` showed the display-precision artifact described
  above, not a genuine mismatch).

Per the guide's tolerance table: this is pure structural reorganization
for all 5 compounds (Archetype 3-minus-peripheral for EE/SPI/TET,
Archetype 3 unchanged for ISO, bespoke single-compartment-with-inflow for
OXO), no Hill-refitting anywhere (every effect term was already a plain
ratio in the original), so bit-exact reproduction is the expected and
achieved result — no floating-point-scale tolerance needed.

## qspserver `/model_manifest` discoverability

Confirmed via a live `/model_manifest` call against the refactored DSL
(220 parameters, 81 output paths total):

- Every renamed PK parameter (`KA_TET`, `F_TET`, `V1_TET`, `CL_TET`,
  `KA_ISO`, `F0_ISO`, `FOODEF_ISO`, `FOOD_ISO`, `LIDOSE_ISO`, `V1_ISO`,
  `CL_ISO`, `Q_ISO`, `V2_ISO`, `KMET_ISO`, `V1_OXO`, `CL_OXO`, `POT_OXO`,
  `KA_SPI`, `F_SPI`, `V1_SPI`, `CL_SPI`, `KA_EE`, `F_EE`, `V1_EE`,
  `CL_EE`) and every renamed/new Hill parameter (`EMAX_TET_AMR`/
  `EC50_TET_AMR`/`GAMMA_TET_AMR`, `EMAX_TET_AI`/`EC50_TET_AI`/
  `GAMMA_TET_AI`, `EMAX_TET_MMP`, `EMAX_ISO`/`EC50_ISO`/`GAMMA_ISO`,
  `EMAX_SPI_LYASE`/`EC50_SPI_LYASE`/`GAMMA_SPI_LYASE`, `KI_SPI_AR`,
  `EMAX_SPI_KAL`/`EC50_SPI_KAL`/`GAMMA_SPI_KAL`, `EMAX_EE_SHBG`/
  `EC50_EE_SHBG`/`GAMMA_EE_SHBG`, `EMAX_EE_LH`/`EC50_EE_LH`/
  `GAMMA_EE_LH`) is present in `parameters` with its original numeric
  default, confirmed present (no `[]` empty-parameters response).
- `C_TET`, `C_ISO`, `C_OXO`, `C_SPI`, `C_EE`, `EFFECT_ISO`,
  `EFFECT_TET_AI`, `EFFECT_TET_AMR`, `EFFECT_EE_SHBG`, `EFFECT_EE_LH`,
  `EFFECT_SPI_LYASE`, `EFFECT_SPI_KAL`, and all 11 renamed compartments
  (`GUT_TET`, `CENT_TET`, `GUT_ISO`, `CENT_ISO`, `PERI_ISO`, `CENT_OXO`,
  `GUT_SPI`, `CENT_SPI`, `GUT_EE`, `CENT_EE`, `CUM_ISO`) are confirmed
  present in `outputPaths` and retrievable via `/run_simulation`'s
  `outputs` selection.
- `model_content` was extracted from the `ACN_CODE <- '...'` R-string
  wrapper (straight text extraction, the same code used to build every
  verification request above); no `.cpp` sibling was left behind —
  extraction was in-memory / scratch-only, discarded after use, per the
  workflow guide.

## No pre-existing build defect

Unlike several other files in this batch, `acn_mrgsolve_model.R` compiles
cleanly as-written under mrgsolve 2.0.1 (confirmed via a live
`/model_manifest` call against the untouched original) — no
`UPSTREAM_ISSUES.md` entry was needed for this file.

## Diff scope confirmation

Beyond the renames, effect-term restructuring (adding `pow(...,GAMMA)`
where the original had an implicit linear ratio, `GAMMA=1` throughout —
verified to reproduce the original exactly, see Results above), and the
`$CAPTURE` additions described above, nothing else in the file differs:
every other compound's `$PARAM`/compartment/`$ODE` lines (benzoyl
peroxide, topical retinoid, clindamycin, azelaic acid, dapsone,
clascoterone), the phenotype presets (`ACN_phenotypes`), the scenario
library (`ACN_scenarios`), and the summary/relapse/resistance helper
functions are byte-identical except for the compartment-name strings the
dosing helpers (`acn_tetracycline`, `acn_isotretinoin`,
`acn_spironolactone`, `acn_coc`) and the summary functions
(`ACN_relapse`, `ACN_cumdose_curve`) reference by name (`"TETD"`->
`"GUT_TET"`, `"ISOD"`->`"GUT_ISO"`, `"SPID"`->`"GUT_SPI"`, `"EED"`->
`"GUT_EE"`, `d$CUMISO`->`d$CUM_ISO`), so the refactored file's own
scenarios and helper functions still dose/read the renamed
compartments correctly.
