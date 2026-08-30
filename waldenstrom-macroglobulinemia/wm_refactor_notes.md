# Refactor notes — `waldenstrom-macroglobulinemia/wm_mrgsolve_model.R`

Four-compound file: **Ibrutinib (IBR)**, **Rituximab (RTX)**, **Venetoclax
(VEN)**, **Zanubrutinib (ZAN)** — all four rows in
`driver-patches/data/compound_perturbation_census.md` were classified
"Redirect concentration (clean single site)". The disease model (BTK
occupancy, NF-κB, LPC/PC tumor burden, IgM, haemoglobin, viscosity, BCL-2,
CD20, apoptosis, NK cells) is shared machinery the four compounds feed
into and is left untouched apart from reading the new `EFFECT_<STEM>`
names in place of the original's local `IBR_effect`/`ZAN_effect`/
`RTX_CD20_inh`/`VEN_BCL2_inhib` variables (same arithmetic, same order of
operations). The two other compounds in this file, **bortezomib**
(`BOR_dose`/`BOR_Cp`/`EC50_BOR`/`kout_Prot`/`kprot_kill`/`Protsm`) and
**bendamustine** (`BENDA_kLPC`/`BENDA_flag`), are out of scope and are
byte-identical to the original.

The original file did not compile at all under mrgsolve 2.0.1, and once
patched to compile, crashed at runtime for every scenario — two
independent, unrelated defects, both fixed syntax-only and disclosed in
the delivered `_refactored.R`, logged as
[`translations/UPSTREAM_ISSUES.md` #73](../translations/UPSTREAM_ISSUES.md).
See "Pre-existing defects" below.

## Archetype per compound

All four compounds' effect terms were already plain, one-line Hill ratios
in the original (`C/(C+EC50)`, implicit `Emax=1`, implicit `gamma=1`) — no
ODE-solved receptor-binding submodel for any of them, despite the file's
own header comment describing rituximab's PK as "2-compartment TMDD
simplified" (it is not; see the Rituximab section below). Every
`EFFECT_<STEM>` in this refactor is therefore a **rename, not a fit**.

### Ibrutinib (IBR) — archetype 3 without a peripheral compartment (depot + central, linear elimination)

PK: `GUT_IBR`→`CENT_IBR` (was `IBR_gut`/`IBR_C`), one-compartment linear
elimination (`CL_IBR`/`V1_IBR`, was `CL_ibr`/`V_ibr`) — archetype 3 with
the peripheral compartment dropped, same allowance the guide gives for a
no-peripheral 2-compartment variant, applied to depot+central here.

`C_IBR` (nM) is the exposed concentration, derived in two steps exactly as
the original did (`C_IBR_NGML` then `/MW_IBR*1000`) rather than collapsed
into one line, so the intermediate ng/mL diagnostic the original's own
`$TABLE`/`$CAPTURE` reported (`Cp_IBR`) is preserved unchanged.

**Hill interface — exact rename, no fit.** The original's `IBR_effect =
Cp_ibr_nM / (Cp_ibr_nM + EC50_ibr_btk)` is already a plain Hill ratio.
`EFFECT_IBR = EMAX_IBR * pow(C_IBR, GAMMA_IBR) / (pow(EC50_IBR, GAMMA_IBR)
+ pow(C_IBR, GAMMA_IBR))`, with `EMAX_IBR = 1` and `GAMMA_IBR = 1` both
new parameters (the original had neither explicit — the ratio already
saturates at 1 with no Hill exponent).

`EFFECT_IBR` feeds into the **shared, disease-level** `BTK_occ`
compartment exactly where `IBR_effect` did (`BTK_input = fmax(EFFECT_IBR,
EFFECT_ZAN)`) — `BTK_occ` itself is left unrenamed, since ibrutinib and
zanubrutinib compete for the same target and neither compound "owns" this
compartment (same category of judgment call as leaving disease-level
states like `NIS`/`ALC` unrenamed in other refactors in this corpus).
This is also exactly the guide's "multiple drugs, one pathway" case: each
compound's `EFFECT_<STEM>` stays a separate, independently-computed named
term, combined only at the point the shared compartment consumes them.

### Zanubrutinib (ZAN) — same archetype as ibrutinib

PK: `GUT_ZAN`→`CENT_ZAN` (was `ZAN_gut`/`ZAN_C`), `CL_ZAN`/`V1_ZAN` (was
`CL_zan`/`V_zan`). `C_ZAN` (nM) via the same two-step ng/mL→nM derivation
as ibrutinib. `EFFECT_ZAN` is the same rename pattern (`EMAX_ZAN = 1`,
`GAMMA_ZAN = 1` both new), feeding the same shared `BTK_occ` compartment
via `fmax(EFFECT_IBR, EFFECT_ZAN)`, unchanged from the original's
`fmax(IBR_effect, ZAN_effect)`.

### Venetoclax (VEN) — same archetype as ibrutinib/zanubrutinib

PK: `GUT_VEN`→`CENT_VEN` (was `VEN_gut`/`VEN_C`), `CL_VEN`/`V1_VEN` (was
`CL_ven`/`V_ven`), `MW_VEN` (was `VEN_MW`). `C_VEN` (µM) is a direct
rename of the original's one-line `Cp_ven_uM` formula (no intermediate
unit needed here — the original's own `$CAPTURE`d unit already matched
what the Hill term reads).

**Hill interface — exact rename, no fit.** `VEN_BCL2_inhib = Cp_ven_uM /
(Cp_ven_uM + EC50_VEN)` → `EFFECT_VEN = EMAX_VEN * pow(C_VEN, GAMMA_VEN) /
(pow(EC50_VEN, GAMMA_VEN) + pow(C_VEN, GAMMA_VEN))`, `EMAX_VEN = 1` and
`GAMMA_VEN = 1` both new (`EC50_VEN` itself already matched the naming
convention, no rename needed). `EFFECT_VEN` is used in exactly the two
places `VEN_BCL2_inhib` was (the `BCL2` turnover equation and the
`apo_drug` composite), same arithmetic, same coefficients (`* 0.5` in
`apo_drug`), same order.

### Rituximab (RTX) — archetype 1 (single compartment, linear elimination) plus a bespoke, disclosed extra loss term

PK: `CENT_RTX` (was `RTX_C`), `CL_RTX` (unchanged name), `V1_RTX` (was
`V_RTX`) — a plain single-compartment linear-elimination PK, **not** the
"2-compartment TMDD simplified" the original's own header comment
(line 80) claims. The original's own `dxdt_RTX_C` line makes this explicit
in its own inline comment: `-(CL_RTX/V_RTX)*RTX_C - k12_RTX*RTX_C; //
(peripheral handled implicitly; simplified 1-cmpt for this model)`.

**`K12_RTX` (was `k12_RTX`) is a bespoke, non-standard term with no slot
in the guide's naming convention**, and it is carried over **exactly
unchanged** rather than renamed to `Q_<STEM>` (which would misleadingly
imply a genuine peripheral compartment exists): it is an *extra*
first-order loss term applied directly to the central compartment, with
no destination compartment and no return flux. The original also declares
`k21_RTX = 0.010` (a return-rate constant implying the author once
intended a real two-compartment model) but **never uses it anywhere in
`$ODE`** — confirmed by grep across the whole file. Both `K12_RTX` (kept,
renamed only for casing) and `k21_RTX` (kept, unrenamed, dead) are
preserved verbatim rather than "fixed" into a proper two-compartment
model or silently dropped — this is a defect in the original's own
pharmacology (a mass-imbalance term dressed as inter-compartmental
clearance, plus one genuinely dead parameter), not something this
refactor is in scope to correct, per the guide's "log what you find,
don't fix it upstream" rule. Not logged as a separate `UPSTREAM_ISSUES.md`
entry (it doesn't block compilation or change reported behavior once
renamed) but flagged here for a reviewer's attention.

`C_RTX` (mg/L) is a direct rename of the original's `Cp_RTX = RTX_C /
V_RTX`.

**Hill interface — exact rename, no fit.** `RTX_CD20_inh = Cp_RTX /
(Cp_RTX + EC50_RTX_CD20)` → `EFFECT_RTX = EMAX_RTX * pow(C_RTX,
GAMMA_RTX) / (pow(EC50_RTX, GAMMA_RTX) + pow(C_RTX, GAMMA_RTX))`,
`EC50_RTX` (was `EC50_RTX_CD20`), `EMAX_RTX = 1`/`GAMMA_RTX = 1` both new.
`EFFECT_RTX` is used in exactly the three places `RTX_CD20_inh` was: the
`CD20` turnover equation, and the ADCC kill terms in both `dxdt_LPC` and
`dxdt_PC` (`kADCC * NK * CD20 * EFFECT_RTX * ...`) — same arithmetic, same
order, same coefficients (including the `* 0.5` scaling on the `PC` ADCC
term).

## Dead parameters preserved, not removed

`IBR_dose`, `ZAN_dose`, `RTX_dose`, `VEN_dose`, and `k_RTX_on` are all
declared in `$PARAM` and never read anywhere in `$MAIN`/`$ODE`/`$TABLE` in
the original (confirmed by grep) — dosing is driven entirely by `ev()`
events in the R driver, not by these placeholder parameters. Same
treatment as `use_IB`/`use_VEN`/`use_OBI` in
`chronic-lymphocytic-leukemia/cll_refactor_notes.md`: kept, unrenamed,
disclosed as dead rather than silently deleted.

## Renaming summary (all parameter values copied verbatim from the original)

| Role | IBR | ZAN | RTX | VEN |
|---|---|---|---|---|
| Depot | `IBR_gut`→`GUT_IBR` | `ZAN_gut`→`GUT_ZAN` | none (IV) | `VEN_gut`→`GUT_VEN` |
| Central | `IBR_C`→`CENT_IBR` | `ZAN_C`→`CENT_ZAN` | `RTX_C`→`CENT_RTX` | `VEN_C`→`CENT_VEN` |
| Absorption/bioavail. | `ka_ibr`→`KA_IBR`, `F_ibr`→`F_IBR` | `ka_zan`→`KA_ZAN`, `F_zan`→`F_ZAN` | n/a | `ka_ven`→`KA_VEN`, `F_ven`→`F_VEN` |
| Clearance/volume | `CL_ibr`→`CL_IBR`, `V_ibr`→`V1_IBR` | `CL_zan`→`CL_ZAN`, `V_zan`→`V1_ZAN` | `CL_RTX` (unchanged), `V_RTX`→`V1_RTX` | `CL_ven`→`CL_VEN`, `V_ven`→`V1_VEN` |
| Bespoke/dead | — | — | `k12_RTX`→`K12_RTX` (bespoke, kept), `k21_RTX` (dead, kept) | — |
| Molecular weight | `MW_ibr`→`MW_IBR` | `MW_zan`→`MW_ZAN` | n/a | `VEN_MW`→`MW_VEN` |
| Affinity | `EC50_ibr_btk`→`EC50_IBR` | `EC50_zan_btk`→`EC50_ZAN` | `EC50_RTX_CD20`→`EC50_RTX` | `EC50_VEN` (unchanged) |
| Hill Emax/gamma | new `EMAX_IBR=1`, `GAMMA_IBR=1` | new `EMAX_ZAN=1`, `GAMMA_ZAN=1` | new `EMAX_RTX=1`, `GAMMA_RTX=1` | new `EMAX_VEN=1`, `GAMMA_VEN=1` |
| Exposed concentration | new `C_IBR` (nM) | new `C_ZAN` (nM) | new `C_RTX` (mg/L) | new `C_VEN` (µM) |
| Compound effect | new `EFFECT_IBR` | new `EFFECT_ZAN` | new `EFFECT_RTX` | new `EFFECT_VEN` |
| Dead dose placeholder | `IBR_dose` (unchanged, dead) | `ZAN_dose` (unchanged, dead) | `RTX_dose`, `k_RTX_on` (unchanged, dead) | `VEN_dose` (unchanged, dead) |

**Left unrenamed, deliberately (disease-level, not compound-specific):**
`BTK_occ`, `NFkB`, `LPC`, `PC`, `IgM`, `Hgb`, `Visc`, `BCL2`, `CD20`,
`Apop`, `NK` (compartments, plus `Protsm` which belongs to
out-of-scope bortezomib); `kout_BTK`, `NFkB_base`, `MYD88_drive`,
`kBTK_NFkB`, `EC50_BOR_NFkB`, `kout_NFkB`, `kel0_LPC`, `kprolif_LPC`,
`NFkB_kLPC`, `LPC0`, `kconv_LPC`, `kel0_PC`, `PC0`, `KMAX_BM`, `ksec_IgM`,
`kel_IgM`, `IgM0`, `Hgb0`, `kprod_Hgb`, `kel_Hgb`, `BMInf_Hgb`,
`Visc_base`, `kIgM_visc`, `Visc_exp`, `BCL2_base`, `kNFkB_BCL2`,
`kout_BCL2`, `kout_CD20`, `kADCC`, `NK0`, `kprod_NK`, `kel_NK` (disease
parameters); `BOR_dose`, `BOR_Cp`, `EC50_BOR`, `kout_Prot`, `kprot_kill`,
`BENDA_kLPC`, `BENDA_flag` (bortezomib/bendamustine, out of scope,
byte-identical to the original).

## Pre-existing defects found (both fixed, syntax-only, in the delivered `_refactored.R`; logged as `translations/UPSTREAM_ISSUES.md` #73)

Neither of these two is related to any of the four compounds' own PK/PD —
both block the file from running at all, original or refactored, until
worked around.

**1. `$ODE` assigns directly to a `$CMT` state (`BMInf`) to hold it at an
algebraic value** (`BMInf = fmin(total_tumor/KMAX_BM, 1.0); dxdt_BMInf =
0;`), which mrgsolve 2.0.1 rejects at compile time
(`assignment of read-only reference 'BMInf'`). Fix: `BMInf` demoted from
a `$CMT` state to a plain `double` local inside `$ODE` (no `dxdt_` line,
since it is no longer a state); `$TABLE`'s `BMInf_pct` recomputed directly
from `LPC`/`PC` with the identical formula. This removes one raw
auto-reported output column (`BMInf`, which mrgsolve auto-includes for
every `$CMT` state regardless of `$CAPTURE`) — but `BMInf_pct` (÷100 of
the same value) was already the only form of this quantity in the
original's own `$CAPTURE` list, so no information is lost.

**2. `_F(CMT) = value;` (used for ibrutinib/zanubrutinib/venetoclax
bioavailability) compiles but crashes `POST /run_simulation` at runtime**
(`munmap_chunk(): invalid pointer`, signal 6) for every scenario tested,
dosed or not — confirmed via a minimal bisected repro (a bare
depot+central 2-compartment model with only `_F(GUT) = F1;` added, no
dosing) reproducing the identical crash on its own. Fix: replaced with
the modern `F_<CMT> = value;` idiom (`F_GUT_IBR = F_IBR;` etc.), confirmed
to apply bioavailability identically (steady-state central-compartment
concentration scales by exactly the `F` value either way).

Both fixes were required simultaneously to get a runnable model through
the API; see `translations/UPSTREAM_ISSUES.md` #73 for the full account,
including the exact compiler error and crash bisection.

## Verification

**Method.** Both DSLs' bare quoted `code <- '...'` blocks were extracted
and run through the qspserver `mrgsolve_api` container (`POST
/model_manifest`, `POST /run_simulation`), spaced ~2–3 s apart per the
container's known concurrency limits. `POST /model_manifest` confirmed
both compile; the refactored model exposes 75 `$PARAM` entries (the
original's 67, plus 8 new: `EMAX_IBR/ZAN/RTX/VEN` and `GAMMA_IBR/ZAN/RTX/
VEN`) and 46 output paths (19 compartments — one fewer than the original's
20 due to `BMInf`'s demotion, see defect 1 — plus 27 `$CAPTURE` entries,
the original's own 19 plus 8 new discoverability diagnostics: `C_IBR`,
`C_ZAN`, `C_RTX`, `C_VEN`, `EFFECT_IBR`, `EFFECT_ZAN`, `EFFECT_RTX`,
`EFFECT_VEN`). Per `POST /model_manifest`, none of the eight new
`C_<STEM>`/`EFFECT_<STEM>` names are listed as `$PARAM` — confirmed by
direct test that adding any of them to `$PARAM` alongside their existing
`$ODE`-local recomputation reproduces the exact same
`assignment of read-only reference` compile error as defect 1 above (a
name cannot simultaneously be a settable `$PARAM` default and a value
freshly recomputed every `$ODE` step in this build). Same confirmed
constraint as `chronic-lymphocytic-leukemia/cll_refactor_notes.md`. All
eight are declared as `double` locals inside `$ODE` and exposed through
`$CAPTURE` only.

Ran all seven of the original's own `run_scenario()` dosing schedules
(amounts, `ii`/`addl`, exactly as computed by the original's own
`ibr_ev()`/`zan_ev()`/`rtx_ev()` helpers and the Venetoclax ramp
construction) through both models, using a **shortened 120-day horizon**
(`delta=24`) rather than the original's own 730-day default — the
original's own scenario functions accept a `days` argument and were run
with `days=120` equivalent dosing, honestly shortened per the guide's
allowance rather than risking an API step-budget timeout on the full
730-day/many-event Venetoclax and Rituximab-cycle schedules; the shorter
window still exercises every compound's full absorption/elimination
cycle, the venetoclax ramp (all 5 dose levels, days 0–27, plus 92 days of
288 mg maintenance), and 5 of rituximab's 6 q28d cycles (day 112 falls
just inside the window; the 6th, at day 140, does not):

1. **Watch_Wait** — no dosing (sanity check untouched code paths)
2. **Ibrutinib** — 420 mg QD × 120 d (`amt=42` into the depot, `F` already
   folded into the dose per the original's own comment, `ii=24, addl=119`)
3. **iR** — ibrutinib QD + rituximab 675 mg at 5 cycle times (0, 672,
   1344, 2016, 2688 h)
4. **Zanubrutinib** — 96 mg BID × 120 d (`ii=12, addl=239`)
5. **R_Benda** — rituximab alone, same 5 cycle times as iR (bendamustine
   itself is out of scope and `BENDA_flag` stays at its default 0 in both
   models, matching the original's own scenario)
6. **BDR** — same rituximab cycles as R_Benda, plus `BOR_Cp=200` set as a
   constant parameter (not a dosing event, matching the original's own
   `extra_p` mechanism; bortezomib itself is out of scope)
7. **Venetoclax** — the 5-step ramp (14.4→36→72→144→288 mg, `*0.72` F
   already folded into each amount per the original), each step dosed QD
   for 7 days except the final maintenance step, which continues to the
   end of the shortened window

**Result: exact match. Maximum absolute deviation = 0.0 across all seven
scenarios**, for every one of the 19 shared `$CAPTURE` outputs
(`Cp_IBR`, `Cp_ZAN`, `Cp_RTX_mgl`, `Cp_VEN_uM`, `BTK_pct`, `NFkB_AU`,
`BCL2_AU`, `LPC_cells`, `PC_cells`, `IgM_gL`, `Hgb_gdL`, `BMInf_pct`,
`Visc_cP`, `NK_cells`, `Protsm_AU`, `CD20_frac`, `Apop_AU`, `HVS_flag`,
`IgM_change_pct`) and all 7 renamed PK compartments checked directly
(`IBR_gut`/`GUT_IBR`, `IBR_C`/`CENT_IBR`, `ZAN_gut`/`GUT_ZAN`,
`ZAN_C`/`CENT_ZAN`, `RTX_C`/`CENT_RTX`, `VEN_gut`/`GUT_VEN`,
`VEN_C`/`CENT_VEN`), across the full 121-point daily time grid for every
scenario. This is the expected outcome per the guide's tolerance table
(archetypes 1/3, pure structural reorganization, no Hill-fitting — every
one of the four compounds' effect terms was already an exact one-line
Hill ratio).

**One self-caught issue during verification, corrected before the exact
match above:** an earlier draft of `$TABLE` recomputed `Cp_IBR`/`Cp_ZAN`/
`Cp_RTX_mgl`/`Cp_VEN_uM` fresh from the state compartments
(`CENT_RTX/V1_RTX`, etc.) rather than aliasing the same `$ODE`-local
variables the way the original's own `$TABLE` did (`double Cp_RTX_mgl =
Cp_RTX;`, where `Cp_RTX` was computed once in `$ODE`). This produced a
one-row mismatch in the `iR`/`R_Benda`/`BDR` scenarios at exactly the
timestamps where a dose coincides with a reporting time (e.g. `Cp_RTX_mgl`
= 150 in the fresh-recompute draft vs. 0 in the original, at `t=0`
immediately after the rituximab bolus): the original's own `$TABLE`
reuses the `$ODE`-local value from the *last integrator step before* the
dose event, so it under-reports the concentration at that one instant,
one row before catching up. Reverting `$TABLE` to alias the same
`$ODE`-locals (`C_IBR_NGML`, `C_ZAN_NGML`, `C_RTX`, `C_VEN`) exactly as
the original did — reproducing this reporting quirk rather than "fixing"
it — restored the exact match reported above. Not logged as a separate
`UPSTREAM_ISSUES.md` entry (it is a minor `$TABLE`-timing nuance specific
to how `$ODE`-local doubles are captured, not a build/compile defect, and
it never affects the underlying integrated state, only one reported row
of one diagnostic column).

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`, four
rows under `waldenstrom-macroglobulinemia`: `IBR` (target BTK, shared
`BTK_occ` compartment with ZAN, `C_IBR`/`EFFECT_IBR`, nM), `RTX` (target
CD20, `C_RTX`/`EFFECT_RTX`, mg/L), `VEN` (target BCL-2,
`C_VEN`/`EFFECT_VEN`, µM), `ZAN` (target BTK, shared `BTK_occ` compartment
with IBR, `C_ZAN`/`EFFECT_ZAN`, nM).
