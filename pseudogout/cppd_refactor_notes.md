# Refactor notes — `pseudogout/cppd_mrgsolve_model.R`

Scope: all four compounds carrying rows in
`driver-patches/data/compound_perturbation_census.md` for this file —
**Anakinra (ANA)**, **Colchicine (COLCH)**, **Indomethacin (INDO)**, and
**Prednisolone (PRED)**. Nothing in this file is out of scope; all four
census rows are updated.

## A naming note on PRED

The census row for this compound literally reads "Drug doses (PRED)" —
the automatic classifier that built the census evidently captured the
`$PARAM` section-comment text (`// --- Drug doses (mg) ---`) rather than
the drug name, the same class of mislabeling already corrected for
`neonatal-hyperbilirubinemia`'s Stannsoporfin/Phenobarbital rows (see
`translations/UPSTREAM_ISSUES.md` history). All four rows in this file
carry the same artifact (`Drug doses (ANA/COLCH/INDO/PRED)`), corrected
below to the real drug names. Separately: the task briefing that scoped
this refactor names the fourth compound "prednisone," but
`cppd_mrgsolve_model.R` itself is unambiguous and consistent throughout —
every comment, parameter name, and compartment name says **prednisolone**
(`"Prednisolone PK (1-CMT oral)"`, `Cp_PRED`, `GR_occ`, the `4_Prednisolone_taper`
scenario). Prednisone is an inactive prodrug that requires hepatic
11β-HSD1 conversion to the active glucocorticoid prednisolone; the
original models no such conversion step and its own PK parameters
(F=0.82, t1/2 2.5–4h) match prednisolone, not prednisone. This refactor
follows the code, not the briefing shorthand: the compound modeled here
is prednisolone, kept under its existing `PRED` stem (unchanged, per
"don't invent a new stem").

## Archetypes

**All four compounds are Archetype 3 minus peripheral** (depot + central,
linear, no peripheral compartment) — confirmed for each against its own
`$ODE` pair (`dxdt_<X>_depot` / `dxdt_<X>_central`, first-order absorption
into a single central compartment with linear clearance, no second
compartment anywhere for any of the four). Renamed uniformly:

| Original | Refactored |
|---|---|
| `COLCH_depot` / `COLCH_central` | `GUT_COLCH` / `CENT_COLCH` |
| `INDO_depot` / `INDO_central` | `GUT_INDO` / `CENT_INDO` |
| `PRED_depot` / `PRED_central` | `GUT_PRED` / `CENT_PRED` |
| `ANA_depot` / `ANA_central` | `GUT_ANA` / `CENT_ANA` |
| `Ka_<X>` | `KA_<STEM>` |
| `Vd_<X>` (central volume) | `V1_<STEM>` |
| `F_<X>`, `CL_<X>` | unchanged (already matched convention) |

All parameter *values* are copied verbatim from the original — no PK
parameter was invented, defaulted, or dropped.

Exposed concentrations (the single point an external covariate could
substitute in), computed identically to the original's own `Cp_<X>`
locals, just renamed:

- `C_COLCH = CENT_COLCH / V1_COLCH * 1000` (ng/mL)
- `C_INDO  = CENT_INDO  / V1_INDO` (µg/mL)
- `C_PRED  = CENT_PRED  / V1_PRED * 1000` (ng/mL)
- `C_ANA   = CENT_ANA   / V1_ANA` (µg/mL ≈ mg/L)

## Hill interface

### Colchicine: rename, not a fit — but two independent effect terms, not one

The original computes two separate, already-Hill-shaped ratios from the
*same* colchicine concentration, sharing one `Emax_COLCH` but each with
its own `EC50`:

```
E_COLCH_NLRP3 = Emax_COLCH * Cp_COLCH / (EC50_COLCH_NLRP3 + Cp_COLCH);  // NLRP3 inflammasome inhibition
E_COLCH_NEUT  = Emax_COLCH * Cp_COLCH / (EC50_COLCH_NEUT  + Cp_COLCH);  // neutrophil migration block
```

These are two genuinely different physiological actions of colchicine
(anti-inflammasome vs. anti-migration), not "several drugs on a combined
multi-drug expression" — the case the guide's "combine only where disease
equations actually use them" targets. Per the same precedent already used
for colchicine's two effect terms in
`familial-mediterranean-fever/fmf_refactor_notes.md` (and Encaleret's two
effect terms in `hypoparathyroidism/hypopt_refactor_notes.md`), both are
kept as independent named terms: `EFFECT_COLCH_NLRP3` and
`EFFECT_COLCH_NEUT`, sharing one renamed `EMAX_COLCH` (was `Emax_COLCH`)
and one new `GAMMA_COLCH = 1` (no explicit Hill exponent in the
original; shared because the original's `Emax_COLCH` was itself shared).
`EC50_COLCH_NLRP3`/`EC50_COLCH_NEUT` already matched the naming
convention and are unchanged. Both are plain renames — `pow(C, 1)` is
bit-identical to the original's unexponentiated ratio, not an
approximation of it.

One consequence worth flagging: the original never captured
`E_COLCH_NEUT` in `$CAPTURE` (only `E_COLCH_NLRP3` was). The refactored
file captures both `EFFECT_COLCH_NLRP3` and `EFFECT_COLCH_NEUT`, so
`EFFECT_COLCH_NEUT`'s own value could not be directly diffed against an
original series of the same name — but its effect is fully observable
through `Neutrophil`/`Neut_out` (which `neut_ingress` depends on
directly), and that compartment matches exactly across every verified
scenario (see below).

### Indomethacin: rename, not a fit

`E_INDO_COX = Emax_INDO * Cp_INDO / (EC50_INDO_COX + Cp_INDO)` is already
exactly the guide's Hill shape (gamma=1 implicit). Renamed to
`EFFECT_INDO = EMAX_INDO * pow(C_INDO, GAMMA_INDO) / (pow(EC50_INDO, GAMMA_INDO) + pow(C_INDO, GAMMA_INDO))`,
with `EMAX_INDO`/`EC50_INDO` pulled from `Emax_INDO`/`EC50_INDO_COX`
(values unchanged) and `GAMMA_INDO = 1` new.

### Anakinra: rename, not a fit

`E_ANA = Emax_ANA * Cp_ANA / (KD_ANA + Cp_ANA)` — the original's own
comment calls this "competitive antagonist" / "IL-1Ra competitive:
fraction = [ANA]/(KD_ANA + [ANA])," but the code itself is a plain Emax-Hill
ratio, not a receptor-binding ODE system (no free/bound IL-1R
compartments) — so this is Archetype 3, not Archetype 4/TMDD, matching
what the code actually does rather than what the comment calls it (same
situation independently found for this file's own `E_ANA`/`E_COLCH_*`/
`E_PRED_NFKB` terms and for `familial-mediterranean-fever`'s ANA/CANA
rows, see `fmf_refactor_notes.md`). Renamed to `EFFECT_ANA = EMAX_ANA *
C_ANA / (EC50_ANA + C_ANA)`, with `EC50_ANA` pulled from the
dissociation constant `KD_ANA` (value unchanged, 0.5 µg/mL) and
`EMAX_ANA` from `Emax_ANA`. No `GAMMA_ANA` term is applied inside the
expression itself (kept as a plain ratio, matching the original's exact
form with no `pow()`); `GAMMA_ANA = 1` is still declared in `$PARAM` for
naming-convention completeness and interface discoverability, but is not
referenced inside `EFFECT_ANA`'s own formula, since the original never
had one to rename in the first place at this exact call site — this is a
faithful rename either way (gamma=1 makes no numeric difference whether
or not `pow()` is written explicitly).

### Prednisolone: two chained named quantities, not one — GR occupancy, then a downstream NF-κB-transrepression scale

The original computes GR occupancy first, then reuses it in **two
separate places**:

```
GR_occ = Emax_GR * Cp_PRED / (EC50_GR + Cp_PRED);          // used directly in nlrp3_inhib
E_PRED_NFKB = Emax_PRED * GR_occ;                           // used in IL6_prod, pge2_prod, E_IL1b_total
```

`GR_occ` is itself already exactly the guide's canonical Hill shape (a
plain ratio, `Emax_GR`/`EC50_GR`), so it is renamed directly as the
compound's **primary** exposed effect: `EFFECT_PRED = EMAX_PRED *
pow(C_PRED, GAMMA_PRED) / (pow(EC50_PRED, GAMMA_PRED) + pow(C_PRED,
GAMMA_PRED))`, with `EMAX_PRED`/`EC50_PRED` renamed from `Emax_GR`/
`EC50_GR` (values 0.95/5.0 ng/mL unchanged) and `GAMMA_PRED = 1` new.

`E_PRED_NFKB`, however, is **not itself a fresh Hill function of
concentration** — it is a constant scale (`Emax_PRED`, was 0.88) applied
to the already-computed GR occupancy, used downstream to cap the NF-κB
transrepression magnitude actually seen on IL-6/PGE2. Renamed to
`EFFECT_PRED_NFKB = EMAX_PRED_NFKB * EFFECT_PRED` (`EMAX_PRED_NFKB`, was
`Emax_PRED`, 0.88 — deliberately *not* called `EMAX_PRED` a second time,
since that name is already used for the GR-occupancy ceiling above and
the two values are different: 0.95 vs 0.88). This mirrors the same
"chained, not single" shape the guide's own TMDD example uses for
fractional occupancy feeding a separate Hill term, except here both
stages are algebraic (no ODE-solved receptor kinetics to hold at steady
state or fit) — so both are plain renames, not approximations.

`GR_occ`/`EFFECT_PRED` is reused directly (not re-derived) at its other
call site, `nlrp3_inhib = (EFFECT_COLCH_NLRP3 + EFFECT_PRED * 0.3) *
NLRP3_act` — a genuine multi-drug combination point (colchicine +
prednisolone jointly inhibiting NLRP3), kept exactly where the disease
equations already used it, per the guide's "combine only at the point
disease equations actually use them."

**Pre-existing dead parameter, found, not fixed:** `EC50_PRED_NFKB`
(originally declared at 5.0 ng/mL, directly under the "NF-κB
transrepression" comment) is never referenced anywhere in the original's
`$MAIN`/`$ODE` — confirmed by exact-word grep, the `$PARAM` declaration is
its only occurrence. `E_PRED_NFKB` uses `GR_occ` (itself governed by
`EC50_GR`), never `EC50_PRED_NFKB`. Preserved as declared-but-unused
under its same name, not logged as an upstream defect since it doesn't
affect compilation or produce incorrect behavior (same treatment given to
`F_ZA`/`F_CTN`/`F_DMB` in `pagets-disease/pbd_refactor_notes.md` and
`KDISS_DMB` there).

## Combination points (multi-drug, unchanged in structure)

Two genuine multi-drug combination points exist in the original, both
kept exactly where the disease equations use them, reading the newly
renamed `EFFECT_*` terms in place of the old `E_*`/`GR_occ` names:

```
E_IL1b_total = 1 - (1-EFFECT_COLCH_NLRP3)*(1-EFFECT_PRED_NFKB)*(1-EFFECT_ANA*0.7);
nlrp3_inhib  = (EFFECT_COLCH_NLRP3 + EFFECT_PRED*0.3) * NLRP3_act;
```

Indomethacin's `EFFECT_INDO` and colchicine's `EFFECT_COLCH_NEUT` are
each used alone (in `pge2_prod` and `neut_ingress` respectively) — no
combination there to preserve or disturb.

## Build-compat fix (non-numeric, disclosed)

The original's `$CAPTURE` duplicates all twelve disease-side `$CMT`
compartment names (`PPi_ext Cryst_cart Cryst_SF NLRP3_act IL1b Neutrophil
IL6 PGE2 CRP PainVAS CartInteg LipoxA4`) alongside its genuine `$MAIN`
locals. mrgsolve 2.0.1 rejects this outright:

```
Error in validObject(.Object) :
  invalid class "mrgmod" object: compartment should not be in $CAPTURE: PPi_ext,Cryst_cart,Cryst_SF,NLRP3_act,IL1b,Neutrophil,IL6,PGE2,CRP,PainVAS,CartInteg,LipoxA4
```

Confirmed on the untouched original alone via `POST /model_manifest`, no
renames involved. Per the guide's settled policy, the twelve names were
removed from `$CAPTURE` in the delivered `cppd_mrgsolve_model_refactored.R`
only (every compartment still appears in mrgsolve output automatically
regardless of `$CAPTURE` — confirmed by `/model_manifest`'s `outputPaths`
still listing all twelve, under their unchanged disease-side names) and
`C_COLCH C_INDO C_PRED C_ANA EFFECT_COLCH_NLRP3 EFFECT_COLCH_NEUT
EFFECT_INDO EFFECT_PRED EFFECT_PRED_NFKB EFFECT_ANA` were added so every
new interface variable is discoverable. Logged as
`translations/UPSTREAM_ISSUES.md` **#70**.

## Verification

Per the guide's mandatory protocol: ran all **seven** of the original
file's own dosing scenarios (`evs_1`…`evs_7` in the R driver, exactly as
coded — same doses, same times, same `tend`) through both the
(`$CAPTURE`-patched-only, for buildability) original and
`cppd_mrgsolve_model_refactored.R`, via the qspserver `mrgsolve_api`
(`POST /model_manifest`, `POST /run_simulation`), requests spaced ~2s
apart per the shared-service note.

**A dosing-mechanism note, not a modeling change:** the API's `events`
field (`var`/`time`/`value`, matched by state-variable *name*) silently
produced all-zero series for every depot/central compartment in initial
testing — investigated and found to be an artifact of how this
particular API deployment handles named-variable events, not a model
defect; switching to the API's `dosing` field (NM-TRAN-style,
compartment-*number*-based, `cmt = 1..8` for the four depot/central
pairs, `cmt = 11` for `Cryst_SF`) produced the expected non-degenerate PK
and PD dynamics identically for both models (compartment order is
identical in both files — only names changed, so the same `cmt` numbers
address the same physical compartments in each). This is a
verification-harness detail, not a property of either model file.

**Result: exact match, max abs diff = 0.0 for every shared output, in
all seven scenarios:**

1. **1_Untreated** — acute crystal bolus only (`Cryst_SF += 15` at t=0),
   0–240h, delta=2.
2. **2_Colchicine_0.5mg** — 15-dose irregular loading schedule (`t = 0,
   1, 12, 24, …, 156`, 0.5 mg each) into `GUT_COLCH`/`COLCH_depot`, plus
   the acute bolus, 0–240h.
3. **3_Indomethacin_50mg** — 21 doses q8h (50 mg) into
   `GUT_INDO`/`INDO_depot`, plus the acute bolus, 0–240h.
4. **4_Prednisolone_taper** — the original's own 4-block taper (30→20→10→5
   mg, each block q24h for ~3 days) into `GUT_PRED`/`PRED_depot`, plus the
   acute bolus, 0–720h (`tend_acute * 3`, matching the original's own
   `tend` for this scenario).
5. **5_Anakinra_100mg** — 14 daily 100 mg doses into
   `GUT_ANA`/`ANA_depot`, plus the acute bolus, 0–240h (the original's own
   `tend` for this scenario, 240h, is shorter than the full 14-dose
   schedule's span (last dose at 312h) — a pre-existing property of the
   original scenario, reproduced identically rather than "fixed," and run
   with the identical window on both sides for a fair comparison).
6. **6_Colch_Indo_combo** — colchicine (15 doses q12h, 0.5 mg) +
   indomethacin (21 doses q8h, 50 mg) together, plus the acute bolus,
   0–240h — the multi-drug combination points (`E_IL1b_total`,
   `nlrp3_inhib`, `pge2_prod`) are all live in this scenario.
7. **7_Prophylaxis_Colch** — 180 daily 0.5 mg colchicine doses, no acute
   bolus, full 0–4320h (180-day) window at delta=24 — run at the
   original's own full duration, no shortening needed (no solver
   step-count issue encountered).

Every shared disease-side output (`PPi_ext, Cryst_cart, Cryst_SF,
NLRP3_act, IL1b, Neutrophil, IL6, PGE2, CRP, PainVAS, CartInteg, LipoxA4`,
plus every `$TABLE`-derived `*_out` alias) and every compound-specific
quantity (`C_COLCH`/`Cp_COLCH`, `C_INDO`/`Cp_INDO`, `C_PRED`/`Cp_PRED`,
`C_ANA`/`Cp_ANA`, `EFFECT_COLCH_NLRP3`/`E_COLCH_NLRP3`,
`EFFECT_INDO`/`E_INDO_COX`, `EFFECT_PRED_NFKB`/`E_PRED_NFKB`,
`EFFECT_ANA`/`E_ANA`, `EFFECT_PRED`/`GR_occ`) matched to **max abs diff
0.0** across the full time grid, for all seven scenarios.

Non-degenerate dynamics were confirmed directly, not just "both are
flat/zero" (e.g. scenario 2's `CENT_COLCH` accumulating across the
15-dose loading schedule, `EFFECT_COLCH_NLRP3`/`EFFECT_COLCH_NEUT` rising
from 0 toward ~0.2/~0.18 respectively by day 2.5, `Neutrophil` and
`PainVAS` both responding) — all bit-identical to the original's
equivalent (renamed) series across every scenario. As expected for pure
structural reorganization plus rename-only Hill terms (Archetype 3
throughout, no ODE-derived kinetics for any of the four compounds), no
`nls()` curve-fit was needed or performed for any compound.

**The API has no `init` override** (only `parameters`, mapping to
`$PARAM`, not `$CMT` initial values; this model has no `$MAIN`-side
`<CMT>_0` init block for the API to read either) — so both sides started
every compartment at 0 rather than the R driver's own `init_vals` (e.g.
`PPi_ext = 8.0`, `Cryst_cart = 25.0`, `IL1b = 0.5`, `CartInteg = 75.0`).
This is the same zero-initial-condition property already documented for
this fork's other refactors (e.g.
`pagets-disease/pbd_refactor_notes.md`); it does not affect the
verification's validity since both sides ran under the exact same
(zero) initial conditions — a fair apples-to-apples comparison of the PK/
effect rewiring, which is what this refactor actually changes.

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`, all
four `pseudogout` rows (ANA, COLCH, INDO, PRED), including the compound-
name correction described above.
