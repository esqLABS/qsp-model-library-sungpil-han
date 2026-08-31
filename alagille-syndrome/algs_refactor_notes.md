# Refactor notes — `alagille-syndrome/algs_mrgsolve_model.R`

## Scope: the census's two listed rows are not drugs; the file actually carries one

`driver-patches/data/compound_perturbation_census.md` carried exactly two rows
for this file, both classified "Normalize duplicate concentration sites, then
redirect":

| Row as listed | What it actually is |
|---|---|
| `BA` | **Not a drug.** Bile acid — the disease's own endogenous, hepatically-synthesised metabolite. Seven `$CMT` compartments (`BAHEP`, `BABIL`, `BADUO`, `BAILE`, `BAPOR`, `BASYS`, `BACOL`) form the enterohepatic circulation the whole model turns on, all produced by the model's own synthesis term (`SYN = SYN0 * CYP`, itself CYP7A1/FXR/FGF19-regulated) and moved by first-order transit/exchange fluxes (`EXCR`, `SPILL`, `RET`, `POV`, `RECAP`, `RENAL`). **No `ev()`/dataset dosing route targets any `BA*` compartment anywhere in the file's 16 scenarios or its R driver code** (confirmed by exact-name grep across the whole file). A drug (the IBAT inhibitor, see below) modulates the bile acid system's *ileal reabsorption fraction*, but bile acid itself is never administered — exactly the same pattern documented in `hypercalcemia-of-malignancy`'s `PTH`/`PTHRP` rows and `graves-disease`'s `T3`/`T4`/`TRAb` rows: an endogenous mediator a drug acts *through*, never itself dosed. |
| `Terminal ileum` | **Not a compound at all.** This is an anatomical site name, lifted verbatim from the `$CMT @annotated` comment text on two different, unrelated compartments — `DILE : IBAT inhibitor, terminal ileum - site of action (ug)` and `BAILE : Terminal ileal luminal bile acid (umol)`. Neither compartment is named "Terminal ileum"; the census's scanner appears to have matched the anatomical phrase inside the comment rather than a compound identifier. Same class of error as `neonatal-hyperbilirubinemia`'s census mislabeling and `essential-tremor`'s "Organ systems (PRM)" row (a process-description phrase mistaken for a drug name). |

Both rows are corrected below (marked "not a drug") rather than force-fitting
a PK block onto a variable with no dosing mechanism, or onto a phrase that
was never a compound identifier in the first place.

**The census had also undercounted this file**, the same failure mode
documented for `prostate-cancer` and `hypercalcemia-of-malignancy`: reading
the model's own `$PARAM`/`$CMT`/`$ODE` blocks and its 16 shipped scenarios
turned up **one** genuine externally-dosed compound the census listed no row
for at all:

| Compound | Stem | Dosed in scenario(s) | Evidence |
|---|---|---|---|
| IBAT (ileal bile acid transporter / ASBT) inhibitor — generic PK/PD shared by maralixibat and odevixibat via a potency multiplier | `IBAT` | sc02 (maralixibat, ICONIC replication), sc03/sc05/sc06/sc06b/sc07/sc09/sc10/sc11/sc12/sc13/sc14 (odevixibat, ASSERT replication and every downstream sweep) | `DEPOT`/`DGUT`/`DILE`/`DSYS` — 4-compartment GI transit chain, driven by a continuous zero-order input (`DOSEKG`, `TSTART`/`TSTOP`), potency-scaled by `POT` (the file's own `ODEV`/`MRX` shorthands: `POT = 1.00`/`0.28` respectively) |

This one compound is refactored below. The two census rows needed a
classification fix (not a drug); no other row existed to correct or confirm.

### What was checked and found NOT to be a PK compound (and why no new row was added for it)

The model's `sc08_drug_panel()` scenario exercises five *other* named
interventions on the same patient (UDCA, cholestyramine, rifampicin,
naltrexone, PEBD), plus an MCT-formula diet arm. Each was checked
specifically because the task flagged this file as a candidate for
under-counted real drugs (maralixibat/odevixibat/UDCA named explicitly).
None of the other five has a PK compartment, a concentration state, or a
timecourse of any kind — each is a bare `0`/`1` on/off flag (or, for
cholestyramine, a fixed 0–1 binding-fraction knob) that switches a
constant-magnitude structural or kinetic multiplier on or off, with no dose
scaling and no `dxdt_` compartment at all:

- **UDCA** (`UDCAF`): the comment says "at 20 mg/kg/day", but no dose
  variable, depot, or concentration exists anywhere for it — `UDCAF` (0/1)
  directly gates two fixed fractional constants (`UDCACH` choleresis gain,
  `UDCAHP` hydrophilicity-cytotoxicity cut) with no dose-response curve.
- **Cholestyramine** (`CHOLB`): a dimensionless 0–1 luminal-binding-fraction
  parameter, set by hand per scenario (`CHOLB = 0.55` in `sc08`/`sc09`); no
  dose, no compartment, no kinetics.
- **Rifampicin** (`RIFF`), **naltrexone** (`NALF`): pure 0/1 flags gating
  fixed fractional effect constants (`RIFITCH`/`RIFHEP`, `NALITCH`).
- **PEBD** (`PEBDF`): a surgical intervention (partial external biliary
  diversion), not a drug at all — 0/1 flag setting a fixed residual-ASBT
  fraction (`PEBDEPS`).
- **MCT formula** (`MCTF`): a dietary intervention, 0/1 flag on a fixed
  energy-recovery constant (`MCTGAIN`).

None of these has "PK" in the sense this guide's naming convention requires
(a compartment + a concentration variable a dose maps into) — there is
nothing to redirect. Consistent with the guide's "don't force a bad fit"
principle, none of the five was refactored or added as a census row; this is
disclosed here rather than silently skipped.

## The original compiles cleanly — no build-compat fix needed

`POST /model_manifest` on the untouched original's own DSL (extracted
verbatim from `algs_code <- '...'`): HTTP 200, 130 `$PARAM @annotated`
parameters and 61 outputs (44 `$CMT` compartments + 17 `$CAPTURE` entries)
listed with no error. No entry was added to `translations/UPSTREAM_ISSUES.md`
— there is no defect to log.

## Archetype — bespoke (none of the guide's four reference shapes fit)

The IBAT inhibitor's PK is a linear, unidirectional three-compartment GI
transit chain (upper-gut depot → mid-gut → terminal ileum), with a small
side branch representing negligible systemic absorption, entirely in raw
mass units (µg) — **there is no volume of distribution or plasma
concentration anywhere in the original for this compound.** Archetypes 2/3
(CL/Q/V) cannot be constructed without inventing a volume parameter the
original never had — the same reasoning `hypercalcemia-of-malignancy`'s
zoledronate (`ZOL`) row gives for keeping `K12`/`K21` as disclosed
micro-rate-constants rather than manufacturing `Q`/`V`. Kept as a disclosed
bespoke transit chain instead:

| Original | Refactored | Role |
|---|---|---|
| `DEPOT` (cmt) | `GUT_IBAT` | upper-gut depot — dosing/entry site (continuous zero-order input, not an `ev()`) |
| `DGUT` (cmt) | `MIDGUT_IBAT` | mid-gut transit — bespoke, no role in the guide's naming table; pure transit, no PD reads it |
| `DILE` (cmt) | `ILE_IBAT` | terminal ileum — **the compound's site of action**, the only compartment PD reads |
| `DSYS` (cmt) | `SYS_IBAT` | systemic — off-target, informational only, never read by PD (same "lung vs plasma"/"tissue vs plasma" pattern as this guide's `copd` ICS/LABA/LAMA rows and `abdominal-aortic-aneurysm`'s doxycycline row) |
| `DOSEKG` | `DOSEKG_IBAT` | dose (µg/kg/day) |
| `POT` | `POT_IBAT` | potency multiplier vs. odevixibat reference (odevixibat 1.00, maralixibat 0.28 — the file's own `ODEV`/`MRX` shorthands) |
| `TSTART`/`TSTOP` | `TSTART_IBAT`/`TSTOP_IBAT` | treatment window (continuous-infusion edge case, per the guide's own "dosing via a zero-order input variable... not an `ev()`" allowance) |
| `KTR1`/`KTR2`/`KTR3` | `KTR1_IBAT`/`KTR2_IBAT`/`KTR3_IBAT` | GI transit rate constants (disclosed micro-rate constants, no volume exists to convert them to `Q`) |
| `FSYS` | `FSYS_IBAT` | fraction diverted to the (PD-inert) systemic branch |
| `KSYSD` | `KSYSD_IBAT` | systemic elimination rate constant |
| `IMAXA` | `EMAX_IBAT` | Hill — maximum achievable fractional ASBT blockade |
| `IC50KG` (a per-kg constant; the actual EC50 is `IC50KG*WT`, WT-scaled at read time) | `IC50KG_IBAT` | Hill — half-maximal ileal drug amount per kg body weight |
| — (none) | `GAMMA_IBAT` (new) | Hill coefficient [math-implied = 1, no explicit exponent in original] |

**Why `EC50_IBAT` is a computed local, not a `$PARAM`:** the original never
declares a fixed-value EC50 parameter — its own `IC50 = IC50KG * WT` is
computed inside `$ODE`/`$TABLE` from `IC50KG` ($PARAM) and `WT` (a
per-patient body-weight covariate, itself a `$PARAM` but one the original
scenarios frequently vary). Turning `EC50_IBAT` into a fixed-value `$PARAM`
would freeze in one patient's body weight and silently break every scenario
that overrides `WT`; kept as a computed local (`EC50_IBAT = IC50KG_IBAT *
WT`), recomputed wherever needed, exactly reproducing the original's own
design choice.

**Why `C_IBAT` is an amount (µg), not a concentration:** the original's own
effect term, `CEFF = POT * DILE`, was already amount-based — no volume term
converts `DILE` (µg) into a concentration anywhere in the file, and `IC50KG`
itself is denominated in µg/kg (an amount per body weight), not a
concentration. Inventing a volume purely to convert `C_IBAT` into a
concentration would add a parameter the original never had, for no
structural gain — so `C_IBAT = POT_IBAT * ILE_IBAT` is kept as the same kind
of quantity the original computed, exposed under the guide's naming
convention rather than silently reinterpreted.

## The Hill interface — exact rename, not a fit

The original's own `INH = IMAXA * CEFF / (IC50 + CEFF)` is already a plain
Emax-style ratio (implicit Hill coefficient = 1, no explicit exponent
anywhere). Refactored to:

```
EFFECT_IBAT = EMAX_IBAT * pow(C_IBAT, GAMMA_IBAT)
              / (pow(EC50_IBAT, GAMMA_IBAT) + pow(C_IBAT, GAMMA_IBAT));
```

with `GAMMA_IBAT = 1.0` (new, explicit) — `pow(x, 1.0)` reproduces `x`
exactly (confirmed bit-identical below, same as every other `GAMMA_* = 1`
rename in this batch). `FREAB_D`/`FREAB`/`FREAB_B` (the disease-side
reabsorption fractions the PEBD surgical arm and the cholestyramine binding
term also feed into) are unchanged in every respect except reading
`EFFECT_IBAT` where the original read `INH` — same values, same shape, same
point of use.

## `$TABLE` recomputation — preserved, not simplified away

This file's own `$TABLE` block carries an explicit design note (defect #6 of
the "NINE DEFECTS FOUND" the original's own header documents): captured
quantities must be *recomputed from compartment state* in `$TABLE`, not read
as leftover locals from the last `$ODE` derivative evaluation, because that
evaluation is not guaranteed to correspond to the reporting timepoint. The
refactor preserves this exactly — `C_IBAT`/`EC50_IBAT`/`EFFECT_IBAT` are
recomputed a second time in `$TABLE`, identically to their `$ODE`
definitions, exactly mirroring how `JCAP`/`SYN`/`EXPO`/etc. are already
handled in the original. Flattening this to a single `$ODE`-only computation
would reintroduce the exact defect the original author found and fixed.

## `$PARAM` vs `$CAPTURE` for `C_IBAT`/`EFFECT_IBAT`

Tried and confirmed working directly (unlike the `$PARAM = 0` route, which
[other refactors in this batch have documented](../pagets-disease/pbd_refactor_notes.md)
fails to compile under mrgsolve 2.0.1 with "assignment of read-only
reference"): `C_IBAT` and `EFFECT_IBAT` are plain `autodec`-typed locals
computed in both `$ODE` and `$TABLE` (same idiom the original already uses
for every other derived quantity in this file), listed bare in `$CAPTURE`
alongside the original's own capture list. Confirmed via `POST
/model_manifest` on the refactored DSL: both appear in `outputPaths`
(61 → 63 total outputs, +`C_IBAT`/+`EFFECT_IBAT`), and all twelve
renamed/new IBAT parameters (`DOSEKG_IBAT, POT_IBAT, TSTART_IBAT,
TSTOP_IBAT, KTR1_IBAT, KTR2_IBAT, KTR3_IBAT, FSYS_IBAT, KSYSD_IBAT,
EMAX_IBAT, IC50KG_IBAT, GAMMA_IBAT`) appear in `parameters` with their
original numeric defaults (130 → 131 total parameters, +`GAMMA_IBAT`).

## Verification

**Method.** Extracted the bare DSL text from both `algs_mrgsolve_model.R`
(`algs_code <- '...'`, unmodified) and `algs_mrgsolve_model_refactored.R`
(same variable name, refactored) and ran them through the local qspserver
`mrgsolve_api` (`http://localhost:8007`), `POST /model_manifest` then `POST
/run_simulation`, requests spaced ~2.2 s apart per the shared-service note.
No `dosing`/`events`/`drivers` fields were needed: both of the file's own
primary trial-replication scenarios dose the IBAT inhibitor as a continuous
zero-order `$ODE` input driven purely by parameters (`DOSEKG`/`POT`/
`TSTART`/`TSTOP`), exactly the guide's "continuous infusion input, not an
event" edge case — so both scenarios were reproduced by passing `parameters`
+ `time: {end, delta}` alone.

**Scenarios run — both of the file's own two primary drug-replication
scenarios, full native windows, no shortening needed (well under the API's
default `maxsteps` at `delta = 1` over 180–438 days):**

1. **`sc02_iconic` — maralixibat 380 µg/kg/day** (`MRX`: `DOSEKG=380,
   POT=0.28`), `TSTART=0, TSTOP=438.3` (ICONIC's own 1.2-year window),
   `PBOMAX=0.7`, `end=438.3, delta=1` (439 points).
2. **`sc03_assert` — odevixibat 120 µg/kg/day, both arms** (`ODEV`:
   `DOSEKG=120, POT=1.00`), `TSTART=0, TSTOP=180`, `PBOMAX=0.8`, `end=180,
   delta=1` (181 points) — run twice, once as the drug arm and once as the
   placebo arm (`DOSEKG=0, TSTART=0, TSTOP=0`), since ASSERT's own design is
   placebo-controlled and the placebo arm exercises the same disease-side
   equations with the IBAT block fully off.

**Result: exact match, max abs diff = 0.0, across all 61 compared outputs,
every scenario** — compared point-by-point across each scenario's full time
grid. The 61 compared quantities per scenario are (a) 57 disease-side
outputs whose names are untouched by this refactor (`BAHEP, BABIL, BADUO,
BAILE, BAPOR, BASYS, BACOL, FGF19, CYP, DPR, DUCTR, HSC, FIB, BILI, GGTX,
ALTX, LPX, XANTH, ATX, SENS, SLEEPD, VITD, VITE, VITK, VITA, IGF1, HTZ, PLT,
CUMBA, HZL, HZQ, HZC, HZV, FECBA, SYNCUM, RETCUM, SBA_OUT, ITCH, ITCH_OBS,
AGE_OUT, EFS, NLS, NLS_POP, EFS_POP, SURV, SURVLIV, TCHOL_OUT, FATABS_O,
FATREL_O, CDUO_OUT, JCAP_OUT, PHI_OUT, INH_OUT, SIG_OUT, EXPO_OUT, PORTP_O,
INR`) and (b) 4 renamed pairs compared old-name-in-original against
new-name-in-refactored (`DEPOT→GUT_IBAT, DGUT→MIDGUT_IBAT, DILE→ILE_IBAT,
DSYS→SYS_IBAT`) — genuinely `0.0` for every point in every scenario, not a
loosened tolerance. This is the expected outcome for a pure structural
reorganization with a `GAMMA_IBAT=1.0` rename (not a fit):
`pow(x,1.0)` reproduced the original's plain ratio bit-identically in every
run, including a direct internal check within the refactored model itself
(`EFFECT_IBAT` vs. the legacy `INH_OUT` diagnostic alias, which the refactor
keeps pointing at `EFFECT_IBAT` for backward-comparable output): max abs
diff 0.0 in all three runs.

`/model_manifest` on the refactored DSL additionally confirmed all twelve
renamed/new IBAT parameters and the six new/renamed IBAT compartments and
covariates (`GUT_IBAT, MIDGUT_IBAT, ILE_IBAT, SYS_IBAT, C_IBAT,
EFFECT_IBAT`) appear correctly (see "`$PARAM` vs `$CAPTURE`" above).

No scratch/debug artifacts were left in the repository — DSL extraction, the
comparison harness, and all intermediate `.cpp`/JSON files ran entirely from
a session-local scratchpad directory outside the repo tree, deleted after
use.

## R driver consistency

Because the IBAT inhibitor's PK parameters (`DOSEKG`, `POT`, `TSTART`,
`TSTOP`, and `IMAXA`) are referenced by name throughout the file's own R
driver (the `ODEV`/`MRX` shorthand lists, all 16 scenario functions, and
`eq_init()`/`sim_eq()`'s parameter-forwarding logic), every one of those
references was renamed in lockstep in the refactored sibling so the shipped
R driver still runs correctly end-to-end (`DOSEKG→DOSEKG_IBAT`,
`POT→POT_IBAT`, `TSTART→TSTART_IBAT`, `TSTOP→TSTOP_IBAT`,
`KTR1/2/3→KTR1_IBAT/KTR2_IBAT/KTR3_IBAT`, `FSYS→FSYS_IBAT`,
`KSYSD→KSYSD_IBAT`, `IMAXA→EMAX_IBAT`, `IC50KG→IC50KG_IBAT`); the one
`ev(..., cmt = "DEPOT", ...)` string in `sc02b_bolus_equivalence()` (a
secondary, non-primary bolus-vs-continuous equivalence check, not one of the
two scenarios used for verification above) was updated to
`cmt = "GUT_IBAT"` for the same reason. A stray data-frame column label
(`IMAXA_assumed` in `sc06_phi_identifiability()`'s output) was also renamed
to `EMAX_IBAT_assumed` for consistency — cosmetic only, not read anywhere
else in the file.

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`: the `BA`
and `Terminal ileum` rows reclassified "not a drug" (endogenous metabolite;
anatomical-site phrase, respectively); one new row added for the IBAT
inhibitor (maralixibat/odevixibat), marked refactored/verified per this
file.
