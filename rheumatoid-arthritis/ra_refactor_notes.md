# Refactor notes — `rheumatoid-arthritis/ra_mrgsolve_model.R`

Calibration run: tocilizumab (TCZ), full TMDD case. Only TCZ's PK + IL-6R
binding block and its interface into the disease equations were touched.
Adalimumab, methotrexate, and baricitinib are byte-for-byte identical to the
original in `ra_mrgsolve_model_refactored.R`.

## Archetype

**Archetype 4 (TMDD)**, confirmed against the actual equations, not assumed
from the file's own "2-compartment TMDD" comment:

- `SC_TCZ`/`CENT_TCZ`/`PERI_TCZ` — depot + central + peripheral PK
  (renamed `GUT_TCZ`/`CENT_TCZ`/`PERI_TCZ`)
- `IL6R_FREE`/`DRUG_RL` — free receptor / drug-receptor complex, each with
  its own ODE (renamed `REC_FREE_TCZ`/`COMPLEX_TCZ`)
- Binding terms (`KON*C_TCZ*IL6R_FREE`, `KOFF*DRUG_RL`) appear symmetrically
  in the central-compartment mass balance and both receptor ODEs — genuine
  receptor-binding kinetics, not an algebraic occupancy shortcut. This
  matches archetype 4 exactly; nothing was flattened.

One structural difference from the guide's archetype-4 template, kept as-is
because the original genuinely has it: the template's receptor turnover is
`KDEG_TCZ*(RTOT_TCZ - REC_FREE_TCZ)` (a restoring term against a *fixed*
pool `RTOT_TCZ`). This file instead uses classic Mager-Jusko-style turnover —
constant zero-order synthesis `KSYN_R` and independent first-order
degradation `KDEG_R` on *both* free receptor and complex. The two are
mathematically related here: at steady state the total pool
`REC_FREE_TCZ + COMPLEX_TCZ` converges to `KSYN_R/KDEG_R = 0.1/0.007 =
14.29`, which is why the original's `RTOT0 = 14.3` baseline parameter and
init condition line up. Renamed to `KSYN_TCZ`/`KDEG_TCZ`, structure
otherwise untouched.

## Renaming (tocilizumab only)

| Original | Refactored |
|---|---|
| `SC_TCZ` | `GUT_TCZ` |
| `IL6R_FREE` | `REC_FREE_TCZ` |
| `DRUG_RL` | `COMPLEX_TCZ` |
| `F_SC` | `F_TCZ` |
| `KON` | `KON_TCZ` |
| `KOFF` | `KOFF_TCZ` |
| `KSYN_R` | `KSYN_TCZ` |
| `KDEG_R` | `KDEG_TCZ` |
| `RTOT0` | `RTOT_TCZ` |
| `OCC_TCZ` | `EFFECT_TCZ` |

`CENT_TCZ`, `PERI_TCZ`, `CL_TCZ`, `V1_TCZ`, `V2_TCZ`, `Q_TCZ`, `KA_TCZ`,
`C_TCZ` already matched the convention and are unchanged. All parameter
*values* are copied verbatim from the original — nothing invented.

Added (new, not present in the original, for the named Hill interface —
see below): `EMAX_TCZ`, `EC50_TCZ`, `GAMMA_TCZ`.

## Hill interface: rename, not a fit — and why

`OCC_TCZ = DRUG_RL/(RTOT+1e-6)` in the original is used directly, verbatim,
as the compound's fractional block on IL-6 signaling and on the combined
`DRUG_ANTI_INFLAM` term — there is no separate Emax/EC50 layer on top of it
in the original. So the real question was not "does a further Hill
transform fit occupancy" but "is occupancy itself, as a function of
concentration, already Hill-shaped." It is, and not approximately:

**Analytic derivation.** At quasi-steady-state, solving the two linear ODEs
for `REC_FREE_TCZ`/`COMPLEX_TCZ` at fixed `C_TCZ` gives
`RTOT(C) = REC_FREE_TCZ + COMPLEX_TCZ = KSYN_TCZ/KDEG_TCZ` (constant,
independent of `C`), and
`Occ(C) = COMPLEX_TCZ/RTOT(C) = C / (C + Kd)`, `Kd = (KOFF_TCZ+KDEG_TCZ)/KON_TCZ`
— an exact Hill function with `gamma = 1`.

**Numerical confirmation.** A standalone probe model (fixed `C_TCZ`, run to
steady state, mrgsolve via the qspserver `mrgsolve_api`) sampled 18
log-spaced concentrations from ~0.001 to 2000 nmol/L — spanning what the
original file's own TCZ IV, TCZ SC, and dose-response (2-10 mg/kg)
scenarios actually produce (~1e-6 to ~1271 nmol/L). `nls()` fit of
`Emax*C^gamma/(EC50^gamma+C^gamma)`:
  - **gamma fixed at 1**: Emax = 1.000000, EC50 = 7.878788 nmol/L,
    R² = 1.0000, max abs residual ≈ 2.8e-10.
  - **gamma free**: converges to gamma = 1.000000 exactly, same Emax/EC50,
    R² = 1.0000.
  - Analytic `Kd = (0.019+0.007)/0.0033 = 7.878788` matches the fitted EC50
    to 6 significant figures.

Given this, the relationship is a **rename, not a fit**. `EFFECT_TCZ` in
the refactored file still computes from the exact ODE state ratio
(`COMPLEX_TCZ/(REC_FREE_TCZ+COMPLEX_TCZ+1e-6)`), not from the closed-form
`Emax*C^gamma/(EC50^gamma+C^gamma)` approximation — this adds zero
approximation error on top of an already-exact relationship, and is why
verification below matches the original exactly rather than "small and
disclosed." `EMAX_TCZ=1`, `EC50_TCZ=7.878788`, `GAMMA_TCZ=1` are still
exposed as named `$PARAM` values (derived, not fit-to-noise) purely so a
future external-covariate substitution has the closed form ready without
needing the receptor compartments.

## `DRUG_ANTI_INFLAM` (shared multi-drug expression)

```
double DRUG_ANTI_INFLAM = 1 - (1-EFFECT_TCZ)*(1-TNF_INH_ADA)*(1-MTX_EFF_INFLAM)*(1-0.6*JAK_INH);
```

Only the first factor (`OCC_TCZ` → `EFFECT_TCZ`) was touched. `TNF_INH_ADA`
(adalimumab), `MTX_EFF_INFLAM` (methotrexate), `JAK_INH` (baricitinib, with
its own `0.6` scaling untouched), and the multiplicative-independence
combination structure itself are all byte-identical to the original. The
same substitution was made everywhere else `OCC_TCZ` was read: `IL6_sig`
and the reported `RO_pct`.

## A pre-existing defect found while verifying (not fixed, logged upstream)

The original file's `$CAPTURE` block lists `CRP INFLAM VDH CART STAT3_P
FLS_ACT` alongside the derived doubles — but those six names are also
`$CMT` compartments. mrgsolve 2.0.1 (the version installed locally and on
the qspserver `mrgsolve_api` container) refuses to build any model whose
`$CAPTURE` repeats a compartment name:

```
Error in validObject(.Object) :
  invalid class "mrgmod" object: compartment should not be in $CAPTURE: CRP,INFLAM,VDH,CART,STAT3_P,FLS_ACT
```

Confirmed on the *original*, untouched file, independently on two mrgsolve
2.0.1 installs (local R and the qspserver container) — this is not
something the refactor introduced. Since compartments are always present in
mrgsolve's output regardless of `$CAPTURE`, removing the duplicate names
from `$CAPTURE` changes nothing about what is reported; it is required
purely so the model builds at all under the currently installed mrgsolve.
Applied identically to both the verification copy of the original and to
`ra_mrgsolve_model_refactored.R` (never to `ra_mrgsolve_model.R` itself).
Logged as a new entry in `translations/UPSTREAM_ISSUES.md`.

## Verification

Per the guide's mandatory protocol: ran the original file's own dosing
scenarios through both the (`$CAPTURE`-patched, for buildability only)
original and the refactored model via the qspserver `mrgsolve_api` service
(`POST /model_manifest`, `POST /run_simulation`), 1-year horizon, daily
(`delta=24`) output grid — matching the original scenarios' own settings.

Scenarios run (dose amounts/timing copied exactly from the original file's
own `ev()` calls, translated to the API's `ii`/`addl` dosing convention):

1. **Scenario 3** — TCZ IV 8 mg/kg q4w (3784 nmol into `CENT_TCZ`, q672h ×14)
2. **Scenario 4** — TCZ SC 162 mg q2w (1095 nmol into `GUT_TCZ`/`SC_TCZ`, q336h ×27)
3. **Scenario 5** — TCZ IV + MTX combo (as above, plus 20 mg into
   `GUT_MTX` q168h ×53) — exercises the shared `DRUG_ANTI_INFLAM` expression
   with a second active compound.

Every `$CAPTURE`d output (`DAS28_est`, `HAQ_DI_est`, `RO_pct`, `C_TCZ`,
`C_ADA`, `C_MTX`, `C_BARI`, `ACR20/50/70`, `REMISSION`, `TNF_INH_ADA`,
`EFFECT_TCZ`/`OCC_TCZ`) plus the raw compartment states
(`CENT_TCZ`, `REC_FREE_TCZ`/`IL6R_FREE`, `COMPLEX_TCZ`/`DRUG_RL`, `CRP`,
`INFLAM`, `VDH`, `CART`, `STAT3_P`, `FLS_ACT`) was compared point-by-point
across the full 367-point time grid for all three scenarios.

**Result: exact match, max abs diff = 0.0 for every output, every
scenario.** This is the expected outcome for a rename (not a fit) per the
guide's own tolerance rule — since `EFFECT_TCZ` is computed from the same
ODE state ratio as `OCC_TCZ` (renamed variables, identical arithmetic,
identical execution order), there is no approximation to introduce error,
so the near-exact match the guide anticipates for a true rename came out as
a literal exact match rather than merely "close."

## Anything else worth flagging

- `RTOT0`/`RTOT_TCZ` is declared as a parameter (documented as "baseline
  membrane IL-6R") but is **not referenced anywhere in the original's
  `$MAIN`/`$ODE`** — the actual receptor pool is whatever
  `IL6R_FREE + DRUG_RL` dynamically works out to via `KSYN_R`/`KDEG_R`
  turnover (which happens to equal `RTOT0` at steady state, see above).
  Preserved as-is (unused) in the refactor; not fixed, not removed, per the
  "don't invent, don't drop the original's parameters" rule.
- All derived PD doubles (`C_TCZ`, `C_ADA`, `C_MTX`, `C_BARI`, `EFFECT_TCZ`,
  `JAK_INH`, `IL6_sig`, `TNF_INH_ADA`, `DRUG_ANTI_INFLAM`, `MTX_EFF_INFLAM`)
  are computed in `$MAIN`, which in mrgsolve is evaluated once per requested
  output/dosing interval using the *start-of-interval* state — this causes
  a one-interval reporting lag relative to the true continuous state,
  clearly visible with the model's own `delta=24` (daily) scenarios (e.g.
  `OCC_TCZ` reads `0` in the very first row after an IV bolus, even though
  the receptor is already ~96% occupied by then; it "catches up" one row
  later). This affects every compound identically and is not specific to
  TCZ, so it was left untouched rather than moved to `$ODE`/`$TABLE` — doing
  so would be a behavioral fix, not a rename, and would break the "exact
  match" verification above. It is exactly why the exact-match result above
  is meaningful evidence of a correct rename rather than a coincidence: the
  refactor reproduces the lag faithfully rather than accidentally fixing it.
  Logged in `translations/UPSTREAM_ISSUES.md`.
- Compartment ordering is unchanged (`GUT_TCZ` is still compartment 1, etc.)
  since renaming was done in place without reordering — so any external
  code addressing compartments by 1-based number, not name, is unaffected.

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`,
`rheumatoid-arthritis | Tocilizumab (TCZ)` row.
