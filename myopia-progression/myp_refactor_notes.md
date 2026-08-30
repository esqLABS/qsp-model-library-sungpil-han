# Refactor notes — `myopia-progression/myp_mrgsolve_model.R`

Scope: **7-methylxanthine ("MX")** only — the compound listed as "Other
treatments (MX)" in `driver-patches/data/compound_perturbation_census.md`.
Only MX's PK block (`MXGUT`/`MXPLA` compartments, `KAMX`/`KELMX`/`VMX`
parameters) and its single point of entry into the disease drive (`EMX`)
were touched. Atropine (11-compartment PK + three-site PD), orthokeratology,
red-light photobiomodulation, and every optical-defocus treatment are
byte-for-byte identical to the original in `myp_mrgsolve_model_refactored.R`
— confirmed by diff (see below).

MX itself is **7-methylxanthine**, an oral xanthine derivative dosed at
400 mg/day in the Danish pilot series (Trier et al. 2008, `myp_references.md`
#62; Cui et al. 2011, #60) and modeled here as a plain oral one-compartment
PK feeding a fixed fractional suppression of the constitutive myopigenic
growth drive — this is a genuinely "clean single site" compound (one PK
block, one downstream read), matching the census row's pre-assigned
"Redirect concentration (clean single site)" classification.

## Archetype determined

**A depot + central variant of Archetype 3, with no peripheral
compartment** (the guide explicitly allows dropping the peripheral
compartment from Archetype 3's template; this drops it and keeps the
depot). The original:

```
dxdt_MXGUT = -KAMX * MXGUT;
dxdt_MXPLA = KAMX * MXGUT / VMX - KELMX * MXPLA;
```

`MXGUT` is an amount compartment (umol); `MXPLA` is unusual in that it is
written directly as a **concentration state** (umol/L) rather than an
amount — the volume division is baked into the absorption term
(`KAMX*MXGUT/VMX`) and elimination is a plain first-order term on the
concentration itself (`-KELMX*MXPLA`). This is mathematically the
change-of-variables `MXPLA ≡ Amt_central / VMX` with `KELMX` playing the
role of a micro-constant `k10 = CL/V`: multiplying the whole `MXPLA` ODE by
`VMX` gives `d(VMX*MXPLA)/dt = KAMX*MXGUT - KELMX*(VMX*MXPLA)`, i.e. a
completely standard amount-based one-compartment-with-depot model with
`CL = KELMX*VMX`. No compartment was added or removed; this is a rename to
the guide's convention (amount-based `CENT_MX` + explicit `CL_MX`/`V1_MX`
+ derived `double C_MX`) rather than a change to the underlying dynamics —
verified below to reproduce the original's concentration trajectory
exactly.

There is no bioavailability parameter in the original (dosing enters
`MXGUT` directly, implicitly `F=1`) — no `F_MX` was invented, per the
guide's "never invent a PK parameter the original didn't have" rule.

## Renaming applied (MX only)

| Original | Refactored | Note |
|---|---|---|
| `MXGUT` (compartment) | `GUT_MX` | amount (umol), unchanged semantics |
| `MXPLA` (compartment) | `CENT_MX` | **now an amount (umol)**, not a concentration — see derivation above; `C_MX = CENT_MX/V1_MX` recovers the original's concentration exactly (verified) |
| `KAMX` | `KA_MX` | value unchanged (6.0 /day) |
| `KELMX`, `VMX` | `CL_MX` = `KELMX*VMX` = 96.0 (L/day), `V1_MX` = `VMX` = 40 (L) | `CL_MX` is a derived value from the original's own two numbers, not invented |
| `EMAXMX` | `EMAX_MX` | value unchanged (0.38) |
| `EC50MX` | `EC50_MX` | value unchanged (1.2 umol/L) |
| — (none; original had no Hill exponent) | `GAMMA_MX` (new, `= 1`, "original had no explicit Hill term") |
| `EMX` (local double in `$ODE`) | `EFFECT_MX` | |
| `MXDOSE` | unchanged | see "anything else flagged" below — not part of the naming convention's role table, and not read anywhere in `$MAIN`/`$ODE` |

`dxdt_MXGUT`/`dxdt_MXPLA` → `dxdt_GUT_MX`/`dxdt_CENT_MX`:

```
dxdt_GUT_MX  = -KA_MX * GUT_MX;
dxdt_CENT_MX =  KA_MX * GUT_MX - CL_MX / V1_MX * CENT_MX;
```

## Hill interface: rename, not a fit

The original's effect term is already the plain Emax/EC50 ratio shape:

```
double EMX = EMAXMX * MXPLA / (EC50MX + MXPLA);
```

Per the guide, this is a rename, not a refit. `GAMMA_MX = 1` reproduces the
ratio exactly for any concentration:

```
double C_MX = CENT_MX / V1_MX;                       // umol/L
double EFFECT_MX = EMAX_MX * pow(C_MX, GAMMA_MX)
                  / (pow(EC50_MX, GAMMA_MX) + pow(C_MX, GAMMA_MX));
```

`EFFECT_MX` is used in exactly the one place `EMX` was used, the
constitutive-drive suppression term:

```
double SIGEFF = SIG * FDA * (1.0 - EATR) * (1.0 - EFFECT_MX) * (1.0 - ERL) * MRUP;
```

Only the `EMX`→`EFFECT_MX` token was touched on this line; the other three
multiplicative factors (`EATR` for atropine, `ERL` for red light, and `MRUP`
for receptor up-regulation) belong to other compounds/mechanisms and are
untouched.

## Diff scope confirmation

`diff -u myp_mrgsolve_model.R myp_mrgsolve_model_refactored.R` touches
exactly 6 hunks: the `KAMX`/`KELMX`/`VMX`/`EMAXMX`/`EC50MX` param lines
(+1 new `GAMMA_MX` line), the `MXGUT`/`MXPLA` compartment declaration
lines, the `EMX`/`C_MX`/`EFFECT_MX` derivation lines in `$ODE`, the
`SIGEFF` line that reads the effect term, the `dxdt_MXGUT`/`dxdt_MXPLA`
lines, and two additive-only documentation changes (two new `$CAPTURE`
lines for qspserver discoverability, and a one-line note on the scenario-28
comment block about the compartment rename). Nothing else in the 797-line
file differs — every other compound's `$PARAM`/`$CMT`/`$ODE`/`$TABLE`/
`$CAPTURE` lines, the optics/defocus/effector-cascade/scleral-creep/choroid/
lens/risk-hazard machinery, and the full 28-scenario list are byte-identical.

## A pre-existing upstream defect found while verifying (not fixed, logged)

The original file does not build under mrgsolve 2.0.1 as written — the
`TRTPD` parameter's annotated description is split across two lines, and
the continuation line has no `name : value` pair, which the annotated-block
parser rejects outright (`Error: improper annotation format`). This
reproduces identically from the untouched original (confirmed via the
qspserver `mrgsolve_api` container) and is unrelated to MX. Logged as
upstream issue **#32** in `translations/UPSTREAM_ISSUES.md`.

Verification below used in-memory-only scratch copies of both the original
and the refactored DSL with the two lines mechanically merged into one
(moving the per-device value list into the existing description's
parentheses — no name, value, or number touched). This identical,
non-numeric workaround was applied to both sides of every comparison.
**Neither the tracked original `myp_mrgsolve_model.R` nor the delivered
`myp_mrgsolve_model_refactored.R` was changed** — both still contain the
original's two-line `TRTPD` annotation exactly as written, and neither
currently builds against mrgsolve 2.0.1 without that same workaround.

## Verification

**Method.** Both the original's and the refactored file's scratch (defect-
patched) copies were POSTed as `model_content` to the qspserver
`mrgsolve_api` service (`http://localhost:8007`, mrgsolve 2.0.1):
`POST /model_manifest` first (confirmed both compile and that every
`KA_MX`/`CL_MX`/`V1_MX`/`EMAX_MX`/`EC50_MX`/`GAMMA_MX` parameter and every
`C_MX`/`EFFECT_MX` output are discoverable in the refactored manifest — see
qspserver compatibility section below), then `POST /run_simulation` for two
scenarios built from the original file's own scenario list, comparing every
shared `$CAPTURE`d/compartment output across the full time grid:

1. **Scenario 28 — 7-methylxanthine 400 mg/day** (`param(MXDOSE=400)` + "a
   daily gut event", per the original's own scenario comment). Reconstructed
   as a daily bolus of `400/194.19*1000 = 2059.838 umol` (MW 194.19 g/mol,
   from `myp_reference_model.py`'s own dosing code) into `MXGUT`/`GUT_MX`
   (compartment index 43 in both files), `ii=1`, `addl=729` (2 years daily),
   `end=730`, `delta=7` — all other parameters at their defaults (no
   atropine, no device).
2. **Untreated baseline** (all defaults, `MXDOSE=0`, no dosing events),
   `end=3650`, `delta=30` — confirms MX-inactive behavior is unaffected.

**Result: exact match.**

- Every output whose underlying computation is untouched by this refactor
  (`ALS`, `ACD`, `LT`, `PLENS`, `CRAD`, `RPRS`, `CHT`, `CBF`, `CHOX`,
  `AGEy`, `SER`, `SERTR`, `AL`, `ALCR`, `EATRo`, `MTRo`, `VIRISK`, `PMMD`,
  `PRD`, `PCNV`, `PGLC`, `PCAT`, `AXCUM`, `SERAUC`) matched with
  **max abs dev = 0.0** in both scenarios.
- `MXGUT` vs `GUT_MX`: **max abs dev = 0.0** (this compartment's ODE and
  units are unchanged).
- `CENT_MX/V1_MX` (the refactored file's derived `C_MX`) vs the original's
  `MXPLA` (which *is* the concentration): **max abs dev = 0.0** in the
  treated scenario, `4.2e-21` (floating-point noise) in the baseline —
  confirms the amount-vs-concentration change of variables is exact, not
  approximate.
- `EFFECT_MX` (refactored, `$CAPTURE`d) vs `EMAXMX*MXPLA/(EC50MX+MXPLA)`
  (recomputed from the original's captured `MXPLA` using the original's own
  constants): differs by up to `4.88e-05` in raw JSON — traced to the
  qspserver API's 4-decimal-place output rounding (`0.33224880883815905`
  reported as `0.3322`; `0.33224880883815905 - 0.3322 = 4.8809e-05`, i.e.
  exactly the rounding residual, not a computational difference). Recomputing
  `EFFECT_MX` from the refactored file's own reported `C_MX` reproduces the
  full-precision value exactly. In the untreated baseline both `C_MX` and
  `EFFECT_MX` are identically `0.0` throughout, as expected.

This is the expected outcome for a rename (not a fit) per the guide's
tolerance rule for Archetypes 1–3: pure structural reorganization, same
math, reorganized — and the one legitimate unit change (`MXPLA`
concentration → `CENT_MX` amount) is proven to be an exact, invertible
change of variables rather than an approximation.

## qspserver compatibility requirements

- **No separate `.cpp` extraction needed.** `myp_mrgsolve_model.R` is
  already written directly in mrgsolve block-marker syntax
  (`$PROB`/`$PARAM`/.../`$CAPTURE`, consumed via `mread()`), not the
  `code <- '...'; mcode(code)` R-string-wrapper pattern the guide's `.cpp`-
  extraction step is aimed at (same situation as `thyroid-eye-disease`).
  Confirmed empirically: the whole file's content, posted as-is (after the
  in-memory `TRTPD` workaround above) to `/model_manifest` and
  `/run_simulation`, compiled and ran correctly — see Verification.
- **`KA_MX`, `CL_MX`, `V1_MX`, `EMAX_MX`, `EC50_MX`, `GAMMA_MX`** are all
  declared in `$PARAM`, so `/model_manifest` lists them with their original
  values (`6`, `96`, `40`, `0.38`, `1.2`, `1` respectively) — confirmed in
  the manifest response.
- **`C_MX` and `EFFECT_MX`** are computed in `$ODE` from compartment state,
  so (as with the TED/RA precedents) they cannot themselves be `$PARAM`
  entries. Both were added to `$CAPTURE` — additive-only, not present in the
  original's `$CAPTURE` — so both are discoverable as **outputs** via
  `/model_manifest`'s `outputPaths` and retrievable via `/run_simulation`'s
  `outputs` list. Confirmed both ways (see manifest excerpt in Verification).
- `$SET`/scenario dosing is untouched and was not relied upon for
  verification — both `/run_simulation` calls above supplied `time`/`dosing`
  directly in the request, matching the original's own scenario 28 dosing
  pattern reconstructed from its comment and `myp_reference_model.py`.

## Anything else worth flagging

- `MXDOSE` (7-MX oral dose, mg/day) is declared in `$PARAM` but **is not
  referenced anywhere in `$MAIN`/`$ODE`** — the original's scenario 28
  achieves MX dosing entirely through an external daily bolus event into
  `MXGUT`, not through this parameter. `MXDOSE` is documentation only in
  the original (its value is never read by the model). Left completely
  untouched in the refactored file per the "don't invent, don't drop the
  original's parameters" rule — this is not a role in the naming
  convention's table (it isn't a `KA`/`CL`/`V1`/`EMAX`/`EC50` etc.), so it
  was not renamed either.
- The scenario-28 comment at the bottom of the file doses into `MXGUT` by
  name; since this refactor renames that compartment to `GUT_MX`, a
  one-line note was added next to that comment in the refactored file
  warning that a dosing event targeting the refactored model needs to use
  `GUT_MX` instead — the comment block itself (all 28 scenarios) is
  otherwise copied verbatim, matching the precedent set for `CEN_TCZ` in
  `thyroid-eye-disease/ted_refactor_notes.md`.
- Compartment ordering is unchanged (`GUT_MX` is still compartment 43,
  `CENT_MX` still compartment 44) since renaming was done in place without
  reordering — confirmed via `/model_manifest`'s `outputPaths`.

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`,
`myopia-progression | Other treatments (MX)` row: target compartment
`CENT_MX` / exposed concentration `C_MX (umol/L)`, status noting the
archetype, the rename-not-fit Hill result, and the exact-match verification.
