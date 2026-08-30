# Refactor notes — `pmr_mrgsolve_model_refactored.R`

Scope: **tocilizumab (TCZ) only** — its PK block and its downstream disease
effect equations. Prednisolone and every disease equation not driven by TCZ
concentration are byte-for-byte copies of `pmr_mrgsolve_model.R`.

## Archetype determination

Read the actual `dxdt_` lines rather than assuming from the comment
("Add simple TMDD"). The original TCZ block is:

- 3 compartments: `DEPOT_TCZ` (SC depot), `CENT_TCZ`, `PERI_TCZ`.
- Linear absorption (`KA_TCZ`, `F_TCZ`) into a two-compartment linear
  distribution model, parameterized with **both** styles at once: `CL_TCZ`/
  `V1_TCZ`/`V2_TCZ`/`Q_TCZ` declared in `$PARAM`, but the `$ODE` block
  re-derives micro-constants from them (`k12_T = Q_TCZ/V1_TCZ`,
  `k21_T = Q_TCZ/V2_TCZ`, `kel_T = CL_TCZ/V1_TCZ`) and writes the equations
  in terms of those micro-constants.
- No `REC_FREE_TCZ`/`COMPLEX_TCZ` compartments anywhere, despite `KON`,
  `KOFF`, and a comment claiming "Add simple TMDD". Grepped the whole file
  for `KON`/`KOFF`: they appear only in the `$PARAM` block and one comment
  (`Kd = KOFF/KON ~ 1 nM`) — **never in an ODE**. They are dead parameters.
  What the comment actually implemented is a single extra term,
  `kel_TMDD = KINT * SIL6R / (SIL6R_BASE + SIL6R)`, added to `CENT_TCZ`'s
  elimination — an *extra clearance rate that is a function of the disease
  state (`SIL6R`)*, not a receptor-binding ODE system.

**Conclusion: Archetype 3 (depot + central + peripheral, linear elimination)**,
expressed in the original as derived micro-constants rather than direct
CL/Q/V — exactly the case the design guide calls out ("Covers both CL/Q/V-
style and micro-constant... originals — always rewrite to CL/Q/V"). This is
**not** Archetype 4: forcing in `REC_FREE_TCZ`/`COMPLEX_TCZ` compartments
that the original never instantiates would be adding structure the original
doesn't have, which the guide explicitly prohibits. The `kel_TMDD` state-
coupling term is real pharmacology that the original models (TCZ's own
clearance accelerates as sIL-6R depletion signals it is complex-bound) and
is **kept, verbatim**, in the refactor — it just doesn't fit any of the five
archetypes' patterns and doesn't need to; it rides along in the canonical
CL/Q/V equation as an added elimination term.

## Parameterization conversion applied

Almost all TCZ parameter names already matched the naming convention
(`KA_TCZ`, `F_TCZ`, `CL_TCZ`, `V1_TCZ`, `V2_TCZ`, `Q_TCZ`) — no renaming
needed there. Changes made:

- `DEPOT_TCZ` → `GUT_TCZ` (naming convention; **same compartment position**,
  4th declared, so the R script's `ev(cmt = 4, ...)` dosing is unaffected —
  verified below).
- `KINT` → `KINT_TCZ`, `KON` → `KON_TCZ`, `KOFF` → `KOFF_TCZ` (adding the
  missing stem suffix per convention; `KON_TCZ`/`KOFF_TCZ` remain unused,
  same as in the original — flagged, not removed, since pruning dead
  parameters wasn't part of the assigned scope).
- ODE rewritten from the micro-constant decomposition
  `-(kel_T + kel_TMDD + k12_T) * CENT_TCZ + k21_T * PERI_TCZ` (with
  `kel_T = CL_TCZ/V1_TCZ`, `k12_T = Q_TCZ/V1_TCZ`, `k21_T = Q_TCZ/V2_TCZ`)
  to the canonical direct form:
  ```
  dxdt_CENT_TCZ =  KA_TCZ * F_TCZ * GUT_TCZ
                   - (CL_TCZ + Q_TCZ) / V1_TCZ * CENT_TCZ
                   + Q_TCZ / V2_TCZ * PERI_TCZ
                   - kel_TMDD * CENT_TCZ;
  dxdt_PERI_TCZ =  Q_TCZ / V1_TCZ * CENT_TCZ - Q_TCZ / V2_TCZ * PERI_TCZ;
  ```
  Algebraic identity: `kel_T + k12_T = CL_TCZ/V1_TCZ + Q_TCZ/V1_TCZ =
  (CL_TCZ+Q_TCZ)/V1_TCZ`, exactly as the design guide's own Archetype 2/3
  examples write it. The two forms differ only in floating-point
  associativity (one division of a sum vs. a sum of two divisions), which is
  the source of the ~1e-7-scale residual reported below.
- `C_TCZ` is the single exposed concentration (nM):
  `C_TCZ = CENT_TCZ / V1_TCZ / 148.0 * 1e6` — the exact same expression and
  operation order the original used to compute `Cp_TCZnM` (via the
  intermediate `Cp_TCZ_h`), just given the canonical name and computed once.
  The redundant intermediate (`Cp_TCZ_h`) was dropped as dead after inlining
  (it was never read except to compute `Cp_TCZnM`).
- The `$MAIN` block's own `Cp_TCZ`/`Cp_TCZ_nM` (a second, independently
  redundant computation — mrgsolve's `$MAIN` and `$ODE` are separate
  functions, so `$MAIN` locals are never visible in `$ODE`, and nothing in
  `$TABLE` reads them either) was **left untouched**, out of scope for this
  task and provably inert.
- `$TABLE`/`$CAPTURE` are untouched — they already computed `Cp_TCZ_out`/
  `Cp_TCZnM_out` independently from `CENT_TCZ`/`V1_TCZ`, so nothing there
  needed to change.

## Hill interface

The original's TCZ→disease interface is already a plain Emax ratio in two
places, so per the guide this is **a rename, not a refit**:

1. **IL-6 signalling blockade** (drives both `inh_IL6_TCZ` and the
   paradoxical sIL-6R rise): original `TCZ_occ_sR = Cp_TCZnM/(0.5+Cp_TCZnM)`,
   `inh_IL6_TCZ = 0.85 * TCZ_occ_sR`. Refactored to:
   ```
   double OCC_TCZ    = (C_TCZ > 0) ? pow(C_TCZ, GAMMA_TCZ) / (pow(EC50_TCZ, GAMMA_TCZ) + pow(C_TCZ, GAMMA_TCZ)) : 0.0;
   double EFFECT_TCZ = EMAX_TCZ * OCC_TCZ;
   double inh_IL6_TCZ = EFFECT_TCZ;
   ```
   with `EC50_TCZ = 0.5` (nM), `EMAX_TCZ = 0.85`, `GAMMA_TCZ = 1` pulled out
   as explicit named parameters (the original had these as bare literals).
   `sIL6R_stim = 1.0 + 0.4 * OCC_TCZ` keeps reading the un-scaled occupancy
   `OCC_TCZ`, exactly as the original's `TCZ_occ_sR` did — this is a second,
   smaller pharmacological consequence of the same receptor occupancy
   (paradoxical sIL-6R rise), not part of the "signalling blockade" Emax
   itself, so it was not folded into `EFFECT_TCZ`.
2. **Direct empirical PMR-AS effect** (a separate dose-response curve, EC50
   50 nM vs. 0.5 nM above — clearly a different, more macroscopic
   relationship, not reducible to the same curve): original
   `eff_TCZ_PMRAS = (EMAX_TCZ_PMRAS * Cp_TCZnM) / (EC50_TCZ_PMRAS + Cp_TCZnM)`.
   Renamed `EFFECT_TCZ_PMRAS`, reading `C_TCZ` and an added explicit
   `GAMMA_TCZ_PMRAS = 1` (the original had no Hill coefficient here either).
   `combined_eff` (which combines TCZ's and prednisolone's separate PMR-AS
   effects at the point of use) now reads `EFFECT_TCZ_PMRAS` in place of
   `eff_TCZ_PMRAS`; the prednisolone side (`eff_GC_PMRAS`) is untouched.

No refit was performed anywhere — both Hill terms reproduce the original
formulas exactly (`GAMMA_TCZ = GAMMA_TCZ_PMRAS = 1` makes `pow(x,1) = x`).

## Verification

**Local `mrgsolve` (2.0.1, `C:\Program Files\R\R-4.6.0\bin\Rscript.exe`)
could not be used** — the original file's `$MAIN` block uses the
`_init_<CMT> = value;` idiom to set initial conditions, which this
mrgsolve build's generated C++ no longer accepts (`error: '_init_CORT' was
not declared in this scope`, etc.). Confirmed with a minimal 1-compartment
reproduction outside this file: `_init_X = A0;` fails the same way, while
`X_0 = A0;` (the modern idiom) compiles and runs. **This is a pre-existing
mrgsolve-version incompatibility in the original file, unrelated to
tocilizumab** — it affects the initials for `CORT`, `IL6`, `SIL6R`, `CRP`,
`ESR`, `BMD`, `PMRAS`, `FLARE`, none of which this task touched. Flagging
here for manual review; **not fixed in either `pmr_mrgsolve_model.R` (never
touched) or `pmr_mrgsolve_model_refactored.R`** (out of the assigned scope,
and the guide's "log what you find, don't fix it upstream" rule extends
naturally to not smuggling unrelated fixes into a scoped refactor's sibling
file either).

Verification instead ran against the **qspserver `mrgsolve_api` container**
(`http://localhost:8007`, confirmed healthy), which compiles server-side on
a Debian/GCC toolchain and does accept the `_init_` idiom directly (its
`/model_manifest` call on the **unpatched, byte-identical extracted DSL**
compiled and ran the original with no changes needed) — no version patch
was applied to the actual comparison; both extracted DSL blocks (original
and refactored `code <- '...'` string, pulled unmodified from the two `.R`
files) were POSTed to `/model_manifest` (both returned 200, confirming both
compile and that `C_TCZ`, `EC50_TCZ`, `EMAX_TCZ`, `GAMMA_TCZ` etc. are
discoverable parameters in the refactored manifest) and then to
`/run_simulation`.

Ran **4 of the original file's own scenarios** (identical dosing built by
reproducing `build_events()`'s exact schedule as `DoseSpec` records, sorted
by time, `cmt=1` for prednisolone / `cmt=4` for TCZ depot — same compartment
position in both models):

| Scenario | Description | Result |
|---|---|---|
| S1_no_treatment | no dosing at all | exact match, 0 deviation on every output |
| S5_TCZ_QW | TCZ 162 mg SC QW + pred 12.5 mg taper | TCZ-driven outputs match to ~1e-7 (see below) |
| S6_TCZ_Q2W | TCZ 162 mg SC Q2W + pred 12.5 mg taper | TCZ-driven outputs match to ~2e-7 |
| S7_TCZ_steroid_free | TCZ 162 mg SC QW, no prednisolone | max relative deviation 2.5e-8 |

For every scenario, every `$CAPTURE`d output was compared across the full
6-hourly time grid (2921 points, t = 0…17520 h). Results, S6 shown in full
(representative):

```
CRP_out        max_rel=0.000e+00
ESR_out        max_rel=0.000e+00
IL6_out        max_rel=0.000e+00
PMRAS_out      max_rel=0.000e+00
CORT_out       max_rel=0.000e+00
BMD_out        max_rel=0.000e+00
Cp_TCZ_out     max_rel=0.000e+00
Cp_TCZnM_out   max_rel=2.236e-07
FLARE_out      max_rel=0.000e+00
SIL6R_out      max_rel=0.000e+00
CRP_norm       max_rel=0.000e+00
```

**Result: verification passed, near-exact as expected for a pure
reparameterization with no Hill-fitting.** The ~2e-7 residual on
`Cp_TCZnM_out` is exactly the floating-point-associativity difference
predicted above (`(CL+Q)/V1` vs `CL/V1 + Q/V1`), well inside "floating-
point-scale" and not a bug.

**One artifact worth flagging, investigated and confirmed benign:** in
scenarios with active prednisolone tapering to 0 (S5, S6), `Cp_PRED_out`
and `Cp_FREE_out` show large *relative* deviations (up to ~14x) at isolated
time points late in the simulation. Inspected the actual values: at every
such point both the original's and the refactored's value are in the
1e-7–1e-9 absolute range (prednisolone has fully washed out by then; peak
concentration earlier in the run is 4.6 μg/mL), i.e. both are numerically
zero and differ only by ODE-solver noise at the ~1e-7 absolute scale. Since
the prednisolone code is byte-for-byte identical between the two files, this
is solver/compilation noise unrelated to the TCZ refactor, not a real
discrepancy — confirmed by absolute-difference inspection, not just
asserted.

## Other findings flagged for manual review

1. **`KON`/`KOFF` (renamed `KON_TCZ`/`KOFF_TCZ`) are dead parameters** in
   the original — declared, documented with a `Kd = KOFF/KON` comment, and
   never used in any ODE. The actual "simple TMDD" is the single
   `kel_TMDD` term described above. Worth a manual decision on whether a
   future pass should either wire up real receptor kinetics (Archetype 4)
   or drop the unused parameters — out of scope here.
2. **mrgsolve-version incompatibility**: the original's `_init_<CMT>`
   initial-condition idiom does not compile under mrgsolve 2.0.1 (local
   Windows R install). It does compile under the qspserver API's mrgsolve
   build. Not fixed anywhere per this task's scope; flagged for whoever next
   touches this file locally.
3. The file's own header comment claims "22 ODEs (Prednisolone PK x5 +
   Tocilizumab PK x4 + ...)"; the actual compartment count is 3 (Prednisolone)
   + 3 (Tocilizumab) + 8 (disease) = 14, not 22. Pre-existing, unrelated to
   TCZ's PK/effect equations, left as-is.
