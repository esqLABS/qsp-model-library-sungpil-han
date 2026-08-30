# Refactor notes — `ppgl_mrgsolve_model.R` (phenoxybenzamine, doxazosin, metyrosine, beta-blocker, sunitinib)

**Scope of this pass.** Per the fork's PK/PD refactor spec
([`FORK_WORKFLOW_GUIDE.md`](../FORK_WORKFLOW_GUIDE.md), Part 2), all 5
compounds this 20-compartment file models were rewritten, per their
existing rows in
[`driver-patches/data/compound_perturbation_census.md`](../driver-patches/data/compound_perturbation_census.md)
(all classified "Redirect concentration (clean single site)"): **BB**
(a beta-blocker, propranolol-class, stem `BB`), **DOX** (doxazosin, stem
`DOX`), **MET** (metyrosine, stem `MET`), **PHE** (phenoxybenzamine, stem
`PHE`), and **SUNIT** (sunitinib, stem `SUNIT`). This file models exactly
these 5 compounds and nothing else — there is no 6th untouched compound to
carve out.

## A structural observation shared by all 5 compounds: the states are already concentrations, not amounts

Every one of this file's PK compartments (`PHE_C`/`PHE_P`, `DOX_C`,
`MET_C`, `BB_C`, `SUNIT_C`/`SUNIT_P`) is written so that **every** term in
its own `dxdt_` equation is already pre-divided by the relevant volume —
e.g. `dxdt_PHE_C = ka_PHE*F_PHE*PHE_gut/V1_PHE - (CL_PHE/V1_PHE)*PHE_C - ...`.
This means the state variable itself tracks concentration directly (a
"unit-volume-reduced" ODE), not amount-that-gets-divided-by-V downstream,
which is the more common convention (and the one the guide's own
archetype templates show, `double C_TCZ = CENT_TCZ / V1_TCZ;`). Both
formulations are mathematically valid and this one is internally
consistent throughout the file (dosing bolus amounts are likewise
pre-scaled to land directly in concentration units, see below), so it was
preserved exactly rather than "corrected" to the more common form — doing
the latter would require re-deriving every dosing amount and would change
nothing this refactor is in scope to change. Consequence for naming: each
compound's `C_<STEM>` is a **direct alias** of its renamed central
compartment (`double C_PHE = CENT_PHE;`), not `CENT_PHE / V1_PHE`.

## A dosing-mechanism observation shared by DOX/MET/BB/SUNIT

The original declares a `DOSE_<STEM>` flag parameter (default `0`) and an
inline `KA_<stem>*F_<stem>*DOSE_<stem>/V_<stem>` absorption-style term
inside each of these four compounds' `dxdt_` lines — but `DOSE_DOX`,
`DOSE_MET`, `DOSE_BB`, and `DOSE_SUNIT` are **never set to anything but
their `0` default anywhere in the original's own `run_scenario()`**
(confirmed by grep: the only assignment is the initial
`params <- c(DOSE_PHE=0, DOSE_DOX=0, DOSE_MET=0, DOSE_BB=0, DOSE_SUNIT=0, SURGERY=0)`).
The actual dosing mechanism for all four is instead a plain `ev()` bolus
dosed directly into the compound's own central compartment
(`cmt = "DOX_C"`, `"MET_C"`, `"BB_C"`, `"SUNIT_C"`), which works because
that compartment is concentration-valued (see above) and mrgsolve adds a
bolus `amt` directly to the compartment's current value. `PHE` is the
only one of the 5 with a real depot/absorption pathway that is actually
used (`PHE_gut`, dosed by its own `ev()`).

**Handling:** the vestigial `KA_<STEM>`/`F_<STEM>`/`DOSE_<STEM>` terms
were kept, renamed, exactly as computed in the original (they always
evaluate to zero either way, so keeping vs. dropping them is numerically
identical) — preserved for maximum fidelity to the original rather than
silently deleting dead code, and called out explicitly here and in the
`$PARAM` comments of the refactored file so a reviewer isn't left
wondering why an unused absorption pathway exists. `kout_PHE`
(phenoxybenzamine) and `IC50_alpha` (cardiovascular section) are two
further genuinely dead parameters in the original — declared, never
referenced anywhere in `$MAIN`/`$ODE`/`$TABLE` — also preserved as-is,
untouched and un-renamed for `IC50_alpha` since it isn't part of any of
the 5 compounds' own `$PARAM` block.

## Archetypes determined

### PHE (phenoxybenzamine) — archetype 3 (depot + central + peripheral)

```
dxdt_GUT_PHE  = -KA_PHE * GUT_PHE;
dxdt_CENT_PHE =  KA_PHE * F_PHE * GUT_PHE / V1_PHE
               - (CL_PHE/V1_PHE) * CENT_PHE
               - (Q_PHE/V1_PHE)  * CENT_PHE
               + (Q_PHE/V2_PHE)  * PERI_PHE;
dxdt_PERI_PHE =  (Q_PHE/V1_PHE)  * CENT_PHE
               - (Q_PHE/V2_PHE)  * PERI_PHE;
```

Exactly the guide's archetype-3 shape (the original's two separate
`-(CL/V1)*C - (Q/V1)*C` terms are algebraically identical to the guide's
combined `-(CL+Q)/V1*C`; left as two terms since that's how the original
wrote it, no behavioral difference). Irreversible α-blockade, a plain
Hill-shaped ratio in the original:

| Original | Refactored | Value |
|---|---|---|
| `PHE_gut`/`PHE_C`/`PHE_P` (compartments) | `GUT_PHE`/`CENT_PHE`/`PERI_PHE` | — |
| `ka_PHE`/`F_PHE`/`CL_PHE`/`V1_PHE`/`V2_PHE`/`Q_PHE` | `KA_PHE`/`F_PHE`/`CL_PHE`/`V1_PHE`/`V2_PHE`/`Q_PHE` | 0.693 (1/h) / 0.27 / 5.8 (L/h) / 45 (L) / 180 (L) / 10 (L/h) |
| `kout_PHE` (dead, unused) | `KOUT_PHE` (dead, unused) | 0.0035 |
| `PHE_C` (state, already concentration) | `C_PHE` | µmol/L, direct alias of `CENT_PHE` |
| `0.012` (inline literal in `PHE_eff`) | `EC50_PHE` | 0.012 (µmol/L), promoted to a named `$PARAM` |
| — (none; implicit linear ratio, no separate Emax) | `EMAX_PHE` (new, `=1.0`), `GAMMA_PHE` (new, `=1.0`) | rename, not a fit |
| `PHE_eff` (local `$ODE` ratio) | `EFFECT_PHE` | used only at `alpha_block`'s Bliss-independence combination with `EFFECT_DOX` (see below) |

### DOX (doxazosin) — archetype 1 (single compartment, direct-bolus dosing)

```
double CL_DOX_eff = 0.693 / T50_DOX * V1_DOX;
dxdt_CENT_DOX = KA_DOX * F_DOX * DOSE_DOX / V1_DOX
             - (CL_DOX_eff / V1_DOX) * CENT_DOX;
```

No depot — dosed by a plain `ev()` bolus directly into `CENT_DOX` (see
the dosing-mechanism note above). One genuine peculiarity, present in the
original and preserved exactly: the `$PARAM CL_DOX = 3.2` is **dead** —
the original instead derives clearance every step from the terminal
half-life, `CL_DOX_calc = 0.693/t50_DOX*V_DOX = 0.693/22*80 ≈ 2.52 L/h`,
which is what actually appears in `dxdt_DOX_C`. This is not a compile
defect (both declarations coexist without error), just an inconsistency
between a declared-but-unused `$PARAM` and the value the model actually
uses — preserved as-is (renamed `CL_DOX_eff`), not "fixed", since altering
which clearance value drives the ODE would be a behavioral change outside
this refactor's scope.

| Original | Refactored | Value |
|---|---|---|
| `DOX_C` (compartment) | `CENT_DOX` | — |
| `ka_DOX`/`F_DOX` (vestigial) | `KA_DOX`/`F_DOX` (vestigial, preserved) | 0.21 (1/h) / 0.65 |
| `CL_DOX` (dead `$PARAM`) | `CL_DOX` (dead `$PARAM`, preserved) | 3.2 (L/h) — never used by either model |
| `V_DOX` | `V1_DOX` | 80 (L) |
| `t50_DOX` | `T50_DOX` | 22 (h) |
| `CL_DOX_calc` (local `$ODE`) | `CL_DOX_eff` | the actual clearance basis, `0.693/T50_DOX*V1_DOX` |
| `DOX_C` (state, already concentration) | `C_DOX` | µmol/L, direct alias of `CENT_DOX` |
| `0.002` (inline literal in `DOX_eff`) | `EC50_DOX` | 0.002 (µmol/L), promoted to a named `$PARAM` |
| — (none; implicit linear ratio, no separate Emax) | `EMAX_DOX` (new, `=1.0`), `GAMMA_DOX` (new, `=1.0`) | rename, not a fit |
| `DOX_eff` (local `$ODE` ratio) | `EFFECT_DOX` | used only at `alpha_block`'s Bliss-independence combination with `EFFECT_PHE` (see below) |

### MET (metyrosine) — archetype 1 (single compartment, direct-bolus dosing)

```
dxdt_CENT_MET = KA_MET * F_MET * DOSE_MET / V1_MET
             - (CL_MET / V1_MET) * CENT_MET;
```

| Original | Refactored | Value |
|---|---|---|
| `MET_C` (compartment) | `CENT_MET` | — |
| `ka_MET`/`F_MET` (vestigial) | `KA_MET`/`F_MET` (vestigial, preserved) | 0.70 (1/h) / 0.85 |
| `CL_MET`/`V_MET` | `CL_MET`/`V1_MET` | 4.5 (L/h) / 40 (L) |
| `MET_C` (state, already concentration) | `C_MET` | µmol/L, direct alias of `CENT_MET` |
| `IC50_MET` | `EC50_MET` | 85 (µmol/L), same value |
| — (none; the original's `0.80 *` max-inhibition multiplier) | `EMAX_MET` (new, `=0.80`) | absorbs the original's own scale, not a fit |
| — (none; implicit linear ratio) | `GAMMA_MET` (new, `=1.0`) | rename, not a fit |
| `MET_inhib` (local `$ODE` ratio) + the inline `0.80 *` scale at `TH_target` | `EFFECT_MET` | `EMAX_MET*C_MET/(EC50_MET+C_MET)` now equals the original's `0.80*MET_inhib` exactly; `TH_target = TH_base*(1-EFFECT_MET)` |

### BB (beta-blocker) — archetype 1 (single compartment, direct-bolus dosing)

```
dxdt_CENT_BB = KA_BB * F_BB * DOSE_BB / V1_BB
            - (CL_BB / V1_BB) * CENT_BB;
```

| Original | Refactored | Value |
|---|---|---|
| `BB_C` (compartment) | `CENT_BB` | — |
| `ka_BB`/`F_BB` (vestigial) | `KA_BB`/`F_BB` (vestigial, preserved) | 0.80 (1/h) / 0.36 |
| `CL_BB`/`V_BB` | `CL_BB`/`V1_BB` | 50 (L/h) / 260 (L) |
| `BB_C` (state, already concentration) | `C_BB` | µmol/L, direct alias of `CENT_BB` |
| `IC50_BB` | `EC50_BB` | 0.022 (µmol/L), same value |
| — (none; the original's `0.40 *` max-blunting multiplier) | `EMAX_BB` (new, `=0.40`) | absorbs the original's own scale, not a fit |
| — (none; implicit linear ratio) | `GAMMA_BB` (new, `=1.0`) | rename, not a fit |
| `beta_block` (local `$ODE` ratio) + the inline `0.40 *` scale at `HR_target` | `EFFECT_BB` | `EMAX_BB*C_BB/(EC50_BB+C_BB)` now equals the original's `0.40*beta_block` exactly; `HR_target *= (1-EFFECT_BB)` |

### SUNIT (sunitinib) — archetype 2 (2-compartment, no depot, direct-bolus dosing)

```
dxdt_CENT_SUNIT = KA_SUNIT * F_SUNIT * DOSE_SUNIT / V1_SUNIT
               - (CL_SUNIT/V1_SUNIT) * CENT_SUNIT
               - (Q_SUNIT/V1_SUNIT)  * CENT_SUNIT
               + (Q_SUNIT/V2_SUNIT)  * PERI_SUNIT;
dxdt_PERI_SUNIT =  (Q_SUNIT/V1_SUNIT)  * CENT_SUNIT
               - (Q_SUNIT/V2_SUNIT)  * PERI_SUNIT;
```

| Original | Refactored | Value |
|---|---|---|
| `SUNIT_C`/`SUNIT_P` (compartments) | `CENT_SUNIT`/`PERI_SUNIT` | — |
| `ka_SUNIT`/`F_SUNIT` (vestigial) | `KA_SUNIT`/`F_SUNIT` (vestigial, preserved) | 0.28 (1/h) / 0.60 |
| `CL_SUNIT`/`V1_SUNIT`/`V2_SUNIT`/`Q_SUNIT` | `CL_SUNIT`/`V1_SUNIT`/`V2_SUNIT`/`Q_SUNIT` | 34 (L/h) / 2230 (L) / 1900 (L) / 6.5 (L/h) |
| `SUNIT_C` (state, already concentration) | `C_SUNIT` | µmol/L, direct alias of `CENT_SUNIT` |
| `IC50_SUNIT_VEGFR` | `EC50_SUNIT` | 0.018 (µmol/L), same value |
| — (none; the original's `0.65 *` max-inhibition multiplier) | `EMAX_SUNIT` (new, `=0.65`) | absorbs the original's own scale, not a fit |
| — (none; implicit linear ratio) | `GAMMA_SUNIT` (new, `=1.0`) | rename, not a fit |
| `SUNIT_inh` (local `$ODE` ratio) + the inline `0.65 *` scale at `kgrowth_eff` | `EFFECT_SUNIT` | `EMAX_SUNIT*C_SUNIT/(EC50_SUNIT+C_SUNIT)` now equals the original's `0.65*SUNIT_inh` exactly; `kgrowth_eff = kgrowth_h*(1-EFFECT_SUNIT)` |

## Why PHE and DOX are NOT collapsed into one shared Hill, and why `alpha_block_max` stays where it is

The original computes the combined alpha-blockade fraction as a
two-drug **Bliss independence** combination of two *separate*,
already-independent saturating ratios, then applies one more,
disease-level scale on top:

```
PHE_eff      = PHE_C / (PHE_C + 0.012);
DOX_eff      = DOX_C / (DOX_C + 0.002);
alpha_block  = 1 - (1 - PHE_eff) * (1 - DOX_eff);
SBP_target  *= (1 - alpha_block_max * alpha_block);   // alpha_block_max = 0.92
```

This is the textbook shape the guide's "combine only at the point of use"
rule describes for two *different* compounds sharing one pathway — so
each keeps its own named `EFFECT_<STEM>` (`EMAX_PHE=1.0`, `EMAX_DOX=1.0`,
i.e. each is the raw, unscaled saturating ratio) and the combination is
left exactly as the original wrote it, just with the renamed inputs:

```
double alpha_block = 1 - (1 - EFFECT_PHE) * (1 - EFFECT_DOX);
```

**`alpha_block_max` (0.92) could not be folded into either drug's own
`EMAX_<STEM>`** the way `EMAX_MET`/`EMAX_BB`/`EMAX_SUNIT` absorb their own
compound-specific scaling factors above. Baking it into, say,
`EMAX_PHE = 0.92` would change `alpha_block`'s own algebra
(`1-(1-0.92·PHE_eff)(1-DOX_eff)` ≠ `0.92·(1-(1-PHE_eff)(1-DOX_eff))` in
general — the two only agree when exactly one drug is present at a time),
so it remains a shared, disease-level constant applied once, after the
Bliss combination, exactly where and how the original applied it. This is
disclosed here rather than silently changed, per the guide's "combine
only at the point of use, never collapse into one shared Hill term" rule.

## Compartment order preserved — dosing by index is unchanged between the two files

The 20 ODE compartments keep exactly the same declaration order as the
original (`GUT_PHE`/`CENT_PHE`/`PERI_PHE`/`CENT_DOX`/`CENT_MET`/`CENT_BB`/
`CENT_SUNIT`/`PERI_SUNIT` then the untouched disease compartments),
confirmed identical via `/model_manifest`'s `outputPaths` for both models
— so a compartment-index dosing record (`cmt=1` for `GUT_PHE`/`PHE_gut`,
`cmt=4` for `CENT_DOX`/`DOX_C`, `cmt=5` for `CENT_MET`/`MET_C`, `cmt=6`
for `CENT_BB`/`BB_C`, `cmt=7` for `CENT_SUNIT`/`SUNIT_C`) targets the same
physical compartment in both files without any translation.

## qspserver `/model_manifest` discoverability

`C_<STEM>`/`EFFECT_<STEM>` are computed as plain `double` locals inside
`$ODE` (not `$PARAM`) — the same pattern used by every prior refactor in
this fork, because a `$PARAM @annotated` (or plain `$PARAM`) declaration
compiles to a read-only reference in mrgsolve 2.0.1, and `$ODE` cannot
assign into it. Discoverability is instead satisfied via a bare
`$CAPTURE` block, which populates `/model_manifest`'s `outputPaths`:

```
$CAPTURE C_PHE C_DOX C_MET C_BB C_SUNIT EFFECT_PHE EFFECT_DOX EFFECT_MET EFFECT_BB EFFECT_SUNIT
```

Confirmed via a live `POST /model_manifest` call against the refactored
DSL (87 parameters, 51 output paths total):

- Every renamed PK parameter (`KA_PHE`, `F_PHE`, `CL_PHE`, `V1_PHE`,
  `V2_PHE`, `Q_PHE`, `KOUT_PHE`, `KA_DOX`, `F_DOX`, `CL_DOX`, `V1_DOX`,
  `T50_DOX`, `KA_MET`, `F_MET`, `CL_MET`, `V1_MET`, `KA_BB`, `F_BB`,
  `CL_BB`, `V1_BB`, `KA_SUNIT`, `F_SUNIT`, `CL_SUNIT`, `V1_SUNIT`,
  `V2_SUNIT`, `Q_SUNIT`) and every renamed/new Hill parameter
  (`EMAX_PHE`/`EC50_PHE`/`GAMMA_PHE`, `EMAX_DOX`/`EC50_DOX`/`GAMMA_DOX`,
  `EMAX_MET`/`EC50_MET`/`GAMMA_MET`, `EMAX_BB`/`EC50_BB`/`GAMMA_BB`,
  `EMAX_SUNIT`/`EC50_SUNIT`/`GAMMA_SUNIT`) is present in `parameters` with
  its original numeric default (confirmed by inspecting the manifest
  response directly — every value matches the original file byte-for-byte
  except for the rename).
- `C_PHE`, `C_DOX`, `C_MET`, `C_BB`, `C_SUNIT`, `EFFECT_PHE`,
  `EFFECT_DOX`, `EFFECT_MET`, `EFFECT_BB`, `EFFECT_SUNIT`, and all 8
  renamed compartments (`GUT_PHE`, `CENT_PHE`, `PERI_PHE`, `CENT_DOX`,
  `CENT_MET`, `CENT_BB`, `CENT_SUNIT`, `PERI_SUNIT`) are confirmed present
  in `outputPaths` and retrievable via `/run_simulation`'s `outputs`
  selection.
- `model_content` was extracted from the `code <- '...'` R-string wrapper
  (straight text extraction) for both the original and the refactored
  file; both were used only in-memory against the local qspserver
  `mrgsolve_api` for verification — no `.cpp` sibling file was left behind
  in the repository, per the task's scratch-cleanup requirement.

## No pre-existing build defect

`ppgl_mrgsolve_model.R` compiles cleanly as-written under mrgsolve 2.0.1
(confirmed via a live `/model_manifest` call against the untouched
original — 40 output paths, all `$PARAM` values present, no error). No
`UPSTREAM_ISSUES.md` entry was needed for this file, and no
syntax-compatibility fix was applied to the refactored sibling.

## Verification

**Method.** Both DSLs (the checked-in original's `code`, extracted
verbatim from its `code <- '...'` R-string, and the refactored file's
equivalent, extracted the same way) were POSTed to the local qspserver
`mrgsolve_api` service at `http://localhost:8007`
(`POST /model_manifest` then `POST /run_simulation`, both sync), which
compiles and runs each DSL block directly with mrgsolve 2.0.1
server-side — no local R/mrgsolve install used. Requests were spaced
~2.2s apart per the task's concurrency note (`max_concurrent_jobs: 2`).
Dosing was submitted as `dosing` records (NM-TRAN-style, `cmt` as a
1-based compartment index — see the "compartment order preserved"
section above for why the same index is valid for both models),
reproducing the original file's own `dose_PHE`/`dose_DOX`/`dose_MET`/
`dose_BB`/`dose_SUNIT` helper functions' amounts/`ii`/`addl` exactly as
`run_scenario()` computes them.

Four of the file's own scenarios were reproduced:

1. **Scenario 3 (Triple therapy: PHE 60 mg/d + Metyrosine 2000 mg/d + BB
   120 mg/d)** — `GUT_PHE` (cmt 1) 60/12/303.8 amt q12h ×27 addl,
   `CENT_MET` (cmt 5) 2000/6/195.2 amt q6h ×55 addl, `CENT_BB` (cmt 6)
   120/8/259.3 amt q8h ×41 addl (all as the original's own `dose_PHE`/
   `dose_MET`/`dose_BB` compute them for `dur_h=336`), `end=720h`,
   `delta=1h` — matching `run_scenario(scenario=3, end_h=720)` exactly
   (the function's own default `delta_h=1`). Exercises PHE, MET, and BB
   together, including the Bliss-independence `alpha_block` combination.
2. **Scenario 2 (Doxazosin monotherapy, 16 mg/d)** — `CENT_DOX` (cmt 4)
   16/24/451.5 amt q24h ×13 addl, `end=720h`, `delta=1h`.
3. **Scenario 4 (Sunitinib, malignant PPGL, 37.5 mg/d)** — `CENT_SUNIT`
   (cmt 7) 37.5/24/532.6 amt q24h ×89 addl, `end=2160h` (matching the
   scenario's own `end_h <- 2160` override), **`delta=24h`** (coarsened
   from the scenario's own `delta_h=1` default purely to keep the output
   grid to a manageable 92 points for this 90-day run — same dosing,
   same duration, only the reporting grid is coarser; disclosed here per
   the guide's honesty requirement even though this was a point-count
   choice, not a solver-step-count failure).
4. **Scenario 5 (Metyrosine monotherapy, 2000 mg/d, full 720h course)** —
   `CENT_MET` (cmt 5) 2000/6/195.2 amt q6h ×119 addl (`dur_h=720` per
   `run_scenario`'s own scenario-5 branch), `end=720h`, `delta=1h`.

**Results.** Compared across all 21 shared `$CAPTURE`d disease/PK outputs
(`PHE_Conc`, `DOX_Conc`, `MET_Conc`, `BB_Conc`, `SUNIT_Conc`,
`TH_activity`, `NE_gran`, `NE_pl`, `EPI_pl`, `NMN_pl`, `MN_pl`,
`SBP_mmHg`, `DBP_mmHg`, `MAP_mmHg`, `HR_bpm`, `Tumor_mL`, `VEGF_pg`,
`Glucose_mM`, `FFA_mM`, `CgA_ng`, `AlphaBlock_frac`) plus all 20
compartments matched positionally (the 8 renamed PK compartments and the
12 untouched disease compartments):

- **Scenario 2 (doxazosin), Scenario 4 (sunitinib), Scenario 5
  (metyrosine): exact match, max abs diff = 0.0** across every shared
  output, every timepoint.
- **Scenario 3 (triple therapy): max abs diff ≈ 1.6×10⁻⁹**, on
  `AlphaBlock_frac` and `PHE_gut`/`GUT_PHE` only — floating-point-scale
  solver noise (both models integrate the identical system with the
  identical dosing; a difference at the 1e-9 level on a quantity of order
  1 is well below the precision either mrgsolve's default ODE solver
  tolerances or double-precision arithmetic can distinguish from zero,
  not a structural or numeric mismatch). Every other of the 41 compared
  output series in this scenario matched at 0.0.

Per the guide's tolerance table: this is pure structural reorganization
for all 5 compounds (archetype 3 for PHE, archetype 1 for DOX/MET/BB,
archetype 2 for SUNIT), no genuine Hill-refitting anywhere (every effect
term was already a plain ratio in the original, with `EMAX`/`EC50`/`GAMMA`
promoted from either an existing `IC50_<STEM>` parameter or an inline
literal, same value either way) — so the near-bit-exact reproduction
achieved (0.0 in 3 of 4 scenarios, ~1e-9 in the 4th) is the expected
result, not a tolerance that needed loosening.

## Diff scope confirmation

Beyond the renames, the vestigial-term/dead-parameter preservation
described above, and the `EFFECT_<STEM>`/`EMAX_<STEM>` promotion of
previously-inline scaling literals (all disclosed and value-preserving),
nothing else in the file differs: the disease-side biosynthesis, tumor,
cardiovascular, metabolic, and biomarker `$ODE` blocks (sections 6–17)
keep their original structure and every non-renamed parameter, and the R
wrapper's `run_scenario()`/summary/plotting code is unchanged except for
the compartment-name strings the dosing helpers
(`dose_PHE`/`dose_DOX`/`dose_MET`/`dose_BB`/`dose_SUNIT`) reference by
name (`"PHE_gut"`→`"GUT_PHE"`, `"DOX_C"`→`"CENT_DOX"`, `"MET_C"`→
`"CENT_MET"`, `"BB_C"`→`"CENT_BB"`, `"SUNIT_C"`→`"CENT_SUNIT"`), so the
refactored file's own scenarios still dose the renamed compartments
correctly.
