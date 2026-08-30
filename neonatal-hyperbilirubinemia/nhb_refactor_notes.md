# Refactor notes — `neonatal-hyperbilirubinemia/nhb_mrgsolve_model.R`

**Scope of this pass.** Per the fork's PK/PD refactor spec
(`FORK_WORKFLOW_GUIDE.md`, Part 2), all **three** rows this file has in
`driver-patches/data/compound_perturbation_census.md` (all classified
"Redirect concentration (clean single site)") were rewritten:
**Stannsoporfin (SNMP)**, **Phenobarbital (PB)**, and **UDCA**
(ursodeoxycholic acid). IVIG (`IGGC`/`IGGP`) is a genuine two-compartment
PK model already following a `CENT`/`PERI` shape, but it is **not** one of
the three census rows for this file and is therefore out of scope — left
byte-identical to the original. Every disease-side equation (haemolysis,
bilirubin production/binding/hepatic handling/enterohepatic shunt,
phototherapy, exchange transfusion, neurotoxicity) is untouched.

## The census mislabeled two of the three compound names — corrected, not skipped

The task instructions specifically flagged the first row, **"Bilirubin
(SNMP)"**, for a sanity check: is this actually the endogenous disease
substrate (bilirubin) mislabeled as a drug, or a real externally-dosed
compound? Checking the actual code:

- `CSNMP`/`ASNMP` (renamed `CENT_SNMP`/`GUT_SNMP`) are a genuine,
  independently-dosed one-compartment-with-depot PK pair: `dose_snmp()`
  builds an `ev()` IM bolus into `ASNMP`, with its own absorption
  (`KASNMP`), volume (`VSNMP`) and elimination (`KESNMP`) parameters, used
  by scenario S10 ("ABO isoimmune + stannsoporfin 4.5 mg/kg IM"). This is
  unambiguously stannsoporfin — a real, exogenously administered HO-1
  inhibitor drug, **not** the endogenous bilirubin pool (`BP`/`BEX`/`TSB`,
  which has no dosing route, is produced from haemoglobin catabolism, and
  is not "a clean single concentration site" in the classifier's own sense
  — it appears in a dozen mass-balance terms across the file, not one).
- The census's "Compound" column text, "Bilirubin (SNMP)", follows the same
  pattern as this file's second row, **"Hepatic handling (PB)"** — a
  process-description noun phrase (what the drug's target process is)
  paired with the drug's own correct PK stem in parentheses. Both rows
  have the right stem (`SNMP`, `PB`) but the wrong "Compound" name; neither
  is actually named after the disease-state variable it modulates.

**Conclusion: both are real, externally-administered drugs, correctly
targeted by their stems — only the census's display name was wrong.**
Corrected in `compound_perturbation_census.md` to "Stannsoporfin (SNMP)"
and "Phenobarbital (PB)" respectively, with a note explaining the original
mislabeling (rather than silently renaming with no trace). Per the task
instructions this is the "correct the classification" branch, not the
"skip, it's not a drug" branch — nothing was skipped.

## Archetype per compound

**All three (Stannsoporfin, Phenobarbital, UDCA): Archetype 3 minus
peripheral (depot + central, linear), with one bespoke deviation from the
literal template.** Each compound in the original has exactly two
compartments — a gut/IM depot and a "central" compartment — first-order
absorption, first-order elimination, no peripheral compartment, no TMDD.
That much is a clean Archetype-3-minus-peripheral match (same shape as
Anakinra in `familial-mediterranean-fever/fmf_refactor_notes.md`).

The one deviation: in the guide's Archetype 1/3 templates, `CENT_<STEM>`
holds an **amount** and `C_<STEM> = CENT_<STEM> / V1_<STEM>`. This file's
own three central compartments do not do that — each one's own `dxdt_`
equation integrates the **concentration directly**:

```
dxdt_CSNMP = KASNMP*ASNMP/(VSNMP*Wc)*10.0 - KESNMP*CSNMP;   // original
```

(`ASNMP/(V*Wc)*10.0` converts the depot's first-order outflow, in mg/h,
straight into a concentration rate, in ug/mL/h — the `*10.0` is a
mg/dL→ug/mL unit conversion, not a bioavailability factor; there is no
explicit `F_<STEM>` in the original for any of the three compounds, so
none was invented). Because `Wc` (body weight) is itself a time-varying
state elsewhere in this file, `d(amount/V)/dt ≠ d(amount)/dt / V` when `V`
changes with time — so converting this to a textbook amount-then-divide
form would **change** the original's own dynamics, not merely rename them.
To honour "pure structural reorganization should match near-exactly" and
"never flatten a model because it's inconvenient," the central
compartments (`CENT_SNMP`/`CENT_PB`/`CENT_UDCA`) were kept as
concentration-integrating states, renamed only, and `C_<STEM>` is defined
as a plain identity of the compartment (`C_SNMP = CENT_SNMP;`, etc.) rather
than a division. This is disclosed here as a bespoke deviation from the
literal Archetype-3 template, per the guide's explicit allowance ("if a
compound's PK genuinely doesn't resemble any archetype... rename to the
convention, isolate it from PD, expose `C_<STEM>`, and note... why").

## Renaming applied

| Original | Refactored | Value | Role |
|---|---|---|---|
| `ASNMP` (cmt) | `GUT_SNMP` | -- | depot (IM) |
| `CSNMP` (cmt) | `CENT_SNMP` | -- | central, concentration-state (bespoke, see above) |
| `KASNMP` | `KA_SNMP` | 0.173 (1/h) | absorption rate |
| `VSNMP` | `V1_SNMP` | 3.00 (dL/kg) | volume of distribution |
| `KESNMP` | `KE_SNMP` | 0.023 (1/h) | elimination rate constant |
| `KISNMP` | `EC50_SNMP` | 8.0 (ug/mL) | Hill EC50 (was named "Ki") |
| -- (none) | `EMAX_SNMP` (new) | 1.0 | Hill Emax [original had none explicit] |
| -- (none) | `GAMMA_SNMP` (new) | 1.0 | Hill exponent [original had none explicit] |
| `CSNMP` (local, identity) | `C_SNMP` | -- | **the exposed concentration** |
| -- (none, was inline in `F_SNMP`) | `EFFECT_SNMP` (new name) | -- | fractional HO-1 inhibition |
| `APB` (cmt) | `GUT_PB` | -- | depot |
| `CPB` (cmt) | `CENT_PB` | -- | central, concentration-state (bespoke, see above) |
| `KAPB` | `KA_PB` | 0.400 (1/h) | absorption rate |
| `VPB` | `V1_PB` | 9.00 (dL/kg) | volume of distribution |
| `KEPB` | `KE_PB` | 0.00693 (1/h) | elimination rate constant |
| `EPB` | `EMAX_PB` | 2.50 | Hill Emax (induction fold-1) |
| `EC50PB` | `EC50_PB` | 20.0 (mg/L) | Hill EC50 |
| -- (none) | `GAMMA_PB` (new) | 1.0 | Hill exponent [original had none explicit] |
| `CPB` (local, identity) | `C_PB` | -- | **the exposed concentration** |
| -- (none, was inline in `F_PB`) | `EFFECT_PB` (new name) | -- | UGT1A1-induction effect |
| `AUDCA` (cmt) | `GUT_UDCA` | -- | depot |
| `CUDCA` (cmt) | `CENT_UDCA` | -- | central, concentration-state (bespoke, see above) |
| `KAUDCA` | `KA_UDCA` | 0.500 (1/h) | absorption rate |
| `VUDCA` | `V1_UDCA` | 3.00 (dL/kg) | volume of distribution |
| `KEUDCA` | `KE_UDCA` | 0.140 (1/h) | elimination rate constant |
| `EMAXU` | `EMAX_UDCA` | 0.80 | Hill Emax |
| `EC50UDCA` | `EC50_UDCA` | 4.00 (umol/L) | Hill EC50 |
| -- (none) | `GAMMA_UDCA` (new) | 1.0 | Hill exponent [original had none explicit] |
| `CUDCA` (local, identity) | `C_UDCA` | -- | **the exposed concentration** |
| `F_UDCA` (was already the plain Hill ratio) | `EFFECT_UDCA` (new name) | -- | enterohepatic-shunt effect |

All parameter *values* are copied verbatim from the original. `F_SNMP`,
`F_PB`, `F_UDCA` are kept as local `$ODE` variables (disease-facing
survival/induction factors, see below) — they are not part of the guide's
naming convention themselves, just the point where `EFFECT_<STEM>` is
algebraically converted into the shape the disease equations already
expect.

## Hill interface: rename, not a fit — for all three compounds

Each original effect term was already exactly the Hill/Emax shape (or its
algebraic complement), so no `nls()` fit was performed for any of the
three compounds:

- **UDCA**: `F_UDCA = EMAXU*CUDCA/(EC50UDCA+CUDCA)` **is already**
  `EFFECT_UDCA` — a direct, one-to-one rename with `GAMMA_UDCA = 1` added.
- **Phenobarbital**: `F_PB = 1.0 + EPB*CPB/(EC50PB+CPB)`. The ratio itself
  is `EFFECT_PB`; the disease equation still reads `F_PB = 1.0 +
  EFFECT_PB`, unchanged in shape.
- **Stannsoporfin**: `F_SNMP = 1.0/(1.0 + CSNMP/KISNMP)` is a competitive
  **inhibitor** term, not a plain increasing ratio. Algebraically,
  `1/(1+C/Ki) = Ki/(Ki+C) = 1 - C/(Ki+C)`, so `EFFECT_SNMP = C/(Ki+C)`
  (the plain Hill ratio, `EMAX=GAMMA=1`) is the named interface, and the
  disease equation reads `F_SNMP = 1.0 - EFFECT_SNMP` — the same algebraic
  complement, not a new shape.

`EMAX_SNMP = 1.0` and all three `GAMMA_<STEM> = 1.0` are new named
parameters (none of the three had an explicit Hill exponent or a
separately-named Emax in the original); their values reproduce the
original's implicit Emax=1/gamma=1 shape exactly.

## A design pitfall found and fixed during this refactor (not an upstream defect)

While the first draft (declaring `double C_SNMP = CENT_SNMP;` inside
`$ODE` and then re-declaring `double C_SNMP = ...;` again inside `$TABLE`
to route it to `$CAPTURE`, mirroring how one might naively read the
guide's Archetype templates) was submitted to `POST /model_manifest`, it
failed to compile:

```
135:10: error: redefinition of 'double {anonymous}::C_SNMP'
54:10: note: 'double {anonymous}::C_SNMP' previously declared here
```

(and the same for `C_PB`, `C_UDCA`, `EFFECT_SNMP`, `EFFECT_PB`,
`EFFECT_UDCA`). This is the same underlying mrgsolve 2.0.1 mechanism as
the pre-existing upstream defects logged in `UPSTREAM_ISSUES.md` #48/#51
(block-local `double` declarations across different DSL blocks are hoisted
into one shared anonymous C++ namespace, so identical names in two blocks
collide even though the DSL author sees them as textually separate) — but
here **the collision was introduced by this refactor's own first draft**,
not present in the checked-in original, so it is **not** logged as an
upstream issue. It is fixed directly in the delivered
`nhb_mrgsolve_model_refactored.R`: `C_SNMP`, `C_PB`, `C_UDCA`,
`EFFECT_SNMP`, `EFFECT_PB`, `EFFECT_UDCA` are declared exactly once, as
file-scope globals (the same pattern the file's own author already used
for `CPNOW`/`CISONOW`/etc. — see the block right before `$MAIN`), assigned
(no `double` keyword) inside `$ODE`, and captured directly in `$TABLE`/
`$CAPTURE` with no re-declaration. Confirmed fixed by the second
`/model_manifest` call succeeding and listing all six names in
`outputPaths` (see below).

The original `nhb_mrgsolve_model.R` itself has **no known build defect** —
its own `code <- '...'` block compiled cleanly on the first
`POST /model_manifest` call, before any refactor changes were made.

## Verification

**Method.** Both files' embedded `nhb_code <- '...'` DSL blocks were
mechanically extracted (regex on the assignment, verbatim quoted text —
23,623 chars original, 27,774 chars refactored) and POSTed to the local
qspserver `mrgsolve_api` service at `http://localhost:8007`
(`/model_manifest` then `/run_simulation`), which compiles and runs each
DSL block directly with mrgsolve 2.0.1 server-side — no local R/mrgsolve
install used. Requests were spaced several seconds apart and run
sequentially (never more than one in flight), respecting the service's
`max_concurrent_jobs: 2` limit and its recent history of crashing under
concurrent load. `POST /run_simulation`'s `dosing` field addresses
compartments by 1-based index, not name; compartment order is identical
between the two files (only names changed), confirmed from
`/model_manifest`'s `outputPaths` (`GUT_SNMP`=23, `GUT_PB`=25,
`GUT_UDCA`=29 in both).

**Scenarios run — the file's own, not invented.** All three are full-
duration (no shortening needed; none of the three approached the API's
default `maxsteps` budget):

1. **Scenario S10, "ABO isoimmune + stannsoporfin 4.5 mg/kg IM"**
   (`AUTOPT=1, IRRSET=30, FBSASET=0.80, PTOFFM=2.0, ABMAT0=0.12, RF=1`;
   15.3 mg bolus at t=24h into `ASNMP`/`GUT_SNMP`, cmt 23; `end=336,
   delta=0.5`, 674 points): **near-exact match** across the 49 shared
   disease-side outputs (all `$CMT` compartments other than the three
   drugs' own PK, plus all pre-existing `$CAPTURE` names) — **max abs diff
   0.0** on 48 of 49, and **1e-4 abs (≈4e-8 relative)** on `BSAOUT` only,
   appearing solely in the last few timepoints (from t≈334.5h). The
   compound's own PK, `CSNMP` vs `CENT_SNMP`, is an **exact match, 0.0**,
   at every timepoint (peak 11.04 ug/mL in both). This tiny, late-only
   deviation is consistent with the guide's "floating-point-scale
   deviation... not a bug" allowance — attributable to the `pow(x,1.0)`
   calls the Hill rewrite introduces (mathematically identical to `x`, not
   bit-identical to the FPU) very slowly compounding through an adaptive
   LSODA stepper over 334h of integration, not a shape or parameter
   change: it does not appear at all in the other two scenarios below,
   which share the same growth/weight ODE but not the SNMP dosing.
2. **Scenario S9, "Crigler-Najjar type II + phenobarbital 5 mg/kg/day"**
   (`AUTOPT=1, IRRSET=30, FBSASET=0.80, PTOFFM=2.0, GENO=0.05`; 17.0 mg
   into `APB`/`GUT_PB`, cmt 25, q24h ×14 doses; `end=336, delta=0.5`, 674
   points): **exact match, max abs diff 0.0** across all 49 shared
   outputs. `CPB` vs `CENT_PB`: exact match, 0.0 (peak 31.73 mg/L in
   both).
3. **The file's own "UDCA 10 mg/kg q12h" arm** (from `ehc_interventions()`
   — dose_udca(3.40), i.e. 34.0 mg into `AUDCA`/`GUT_UDCA`, cmt 29, q12h
   ×28 doses, all other parameters at baseline; `end=240, delta=6`, 41
   points, matching that function's own call): **exact match, max abs
   diff 0.0** across all 49 shared outputs. `CUDCA` vs `CENT_UDCA`: exact
   match, 0.0 (peak 22.94 umol/L in both).

Two of three compounds (Phenobarbital, UDCA) verify bit-for-bit exact,
consistent with the guide's tolerance table for Archetypes 1–3 ("pure
structural reorganization... expect a near-exact match"). Stannsoporfin
verifies to floating-point-scale precision (4e-8 relative) rather than
bit-exact, isolated to a single late-appearing output and explained above;
this is disclosed rather than smoothed over, per "a gate failure is a
question, not a verdict."

## qspserver `/model_manifest` discoverability

Confirmed via `POST /model_manifest` against the extracted, renamed DSL:
`KA_SNMP`, `V1_SNMP`, `KE_SNMP`, `EC50_SNMP`, `EMAX_SNMP`, `GAMMA_SNMP`,
`KA_PB`, `V1_PB`, `KE_PB`, `EC50_PB`, `EMAX_PB`, `GAMMA_PB`, `KA_UDCA`,
`V1_UDCA`, `KE_UDCA`, `EC50_UDCA`, `EMAX_UDCA`, `GAMMA_UDCA` all appear in
the manifest's `parameters` with their original numeric defaults.
`C_SNMP`, `C_PB`, `C_UDCA`, `EFFECT_SNMP`, `EFFECT_PB`, `EFFECT_UDCA` are
state-derived (computed in `$ODE` from the compartment concentrations),
so — per the same reasoning established for `C_COL`/`EFFECT_ANA`/etc. in
`familial-mediterranean-fever/fmf_refactor_notes.md` — they cannot also be
`$PARAM` entries; all six appear in the manifest's `outputPaths` via the
file's own (extended) `$CAPTURE` list, confirmed discoverable. `GUT_SNMP`,
`CENT_SNMP`, `GUT_PB`, `CENT_PB`, `GUT_UDCA`, `CENT_UDCA` also appear in
`outputPaths` as ordinary compartments (indices 23–26, 29–30), unchanged
from the original's compartment order.

No `.cpp` extraction file was left behind — extraction was in-memory only,
used to build the verification requests above and then discarded, per the
workflow guide.

## Anything else flagged

- No compound other than Stannsoporfin, Phenobarbital, and UDCA was
  touched. IVIG's parameters (`VIGG`, `KIGGCP`, `KIGGPC`, `KIGGEL`, `IGG0`,
  `IMAXIVIG`, `IC50IVIG`) and its two compartments (`IGGC`, `IGGP`) are
  byte-identical to the original — it is not one of this file's three
  census rows.
- The R-side scenario list (`scenarios`), `dose_snmp()`/`dose_pheno()`/
  `dose_udca()`, and `ehc_interventions()` were updated only where they
  name a renamed compartment (`"ASNMP"→"GUT_SNMP"`, `"APB"→"GUT_PB"`,
  `"AUDCA"→"GUT_UDCA"` inside the three `ev()` calls) — same dosing
  amounts, same timing throughout. No other R-side code references any of
  the renamed parameters or compartments by name.
- All ten shipped scenarios (`S1`–`S10`) and all eight analysis functions
  (`iso_bf_table`, `pt_dose_response`, `photoisomer_run`, `exchange_run`,
  `cn_ceiling`, `gene_therapy_run`, `genotype_feeding`,
  `ehc_interventions`) still reference the model object the same way; none
  needed changes beyond the three `ev()` compartment-name updates above.

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`, the
three `neonatal-hyperbilirubinemia` rows — the first two also corrected
from their original, misleading "Bilirubin (SNMP)" / "Hepatic handling
(PB)" compound names to "Stannsoporfin (SNMP)" / "Phenobarbital (PB)" (see
above).
