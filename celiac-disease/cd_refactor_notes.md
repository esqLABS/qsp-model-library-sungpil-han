# Refactor notes — `celiac-disease/cd_mrgsolve_model.R`

**Scope of this pass.** This file has exactly **one** row in
`driver-patches/data/compound_perturbation_census.md`, generically labeled
**"DRUG"** and classified "Normalize duplicate concentration sites, then
redirect." Every disease-side equation (gluten processing, tTG2, intestinal
permeability, innate IL-15, IEL/CD4T/Th1/Th17/B-cell/serology dynamics,
Marsh histopathology, nutritional sequelae) is untouched.

## What "DRUG" turned out to be

The file's own `$PARAM` block resolves the ambiguity directly:

```
// --- Drug PK parameters (Larazotide example) ---
F_oral    : 0.5  : Bioavailability (gut-local drug; systemic F~0.01)
ka_drug   : 0.6  : h-1, first-order absorption rate
CL_drug_L_h : 8.0 : L/h, total drug clearance
Vd_drug_L : 6.0  : L, volume of distribution
MW_drug   : 406  : g/mol, molecular weight (larazotide)
```

The PK block's own comments name the compound: **larazotide (AT-1001)**.
This is the real identity behind the generic "DRUG" census label, and the
stem used throughout this refactor is **`LARA`**.

## A real quirk in the original: one PK block, three selectable PD targets

`Drug_type` (0/1/2/3) is a disease-model switch, not a second or third
compound. The original gives **only one** PK compartment pair
(`DrugGut`/`DrugPlasma`), and its single concentration (`Cp`) feeds three
mutually-exclusive Hill terms depending on which "drug" is selected for a
scenario:

```r
double Cp = DrugPlasma;
double E_lara = (Drug_type == 1) ? Emax_drug * Cp / (Cp + EC50_lara_IP)  : 0.0;
double E_ZED  = (Drug_type == 2) ? Emax_drug * Cp / (Cp + EC50_ZED_tTG)  : 0.0;
double E_AMG  = (Drug_type == 3) ? Emax_drug * Cp / (Cp + EC50_AMG_IL15) : 0.0;
```

ZED1227 and AMG714 never get their own absorption/clearance/volume
parameters anywhere in the file — dosing "ZED1227" or "AMG714" (scenarios 5
and 6) really means "dose the same generic larazotide-parameterized PK
compartment, then read out the tTG2- or IL-15-targeted branch instead of
the IP-targeted one." This is a real simplification in the original, not
something introduced by this refactor, and it is disclosed rather than
"fixed": per the guide, parameter values are never invented, so no
separate PK was manufactured for ZED1227/AMG714 to make the file look more
mechanistically distinct than it actually is.

**Scope consequence.** Only the larazotide identity (`LARA`) is this
census row's compound and is renamed to the naming convention. `E_ZED` and
`E_AMG` are **not** renamed to `EFFECT_ZED`/`EFFECT_AMG` — they have no
census row and no PK of their own to redirect. They are left with their
original names and exact original formulas, touched only in the single,
forced, mechanical way that removing the duplicate `Cp`/`Emax_drug` names
requires: they now read `C_LARA`/`EMAX_LARA` instead of `Cp`/`Emax_drug` —
literally the same shared variables, renamed once, not a second compound's
parameters being renamed.

## The census classification, applied literally

"Normalize duplicate concentration sites, then redirect" describes exactly
two pre-existing sites that both re-read `DrugPlasma` under different
names:

1. `double Cp = DrugPlasma;` inside `$ODE` (feeds `E_lara`/`E_ZED`/`E_AMG`)
2. `capture Drug_Cp = DrugPlasma;` inside `$TABLE` (a second, independent
   re-reading, exposed only for reporting)

These are normalized into one file-scope variable, `C_LARA`, declared once
in a new `$GLOBAL` block (the original had no `$GLOBAL` block at all) and
assigned once inside `$ODE` — the same collision-avoidance pattern already
established in `neonatal-hyperbilirubinemia/nhb_refactor_notes.md` (mrgsolve
2.0.1 hoists every block-local `double NAME` into one shared anonymous
C++ namespace, so declaring the same name as a local in two different DSL
blocks collides even though each block is textually separate). `Drug_Cp`
is removed; `C_LARA` is captured directly instead — a genuine
single-site concentration definition where the original had two.

## Archetype

**Archetype 3 minus peripheral (depot + central, linear), with one bespoke
deviation** — the same deviation already documented for
Stannsoporfin/Phenobarbital/UDCA in
`neonatal-hyperbilirubinemia/nhb_refactor_notes.md`: the original's
"central" compartment integrates **concentration directly**, not an amount
later divided by volume:

```
dxdt_DrugPlasma = F_oral * ka_drug * DrugGut * 1000.0 / Vd_drug_L
    - (CL_drug_L_h / Vd_drug_L) * DrugPlasma;
```

`DrugGut` (mg, amount) feeds this rate as mg/h, converted through
`*1000.0/Vd_drug_L` straight into ng/mL/h — `DrugPlasma` is declared and
used throughout as ng/mL, never divided by volume elsewhere. Renamed
`CENT_LARA` keeps this exact structure; `C_LARA = CENT_LARA` is an
identity alias, not a division, matching the original's own math exactly
rather than flattening it into the literal Archetype-3 amount/volume
template.

## Renaming applied

| Original | Refactored | Value | Role |
|---|---|---|---|
| `DrugGut` (cmt) | `GUT_LARA` | -- | depot (oral) |
| `DrugPlasma` (cmt) | `CENT_LARA` | -- | central, concentration-state (bespoke, see above) |
| `F_oral` | `F_LARA` | 0.5 | bioavailability |
| `ka_drug` | `KA_LARA` | 0.6 (1/h) | absorption rate |
| `CL_drug_L_h` | `CL_LARA` | 8.0 (L/h) | clearance |
| `Vd_drug_L` | `V1_LARA` | 6.0 (L) | volume of distribution |
| `MW_drug` | `MW_LARA` | 406 (g/mol) | molecular weight (informational; never used in `$ODE`) |
| `EC50_lara_IP` | `EC50_LARA` | 50 (ng/mL) | Hill EC50 |
| `Emax_drug` | `EMAX_LARA` | 0.85 | Hill Emax — **shared verbatim** by the original across all three `Drug_type` branches (see above); renaming it does not change what `E_ZED`/`E_AMG` compute |
| -- (none) | `GAMMA_LARA` (new) | 1 | Hill exponent [implicit in the original] |
| `Cp` / `Drug_Cp` (two duplicate sites) | `C_LARA` | -- | **the exposed concentration** (single, normalized site) |
| `E_lara` | `EFFECT_LARA` | -- | rename, not a fit — this census row's own named Hill interface |
| `E_ZED` | `E_ZED` (unchanged name) | -- | **not in scope** — no separate PK/census row; only its `Cp`/`Emax_drug` reads were redirected to `C_LARA`/`EMAX_LARA` |
| `E_AMG` | `E_AMG` (unchanged name) | -- | **not in scope** — same as above |
| `EC50_ZED_tTG`, `EC50_AMG_IL15` | unchanged | 80, 15 | not this census row's parameters |

All parameter *values* are copied verbatim from the original. `Drug_type`
is untouched — it is the disease model's own selector logic, not a PK
identifier.

## Hill interface: rename, not a fit

`E_lara = Emax_drug * Cp / (Cp + EC50_lara_IP)` **is already** the
canonical `EMAX*C^gamma/(EC50^gamma+C^gamma)` shape with implicit
`gamma=1`. `EFFECT_LARA` is a one-to-one rename with `GAMMA_LARA = 1`
added explicitly (`pow(C_LARA,1)` is mathematically identical to
`C_LARA`, not a behavioral change). No `nls()` fit was needed or
performed.

## When the original doesn't compile at all

The untouched original fails to compile under mrgsolve 2.0.1:

```
Error in validObject(.Object) :
  invalid class "mrgmod" object: compartment should not be in $CAPTURE:
  AntiTTG,VH,CrD,AbsArea
```

`AntiTTG`, `VH`, `CrD`, `AbsArea` are declared as compartments in `$INIT
@annotated` and separately re-listed in `$CAPTURE` — mrgsolve 2.0.1
rejects a compartment name that also appears in `$CAPTURE`. This is
unrelated to the drug PK/PD block in scope for this refactor (all four
names belong to the disease side of the model: serology and
histopathology state variables). Logged as
`translations/UPSTREAM_ISSUES.md` **#106**; the checked-in
`cd_mrgsolve_model.R` is completely untouched and still carries the defect
exactly as written.

**Fix applied directly to the delivered `cd_mrgsolve_model_refactored.R`**
(per the guide's settled policy — syntax-only, non-numeric): the four
compartment names are dropped from the `$CAPTURE` line. Compartments are
always present in mrgsolve's output regardless of `$CAPTURE`, so nothing
about what is reported changes; `AntiTTG`/`VH`/`CrD`/`AbsArea` remain
fully readable in every simulation output as ordinary compartments.

## Verification

**Method.** Both files' embedded `cd_model_code <- '...'` DSL blocks were
mechanically extracted (find the quoted-string assignment, pull the
contents verbatim — 11,102 chars original, 13,822 chars refactored) and
POSTed to the local qspserver `mrgsolve_api` service at
`http://localhost:8007` (`/model_manifest` then `/run_simulation`), which
compiles and runs each DSL block directly with mrgsolve 2.0.1 server-side
— no local R/mrgsolve install used. For the "original" side of the
comparison, the same single, disclosed `$CAPTURE` build-compat fix
described above was applied to a scratch copy only (never to the
checked-in `cd_mrgsolve_model.R`) so that a compiled baseline exists to
verify against at all. Requests were spaced ~2.5s apart and run
sequentially (never more than one in flight), respecting the service's
`max_concurrent_jobs: 2` limit and its history of crashing under
concurrent load. `POST /run_simulation`'s `dosing` field addresses
compartments by 1-based index; confirmed identical between both files via
`/model_manifest`'s `outputPaths` (`GUT_LARA`/`DrugGut` both index 19,
`CENT_LARA`/`DrugPlasma` both index 20).

**Scenarios run — all six of the file's own, not invented, full 1-year
duration (`end=8760h, delta=24h`, 366 points; no shortening needed, none
of the six approached the API's default `maxsteps` budget):**

1. `1_Untreated_Normal_Diet` (GFD=0, no drug)
2. `2_Strict_GFD` (GFD=1, no leak, no drug)
3. `3_Partial_GFD_5pct_leak` (GFD=1, 5% leak, no drug)
4. `4_GFD_plus_Larazotide` (GFD=1, 10% leak, Drug_type=1, 2 mg q8h — 1095
   doses, `addl=1094`)
5. `5_GFD_plus_ZED1227` (GFD=1, 10% leak, Drug_type=2, 300 mg q24h — 365
   doses, `addl=364`)
6. `6_GFD_plus_AMG714_RCD` (GFD=1, 20% leak, Drug_type=3, 150 mg q168h —
   52 doses, `addl=51`)

All dosing amounts/intervals/durations reproduce the file's own
`scenarios` list and `run_scenario()` dosing logic exactly (the original's
R-side event object also sets `rate=-2`, an mrgsolve-local infusion-duration
idiom relevant only to a local `mrgsim()`/`ev()` call; the qspserver API's
own `DoseSpec` has no equivalent field and does not need one — ordinary
bolus dosing at the same amounts/times/compartment reproduces the same
scenario).

**Result: exact match, max abs diff 0.0**, across all 26 shared outputs
(18 disease-state compartments + 8 `$CAPTURE`d derived outputs) plus the
compound's own PK pair (`DrugGut`/`GUT_LARA`, `DrugPlasma`/`CENT_LARA`),
at every one of 366 timepoints, in **all six scenarios** — including both
scenarios (5, 6) that exercise the untouched `E_ZED`/`E_AMG` branches now
reading the renamed `C_LARA`/`EMAX_LARA`. `C_LARA` and `EFFECT_LARA` (new,
not present in the original) were sanity-checked directly: `C_LARA` is
nonzero in all three drug scenarios (peak 1.13 ng/mL under larazotide
dosing, smaller residual peaks of 0.011/0.006 ng/mL under the ZED1227/
AMG714 dosing amounts, since all three share the same PK), while
`EFFECT_LARA` is nonzero only when `Drug_type==1` (peak 0.019, scenario 4)
and exactly zero in every other scenario — confirming the Hill interface
fires only for its own compound, per the guide's naming-convention item 4.

This is a pure structural reorganization (rename + duplicate-site
normalization, no Hill-fitting), consistent with the guide's tolerance
table for Archetypes 1-3: "expect a near-exact match... anything beyond
floating-point-scale deviation means a bug." Here the match is bit-exact,
not merely floating-point-scale.

## qspserver `/model_manifest` discoverability

Confirmed via `POST /model_manifest` against the extracted, renamed DSL:
`F_LARA`, `KA_LARA`, `CL_LARA`, `V1_LARA`, `MW_LARA`, `EC50_LARA`,
`EMAX_LARA`, `GAMMA_LARA`, `EC50_ZED_tTG`, `EC50_AMG_IL15`, `Drug_type` all
appear in the manifest's `parameters` with their original numeric
defaults. `C_LARA` and `EFFECT_LARA` are state-derived (computed in `$ODE`
from `CENT_LARA` each step), so — per the same reasoning established for
`C_SNMP`/`EFFECT_SNMP` etc. in
`neonatal-hyperbilirubinemia/nhb_refactor_notes.md` — they cannot also be
`$PARAM` entries (mrgsolve 2.0.1 treats `$PARAM` members as read-only in
`$ODE`); both appear in the manifest's `outputPaths` via `$CAPTURE`,
confirmed discoverable. `GUT_LARA`/`CENT_LARA` also appear in
`outputPaths` as ordinary compartments (indices 19-20), unchanged position
from the original's compartment order.

No `.cpp` extraction file was left behind — extraction was scratch-only,
used to build the verification requests above and then discarded.

## Anything else flagged

- No compound other than larazotide (`LARA`) has a census row in this
  file; `E_ZED`/`E_AMG` and their `EC50_ZED_tTG`/`EC50_AMG_IL15` parameters
  are left with their original names and exact original formulas (see
  "Scope consequence" above).
- The R-side `scenarios` list and `run_scenario()` were updated only where
  they name the renamed compartment (`"DrugGut"` -> `"GUT_LARA"`, two
  occurrences in the `ev()` calls). Same dosing amounts, same timing,
  same six scenarios throughout. No other R-side code references any
  renamed parameter or compartment by name.
- `cd_shiny_app.R` and `cd_qsp_model.dot`/`.svg`/`.png` are untouched —
  out of scope for this refactor (only the mrgsolve model file and its
  siblings named in the task are deliverables).

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`: the
`celiac-disease` row's "Compound" column corrected from generic "DRUG" to
"Larazotide (LARA)", target/pathway filled in (paracellular tight-junction
permeability / zonulin pathway), and the notes column filled in per the
above.

## Discoverability fix

A corpus-wide discoverability audit found `C_LARA` was never written as a
single contiguous `double C_<STEM> = <expr>;` statement anywhere in the
file — it was `$GLOBAL`-forward-declared (`double C_LARA, EFFECT_LARA;`)
and bare-assigned exactly once, in `$ODE` (`C_LARA = CENT_LARA;`, line
~212). Correct, working code (not a bug — this file's own `$GLOBAL` comment
already explains this pattern exists specifically to avoid an $ODE-vs-
$TABLE duplicate-declaration collision), but not literal-text-discoverable
by tooling that regexes for `double C_<STEM> = ...;`.

Given the exact collision this file's own `$GLOBAL` comment warns about,
a naive "just add a new `double C_LARA = CENT_LARA;` line in `$TABLE`,
leave the `$GLOBAL` forward-declare untouched" was not attempted blind —
the identical failure mode was hit and diagnosed first on
`breast-cancer/bc_mrgsolve_model_refactored.R` in this same batch (see
that file's own refactor notes): mrgsolve auto-declares a persistent,
`$CAPTURE`-visible class member for *every* `double NAME = ...;`
initializing statement found anywhere in the block, so adding one in
`$TABLE` while a bare `double C_LARA;` forward-declare still exists in
`$GLOBAL` produces two competing declarations of the same member
("reference to 'C_LARA' is ambiguous").

**Fix applied:** removed `C_LARA` from the `$GLOBAL` bare forward-declare
(`double C_LARA, EFFECT_LARA;` → `double EFFECT_LARA;`), so the new
`$TABLE` line becomes its sole declaration — the same mechanism this file
already uses for its `capture NAME = expr;` outputs (`VH_CD_ratio`,
`Marsh_score`, etc.), just spelled with an explicit `double` instead of
`capture` (both trigger the identical auto-declare, so mixing them for the
same name is what causes the collision — using only one form per name is
what avoids it). `EFFECT_LARA` is untouched (out of scope for this fix;
still forward-declared in `$GLOBAL` and bare-assigned once in `$ODE`,
exactly as before). The `$ODE` bare assignment (`C_LARA = CENT_LARA;`) and
every downstream read (`EFFECT_LARA`, `E_ZED`, `E_AMG`) are untouched —
they now simply target the auto-declared member instead of the
manually-declared one, with identical storage semantics.

Added line, immediately before `$CAPTURE`:

```
double C_LARA = CENT_LARA;
```

**Verification:** `Rscript -e 'parse(...)'` succeeds with no error. The
DSL string was extracted and posted to qspserver's `mrgsolve_api` at
`localhost:8007`: `/model_manifest` compiles cleanly and lists `C_LARA` in
`outputPaths` alongside `EC50_LARA` in `parameters`. `/run_simulation` was
run for the file's own Scenario 4 ("GFD + Larazotide", `Drug_type=1`,
`GFD=1`, `GFD_leak=0.10`, 2 mg TID into `GUT_LARA` — reproduced as
`ii=8, addl=1094` bolus dosing over the file's own 8760 h / daily-sampling
horizon) against both the pre-fix original DSL (`git show HEAD:...`) and
the fixed DSL, identical dosing and parameters. `C_LARA`, `EFFECT_LARA`,
and `VH_CD_ratio` (367 time points) were numerically identical between the
two runs (max abs diff = 0 for every captured column) — confirming the fix
changes nothing numeric. Grep-confirmed `double C_LARA = ` and `EC50_LARA`
both now appear in the file.
