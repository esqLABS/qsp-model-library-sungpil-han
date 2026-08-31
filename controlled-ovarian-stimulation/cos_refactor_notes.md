# Refactor notes — `controlled-ovarian-stimulation/cos_mrgsolve_model.R`

## Scope: 5 compounds refactored (the 3 assigned + 2 found during audit); 3 more found and deliberately left out of scope

The census carried three rows for this file, all classified "Normalize
duplicate concentration sites, then redirect": `AGO`, `ANT`, `HCG`, target
unknown (`?`). Checked against the actual code:

| Row | Real identity | Evidence |
|---|---|---|
| `AGO` | **Triptorelin** (a GnRH agonist) | The original's own `$PARAM` comments name it explicitly: `KAAGO : ... triptorelin absorption`, `KELAGO : ... triptorelin elimination`, `VAGO : ... triptorelin V/F`. The stem `AGO` is a class abbreviation the code itself already disambiguates. |
| `ANT` | **Ganirelix** (a GnRH antagonist) | The header's own calibration section is explicit and cites a real PMID: "ganirelix 0.25 mg s.c.: Cmax 11.2 ng/mL -> 11.4, t1/2 13 h, LH suppression 70-80% -> 77% (PMID 10593372)". The protocol builder's own comment (`ant_dose = 0.25, # mg/d ganirelix/cetrorelix`) allows cetrorelix as an untested alternative label, but the calibration target — the only PK data actually fit — is ganirelix. |
| `HCG` | **Human chorionic gonadotropin** (the trigger shot) | Already a real compound name, not a class placeholder. Also carries hp-hMG's dosed LH activity in the model's own `hmg_lh` protocol arm, added into the same `HCGD`/`HCG` compartment (not a second compound — same receptor, same PK, just a second dosing route into the existing pool). |

The task's own brief flagged "an FSH/gonadotropin stimulation drug" as a
plausible corpus-census undercount, the same failure mode already
documented for `prostate-cancer`, `hypercalcemia-of-malignancy`, and
`alagille-syndrome` in this batch. Reading the model's own `$PARAM`/`$CMT`/
`$ODE` blocks confirmed this: the file has **eight** externally-dosed
compounds total, not three.

| Compound | Stem | Dosed in the file's own scenarios | In census? |
|---|---|---|---|
| Recombinant FSH | `FSH` | s01-s21 (daily stimulation, most scenarios) | **No — added below** |
| Corifollitropin alfa | `CORI` | s10 (`cori = 150`) | **No — added below** |
| Ganirelix | `ANT` | s01-s10, s12-s21 (antagonist protocol) | Yes (was `?`) |
| hCG | `HCG` | s01, s03, s05-s09, s12, s14-s18, s20-s21 (trigger, dual trigger, hp-hMG) | Yes (was `?`) |
| Triptorelin | `AGO` | s04, s12, s13, s17, s19 (agonist/dual trigger) | Yes (was `?`) |
| Cabergoline | `CAB` | s06, s20 (OHSS prevention) | No — **out of scope, see below** |
| Letrozole | `LET` | none of the file's 22 shipped scenarios sets `letro = TRUE` (the protocol builder supports it; no example exercises it) | No — **out of scope, see below** |
| Exogenous vaginal/IM progesterone | `P4D`/`P4` | s01-s10, s12-s21 (`luteal_p4 = 600`) | No — **out of scope, see below** |

`FSH` and `CORI` are refactored below, alongside the three assigned
compounds — five total. `CAB`, `LET`, and the exogenous-progesterone depot
are real, separately-dosed compounds this file also models, found during
the same audit, but are **deliberately left out of scope** for this
refactor (same-in-kind decision as `clostridioides-difficile-infection`
scoping fidaxomicin/metronidazole/rifaximin/ridinilazole/the index
antibiotic/the live biotherapeutic out of its own refactor): the task's own
brief bounded this refactor to the three assigned compounds plus checking
specifically for an undercounted FSH/gonadotropin drug, and the corpus
convention (confirmed by that cdi precedent) is that finding a real
compound during an audit does not, by itself, obligate refactoring it —
disclosure is the requirement, not exhaustive coverage. No census row is
added for `CAB`/`LET`/`P4D` (matching that same precedent: out-of-scope
compounds are disclosed here, not given half-finished census entries).
`CAB`/`LET`/`P4D` are left completely untouched: same compartment names
(`CAB`, `LET`, `P4D`), same parameter names (`KELCAB`, `EMAXCAB`,
`EC50CAB`, `KELLET`, `EMAXLET`, `EC50LET`, `KAP4`, `FP4`), same equations,
byte-for-byte, in the delivered `_refactored.R`.

## The original does not compile under mrgsolve 2.0.1 — two pre-existing, unrelated defects

`POST /model_manifest` on the untouched original's own extracted DSL
(`cos_code <- '...'`) fails before any of the five in-scope compounds'
code is even reached. Logged as **`translations/UPSTREAM_ISSUES.md`
entry #132**:

1. **Loop-counter hoisting collision.** Four separate, non-overlapping
   `for (int i = 0; i < 10; ++i) { ... }` loops (one in `$MAIN`, two in
   `$ODE`, one in `$TABLE`) each declare `i` inside their own `for`
   statement, but mrgsolve's variable-hoisting preprocessor lifts every
   `int i` into one shared anonymous-namespace scope without
   deduplicating: `error: redefinition of 'int {anonymous}::i'`. Same
   defect class as `kidney-transplant-rejection` #27 and
   `hereditary-spherocytosis`, just with `i` instead of `k`.
2. **Multi-declarator hoisting defect**, surfacing only once (1) is fixed.
   Three lines mix an initializer with a following comma-separated
   declarator — `double MG = 0.0, MGA = 0.0, NSM = 0.0, SMALL = 0.0;`,
   `double NASPF = 0.0, NMIIF = 0.0, CLASP = 0.0, CLRUP = 0.0, DMAXF =
   0.0;` (both `$ODE`), and `double NF11 = 0.0, NF14 = 0.0, NF17 = 0.0,
   MGOUT = 0.0, OVOL = 0.0;` (`$TABLE`). mrgsolve's hoister keeps only the
   first name on each line; gcc reports every later name as undeclared.
   Same defect class as `diabetic-ketoacidosis` #37.

Neither defect is inside any of the five in-scope compounds' own blocks —
both are pre-existing, incidental to the disease-side follicle/scan code.
Per the guide's settled policy, **not fixed in the checked-in original**;
the delivered `cos_mrgsolve_model_refactored.R` carries the syntax-only fix
forward directly:

- The four colliding `i` loop counters renamed to four distinct names
  (`$MAIN`'s threshold-setup loop -> `iTH`, `$ODE`'s follicle-state-gather
  loop -> `iG`, `$ODE`'s main follicle-dynamics loop -> `iF`, `$TABLE`'s
  follicle-count scan -> `iC`), every `[i]` array index inside each loop
  body renamed in lockstep. No loop bound, increment, or body logic
  changed.
- The three multi-declarator lines each split into one
  `double NAME = value;` statement per variable, values unchanged.

Both changes are confirmed non-numeric below (exact-match verification).

## Archetypes

### FSH (recombinant FSH) — Archetype 3 minus peripheral (depot+central)

```
GUT_FSH  -> $ODE:  dxdt_GUT_FSH  = -KA_FSH*GUT_FSH;
CENT_FSH -> $ODE:  dxdt_CENT_FSH =  KA_FSH*GUT_FSH - CL_FSH*CENT_FSH;
C_FSH = CENT_FSH / V1_FSH
```
| Original | Refactored |
|---|---|
| `FSHDEP`/`FSHC` (cmt) | `GUT_FSH`/`CENT_FSH` |
| `KAF`/`KELFX`/`VF` | `KA_FSH`/`CL_FSH`/`V1_FSH` |

No `F_FSH` term — bioavailability (`FB$FSH = 0.80`) is pre-applied to the
dosed amount in the R driver's own `cos_events()`, exactly as the original
does; the refactor keeps this rather than inventing an in-DSL `F_FSH`
multiplier the original never had.

**Why there is no `EFFECT_FSH`:** the original has no scalar FSH effect
term at all. Exogenous FSH is summed with endogenous FSH and
corifollitropin's FSH-equivalent contribution into one pooled variable,
`FSHTOT = FSHE + C_FSH + C_CORI`, which is then read *inside a 10-iteration
per-follicle loop* against a per-slot threshold (`TEFF`) that is itself
dynamic — falling as the follicle enlarges, rising with AMH. There is no
single EC50 to name; the "Hill" is really ten simultaneous Hills against
ten different, state-dependent thresholds. Manufacturing a single
`EFFECT_FSH` would misrepresent this — same "don't force a bad fit"
reasoning as `hypercalcemia-of-malignancy`'s denosumab row (no scalar
effect term exists in the original either). `C_FSH` alone satisfies the
interface requirement (PD reads exactly one named concentration).

### Corifollitropin alfa (CORI) — bespoke, shares FSH's volume

```
GUT_CORI  -> $ODE:  dxdt_GUT_CORI  = -KA_CORI*GUT_CORI;
CENT_CORI -> $ODE:  dxdt_CENT_CORI =  KA_CORI*GUT_CORI - CL_CORI*CENT_CORI;
C_CORI = POT_CORI * CENT_CORI / V1_FSH        // FSH-equivalent, IU/L
```
| Original | Refactored |
|---|---|
| `CORID`/`CORIC` (cmt) | `GUT_CORI`/`CENT_CORI` |
| `KACO`/`KELCO`/`POTCO` | `KA_CORI`/`CL_CORI`/`POT_CORI` |

**Why `C_CORI` reuses `V1_FSH` rather than getting its own `V1_CORI`:**
the original's own pooling term, `CFSHEX = FSHC/VF + POTCO*CORIC/VF`,
divides *both* drugs' central amounts by the same volume, `VF` — the
original never gives corifollitropin an independent volume of
distribution. Archetype 3's worked example assumes each compound gets its
own `V1_<STEM>`; inventing `V1_CORI` here would add a parameter the
original never had, for no structural gain, so the refactor keeps the
original's own shared-volume design and discloses it rather than silently
"completing" it. `POT_CORI` (was `POTCO`) is the potency multiplier that
converts corifollitropin's own central amount into FSH-equivalent units on
that shared volume, unchanged from the original.

**No independent `EFFECT_CORI`**, same reasoning as FSH: corifollitropin
feeds the identical combined `FSHTOT` pool, read by the same per-follicle
dynamic-threshold loop — same "shared concentration, no independent
effect" pattern this fork has already used for acne-vulgaris's ISO/OXO
parent+metabolite pair and achondroplasia's Vosoritide/TransCon-CNP shared
receptor pair.

### Ganirelix (ANT) and Triptorelin (AGO) — Archetype 3 minus peripheral, plus a genuinely shared competitive-occupancy Hill

Both are plain depot+central PK, structurally identical to FSH/HCG above:

| Original | Refactored |
|---|---|
| `ANTD`/`ANTC` (cmt) | `GUT_ANT`/`CENT_ANT` |
| `KAANT`/`KELANT`/`VANT`/`IC50ANT` | `KA_ANT`/`CL_ANT`/`V1_ANT`/`IC50_ANT` |
| `AGOD`/`AGOC` (cmt) | `GUT_AGO`/`CENT_AGO` |
| `KAAGO`/`KELAGO`/`VAGO`/`EC50AGO` | `KA_AGO`/`CL_AGO`/`V1_AGO`/`EC50_AGO` |

No `F_ANT`/`F_AGO` — same "bioavailability pre-applied to the dosed
amount" pattern as FSH.

**The Hill interface is genuinely bespoke, not a plain per-compound
rename.** The original's own GnRH-receptor block computes agonist and
antagonist occupancy *together*, because they competitively bind the same
receptor site:

```c
double XA = CAGO / EC50AGO;
double XN = CANT / IC50ANT;
double OCCAGO = XA / (1.0 + XA + XN);
double OCCANT = XN / (1.0 + XA + XN);
```

This is real, disclosed competitive-binding pharmacology — each
compound's own occupancy formula depends on *both* drugs' concentrations,
because that is what competing for one receptor site means. The guide's
"Multiple drugs, one pathway" rule against collapsing several drugs into
one shared `EFFECT_<STEM>` does not apply here in its usual sense: this
refactor keeps `EFFECT_AGO` and `EFFECT_ANT` as **two separate named
interfaces** (satisfying independent driveability — each has its own name,
its own `$CAPTURE` entry, and could be redirected independently), it is
only that each one's own defining formula legitimately reads both `C_AGO`
and `C_ANT`, because the original's own mechanism does too. Renamed via
`$GLOBAL` macros:

```c
#define XOCC_AGO   (C_AGO / EC50_AGO)
#define XOCC_ANT   (C_ANT / IC50_ANT)
#define EFFECT_AGO (XOCC_AGO / (1.0 + XOCC_AGO + XOCC_ANT))
#define EFFECT_ANT (XOCC_ANT / (1.0 + XOCC_AGO + XOCC_ANT))
```

Exact rename, not a fit — same algebraic shape as the original's
`OCCAGO`/`OCCANT`, gamma implicitly 1 (no exponent in the original's own
competitive-binding formula, so none was added). Downstream, `FANTF`/
`FANTL` (the antagonist's maximal fractional suppression of FSH/LH
release) and `AMPA`/`KGRDOWN` (the agonist's trigger-drive amplitude and
GnRH-receptor down-regulation rate) are kept **unrenamed** — they are
disease-side ceiling/magnitude parameters multiplying `EFFECT_ANT`/
`EFFECT_AGO`, not part of either compound's own PK/Hill role table, the
same "downstream flux parameter, not part of the compound's own block"
treatment `hypercalcemia-of-malignancy` gives `FUR_DIU`/`FUR_NA`.

### hCG — Archetype 3 minus peripheral, no independent effect term

```
GUT_HCG  -> $ODE:  dxdt_GUT_HCG  = -KA_HCG*GUT_HCG;
CENT_HCG -> $ODE:  dxdt_CENT_HCG =  KA_HCG*GUT_HCG - CL_HCG*CENT_HCG;
C_HCG = CENT_HCG / V1_HCG
```
| Original | Refactored |
|---|---|
| `HCGD`/`HCG` (cmt) | `GUT_HCG`/`CENT_HCG` (renamed to avoid the cmt name colliding with the stem `HCG`) |
| `KAHCG`/`KELHCG`/`VHCG` | `KA_HCG`/`CL_HCG`/`V1_HCG` |

**Why there is no `EFFECT_HCG`:** hCG feeds the same pooled ligand
variable as endogenous LH and pregnancy-derived hCG (`LHEQ = LH + C_HCG +
PHCG`), which is then read through **three different downstream kernels**
(luteal/theca support, meiotic commitment, VEGF/permeability) — the
model's own central thesis ("one ligand pool, three response kernels").
There is no private hCG-only effect term to rename; manufacturing one
would misrepresent the shared-ligand mechanism the model exists to
demonstrate. Same "no scalar EFFECT term, purely structural" pattern as
`hypercalcemia-of-malignancy`'s denosumab row. `C_HCG` alone satisfies the
interface requirement.

## `$GLOBAL` macros — normalizing the "duplicate concentration sites" the census flagged

The original computed each of these five compounds' own concentration
**twice**: once as a local `double` inside `$ODE` (`CFSHEX`/`CANT`/`CAGO`/
`CHCG`), and again, independently, inside `$TABLE` (`FSHTOTAL`/`CANTOUT`/
`CHCGOUT`/`CAGOOUT`) — exactly the "duplicate concentration sites" pattern
this file's census rows were classified for. Normalized into one `$GLOBAL`
macro per compound (`C_FSH`, `C_CORI`, `C_ANT`, `C_HCG`, `C_AGO`,
`EFFECT_AGO`, `EFFECT_ANT`), invoked identically wherever `$ODE` or
`$TABLE` needs the value, so there is exactly one textual definition per
quantity rather than two that could silently drift. This also sidesteps
the dose-instant stale-value reporting artifact described in the fork
workflow guide — the same fix `clostridioides-difficile-infection`'s
refactor used, chosen here because this file already has a `$GLOBAL`-macro
precedent to match (`#define HL(x, k, n) ...`).

`$CAPTURE` lists the macro names directly (`C_FSH C_CORI C_ANT C_HCG C_AGO
EFFECT_AGO EFFECT_ANT`), the same idiom `clostridioides-difficile-infection`
uses and confirmed compiling here too — this is not a `double NAME = NAME;`
self-referential declaration (which would fail), just a bare token in
`$CAPTURE` that the preprocessor expands at its own use site. The original's
own `CANTOUT`/`CHCGOUT`/`CAGOOUT` (its duplicate re-derivations) are removed
as redundant now that `C_ANT`/`C_HCG`/`C_AGO` are captured directly under
their own names; `FSHTOTAL`/`LHEQOUT` (already-combined, disease-side pooled
names, not per-compound) are kept, recomputed via the macros instead of a
second independent derivation.

## `$PARAM` vs `$CAPTURE`

Not attempted as `$PARAM = 0` declarations — the read-only-reference
compile failure this produces under mrgsolve 2.0.1 is already documented
in `pagets-disease/pbd_refactor_notes.md` and every refactor in this batch
since. `C_FSH`/`C_CORI`/`C_ANT`/`C_HCG`/`C_AGO`/`EFFECT_AGO`/`EFFECT_ANT`
are `$GLOBAL` macros, listed bare in `$CAPTURE`. Confirmed via `POST
/model_manifest` on the refactored DSL: all seven appear in `outputPaths`
(81 -> 85 total outputs: +7 new `C_FSH/C_CORI/C_ANT/C_HCG/C_AGO/
EFFECT_AGO/EFFECT_ANT`, -3 removed duplicates `CANTOUT`/`CHCGOUT`/
`CAGOOUT`, net +4), and all seventeen
renamed PK parameters (`KA_FSH, CL_FSH, V1_FSH, KA_CORI, CL_CORI,
POT_CORI, KA_ANT, CL_ANT, V1_ANT, IC50_ANT, KA_HCG, CL_HCG, V1_HCG,
KA_AGO, CL_AGO, V1_AGO, EC50_AGO`) appear in `parameters` with their
original numeric defaults unchanged — 142 parameters both sides (pure
renames, no count change).

## Verification

**Method.** Extracted the bare DSL text from both `cos_mrgsolve_model.R`
(`cos_code <- '...'`) and `cos_mrgsolve_model_refactored.R` (same variable
name). Since the untouched original does not compile (see above), an
**in-memory-only** scratch copy with the identical syntax-only fix was used
for the original side of the comparison — the checked-in original was
never modified. Both DSLs run through the local qspserver `mrgsolve_api`
(`http://localhost:8007`), `POST /model_manifest` then `POST
/run_simulation`, requests spaced ~2s apart per the shared-service note.

**Dosing.** `/run_simulation`'s `events` field (named-variable, `method:
"add"`) silently no-ops against this model — doses landed but never
changed any compartment's trajectory (confirmed: max value 0.0 across the
whole run on both sides, a harness-methodology finding, not a model
defect). Switched to the `dosing` field (numeric 1-based `cmt`, mrgsolve's
own `as.ev()` event convention) with compartment indices read from each
side's own `POST /model_manifest` `compartments` list — this produced real
non-zero PK trajectories (e.g. `s04`'s peak `CENT_ANT` = 247.1, `CENT_AGO`
= 68.8) and is the request shape used for the results below. Doses were
computed by the same arithmetic as the original's own `cos_events()`/
`cos_protocol()` R functions (bioavailability fractions, daily/single-dose
schedules, trigger timing) — pure arithmetic, no R or mrgsolve invoked
locally — using each scenario's own documented reference trigger day
(`stim_days`, from the file's own header/README/`cos_reference_output.txt`)
rather than re-running the original's two-pass `cos_trigger_day()`
simulation; since the identical dosing table is applied to both sides, this
does not affect the validity of the original-vs-refactored comparison.

**Scenarios run — 2 of the file's own 22, chosen to exercise all five
in-scope compounds at least once:**

1. **`s04_hi_ago_freeze`** — PCOS patient (`AFC=25, T50=9.0, TONE=1.8,
   AGE=31`), rFSH 150 IU/d days 2-14, ganirelix 0.25 mg/d days 6-16,
   triptorelin 0.2 mg agonist trigger at day 14.4, no luteal support
   (`fresh=FALSE`). Exercises **FSH, ANT, AGO**. `end=19, delta=0.2`
   (96 points before dosing, 121 with dose-time augmentation).
2. **`s10_corifollitropin`** — normal patient (`AFC=12, T50=9.0, TONE=1.0,
   AGE=32`), corifollitropin alfa 150 µg at day 2, rFSH 150 IU/d days 9-13
   (ENGAGE-style day-8-on top-up), ganirelix 0.25 mg/d days 6-15, hCG
   10 000 IU trigger at day 13.4, luteal support. Exercises **FSH, CORI,
   ANT, HCG**. `end=18, delta=0.2` (91 points before augmentation, 108
   with).

Neither window needed shortening for a solver-step-count reason — both ran
to their full requested horizon on both sides without a `maxsteps` error.

**Result: exact match, max abs diff = 0.0, across all 28 shared/renamed
outputs, both scenarios** — compared point-by-point across each scenario's
full time grid (121 and 108 points respectively). The 28 compared
quantities per scenario are (a) 20 disease-side outputs whose names are
untouched by this refactor (`LH, FSHE, E2, P4, OOC, MIIC, GRR, FSHTOTAL,
LHEQOUT, HCT, NF11, NF14, NF17, MGOUT, OVOL, OHSSG, TWOPN, BLAST, EUPL,
CLBR`) and (b) 8 renamed pairs compared old-name-in-original against
new-name-in-refactored (`ANTC->CENT_ANT, HCG->CENT_HCG, AGOC->CENT_AGO,
CANTOUT->C_ANT, CHCGOUT->C_HCG, CAGOOUT->C_AGO, FSHC->CENT_FSH,
CORIC->CENT_CORI`) — genuinely `0.0` for every point in every scenario,
the expected outcome for a pure structural reorganization (renames plus
one `$GLOBAL`-macro normalization, no Hill-fitting anywhere: `AGO`/`ANT`'s
competitive-occupancy formula is an exact algebraic rename of the
original's own `OCCAGO`/`OCCANT`).

**`EFFECT_AGO`/`EFFECT_ANT` spot-check** (not comparable against the
original, which never captured `OCCAGO`/`OCCANT`): at `s04`, t=14.6 (0.2 d
post-agonist-trigger), `C_AGO=2.2936, C_ANT=7.4044` gives
`XOCC_AGO=45.872, XOCC_ANT=16.828`, so `EFFECT_AGO =
45.872/(1+45.872+16.828) = 0.72012` and `EFFECT_ANT = 16.828/63.700 =
0.26418` — matches the API's reported `EFFECT_AGO=0.7201`/
`EFFECT_ANT=0.2642` to within the API's own 4-decimal output rounding
(diff ~2.3e-5, consistent with rounding, not a formula error).

`/model_manifest` on the refactored DSL additionally confirmed all
seventeen renamed PK parameters and the ten new/renamed compartments and
covariates (`GUT_FSH, CENT_FSH, GUT_CORI, CENT_CORI, GUT_ANT, CENT_ANT,
GUT_HCG, CENT_HCG, GUT_AGO, CENT_AGO, C_FSH, C_CORI, C_ANT, C_HCG, C_AGO,
EFFECT_AGO, EFFECT_ANT`) appear correctly (see "`$PARAM` vs `$CAPTURE`"
above).

No scratch/debug artifacts were left in the repository — DSL extraction,
the two-scenario comparison harness, the in-memory-only original-side
build-defect workaround, and all intermediate `.cpp`/JSON files ran
entirely from a session-local scratchpad directory outside the repo tree,
deleted after use.

## R driver consistency

Because five compartments were renamed, every reference to their old names
in this file's own embedded R driver (`cos_events()`) was updated in
lockstep so the shipped R driver still runs correctly end-to-end:
`"FSHDEP"->"GUT_FSH"` (2 call sites), `"CORID"->"GUT_CORI"` (1),
`"ANTD"->"GUT_ANT"` (1), `"HCGD"->"GUT_HCG"` (2),
`"AGOD"->"GUT_AGO"` (1) — confirmed by exact-count grep against the
original before and after. `"CAB"`/`"LET"`/`"P4D"` (out of scope) are
unchanged. The `FB` bioavailability list's own keys (`FSH`, `CORI`, `ANT`,
`HCG`, `AGO`, `HMG`) are plain R list keys, not compartment references —
left unchanged, same as the guide's precedent for internal driver-variable
names that never reach the DSL.

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`: the
`AGO`, `ANT`, `HCG` rows filled in with real identities (triptorelin,
ganirelix, hCG), target/pathway, redirect site, and this refactor's
verification summary; two new rows added for `CORI` (corifollitropin
alfa) and `FSH` (recombinant FSH), the corpus-census undercount the task
flagged as likely. No rows added for `CAB`/`LET`/`P4D` — real compounds
found during the audit, disclosed above, deliberately out of scope for
this refactor.
