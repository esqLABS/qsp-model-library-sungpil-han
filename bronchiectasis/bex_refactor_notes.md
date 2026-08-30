# Refactor notes — `bronchiectasis/bex_mrgsolve_model.R`

Scope: all **three** compounds this file gives their own PK compartments and
that have existing rows in `driver-patches/data/compound_perturbation_census.md`
— **Azithromycin (AZM)**, **Ciprofloxacin (CIPRO)**, and **Tobramycin
(TOBRA)** — all classified "Redirect concentration (clean single site)".
Dornase Alfa (the file's fourth pharmacological input, an inhaled mucoactive
with no census row) is completely untouched: its PK compartment
(`LUNG_DNase`), parameters (`Kd_DNase`, `Emax_DNase`, `EC50_DNase`), and
effect (`E_DNase_MCC`) are byte-identical to the original everywhere they
appear.

## Headline finding: the original does not compile at all — unrelated to any of the three compounds' own PK

Before any renaming could even be attempted, `POST /model_manifest` on the
untouched original (via the qspserver `mrgsolve_api` container,
`http://localhost:8007`) failed to build. Root cause, logged as
`translations/UPSTREAM_ISSUES.md` #46: the `$TABLE` block ends in three bare
lowercase `capture ...` lines with **no `$CAPTURE` (or any) block header at
all** — not valid mrgsolve DSL, so mrgsolve tries to compile `capture(...)`
as a C++ function call and fails outright. Once patched with a real
`$CAPTURE` header, a second, already-familiar defect (same class as issues
#30/#34/#41/#45) is revealed: the first `$CAPTURE` line also lists six names
(`IL8, NE, MCC, AD, BIOFILM, EXAC`) that duplicate `$CMT` compartment names,
which mrgsolve 2.0.1 also refuses to build.

Per the fork guide's settled policy for this situation
(`FORK_WORKFLOW_GUIDE.md`, "When the original doesn't compile at all"), both
fixes are **syntax-only, non-numeric, non-behavioral**, and — because a
"verified" deliverable nobody can actually run defeats the point of this
refactor — both are applied directly in the delivered `_refactored.R`, not
just a throwaway scratch copy:

1. Each bare `capture ...` line becomes its own `$CAPTURE` header line
   (`capture` → `$CAPTURE`; no symbol added, removed, or reordered).
2. The six compartment names are dropped from the first `$CAPTURE` line
   (`$CAPTURE FEV1 FEV1_pct log10_BACT`) — compartments are always present
   in the output regardless of `$CAPTURE`, so this changes nothing about
   what is reported.

Neither change touches a single equation or parameter value. The tracked
`bex_mrgsolve_model.R` is completely untouched and still carries both
defects exactly as written; the workaround lives only in
`bex_mrgsolve_model_refactored.R` and is called out again, inline, at the
`$CAPTURE` block itself. See `translations/UPSTREAM_ISSUES.md` #46 for the
full defect writeup (also notes, in passing and not fixed, a pre-existing
day/hour unit-labeling mismatch in the three refactored compounds' PK rate
constants, the same defect class as `abdominal-aortic-aneurysm` issue #44 —
it does not visibly break any scenario within the windows verified here, so
it is not the subject of this refactor).

## Archetype per compound

**Azithromycin (AZM): Archetype 3 (depot + central + peripheral, linear).**
The original already used `GUT_AZM`/`CENT_AZM`/`PERI_AZM` compartment names
matching this fork's convention almost exactly — only the parameter names
needed renaming (`Ka_AZM`→`KA_AZM`, `Vd_AZM`→`V1_AZM`, `Vp_AZM`→`V2_AZM`,
`CL_AZM`/`F_AZM`/`Q_AZM` unchanged). The ODEs are a pure identifier rename;
`(CL_AZM/V1_AZM) + (Q_AZM/V1_AZM)` is arithmetically the combined
`(CL_AZM+Q_AZM)/V1_AZM` the guide's template shows, not a behavioral change.

**Tobramycin (TOBRA): bespoke, 2-compartment, one-way, no back-flux** — does
not fit archetype 2 (which is bidirectional `Q/V1`/`Q/V2` exchange) or
archetype 3 (which is depot→central with reversible peripheral exchange).
The original's `LUNG_Tobra` is dosed directly (inhaled), has its **own**
local elimination, *and* leaks a fraction of its content one-way into
`CENT_Tobra` (systemic), which never flows back. This is the same class of
non-standard structure the `mucopolysaccharidosis-type-1` (laronidase) and
`age-related-macular-degeneration` refactors already handled as bespoke —
renamed to the convention, isolated from PD, single exposed concentration,
but not forced into a bidirectional-exchange shape. See "Naming" below for
which compartment plays which (non-table) role.

**Ciprofloxacin (CIPRO): Archetype 3 without a peripheral compartment**
(depot + central only) — explicitly one of the guide's own named variants
("Drop `GUT_TCZ`/`KA_TCZ`/`F_TCZ` for a 2-compartment-no-depot variant" runs
the same idea in reverse: keep the depot, drop the peripheral). Pure
identifier rename, no structural change.

## Naming

| Role | AZM | TOBRA | CIPRO |
|---|---|---|---|
| Depot (dosed) | `GUT_AZM` (unchanged) | `LUNG_TOBRA` (was `LUNG_Tobra`) — bespoke role: dosed, local, **PD-exposed** depot, not a plain absorption `GUT_` | `GUT_CIPRO` (was `GUT_Cipro`) |
| Central / systemic | `CENT_AZM` (unchanged) | `CENT_TOBRA` (was `CENT_Tobra`) — bespoke role: one-way systemic sink, **informational only**, never read by PD | `CENT_CIPRO` (was `CENT_Cipro`) |
| Peripheral | `PERI_AZM` (unchanged) | — (none in the original) | — (none in the original) |
| `KA_/F_/CL_/V1_/V2_/Q_` | `KA_AZM`(was `Ka_AZM`), `F_AZM`(unchanged), `CL_AZM`(unchanged), `V1_AZM`(was `Vd_AZM`), `V2_AZM`(was `Vp_AZM`), `Q_AZM`(unchanged) | `F_TOBRA`(was `F_Tobra_inh`), `V1_TOBRA`(was `Vd_Tobra`), `CL_TOBRA`(was `CL_Tobra`) — systemic-sink params only | `KA_CIPRO`(was `Ka_Cipro`), `F_CIPRO`(was `F_Cipro`), `V1_CIPRO`(was `Vd_Cipro`), `CL_CIPRO`(was `CL_Cipro`) |
| Local/bespoke rate constant | — | `KE_TOBRA_LUNG` (was `ke_Tobra_lung`) — local lung elimination rate; no table slot fits a dosed depot with its own elimination, so this stays a bespoke, stem-anchored name | — |
| ELF penetration ratio | — | — | `ELF_RATIO_CIPRO` (was `ELF_ratio`) |
| Exposed concentration | `C_AZM` (was the value carried by `C_AZM_lung`) | `C_TOBRA` (was `C_Tobra_lung`) | `C_CIPRO` (was `C_Cipro_ELF`) |
| Informational-only concentration (not exposed, not read by PD) | `C_AZM_PLASMA` (was `C_AZM_plasma`) | `C_TOBRA_SYS` (was `C_Tobra_sys` — dead in the original, still dead here) | `C_CIPRO_PLASMA` (was `C_Cipro_plasma`) |
| Hill: Emax/EC50/Gamma | `EMAX_AZM_IL8`/`EC50_AZM_IL8`/`GAMMA_AZM_IL8`, `EMAX_AZM_QS`/`EC50_AZM_QS`/`GAMMA_AZM_QS` | `EMAX_TOBRA`/`EC50_TOBRA`/`GAMMA_TOBRA` | `EMAX_CIPRO`/`EC50_CIPRO`/`GAMMA_CIPRO` |
| Effect(s) on disease | `EFFECT_AZM_IL8` (was `E_AZM_IL8`), `EFFECT_AZM_QS` (was `E_AZM_QS`) | `EFFECT_TOBRA` (was `E_Tobra_kill`) | `EFFECT_CIPRO` (was `E_Cipro_kill`) |
| Decorative, unused elsewhere | — | `MIC_PA_TOBRA` (was `MIC_PA_Tobra`) | `MIC_PA_CIPRO` (was `MIC_PA_Cipro`) |

All parameter *values* are copied verbatim from the original — nothing
invented, nothing defaulted, nothing dropped.

### AZM has two independent downstream effects, kept separate

Per the guide's "never collapse several drugs into one shared Hill term"
rule (which also applies within one compound with more than one target):
`EFFECT_AZM_IL8` (IL-8 suppression, feeds `IL8_prod`) and `EFFECT_AZM_QS`
(quorum-sensing inhibition, feeds `Biofilm_growth`) are two separate named
effects, each with its own `EMAX_AZM_*`/`EC50_AZM_*`/`GAMMA_AZM_*` triplet,
exactly mirroring how the original kept `E_AZM_IL8`/`E_AZM_QS` separate.
Both read the same exposed `C_AZM`.

### A pre-existing internal inconsistency in AZM's own concentration math (found, not fixed)

The original's `C_AZM_plasma = CENT_AZM / (Vd_AZM * 70)` divides by
`Vd_AZM*70` (i.e. treats `Vd_AZM=31` as an L/kg value scaled to a 70 kg
patient, ≈2170 L), while the **elimination/inter-compartmental terms in the
ODE** (`CL_AZM/Vd_AZM`, `Q_AZM/Vd_AZM`) use `Vd_AZM` alone (31, unscaled).
So the same parameter is used two different ways ~70-fold apart within the
same compound's own block. This looks like a genuine authoring
inconsistency, not something introduced by the refactor: it is reproduced
exactly, unchanged, in `C_AZM_PLASMA`/`V1_AZM` here (only renamed), and the
verification below confirms the refactored model's `C_AZM`/`AZM_lung_conc`
values match the original's `C_AZM_lung` values exactly, meaning this
inconsistency (whatever its clinical realism) survived the rename bit-for-
bit. Not logged as its own `UPSTREAM_ISSUES.md` entry since it does not
block building or running (unlike the `$CAPTURE` defect above) — noted here
for a future reviewer instead. Ciprofloxacin's equivalent math
(`C_Cipro_plasma = CENT_Cipro/(Vd_Cipro*70)`, and the ODE's own elimination
term `CL_Cipro/(Vd_Cipro*70)`) is internally consistent by contrast — both
use the same `*70` scaling — so this is specific to AZM, not a repo-wide
pattern in this file.

### Tobramycin's `CENT_TOBRA` (was `CENT_Tobra`) is dead output, preserved as such

`C_Tobra_sys` (informational plasma-equivalent concentration, now
`C_TOBRA_SYS`) is computed in the original's `$MAIN` but never read by any
`$ODE`/`$TABLE` expression, never captured, and never plotted — the
compartment `CENT_Tobra`/`CENT_TOBRA` itself only accumulates and clears.
Preserved exactly as dead/unused, per the same "don't invent, don't drop
the original's own quantities" treatment given to unused parameters in
other refactors (e.g. `RTOT_TCZ` in the rheumatoid-arthritis precedent).

## The Hill interface: two renames-only, one already-fitted power law

- **AZM** (`EFFECT_AZM_IL8`, `EFFECT_AZM_QS`): the original was already a
  plain `Emax*C/(EC50+C)` ratio for both effects, no Hill exponent. Rename
  only — `GAMMA_AZM_IL8 = 1` and `GAMMA_AZM_QS = 1` are added as explicit
  named parameters (none existed before), and `pow(x,1)==x` makes the new
  formula arithmetically identical to the original ratio.
- **CIPRO** (`EFFECT_CIPRO`): same situation — the original's
  `Emax_Cipro * C_Cipro_ELF / (EC50_Cipro + C_Cipro_ELF)` had no exponent.
  `GAMMA_CIPRO = 1` added, rename only.
- **TOBRA** (`EFFECT_TOBRA`): the original was **already** exactly the
  guide's Hill formula, `Emax_Tobra * pow(C, Hill_Tobra) / (pow(EC50_Tobra,
  Hill_Tobra) + pow(C, Hill_Tobra))`, with an explicit, previously-fitted
  Hill coefficient (`Hill_Tobra = 1.5`). This is a pure rename
  (`Hill_Tobra` → `GAMMA_TOBRA`, value unchanged) — no refitting was done or
  needed anywhere in this file.

The combined bactericidal term (`E_bact_kill = 1 - (1-EFFECT_TOBRA*TIP_on)
* (1-EFFECT_CIPRO*Cipro_on)`, multiplicative independence between the two
antibiotics) is untouched except for the two renamed inputs — this already
matched the guide's "combine only at the point disease equations actually
use them" rule before the refactor.

## qspserver compatibility

Followed the precedent set by the 21 refactors before this one (see e.g.
`abdominal-aortic-aneurysm/aaa_refactor_notes.md`): `C_<STEM>` and
`EFFECT_<STEM>` are plain `double` locals computed in `$MAIN`/`$ODE`
(assigning to a `$PARAM`-declared symbol from those blocks is a compile
error in this mrgsolve build), and discoverability is satisfied instead by
listing every one of them in `$CAPTURE`. Confirmed via `POST
/model_manifest` on the qspserver `mrgsolve_api` container that:
- `outputPaths` includes all three exposed concentrations (`C_AZM`,
  `C_TOBRA`, `C_CIPRO`), all four effect terms (`EFFECT_AZM_IL8`,
  `EFFECT_AZM_QS`, `EFFECT_TOBRA`, `EFFECT_CIPRO`), and the original's own
  output names (`AZM_lung_conc`, `Tobra_lung_conc`, `Cipro_ELF_conc`,
  `FEV1`, `FEV1_pct`, `log10_BACT`, ...), kept for backward output
  comparability.
- `parameters` includes every renamed Hill parameter
  (`EMAX_AZM_IL8`/`EC50_AZM_IL8`/`GAMMA_AZM_IL8`, `EMAX_AZM_QS`/
  `EC50_AZM_QS`/`GAMMA_AZM_QS`, `EMAX_TOBRA`/`EC50_TOBRA`/`GAMMA_TOBRA`,
  `EMAX_CIPRO`/`EC50_CIPRO`/`GAMMA_CIPRO`) and every renamed PK parameter
  for all three compounds.
- `compartments` shows the renamed 1-based compartment order is unchanged
  from the original (`GUT_AZM, CENT_AZM, PERI_AZM, LUNG_TOBRA, CENT_TOBRA,
  GUT_CIPRO, CENT_CIPRO, LUNG_DNase, BACT, ...`), so dosing by compartment
  number (as `POST /run_simulation`'s `dosing` field uses) is unaffected.

## Verification

Per the guide's mandatory protocol: ran the original file's own dosing
regimens (amount, interval, target compartment) through both the
`$CAPTURE`-patched original and the refactored model via the qspserver
`mrgsolve_api` service (`POST /run_simulation`), comparing every shared
output point-by-point. The API container showed intermittent, non-
deterministic crashes ("worker killed by signal 6", `double free or
corruption`) on **both** models under repeated load — confirmed flaky by
re-running the identical request several times and getting a mix of
successes and failures on each side independently; this is a container/
solver-runtime instability, not a property of either model, and every
comparison below is from a run where both sides completed successfully.

**AZM — Scenario 2's own regimen (250 mg, 3×/week, `ii=7/3` days), 60-day
window** (shortened from the original R script's 2-year horizon per the
guide's solver-budget allowance — both models hit intermittent solver
stress near t≈40–58d under the API's default `maxsteps`, unrelated to
which model is "right"; 60 days was comfortably reproducible on repeated
runs and is described honestly here as shortened): `GUT_AZM`, `CENT_AZM`,
`PERI_AZM`, `FEV1_pct`, `log10_BACT`, `IL8`, `BIOFILM`, `AZM_lung_conc` —
**max abs diff = 0 for every output, every time point (32 points).**

**TOBRA — Scenario 3's own regimen (300 mg BID, 28-day-on/28-day-off
cycle), first 60 days of the cycle** (explicit dose times reproduced from
the original R script's own `on_days`/`t_tip_on` construction, truncated to
the first on/off cycle rather than the full 2-year horizon): `LUNG_Tobra`↔
`LUNG_TOBRA`, `CENT_Tobra`↔`CENT_TOBRA`, `FEV1_pct`, `log10_BACT`,
`Tobra_lung_conc` — **max abs diff = 0 for every output, every time point
(31 points).**

**CIPRO — Scenario 6's own regimen (exacerbation initial condition: `BACT=
2.0, BIOFILM=0.50, NEUT=12.0, IL8=60, NE=6.0, MCC=0.25, AD=0.60, EXAC=
0.85`, reproduced via `replace`-method events at t=0 since the API's
`SimRequest` has no direct initial-condition override field; then
Ciprofloxacin 500 mg BID × 14 days)**: `GUT_Cipro`↔`GUT_CIPRO`,
`CENT_Cipro`↔`CENT_CIPRO`, `FEV1_pct`, `log10_BACT`, `Cipro_ELF_conc`,
`EXAC` — **max abs diff = 0 for every output, every time point (30
points).**

**Sanity check — no-treatment natural history, 60 days, default initial
conditions:** `BACT`, `FEV1_pct`, `log10_BACT`, `IL8`, `BIOFILM`, `MCC`,
`AD` — **max abs diff = 0**, confirming the disease-PD core and the
untouched DNase block are unaffected by the refactor.

**Result: exact match (max abs diff = 0) for every output, every scenario,
every compound.** This is the expected outcome for archetype 3 (AZM,
CIPRO) and the bespoke one-way structure (TOBRA) — pure PK reorganization
plus GAMMA=1/GAMMA=1.5 Hill renames, with no fitting anywhere in this file.

## Anything else worth flagging

- The three compounds' PK rate constants are commented as h⁻¹ (e.g.
  `Ka_AZM = 0.5, // Absorption rate (h-1)`) while every dosing interval and
  `tgrid` in the R script runs on a day-valued clock — the same day/hour
  labeling mismatch class as `abdominal-aortic-aneurysm` issue #44, logged
  in `UPSTREAM_ISSUES.md` #46 as an observation. Unlike #44, none of the
  scenarios verified here (natural history, or each compound's own regimen
  out to 60 days) diverges or blows up because of it, so it is not treated
  as a build-blocking defect — just flagged for whoever eventually
  recalibrates this file's absorption/elimination half-lives against the
  cited clinical trials.
- Compartment ordering is unchanged (`GUT_AZM` is still compartment 1,
  `LUNG_TOBRA` (was `LUNG_Tobra`) is still compartment 4, `GUT_CIPRO` (was
  `GUT_Cipro`) is still compartment 6, etc.), so 1-based compartment-number
  dosing is unaffected between the original and the refactored file.
- Dornase Alfa (`LUNG_DNase`, `Kd_DNase`, `Emax_DNase`, `EC50_DNase`,
  `E_DNase_MCC`, `DNase_on`) is completely untouched — no census row, not
  one of the three scoped compounds.
- Outside the DSL block, the surrounding R script needed exactly two
  string updates for the renamed compartments it doses into by name:
  `cmt = "LUNG_Tobra"` → `cmt = "LUNG_TOBRA"` (in `make_events()` and
  Scenarios 3/4) and `cmt = "GUT_Cipro"` → `cmt = "GUT_CIPRO"` (in
  `make_events()`). `GUT_AZM` needed no change (already matched the
  convention). No other line of the R wrapper (scenarios, plotting,
  summary table) differs from the original.

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`, the
`bronchiectasis | AZM`, `bronchiectasis | CIPRO`, and `bronchiectasis |
TOBRA` rows.
