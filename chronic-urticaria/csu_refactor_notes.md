# CSU (chronic-urticaria) PK/PD refactor notes

Refactored sibling: `csu_mrgsolve_model_refactored.R` (+ extracted
`csu_mrgsolve_model_refactored.cpp`, byte-identical to the quoted DSL block
in the `.R`, for qspserver `model_content` calls).
Original (never edited): `csu_mrgsolve_model.R`.

## Confirmed compound identities (not the census's generic classifier labels)

The census file listed two of this model's four compounds under generic
classifier labels ("BTK inhibitor", "H1-antihistamine (AH)"). Both resolve
to specific, named real drugs once the code and its own dosing-scenario
labels are read, not left as placeholders:

| Census label | Confirmed identity | Evidence in the code |
|---|---|---|
| H1-antihistamine (AH) | **Cetirizine** | `$PARAM` comment "H1-antihistamine PK (1-compartment, cetirizine prototype)"; the model's own dosing scenario is literally named `"Cetirizine 10 mg QD"` (`scenarios` list, id 2) |
| BTK inhibitor | **Remibrutinib** | `$PARAM` comment "BTK inhibitor PK (1-compartment oral, remibrutinib prototype)"; `dose_BTK` variable and scenario comment both say "remibrutinib prototype" |
| Dupilumab | Dupilumab (already named) | `$PROB`/`$PARAM` comments, LIBERTY-CSU CUPID A/B citation |
| Omalizumab | Omalizumab (already named) | `$PROB`/`$PARAM` comments, GLACIAL/ASTERIA citations |

Neither AH nor BTK is a deliberately generic class-level placeholder with
made-up parameters — the code names a specific real compound for each, so
the refactor keeps the existing code-level stems (`AH`, `BTK`) per the
guide's "don't invent a new stem" rule, while recording the confirmed real
identity here and in the compound census.

## Archetype per compound

**AH (cetirizine) — Archetype 3 minus peripheral** (depot + central,
linear elimination, no tissue distribution). `GUT_AH`/`CENT_AH`, renamed
from `AH_GI`/`AH_plasma`. `EFFECT_AH` renamed from `E_AH`; the original's
effect term was already a Hill ratio (`pow(C_AH,n_AH)/(pow(EC50_AH,n_AH)+
pow(C_AH,n_AH))`), so `EMAX_AH`/`EC50_AH`/`GAMMA_AH` are a straight rename
of `kinh_AH`/`EC50_AH`/`n_AH` — a rename, not a refit.

**OMA (omalizumab) — Archetype 3 (2-cpt SC) PK, plus a bespoke drug-target
binding system.** `GUT_OMA`/`CENT_OMA`/`PERI_OMA`, renamed from
`OMA_depot`/`OMA_c`/`OMA_p`. The binding rate constants `kbind_OMA`/
`kdis_OMA` are renamed to `KON_OMA`/`KOFF_OMA` per the guide's TMDD
on/off-rate naming — but the binding *system itself* is kept bespoke
rather than forced into the full Archetype 4 shape, for two reasons found
by reading the original:

1. `IgE_free` is not a private receptor pool belonging only to omalizumab
   — it is the disease's own free-IgE state, read directly by the
   FcεRI-arming equation (`dxdt_MC_primed`) independent of any drug. Archetype
   4's `REC_FREE_<STEM>` is meant for a compound-private receptor; renaming
   `IgE_free`/`IgE_OMA` to `REC_FREE_OMA`/`COMPLEX_OMA` would suggest a
   drug-private pool that isn't what this is, so the original compartment
   names are kept unchanged (only the PK compartments proper — depot/
   central/peripheral — were renamed).
2. The complex (`IgE_OMA`) has **no clearance term** in the original
   (`dxdt_IgE_OMA = KON_OMA*C_OMA*IgE_free - KOFF_OMA*IgE_OMA` only) and
   `RTOT` is not fixed (IgE has its own independent synthesis/degradation,
   `ksyn_IgE`/`kdeg_IgE`), so the assumptions the Archetype 4 template
   makes (a conserved total receptor pool, a complex that degrades) do not
   hold here. This was tested directly: at true mass-balance steady state,
   `IgE_free_ss = ksyn_IgE/kdeg_IgE` **independent of `C_OMA`**, because
   nothing removes the accumulated complex — so a steady-state Hill fit of
   "occupancy vs. `C_OMA`" (the guide's prescribed approach when an effect
   emerges from ODE-solved kinetics) is degenerate: every fit point
   converges toward the same zero-drug-effect asymptote regardless of
   dose, for a reason unrelated to the drug's actual potency. No Hill fit
   was attempted for this reason; per the guide, a poor/inapplicable fit is
   a finding, not something to force. `EC50_OMA`/`EMAX_OMA`/`GAMMA_OMA` do
   not exist (neither did they in the original), so Omalizumab is not part
   of the auto-discoverable `C_<STEM>`+`EC50_<STEM>` Hill-interface pair —
   `C_OMA` is still exposed and $PARAM-manifest-visible, just without a
   Hill potency parameter, which the guide explicitly allows ("plenty of
   disease models... use [PK] for constructs that don't fit the Hill
   interface... not a bug").

`IgE_suppression` (`(1 - IgE_free/IgE_norm)*100`, already in the
original's `$TABLE`) is kept as Omalizumab's disease-facing readout, since
that already is "one named function of [IgE], not buried inside a combined
multi-drug expression" — it just isn't a Hill function of `C_OMA` directly,
because the drug's action here is a genuine multi-state binding-kinetics
outcome, not a memoryless concentration→effect map.

**DUP (dupilumab) — Archetype 3 (2-cpt SC), linear.**
`GUT_DUP`/`CENT_DUP`/`PERI_DUP`, renamed from `DUP_depot`/`DUP_c`/`DUP_p`.
`EFFECT_DUP_IL4`/`EFFECT_DUP_IL13` renamed from `E_DUP_IL4`/`E_DUP_IL13`;
`GAMMA_DUP=1.0` added (the original had no Hill exponent for either).
**Major finding, preserved verbatim, not fixed** — see below.

**BTK (remibrutinib) — Archetype 3 minus peripheral** (depot + central,
linear). `GUT_BTK`/`CENT_BTK`, renamed from `BTK_GI`/`BTK_plasma`.
`EFFECT_BTK` renamed from `E_BTK`; the original's effect term was a plain
ratio (`kinh_BTK*C_BTK/(EC50_BTK+C_BTK)`, no Hill exponent), so
`GAMMA_BTK=1.0` was added — a rename, not a refit.

## Structural / discoverability implementation notes

- `C_<STEM>` and `EFFECT_<STEM>` are **not** `$PARAM` members. A `$PARAM`
  compiles as a read-only reference inside `$ODE` in this mrgsolve build,
  so a value recomputed every timestep cannot be one — this is the same
  constraint already validated and documented in
  `copd/copd_refactor_notes.md` and `breast-cancer/bc_refactor_notes.md`.
  `EFFECT_<STEM>` values are predeclared (bare, no initializer) in
  `$GLOBAL` and bare-assigned in `$ODE`; `C_<STEM>` values get their sole
  literal `double C_<STEM> = <expr>;` declaration in `$TABLE` (identical
  formula to the bare `$ODE` recompute), matching the pattern already
  validated in `copd_mrgsolve_model_refactored.R` — mrgsolve does not carry
  an `$ODE`-scoped local into `$TABLE`, so this is what actually makes
  `C_<STEM>` available for reporting at every `$TABLE`-evaluated row, and
  it is also literal-text-discoverable by downstream tooling. This is the
  actual fix for the "duplicate concentration site" pattern the compound
  census flagged for all four compounds: the original computed the same
  concentration twice under **two different names** (e.g. `C_AH` in
  `$ODE`, `CONC_AH` in `$TABLE`); the refactor uses one canonical name,
  `C_<STEM>`, at both sites.
- `EC50_<STEM>` exists for AH, DUP, and BTK (so all three are discoverable
  by the `double C_<STEM> = ...;` + `EC50_<STEM>` pattern); it does not
  exist for OMA (see Archetype/bespoke discussion above).
- `$CAPTURE @annotated` lists only genuinely-derived `$TABLE`/effect
  quantities, never a `$CMT` compartment name — the original's own
  attempt to do the latter (via bare `capture NAME NAME NAME` lines) is
  exactly what made it fail to compile; see the upstream defect below.

## Upstream defect found: original does not compile at all (UPSTREAM_ISSUES.md #147)

The original's `$TABLE` block ends with four lines of the form `capture
CONC_AH CONC_OMA CONC_DUP CONC_BTK` (multiple space-separated names on one
`capture` statement, no `$CAPTURE` block header anywhere in the file).
Confirmed via the qspserver `mrgsolve_api`, posting the **untouched**
original's own extracted DSL block: this fails to compile —

```
490:19: error: expected initializer before 'CONC_OMA'
```

— because this build's `capture` directive takes exactly one name per
statement, not a space-separated list. Wrapping the same fifteen names in
a proper `$CAPTURE` block surfaces a second, independent defect: seven of
them (`IgE_free`, `MC_primed`, `MC_act`, `Hist_skin`, `Hist_plasm`,
`IL31_skin`, `IL33_skin`) are already `$CMT` compartments, and mrgsolve
rejects a compartment name inside `$CAPTURE`:

```
Error in validObject(.Object) :
  invalid class "mrgmod" object: compartment should not be in $CAPTURE: ...
```

**Fix applied in the refactored file only** (never in the original):
a proper `$CAPTURE @annotated` block, listing only the eight genuinely-
derived `$TABLE` quantities plus the new `EFFECT_*`/`C_*` names — no
`$CMT` compartment re-listed. This is a syntax-only, non-numeric fix (it
only changes which output block a name lives under, not any equation) —
confirmed by the verification below matching the original bit-for-bit.
Logged as `translations/UPSTREAM_ISSUES.md` entry #147.

## Major finding: Dupilumab has zero simulated effect (UPSTREAM_ISSUES.md #148)

`E_DUP_IL4`/`E_DUP_IL13` (renamed `EFFECT_DUP_IL4`/`EFFECT_DUP_IL13`) are
computed in `$ODE` exactly as the original computed them, but — as in the
original — **neither is referenced by any `dxdt_` equation**. Dupilumab's
own PK reaches a real, non-trivial concentration (`CONC_DUP`/`C_DUP` peaks
at ~81 mg/L under the file's own "Dupilumab 300 mg q2wk" scenario), but
that concentration never reaches the disease system (`MC_primed`,
`MC_act`, `IL31_skin`, `IL33_skin` all read only `EFFECT_AH`/`EFFECT_BTK`
or `MC_act` itself).

Confirmed empirically, not just by reading the code: running the
original's own "Dupilumab 300 mg q2wk" scenario and its own "No treatment"
scenario side by side and comparing `UAS7` (and every other disease-state
output) at every shared time point gives a **maximum absolute difference
of 0.0** — bit-for-bit identical trajectories with and without an
81 mg/L dose of the drug on board, despite the model's own `$PROB` header
and calibration comments citing LIBERTY-CSU CUPID A/B (Simpson 2023 NEJM)
as the basis for this arm.

**Not fixed** — wiring `EFFECT_DUP_IL4`/`EFFECT_DUP_IL13` into the disease
ODEs would be new pharmacology the original never had, which is out of
scope for a rename-and-reorganize refactor. Preserved verbatim in the
refactored file (computed, unused, exactly as in the original) and logged
as `translations/UPSTREAM_ISSUES.md` entry #148.

## Finding: antihistamine-dosed scenarios blow up to NaN (UPSTREAM_ISSUES.md #149)

Three of the file's own seven scenarios (`"Cetirizine 10 mg QD"`,
`"High-dose AH 40 mg/day"`, `"Omalizumab 300 mg + AH"`) run cleanly for a
while and then every captured output simultaneously becomes `NaN` for the
rest of the 24-week horizon (onset at t=936h/2280h/1152h respectively —
the higher, longer-dosed scenario fails *later*). The remaining four
scenarios (no treatment, omalizumab alone, dupilumab, BTKi) run to
completion with no `NaN` at all. Most likely mechanism (circumstantial,
not a line-by-line solver trace): `EFFECT_AH` uses `pow(C_AH, GAMMA_AH)`
with `GAMMA_AH=1.5`, a fractional exponent; `pow()` of a negative base to
a non-integer power is `NaN`, and `C_AH` (unfloored, first-order decay) can
undershoot fractionally below zero during a long, quiet post-dose decay
tail under an adaptive-step solver. AH is the only compound with a
fractional (non-1, non-implicit) Hill exponent among the four, and the
only scenario family that never hits this (omalizumab-only) never touches
`C_AH` at all (stays exactly 0). Reproduced **identically** (same onset
time index, same values before onset) on the original and the refactored
file — confirming this is a pre-existing property of the original's own
ODE system, unrelated to the refactor. Preserved verbatim, not fixed;
logged as `translations/UPSTREAM_ISSUES.md` entry #149.

## Finding: `IgE_free`'s initial condition is not at the model's own equilibrium (UPSTREAM_ISSUES.md #150)

`IgE_free_0 = IgE0 = 300` nM, but the model's own no-drug dynamics
(`dxdt_IgE_free = ksyn_IgE - kdeg_IgE*IgE_free`) equilibrate at
`ksyn_IgE/kdeg_IgE = 0.3` nM — three orders of magnitude lower. So even
the file's own "No treatment" scenario shows a large, purely numerical
IgE decay (300 to 266.1 to 236.1 nM over the first two weeks, continuing
toward ~0.3 nM), unrelated to any drug. Reproduced identically in the
refactor (initial conditions and parameter values are unchanged by a
PK-reorganization refactor). Not fixed; logged as
`translations/UPSTREAM_ISSUES.md` entry #150.

## Minor finding (not separately logged upstream): several disease-network parameters are declared but never used

`ksyn_IL4`, `kdeg_IL4`, `ksyn_IL13`, `kdeg_IL13`, `keo_in`, `keo_out`,
`Eo0`, `kUAS_H`, `kUAS_IL31` are all declared in the original's `$PARAM`
but never appear in any `$ODE`/`$TABLE` equation (no `IL4`/`IL13`/`Eo`
compartments exist at all, and the `UAS7` formula uses hardcoded
coefficients `0.5`/`0.3`/`0.2` rather than `kUAS_H`/`kUAS_IL31`) — likely
vestigial from an earlier model iteration. Carried over unchanged in the
refactor (removing an unused parameter is a cleanup decision out of scope
for this refactor, not a rename), each annotated `-- unused, see notes` in
`$PARAM`.

## Verification

All seven of the original's own dosing scenarios (`scenarios` list in the
`.R` driver code: No treatment, Cetirizine 10 mg QD, High-dose AH 40
mg/day, Omalizumab 300 mg q4wk, Omalizumab 300 mg + AH, Dupilumab 300 mg
q2wk, BTKi 25 mg QD) were run through the qspserver `mrgsolve_api`
(`/model_manifest` then `/run_simulation`, `POST http://localhost:8007`),
identical dosing (same amounts, `ii`/`addl`, and 1-based compartment
indices — `$CMT` declaration order is unchanged between the original and
the refactored file, so dosing into compartment N hits the same physical
entity in both), full `end=4032h, delta=24h` time grid, for both the
original (via the syntax-only fix in entry #147 above — the untouched
original does not compile on its own) and the refactored file.

Every `$CAPTURE`d output was compared at every shared time point
(`CONC_AH`/`C_AH`, `CONC_OMA`/`C_OMA`, `CONC_DUP`/`C_DUP`,
`CONC_BTK`/`C_BTK`, `IgE_free`, `IgE_suppression`, `MC_primed`, `MC_act`,
`Hist_skin`, `Hist_plasm`, `IL31_skin`, `IL33_skin`, `UAS7`, `WCU`, `CR`).

**Result: exact match, all 7 scenarios, all 15 compared outputs — max
absolute difference 0.0, max relative difference 0.0.** This includes the
three scenarios that hit the pre-existing NaN blowup (entry #149): the
NaN onset happens at the identical time index in both files, and every
numeric value before onset matches bit-for-bit. This is the expected
result for archetypes 1-3 and the bespoke OMA structure (pure
reorganization/renaming, no Hill-fitting) per the guide's tolerance
section — nothing beyond floating-point-scale deviation was found.

`model_manifest` was also called once on the finished refactored file to
confirm every renamed parameter and covariate is declared correctly (63
parameters, 30 output paths including all four `C_<STEM>` and four
`EFFECT_<STEM>` values, all 18 `$CMT` compartments in the original
declaration order). `Rscript -e 'parse("csu_mrgsolve_model_refactored.R")'`
succeeds with no error. The extracted `.cpp` was confirmed byte-identical
to the quoted string body in the `.R` file.

## qspserver infrastructure note (not a defect in either model)

During verification, `/model_manifest` and `/run_simulation` calls to
`http://localhost:8007` intermittently failed for **several minutes**,
for both the original and the refactored model, and even for a
deliberately invalid `model_content` payload, with:

```
Error in readRDS(cache_file) :
  embedded nul in string: '//  h^-1\0\004\0\to 0'
```

— a corrupted server-side compile cache, unrelated to this model's
content (confirmed by the failure persisting across different
`model_content` strings and across a `/model_manifest` call with garbage
input). Recovered on its own after retrying with a delay; all scenarios
above were eventually completed successfully once it cleared. Consistent
with the task brief's warning that this shared API "has crashed under
heavier load before" from concurrent agents.
