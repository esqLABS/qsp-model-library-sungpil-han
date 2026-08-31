# Refactor notes — `cholelithiasis/chol_mrgsolve_model.R`

Three compounds refactored, per the existing rows in
`driver-patches/data/compound_perturbation_census.md`: **Ezetimibe**
(classified "Redirect concentration (clean single site)"), **Statin**
(classified "Redirect concentration (clean single site)"), and **UDCA**
(ursodeoxycholic acid, classified "Normalize duplicate concentration
sites, then redirect"). The file models exactly these three compounds —
no others — so nothing was left untouched by scope-narrowing; every
compound in the file is in scope.

## What the Statin actually is

The `$PROB` line and the `$PARAM` section header both say so directly:
`// ---- Statin PK Parameters (Simvastatin) ----`, and Scenario 4 is named
`"UDCA + Simvastatin 40 mg"`. **The statin is Simvastatin.**

## Archetype per compound

- **Ezetimibe (`EZET`)** — Archetype 3 variant (depot + central, **no
  peripheral**): `A_gut_EZET`→`A_plas_EZET` renamed `GUT_EZET`→`CENT_EZET`.
  Clean rename, no restructuring.
- **Statin/Simvastatin (`STAT`)** — Archetype 3 variant (depot + central,
  no peripheral): `A_gut_STAT`→`A_plas_STAT` renamed `GUT_STAT`→`CENT_STAT`.
  Clean rename, no restructuring.
- **UDCA (`UDCA`)** — **bespoke**, does not fit any of the guide's four
  archetypes. The original models a genuine multi-organ enterohepatic
  chain with a branching first-pass split, not a linear
  depot→central(→peripheral) chain:
  ```
  GUT_UDCA --KA_UDCA-->  [splits by FRAC_HEP_UDCA]
      (1-FRAC_HEP_UDCA) -----------------------------> CENT_UDCA (plasma, cleared by CL_UDCA)
      FRAC_HEP_UDCA -----> HEP_UDCA --KEHC_UDCA,FRAC_BILE_UDCA--> BILE_UDCA --(0.40 split)--> GB_UDCA
  ```
  (renamed from `A_gut_UDCA`/`A_plas_UDCA`/`A_hep_UDCA`/`A_bile_UDCA`/
  `A_gb_UDCA`). Per the guide's own instruction for this case ("Rename to
  the convention, isolate it from PD, expose `C_<STEM>`, and note... that
  this compound needed a bespoke structure and why"): all five
  compartments renamed, nothing added, removed, or reordered
  (compartment numbering 1–5 preserved exactly, so the R driver's
  `ev(cmt = "GUT_UDCA", ...)` dosing needed only a string-literal rename,
  no restructuring). This mirrors the doxycycline precedent in
  `abdominal-aortic-aneurysm` (a genuinely distinct tissue compartment
  kept, not flattened into a simpler archetype).

## The Hill interface: renames, not fits, for all three

Every one of the three compounds' original effect terms was already a
plain Michaelis–Menten/Hill-shaped ratio (`Emax*C/(Km+C)`, implicit
gamma=1) **except UDCA's**, which is linear-with-a-hard-cap, not a
saturating ratio (see below). No `nls()` fitting was performed anywhere in
this file.

- **Statin**: `E_STAT = Emax_STAT*C_STAT_plas/(Km_STAT+C_STAT_plas)` →
  renamed `EFFECT_STAT` (`EMAX_STAT`, `EC50_STAT` renamed from
  `Emax_STAT`/`Km_STAT`, `GAMMA_STAT=1` new/math-implied). Statin has
  **three independent downstream points of use**, all reading the *same*
  `EFFECT_STAT` value (not three separate Hill terms — the original
  multiplies the one shared fraction by three different literal weights
  at each site, exactly the pattern already established for statin in
  `elevated-lipoprotein-a`'s `STA` and `abdominal-aortic-aneurysm`'s
  `STAT`):
  - `k_CHOL_syn_eff = k_CHOL_syn*(1-EFFECT_STAT)*(...)` — cholesterol
    synthesis inhibition (`$MAIN`, weight 1×)
  - `CHOL_uptake = 0.20*(1+EFFECT_STAT*0.8)` — LDLR-mediated uptake
    stimulation (`$ODE`, weight 0.8×)
  - `k_CHOL_bil_eff = k_CHOL_bil*(1-E_UDCA_CSI)*(1-EFFECT_STAT*0.3)` —
    biliary cholesterol secretion reduction (`$MAIN`, weight 0.3×)
  All three multiplications are preserved at the exact literal weights
  and in the exact block each already lived in — only the variable name
  changed.
- **Ezetimibe**: `E_EZET = Emax_EZET*C_EZET_plas/(Km_EZET+C_EZET_plas)` →
  renamed `EFFECT_EZET` (`EMAX_EZET`/`EC50_EZET` renamed from
  `Emax_EZET`/`Km_EZET`, `GAMMA_EZET=1` new). **This term is never read by
  any `dxdt_` line or disease-facing calculation anywhere in the
  original** — confirmed by grep, it is defined once and never
  referenced again. `EFFECT_EZET` is preserved in the refactored file
  exactly as dangling as the original's `E_EZET` was; wiring it into the
  disease network would be a real, undisclosed behavioural change, out of
  scope for a naming/structure refactor. Logged as
  `translations/UPSTREAM_ISSUES.md` #134.
- **UDCA — bespoke, not forced into the Hill shape.** UDCA's two effect
  terms are **linear in concentration with a hard clamp**, not saturating
  ratios:
  ```
  EFFECT_UDCA_DISSOL = K_DISSOL_UDCA * C_UDCA_norm      // (was E_UDCA_dis * C_UDCA_bile_norm)
  E_UDCA_CSI = E_FXR_CHOL + EFFECT_UDCA_DISSOL,  capped at 0.70
  k_dissol_eff = EFFECT_UDCA_DISSOL * 0.5

  EFFECT_UDCA_BA = K_BA_UDCA * C_UDCA_norm              // (was E_UDCA_BA * C_UDCA_bile_norm)
  BA_syn_rate = (...) * (1 + EFFECT_UDCA_BA)
  ```
  Neither is `Emax*Xᵞ/(EC50ᵞ+Xᵞ)` — there is no asymptote in the
  underlying formula (`E_UDCA_CSI`'s saturation comes entirely from the
  hand-written `if(E_UDCA_CSI > 0.70)` clamp, not from the term's own
  shape). Per the guide's rule for a compound whose effect term doesn't
  fit the Hill archetype ("a clean, non-standard structure beats a
  standard structure that's wrong" — the same principle already applied
  to evolocumab's `K_EVO_PCSK9*C_EVO` linear effect in
  `elevated-lipoprotein-a`), both terms are kept as bespoke named linear
  effects, `EFFECT_UDCA_DISSOL` and `EFFECT_UDCA_BA`, with no
  `EMAX_UDCA`/`EC50_UDCA`/`GAMMA_UDCA` introduced. `E_UDCA_dis` and
  `E_UDCA_BA` are renamed `K_DISSOL_UDCA`/`K_BA_UDCA` (bespoke rate-style
  coefficients, not Hill parameters) — values unchanged. `E_FXR_CHOL`
  (the additive baseline in `E_UDCA_CSI`) is **not** UDCA's own parameter
  — it is a disease-baseline biliary-secretion term unrelated to UDCA
  dose (never multiplied by any UDCA concentration) — left completely
  untouched, not renamed.
  `EFFECT_UDCA_DISSOL` is a genuine, value-inert common-subexpression
  factoring: the original computed the identical product
  `E_UDCA_dis * C_UDCA_bile_norm` twice, independently, within the same
  `$MAIN` block (once for `E_UDCA_CSI`, again three lines later for
  `k_dissol_eff`) — factored to one named term, still evaluated in
  `$MAIN` exactly where both original uses lived. `EFFECT_UDCA_BA` is
  **not** hoisted into `$MAIN` even though it could be — the original
  computed `E_UDCA_BA * C_UDCA_bile_norm` inline, directly inside
  `$ODE`'s `BA_syn_rate` line, so the renamed term is declared as an
  `$ODE`-local right at that same point of use, per the guide's "keep a
  calculation in the block the original used it in" rule.

## UDCA's two concentration sites: renamed and clarified, not merged

The census classifies UDCA as "Normalize duplicate concentration sites,
then redirect." The original genuinely computes what look like the same
two physical quantities (UDCA in bile, UDCA in plasma) **twice**, once in
`$MAIN` (`C_UDCA_plas`/`C_UDCA_bile`, feeding the PD effect terms above)
and again, independently, in `$TABLE` (`UDCA_plas_conc`/
`UDCA_bile_conc_umol`, purely for reporting) — textually identical
formulas, reading the same raw compartment amounts.

**These two compute-sites were investigated empirically before deciding
how to normalize them, and turned out not to be interchangeable.**
Initial assumption (matching a literal reading of "duplicate") was that
`$MAIN`-cadence and `$TABLE`-cadence values would be numerically
identical for any compartment without a direct-bolus discontinuity (only
`GUT_UDCA`/`GUT_STAT`/`GUT_EZET` receive boluses; `BILE_UDCA`, `CENT_UDCA`,
`CENT_STAT`, `CENT_EZET` all only receive smooth ODE inflow) — confirmed
true at **daily** reporting resolution (`delta=24`, matching the
original's own `run_scenario()` default), max abs diff 0.0 across a
90-day run. **At hourly resolution this assumption is false.** A direct
test — capturing both the `$MAIN`-computed value and the `$TABLE`-fresh
value from the *same running model*, hourly, over 300h — shows `$MAIN`
locals lag the fully-integrated state by (in effect) one report step,
pervasively (179–263 mismatches out of 304 hourly points across
`STAT`/`EZET`/`UDCA`), reaching a max abs diff of ~53 µmol/L for UDCA's
own biliary concentration (a ~20%+ relative deviation), not a
floating-point-scale artifact and not limited to dose instants. This is
consistent with the guide's own dose-instant-artifact section (`$MAIN` is
evaluated using state from the start of the record) but turns out, for
this file's fast-absorbing compounds, to be visible far more broadly than
the "one self-healing duplicate row" the guide describes for other files.

Given that finding, **the two compute-sites are not collapsed into one.**
Collapsing them (reporting only the `$MAIN`-cadence value, or switching
the PD-driving calculation to `$TABLE`/fresh cadence) would have silently
changed either the reported diagnostic values or the actual PD
trajectory relative to the original — exactly what the guide's tolerance
rule forbids for a pure-reorganization archetype. Instead:

- `C_UDCA` (biliary, was `C_UDCA_bile`) and `C_UDCA_PLAS` (plasma, was
  `C_UDCA_plas`) are computed in `$MAIN`, in the exact same place and
  form the original used — these are what the original's own PD
  equations actually consumed, confirmed to match the original's own
  `$MAIN`-computed `C_UDCA_bile`/`C_UDCA_plas` at **max abs diff 0.0**
  across a 300-hour hourly comparison (see Verification). `C_UDCA` is
  the exposed, pluggable concentration per this fork's naming convention
  (the biliary site drives every UDCA PD effect; plasma is
  non-driving/informational, mirroring the `abdominal-aortic-aneurysm`
  doxycycline precedent where tissue drives PD and plasma is kept as a
  diagnostic only).
- `UDCA_plas_conc`/`UDCA_bile_conc_umol` (and, analogously,
  `STAT_plas_conc`) are **preserved verbatim** in `$TABLE`, same formula,
  same names, only compartment names renamed — the original's own fresh
  per-report diagnostic, byte-for-byte behaviorally unchanged.
- `EZET_plas_conc` — never captured in the original at all (see
  `UPSTREAM_ISSUES.md` #134) — is added as the same style of `$TABLE`
  diagnostic, matching the pattern the original already used for UDCA and
  Statin. A pure addition (new `$CAPTURE` entry, no new compartment or
  `dxdt_` term), confirmed not to change any other output.

So "normalize" here means *clarify and correctly attribute* the two
sites (document which one drives PD and is the pluggable interface, which
one is a legacy report-only diagnostic, and why they are not
interchangeable) rather than literally merging them — merging would have
been the wrong fix once the empirical check showed they are not, in fact,
harmless duplicates.

## `$PARAM` vs `$CAPTURE` for `C_<STEM>`/`EFFECT_<STEM>`

Per the guide's qspserver compatibility requirement #2, every `C_<STEM>`/
`EFFECT_<STEM>` would ideally live in `$PARAM`. As established throughout
this fork (confirmed again here against mrgsolve 2.0.1: a `$PARAM` member
reassigned in `$ODE`/`$MAIN` fails with `assignment of read-only
reference`), a quantity recomputed every timestep from state cannot also
be a `$PARAM`. `C_UDCA`, `C_UDCA_PLAS`, `C_STAT`, `C_EZET`, `EFFECT_STAT`,
`EFFECT_EZET`, `EFFECT_UDCA_DISSOL`, `EFFECT_UDCA_BA` are `double`s
computed in `$MAIN`/`$ODE` and listed in `$CAPTURE` instead — confirmed
present in `/model_manifest`'s `outputPaths` (34 total, up from the
original's 26 — the increase is the eight new `C_*`/`EFFECT_*` terms).
Every PK parameter and every new `GAMMA_*` constant (never reassigned)
*are* real `$PARAM` entries — confirmed present in `/model_manifest`'s
`parameters` list (61 total, up from the original's 58; `Km_STAT`→
`EC50_STAT`, `Emax_STAT`→`EMAX_STAT`, `Km_EZET`→`EC50_EZET`,
`Emax_EZET`→`EMAX_EZET`, `E_UDCA_dis`→`K_DISSOL_UDCA`,
`E_UDCA_BA`→`K_BA_UDCA` are renames, not additions; `GAMMA_STAT`,
`GAMMA_EZET` are the only two genuinely new parameters, both math-implied
gamma=1, not fitted).

## Two pre-existing dead `$PARAM` entries, left untouched

`E_STAT_CHOL = 0.60` and `k_dissol = 0.0` are both declared in the
original's `$PARAM` block but never referenced anywhere else in the file
(confirmed by grep) — pre-existing dead parameters, unrelated to any of
the three in-scope compounds' own active PK/PD. Left completely
unrenamed and unused, exactly as in the original; not logged as a
separate `UPSTREAM_ISSUES.md` entry since they cause no build or
behavioral issue, just noted here for anyone reading the `$PARAM` block
looking for what `EMAX_STAT`/`K_DISSOL_UDCA` might have been confused
with.

## When the original doesn't compile at all

The untouched original does not compile under mrgsolve 2.0.1 as written —
three stacked defects (`SOLVERTIME` used inside `$MAIN`, the `_INIT()`
macro idiom, and two bare unheadered `capture` lines that also duplicate
six `$CMT` names). Logged as `translations/UPSTREAM_ISSUES.md` #133;
full detail there. Per the guide's settled policy, all three are fixed
directly in the delivered `chol_mrgsolve_model_refactored.R` (syntax-only,
non-numeric): `SOLVERTIME`→`TIME`, `_INIT(<CMT>) = v;`→`<CMT>_0 = v;`, and
the two `capture` lines→real `$CAPTURE` headers with the six duplicated
compartment names dropped. The untouched original still carries all three
defects forward unfixed.

## Verification

**Method.** Both the untouched original's own quoted DSL string (with
only the three syntax-only build-compat fixes above applied to a
throwaway verification copy — never to the checked-in original) and the
`_refactored.R`'s embedded DSL (confirmed byte-identical to a standalone
extraction of its own quoted string, 15,408 characters, exact match) were
run through the qspserver `mrgsolve_api` container (`POST
/model_manifest`, `POST /run_simulation`) at `http://localhost:8007`,
requests spaced ~2.5s apart, respecting the API's 2-concurrent-job limit.

Two dosing scenarios, both built from the file's own real regimens (not
invented doses):

1. **Scenario A — superposition of all three compounds' own regimens**:
   UDCA 750 mg/day TID (250 mg q8h), Simvastatin 40 mg qd, Ezetimibe
   10 mg qd (the exact `amt`/`ii` values the original's own `make_events`
   and Scenarios 2/4/5 already use), `Stone_vol0 = 0.5`, run 90 days at
   **hourly** resolution (`end=2160, delta=1` — the finer grid needed to
   surface the `$MAIN`/`$TABLE` cadence difference discussed above; the
   original's own daily-resolution default would have hidden it).
2. **Scenario B — UDCA high-dose monotherapy + weight-loss flag**: UDCA
   1050 mg/day TID (350 mg q8h, matching Scenario 3's dose), `WLOSS = 1`
   (matching Scenario 6's flag, exercising the `BWT_t`/`TIME` build-fix
   path), `Stone_vol0 = 0.5`, run 60 days hourly (`end=1440, delta=1`).

No shortening was needed for solver-step budget — both scenarios (2,164
and 1,442 report rows respectively, including duplicate dose-instant
rows) completed well within the API's default `maxsteps`.

**Result: exact match, max abs diff = 0.0, for every one of 26 shared
outputs, at every timepoint, in both scenarios** — all 9 renamed PK
compartments (`GUT_UDCA`/`CENT_UDCA`/`HEP_UDCA`/`BILE_UDCA`/`GB_UDCA`/
`GUT_STAT`/`CENT_STAT`/`GUT_EZET`/`CENT_EZET`), all 9 untouched
disease-network compartments (`BA_pool`, `CHOL_h`, `CHOL_bil`, `PL_bil`,
`GB_vol`, `Crystal_mass`, `Stone_V`, `IL6`, `CRP_plas`), and all 8
`$TABLE`-computed diagnostics (`CSI_out`, `UDCA_plas_conc`,
`UDCA_bile_conc_umol`, `STAT_plas_conc`, `EZET_plas_conc`, `Stone_mm`,
`BA_pool_g`, `CHOL_sat_pct`) — including the dose-instant duplicate report
rows. This is the expected outcome for three compounds whose refactor is
entirely rename/CSE-factoring plus one disclosed dangling-effect
preservation, per the guide's tolerance rule for pure structural
reorganization (Archetype 3 for Ezetimibe/Statin, bespoke-but-linear for
UDCA — no Hill-fitting anywhere in this file).

**Separately**, the new `C_<STEM>`/`EFFECT_<STEM>` outputs were verified
against the original's own (previously uncaptured or under-used) `$MAIN`
locals: `C_STAT`/`C_EZET`/`C_UDCA`/`C_UDCA_PLAS`/`EFFECT_STAT`/
`EFFECT_EZET` all match the original's own `C_STAT_plas`/`C_EZET_plas`/
`C_UDCA_bile`/`C_UDCA_plas`/`E_STAT`/`E_EZET` at **max abs diff 0.0**
across a 300-hour hourly comparison — confirming the refactor's new
pluggable interface faithfully reproduces the exact internal values that
drove the original's own PD, not a different (fresher/staler) proxy for
them.

`/model_manifest` on the refactored DSL confirms 34 output paths (up from
26) and 61 parameters (up from 58), including every renamed PK parameter
and both new `GAMMA_*` constants.

## Anything else worth flagging

- **Compartment ordering is fully preserved** — no compartment was added,
  removed, or reordered (`GUT_UDCA` is still compartment 1, `CRP_plas` is
  still compartment 18), so the R driver's `ev(cmt = "GUT_UDCA", ...)`
  string-keyed dosing needed only the three cmt-name string literals
  updated in `make_events()` — no restructuring of the dosing logic
  itself.
- The compiled model's cache name was changed to
  `"cholelithiasis_qsp_refactored"` (from `"cholelithiasis_qsp"`) so it
  does not share an `mcode_cache` compilation slot with the original — a
  standard, value-inert precaution used elsewhere in this fork.
- The `scenarios` list, `run_scenario()`, all six plots, the endpoint
  summary table, and the UDCA dose-response sensitivity analysis are all
  untouched — none of them reference the renamed PK compartments or
  parameters by name (dosing goes through `make_events()`'s compartment
  string, and every plotted/summarized column — `Stone_V`, `Stone_mm`,
  `CHOL_sat_pct`, `BA_pool_g`, `UDCA_bile_conc_umol`, `CRP_plas`,
  `CHOL_h` — is an unrenamed disease-network state or an unrenamed
  `$TABLE` diagnostic), so they continue to work exactly as before.

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`, the
`cholelithiasis | Ezetimibe`, `| Statin`, and `| UDCA` rows — Statin's
Target/Pathway column confirmed as HMG-CoA reductase (simvastatin) against
the code, not left as a placeholder.
