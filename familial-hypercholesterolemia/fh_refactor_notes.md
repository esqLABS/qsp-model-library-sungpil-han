# Refactor notes — `familial-hypercholesterolemia/fh_mrgsolve_model.R`

Scope: all **three** compounds tracked for this file in
`driver-patches/data/compound_perturbation_census.md` — **Ezetimibe**,
**PCSK9 inhibitor**, and **Statin**. Lomitapide and Bempedoic acid are
flag-driven (`DOSE_LOMT`, `DOSE_BEMP`; no PK compartments of their own) and
are byte-for-byte identical to the original everywhere they appear —
neither is in the census, and neither has anything to isolate/rename.
Inclisiran is mentioned in the file's header comment and README, but the
code contains only a dead placeholder variable (`PCSK9_inclisiran_red =
1.0;`, never read anywhere) — no real PK/PD to refactor; left untouched.

## Confirmed compound identities (not generic census labels)

The census rows list "PCSK9 inhibitor" and "Statin" as generic placeholders.
The code and `README_en.md` both confirm specific real drugs:

- **Statin → Rosuvastatin.** `$PARAM` comment: "Statin PK (Rosuvastatin-like,
  oral, mg → μg/mL)"; `README_en.md` scenario table names "Rosuvastatin 40
  mg/d" explicitly (lines 117–121); the R script's own scenario labels say
  "Rosuvastatin 40 mg" throughout.
- **PCSK9 inhibitor → Evolocumab.** `$PARAM` comment: "PCSK9 inhibitor PK
  (Evolocumab-like, subcutaneous mAb)"; dose 420 mg q4w and MW 144,000 Da
  both match evolocumab exactly (real evolocumab MW ≈144 kDa, FOURIER-trial
  dosing regimen); `README_en.md` names "Evolocumab 420 mg q4w" and cites
  the FOURIER trial (Sabatine 2017) as its calibration source.
- **Ezetimibe** is already a specific real drug in the original — no
  disambiguation needed. Target NPC1L1, 10 mg/day, matches IMPROVE-IT
  (Cannon 2015) per the README's own reference table.

## Archetype per compound

**Statin (Rosuvastatin): bespoke, not a listed archetype.** `GUT_STATIN`
feeds **both** `CENT_STATIN` (systemic) and `LIV_STATIN` (hepatic site of
action) as two independent first-pass splits of the same absorbed dose —
`dxdt_LIV_STATIN` is driven directly from `GUT_STATIN` (times a fixed 0.80
first-pass-extraction fraction), not from `CENT_STATIN` via a Q-mediated
exchange. This is *not* Archetype 3's central/peripheral pair (no
inter-compartmental clearance between `CENT_STATIN` and `LIV_STATIN` at
all) — it is the original's own structure, preserved and renamed exactly,
per the guide's "none of these fit" carve-out. Two concentration sites
genuinely matter here (the guide's own exception clause): `C_STATIN`
(hepatic, μg/mL) is the single PD-facing, canonical concentration — it is
the *only* one any effect equation reads (`EFFECT_STATIN`, i.e. HMGCR
inhibition). `C_STATIN_SYS` (systemic, μg/mL) is retained purely for
reporting, matching what the original's `C_statin_sys`/`C_sys_S` already
were — read by nothing.

**PCSK9 inhibitor (Evolocumab): Archetype 4 (TMDD).** `GUT_PCSK9I` (was
`SC_PCSK9I`) → `CENT_PCSK9I` → `PERI_PCSK9I` is a standard depot+
central+peripheral PK chain; `COMPLEX_PCSK9I` (was `COMP_PK9`) is the
drug–target complex, formed from `PCSK9I_nM` (drug) binding `PCSK9_pl`
(target). This *is* real receptor-binding kinetics as its own dynamic
system (`KON_PCSK9I`/`KOFF_PCSK9I` binding ODEs, not an algebraic ratio) —
kept in full, not flattened. One deliberate deviation from the guide's
worked TCZ example: **the TMDD "free target" (`PCSK9_pl`) is not renamed
to a `PCSK9I_`-stemmed name.** `PCSK9_pl` is declared under "PD
compartments" in the *original's own* `$CMT` block (not under the PK
section with the other PCSK9i compartments), is read and dosed
independently of any drug (`kin_PK9`/`kout_PK9` turnover, statin feedback,
its own scenario-level `init_vals["PCSK9_pl"] <- 400` override), and drives
LDLR degradation (`PCSK9_effect_LDLR`) regardless of whether any PCSK9i is
present at all. It is disease physiology that the drug happens to bind,
not the drug's own receptor pool — renaming it to a compound-stemmed name
would misrepresent what it is. Only `COMPLEX_PCSK9I` (the drug-specific
bound fraction) is renamed to the compound's stem.

**Ezetimibe: Archetype 3 without a peripheral compartment** (depot +
central only, exactly as the original already models it — nothing added or
removed).

## Naming

| Role | Statin | PCSK9 inhibitor | Ezetimibe |
|---|---|---|---|
| Depot | `GUT_STATIN` (was `GUT_S`) | `GUT_PCSK9I` (was `SC_PCSK9I`) | `GUT_EZE` (unchanged) |
| Central | `CENT_STATIN` (was `CENT_S`) | `CENT_PCSK9I` (unchanged) | `CENT_EZE` (unchanged) |
| Peripheral / 2nd site | `LIV_STATIN` (was `LIV_S`; hepatic, bespoke — see above) | `PERI_PCSK9I` (unchanged) | — |
| TMDD complex | — | `COMPLEX_PCSK9I` (was `COMP_PK9`) | — |
| TMDD free target | — | `PCSK9_pl` (unchanged; disease-side, deliberately not stem-renamed — see above) | — |
| `KA_/F_/CL_/V1_/V2_/Q_` | `KA_STATIN`/`F_STATIN`/`V1_STATIN`(per-kg)/`CL_STATIN`/`CL_HEP_STATIN` (was `ka_S`/`F_S`/`Vd_S`/`CL_S`/`CL_H`) | `KA_PCSK9I`/`F_PCSK9I`/`V1_PCSK9I`/`V2_PCSK9I`/`CL_PCSK9I`/`Q_PCSK9I`/`KTMDD_PCSK9I` (was `ka_P`/`F_P`/`Vc_P`/`Vp_P`/`CL_P`/`Q_P`/`kTMDD`) | `KA_EZE`/`F_EZE`/`V1_EZE`/`CL_EZE` (was `ka_EZE`/`F_EZE`/`Vd_EZE`/`CL_EZE`) |
| Binding on/off | — | `KON_PCSK9I`/`KOFF_PCSK9I` (was `kon_PK9`/`koff_PK9`) | — |
| Exposed concentration | `C_STATIN` (hepatic, μg/mL; was `C_liv_S`/`C_statin_liv`, two names for one quantity) | `C_PCSK9I` (mg/L; was `C_pcsk9i_mgL`, $TABLE-only, now the single site both PK and reporting read) | `C_EZE` (mg/L; was `C_EZE`/`C_eze_mgL`, two names for one quantity) |

`V1_STATIN` and `V1_EZE` are **per-kg** apparent volumes (`Vd_S`/`Vd_EZE`,
combined with the shared `BW` parameter at every use site, exactly as the
original did) — a deliberate deviation from the naming table's flat-volume
default, disclosed here rather than silently reinterpreted.

`DOSE_STATIN`/`DOSE_PCSK9I`/`DOSE_EZE` (was `DOSE_S`/`DOSE_P`/`DOSE_EZE`)
are **not read by any `$ODE`/`$TABLE` expression** in either the original
or the refactor — dosing happens entirely through compartment events
(`GUT_STATIN`/`GUT_PCSK9I`/`GUT_EZE`). These three parameters are pure
R-side bookkeeping (used only to build `params_update` for record-keeping/
display), same treatment as the `RTOT_TCZ`/`DOXY_ON` precedents elsewhere
in this corpus: preserved as-is, not wired up, not removed.

All parameter *values* are copied verbatim from the original
(`ka_S=0.45` → `KA_STATIN=0.45`, etc.) — nothing invented, nothing
defaulted, nothing dropped.

## The Hill interface

**Statin → `EFFECT_STATIN`.** The original's `Inh_S = Emax_SI * C_liv_S /
(EC50_SI + C_liv_S)` was already exactly the Hill shape with no exponent.
Pure rename/promotion: `Emax_SI`→`EMAX_STATIN`, `EC50_SI`→`EC50_STATIN`,
`GAMMA_STATIN=1` added (no original Hill coefficient). `EFFECT_STATIN =
EMAX_STATIN*pow(C_STATIN,GAMMA_STATIN)/(pow(EC50_STATIN,GAMMA_STATIN)+
pow(C_STATIN,GAMMA_STATIN))` — with `GAMMA=1`, `pow(x,1)==x`, so this is
arithmetically identical to the original (confirmed by the exact-match
verification below, same empirical precedent already established for
`GAMMA=1` renames elsewhere in this corpus, e.g. `abdominal-aortic-
aneurysm`).

**Ezetimibe → `EFFECT_EZE`.** The original's `EZE_effect_VLDL = 1.0 -
fEZE_abs*C_EZE/(0.05+C_EZE)` buried its EC50 as an **embedded numeric
literal (`0.05`)** with no named parameter at all — the only one of the
three compounds with this pattern. Promoted (not fit): `fEZE_abs` →
`EMAX_EZE`, the literal `0.05` → `EC50_EZE=0.05` (new named `$PARAM`),
`GAMMA_EZE=1` added. `EFFECT_EZE` computed the same Hill-formula way as
`EFFECT_STATIN`; `EZE_effect_VLDL = 1.0 - EFFECT_EZE` — same value, same
shape, now discoverable as a parameter instead of a literal buried in an
expression.

**PCSK9 inhibitor → `EFFECT_PCSK9I` (TMDD, occupancy — not a fit).**
Per the guide's Archetype-4 pattern (`COMPLEX_TCZ/RTOT_TCZ` in the
`rheumatoid-arthritis` precedent), the Hill input for a TMDD compound is
fractional target engagement, not raw concentration. This model has no
fixed `RTOT` (unlike a classic membrane receptor pool with homeostatic
turnover to a constant total) — `PCSK9_pl` undergoes its own continuous
synthesis (`kin_PK9`) and clearance (`kout_PK9`), with the *current* total
circulating PCSK9 (free + bound) varying over time and with statin
feedback. So the occupancy denominator is defined dynamically:

```c
double OCC_PCSK9I    = (COMPLEX_PCSK9I * 1000.0) / (PCSK9_pl + COMPLEX_PCSK9I * 1000.0);
double EFFECT_PCSK9I = OCC_PCSK9I;
```

(the `*1000.0` factor matches the unit-equivalence convention the
*original* already uses everywhere it combines `PCSK9_pl` [ng/mL] and
`COMP_PK9`/`COMPLEX_PCSK9I` [nmol/L-equiv], e.g. `PCSK9_free = PCSK9_pl -
COMP_PK9*1000.0`).

**This is a diagnostic addition, not a rewiring of the disease
mechanism.** Unlike the RA/tocilizumab precedent — where `EFFECT_TCZ` *is*
literally the same quantity the original's own `IL6_sig` term already
consumed, so exposing it under the new name changed nothing about how the
disease equations are driven — this file's disease side does **not** read
occupancy at all in either the original or the refactor. The original's
actual PCSK9i-disease coupling is pure mass action, spread across two
ODEs: `COMPLEX_PCSK9I` accumulates bound PCSK9 (removing it from
`PCSK9_pl`'s availability), and `PCSK9_free = PCSK9_pl -
COMPLEX_PCSK9I*1000` (unchanged, disease-side) is what actually drives
LDLR degradation via the pre-existing `PCSK9_effect_LDLR` Hill term (using
`EC50_PK9_LR`, a disease parameter, not a drug parameter). `EFFECT_PCSK9I`/
`OCC_PCSK9I` is a **new, purely reported** quantity — it does not feed
back into `$ODE` anywhere, so its addition cannot change any simulated
trajectory (confirmed below: it doesn't move a single existing output by
more than floating-point noise). It exists so a downstream tool can
discover "how neutralized is this compound's target" as a bounded [0,1]
number, per the guide's driveability intent, without inventing a
mechanism the original doesn't have.

`EC50_PCSK9I = (KOFF_PCSK9I + kout_PK9) / KON_PCSK9I ≈ 2.0008` nM is the
analytically-derived quasi-steady-state Kd for this binding system — the
same closed-form pattern as `EC50_TCZ = (KOFF_TCZ+KDEG_TCZ)/KON_TCZ` in the
RA precedent (`kout_PK9` plays the same role there as `KDEG_TCZ`: the
first-order clearance rate of the free/complex pool). **Not consumed by
`EFFECT_PCSK9I`'s formula above** (which is the exact ODE-state ratio, no
fit needed) — exposed purely for discoverability/future closed-form use,
identical treatment to `EMAX_TCZ`/`EC50_TCZ`/`GAMMA_TCZ` in the RA file.
No `nls()` fit was performed or needed here, since nothing about
`EFFECT_PCSK9I` is being approximated — it is read directly off ODE state,
not fit to it.

## Five pre-existing build defects found and fixed (syntax-only)

None of these five is inside any of the three compounds' own PK/PD block
(defect 5 touches the disease-side `PCSK9_free` line, not compound PK,
though the rename of `COMP_PK9`→`COMPLEX_PCSK9I` that this refactor already
needed to make happens to sit right next to it). Logged as
`translations/UPSTREAM_ISSUES.md` entry **#146**, with the full reproduction
detail there. Summary:

1. **`$PLUGIN autodiff` unavailable** in this mrgsolve/qspserver build
   (`Error: plugin autodiff could not be found.`). Unused elsewhere in the
   DSL — dropped.
2. **`$INIT` evaluated as one `list(...)` call with no access to `$PARAM`
   values** — fails on any line referencing a sibling `$PARAM` (six of the
   original's seventeen `$INIT` lines do: `PCSK9_pl = PCSK9_0`, `VLDL_C =
   VLDL0`, `IDL_C = IDL0`, `LDL_C = LDL0_het`, `HDL_C = HDL0`, `TG_C =
   TG0`).
3. **`$CMT` + a separate `$INIT` jointly declaring the same compartments is
   rejected outright** (`Duplicated model names: ...`), independent of
   defect 2.
   - **Fix for both 2 and 3**: deleted the `$INIT` block entirely; every
     initial value (including the six param-referencing ones) is now set
     in `$MAIN` via the `<cmt>_0 = <expr>;` idiom, e.g. `PCSK9_pl_0 =
     PCSK9_0;`. Same values, same defaults, syntax-only.
4. **`$CAPTURE` repeating nine compartment names already in `$CMT`**
   (`HMGCR_rel, LDLR_rel, PCSK9_pl, COMP_PK9→COMPLEX_PCSK9I, VLDL_C, IDL_C,
   LDL_C, HDL_C, TG_C`) is rejected (`compartment should not be in
   $CAPTURE`). Fix: removed these nine from `$CAPTURE` — they still appear
   in every output row automatically as compartment states, so nothing is
   lost.
5. **`PCSK9_free` declared twice** (`$ODE` and `$TABLE`, same name,
   equivalent formula) — a genuine C++ redefinition error, confirming
   `$ODE`/`$TABLE` share one variable scope in this build. Fix: kept the
   single `$ODE` declaration (needed there — it's read inside the
   differential system every solver substep, not just at reporting time)
   and `$CAPTURE`d it directly; removed the `$TABLE` duplicate.

**Each fix was isolated one at a time**, resubmitting to
`/model_manifest` after each patch and confirming the *next* reported
error was independent of the one just fixed — five distinct compile
errors in sequence, not one defect surfacing differently each time.

## Normalized duplicate concentration sites (the census's own ask, and separately confirmed necessary)

Beyond the five compile-blocking defects, the original also computes each
of Statin's and Ezetimibe's concentrations **twice under two different
names** — `C_liv_S`/`C_sys_S` in `$ODE` (needed there, feeding `Inh_S`
every substep) vs. `C_statin_liv`/`C_statin_sys` recomputed with the
identical formula in `$TABLE`; `C_EZE` in `$ODE` vs. `C_eze_mgL`
recomputed in `$TABLE`. Because the two names differ, this does **not**
block compilation (unlike `PCSK9_free`'s exact-name collision above) — but
it is precisely what the census's "Normalize duplicate concentration
sites, then redirect" classification for Statin and Ezetimibe calls for,
and it is also what the guide's "discoverability" contract needs (a single
contiguous `double C_<STEM> = <expr>;` statement per compound). Fixed by
declaring `C_STATIN`/`C_STATIN_SYS`/`C_EZE` once each in `$ODE` (kept
exactly where the original computed the PD-feeding one, per the guide's
"keep a calculation in the block the original used it in" rule) and
`$CAPTURE`-ing them directly, with no second `$TABLE` recomputation.
PCSK9 inhibitor's own concentration was already a single, clean site
(`C_pcsk9i_mgL`, `$TABLE`-only, matching the census's "clean single site"
classification for it) — moved to `$ODE` as `C_PCSK9I` so `PCSK9I_nM` (used
in the TMDD binding ODEs) derives from it instead of being recomputed
independently from `CENT_PCSK9I` a second time.

**Disclosed cosmetic side-effect**: this file's own tests found no
"dose-instant reporting artifact" for these specific variables in the
verification run below (see the exact-match result), but per the guide's
documented pattern, an `$ODE`-declared `double` can in principle read a
stale value on the reporting row generated exactly at a dose instant. If
this shows up under other dosing patterns, the fix (matching
`clostridioides-difficile-infection`'s precedent) would be to declare
`C_STATIN`/`C_PCSK9I`/`C_EZE` as `$GLOBAL` macros instead — not applied
here since the original had no pre-existing `$GLOBAL`-macro convention to
match, and the verification below shows no artifact at any of the shared
dosing scenarios' reporting times.

## Verification

Per the guide's mandatory protocol: ran every one of the original file's
own six R-script scenarios (`1. HetFH — No Treatment`, `2. HetFH —
Rosuvastatin 40 mg/d`, `3. HetFH — Rosuvastatin 40 mg + Ezetimibe 10 mg`,
`4. HetFH — Evolocumab 420 mg q4w`, `5. HetFH — Rosuvastatin 40 mg +
Evolocumab 420 mg q4w`, `6. HomFH — Lomitapide + Statin + Evolocumab`)
through both models via the qspserver `mrgsolve_api` service
(`POST /model_manifest` then `POST /run_simulation`), identical dosing
(`ii`/`addl` reproducing each scenario's own daily/q4w regimen), identical
parameter overrides (`LDLR_fxn`, `LDL0_het`, `DOSE_*`, `DOSE_LOMT`,
`DOSE_BEMP`), identical time grid (`end=365*24=8760h`, `delta=12h`,
matching the original R script's own `mrgsim(..., end=365*24, delta=12)`
exactly — no shortening needed, all runs completed within the API's
default solver step budget).

Both the original (patched only with the five syntax-only build fixes
above, applied in-memory to a scratch copy — never to the checked-in file)
and the refactored model were run through the *same* two-in-a-row
`/run_simulation` request pattern, then every shared output compared
point-by-point across the full ~731-733-point time grid per scenario:
all 17 disease/PD outputs common to both (`HMGCR_rel, LDLR_rel, PCSK9_pl,
VLDL_C, IDL_C, LDL_C, HDL_C, TG_C, LDLR_pct, PCSK9_free, NonHDL_C, TC,
LDL_reduction_pct, CVD_risk_10yr, LDL_goal_55, LDL_goal_70`, plus
`GUT_EZE`/`CENT_EZE`) and all 11 renamed PK/concentration pairs
(`GUT_S`↔`GUT_STATIN`, `CENT_S`↔`CENT_STATIN`, `LIV_S`↔`LIV_STATIN`,
`SC_PCSK9I`↔`GUT_PCSK9I`, `CENT_PCSK9I`↔`CENT_PCSK9I`,
`PERI_PCSK9I`↔`PERI_PCSK9I`, `COMP_PK9`↔`COMPLEX_PCSK9I`,
`C_statin_liv`↔`C_STATIN`, `C_statin_sys`↔`C_STATIN_SYS`,
`C_pcsk9i_mgL`↔`C_PCSK9I`, `C_eze_mgL`↔`C_EZE`).

**Result: exact match everywhere except one point.** Across all 6
scenarios × 29 compared series (≈4,400 individual value comparisons):
max absolute difference **0.0** in 5 of 6 scenarios; scenario 4
(Evolocumab 420 mg q4w monotherapy) showed a single point of difference —
`CENT_PCSK9I` at t=8100h, 269.3825 vs. 269.3824 (abs diff `1.0e-4`, relative
diff `3.7e-7`) — one value out of ~730 time points for that variable,
consistent with floating-point/serialization noise, not a structural
divergence (every other point in that same series, and every other
variable in that same scenario, matched exactly). This is the exact-match
result the guide anticipates for Archetypes 3/4 plus `gamma=1` Hill
renames (no fitting was needed or performed anywhere in this file).

`EFFECT_STATIN`/`EFFECT_EZE`/`EFFECT_PCSK9I`/`OCC_PCSK9I` (new,
discoverability-only outputs that don't exist in the original) were
separately sanity-checked via `/run_simulation` on the triple-therapy
scenario: all bounded in `[0,1]` and behaving as expected (`EFFECT_STATIN`
→0.81 at steady exposure matching `EMAX_STATIN=0.92` asymptote,
`EFFECT_EZE`→0.52 approaching `EMAX_EZE=0.55`, `EFFECT_PCSK9I`/
`OCC_PCSK9I`≈0.99, consistent with evolocumab's near-irreversible binding
— `KOFF_PCSK9I=1e-5` is tiny relative to `KON_PCSK9I=0.012`).

Confirmed via `POST /model_manifest` that the refactored DSL compiles and
`outputPaths`/`parameters`/`compartments` list every renamed and newly
added name (`C_STATIN`, `C_STATIN_SYS`, `C_PCSK9I`, `C_EZE`,
`EFFECT_STATIN`, `EFFECT_PCSK9I`, `EFFECT_EZE`, `OCC_PCSK9I`, and all
`EC50_*`/`EMAX_*`/`GAMMA_*` parameters).

`Rscript -e 'parse("fh_mrgsolve_model_refactored.R")'` succeeds with no
error (31 top-level expressions parsed).

## qspserver infrastructure note (not specific to this file)

The `mrgsolve_api`'s shared `inline` project cache
(`/soloc/mrgsolve-so-2.0.1-x86_64-pc-linux-gnu/inline/mrgmod_cache.RDS`)
was corrupted at the start of this session (`readRDS`: "embedded nul in
string") — every `model_content` request goes through this one shared
cache bucket (`mcode_cache(model="inline", ...)` in `_runner.R`), so this
blocked *any* inline submission, not just this file's. Matches the known
issue already documented in `driver-patches/HANDOFF.md`. Cleared with
`docker exec qspserver-mrgsolve_api-1 rm -f
/soloc/mrgsolve-so-2.0.1-x86_64-pc-linux-gnu/inline/mrgmod_cache.RDS`,
confirmed fixed by a subsequent successful `/model_manifest` call. Flagged
here in case it recurs for the next agent picking up this queue.

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`, the
`familial-hypercholesterolemia | Ezetimibe`, `familial-hypercholesterolemia
| PCSK9 inhibitor`, and `familial-hypercholesterolemia | Statin` rows —
compound identities confirmed as Ezetimibe, Evolocumab, and Rosuvastatin
respectively (see above).
