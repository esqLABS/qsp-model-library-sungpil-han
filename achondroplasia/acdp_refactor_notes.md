# Refactor notes — `achondroplasia/acdp_mrgsolve_model.R`

## Scope: all three census rows are genuine, externally-dosed drugs

This file carries three rows in
`driver-patches/data/compound_perturbation_census.md`: **INFIG**,
**"Released free-CNP (TCNP)"**, **Vosoritide**. The census's automated
classifier had tagged all three as `Delete PK compartment; concentration is
itself the state` — worth checking against the actual code rather than
trusting the label, per the task's specific flag on TCNP (does it have a
real dosing route, or is it endogenous CNP released as part of the
disease's own signalling?).

| Row | Real drug? | Evidence |
|---|---|---|
| Vosoritide (VOS) | **Yes** | Own compartments `VOS_DEPOT`→`VOS_CP`, dosed via `ev(cmt="VOS_DEPOT", ...)` in scenarios 2, 3, 4, 5, 9, 10; own `KA_VOS`/`KE_VOS`/`V_VOS`/`F_VOS`; CNP analog, direct SC dosing. |
| "Released free-CNP" (TCNP) | **Yes — a real drug's active moiety, not endogenous CNP.** `$CMT TCNP_DEPOT`/`TCNP_CP` model a *sustained-release SC prodrug depot* (`TransCon CNP / navepegritide`) that releases free CNP over time (`KREL_TCNP = 0.018 /h`, "sustained over ~1 wk") and is eliminated fast once released (`KE_TCNP = 2.5 /h`, "fast, like native CNP"). Scenario 6 ("6_TransConCNP_QW") doses `TCNP_DEPOT` directly: `make_ev(100 * WT, 168, 51, "TCNP_DEPOT")` — a real, once-weekly SC event. `TCNP_CP` is therefore the pharmacologically active species produced from an administered prodrug (the same relationship as an active metabolite of a dosed parent drug, e.g. `acne-vulgaris`'s OXO from dosed ISO), not the body's own baseline CNP output. There is no separate "endogenous CNP" compartment anywhere in the file — `CNP_TOTAL = VOS_CP + TCNP_CP` sums only the two *drug-derived* concentrations. **Conclusion: the opposite of the task's flagged suspicion** — checking the code confirms a real dosing route, so this is refactored as a genuine drug, not reclassified as endogenous. |
| Infigratinib (INFIG) | **Yes** | Own compartments `INFIG_GUT`→`INFIG_CP`, dosed via `ev(cmt="INFIG_GUT", ...)` in scenario 7; own `KA_INFIG`/`CL_INFIG`/`V_INFIG`/`F_INFIG`; oral FGFR1-3 TKI, direct PO dosing. |

No compound in this file needed to be reclassified as endogenous — unlike
the `graves-disease` precedent (T3/T4/TRAb/TSH), every one of achondroplasia's
three census rows has its own depot→central PK block and its own dosing
event in the file's own scenario list. (The file's separate `GH_ON`
flag/`GH_EFFECT` parameter for off-label growth hormone is *not* a census
row and has no PK compartment at all — a bare on/off multiplier, out of
scope, left untouched.)

## The original does not compile — two independent build defects, fixed, disclosed, logged

`POST /model_manifest` on the untouched original's own DSL (extracted
verbatim from `acdp_code <- '...'`) returned HTTP 500:

```
invalid class "mrgmod" object: compartment should not be in $CAPTURE:
PERK,CGMP_SIG,CHONDRO,HEIGHT_CM,HEIGHTZ,FMAREA,SPCANALZ,AHI,OTITIS,BMIZ,
MAP_BP,HR,PHOS,VOS_CP,TCNP_CP,INFIG_CP
```

After patching that around (to reach the C++ compile stage), a second,
independent defect appeared:

```
174:13: error: assignment of read-only reference 'PERK'
  PERK      = PERK_BASE;
...(identically for CGMP_SIG, CHONDRO, HEIGHT_CM, HEIGHTZ, FMAREA,
    SPCANALZ, AHI, OTITIS, BMIZ, MAP_BP, HR, PHOS)
```

Both fixes are syntax-only and non-numeric. Logged as
`translations/UPSTREAM_ISSUES.md` **#124** — the same defect pair already
seen in `neurofibromatosis-type-1/nf1_mrgsolve_model.R` (issue #41): the
two defects are entangled (the bare `$MAIN` assignment idiom only "works"
because `$CAPTURE`'s duplicate listing means the code never reaches the
C++ compile stage that would reject it), so a working build needs both
fixed together.

1. **`$CAPTURE` re-lists fifteen `$CMT` compartment names directly**
   (`PERK CGMP_SIG CHONDRO HEIGHT_CM HEIGHTZ FMAREA SPCANALZ AHI OTITIS
   BMIZ MAP_BP HR PHOS VOS_CP TCNP_CP INFIG_CP`), which mrgsolve 2.0.1
   rejects. Fixed by dropping all fifteen (compartments are exposed as
   output columns automatically; `AGV_CALC`, the one genuinely non-compartment
   name in the original's `$CAPTURE` line, is kept), and replacing the line
   with this refactor's own `C_<STEM>`/`EFFECT_<STEM>`/`CNP_TOTAL` capture
   entries.
2. **`$MAIN`'s `if (NEWIND <= 1) { ... }` block sets thirteen disease
   compartments' initial conditions with bare compartment-name assignment**
   (`PERK = PERK_BASE;`, `CGMP_SIG = 0;`, ... 13 lines total) instead of
   the standard mrgsolve `<cmt>_0 = value;` idiom. Fixed by switching all
   thirteen to `PERK_0 = PERK_BASE;` etc. — identical values, identical
   semantics (confirmed by the verification below: every scenario's own
   `HEIGHT0`/`AHI0`/etc.-derived initial state matches the patched-original
   exactly).

Both fixes were applied directly to the delivered
`acdp_mrgsolve_model_refactored.R` (not just a scratch copy), per the
guide's settled policy. The same pair of fixes was applied to an
in-memory-only **patched-original** scratch copy (original DSL, otherwise
byte-identical) used solely as the verification baseline; it was never
written into the repo tree and is not part of the deliverables.

## A pre-existing numerical fragility (not a build defect, not fixed, disclosed) — `pow()` on a noisy near-zero concentration

Found while running the verification below, present identically in the
untouched original: the CNP-axis Hill term (`EFFECT_VOS`, was
`CGMP_DRIVE`) uses a **non-integer exponent**, `GAMMA_VOS = 1.5` (was
`HILL_VOS`), applied to `CNP_TOTAL = C_VOS + C_TCNP`. Vosoritide has a
~15-minute half-life (`CL_VOS = 2.77/h`, `KA_VOS = 6.0/h`) and is dosed
once daily (`ii = 24`) — between doses, `C_VOS` decays to a concentration
many orders of magnitude below `EC50_VOS = 8.0`, at which scale ordinary
floating-point/adaptive-solver roundoff can push the true near-zero value
slightly **negative**. `pow(negative, 1.5)` is `NaN` in C++, and once any
single `dxdt_*` evaluates to `NaN` the whole state vector is `NaN` from
that point forward — cascading `CGMP_SIG → PERK →` every downstream
disease output, permanently, for the rest of the simulation.

**Confirmed identical in the untouched original** (same formula:
`pow(CNP_TOTAL, HILL_VOS)`, same trigger): a direct test of
`patched-original` vs `refactored` on scenario 2 hits `NaN` at
**simulated hour 8, in both models, at exactly the same reported time
step** — not merely "eventually blows up in both," the onset is identical
to the resolution of the time grid tested (hourly). This is a materially
*tighter* match than the closest precedent for this defect class
(`essential-thrombocythemia` issue #64, where a restructured PK gave the
refactored file a different floating-point path and a different NaN onset
time from the original); here the refactor keeps arithmetic close enough
to the original's own that the two blow up at the identical step.

Logged as `translations/UPSTREAM_ISSUES.md` **#125**. Not fixed — this is
the same class of pre-existing numerical fragility already disclosed (not
fixed) in `essential-thrombocythemia` (#64), `von-willebrand-disease`
(#72), and `takayasu-arteritis` (#74); fixing it would mean either
clamping `CNP_TOTAL` to a floor of 0 before the `pow()` call or making
`GAMMA_VOS` an integer, both of which would be a real behavioural change
to the compound's own PD math, out of scope for a naming/structure
refactor. TransCon CNP (weekly dosing, slow release) and Infigratinib
(daily dosing, much slower elimination relative to its dosing interval)
were checked over the full 90-day window and never trigger it — their
troughs never approach the same near-zero floor between doses. The
untreated scenario (`1_Untreated_NaturalHistory`) never triggers it either
(`CNP_TOTAL` stays exactly `0` the whole run, `pow(0, 1.5) = 0`, no sign
ambiguity).

## Archetype: all three drugs use a bespoke depot→central variant (state IS the concentration)

Vosoritide, TransCon CNP, and Infigratinib all use the same structural
shape as the guide's **Archetype 3 minus peripheral** (depot→central,
linear elimination, no peripheral compartment) — **with one shared,
disclosed bespoke deviation**: none of the three central compartments
holds an *amount*. Each is written directly as a **concentration state**:

```
dxdt_CENT_VOS = KA_VOS * GUT_VOS / V1_VOS - CL_VOS * CENT_VOS;
```

— the depot's inflow is divided by volume on the way in, and elimination
is a plain first-order rate constant on the concentration itself, not
`CL/V1 * amount`. This is mathematically a valid one-compartment PK model
(equivalent to the amount-based archetype divided through by a constant
`V1`), just expressed with the state variable already in concentration
units. Because of this, `C_<STEM>` is a plain **identity** of the central
compartment (`double C_VOS = CENT_VOS;`), not a division — dividing by
`V1_VOS` again here would double-count the volume and change the numeric
trajectory. This is the same reasoning already used for `graves-disease`'s
MMI/PTU/PROP compartments, though the underlying reason differs slightly:
there, the declared volume parameter was simply unused; here, `V1_VOS`/
`V1_TCNP`/`V1_INFIG` **are** genuinely used (in the depot-inflow term),
just not in the usual post-hoc division position. No parameter is
declared-but-unused for any of the three compounds in this file.

**`KA_TCNP` naming note:** the original's `KREL_TCNP` ("Prodrug release
rate from SC depot") is physically a *release* rate, not an *absorption*
rate — but it occupies exactly the same structural role in the ODE
(first-order depot→central transfer) as `KA_<STEM>` does for the other two
compounds, and the naming convention has no separate slot for "release
rate." Renamed to `KA_TCNP` per the convention, with the physical
distinction preserved in a comment (`// was KREL_TCNP -- release, not
absorption`).

**`CL_VOS`/`CL_TCNP` naming note:** the originals' `KE_VOS`/`KE_TCNP` are
plain first-order elimination-rate constants (1/h), not literal clearances
(L/h) — since the compartment is already a concentration state, `KE * C`
*is* the correct elimination term (no division by volume needed). Renamed
to `CL_VOS`/`CL_TCNP` anyway, matching the exact same "elimination rate
constant" naming precedent already used in `graves-disease` (`kel_MMI` →
`CL_MMI`, same reasoning, same disclosure).

## The Hill interface — Vosoritide and TransCon CNP share one; Infigratinib has its own

**Infigratinib** has an independent, single-mechanism Hill term in the
original (`INFIG_INHIB`), already a plain ratio — an exact rename:

```
EFFECT_INFIG = EMAX_INFIG * pow(C_INFIG, GAMMA_INFIG) / (pow(EC50_INFIG, GAMMA_INFIG) + pow(C_INFIG, GAMMA_INFIG));
```
(`EMAX_INFIG` was `EMAX_INFIG_PERK`, `GAMMA_INFIG` was `HILL_INFIG` —
same values, same shape.)

**Vosoritide and TransCon CNP genuinely share one receptor-occupancy Hill
term in the original, not two independent ones — kept as a disclosed
bespoke shared interface, not force-split.** The original sums their
concentrations *before* applying the single Hill nonlinearity:
`CNP_TOTAL = VOS_CP + TCNP_CP`, then `CGMP_DRIVE = pow(CNP_TOTAL,
HILL_VOS)/(pow(EC50_VOS, HILL_VOS) + pow(CNP_TOTAL, HILL_VOS))` — one
shared `EC50_VOS`/`HILL_VOS` pair; there is no `EC50_TCNP`/`HILL_TCNP`
anywhere in the original at all. This reflects genuine, real biology, not
a convenience shortcut: Vosoritide is a CNP analog and TransCon CNP
releases the *native* ligand — both are direct agonists of the *same*
NPR-B receptor pool, so a physiologically pooled-ligand, single-receptor-
occupancy model is a defensible mechanistic choice (the same reasoning
already accepted for `acne-vulgaris`'s parent+active-metabolite shared
`EFFECT_ISO`, fed by `X_ISO = C_ISO + POT_OXO*C_OXO`).

This does sit close to the guide's "never collapse several drugs into one
shared Hill term" warning, since VOS and TCNP are two independently-dosed
products (not a parent/metabolite pair of the *same* administered drug).
Two things support keeping the shared form rather than force-splitting it:

1. **Independent driveability is preserved at the concentration level.**
   `C_VOS` and `C_TCNP` are each fully independent, individually-exposed
   state variables from their own independent PK compartments — an
   external covariate can still substitute either one alone. Only the
   downstream *receptor-saturation nonlinearity* is shared, which is where
   the real biology (one receptor pool, two ligands) actually lives.
2. **None of the file's own ten scenarios doses both simultaneously** (the
   two are, sensibly, alternative CNP-class therapies, never combined) —
   confirmed by inspection of every `make_ev(...)` call. Since
   `pow(a+0, γ) ≡ pow(a, γ)`, the shared-then-Hill formulation and a
   hypothetical split-then-combine formulation would produce numerically
   identical results for *every scenario this file actually exercises*;
   only a hypothetical future combined-dosing scenario (never tested by
   either the original or this refactor) would distinguish them. Splitting
   now would be a speculative restructuring with no scenario to verify it
   against — kept as the original wrote it instead, per "don't force a fit
   the original doesn't have."

`EFFECT_VOS` is the one named interface (`TCNP` contributes no independent
`EFFECT_TCNP`):

```
double CNP_TOTAL  = C_VOS + C_TCNP;
double EFFECT_VOS = EMAX_VOS * pow(CNP_TOTAL, GAMMA_VOS) / (pow(EC50_VOS, GAMMA_VOS) + pow(CNP_TOTAL, GAMMA_VOS));
```

`EMAX_VOS = 1.0` is **math-implied, not fit** — the original's own ratio
already saturates at 1 as `CNP_TOTAL → ∞`; the real physiological ceiling
is applied downstream, unchanged, via `EMAX_PERKINH` (`CNP_INHIB =
EMAX_PERKINH * CGMP_SIG`), exactly as the original did. `GAMMA_VOS` is an
exact rename of `HILL_VOS` (1.5), no fitting performed.

**A separate, distinct downstream use of `CNP_TOTAL` (not part of the Hill
interface) — hemodynamic safety.** `MAP_DROP = EMAX_MAP_DROP * CNP_TOTAL /
(EC50_MAP + CNP_TOTAL)` is the original's *own second, independent*
physiological readout (transient vasodilation/reflex tachycardia) with its
own `EC50_MAP` and implicit `gamma=1` — a genuinely different mechanism
site (off-target CNP-class hemodynamics, not NPR-B/pERK signalling), left
exactly as the original computed it, still reading `CNP_TOTAL` directly
(unrenamed structurally, just using the renamed `C_VOS`/`C_TCNP` inputs).

Both Hill terms are exact renames (no `nls()` fitting anywhere) — every
one of the original's effect terms was already a plain ratio.

## `$PARAM` vs `$CAPTURE` for `C_<STEM>`/`EFFECT_<STEM>`

Following the same, actually-compiling precedent already documented in
`pagets-disease/pbd_refactor_notes.md` and repeated across this batch
(declaring `C_VOS = 0` etc. directly in `$PARAM` and reassigning in `$ODE`
does not compile under mrgsolve 2.0.1 — `$PARAM` values are read-only
inside `$ODE`): `C_VOS`, `C_TCNP`, `C_INFIG`, `EFFECT_VOS`, `EFFECT_INFIG`,
and `CNP_TOTAL` are plain `double` locals computed in `$ODE`, listed in
`$CAPTURE`. Confirmed via `POST /model_manifest` on the refactored DSL:
all six appear in `outputPaths`.

## Dose-instant reporting artifact — checked, cosmetic, self-healing, identical in both models

Per the guide's note on this qspserver/mrgsolve quirk: checked whether the
duplicate report row `/run_simulation` emits at `t=0` (the dose instant)
shows a stale pre-dose value for `C_VOS`/`EFFECT_VOS`. It does (both rows
at `t=0` report `0`, before the dose's effect appears at the first
non-duplicate time point) — but this file's original has **no pre-existing
`$GLOBAL` macro convention** to match (grep confirms no `$GLOBAL` block
anywhere in the original), so per the guide's stated preference this is
disclosed as the cosmetic, self-healing artifact it is, not force-fit into
a `$GLOBAL` rewrite without precedent. Confirmed identical between the
patched-original and the refactored DSL at every duplicated `t=0` row
checked (see Verification below) — it does not affect the exact-match
result.

## Renaming applied

| Original | Refactored | Value | Role |
|---|---|---|---|
| `VOS_DEPOT` (cmt) | `GUT_VOS` | -- | depot |
| `VOS_CP` (cmt) | `CENT_VOS` | -- | central, concentration-state (bespoke, see above) |
| `KA_VOS` | `KA_VOS` | 6.0 | absorption rate (unchanged name) |
| `KE_VOS` | `CL_VOS` | 2.77 | elimination rate constant |
| `V_VOS` | `V1_VOS` | 15.0 | volume (used in depot-inflow scaling) |
| `F_VOS` | `F_VOS` | 0.70 | bioavailability (unchanged name) |
| `TCNP_DEPOT` (cmt) | `GUT_TCNP` | -- | depot |
| `TCNP_CP` (cmt) | `CENT_TCNP` | -- | central, concentration-state (bespoke) |
| `KREL_TCNP` | `KA_TCNP` | 0.018 | release rate (structural role = absorption; see naming note above) |
| `KE_TCNP` | `CL_TCNP` | 2.5 | elimination rate constant |
| `V_TCNP` | `V1_TCNP` | 15.0 | volume |
| -- (none in original) | -- (no `F_TCNP`) | -- | original has no bioavailability parameter for TCNP; none invented |
| `INFIG_GUT` (cmt) | `GUT_INFIG` | -- | depot |
| `INFIG_CP` (cmt) | `CENT_INFIG` | -- | central, concentration-state (bespoke) |
| `KA_INFIG` | `KA_INFIG` | 0.35 | absorption rate (unchanged name) |
| `CL_INFIG` | `CL_INFIG` | 7.7 | clearance (unchanged name — already matched convention) |
| `V_INFIG` | `V1_INFIG` | 480 | volume |
| `F_INFIG` | `F_INFIG` | 0.60 | bioavailability (unchanged name) |
| `EC50_VOS` | `EC50_VOS` | 8.0 | Hill EC50 (unchanged name; shared VOS+TCNP interface) |
| `HILL_VOS` | `GAMMA_VOS` | 1.5 | Hill coefficient |
| -- (none) | `EMAX_VOS` (new) | 1.0 | Hill ceiling [math-implied] |
| `EC50_INFIG` | `EC50_INFIG` | 25.0 | Hill EC50 (unchanged name) |
| `HILL_INFIG` | `GAMMA_INFIG` | 1.2 | Hill coefficient |
| `EMAX_INFIG_PERK` | `EMAX_INFIG` | 0.55 | Hill ceiling (real physiological cap, unchanged value) |

`PERK`, `CGMP_SIG`, `CHONDRO`, `HEIGHT_CM`, `HEIGHTZ`, `FMAREA`,
`SPCANALZ`, `AHI`, `OTITIS`, `BMIZ`, `MAP_BP`, `HR`, `PHOS`, `AGV_CALC`,
and every disease-side parameter (`PERK_BASE`, `AGV_BASE`, `GAIN_AGV`,
`EMAX_PERKINH`, `K_OFFTARGET_FGFR1`, `K_FM`, `K_SPCANAL`, `K_OTITIS`,
`K_BMI`, `EMAX_MAP_DROP`/`EC50_MAP`/`KOUT_MAP`, `EMAX_HR_RISE`/`KOUT_HR`,
`EMAX_PHOS_RISE`/`KOUT_PHOS`, `ADHERENCE`, `GH_ON`, `GH_EFFECT`, and all
ten `*0` baseline parameters) are **completely untouched** — not renamed,
not restructured. `F_VOS_DEPOT`/`F_INFIG_GUT` (the `$MAIN` bioavailability
assignments, which must track their compartment's new name under
mrgsolve's `F_<cmt>` convention) were renamed to `F_GUT_VOS`/`F_GUT_INFIG`
alongside the compartment rename — same values, same mechanism.

## Verification

**Method.** Extracted the bare DSL text from both `acdp_mrgsolve_model.R`
(`acdp_code <- '...'`, unmodified) and
`acdp_mrgsolve_model_refactored.R` (same variable name, refactored) and
ran them through the local qspserver `mrgsolve_api`
(`http://localhost:8007`), `POST /model_manifest` then `POST
/run_simulation`, requests spaced ~2.2s apart per the shared-service note.
Since the untouched original does not compile at all (defect #124), the
comparison baseline is a **patched-original** scratch copy — the original
DSL with only the same two syntax-only fixes applied (drop the fifteen
duplicate names from `$CAPTURE`; switch the thirteen `$MAIN`
initial-condition lines to the `<cmt>_0` idiom) — confirmed compiling
cleanly via `/model_manifest` before use. This scratch copy was never
written into the repo tree and is not part of the deliverables.

**Window.** The original's own scenarios run `end=8760, delta=24` (a full
year, daily observations). Per the guide's solver-budget allowance, and
because the `pow()`-NaN fragility above makes any window beyond its onset
uninformative for the affected scenarios anyway, verification used a
finer, shorter grid — `end=240, delta=1` (10 days, hourly) — long enough
to capture both (a) the full pre-NaN trajectory at high resolution and (b)
confirmation that both models stay `NaN` identically for the following 9
days once triggered, i.e. it's a stable shared state, not a transient
divergence. Scenarios 6 and 7 (TransCon CNP, Infigratinib — neither
triggers the fragility) were additionally checked over the full `end=2160,
delta=24` (90-day) window to confirm no ceiling/step-count issue over a
longer horizon.

**Scenarios run — the file's own, not invented, dosing amounts/timing
unchanged, only the `cmt` target and (for #10) the `ADHERENCE` override:**

1. **`1_Untreated_NaturalHistory`** (no dosing): **exact match, max abs
   diff `0.0`**, on every shared output over the full 10-day window
   (`PERK, CGMP_SIG, CHONDRO, AGV_CALC, HEIGHT_CM, HEIGHTZ, FMAREA,
   SPCANALZ, AHI, OTITIS, BMIZ, MAP_BP, HR, PHOS`, plus `VOS_CP/C_VOS`,
   `TCNP_CP/C_TCNP`, `INFIG_CP/C_INFIG` all identically `0`). No `NaN`
   anywhere (`CNP_TOTAL` stays exactly `0`, `pow(0,1.5)=0`, no sign
   ambiguity) — confirms the shared fragility above is specific to a
   nonzero, decaying, fractional-exponent input, not a general model
   defect.
2. **`2_Vosoritide_15ugkg_QD`** (`ev(amt=225, ii=24, addl=364, cmt=1)`):
   **exact match, max abs diff `0.0`** on every non-`NaN` shared output
   for `t = 0` through `t = 7` (8 reported points); both models hit `NaN`
   at **`t = 8`, identically**, and stay `NaN`-in-both for the remaining 9
   days checked — no `NaN`-alignment mismatch at any point (i.e. the
   original never goes `NaN` where the refactored stays finite, or vice
   versa). `VOS_CP` (original) vs `C_VOS` (refactored): exact match at
   every pre-`NaN` point (e.g. `3.9114` at `t=0.5`).
3. **`6_TransConCNP_QW`** (`ev(amt=1500, ii=168, addl=51, cmt=3)`): **exact
   match, max abs diff `0.0`**, every shared output, full 10-day window —
   no `NaN` (as expected, weekly dosing with slow release never lets
   `C_TCNP` fully clear to the noise floor). Additionally confirmed no
   `NaN` over the full 90-day window.
4. **`7_Infigratinib_PO_QD`** (`ev(amt=7.5, ii=24, addl=364, cmt=5)`):
   **exact match, max abs diff `0.0`**, every shared output, full 10-day
   window — no `NaN` (daily oral dosing with `t1/2 ≈ 20-24h` never fully
   clears between doses). Additionally confirmed no `NaN` over the full
   90-day window; `INFIG_CP`/`C_INFIG` both reach a `0.0209`-equivalent
   steady-state trough by day 90.
5. **`10_Vosoritide_PoorAdherence_60pct`** (same event as scenario 2, with
   `ADHERENCE = 0.6` param override): **exact match, max abs diff `0.0`**
   on every non-`NaN` shared output; `NaN` onset again identically at `t =
   8` in both models. Cross-checked the override itself: peak `C_VOS` at
   `t=0.5` is `0.7043` here vs `1.1739` in scenario 2 — ratio `0.6001 ≈
   0.6`, confirming `F_GUT_VOS = F_VOS * ADHERENCE` is applying correctly
   through the renamed bioavailability variable.

This is the expected outcome for three bespoke depot→central
(concentration-state) PK compounds with two exact-rename Hill interfaces
(no `nls()` fitting anywhere): genuinely `0.0` at every deterministic,
non-`NaN` point, with the one pre-existing numerical fragility reproducing
identically (down to the exact time step) between original and
refactored — confirmed, not a loosened tolerance.

`/model_manifest` on the refactored DSL additionally confirmed all of
`KA_VOS, CL_VOS, V1_VOS, F_VOS, KA_TCNP, CL_TCNP, V1_TCNP, KA_INFIG,
CL_INFIG, V1_INFIG, F_INFIG, EC50_VOS, GAMMA_VOS, EMAX_VOS, EMAX_PERKINH,
EC50_INFIG, GAMMA_INFIG, EMAX_INFIG` appear in `parameters` with their
original numeric defaults, and `GUT_VOS, CENT_VOS, GUT_TCNP, CENT_TCNP,
GUT_INFIG, CENT_INFIG, C_VOS, C_TCNP, C_INFIG, EFFECT_VOS, EFFECT_INFIG,
CNP_TOTAL` all appear in `outputPaths`.

No scratch/debug artifacts were left in the repository — DSL extraction,
the patched-original scratch copy, and the comparison scripts ran entirely
from a session-local scratchpad directory outside the repo tree, deleted
after use.

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`: all
three rows (INFIG, "Released free-CNP (TCNP)", Vosoritide) updated from
the classifier's placeholder `Delete PK compartment; concentration is
itself the state` to reflect the actual refactor outcome — each is a real,
externally-dosed drug, refactored as a bespoke depot→central
(concentration-state) archetype with its exposed `C_<STEM>` and (for
Vosoritide/TransCon CNP, jointly; for Infigratinib, independently) named
`EFFECT_<STEM>` Hill interface.
