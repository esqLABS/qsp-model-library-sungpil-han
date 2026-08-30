# Refactor notes — `takotsubo-syndrome/tts_mrgsolve_model.R`

Compounds refactored: **Apixaban (APX)**, **Carvedilol (CAR)**, **Dobutamine
(DOB)**, **Esmolol (ESM)**, **Furosemide (FUR)**, **Levosimendan (LEV)**,
**Metoprolol (MET)**, **Milrinone (MIL)**, **OR-1896 (OR)**, **Phenylephrine
(PHE)**, and **Ramipril/ramiprilat (RAM)** — the eleven rows in
`driver-patches/data/compound_perturbation_census.md` classified
`takotsubo-syndrome | APX/CAR/DOB/ESM/FUR/LEV/MET/MIL/OR/PHE/RAM`, all
"Redirect concentration inside #define macro". These are the file's entire
set of exogenous compounds. Every disease-side compartment and equation
(the three-segment beta-AR/Gs/Gi switch, cAMP/PKA/calcium, energetics and
injury, mechanics/LVOT, haemodynamics, neurohormonal axis, electrophysiology,
thrombus, recurrence, and the four counterfactual integrators) is
byte-for-byte identical to `tts_mrgsolve_model.R` — confirmed programmatically
by diffing the extracted DSL text token-by-token against the original outside
the touched lines. Only the eleven compounds' own PK compartments, their
exposed concentrations, and their effect terms were renamed and pulled out of
the `TTS_ALGEBRA` `#define` macro's un-renamed identifiers.

## Drug identity (per the task's request to confirm from code, not guess)

- **OR** = **OR-1896**, levosimendan's long-half-life active metabolite —
  confirmed from the original's own parameter names and comments (`V_OR`,
  `CL_OR` "OR-1896 clearance (L/h ; t1/2 ~77 h)", `FM_LEV` "fraction of
  levosimendan converted to OR-1896", `POT_OR` "OR-1896 potency relative to
  parent") and from its own dynamics (`dxdt_ORM = FM_LEV*(CL_LEV/V_LEV)*LEVC
  - (CL_OR/V_OR)*ORM` — a metabolite formed from levosimendan's own
  clearance, with no independent dose route anywhere in the file). **Not**
  an oral beta-blocker (the task brief flagged this as a possibility to
  check): OR-1896 acts through the same myofilament Ca²⁺-sensitisation
  mechanism as its parent, not the beta-adrenergic receptor, and it has no
  `KA_OR`/depot compartment (i.e. no oral route at all — it only ever
  exists as levosimendan's own metabolite, formed continuously during an IV
  levosimendan infusion).
- **APX** = Apixaban (oral factor Xa inhibitor, thromboprophylaxis).
- **CAR** = Carvedilol (oral non-selective beta-blocker + alpha1-blocker).
- **DOB** = Dobutamine (IV inotrope, beta1-selective partial agonist with
  minor beta2 activity).
- **ESM** = Esmolol (IV ultra-short-acting beta-blocker).
- **FUR** = Furosemide (oral/IV loop diuretic).
- **LEV** = Levosimendan (IV calcium sensitiser / K_ATP-channel opener).
- **MET** = Metoprolol (oral beta-blocker).
- **MIL** = Milrinone (IV PDE3 inhibitor inotrope).
- **PHE** = Phenylephrine (IV alpha1 agonist pressor).
- **RAM** = Ramipril, dosed as the prodrug and modelled as its active
  moiety ramiprilat (oral ACE inhibitor).

## Why "inside #define macro" needed more than a rename

All eleven compounds' concentrations were computed as assignments inside
`TTS_ALGEBRA`, a single `#define` macro expanded identically inside both
`$ODE` and `$TABLE` (the file's own stated reason: "a `$TABLE` block which
recomputes a flux with even slightly different gating reports a trajectory
that is not the one being integrated. A macro cannot drift"). Unlike the
`acute-pancreatitis`/`prurigo-nodularis` precedent files, this file's
authors had *already* worked around the classic `$ODE`-local + `$CAPTURE`
mrgsolve collision themselves: every quantity the macro touches (`C_DOB`,
`C_MIL`, `C_LEV`, `C_OR`, `C_PHE`, `C_RAM`, `C_FUR`, `C_APX`, `C_ESM`,
`C_MET`, `C_CAR`, plus what were `DOBEF`/`PDE3INH`/`CASENS`/`ACEINH`/
`FUREF`/`ANTIC`/`INH1`/`INH2`) is `double`-declared once at `$GLOBAL` scope,
and the macro only *assigns* them — never re-declares — so the same macro
text can expand in both `$ODE` and `$TABLE` without a duplicate-declaration
error. This refactor followed that same, already-working pattern for every
new identifier it introduces (see below), rather than moving anything to
`$ODE`-local scope.

Six of the eleven compounds' concentration variables were **already** named
exactly `C_<STEM>` per this fork's convention (`C_DOB`, `C_MIL`, `C_LEV`,
`C_OR`, `C_PHE`, `C_RAM`, `C_FUR`, `C_APX`, `C_ESM`, `C_MET`, `C_CAR` — all
eleven, in fact) — so, unlike `acute-pancreatitis`/`prurigo-nodularis`, no
"amount compartment was itself called `C_<STEM>`" collision existed here.
What still needed doing:

1. **Amount compartments renamed** to the `CENT_<STEM>`/`GUT_<STEM>`/
   `PERI_<STEM>` convention (`DOB`→`CENT_DOB`, `MIL`→`CENT_MIL`,
   `LEVC`/`LEVP`→`CENT_LEV`/`PERI_LEV`, `ORM`→`CENT_OR`, `ESM`→`CENT_ESM`,
   `METD`/`METC`→`GUT_MET`/`CENT_MET`, `CARD`/`CARC`→`GUT_CAR`/`CENT_CAR`,
   `PHE`→`CENT_PHE`, `RAMD`/`RAMC`→`GUT_RAM`/`CENT_RAM`, `FURD`/`FURC`→
   `GUT_FUR`/`CENT_FUR`, `APXD`/`APXC`→`GUT_APX`/`CENT_APX`), preserving the
   `$INIT` declaration order exactly (mrgsolve infers compartment numbers
   from `$INIT` order in this file — there is no separate `$CMT` block —
   so a pure 1-for-1, order-preserving rename keeps every compound's
   1-based compartment index identical to the original; confirmed directly
   via `/model_manifest`'s `outputPaths`, see Verification).
2. **Volume/clearance parameters renamed** to the `V1_<STEM>`/`V2_<STEM>`
   convention (`V_<STEM>`→`V1_<STEM>` for all eleven; `VP_LEV`→`V2_LEV`).
   `CL_<STEM>`, `KA_<STEM>`, `F_<STEM>`, `Q_LEV` were already convention-
   compliant and untouched.
3. **Effect terms pulled out of the macro's un-named arithmetic** into
   explicit `EFFECT_<STEM>` variables (see per-compound table below),
   declared at `$GLOBAL` next to the existing declarations, following the
   file's own established macro-safe pattern.
4. **`$CAPTURE`d**: none of `C_DOB`, `C_MIL`, `C_LEV`, `C_OR`, `C_ESM`,
   `C_MET`, `C_CAR`, `C_PHE`, `C_RAM`, `C_FUR`, `C_APX`, or any of the old
   `DOBEF`/`PDE3INH`/`INH1`/`INH2`/`INHA1` were in the original's
   `$CAPTURE` block (only `CASENS`, `ACEINH`, `ANTIC`, `FFLOSS` were, of
   which the first three are renamed below). Per the guide's
   qspserver-compatibility requirement 4, all eleven `C_<STEM>` plus every
   `EFFECT_<STEM>` (16 new capture entries in total, since CAR/ESM/MET each
   expose two or three named effect terms — see below) were added to
   `$CAPTURE`, confirmed discoverable via `/model_manifest`'s `outputPaths`
   (see Verification).

No compound in this file needed reordering `$ODE` (unlike
`acute-pancreatitis`, where the PK block sat after the disease block that
read it) — `TTS_ALGEBRA` already computes every compound's concentration
and effect *before* any disease-side algebra reads it, in one pass, and the
`dxdt_` PK lines already sit in their own contiguous block inside `$ODE`.

## Archetype per compound

- **Dobutamine (DOB)**: **archetype 1** (single compartment, linear
  elimination, no depot). `CENT_DOB` (was bare `DOB`), `CL_DOB`/`V1_DOB`
  (was `V_DOB`) unchanged in value. **Continuous zero-order infusion**
  (`RATE_DOB`, a `$PARAM` added directly into `dxdt_CENT_DOB`, not an
  `ev()`/dosing event) — the guide's edge case, preserved exactly; the
  original's own scenarios (S08, S20) drive it by splicing `RATE_DOB` on
  and off inside R, which this refactor does not change.
- **Milrinone (MIL)**: **archetype 1**. `CENT_MIL` (was `MIL`),
  `CL_MIL`/`V1_MIL` unchanged. Continuous zero-order infusion (`RATE_MIL`),
  same edge case.
- **Levosimendan (LEV)**: **archetype 2** (central + peripheral, no depot;
  `CL_LEV`/`V1_LEV`/`V2_LEV` (was `V_LEV`/`VP_LEV`)/`Q_LEV` unchanged).
  Continuous zero-order infusion (`RATE_LEV`) directly into `CENT_LEV`, same
  edge case. The original writes the CL and Q outflows as two separate
  terms (`- (CL_LEV/V_LEV)*LEVC - (Q_LEV/V_LEV)*LEVC + ...`) rather than
  the archetype template's single combined `-(CL+Q)/V1*CENT` term; this
  refactor keeps the original's own two-term decomposition verbatim
  (algebraically identical, a pure rename) rather than consolidating it, to
  minimise the risk of a transcription error in a heavily-load-bearing line.
- **OR-1896 (OR)**: **bespoke, archetype-1-like**. `CENT_OR` (was `ORM`) is
  a single compartment with linear elimination (`CL_OR`/`V1_OR`, was
  `V_OR`), but its *input* is not its own dose/depot — it is fed entirely by
  a fixed fraction (`FM_LEV`) of levosimendan's own clearance
  (`dxdt_CENT_OR = FM_LEV*(CL_LEV/V1_LEV)*CENT_LEV - (CL_OR/V1_OR)*CENT_OR`),
  exactly as the original models it. No `GUT_OR`/`KA_OR` exists because none
  existed in the original — this is disclosed as bespoke rather than forced
  into a depot-and-central shape, following the same precedent as
  `acne-vulgaris`'s OXO (4-oxo-isotretinoin, ISO's own active metabolite,
  fed from `KMET_ISO` rather than its own dose).
- **Esmolol (ESM)**: **archetype 1**. `CENT_ESM` (was `ESM`),
  `CL_ESM`/`V1_ESM` unchanged. Continuous zero-order infusion (`RATE_ESM`).
- **Metoprolol (MET)**: **archetype 3 minus peripheral** (depot + central,
  linear, oral). `GUT_MET`/`CENT_MET` (was `METD`/`METC`), `KA_MET`,
  `F_MET`, `CL_MET`/`V1_MET` (was `V_MET`) unchanged.
- **Carvedilol (CAR)**: **archetype 3 minus peripheral** (depot + central,
  linear, oral). `GUT_CAR`/`CENT_CAR` (was `CARD`/`CARC`), `KA_CAR`,
  `F_CAR`, `CL_CAR`/`V1_CAR` (was `V_CAR`) unchanged.
- **Phenylephrine (PHE)**: **archetype 1**. `CENT_PHE` (was bare `PHE` —
  see naming caution below), `CL_PHE`/`V1_PHE` (was `V_PHE`) unchanged.
  Continuous zero-order infusion (`RATE_PHE`).
- **Ramipril/ramiprilat (RAM)**: **archetype 3 minus peripheral** (depot +
  central, linear, oral). `GUT_RAM`/`CENT_RAM` (was `RAMD`/`RAMC`),
  `KA_RAM`, `F_RAM` (folds the prodrug-to-ramiprilat conversion efficiency
  into oral bioavailability, exactly as the original already does — no
  separate conversion compartment existed to preserve), `CL_RAM`/`V1_RAM`
  (was `V_RAM`) unchanged.
- **Furosemide (FUR)**: **archetype 3 minus peripheral** (depot + central,
  linear, oral). `GUT_FUR`/`CENT_FUR` (was `FURD`/`FURC`), `KA_FUR`,
  `F_FUR`, `CL_FUR`/`V1_FUR` (was `V_FUR`) unchanged.
- **Apixaban (APX)**: **archetype 3 minus peripheral** (depot + central,
  linear, oral). `GUT_APX`/`CENT_APX` (was `APXD`/`APXC`), `KA_APX`,
  `F_APX`, `CL_APX`/`V1_APX` (was `V_APX`) unchanged.

No compound needed TMDD (archetype 4) — there is no receptor-binding ODE
system anywhere in this file for any of the eleven; every beta/alpha
receptor interaction is an algebraic occupancy ratio, not an integrated
complex.

### A short-token rename hazard, avoided deliberately (`PHE`, `DOB`, `MIL`, `ESM`)

Several of this file's bare compartment names are short tokens that are
also substrings of their *own* parameter names (`PHE` inside `RATE_PHE`,
`C_PHE`, `V_PHE`, `CL_PHE`, `POT_PHE`; similarly `DOB`, `MIL`, `ESM` each
inside `RATE_<STEM>`, `C_<STEM>`, `V_<STEM>`, `CL_<STEM>`, and — for
DOB specifically — `EC50_DOB`/`EMAX_DOB`/`FB2_DOB`/`DOBEF`). A blind
find-and-replace of the bare token (`PHE`→`CENT_PHE`, etc.) across the
whole file would have corrupted every one of those parameter names into
nonsense (`V_PHE`→`V_CENT_PHE`, `C_DOB`→`C_CENT_DOB`, ...). This refactor
avoided that entirely by editing full, exact, contextually-unique
lines/expressions (the `$PARAM` block, the `$INIT` block, and each
specific `TTS_ALGEBRA`/`$ODE` line), never a blind bare-token substitution
— confirmed by grepping the delivered file for every old bare token
(`DOB`, `MIL`, `LEVC`, `ORM`, `PHE`, `ESM`, `METD`, `METC`, `CARD`, `CARC`,
`RAMD`, `RAMC`, `FURD`, `FURC`, `APXD`, `APXC`) and finding zero remaining
occurrences outside the intended renamed forms.

## Hill interface

### Rename, not a fit (six compounds: DOB, MIL, LEV, RAM, FUR, APX)

Each of these six had an effect term already shaped exactly
`Emax*C/(C+EC50)` (or, for four of them, an implicit-`Emax=1` plain ratio
`C/(C+EC50)`) inside the macro — a rename, per the guide's rule:

| Compound | Old name | New name | Notes |
|---|---|---|---|
| DOB | `DOBEF` (`EMAX_DOB*C_DOB/(C_DOB+EC50_DOB)`) | `EFFECT_DOB` | `EMAX_DOB`/`EC50_DOB` already named; `GAMMA_DOB=1` pulled out explicit |
| MIL | `PDE3INH` (`IMAX_MIL*C_MIL/(C_MIL+IC50_MIL)`) | `EFFECT_MIL` | `IMAX_MIL`→`EMAX_MIL`, `IC50_MIL`→`EC50_MIL`, `GAMMA_MIL=1` new |
| LEV (+OR) | `CASENS` (`EMAX_LEV*LEVEQ/(LEVEQ+EC50_LEV)`) | `EFFECT_LEV` (shared, see below) | `LEVEQ`→`X_LEV`; `EMAX_LEV`/`EC50_LEV` already named; `GAMMA_LEV=1` new |
| RAM | `ACEINH` (`IMAX_RAM*C_RAM/(C_RAM+IC50_RAM)`) | `EFFECT_RAM` | `IMAX_RAM`→`EMAX_RAM`, `IC50_RAM`→`EC50_RAM`, `GAMMA_RAM=1` new |
| FUR | `FUREF` (`C_FUR/(C_FUR+EC50_FUR)`) | `EFFECT_FUR` | implicit `Emax=1`→`EMAX_FUR=1` new, `EC50_FUR` already named, `GAMMA_FUR=1` new |
| APX | `ANTIC` (`C_APX/(C_APX+EC50_APX)`) | `EFFECT_APX` | implicit `Emax=1`→`EMAX_APX=1` new, `EC50_APX` already named, `GAMMA_APX=1` new |

Each `EFFECT_<STEM>` is now written with explicit `pow(X,GAMMA_<STEM>)` on
both numerator and denominator (`GAMMA_<STEM>=1` in every case, matching
the original's implicit shape exactly), and replaces its old macro-internal
name verbatim at every downstream call site (`dxdt_ANGII`, `dxdt_HREC`,
`dxdt_VOL`, `dxdt_KPL`, `dxdt_THR`, `dxdt_HEMB`, `CONT_AP/MD/BS`,
`CTN_AP/MD/BS`, `CO2_AP/MD/BS`, `CO3_AP/MD/BS`, `SVRT` — a straightforward
token rename, no coefficient touched).

### Disclosed bespoke #1 — shared metabolite Hill (LEV + OR)

`EFFECT_LEV` is genuinely **shared** between levosimendan and its own
active metabolite OR-1896, over a potency-weighted composite exposure
`X_LEV = C_LEV + POT_OR*C_OR` (was `LEVEQ`) — this is the *same* molecule's
own metabolite acting through the *same* mechanism (myofilament Ca²⁺
sensitisation), not two independent drugs collapsed into one term against
the guide's rule. No independent `EFFECT_OR` exists, because none existed
in the original (`CASENS` was already computed from the combined `LEVEQ`).
This exact situation — a parent and its active metabolite sharing one Hill
term over a composite exposure — has a direct precedent in this fork:
`acne-vulgaris`'s ISO/OXO row (`EFFECT_ISO` shared with OXO via
`X_ISO = C_ISO + POT_OXO*C_OXO`), disclosed there the same way.

### Disclosed bespoke #2 — linear competitive-occupancy terms (ESM, MET, CAR)

Esmolol, metoprolol, and carvedilol are three simultaneous competitive
antagonists at the same receptor pools. The original computes shared
receptor-availability fractions directly from all three drugs' raw
concentrations in one line each:

```
INH1 = 1.0 / (1.0 + C_ESM/KI1_ESM + C_MET/KI1_MET + C_CAR/KI1_CAR);   // beta1 available
INH2 = 1.0 / (1.0 + C_ESM/KI2_ESM + C_MET/KI2_MET + C_CAR/KI2_CAR);   // beta2 available
INHA1 = 1.0/(1.0 + C_CAR/KIA_CAR);                                     // alpha1 available (carvedilol only)
```

This is **not** a classic saturating `Emax*C/(C+EC50)` Hill ratio for any
one compound in isolation — each term is *linear* in that compound's own
concentration, and the saturation only emerges from the *sum* in the
shared denominator. Forcing an `EMAX_<STEM>`/`EC50_<STEM>`/`GAMMA_<STEM>`
Hill shape onto any one of these three in isolation would misrepresent
what the equation actually is. Per the guide ("keep each compound's
`EFFECT_<STEM>` separate; combine them only at the point the disease
equations actually use them"), each compound's own linear competitive term
was pulled out and named independently:

```
EFFECT_ESM_B1 = C_ESM/KI1_ESM;   EFFECT_ESM_B2 = C_ESM/KI2_ESM;
EFFECT_MET_B1 = C_MET/KI1_MET;   EFFECT_MET_B2 = C_MET/KI2_MET;
EFFECT_CAR_B1 = C_CAR/KI1_CAR;   EFFECT_CAR_B2 = C_CAR/KI2_CAR;   EFFECT_CAR_A1 = C_CAR/KIA_CAR;
INH1  = 1.0 / (1.0 + EFFECT_ESM_B1 + EFFECT_MET_B1 + EFFECT_CAR_B1);
INH2  = 1.0 / (1.0 + EFFECT_ESM_B2 + EFFECT_MET_B2 + EFFECT_CAR_B2);
INHA1 = 1.0 / (1.0 + EFFECT_CAR_A1);
```

`INH1`/`INH2`/`INHA1` themselves are left as-is (renamed nothing) — they
are genuinely combined, disease-side receptor-availability quantities, not
any one compound's own effect, and the guide's combine-at-point-of-use rule
is satisfied exactly by keeping them as the single place the three
independent `EFFECT_<STEM>` terms are summed. This is algebraically
identical to the original (`C_ESM/KI1_ESM` computed once either way), so
no numeric change results — confirmed by the multi-drug stress test below.

### Disclosed bespoke #3 — linear agonist-drive term (PHE)

Phenylephrine's contribution to alpha1 tone is `POT_PHE*C_PHE`, added
directly into the shared agonist-drive sum `AG_A1 = NE + AE_A1*EPI +
POT_PHE*C_PHE` — again linear in phenylephrine's own concentration, with
saturation occurring only downstream, in the shared `OA1 = INHA1*AG_A1/
(AG_A1+KD_A1)` occupancy-binding equation. This mirrors exactly how NE and
EPI already enter that same sum (also unscaled by any Emax/EC50 of their
own). Named `EFFECT_PHE = POT_PHE*C_PHE` and substituted at its one call
site — a rename of an already-linear term, not a fit, and not forced into
an Emax/EC50 shape it does not have.

## qspserver `/model_manifest` discoverability

Confirmed via `POST /model_manifest` against the extracted refactored DSL
(`http://localhost:8007`), diffed parameter-by-parameter against the
original's own manifest:

- Every renamed parameter (`V1_DOB`, `V1_MIL`, `V1_LEV`, `V2_LEV`, `V1_OR`,
  `V1_ESM`, `V1_MET`, `V1_CAR`, `V1_PHE`, `V1_RAM`, `V1_FUR`, `V1_APX`,
  `EC50_MIL`, `EMAX_MIL`, `EC50_RAM`, `EMAX_RAM`) carries **exactly** the
  original's own value under its new name — confirmed by an automated
  parameter-value diff across all 239 original parameters (0 mismatches).
- The only parameters present in the refactored manifest but absent from
  the original are the 8 explicitly-disclosed new ones: `GAMMA_DOB`,
  `GAMMA_MIL`, `GAMMA_LEV`, `GAMMA_RAM`, `GAMMA_FUR`, `GAMMA_APX`,
  `EMAX_FUR`, `EMAX_APX` — all `=1` (implicit in the original, made
  explicit per the guide), no other new parameter exists.
- All eleven `C_<STEM>` (`C_DOB`, `C_MIL`, `C_LEV`, `C_OR`, `C_ESM`,
  `C_MET`, `C_CAR`, `C_PHE`, `C_RAM`, `C_FUR`, `C_APX`) plus `X_LEV` and all
  sixteen `EFFECT_<STEM>` terms (`EFFECT_DOB`, `EFFECT_MIL`, `EFFECT_LEV`,
  `EFFECT_ESM_B1`, `EFFECT_ESM_B2`, `EFFECT_MET_B1`, `EFFECT_MET_B2`,
  `EFFECT_CAR_B1`, `EFFECT_CAR_B2`, `EFFECT_CAR_A1`, `EFFECT_PHE`,
  `EFFECT_RAM`, `EFFECT_FUR`, `EFFECT_APX`) are state-dependent (not
  covariates a client sets directly), so — per this fork's established
  precedent (`acute-pancreatitis`, `prurigo-nodularis`) — they cannot be
  `$PARAM` defaults; all appear in the manifest's `outputPaths` instead
  (128 `outputPaths` total: 64 `$CMT` states + 64 `$CAPTURE`d quantities,
  vs. the original's 64 `$CMT` + 40 `$CAPTURE`).
- Compartment order/index confirmed **identical** between original and
  refactored manifests at every one of the 64 positions (e.g. position 4
  is `DOB` in the original and `CENT_DOB` in the refactored, position 10 is
  `METD`/`GUT_MET`, etc.) — the pure order-preserving rename means every
  compound's own numeric `cmt=` dosing target is unchanged.
- No separate `.cpp` extraction was left in the repository; the quoted DSL
  block was extracted verbatim to a scratchpad file for every API call
  below, confirmed byte-identical (Python string equality, both 61,498
  characters) to the `tts_code_refactored <- r"---(...)---"` string in the
  delivered `tts_mrgsolve_model_refactored.R`, then the scratch copy was
  deleted after verification.

## Upstream build status

`tts_mrgsolve_model.R` **compiles cleanly** under mrgsolve 2.0.1 via the
qspserver API — no pre-existing build defect was found for this file (its
own author had already worked around the `$ODE`-local/`$CAPTURE` and
`$GLOBAL`-declare-once patterns that have tripped up other files in this
corpus; see "Why #define macro needed more than a rename" above). No
`UPSTREAM_ISSUES.md` entry was needed for this model.

## Verification

Per the guide's mandatory protocol, via the qspserver `mrgsolve_api`
service (`POST /model_manifest`, `POST /run_simulation`,
`http://localhost:8007`, requests spaced ~2.2 s apart), comparing the
untouched original against `tts_mrgsolve_model_refactored.R`'s extracted
DSL, using **each compound's own dosing scenario already defined in the
original file's own R code** (`scenarios`/`oral_events`), same dosing, same
default parameters (the "emo" phenotype's `AMP_TOT=250`/`FRAC_E=0.75`/
`TAU_SUR=1.6`/`E2=15`/`SEPT=0.35` are already this model's own `$PARAM`
defaults, so no phenotype override was needed to reproduce S02's trigger).
Continuous-infusion compounds (DOB, MIL, LEV, ESM, PHE) were driven through
the API using mrgsolve's own NM-TRAN-style zero-order-infusion `DoseSpec`
(`amt = rate*duration`, `rate = <the scenario's own RATE_<STEM> value>`,
`time = <the scenario's own onset>`) directly into each compound's own
compartment — mathematically identical to the R harness's own
`RATE_<STEM>` on/off splicing (`dxdt_CENT_<STEM> = rate - CL/V*CENT_<STEM>`
during the infusion window is exactly what a zero-order NM-TRAN infusion of
that same rate produces), not a different dosing scheme; this is a
verification-harness convenience, not a change to the model's own
`RATE_<STEM>`-based dosing convention (which is unchanged and still what
`tts_mrgsolve_model_refactored.R`'s own R scenarios use).

| Compound | Scenario reproduced | Window | Points | Max abs diff |
|---|---|---|---|---|
| DOB | S08 (dobutamine 5 ug/kg/min, h4–52) | 0–100 h | 101 | **0.0** |
| MIL | S10 (milrinone 0.5 ug/kg/min, h4–52, falsification-3) | 0–100 h | 101 | **0.0** |
| LEV (+OR) | S09 (levosimendan 0.1 ug/kg/min x24h from h4) | 0–100 h | 101 | **0.0** |
| ESM | S14 (esmolol, non-obstructive) | 0–100 h | 101 | **0.0** |
| PHE | S13 (phenylephrine 30 mg/h, h4–52) | 0–100 h | 101 | **0.0** |
| FUR | `oral_events$furosemide` (40 mg q12h x13, full schedule) | 0–180 h | 181 | **0.0** |
| MET | S27 (metoprolol 50 mg bd, **full** 365 d schedule, 730 doses) | 0–8760 h | 732 | **0.0** |
| CAR | `oral_events$carvedilol` (12.5 mg q12h, **full** 365 d schedule) | 0–8760 h | 732 | **0.0** |
| RAM | S26 (ramipril 5 mg daily, **full** 365 d schedule) | 0–8760 h | 732 | **0.0** |
| APX | S28 (apixaban 5 mg bd, **full** 90 d schedule) | 0–2160 h | 362 | **0.0** |

For every row, every compared output was exact (`LVEF`, `MAP`, `CO`, plus
compound-specific readouts — `O2_AP`/`O2_BS`/`GI_AP` for the beta-agonists,
`CONT_AP`/`CQ_AP` for MIL, `CENT_LEV`/`PERI_LEV`/`CENT_OR` states for
LEV/OR, `VOL`/`KPL` for FUR, `HR`/`O2_AP` for MET/CAR, `ANGII`/`ALDO`/`PREC`
for RAM, `THR`/`PEMB` for APX), at every reported time point, with **zero
shortening needed for any scenario** — none hit the API's default solver
step-count limit, including the four **full 365-day / 730-dose** oral
schedules run at their original length (this file's own `$SET
end=2160/delta=0.25` and `mcode_cache(... maxsteps=2000000)` build options
are R-side and do not carry over to the bare `model_content` compiled
server-side, but this system evidently does not need the larger budget at
these output resolutions — no `maxsteps`/solver error was ever raised).

**Additional stress test (not a named original scenario, but exercising the
disclosed bespoke ESM/MET/CAR decomposition together)**: esmolol continuous
infusion (S14's own rate) simultaneously with metoprolol's own oral
schedule, in the obstructive phenotype (`SEPT=1.0`, matching S15/S16's own
parameter), 0–240 h, 243 points — `LVEF`, `MAP`, `CO`, `GRAD`, `O2_AP`,
`O2_BS` all exact (max abs diff **0.0**), confirming the three-drug
competitive-occupancy decomposition reproduces the original's combined
`INH1`/`INH2` exactly under simultaneous multi-drug exposure, not just each
drug alone.

This is the expected outcome per the guide's tolerance rule: every one of
this file's eleven compounds is a pure structural reorganisation (rename +
macro-to-named-variable extraction, or a disclosed-but-still-exact linear
decomposition for ESM/MET/CAR/PHE, or a disclosed-but-still-exact shared
Hill for LEV/OR) — no ODE-derived Hill-fitting was needed anywhere in this
file, so exact (0.0) agreement is the correct result, not a loosened
tolerance.

## Anything else flagged

- No compound other than APX, CAR, DOB, ESM, FUR, LEV, MET, MIL, OR, PHE,
  RAM was touched.
- The R-side dosing helper `oral_events` (`ramipril`, `metoprolol`,
  `apixaban`, `furosemide`, `carvedilol`) was updated to dose into the
  renamed depot compartments (`GUT_RAM`, `GUT_MET`, `GUT_APX`, `GUT_FUR`,
  `GUT_CAR`) so the refactored sibling's own R scenario-driver code stays
  internally consistent — a rename of the dosing target only (same
  compartment, same 1-based index, same dose amounts/timing/intervals), not
  a behavioural change. `recurrence_1y()`'s one direct reference to the
  old capture name (`o$ACEINH`) was updated to `o$EFFECT_RAM`.
  `mcode_cache()`'s model name was changed to `"tts_qsp_refactored"` to
  avoid any compiled-model cache collision with the original.
- All scratch/debug files used for API verification (`.cpp` DSL
  extractions, transformation scripts, request/response JSON, the
  comparison harness) were created under the session scratchpad directory
  only and deleted after use — none were left in `takotsubo-syndrome/` or
  anywhere else in the repo.

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`, the
`takotsubo-syndrome | APX`, `takotsubo-syndrome | CAR`,
`takotsubo-syndrome | DOB`, `takotsubo-syndrome | ESM`,
`takotsubo-syndrome | FUR`, `takotsubo-syndrome | LEV`,
`takotsubo-syndrome | MET`, `takotsubo-syndrome | MIL`,
`takotsubo-syndrome | OR`, `takotsubo-syndrome | PHE`, and
`takotsubo-syndrome | RAM` rows (with OR's real identity, OR-1896, filled
in per the task's request).
