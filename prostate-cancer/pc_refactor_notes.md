# Refactor notes — `prostate-cancer/pc_mrgsolve_model.R`

## Scope of this pass — the census undercounted this file by 7 compounds

The census (`driver-patches/data/compound_perturbation_census.md`) had
exactly **one** row for this file: `AR Signaling (DEG)`, classified
"Redirect concentration (clean single site)". The task instructions
flagged this for a sanity check, following the exact precedent in
`neonatal-hyperbilirubinemia/nhb_refactor_notes.md` (where "Bilirubin
(SNMP)" and "Hepatic handling (PB)" were process-description phrases
paired with the *right* drug stem under the *wrong* display name).

Reading the actual code: this file contains **eight** distinct, real,
independently-dosed compounds, each with its own PK block and its own
dosing route — Leuprolide, Degarelix, Relugolix, Enzalutamide,
Abiraterone, Docetaxel, Olaparib, Denosumab. "DEG" is not a degradation-
rate parameter (the alternative the task instructions raised as a
possibility) — it is Degarelix's own stem, and `Deg_c` genuinely is "a
clean single concentration site": `Deg_blockade = Deg_Emax * Deg_c /
(Deg_c + Deg_EC50)`, feeding directly and only into `GnRH_total`. So the
census's stem was right and its display name was wrong, exactly the nhb
pattern — except here the classifier also missed **seven other real
compounds** in the same file entirely. This is disclosed as a severe,
file-specific undercount (the census's own limitations section already
warns detection rate varies materially across independent
reimplementations of the classifier), not something to skip past.

**All eight compounds were refactored** and all eight now have their own
row in the census (see the Census section at the end of this file).

## Real compound identities (per the task's two questions)

- (a) The real, externally-dosed drugs are: **Leuprolide** (GnRH agonist,
  7.5 mg IM monthly depot), **Degarelix** (GnRH antagonist, 240 mg SC
  loading), **Relugolix** (oral GnRH antagonist, 120 mg QD), **Enzalutamide**
  (oral AR pathway inhibitor, 160 mg QD), **Abiraterone** (oral CYP17A1
  inhibitor, 1000 mg QD), **Docetaxel** (IV chemotherapy, 75 mg/m² q3w),
  **Olaparib** (oral PARP inhibitor, 300 mg BID), **Denosumab** (SC RANKL
  antibody, 120 mg q4w).
- (b) "DEG" refers to **Degarelix**, not a degradation-rate parameter. The
  file does separately have many genuine degradation-rate parameters named
  with a `_deg` suffix (`kLH_deg`, `kT_deg`, `kDHT_deg`, `kAR_deg`,
  `kAR_nuc_deg`, `kPSA_deg`, `kOC_deg`, `kOB_deg` — none of these are
  compounds), which is plausibly what made the census's classifier
  mis-derive "AR Signaling" as a display label (`kAR_deg` sits inside the
  `// ---- AR Signaling ----` `$PARAM` section) even though the actual
  concentration site it flagged, `Deg_c`, belongs to a completely different
  block (`// ---- Degarelix PK (240 mg SC loading) ----`) with nothing to
  do with AR signaling at all.

## Archetype per compound

All eight are **Archetype 3 minus peripheral** (depot + central, linear,
no peripheral compartment) **except Docetaxel**, which is bespoke:

- **LEUP, DEG, REL, ENZ, ABI, OLA, DEN** (7 of 8): depot compartment with
  first-order release/absorption into a central compartment with
  first-order elimination. Faithful to a template already used elsewhere
  in this fork (e.g. LETRO/OLAP in `breast-cancer/bc_refactor_notes.md`).
  REL, ENZ, ABI, OLA have an explicit oral bioavailability (`F_<STEM>`);
  LEUP, DEG, DEN do not (the original never declared one for these three
  either — an IM/SC depot releasing 1:1 into the central compartment).
- **DOC** (Docetaxel): **bespoke 2-compartment**, kept exactly as the
  original built it — `dxdt_CENT_DOC = -(KE1_DOC+K12_DOC)*CENT_DOC +
  K21_DOC*PERI_DOC`, `dxdt_PERI_DOC = K12_DOC*CENT_DOC -
  (K21_DOC+KE2_DOC)*PERI_DOC`. This does **not** fit the guide's Archetype
  2 template (which eliminates only from the central compartment) because
  the original's own peripheral compartment *also* eliminates
  (`KE2_DOC`, was `kDoc_elim2`) — a genuine dual-elimination biexponential
  design, not a textbook CL/Q/V system. Per the guide's explicit allowance
  ("if a compound's PK genuinely doesn't resemble any archetype... rename
  to the convention... note why"), the micro-constants were renamed in
  place (`KE1_DOC`/`K12_DOC`/`K21_DOC`/`KE2_DOC`) rather than forced into
  CL/Q/V. Docetaxel is also dosed directly into `CENT_DOC` as an IV bolus
  in already-converted µM units (`mg_to_uM_docetaxel()` in the R script,
  unchanged) — `CENT_DOC` is a concentration-state, not an amount depot,
  exactly matching the original's own `Doc_c` dosing convention (the "edge
  case" the guide anticipates: preserve non-event-style dosing as-is).

## The Leuprolide flare term is a separate PD sub-effect, not folded into EFFECT_LEUP

The original computes `GnRH_agonist_effect = (1 + GnRH_flare*Flare_eff) *
(1 - 0.97*Leup_suppress)` — two independent multiplicative terms: a
transient agonist "flare" (its own state, `Flare_eff`/`FLARE_LEUP`) and a
concentration-driven desensitization ratio. Only the second term is a
plain Hill ratio in `C_LEUP`; the first depends on a *different* state
variable entirely. Per the guide ("keep each compound's `EFFECT_<STEM>`
separate; combine them only at the point disease equations actually use
them"), `EFFECT_LEUP = EMAX_LEUP*C_LEUP/(EC50_LEUP+C_LEUP)` captures only
the concentration-dependent desensitization; the flare multiplier is kept
as a separate disclosed factor in `GnRH_agonist_effect =
(1+FLARE_MULT_LEUP*FLARE_LEUP)*(1-EFFECT_LEUP)`. `EMAX_LEUP=0.97` and
`EC50_LEUP=2.0` are new named parameters (the original hardcoded these two
literals directly in the ternary expression; comment already stated
"EC50=2 ng/mL").

**The flare mechanism is dead code in every one of the file's own 7
scenarios.** `Flare_eff`/`FLARE_LEUP` starts at 0 and is never dosed or set
nonzero by any scenario — the first draft of Scenario 2's event table
(`sc2_events <- ev(...) %>% mutate(Flare_eff = ifelse(time==0, 1.0, 0.0))`)
is immediately overwritten by a second, simplified `sc2_events <-
ev(data.frame(...))` two lines later that drops the flare entirely (see
lines ~519–536 of both the original and this refactored file — preserved
exactly, unchanged, including the overwrite). `GnRH_flare`/`flare_decay`
(now `FLARE_MULT_LEUP`/`FLARE_DECAY_LEUP`) and the `FLARE_LEUP` state
therefore have zero effect on any of the 7 shipped scenarios. Kept exactly
as the original had it (renamed only) — this is a pre-existing modelling
quirk, not something introduced or fixed by this refactor.

## Olaparib: HRR_def gates the effect, is not part of its own Hill term

`Ola_kill = HRR_def * Ola_Emax * Ola_c/(Ola_c+Ola_EC50)` combines
Olaparib's own concentration-driven Hill ratio with a *disease/genotype
covariate* (`HRR_def`, 0/1, BRCA2-mutation status) that has nothing to do
with Olaparib's own PK. `EFFECT_OLA` captures only the Hill ratio;
`HRR_def * EFFECT_OLA` is written explicitly at the one place it's used
(`k_death_eff`), preserving the original's exact combination.

## Bespoke deviation, disclosed: Degarelix's own volume was wrong in the original

`$PARAM` declares `V_Deg = 1000.0` ("Volume of distribution (L)"), but
`dxdt_Deg_c = kDeg_abs*Deg_sc/V_Den - kDeg_elim*Deg_c` divides by `V_Den`
(Denosumab's own volume, 3.0 L) instead — `V_Deg` is declared but **never
referenced anywhere in `$ODE`**. This is an apparent copy-paste error in
the original (the file has a dozen similarly-named `V_*` parameters
declared close together; confirmed by grepping the whole file for `V_Deg`
— it appears exactly twice, both in `$PARAM`/`$CMT` comments, never in
`$ODE`). The practical effect: Degarelix's simulated plasma concentration
is ~333x higher than a correctly-parameterized 1000 L volume would give,
and Degarelix's kinetics are silently tied to whatever `V_Den` happens to
be at run time (Denosumab's own volume).

Per "log what you find, don't fix it upstream" and "parameter values
always come from the original," this is preserved **numerically** —
`V1_DEG = 3.0`, the value actually in effect, not the declared-but-dead
`1000.0` — while being fixed **structurally**: Degarelix now has its own
independent `V1_DEG` parameter, no longer silently sharing Denosumab's
`V1_DEN` at run time. This satisfies the guide's "own compartments... not
interleaved" plumbing requirement (an override to one compound's volume no
longer silently perturbs the other) while reproducing the checked-in
original's actual trajectory exactly. Logged as
[`UPSTREAM_ISSUES.md` #67](../translations/UPSTREAM_ISSUES.md).

## Enzalutamide/Olaparib concentration unit label vs. actual computation

Both `Enz_c` and `Ola_c` are commented/labelled "(µM)" but their `$ODE`
equations never apply a molecular-weight conversion (unlike Abiraterone,
which explicitly divides by `Abi_MW`) — they compute a plain mg/L number
and call it µM. Preserved exactly (renamed only, `C_ENZ`/`C_OLA`, both
`$CAPTURE`-annotated "(uM, nominal)" to flag this for a reader). Not
fixed — this is a labelling/units inconsistency in the original, not a
build defect, and changing it would change the compound's numeric
potency relative to its own `EC50`.

## Renaming applied (values unchanged from the original except where noted)

| Original | Refactored | Value | Role |
|---|---|---|---|
| `Leup_depot`/`Leup_c`/`Flare_eff` (cmt) | `GUT_LEUP`/`CENT_LEUP`/`FLARE_LEUP` | -- | LEUP depot/central/flare |
| `kLeup_rel`/`kLeup_elim`/`V_Leup` | `KA_LEUP`/`KE_LEUP`/`V1_LEUP` | 0.033/0.693/40.0 | LEUP PK |
| `GnRH_flare`/`flare_decay` | `FLARE_MULT_LEUP`/`FLARE_DECAY_LEUP` | 3.0/0.5 | LEUP flare (dead code, see above) |
| (hardcoded `0.97`, `2.0`) | `EMAX_LEUP`/`EC50_LEUP` (new) | 0.97/2.0 | LEUP Hill |
| -- | `GAMMA_LEUP` (new) | 1.0 | LEUP Hill exponent |
| `Deg_sc`/`Deg_c` (cmt) | `GUT_DEG`/`CENT_DEG` | -- | DEG depot/central |
| `kDeg_abs`/`kDeg_elim` | `KA_DEG`/`KE_DEG` | 0.15/0.023 | DEG PK |
| `V_Deg` (dead) / actual `V_Den` | `V1_DEG` | 3.0 | DEG volume — see bespoke note above |
| `Deg_Emax`/`Deg_EC50` | `EMAX_DEG`/`EC50_DEG` | 0.98/0.001 | DEG Hill |
| -- | `GAMMA_DEG` (new) | 1.0 | DEG Hill exponent |
| `Rel_gut`/`Rel_c` (cmt) | `GUT_REL`/`CENT_REL` | -- | REL depot/central |
| `kRel_abs`/`F_Rel`/`kRel_elim`/`V_Rel` | `KA_REL`/`F_REL`/`KE_REL`/`V1_REL` | 1.4/0.12/1.0/2800.0 | REL PK |
| (hardcoded `0.98`) / `Rel_EC50` | `EMAX_REL` (new)/`EC50_REL` | 0.98/0.005 | REL Hill |
| -- | `GAMMA_REL` (new) | 1.0 | REL Hill exponent |
| `Enz_gut`/`Enz_c` (cmt) | `GUT_ENZ`/`CENT_ENZ` | -- | ENZ depot/central |
| `kEnz_abs`/`F_Enz`/`kEnz_elim`/`V_Enz` | `KA_ENZ`/`F_ENZ`/`KE_ENZ`/`V1_ENZ` | 1.5/0.84/0.114/110.0 | ENZ PK |
| `Enz_Emax`/`Enz_EC50` | `EMAX_ENZ`/`EC50_ENZ` | 0.95/3.0 | ENZ Hill |
| -- | `GAMMA_ENZ` (new) | 1.0 | ENZ Hill exponent |
| `Abi_gut`/`Abi_c` (cmt) | `GUT_ABI`/`CENT_ABI` | -- | ABI depot/central |
| `kAbi_abs`/`F_Abi`/`kAbi_elim`/`V_Abi` | `KA_ABI`/`F_ABI`/`KE_ABI`/`V1_ABI` | 0.8/0.10/1.7/19669.0 | ABI PK |
| `Abi_MW` | `MW_ABI` | 391.6 | ABI molecular weight |
| `Abi_Emax`/`Abi_EC50` | `EMAX_ABI`/`EC50_ABI` | 0.95/0.05 | ABI Hill |
| -- | `GAMMA_ABI` (new) | 1.0 | ABI Hill exponent |
| `Doc_c`/`Doc_p` (cmt) | `CENT_DOC`/`PERI_DOC` | -- | DOC central/peripheral (bespoke) |
| `kDoc_elim1`/`kDoc_k12`/`kDoc_k21`/`kDoc_elim2` | `KE1_DOC`/`K12_DOC`/`K21_DOC`/`KE2_DOC` | 3.94/1.5/0.8/0.231 | DOC PK (bespoke) |
| `Doc_Emax`/`Doc_EC50` | `EMAX_DOC`/`EC50_DOC` | 0.90/0.05 | DOC Hill |
| -- | `GAMMA_DOC` (new) | 1.0 | DOC Hill exponent |
| `Ola_gut`/`Ola_c` (cmt) | `GUT_OLA`/`CENT_OLA` | -- | OLA depot/central |
| `kOla_abs`/`F_Ola`/`kOla_elim`/`V_Ola` | `KA_OLA`/`F_OLA`/`KE_OLA`/`V1_OLA` | 1.4/0.66/1.6/167.0 | OLA PK |
| `Ola_Emax`/`Ola_EC50` | `EMAX_OLA`/`EC50_OLA` | 0.80/0.1 | OLA Hill |
| -- | `GAMMA_OLA` (new) | 1.0 | OLA Hill exponent |
| `Den_sc`/`Den_c` (cmt) | `GUT_DEN`/`CENT_DEN` | -- | DEN depot/central |
| `kDen_abs`/`kDen_elim`/`V_Den` | `KA_DEN`/`KE_DEN`/`V1_DEN` | 0.062/0.023/3.0 | DEN PK |
| `Den_Emax`/`Den_EC50` | `EMAX_DEN`/`EC50_DEN` | 0.95/0.5 | DEN Hill |
| -- | `GAMMA_DEN` (new) | 1.0 | DEN Hill exponent |
| `Leup_suppress` (partial)/`Deg_blockade`/`Rel_blockade`/`Enz_AR_inh`/`Abi_T_inh`/`Doc_kill`/`Ola_kill` (minus `HRR_def`)/`Den_RANKL_inh` | `EFFECT_LEUP`/`EFFECT_DEG`/`EFFECT_REL`/`EFFECT_ENZ`/`EFFECT_ABI`/`EFFECT_DOC`/`EFFECT_OLA`/`EFFECT_DEN` | -- | the 8 named Hill effects |

`HRR_def`, `PTEN_loss`, and all disease-side parameters/compartments
(HPG axis, AR signaling, tumor kinetics, PI3K/AKT, bone metastasis) are
unchanged, same names, same values.

## Fix applied for the original's own build defects (disclosed, non-numeric)

The original does not compile under mrgsolve 2.0.1 at all, for four
reasons unrelated to any compound's own archetype (three fully generic,
one specific to Leuprolide's own effect block) — full detail in
[`UPSTREAM_ISSUES.md` #67](../translations/UPSTREAM_ISSUES.md):

1. 16 of 33 `$INIT @annotated` lines were missing the required third
   (description) field — `parse annotated init block (INIT)` error.
2. Once fixed, `$CMT` + `$INIT` jointly redeclare all 33 compartments —
   `Duplicated model names` error (same family as #61/#63/#64/#65/#66).
3. `$CAPTURE` duplicated all 18 of its own entries against `$CMT` —
   `compartment should not be in $CAPTURE` error.
4. `self.trt_leup` (undefined `databox` member) inside Leuprolide's own
   GnRH-agonist-effect guard — `class databox has no member named
   trt_leup`.

Fixes 1–3 are applied directly to the delivered
`pc_mrgsolve_model_refactored.R`, syntax-only, non-numeric: the incomplete
`$INIT` lines (and then the whole `$INIT` block) are replaced by
`<CMT>_0 = value;` assignments in `$MAIN` (identical values, identical
compartments); `$CAPTURE` no longer repeats any `$CMT` name (all 16 of its
entries are now the new `C_<STEM>`/`EFFECT_<STEM>` derived quantities —
mrgsolve reports every compartment's own state via `/model_manifest`'s
`outputPaths` regardless of `$CAPTURE`, confirmed).

Fix 4 is treated as an **in-scope Leuprolide design decision**, not a
generic defect (per the guide's point 4: "if a defect is genuinely inside
the scope compound's own block... fix it as part of the refactor itself
and say so"): the ternary
`(NEWIND <= 1 || self.trt_leup == 0) ? 1.0 : (formula)` is replaced by the
unconditional formula. Because the file never compiled as originally
written, there is no prior runtime behaviour for this branch to preserve —
but the unconditional formula is (i) what the comment directly above it
describes ("leuprolide: flare then desensitize"), and (ii) numerically
neutral whenever Leuprolide is absent, since `Leup_c = Flare_eff = 0`
already makes the formula evaluate to exactly `1.0`, the same value the
guarded branch would have returned. Confirmed by this file's own untreated
scenario (below): LH/T converge to the same values whether or not the
guard is present, since no Leuprolide is dosed.

## A design pitfall found and fixed during this refactor (not an upstream defect)

The first draft computed `C_<STEM>`/`EFFECT_<STEM>` and all disease-facing
derived quantities (`GnRH_total`, `AKT_ss`, `AR_nuc_eff`, `k_death_eff`,
`RANKL_eff`, etc.) inside **`$ODE`** rather than `$MAIN`, on the reasoning
that `$ODE` is evaluated continuously by the solver and would avoid any
capture-timing artifact. This was wrong for this file specifically: the
**original's own** analogous quantities (`Leup_suppress`, `Enz_AR_inh`,
`Doc_kill`, `AKT_ss`, etc.) are computed inside **`$MAIN`**, which mrgsolve
evaluates once per output record, *before* that record's own `$ODE`
integration runs — making them piecewise-constant across each output
interval, not continuously updated. Verified empirically: with the
derived-quantity block moved into `$ODE`, Scenario 2 (ADT alone) showed a
**substantial**, non-floating-point divergence in LH at day 1 (2.73
absolute, ~2.7/5 ≈ 54% of the baseline value) that shrank back down only
after ~7 days — because the refactored model was now applying Leuprolide's
suppressive effect on GnRH continuously, while the original applies it
only once per day (per output interval), a real behavioural difference,
not a display artifact. Moving the block back into `$MAIN` (matching the
original's own placement exactly) reproduced the original's actual
once-per-interval update cadence and eliminated this divergence (see
Verification below). This is disclosed here as a genuine modelling
subtlety of *this file specifically* (not present in, e.g., `nhb` or `bc`,
whose analogous derived quantities were already computed in `$ODE` in
their own originals) — reproducing the original's own block placement, not
"improving" on it, is what the guide's verification mandate requires.

**Consequence for `$CAPTURE`d `C_<STEM>`/`EFFECT_<STEM>` outputs.**
Because these are computed in `$MAIN` (once per output record, before that
record's own integration), a captured `C_<STEM>` value is **one output
interval "behind"** the true current `CENT_<STEM>` state at that moment —
confirmed exactly: `C_LEUP[i] == Leup_c_original[i-1]` at every index,
bit-for-bit. This is not a numerical bug in the rename; it is an authentic
reproduction of the original's own behaviour (the original's `Leup_c` is a
genuine `$CMT`, always current, but its *own* downstream effect on GnRH
was *already* one-interval-delayed in the checked-in original, since
`Leup_suppress` is computed the same way, in `$MAIN`, from that same
state). The fair, apples-to-apples comparison used below is therefore
`CENT_<STEM>` (a raw, unlagged `$CMT` state, captured directly) divided by
its own volume outside the model, compared against the original's raw
`$CMT` concentration — which matches exactly (see table below) — rather
than the `$MAIN`-computed `C_<STEM>` capture, which is provided purely for
qspserver discoverability (per the compat requirement) and disclosed here
as one-interval-delayed by design.

## Verification

**Method.** Both the original's own model code (with the four build-compat
fixes above applied to an in-memory-only scratch copy, so it would compile
at all — never applied to the checked-in `pc_mrgsolve_model.R`) and
`pc_mrgsolve_model_refactored.R`'s embedded DSL were extracted as bare
mrgsolve DSL text and run through the qspserver `mrgsolve_api` container
(`POST /model_manifest`, `POST /run_simulation`) at `http://localhost:8007`,
requests spaced ~2s apart, never more than one in flight (respecting the
service's `max_concurrent_jobs: 2` limit and its history of crashing under
concurrent load). The extracted refactored DSL was confirmed byte-identical
to the `_refactored.R`'s embedded quoted string except for a small number
of possessive apostrophes in comments, which cannot appear inside an R
single-quoted string literal and were rewritten without changing meaning
(`"the original's"` → `"the originals"`); no arithmetic, parameter, or
compartment content differs. Compartment order (and therefore every
1-based dosing index) is **identical** between the two models — the
refactor renamed every compartment in place without adding, removing, or
reordering any of the 33 — confirmed from `/model_manifest`'s
`outputPaths` (e.g. `GUT_LEUP`=17, `CENT_DEG`=21, `CENT_DOC`=28,
`CENT_DEN`=33 in both).

**Scenarios run.** All 7 of the original's own scenarios (dose amounts,
`ii`/timing, and the `HRR_def=1.0` override for Scenario 6 all copied
verbatim from the original's own R code), full 3-year horizon (`end=1095,
delta=1`, no shortening needed — none of the scenarios approached the
API's default `maxsteps` budget), plus 3 additional single-dose checks
constructed for Degarelix, Relugolix, and Denosumab — **none of the file's
own 7 scenarios doses these three drugs at all** (only Leuprolide is ever
used as the "ADT backbone"; Degarelix/Relugolix/Denosumab have full PK
blocks but are never exercised by any shipped scenario), so a scenario had
to be constructed to verify their rename at all. Each constructed scenario
uses that compound's own declared dose amount (240 mg SC bolus for
Degarelix; 120 mg PO QD for Relugolix; 120 mg SC q4w x4 for Denosumab)
over a shorter window (60–180 days) — disclosed here as deviating from
"the file's own scenarios" only because there is nothing else to run for
these three compounds.

For each scenario, all 16 disease-side outputs shared between the two
models (`LH`, `T`, `DHT`, `AR_free`, `AR_DHT`, `AR_nuc`, `PSA`, `TC_p`,
`TC_q`, `CRPC_frac`, `ARv7_frac`, `AKT_act`, `OC`, `OB`, `BMD`,
`BoneMets`) were compared point-by-point, plus each scenario's own
dosed compound's raw PK (`GUT_<STEM>` vs. the original's depot compartment,
and `CENT_<STEM>`/`V1_<STEM>` reconstructed vs. the original's raw central
concentration `$CMT` — see the note above on why the `$MAIN`-computed
`C_<STEM>` capture is compared this way rather than directly).

| Scenario | Dosing | n pts | Shared-output max abs diff | PK reconstruction max abs diff |
|---|---|---|---|---|
| S1 Untreated | none | 1096 | 2.1e-3 (`TC_p`, scale ~420) | — |
| S2 ADT (Leuprolide) | Leup 7.5mg q28d x39 | 1135 | 1.97e-2 (`LH`, scale ~2-8) | `GUT_LEUP` 0.0 exact; `CENT_LEUP/V1*1000` vs `Leup_c` 1.2e-3 |
| S3 ADT+Enzalutamide | + Enz 160mg QD x1095 | 2230 | **0.0 exact** | `GUT_ENZ` 0.0 exact; recon vs `Enz_c` 5.0e-5 |
| S4 ADT+Abiraterone | + Abi 1000mg QD x1095 | 2230 | **0.0 exact** | `GUT_ABI` 0.0 exact; recon vs `Abi_c` 4.7e-5 |
| S5 ADT+Docetaxel x6 | + Doc IV q21d x6 | 1141 | 3.1e-2 (`LH`) | `CENT_DOC` vs `Doc_c` ~1e-91 (bit-exact); `PERI_DOC` vs `Doc_p` ~1e-59 (bit-exact) |
| S6 ADT+Olaparib (HRR_def=1) | + Ola 300mg BID x2190 | 3325 | **0.0 exact** | `GUT_OLA` 1.0e-4; recon vs `Ola_c` 4.6e-5 |
| S7 Sequential ADT→ARPI→Doc | Leup + Enz d180-449 + Doc d450+ | 1411 | 5.97e-2 (`TC_p`, scale ~400-1000) | — |
| Xc Degarelix (constructed) | 240mg SC bolus, 180d | 182 | 1.0e-4 (`TC_p`) | `GUT_DEG` 0.0 exact; recon vs `Deg_c` 3.3e-5 |
| Xc Relugolix (constructed) | 120mg PO QD, 60d | 121 | **0.0 exact** | `GUT_REL` 0.0 exact; recon vs `Rel_c` 1.72e-2 (scale: `C_REL` reaches several hundred ng/mL) |
| Xc Denosumab (constructed) | 120mg SC q4w x4, 180d | 188 | 1.0e-4 (`TC_p`) | `GUT_DEN` 0.0 exact; recon vs `Den_c` 3.3e-5 |

**Result: floating-point-scale match on every scenario** — 3 of 10 are
bit-exact (`0.0`) on every shared output; the remaining 7 range from
`1.0e-4` to `5.97e-2` absolute, all against output magnitudes of order
1–1000 (worst-case relative deviation ≈0.5% on `LH`, S2). This is
consistent with the guide's tolerance for Archetypes 1–3 ("pure structural
reorganization... expect a near-exact match... anything beyond
floating-point-scale deviation means a bug") and with this fork's own
documented precedent for how such reorganization compounds through an
adaptive-step ODE solver over long horizons (see
`essential-thrombocythemia/et_refactor_notes.md` and
`translations/UPSTREAM_ISSUES.md`'s entry for that file): the same
algebra, differently ordered, produces different floating-point rounding
at each solver step, which very slowly compounds over hundreds of days —
S1 (fully untreated, no compound active at all) already shows a
non-zero 2.1e-3 baseline from this effect alone, purely from renaming
variables with no behavioural change whatsoever. No scenario showed a
qualitative divergence (wrong sign, wrong magnitude, or unbounded growth);
every raw PK compartment (`GUT_<STEM>`, and `CENT_<STEM>` reconstructed to
a concentration) matches to the same floating-point-scale precision.

**No Hill fitting was needed for any of the 8 compounds** — every one of
the original's effect terms was already a plain concentration ratio
(`C/(C+EC50)`, `Emax*C/(EC50+C)`, or their algebraic complement), so every
`EFFECT_<STEM>` is a rename, not an approximation, consistent with the
exact/near-exact verification results above.

## qspserver `/model_manifest` discoverability

Confirmed via `POST /model_manifest` against the extracted, renamed DSL:
all 104 `$PARAM` entries are listed, including every `KA_<STEM>`,
`KE_<STEM>`/`CL_<STEM>`-equivalent, `F_<STEM>`, `V1_<STEM>`,
`EMAX_<STEM>`, `EC50_<STEM>`, `GAMMA_<STEM>` for all 8 compounds (e.g.
`KA_LEUP, KE_LEUP, V1_LEUP, FLARE_MULT_LEUP, FLARE_DECAY_LEUP, EMAX_LEUP,
EC50_LEUP, GAMMA_LEUP`; `KA_DEG, KE_DEG, V1_DEG, EMAX_DEG, EC50_DEG,
GAMMA_DEG`; and so on for REL/ENZ/ABI/DOC/OLA/DEN). `C_<STEM>` and
`EFFECT_<STEM>` (all 16) are state-derived (computed in `$MAIN` from
compartment values — see the design-pitfall note above for why `$MAIN` and
not `$ODE`), so — per the same reasoning already established in
`breast-cancer/bc_refactor_notes.md` and
`neonatal-hyperbilirubinemia/nhb_refactor_notes.md` (mrgsolve 2.0.1
compiles `$PARAM` members as read-only references inside `$ODE`/`$MAIN`,
so a value recomputed every timestep from state cannot also live in
`$PARAM`) — they cannot also be `$PARAM` entries; all 16 appear in the
manifest's `outputPaths` via `$CAPTURE`, confirmed discoverable. All 33
compartments (`GUT_LEUP`, `CENT_LEUP`, `FLARE_LEUP`, ... `GUT_DEN`,
`CENT_DEN`) also appear in `outputPaths` as ordinary compartments, in the
same order/index as the original.

No `.cpp` extraction file was left behind — extraction was in-memory only,
used to build the verification requests above and then discarded, per the
workflow guide.

## Anything else flagged

- The R-side scenario definitions (`sc2_events`...`sc7_events`,
  `enz_daily`, `abi_daily`, `doc_cycles`, `ola_bid`, `enz_from180`,
  `doc_from450`) were updated **only** where they name a renamed
  compartment inside an `ev()` call (`"Leup_depot"→"GUT_LEUP"`,
  `"Enz_gut"→"GUT_ENZ"`, `"Abi_gut"→"GUT_ABI"`, `"Doc_c"→"CENT_DOC"`,
  `"Ola_gut"→"GUT_OLA"`). Every dose amount, `ii`, `time` sequence, and
  `param_override` (`HRR_def=1.0` for Scenario 6, `PTEN_loss`/`HRR_def` in
  the sensitivity grid) is identical to the original. No plot
  (`p1`–`p7`) references any renamed drug-specific column — all seven
  plot only disease-side outputs (`PSA`, `T`, `TC_p`, `TC_q`, `CRPC_frac`,
  `ARv7_frac`, `BMD`, `BoneMets`, `AR_nuc`), unchanged names — so no plot
  code needed updating beyond the model-name string passed to `mcode()`.
- Compartment order and 1-based indices are unchanged from the original
  for all 33 compartments (nothing was added, removed, or reordered), so
  no external code addressing this model by compartment number is
  affected.
- Every other disease-side computation (HPG axis, AR signaling, tumor
  kinetics, PI3K/AKT, bone metastasis) is untouched apart from reading the
  renamed drug-effect variables it already read before (e.g. `Enz_AR_inh`
  → `EFFECT_ENZ` inside `AR_bind`/`AR_nuc_eff`/`CRPC_growth`/`ARv7_growth`).

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`: the
single, mislabeled `prostate-cancer | AR Signaling (DEG)` row was replaced
with eight rows, one per real compound identified above (Leuprolide,
Degarelix, Relugolix, Enzalutamide, Abiraterone, Docetaxel, Olaparib,
Denosumab), each classified "Redirect concentration (clean single site)"
— every one of the eight exposes exactly one concentration variable
feeding exactly one Hill-shaped effect term, the cleanest of the census's
own categories.
