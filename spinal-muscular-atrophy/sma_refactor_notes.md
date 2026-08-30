# Refactor notes — `spinal-muscular-atrophy/sma_mrgsolve_model.R`

## The census entry was mislabeled — corrected, not skipped

The task instructions flagged the file's only
`driver-patches/data/compound_perturbation_census.md` row, **"SMN2
Splicing (RIS)"**, for a sanity check: is this a process-description
phrase mislabeled as a drug (the same mistake already found and
corrected once in `neonatal-hyperbilirubinemia/nhb_refactor_notes.md`),
or a real compound? Reading the actual code confirms it is the same
pattern again:

- "SMN2 Splicing" is the *process* the drugs act on (the disease
  mechanism), not a compound name — the census classifier appears to
  have grabbed the nearest doc-comment phrase again.
- The parenthetical **"(RIS)"** is the real PK stem in the file:
  `Emax_RIS`/`EC50_RIS`/`hill_RIS`, `ka_RIS`/`F_RIS`/`Vd_RIS`/`CL_RIS`/
  `Kp_brain`/`fu_RIS`, and three dosed compartments (`A_gut_RIS`,
  `A_plasma_RIS`, `A_CNS_RIS`) with their own `ev_risdiplam()`/
  `ev_risdiplam_peds()` dosing builders (oral, 5 mg/day adult or
  0.2 mg/kg/day pediatric) driving Scenarios 3 and 6. This is
  unambiguously **risdiplam**, the oral SMN2-splicing-modifier small
  molecule (Evrysdi), a real, exogenously administered drug — exactly
  as the task's hypothesis anticipated.
- **The file models two further real, independently-dosed compounds the
  census never listed at all**: **nusinersen** (`Emax_NUS`/`EC50_NUS`/
  `hill_NUS`, intrathecal CSF PK `A_CSF_L`/`A_CSF_C`/`A_CNS_NUS`,
  `ev_nusinersen()`, Scenarios 2 and 5 — ENDEAR/CHERISH dosing) and
  **onasemnogene abeparvovec / Zolgensma** (`k_vg_clear`/`k_MN_trans`/
  `k_tg_txn`/`k_tg_mRNA_deg`/`k_tg_prot`/`eff_tg`/`AAV9_Ab_block`, IV
  AAV9 vector PK `A_plasma_ZOL`/`A_MN_ZOL`/`A_tg_mRNA`, `ev_zolgensma()`,
  Scenario 4 — a single 1.1×10¹⁴ vg/kg-style IV dose). Both are named
  explicitly in the file's own header comment and calibration references
  (Darras 2019 NEJM ENDEAR, Mercuri 2018 NEJM CHERISH, Day 2021 NEJM
  STR1VE) — this is not a stretch reading.

**Conclusion: all three — nusinersen, risdiplam, and onasemnogene
abeparvovec (Zolgensma) — are real, externally-dosed drugs with their
own PK, and all three were refactored.** The census is corrected below:
the existing row is renamed to risdiplam and two new rows are added for
nusinersen and Zolgensma, all recorded against this same file.

## Archetype per compound

**Risdiplam (RIS) — Archetype 3 minus peripheral (depot + central,
linear) for the plasma PK that PD actually reads**, same shape as
Anakinra in `familial-mediterranean-fever/fmf_refactor_notes.md` and the
three compounds in `neonatal-hyperbilirubinemia/nhb_refactor_notes.md`.
`GUT_RIS` (depot, was `A_gut_RIS`) → `CENT_RIS` (was `A_plasma_RIS`),
first-order absorption (`KA_RIS`, was `ka_RIS`, converted /h→/day),
first-order elimination (`CL_RIS`/`V1_RIS`, was `CL_RIS`/`Vd_RIS`), no
peripheral compartment. `EFFECT_RIS`'s Hill formula reads `C_RIS`,
defined as `(CENT_RIS/V1_RIS)*1000.0` (plasma, ng/mL) — a pure rename of
the original's `C_plasma_ng`.

One bespoke, disclosed deviation: the original also carries a third
compartment, `A_CNS_RIS` (renamed `CNSTISSUE_RIS`), fed directly from
the gut depot scaled by `Kp_brain` (not from plasma) and eliminated at
`CL_day*Kp_brain*(A_CNS_RIS/Vd_RIS)` — a non-standard construction, not a
real two-way peripheral exchange. **It is preserved unchanged and
renamed only, not removed**, per "don't flatten a mechanistically rich
model" — but it is genuinely PD-silent in the original: the Hill effect
term reads `C_plasma_ng` only, never `C_CNS_ng` (which the original computes
in `$ODE` and then never uses again either). This is disclosed rather
than quietly dropped; keeping it is the more conservative choice even
though it changes no output.

**Nusinersen (NUS) — bespoke, does not fit any of the guide's four
archetypes.** The original's own structure is: `A_CSF_L` (dosing site,
renamed `CENT_NUS`) --one-way bulk CSF flow (`Q_CSF`, renamed
`Q_NUS`)--> `A_CSF_C` (renamed `PERI_NUS`), plus a **separate, saturable
Michaelis-Menten uptake** (`Vmax_NUS`/`Km_NUS`, renamed `VMAX_NUS`/
`KM_NUS`) from `CENT_NUS` into `A_CNS_NUS` (renamed `TISSUE_NUS`, the
compartment the Hill effect actually reads). This is neither a linear
1-3-compartment archetype (the CSF flow is one-way, not a symmetric
`Q`-exchange, and the CNS uptake is nonlinear) nor Archetype 4/TMDD
(there is no receptor pool, on/off-rate, or complex-formation dynamic —
just a saturable *transport* flux). Per the guide's explicit allowance
("if a compound's PK genuinely doesn't resemble any archetype... rename
to the convention, isolate it from PD, expose `C_<STEM>`, and note...
why"), this is handled as bespoke: renamed to the `_NUS` convention,
isolated in its own `$ODE` block, with `C_NUS = TISSUE_NUS/V_CNS_NUS`
as the single PD-reading concentration. One parameter, `ka_CNS_NUS`
(renamed `KA_CNS_NUS`), is declared in the original but never referenced
in any `dxdt_`/`$TABLE` expression — a pre-existing unused/vestigial
parameter, preserved and flagged rather than silently dropped.

**Onasemnogene abeparvovec / Zolgensma (ZOL) — bespoke, both in PK shape
and in how it reaches the disease system.** PK: `A_plasma_ZOL` (renamed
`CENT_ZOL`) has pure linear elimination, `dxdt = -k_vg_clear*A_plasma_ZOL`
— structurally identical to Archetype 1's single-compartment shape —
but the original's own compartment comment ("vg/mL * Vdist") makes clear
the state is already concentration-valued, with no separate volume term
to divide by. This is the same "state IS the concentration" deviation
disclosed for Stannsoporfin/Phenobarbital/UDCA in
`neonatal-hyperbilirubinemia/nhb_refactor_notes.md` — `C_ZOL = CENT_ZOL`
is an identity, not a division, preserved rather than forced into a
textbook amount/volume split that would change the original's dynamics.
Downstream, `CENT_ZOL` drives an **irreversible transduction →
transcription cascade** (`A_MN_ZOL`→`A_tg_mRNA`, renamed
`MNLOAD_ZOL`→`MRNA_ZOL`) that feeds `SMN_synthesis` directly — this is
not a concentration-response relationship of any kind (no equilibrium,
no saturation curve, nothing an `nls()` Hill fit would be fitting
against), so per the guide's Hill-interface section ("if you can't write
it down as one line... fit an approximation" only applies when there
*is* a concentration-response relationship to approximate) no fit was
attempted. Instead, `EFFECT_ZOL` is declared as an **identity of
`MRNA_ZOL`**, the exact quantity the original fed straight into
`SMN_synthesis` — this names the point where the disease reads ZOL's
contribution without inventing pharmacology the original doesn't have.
Disclosed here as bespoke on both counts (PK shape and effect
interface), per "if a compound's structure genuinely doesn't fit any of
the guide's four archetypes, say so and handle it as bespoke."

## Renaming applied (parameter values are all copied verbatim from the original)

| Original | Refactored | Role |
|---|---|---|
| `A_CSF_L` (cmt 1) | `CENT_NUS` | NUS dosing/"central" CSF site |
| `A_CSF_C` (cmt 2) | `PERI_NUS` | NUS "peripheral" CSF site |
| `A_CNS_NUS` (cmt 3) | `TISSUE_NUS` | NUS CNS tissue — **the exposed PD site** |
| `V_CSF_L` | `V1_NUS` | NUS central volume |
| `V_CSF_C` | `V2_NUS` | NUS peripheral volume |
| `V_CNS` | `V_CNS_NUS` | NUS tissue volume (bespoke, no table slot) |
| `Q_CSF` | `Q_NUS` | NUS one-way CSF bulk flow |
| `CL_ASO_CSF` | `CL_NUS` | NUS clearance from CSF pools |
| `CL_ASO_CNS` | `CL_CNS_NUS` | NUS elimination from CNS tissue (bespoke) |
| `ka_CNS_NUS` | `KA_CNS_NUS` | **unused in original — preserved, flagged** |
| `Vmax_NUS` | `VMAX_NUS` | NUS saturable-uptake Vmax |
| `Km_NUS` | `KM_NUS` | NUS saturable-uptake Km |
| `Emax_NUS` | `EMAX_NUS` | NUS Hill Emax |
| `EC50_NUS` | `EC50_NUS` | NUS Hill EC50 (name unchanged) |
| `hill_NUS` | `GAMMA_NUS` | NUS Hill coefficient |
| `C_CNS_N` (local) | `C_NUS` | **the exposed concentration** |
| `Emax_eff_NUS` (local) | `EFFECT_NUS` | NUS splicing effect |
| `A_gut_RIS` (cmt 4) | `GUT_RIS` | RIS depot |
| `A_plasma_RIS` (cmt 5) | `CENT_RIS` | RIS plasma — **the exposed PD site** |
| `A_CNS_RIS` (cmt 6) | `CNSTISSUE_RIS` | RIS CNS partition — PD-silent, preserved |
| `ka_RIS` | `KA_RIS` | RIS absorption rate |
| `F_RIS` | `F_RIS` | RIS bioavailability (name unchanged) |
| `Vd_RIS` | `V1_RIS` | RIS volume of distribution |
| `CL_RIS` | `CL_RIS` | RIS clearance (name unchanged) |
| `Kp_brain` | `KP_RIS` | RIS brain:plasma partition (bespoke) |
| `fu_RIS` | `FU_RIS` | **unused in original — preserved, flagged** |
| `Emax_RIS` | `EMAX_RIS` | RIS Hill Emax |
| `EC50_RIS` | `EC50_RIS` | RIS Hill EC50 (name unchanged) |
| `hill_RIS` | `GAMMA_RIS` | RIS Hill coefficient |
| `C_plasma_ng` (local) | `C_RIS` | **the exposed concentration** |
| `Emax_eff_RIS` (local) | `EFFECT_RIS` | RIS splicing effect |
| `A_plasma_ZOL` (cmt 7) | `CENT_ZOL` | ZOL plasma, concentration-valued (bespoke) |
| `A_MN_ZOL` (cmt 8) | `MNLOAD_ZOL` | ZOL transduced MN vg load (bespoke) |
| `A_tg_mRNA` (cmt 9) | `MRNA_ZOL` | ZOL transgene mRNA — **the PD-reading quantity** |
| `k_vg_clear` | `CL_ZOL` | ZOL plasma clearance |
| `k_MN_trans` | `K_TRANS_ZOL` | ZOL transduction rate (bespoke) |
| `k_tg_txn` | `K_TXN_ZOL` | ZOL transgene transcription rate |
| `k_tg_mRNA_deg` | `K_MRNA_DEG_ZOL` | ZOL transgene mRNA degradation |
| `k_tg_prot` | `K_PROT_ZOL` | **feeds only a dead-code toggle in the original — preserved, flagged** |
| `eff_tg` | `EFF_TG_ZOL` | ZOL transgene expression efficiency |
| `AAV9_Ab_block` | `AAV9_AB_BLOCK_ZOL` | ZOL anti-vector antibody block |
| `A_tg_mRNA` (direct read in `SMN_synthesis`) | `EFFECT_ZOL` (identity of `MRNA_ZOL`) | **the disease-facing ZOL interface — not a Hill function, see above** |

Compartment order (and therefore the 1-based `cmt` indices `ev_nusinersen()`/
`ev_risdiplam()`/`ev_zolgensma()` dose into: 1, 4, 7) is unchanged from the
original — only names changed.

## Hill interface: rename, not a fit — for nusinersen and risdiplam

Both original effect terms were already exactly the Hill/Emax shape, so
no `nls()` fit was performed for either:

- **Nusinersen**: `Emax_eff_NUS = Emax_NUS*pow(C_CNS_N,hill_NUS)/(pow(EC50_NUS,hill_NUS)+pow(C_CNS_N,hill_NUS))`
  is exactly `EFFECT_NUS`, a one-to-one rename.
- **Risdiplam**: same shape, `EFFECT_RIS` is a one-to-one rename of
  `Emax_eff_RIS`.
- The combination rule in `$ODE` — `E7I_base + (E7I_max-E7I_base)*
  (EFFECT_NUS + EFFECT_RIS - EFFECT_NUS*EFFECT_RIS)`, then clamped and
  scaled by `SMN2_copies/2` — is unchanged from the original, and keeps
  each compound's `EFFECT_<STEM>` separate until the single point disease
  equations combine them, per the guide's "never collapse several drugs
  into one shared Hill term" rule.
- **Zolgensma** has no Hill/occupancy relationship in the original (see
  archetype discussion above) — `EFFECT_ZOL` is an identity of the
  transgene mRNA pool, not a fitted curve. No R² is reported because
  there is no curve being approximated; this is disclosed as the correct
  bespoke handling per the guide, not a skipped fit.

## A build defect in the original, fixed syntax-only in the delivered file

The untouched original does not compile under mrgsolve 2.0.1:
`$TABLE` declares `double CMAP = CMAP_max*MN_pool*NMJ_score;` and then
writes `capture CMAP = CMAP;` — a self-referential capture that makes
mrgsolve re-declare `CMAP` a second time, colliding with the `double`
already in scope. Confirmed via `POST /model_manifest` on the untouched
original alone:

```
65:11: error: redefinition of 'capture {anonymous}::CMAP'
56:10: note: 'double {anonymous}::CMAP' previously declared here
```

...and identically for `HFMSE` and `RULM` (both follow the same
`double X = ...; capture X = X;` pattern). `FVC`/`CHOP_INTEND` do **not**
hit this, because their own local variable names (`FVC_pct`, `CHOP`)
already differ from their capture names — confirming the mechanism is
purely the self-referential name match, not something structural about
these three variables. Logged as `UPSTREAM_ISSUES.md` #66, together with
two additional, non-build-blocking findings (below).

**Fix applied directly to the delivered `sma_mrgsolve_model_refactored.R`**,
per the guide's settled policy for a non-compiling original: the three
colliding locals were renamed `CMAP_val`/`HFMSE_val`/`RULM_val` (the same
pattern the original already used successfully for `FVC_pct`/`CHOP`),
and their `capture` lines updated to read from the renamed local instead
of re-declaring the capture name. This changes nothing numeric — see
verification below. The checked-in original (`sma_mrgsolve_model.R`) is
untouched and still carries the defect exactly as written.

A second, related design point (not an upstream defect, a decision made
during this refactor): `C_NUS`/`C_RIS`/`C_ZOL`/`EFFECT_NUS`/`EFFECT_RIS`/
`EFFECT_ZOL` are declared once as `$GLOBAL` doubles and assigned (no
`double` keyword) in `$ODE`, then listed **bare** in a `$CAPTURE` block
(`C_NUS C_RIS C_ZOL EFFECT_NUS EFFECT_RIS EFFECT_ZOL`, no `= expr`) —
never as `capture NAME = NAME;` — for exactly the same reason the three
`$TABLE` locals above needed fixing: writing `capture C_NUS = C_NUS;`
would re-collide, this time between `$GLOBAL` and `$TABLE` scope. This
is the same pattern established in
`neonatal-hyperbilirubinemia/nhb_refactor_notes.md` for its
`C_SNMP`/`C_PB`/`C_UDCA`/`EFFECT_SNMP`/`EFFECT_PB`/`EFFECT_UDCA`.

## A verification-time bug found and fixed in this refactor's own first draft (not an upstream defect)

The first draft reused an `$ODE`-block-local `double conc_CENT_NUS =
CENT_NUS/V1_NUS;` (needed for the CSF-flow/uptake flux math) directly in
`$TABLE` for the informational `C_CSF_lumbar` capture. This compiled
and ran, but comparison against the original showed a **400 ng/mL**
discrepancy (`C_CSF_lumbar`, refactored vs. original) at exactly the
timepoints coinciding with a nusinersen dose — 400 = 12000 ng dose /
30 mL (`V1_NUS`), i.e. the entire size of one bolus. The cause: mrgsolve
can call `$TABLE` for a reported row without having just re-evaluated
`$ODE` at that exact instant, so a block-local `double` computed inside
`$ODE` can still hold the value from the solver's last internal
derivative evaluation rather than the state exactly at this reported
time — invisible for a continuously-varying compartment (which is why
`C_NUS`/`EFFECT_NUS`/`E7_inclusion`, all downstream of the *continuous*
`TISSUE_NUS`, were unaffected), but a full dose-sized artifact exactly
at a row where a bolus lands on the *directly-dosed* compartment being
read. **Fixed** by recomputing every exposed/reported quantity directly
from state *inside* `$TABLE` (`double CENT_NUS_conc_out = CENT_NUS /
V1_NUS;` for the informational output, and reassigning — not
re-declaring — the six `$GLOBAL` exposed variables from state at the top
of `$TABLE`), exactly mirroring how the original always computed every
`$TABLE` quantity fresh from state rather than reusing an `$ODE` local.
This is the same general risk for any compound dosed directly into its
own exposed/"central" compartment (nusinersen and Zolgensma both are,
here) rather than through a depot — worth flagging for any future
refactor of a similar bespoke, depot-less dosing structure.

## Two additional findings, disclosed, not fixed (logged in `UPSTREAM_ISSUES.md` #66)

1. **`E7_inclusion` (the captured output) is not the same quantity that
   drives disease dynamics.** `$ODE`'s `E7I_current` combines
   `EFFECT_NUS`/`EFFECT_RIS` with a product-complement term and scales
   by `SMN2_copies/2`; `$TABLE`'s own `E7I_out` (captured as
   `E7_inclusion`) independently recomputes the same two Hill terms from
   the same underlying concentrations but sums them plainly, with no
   complement term and no `SMN2_copies` scaling. This inconsistency
   pre-exists in the original (confirmed algebraically — the two `C`
   values used in each are identical) and is preserved as-is here
   (renamed variables only, same formula in each place); it does not
   affect verification since both the refactored and original files
   compute the same (inconsistent) `E7I_out` formula identically.
2. **`k_tg_prot`/`SMN_from_ZOL` is dead code.** `double SMN_from_ZOL =
   k_tg_prot > 0 ? A_tg_mRNA : 0.0;` is computed but never read anywhere
   else — `SMN_synthesis` reads `A_tg_mRNA` (renamed `MRNA_ZOL`, i.e.
   `EFFECT_ZOL`) directly. `k_tg_prot`'s value has no effect on any
   model output through this expression. Preserved verbatim (renamed
   `K_PROT_ZOL`), not fixed or removed, per "log what you find, don't
   fix" — this is inert either way.

## Verification

**Method.** Both files' embedded DSL blocks (`code <- '...'` equivalent
— here `mrgsolve::mcode("sma_qsp"/"sma_qsp_refactored", '...')`) were
mechanically extracted (regex on the `mcode(...)` call, verbatim quoted
text) and POSTed to the local qspserver `mrgsolve_api` service at
`http://localhost:8007` (`/model_manifest` then `/run_simulation`),
which compiles and runs each DSL block directly with mrgsolve 2.0.1
server-side — no local R/mrgsolve install used. Requests were spaced
~2 seconds apart and run sequentially (never more than one in flight),
respecting the service's `max_concurrent_jobs: 2` limit and its history
of crashing under concurrent load (one transient `dyn.load`/"file too
short" cache error was hit and cleared on retry, unrelated to either
model's own content). The refactored DSL, once embedded inside R's
single-quoted string literal in the delivered `.R` file, requires
escaping the 11 literal apostrophes in its comments as `\'` for valid R
syntax; confirmed these are purely a source-syntax artifact — applying
R's own single-quote unescaping (`\'` → `'`) to the embedded text
reproduces the exact 19,276-character DSL string that was compiled and
run against qspserver, byte-for-byte identical.

**Fix used for verification.** Because the untouched original does not
compile (see above), the identical syntax-only `CMAP_val`/`HFMSE_val`/
`RULM_val` rename was applied to an **in-memory-only scratch copy** of
the original (never to the checked-in `sma_mrgsolve_model.R`) so both
sides could actually build and run. No `.cpp` extraction file was left
behind — extraction was in-memory only, used to build the verification
requests and then discarded.

**Scenarios run — all six of the file's own `run_scenarios()`, not
invented, full 730-day duration (no shortening needed; none approached
the API's default `maxsteps` budget):**

1. **Scenario 1, "SMA Type I — No Treatment"** (`SMN2_copies=2, MN0=1.0,
   k_MN_death=0.004`, no dosing): exact match, max abs diff 0.0, across
   all 13 shared `$CAPTURE` outputs, 732 rows.
2. **Scenario 2, "SMA Type I — Nusinersen"** (same params; intrathecal
   12,000 ng boluses at t=0,14,28,63,183,303,423,543,663 into cmt 1):
   exact match, max abs diff 0.0, 740 rows (including the fix described
   above for `C_CSF_lumbar` — confirmed 0.0 only after that fix; the
   first draft showed a 400 ng/mL artifact at every dose row before it).
3. **Scenario 3, "SMA Type II — Risdiplam"** (`SMN2_copies=3, MN0=1.0,
   k_MN_death=0.002`; 5 mg/day oral into cmt 4, `ii=1, addl=730`): exact
   match, max abs diff 0.0, 732 rows.
4. **Scenario 4, "Presymptomatic SMA — Zolgensma"** (`SMN2_copies=2,
   MN0=0.95, k_MN_death=0.004`; single 16.5-unit IV bolus into cmt 7):
   exact match, max abs diff 0.0, 732 rows.
5. **Scenario 5, "SMA Type II — Nusinersen Late Start"**
   (`SMN2_copies=3, MN0=1.0, k_MN_death=0.002`; boluses at
   t=365,379,393,428,548,668): exact match, max abs diff 0.0, 737 rows
   (same `C_CSF_lumbar` fix confirmed here too — this scenario is what
   originally exposed the bug, at the dose nearest t≈366).
6. **Scenario 6, "SMA Type II — Risdiplam Pediatric (15 kg)"**
   (`SMN2_copies=3, MN0=1.0, k_MN_death=0.002`; 3 mg/day into cmt 4,
   `ii=1, addl=730`): exact match, max abs diff 0.0, 732 rows.

All six scenarios verify **bit-for-bit exact (max abs diff 0.0)** across
all 13 of the original's own `$CAPTURE` outputs (`CMAP`, `FVC`, `HFMSE`,
`CHOP_INTEND`, `RULM`, `E7_inclusion`, `SMN_protein`, `MN_fraction`,
`NMJ_maturity`, `C_CSF_lumbar`, `C_CNS_nusinersen`, `C_plasma_risdiplam`,
`Transgene_mRNA`), consistent with the guide's tolerance table for
"pure structural reorganization" (Archetypes 1-3 and, here, the bespoke
NUS/ZOL structures too, since none of them required a Hill-fit
approximation).

## qspserver `/model_manifest` discoverability

Confirmed via `POST /model_manifest` against the extracted, renamed DSL:
all 55 renamed `$PARAM` entries (29 disease-side plus 11 NUS + 9 RIS +
7 ZOL PK/Hill parameters) appear in the manifest's `parameters` with
their original numeric defaults, including the two preserved-but-unused
ones (`KA_CNS_NUS`, `FU_RIS`) and the dead-code-feeding one
(`K_PROT_ZOL`). `C_NUS`/`C_RIS`/`C_ZOL`/`EFFECT_NUS`/`EFFECT_RIS`/
`EFFECT_ZOL` are state-derived (computed in `$ODE`, refreshed in
`$TABLE`), so — per the same reasoning established for
`C_SNMP`/`EFFECT_SNMP`/etc. in
`neonatal-hyperbilirubinemia/nhb_refactor_notes.md` — they cannot also be
`$PARAM` entries; all six appear in the manifest's `outputPaths` via the
bare `$CAPTURE` list, confirmed discoverable, alongside the file's own
13 renamed capture names and all 17 `$CMT` compartments (36
`outputPaths` total).

## Anything else flagged

- No compound other than nusinersen, risdiplam, and Zolgensma exists in
  this file — all three were in scope (the census previously named only
  one, mislabeled; see above).
- The R-side dosing builders (`ev_nusinersen()`, `ev_risdiplam()`,
  `ev_risdiplam_peds()`, `ev_zolgensma()`) needed no changes beyond
  updated comments — they address compartments by 1-based `cmt` index,
  unchanged from the original, not by name.
- `virtual_population()`'s `idata` data frame referenced `Emax_NUS` by
  name (renamed `EMAX_NUS` by this refactor's Hill-parameter naming
  convention) — updated to match; `EC50_NUS`/`k_prot_deg`/`k_MN_death`/
  `SMN2_copies` were already spelled identically in the original and
  needed no change.
- `run_scenarios()`, `plot_results()`, `sensitivity_analysis()` reference
  only `$CAPTURE` output names (all unchanged) and disease-side
  parameters (all unchanged) — no other updates needed.

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`: the
existing `spinal-muscular-atrophy` row corrected from "SMN2 Splicing
(RIS)" to "Risdiplam (RIS)", and two new rows added for "Nusinersen
(NUS)" and "Onasemnogene abeparvovec / Zolgensma (ZOL)" — all three
against this same file, per the findings above.
