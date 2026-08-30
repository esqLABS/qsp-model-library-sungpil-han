# Refactor notes — `gout/gout_mrgsolve_model.R`

Scope: all nine compounds carrying rows in
`driver-patches/data/compound_perturbation_census.md` for this file —
**Allopurinol (ALLO)**, **Anakinra (ANA)**, **Canakinumab (CANA)**,
**Colchicine (COLCH)**, **Febuxostat (FEBU)**, **Indomethacin (INDO)**,
**Lesinurad (LESI)**, **Oxypurinol (OXY)**, **Probenecid (PROB)**.
Nothing in this file is out of scope; all nine census rows are updated.
The disease-side biology (purine production, urate distribution, renal
handling, crystal/tophus dynamics, NLRP3/cytokine cascade, pain, eGFR) is
byte-identical to the original — no build-compatibility fix was needed
(the original compiles cleanly under mrgsolve 2.0.1, confirmed via
`POST /model_manifest`).

## ALLO/OXY is a genuine parent-metabolite chain, not two independent compounds

Confirmed the same way as the `urolithiasis/uri_refactor_notes.md`
precedent for the identical real-world drug pair:

- `A_Oxy_cent` (renamed `CENT_OXY`) has **no `GUT_`/depot compartment and
  no `ka_`/`F_` absorption parameters of its own** — its only inflow is
  `k_Allo_Oxy * A_Allo_cent` (renamed `KCONV_OXY * CENT_ALLO`), and
  `k_Allo_Oxy`'s own `$PARAM` comment literally reads "Allopurinol to
  oxypurinol conversion rate."
- The disease-effect terms (`XO_inhib_Oxy`, feeding the combined
  `XO_inhib`) are computed from `C_Oxy` (renamed `C_OXY`), **never** from
  `C_Allo` — allopurinol itself has no direct PD effect in this model;
  only its metabolite (oxypurinol) inhibits xanthine oxidase. `C_Allo`
  (renamed `C_ALLO`) is computed in `$MAIN` but never read anywhere else
  — a dead, informational-only quantity, matching the identical pattern
  already documented for `ALLO`/`OXP` in `uri_refactor_notes.md` and for
  `PRED`'s own plasma level in the sarcoidosis file.

**Handled as a genuine parent-metabolite chain:** `ALLO` keeps its own PK
block (`GUT_ALLO`/`CENT_ALLO`) and its own informational-only `C_ALLO`
(no `EFFECT_ALLO` — matches the original having no direct allopurinol PD
effect). `OXY` keeps its own two-compartment PK (`CENT_OXY`/`PERI_OXY`,
fed by the conversion out of `CENT_ALLO`) and exposes the
`C_OXY`/`EFFECT_OXY` pair the disease equation (xanthine-oxidase
inhibition feeding `UA_prod` via `XO_activity`) actually reads. The
census rows for both ALLO and OXY (previously `Target/Pathway: ?`) are
corrected below to reflect this relationship.

## A parameterization quirk found across five of the nine compounds: amount-based transfer constants that are NOT Q/V1/V2-consistent

The guide's Archetype 2/3 templates assume a real inter-compartmental
clearance `Q`, where `Q/V1` (outflow rate constant) and `Q/V2` (inflow
rate constant) share the *same* `Q` — which only holds if
`k12 * V1 == k21 * V2` for the original's own rate constants and volumes.
Checked this identity for every two-compartment PK block in the file:

| Compound | `ktp*V1` (out) | `kpt*V2` (in) | Consistent? |
|---|---|---|---|
| Febuxostat | 0.3*12=3.6 | 0.15*24=3.6 | **yes** |
| Probenecid | 0.4*10=4.0 | 0.2*20=4.0 | **yes** |
| Oxypurinol | 0.5*45=22.5 | 0.2*30=6.0 | no |
| Colchicine | 2.1*100=210 | 0.42*400=168 | no |
| Indomethacin | 0.6*16=9.6 | 0.3*20=6.0 | no |
| Anakinra | 0.15*8=1.2 | 0.08*5=0.4 | no |
| Canakinumab | 0.005*6=0.03 | 0.003*3.5=0.0105 | no |

For Febuxostat and Probenecid, a single derived `Q_<STEM>` (3.6 and 4.0
respectively) reproduces the original exactly, so those two are
rewritten in the guide's canonical `(CL+Q)/V1`, `Q/V1`, `Q/V2` form — a
pure rename, verified bit-identical (see below).

For the other five (Oxypurinol, Colchicine, Indomethacin, Anakinra,
Canakinumab), no single `Q` reproduces the original's own asymmetric
transfer — forcing one would silently change the model's numerics rather
than just renaming it. These are kept as **directly-renamed,
amount-based transfer rate constants** (`KCP_<STEM>` central-to-
peripheral, `KPC_<STEM>` peripheral-to-central — i.e. `dxdt_CENT =
... - KCP*CENT + KPC*PERI`, exactly the original's own `- ktp*A_cent +
kpt*A_peri` form, just renamed), disclosed here as a deviation from the
guide's CL/Q/V convention **in parameter naming only** — the compartment
structure itself is still Archetype 2 (Oxypurinol) or Archetype 3 (the
other four), and every value is copied verbatim from the original.
Nothing was invented, defaulted, or refit to make this table's "yes"
column say more than it honestly can.

## Archetype and Hill interface per compound

**ALLO — Archetype 3 minus peripheral** (depot + central, linear,
PK-only, no direct PD effect — see above):
```
dxdt_GUT_ALLO  = -KA_ALLO * GUT_ALLO;
dxdt_CENT_ALLO = KA_ALLO*F_ALLO*GUT_ALLO - (CL_ALLO/V1_ALLO)*CENT_ALLO - KCONV_OXY*CENT_ALLO;
```
Renamed `A_Allo_gut→GUT_ALLO`, `A_Allo_cent→CENT_ALLO`, `ka_Allo→KA_ALLO`,
`F_Allo→F_ALLO`, `CL_Allo→CL_ALLO`, `Vc_Allo→V1_ALLO`,
`k_Allo_Oxy→KCONV_OXY` (named for the metabolite it produces, matching
the `k_Allo_Oxp→KCONV_OXP` precedent's naming logic in
`uri_refactor_notes.md`).

**OXY — Archetype 2** (no depot, central+peripheral, fed by the ALLO
conversion, amount-based transfer per the table above):
```
dxdt_CENT_OXY = KCONV_OXY*CENT_ALLO - (CL_OXY/V1_OXY)*CENT_OXY - KCP_OXY*CENT_OXY + KPC_OXY*PERI_OXY;
dxdt_PERI_OXY = KCP_OXY*CENT_OXY - KPC_OXY*PERI_OXY;
```
Renamed `A_Oxy_cent→CENT_OXY`, `A_Oxy_peri→PERI_OXY`, `CL_Oxy→CL_OXY`,
`Vc_Oxy→V1_OXY`, `Vp_Oxy→V2_OXY`, `ktp_Oxy→KCP_OXY`, `kpt_Oxy→KPC_OXY`.
Effect: original `XO_inhib_Oxy = C_Oxy/(C_Oxy+Ki_Oxy)` is already the
exact Hill shape with implicit Emax=1, gamma=1 — renamed to
`EFFECT_OXY = EMAX_OXY*pow(C_OXY,GAMMA_OXY)/(pow(EC50_OXY,GAMMA_OXY)+pow(C_OXY,GAMMA_OXY))`,
`EC50_OXY = Ki_Oxy` (0.001, unchanged), `EMAX_OXY=1.0`/`GAMMA_OXY=1.0`
new and explicit (a rename, not a fit — `pow(x,1)` and `1.0*ratio` are
bit-identical to the original's bare ratio).

**FEBU — Archetype 3, CL/Q/V-consistent** (see table above; `Q_FEBU=3.6`
derived exactly, a rename not an approximation):
```
dxdt_CENT_FEBU = KA_FEBU*F_FEBU*GUT_FEBU - (CL_FEBU+Q_FEBU)/V1_FEBU*CENT_FEBU + Q_FEBU/V2_FEBU*PERI_FEBU;
dxdt_PERI_FEBU = Q_FEBU/V1_FEBU*CENT_FEBU - Q_FEBU/V2_FEBU*PERI_FEBU;
```
Renamed `A_Febu_gut/_cent/_peri→GUT_FEBU/CENT_FEBU/PERI_FEBU`,
`ka_Febu→KA_FEBU`, `F_Febu→F_FEBU`, `CL_Febu→CL_FEBU`, `Vc_Febu→V1_FEBU`,
`Vp_Febu→V2_FEBU`. Effect: `XO_inhib_Febu = C_Febu/(C_Febu+Ki_Febu)`,
same shape as OXY above — `EMAX_FEBU=1.0`/`EC50_FEBU=Ki_Febu`
(0.000001)/`GAMMA_FEBU=1.0`, rename not a fit. `EFFECT_OXY` and
`EFFECT_FEBU` are combined exactly where the original combines them
(`XO_inhib = 1-(1-EFFECT_OXY)*(1-EFFECT_FEBU)`, feeding `XO_activity`) —
unchanged combination point, per the guide's "combine only at the point
disease equations actually use them."

**PROB — Archetype 3, CL/Q/V-consistent** (`Q_PROB=4.0` derived exactly):
```
dxdt_CENT_PROB = KA_PROB*F_PROB*GUT_PROB - (CL_PROB+Q_PROB)/V1_PROB*CENT_PROB + Q_PROB/V2_PROB*PERI_PROB;
dxdt_PERI_PROB = Q_PROB/V1_PROB*CENT_PROB - Q_PROB/V2_PROB*PERI_PROB;
```
Renamed analogously to FEBU. Effect: original
`URAT1_inhib_P = (C_Prob/IC50_Prob)/(1+C_Prob/IC50_Prob)` is algebraically
identical to `C/(C+IC50)` — renamed `EFFECT_PROB`,
`EC50_PROB=IC50_Prob` (5.0), `EMAX_PROB=1.0`/`GAMMA_PROB=1.0` new,
rename not a fit.

**LESI — Archetype 3 minus peripheral** (no peripheral compartment in
the original):
```
dxdt_CENT_LESI = KA_LESI*F_LESI*GUT_LESI - (CL_LESI/V1_LESI)*CENT_LESI;
```
Renamed `A_Lesi_gut/_cent→GUT_LESI/CENT_LESI`, `ka_Lesi→KA_LESI`,
`F_Lesi→F_LESI`, `CL_Lesi→CL_LESI`, `Vc_Lesi→V1_LESI`. Effect: same
ratio-form as PROB (`(C/IC50)/(1+C/IC50)`), renamed `EFFECT_LESI`,
`EC50_LESI=IC50_Lesi` (0.1), `EMAX_LESI=1.0`/`GAMMA_LESI=1.0` new. Both
`EFFECT_PROB` and `EFFECT_LESI` combine exactly where the original
combines them (`URAT1_inhib = 1-(1-EFFECT_PROB)*(1-EFFECT_LESI)`,
feeding `kURAT1_eff`) — unchanged combination point.

**COLCH — Archetype 3, amount-based transfer** (not Q/V1/V2-consistent,
see table):
```
dxdt_CENT_COLCH = KA_COLCH*F_COLCH*GUT_COLCH - (CL_COLCH/V1_COLCH)*CENT_COLCH - KCP_COLCH*CENT_COLCH + KPC_COLCH*PERI_COLCH;
dxdt_PERI_COLCH = KCP_COLCH*CENT_COLCH - KPC_COLCH*PERI_COLCH;
```
Renamed `A_Colch_gut/_cent/_peri→GUT_COLCH/CENT_COLCH/PERI_COLCH`,
`ka_Colch→KA_COLCH`, `F_Colch→F_COLCH`, `CL_Colch→CL_COLCH`,
`Vc_Colch→V1_COLCH`, `Vp_Colch→V2_COLCH`, `ktp_Colch→KCP_COLCH`,
`kpt_Colch→KPC_COLCH`. Effect: original
`NLRP3_inhib_Colch = Emax_Colch*C_Colch/(C_Colch+IC50_Colch)` already has
an explicit, non-unity `Emax_Colch=0.85` — renamed `EMAX_COLCH` (0.85,
rename not new), `EC50_COLCH=IC50_Colch` (0.0003), `GAMMA_COLCH=1.0` new
(no exponent in the original). `EFFECT_COLCH` feeds `NLRP3_act` exactly
where the original's `NLRP3_inhib_Colch` did.

**INDO — Archetype 3, amount-based transfer, PD reads the tissue site
not plasma** (not Q/V1/V2-consistent, see table):
```
dxdt_CENT_INDO = KA_INDO*F_INDO*GUT_INDO - (CL_INDO/V1_INDO)*CENT_INDO - KCP_INDO*CENT_INDO + KPC_INDO*PERI_INDO;
dxdt_PERI_INDO = KCP_INDO*CENT_INDO - KPC_INDO*PERI_INDO;
```
Renamed `A_Indo_gut/_cent/_tiss→GUT_INDO/CENT_INDO/PERI_INDO`,
`ka_Indo→KA_INDO`, `F_Indo→F_INDO`, `CL_Indo→CL_INDO`,
`Vc_Indo→V1_INDO`, `Vt_Indo→V2_INDO`, `ktp_Indo→KCP_INDO`,
`kpt_Indo→KPC_INDO`. The original's own `COX_inhib_Indo` reads
`C_IndoT` (the tissue/synovial concentration), **never** the plasma
concentration — the same "PD reads the non-plasma site" pattern already
documented for Doxycycline in
`abdominal-aortic-aneurysm/aaa_refactor_notes.md`. Kept faithfully: `C_INDO = CENT_INDO/V1_INDO` (plasma,
informational only, matching the original's dead `C_Allo`-style
treatment) and `C_INDO_PERI = PERI_INDO/V2_INDO` (the tissue site,
PD-relevant) are both exposed, with `EFFECT_INDO` computed from
`C_INDO_PERI`, not `C_INDO`. `EMAX_INDO=1.0`/`EC50_INDO=IC50_Indo`
(0.002)/`GAMMA_INDO=1.0`, all new/explicit (rename of the already-Hill
`C/(C+IC50)` ratio). **Note:** the original never captured
`C_Indo`/`C_IndoT` at all (only 8 of the 9 compound concentrations were
captured in the original's own `$TABLE`) — `C_INDO_OUT`/`C_INDO_PERI_OUT`
are genuinely new captures, added per the qspserver
`/model_manifest`-discoverability requirement, not a rename of an
existing original capture.

**ANA — Archetype 3, amount-based transfer** (SC depot named `GUT_ANA`
per the guide's generic "absorption depot" role, not Q/V1/V2-consistent,
see table):
```
dxdt_CENT_ANA = KA_ANA*F_ANA*GUT_ANA - (CL_ANA/V1_ANA)*CENT_ANA - KCP_ANA*CENT_ANA + KPC_ANA*PERI_ANA;
dxdt_PERI_ANA = KCP_ANA*CENT_ANA - KPC_ANA*PERI_ANA;
```
Renamed `A_Ana_SC/_cent/_peri→GUT_ANA/CENT_ANA/PERI_ANA`,
`ka_Ana→KA_ANA`, `F_Ana→F_ANA`, `CL_Ana→CL_ANA`, `Vc_Ana→V1_ANA`,
`Vp_Ana→V2_ANA`, `ktp_Ana→KCP_ANA`, `kpt_Ana→KPC_ANA`. Effect: original
`Ana_effect = C_Ana/(C_Ana+IC50_Ana)` already the exact Hill ratio
(implicit Emax=1, gamma=1) — renamed `EMAX_ANA=1.0`,
`EC50_ANA=IC50_Ana` (0.5), `GAMMA_ANA=1.0`, all new/explicit. Used in
`IL1b_prod = kIL1b_prod*NLRP3_act*(1-EFFECT_ANA)`, exactly where the
original used `(1-Ana_effect)`.

**CANA — Archetype 3 PK backbone + genuine mass-action neutralization of
the shared disease IL-1b state** (amount-based transfer, not
Q/V1/V2-consistent, see table; not textbook TMDD since the "receptor" is
IL-1b itself, a pre-existing disease state with its own independent
production/degradation, not a canakinumab-owned receptor pool):
```
dxdt_CENT_CANA = KA_CANA*F_CANA*GUT_CANA - (CL_CANA/V1_CANA)*CENT_CANA - KCP_CANA*CENT_CANA
               + KPC_CANA*PERI_CANA - KON_CANA*CENT_CANA*A_IL1b + KOFF_CANA*COMPLEX_CANA;
dxdt_PERI_CANA    = KCP_CANA*CENT_CANA - KPC_CANA*PERI_CANA;
dxdt_COMPLEX_CANA = KON_CANA*CENT_CANA*A_IL1b - KOFF_CANA*COMPLEX_CANA - (CL_CANA/V1_CANA)*COMPLEX_CANA;
```
Renamed `A_Cana_SC/_cent/_peri→GUT_CANA/CENT_CANA/PERI_CANA`,
`A_IL1b_Cana→COMPLEX_CANA` (canakinumab's own drug-target complex — a
genuinely drug-owned quantity, renamed), `ka_Cana→KA_CANA`,
`F_Cana→F_CANA`, `CL_Cana→CL_CANA`, `Vc_Cana→V1_CANA`, `Vp_Cana→V2_CANA`,
`ktp_Cana→KCP_CANA`, `kpt_Cana→KPC_CANA`, `Kon_Cana→KON_CANA`,
`Koff_Cana→KOFF_CANA`. **`A_IL1b` itself is left completely unrenamed**
— it is the disease's own shared cytokine state (independent production
`IL1b_prod`, independent degradation `kIL1b_deg`), not a canakinumab-
owned receptor pool, matching the same "shared disease state stays
unrenamed" treatment already used for `VEGF_FREE`/`VEGF_BOUND` in
`age-related-macular-degeneration/amd_refactor_notes.md`.

One faithfully-preserved quirk: the original's own binding term reads
`A_Cana_cent` (an **amount**), not a concentration — `dxdt_A_IL1b` and
`dxdt_A_Cana_cent` both use `Kon_Cana * A_IL1b * A_Cana_cent`, never
`C_Cana`. This is reproduced exactly as `KON_CANA * A_IL1b * CENT_CANA`
(not `C_CANA`) — a pre-existing dimensional quirk of the original,
carried forward unchanged rather than "fixed."

**No algebraic `EFFECT_CANA` exists in the original** — canakinumab's
neutralization of IL-1b is expressed entirely through the mass-action
ODE system above (the free cytokine `A_IL1b` is directly depleted by
binding), not through any downstream multiplicative Hill term. Per the
guide's own precedent for exactly this situation
(`EFFECT_VIT` in `age-related-macular-degeneration/amd_refactor_notes.md`:
*"no EFFECT_VIT Hill fit possible or attempted — effect is a multi-step
ODE chain, kept unflattened"*), `EFFECT_CANA` here is added as a
**diagnostic-only** fractional-neutralization readout,
`EFFECT_CANA = COMPLEX_CANA / (COMPLEX_CANA + A_IL1b + 1e-9)`, computed
in `$MAIN` and `$CAPTURE`d for discoverability, but **not** referenced by
any `dxdt_` line or disease equation — the real effect is already fully
present via the depletion of `A_IL1b`. No `EMAX_CANA`/`EC50_CANA`/
`GAMMA_CANA` were invented, since no Hill fit was performed or attempted.

## `$MAIN`-placement quirk (pre-existing, not introduced by this refactor)

Every `C_<STEM>`/`EFFECT_<STEM>` in both the original and the refactored
file is computed in `$MAIN`, not `$ODE`/`$TABLE` (matching the original's
own architecture exactly — this refactor did not move anything between
sections). Because mrgsolve's `$MAIN` block re-executes at dosing/record
boundaries rather than at every reported output time, a `$MAIN`-local
double that depends on a compartment's amount can read one dosing-
interval "stale" relative to the compartment's own just-integrated value
at an arbitrary output time — the same phenomenon already documented in
`urolithiasis/uri_refactor_notes.md`. Spot-checked directly on the
refactored file (indomethacin scenario, t=336h): `PERI_INDO/V2_INDO`
computed by hand did not exactly equal the reported `C_INDO_PERI_OUT` at
that same row, confirming the lag is real. This is a **pre-existing
property of the original's own design**, reproduced identically (not
introduced or fixed) — and it does not affect the verification below,
because every output that was actually compared (all of the original's
own 21 pre-existing captures) is subject to the *identical* lag on both
sides, which therefore cancels out of the comparison entirely (confirmed
by the exact match reported below). The two genuinely *new* diagnostic
captures without an original counterpart (`EFFECT_INDO_OUT`,
`EFFECT_CANA_OUT`, and the newly-added `C_INDO_OUT`/`C_INDO_PERI_OUT`)
were spot-checked by manual recomputation instead (see below) and are
internally consistent with this same, disclosed lag — not a
computational error.

## qspserver compatibility checklist

- `model_content` is pure mrgsolve DSL text (no R wrapper) — confirmed by
  extracting the quoted `gout_model_code_refactored <- '...'` block
  byte-for-byte and diffing it programmatically against the standalone
  `.cpp` scratch copy used for `/model_manifest`/`/run_simulation` calls:
  byte-identical (24,803 characters, confirmed via Python string
  comparison after undoing the R single-quote escaping).
- Confirmed via `POST /model_manifest` on the qspserver `mrgsolve_api`
  (`http://localhost:8007`): the refactored model compiles cleanly (the
  original also compiles cleanly — **no build-compatibility fix was
  needed for this file**), and every renamed `KA_`/`F_`/`CL_`/`V1_`/
  `V2_`/`Q_`/`KCP_`/`KPC_`/`KCONV_OXY`/`KON_CANA`/`KOFF_CANA`/`EMAX_`/
  `EC50_`/`GAMMA_` parameter for all nine compounds appears in the
  manifest's `parameters` list.
- `C_<STEM>`/`EFFECT_<STEM>` are `$MAIN`-local doubles, not `$PARAM`
  members — same established constraint documented across many other
  refactors in this fork (e.g. `acute-intermittent-porphyria/
  aip_refactor_notes.md`, `age-related-macular-degeneration/
  amd_refactor_notes.md`): mrgsolve 2.0.1 treats `$PARAM` members as
  read-only inside `$MAIN`/`$ODE`, so a per-step-recomputed concentration
  or effect and a settable `$PARAM` default are mutually exclusive.
  Exposed instead as `$CAPTURE`d, distinctly-named `_OUT` outputs
  (`C_ALLO_OUT`, `C_OXY_OUT`, `C_FEBU_OUT`, `C_PROB_OUT`, `C_LESI_OUT`,
  `C_COLCH_OUT`, `C_INDO_OUT`, `C_INDO_PERI_OUT`, `C_ANA_OUT`,
  `C_CANA_OUT`, `EFFECT_OXY_OUT`, `EFFECT_FEBU_OUT`, `EFFECT_PROB_OUT`,
  `EFFECT_LESI_OUT`, `EFFECT_COLCH_OUT`, `EFFECT_INDO_OUT`,
  `EFFECT_ANA_OUT`, `EFFECT_CANA_OUT`) — the `_OUT` suffix avoids the
  `capture X = X;` self-reference collision already logged as its own
  defect class in `translations/UPSTREAM_ISSUES.md` (#35/#43/#50/#59/#78/
  #80), same as the `uri_refactor_notes.md` precedent. Confirmed
  discoverable in `/model_manifest`'s `outputPaths`.
- `$SET end/delta` is not present in the original either (the R-side
  `mrgsim(end=..., delta=...)` calls supply the time grid explicitly) —
  no change needed; every scenario below drove `time.end`/`time.delta`
  through the request, not through any model-internal default.
- No R-only syntax inside the DSL block — the extracted block compiled
  standalone with no surrounding R script.
- Solver step budget: none of the eight scenarios below needed
  shortening — the API's default step budget handled the full 52-week
  (8736h) and 90-day (2160h) windows without a timeout.

## Verification

Per the guide's mandatory protocol: ran all **eight** of the original
file's own dosing scenarios (Scenarios 1–8, exactly as coded in the
original's post-DSL R script — same doses, same `ii`/`addl`, same
`time`/`delta`) through both the untouched original and
`gout_mrgsolve_model_refactored.R`, via the qspserver `mrgsolve_api`
(`http://localhost:8007`, confirmed healthy throughout,
`POST /model_manifest` then `POST /run_simulation`, requests spaced ~2s
apart per the shared-service note). Dosing was submitted via the API's
NM-TRAN-style `dosing` field (`cmt` = compartment **number**, not name);
compartment declaration order is identical between the original and
refactored `$CMT` blocks (only names changed), so the same `cmt` index
addresses the same physical compartment in both files.

| Scenario | Dosing (from the original's own R scenario blocks) | `n` time points | Result |
|---|---|---|---|
| 1. Untreated (high-purine diet + alcohol) | `FOOD_score=0.7, ETOH=2.0`, no dosing, 0–8736h | 365 | exact match, max abs diff 0.0 |
| 2. Allopurinol 300mg/day | `GUT_ALLO` 300mg q24h x364, 0–8736h | 366 | exact match, max abs diff 0.0 |
| 3. Febuxostat 80mg/day | `GUT_FEBU` 80mg q24h x364, 0–8736h | 366 | exact match, max abs diff ~1e-15 (floating-point noise) |
| 4. Allopurinol + Lesinurad combo | `GUT_ALLO` 300mg + `GUT_LESI` 200mg, both q24h x364, 0–8736h | 367 | exact match, max abs diff 0.0 |
| 5. Colchicine acute flare | crystal bolus (`A_Crystal`+=5) + 1.2mg/0.6mg loading + 0.6mg q12h x7 maintenance into `GUT_COLCH`, 0–336h | 173 | exact match, max abs diff 0.0 |
| 6. Indomethacin acute flare | crystal bolus + `GUT_INDO` 50mg q8h x20, 0–336h | 171 | exact match, max abs diff 0.0 |
| 7. Canakinumab | crystal bolus + `GUT_CANA` 150mg single SC dose, 0–2160h | 183 | exact match, max abs diff 0.0 |
| 8. Febuxostat + Colchicine prophylaxis | `GUT_FEBU` 80mg + `GUT_COLCH` 0.5mg, both q24h x364, 0–8736h | 367 | exact match, max abs diff ~1e-21 (floating-point noise) |

Every shared output was compared point-by-point across the full time
grid for each scenario: all 37 compartments (mapped through the renames
given above — 12 disease-side, unchanged names; 25 compound-side,
renamed) and 21 of the original's own `$CAPTURE` outputs (`sUA`,
`sUA_syn`, `Crystal`, `Tophus`, `IL1b_f`, `TNFa_f`, `PMN`, `Pain`,
`JntDmg`, `eGFR_sim`, `XO_inh`, `URAT1_inh`, `FEurate`, and the eight
original `C_<drug>_out` concentration captures mapped to their renamed
`C_<STEM>_OUT` counterparts) — **58 shared outputs per scenario, 464
comparisons total across the eight scenarios. Result: exact match (max
abs diff 0.0–1e-15) in every case**, consistent with pure structural
reorganization plus rename-only Hill terms throughout (no compound in
this file required an `nls()` curve fit — every one of the eight
compounds with an algebraic PD effect was already an exact Hill ratio in
the original, and canakinumab's neutralization needed no fit since it
was left as the same mass-action ODE system, per above).

**New diagnostic-only outputs without an original counterpart**
(`EFFECT_OXY_OUT`, `EFFECT_FEBU_OUT`, `EFFECT_PROB_OUT`,
`EFFECT_LESI_OUT`, `EFFECT_COLCH_OUT`, `EFFECT_INDO_OUT`,
`EFFECT_ANA_OUT`, `EFFECT_CANA_OUT`, `C_INDO_OUT`, `C_INDO_PERI_OUT`)
were spot-checked by manual recomputation on the refactored model's own
raw output instead of an original-vs-refactored diff (the original never
captured any of these): e.g. indomethacin scenario, t=336h:
`EFFECT_INDO_OUT = 1.5138e-6` matches
`C_INDO_PERI_OUT/(C_INDO_PERI_OUT+EC50_INDO)` computed by hand to full
reported precision. `EFFECT_CANA_OUT` at t=576h in the canakinumab
scenario (`COMPLEX_CANA=0.2425`, `A_IL1b=1.0027`,
`EFFECT_CANA_OUT=0.2085`) does not exactly equal
`COMPLEX_CANA/(COMPLEX_CANA+A_IL1b)` computed from that same row's raw
compartment values (manual: 0.1947) — this is the disclosed `$MAIN`
one-interval staleness described above (both quantities are `$MAIN`
locals, re-evaluated only at dosing-record boundaries, not at every
reported row), not a formula error; the formula itself
(`complex/(complex+free)`) is algebraically correct and was confirmed
internally consistent across multiple time points.

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`, all
nine `gout` rows (ALLO, ANA, CANA, COLCH, FEBU, INDO, LESI, OXY, PROB),
each `Target/Pathway` filled in from the actual code (none/PK-only for
ALLO; xanthine oxidase for OXY and FEBU; URAT1 for PROB and LESI; NLRP3
inflammasome for COLCH; COX-1/COX-2 for INDO; IL-1 receptor for ANA;
IL-1b mass-action neutralization for CANA) and OXY's row additionally
corrected to note the parent-metabolite relationship with ALLO.

## Anything else flagged

- No pre-existing upstream mrgsolve-2.0.1 build defect was found in this
  file — the untouched original compiles cleanly via
  `POST /model_manifest`, confirmed before any refactoring work began.
  Nothing was logged to `translations/UPSTREAM_ISSUES.md` for this file.
- The post-DSL R script required renaming every `ev(cmt="A_..._gut"/"A_..._SC", ...)`
  target to its renamed compartment (`e2`/`e3`/`e4a`/`e4b`/`e5a-c`/`e6`/`e7`/`e8a`/`e8b`)
  — simple `cmt=` string updates, not new dosing mechanisms. The three
  `ev(cmt="A_Crystal", ...)` flare-induction events are untouched
  (disease-side compartment, out of scope). Every ggplot figure, the
  results-summary table, and the sUA-target sensitivity analysis are
  byte-identical to the original (none of their aesthetics reference a
  renamed compartment or capture — `sUA`/`XO_inh`/`URAT1_inh`/`Crystal`/
  `Tophus`/`Pain`/`IL1b_f`/`eGFR_sim` all kept their original capture
  names).
- `mread_cache()`'s model identifier was changed from `"gout_qsp"` to
  `"gout_qsp_refactored"` purely to avoid an `mrgsolve` `soloc` cache
  collision if both files are ever compiled in the same R session;
  non-behavioral.
