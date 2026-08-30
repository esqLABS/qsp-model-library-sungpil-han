# Refactor notes — `bile-acid-diarrhea/bam_mrgsolve_model.R`

## Scope: the census undercounted this file

Per the fork's PK/PD refactor spec (`FORK_WORKFLOW_GUIDE.md`, Part 2), the
census (`driver-patches/data/compound_perturbation_census.md`) carried five
rows for this file, and three of the five labels were process-description
phrases rather than drug names: `Drugs (LOP)`, `Drugs (OND)`, `Drugs (RIF)`,
`FXR / FGF19 (OCA)`, `FXR / FGF19 (TRO)`. Reading the actual code (per the
guide's "README/census is context, not ground truth" rule) confirms:

- **LOP = loperamide** (opioid mu-agonist, slows colonic transit)
- **OND = ondansetron** (5-HT3 receptor antagonist)
- **RIF = rifaximin** (non-absorbable antibiotic; census already had the
  right target, "Bacterial RNA polymerase")
- **OCA = obeticholic acid** (FXR agonist)
- **TRO = tropifexor** (non-steroidal FXR agonist)

The file also models **two further compounds with no census row at all**:
**elobixibat (ELO)**, a pharmacological ASBT inhibitor, and **colesevelam
(SEQ)**, a bile-acid sequestrant — the same undercounting pattern flagged in
`prostate-cancer` (8 real drugs vs. 1 listed). All **seven** compounds were
refactored in this pass; the census corrections record all seven.

## Archetype per compound

| Compound | Archetype | Notes |
|---|---|---|
| Rifaximin (RIF) | 1 variant, no volume | Single compartment, no depot, linear elimination. The original never divides `RIFg` by a volume before using it in the effect term (IC50 is in mg, an amount) — kept exactly as the raw amount, same treatment `BEZ_plasma`/`C_BEZ` got in `primary-sclerosing-cholangitis/psc_refactor_notes.md`. |
| Ondansetron (OND) | 1 variant, no volume | Same shape as RIF: single compartment, no depot, raw-amount IC50 (mg), no volume in the original. |
| Loperamide (LOP) | 3 variant, no volume | Depot (`LOPg`) + effect-site compartment (`LOPe`), but again no volume anywhere — `IC50_lop` is in mg. Renamed `GUT_LOP`/`CENT_LOP`, no `V1_LOP` invented. |
| Elobixibat (ELO) | 1 variant, shared volume | Single compartment (no depot — dosed directly into the ileal lumen), whose "concentration" (`C_ELO`) legitimately reuses the disease's own shared luminal volume `V_ile` (`CENT_ELO/(3*V_ile)*1000`), exactly as the original did — not an invented drug-specific volume. |
| Obeticholic acid (OCA) | 3 variant, bespoke absorption | Depot (`GUT_OCA`) + enterocyte tissue site (`ENT_OCA`, the actual PD-facing site) + systemic plasma (`CENT_OCA`, tracked but never read by any PD equation). **Bespoke**: the depot-to-enterocyte rate `KA_OCA` is not a constant — it is computed each step from `phi`/`ASBTx`/`ENTMASS` (`KA_OCA = KA0_OCA*(FRACPASS_OCA + (1-FRACPASS_OCA)*phi*ASBTx*ENTMASS)`), because part of OCA's own gut uptake genuinely rides the same ASBT transporter the disease's own bile-acid lesion (`phi`) acts on. Kept exactly as the original computed it — flattening this into a constant `KA_OCA` would remove real, intentional pharmacology (OCA absorption itself degrades as the ileal lesion worsens). |
| Tropifexor (TRO) | 3 variant minus peripheral | Depot (`GUT_TRO`) + central (`CENT_TRO`), converted to the guide's preferred `CL/V1` convention: `CL_TRO = 3.3` (promoted from `ke_tro * V_trop = 0.06 * 55`, exact), `V1_TRO = 55` (renamed from `V_trop`). The exposed `C_TRO` is the ileal-partitioned concentration (`FPART_TRO*CENT_TRO/V1_TRO`, renamed from `ftro_ile`), not raw plasma — plasma alone is never read by any PD equation, matching the same "only the PD-facing site is exposed" precedent as OCA and as `BILE_UDCA`/`BILE_OCA` in the PSC refactor. |
| Colesevelam (SEQ) | **Bespoke — none of the four archetypes fit** | A three-stage lumenal transit chain (`GUT_SEQ` stomach → `SI_SEQ` small intestine → `COL_SEQ` colon → out) coupled to disease bile-acid pools by **irreversible mass-action binding against a finite, consumable capacity** (`BND_SEQ`, the bound-bile-acid mass), not receptor occupancy. There is no single "concentration" the way a Hill interface expects; the genuine single coupling point is the free binding capacity, `EFFECT_SEQ_FREECAP = CAP_SEQ*(SI_SEQ+COL_SEQ) - BND_SEQ` (renamed from `cap_free`, `seq_cap`), which both binding fluxes (`Jbind_ile`, `Jbind_col`) draw from. Binding happens at two genuinely distinct physical sites (small intestine, colon) with independent flux equations, so **two** exposed masses are used — `C_SEQ_SI` (`SI_SEQ`) and `C_SEQ_COL` (`COL_SEQ`) — matching the guide's explicit "two only when a genuinely different tissue site matters" exception. No `EMAX_SEQ`/`EC50_SEQ`/`GAMMA_SEQ` were invented: forcing a Hill-shaped effect term onto a consumable-capacity mass-action system would misrepresent it. |

## Renaming applied (drug identifiers only; every value copied verbatim)

| Original | Refactored | Value | Role |
|---|---|---|---|
| `seq_cap` | `CAP_SEQ` | 0.22 | binding capacity per mg resin |
| `kbind` | `KBIND_SEQ` | 0.020 | second-order binding rate |
| `kseq12`/`kseq23`/`kseq3out` | `K1_SEQ`/`K2_SEQ`/`K3_SEQ` | 1.30/0.55/0.05 | transit-chain rates (bespoke role, no depot/peri/cent split fits a 3-stage lumenal chain) |
| `SEQ1`/`SEQ2`/`SEQ3` (cmt) | `GUT_SEQ`/`SI_SEQ`/`COL_SEQ` | — | transit-chain compartments |
| `BND` (cmt) | `BND_SEQ` | — | bound bile acid mass |
| (derived) | `C_SEQ_SI`/`C_SEQ_COL` (new) | — | **the two exposed masses** (small intestine, colon binding sites) |
| (derived) | `EFFECT_SEQ_FREECAP` (new) | — | **the single coupling variable** (free binding capacity) |
| `ka_oca` | `KA0_OCA` | 0.55 | absorption base rate |
| `foca_pass` | `FRACPASS_OCA` | 0.40 | ASBT-independent uptake fraction |
| `kOCAep` | `KEP_OCA` | 1.10 | enterocyte -> plasma |
| `kOCAout` | `KE_OCA` | 0.055 | systemic elimination |
| `kOCAgo` | `KGO_OCA` | 0.15 | luminal transit loss |
| `OCAg`/`OCAe`/`OCAp` (cmt) | `GUT_OCA`/`ENT_OCA`/`CENT_OCA` | — | depot / enterocyte (PD-facing) / plasma (untouched by PD) |
| `EC50_oca` | `EC50_OCA_FXR` | 0.16 (uM) | Hill EC50 |
| (implicit, ratio is linear) | `EMAX_OCA_FXR` — **not added** | — | see "Hill interface" section below: OCA/TRO occupancy is linear, not Emax-saturating per compound |
| (derived) | `C_OCA` (new) | — | **the exposed concentration** (`ENT_OCA/V_ent*1000`, enterocyte site) |
| (derived) | `EFFECT_OCA_FXR` (new) | — | OCA linear contribution to FXR occupancy |
| `ka_tro` | `KA_TRO` | 0.75 | absorption rate |
| `V_trop` | `V1_TRO` | 55.0 (L) | central volume |
| `ke_tro` | -- (removed) | 0.06 | replaced by `CL_TRO`/`V1_TRO` per the guide's CL/V convention |
| -- (promoted) | `CL_TRO` (new) | 3.3 (L/h) | `CL_TRO/V1_TRO = 0.06`, exact — `0.06 * 55 = 3.3` |
| `ftro_ile` | `FPART_TRO` | 0.55 | plasma -> ileal partition |
| `TROg`/`TROp` (cmt) | `GUT_TRO`/`CENT_TRO` | — | depot / central |
| `EC50_tro` | `EC50_TRO_FXR` | 0.0003 (uM) | Hill EC50 |
| (derived) | `C_TRO` (new) | — | **the exposed concentration** (ileal-partitioned, `FPART_TRO*CENT_TRO/V1_TRO`) |
| (derived) | `EFFECT_TRO_FXR` (new) | — | TRO linear contribution to FXR occupancy |
| `IC50_elo` | `EC50_ELO` | 0.35 (uM) | Hill EC50 |
| -- (implicit, ratio max = 1) | `EMAX_ELO` (new) | 1.0 | Hill Emax, "implicit in original ratio `1/(1+C/IC50)`" |
| -- (none) | `GAMMA_ELO` (new) | 1.0 | Hill exponent, "no explicit Hill term in original" |
| `kelo_out` | `KE_ELO` | 0.55 | luminal clearance |
| `ELOg` (cmt) | `CENT_ELO` | — | single compartment, no depot |
| (derived) | `C_ELO` (new) | — | **the exposed concentration** (`CENT_ELO/(3*V_ile)*1000`, shared luminal volume) |
| (derived) | `EFFECT_ELO` (new) | — | fractional ASBT inhibition (reporting only — see below) |
| `ka_lop` | `KA_LOP` | 0.90 | absorption rate |
| `ke_lop` | `KE_LOP` | 0.14 | effect-site loss (no volume in original) |
| `IC50_lop` | `EC50_LOP` | 2.50 (mg) | Hill EC50 |
| -- (implicit) | `EMAX_LOP` (new) | 1.0 | Hill Emax |
| -- (none) | `GAMMA_LOP` (new) | 1.0 | Hill exponent |
| `LOPg`/`LOPe` (cmt) | `GUT_LOP`/`CENT_LOP` | — | depot / effect site |
| (derived) | `C_LOP` (new) | — | **the exposed concentration** (raw amount, mg) |
| (derived) | `EFFECT_LOP` (new) | — | fractional motility inhibition (reporting only) |
| `ke_ond` | `KE_OND` | 0.19 | elimination |
| `IC50_ond` | `EC50_OND` | 3.2 (mg) | Hill EC50 |
| `fOND` | `EMAX_OND` | 0.45 | Hill Emax (already a named max-fraction parameter) |
| -- (none) | `GAMMA_OND` (new) | 1.0 | Hill exponent |
| `ONDc` (cmt) | `CENT_OND` | — | single compartment, no depot |
| (derived) | `C_OND` (new) | — | **the exposed concentration** (raw amount, mg) |
| (derived) | `EFFECT_OND` (new) | — | fractional 5-HT3 blockade |
| `ke_rif` | `KE_RIF` | 0.35 | colonic clearance |
| `IC50_rif` | `EC50_RIF` | 900.0 (mg) | Hill EC50 |
| `fRIF` | `EMAX_RIF` | 0.60 | Hill Emax (already a named max-fraction parameter) |
| -- (none) | `GAMMA_RIF` (new) | 1.0 | Hill exponent |
| `RIFg` (cmt) | `CENT_RIF` | — | single compartment, no depot |
| (derived) | `C_RIF` (new) | — | **the exposed concentration** (raw amount, mg) |
| (derived) | `EFFECT_RIF` (new) | — | fractional suppression of 7alpha-dehydroxylation |

All values are copied verbatim from the original, except `CL_TRO` (a
promoted arithmetic product of two original values, algebraically exact).

## Hill interface

**RIF, OND, LOP, ELO: rename, not a fit.** All four effect terms were
already `Emax*C/(IC50+C)`-shaped (or the algebraically identical
`1/(1+C/IC50)` for ELO/LOP), a plain ratio with an implicit exponent of 1.
`EMAX_<STEM>`/`EC50_<STEM>`/`GAMMA_<STEM>=1` were pulled out with the
original's own values, and `EFFECT_<STEM>` is captured for discoverability.

**OCA and TRO: bespoke — not a per-compound Hill formula.** The original
computes FXR occupancy as `occ = fu_ent*C_ba/EC50_ba + fu_ent*C_oca/EC50_oca
+ C_tro/EC50_tro; FXRi = occ/(1+occ)`. Each compound's own contribution
(`EFFECT_OCA_FXR`, `EFFECT_TRO_FXR`) is a **plain linear term** (`C/EC50`,
no `+C` in the denominator) — the saturating nonlinearity is shared and
competitive across all three FXR ligands (endogenous bile acid, OCA, TRO)
and is applied exactly once, at the combined `occ/(1+occ)` step, not per
compound. This is a genuinely different structure from the guide's
canonical `EMAX*C/(EC50+C)` Hill ratio — forcing that shape onto each
compound individually would either double-count the saturation or silently
remove the competitive-occupancy physiology the original intentionally
modeled (OCA and endogenous bile acid share the same buffered pool via
`fu_ent`; tropifexor reaches FXR unbuffered from plasma, which is exactly
why its effect is phi-independent). Each compound still gets its own named
`EFFECT_<STEM>_FXR` term, and the guide's own "keep each compound's
`EFFECT_<STEM>` separate; combine only where disease equations actually use
them" rule is followed exactly — the combination point (`occ = ... +
EFFECT_OCA_FXR + EFFECT_TRO_FXR`) is the same point the original already
combined them.

**Colesevelam (SEQ): no Hill interface at all.** See the archetype table
above — this is a mass-action consumable-capacity system, not receptor
occupancy, and no `EMAX`/`EC50`/`GAMMA` were invented for it.

## A bisected floating-point subtlety: algebraically-equivalent rearrangements are not bit-identical

The first verification pass found small (~1e-4 absolute, on values in the
hundreds-to-tens-of-thousands range, i.e. ~1e-7 to ~1e-9 relative), bounded,
non-growing mismatches across most of the 56-day trajectory in **every**
scenario tested, including the fully undosed "healthy control" scenario.
Per the guide's "a gate failure is a question, not a verdict" principle,
this was run down rather than accepted or hidden.

**Method.** A dummy no-op identifier rename (`kappa` -> `kappaXYZ`,
touching nothing computational) reproduced **zero** mismatch, ruling out
generic compiler/codegen non-reproducibility across two separate
`mcode_cache` compiles. A cumulative bisection — reapplying this refactor's
own edits one block at a time and re-running the undosed scenario after
each — isolated the mismatch to exactly one class of edit: rewriting
`Ielo = 1.0/(1.0 + C/EC50)` as `Ielo = 1.0 - EFFECT` where `EFFECT =
EMAX*C/(EC50+C)`. These two forms are algebraically identical for all `C`,
but **not IEEE754 bit-identical** — and `CENT_ELO`, though never dosed in
this scenario, is not held at an exact analytic `0.0` by the adaptive ODE
solver (a pure-decay state started at 0 can still carry solver-level
floating-point noise at the ~1e-16 scale). The sub-ULP difference between
the two rearrangements, applied to that near-zero noise, is tiny on its
own, but this model's own feedback loop (FGF19 -> CYP7A1 -> bile-acid pool,
exactly the "loop gain" the model's own organising idea is about) amplifies
it — boundedly, never divergently, confirmed by checking that the maximum
mismatch never grows across the full 1344 h / 5377-point trajectory — into
a visible last-digit mismatch on `POOLtot`/`FGF19`/`GB`/`FEC`/etc. The same
mechanism, from the same class of rearrangement, was independently found
and confirmed for tropifexor (`FPART_TRO*CENT_TRO/V1_TRO` vs.
`(CENT_TRO/V1_TRO)*FPART_TRO` — reordering a multiply and a divide is not
bit-identical either) once ELO's fix was in place and a standalone
tropifexor-dosed scenario was tested.

**Fix.** For every compound whose original effect formula was not already
exactly `EMAX*C/(EC50+C)` (ELO, LOP: originally `1/(1+C/EC50)`; TRO:
originally `ftro_ile*TROp/V_trop`, a specific multiply-then-divide order),
the value that actually drives the disease ODEs (`Ielo`, `LOP`, `C_TRO`)
keeps the **exact expression and operation order the original already
used**, with only identifiers renamed. The canonical, named
`EMAX_<STEM>/EC50_<STEM>/GAMMA_<STEM>`-based `EFFECT_<STEM>` is still
computed and captured for discoverability, but as a **separate,
non-load-bearing** quantity that never feeds back into any `dxdt_`. For
RIF and OND, whose original formula was already exactly `EMAX*C/(EC50+C)`,
no rearrangement was ever needed — the load-bearing expression is a pure
identifier rename. This is disclosed here in full because it is exactly
the kind of subtlety the guide's tolerance section anticipates, not a
tolerance loosened to force a pass: the fix was verified, not assumed, by
re-running every scenario below after each change.

## A pre-existing upstream build defect (autodec multi-declarator), unrelated to any of the seven compounds

`bam_mrgsolve_model.R` does not compile under mrgsolve 2.0.1 as written.
Same defect class as issues #37/#61/#66 (multi-declarator `double`
statements): `$PLUGIN autodec`'s `$ODE` preprocessing strips `double` from
every comma-separated declarator but the first, so

```
double J12 = kt*ILE1, J23 = kt*ILE2, J3c = kile*TRANS*ILE3;
double fb1 = ILE1/ILEtot, fb2 = ILE2/ILEtot, fb3 = ILE3/ILEtot;
```
and
```
double fc1 = CCA/tot_col, fc2 = CCDCA/tot_col;
double fc3 = CDCA/tot_col, fc4 = CLCA/tot_col;
```

give `'J23' was not declared in this scope`, `'fb2' was not declared`,
`'fc2' was not declared`, `'fc4' was not declared`. All four lines sit in
core disease-side scaffolding (ileal transit fluxes, colonic pool
fractions) — none is inside any of the seven refactored compounds' own
blocks, so per the guide's policy this is a syntax-only, non-numeric
build-compatibility fix, logged as **`UPSTREAM_ISSUES.md` #75**, and
applied identically to the delivered `bam_mrgsolve_model_refactored.R`:
each multi-declarator line was split into separate `double` statements,
same values, same order, no behavior changed. Confirmed by the exact-match
verification below (both the untouched original, patched in-memory with
only this fix for the verification call, and the delivered refactored
file carry this same fix).

## Verification

**Method.** Both files' embedded mrgsolve `code <- '...'` DSL blocks were
mechanically extracted (regex on the assignment, verbatim quoted text; no
separate `.cpp` file left behind — extraction was in-memory/scratch-only)
and POSTed to the local qspserver `mrgsolve_api` service at
`http://localhost:8007` (`/model_manifest` then `/run_simulation`), which
compiles and runs each DSL block directly with mrgsolve 2.0.1 server-side —
no local R/mrgsolve install used. Requests were spaced >=2 s apart,
respecting the service's stated `max_concurrent_jobs: 2` limit and its
recent-crash history. The original's copy was additionally patched
in-memory with the autodec fix above (never written back to the tracked
`bam_mrgsolve_model.R`).

**Results: exact match (max abs deviation 0.0) on every scenario tested**,
full 56-day window (`end=1344h`, `delta=0.25h`, matching the original's own
`sim_bam()` defaults — no shortening needed, well within the API's default
`maxsteps`), across all 44 shared disease-side `$CMT`/`$CAPTURE` outputs
plus all 13 renamed drug compartments:

1. All 12 of `bam_scenarios()`'s own named scenarios (`01`–`14`, including
   the central FEEDBACK-live-vs-frozen experiment, all four aetiological
   types, colesevelam/OCA/loperamide/elobixibat dosing, and their
   combination).
2. Standalone ondansetron (4 mg BID), rifaximin (550 mg TID), and
   tropifexor (0.090 mg QD) at each compound's own `ev_ond()`/`ev_rif()`/
   `ev_tro()` default dose — **the original's own `bam_scenarios()` never
   actually doses these three compounds**, despite defining working helper
   functions for all three; this is the only way to verify OND/RIF/TRO's
   renamed PK+effect blocks against real, non-zero dosing at all, so they
   were exercised here for completeness on every compound the census now
   tracks.

Every one of the above 17 comparisons returned **max abs diff = 0.0**
across every shared output — exact, not merely close — once the bisected
floating-point issue above was fixed. Compartment indices are identical
between the original and the refactored model (compartment order in `$CMT`
is unchanged, only names changed), so the same 1-based `cmt=` dosing
indices work for both.

## qspserver `/model_manifest` discoverability

Confirmed via `POST /model_manifest` against the extracted, renamed DSL:
all renamed PK parameters (`CAP_SEQ`, `KBIND_SEQ`, `K1_SEQ`, `K2_SEQ`,
`K3_SEQ`, `KA0_OCA`, `FRACPASS_OCA`, `KEP_OCA`, `KE_OCA`, `KGO_OCA`,
`EC50_ELO`, `EMAX_ELO`, `GAMMA_ELO`, `KE_ELO`, `KA_TRO`, `V1_TRO`,
`CL_TRO`, `FPART_TRO`, `KA_LOP`, `KE_LOP`, `EC50_LOP`, `EMAX_LOP`,
`GAMMA_LOP`, `KE_OND`, `EC50_OND`, `EMAX_OND`, `GAMMA_OND`, `KE_RIF`,
`EC50_RIF`, `EMAX_RIF`, `GAMMA_RIF`, `EC50_OCA_FXR`, `EC50_TRO_FXR`) appear
in the manifest's `parameters` with their original numeric defaults (119
parameters total). `C_<STEM>` and `EFFECT_<STEM>` for every compound are
state-derived (computed in `$ODE` from compartment values, several
deliberately kept as non-load-bearing reporting-only quantities per the
floating-point fix above) — per the same reasoning established for
`EFFECT_COL_NEU`/`EFFECT_ANA`/`EFFECT_CANA` in
`familial-mediterranean-fever/fmf_refactor_notes.md` and
`C_OCA`/`C_UDCA`/`C_BEZ` in `primary-sclerosing-cholangitis/psc_refactor_notes.md`,
they cannot also be `$PARAM` entries; all 15 (`C_SEQ_SI`, `C_SEQ_COL`,
`EFFECT_SEQ_FREECAP`, `C_OCA`, `EFFECT_OCA_FXR`, `C_TRO`, `EFFECT_TRO_FXR`,
`C_ELO`, `EFFECT_ELO`, `C_LOP`, `EFFECT_LOP`, `C_OND`, `EFFECT_OND`,
`C_RIF`, `EFFECT_RIF`) appear in the manifest's `outputPaths` via the
file's own `$CAPTURE @annotated` list, confirmed discoverable, alongside
all 13 renamed drug compartments and 31 unrenamed disease-side
compartments/captures.

## Anything else flagged

- The R-side dosing helpers (`ev_seq`, `ev_oca`, `ev_tro`, `ev_elo`,
  `ev_lop`, `ev_ond`, `ev_rif`) and `sim_bam()`'s placebo event were updated
  to the renamed `cmt=` targets (`"SEQ1"` -> `"GUT_SEQ"`, `"OCAg"` ->
  `"GUT_OCA"`, `"TROg"` -> `"GUT_TRO"`, `"ELOg"` -> `"CENT_ELO"`, `"LOPg"`
  -> `"GUT_LOP"`, `"ONDc"` -> `"CENT_OND"`, `"RIFg"` -> `"CENT_RIF"`). Same
  dosing amounts, same timing, same scenario list, same `bam_readout()`
  post-processing throughout — checked by grepping the full original file
  for every old PK identifier outside the DSL block.
- `mcode_cache("bam_qsp", ...)` was renamed to `mcode_cache("bam_qsp_refactored", ...)`
  so the two files never collide in mrgsolve's on-disk model cache.
- No disease-side equation, parameter, or compartment was touched beyond
  the autodec syntax fix above; `phi_from_resection()`, `bam_readout()`,
  and every calibration-provenance comment are unchanged.

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`: the five
existing `bile-acid-diarrhea` rows were corrected to real drug names
(`Loperamide (LOP)`, `Ondansetron (OND)`, `Rifaximin (RIF)`, `Obeticholic
acid (OCA)`, `Tropifexor (TRO)`), and two new rows were added for the
previously-uncounted `Elobixibat (ELO)` and `Colesevelam (SEQ)`.
