# Refactor notes — `essential-tremor/et_mrgsolve_model.R`

**Scope of this pass.** Per the fork's PK/PD refactor spec
(`FORK_WORKFLOW_GUIDE.md`, Part 2) and the task instructions, all eleven
of this file's rows in `driver-patches/data/compound_perturbation_census.md`
were refactored: **propranolol (PRP)**, **atenolol (ATN)**, **nadolol
(NAD)**, **primidone (PRM)**, **phenobarbital (PB)**, **PEMA (PEM)**,
**topiramate (TOP)**, **gabapentin (GBP)**, **ethanol (ETH)**,
**1-octanol (OCT)**, and a **T-type/Cav3.1 calcium-channel blocker
(TTB)**. Two of the census rows were process-description phrases rather
than drug names — corrected below. Everything else in the file (the
oscillator core, biomechanics, disease/adaptation states, Vim
lesion/DBS, botulinum toxin `A_BTXT`/`A_BTXG`/`SNAPT`/`SNAPG`, and every
organ-system readout) is untouched apart from the three build-compat
renames disclosed under "Upstream defects" below.

## The two mislabeled census rows

The task flagged **"Organ systems (PRM)"** and **"Organ systems (PRP)"**
as suspicious — process-description phrases, not drug names, paired with
a correct stem in parentheses (the same pattern documented for
`neonatal-hyperbilirubinemia`'s "Bilirubin (SNMP)" / "Hepatic handling
(PB)" rows). Checking the code:

- **PRM is primidone.** `MW_PRM = 218.25` matches primidone's molecular
  weight; it has its own oral depot/absorption (`KA_PRM`, `F_PRM`, gut
  compartment `A_PRMG`) and is the compound the file's own scenarios S8/S9
  dose directly (`reg("A_PRMG", 250/3, ii=8)`, "primidone: the parent
  molecule is the active moiety").
- **PRP is propranolol.** `MW_PRP = 259.3` matches propranolol's
  molecular weight; it has beta1/beta2 receptor-affinity parameters
  (`KI_PRP_B1`, `KI_PRP_B2`), a dose-dependent first-pass-escape
  bioavailability term (`F0_PRP + FMX_PRP*DPRP/(FD50_PRP+DPRP)`, matching
  propranolol's well-known saturable presystemic CYP2D6/1A2 metabolism),
  and is the compound dosed in scenario S1 ("beta-blockade, and why
  beta1-selectivity fails") and named explicitly throughout the file's own
  comments ("Propranolol acts on beta2 receptors on the muscle SPINDLE").

**Conclusion: both are real, externally-dosed drugs, correctly targeted by
their stems — only the census's display name was a process description.**
Corrected to "Primidone (PRM)" and "Propranolol (PRP)" in the census (see
below).

## Other identities verified against the code (not assumed)

- **ATN = atenolol** (`MW_ATN = 266.3`; `KI_ATN_B2 = 1000` nM, beta1-
  selective; no brain-biophase compartment anywhere for it, consistent
  with atenolol's known poor CNS penetration — a real pharmacological fact
  the model already encodes structurally, not something this refactor
  added).
- **NAD = nadolol** (`MW_NAD = 309.4`; non-selective, peripherally
  restricted — same "no brain compartment" pattern as atenolol).
- **ETH = ethanol** — genuinely dosed, not a non-drivable trigger. It has
  its own oral depot (`A_ETHG`, g), a saturable (Michaelis-Menten) central
  elimination, and a dynamic brain-biophase compartment; the file's own
  scenarios S12-S15 dose it directly via `ev(amt=.., cmt="A_ETHG")`. Kept
  as a fully-dosed compound archetype, not corrected to a trigger-only
  row.
- **OCT = 1-octanol** (`MW_OCT = 130.2`, matching 1-octanol's C8H18O
  formula weight; the file's own comments name it explicitly, "1-octanol:
  tremor benefit per unit intoxication", and reference its real clinical
  trial history for essential tremor).
- **PEM = PEMA** (phenylethylmalonamide), primidone's second metabolite —
  see "Bespoke: primidone's metabolite pair" below.
- **TTB = a T-type/Cav3.1 calcium-channel blocker.** The file's own
  scenario comments name the real compound this represents explicitly:
  "The species discordance that killed CX-8998's rationale is ONE
  parameter" (CX-8998, later PRAX-944, is Cavion/Praxis Precision
  Medicines' Cav3.1-selective T-type calcium-channel blocker studied for
  essential tremor). The in-code stem `TTB` is kept (per the guide,
  "keyed to its existing abbreviation stem... don't invent a new stem").

## Archetype per compound

| Compound | Archetype | Notes |
|---|---|---|
| **PRP** (propranolol) | Archetype 3 (depot+central+peripheral) **plus a bespoke dynamic biophase (brain-exposure) compartment** | `GUT_PRP`→`CENT_PRP`↔`PERI_PRP` is textbook Archetype 3; `C_PRPB` (renamed `CE_PRP`) is a genuine first-order Sheiner-type effect-site compartment (`dxdt_CE_PRP = KE0_PRP*(KPUU_PRP*CF_PRP - CE_PRP)`), not covered by the guide's four-slot table — same "extra role beyond the guide's four-slot table" situation documented for `LIV_GIV`/`LIV_HEM` (`acute-intermittent-porphyria`) and brain-exposure compartments in `huntingtons-disease`. Named `CE_<STEM>` throughout (see naming table). |
| **ATN** (atenolol) | Archetype 3 minus peripheral (depot+central) | No brain compartment (real atenolol CNS penetration is poor) — nothing bespoke. |
| **NAD** (nadolol) | Archetype 3 minus peripheral (depot+central) | Same shape as ATN. |
| **PRM** (primidone) | Archetype 3 minus peripheral (depot+central) **plus the same bespoke dynamic biophase pattern as PRP**, **plus is the source of a genuine metabolite pair (PB, PEM)** | `GUT_PRM`→`CENT_PRM` plus `CE_PRM` (dynamic biophase). Primidone's own elimination rate (`KE_PRM = CL_PRM/V1_PRM`) is split by fixed fractions (`FM_PB`, `FM_PEM`) into phenobarbital's and PEMA's own central compartments — see below. |
| **PB** (phenobarbital) | **Bespoke**: Archetype-1-like single compartment (`CENT_PB`) fed entirely by metabolic conversion from `CENT_PRM`, no `GUT_PB`, no independent dosing route in this file — **plus its own dynamic biophase compartment** (`CE_PB`) | Same pattern as `OXO` (4-oxo-isotretinoin, `acne-vulgaris`) and `HTBZ` (huntingtons-disease): a real metabolite with its own genuine elimination (`CL_PB`/`V1_PB`), not merely a scaled copy of the parent. Clinically phenobarbital *can* be dosed independently, but this file only ever produces it via primidone metabolism — disclosed as a modelling choice of the original, not "fixed" to add independent dosing (that would invent PK the original never had). |
| **PEM** (PEMA) | **Bespoke**: same pattern as PB — Archetype-1-like single compartment (`CENT_PEM`) fed by metabolic conversion from `CENT_PRM`, no `GUT_PEM`, no independent dosing route | No dynamic biophase compartment for PEM (its GABA-A effect reads the unit-converted plasma concentration `CU_PEM` directly, an algebraic, not dynamic, quantity). |
| **TOP** (topiramate) | Archetype 3 minus peripheral (depot+central) | Brain-equivalent concentration `CB_TOP` is an **algebraic** (instantaneous) partition (`KP_TOP*CU_TOP`), not a dynamic compartment — distinct from PRP/PRM/PB/TTB's genuine KE0-driven biophase kinetics. Two independent named effects (cortical, GABA — see Hill interface below), matching the `STAT` precedent (`abdominal-aortic-aneurysm`) for a single compound with multiple independently-named downstream effects. |
| **GBP** (gabapentin) | Archetype 3 minus peripheral (depot+central) | Same algebraic-partition pattern as TOP (`CB_GBP = KP_GBP*CU_GBP`); one named cortical effect. |
| **ETH** (ethanol) | **Bespoke**: depot+central with **saturable (Michaelis-Menten) central elimination**, not first-order CL — **plus a dynamic biophase compartment** | `dxdt_CENT_ETH = KA_ETH*GUT_ETH - VMAX_ETH*C_ETH/(KM_ETH+C_ETH)` does not fit any of the guide's four linear archetypes; kept exactly as the original computed it (guide: "don't add or remove compartments to make it fit one"). Three independent named effects (GABA, olivary, cortical — see below). |
| **OCT** (1-octanol) | Archetype 3 minus peripheral (depot+central) | Same algebraic-partition pattern as TOP/GBP (`CB_OCT = KP_OCT*CU_OCT`). Two independent named effects (GABA, olivary). |
| **TTB** (T-type/Cav3.1 blocker) | Archetype 3 minus peripheral (depot+central) **plus the same bespoke dynamic biophase pattern as PRP/PRM/PB** | `GUT_TTB`→`CENT_TTB` plus `CE_TTB` (dynamic biophase, `KE0_TTB`-driven). One named effect (olivary Cav3.1 block). |

Four compounds (PRP, PRM, PB, TTB) share the exact same first-order,
KE0-driven "plasma → brain-free" biophase-equilibration structure — a
recurring motif in this file, not a one-off. None of the guide's four
named roles (`GUT_`/`CENT_`/`PERI_`/`REC_FREE_`) cover it, so it is named
`CE_<STEM>` (concentration, effect-site) throughout, distinct from the
true bidirectional, mass-conserving `PERI_PRP` propranolol already has for
its ordinary peripheral tissue. `huntingtons-disease` used `PERI_<STEM>`
for an analogous one-way brain-exposure role, but that file had no
competing true peripheral compartment to collide with; here PRP has both
a genuine `PERI_PRP` (Q/V2 exchange) and this separate biophase link, so
reusing `PERI_` would have been ambiguous — `CE_<STEM>` is used instead
and disclosed here as the bespoke naming choice.

## The beta-adrenoceptor occupancy system (PRP + ATN + NAD): a bespoke,
## non-Hill shared-receptor structure

Propranolol, atenolol, and nadolol are all competitive antagonists at the
**same** peripheral beta-adrenoceptor pool. The original computes this as
a Gaddum-equation competitive equilibrium: each drug's own free
concentration divided by its own receptor affinity (`Ki`) sums into one
shared occupancy fraction (`OCCB2`, `OCCB1`), which every downstream
tremor/HR/SBP/FEV1 effect reads. This is **not** the canonical
`EMAX*Xᵞ/(EC50ᵞ+Xᵞ)` Hill shape — it's a linear competitive quotient
combined *before* the shared saturation, structurally closer to TMDD's
"build the interface from occupancy, not raw concentration" guidance than
to a plain per-drug Hill term.

Per the guide's "combine them only at the point the disease equations
actually use them" (never collapse several drugs into one shared term
prematurely), each drug's own quotient is now its own named variable —
`EFFECT_PRP_B2 = CF_PRP/KI_PRP_B2`, `EFFECT_ATN_B2 = CF_ATN/KI_ATN_B2`,
`EFFECT_NAD_B2 = CF_NAD/KI_NAD_B2` (and the `_B1` equivalents) — summed
into `IQ_B2`/`IQ_B1` only at that point, exactly where the original
already summed them (a rename/reorganization, not a new algebraic shape).
No `GAMMA_<STEM>` was added to these six terms — they are disclosed as a
bespoke, non-Hill "occupancy quotient" role, not forced into the Hill
template it doesn't fit.

## The Hill interface: rename (gamma=1), not a fit, for every compound

Every compound's disease-relevant effect term in this file was already
exactly the plain Hill/Emax ratio shape (`C/(EC50+C)`, no explicit
exponent) — so every one is a rename, per the guide's explicit
instruction, with an explicit `GAMMA_<STEM>[_role] = 1.0` parameter added
(none of the originals had one) and the canonical
`EMAX*pow(X,GAMMA)/(pow(EC50,GAMMA)+pow(X,GAMMA))` form used throughout
(mathematically identical to the plain ratio at gamma=1). No `nls()` fit
was needed for any of the eleven compounds.

Several compounds have **more than one** independently-named effect —
kept separate per the guide ("Multiple drugs, one pathway: keep each
compound's `EFFECT_<STEM>` separate; combine them only at the point the
disease equations actually use them"), following the same pattern already
used for STAT's three independent effects (`abdominal-aortic-aneurysm`):

| Compound | Named effect(s) | Feeds |
|---|---|---|
| PRP | `EFFECT_PRP_CBL` | cerebello-thalamic branch (`PHI_CBL`) |
| PRM | `EFFECT_PRM_GABA`, `EFFECT_PRM_NACH` | GABA-A potentiation (`P_RAW`), cortical Na-channel term (`PHI_CTX`) |
| PB | `EFFECT_PB_GABA` | GABA-A potentiation |
| PEM | `EFFECT_PEM_GABA` | GABA-A potentiation |
| TOP | `EFFECT_TOP_GABA`, `EFFECT_TOP_CTX` | GABA-A potentiation, cortical term |
| GBP | `EFFECT_GBP` | cortical term |
| ETH | `EFFECT_ETH_GABA`, `EFFECT_ETH_OL`, `EFFECT_ETH_CTX` | GABA-A potentiation, olivary branch, cortical term |
| OCT | `EFFECT_OCT_GABA`, `EFFECT_OCT_OL` | GABA-A potentiation, olivary branch |
| TTB | `EFFECT_TTB_OL` | olivary branch (Cav3.1 block) |

**The one place the original genuinely combined two drugs into a shared
term (`ETH_OL`) was reorganized, not force-split apart.** The original
computed `ETH_OL = EMAX_ETHO*C_ETHB/(EC50_ETHO+C_ETHB) +
EMAX_OCTO*CB_OCT/(EC50_OCTO+CB_OCT)`, then capped the *sum* at 0.95 — a
genuine shared physiological ceiling (both act on the same olivary
Cav3.1-adjacent mechanism and cannot together exceed near-total
suppression), not a shortcut. Per the guide, this is now `EFFECT_ETH_OL`
and `EFFECT_OCT_OL` as two independent named terms, summed and capped
**only at the point PHI_OL is computed** — same numeric result, but each
drug's own term is now independently driveable before the shared cap is
applied.

## Magic-number promotions (same values, now named/discoverable)

Six inline literal EC50-like constants that were never named parameters
in the original were promoted to explicit `$PARAM` entries, following the
precedent set for STAT's NF-kB/ROS terms (`abdominal-aortic-aneurysm`,
"hardcoded magic-number ratios... were promoted to named Emax/EC50
params, same values, same shape"):

| Original inline literal | New parameter | Value | Used in |
|---|---|---|---|
| `300.0` (propranolol term in `SED_RAW`) | `EC50_SED_PRP` | 300.0 | Sedation index |
| `4.0` (gabapentin term in `SED_RAW`) | `EC50_SED_GBP` | 4.0 | Sedation index |
| `30.0` (topiramate term in `SED_RAW`) | `EC50_SED_TOP` | 30.0 | Sedation index |
| `0.90` (ethanol term in `SED_RAW`) | `EC50_SED_ETH` | 0.90 | Sedation index |
| `120.0` (phenobarbital term in `ATAXd`) | `EC50_ATAX_PB` | 120.0 | Ataxia index |
| `1.20` (ethanol term in `ATAXd`) | `EC50_ATAX_ETH` | 1.20 | Ataxia index |

These feed clinical-endpoint biomarkers (sedation, ataxia), not the
tremor-gain Hill interface itself, so they were promoted to named
parameters (for manifest discoverability) but **not** converted to the
`pow()`-based canonical Hill form — same as the original's own already-
named `EC50_SED_PB`/`EC50_SED_PRM` neighbours in the same expressions,
which were left exactly as they were.

## Full renaming table

| Role | Original | Refactored |
|---|---|---|
| PRP depot/central/peripheral | `A_PRPG`/`A_PRPC`/`A_PRPP` | `GUT_PRP`/`CENT_PRP`/`PERI_PRP` |
| PRP biophase (bespoke) | `C_PRPB` | `CE_PRP` |
| ATN depot/central | `A_ATNG`/`A_ATNC` | `GUT_ATN`/`CENT_ATN` |
| ATN volume | `V_ATN` | `V1_ATN` |
| NAD depot/central | `A_NADG`/`A_NADC` | `GUT_NAD`/`CENT_NAD` |
| NAD volume | `V_NAD` | `V1_NAD` |
| PRM depot/central | `A_PRMG`/`A_PRMC` | `GUT_PRM`/`CENT_PRM` |
| PRM volume, elimination | `V_PRM`, `kel_prm` | `V1_PRM`, `KE_PRM` |
| PRM biophase (bespoke) | `C_PRMB` | `CE_PRM` |
| PB central (bespoke, metabolite) | `A_PBC` | `CENT_PB` |
| PB volume | `V_PB` | `V1_PB` |
| PB biophase (bespoke) | `C_PBB` | `CE_PB` |
| PEM central (bespoke, metabolite) | `A_PEMC` | `CENT_PEM` |
| PEM volume | `V_PEM` | `V1_PEM` |
| TOP depot/central | `A_TOPG`/`A_TOPC` | `GUT_TOP`/`CENT_TOP` |
| TOP volume | `V_TOP` | `V1_TOP` |
| GBP depot/central | `A_GBPG`/`A_GBPC` | `GUT_GBP`/`CENT_GBP` |
| GBP volume | `V_GBP` | `V1_GBP` |
| ETH depot/central (bespoke, saturable) | `A_ETHG`/`A_ETHC` | `GUT_ETH`/`CENT_ETH` |
| ETH volume | `V_ETH` | `V1_ETH` |
| ETH biophase (bespoke) | `C_ETHB` | `CE_ETH` |
| OCT depot/central | `A_OCTG`/`A_OCTC` | `GUT_OCT`/`CENT_OCT` |
| OCT volume | `V_OCT` | `V1_OCT` |
| TTB depot/central | `A_TTBG`/`A_TTBC` | `GUT_TTB`/`CENT_TTB` |
| TTB volume | `V_TTB` | `V1_TTB` |
| TTB biophase (bespoke) | `C_TTBB` | `CE_TTB` |
| PRP cerebellar Emax/EC50 | `EMAX_PRPC`/`EC50_PRPC` | `EMAX_PRP_CBL`/`EC50_PRP_CBL` |
| PRM GABA Emax/EC50 | `EMAX_PRM`/`EC50_PRM` | `EMAX_PRM_GABA`/`EC50_PRM_GABA` |
| PRM cortical (Na-channel) Emax/EC50 | `EMAX_NACH`/`EC50_NACH` | `EMAX_PRM_NACH`/`EC50_PRM_NACH` |
| PB GABA Emax/EC50 | `EMAX_PB`/`EC50_PB` | `EMAX_PB_GABA`/`EC50_PB_GABA` |
| PEM GABA Emax/EC50 | `EMAX_PEM`/`EC50_PEM` | `EMAX_PEM_GABA`/`EC50_PEM_GABA` |
| TOP cortical Emax/EC50 | `EMAX_TOPC`/`EC50_TOPC` | `EMAX_TOP_CTX`/`EC50_TOP_CTX` |
| TOP GABA Emax/EC50 | `EMAX_TOPG`/`EC50_TOPG` | `EMAX_TOP_GABA`/`EC50_TOP_GABA` |
| ETH GABA Emax/EC50 | `EMAX_ETHG`/`EC50_ETHG` | `EMAX_ETH_GABA`/`EC50_ETH_GABA` |
| ETH olivary Emax/EC50 | `EMAX_ETHO`/`EC50_ETHO` | `EMAX_ETH_OL`/`EC50_ETH_OL` |
| ETH cortical Emax/EC50 | `EMAX_ETHC`/`EC50_ETHC` | `EMAX_ETH_CTX`/`EC50_ETH_CTX` |
| ETH intoxication EC50 | `EC50_INT` | `EC50_INTOX_ETH` |
| OCT GABA Emax/EC50 | `EMAX_OCTG`/`EC50_OCTG` | `EMAX_OCT_GABA`/`EC50_OCT_GABA` |
| OCT olivary Emax/EC50 | `EMAX_OCTO`/`EC50_OCTO` | `EMAX_OCT_OL`/`EC50_OCT_OL` |
| OCT intoxication EC50 | `EC50_INTO` | `EC50_INTOX_OCT` |
| TTB olivary Emax/EC50 | `EMAX_TT`/`IC50_TT` | `EMAX_TTB_OL`/`EC50_TTB_OL` (IC50→EC50 per convention, same as SNMP in `neonatal-hyperbilirubinemia`) |
| GBP bioavailability override | `F_A_PRPG` etc. | `F_GUT_PRP` etc. (five compounds: PRP, ATN, NAD, PRM, TOP — the only five with an explicit bioavailability override in `$MAIN`) |
| Effect terms (rename) | `P_PRM`, `P_PB`, `P_PEM`, `P_ETH`, `P_OCT`, `P_TOP` | `EFFECT_PRM_GABA`, `EFFECT_PB_GABA`, `EFFECT_PEM_GABA`, `EFFECT_ETH_GABA`, `EFFECT_OCT_GABA`, `EFFECT_TOP_GABA` |
| Cav3.1 block term | `BLK_TT` | `EFFECT_TTB_OL` |

All parameter **values** are copied verbatim from the original — nothing
was invented, defaulted, or refit.

## `F_TTB` — a pre-existing dead parameter, disclosed but not "fixed"

`F_TTB = 0.70` is declared in `$PARAM` but **never used anywhere** in the
original — there is no `F_A_TTBG` (now `F_GUT_TTB`) override in `$MAIN`
the way there is for propranolol/atenolol/nadolol/primidone/topiramate,
so mrgsolve's default bioavailability (1.0) silently applies to every
T-type-blocker dose regardless of this parameter's value. This is a real
defect in the original (a stated 70% oral bioavailability that the model
never actually applies), but per the guide's near-exact-match requirement
for a linear archetype, **adding `F_GUT_TTB = F_TTB;` to `$MAIN` would
change TTB's numeric behaviour** (introducing a 30% bioavailability cut
that was never actually in effect), which is out of scope for a rename/
reorganization refactor. `F_TTB` is therefore carried forward, still
declared, still unused, exactly as in the original — disclosed here
rather than silently fixed or silently dropped.

## Upstream defects (build-compat, logged as `UPSTREAM_ISSUES.md` #105)

The untouched original does not compile under mrgsolve 2.0.1. Three
independent defects, each surfacing only after the previous was fixed:

1. **`$PARAM @annotated` line `EMAX_DBS: 0.90`** has two colon-separated
   fields where the parser requires three (name : value : description).
   Fixed by adding the missing description (`: DBS max block fraction`).
   Unrelated to any of the eleven compounds (DBS/lesion section).
2. **`EPS` is a reserved word** under mrgsolve 2.0.1 (`"Reserved words in
   model names: EPS"`). Renamed `EPS_NOISE` everywhere (one declaration,
   three uses in the oscillator's noise-seeding term). Unrelated to any
   compound (the tremor oscillator's own noise floor).
3. **`F_PRP` is used in `$MAIN` but never declared**, a genuine C++
   compile failure (`'F_PRP' was not declared in this scope`). Fixed by
   adding it to the file's existing `$GLOBAL` double-declaration list.
   This one **is** inside propranolol's own block — per the guide's
   allowance for a defect genuinely inside the scope compound's own
   code, it is fixed as part of this refactor's own archetype work for
   PRP (not treated as a foreign incidental patch), while still logged
   as a real pre-existing defect in the checked-in original.

All three are syntax-only, non-numeric, name/declaration-only changes.
Full detail, reproduction steps, and the exact fixes: `UPSTREAM_ISSUES.md`
entry #105. The tracked `et_mrgsolve_model.R` is completely untouched and
still carries all three defects exactly as written.

## Verification

**Method.** Both DSL blocks (the original's own, with the three build-
compat fixes above applied to make it compile at all — no compound-level
change — and the refactored sibling's) were extracted from their
respective `code <- '...'` R strings and POSTed to the local qspserver
`mrgsolve_api` service at `http://localhost:8007` (`/model_manifest` then
`/run_simulation`), which compiles and runs each DSL block directly with
mrgsolve 2.0.1 server-side. Requests were spaced roughly two seconds
apart, one in flight at a time, respecting the service's
`max_concurrent_jobs: 2` limit. Every step of the refactor (pure renaming
alone, then the structural reorganization on top of that) was verified
independently before being combined, and the final artifact was verified
once more by re-extracting its `code <- '...'` block from the finished
`et_mrgsolve_model_refactored.R` itself (not a separately-tracked copy).

**Scenarios run.** Nine of the file's own native dosing scenarios (S1,
S8, S12, S14, S15 ×2, S16), spanning all eleven compounds except
topiramate and gabapentin — full duration exactly as coded (no
shortening needed; none approached the API's default `maxsteps` budget):

1. **S1 propranolol** — 160 mg qd × 24 weeks (`ii=24, addl=167`,
   `end=4032h`): dosed into `GUT_PRP` (index 1; `A_PRPG` in the original).
2. **S1 atenolol** — 100 mg qd × 24 weeks: `GUT_ATN` (index 5).
3. **S1 nadolol** — 120 mg qd × 24 weeks: `GUT_NAD` (index 7).
4. **S8 primidone** — 83.33 mg (250/3) q8h × 24 weeks (`ii=8, addl=503`):
   `GUT_PRM` (index 9) — exercises the full PRM→PB/PEM metabolite chain.
5. **S12 ethanol** — 28 g single bolus, `end=40h`: `GUT_ETH` (index 19).
6. **S14 daily ethanol** — 42 g qd × 90 doses starting at t=19h
   (`ii=24, addl=89`), `end=2184h` (91 days): exercises the chronic
   adaptation/rebound (`ADAPTF`/`ADAPTS`) pathway.
7. **S15 1-octanol** — 560 mg and 1120 mg single bolus, `end=24h` each:
   `GUT_OCT` (index 22).
8. **S16 T-type blocker** — 100 mg qd × 12 weeks (`ii=24, addl=83`),
   `end=2016h`: `GUT_TTB` (index 24).

All eleven output names captured in both files (`GTOT`, `MU`, `PHI_C`,
`PHI_P`, `PHI_OL`, `PHI_CBL`, `PHI_TH`, `PHI_SPIN`, `PHI_NMJ`, `A_UL`,
`A_LC`, `A_PHYS`, `A_HD`, `A_VX`, `FNEUR`, `F0M`, `F_OBS`, `T_R`,
`TETRAS_PS`, `TETRAS_ADL`, `FTM`, `SPIRAL`, `BAINF`, `QUEST`, `OCCB2`,
`OCCB1`, `P_RAW`, `P_EFF`, the primidone/phenobarbital GABA terms, `REB`,
the T-type block term, `LES_EFF`, `DBSB`, `GRIP`, `HR`, `SBP`, `SED`,
`ATAX`, `INTOX`, and all five plasma concentrations the original already
captured) matched with **max abs diff 0.0** (bit-exact) in every scenario.

**Topiramate and gabapentin have no native scenario.** Neither the
original file nor any of its 25 shipped scenarios ever doses `A_TOPG` or
`A_GBPC` — both compounds have full depot+central PK and named Hill
effects wired into the disease equations, but are never exercised. Per
the guide's verification requirement ("run every dosing scenario already
defined in the original's own file's own R code — not invented ones"),
there is no native scenario to run for either. As an additional,
disclosed-as-non-native sanity check (not a substitute for a native
scenario, since none exists), both were dosed identically into the
original and refactored DSL (topiramate 100 mg qd × 10 days into
`GUT_TOP`/index 15; gabapentin 300 mg q8h × 10 days into `GUT_GBP`/index
17) — **max abs diff 0.0** in both, confirming the renamed/restructured
cortical and GABA-A terms compute identically to the original's inline
formulas when actually exercised.

**Tolerance.** All eleven compounds are pure structural reorganization
(rename + reorganize; no compound's effect term required an `nls()` fit,
since every one was already the plain Hill ratio shape). Per the guide's
tolerance table for Archetypes 1–3 (and the file's bespoke variants
thereof), an exact match was expected — and obtained: **0.0 max absolute
difference across all eleven verification scenarios**, including after
introducing `pow(x, 1.0)` calls for the canonical Hill form (this
particular build's `pow()` returns bit-identical results for exponent
1.0, so — unlike the floating-point-scale deviation disclosed for a
similar rewrite in `neonatal-hyperbilirubinemia`'s Stannsoporfin term —
no residual deviation appears here at all).

## qspserver `/model_manifest` discoverability

Confirmed via `POST /model_manifest` against the extracted, renamed DSL
(248 parameters, 119 output paths, up from the original's 228 parameters
after the three build-compat fixes and 45 `$CAPTURE` entries):

- Every `C_<STEM>` (`C_PRP`, `C_ATN`, `C_NAD`, `C_PRM`, `C_PB`, `C_PEM`,
  `C_TOP`, `C_GBP`, `C_OCT`, `C_TTB`) and the three algebraic
  brain-equivalent concentrations (`CB_TOP`, `CB_GBP`, `CB_OCT`) appear in
  `outputPaths` via the extended `$CAPTURE` list (state-derived, so per
  the same reasoning established in `neonatal-hyperbilirubinemia`'s and
  `familial-mediterranean-fever`'s notes, they cannot also be `$PARAM`
  entries).
- Every `EFFECT_<STEM>[_role]` (all seventeen: the six beta-receptor
  quotients, the six GABA terms, the three cortical terms, the two
  olivary terms, the cerebellar term, the Cav3.1 term) appears in
  `outputPaths`.
- Every new `GAMMA_<STEM>[_role]` (fourteen, all = 1.0) and every promoted
  `EC50_SED_<STEM>`/`EC50_ATAX_<STEM>` (six, same values as the original's
  inline literals) appears in `parameters` with the correct default,
  spot-checked directly against the manifest response.
- `GUT_<STEM>`/`CENT_<STEM>`/`PERI_PRP`/`CE_<STEM>` all appear in
  `compartments`, in the same order as the original's `A_*`/`C_*B`
  compartments (only names changed — confirmed by index-based dosing
  against both files producing identical results throughout verification).

No `.cpp` extraction file was left behind — extraction was in-memory only,
used to build the verification requests above and then discarded, per the
workflow guide.

## Anything else flagged

- No compound other than the eleven census rows was touched. Botulinum
  toxin (`A_BTXT`, `A_BTXG`, `SNAPT`, `SNAPG`, `KDEG_BTX`, `KCL_BTX`,
  `KR_SNAP`, `FSPILL`) and benzodiazepine (`EMAX_BZD`, `BZDLEV`, no PK
  model) are byte-identical to the original — neither is one of this
  file's eleven census rows.
- The R-side scenario code (`reg()`/`ev()` calls, `select()`/`mutate()`
  column references, and the one `param(mod, EMAX_TT=..., IC50_TT=...)`
  override in S16) was updated only where it names a renamed compartment,
  parameter, or `$CAPTURE` variable — same dosing amounts, same timing,
  same scenario logic throughout. All 25 of the file's own scenarios
  (`S0`-`S25`) still run against the renamed model.

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`: all
eleven `essential-tremor` rows updated with target/pathway and the
outcome of this refactor; the two mislabeled rows ("Organ systems (PRM)",
"Organ systems (PRP)") corrected to "Primidone (PRM)" and "Propranolol
(PRP)" with a note explaining the original mislabeling (per the same
pattern documented in `neonatal-hyperbilirubinemia/nhb_refactor_notes.md`).
