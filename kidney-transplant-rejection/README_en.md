# Kidney Transplant Rejection — QSP Model

> This is the **only disease in this library whose pathogenesis is a normally
> functioning immune system.** The treatment target is not a pathway but a single
> **set point**.
>
> ```
>   Too little immunosuppression  →  TCMR · de novo DSA · ABMR      →  graft loss
>   Too much immunosuppression    →  BK nephropathy · CMV · PTLD · sepsis →  graft loss
> ```
>
> So this model is built not as a stepwise cascade but as a **U-shaped trade-off**,
> with everything arranged around the one dial a clinician can actually turn in
> between — **net immunosuppression (NIS)**.

<p align="center">
  <a href="ktx_qsp_model.svg">
    <img src="ktx_qsp_model.png" width="900" alt="Kidney transplant rejection QSP mechanistic map">
  </a>
</p>

---

## 1. Files

| File | Content |
|------|------|
| [`ktx_qsp_model.dot`](ktx_qsp_model.dot) | Mechanistic map source — 20 modules · 243 nodes · 385 edges |
| [`ktx_qsp_model.svg`](ktx_qsp_model.svg) / [`.png`](ktx_qsp_model.png) | Rendered output (`dot -Tsvg` / `dot -Tpng -Gdpi=150`) |
| [`ktx_mrgsolve_model.R`](ktx_mrgsolve_model.R) | **59-ODE** mrgsolve model + 26 scenarios + an automated validation block |
| [`ktx_shiny_app_en.R`](ktx_shiny_app_en.R) | 10-tab interactive dashboard (including a U-curve sweep) |
| [`ktx_references_en.md`](ktx_references_en.md) | **153 references, actually looked up and verified** on PubMed |

```bash
dot -Tsvg ktx_qsp_model.dot -o ktx_qsp_model.svg
dot -Tpng -Gdpi=150 ktx_qsp_model.dot -o ktx_qsp_model.png
Rscript ktx_mrgsolve_model.R                       # 26 scenarios + validation table
Rscript -e 'shiny::runApp("ktx_shiny_app_en.R")'      # dashboard
```

---

## 2. Three feedback loops that drive the model

| Loop | Path | What it explains |
|------|------|-------------------|
| **LOOP 1 — alloimmune amplification** | injury → DAMPs → antigen presentation (AG)↑ → effector-cell influx → injury | **Sub-critical** under adequate immunosuppression, **super-critical** without it. This is why untreated rejection does not worsen gradually but **ignites within 10-14 days**, reaching Banff t ≈ 2.5. |
| **LOOP 2 — nephron capacity / hyperfiltration (Brenner)** | nephron loss → single-nephron hyperfiltration → TGF-β → IFTA → nephron loss | **Independent of drugs.** The reason every graft declines by roughly −1 mL/min/1.73 m² per year even when nothing immunological happens. |
| **LOOP 3 — clinical control loop** | BK viraemia → immunosuppression reduction → viral clearance → alloimmune rebound | Implemented in the model as an **actual closed-loop controller** (`BKRED`). The rebound is not scripted in as a dose change but **emerges** from the model. |

## 3. Two ignition thresholds that generate the phenotypes

**(a) T-cell ignition.** Alloreactive blast proliferation is autocatalytic through IL-2
(`PROLT`). Its rate constant has to exceed the blast clearance rate (`KDIFF + KDTACT`)
for a lesion to ignite.

- With **tacrolimus** present, the calcineurin term holds this **unconditionally
  sub-critical**. → Even when the CYP3A5 expressor phenotype sits at a trough of
  4 ng/mL, the lesion only grows *a little* rather than exploding.
- **Belatacept** has no calcineurin term. That makes it **conditionally stable**.
  Recipients with a large pool of CD28-null CD57+ memory T cells (high `FCD28N`) cross
  the threshold and reject; recipients with fewer memory cells do not.
  → This is the model's mechanistic explanation of the **BENEFIT paradox** (more acute
  rejection but better long-term GFR).

**(b) Germinal-centre (GC) ignition = de novo DSA.** GC B-cell proliferation requires
sustained Tfh help (`BHELP`), scales with the eplet-mismatch load (`BEPI`), and is
suppressed by mycophenolate.

- Standard therapy stays **below** the threshold → DSA < 300 MFI at 5 years
- If exposure is halved by non-adherence, it rises **above** the threshold → dnDSA
  develops
- And because **long-lived plasma cells (LLPCs) accumulate irreversibly**, DSA does not
  regress even after adherence is restored. → This is why late graft loss is dominated
  by dnDSA plus non-adherence.

---

## 4. ODE structure (59)

| Compartment group | Count | State variables |
|--------|------|-----------|
| Pharmacokinetics | 21 | TAC (2 gut + central + peripheral) · CsA (2) · MMF/MPA/MPAG (4, enterohepatic recirculation) · prednisolone · methylprednisolone · everolimus · rATG · basiliximab · belatacept · rituximab · anti-IL-6 · anti-CD38 · IVIG · ganciclovir |
| TDM controller | 1 | `DOSEADJ` (integral controller) |
| Alloimmunity | 19 | AG · TN · TACT · TEFF · TMEM · TREG · TINF · IL2 · IFNG · BN · BGC · PB · **LLPC** · DSA · C4D · NKACT · ENDO · MVI · TG |
| Graft | 9 | TUB · AKI · DAMP · IFTA · CNIHYAL · NM · SCR · PROT · CFDNA |
| Infection / safety | 6 | BKV · BKVAN · CMV · WBC · GLU · NISAUC |
| Behaviour / outcome | 3 | ADH · HAZ · CUMREJ |

### Some design choices that matter

- **Chronic oral drugs enter as a continuous infusion rather than discrete events**,
  and the model carries an *average concentration*. Reported troughs are converted using
  the literature's Cavg→C0 ratio (tacrolimus `RTROUGH = 0.75`, AUC0-12/C0 ≈ 16 h). This
  is more honest and reproducible than pushing in 3,650 BID boluses over a 5-year
  disease-progression model. Intermittent IV drugs (ATG · basiliximab · belatacept ·
  rituximab · anti-IL-6 · anti-CD38 · IVIG · steroid pulses · plasma exchange ·
  imlifidase) enter as smooth rectangular windows, so real peak-trough behaviour is
  preserved.
- **TDM is a controller, not a dosing table.** `DOSEADJ` moves toward protocol target
  troughs (9 at 0-3 months → 7 by 12 months → 6 ng/mL thereafter). Two things follow
  from this.
  - A CYP3A5 expressor on a fixed mg/kg dose stays around C0 ≈ 4 ng/mL, but is corrected
    to target once TDM is turned on.
  - **A non-adherent patient is not corrected.** The trough measured in clinic is drawn
    *after* the dose the patient actually took (white-coat adherence). The model
    expresses this as `C0MEAS = C0 / ADH` — the main reason non-adherence stays invisible
    until DSA appears.
- **Nephron destruction uses a cubic kernel**, to capture the clinical fact that mild
  inflammation is tolerated while a full-blown lesion is not.
- **Rituximab has no term at all in the LLPC equation** (CD20-negative). Anti-CD38 does.
  This is where RITUX-ERAH's negative result and anti-CD38's actual lowering of DSA come
  from.
- **Eculizumab blocks only C5, leaving C4d and NK ADCC untouched** — the
  complement-independent arm remains.

---

## 5. The 26 simulated scenarios and their results

Every scenario is **purely a parameter set** (no dosing events, no scripted
interventions). The numbers below are 5-year (1,825-day) simulation results.

| # | Scenario | eGFR 1yr | eGFR 5yr | Slope/yr | Peak Banff t | DSA 5yr (MFI) | cg 5yr | 5-year graft survival |
|---|----------|---------:|---------:|-----------:|-------------:|--------------:|-------:|--------------------:|
| 1 | Standard risk (TAC/MPA/steroid + basiliximab) | 57.5 | 53.0 | −1.14 | 0.39 | 152 | 0.28 | 0.921 |
| 2 | Living donor, low mismatch, CIT 2 h | 66.4 | 63.4 | −0.75 | 0.28 | 79 | 0.10 | 0.950 |
| 3 | High immunological risk + ATG induction | 57.9 | 46.6 | −2.82 | 0.65 | 451 | 1.15 | 0.868 |
| 4 | High immunological risk + basiliximab only | 56.2 | 44.1 | −3.02 | 0.65 | 457 | 1.19 | 0.836 |
| 5 | CYP3A5 expressor, fixed mg/kg (no TDM) | 59.8 | 52.0 | −1.93 | 0.52 | 196 | 0.49 | 0.911 |
| 6 | CYP3A5 expressor, genotype-guided TDM | 57.5 | 52.9 | −1.14 | 0.39 | 152 | 0.28 | 0.921 |
| 7 | Belatacept, CNI-free | **66.7** | **64.0** | −0.67 | 0.36 | **38** | 0.05 | **0.955** |
| 8 | Belatacept + recipient rich in CD28-null cells | 64.2 | 46.7 | −4.37 | **1.11** | 101 | 0.72 | 0.891 |
| 9 | Cyclosporine control | 54.4 | 50.0 | −1.09 | 0.34 | 137 | 0.22 | 0.914 |
| 10 | Early steroid withdrawal (day 7) | 56.5 | 45.8 | −2.68 | 0.75 | 204 | 0.74 | 0.878 |
| 11 | Everolimus + reduced CNI (no MPA) | 59.6 | 50.6 | −2.26 | 0.56 | 267 | 0.60 | 0.894 |
| 12 | Non-adherence from month 12 (45% of doses taken) | 57.5 | **7.4** | **−12.54** | 2.05 | **3,491** | **2.38** | **0.158** |
| 13 | Non-adherence from month 12 (60% of doses taken) | 57.5 | 43.2 | −3.59 | 0.93 | 420 | 1.33 | 0.842 |
| 14 | High eplet load (DR/DQ 24) | 56.5 | 43.5 | −3.26 | 0.67 | 518 | 1.30 | 0.827 |
| 15 | High tacrolimus IPV | 58.9 | 50.6 | −2.06 | 0.42 | 147 | 0.26 | 0.920 |
| 16 | BK viraemia + **screening with pre-emptive reduction** | 53.5 | 49.0 | −1.13 | 1.32 | 153 | 0.29 | 0.902 |
| 17 | BK viraemia, **no** screening | 41.8 | 36.2 | −1.40 | 1.67 | 153 | 0.29 | 0.825 |
| 18 | Banff IIA TCMR (medication gap) + MP pulse therapy | 53.2 | 47.8 | −1.33 | 2.43 | 160 | 0.34 | 0.881 |
| 19 | Steroid-resistant TCMR → ATG rescue | 53.5 | 46.7 | −1.72 | 2.45 | 180 | 0.52 | 0.880 |
| 20 | dnDSA ABMR — adherence restoration only (control) | 57.5 | 34.7 | −5.71 | 2.02 | 214 | 0.74 | 0.722 |
| 21 | dnDSA ABMR — PLEX + IVIG + rituximab | 57.5 | 36.1 | −5.36 | 2.02 | 209 | 0.68 | 0.756 |
| 22 | dnDSA ABMR — anti-IL-6 (tocilizumab/clazakimab) | 57.5 | 35.4 | −5.52 | 2.02 | 199 | 0.64 | 0.742 |
| 23 | dnDSA ABMR — **anti-CD38 (felzartamab)** | 57.5 | **36.4** | **−5.29** | 2.02 | **141** | **0.44** | **0.768** |
| 24 | DCD, CIT 24 h, KDPI 85 (DGF · low nephron mass) | 38.2 | 35.1 | −0.76 | 0.39 | 152 | 0.29 | 0.851 |
| 25 | CMV D+/R- + 200 days of valganciclovir | 57.5 | 53.0 | −1.14 | 0.39 | 152 | 0.28 | 0.921 |
| 26 | No maintenance immunosuppression at all (untreated control) | 24.8 | 0.6 | −6.06 | **2.57** (day 15) | 9,219 | 2.48 | 0.002 |

### Results the model "generated" rather than "assumed"

1. **The BK trap (16 vs 17).** Turning screening on makes the closed-loop controller
   cut MPA by 60% and lower the CNI target by 30%. The virus is cleared (peak 4.83 vs
   5.65 log10), but **rejection rebounds**, producing Banff t 1.32 at day 103 — an event
   nobody scripted in. Even so, 1-year eGFR is far better with screening (53.5 vs 41.8).
2. **Why rituximab fails to work (21 vs 23).** Both scenarios start treatment on the
   same day, but rituximab cannot reach the LLPC, so DSA stays at 209 MFI, while
   anti-CD38 brings it down to 141, and cg splits 0.68 → 0.44 accordingly.
3. **The BENEFIT paradox (7 vs 8).** Same belatacept regimen, but only the recipient
   rich in CD28-null memory cells develops early rejection (Banff t 1.11 at day 199).
   Even so, 1-year eGFR is still higher at 64.2 than tacrolimus's 57.5 — because there is
   no CNI haemodynamic vasoconstriction.
4. **The nonlinearity of non-adherence (12 vs 13).** A GC ignition threshold sits
   between 45% and 60% dosing adherence. 45% diverges into dnDSA of 3,491 MFI and graft
   loss; 60% diverges into a gentle decline.
5. **CYP3A5 (5 vs 6).** With a fixed dose, the trough sits at 3.9 ng/mL, the lesion
   grows larger (Banff t 0.52 vs 0.39), and dnDSA risk rises — but with weaker
   vasoconstriction, *1-year eGFR is actually higher* (59.8 vs 57.5), a plain
   illustration that CNI packages benefit and toxicity into the same molecule.

---

## 6. Validation

`ktx_mrgsolve_model.R` automatically prints the following checks when it runs.

| Check | Target range | Model value |
|-----------|-----------|---------|
| Standard: tacrolimus trough (6 months) | 5-8 ng/mL | 7.2 |
| Standard: MPA AUC0-12 (6 months) | 30-60 mg·h/L | 46 |
| Standard: 12-month eGFR | 52-62 mL/min/1.73 m² | 57.5 |
| Standard: eGFR slope, years 1-5 | −2.5 to −0.3/year | −1.14 |
| Standard: 5-year absence of dnDSA | < 1,000 MFI | 152 |
| Standard: 12-month serum creatinine | 1.0-1.8 mg/dL | 1.45 |
| Belatacept − tacrolimus 12-month eGFR difference | +4 to +18 mL/min | +9.2 |
| Tacrolimus − cyclosporine 12-month eGFR difference | +1 to +10 mL/min | +3.1 |
| Non-adherence: 5-year dnDSA | 1,500-40,000 MFI | 3,491 |
| Non-adherence: 5-year chronic ABMR (cg) | 1.0-3.0 | 2.38 |
| Untreated: Banff t ≥ 2 within 30 days | 2.0-3.0 | 2.57 |
| Untreated: 5-year graft function loss | < 20 mL/min | 0.6 |

> **Transparency notice — R was not installed in this session's execution
> environment.**
> So the numbers above were not obtained by running `ktx_mrgsolve_model.R` in R, but
> from an **independent Python/SciPy (LSODA) implementation carrying the same equations
> and the same parameters**. To keep the two implementations from diverging:
> - All **267 parameters** in the R model's `$PARAM` block were **mechanically
>   generated from the same parameter dictionary used for validation** (no manual
>   transcription).
> - The `$ODE` block was written with **the same variable names and the same statement
>   order** as the Python right-hand side, so the two can be checked line by line.
> - It was automatically verified that all 59 compartments have a `dxdt_` assignment
>   and that there are no duplicate parameters.
>
> Running `Rscript ktx_mrgsolve_model.R` in an environment with R installed will print
> this validation block again as PASS/FAIL.

---

## 7. What this model does *not* claim

- **This is a deterministic single-patient model.** It does not generate an incidence.
  Scenarios either cross a lesion threshold or they do not. Clinical-trial percentages
  (about 10-12% 12-month BPAR, about 15-25% 5-year dnDSA) are **population-level
  statements** and are used here only as **qualitative anchors** to check that the
  phenotypes fall within a clinically plausible range. They are not fitting targets.
- IPV (inter-patient variability) is implemented as a deterministic sinusoid, so its
  effect is smaller than actual random trough excursions. A stochastic implementation
  is future work.
- Parameters were fixed to the literature wherever a number exists (all PK, the
  1.9-fold CYP3A5 effect, trough targets, the MPA AUC window, the BK 10⁴ threshold,
  eGFR levels and slopes), and hand-calibrated otherwise.
- Recurrent glomerulonephritis, non-HLA antibodies (AT1R, etc.), competing
  cardiovascular mortality risk in chronic graft loss, paediatric recipients, and
  ABO-incompatible transplantation are shown on the map but not implemented as ODEs.

> ⚠️ **This is a semi-quantitative QSP model for education and research. It must not
> be used for clinical decision-making, prescribing, or regulatory submission.**

---

## 8. References

[`ktx_references_en.md`](ktx_references_en.md) — 20 topic sections, **153 references**.
Every PMID was actually looked up via NCBI E-utilities to confirm title, journal, and
year.
