# Nontuberculous Mycobacterial Lung Disease (MAC-PD) — *Mycobacterium avium* Complex Pulmonary Disease QSP Model

<p align="center">
  <a href="ntm_qsp_model.svg"><img src="ntm_qsp_model.png" width="900" alt="MAC-PD QSP mechanistic map"></a>
</p>

---

## The question this model tries to answer

The standard therapy for MAC lung disease is **a macrolide + ethambutol + rifampin for 12 months or longer**.
Yet this regimen cures 60–70% of patients with the nodular bronchiectatic form and fails to cure even half of
those with the fibrocavitary form. Why does the same drug combination behave so differently in different
patients against the same organism — and why does **intravenous amikacin barely work, while the very same
molecule packaged in a liposome and inhaled (ALIS) does**?

This model was designed not to "assume" the answer but to make it **fall out of the arithmetic**. Doing so
required two structural decisions.

---

## Structural decision 1 — split the bacteria not by phenotype but by **physical niche**

MAC in the lung is not "one population, one drug exposure." It occupies **four physically separate sites**,
and each site differs completely in (a) its growth rate and (b) which drugs can even reach it.

| Niche | Location | Growth | Drugs that can reach it |
|---|---|---|---|
| **B_E** | Airway lumen · ELF, free-floating bacteria | Fast (doubling time ≈ 20 h) | All |
| **B_B** | Biofilm in bronchiectatic segments · mucus plugs | Slow, phenotypic resistance | Only those that cross the EPS diffusion barrier |
| **B_I** | Macrophage phagosome (pH 5.2) | Doubling time 48–72 h | Only those that cross the cell membrane |
| **B_C** | Cavity wall · caseous necrosis | Almost non-replicating | Lipophilic drugs only |

So in this model **every drug has not one but four effective concentrations.** The ELF concentration is
multiplied by a niche-specific penetration coefficient (`PBF*`, `PCS*`), while the intracellular compartment is
computed with its own separate transport equation.

The result is that **the cavitary and nodular bronchiectatic phenotypes are not different models but different
initial conditions of the same model.** `CAVFLAG` changes only a single initial value, that of `B_C`
(3×10⁵ → 1×10⁸). Because essentially no drug reaches the caseum during oral triple therapy (macrolide
penetration 0.15, amikacin 0.02, rifampin 0.30), the gap in culture conversion between the two phenotypes
comes not from a parameter but from **subtraction**.

---

## Structural decision 2 — enter phagosomal pH **once** and let it produce two opposite outcomes simultaneously

`PHPHAG = 5.2`. This single value acts in two directions.

**(i) Accumulation — Henderson-Hasselbalch ion trapping.**
Macrolides are weak bases (pKa ≈ 8.7), so inside the acidic phagosome they become protonated and cannot escape.

```
R_trap = (1 + 10^(pKa − pH_in)) / (1 + 10^(pKa − pH_out))
       = (1 + 10^3.5) / (1 + 10^1.3)
       = 151-fold
```

Azithromycin's famous "100–1000-fold accumulation in pulmonary macrophages" comes out of this.
**In this model that number is not entered as a parameter — it is calculated.**

**(ii) Loss of activity — the same acidic pH raises the MIC.**
Only the un-ionised molecule can bind the 50S peptidyl transferase centre, so

```
MIC_M(pH) = MIC_M0 × 10^(GM × (7.4 − pH_in)) = 12.6-fold
```

**Conclusion: net gain in potency = 151 / 12.6 ≈ 12-fold — not 151-fold.**
**The drug that reaches the phagosome best is exactly the drug that is inactivated most there.**

And because this residual sub-MIC pressure falls on **the compartment the companion drugs cannot reach**, the
rrl A2058G resistance mutation is born nowhere else but **inside the phagosome**. In the model, `R_I` (the
intracellular resistant-bacteria compartment) is not switched on — it **emerges** from this asymmetry.

Amikacin fails at the same pH for **the opposite biophysical reason** — aminoglycoside uptake into the
bacterium is proton-motive-force (PMF) dependent and collapses below pH 6 (`GK` 0.85 > `GM` 0.50, a steeper
penalty). Being a polycation on top of that, **free amikacin cannot even cross the macrophage plasma membrane
at all.**

This is the entire argument for ALIS. Because liposomes are **phagocytosed whole**, they deposit directly in
the lung (bypassing barrier 1) and are delivered all the way into the cell (bypassing barrier 3) — while
routing around the plasma compartment that carries the ototoxicity.

In this model, ALIS's **efficacy** is driven by `KMAC`/`KELF` (lung), while its **ototoxicity** is driven by
`KPERI` (perilymph, which tracks plasma). Because the two pathways are structurally separated, "high in the
lung, low in the ear" is **a result, not an assumption.**

---

## Verification

The R model in this repository was **independently reimplemented and numerically verified**. All 47 ODEs were
separately transcribed into Python/scipy (LSODA) and integrated for 12 scenarios over 540 days, to confirm
that the values in the table below actually arise from the model structure rather than from parameter tuning.

The verification process actually uncovered three structural defects, all of which were fixed — ① the
biofilm/caseum compartments were replicating without limit while unreached by any drug, neutralising every
regimen; ② amikacin's plasma↔ELF partition term, because the ELF volume (25 mL) is small, was returning more
than 3 g/day of inhaled drug back into plasma, destroying ALIS's lung/plasma separation; ③ there was no
resistant biofilm compartment, so resistant organisms had structurally nowhere to establish themselves. Had
this not been verified, all three would have remained.

### Derived values (not parameters)

| Derived quantity | Calculated value | Literature range |
|---|---|---|
| Ion-trapping accumulation ratio `R_trap` | **151×** | Azithromycin macrophage:plasma 10²–10³ (Olsen 1996, Rodvold 1997) |
| Rise in intracellular macrolide MIC | **12.6×** | 4–16-fold at acidic pH (Tulkens 1991, Lemaire 2005) |
| Net intracellular potency gain | **12.0×** | — (this model's prediction) |
| Rise in intracellular amikacin MIC | **74×** | Aminoglycosides are essentially inactive below pH 6 (Maurin 2001) |

### Scenario results (verification run output, sputum log10 CFU/mL)

| # | Scenario | 6 mo | 12 mo | 18 mo<br>(after end of therapy) | 12-mo<br>culture conversion | Hearing loss | 12-mo<br>resistant fraction |
|---|---|---|---|---|---|---|---|
| 1 | Watchful waiting (nodular bronchiectatic) | 5.89 | 5.89 | 5.89 | ✗ | 0 dB | 2.9e-08 |
| 2 | AZM+EMB+RIF 3×/week (nodular) | **−0.33** | **−5.30** | −6.00 | ✓ | 0 dB | 5.4e-15 |
| 3 | AZM+EMB+RIF daily (cavitary) | 3.22 | 0.84 | 5.88 (relapse) | ✓ delayed | 0 dB | 3.4e-09 |
| 4 | + IV amikacin for 3 months | 3.37 | 1.02 | 5.88 (relapse) | ✗ | **15.5 dB** | 3.3e-09 |
| 5 | + Inhaled ALIS 590 mg/day | 2.42 | −0.30 | 5.87 (relapse) | ✓ | **0.5 dB** | 1.1e-10 |
| 6 | Azithromycin monotherapy (nodular) | 4.15 | 3.48 | 5.88 | ✗ | 0 dB | 1.7e-08 |
| 7 | CLR+EMB+RIF (CYP3A interaction) | 3.66 | 2.48 | 5.89 | ✗ | 0 dB | 7.8e-10 |
| 8 | AZM+EMB+CFZ+ALIS (rifampin excluded) | 1.43 | **−3.66** | **−6.00 (sustained)** | ✓ | 0.5 dB | 4.1e-13 |
| 9 | #8 + airway clearance + nutrition | 1.41 | −3.67 | −6.00 (sustained) | ✓ | 0.5 dB | 3.3e-13 |
| 10 | Anti-IFN-γ host + ALIS | 2.42 | −0.30 | 5.99 (relapse) | ✓ | 0.5 dB | 1.1e-10 |
| 11 | Cavitary AZM monotherapy (normal host) | 5.05 | 4.64 | 5.89 | ✗ | 0 dB | 1.8e-08 |
| 12 | Cavitary monotherapy + anti-IFN-γ | 5.95 | 5.94 | 5.96 | ✗ | 0 dB | **1.00 (complete sweep)** |

### Three things to read out of this table

**① Scenario 2 vs. 3 — the gap comes not from a parameter but from a single initial value.**
The only difference between the two scenarios is the initial value of `B_C` (the caseum compartment): 3×10⁵
versus 1×10⁸. Not a single other parameter differs. The nodular bronchiectatic form converts by 6 months,
while the cavitary form barely reaches conversion at 12 months and then relapses once therapy stops.

**② Scenario 4 vs. 5 — the same molecule, a different address.**
IV amikacin gives no benefit (12-month value 1.02 versus 0.84 — actually slightly worse) while taking away
**15.5 dB of hearing**, because the intracellular concentration is zero. Inhaling the same amikacin as a
liposome gives an intracellular concentration of 897 mg/L with only **0.5 dB** of hearing loss — a 31-fold
separation. This is an arithmetic consequence of `KMAC` (lung) and `KPERI` (perilymph, which tracks plasma)
sitting on different pathways.

**③ Scenario 11 vs. 12 — selection pressure alone does not produce resistance.**
Both scenarios are macrolide monotherapy, and in both, Φ → 1 in every niche. Yet when the host is normal, the
resistant fraction stays at 10⁻⁸; once anti-IFN-γ autoantibodies collapse the Th1 axis, it **sweeps all the
way to 100%.** In other words, emergence of resistance requires **two conditions at once** — ① selection
pressure, and ② a host state in which the net growth rate of the resistant organisms turns positive (high
bacterial burden plus impaired clearance). This agrees with the clinical observation (Griffith 2006, Morimoto
2016) that acquired macrolide resistance is not "a property of the drug" but the product of high bacterial
burden, functional monotherapy, and impaired clearance occurring together.

> **On the scope of the verification.** The numbers above come from an **independent reimplementation** of
> all 47 ODEs in Python/scipy (LSODA). This verifies the model's *equations and dynamics*; whether the
> mrgsolve (R) code carried in this repository itself compiles and runs could not be confirmed, because this
> environment has no R. The two implementations were mapped 1:1 to share the same equations and the same
> parameters.

---

## File Layout

| File | Content |
|---|---|
| [`ntm_qsp_model.dot`](ntm_qsp_model.dot) | Mechanistic map source — **194 nodes, 15 subgraph clusters** |
| [`ntm_qsp_model.svg`](ntm_qsp_model.svg) | Vector map (zoom in to view) |
| [`ntm_qsp_model.png`](ntm_qsp_model.png) | Raster map (150 dpi) |
| [`ntm_mrgsolve_model.R`](ntm_mrgsolve_model.R) | **47-ODE** mrgsolve model + 12 treatment scenarios + sensitivity analysis |
| [`ntm_shiny_app_en.R`](ntm_shiny_app_en.R) | **11-tab** interactive dashboard |
| [`ntm_references_en.md`](ntm_references_en.md) | **76 references**, 12 sections |

### Map clusters (15)

Host susceptibility phenotype · exposure and airway defence · **four bacterial niches** · resistance evolution
· innate immunity and the phagosome · adaptive immunity (IL-12/IFN-γ) · tissue destruction (Cole vicious
cycle · cavitation) · macrolide PK/PD · ethambutol · rifamycin and CYP3A interaction · **IV amikacin versus
ALIS** · salvage therapy · clinical endpoints · diagnosis and monitoring · legend

### ODE compartments (47)

- **Drug PK (22)** — macrolide 5 (gut · plasma · peripheral · ELF · intracellular), ethambutol 3, rifampin 3 +
  1 CYP3A enzyme compartment, amikacin 7 (plasma · peripheral · ELF · liposomal lung depot · intracellular ·
  perilymph · renal cortex), clofazimine 3
- **Bacteria (7)** — `B_E` `B_B` `B_I` `B_C` + resistant `R_E` `R_B` `R_I`
  (the initial value of each resistant compartment = mutation rate × the burden of the corresponding
  compartment → if the burden is low, the mutant simply does not exist to begin with)
- **Host (11)** — macrophages · IFN-γ · TNF-α · neutrophils · MMP · mucus · mucociliary clearance ·
  bronchiectasis · cavity · symptoms · body weight
- **Toxicity and tracking (7)** — hearing · optic nerve · renal tubule · liver · QTc effect compartment ·
  cumulative days of conversion · cumulative exposure

---

## Treatment Scenarios (12)

| # | Scenario | What it shows |
|---|---|---|
| 1 | Watchful waiting | Untreated natural history — the option guidelines allow |
| 2 | AZM+EMB+RIF 3×/week, nodular bronchiectatic | The case where standard therapy works well |
| 3 | AZM+EMB+RIF daily, cavitary | Same drugs, only `B_C`'s initial value differs → the arithmetic origin of the gap |
| 4 | + IV amikacin for 3 months | Plasma and perilymph rise while the intracellular compartment stays at **zero** |
| 5 | + Inhaled ALIS 590 mg/day | The CONVERT trial design — the same molecule, a different address |
| 6 | Azithromycin **monotherapy** (nodular) | Φ → 1, yet resistance cannot establish itself if the host is normal |
| 7 | CLR+EMB+RIF | CYP3A induction means **the third drug eats the first** |
| 8 | AZM+EMB+CFZ+ALIS | Rifampin excluded · replaced with a lipophilic drug that can reach the caseum |
| 9 | #8 + airway clearance therapy + nutritional support | Non-drug levers (the MCC and TNF-depletion loop) |
| 10 | Host with anti-IFN-γ autoantibodies + ALIS | Why drugs alone cannot rescue this patient |
| 11 | Cavitary AZM monotherapy, normal host | Selection pressure alone does not produce resistance |
| 12 | #11 + anti-IFN-γ autoantibodies | Resistance only sweeps once the two conditions coincide |

---

## How to Run

```r
# Model + 12 scenarios + sensitivity analysis
source("ntm_mrgsolve_model.R")

# Interactive dashboard (11 tabs)
shiny::runApp("ntm_shiny_app_en.R")
```

```bash
# Re-render the maps
dot -Tsvg ntm_qsp_model.dot -o ntm_qsp_model.svg
dot -Tpng -Gdpi=150 -Gsize="30,18" ntm_qsp_model.dot -o ntm_qsp_model.png
```

Required packages: `mrgsolve`, `dplyr`, `tidyr`, `ggplot2`, `shiny`, `DT`

---

## Predictions This Model Makes (in falsifiable form)

1. **Drop ethambutol and the regimen quietly becomes macrolide monotherapy.** In the model, EMB's main
   contribution is not direct killing but enhancing cell-wall permeability (the `PERM` term), so removing EMB
   drops the effective concentrations of both the macrolide and amikacin **simultaneously**, raising the
   resistance-selection gate Φ toward 1.
2. **Combining rifampin with clarithromycin cuts the anchor drug's exposure by more than 60%.**
   Azithromycin escapes this trap because only 0.05 of its clearance is CYP3A-dependent — in the model this
   comes out arithmetically from a single difference in the `FCYP3A` parameter.
3. **IV amikacin's intracellular concentration is zero.** IV amikacin therefore pays the full ototoxicity cost
   while contributing nothing at all to one of the three niches (`B_I`).
4. **The caseum compartment is the reservoir of relapse.** The curve that climbs back up after the end of
   therapy always starts from `B_C`. This is why the guideline criterion of "12 months of sustained culture
   conversion" exists.
5. **Resistance does not arise from selection pressure alone.** Macrolide monotherapy drives Φ → 1 in every
   niche, but as long as the IFN-γ axis and mucociliary clearance remain intact, mutants are cleared faster
   than they replicate and stay at 10⁻⁸. For resistance to sweep, **selection pressure and impaired clearance
   must coincide** (scenario 11 vs. 12). This is a conclusion the model produced on its own, one that was not
   anticipated before verification.

---

## Limitations

- This is a **deterministic single-patient model**. It has no inter-individual variability and no stochastic
  extinction, so culture conversion is represented as crossing a threshold rather than as a probability. It
  cannot be compared directly against a clinical trial's conversion percentage — only the direction and order
  of magnitude should be compared.
- **It does not distinguish relapse from reinfection** (a distinction that, in practice, requires genotyping).
- **Parameters are literature-based approximations** and have not been fitted to patient data.
- Rifampin's contribution coming out near zero reflects its short half-life, which makes the time-averaged
  Hill term based on instantaneous concentration small; this is consistent with MAC's inherently low rpoB
  affinity, but it is also **a consequence of a modelling choice**. It should not be read as a judgement on
  rifampin's clinical value.

> ⚠️ **This is a QSP model for education and research purposes. It must not be used for clinical
> decision-making, prescribing, or regulatory submission.**

---

*Claude Code Routine · QSP Disease Model Library · 2026-08-01*
