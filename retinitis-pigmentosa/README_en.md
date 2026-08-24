# Retinitis Pigmentosa (RP) — QSP Model

> Integrated Quantitative Systems Pharmacology model linking RHO/RPGR/USH2A/
> PDE6/RPE65 genotype-dependent primary rod photoreceptor apoptosis
> (rhodopsin misfolding/ER stress, PDE6-cGMP-Ca2+ excitotoxicity, RPGR-
> ciliopathy transport failure) to secondary cone death (RdCVF loss,
> oxidative stress) and microglial/gliotic amplification — producing
> progressive night blindness, visual field constriction, and central
> vision loss — coupled to voretigene neparvovec (AAV2-RPE65 subretinal
> gene therapy), investigational RPGR gene augmentation (AAV8/AAV5-RPGR),
> MCO-010 optogenetic gene therapy (AAV2 multi-characteristic opsin,
> intravitreal), CNTF encapsulated-cell neuroprotection, N-acetylcysteine
> (antioxidant), and vitamin A palmitate/DHA supplementation PK/PD.

| Deliverable | File |
|---|---|
| 🗺️ Mechanistic map (DOT) | [`rp_qsp_model_en.dot`](rp_qsp_model_en.dot) |
| 🖼️ Map (SVG)             | [`rp_qsp_model.svg`](rp_qsp_model.svg) |
| 🖼️ Map (PNG, 150 dpi)    | [`rp_qsp_model.png`](rp_qsp_model.png) |
| ⚙️ mrgsolve ODE model     | [`rp_mrgsolve_model.R`](rp_mrgsolve_model.R) |
| 📊 Shiny dashboard        | [`rp_shiny_app_en.R`](rp_shiny_app_en.R) |
| 📚 References             | [`rp_references_en.md`](rp_references_en.md) |

---

## 1. Disease in One Paragraph

Retinitis pigmentosa is the most common inherited retinal dystrophy, caused by mutations in over 90
genes, and is a progressive disease in which primary death of rod photoreceptors leads to secondary
death of cone photoreceptors and ultimately blindness. The major causes include autosomal dominant
**RHO** (rhodopsin misfolding → ER stress/UPR), X-linked **RPGR** (ciliary transport defect),
autosomal recessive **USH2A** (may be accompanied by Usher syndrome), **PDE6A/B** (failed cGMP
hydrolysis → CNG channel hyperactivation → Ca2+ toxicity), and **RPE65** (deficiency of the
visual-cycle isomerase, seen in Leber congenital amaurosis/early-onset severe retinal dystrophy). Once
rod photoreceptors are lost, retinal oxygen consumption falls, exposing cones to a relatively
hyperoxic environment, and depletion of rod-derived cone viability factor (RdCVF/NXNL1) triggers
oxidative-stress-mediated secondary cone death. Microglial activation and Müller cell reactive gliosis
amplify this process, and the clinical course progresses from night blindness (the first symptom) →
progressive peripheral visual field loss (ring scotoma → tunnel vision) → eventual central vision
loss. For patients with biallelic RPE65 mutations, **voretigene neparvovec** (AAV2-RPE65 subretinal
gene therapy, FDA-approved based on the 2017 phase 3 Russell Lancet study) restores the visual cycle,
and AAV8/AAV5-RPGR gene therapy is in clinical development for X-linked RP (Cehajic-Kapetanovic 2020
Nat Med, XIRIUS phase 3). **MCO-010 optogenetic therapy** (a multi-characteristic opsin expressed in
surviving bipolar/ganglion cells), which works independently of photoreceptor survival, can be applied
to patients regardless of disease stage or genotype; CNTF sustained-release implants (neuroprotective,
with a paradoxical accompanying ERG suppression), N-acetylcysteine (antioxidant), and vitamin A
palmitate + DHA (modest slowing of progression, Berson 1993/2004) are used as
adjunctive/neuroprotective therapies. Patients with end-stage disease are candidates for the Argus II
retinal prosthesis (an electrical-stimulation bypass pathway, commercially discontinued in 2020) or
low-vision aids.

## 2. Mechanistic Map Clusters (12 Clusters, 136 Nodes)

1. Genetic etiology (RHO dominant · RPGR/RP2 X-linked · USH2A/PDE6A/B/CRB1/NR2E3/RP1/EYS
   recessive · PRPF splicing-factor dominant · BBS ciliopathy · Usher types I/II/III)
2. Molecular pathophysiology: photoreceptor protein abnormalities (rhodopsin misfolding/UPR, PDE6
   loss of function → cGMP accumulation → CNG channel hyperactivation → Ca2+ toxicity →
   calpain/caspase, RPGR-ORF15/BBSome ciliary transport defects, PRPF splicing abnormalities,
   PARP1/AIF cell death)
3. Rod-cone interaction and secondary death (loss of RdCVF/NXNL1, decreased retinal oxygen
   consumption → hyperoxic exposure → oxidative stress → cone death, disrupted glucose/lactate
   shuttle, mTOR/autophagy abnormalities, relative preservation of foveal cones)
4. Retinal pigment epithelium/visual cycle (RPE65 isomerase, LRAT, 11-cis-retinal regeneration,
   apo-opsin accumulation/constitutive signalling, RPE phagocytosis, pigment migration/bone-spicule
   pigmentation, vascular attenuation, optic atrophy)
5. Neuroinflammation/retinal remodelling (Müller cell reactive gliosis, microglial activation,
   complement activation, bipolar/horizontal cell dendrite retraction, aberrant neurite sprouting,
   relative preservation of ganglion cells)
6. Clinical phenotype and staging (night blindness → ring scotoma → tunnel vision → central vision
   loss, Usher-associated hearing loss, Bardet-Biedl obesity/polydactyly/renal disease, posterior
   subcapsular cataract, genotype-specific progression rate)
7. Clinical endpoints and biomarkers (ERG rod/cone amplitude, Goldmann visual field, BCVA, OCT EZ
   width, fundus autofluorescence ring, FST, MLMT, CST, NEI-VFQ-25)
8. Gene therapy PK/PD (voretigene neparvovec subretinal AAV2-RPE65, botaretigene
   sparoparvovec/cotoretigene toliparvovec AAV8/AAV5-RPGR, CRISPR editing)
9. Optogenetic therapy PK/PD (MCO-010 intravitreal AAV2 multi-characteristic opsin, a bypass
   pathway expressed in surviving bipolar/ganglion cells)
10. Neuroprotective/antioxidant therapy PK/PD (CNTF encapsulated-cell implant, N-acetylcysteine,
    vitamin A palmitate + DHA, investigational recombinant RdCVF, historical valproic acid)
11. Complication management (cystoid macular edema — dorzolamide/acetazolamide, posterior
    subcapsular cataract surgery, low-vision rehabilitation)
12. Retinal prosthesis/assistive devices (Argus II epiretinal electrode array, low-vision aids)

## 3. mrgsolve Model (23 ODE Compartments)

* **Drug/vector PK (10 compartments)** — voretigene neparvovec (subretinal vector genome/RPE65
  expression, 2 compartments), investigational RPGR gene therapy (vector/expression, 2 compartments),
  MCO-010 (vector/opsin expression, 2 compartments), CNTF implant (tissue concentration, 1
  compartment), N-acetylcysteine (oral depot/plasma, 2 compartments), vitamin A/DHA depot (1
  compartment).
* **Disease/PD (6 compartments)** — rod survival fraction (ROD_FRAC), cone survival fraction
  (CONE_FRAC), oxidative-stress index (ROS), microglial activation index (MICROGLIA), retinal
  ganglion cell survival fraction (RGC_FRAC).
* **Clinical endpoints (7 compartments, delayed-relaxation model)** — ERG rod/cone amplitude, visual
  field area, BCVA, FST, MLMT, cystoid macular edema central subfield thickness (CME_CST).
* In the RPE65 genotype, lower effective visual-cycle flux increases apo-opsin toxicity, worsening
  the rod death rate, and voretigene neparvovec expression reverses this via an Emax model. In the
  XLRP genotype, RPGR gene-therapy expression competitively corrects the contribution of the ciliary
  transport defect. MCO-010 works independently of photoreceptor survival, but provides a bypass
  signal proportional to the ganglion-cell survival fraction. CNTF was modelled to mitigate rod/cone
  death while also causing a reversible suppression of ERG amplitude (the Birch 2013 paradox).

### 10 Scenarios

| # | Scenario | Calibration Basis |
|---|---|---|
| 1 | Natural history — RHO-adRP (SEVERITY=1.0) | Berson 1985 Am J Ophthalmol |
| 2 | Natural history — RPGR-XLRP (SEVERITY=1.6) | Hartong 2006 Lancet |
| 3 | Natural history — RPE65-LCA/EOSRD (SEVERITY=2.5) | Cideciyan 2013 PNAS |
| 4 | Voretigene neparvovec (biallelic RPE65) | Russell 2017 Lancet, Maguire 2019 Ophthalmology |
| 5 | Investigational RPGR gene therapy (XLRP) | Cehajic-Kapetanovic 2020 Nat Med, XIRIUS 2024 |
| 6 | MCO-010 optogenetic therapy (end-stage, genotype-independent) | Busskamp 2010 Science, Sahel 2021 Nat Med |
| 7 | CNTF encapsulated-cell implant | Sieving 2006 PNAS, Birch 2013 Am J Ophthalmol |
| 8 | Oral N-acetylcysteine | Campochiaro 2020 J Clin Invest |
| 9 | Oral vitamin A palmitate + DHA | Berson 1993/2004 Arch Ophthalmol |
| 10 | Voretigene neparvovec + NAC combination | Combined neuroprotection hypothesis (extrapolation) |

## 4. Shiny Dashboard (8 Tabs)

1. **Patient profile** — genotype (severity) selection, age, simulation duration.
2. **PK** — gene-therapy expression (RPE65/RPGR/MCO opsin), CNTF tissue concentration, NAC plasma
   concentration, vitamin A depot.
3. **Key PD metrics** — rod/cone survival fraction, oxidative stress/microglial activation.
4. **Clinical endpoints** — ERG, visual field area/BCVA, FST/MLMT.
5. **Scenario comparison** — overlaid comparison of multiple scenarios plus a summary table.
6. **Biomarkers** — ganglion cell survival fraction, cystoid macular edema central subfield
   thickness (including a CAI co-administration option).
7. **Gene therapy · optogenetics** — visual-cycle restoration index (RPE65), optogenetic bypass
   signal (MCO-010).
8. **References** — the full reference list.

## 5. Usage

```bash
# 1) Render the mechanistic map
dot -Tsvg rp_qsp_model_en.dot -o rp_qsp_model.svg
dot -Tpng -Gdpi=150 rp_qsp_model_en.dot -o rp_qsp_model.png
```

```r
# 2) R/mrgsolve simulation
install.packages(c("mrgsolve","dplyr","tidyr","ggplot2","shiny","DT"))
library(mrgsolve)
mod <- mread("rp_mrgsolve_model.R") %>% param(SEVERITY = 2.5, IS_RPE65 = 1)
e_gt <- ev(amt = 100, cmt = "GT65_VG", time = 0)  # voretigene neparvovec, single subretinal dose
out <- mod %>% ev(e_gt) %>% mrgsim(end = 24*365*10, delta = 24)  # 10-year follow-up
plot(out, c("ROD_FRAC", "CONE_FRAC", "ERG_ROD", "MLMT"))

# 3) Launch the Shiny dashboard
shiny::runApp("rp_shiny_app_en.R")
```

## 6. Key Clinical Calibration Evidence

| Endpoint | Comparator | Evidence |
|---|---|---|
| Natural-history ERG amplitude decline rate (~15-20%/year) | 3-year natural-history follow-up | Berson 1985 Am J Ophthalmol |
| Voretigene neparvovec MLMT/FST improvement | 1-year phase 3, 4-year durability | Russell 2017 Lancet, Maguire 2019 Ophthalmology |
| RPGR gene therapy, dose-dependent retinal atrophy | Phase 1/2, high-dose group | Cehajic-Kapetanovic 2020 Nat Med |
| Optogenetic partial restoration of visual function | Single-patient case report | Sahel 2021 Nat Med |
| CNTF paradoxical ERG suppression | Randomised phase 3 | Birch 2013 Am J Ophthalmol |
| Oral vitamin A supplementation, attenuated ERG decline (~20%) | DBA randomised trial | Berson 1993 Arch Ophthalmol |
| N-acetylcysteine, improved cone function | Phase 1 dose escalation | Campochiaro 2020 J Clin Invest |
| CAI (dorzolamide/acetazolamide) macular edema response rate | Open-label/retrospective study | Fishman 1989, Grover 1997 |

## 7. Model Validation Status

This container does not have an R/mrgsolve execution environment installed (no `Rscript`), so the
mrgsolve model was completed only through the stage of **literature-based parameter design and a
self-review of the code (compartment/parameter consistency, bracket/dimension checks)**; the numbers
were not verified by actually compiling and integrating the model. The `.dot` file was actually
rendered into SVG/PNG and checked using a local installation of Graphviz `dot` (136 nodes, 12
clusters). Each of the 45 PubMed references had its PMID individually verified against the actual
PubMed page (the MCO-010 RESTORE trial had no peer-reviewed paper at the time of writing, so it is
mentioned only at the level of a conference abstract and is not listed as a separate reference). It is
recommended that anyone with an mrgsolve/R environment run the model as described in "Usage" above to
confirm the numerical integration results.
