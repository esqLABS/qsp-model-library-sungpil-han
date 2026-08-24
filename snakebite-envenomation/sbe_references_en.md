# Venomous snakebite (envenoming) + antivenom QSP model — references
# Snakebite Envenoming and Antivenom QSP Model — References

This list is the evidence for the structural choices and parameters that went into `sbe_qsp_model_en.dot` (the mechanistic map · 145 nodes · 21 clusters),
`sbe_mrgsolve_model_en.R` (the 50-ODE model), and `sbe_reference_model_en.py` (an independent verification
implementation).

**Every entry is a real record looked up directly with the PubMed E-utilities (esearch +
esummary)** — the PMID, author, year, journal, and title are transferred exactly as returned, and not one
citation was written from memory.

Every entry below is a real record retrieved directly from PubMed via the NCBI
E-utilities API. PMIDs, first authors, years, journals and titles are
transcribed from that retrieval rather than written from memory.

---

## What this model took from the literature, and what it did not
## What is sourced, and what is not

Distinguishing the two honestly is half of what this file is for. The most dangerous thing in a QSP model is
a calibration parameter that reads like an observation.

**Sourced from the literature**

- Terminal half-life by antivenom fragment — ovine Fab 12–23 h, equine F(ab')₂ ≈ 130 h,
  equine whole IgG 90–200 h. The two-compartment parameters of the model were chosen to reproduce these
  values, and `report_product_pk()` actually computes them to confirm it.
- The **labelled potency** of Indian polyvalent antivenom (ASV) — neutralisation of 0.6 mg of *Daboia russelii*
  venom per mL, that is 6 mgNE per 10 mL vial. The antivenom unit of the model (mgNE) is that regulatory
  figure itself and not a fitted value.
- Mass fraction by venom family (venomics/proteomics: Echis 68% SVMP, Naja 55% 3FTx,
  Bungarus 66% PLA2, and so on).
- The diagnostic threshold at which the 20-minute whole blood clotting test (20WBCT) turns positive, a fibrinogen
  of about 0.5 g/L.
- A normal fibrinogen of 2.8 g/L, a half-life of 96 h, and the range of acute-phase upregulation of synthesis.
- The distinction between the quasi-reversible dissociation of the α-neurotoxins (tens of hours) and the
  **irreversible** destruction of the motor end plate by the β-neurotoxins (PLA2), together with the clinical fact
  that anticholinesterases are of no use against the latter.
- The frequency of dry bites, the frequency of acute antivenom reactions and of serum sickness and how it differs by product,
  and the randomised trial result that heparin is of no use in VICC.
- The clinical time course of Russell's viper bite — the timing of the onset and recovery of VICC, the signs of neurotoxicity,
  and the 3–5 day peak of acute kidney injury.

**Model-level hypotheses, NOT observations**

- **The ε vector (neutralisation efficiency by family) = SVMP 1.20 / SVSP 1.00 / PLA2 0.60 / 3FTx 0.35.**
  That antivenom covers the large enzymes well and the small, poorly immunogenic toxins badly is
  a qualitative conclusion of the antivenomics literature, but **these four numbers are not quantitatively reported
  values; they are model parameters.** They were chosen so that the arithmetic reproduces the clinical practice
  of needing 20–30 vials for a cobra bite.
- **The slow reservoir (f_slow, ka_s).** The observation that venom antigen is detectable for days is real, but the
  compartmentalisation in which "35% of the injected amount is sequestered in a deep compartment and released with a t½ of 45 h" is a
  structural hypothesis of the model. The sensitivity analysis of §E tests the limits of this hypothesis explicitly.
- The exact size of the local neutralisation penalty at the bite site (k_b,loc/k_b0 = 0.038).
- The neuromuscular safety factor SF₀ = 4.0, the ptosis threshold SF < 2.2, and the ventilation threshold
  VC_frac < 0.15 — the order and the direction are robust but the absolute values are calibrated.
- Every rate constant in the renal loss terms (k_n,fib · k_n,cast · k_n,isch · k_n,dir),
  and the assumption that 10% of the loss becomes permanent scarring.
- The baseline hazard of the death and haemorrhage hazard functions.
- The dissociation rate of the venom-antivenom complex, k_off = 0 (the default). §E sweeps this value to
  compute what it would have to be in order to explain clinical recurrent coagulopathy.

The **ratios, orderings, and time windows** that are reported are far more robust than the absolute values. All four
claims of this model are about ratios.

---

## 1. Burden of disease · epidemiology · global health strategy
## Burden, epidemiology and the global health strategy

Why this is a pharmacology problem and not only a first-aid problem: the deaths
are concentrated where the drug supply and the ventilator are not.

- GBD 2019 Snakebite Envenomation Collaborators. **Global mortality of snakebite envenoming between 1990 and 2019.** *Nat Commun*. 2022. [PMID 36284094](https://pubmed.ncbi.nlm.nih.gov/36284094/)
- Afroz A, et al. **Snakebite envenoming: A systematic review and meta-analysis of global morbidity and mortality.** *PLoS Negl Trop Dis*. 2024. [PMID 38574167](https://pubmed.ncbi.nlm.nih.gov/38574167/)
- Longbottom J, et al. **Vulnerability to snakebite envenoming: a global mapping of hotspots.** *Lancet*. 2018. [PMID 30017551](https://pubmed.ncbi.nlm.nih.gov/30017551/)
- Iliyasu G, et al. **Case fatality rate and burden of snakebite envenoming in children — A systematic review and meta-analysis.** *Toxicon*. 2023. [PMID 37739273](https://pubmed.ncbi.nlm.nih.gov/37739273/)
- Gololo AA, et al. **Epidemiological models to estimate the burden of snakebite envenoming: A systematic review.** *Trop Med Int Health*. 2025. [PMID 39743841](https://pubmed.ncbi.nlm.nih.gov/39743841/)
- Menon JC, et al. **ICMR task force project — survey of the incidence, mortality, morbidity and socio-economic burden of snakebite in India: A study protocol.** *PLoS One*. 2022. [PMID 35994445](https://pubmed.ncbi.nlm.nih.gov/35994445/)
- Alcoba G, et al. **Snakebite epidemiology and health-seeking behavior in Akonolinga health district, Cameroon: Cross-sectional study.** *PLoS Negl Trop Dis*. 2020. [PMID 32584806](https://pubmed.ncbi.nlm.nih.gov/32584806/)
- Pach S, et al. **Paediatric snakebite envenoming: the world's most neglected 'Neglected Tropical Disease'?** *Arch Dis Child*. 2020. [PMID 32998874](https://pubmed.ncbi.nlm.nih.gov/32998874/)
- Minghui R, et al. **WHO's Snakebite Envenoming Strategy for prevention and control.** *Lancet Glob Health*. 2019. [PMID 31129124](https://pubmed.ncbi.nlm.nih.gov/31129124/)
- Williams DJ, et al. **Strategy for a globally coordinated response to a priority neglected tropical disease: Snakebite envenoming.** *PLoS Negl Trop Dis*. 2019. [PMID 30789906](https://pubmed.ncbi.nlm.nih.gov/30789906/)
- Chippaux JP. **The WHO strategy for prevention and control of snakebite envenoming: a sub-Saharan Africa plan.** *J Venom Anim Toxins Incl Trop Dis*. 2019. [PMID 31839803](https://pubmed.ncbi.nlm.nih.gov/31839803/)
- Harrison RA, et al. **The time is now: a call for action to translate recent momentum on tackling tropical snakebite into sustained benefit for victims.** *Trans R Soc Trop Med Hyg*. 2019. [PMID 30668842](https://pubmed.ncbi.nlm.nih.gov/30668842/)

## 2. Integrative reviews (the overall skeleton of the model)
## Integrative reviews — the skeleton the map is built on

- Gutiérrez JM, Calvete JJ, Habib AG, Harrison RA, Williams DJ, Warrell DA. **Snakebite envenoming.** *Nat Rev Dis Primers*. 2017. [PMID 28905944](https://pubmed.ncbi.nlm.nih.gov/28905944/)
- Gutiérrez JM, et al. **Understanding and confronting snakebite envenoming: The harvest of cooperation.** *Toxicon*. 2016. [PMID 26615826](https://pubmed.ncbi.nlm.nih.gov/26615826/)
- Clare RH, et al. **Small Molecule Drug Discovery for Neglected Tropical Snakebite.** *Trends Pharmacol Sci*. 2021. [PMID 33773806](https://pubmed.ncbi.nlm.nih.gov/33773806/)

## 3. The venom proteome — the basis for dividing it into four families
## The venom proteome — the basis for four toxin classes

The mass fractions in `SNAKES` come from this literature. They are the reason
*Echis* is a pure coagulopathy, *Naja* a pure paralysis, and *Daboia* both.

- Calvete JJ, et al. **Exploring the venom proteome of the western diamondback rattlesnake, Crotalus atrox, via snake venomics and combinatorial peptide ligand library approaches.** *J Proteome Res*. 2009. [PMID 19371136](https://pubmed.ncbi.nlm.nih.gov/19371136/)
- Serrano SM, et al. **A multifaceted analysis of viperid snake venoms by two-dimensional gel electrophoresis: an approach to understanding venom proteomics.** *Proteomics*. 2005. [PMID 15627971](https://pubmed.ncbi.nlm.nih.gov/15627971/)
- Slagboom J, et al. **High-Throughput Venomics.** *J Proteome Res*. 2023. [PMID 37010854](https://pubmed.ncbi.nlm.nih.gov/37010854/)
- Oh AMF, et al. **Venomics of Bungarus caeruleus (Indian krait): Comparable venom profiles, variable immunoreactivities among specimens from Sri Lanka, India and Pakistan.** *J Proteomics*. 2017. [PMID 28476572](https://pubmed.ncbi.nlm.nih.gov/28476572/)
- Hia YL, et al. **Comparative venom proteomics of banded krait (Bungarus fasciatus) from five geographical locales.** *Acta Trop*. 2020. [PMID 32278639](https://pubmed.ncbi.nlm.nih.gov/32278639/)
- Bittenbinder MA, et al. **Coagulotoxic Cobras: Clinical Implications of Strong Anticoagulant Actions of African Spitting Naja Venoms That Are Not Neutralised by Antivenom but Are by LY315920 (Varespladib).** *Toxins (Basel)*. 2018. [PMID 30518149](https://pubmed.ncbi.nlm.nih.gov/30518149/)

## 4. SVMP — the metalloproteinases: haemorrhagins · prothrombin/FX activation · digestion of the basement membrane
## Snake venom metalloproteinases

- Kamiguti AS, Hay CR, Theakston RD, Zuzel M. **Insights into the mechanism of haemorrhage caused by snake venom metalloproteinases.** *Toxicon*. 1996. [PMID 8817809](https://pubmed.ncbi.nlm.nih.gov/8817809/)
- Moura-da-Silva AM, et al. **Jararhagin, a hemorrhagic snake venom metalloproteinase from Bothrops jararaca.** *Toxicon*. 2012. [PMID 22534074](https://pubmed.ncbi.nlm.nih.gov/22534074/)
- Valente RH, et al. **BJ46a, a snake venom metalloproteinase inhibitor. Isolation, characterization, cloning and insights into its mechanism of action.** *Eur J Biochem*. 2001. [PMID 11358523](https://pubmed.ncbi.nlm.nih.gov/11358523/)
- Bickler PE. **Amplification of Snake Venom Toxicity by Endogenous Signaling Pathways.** *Toxins (Basel)*. 2020. [PMID 31979014](https://pubmed.ncbi.nlm.nih.gov/31979014/)

## 5. SVSP — the serine proteinases: thrombin-like enzymes and the consumption of fibrinogen
## Snake venom serine proteases and thrombin-like enzymes

The `kcat_svsp_fg` term and the Michaelis–Menten form of fibrinogen consumption
come from this pharmacology: these enzymes cleave fibrinopeptide A only, so the
fibrin they make is friable and the fibrinogen they consume is simply gone.

- Stocker K, Barlow GH. **Thrombin-like snake venom proteinases.** *Toxicon*. 1982. [PMID 7043783](https://pubmed.ncbi.nlm.nih.gov/7043783/)
- Aronson DL. **Comparison of the actions of thrombin and the thrombin-like venom enzymes ancrod and batroxobin.** *Thromb Haemost*. 1976. [PMID 1036831](https://pubmed.ncbi.nlm.nih.gov/1036831/)
- Hahn BS, et al. **Purification and molecular cloning of calobin, a thrombin-like enzyme from Agkistrodon caliginosus (Korean viper).** *J Biochem*. 1996. [PMID 8797081](https://pubmed.ncbi.nlm.nih.gov/8797081/)
- Cho SY, et al. **Purification and characterization of calobin II, a second type of thrombin-like enzyme from Agkistrodon caliginosus (Korean viper).** *Toxicon*. 2001. [PMID 11024490](https://pubmed.ncbi.nlm.nih.gov/11024490/)

## 6. svPLA2 — myotoxicity · presynaptic neurotoxicity · oedema
## Secreted phospholipase A₂: myotoxicity, presynaptic neurotoxicity, oedema

- Gutiérrez JM, Rucavado A, Chaves F, Díaz C, Escalante T. **Experimental pathology of local tissue damage induced by Bothrops asper snake venom.** *Toxicon*. 2009. [PMID 19303033](https://pubmed.ncbi.nlm.nih.gov/19303033/)
- Gutiérrez JM, León G, Rojas G, Lomonte B, Rucavado A, Chaves F. **Neutralization of local tissue damage induced by Bothrops asper (terciopelo) snake venom.** *Toxicon*. 1998. [PMID 9792169](https://pubmed.ncbi.nlm.nih.gov/9792169/)
- Chaves F, et al. **Histopathological and biochemical alterations induced by intramuscular injection of Bothrops asper (terciopelo) venom in mice.** *Toxicon*. 1989. [PMID 2815106](https://pubmed.ncbi.nlm.nih.gov/2815106/)
- Rucavado A, et al. **Analysis of wound exudates reveals differences in the patterns of tissue damage and inflammation induced by the venoms of Daboia russelii and Bothrops asper.** *Toxicon*. 2020. [PMID 32781076](https://pubmed.ncbi.nlm.nih.gov/32781076/)
- Rucavado A, et al. **Increments in cytokines and matrix metalloproteinases in skeletal muscle after injection of tissue-damaging toxins from the venom of the snake Bothrops asper.** *Mediators Inflamm*. 2002. [PMID 12061424](https://pubmed.ncbi.nlm.nih.gov/12061424/)
- Silva A, et al. **Clinical and Pharmacological Investigation of Myotoxicity in Sri Lankan Russell's Viper (Daboia russelii) Envenoming.** *PLoS Negl Trop Dis*. 2016. [PMID 27911900](https://pubmed.ncbi.nlm.nih.gov/27911900/)
- Sanhajariya S, et al. **Investigating myotoxicity following Australian red-bellied black snake (Pseudechis porphyriacus) envenomation.** *PLoS One*. 2021. [PMID 34506531](https://pubmed.ncbi.nlm.nih.gov/34506531/)
- Silva MDS, et al. **NLRP3 inflammasome activation in human peripheral blood mononuclear cells induced by venoms secreted PLA₂s.** *Int J Biol Macromol*. 2022. [PMID 35074331](https://pubmed.ncbi.nlm.nih.gov/35074331/)

## 7. 3FTx — the three-finger toxins: postsynaptic α-neurotoxins and the cytotoxins
## Three-finger toxins: postsynaptic α-neurotoxins and cytotoxins

The apparent perijunctional K_d and the slow off-rate (t½ ≈ 35 h, and ~230 h for
α-bungarotoxin) come from this literature, and they are why an occupancy decays
on its own clock rather than being stripped off by antivenom.

- Nirthanan S. **Snake three-finger α-neurotoxins and nicotinic acetylcholine receptors: molecules, mechanisms and medicine.** *Biochem Pharmacol*. 2020. [PMID 32710970](https://pubmed.ncbi.nlm.nih.gov/32710970/)
- Nirthanan S, Gwee MC. **Three-finger alpha-neurotoxins and the nicotinic acetylcholine receptor, forty years on.** *J Pharmacol Sci*. 2004. [PMID 14745112](https://pubmed.ncbi.nlm.nih.gov/14745112/)
- Shenkarev ZO, et al. **Membrane-mediated interaction of non-conventional snake three-finger toxins with nicotinic acetylcholine receptors.** *Commun Biol*. 2022. [PMID 36477694](https://pubmed.ncbi.nlm.nih.gov/36477694/)
- Bekbossynova A, et al. **Venom-Derived Neurotoxins Targeting Nicotinic Acetylcholine Receptors.** *Molecules*. 2021. [PMID 34204855](https://pubmed.ncbi.nlm.nih.gov/34204855/)
- Choudhury M, et al. **Snake Venom Three-Finger Neurotoxins and Neurotoxin-Like Proteins: Insights Into Their Structural and Functional Aspects.** *Basic Clin Pharmacol Toxicol*. 2025. [PMID 41230907](https://pubmed.ncbi.nlm.nih.gov/41230907/)
- Tamiya N, Sato A. **Studies on sea snake venom.** *Proc Jpn Acad Ser B Phys Biol Sci*. 2011. [PMID 21422738](https://pubmed.ncbi.nlm.nih.gov/21422738/)

## 8. β-neurotoxins — presynaptic destruction: the injury antivenom cannot reverse
## β-neurotoxins: the presynaptic destruction antivenom cannot undo

This section is the evidentiary basis of the model's claim 3. The `dTERM/dt`
equation contains a destruction term and a 5-day regeneration term, and neither
antivenom nor neostigmine appears in it.

- Galappaththige J, et al. **Neurotoxicity of Sri Lankan Krait (Bungarus ceylonicus) and Common Krait (Bungarus caeruleus) Venoms and Their Neutralisation by Commercial Antivenoms.** *Toxins (Basel)*. 2025. [PMID 41003503](https://pubmed.ncbi.nlm.nih.gov/41003503/)
- Liang Q, et al. **In Vitro Neurotoxicity of Chinese Krait (Bungarus multicinctus) Venom and Neutralization by Antivenoms.** *Toxins (Basel)*. 2021. [PMID 33440641](https://pubmed.ncbi.nlm.nih.gov/33440641/)
- Herkert M, et al. **Beta-bungarotoxin is a potent inducer of apoptosis in cultured rat neurons by receptor-mediated internalization.** *Eur J Neurosci*. 2001. [PMID 11576186](https://pubmed.ncbi.nlm.nih.gov/11576186/)
- Othman IB, Spokes JW, Dolly JO. **Preparation of neurotoxic ³H-beta-bungarotoxin: demonstration of saturable binding to brain synapses and its inhibition by toxin I.** *Eur J Biochem*. 1982. [PMID 7173209](https://pubmed.ncbi.nlm.nih.gov/7173209/)

## 9. The clinical picture of the neurotoxic bite — kraits and cobras
## Clinical neurotoxic envenoming — krait and cobra

The ventilation durations (4–10 days), the delayed krait onset, and the
observation that anticholinesterases help the cobra and not the krait, are all
read off these series.

- Kularatne SA. **Common krait (Bungarus caeruleus) bite in Anuradhapura, Sri Lanka: a prospective clinical study, 1996-98.** *Postgrad Med J*. 2002. [PMID 12151569](https://pubmed.ncbi.nlm.nih.gov/12151569/)
- Bawaskar HS, Bawaskar PH. **Envenoming by the common krait (Bungarus caeruleus) and Asian cobra (Naja naja): clinical manifestations and their management in a rural setting.** *Wilderness Environ Med*. 2004. [PMID 15636376](https://pubmed.ncbi.nlm.nih.gov/15636376/)
- Theakston RD, et al. **Envenoming by the common krait (Bungarus caeruleus) and Sri Lankan cobra (Naja naja naja): efficacy and complications of therapy with Haffkine antivenom.** *Trans R Soc Trop Med Hyg*. 1990. [PMID 2389328](https://pubmed.ncbi.nlm.nih.gov/2389328/)
- Pannu AK, et al. **Efficacy of a low dose of antivenom for severe neuroparalysis in Bungarus caeruleus (common krait) envenomation: a pilot study.** *Toxicol Res (Camb)*. 2024. [PMID 38450179](https://pubmed.ncbi.nlm.nih.gov/38450179/)
- Gupta A, et al. **Unusually prolonged neuromuscular weakness caused by krait (Bungarus caeruleus) bite: Two case reports.** *Toxicon*. 2021. [PMID 33497743](https://pubmed.ncbi.nlm.nih.gov/33497743/)
- Faiz MA, et al. **Bites by the Monocled Cobra, Naja kaouthia, in Chittagong Division, Bangladesh: Epidemiology, Clinical Features of Envenoming and Management of 70 Identified Cases.** *Am J Trop Med Hyg*. 2017. [PMID 28138054](https://pubmed.ncbi.nlm.nih.gov/28138054/)
- Sarin K, et al. **Clinical profile & complications of neurotoxic snake bite & comparison of two regimens of polyvalent anti-snake venom in its treatment.** *Indian J Med Res*. 2017. [PMID 28574015](https://pubmed.ncbi.nlm.nih.gov/28574015/)
- Silva A, et al. **Neurotoxicity in Russell's viper (Daboia russelii) envenoming in Sri Lanka: a clinical and neurophysiological study.** *Clin Toxicol (Phila)*. 2016. [PMID 26923566](https://pubmed.ncbi.nlm.nih.gov/26923566/)

## 10. Venom-induced consumption coagulopathy (VICC)
## Venom-induced consumption coagulopathy

- Maduwage K, Isbister GK. **Current treatment for venom-induced consumption coagulopathy resulting from snakebite.** *PLoS Negl Trop Dis*. 2014. [PMID 25340841](https://pubmed.ncbi.nlm.nih.gov/25340841/)
- Maduwage K, Buckley NA, de Silva HJ, Lalloo DG, Isbister GK. **Snake antivenom for snake venom induced consumption coagulopathy.** *Cochrane Database Syst Rev*. 2015. [PMID 26058967](https://pubmed.ncbi.nlm.nih.gov/26058967/)
- Wedasingha S, et al. **Bedside Coagulation Tests in Diagnosing Venom-Induced Consumption Coagulopathy in Snakebite.** *Toxins (Basel)*. 2020. [PMID 32927702](https://pubmed.ncbi.nlm.nih.gov/32927702/)
- Valenta J, et al. **Fibrinogenolysis in Venom-Induced Consumption Coagulopathy after Viperidae Snakebites: A Pilot Study.** *Toxins (Basel)*. 2022. [PMID 36006200](https://pubmed.ncbi.nlm.nih.gov/36006200/)
- Little M, Isbister GK. **Is D-dimer the new test for venom-induced consumption coagulopathy after snakebite?** *Med J Aust*. 2022. [PMID 35843724](https://pubmed.ncbi.nlm.nih.gov/35843724/)
- Rajkumar B, et al. **Venom induced consumption coagulopathy and performance of 20-min whole blood clotting test for its detection in viperid envenomation.** *J R Coll Physicians Edinb*. 2022. [PMID 36300884](https://pubmed.ncbi.nlm.nih.gov/36300884/)
- Silva A, et al. **Indian Polyvalent Antivenom Accelerates Recovery From Venom-Induced Consumption Coagulopathy (VICC) in Sri Lankan Russell's Viper (Daboia russelii) Envenoming.** *Front Med (Lausanne)*. 2022. [PMID 35321467](https://pubmed.ncbi.nlm.nih.gov/35321467/)

**This section is the basis for claim 1.** Brown 2009 in particular deals directly, in the clinic, with the structure of this
model in which factor replacement changes the 'depth' of the recovery rather than its 'rate'.

- Brown SG, Caruso N, Borland ML, McCoubrie DL, Celenza A, Isbister GK. **Clotting factor replacement and recovery from snake venom-induced consumptive coagulopathy.** *Intensive Care Med*. 2009. [PMID 19547954](https://pubmed.ncbi.nlm.nih.gov/19547954/)

## 11. The 20-minute whole blood clotting test — the gate that costs nothing
## The 20-minute whole blood clotting test

The model's `FG_unclot = 0.5 g/L` threshold and the "dry bite gate" argument in
D6 rest on this diagnostic literature.

- Lamb T, et al. **The 20-minute whole blood clotting test (20WBCT) for snakebite coagulopathy — A systematic review and meta-analysis of diagnostic test accuracy.** *PLoS Negl Trop Dis*. 2021. [PMID 34375338](https://pubmed.ncbi.nlm.nih.gov/34375338/)
- Lamb T, et al. **Correction: The 20-minute whole blood clotting test (20WBCT) for snakebite coagulopathy.** *PLoS Negl Trop Dis*. 2023. [PMID 36693087](https://pubmed.ncbi.nlm.nih.gov/36693087/)
- Thongtonyong N, Chinthammitr Y. **Sensitivity and specificity of 20-minute whole blood clotting test, prothrombin time, activated partial thromboplastin time tests in diagnosis of defibrination syndrome.** *Toxicon*. 2020. [PMID 32712023](https://pubmed.ncbi.nlm.nih.gov/32712023/)
- Tianyi FL, et al. **Diagnostic characteristics of the 20-minute whole blood clotting test in detecting venom-induced consumptive coagulopathy following carpet viper envenoming.** *PLoS Negl Trop Dis*. 2023. [PMID 37363905](https://pubmed.ncbi.nlm.nih.gov/37363905/)
- Hamza M, et al. **Performance of the 20 minutes Whole Blood Clotting Test in detection, monitoring and antivenom therapy of West African Carpet viper (Echis romani) envenoming.** *Toxicon*. 2023. [PMID 36640811](https://pubmed.ncbi.nlm.nih.gov/36640811/)
- Suseel A, et al. **Comparing modified Lee and White method against 20-minute whole blood clotting test as bedside coagulation screening test in snake envenomation victims.** *J Venom Anim Toxins Incl Trop Dis*. 2023. [PMID 37342654](https://pubmed.ncbi.nlm.nih.gov/37342654/)

## 12. Thrombocytopenia — the C-type lectins and GPIb
## Venom-induced thrombocytopenia

- Lu Q, Clemetson JM, Clemetson KJ. **Snake venom C-type lectins interacting with platelet receptors. Structure-function relationships and effects on haemostasis.** *Toxicon*. 2005. [PMID 15876445](https://pubmed.ncbi.nlm.nih.gov/15876445/)
- Long C, et al. **Potential Role of Platelet-Activating C-Type Lectin-Like Proteins in Viper Envenomation Induced Thrombotic Microangiopathy Symptom.** *Toxins (Basel)*. 2020. [PMID 33260875](https://pubmed.ncbi.nlm.nih.gov/33260875/)
- Shen C, et al. **Viper venoms drive the macrophages and hepatocytes to sequester and clear platelets: novel mechanism and therapeutic strategy for venom-induced thrombocytopenia.** *Arch Toxicol*. 2021. [PMID 34519865](https://pubmed.ncbi.nlm.nih.gov/34519865/)
- Thomazini CM, et al. **Involvement of von Willebrand factor and botrocetin in the thrombocytopenia induced by Bothrops jararaca snake venom.** *PLoS Negl Trop Dis*. 2021. [PMID 34478462](https://pubmed.ncbi.nlm.nih.gov/34478462/)
- Xu G, et al. **How does agkicetin-C bind on platelet glycoprotein Ibalpha and achieve its platelet effects?** *Toxicon*. 2005. [PMID 15777951](https://pubmed.ncbi.nlm.nih.gov/15777951/)

## 13. The pharmacokinetics of venom and of antivenom — the PK skeleton of this model
## Venom and antivenom pharmacokinetics

Isbister 2015 is the anchor for the antivenom disposition parameters; Morris 2022
is the closest prior computational model of the same system and the reason this
model was written as two clocks rather than one.

- Isbister GK, et al. **Population Pharmacokinetics of an Indian F(ab')₂ Snake Antivenom in Patients with Russell's Viper (Daboia russelii) Bites.** *PLoS Negl Trop Dis*. 2015. [PMID 26135318](https://pubmed.ncbi.nlm.nih.gov/26135318/)
- Sanhajariya S, et al. **Population pharmacokinetics of Pseudechis porphyriacus (red-bellied black snake) venom in snakebite patients.** *Clin Toxicol (Phila)*. 2021. [PMID 33832399](https://pubmed.ncbi.nlm.nih.gov/33832399/)
- Morris NM, et al. **Developing a computational pharmacokinetic model of systemic snakebite envenomation and antivenom treatment.** *Toxicon*. 2022. [PMID 35716719](https://pubmed.ncbi.nlm.nih.gov/35716719/)
- Audebert F, et al. **Viper bites in France: clinical and biological evaluation; kinetics of envenomations.** *Hum Exp Toxicol*. 1994. [PMID 7826686](https://pubmed.ncbi.nlm.nih.gov/7826686/)
- Theakston RD. **Snake venoms in science and clinical medicine. 2. Applied immunology in snake venom research.** *Trans R Soc Trop Med Hyg*. 1989. [PMID 2617643](https://pubmed.ncbi.nlm.nih.gov/2617643/)
- Chavez-Olortegui C, et al. **An enzyme linked immunosorbent assay (ELISA) that discriminates between Bothrops atrox and Lachesis muta muta venoms.** *Toxicon*. 1993. [PMID 8503131](https://pubmed.ncbi.nlm.nih.gov/8503131/)

## 14. Recurrence — a short fragment missing a long tail
## Recurrence after a short-lived fragment

**The clinical basis for claim 2.** Recurrent coagulopathy, thrombocytopenia, and antigenaemia after ovine Fab have been
reported repeatedly, and the maintenance dose schedule on the CroFab label (2 vials q6h ×3) is the
consequence. §E states explicitly the negative result that the *severe* recurrence this literature reports is not
explained by reservoir release alone in this model.

- Ruha AM, Curry SC, Albrecht C, Riley B, Pizon A. **Late hematologic toxicity following treatment of rattlesnake envenomation with crotalidae polyvalent immune Fab antivenom.** *Toxicon*. 2011. [PMID 20920516](https://pubmed.ncbi.nlm.nih.gov/20920516/)
- Miller AD, et al. **Recurrent coagulopathy and thrombocytopenia in children treated with crotalidae polyvalent immune fab: a case series.** *Pediatr Emerg Care*. 2010. [PMID 20693856](https://pubmed.ncbi.nlm.nih.gov/20693856/)
- Fazelat J, Teperman SH, Touger M. **Recurrent hemorrhage after western diamondback rattlesnake envenomation treated with crotalidae polyvalent immune fab (ovine).** *Clin Toxicol (Phila)*. 2008. [PMID 18608290](https://pubmed.ncbi.nlm.nih.gov/18608290/)
- Lavonas EJ, et al. **Initial experience with Crotalidae polyvalent immune Fab (ovine) antivenom in the treatment of copperhead snakebite.** *Ann Emerg Med*. 2004. [PMID 14747809](https://pubmed.ncbi.nlm.nih.gov/14747809/)
- Keating GM. **Crotalidae polyvalent immune Fab: in patients with North American crotaline envenomation.** *BioDrugs*. 2011. [PMID 21443271](https://pubmed.ncbi.nlm.nih.gov/21443271/)
- Keating GM. **Crotalidae polyvalent immune Fab: a guide to its use in North American crotaline envenomation.** *Clin Drug Investig*. 2012. [PMID 22765769](https://pubmed.ncbi.nlm.nih.gov/22765769/)
- Boyer L, et al. **Safety of intravenous equine F(ab')₂: insights following clinical trials involving 1534 recipients of scorpion antivenom.** *Toxicon*. 2013. [PMID 23916602](https://pubmed.ncbi.nlm.nih.gov/23916602/)

## 15. Antivenom dose · potency · clinical trials
## Antivenom dose, potency and trials

The labelled-potency unit (mgNE) and the "10 vials is a 0.98× margin" arithmetic
in section B come from this literature.

- Abubakar IS, et al. **Randomised controlled double-blind non-inferiority trial of two antivenoms for saw-scaled or carpet viper (Echis ocellatus) envenoming in Nigeria.** *PLoS Negl Trop Dis*. 2010. [PMID 20668549](https://pubmed.ncbi.nlm.nih.gov/20668549/)
- Abubakar SB, et al. **Pre-clinical and preliminary dose-finding and safety studies to identify candidate antivenoms for treatment of envenoming by saw-scaled or carpet vipers (Echis ocellatus) in Nigeria.** *Toxicon*. 2010. [PMID 19874841](https://pubmed.ncbi.nlm.nih.gov/19874841/)
- Meyer WP, et al. **First clinical experiences with a new ovine Fab Echis ocellatus snake bite antivenom in Nigeria: randomized comparative trial with Institute Pasteur Serum.** *Am J Trop Med Hyg*. 1997. [PMID 9129531](https://pubmed.ncbi.nlm.nih.gov/9129531/)
- Sagar P, et al. **Comparison of two Anti Snake Venom protocols in hemotoxic snake bite: A randomized trial.** *J Forensic Leg Med*. 2020. [PMID 32658754](https://pubmed.ncbi.nlm.nih.gov/32658754/)
- Patra A, Mukherjee AK. **Assessment of quality and pre-clinical efficacy of a newly developed polyvalent antivenom against the medically important snakes of Sri Lanka.** *Sci Rep*. 2021. [PMID 34521877](https://pubmed.ncbi.nlm.nih.gov/34521877/)
- Lim ASS, et al. **Immunoreactivity and neutralization efficacy of Pakistani Viper Antivenom (PVAV) against venoms of Saw-scaled Vipers (Echis carinatus subspp.) and Western Russell's Viper.** *Acta Trop*. 2024. [PMID 38097152](https://pubmed.ncbi.nlm.nih.gov/38097152/)
- Chippaux JP. **[Guidelines for the production, control and regulation of snake antivenom immunoglobulins].** *Biol Aujourdhui*. 2010. [PMID 20950580](https://pubmed.ncbi.nlm.nih.gov/20950580/)
- World Health Organization. **WHO Expert Committee on Biological Standardization.** *World Health Organ Tech Rep Ser*. 2012. [PMID 22900409](https://pubmed.ncbi.nlm.nih.gov/22900409/)
- Scheske L, Ruitenberg J, Bissumbhar B. **Needs and availability of snake antivenoms: relevance and application of international guidelines.** *Int J Health Policy Manag*. 2015. [PMID 26188809](https://pubmed.ncbi.nlm.nih.gov/26188809/)

## 16. Antivenom adverse reactions — the cost side of the ledger
## Antivenom adverse reactions — the cost side of the ledger

The rate-driven anaphylactoid term (`rho * AVRATE`) and the immune-complex →
serum-sickness cascade are parameterised against these product-specific rates.

- Habib AG. **Effect of pre-medication on early adverse reactions following antivenom use in snakebite: a systematic review and meta-analysis.** *Drug Saf*. 2011. [PMID 21879781](https://pubmed.ncbi.nlm.nih.gov/21879781/)
- Morais V. **Antivenom therapy: efficacy of premedication for the prevention of adverse reactions.** *J Venom Anim Toxins Incl Trop Dis*. 2018. [PMID 29507580](https://pubmed.ncbi.nlm.nih.gov/29507580/)
- Ryan NM, Kearney RT, Brown SG, Isbister GK. **Incidence of serum sickness after the administration of Australian snake antivenom (ASP-22).** *Clin Toxicol (Phila)*. 2016. [PMID 26490786](https://pubmed.ncbi.nlm.nih.gov/26490786/)
- Tu M, et al. **Safety profile of antivenom in a cohort of patients envenomed by Deinagkistrodon acutus in Hangzhou, Zhejiang Province, Southeast China.** *Clin Toxicol (Phila)*. 2025. [PMID 40106271](https://pubmed.ncbi.nlm.nih.gov/40106271/)
- LoVecchio F, et al. **Incidence of immediate and delayed hypersensitivity to Centruroides antivenom.** *Ann Emerg Med*. 1999. [PMID 10533009](https://pubmed.ncbi.nlm.nih.gov/10533009/)
- Sheikh A. **Glucocorticosteroids for the treatment and prevention of anaphylaxis.** *Curr Opin Allergy Clin Immunol*. 2013. [PMID 23507835](https://pubmed.ncbi.nlm.nih.gov/23507835/)

## 17. Venom-induced acute kidney injury — the kidney as an integral
## Venom-induced acute kidney injury

Claim 4. The day-3-to-5 creatinine peak, the fibrin-microthrombus and pigment
cast pathways, and bilateral cortical necrosis as the irreversible end of the
same integral, come from this literature.

- Sitprija V. **Snakebite nephropathy.** *Nephrology (Carlton)*. 2006. [PMID 17014559](https://pubmed.ncbi.nlm.nih.gov/17014559/)
- Chugh KS. **Snake-bite-induced acute renal failure in India.** *Kidney Int*. 1989. [PMID 2651763](https://pubmed.ncbi.nlm.nih.gov/2651763/)
- Vikrant S, Jaryal A, Parashar A. **Clinicopathological spectrum of snake bite-induced acute kidney injury from India.** *World J Nephrol*. 2017. [PMID 28540205](https://pubmed.ncbi.nlm.nih.gov/28540205/)
- Rao PSK, et al. **Snakebite envenomation-associated acute kidney injury: a South-Asian perspective.** *Trans R Soc Trop Med Hyg*. 2025. [PMID 39749490](https://pubmed.ncbi.nlm.nih.gov/39749490/)
- Alvitigala BY, et al. **Snakebite-associated acute kidney injury in South Asia: narrative review on epidemiology, pathogenesis and management.** *Trans R Soc Trop Med Hyg*. 2025. [PMID 39749470](https://pubmed.ncbi.nlm.nih.gov/39749470/)
- Meena P, et al. **The kidney histopathological spectrum of patients with kidney injury following snakebite envenomation in India: scoping review of five decades.** *BMC Nephrol*. 2024. [PMID 38515042](https://pubmed.ncbi.nlm.nih.gov/38515042/)
- Prakash J, Singh TB, Ghosh B, et al. **Changing picture of renal cortical necrosis in acute kidney injury in developing country.** *World J Nephrol*. 2015. [PMID 26558184](https://pubmed.ncbi.nlm.nih.gov/26558184/)

## 18. Local injury · compartment syndrome · fasciotomy
## Local injury, compartment syndrome and fasciotomy

Claim 5: the compartment intravenous antibody cannot reach.

- Spyres MB, et al. **Compartment Syndrome after Crotalid Envenomation in the United States: A Review of the North American Snakebite Registry from 2013 to 2021.** *Wilderness Environ Med*. 2023. [PMID 37474357](https://pubmed.ncbi.nlm.nih.gov/37474357/)
- Kim YH, et al. **Fasciotomy in compartment syndrome from snakebite.** *Arch Plast Surg*. 2019. [PMID 30685944](https://pubmed.ncbi.nlm.nih.gov/30685944/)
- Chen C, et al. **Severe compartment syndrome following snakebite in a man: A case report.** *Asian J Surg*. 2024. [PMID 39237419](https://pubmed.ncbi.nlm.nih.gov/39237419/)
- Sassoè-Pognetto M, et al. **Acute compartment syndrome and fasciotomy after a viper bite in Italy: a case report.** *Ital J Pediatr*. 2024. [PMID 38627836](https://pubmed.ncbi.nlm.nih.gov/38627836/)
- Hou YT, et al. **Prediction of Compartment Syndrome after Protobothrops mucrosquamatus Snakebite by Diastolic Retrograde Arterial Flow: A Case Report.** *Medicina (Kaunas)*. 2022. [PMID 35893111](https://pubmed.ncbi.nlm.nih.gov/35893111/)

## 19. Russell's viper — the reference species of this model
## Daboia russelii — the model's reference species

- Kularatne SA. **Epidemiology and clinical picture of the Russell's viper (Daboia russelii russelii) bite in Anuradhapura, Sri Lanka: a prospective study of 336 patients.** *Southeast Asian J Trop Med Public Health*. 2003. [PMID 15115100](https://pubmed.ncbi.nlm.nih.gov/15115100/)
- Kularatne SA, et al. **Revisiting Russell's viper (Daboia russelii) bite in Sri Lanka: is abdominal pain an early feature of systemic envenoming?** *PLoS One*. 2014. [PMID 24587278](https://pubmed.ncbi.nlm.nih.gov/24587278/)

## 20. What does not work — heparin and tranexamic acid
## What does NOT work: heparin and tranexamic acid

**These are the clinical negative results that this model reproduces by calculation in D1.** Heparin cannot
regulate thrombin that is generated outside the cascade, and tranexamic acid has only 5% of the problem it can touch,
because 95% of the fibrinogen loss is direct enzymatic cleavage.

- Tin Na Swe, et al. **Heparin therapy in Russell's viper bite victims with disseminated intravascular coagulation: a controlled trial.** *Southeast Asian J Trop Med Public Health*. 1992. [PMID 1345132](https://pubmed.ncbi.nlm.nih.gov/1345132/)
- Myint-Lwin, et al. **Heparin therapy in Russell's viper bite victims with impending DIC (a controlled trial).** *Southeast Asian J Trop Med Public Health*. 1989. [PMID 2532790](https://pubmed.ncbi.nlm.nih.gov/2532790/)
- Paul V, et al. **Trial of heparin in viper bites.** *J Assoc Physicians India*. 2003. [PMID 12725259](https://pubmed.ncbi.nlm.nih.gov/12725259/)
- Paul V, et al. **Trial of low molecular weight heparin in the treatment of viper bites.** *J Assoc Physicians India*. 2007. [PMID 17844693](https://pubmed.ncbi.nlm.nih.gov/17844693/)
- Huang YN, et al. **Effects of heparin on venom-induced consumption coagulopathy: a meta-analysis of randomized controlled trials.** *Trans R Soc Trop Med Hyg*. 2025. [PMID 39749494](https://pubmed.ncbi.nlm.nih.gov/39749494/)
- Larréché S, et al. **[Evaluation of the Efficacy and Tolerance of Tranexamic Acid Combined with Inoserp™ PAN-AFRICA Antivenom in the Treatment of Hemorrhagic Syndrome].** *Med Trop Sante Int*. 2026. [PMID 42221451](https://pubmed.ncbi.nlm.nih.gov/42221451/)
- Yamamoto A, et al. **Therapeutic Effects of Single and Combined Anti-Disseminated Intravascular Coagulation (DIC) Drugs in a Rat Venom-Induced Consumption Coagulopathy Model.** *Toxins (Basel)*. 2026. [PMID 41893574](https://pubmed.ncbi.nlm.nih.gov/41893574/)
- Abraham SV, et al. **Hematotoxic Snakebite Victim with Trauma: The Role of Guided Transfusion, Rotational Thromboelastometry, and Tranexamic Acid.** *Wilderness Environ Med*. 2020. [PMID 33162320](https://pubmed.ncbi.nlm.nih.gov/33162320/)

## 21. Small molecule inhibitors — the drug that goes where antivenom cannot
## Small-molecule inhibitors — the drug that goes where antibody cannot

The `VAR_*` compartments and `INH_PLA2` factor, and the D5 result, are built on
this programme. The BRAVO trial is the reason varespladib is in the ODEs at all
rather than only in the map.

- Gerardo CJ, et al. **Oral varespladib for the treatment of snakebite envenoming in India and the USA (BRAVO): a phase II randomised clinical trial.** *BMJ Glob Health*. 2024. [PMID 39442939](https://pubmed.ncbi.nlm.nih.gov/39442939/)
- Carter RW, et al. **The BRAVO Clinical Study Protocol: Oral Varespladib for Inhibition of Secretory Phospholipase A2 in the Treatment of Snakebite Envenoming.** *Toxins (Basel)*. 2022. [PMID 36668842](https://pubmed.ncbi.nlm.nih.gov/36668842/)
- Quiroz S, et al. **Inhibitory Effects of Varespladib, CP471474, and Their Potential Synergistic Activity on Bothrops asper and Crotalus durissus cumanensis Venoms.** *Molecules*. 2022. [PMID 36500682](https://pubmed.ncbi.nlm.nih.gov/36500682/)
- Woliver C, et al. **Phospholipase A2 inhibitor may shorten the duration of clinical signs in the treatment of neurotoxicity caused by eastern coral snake (Micrurus fulvius) envenomation.** *J Am Vet Med Assoc*. 2025. [PMID 40738154](https://pubmed.ncbi.nlm.nih.gov/40738154/)

## 22. Next-generation antivenoms — recombinant antibodies and nanobodies
## Next-generation antivenom

Where the ε vector stops being a fact of nature and becomes a design variable.

- Kini RM, et al. **Biosynthetic Oligoclonal Antivenom (BOA) for Snakebite and Next-Generation Treatments for Snakebite Victims.** *Toxins (Basel)*. 2018. [PMID 30551565](https://pubmed.ncbi.nlm.nih.gov/30551565/)
- Knudsen C, Laustsen AH. **Engineering and design considerations for next-generation snakebite antivenoms.** *Toxicon*. 2019. [PMID 31173790](https://pubmed.ncbi.nlm.nih.gov/31173790/)
- Fernandes CFC, et al. **Engineering of single-domain antibodies for next-generation snakebite antivenoms.** *Int J Biol Macromol*. 2021. [PMID 34118288](https://pubmed.ncbi.nlm.nih.gov/34118288/)
- Jenkins TP, Laustsen AH. **Cost of Manufacturing for Recombinant Snakebite Antivenoms.** *Front Bioeng Biotechnol*. 2020. [PMID 32766215](https://pubmed.ncbi.nlm.nih.gov/32766215/)
- Jenkins TP, et al. **Toxin Neutralization Using Alternative Binding Proteins.** *Toxins (Basel)*. 2019. [PMID 30658491](https://pubmed.ncbi.nlm.nih.gov/30658491/)
- Knudsen C, et al. **Novel Snakebite Therapeutics Must Be Tested in Appropriate Rescue Models to Robustly Assess Their Preclinical Efficacy.** *Toxins (Basel)*. 2020. [PMID 32824899](https://pubmed.ncbi.nlm.nih.gov/32824899/)
- Tabares Vélez S, et al. **Standard Quality Characteristics and Efficacy of a New Third-Generation Antivenom Developed in Colombia Covering Micrurus spp. Venoms.** *Toxins (Basel)*. 2024. [PMID 38668608](https://pubmed.ncbi.nlm.nih.gov/38668608/)

## 23. First aid and lymphatic flow — an attempt to control the reservoir
## First aid and lymphatic flow — trying to control the depot

The pressure-immobilisation node in cluster 2, and the model's ka being a
lymphatic rate constant rather than a capillary one, come from here.

- Parker-Cote J, Meggs WJ. **First Aid and Pre-Hospital Management of Venomous Snakebites.** *Trop Med Infect Dis*. 2018. [PMID 30274441](https://pubmed.ncbi.nlm.nih.gov/30274441/)
- van Helden DF, et al. **Pharmacological approaches that slow lymphatic flow as a snakebite first aid.** *PLoS Negl Trop Dis*. 2014. [PMID 24587472](https://pubmed.ncbi.nlm.nih.gov/24587472/)
- Howarth DM, Southee AE, Whyte IM. **Lymphatic flow rates and first-aid in simulated peripheral snake or spider envenomation.** *Med J Aust*. 1994. [PMID 7830641](https://pubmed.ncbi.nlm.nih.gov/7830641/)

## 24. Supply · cost · access — the variables that decide the real outcome
## Supply, cost and access — the variable that actually decides outcome

A QSP model can compute the optimal vial count; it cannot compute whether the
vial is in the building. These papers are why the "no ventilator" and "no
dialysis" switches exist in the model at all.

- Hamza M, et al. **Cost-Effectiveness of Antivenoms for Snakebite Envenoming in 16 Countries in West Africa.** *PLoS Negl Trop Dis*. 2016. [PMID 27027633](https://pubmed.ncbi.nlm.nih.gov/27027633/)
- Habib AG, et al. **Cost-effectiveness of antivenoms for snakebite envenoming in Nigeria.** *PLoS Negl Trop Dis*. 2015. [PMID 25569252](https://pubmed.ncbi.nlm.nih.gov/25569252/)
- Ooms GI, et al. **Availability, affordability and stock-outs of commodities for the treatment of snakebite in Kenya.** *PLoS Negl Trop Dis*. 2021. [PMID 34398889](https://pubmed.ncbi.nlm.nih.gov/34398889/)
- Gampini S, et al. **Retrospective study on the incidence of envenomation and accessibility to antivenom in Burkina Faso.** *J Venom Anim Toxins Incl Trop Dis*. 2016. [PMID 26985188](https://pubmed.ncbi.nlm.nih.gov/26985188/)

## 25. Long-term sequelae — the burden mortality misses
## Long-term sequelae — the burden mortality statistics miss

The model's `FIBR` (permanent renal scar) and `NEC` (19-day healing half-life)
states exist because of this literature: the patient who survives is not the same as
the patient who recovers.

- Waiddyanatha S, Silva A, Siribaddana S, Isbister GK. **Long-term Effects of Snake Envenoming.** *Toxins (Basel)*. 2019. [PMID 30935096](https://pubmed.ncbi.nlm.nih.gov/30935096/)
- Jayawardana S, Gnanathasan A, Arambepola C, Chang T. **Long-term health complications following snake envenoming.** *J Multidiscip Healthc*. 2018. [PMID 29983571](https://pubmed.ncbi.nlm.nih.gov/29983571/)

## 26. Methodology — mrgsolve and QSP
## Methodology — mrgsolve and QSP

- Elmokadem A, Riggs MM, Baron KT. **Quantitative Systems Pharmacology and Physiologically-Based Pharmacokinetic Modeling With mrgsolve: A Hands-On Tutorial.** *CPT Pharmacometrics Syst Pharmacol*. 2019. [PMID 31652028](https://pubmed.ncbi.nlm.nih.gov/31652028/)

---

## Verification

Every PMID was obtained by the following procedure.

1. About 50 topic-based queries were run through NCBI `esearch.fcgi` (db=pubmed, sort=relevance)
2. The returned PMIDs were looked up with `esummary.fcgi` to obtain the author, year, journal, and title
3. **Anything not in those results was not put into this file.** There is not one citation written from memory,
   not one guessed PMID, and not one entry of the "this is probably the paper" kind.
4. Each PMID link is of the form `https://pubmed.ncbi.nlm.nih.gov/<PMID>/`, using only the numbers
   confirmed in the lookups above.

What this file does **not** support is set out in the "what it did not take" section above.
The commonest dishonesty in a QSP model is not pretending to evidence that does not exist; it is placing a
calibration parameter as though it were an observation.
