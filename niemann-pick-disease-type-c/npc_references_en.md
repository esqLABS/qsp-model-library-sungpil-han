# Niemann-Pick disease type C (NPC) — References
# Niemann-Pick Disease Type C — Annotated References

The literature underlying every structural assumption and every parameter of this model
(`npc_qsp_model_en.dot` · `npc_mrgsolve_model_en.R` · `npc_reference_model_en.py` · `npc_shiny_app_en.R`).

**Every PMID was queried individually through the NCBI E-utilities and its title, authors, year and journal confirmed.**
Link format: `https://pubmed.ncbi.nlm.nih.gov/<PMID>/`

References used directly as a **quantitative calibration target** for the model are marked 🎯; those
not used in calibration but only for **validation of the model's predictions** are marked 🔍.

---

## Contents

| Section | Topic | Entries |
|----|------|------|
| [1](#1-genetics-and-molecular-basis) | Genetics · molecular basis | 10 |
| [2](#2-the-mechanism-of-lysosomal-cholesterol-egress) | The mechanism of lysosomal cholesterol egress | 13 |
| [3](#3-secondary-lipid-storage--lysosomal-calcium) | Secondary lipid storage · lysosomal calcium | 8 |
| [4](#4-autophagy--mtorc1--tfeb) | Autophagy · mTORC1 · TFEB | 8 |
| [5](#5-neuropathology--neuroinflammation) | Neuropathology · neuroinflammation | 13 |
| [6](#6-biomarkers) | Biomarkers | 14 |
| [7](#7-natural-history--severity-scales) | Natural history · severity scales | 19 |
| [8](#8-miglustat) | Miglustat | 13 |
| [9](#9-arimoclomol) | Arimoclomol | 7 |
| [10](#10-levacetylleucine--n-acetyl-l-leucine) | Levacetylleucine | 11 |
| [11](#11-cyclodextrin-2-hpβcd--adrabetadex) | Cyclodextrin | 12 |
| [12](#12-approaches-in-development--experimental-approaches) | Approaches in development | 8 |
| [13](#13-qsp-methodology) | QSP methodology | 8 |
| | **Total** | **144** |

---

## 1. Genetics and molecular basis

The basis for cluster 2 of the model (`cluster_gene`), the state variables `NPC1_ER`·`NPC1_L`, and the
genotype table `GENOTYPES`.

1. Carstea ED et al. **Niemann-Pick C1 disease gene: homology to mediators of cholesterol homeostasis.** *Science* 1997. — The first cloning of NPC1; a 13-TM membrane protein with homology to a sterol-sensing domain. [PMID 9211849](https://pubmed.ncbi.nlm.nih.gov/9211849/)
2. Naureckiene S et al. **Identification of HE1 as the second gene of Niemann-Pick C disease.** *Science* 2000. — Identification of NPC2 (= HE1). Represented in the model by `f_npc2`. [PMID 11125141](https://pubmed.ncbi.nlm.nih.gov/11125141/)
3. Vanier MT. **Niemann-Pick disease type C.** *Orphanet J Rare Dis* 2010. — The standard review. The basis for the classification by age at onset (perinatal · early infantile · late infantile · juvenile · adult). [PMID 20525256](https://pubmed.ncbi.nlm.nih.gov/20525256/)
4. Wassif CA et al. **High incidence of unrecognized visceral/neurological late-onset Niemann-Pick disease, type C1, predicted by analysis of massively parallel sequencing data sets.** *Genet Med* 2016. — Under-diagnosis of the adult form, and carrier frequency. The basis for the `mild/mild` genotype. [PMID 25764212](https://pubmed.ncbi.nlm.nih.gov/25764212/)
5. Nakasone N et al. **Endoplasmic reticulum-associated degradation of Niemann-Pick C1: evidence for the role of heat shock proteins and identification of lysine residues that accept ubiquitin.** *J Biol Chem* 2014. — 🎯 ERAD of mutant NPC1 and **HSP-dependent rescue of folding**. The direct basis for the model's `kerad`·`Emax_fold`·`theta_eff` structure. [PMID 24891511](https://pubmed.ncbi.nlm.nih.gov/24891511/)
6. Burton BK et al. **Estimating the prevalence of Niemann-Pick disease type C (NPC) in the United States.** *Mol Genet Metab* 2021. [PMID 34304992](https://pubmed.ncbi.nlm.nih.gov/34304992/)
7. Labrecque M et al. **Estimated prevalence of Niemann-Pick type C disease in Quebec.** *Sci Rep* 2021. [PMID 34799641](https://pubmed.ncbi.nlm.nih.gov/34799641/)
8. Geberhiwot T et al. **Consensus clinical management guidelines for Niemann-Pick disease type C.** *Orphanet J Rare Dis* 2018. — The diagnostic and treatment standard. [PMID 29625568](https://pubmed.ncbi.nlm.nih.gov/29625568/)
9. Yoon HJ et al. **The point mutation of the cholesterol trafficking membrane protein NPC1 may affect its proper function.** *Comput Biol Chem* 2022. [PMID 35850050](https://pubmed.ncbi.nlm.nih.gov/35850050/)
10. Elghobashi-Meinhardt N. **Cholesterol Transport in Wild-Type NPC1 and P691S: Molecular Dynamics Simulations Reveal Changes in Dynamical Behavior.** *Int J Mol Sci* 2020. — Differences in residual function between variants. [PMID 32331453](https://pubmed.ncbi.nlm.nih.gov/32331453/)

---

## 2. The mechanism of lysosomal cholesterol egress

The basis for cluster 3 of the model (`cluster_egress`), the state variables `CHOL_V`·`CHOL_C`, and the
Michaelis-Menten egress term `Vmax * f_eg * CHOL/(Km + CHOL)`.

11. Pentchev PG et al. **A defect in cholesterol esterification in Niemann-Pick disease (type C) patients.** *Proc Natl Acad Sci U S A* 1985. — The first biochemical description of the primary defect. [PMID 3865225](https://pubmed.ncbi.nlm.nih.gov/3865225/)
12. Liscum L, Faust JR. **Low density lipoprotein (LDL)-mediated suppression of cholesterol synthesis and LDL uptake is defective in Niemann-Pick type C fibroblasts.** *J Biol Chem* 1987. — 🎯 The direct basis for the SREBP2 paradox (cluster 4), in which **storage overflows while the ER starves**. [PMID 3680287](https://pubmed.ncbi.nlm.nih.gov/3680287/)
13. Slotte JP et al. **Intracellular transport of cholesterol in type C Niemann-Pick fibroblasts.** *Biochim Biophys Acta* 1989. [PMID 2804059](https://pubmed.ncbi.nlm.nih.gov/2804059/)
14. Infante RE et al. **NPC2 facilitates bidirectional transfer of cholesterol between NPC1 and lipid bilayers, a step in cholesterol egress from lysosomes.** *Proc Natl Acad Sci U S A* 2008. — NPC2→NPC1 transfer in series. Represented in the model as `f_eg = f_NPC1 × f_npc2` (a product, hence in series). [PMID 18772377](https://pubmed.ncbi.nlm.nih.gov/18772377/)
15. Wang ML et al. **Identification of surface residues on Niemann-Pick C2 essential for hydrophobic handoff of cholesterol to NPC1 in lysosomes.** *Cell Metab* 2010. — The hydrophobic handoff. [PMID 20674861](https://pubmed.ncbi.nlm.nih.gov/20674861/)
16. Li X et al. **Clues to the mechanism of cholesterol transfer from the structure of NPC1 middle lumenal domain bound to NPC2.** *Proc Natl Acad Sci U S A* 2016. [PMID 27551080](https://pubmed.ncbi.nlm.nih.gov/27551080/)
17. Qian H et al. **Structural Basis of Low-pH-Dependent Lysosomal Cholesterol Egress by NPC1 and NPC2.** *Cell* 2020. — 🎯 The structural basis for egress being **dependent on low pH**. The basis for the feedback in which a rise in lysosomal pH (`Kph`) impairs not only hydrolase activity but egress itself. [PMID 32544384](https://pubmed.ncbi.nlm.nih.gov/32544384/)
18. Pfeffer SR. **NPC intracellular cholesterol transporter 1 (NPC1)-mediated cholesterol export from lysosomes.** *J Biol Chem* 2019. [PMID 30710017](https://pubmed.ncbi.nlm.nih.gov/30710017/)
19. Sandhu J et al. **Aster Proteins Facilitate Nonvesicular Plasma Membrane to ER Cholesterol Transport in Mammalian Cells.** *Cell* 2018. — The `EGRESS_PM` route of cluster 3. [PMID 30220461](https://pubmed.ncbi.nlm.nih.gov/30220461/)
20. Naito T, Saheki Y. **GRAMD1-mediated accessible cholesterol sensing and transport.** *Biochim Biophys Acta Mol Cell Biol Lipids* 2021. [PMID 33932585](https://pubmed.ncbi.nlm.nih.gov/33932585/)
21. Ferrari A et al. **Aster Proteins Regulate the Accessible Cholesterol Pool in the Plasma Membrane.** *Mol Cell Biol* 2020. [PMID 32719109](https://pubmed.ncbi.nlm.nih.gov/32719109/)
22. Long T et al. **Structural basis for itraconazole-mediated NPC1 inhibition.** *Nat Commun* 2020. — Pharmacological inhibition of NPC1 (verification in the reverse direction). [PMID 31919352](https://pubmed.ncbi.nlm.nih.gov/31919352/)
23. Elghobashi-Meinhardt N. **Niemann-Pick type C disease: a QM/MM study of conformational changes in cholesterol in the NPC1(NTD) and NPC2 binding pockets.** *Biochemistry* 2014. [PMID 25251378](https://pubmed.ncbi.nlm.nih.gov/25251378/)

---

## 3. Secondary lipid storage · lysosomal calcium

The basis for clusters 5 and 6 of the model and the state variables `SPH`·`CA_LY`·`GSL_V`·`GSL_C`.

24. Lloyd-Evans E et al. **Niemann-Pick disease type C1 is a sphingosine storage disease that causes deregulation of lysosomal calcium.** *Nat Med* 2008. — 🎯 The key paper, establishing that **sphingosine is the first lipid to accumulate and that it blocks refilling of the acidic Ca²⁺ store**. The basis for the model's `SPH → CA_LY` inhibition term (`Ksph`) and for the structure in which `SPH` has an `f_NPC1`-dependent egress. [PMID 18953351](https://pubmed.ncbi.nlm.nih.gov/18953351/)
25. Lloyd-Evans E, Platt FM. **Lysosomal Ca²⁺ homeostasis: role in pathogenesis of lysosomal storage diseases.** *Cell Calcium* 2011. [PMID 21724254](https://pubmed.ncbi.nlm.nih.gov/21724254/)
26. Zervas M et al. **Critical role for glycosphingolipids in Niemann-Pick disease type C.** *Curr Biol* 2001. — Reducing GSL synthesis improves the mouse phenotype → the rationale for miglustat. [PMID 11525744](https://pubmed.ncbi.nlm.nih.gov/11525744/)
27. Vanier MT. **Lipid changes in Niemann-Pick disease type C brain: personal experience and review of the literature.** *Neurochem Res* 1999. — 🔍 The fold elevation of GM2/GM3 in the NPC brain. The basis for validating the model's CNS GSL accumulation (model prediction ~2.6-fold). [PMID 10227680](https://pubmed.ncbi.nlm.nih.gov/10227680/)
28. Shen D et al. **Lipid storage disorders block lysosomal trafficking by inhibiting a TRP channel and lysosomal calcium release.** *Nat Commun* 2012. — Inhibition of TRPML1 by cholesterol (`Ktr_chol`). [PMID 22415822](https://pubmed.ncbi.nlm.nih.gov/22415822/)
29. Pagano RE. **Endocytic trafficking of glycosphingolipids in sphingolipid storage diseases.** *Philos Trans R Soc Lond B Biol Sci* 2003. [PMID 12803922](https://pubmed.ncbi.nlm.nih.gov/12803922/)
30. Walkley SU et al. **Initiation and growth of ectopic neurites and meganeurites during postnatal cortical development in ganglioside storage disease.** *Brain Res Dev Brain Res* 1990. — GM2-driven meganeurites (cluster 11, `MEGANEURITE`). [PMID 2108821](https://pubmed.ncbi.nlm.nih.gov/2108821/)
31. Kuech EM et al. **Alterations in membrane trafficking and pathophysiological implications in lysosomal storage disorders.** *Biochimie* 2016. [PMID 27664461](https://pubmed.ncbi.nlm.nih.gov/27664461/)

---

## 4. Autophagy · mTORC1 · TFEB

The basis for cluster 7 of the model and the state variables `AUTOPH`·`HYD`.

32. Elrick MJ et al. **Impaired proteolysis underlies autophagic dysfunction in Niemann-Pick type C disease.** *Hum Mol Genet* 2012. — 🎯 That autophagic **induction is normal and it is degradation that fails**. The basis on which the model writes `AUTOPH` with a constant influx and a reduced degradation. [PMID 22872701](https://pubmed.ncbi.nlm.nih.gov/22872701/)
33. Settembre C et al. **TFEB links autophagy to lysosomal biogenesis.** *Science* 2011. [PMID 21617040](https://pubmed.ncbi.nlm.nih.gov/21617040/)
34. Castellano BM et al. **Lysosomal cholesterol activates mTORC1 via an SLC38A9-Niemann-Pick C1 signaling complex.** *Science* 2017. — Lysosomal cholesterol → mTORC1. [PMID 28336668](https://pubmed.ncbi.nlm.nih.gov/28336668/)
35. Lim CY et al. **ER-lysosome contacts enable cholesterol sensing by mTORC1 and drive aberrant growth signalling in Niemann-Pick type C.** *Nat Cell Biol* 2019. [PMID 31548609](https://pubmed.ncbi.nlm.nih.gov/31548609/)
36. Davis OB et al. **NPC1-mTORC1 Signaling Couples Cholesterol Sensing to Organelle Homeostasis and Is a Targetable Pathway in Niemann-Pick Type C.** *Dev Cell* 2021. [PMID 33308480](https://pubmed.ncbi.nlm.nih.gov/33308480/)
37. Kataura T et al. **Targeting the autophagy-NAD axis protects against cell death in Niemann-Pick type C1 disease models.** *Cell Death Dis* 2024. [PMID 38821960](https://pubmed.ncbi.nlm.nih.gov/38821960/)
38. Dai S et al. **Methyl-β-cyclodextrin restores impaired autophagy flux in Niemann-Pick C1-deficient cells through activation of AMPK.** *Autophagy* 2017. [PMID 28613987](https://pubmed.ncbi.nlm.nih.gov/28613987/)
39. Lee H et al. **Pathological roles of the VEGF/SphK pathway in Niemann-Pick type C neurons.** *Nat Commun* 2014. [PMID 25417698](https://pubmed.ncbi.nlm.nih.gov/25417698/)

---

## 5. Neuropathology · neuroinflammation

The basis for clusters 11, 12 and 13 of the model and the state variables `PC`·`PC_S`·`PC_LOST`·`INFL`·`SYN`·`CBL`.

40. Higashi Y et al. **Cerebellar degeneration in the Niemann-Pick type C mouse.** *Acta Neuropathol* 1993. — The first description of Purkinje cell loss. [PMID 8382896](https://pubmed.ncbi.nlm.nih.gov/8382896/)
41. Sarna JR et al. **Patterned Purkinje cell degeneration in mouse models of Niemann-Pick type C disease.** *J Comp Neurol* 2003. — The anterior→posterior lobe gradient. [PMID 12528192](https://pubmed.ncbi.nlm.nih.gov/12528192/)
42. German DC et al. **Neurodegeneration in the Niemann-Pick C mouse: glial involvement.** *Neuroscience* 2002. [PMID 11823057](https://pubmed.ncbi.nlm.nih.gov/11823057/)
43. Elrick MJ et al. **Conditional Niemann-Pick C mice demonstrate cell autonomous Purkinje cell neurodegeneration.** *Hum Mol Genet* 2010. — 🎯 **Cell-autonomous death**. The basis on which the model treats `CHOL_CNS → PC_S` as a cell-autonomous route and uses neuroinflammation only as an *amplifier*. [PMID 20007718](https://pubmed.ncbi.nlm.nih.gov/20007718/)
44. Dinkel L et al. **Myeloid cell-specific loss of NPC1 in mice recapitulates microgliosis and neurodegeneration in patients.** *Sci Transl Med* 2024. — Loss confined to microglia also produces neurodegeneration → the basis for the `INFL → PC_S` amplification loop. [PMID 39630885](https://pubmed.ncbi.nlm.nih.gov/39630885/)
45. Love S et al. **Neurofibrillary tangles in Niemann-Pick disease type C.** *Brain* 1995. — NFTs of the Alzheimer type. [PMID 7894998](https://pubmed.ncbi.nlm.nih.gov/7894998/)
46. Bu B et al. **Niemann-Pick disease type C yields possible clue for why cerebellar neurons do not form neurofibrillary tangles.** *Neurobiol Dis* 2002. [PMID 12505421](https://pubmed.ncbi.nlm.nih.gov/12505421/)
47. Malnar M et al. **Bidirectional links between Alzheimer's disease and Niemann-Pick type C disease.** *Neurobiol Dis* 2014. [PMID 24907492](https://pubmed.ncbi.nlm.nih.gov/24907492/)
48. Mattsson N et al. **Gamma-secretase-dependent amyloid-beta is increased in Niemann-Pick type C: a cross-sectional study.** *Neurology* 2011. [PMID 21205675](https://pubmed.ncbi.nlm.nih.gov/21205675/)
49. Woś M et al. **Mitochondrial dysfunction in fibroblasts derived from patients with Niemann-Pick type C disease.** *Arch Biochem Biophys* 2016. — The model's `MITO`·`ROS`. [PMID 26869201](https://pubmed.ncbi.nlm.nih.gov/26869201/)
50. Takikita S et al. **Perturbed myelination process of premyelinating oligodendrocyte in Niemann-Pick type C mouse.** *J Neuropathol Exp Neurol* 2004. — 🎯 The basis for the developmental vulnerability (`v_dev`·`tau_dev`): the same biochemical insult exacts a greater cost during the period of active myelination. [PMID 15217094](https://pubmed.ncbi.nlm.nih.gov/15217094/)
51. Burbulla LF et al. **Modeling Brain Pathology of Niemann-Pick Disease Type C Using Patient-Derived Neurons.** *Mov Disord* 2021. [PMID 33438272](https://pubmed.ncbi.nlm.nih.gov/33438272/)
52. Wheeler S, Sillence DJ. **Niemann-Pick type C disease: cellular pathology and pharmacotherapy.** *J Neurochem* 2020. — A review of the cellular pathology and the pharmacotherapy. [PMID 31608980](https://pubmed.ncbi.nlm.nih.gov/31608980/)

---

## 6. Biomarkers

The basis for cluster 9 of the model and the state variables `TRIOL`·`PPCS`·`TCG`·`NFL`·`CALB`.
**The model's first structural claim (the compartment claim) hangs on this section.**

53. Porter FD et al. **Cholesterol oxidation products are sensitive and specific blood-based biomarkers for Niemann-Pick C1 disease.** *Sci Transl Med* 2010. — The first establishment of 7-KC · C-triol. [PMID 21048217](https://pubmed.ncbi.nlm.nih.gov/21048217/)
54. Jiang X et al. **A sensitive and specific LC-MS/MS method for rapid diagnosis of Niemann-Pick C1 disease from human plasma.** *J Lipid Res* 2011. [PMID 21518695](https://pubmed.ncbi.nlm.nih.gov/21518695/)
55. Kuchar L et al. **Quantitation of plasmatic lysosphingomyelin and lysosphingomyelin-509 for differential screening of Niemann-Pick A/B and C diseases.** *Anal Biochem* 2017. [PMID 28259515](https://pubmed.ncbi.nlm.nih.gov/28259515/)
56. Sidhu R et al. **N-acyl-O-phosphocholineserines: structures of a novel class of lipids that are biomarkers for Niemann-Pick C1 disease.** *J Lipid Res* 2019. — The actual structure of lysoSM-509 is PPCS. The model's state variable `PPCS`. [PMID 31201291](https://pubmed.ncbi.nlm.nih.gov/31201291/)
57. Jiang X et al. **Development of a bile acid-based newborn screen for Niemann-Pick disease type C.** *Sci Transl Med* 2016. — Bile acid B (TCG). The model's state variable `TCG`. [PMID 27147587](https://pubmed.ncbi.nlm.nih.gov/27147587/)
58. Jiang X et al. **Diagnosis of Niemann-Pick C1 by measurement of bile acid biomarkers in archived newborn dried blood spots.** *Mol Genet Metab* 2019. — Already elevated from the neonatal period → the basis on which the model has the marker rise immediately after birth. [PMID 30172462](https://pubmed.ncbi.nlm.nih.gov/30172462/)
59. Mazzacuva F et al. **Identification of novel bile acids as biomarkers for the early diagnosis of Niemann-Pick C disease.** *FEBS Lett* 2016. [PMID 27139891](https://pubmed.ncbi.nlm.nih.gov/27139891/)
60. Bradbury A et al. **Cerebrospinal Fluid Calbindin D Concentration as a Biomarker of Cerebellar Disease Progression in Niemann-Pick Type C1 Disease.** *J Pharmacol Exp Ther* 2016. — 🎯 CSF calbindin is a read-out of the **Purkinje death flux**. The basis on which the model makes `CALB` proportional to the rate of death (`die`) — the decisive point being that it is the *rate of death*, not the amount stored. [PMID 27307499](https://pubmed.ncbi.nlm.nih.gov/27307499/)
61. Alam MS et al. **Plasma signature of neurological disease in the monogenetic disorder Niemann-Pick Type C.** *J Biol Chem* 2014. [PMID 24488491](https://pubmed.ncbi.nlm.nih.gov/24488491/)
62. Sidhu R et al. **A validated LC-MS/MS assay for quantification of 24(S)-hydroxycholesterol in plasma and cerebrospinal fluid.** *J Lipid Res* 2015. — A brain-derived oxysterol (the model's `OHC24` node). [PMID 25866316](https://pubmed.ncbi.nlm.nih.gov/25866316/)
63. Agrawal N et al. **Neurofilament light chain in cerebrospinal fluid as a novel biomarker in evaluating both clinical severity and therapeutic response in Niemann-Pick disease type C1.** *Genet Med* 2023. — The model's state variable `NFL`. [PMID 36470574](https://pubmed.ncbi.nlm.nih.gov/36470574/)
64. Eratne D et al. **Plasma neurofilament light chain is increased in Niemann-Pick Type C but glial fibrillary acidic protein is not.** *Acta Neuropsychiatr* 2024. [PMID 38533577](https://pubmed.ncbi.nlm.nih.gov/38533577/)
65. Stern S et al. **Evaluation of the landscape of pharmacodynamic biomarkers in Niemann-Pick Disease Type C (NPC).** *Orphanet J Rare Dis* 2024. — 🔍 A survey of the present limits of the marker-to-clinic link. Corresponds directly to the model's compartment claim. [PMID 39061081](https://pubmed.ncbi.nlm.nih.gov/39061081/)
66. Pataj Z et al. **Quantification of oxysterols in human plasma and red blood cells by liquid chromatography high-resolution tandem mass spectrometry.** *J Chromatogr A* 2016. [PMID 26607314](https://pubmed.ncbi.nlm.nih.gov/26607314/)

---

## 7. Natural history · severity scales

The basis for the model's clinical scale mappings (`SARA`·`NPCCSS5`·`NPCCSS4`·`NPCCSS17`) and for
**the second structural claim (the reserve claim)**.

67. Yanjanin NM et al. **Linear clinical progression, independent of age of onset, in Niemann-Pick disease, type C.** *Am J Med Genet B Neuropsychiatr Genet* 2010. — 🎯🔍 **The paper that determined the model's structure.** Because progression is (a) linear and (b) independent of age at onset, damage was written as an *integral* and the rate of death as a *saturating gate*. The linearity was used for calibration, the independence of age at onset for validation. [PMID 19415691](https://pubmed.ncbi.nlm.nih.gov/19415691/)
68. Patterson MC et al. **Validation of the 5-domain Niemann-Pick type C Clinical Severity Scale.** *Orphanet J Rare Dis* 2021. — The 5-domain scale. [PMID 33579322](https://pubmed.ncbi.nlm.nih.gov/33579322/)
69. Mengel E et al. **Clinical disease progression and biomarkers in Niemann-Pick disease type C: a prospective cohort study.** *Orphanet J Rare Dis* 2020. — 🎯 **The source of four of the model's quantitative targets**: 1.5 points/year on the 5-domain scale, ~2.9 points/year on the 17-domain scale, plasma triol 88.31 in patients vs 5.97 ng/mL in controls, and a triol–5-domain Spearman ρ = 0.265. [PMID 33228797](https://pubmed.ncbi.nlm.nih.gov/33228797/)
70. Mengel E et al. **Correction to: Clinical disease progression and biomarkers in Niemann-Pick disease type C.** *Orphanet J Rare Dis* 2021. [PMID 34074315](https://pubmed.ncbi.nlm.nih.gov/34074315/)
71. Cortina-Borja M et al. **Annual severity increment score as a tool for stratifying patients with Niemann-Pick disease type C and for recruitment to clinical trials.** *Orphanet J Rare Dis* 2018. [PMID 30115089](https://pubmed.ncbi.nlm.nih.gov/30115089/)
72. Imrie J et al. **The natural history of Niemann-Pick disease type C in the UK.** *J Inherit Metab Dis* 2007. [PMID 17160617](https://pubmed.ncbi.nlm.nih.gov/17160617/)
73. Mengel E et al. **Niemann-Pick disease type C symptomatology: an expert-based clinical description.** *Orphanet J Rare Dis* 2013. [PMID 24135395](https://pubmed.ncbi.nlm.nih.gov/24135395/)
74. Stampfer M et al. **Niemann-Pick disease type C clinical database: cognitive and coordination deficits are early disease indicators.** *Orphanet J Rare Dis* 2013. [PMID 23433426](https://pubmed.ncbi.nlm.nih.gov/23433426/)
75. Walterfang M et al. **Dysphagia as a risk factor for mortality in Niemann-Pick disease type C: systematic literature review and evidence from studies with miglustat.** *Orphanet J Rare Dis* 2012. — 🎯 The basis for making the model's survival hazard proportional to the **square of swallowing function** (`h_swal`). [PMID 23039766](https://pubmed.ncbi.nlm.nih.gov/23039766/)
76. Bianconi SE et al. **Evaluation of age of death in Niemann-Pick disease, type C: Utility of disease support group websites to understand natural history.** *Mol Genet Metab* 2019. — 🔍 The distribution of age at death. [PMID 30850267](https://pubmed.ncbi.nlm.nih.gov/30850267/)
77. Gardin A et al. **A Retrospective Multicentric Study of 34 Patients with Niemann-Pick Type C Disease and Early Liver Involvement.** *J Pediatr* 2023. — The model's `h_liver` and the perinatal form. [PMID 36265573](https://pubmed.ncbi.nlm.nih.gov/36265573/)
78. Schmitz-Hübsch T et al. **Scale for the assessment and rating of ataxia: development of a new clinical scale.** *Neurology* 2006. — The original SARA paper (0-40 points). [PMID 16769946](https://pubmed.ncbi.nlm.nih.gov/16769946/)
79. Bremova-Ertl T et al. **A cross-sectional, prospective ocular motor study in 72 patients with Niemann-Pick disease type C.** *Eur J Neurol* 2021. — Quantification of saccade velocity. [PMID 34096670](https://pubmed.ncbi.nlm.nih.gov/34096670/)
80. Hopf S et al. **Vertical saccadic palsy and foveal retinal thinning in Niemann-Pick disease type C.** *PLoS One* 2021. [PMID 34086834](https://pubmed.ncbi.nlm.nih.gov/34086834/)
81. Grillini A et al. **Measuring saccades in patients with Niemann-Pick type C: A comparison between video-oculography and a novel method.** *Clin Park Relat Disord* 2022. [PMID 36338825](https://pubmed.ncbi.nlm.nih.gov/36338825/)
82. King KA et al. **Auditory phenotype of Niemann-Pick disease, type C1.** *Ear Hear* 2014. — 🎯 The baseline hearing loss of the disease itself. Required in order to **separate** it from cyclodextrin ototoxicity. [PMID 24225652](https://pubmed.ncbi.nlm.nih.gov/24225652/)
83. King KA et al. **Hearing loss is an early consequence of Npc1 gene deletion in the mouse model of Niemann-Pick disease, type C.** *J Assoc Res Otolaryngol* 2014. [PMID 24839095](https://pubmed.ncbi.nlm.nih.gov/24839095/)
84. Ong LT et al. **Psychosis symptoms associated with Niemann-Pick disease type C.** *Psychiatr Genet* 2021. — Psychiatric symptoms as the first presentation of the adult form. [PMID 34133410](https://pubmed.ncbi.nlm.nih.gov/34133410/)
85. Patterson MC et al. **Recommendations for the diagnosis and management of Niemann-Pick disease type C: an update.** *Mol Genet Metab* 2012. [PMID 22572546](https://pubmed.ncbi.nlm.nih.gov/22572546/)

---

## 8. Miglustat

The basis for cluster 16 of the model and the `mig_*` parameters.

86. Patterson MC et al. **Miglustat for treatment of Niemann-Pick C disease: a randomised controlled study.** *Lancet Neurol* 2007. — 🎯 The primary endpoint was horizontal saccadic eye movement velocity (HSEM). Improved swallowing, stable hearing, delayed deterioration of gait. Corresponds to the model's prediction of the miglustat effect as a **visceral > CNS asymmetry**. [PMID 17689147](https://pubmed.ncbi.nlm.nih.gov/17689147/)
87. Wraith JE et al. **Miglustat in adult and juvenile patients with Niemann-Pick disease type C: long-term data from a clinical trial.** *Mol Genet Metab* 2010. [PMID 20045366](https://pubmed.ncbi.nlm.nih.gov/20045366/)
88. Patterson MC et al. **Long-term miglustat therapy in children with Niemann-Pick disease type C.** *J Child Neurol* 2010. [PMID 19822772](https://pubmed.ncbi.nlm.nih.gov/19822772/)
89. Patterson MC et al. **Stable or improved neurological manifestations during miglustat therapy in patients from the international disease registry for Niemann-Pick disease type C.** *Orphanet J Rare Dis* 2015. [PMID 26017010](https://pubmed.ncbi.nlm.nih.gov/26017010/)
90. Patterson MC et al. **Treatment outcomes following continuous miglustat therapy in patients with Niemann-Pick disease Type C.** *Orphanet J Rare Dis* 2020. [PMID 32334605](https://pubmed.ncbi.nlm.nih.gov/32334605/)
91. Pineda M et al. **Miglustat in Niemann-Pick disease type C patients: a review.** *Orphanet J Rare Dis* 2018. [PMID 30111334](https://pubmed.ncbi.nlm.nih.gov/30111334/)
92. Solomon BI et al. **Association of Miglustat With Swallowing Outcomes in Niemann-Pick Disease, Type C1.** *JAMA Neurol* 2020. — 🔍 The most quantitative evidence on swallowing outcomes. Validation of the model's `SWALLOW` route. [PMID 32897301](https://pubmed.ncbi.nlm.nih.gov/32897301/)
93. Platt FM et al. **Prevention of lysosomal storage in Tay-Sachs mice treated with N-butyldeoxynojirimycin.** *Science* 1997. — The principle of substrate reduction therapy. [PMID 9103204](https://pubmed.ncbi.nlm.nih.gov/9103204/)
94. Butters TD et al. **Imino sugar inhibitors for treating the lysosomal glycosphingolipidoses.** *Glycobiology* 2005. — 🎯 The range of IC₅₀ for GCS inhibition (`mig_IC50`). [PMID 15901676](https://pubmed.ncbi.nlm.nih.gov/15901676/)
95. Andersson U et al. **N-butyldeoxygalactonojirimycin: a more selective inhibitor of glycosphingolipid biosynthesis than N-butyldeoxynojirimycin, in vitro and in vivo.** *Biochem Pharmacol* 2000. — The mechanistic basis for the inhibition of the intestinal disaccharidases (diarrhoea). [PMID 10718340](https://pubmed.ncbi.nlm.nih.gov/10718340/)
96. Shayman JA. **The development and use of small molecule inhibitors of glycosphingolipid metabolism for lysosomal storage diseases.** *J Lipid Res* 2014. [PMID 24534703](https://pubmed.ncbi.nlm.nih.gov/24534703/)
97. Maegawa GH et al. **Pharmacokinetics, safety and tolerability of miglustat in the treatment of pediatric patients with GM2 gangliosidosis.** *Mol Genet Metab* 2009. — 🎯 Paediatric PK. The model's `mig_V`·`mig_CL` and its body-surface-area-based dosing. [PMID 19447653](https://pubmed.ncbi.nlm.nih.gov/19447653/)
98. Belmatoug N et al. **Gastrointestinal disturbances and their management in miglustat-treated patients.** *J Inherit Metab Dis* 2011. — The model's `MIG_GI_AE` → dose-reduction route. [PMID 21779792](https://pubmed.ncbi.nlm.nih.gov/21779792/)

---

## 9. Arimoclomol

The basis for cluster 17 of the model and the `ari_*`·`Emax_fold` parameters.

99. Mengel E et al. **Efficacy and safety of arimoclomol in Niemann-Pick disease type C: Results from a double-blind, randomised, placebo-controlled, multinational phase 2/3 trial of a novel treatment.** *J Inherit Metab Dis* 2021. — 🎯 A 12-month difference on the 5-domain NPCCSS of −1.40 (95% CI −2.76, −0.03; p = 0.046); progression in the placebo arm 2.15; **subgroup on concomitant miglustat −2.06 (p = 0.006)**. The calibration target for the model's `Emax_fold`, and the test of its prediction of a multiplicative interaction. [PMID 34418116](https://pubmed.ncbi.nlm.nih.gov/34418116/)
100. Mengel E et al. **Efficacy results from a 12-month double-blind randomized trial of arimoclomol for treatment of Niemann-Pick disease type C (NPC): Presenting a rescored 4-domain NPC Clinical Severity Scale.** *Mol Genet Metab Rep* 2025. — 🔍 An R4DNPCCSS difference of −1.70 (95% CI −3.05, −0.34; p = 0.016). The reasons for excluding the cognitive domain and rescoring the swallowing domain. Validation of the model's `n4_frac`·`n4_gain` predictions. [PMID 40520915](https://pubmed.ncbi.nlm.nih.gov/40520915/)
101. Kirkegaard T et al. **Hsp70 stabilizes lysosomes and reverts Niemann-Pick disease-associated lysosomal pathology.** *Nature* 2010. — 🎯 Stabilisation of the lysosome by HSP70 (`e_hyd_hsp`). [PMID 20111001](https://pubmed.ncbi.nlm.nih.gov/20111001/)
102. Petersen NH, Kirkegaard T. **Connecting Hsp70, sphingolipid metabolism and lysosomal stability.** *Cell Cycle* 2010. [PMID 20519957](https://pubmed.ncbi.nlm.nih.gov/20519957/)
103. Gray J et al. **Heat shock protein amplification improves cerebellar myelination in the Npc1^nih mouse model.** *EBioMedicine* 2022. [PMID 36455410](https://pubmed.ncbi.nlm.nih.gov/36455410/)
104. Keam SJ. **Arimoclomol: First Approval.** *Drugs* 2025. — 🎯 The approved doses (47/62/93/124 mg tid by body weight) and the PK. [PMID 39715913](https://pubmed.ncbi.nlm.nih.gov/39715913/)
105. Cudkowicz ME et al. **Arimoclomol at dosages up to 300 mg/day is well tolerated and safe in amyotrophic lateral sclerosis.** *Muscle Nerve* 2008. — 🎯 Linear PK, t½ ~4 h, and a **dose-dependent rise in CSF concentration**. The model's `ari_Kp_csf`. [PMID 18551622](https://pubmed.ncbi.nlm.nih.gov/18551622/)

---

## 10. Levacetylleucine / N-acetyl-L-leucine

The basis for cluster 18 of the model and the `nal_*` parameters.

106. Bremova-Ertl T et al. **Trial of N-Acetyl-l-Leucine in Niemann-Pick Disease Type C.** *N Engl J Med* 2024. — 🎯 IB1001-301. A 12-week change in SARA of −1.97 ± 2.43 (drug) vs −0.60 ± 2.39 (placebo), LS mean difference **−1.28 (95% CI −1.91, −0.65; p < 0.001)**; baseline SARA 15.91 ± 7.65. The calibration target for the model's `nal_Emax_sym`. [PMID 38294974](https://pubmed.ncbi.nlm.nih.gov/38294974/)
107. Bremova-Ertl T et al. **Efficacy and safety of N-acetyl-L-leucine in Niemann-Pick disease type C.** *J Neurol* 2022. [PMID 34387740](https://pubmed.ncbi.nlm.nih.gov/34387740/)
108. Patterson MC et al. **Disease-Modifying, Neuroprotective Effect of N-Acetyl-l-Leucine in Adult and Pediatric Patients With Niemann-Pick Disease Type C.** *Neurology* 2025. — 🔍 The long-term extension. Change in SARA from baseline of −1.88 ± 2.89 at 12 months and −1.64 ± 3.24 at 18 months. **This is an item the model does not reproduce**: fitting `nal_Emax_dm` to these values would require blocking 247% of the CNS lipid influx (physically impossible), so they were not used for calibration; the 8% supported by the mouse data was used instead and the discrepancy reported as it stands. See section 7 ② of the README. [PMID 40513057](https://pubmed.ncbi.nlm.nih.gov/40513057/)
109. Fields T et al. **N-acetyl-L-leucine for Niemann-Pick type C: a multinational double-blind randomized placebo-controlled crossover study.** *Trials* 2023. — 🎯 **A crossover design, a 12-week period, a washout**. The direct basis for the model's third structural claim (the design claim). [PMID 37248494](https://pubmed.ncbi.nlm.nih.gov/37248494/)
110. Fields T et al. **A master protocol to investigate a novel therapy acetyl-L-leucine for three ultra-rare neurodegenerative diseases.** *Trials* 2021. [PMID 33482890](https://pubmed.ncbi.nlm.nih.gov/33482890/)
111. Bremova T et al. **Acetyl-dl-leucine in Niemann-Pick type C: A case series.** *Neurology* 2015. — 🎯 **Onset of effect within days to weeks, and deterioration again on withdrawal**. The basis on which the model sets the effect-site `ke0` at 0.15/d (t½ ~4.6 days) and hangs the effect on `D_rev`. [PMID 26400580](https://pubmed.ncbi.nlm.nih.gov/26400580/)
112. Kaya E et al. **Acetyl-leucine slows disease progression in lysosomal storage disorders.** *Brain Commun* 2021. — 🎯 Slowed progression in the mouse → the basis for the existence of `nal_Emax_dm` (the disease-modifying component). [PMID 33738443](https://pubmed.ncbi.nlm.nih.gov/33738443/)
113. Günther L et al. **N-acetyl-L-leucine accelerates vestibular compensation after unilateral labyrinthectomy by action in the cerebellum and thalamus.** *PLoS One* 2015. — The symptomatic mechanism (`NAL_VEST`). [PMID 25803613](https://pubmed.ncbi.nlm.nih.gov/25803613/)
114. Vibert N, Vidal PP. **In vitro effects of acetyl-DL-leucine (tanganil) on central vestibular neurons and vestibulo-ocular networks of the guinea-pig.** *Eur J Neurosci* 2001. [PMID 11207808](https://pubmed.ncbi.nlm.nih.gov/11207808/)
115. Schniepp R et al. **Acetyl-DL-leucine improves gait variability in patients with cerebellar ataxia - a case series.** *Cerebellum Ataxias* 2016. [PMID 27073690](https://pubmed.ncbi.nlm.nih.gov/27073690/)
116. Martakis K et al. **Safety and efficacy of levacetylleucine in ataxia-telangiectasia: a phase 3, randomised, double-blind trial.** *Lancet Neurol* 2026. — 🔍 The same design and the same size of effect in a different disease. External validation of the symptomatic-mechanism interpretation. [PMID 42309084](https://pubmed.ncbi.nlm.nih.gov/42309084/)

---

## 11. Cyclodextrin (2-HPβCD / adrabetadex)

The basis for cluster 19 of the model and the `cd_*` parameters.
**The evidence for the mechanistically inseparable toxicity (efficacy and ototoxicity being one and the same mechanism) is in this section.**

117. Ory DS et al. **Intrathecal 2-hydroxypropyl-β-cyclodextrin decreases neurological disease progression in Niemann-Pick disease, type C1: a non-randomised, open-label, phase 1-2 trial.** *Lancet* 2017. — 🎯 Open label. Progression 21/21 (historical control) vs 7/14 (treated). **Scored as "NSS minus hearing"** — that hearing had to be taken out of the endpoint because the drug damages hearing is the direct basis for the model's toxicity term. [PMID 28803710](https://pubmed.ncbi.nlm.nih.gov/28803710/)
118. Liu B et al. **Reversal of defective lysosomal transport in NPC disease ameliorates liver dysfunction and neurodegeneration in the npc1-/- mouse.** *Proc Natl Acad Sci U S A* 2009. — Storage is reduced even by a single dose. [PMID 19171898](https://pubmed.ncbi.nlm.nih.gov/19171898/)
119. Vite CH et al. **Intracisternal cyclodextrin prevents cerebellar dysfunction and Purkinje cell death in feline Niemann-Pick type C1 disease.** *Sci Transl Med* 2015. — 🎯 Prevention of Purkinje death in the feline model. The basis for the term in which cyclodextrin lowers `CHOL_C` directly and independently of NPC1. [PMID 25717099](https://pubmed.ncbi.nlm.nih.gov/25717099/)
120. Vance JE, Peake KB. **Function of the Niemann-Pick type C proteins and their bypass by cyclodextrin.** *Curr Opin Lipidol* 2011. — 🎯 **Bypassing NPC1**. The basis for the model's genotype-independent prediction (the effect is retained even in null/null and in NPC2). [PMID 21412152](https://pubmed.ncbi.nlm.nih.gov/21412152/)
121. Crumling MA et al. **Hearing loss and hair cell death in mice given the cholesterol-chelating agent hydroxypropyl-β-cyclodextrin.** *PLoS One* 2012. — 🎯 Death of the outer hair cells. The model's `cd_koto`·`OHC`. [PMID 23285273](https://pubmed.ncbi.nlm.nih.gov/23285273/)
122. Crumling MA et al. **Cyclodextrins and Iatrogenic Hearing Loss: New Drugs with Significant Risk.** *Front Cell Neurosci* 2017. [PMID 29163061](https://pubmed.ncbi.nlm.nih.gov/29163061/)
123. Ward S et al. **2-hydroxypropyl-beta-cyclodextrin raises hearing threshold in normal cats and in cats with Niemann-Pick type C disease.** *Pediatr Res* 2010. — 🎯 It occurs in normal animals as well → a drug toxicity unrelated to the disease. [PMID 20357695](https://pubmed.ncbi.nlm.nih.gov/20357695/)
124. Cronin S et al. **Hearing Loss and Otopathology Following Systemic and Intracerebroventricular Delivery of 2-Hydroxypropyl-Beta-Cyclodextrin.** *J Assoc Res Otolaryngol* 2015. — 🎯 **Differences by route of administration**. The basis on which the model hangs the ototoxicity on CSF concentration (by way of the cochlear aqueduct) and separates it from the intravenous route. [PMID 26055150](https://pubmed.ncbi.nlm.nih.gov/26055150/)
125. Takahashi S et al. **Susceptibility of outer hair cells to cholesterol chelator 2-hydroxypropyl-β-cyclodextrine is prestin-dependent.** *Sci Rep* 2016. — 🎯 Prestin dependence → **efficacy (extraction of membrane cholesterol) and toxicity are the same mechanism**. [PMID 26903308](https://pubmed.ncbi.nlm.nih.gov/26903308/)
126. Zhou Y et al. **The susceptibility of cochlear outer hair cells to cyclodextrin is not related to their electromotile activity.** *Acta Neuropathol Commun* 2018. — 🔍 A result contradicting the paper above. The model states explicitly that the mechanistic detail is unresolved. [PMID 30249300](https://pubmed.ncbi.nlm.nih.gov/30249300/)
127. Matsuo M et al. **Effects of intracerebroventricular administration of 2-hydroxypropyl-β-cyclodextrin in a patient with Niemann-Pick Type C disease.** *Mol Genet Metab Rep* 2014. [PMID 27896112](https://pubmed.ncbi.nlm.nih.gov/27896112/)
128. El-Darzi N et al. **2-Hydroxypropyl-β-cyclodextrin reduces retinal cholesterol in wild-type and Cyp27a1−/−Cyp46a1−/− mice.** *Br J Pharmacol* 2021. — It extracts cholesterol from off-target tissues as well. [PMID 32698250](https://pubmed.ncbi.nlm.nih.gov/32698250/)

---

## 12. Approaches in development · experimental approaches

The basis for cluster 20 of the model.

129. Chandler RJ et al. **Systemic AAV9 gene therapy improves the lifespan of mice with Niemann-Pick disease, type C1.** *Hum Mol Genet* 2017. [PMID 27798114](https://pubmed.ncbi.nlm.nih.gov/27798114/)
130. Pipalia NH et al. **Histone deacetylase inhibitor treatment dramatically reduces cholesterol accumulation in Niemann-Pick type C1 mutant human fibroblasts.** *Proc Natl Acad Sci U S A* 2011. [PMID 21436030](https://pubmed.ncbi.nlm.nih.gov/21436030/)
131. Pipalia NH et al. **Histone deacetylase inhibitors correct the cholesterol storage defect in most Niemann-Pick C1 mutant cells.** *J Lipid Res* 2017. — 🔍 **Differences in response between variants** → support for the model's genotype × mechanism interaction. [PMID 28193631](https://pubmed.ncbi.nlm.nih.gov/28193631/)
132. Pipalia NH et al. **HSP90 inhibitors reduce cholesterol storage in Niemann-Pick type C1 mutant fibroblasts.** *J Lipid Res* 2021. — An alternative approach to the chaperone route. [PMID 34481829](https://pubmed.ncbi.nlm.nih.gov/34481829/)
133. Maceyka M, Spiegel S. **The potential of histone deacetylase inhibitors in Niemann-Pick type C disease.** *FEBS J* 2013. [PMID 23992240](https://pubmed.ncbi.nlm.nih.gov/23992240/)
134. Hovakimyan M et al. **Combined therapy with cyclodextrin/allopregnanolone and miglustat improves motor but not cognitive functions in Niemann-Pick Type C1 mice.** *Neuroscience* 2013. — 🔍 Whether the combination is additive or synergistic. [PMID 23948640](https://pubmed.ncbi.nlm.nih.gov/23948640/)
135. Maass F et al. **Reduced cerebellar neurodegeneration after combined therapy with cyclodextrin/allopregnanolone and miglustat in NPC1 mutant mice.** *J Neurosci Res* 2015. [PMID 25400034](https://pubmed.ncbi.nlm.nih.gov/25400034/)
136. Ebner L et al. **Evaluation of Two Liver Treatment Strategies in a Mouse Model of Niemann-Pick-Disease Type C1.** *Int J Mol Sci* 2018. [PMID 29587349](https://pubmed.ncbi.nlm.nih.gov/29587349/)

---

## 13. QSP methodology

137. Elmokadem A et al. **Quantitative Systems Pharmacology and Physiologically-Based Pharmacokinetic Modeling With mrgsolve: A Hands-On Tutorial.** *CPT Pharmacometrics Syst Pharmacol* 2019. — The standard citation for mrgsolve. [PMID 31652028](https://pubmed.ncbi.nlm.nih.gov/31652028/)
138. Helmlinger G et al. **Quantitative Systems Pharmacology: An Exemplar Model-Building Workflow With Applications in Cardiovascular, Metabolic, and Oncology Drug Development.** *CPT Pharmacometrics Syst Pharmacol* 2019. [PMID 31087533](https://pubmed.ncbi.nlm.nih.gov/31087533/)
139. Allen RJ et al. **Efficient Generation and Selection of Virtual Populations in Quantitative Systems Pharmacology Models.** *CPT Pharmacometrics Syst Pharmacol* 2016. — How the Shiny app generates its virtual patient population. [PMID 27069777](https://pubmed.ncbi.nlm.nih.gov/27069777/)
140. Rieger TR et al. **Improving the generation and selection of virtual populations in quantitative systems pharmacology models.** *Prog Biophys Mol Biol* 2018. [PMID 29902482](https://pubmed.ncbi.nlm.nih.gov/29902482/)
141. Bai JPF et al. **FDA-Industry Scientific Exchange on assessing quantitative systems pharmacology models in clinical drug development.** *AAPS J* 2021. [PMID 33931790](https://pubmed.ncbi.nlm.nih.gov/33931790/)
142. Bai JPF et al. **Translational Quantitative Systems Pharmacology in Drug Development: from Current Landscape to Good Practices.** *AAPS J* 2019. [PMID 31161268](https://pubmed.ncbi.nlm.nih.gov/31161268/)
143. Braniff N et al. **An integrated quantitative systems pharmacology virtual population approach for calibration with oncology efficacy data.** *CPT Pharmacometrics Syst Pharmacol* 2025. [PMID 39508122](https://pubmed.ncbi.nlm.nih.gov/39508122/)
144. Traynard P et al. **Logic Modeling in Quantitative Systems Pharmacology.** *CPT Pharmacometrics Syst Pharmacol* 2017. [PMID 28681552](https://pubmed.ncbi.nlm.nih.gov/28681552/)

---

## Summary of the quantitative calibration targets

| # | Target | Source | Value | Use |
|---|------|------|-----|------|
| T1 | Plasma C-triol, patients | PMID 33228797 | 88.31 ng/mL | 🎯 Calibration (`ktri`, `Ktri`) |
| T2 | Plasma C-triol, controls | PMID 33228797 | 5.97 ng/mL | 🎯 Calibration (`ktri`, `Ktri`) |
| T3 | 5-domain NPCCSS progression | PMID 33228797 | 1.5 points/year | 🎯 Calibration (`kdie`) |
| T4 | 17-domain NPCCSS progression | PMID 33228797 | 2.7–2.9 points/year | 🔍 Prediction validation |
| T5 | Arimoclomol, 12-month 5-domain difference | PMID 34418116 | −1.40 points | 🎯 Calibration (`Emax_fold`) |
| T6 | Arimoclomol, 12-month 4-domain difference | PMID 40520915 | −1.70 points | 🔍 Prediction validation |
| T7 | Arimoclomol, progression in the placebo arm | PMID 34418116 | +2.11–2.15 points/year | 🔍 Prediction validation |
| T8 | Levacetylleucine, 12-week SARA difference | PMID 38294974 | −1.28 points | 🎯 Calibration (`nal_Emax_sym`) |
| T9 | C-triol vs 5-domain Spearman ρ | PMID 33228797 | 0.265 | 🔍 Prediction validation |
| T10 | Adrabetadex, open-label progression | PMID 28803710 | 7/14 vs 21/21 | 🔍 Prediction validation |
| T11 | Baseline SARA (IB1001-301) | PMID 38294974 | 15.91 ± 7.65 | 🎯 Calibration (`sara_k`) |
| T12 | Levacetylleucine, 18-month SARA change | PMID 40513057 | −1.64 points | 🔍 Validation — **discrepant** |
| T13 | Levacetylleucine, 12-month SARA change | PMID 40513057 | −1.88 points | 🔍 Validation — **discrepant** |
| T14 | Linearity of progression, independent of age at onset | PMID 19415691 | — | 🔍 Structural validation |
| T15 | Arimoclomol, miglustat subgroup | PMID 34418116 | −2.06 points | 🔍 Prediction validation |

Seven targets were used in calibration (T1·T2·T3·T8·T11, together with the natural-history anchors combined with them) and eight were used for validation only. Of the validation targets, T5·T6·T12·T13·T15 were **not reproduced**, and are reported as they stand without adjustment. The results are
in the "model–data comparison table" of [README_en.md](README_en.md), and **items that did not
match are likewise reported as they stand, without adjustment**.

---

## Disclaimer

This reference list and the model are for educational and research purposes. They must not be used
for diagnostic or treatment decisions in an individual patient. The figures cited are values that can
be verified in the source articles, and the model parameters are approximations tuned to reproduce them.
