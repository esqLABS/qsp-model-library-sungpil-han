# Multiple System Atrophy (MSA) — references

> **PMID verification.** Every entry below was queried programmatically through the NCBI E-utilities
> (`esearch` + `esummary`) to **confirm that the PMID · first author ·
> journal · year · title match the actual PubMed record**.
> Citations that could not be verified were not included in this file. The one exception is
> the verdiperstat phase 3 trial (M-STAR) in §12, which is not a paper indexed in PubMed and is
> given by its clinical trial registration number instead — a fact stated explicitly below.
>
> Link format: `https://pubmed.ncbi.nlm.nih.gov/<PMID>/`

---

## 1. Diagnostic criteria and clinical overview

The clinical endpoints of the MSA model (the §22 cluster) and its phenotype classification (MSA-P / MSA-C)
follow the 2022 MDS criteria. That cognitive impairment is no longer an exclusion criterion,
and that a "prodromal MSA" category has been newly introduced, are the basis for the model's early-intervention scenario
(scenario 19, anti-α-synuclein antibody given early vs late).

1. Wenning GK, et al. **The Movement Disorder Society Criteria for the Diagnosis of Multiple System Atrophy.** *Mov Disord* 2022. [PMID 35445419](https://pubmed.ncbi.nlm.nih.gov/35445419/)
2. Gilman S, et al. **Second consensus statement on the diagnosis of multiple system atrophy.** *Neurology* 2008. [PMID 18725592](https://pubmed.ncbi.nlm.nih.gov/18725592/)
3. Sekiya H, et al. **Validation Study of the MDS Criteria for the Diagnosis of Multiple System Atrophy in the Mayo Clinic Brain Bank.** *Neurology* 2023. [PMID 37816641](https://pubmed.ncbi.nlm.nih.gov/37816641/)
4. Sun Y, et al. **Comparison of the second consensus statement with the movement disorder society criteria for multiple system atrophy.** *Parkinsonism Relat Disord* 2023. [PMID 36529110](https://pubmed.ncbi.nlm.nih.gov/36529110/)
5. Lamotte G, et al. **Movement disorder society criteria for the diagnosis of multiple system atrophy — what's new?** *Clin Auton Res* 2022. [PMID 35633428](https://pubmed.ncbi.nlm.nih.gov/35633428/)
6. Stankovic I, et al. **A Review on the Clinical Diagnosis of Multiple System Atrophy.** *Cerebellum* 2023. [PMID 35986227](https://pubmed.ncbi.nlm.nih.gov/35986227/)
7. Wenning GK, et al. **Multiple system atrophy.** *Handb Clin Neurol* 2013. [PMID 24095129](https://pubmed.ncbi.nlm.nih.gov/24095129/)
8. Peeraully T. **Multiple system atrophy.** *Semin Neurol* 2014. [PMID 24963676](https://pubmed.ncbi.nlm.nih.gov/24963676/)
9. Wenning GK, et al. **Multiple system atrophy: a review of 203 pathologically proven cases.** *Mov Disord* 1997. [PMID 9087971](https://pubmed.ncbi.nlm.nih.gov/9087971/)
10. Lin DJ, et al. **The Diagnosis and Natural History of Multiple System Atrophy, Cerebellar Type.** *Cerebellum* 2016. [PMID 26467153](https://pubmed.ncbi.nlm.nih.gov/26467153/)
11. Liu M, et al. **Multiple system atrophy: an update and emerging directions of biomarkers and clinical trials.** *J Neurol* 2024. [PMID 38483626](https://pubmed.ncbi.nlm.nih.gov/38483626/)

---

## 2. Natural history · progression rate · survival

**Calibration targets.** The `KND` of the mrgsolve model
(the neurodegeneration driving coefficient) and the survival hazard function (`H0`, `B_U2`, `B_STR`, `B_RESP`) were
calibrated to two figures reported in the cohorts below:
**a UMSARS-II increase of 5–8 points per year** and **a median survival of 6–10 years from onset**.
The model's own self-test forces the slope over the 4→8 year interval and the median survival
into the ranges 3–11 points/year and 6–11 years respectively.

12. Wenning GK, et al. **The natural history of multiple system atrophy: a prospective European cohort study.** *Lancet Neurol* 2013. [PMID 23391524](https://pubmed.ncbi.nlm.nih.gov/23391524/)
13. Low PA, et al. **Natural history of multiple system atrophy in the USA: a prospective cohort study.** *Lancet Neurol* 2015. [PMID 26025783](https://pubmed.ncbi.nlm.nih.gov/26025783/)
14. Goldstein DS, et al. **Survival in synucleinopathies: A prospective cohort study.** *Neurology* 2015. [PMID 26432848](https://pubmed.ncbi.nlm.nih.gov/26432848/)
15. Wenning GK, et al. **Development and validation of the Unified Multiple System Atrophy Rating Scale (UMSARS).** *Mov Disord* 2004. [PMID 15452868](https://pubmed.ncbi.nlm.nih.gov/15452868/)
16. Palma JA, et al. **Limitations of the Unified Multiple System Atrophy Rating Scale as outcome measure for clinical trials and a roadmap for improvement.** *Clin Auton Res* 2021. [PMID 33554315](https://pubmed.ncbi.nlm.nih.gov/33554315/)
17. Geser F, et al. **The European Multiple System Atrophy-Study Group (EMSA-SG).** *J Neural Transm (Vienna)* 2005. [PMID 16049636](https://pubmed.ncbi.nlm.nih.gov/16049636/)
18. Xiao Y, et al. **Modified version of unified multiple system atrophy rating scale for remote video-based assessments.** *NPJ Parkinsons Dis* 2023. [PMID 37891215](https://pubmed.ncbi.nlm.nih.gov/37891215/)
19. Kaufmann H, et al. **Multiple System Atrophy Combined Outcome Assessment (MuSyCA): process, format, and validation plan.** *Clin Auton Res* 2026. [PMID 41762390](https://pubmed.ncbi.nlm.nih.gov/41762390/)
20. Feng T, et al. **Natural history and 12-month progression of multiple system atrophy in a Chinese cohort.** *BMC Neurol* 2026. [PMID 42286509](https://pubmed.ncbi.nlm.nih.gov/42286509/)

---

## 3. Oligodendroglial α-synuclein pathology / GCI

**The basis of the model's fourth structural claim.** MSA is a disease in which α-synuclein
aggregates in **oligodendrocytes**, not in neurons. The model's
`ASYNM → ASYNO → GCI` cascade is catalysed solely by `P25A` (the redistribution of p25α/TPPP from myelin
to the cell body), and that `P25A` is in turn generated by a fall in `MYE` (myelin
integrity) — that is, a self-amplifying loop. Neuronal death occurs **downstream**,
through the `WTROPH` term (lost oligodendroglial trophic support), so a treatment that
targets neurons alone is in principle aimed at the wrong compartment.

21. Wakabayashi K, et al. **Alpha-synuclein immunoreactivity in glial cytoplasmic inclusions in multiple system atrophy.** *Neurosci Lett* 1998. [PMID 9682846](https://pubmed.ncbi.nlm.nih.gov/9682846/)
22. Ndayisaba A, et al. **Multiple System Atrophy: Pathology, Pathogenesis, and Path Forward.** *Annu Rev Pathol* 2025. [PMID 39405585](https://pubmed.ncbi.nlm.nih.gov/39405585/)
23. Reddy K, et al. **Multiple system atrophy: α-Synuclein strains at the neuron-oligodendrocyte crossroad.** *Mol Neurodegener* 2022. [PMID 36435784](https://pubmed.ncbi.nlm.nih.gov/36435784/)
24. Hoffmann A, et al. **Oligodendroglial α-synucleinopathy-driven neuroinflammation in multiple system atrophy.** *Brain Pathol* 2019. [PMID 30444295](https://pubmed.ncbi.nlm.nih.gov/30444295/)
25. Wiseman JA, et al. **Neuronal α-synuclein toxicity is the key driver of neurodegeneration in multiple system atrophy.** *Brain* 2025. [PMID 39908177](https://pubmed.ncbi.nlm.nih.gov/39908177/)
26. Mavroeidi P, et al. **Endogenous oligodendroglial alpha-synuclein and TPPP/p25α orchestrate alpha-synuclein pathology in experimental multiple system atrophy models.** *Acta Neuropathol* 2019. [PMID 31011860](https://pubmed.ncbi.nlm.nih.gov/31011860/)
27. Lindersson E, et al. **p25alpha Stimulates alpha-synuclein aggregation and is co-localized with aggregated alpha-synuclein in alpha-synucleinopathies.** *J Biol Chem* 2005. [PMID 15590652](https://pubmed.ncbi.nlm.nih.gov/15590652/)
28. Ferreira N, et al. **Multiple system atrophy-associated oligodendroglial protein p25α stimulates formation of novel α-synuclein strain with enhanced neurodegenerative potential.** *Acta Neuropathol* 2021. [PMID 33978813](https://pubmed.ncbi.nlm.nih.gov/33978813/)
29. Mavroeidi P, et al. **Autophagy mediates the clearance of oligodendroglial SNCA/alpha-synuclein and TPPP/p25A in multiple system atrophy models.** *Autophagy* 2022. [PMID 35000546](https://pubmed.ncbi.nlm.nih.gov/35000546/)
30. Kragh CL, et al. **FAS-dependent cell death in α-synuclein transgenic oligodendrocyte models of multiple system atrophy.** *PLoS One* 2013. [PMID 23372841](https://pubmed.ncbi.nlm.nih.gov/23372841/)

---

## 4. Selective neuronal loss — why the preganglionic side collapses

**The anatomical basis of the model's first structural claim.** What disappears in MSA is
the **preganglionic** sympathetic neurons of the spinal intermediolateral nucleus (IML) and the brainstem cardiovascular neurons, while
the **postganglionic** terminals are relatively preserved. The model expresses this asymmetry with
a single parameter, `VULN_PG = 0.08` (a tenth or less of the other populations), and
separates `G_CENT` (central gain) from `POSTG` (postganglionic integrity) as multiplicative terms.
Raising `VULN_PG` alone in the same equations gives the pure autonomic failure (PAF) phenotype,
and that alone makes the response to atomoxetine disappear (scenarios 7/8).

31. Oppenheimer DR. **Lateral horn cells in progressive autonomic failure.** *J Neurol Sci* 1980. [PMID 6247458](https://pubmed.ncbi.nlm.nih.gov/6247458/)
32. Terao S, et al. **Disease-specific patterns of neuronal loss in the spinal ventral horn in amyotrophic lateral sclerosis, multiple system atrophy and Werdnig-Hoffmann disease.** *J Neurol* 1994. [PMID 8195817](https://pubmed.ncbi.nlm.nih.gov/8195817/)
33. Cortelli P, et al. **Autonomic blood pressure control.** *Handb Clin Neurol* 2026. [PMID 41896018](https://pubmed.ncbi.nlm.nih.gov/41896018/)

---

## 5. Prion-like propagation — seeding and strain biology

The model's `SEED` compartment (extracellular seeding-competent α-synuclein) and
the `SEED → ASYN_UPTAKE` positive feedback rest on the propagation experiments below. This
compartment is the **only** target an anti-α-synuclein antibody (`MAB`) can bind, and
for that reason the antibody's effect in the model is gated by the number of neurons still remaining.

34. Prusiner SB, et al. **Evidence for α-synuclein prions causing multiple system atrophy in humans with parkinsonism.** *Proc Natl Acad Sci U S A* 2015. [PMID 26324905](https://pubmed.ncbi.nlm.nih.gov/26324905/)
35. Holec SAM, et al. **Multiple system atrophy prions transmit neurological disease to mice expressing wild-type human α-synuclein.** *Acta Neuropathol* 2022. [PMID 36018376](https://pubmed.ncbi.nlm.nih.gov/36018376/)
36. Holec SAM, et al. **α-synuclein prion strains differentially adapt after passage in mice.** *PLoS Pathog* 2024. [PMID 39642110](https://pubmed.ncbi.nlm.nih.gov/39642110/)
37. Dhillon JS, et al. **Dissecting α-synuclein inclusion pathology diversity in multiple system atrophy: implications for the prion-like transmission hypothesis.** *Lab Invest* 2019. [PMID 30737468](https://pubmed.ncbi.nlm.nih.gov/30737468/)
38. Jellinger KA, et al. **Is Multiple System Atrophy a Prion-like Disorder?** *Int J Mol Sci* 2021. [PMID 34576255](https://pubmed.ncbi.nlm.nih.gov/34576255/)

---

## 6. Genetics · the coenzyme Q10 axis

The model parameter `COQ2F` (CoQ10 biosynthetic capacity) represents a carrier of a COQ2 variant, and
in the ubiquinol scenario (scenario 17) the difference in response between the wild type and the COQ2-deficient form
is **not coded in** but emerges from the `CQ → OXS → NDRIVE` route.

39. Multiple-System Atrophy Research Collaboration. **Mutations in COQ2 in familial and sporadic multiple-system atrophy.** *N Engl J Med* 2013. [PMID 23758206](https://pubmed.ncbi.nlm.nih.gov/23758206/)
40. Ogaki K, et al. **Analysis of COQ2 gene in multiple system atrophy.** *Mol Neurodegener* 2014. [PMID 25373618](https://pubmed.ncbi.nlm.nih.gov/25373618/)
41. Porto KJ, et al. **COQ2 V393A confers high risk susceptibility for multiple system atrophy in East Asian population.** *J Neurol Sci* 2021. [PMID 34455210](https://pubmed.ncbi.nlm.nih.gov/34455210/)
42. Procopio R, et al. **Genetic mutation analysis of the COQ2 gene in Italian patients with multiple system atrophy.** *Gene* 2019. [PMID 31398377](https://pubmed.ncbi.nlm.nih.gov/31398377/)
43. Scholz SW, et al. **SNCA variants are associated with increased risk for multiple system atrophy.** *Ann Neurol* 2009. [PMID 19475667](https://pubmed.ncbi.nlm.nih.gov/19475667/)
44. Federoff M, et al. **Multiple system atrophy: the application of genetics in understanding etiology.** *Clin Auton Res* 2015. [PMID 25687905](https://pubmed.ncbi.nlm.nih.gov/25687905/)
45. Bougea A, et al. **Genetics of Multiple System Atrophy and Progressive Supranuclear Palsy: A Systemized Review of the Literature.** *Int J Mol Sci* 2023. [PMID 36982356](https://pubmed.ncbi.nlm.nih.gov/36982356/)
46. Li XY, et al. **Genetic profiles of multiple system atrophy revealed by exome sequencing, long-read sequencing and spinocerebellar ataxia repeat expansion analysis.** *Eur J Neurol* 2024. [PMID 39152783](https://pubmed.ncbi.nlm.nih.gov/39152783/)

---

## 7. Neuroinflammation · myeloperoxidase

**The basis for why scenario 18 has to be null.** In MSA transgenic mice,
MPO inhibition was neuroprotective when started early (#47), but with **delayed initiation** it
failed to protect neurons even though it clearly suppressed microglia (#48). The model
reproduces this distinction with `WINDOW` logic — that is, MPO inhibition blocks only the `MGL → MPO → OXS`
route, and by the time of diagnosis the `GCI` and `1-MYE` terms already dominate `NDRIVE`,
so the UMSARS-II curve barely moves.

47. Stefanova N, et al. **Myeloperoxidase inhibition ameliorates multiple system atrophy-like degeneration in a transgenic mouse model.** *Neurotox Res* 2012. [PMID 22161470](https://pubmed.ncbi.nlm.nih.gov/22161470/)
48. Kaindlstorfer C, et al. **Failure of Neuroprotection Despite Microglial Suppression by Delayed-Start Myeloperoxidase Inhibition in a Model of Advanced Multiple System Atrophy.** *Neurotox Res* 2015. [PMID 26194617](https://pubmed.ncbi.nlm.nih.gov/26194617/)

---

## 8. Orthostatic hypotension (nOH) · supine hypertension · nocturnal polyuria — the nocturnal loop

**The basis of the model's third structural claim, and the most important section in this file.**
Goldstein (#49) showed that in autonomic failure **supine hypertension and orthostatic hypotension are
linked and appear together**; Okamoto (#50) showed the loss of nocturnal blood pressure
dipping; and Shibao (#51) showed that supine hypertension causes nocturnal volume loss through
**pressure natriuresis**. Mathias (#54) ran the experiment that closes
the loop — at bedtime, desmopressin reduced nocturnal polyuria and overnight weight loss
and **improved the orthostatic hypotension of the following
morning**.

In the model both `UNAV` and `UVOL` carry an `exp(K*(MAP - MAPNAT))` term, so the blood pressure
that rises overnight while the patient lies down spends in advance the volume needed that morning.
`KPN` (the pressure term) is left **deliberately shallow**, so that a 25 mmHg rise gives about a 1.5–2-fold
nocturnal natriuresis, while long-term volume homeostasis is handled by a separate, steep
`KVOLN` (volume error) term. Without this separation the loop either runs away unrealistically
(a steep pressure term) or total sodium diverges (with the pressure term alone).

49. Goldstein DS, et al. **Association between supine hypertension and orthostatic hypotension in autonomic failure.** *Hypertension* 2003. [PMID 12835329](https://pubmed.ncbi.nlm.nih.gov/12835329/)
50. Okamoto LE, et al. **Nocturnal blood pressure dipping in the hypertension of autonomic failure.** *Hypertension* 2009. [PMID 19047577](https://pubmed.ncbi.nlm.nih.gov/19047577/)
51. Shibao C, et al. **Clonidine for the treatment of supine hypertension and pressure natriuresis in autonomic failure.** *Hypertension* 2006. [PMID 16391172](https://pubmed.ncbi.nlm.nih.gov/16391172/)
52. Fanciulli A, et al. **Consensus statement on the definition of neurogenic supine hypertension in cardiovascular autonomic failure by the American Autonomic Society and the European Federation of Autonomic Societies.** *Clin Auton Res* 2018. [PMID 29766366](https://pubmed.ncbi.nlm.nih.gov/29766366/)
53. Park JW, et al. **Advances in the Pathophysiology and Management of Supine Hypertension in Patients with Neurogenic Orthostatic Hypotension.** *Curr Hypertens Rep* 2022. [PMID 35230654](https://pubmed.ncbi.nlm.nih.gov/35230654/)
54. Mathias CJ, et al. **The effect of desmopressin on nocturnal polyuria, overnight weight loss, and morning postural hypotension in patients with autonomic failure.** *Br Med J (Clin Res Ed)* 1986. [PMID 3089519](https://pubmed.ncbi.nlm.nih.gov/3089519/)
55. Umbertini E, et al. **Understanding nocturnal polyuria in cardiovascular autonomic failure: Pathophysiological mechanisms and clinical implications.** *Auton Neurosci* 2026. [PMID 42241932](https://pubmed.ncbi.nlm.nih.gov/42241932/)
56. Norcliffe-Kaufmann L, et al. **Orthostatic heart rate changes in patients with autonomic failure caused by neurodegenerative synucleinopathies.** *Ann Neurol* 2018. [PMID 29405350](https://pubmed.ncbi.nlm.nih.gov/29405350/)
57. Pavy-Le Traon A, et al. **New insights into orthostatic hypotension in multiple system atrophy: a European multicentre cohort study.** *J Neurol Neurosurg Psychiatry* 2016. [PMID 25977316](https://pubmed.ncbi.nlm.nih.gov/25977316/)
58. Jiang Q, et al. **Orthostatic Hypotension in Multiple System Atrophy: Related Factors and Disease Prognosis.** *J Parkinsons Dis* 2023. [PMID 38143372](https://pubmed.ncbi.nlm.nih.gov/38143372/)
59. Idiaquez JF, et al. **Neurogenic Orthostatic Hypotension. Lessons From Synucleinopathies.** *Am J Hypertens* 2021. [PMID 33705537](https://pubmed.ncbi.nlm.nih.gov/33705537/)

---

## 9. Dissociation of the neurohormonal responses: NE, renin, AVP

There are three parts to the neurohormonal fingerprint of MSA that the model has to reproduce:
**(i) supine plasma NE is normal** (because the postganglionic terminals are alive),
**(ii) the rise in NE on standing is blunted** (because there is no central drive),
**(iii) renin does not rise even on standing** (because β1 sympathetic stimulation of the
juxtaglomerular apparatus has gone). In the model (i)–(iii) follow automatically from the `POSTG`, `G_CENT` and
`KBRENIN*SNA*POSTG` terms respectively. Because the baroreceptor arm of AVP
(`KBAVP*G_CENT*...`) is gated by `G_CENT`, the fact that hypotension can no longer
release AVP comes out of the same structure.

60. Biaggioni I, et al. **Hyporeninemic normoaldosteronism in severe autonomic failure.** *J Clin Endocrinol Metab* 1993. [PMID 7680352](https://pubmed.ncbi.nlm.nih.gov/7680352/)
61. Giza RJ, et al. **Clinical and neurohormonal characteristics in African Americans with neurogenic orthostatic hypotension.** *Clin Auton Res* 2021. [PMID 33502643](https://pubmed.ncbi.nlm.nih.gov/33502643/)
62. Mendoza-Velásquez JJ, et al. **Autonomic Dysfunction in α-Synucleinopathies.** *Front Neurol* 2019. [PMID 31031694](https://pubmed.ncbi.nlm.nih.gov/31031694/)

---

## 10. Cardiac sympathetic imaging and skin biopsy — the fork between MSA and PD/PAF

The model's central assumption that `POSTG` is preserved in MSA and lost in PD/PAF is
**directly measurable clinically**: myocardial ¹²³I-MIBG uptake is
normal in MSA and reduced in PD/PAF. The model expresses this as the `MIBG` node and,
through `POSTG`, ties it to the **same state variable** as the atomoxetine response — that is,
MIBG becomes a predictive biomarker in the model.

63. King AE, et al. **Meta-analysis of 123I-MIBG cardiac scintigraphy for the diagnosis of Lewy body-related disorders.** *Mov Disord* 2011. [PMID 21480373](https://pubmed.ncbi.nlm.nih.gov/21480373/)
64. Alves Do Rego C, et al. **Prospective study of relevance of (123)I-MIBG myocardial scintigraphy and clonidine GH test to distinguish Parkinson's disease and multiple system atrophy.** *J Neurol* 2018. [PMID 29956027](https://pubmed.ncbi.nlm.nih.gov/29956027/)
65. Catalan M, et al. **(123)I-Metaiodobenzylguanidine Myocardial Scintigraphy in Discriminating Degenerative Parkinsonisms.** *Mov Disord Clin Pract* 2021. [PMID 34295947](https://pubmed.ncbi.nlm.nih.gov/34295947/)
66. Yang T, et al. **(131)I-MIBG myocardial scintigraphy for differentiation of Parkinson's disease from multiple system atrophy or essential tremor.** *J Neurol Sci* 2017. [PMID 28131225](https://pubmed.ncbi.nlm.nih.gov/28131225/)
67. Donadio V, et al. **Skin sympathetic fiber α-synuclein deposits: a potential biomarker for pure autonomic failure.** *Neurology* 2013. [PMID 23390175](https://pubmed.ncbi.nlm.nih.gov/23390175/)

---

## 11. The postsynaptic origin of levodopa failure

**The decisive basis of the model's second structural claim.** Churchyard (#68) showed directly that the dopa
resistance of MSA is a **loss of postsynaptic D2 receptors**, and
Sawle (#69) separated the pre- from the postsynaptic change by imaging. The model
writes striatal output as `STRIAT = G_POST × DA/(EC50DA + DA)` and, in MSA,
lets `G_POST = NMSN^GEXP_MSN` fall. Levodopa can raise only `DA`,
so **at the same brain exposure PD is almost normalised while MSA-P responds less and
less** — the two red flags of loss within 1–2 years after an initial partial response, and the absence
of dyskinesia, both come out of this one line.

68. Churchyard A, et al. **Dopa resistance in multiple-system atrophy: loss of postsynaptic D2 receptors.** *Ann Neurol* 1993. [PMID 8338346](https://pubmed.ncbi.nlm.nih.gov/8338346/)
69. Sawle GV, et al. **Asymmetrical pre-synaptic and post-synaptic changes in the striatal dopamine projection in dopa naïve parkinsonism.** *Brain* 1993. [PMID 8353712](https://pubmed.ncbi.nlm.nih.gov/8353712/)
70. Booij J, et al. **The clinical benefit of imaging striatal dopamine transporters with [123I]FP-CIT SPET in differentiating patients with presynaptic parkinsonism from those with other forms of parkinsonism.** *Eur J Nucl Med* 2001. [PMID 11315592](https://pubmed.ncbi.nlm.nih.gov/11315592/)

---

## 12. Pressor pharmacology — selectivity by lesion site

### 12.1 Midodrine — acts downstream of every lesion
In the model the active metabolite desglymidodrine is added **directly**
to the α1 agonist term `AGON`, so it works irrespective of whether the lesion is central or postganglionic. And because
`A1R` (denervation supersensitivity) multiplies it, the clinical observation that **the same mg produces a larger
ΔSBP as the disease progresses** is derived in the model (a self-test item).

71. Low PA, et al. **Efficacy of midodrine vs placebo in neurogenic orthostatic hypotension. A randomized, double-blind multicenter study.** *JAMA* 1997. [PMID 9091692](https://pubmed.ncbi.nlm.nih.gov/9091692/)
72. Jankovic J, et al. **Neurogenic orthostatic hypotension: a double-blind, placebo-controlled study with midodrine.** *Am J Med* 1993. [PMID 7687093](https://pubmed.ncbi.nlm.nih.gov/7687093/)
73. Wright RA, et al. **A double-blind, dose-response study of midodrine in neurogenic orthostatic hypotension.** *Neurology* 1998. [PMID 9674789](https://pubmed.ncbi.nlm.nih.gov/9674789/)
74. Fouad-Tarazi FM, et al. **Alpha sympathomimetic treatment of autonomic insufficiency with orthostatic hypotension.** *Am J Med* 1995. [PMID 7503082](https://pubmed.ncbi.nlm.nih.gov/7503082/)

### 12.2 Droxidopa — **requires** postganglionic AADC
The model: `DROXNE = KAADC × CDRX × POSTG × (1 − CARBI)`. This one line
produces, at the same time, both (a) that droxidopa works only if the terminals are alive and (b) that **the carbidopa
co-administered for the parkinsonism occupies the same enzyme and so the two antagonise each other**
(scenario 6).

75. Kaufmann H, et al. **Droxidopa for neurogenic orthostatic hypotension: a randomized, placebo-controlled, phase 3 trial.** *Neurology* 2014. [PMID 24944260](https://pubmed.ncbi.nlm.nih.gov/24944260/)
76. Biaggioni I, et al. **Randomized withdrawal study of patients with symptomatic neurogenic orthostatic hypotension responsive to droxidopa.** *Hypertension* 2015. [PMID 25350981](https://pubmed.ncbi.nlm.nih.gov/25350981/)
77. Elgebaly A, et al. **Meta-analysis of the safety and efficacy of droxidopa for neurogenic orthostatic hypotension.** *Clin Auton Res* 2016. [PMID 26951135](https://pubmed.ncbi.nlm.nih.gov/26951135/)
78. Strassheim V, et al. **Droxidopa for orthostatic hypotension: a systematic review and meta-analysis.** *J Hypertens* 2016. [PMID 27442791](https://pubmed.ncbi.nlm.nih.gov/27442791/)
79. Isaacson S, et al. **Long-term safety of droxidopa in patients with symptomatic neurogenic orthostatic hypotension.** *J Am Soc Hypertens* 2016. [PMID 27614923](https://pubmed.ncbi.nlm.nih.gov/27614923/)
80. Chen JJ, et al. **Standing and Supine Blood Pressure Outcomes Associated With Droxidopa and Midodrine in Patients With Neurogenic Orthostatic Hypotension.** *Ann Pharmacother* 2018. [PMID 29972032](https://pubmed.ncbi.nlm.nih.gov/29972032/)

### 12.3 NET inhibitors (atomoxetine, ampreloxetine) — experimental evidence of MSA selectivity
NET blockade amplifies NE release that is **already happening**, so it is large with a central lesion +
intact postganglionic terminals (= MSA) and small with a postganglionic lesion (= PAF). In the model
this selectivity is nothing but the structure in which NET blockade multiplies into
`NEREL = SNA × NEVES × POSTG`, and nowhere is it written that "atomoxetine is selective for
MSA" (scenarios 7/8, forced in the self-test).

81. Byun JI, et al. **Efficacy of atomoxetine versus midodrine for neurogenic orthostatic hypotension.** *Ann Clin Transl Neurol* 2020. [PMID 31856425](https://pubmed.ncbi.nlm.nih.gov/31856425/)
82. Mwesigwa N, et al. **Atomoxetine on neurogenic orthostatic hypotension: a randomized, double-blind, placebo-controlled crossover trial.** *Clin Auton Res* 2024. [PMID 39294522](https://pubmed.ncbi.nlm.nih.gov/39294522/)
83. Okamoto LE, et al. **Synergistic effect of norepinephrine transporter blockade and α-2 antagonism on blood pressure in autonomic failure.** *Hypertension* 2012. [PMID 22311903](https://pubmed.ncbi.nlm.nih.gov/22311903/)
84. Lo A, et al. **Pharmacokinetics and pharmacodynamics of ampreloxetine, a novel, selective norepinephrine reuptake inhibitor, in symptomatic neurogenic orthostatic hypotension.** *Clin Auton Res* 2021. [PMID 33782836](https://pubmed.ncbi.nlm.nih.gov/33782836/)
85. Kaufmann H, et al. **Safety and efficacy of ampreloxetine in symptomatic neurogenic orthostatic hypotension: a phase 2 trial.** *Clin Auton Res* 2021. [PMID 34657222](https://pubmed.ncbi.nlm.nih.gov/34657222/)
86. Hoxhaj P, et al. **Ampreloxetine Versus Droxidopa in Neurogenic Orthostatic Hypotension: A Comparative Review.** *Cureus* 2023. [PMID 37303338](https://pubmed.ncbi.nlm.nih.gov/37303338/)

### 12.4 α2 antagonists (yohimbine) — require residual central drive
87. Onrot J, et al. **Oral yohimbine in human autonomic failure.** *Neurology* 1987. [PMID 3808301](https://pubmed.ncbi.nlm.nih.gov/3808301/)
88. Biaggioni I, et al. **Manipulation of norepinephrine metabolism with yohimbine in the treatment of autonomic failure.** *J Clin Pharmacol* 1994. [PMID 8089252](https://pubmed.ncbi.nlm.nih.gov/8089252/)

### 12.5 Pyridostigmine — amplifies only reflex neurotransmission
Ganglionic AChE inhibition amplifies the preganglionic transmission that is **driven by the baroreceptors**, so
it does almost nothing while the patient is lying down (reflex error ≈ 0). The model
expresses this by multiplying `PYRAMP` into the **baroreceptor term of `SNA_ss` only**, and the result is that
the observation "it raises standing blood pressure while worsening supine hypertension relatively less"
is derived (scenario 13, forced in the self-test).

89. Singer W, et al. **Pyridostigmine treatment trial in neurogenic orthostatic hypotension.** *Arch Neurol* 2006. [PMID 16476804](https://pubmed.ncbi.nlm.nih.gov/16476804/)
90. Okamoto LE, et al. **Clinical Correlates of Efficacy of Pyridostigmine in the Treatment of Orthostatic Hypotension.** *Hypertension* 2025. [PMID 39727053](https://pubmed.ncbi.nlm.nih.gov/39727053/)
91. Holder AC, et al. **Pyridostigmine for the Management of Neurogenic Orthostatic Hypotension: A Systemic Review.** *J Geriatr Psychiatry Neurol* 2025. [PMID 39043171](https://pubmed.ncbi.nlm.nih.gov/39043171/)

### 12.6 Volume expansion · non-pharmacological measures
92. Veazie S, et al. **Fludrocortisone for orthostatic hypotension.** *Cochrane Database Syst Rev* 2021. [PMID 34000076](https://pubmed.ncbi.nlm.nih.gov/34000076/)
93. van Lieshout JJ, et al. **Fludrocortisone and sleeping in the head-up position limit the postural decrease in cardiac output in autonomic failure.** *Clin Auton Res* 2000. [PMID 10750642](https://pubmed.ncbi.nlm.nih.gov/10750642/)
94. May M, et al. **The osmopressor response to water drinking.** *Am J Physiol Regul Integr Comp Physiol* 2011. [PMID 21048076](https://pubmed.ncbi.nlm.nih.gov/21048076/)
95. Okamoto LE, et al. **Efficacy of Servo-Controlled Splanchnic Venous Compression in the Treatment of Orthostatic Hypotension: A Randomized Comparison With Midodrine.** *Hypertension* 2016. [PMID 27271310](https://pubmed.ncbi.nlm.nih.gov/27271310/)
96. van der Stam AH, et al. **The Impact of Head-Up Tilt Sleeping on Orthostatic Tolerance: A Scoping Review.** *Biology (Basel)* 2023. [PMID 37626994](https://pubmed.ncbi.nlm.nih.gov/37626994/)
97. van der Stam AH, et al. **Tolerability and efficacy of full-body head-up tilt sleeping in Parkinson's disease and multiple system atrophy.** *NPJ Parkinsons Dis* 2026. [PMID 42143029](https://pubmed.ncbi.nlm.nih.gov/42143029/)
98. van der Stam AH, et al. **Study protocol for the Heads-Up trial: a phase II randomized controlled trial investigating head-up tilt sleeping to alleviate orthostatic intolerance.** *BMC Neurol* 2024. [PMID 38166676](https://pubmed.ncbi.nlm.nih.gov/38166676/)

### 12.7 Post-prandial hypotension and octreotide
99. Armstrong E, et al. **The effects of the somatostatin analogue, octreotide, on postural hypotension, before and after food ingestion, in primary autonomic failure.** *Clin Auton Res* 1991. [PMID 1822761](https://pubmed.ncbi.nlm.nih.gov/1822761/)
100. Alam M, et al. **Effects of the peptide release inhibitor, octreotide, on daytime hypotension and on nocturnal hypertension in primary autonomic failure.** *J Hypertens* 1995. [PMID 8903629](https://pubmed.ncbi.nlm.nih.gov/8903629/)
101. Smith GD, et al. **Effect of the somatostatin analogue, octreotide, on exercise-induced hypotension in human subjects with chronic sympathetic failure.** *Clin Sci (Lond)* 1995. [PMID 7493436](https://pubmed.ncbi.nlm.nih.gov/7493436/)
102. Jansen RW, et al. **Postprandial hypotension: epidemiology, pathophysiology, and clinical management.** *Ann Intern Med* 1995. [PMID 7825766](https://pubmed.ncbi.nlm.nih.gov/7825766/)
103. Chaudhuri KR, et al. **Alcohol ingestion lowers supine blood pressure, causes splanchnic vasodilatation and worsens postural hypotension in primary autonomic failure.** *J Neurol* 1994. [PMID 8164016](https://pubmed.ncbi.nlm.nih.gov/8164016/)

### 12.8 Treatment guidelines · comprehensive reviews
104. Park JW, et al. **Pharmacologic treatment of orthostatic hypotension.** *Auton Neurosci* 2020. [PMID 32979782](https://pubmed.ncbi.nlm.nih.gov/32979782/)
105. Eschlböck S, et al. **Evidence-based treatment of neurogenic orthostatic hypotension and related symptoms.** *J Neural Transm (Vienna)* 2017. [PMID 29058089](https://pubmed.ncbi.nlm.nih.gov/29058089/)
106. Palma JA, et al. **Management of Orthostatic Hypotension.** *Continuum (Minneap Minn)* 2020. [PMID 31996627](https://pubmed.ncbi.nlm.nih.gov/31996627/)
107. Shibao C, et al. **Pharmacotherapy of autonomic failure.** *Pharmacol Ther* 2012. [PMID 21664375](https://pubmed.ncbi.nlm.nih.gov/21664375/)
108. Chen B, et al. **Non-pharmacological and drug treatment of autonomic dysfunction in multiple system atrophy: current status and future directions.** *J Neurol* 2023. [PMID 37477834](https://pubmed.ncbi.nlm.nih.gov/37477834/)
109. Arbique D, et al. **Management of neurogenic orthostatic hypotension.** *J Am Med Dir Assoc* 2014. [PMID 24388946](https://pubmed.ncbi.nlm.nih.gov/24388946/)
110. Vidal-Petiot E, et al. **Orthostatic hypotension: Review and expert position statement.** *Rev Neurol (Paris)* 2024. [PMID 38123372](https://pubmed.ncbi.nlm.nih.gov/38123372/)
111. Fanciulli A, et al. **Management of Orthostatic Hypotension in Parkinson's Disease.** *J Parkinsons Dis* 2020. [PMID 32716319](https://pubmed.ncbi.nlm.nih.gov/32716319/)

---

## 13. Attempts at disease-modifying treatment — mostly negative

**⚠️ A caution on the verdiperstat phase 3 trial (M-STAR).** The failure of the primary endpoint of this trial
(NCT03952806, Biohaven) **could not be confirmed against a paper indexed in PubMed.**
This file therefore does not present that result as an academic citation but gives only the registration
number, and it places the literature basis of the model's scenario 18 (verdiperstat null) in the
preclinical delayed-initiation experiment of §7 (#48) and the rifampicin failure below (#112).
Scenario 18 should be read not as "a reproduction of the phase 3 result" but as the structural
prediction that **"single inhibition of a downstream inflammatory target, begun at the time of
diagnosis, has to be null in this model structure"**.

112. Low PA, et al. **Efficacy and safety of rifampicin for multiple system atrophy: a randomised, double-blind, placebo-controlled trial.** *Lancet Neurol* 2014. [PMID 24507091](https://pubmed.ncbi.nlm.nih.gov/24507091/)
113. Singer W, et al. **Optimizing clinical trial design for multiple system atrophy: lessons from the rifampicin study.** *Clin Auton Res* 2015. [PMID 25763826](https://pubmed.ncbi.nlm.nih.gov/25763826/)
114. Mitsui J, et al. **High-dose ubiquinol supplementation in multiple-system atrophy: a multicentre, randomised, double-blinded, placebo-controlled phase 2 trial.** *EClinicalMedicine* 2023. [PMID 37256098](https://pubmed.ncbi.nlm.nih.gov/37256098/)
115. Mitsui J, et al. **Three-Year Follow-Up of High-Dose Ubiquinol Supplementation in a Case of Familial Multiple System Atrophy with Compound Heterozygous COQ2 Mutations.** *Cerebellum* 2017. [PMID 28150130](https://pubmed.ncbi.nlm.nih.gov/28150130/)
116. Singer W, et al. **Intrathecal administration of autologous mesenchymal stem cells in multiple system atrophy.** *Neurology* 2019. [PMID 31152011](https://pubmed.ncbi.nlm.nih.gov/31152011/)
117. Poewe W, et al. **Therapeutic advances in multiple system atrophy and progressive supranuclear palsy.** *Mov Disord* 2015. [PMID 26227071](https://pubmed.ncbi.nlm.nih.gov/26227071/)
118. Eschlböck S, et al. **Interventional trials in atypical parkinsonism.** *Parkinsonism Relat Disord* 2016. [PMID 26421389](https://pubmed.ncbi.nlm.nih.gov/26421389/)
119. Jeong SH, et al. **Drug repurposing for disease-modifying effects in multiple system atrophy.** *Transl Neurodegener* 2026. [PMID 42010648](https://pubmed.ncbi.nlm.nih.gov/42010648/)

---

## 14. Biomarkers: NfL, seed amplification, MRI

The model's `NFL` compartment is released in proportion to the **rate** of neuronal loss and is lost
first-order, so in the natural history it peaks in mid-course (about 6 years) and comes down slowly
as the disease exhausts itself. This fits the observation that NfL is **an index of the rate of progression
rather than of cross-sectional severity** (#120–#123), and at the same time keeps the model from
overestimating NfL late in the disease.

120. Chelban V, et al. **Neurofilament light levels predict clinical progression and death in multiple system atrophy.** *Brain* 2022. [PMID 35903017](https://pubmed.ncbi.nlm.nih.gov/35903017/)
121. Zhang L, et al. **Neurofilament Light Chain Predicts Disease Severity and Progression in Multiple System Atrophy.** *Mov Disord* 2022. [PMID 34719813](https://pubmed.ncbi.nlm.nih.gov/34719813/)
122. Singer W, et al. **Neurofilament light chain in spinal fluid and plasma in multiple system atrophy: a prospective, longitudinal biomarker study.** *Clin Auton Res* 2023. [PMID 37603107](https://pubmed.ncbi.nlm.nih.gov/37603107/)
123. Huang J, et al. **Plasma Neurofilament Light Chain as a Biomarker for Motor Progression and Disease Milestones in Multiple System Atrophy: An Update.** *Mov Disord* 2026. [PMID 41912379](https://pubmed.ncbi.nlm.nih.gov/41912379/)
124. Jin B, et al. **Plasma neurofilament light chain as a predictor of multiple system atrophy in idiopathic REM sleep behavior disorder.** *J Neurol* 2025. [PMID 41359208](https://pubmed.ncbi.nlm.nih.gov/41359208/)
125. Fernandes Gomes B, et al. **α-Synuclein seed amplification assay as a diagnostic tool for parkinsonian disorders.** *Parkinsonism Relat Disord* 2023. [PMID 37591709](https://pubmed.ncbi.nlm.nih.gov/37591709/)
126. Rossi M, et al. **Comparison of Two α-Synuclein Seed Amplification Assays for Discrimination of Parkinson Disease and Atypical Parkinsonism.** *Mov Disord* 2025. [PMID 40879244](https://pubmed.ncbi.nlm.nih.gov/40879244/)
127. Grossauer A, et al. **α-Synuclein Seed Amplification Assays in the Diagnosis of Synucleinopathies Using Cerebrospinal Fluid — A Systematic Review and Meta-Analysis.** *Mov Disord Clin Pract* 2023. [PMID 37205253](https://pubmed.ncbi.nlm.nih.gov/37205253/)
128. Sugiyama A, et al. **Revisiting 'hot cross bun' sign: a multicentre MRI study of 97 patients with autopsy-confirmed multiple system atrophy.** *J Neurol Neurosurg Psychiatry* 2026. [PMID 41083252](https://pubmed.ncbi.nlm.nih.gov/41083252/)
129. Portet M, et al. **Hot cross bun sign.** *J Neurol* 2019. [PMID 31254063](https://pubmed.ncbi.nlm.nih.gov/31254063/)

---

## 15. Sleep · stridor · sudden death

The model's survival hazard function gives independent weights to `B_STR` (stridor) and
`B_RESP` (loss of respiratory neurons). This reflects the observation that stridor predicts
survival **separately** from the UMSARS motor score and that tracheostomy prolongs survival
(#131–#133).

130. Silber MH, et al. **Stridor and death in multiple system atrophy.** *Mov Disord* 2000. [PMID 10928581](https://pubmed.ncbi.nlm.nih.gov/10928581/)
131. Cortelli P, et al. **Stridor in multiple system atrophy: Consensus statement on diagnosis, prognosis, and treatment.** *Neurology* 2019. [PMID 31570638](https://pubmed.ncbi.nlm.nih.gov/31570638/)
132. Giannini G, et al. **Early stridor onset and stridor treatment predict survival in 136 patients with MSA.** *Neurology* 2016. [PMID 27566741](https://pubmed.ncbi.nlm.nih.gov/27566741/)
133. Giannini G, et al. **Tracheostomy is associated with increased survival in Multiple System Atrophy patients with stridor.** *Eur J Neurol* 2022. [PMID 35384153](https://pubmed.ncbi.nlm.nih.gov/35384153/)
134. Laga A, et al. **A strategic approach of the management of sleep-disordered breathing in multiple system atrophy.** *J Clin Sleep Med* 2025. [PMID 39539061](https://pubmed.ncbi.nlm.nih.gov/39539061/)
135. Giannini G, et al. **REM Sleep Behaviour Disorder in Multiple System Atrophy: From Prodromal to Progression of Disease.** *Front Neurol* 2021. [PMID 34194385](https://pubmed.ncbi.nlm.nih.gov/34194385/)
136. Postuma RB, et al. **Evolution of Prodromal Multiple System Atrophy from REM Sleep Behavior Disorder: A Descriptive Study.** *J Parkinsons Dis* 2022. [PMID 35094998](https://pubmed.ncbi.nlm.nih.gov/35094998/)
137. Postuma RB, et al. **Risk and predictors of dementia and parkinsonism in idiopathic REM sleep behaviour disorder: a multicentre study.** *Brain* 2019. [PMID 30789229](https://pubmed.ncbi.nlm.nih.gov/30789229/)

---

## 16. Urogenital autonomic failure — Onuf's nucleus

In the model `NONUF` (Onuf's nucleus) is set as the **most vulnerable** population (`VULN_ONUF = 1.05`),
so post-void residual volume (PVR) and erectile dysfunction appear before the motor symptoms. This reflects
the prospective observation that bladder dysfunction can be the **first symptom** of MSA (#139),
and the report that erectile dysfunction is the earliest symptom in men (#138).

138. Papatsoris AG, et al. **Urinary and erectile dysfunction in multiple system atrophy (MSA).** *Neurourol Urodyn* 2008. [PMID 17563111](https://pubmed.ncbi.nlm.nih.gov/17563111/)
139. Sakakibara R, et al. **Bladder dysfunction as the initial presentation of multiple system atrophy: a prospective cohort study.** *Clin Auton Res* 2019. [PMID 30043182](https://pubmed.ncbi.nlm.nih.gov/30043182/)

---

## 17. Modelling tools

- **mrgsolve** — ODE-based PK/PD and QSP simulation: <https://mrgsolve.org/>
- **Graphviz** — rendering the mechanistic map: <https://graphviz.org/>
- **Shiny** — interactive dashboard: <https://shiny.posit.co/>
- gPKPDviz (a Shiny simulation tool built on mrgsolve), paper:
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC10941578/>

---

## ⚠️ Disclaimer

This reference list and the QSP model linked to it were written for **educational and research purposes**.
The parameters are **illustrative approximations** calibrated to the ranges of the values reported in the
literature above; they have not been fitted to or validated against individual patient data.
**They must not be used directly for actual clinical decision-making, prescribing or regulatory
submission.**
